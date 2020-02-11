#include "Executor.hpp"
#include "Functor.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
namespace {
    #define STRINGTIFY(x) _MAKE_STR(x)
    #define _MAKE_STR(x) #x
    #define SHUT_DOWN_CODE -1
}

void SingleProcessFunctorExecutor::PostFunctorAsync(std::unique_ptr<FunctorBase>&& in_functor) 
{
    // Note: for single process executor, there is no differences
    // b/w sync and async (everything is synchronous).
    in_functor->execute();
}


void SingleProcessFunctorExecutor::CallFunctorSync(std::unique_ptr<FunctorBase>&& in_functor, SerializationType& out_Result) 
{
    in_functor->execute(&out_Result);
}


CommSpawnFunctorExecutor::CommSpawnFunctorExecutor(int in_nbProcs)
{
    int world_size;
    int universe_size = in_nbProcs + 1;
    int rank;
   
    // Singleton init
    MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
    if (world_size != 1) 
    {
		printf("FAILURE: Started %d executor processes. Please only start 1 executor process.\n", world_size);
		exit(1);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) 
    {
		printf("FAILURE: Invalid Rank (%d) for the Executor. Expect 0.\n", rank);
		exit(1);
	}

    std::string exeDir;
    // We should get this defined from CMake
#ifdef NODE_SERVICE_EXE_DIR
    exeDir = STRINGTIFY(NODE_SERVICE_EXE_DIR);
#endif
    // If none, use current directory.
    // The MPI_Comm_spawn will fail eventually if it cannot locate the executable.
    if (exeDir.empty())
    {
        exeDir = ".";
    }
    
    const auto exeProcess = exeDir + "/PlatformNode";
	
    // Start Platform Node process
	const auto rc = MPI_Comm_spawn(exeProcess.c_str(), MPI_ARGV_NULL, universe_size - 1,  
			    MPI_INFO_NULL, 0, MPI_COMM_SELF, &m_intercomm,  
			    MPI_ERRCODES_IGNORE);
	if (rc != MPI_SUCCESS) 
    {
		printf("FAILURE: MPI_Comm_spawn(): %d\n", rc);
		exit(1);
	}

    puts("Successfully launch Node service!");
}

void CommSpawnFunctorExecutor::PostFunctorAsync(std::unique_ptr<FunctorBase>&& in_functor) 
{
    auto strBuffer = serializeObject(in_functor);
    std::cout << " Serialize functor '" << in_functor->name() << "' to " << strBuffer.size() << " bytes \n";
    auto bufferSize = strBuffer.size() + 1;
	MPI_Bcast(&bufferSize, 1, MPI_INT, MPI_ROOT, m_intercomm);
	MPI_Bcast(&strBuffer[0], bufferSize, MPI_CHAR, MPI_ROOT, m_intercomm);
	MPI_Barrier(m_intercomm);
}

void CommSpawnFunctorExecutor::CallFunctorSync(std::unique_ptr<FunctorBase>&& in_functor, SerializationType& outResult) 
{
    // TODO: we need to get the result from this one
    auto strBuffer = serializeObject(in_functor);
    std::cout << " Serialize functor '" << in_functor->name() << "' to " << strBuffer.size() << " bytes \n";
    auto bufferSize = strBuffer.size() + 1;
	MPI_Bcast(&bufferSize, 1, MPI_INT, MPI_ROOT, m_intercomm);
	MPI_Bcast(&strBuffer[0], bufferSize, MPI_CHAR, MPI_ROOT, m_intercomm);
	MPI_Barrier(m_intercomm);
}

CommSpawnFunctorExecutor::~CommSpawnFunctorExecutor()
{
    // Send the shut down code to the child processes
    int token = SHUT_DOWN_CODE;
	MPI_Bcast(&token, 1, MPI_INT, MPI_ROOT, m_intercomm);
	puts("CommSpawnFunctorExecutor Exit!");
	MPI_Finalize();
}
