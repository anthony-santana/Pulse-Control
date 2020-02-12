#include "PlatformNode.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include "Functor.hpp"

#define SHUT_DOWN_CODE -1

namespace MPI {
     PlatformNode::PlatformNode(MPI_Comm in_host, int in_rank):
        m_host(in_host),
        m_shutdownStarted(false),
        m_rank(in_rank)
    {
		// printf("PlatformNode rank %d created.\n", m_rank);
	}

    void PlatformNode::start()
    {
        while (true)
        {
			int messageSize;
			MPI_Bcast(&messageSize, 1, MPI_INT, 0, m_host);		
            // Use a special code to indicate that the host want us to shut down.
            m_shutdownStarted = (messageSize == SHUT_DOWN_CODE);
			if (m_shutdownStarted)
			{
				break;
			}
            
            // Otherwise, receive the buffer
			char buffer[messageSize];
			MPI_Bcast(&buffer, messageSize, MPI_CHAR, 0, m_host);
            // printf("PlatformNode rank %d received a buffer of length %d.\n", m_rank, messageSize);
            std::unique_ptr<FunctorBase> receivedFunctor;
            deserializeBuffer(buffer, messageSize, receivedFunctor);
            SerializationType serializedResult;
            processFunctor(receivedFunctor, serializedResult);
            MPI_Barrier(m_host);
            
            // Got some results from the functor execution, return to host
            if (!serializedResult.str().empty() && m_rank == 0)
            { 
                auto strBuffer = serializedResult.str();                
                auto bufferSize = strBuffer.size() + 1;
                // printf("PlatformNode rank %d is sending a buffer of length %ld to host.\n", m_rank, bufferSize);
                MPI_Send(&bufferSize, 1, MPI_INT, 0, 0, m_host);
                MPI_Send(&strBuffer[0], bufferSize, MPI_CHAR, 0, 0, m_host);
            }
        }
        
        // Exit the loop
        if (m_rank == 0)
        {
            puts("Shut-down Platform Node!");
        }
    }

    void PlatformNode::processFunctor(std::unique_ptr<FunctorBase>& in_functor, SerializationType& out_result)
    {
        // printf("PlatformNode rank %d processes a '%s' functor.\n", m_rank, in_functor->name().c_str());
        in_functor->execute(&out_result);
    }
}

int main(int argc, char *argv[]) 
{ 
	int returnCode = 0;
  	MPI_Comm parent;
    struct utsname uts;
    int size, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_get_parent(&parent);
	
    const auto finish = [&returnCode](){
        MPI_Finalize(); 
	    exit(returnCode);
    };

    const auto isRank0 = [&rank](){
        return rank == 0;
    };
    
    // This is a slave process, hence must be started via MPI_Comm_spawn,
    // i.e. MPI_Comm_get_parent must return a valid parent.
    if (parent == MPI_COMM_NULL) 
    {
		puts("FATAL: No parent process! This process must be launched by using MPI_Comm_spawn.");
		returnCode = 1;
		finish();
	}

	MPI_Comm_remote_size(parent, &size);
	
    if (size != 1) 
    {
		puts("FATAL: Something's wrong with the parent.");
		returnCode = 2;
		finish();
	}

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
    if (isRank0())
    {
       	uname(&uts);
        printf("Successfully start MPI platform node on '%s' with %d processes.\n", uts.nodename, size);
    }

    MPI::PlatformNode nodeService(parent, rank);
	// Start the service
    nodeService.start();

    finish();
}  