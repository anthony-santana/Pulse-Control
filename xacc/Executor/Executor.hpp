#pragma once
#include <memory>
#include <mpi.h>
#include "Serialization.hpp"

class FunctorBase;

class FunctorExecutorBase
{
public:
    virtual void PostFunctorAsync(std::unique_ptr<FunctorBase>&& in_functor) = 0;
    virtual void CallFunctorSync(std::unique_ptr<FunctorBase>&& in_functor, SerializationType& out_Result) = 0;
};

// Non-MPI executor
class SingleProcessFunctorExecutor: public FunctorExecutorBase
{
public:
    virtual void PostFunctorAsync(std::unique_ptr<FunctorBase>&& in_functor) override;
    virtual void CallFunctorSync(std::unique_ptr<FunctorBase>&& in_functor, SerializationType& outResult) override;
};

// MPI Process Spawning executor:
// This will spawn MPI processes to execute QuaC functors.
// Intended use case: laptop/workstation where we have complete control of the resources.
// This can also be used on clusters provided that they have supports for MPI_Comm_Spawn.
class CommSpawnFunctorExecutor: public FunctorExecutorBase
{
public:
    CommSpawnFunctorExecutor(int in_nbProcs);
    virtual void PostFunctorAsync(std::unique_ptr<FunctorBase>&& in_functor) override;
    virtual void CallFunctorSync(std::unique_ptr<FunctorBase>&& in_functor, SerializationType& outResult) override;
    ~CommSpawnFunctorExecutor();
private:
    MPI_Comm m_intercomm;
};