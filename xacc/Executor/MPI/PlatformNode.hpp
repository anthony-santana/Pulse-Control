#pragma once

#include <mpi.h>
#include <memory>
#include "Serialization.hpp"

class FunctorBase;

namespace MPI {
class PlatformNode
{
public:
    // Constructor
    PlatformNode(MPI_Comm in_host, int in_rank);

    void start();
private:
    void processFunctor(std::unique_ptr<FunctorBase>& in_functor, SerializationType& out_result);

private:
    MPI_Comm m_host;
    bool m_shutdownStarted;
    int m_rank;
};
}

