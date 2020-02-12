#pragma once

#include "TearDown.hpp"

// Function call to finalize/clean-up the executor.
extern void finalizeQuacExecutor();

// Impl QuaC TearDown interface to shut-down QuaC Executor when we are done (XACC::Finalize)
class QuacTearDown : public xacc::TearDown  
{
public:
    virtual void tearDown() override 
    {
        finalizeQuacExecutor();
    } 
};
