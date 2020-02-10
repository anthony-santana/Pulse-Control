#pragma once
#include <memory>
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