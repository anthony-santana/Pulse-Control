#include "Executor.hpp"
#include "Functor.hpp"

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