#include "QuaC_Accelerator.hpp"
#include "QuaC_Circuit_Visitor.hpp"
namespace QuaC {
    void QuaC_Accelerator::initialize(const HeterogeneousMap& params)  
    {
       // DEBUG
        std::cout << ">> DEBUG: QuaC_Accelerator initialize ...\n";
        m_visitor = std::make_shared<CircuitVisitor>();      
        
        // TODO
        // Initialize noise, gate time params etc. since QuaC is able to simulate noisy circuits
    }

    void QuaC_Accelerator::execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::shared_ptr<CompositeInstruction> compositeInstruction)  
    {
        m_visitor->initialize(buffer);
        // Walk the IR tree, and visit each node
        InstructionIterator it(compositeInstruction);
        while (it.hasNext()) 
        {
            auto nextInst = it.next();
            if (nextInst->isEnabled()) 
            {
                nextInst->accept(m_visitor);
            }
        }

        m_visitor->finalize();
    }
    void QuaC_Accelerator::execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::vector<std::shared_ptr<CompositeInstruction>> compositeInstructions)  
    {
       // TODO
    }
}

