#include "QuaC_Circuit_Visitor.hpp"
#include "xacc.hpp"

extern "C" {
#include "interface_xacc_ir.h"
}

namespace QuaC {
    void CircuitVisitor::initialize(std::shared_ptr<AcceleratorBuffer> buffer) 
    {
        XACC_QuaC_Initialize(buffer->size());
    }

    void CircuitVisitor::finalize() 
    {
        const auto result = XACC_QuaC_ExecuteCircuit(0, nullptr);
        XACC_QuaC_Finalize();
    }

    void CircuitVisitor::visit(Hadamard& h)  
    {
        int qbitOperands[1] =  { h.bits()[0] };
        XACC_QuaC_AddInstruction(h.name().c_str(), qbitOperands, 1, 0, nullptr);
    }
    
    void CircuitVisitor::visit(CNOT& cnot)  
    {
        int qbitOperands[2] =  { cnot.bits()[0], cnot.bits()[1] };
        XACC_QuaC_AddInstruction(cnot.name().c_str(), qbitOperands, 2, 0, nullptr);
    }
    
    void CircuitVisitor::visit(Rz& rz)  
    {
        // TODO
    }
    
    void CircuitVisitor::visit(Ry& ry)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(Rx& rx)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(X& x)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(Y& y)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(Z& z)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(CY& cy)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(CZ& cz)  
    {
      // TODO
    }
    
    void CircuitVisitor::visit(Swap& s)  
    {
      // TODO
    }
    
    void CircuitVisitor::visit(CRZ& crz)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(CH& ch)  
    {
       // TODO
    }
    
    void CircuitVisitor::visit(S& s)  
    {
        // TODO
    }
    
    void CircuitVisitor::visit(Sdg& sdg)  
    {
      // TODO
    }
    
    void CircuitVisitor::visit(T& t)  
    {
        // TODO
    }
    
    void CircuitVisitor::visit(Tdg& tdg)  
    {
        // TODO
    }
    
    void CircuitVisitor::visit(CPhase& cphase)  
    {
       // TODO
    }
   
    void CircuitVisitor::visit(Identity& i)  
    {
        // TODO
    }
    
    void CircuitVisitor::visit(U& u)  
    {
       // TODO
    }

    void CircuitVisitor::visit(Measure& measure)  
    {
       // TODO
    }
}