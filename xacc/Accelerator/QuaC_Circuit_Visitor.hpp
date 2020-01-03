#pragma one

#include "Identifiable.hpp"
#include "AllGateVisitor.hpp"
#include "AcceleratorBuffer.hpp"
#include "OptionsProvider.hpp"

using namespace xacc;
using namespace xacc::quantum;

namespace QuaC {
    // Circuit (gate-by-gate) visitor:
    // NOTE: we will create a `Pulse`-style visitor backend based on QuaC.
    // For initial set-up, we have a gate-based backend.
    class CircuitVisitor : public AllGateVisitor, public OptionsProvider, public Cloneable<CircuitVisitor>
    {
    public:
        void initialize(std::shared_ptr<AcceleratorBuffer> buffer);
        void finalize();

        void visit(Hadamard& h) override;
        void visit(CNOT& cnot) override;
        void visit(Rz& rz) override;
        void visit(Ry& ry) override;
        void visit(Rx& rx) override;
        void visit(X& x) override;
        void visit(Y& y) override;
        void visit(Z& z) override;
        void visit(CY& cy) override;
        void visit(CZ& cz) override;
        void visit(Swap& s) override;
        void visit(CRZ& crz) override;
        void visit(CH& ch) override;
        void visit(S& s) override;
        void visit(Sdg& sdg) override;
        void visit(T& t) override;
        void visit(Tdg& tdg) override;
        void visit(CPhase& cphase) override;
        void visit(Measure& measure) override;
        void visit(Identity& i) override;
        void visit(U& u) override;

        virtual std::shared_ptr<CircuitVisitor> clone() { return std::make_shared<CircuitVisitor>(); }
    };    
}