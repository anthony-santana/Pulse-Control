#pragma one

#include "Identifiable.hpp"
#include "OptionsProvider.hpp"
#include "AcceleratorBuffer.hpp"
#include "Cloneable.hpp"
#include "PulseChannelController.hpp"
#include "AllGateVisitor.hpp"
#include "Hamiltonian.hpp"
#include "Pulse.hpp"

using namespace xacc;
using namespace xacc::quantum;

namespace QuaC {    
    class PulseSystemModel;
    class PulseVisitor : public AllGateVisitor, public InstructionVisitor<Pulse>, public OptionsProvider, public Cloneable<PulseVisitor>, public IChannelNameResolver
    {
    public:
        void initialize(std::shared_ptr<AcceleratorBuffer> buffer, PulseSystemModel* in_systemModel, const HeterogeneousMap& in_params = {});
        void solve();
        void finalize();
        virtual std::shared_ptr<PulseVisitor> clone() { return std::make_shared<PulseVisitor>(); }
        virtual int GetChannelId(const std::string& in_channelName) override;
        // Gate visit
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
        void visit(Pulse& p) override;

    private:
        void schedulePulses(const std::shared_ptr<CompositeInstruction>& in_pulseInstruction);
    private:
        std::unique_ptr<PulseChannelController> m_pulseChannelController;
        std::shared_ptr<CompositeInstruction> m_pulseComposite;
        PulseSystemModel* m_systemModel; 
        std::shared_ptr<AcceleratorBuffer> m_buffer;
    };    
}