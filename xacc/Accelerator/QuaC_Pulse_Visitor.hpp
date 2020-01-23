#pragma one

#include "Identifiable.hpp"
#include "OptionsProvider.hpp"
#include "AcceleratorBuffer.hpp"
#include "Cloneable.hpp"
#include "PulseChannelController.hpp"
#include "AllGateVisitor.hpp"

using namespace xacc;

namespace QuaC {    
    class PulseVisitor : public OptionsProvider, public Cloneable<PulseVisitor>
    {
    public:
        void initialize(std::shared_ptr<AcceleratorBuffer> buffer, const HeterogeneousMap& in_params = {}, const PulseLib& in_importedPulses = {});
        void solve(const std::shared_ptr<CompositeInstruction>& in_pulseInstruction);
        void finalize();
        virtual std::shared_ptr<PulseVisitor> clone() { return std::make_shared<PulseVisitor>(); }
        static std::vector<std::complex<double>> PulseSamplesToComplexVec(const std::vector<std::vector<double>>& in_samples);
    private:
        void schedulePulses(const std::shared_ptr<CompositeInstruction>& in_pulseInstruction);
    private:
        std::unique_ptr<PulseChannelController> m_pulseChannelController;
    };    
}