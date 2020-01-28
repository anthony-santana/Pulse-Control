#pragma once

#include "Identifiable.hpp"
#include "AllGateVisitor.hpp"
#include "xacc.hpp"

using namespace xacc;
using namespace xacc::quantum;

namespace QuaC {
    class PulseVisitor;
    class QuaC_Accelerator : public Accelerator {

    public:
        // Identifiable interface impls
        virtual const std::string name() const override { return "QuaC"; }
        virtual const std::string description() const override { return "XACC Simulation Accelerator based on QuaC library."; }
        
        // Accelerator interface impls
        virtual void initialize(const HeterogeneousMap& params = {}) override;
        virtual void updateConfiguration(const HeterogeneousMap& config) override {};

        virtual void contributeInstructions(const std::string& custom_json_config) override;

        virtual const std::vector<std::string> configurationKeys() override { return {}; }
        virtual void execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::shared_ptr<CompositeInstruction> compositeInstruction) override;
        virtual void execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::vector<std::shared_ptr<CompositeInstruction>> compositeInstructions) override;
    private:    
        HeterogeneousMap m_params;
        // TODO: unify the Pulse backend API
        std::shared_ptr<PulseVisitor> m_pulseVisitor;
        std::unordered_map<std::string, std::vector<std::complex<double>>> m_importedPulses;
    };
}
