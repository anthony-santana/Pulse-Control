
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

// ============================================================================
// Simulate gaussian pulse to implement an X gate
// that was optimized by GOAT.
// ============================================================================

namespace {
    // Note: this is Hamiltonian is in the rotating frame
    const std::string singleQubitHamiltonianJson = R"#(
        {
            "description": "One-qubit Hamiltonian (used by GOAT)",
            "h_latex": "",
            "h_str": ["omega0*Z0", "omegaa*X0||D0"],
            "osc": {},
            "qub": {
                "0": 2
            },
            "vars": {
                "omega0": 0.0,
                "omegaa": 0.062832
            }
        }
    )#";
}

int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    if (!systemModel->loadHamiltonianJson(singleQubitHamiltonianJson))
    {
        return -1;
    }
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 0.5;
    // Rotating frame
    // LO freq = Cavity Freq = 0.0
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    const size_t nbSamples = 200;
    const double gaussSigma = 19.946950;
    channelConfigs.addOrReplacePulse("gaussian", QuaC::GaussianPulse(nbSamples, gaussSigma, channelConfigs.dt));    
    systemModel->setChannelConfigs(channelConfigs);
    // Get QuaC accelerator
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

    // Create a simple pulse program using IR
    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("gaussian", "d0");
    // Add those pulse IR's to the composite and execute.
    compositeInst->addInstruction(pulseInst);
  
    auto qubitReg = xacc::qalloc(1);    
    quaC->execute(qubitReg, compositeInst);
    qubitReg->print();
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}