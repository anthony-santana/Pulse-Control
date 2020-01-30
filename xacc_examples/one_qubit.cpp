
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    auto systemModel = new QuaC::PulseSystemModel();

    // An example one-qubit Hamiltonian
    const std::string hamiltonianJson = R"(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["-0.5*omega0*Z0", "0.5*omegaa*X0||D0"],
            "osc": {},
            "qub": {
                "0": 2
            },
            "vars": {
                "omega0": 6.2831853,
                "omegaa": 0.031459
            }
        }
    )";

    if (!systemModel->loadHamiltonianJson(hamiltonianJson))
    {
        return -1;
    }
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;
    channelConfigs.loFregs_dChannels.emplace_back(1.0);
    
    const size_t nbSamples = 100;
    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));
    systemModel->setChannelConfigs(channelConfigs);


    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    compositeInst->addInstruction(pulseInst);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);

    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}