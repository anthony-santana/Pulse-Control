
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();

    // An example one-qubit Hamiltonian
    // H = epsilon / 2.0 * sigmaz() + delta / 2.0 * sigmax()
    // Ref: https://nbviewer.jupyter.org/github/qutip/qutip-notebooks/blob/master/examples/qubit-dynamics.ipynb
    // epsilon = cavity frequency = 0.0
    // delta = atom frequency = 2*pi
    const std::string hamiltonianJson = R"(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["0.5*epsilon*Z0", "0.5*delta*X0||D0"],
            "osc": {},
            "qub": {
                "0": 2
            },
            "vars": {
                "epsilon": 0.0,
                "delta": 6.2831853
            }
        }
    )";

    if (!systemModel->loadHamiltonianJson(hamiltonianJson))
    {
        return -1;
    }
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 0.05;
    // LO freq = Cavity Freq = 0.0
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    // 100 samples * dt (0.05) = 5
    const size_t nbSamples = 100;

    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));
    
    systemModel->setChannelConfigs(channelConfigs);
    // Qubit decay gamma
    const double gamma = 0.15;
    systemModel->setQubitT1(0, 1.0/gamma);

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