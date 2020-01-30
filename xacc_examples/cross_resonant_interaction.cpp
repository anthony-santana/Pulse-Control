
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"


// Based on Eugene's Cross_Resonant_Interaction IPython Notebook
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();

    // Hamiltonian: (w_0*n(N,D,0) + d*n(N,D,0)*(n(N,D,0)-1) + w_1*n(N,D,1) + d*n(N,D,1)*(n(N,D,1)-1)) 
    // + J*(a(N,D,0).dag()*a(N,D,1) + a(N,D,0)*a(N,D,1).dag()) + O*(a(N,D,1).dag()+a(N,D,1))*cos(w_C * t)
    // Note: we use IBM notation O/P operators are the number operator (N)
    // "h_str": ["w_0*O0", "d*O0*(O0-I0)", "w_1*O1", "d*O1*(O1-I1)", "J*(SP0*SM1 + SM0*SP1)", "O*(SM0 + SP0)||D0", "O*(SM1 + SP1)||D1"],
    const std::string hamiltonianJson = R"#(
        {
            "description": "Two-qutrit Hamiltonian.",
            "h_latex": "",
            "h_str": ["w_0*O0", "d*O0*(O0-I0)", "w_1*O1", "d*O1*(O1-I1)", "J*(SP0*SM1 + SM0*SP1)", "O*(SM0 + SP0)||D0"],
            "osc": {},
            "qub": {
                "0": 3,
                "1": 3
            },
            "vars": {
                "w_0": 5.114,
                "w_1": 4.914,
                "d": -0.33,
                "J": 0.004,
                "O": 0.060
            }
        }
    )#";

    if (!systemModel->loadHamiltonianJson(hamiltonianJson))
    {
        return -1;
    }
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 5.0;
    // LO freq = [ w_1/(2*pi), w_0/(2*pi)]
    channelConfigs.loFregs_dChannels.emplace_back(4.914*M_1_PI/2.0);
    channelConfigs.loFregs_dChannels.emplace_back(5.114*M_1_PI/2.0);

    const size_t nbSamples = 1000;
    // Constant drive
    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));
    
    systemModel->setChannelConfigs(channelConfigs);
    
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

    auto qubitReg = xacc::qalloc(2);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInstD0 = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    auto pulseInstD1 = std::make_shared<xacc::quantum::Pulse>("square", "d1");

    compositeInst->addInstruction(pulseInstD0);
    compositeInst->addInstruction(pulseInstD1);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);

    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}