#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    
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
                "omegaa": 0.02
            }
        }
    )#";

    if (!systemModel->loadHamiltonianJson(singleQubitHamiltonianJson))
    {
        return -1;
    }
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 0.1;
    const size_t nbSamples = 1000;
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    const std::string fourierSeries = R"#(
        0.726924 + 
        0.065903*cos(1*0.1*t) + 0.128627*sin(1*0.1*t) +
        0.079360*cos(2*0.1*t) + 0.111686*sin(2*0.1*t) + 
        0.096717*cos(3*0.1*t) + 0.096822*sin(3*0.1*t) + 
        0.106937*cos(4*0.1*t) + 0.092216*sin(4*0.1*t) + 
        0.215306*cos(5*0.1*t) + 0.118562*sin(5*0.1*t) +
        0.117682*cos(6*0.1*t) + 0.126134*sin(6*0.1*t) + 
        0.100447*cos(7*0.1*t) + 0.120409*sin(7*0.1*t) + 
        0.103292*cos(8*0.1*t) + 0.108712*sin(8*0.1*t))#";
    
    channelConfigs.addOrReplacePulse("fourier", QuaC::PulseFunc(fourierSeries, nbSamples, channelConfigs.dt));
    
    systemModel->setChannelConfigs(channelConfigs);
    

    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("fourier", "d0");
    compositeInst->addInstruction(pulseInst);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    qubitReg->print();
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}