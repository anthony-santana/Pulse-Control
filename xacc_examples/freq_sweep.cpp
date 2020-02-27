#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "xacc_service.hpp"
#include "PulseSystemModel.hpp"

// ============================================================================
// Pulse Demo 1: Finding the qubit frequency using a frequency sweep
// ============================================================================
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    const std::string hamiltonianJson = R"(
    {
        "description": "Qubits are modelled as a two level system. System of 2 qubits.\n",
        "h_str": ["_SUM[i,0,1,wq{i}/2*Z{i}]", "_SUM[i,0,1,omegad{i}*X{i}||D{i}]", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
        "osc": {},
        "qub": {
            "0": 2,
            "1": 2
        },
        "vars": {
            "omegad0": 1.303125,
            "omegad1": 0.97, 
            "wq0": 30.91270129264568,
            "wq1": 30.36010168900955,
            "jq0q1": 0.04
        }
    })";
    
    // Sweeping Q1 freq
    const auto freq_q1 = xacc::linspace(30.7, 31.0, 50);
    std::vector<double> prob1VecQ1;
    
    // Sweeping Q2 freq
    const auto freq_q2 = xacc::linspace(30.15, 30.45, 50);
    std::vector<double> prob1VecQ2;


    bool loadPulseCmdDef = false;
    auto xasmCompiler = xacc::getCompiler("xasm");    
    
    // X gate on the first qubit (sweeping LO freq)
    auto ir1 = xasmCompiler->compile(R"(__qpu__ void testQ1(qbit q) {
        X(q[0]);
        Measure(q[0]);
    })");    
    auto program1 = ir1->getComposite("testQ1");

    // X gate on the second qubit (sweeping LO freq)
    auto ir2 = xasmCompiler->compile(R"(__qpu__ void testQ2(qbit q) {
        X(q[1]);
        Measure(q[1]);
    })");    
    auto program2 = ir2->getComposite("testQ2");
    
    // First qubit
    for (const auto& freq : freq_q1)
    {
        std::vector<double> loFreqs {
            freq,
            // We don't drive second qubit, hence this freq doesn't matter
            30.36010168900955
        };
        
        // We don't use U channels (cross-resonance) at the moment.
        std::vector<double> loFreqs_uChannels(20, 0.0);

        for (auto& freq : loFreqs)
        {
            freq = freq/(2*M_PI);
        }

        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        const bool loadOk = systemModel->loadHamiltonianJson(hamiltonianJson);
        assert(loadOk);
        
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 3.5555555555555554;
        channelConfigs.loFregs_dChannels = loFreqs;
        channelConfigs.loFregs_uChannels = loFreqs_uChannels;
        systemModel->setChannelConfigs(channelConfigs);    

        const int NB_SHOTS = 10000;
        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), 
                                                    std::make_pair("shots", NB_SHOTS) });
        
        
        if (!loadPulseCmdDef)
        {
            // Contribute all cmd-defs in the backend Json as XACC instruction.
            // This will activate Gates -> Pulses decomposition when simulating the circuit.
            const std::string jsonConfigFile = std::string(GATEIR_TEST_FILE_DIR) + "/test_backends.json";
            std::ifstream backendFile(jsonConfigFile);
            std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
            quaC->contributeInstructions(jjson);   
            loadPulseCmdDef = true;  
        }
    
        // Note: this may take a long time to run.
        auto qubitReg = xacc::qalloc(2);
    
        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(qubitReg, program1);
        
        prob1VecQ1.emplace_back(qubitReg->computeMeasurementProbability("1"));
    }

    // Second qubit
    for (const auto& freq : freq_q2)
    {
        std::vector<double> loFreqs {
            // We don't drive first qubit, hence this freq doesn't matter
            30.91270129264568,
            freq
        };
        
        // We don't use U channels (cross-resonance) at the moment.
        std::vector<double> loFreqs_uChannels(20, 0.0);

        for (auto& freq : loFreqs)
        {
            freq = freq/(2*M_PI);
        }

        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        const bool loadOk = systemModel->loadHamiltonianJson(hamiltonianJson);
        assert(loadOk);
        
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 3.5555555555555554;
        channelConfigs.loFregs_dChannels = loFreqs;
        channelConfigs.loFregs_uChannels = loFreqs_uChannels;
        systemModel->setChannelConfigs(channelConfigs);    

        const int NB_SHOTS = 10000;
        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), 
                                                    std::make_pair("shots", NB_SHOTS) });
        
        // Note: this may take a long time to run.
        auto qubitReg = xacc::qalloc(2);
    
        // Run the Pulse simulation for second program (X on q[1])
        quaC->execute(qubitReg, program2);
        
        prob1VecQ2.emplace_back(qubitReg->computeMeasurementProbability("1"));
    }

    // First qubit's sweep result
    for (int i = 0; i < prob1VecQ1.size(); ++i)
    {
        std::cout << " f = " << freq_q1[i] << "; P(1) = " << prob1VecQ1[i] << "\n";
    }
    
    // Second qubit's sweep result
    for (int i = 0; i < prob1VecQ2.size(); ++i)
    {
        std::cout << " f = " << freq_q2[i] << "; P(1) = " << prob1VecQ2[i] << "\n";
    }
    
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}