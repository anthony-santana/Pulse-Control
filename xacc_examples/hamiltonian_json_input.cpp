#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "xacc_service.hpp"
#include "PulseSystemModel.hpp"

// ============================================================================
// Using cmd-def's of single-qubit gates from IBM Poughkeepsie backend 
// ============================================================================
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    // NOTES: 
    // (1) Although we will load all cmd-def's from the IBM Json,
    // we need to *modify* the Hamiltonian to make it work.
    // This is because all the driving terms in the IBM Json have ZERO coefficients
    // which doesn't make any sense (we cannot drive the qubits).
    // For this example, we use a trial-and-error approach to set these coefficients
    // based on the cmd-def for X gates.
    // i.e. we do X(q[i]) which will look up the cmd-def and associated pulses from the 
    // IBM Json, drive the qubits with those pulses, and tune the coefficients to 
    // make the X gate.  
    // Hence, the omegad0, omegad1, etc. in the below JSON are approximate that we 
    // have calibrated (very roughly, upto 10% error rate for X gates) to match
    // those pulses defined in the JSON files.
    // (2) We only use 10-qubit for demonstration, 
    //  i.e. we will only ever use a subset of cmd-defs defined.
    const std::string hamiltonianJson = R"(
    {
        "description": "Qubits are modelled as a two level system.\n",
        "h_str": ["_SUM[i,0,9,wq{i}/2*Z{i}]", "_SUM[i,0,9,omegad{i}*X{i}||D{i}]"],
        "osc": {},
        "qub": {
            "0": 2,
            "1": 2,
            "2": 2,
            "3": 2,
            "4": 2,
            "5": 2,
            "6": 2,
            "7": 2,
            "8": 2,
            "9": 2
        },
        "vars": {
            "omegad0": 1.25,
            "omegad1": 0.97, 
            "omegad2": 2.5,
            "omegad3": 1.05,
            "omegad4": 0.845,
            "omegad5": 1.25,
            "omegad6": 1.25,
            "omegad7": 0.97,
            "omegad8": 1.33,
            "omegad9": 1.45,
            "wq0": 30.91270129264568,
            "wq1": 30.36010168900955,
            "wq2": 31.041771660759178,
            "wq3": 28.36701429077905,
            "wq4": 29.298278199939336,
            "wq5": 31.14757431485366,
            "wq6": 31.387741224914162,
            "wq7": 30.232349262897678,
            "wq8": 31.502130591468386,
            "wq9": 31.769632280927663
        }
    })";
    
    std::vector<double> loFreqs {
        30.91270129264568,
        30.36010168900955,
        31.041771660759178,
        28.36701429077905,
        29.298278199939336,
        31.14757431485366,
        31.387741224914162,
        30.232349262897678,
        31.502130591468386,
        31.769632280927663
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
    
    
    // Contribute all cmd-defs in the backend Json as XACC instruction.
    // This will activate Gates -> Pulses decomposition when simulating the circuit.
    const std::string jsonConfigFile = std::string(GATEIR_TEST_FILE_DIR) + "/test_backends.json";
    std::ifstream backendFile(jsonConfigFile);
    std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
    quaC->contributeInstructions(jjson);     
    
    // Note: this may take a long time to run.
    auto qubitReg = xacc::qalloc(10);
    auto provider = xacc::getIRProvider("quantum");
    auto xasmCompiler = xacc::getCompiler("xasm");
    
    // Hadamard by pulses:
    // Note: we only calibrate the Hamiltonian by the X gate (pulse named Xp_....)
    // we expect the calibrated Hamiltonian coefficients should also work reasonably well
    // for other pulses defined in this backend library.
    // e.g the Hadamard gate will be decomposed into X90p... and X90m... pulses along with a few frame changes.  
    auto ir = xasmCompiler->compile(R"(__qpu__ void test(qbit q) {
      H(q[0]);
      Measure(q[0]);
    })", quaC);
    
    auto program = ir->getComposite("test");
  
    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, program);

    qubitReg->print();
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}