#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include <chrono>
#include <random>

// This test will also validate our Hamiltonian parsing (SUM terms)
// Run the 10-qubit Heisenberg spin chain simulation with MPI
TEST(QuaCMPITester, checkSpinChain) 
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();

    // Ref https://nbviewer.jupyter.org/github/qutip/qutip-notebooks/blob/master/examples/spin-chain.ipynb
    // 10-qubit spin chain
    // Note: the sum term format is that of IBM
    const std::string hamiltonianJson = R"(
        {
            "description": "Heisenberg spin-chain Hamiltonian.",
            "h_latex": "",
            "h_str": ["_SUM[i,0,9,-0.5*h{i}*Z{i}]", "_SUM[i,0,8,-0.5*Jx{i}*X{i}*X{i+1}]", "_SUM[i,0,8,-0.5*Jy{i}*Y{i}*Y{i+1}]", "_SUM[i,0,8,-0.5*Jz{i}*Z{i}*Z{i+1}]"],
            "osc": {},
            "qub": {
                "0": 2
            },
            "vars": {
                "h0": 6.2831853,
                "h1": 6.2831853,
                "h2": 6.2831853,
                "h3": 6.2831853,
                "h4": 6.2831853,
                "h5": 6.2831853,
                "h6": 6.2831853,
                "h7": 6.2831853,
                "h8": 6.2831853,
                "h9": 6.2831853,
                "Jx0": 0.62831853,
                "Jx1": 0.62831853,
                "Jx2": 0.62831853,
                "Jx3": 0.62831853,
                "Jx4": 0.62831853,
                "Jx5": 0.62831853,
                "Jx6": 0.62831853,
                "Jx7": 0.62831853,
                "Jx8": 0.62831853,
                "Jy0": 0.62831853,
                "Jy1": 0.62831853,
                "Jy2": 0.62831853,
                "Jy3": 0.62831853,
                "Jy4": 0.62831853,
                "Jy5": 0.62831853,
                "Jy6": 0.62831853,
                "Jy7": 0.62831853,
                "Jy8": 0.62831853,
                "Jz0": 0.62831853,
                "Jz1": 0.62831853,
                "Jz2": 0.62831853,
                "Jz3": 0.62831853,
                "Jz4": 0.62831853,
                "Jz5": 0.62831853,
                "Jz6": 0.62831853,
                "Jz7": 0.62831853,
                "Jz8": 0.62831853
            }
        }
    )";

    const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
    EXPECT_TRUE(loadOK);
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 0.5;
    // LO freq = Cavity Freq = 0.0
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    const size_t nbSamples = 100;

    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));
    
    systemModel->setChannelConfigs(channelConfigs);
    
    // intial state, first spin in state |1>, the rest in state |0>
    systemModel->setQubitInitialPopulation(0, 1.0);
    // dephasing rate
    const double gamma = 0.01;
    for (int i = 0; i < 10; ++i)
    {
        systemModel->setQubitT1(i, 1.0/gamma);
    }

    auto quaC = xacc::getAccelerator("QuaC", 
            {   std::make_pair("system-model", systemModel), 
                // Request MPI executor!
                std::make_pair("execution-mode", "MPI") 
            });  

    // 10 spin qubits
    auto qubitReg = xacc::qalloc(10);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    compositeInst->addInstruction(pulseInst);
    // Time the execution to compare with the same test in SimpleTester (single process)
    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "[MPI] Heisenberg spin chain elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
    EXPECT_EQ(finalPopulations.size(), 10);
    // All populations should be pretty small (corresponding to |0> state, +1 expectation-Z)
    qubitReg->print();
    for (const auto& data: finalPopulations)
    {
        EXPECT_TRUE(data < 0.25);
    }
}

// Test the X gate for IBM backend
// To reduce test time, just use 5 qubits of the 20 qubits
TEST(QuaCMPITester, testXGateIBMLib) 
{
    const std::string hamiltonianJson = R"(
    {
        "description": "Qubits are modelled as a two level system.\n",
        "h_str": ["_SUM[i,0,4,wq{i}/2*Z{i}]", "_SUM[i,0,4,omegad{i}*X{i}||D{i}]"],
        "osc": {},
        "qub": {
            "0": 2,
            "1": 2,
            "2": 2,
            "3": 2,
            "4": 2
        },
        "vars": {
            "omegad0": 1.25,
            "omegad1": 0.97, 
            "omegad2": 2.5,
            "omegad3": 1.05,
            "omegad4": 0.845,
            "wq0": 30.91270129264568,
            "wq1": 30.36010168900955,
            "wq2": 31.041771660759178,
            "wq3": 28.36701429077905,
            "wq4": 29.298278199939336
        }
    })";
    
    std::vector<double> loFreqs {
        30.91270129264568,
        30.36010168900955,
        31.041771660759178,
        28.36701429077905,
        29.298278199939336
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
    const int NB_QUBITS = 5;

    // Random select a qubit to run
    std::random_device dev;
    std::mt19937 rng(dev());
    // Select a qubit from 0->4
    std::uniform_int_distribution<std::mt19937::result_type> dist(0, NB_QUBITS - 1);
    const int qubitIdx = dist(rng);
    std::cout << "Perform X[q" << qubitIdx << "]\n";
    
    // Template for the XASM (randomly assign a qubit to the X gate)
    const std::string xasmTmpl = 
        R"(__qpu__ void testX(qbit q) {
            X(q[%%Idx%%%]);  
            Measure(q[0]);
            Measure(q[1]);
            Measure(q[2]);
            Measure(q[3]);
            Measure(q[4]);
        })";

    // Randomly pick a qubit
    const std::string qubitIdxPlaceholder = "%%Idx%%%";
    std::string randomXasm = xasmTmpl;
    randomXasm.replace(randomXasm.find(qubitIdxPlaceholder), qubitIdxPlaceholder.length(), std::to_string(qubitIdx));

    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), 
                                                std::make_pair("execution-mode", "MPI"),
                                                std::make_pair("shots", NB_SHOTS) });
    
    // Load cmd-def from IBM Json
    const std::string jsonConfigFile = std::string(GATEIR_TEST_FILE_DIR) + "/test_backends.json";
    std::ifstream backendFile(jsonConfigFile);
    std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
    quaC->contributeInstructions(jjson);     
   
    auto qubitReg = xacc::qalloc(NB_QUBITS);
    auto provider = xacc::getIRProvider("quantum");
    auto xasmCompiler = xacc::getCompiler("xasm");
    
    // Run the XASM circuit (randomly pick a qubit for X gate)
    auto ir = xasmCompiler->compile(randomXasm, quaC);    
    auto program = ir->getComposite("testX"); 
    quaC->execute(qubitReg, program);

    // Construct the expected bit string: '1' for the qubit that was picked, 
    // '0' for all others.
    std::string expectedBitString;
    for (int i = 0; i < NB_QUBITS; ++i)
    {
        expectedBitString.push_back(i == qubitIdx ? '1' : '0');
    }

    const double probResult = qubitReg->computeMeasurementProbability(expectedBitString);
    
    // The probability of the expected bitstring should be close to 1.0
    EXPECT_NEAR(probResult, 1.0, 0.1);
    // Debug:
    qubitReg->print();
    // Finalize the XACC Framework
	xacc::Finalize();
}


int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
