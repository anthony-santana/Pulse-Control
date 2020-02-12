#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include <chrono>

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


int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
