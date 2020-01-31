#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

// NOTE: these tests are designed to validate time-independent evolution,
// we will have tests for cases where we have time-dependent terms.

TEST(SimpleTester, checkOneQubitX) 
{
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

    const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
    EXPECT_TRUE(loadOK);
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 0.005;
    // LO freq = Cavity Freq = 0.0
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    // delta = 2*pi => period = 1 => PI pulse (X) = 0.5
    // 100 samples * dt (0.005) = 0.5 => PI pulse
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
    const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
    EXPECT_EQ(finalPopulations.size(), 1);
    // Final population should be very close to 1 (population inversion)
    EXPECT_NEAR(finalPopulations[0], 1.0, 0.05);
}

TEST(SimpleTester, checkOneQubitXscaling) 
{
    const std::string hamiltonianJsonTemplate = R"(
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
                "delta": {{delta}}
            }
        }
    )";

    const std::string deltaVar = "{{delta}}";
    // Scale the freq by these factors
    std::vector<double> scales { 2.0, 0.1234, 100 };
    
    for (const auto& scaleFactor : scales)
    {
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        auto hamiltonianJson = hamiltonianJsonTemplate;
        // frequency divide by scale factor
        const double deltaVal = 2.0 * M_PI / scaleFactor;
        hamiltonianJson.replace(hamiltonianJson.find(deltaVar), deltaVar.length(), std::to_string(deltaVal));
        const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
        EXPECT_TRUE(loadOK);
        BackendChannelConfigs channelConfigs;
        // dt times scale factor => should be unchanged.
        channelConfigs.dt = 0.005 * scaleFactor;
        // LO freq = Cavity Freq = 0.0
        channelConfigs.loFregs_dChannels.emplace_back(0.0);
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
        const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
        
        // Final population should be very close to 1 (population inversion)
        EXPECT_EQ(finalPopulations.size(), 1);
        EXPECT_NEAR(finalPopulations[0], 1.0, 0.01);
    }
}

TEST(SimpleTester, checkTwoBodyHamiltonian) 
{
    const std::vector<double> durations { 5.0, 10.0, 15.0, 20.0, 25.0 };
    double totalPopulation = 1.0;
    for (const auto& duration : durations)
    {
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        // Jaynes-Cumming model atom + cavity coupling
        // https://nbviewer.jupyter.org/github/qutip/qutip-notebooks/blob/master/examples/rabi-oscillations.ipynb
        // Cavity (Q0) is 15-level
        // Atom (Q1) is 2-level
        const std::string hamiltonianJson = R"#(
            {
                "description": "Jaynes-Cumming Hamiltonian.",
                "h_latex": "",
                "h_str": ["wc*O0", "wa*O1", "g*(SP0*SM1 + SM0*SP1)"],
                "osc": {},
                "qub": {
                    "0": 15,
                    "1": 2
                },
                "vars": {
                    "wc": 6.2831853,
                    "wa": 6.2831853,
                    "g": 0.314159
                }
            }
        )#";

        const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
        EXPECT_TRUE(loadOK);
        // cavity dissipation rate
        const double kappa = 0.005;
        systemModel->setQubitT1(0, 1.0/kappa);
        // atom dissipation rate
        const double gamma = 0.05;           
        systemModel->setQubitT1(1, 1.0/gamma);
        // Set the atom to the excited state:
        systemModel->setQubitInitialPopulation(1, 1.0);
        
        BackendChannelConfigs channelConfigs;
        // We run the simulation to the fixed duration
        const size_t nbSamples = 100;        
        channelConfigs.dt = duration/nbSamples;
        // LO freq = Cavity Freq = 0.0
        channelConfigs.loFregs_dChannels.emplace_back(0.0);
        channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));    
        systemModel->setChannelConfigs(channelConfigs);
        
        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    
        // 2 'qubit' system (the cavity is many-level system)
        auto qubitReg = xacc::qalloc(2);    

        auto provider = xacc::getIRProvider("quantum");
        auto compositeInst = provider->createComposite("test_pulse");
        auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
        compositeInst->addInstruction(pulseInst);

        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(qubitReg, compositeInst);
        const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
        EXPECT_EQ(finalPopulations.size(), 2);
        // Current atom + cavity population
        const double totalPop = finalPopulations[0] + finalPopulations[1];
        // It must decay (less than that of the previous run for a shorter duration)
        EXPECT_TRUE(totalPopulation > totalPop);
        totalPopulation = totalPop;
    }
}

// This test will also validate our Hamiltonian parsing (SUM terms)
TEST(SimpleTester, checkSpinChain) 
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

    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    
    // 10 spin qubits
    auto qubitReg = xacc::qalloc(10);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    compositeInst->addInstruction(pulseInst);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
    EXPECT_EQ(finalPopulations.size(), 10);
    // All populations should be pretty small (corresponding to |0> state, +1 expectation-Z)
    for (const auto& data: finalPopulations)
    {
        EXPECT_TRUE(data < 0.5);
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
