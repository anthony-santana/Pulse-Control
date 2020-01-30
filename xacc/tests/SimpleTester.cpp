#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

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



int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
