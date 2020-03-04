#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include "CommonGates.hpp"

namespace {
    // Three-level single-qubit Hamiltonian
    // This is the same as a regular 2-level Hamiltonian but Pauli X and Z operators
    // are converted to ladder operators.
    const std::string singleQubitHamiltonianJsonTemplate = R"(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["-0.5*omega0*I0", "omega0*O0", "omegaa*(SM0+SP0)||D0"],
            "osc": {},
            "qub": {
                "0": 3
            },
            "vars": {
                "omega0": {{omega0}},
                "omegaa": {{omegaa}}
            }
        }
    )";

   
    std::string createHamiltonianJson1Q(double in_omega0, double in_omegaA)
    {
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string omegaaPlaceHolder = "{{omegaa}}";
        std::string hamiltonianJson = singleQubitHamiltonianJsonTemplate;
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(in_omega0));
        hamiltonianJson.replace(hamiltonianJson.find(omegaaPlaceHolder), omegaaPlaceHolder.length(), std::to_string(in_omegaA));
        return hamiltonianJson;
    }
}


TEST(MultiLevelTester, testThreeLevel)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    // Params:
    const int total_samples = 100;
    const double omega_0 = 2 * M_PI;
    const double omega_d0 = omega_0;
    
    // Max |1> population: 
    // Expected: Pop(0) = 4/9; Pop(1) = 1/3; Pop(2) = 2/9 
    const double omega_a = M_PI / (std::sqrt(3.0) * total_samples);

    const std::string hamiltonianJson = createHamiltonianJson1Q(omega_0, omega_a);

    const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
    EXPECT_TRUE(loadOK);
     BackendChannelConfigs channelConfigs;
    
    channelConfigs.dt = 1.0;
    // LO freq: this is non-zero, hence the drive signal will be modulated. 
    channelConfigs.loFregs_dChannels.emplace_back(omega_d0/(2*M_PI));
      
    // Add the square pulse
    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(total_samples));    
    systemModel->setChannelConfigs(channelConfigs);
    const int NB_SHOTS = 1000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS) });    

    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    compositeInst->addInstruction(pulseInst);
       
    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    qubitReg->print();
    const auto finalPopulations = qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
    // We expect to have 3 numbers
    EXPECT_EQ(finalPopulations.size(), 3);
     // Expected: Pop(0) = 4/9; Pop(1) = 1/3; Pop(2) = 2/9 
    EXPECT_NEAR(finalPopulations[0], 4.0/9.0, 0.01);
    EXPECT_NEAR(finalPopulations[1], 1.0/3.0, 0.01);
    EXPECT_NEAR(finalPopulations[2], 2.0/9.0, 0.01);
}


int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}