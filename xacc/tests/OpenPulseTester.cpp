#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include "CommonGates.hpp"

namespace
{
    const int NB_SHOTS = 256;
    
    // Expected state vector (2 complex numbers) when driving a single qubit with
    // a Gaussian pulse. 
    // Used for test validation.
    std::vector<std::complex<double>> expectedStateVectorGaussianDrive(double in_time, double in_sigma, double in_driveStrength)
    {
        // Analytical result:
        // [cos(x), -isin(x)] with
        // x = 0.5* sqrt(pi/2) * sigma * drive_strength * ERF(t / (sqrt(2) * sigma))
        const double x =  0.5 * sqrt(M_PI/2) * in_sigma * in_driveStrength * erf(in_time/ (sqrt(2) * in_sigma));
        const std::complex<double> state0 { cos(x), 0.0};
        const std::complex<double> state1 { 0.0, -sin(x)};
        return { state0, state1 };
    }

    const std::string singleQubitHamiltonianJsonTemplate = R"(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["-0.5*omega0*Z0", "omegaa*X0||D0"],
            "osc": {},
            "qub": {
                "0": {{dim}}
            },
            "vars": {
                "omega0": {{omega0}},
                "omegaa": {{omegaa}}
            }
        }
    )";

    // Create a Hamiltonian Json:  
    // in_omega0: frequency of qubit (rad)
    // in_omegaA: drive strength (rad)
    // in_dim(default = 2): dimension of qubit/qutrit
    std::string createHamiltonianJson1Q(double in_omega0, double in_omegaA, int in_dim = 2)
    {
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string omegaaPlaceHolder = "{{omegaa}}";
        const std::string dimlaceHolder = "{{dim}}";
        std::string hamiltonianJson = singleQubitHamiltonianJsonTemplate;
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(in_omega0));
        hamiltonianJson.replace(hamiltonianJson.find(omegaaPlaceHolder), omegaaPlaceHolder.length(), std::to_string(in_omegaA));
        hamiltonianJson.replace(hamiltonianJson.find(dimlaceHolder), dimlaceHolder.length(), std::to_string(in_dim));
        return hamiltonianJson;
    }
}

// Test Open Pulse: i.e. quantum gates by pulses
// Test single qubit gates (square drive)
TEST(SimpleTester, testXGate) 
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    // Params:
    const int total_samples = 100;
    const double omega_0 = 2 * M_PI;
    const double omega_d0 = omega_0;
    const double omega_a = M_PI / total_samples;

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
    
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS) });    

    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    compositeInst->addInstruction(pulseInst);
    
    // Add a measure to get bitstrings
  	auto meas = std::make_shared<xacc::quantum::Measure>(0);
    compositeInst->addInstruction(meas);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
    EXPECT_EQ(finalPopulations.size(), 1);
    EXPECT_NEAR(finalPopulations[0], 1.0, 0.01);

    // Check via bitstrings
    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    EXPECT_NEAR(prob0, 0.0, 0.1);
    EXPECT_NEAR(prob1, 1.0, 0.1);
}

TEST(SimpleTester, testHadamardGate) 
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    // Params:
    const int total_samples = 100;
    const double omega_0 = 2 * M_PI;
    const double omega_d0 = omega_0;
    // Hadamard is rotation by PI/2 (half that of X gate)
    const double omega_a = M_PI / (2 * total_samples);
    // Frame change
    const double phi = -M_PI / 2;
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
    
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS) });    

    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_hadamard_pulse");
    
    // Frame-change by phi = -pi/2    
    auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
    xacc::InstructionParameter fcParameter(phi);
    fcInst->setParameter(0, fcParameter);

    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    
    compositeInst->addInstruction(fcInst);
    compositeInst->addInstruction(pulseInst);
    
    // Request a measurement via 'acquire' pulse instructions 
    auto aqInst = std::make_shared<xacc::quantum::Pulse>("acquire");
    aqInst->setBits({0});
    compositeInst->addInstruction(aqInst);
    
    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
    EXPECT_EQ(finalPopulations.size(), 1);
    // Population ~ 0.5 (i.e. if doing shots, we get half/half probability)
    EXPECT_NEAR(finalPopulations[0], 0.5, 0.01);

    // Check via bitstrings
    qubitReg->print();
    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    // Use a relaxed limit due to limited shot count
    EXPECT_NEAR(prob0, 0.5, 0.1);
    EXPECT_NEAR(prob1, 0.5, 0.1);
}

// Test Gaussian drive
TEST(SimpleTester, testGaussianDrive) 
{
    // Params:
    const int total_samples = 100;
    const double omega_0 = 2 * M_PI;
    const double omega_d0 = omega_0;
    const double omega_a = M_PI / total_samples;  
    // Test cases: different pulse width (Gaussian sigma param)
    const std::vector<int> gaussSigmaSampleCounts { total_samples / 6, total_samples / 3, total_samples };
    
    for (const auto& gaussSigmaSampleCount : gaussSigmaSampleCounts)
    {
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        const std::string hamiltonianJson = createHamiltonianJson1Q(omega_0, omega_a);

        const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
        EXPECT_TRUE(loadOK);    
        BackendChannelConfigs channelConfigs;
    
        channelConfigs.dt = 1.0;
        // LO freq: this is non-zero, hence the drive signal will be modulated. 
        channelConfigs.loFregs_dChannels.emplace_back(omega_d0/(2*M_PI));        
        
        const double  gaussSigma = gaussSigmaSampleCount * channelConfigs.dt;

        // Add the Gaussian pulse
        channelConfigs.addOrReplacePulse("gaussian", QuaC::GaussianPulse(total_samples, gaussSigma));    
        systemModel->setChannelConfigs(channelConfigs);
        
        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

        auto qubitReg = xacc::qalloc(1);    

        auto provider = xacc::getIRProvider("quantum");
        auto compositeInst = provider->createComposite("test_gaussian_pulse");
        auto pulseInst = std::make_shared<xacc::quantum::Pulse>("gaussian", "d0");
        compositeInst->addInstruction(pulseInst);


        const auto expectedStateVector = expectedStateVectorGaussianDrive(total_samples*channelConfigs.dt, gaussSigma, omega_a);
        const auto expectedProb0 = std::pow(std::abs(expectedStateVector[0]), 2.0);
        const auto expectedProb1 = std::pow(std::abs(expectedStateVector[1]), 2.0);
        // Check our analytical results as well
        EXPECT_NEAR(expectedProb0 + expectedProb1, 1.0, 1e-9);

        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(qubitReg, compositeInst);
        const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
        EXPECT_EQ(finalPopulations.size(), 1);
        
        // The simulation result is within 10% error bar (compared to analytical)
        EXPECT_NEAR(finalPopulations[0], expectedProb1, 0.1*expectedProb1);
    }
}

// Do a frame change by PI, confirm that the gate is reverse.
// i.e. X(pi/2) FC(pi) X(pi/2) = Identity (not a full X gate)
TEST(SimpleTester, testFrameChangeCancel) 
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    // Params:
    const int total_samples = 100;
    const double omega_0 = 2 * M_PI;
    const double omega_d0 = omega_0;
    // Hadamard is rotation by PI/2 (half that of X gate)
    const double omega_a = M_PI / (2 * total_samples);
    // Frame change by PI
    const double phi = M_PI;
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
    
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_fc_pulse");
    
    // Frame-change by phi = Pi    
    auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
    xacc::InstructionParameter fcParameter(phi);
    fcInst->setParameter(0, fcParameter);
    // Same square pulse (X(pi/2))
    auto pulseInst1 = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    auto pulseInst2 = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    
    compositeInst->addInstruction(pulseInst1);
    compositeInst->addInstruction(fcInst);
    compositeInst->addInstruction(pulseInst2);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
    EXPECT_EQ(finalPopulations.size(), 1);
    // This should be very close to 0.0 because of the FC
    EXPECT_NEAR(finalPopulations[0], 0.0, 1e-3);
}

// Doing a pi/4 pulse, then pi phase change, then do pi/8 pulse.
// => expect net rotation = pi/4 - pi/8 = pi/8 
TEST(SimpleTester, testFrameChangeNet) 
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    // Params:
    const int duration1 = 100; // Pi/4
    const int duration2 = duration1 / 2; // Pi/8
    // Sample count for Pi rotation (for drive strength normalization)
    const int piSamples = 4*duration1;
    const double omega_0 = 2 * M_PI;
    const double omega_d0 = omega_0;
    const double omega_a = M_PI / piSamples;
    // Frame change by PI
    const double phi = M_PI;
    const std::string hamiltonianJson = createHamiltonianJson1Q(omega_0, omega_a);

    const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
    EXPECT_TRUE(loadOK);
    
    BackendChannelConfigs channelConfigs;
    
    channelConfigs.dt = 1.0;
    // LO freq: this is non-zero, hence the drive signal will be modulated. 
    channelConfigs.loFregs_dChannels.emplace_back(omega_d0/(2*M_PI));
      
    // Add the square pulses
    channelConfigs.addOrReplacePulse("square1", QuaC::SquarePulse(duration1));    
    channelConfigs.addOrReplacePulse("square2", QuaC::SquarePulse(duration2));    
    systemModel->setChannelConfigs(channelConfigs);
    
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_fc_pulse_net");
    
    // Frame-change by phi = Pi    
    auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
    xacc::InstructionParameter fcParameter(phi);
    fcInst->setParameter(0, fcParameter);
    // Two square pulses
    auto pulseInst1 = std::make_shared<xacc::quantum::Pulse>("square1", "d0");
    auto pulseInst2 = std::make_shared<xacc::quantum::Pulse>("square2", "d0");
    
    compositeInst->addInstruction(pulseInst1);
    compositeInst->addInstruction(fcInst);
    compositeInst->addInstruction(pulseInst2);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
    const auto expectedVal = 1.0 - std::pow(std::cos((M_PI/4.0 - M_PI / 8.0) / 2.0), 2.0);
    EXPECT_EQ(finalPopulations.size(), 1);
    const double error = std::abs(expectedVal  - finalPopulations[0])/expectedVal;
    // relative error must be less than 10% from the analytical expectation
    EXPECT_LT(error, 0.1);
}


int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
