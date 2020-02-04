#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include "CommonGates.hpp"

namespace {
    const int NB_SHOTS = 256;
    
    const std::string twoQubitHamiltonianJsonTemplate = R"#(
        {
            "description": "Two-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["-0.5*omega0*Z0", "omegaa*X0||D0", "omegai*(SP0*SM1+SM0*SP1)||U1"],
            "osc": {},
            "qub": {
                "0": {{dim}},
                "1": {{dim}}
            },
            "vars": {
                "omega0": {{omega0}},
                "omegaa": {{omegaa}},
                "omegai": {{omegai}}
            }
        }
    )#";

    // Create a 2-qubit Hamiltonian Json:  
    // in_omega0: frequency of qubit (rad)
    // in_omegaA: drive strength (rad)
    // in_omegaI: control (U) strength (rad)
    // in_dim(default = 2): dimension of qubit/qutrit
    std::string createHamiltonianJson2Q(double in_omega0, double in_omegaA, double in_omegaI, int in_dim = 2)
    {
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string omegaaPlaceHolder = "{{omegaa}}";
        const std::string omegaiPlaceHolder = "{{omegai}}";
        const std::string dimlaceHolder = "{{dim}}";
        std::string hamiltonianJson = twoQubitHamiltonianJsonTemplate;
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(in_omega0));
        hamiltonianJson.replace(hamiltonianJson.find(omegaaPlaceHolder), omegaaPlaceHolder.length(), std::to_string(in_omegaA));
        hamiltonianJson.replace(hamiltonianJson.find(omegaiPlaceHolder), omegaiPlaceHolder.length(), std::to_string(in_omegaI));
        while (hamiltonianJson.find(dimlaceHolder) != std::string::npos)
        {
            hamiltonianJson.replace(hamiltonianJson.find(dimlaceHolder), dimlaceHolder.length(), std::to_string(in_dim));
        }
        return hamiltonianJson;
    }

    // Two-qubit cross-resonance model (qubit dimension = 2)
    // Ref: Rigetti and Devoret: PRB 81, 134507 (2010) 
    // Hamiltonian is from Eq. (1)
    const std::string twoQubitCrossResonanceHamiltonianJsonTemplate = R"#(
        {
            "description": "Two-qubit cross-resonance Hamiltonian.",
            "h_latex": "",
            "h_str": ["0.5*omega0*Z0", "Omega0*X0||D0", "0.5*omega1*Z1", "Omega1*X1||D1", "0.5*omegaxx*X0*X1"],
            "osc": {},
            "qub": {
                "0": 2,
                "1": 2
            },
            "vars": {
                "omega0": {{omega0}},
                "Omega0": {{Omega0}},
                "omega1": {{omega1}},
                "Omega1": {{Omega1}},
                "omegaxx": {{omegaxx}}
            }
        }
    )#";

    std::string createHamiltonianJsonCrossResonance(double in_omega0, double in_Omega0, double in_omega1, double in_Omega1, double in_omegaxx)
    {
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string Omega0PlaceHolder = "{{Omega0}}";
        const std::string omega1PlaceHolder = "{{omega1}}";
        const std::string Omega1PlaceHolder = "{{Omega1}}";
        const std::string omegaxxPlaceHolder = "{{omegaxx}}";

        std::string hamiltonianJson = twoQubitCrossResonanceHamiltonianJsonTemplate;
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(in_omega0));
        hamiltonianJson.replace(hamiltonianJson.find(Omega0PlaceHolder), Omega0PlaceHolder.length(), std::to_string(in_Omega0));

        hamiltonianJson.replace(hamiltonianJson.find(omega1PlaceHolder), omega1PlaceHolder.length(), std::to_string(in_omega1));
        hamiltonianJson.replace(hamiltonianJson.find(Omega1PlaceHolder), Omega1PlaceHolder.length(), std::to_string(in_Omega1));

        hamiltonianJson.replace(hamiltonianJson.find(omegaxxPlaceHolder), omegaxxPlaceHolder.length(), std::to_string(in_omegaxx));
        
        return hamiltonianJson;
    }
}

// Test 2-qubit interaction via a SWAP gate
TEST(CircuitTester, testInteraction)
{
    // Run the SWAP test, if shouldSwap = true then we add the SWAP pulse,
    // otherwise, we don't add the SWAP pulse (to verify that the SWAP pulse is indeed effective)
    const auto runSwapTest = [](bool shouldSwap){
        // Params
        const int total_samples = 100;
        const double omegaiSwap = M_PI / total_samples;
        const double omega0 = 2 * M_PI;
        // We want a Pi pulse on Q0 then swap to verify state swaps from 01 -> 10 
        const double omegaa = M_PI / total_samples;

        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        const std::string hamiltonianJson = createHamiltonianJson2Q(omega0, omegaa, omegaiSwap);

        const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
        EXPECT_TRUE(loadOK);
        
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 1.0;
        // D0 Channel: drive on resonance
        const auto omega_d0 = omega0;
        const auto omega_d1 = 0.0;
        channelConfigs.loFregs_dChannels.emplace_back(omega_d0/(2*M_PI));
        channelConfigs.loFregs_dChannels.emplace_back(omega_d1/(2*M_PI));

        // U Channels:
        channelConfigs.loFregs_uChannels.emplace_back(omega_d0/(2*M_PI));
        channelConfigs.loFregs_uChannels.emplace_back((omega_d1 - omega_d0)/(2*M_PI));
        // Add the square pulse
        channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(total_samples));    
        systemModel->setChannelConfigs(channelConfigs);
        //systemModel->setQubitInitialPopulation(0, 1.0);

        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS) });    
        // Need 2 qubits for this
        auto qubitReg = xacc::qalloc(2);    

        auto provider = xacc::getIRProvider("quantum");
        auto compositeInst = provider->createComposite("test_swap" + std::to_string(shouldSwap));
        // Note: we can use the same square pulse on both D0 and U1 channels
        // because the drive strength has been calculated to produce the correct amount of rotation. 
        auto xInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
        auto interactionInst = std::make_shared<xacc::quantum::Pulse>("square", "u1");
        interactionInst->setStart(total_samples);
        // Create a sequence of pulses for testing two qubit interaction. 
        // - do a PI pulse on Q0 (X gate)
        // - drive U1 channel to start Q0-Q1 SWAP interaction 
        compositeInst->addInstruction(xInst);
        // Only add SWAP pulse if requested
        if (shouldSwap)
        {
            compositeInst->addInstruction(interactionInst);
        }
        
        // Measure both qubits
        auto meas0 = std::make_shared<xacc::quantum::Measure>(0);
        auto meas1 = std::make_shared<xacc::quantum::Measure>(1);

        compositeInst->addInstruction(meas0);
        compositeInst->addInstruction(meas1);

        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(qubitReg, compositeInst);
        
        // Check via bitstrings
        const double prob01 = qubitReg->computeMeasurementProbability("01");
        const double prob10 = qubitReg->computeMeasurementProbability("10");
        EXPECT_NEAR(prob01 + prob10, 1.0, 1e-3);
        if (shouldSwap)
        {
            // SWAP was run
            // The final result should be Q0 = 0; Q1 = 1 after the SWAP, hence bitstring = '01'
            // Note: bit string is constructed by 'push_back' hence Q0 will be the first bit in the string.
            EXPECT_NEAR(prob01, 1.0, 1e-3);
        }
        else
        {
            // If no SWAP, we should get an '10' bitstring (just a PI pulse on Q0)
            EXPECT_NEAR(prob10, 1.0, 1e-3);
        }
    };
    

    // Run SWAP and non-SWAP tests
    runSwapTest(true);
    runSwapTest(false);
} 

// 2-qubit cross resonance test
TEST(CircuitTester, testCrossResonance)
{
    // Simulate the entanglement concurrence
    // see inset of Fig. (3) of PRB 81, 134507 (2010) 
    const std::vector<int> sampleCounts {15, 30, 50};    
    const std::vector<double> expectedResults { 0.6, 0.9, 0.4};  
    
    for (int i = 0; i < sampleCounts.size(); ++i)
    {
        // Vary the simulation endtime
        const int total_samples = sampleCounts[i];        
        // Coupling
        const double omegaXX = 0.13;
        // Qubit 1: 
        const double omega1 = 12.8;
        // Qubit 2: 
        const double omega2 = 16.1;
        
        // Drive strength
        const double Omega = 0.5;
        
        const double Omega1 = Omega; 
        // We only drive Q0
        const double Omega2 = 0.0;

        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        const std::string hamiltonianJson = createHamiltonianJsonCrossResonance(omega1, Omega1, omega2, Omega2, omegaXX);

        const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
        EXPECT_TRUE(loadOK);

        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 1.0;
        
        // D Channels: drive using cross-resonance freqs.
        const auto omega_d0 = omega2;
        const auto omega_d1 = omega1;
        channelConfigs.loFregs_dChannels.emplace_back(omega_d0/(2*M_PI));
        channelConfigs.loFregs_dChannels.emplace_back(omega_d1/(2*M_PI));
    
        // Add the square pulse
        channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(total_samples));    
        systemModel->setChannelConfigs(channelConfigs);

        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS) });    
        
        // Need 2 qubits for this
        auto qubitReg = xacc::qalloc(2);    

        auto provider = xacc::getIRProvider("quantum");
        auto compositeInst = provider->createComposite("test_cross_resonance" + std::to_string(total_samples));
    
        // Drive on D0 channel
        auto inst1 = std::make_shared<xacc::quantum::Pulse>("square", "d0");
        compositeInst->addInstruction(inst1);

        // Measure both qubits
        auto meas0 = std::make_shared<xacc::quantum::Measure>(0);
        auto meas1 = std::make_shared<xacc::quantum::Measure>(1);
        compositeInst->addInstruction(meas0);
        compositeInst->addInstruction(meas1);
        
        // Run the Pulse simulation with the Hamiltonian provided
        // HACKERY: The Pulse simulator can compute more properties than just occupation expectation,
        // we will just use the existing ExtraInfo data of the buffer to request the visitor to fill in those 
        // results when it completes the simulation.
        // This will just be a map from string (pair of qubit indices) to double (result).
        // We need this, especially for unit testing, to quickly verify the correctness of the simulator.
        std::map<std::string, double> concurrentInfoRequest;
        concurrentInfoRequest.emplace("0_1", 0.0);
        qubitReg->addExtraInfo("Concurrence", concurrentInfoRequest);

        quaC->execute(qubitReg, compositeInst);
        // Debug
        qubitReg->print();
        // Verify the concurrence result for this end time
        auto concurrenceResult = (*qubitReg)["Concurrence"].as<std::map<std::string, double>>();
        const double result = concurrenceResult["0_1"];
        EXPECT_NEAR(result, expectedResults[i], 0.1);
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
