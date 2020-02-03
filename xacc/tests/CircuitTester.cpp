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

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
