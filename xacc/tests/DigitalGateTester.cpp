#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include "CommonGates.hpp"

namespace {
    const std::string singleQubitHamiltonianJsonTemplate = R"#(
        {
            "description": "One-qubit transmon qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["omega0*O0", "g*(SP0 + SM0)||D0", "0.5*delta*O0*(I0-O0)"],
            "osc": {},
            "qub": {
                "0": {{dim}}
            },
            "vars": {
                "omega0": {{omega0}},
                "g": {{g}},
                "delta": {{delta}}
            }
        }
    )#";
    
    std::string createTransmonHamiltonianJson()
    {
        const double omega0 = 5.35*2*M_PI; 
        const double delta = 0.35*2*M_PI; 
        // Driving field coefficient
        const double g = 2.0*M_PI/100;        
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string gPlaceHolder = "{{g}}";
        const std::string deltaPlaceHolder = "{{delta}}";
        const std::string dimPlaceHolder = "{{dim}}";

        std::string hamiltonianJson = singleQubitHamiltonianJsonTemplate;
        // Substitute the value to the Hamiltonian
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(omega0));
        hamiltonianJson.replace(hamiltonianJson.find(gPlaceHolder), gPlaceHolder.length(), std::to_string(g));
        hamiltonianJson.replace(hamiltonianJson.find(deltaPlaceHolder), deltaPlaceHolder.length(), std::to_string(delta));
        hamiltonianJson.replace(hamiltonianJson.find(dimPlaceHolder), dimPlaceHolder.length(), "2");
        return hamiltonianJson;
    }
}


TEST(DigitalGateTester, testHadamard)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    
    // We run 10000 shots
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto qubitReg = xacc::qalloc(1);    
    auto xasmCompiler = xacc::getCompiler("xasm");
    auto ir = xasmCompiler->compile(R"(__qpu__ void test1(qbit q) {
      H(q[0]);
      Measure(q[0]);
    })", quaC);

    auto program = ir->getComposite("test1");
    quaC->execute(qubitReg, program);

    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    EXPECT_NEAR(prob0, 0.5, 0.1);
    EXPECT_NEAR(prob1, 0.5, 0.1);
}


TEST(DigitalGateTester, testRx)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto xasmCompiler = xacc::getCompiler("xasm");
    
    // Rotate Rx(theta) => expected result: cos(theta/2) |0> + sin(theta/2) |1>
    auto ir = xasmCompiler->compile(R"(__qpu__ void test2(qbit q, double t) {
      Rx(q[0], t);
      Measure(q[0]);
    })", quaC);
    auto program = ir->getComposite("test2");
    
    const auto angles = xacc::linspace(-xacc::constants::pi, xacc::constants::pi, 20);
    for (const auto& angle : angles)
    {
        auto qubitReg = xacc::qalloc(1);    
        auto evaled = program->operator()({ angle });
        quaC->execute(qubitReg, evaled);
        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        // cos(theta/2) |0> + sin(theta/2) |1>
        const double expectedProb0 = cos(angle/2)*cos(angle/2);
        const double expectedProb1 = 1.0 - expectedProb0;
        EXPECT_NEAR(prob0, expectedProb0, 0.1);
        EXPECT_NEAR(prob1, expectedProb1, 0.1);
    }
}

TEST(DigitalGateTester, testMultipleGates)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto xasmCompiler = xacc::getCompiler("xasm");
    
    auto ir = xasmCompiler->compile(R"(__qpu__ void test3(qbit q) {
      H(q[0]);
      Z(q[0]);
      H(q[0]);
      Measure(q[0]);
    })", quaC);
    auto program = ir->getComposite("test3");
    auto qubitReg = xacc::qalloc(1);    
    quaC->execute(qubitReg, program);

    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    // Should be in |1> state
    EXPECT_NEAR(prob0, 0.0, 0.01);
    EXPECT_NEAR(prob1, 1.0, 0.01);
}

TEST(DigitalGateTester, testMixPulseDigital)
{
    // Similar to the above test: H-Z-H is equiv to an X gate,
    // but use pulses for the Hadamard gate instead.
    
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);   
    const int total_samples = 100; 
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(5.35);
    
    // Optimal pulse for X(pi/2) of this system    
    const double  gaussSigma = 27.7612; //ns 
    // Add the Gaussian pulse
    channelConfigs.addOrReplacePulse("gaussian", QuaC::GaussianPulse(total_samples, gaussSigma));    
    systemModel->setChannelConfigs(channelConfigs);
    
    // We run 10000 shots
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto qubitReg = xacc::qalloc(1);    
    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_hadamard_gaussian_pulse");
    
    // First Hadamard (Pulse)
    {
        // We need to frame change in order to create a Hadarmard
        {
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
        {
            // Add X(pi/2) pulse (Gaussian with optimal pulse width)
            auto pulseInst = std::make_shared<xacc::quantum::Pulse>("gaussian", "d0");
            pulseInst->setBits({0});
            // Add the Gaussian pulse
            compositeInst->addInstruction(pulseInst);
        }
        {
            // Another frame-change
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
    }
        
    // We should have applied a Hadamard gate by *pulse*
    // Now apply a digital Z gate:
    {
        auto zGate = std::make_shared<xacc::quantum::Z>(0);
        compositeInst->addInstruction(zGate);
    }
    
    // Next: another Hadamard gate: by pulse
    {
        {
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
        {
            // Add X(pi/2) pulse (Gaussian with optimal pulse width)
            auto pulseInst = std::make_shared<xacc::quantum::Pulse>("gaussian", "d0");
            pulseInst->setBits({0});
            // Add the Gaussian pulse
            compositeInst->addInstruction(pulseInst);
        }
        {
            // Another frame-change
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
    }
   
    // Measure qubit after the pulse
    auto meas = std::make_shared<xacc::quantum::Measure>(0);
    compositeInst->addInstruction(meas);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    
    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    // Should be in |1> state
    EXPECT_NEAR(prob0, 0.0, 0.1);
    EXPECT_NEAR(prob1, 1.0, 0.1);  
}

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}