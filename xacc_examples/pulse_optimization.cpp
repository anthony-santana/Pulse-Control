#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include <cassert>
#include "CommonGates.hpp"
#include "Optimizer.hpp"

// =======================================================================================
// This examples demonstrates single-qubit pulse optimization using 
// existing facilities from the XACC framwork.
// Introduction/Motivation:
// To realize single-qubit rotations, we need to specify a certain pulse shape 
// for driving channels. The areas under these functions determine the angle of rotation.
// Pulses can be shaped to reduce unwanted effects such as, e.g., leakage due to weak anharmonicity. 
// XACC-QuaC pulse-level simulator provides an emulator backend for pulse optimization design
// e.g. optimizing the design itself before deploying/testing on physical quantum hardware. 
// =======================================================================================

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
    
    // Create a transmon qubit Hamiltonian:
    // - omega: qubit transition frequency 
    // - driving strength: the coefficient to be multiplied with the driving channel.
    // - anharmonicity (delta) of the transmon qubit
    // - dimension cutoff of the qubit 
    std::string createTransmonHamiltonianJson(double in_omega0, double in_drivingStrength, double in_anharmonicity, int in_dim = 2)
    {
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string gPlaceHolder = "{{g}}";
        const std::string deltaPlaceHolder = "{{delta}}";
        const std::string dimPlaceHolder = "{{dim}}";

        std::string hamiltonianJson = singleQubitHamiltonianJsonTemplate;
        // Substitute the value to the Hamiltonian
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(in_omega0));
        hamiltonianJson.replace(hamiltonianJson.find(gPlaceHolder), gPlaceHolder.length(), std::to_string(in_drivingStrength));
        hamiltonianJson.replace(hamiltonianJson.find(deltaPlaceHolder), deltaPlaceHolder.length(), std::to_string(in_anharmonicity));
        hamiltonianJson.replace(hamiltonianJson.find(dimPlaceHolder), dimPlaceHolder.length(), std::to_string(in_dim));
        return hamiltonianJson;
    }
}


int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    // Hamiltonian of the transmon qubit system
    // Transmon qubit transition freq.
    const double omega0 = 5.35*2*M_PI; // 5.35 GHz
    
    // Anharmonicity
    const double delta = 0.35*2*M_PI; // 350 MHz
    const int total_samples = 100;
    
    // Driving field coefficient
    const double g = 2.0*M_PI/total_samples;

    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const std::string hamiltonianJson = createTransmonHamiltonianJson(omega0, g, delta);
    const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0; //ns
    // D0 Channel: drive on resonance
    const auto omega_d0 = omega0;
    channelConfigs.loFregs_dChannels.emplace_back(omega_d0/(2*M_PI));

    xacc::OptFunction f(
        [&](const std::vector<double>& x, std::vector<double>& dx) {
            // Problem statement: 
            // We will try to select an optimal Gaussian pulse width to implement an
            // X(pi/2) gate which produces a 50/50 probability distribution.
            // Note: X(pi/2) is used to construct H gate by an additional frame change.
            // For this simple optimization example, we just optimize a single param that is the the pulse width.            
            if (x[0] <= 0.0 || x[0] > 1.0)
            {
                // The param here is expected to be within 0.0 and 1.0
                return 999999.999;
            }

            const double  gaussSigma = x[0] * total_samples * channelConfigs.dt; 
            // Add the Gaussian pulse
            channelConfigs.addOrReplacePulse("gaussian", QuaC::GaussianPulse(total_samples, gaussSigma));    
            systemModel->setChannelConfigs(channelConfigs);
            
            // We run 10000 shots
            const int NB_SHOTS = 10000;
            auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

            auto qubitReg = xacc::qalloc(1);    

            auto provider = xacc::getIRProvider("quantum");
            auto compositeInst = provider->createComposite("test_gaussian_pulse");
            auto pulseInst = std::make_shared<xacc::quantum::Pulse>("gaussian", "d0");
            // Add the Gaussian pulse
            compositeInst->addInstruction(pulseInst);
            
            // Mesure qubit after the pulse
            auto meas = std::make_shared<xacc::quantum::Measure>(0);
            compositeInst->addInstruction(meas);

            // Run the Pulse simulation with the Hamiltonian provided
            quaC->execute(qubitReg, compositeInst);
            // Debug:
            // qubitReg->print();
            
            const double prob0 = qubitReg->computeMeasurementProbability("0");
            const double prob1 = qubitReg->computeMeasurementProbability("1");
            // The optimizer will try to minimize the different between '0' and '1' probabilities
            // i.e. to achieve a perfect 50/50 distribution.
            // Due to random sampling error, we just zero out small differences so that it converges a bit faster.
            // Note: we can easily extend this *Cost function* to include things like gate time, drive strength constrainsts, etc.
            // by adding drive amplitude as an optimization parameter as well.
            const auto probDiff = std::abs(prob0 - prob1);
            return probDiff < 0.01 ? 0.0 : probDiff;
        },
        // We only have one variable
        1
    );
    
    // Get the XACC nlopt service
    // For Pulse optimization, normally, we should already have a guess about the initial values, hence should feed it to the optimizer.
    auto optimizer = xacc::getOptimizer("nlopt", { std::make_pair("nlopt-maxeval", 20), std::make_pair("initial-parameters", std::vector<double>{ 0.25 }) });
    auto result = optimizer->optimize(f);
    
    std::cout << "Optimal Gaussian pulse width (sigma): " << result.second[0] * total_samples * channelConfigs.dt << "ns.\n";

    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}