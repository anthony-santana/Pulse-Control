#include "xacc.hpp"
#include "Pulse.hpp"

// ============================================================================
// This example demonstrate the use of pulse level IR:
// We use pulses from Q-CTRL (https://q-ctrl.com/) Black Opal library to perform 
// X gates (rather than the default Gaussian pulse)
// ============================================================================

int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    // Use our fake 1Q system backend
    // TODO: eventually, we need to be able to specify a backend here (with its Hamiltonian, device dt, and LO freqs, etc.)
    const std::string backendName = "Fake1Q";
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("backend", backendName) });    


    // Create a simple pulse program using IR
    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    // Note: we need to scale the pulse to match out Fake1Q backend
    const auto scalePulseSamples = [](std::vector<std::vector<double>>& io_samples){
        const double scaleFactor = 0.725;
        for (auto& elementVec : io_samples)
        {
            for (auto& element: elementVec)
            {
                element *= scaleFactor;
            }
        }
    };
   
    // Walsh-Gaussian pulse from Black Opal library
    auto pulseWalshGaussian = std::make_shared<xacc::quantum::Pulse>("BlackOpalWalshGaussian", "d0");
    {
        // Pulse reference: https://arxiv.org/abs/1809.03452 (sec 7.1)
        std::vector<std::vector<double>> samples = {
            {0.0, 0.0}, {0.013434, 0.0}, {0.058597, 0.0}, {0.146446, 0.0},  
            {0.229930, 0.0}, {0.229930, 0.0}, {0.1464464, 0.0},  
            {0.058597, 0.0}, {0.013434, 0.0}, {0.0, 0.0}, {0.0, 0.0},  
            {0.009035, 0.0}, {0.039411, 0.0}, {0.098498, 0.0},  
            {0.154648, 0.0}, {0.154648, 0.0}, {0.098498, 0.0},  
            {0.039411, 0.0}, {0.009035, 0.0}, {0.0, 0.0}, {0.0, 0.0},  
            {0.009035, 0.0}, {0.039411, 0.0}, {0.098498, 0.0},  
            {0.154648, 0.0}, {0.154648, 0.0}, {0.098498, 0.0},  
            {0.039411, 0.0}, {0.009035, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 
            {0.013434, 0.0}, {0.058597, 0.0}, {0.146446, 0.0},  
            {0.229930, 0.0}, {0.229930, 0.0}, {0.146446, 0.0},  
            {0.058597, 0.0}, {0.013434, 0.0}, {0.0, 0.0}
        };   

        scalePulseSamples(samples);
        pulseWalshGaussian->setSamples(samples);
    }
    
    //  BB1 Gaussian pulse from Black Opal library
    auto pulseBB1Gaussian = std::make_shared<xacc::quantum::Pulse>("BlackOpalBB1Gaussian", "d0");
    {
        // Pulse reference: https://arxiv.org/abs/1809.03452 (sec 7.1)
        std::vector<std::vector<double>> samples = {
            {0.0, 0.0}, {0.014980, 0.0}, {0.065339, 0.0}, {0.163296, 0.0}, 
            {0.256385, 0.0}, {0.256385, 0.0}, {0.163296, 0.0}, 
            {0.065339, 0.0}, {0.014980, 0.0}, {0.0, 0.0}, {0.0, 0.0},  
            {-0.003745, 0.014504}, {-0.016335, 0.063264}, 
            {-0.040824, 0.158111}, {-0.064096, 0.248244},  
            {-0.064096, 0.248244}, {-0.040824, 0.158111}, 
            {-0.016335, 0.063264}, {-0.003745, 0.014504}, {0.0, 0.0}, 
            {0.0, 0.0}, {0.020597, -0.021756}, {0.089841, -0.094896}, 
            {0.224532, -0.237166}, {0.352530, -0.372366},  
            {0.352530, -0.372366}, {0.224532, -0.237166}, 
            {0.089841, -0.094896}, {0.020597, -0.021756}, {0.0, 0.0}, 
            {0.0, 0.0}, {-0.003745, 0.014504}, {-0.016335, 0.063264},  
            {-0.040824, 0.158111}, {-0.064096, 0.248244}, 
            {-0.064096, 0.248244}, {-0.040824, 0.158111},  
            {-0.016335, 0.063264}, {-0.003745, 0.014504}, {0.0, 0.0}
        };

        scalePulseSamples(samples);
        pulseBB1Gaussian->setSamples(samples);
    }
    
    // Add the pulse and execute (use one of the two)
    compositeInst->addInstruction(pulseWalshGaussian);
    // compositeInst->addInstruction(pulseBB1Gaussian);

    auto qubitReg = xacc::qalloc(1);    
    quaC->execute(qubitReg, compositeInst);
    
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}