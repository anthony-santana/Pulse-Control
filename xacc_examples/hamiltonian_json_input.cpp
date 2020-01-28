
#include "xacc.hpp"
#include "Pulse.hpp"

int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    // Use the mock 2Q backend
    const std::string backendName = "Fake2Q";
   
    // An example two-qubit Hamiltonian (adapted from Qiskit AER pulse_sim)
    const std::string hamiltonianJson = R"(
        {
            "description": "Two-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["pi*(2*v0-alpha0)*O0", "pi*alpha0*O0*O0", "2*pi*r*X0||D0", "2*pi*r*X0||U1", "2*pi*r*X1||U0", "pi*(2*v1-alpha1)*O1", "pi*alpha1*O1*O1", "2*pi*r*X0||D1", "2*pi*j*SP0*SM1", "2*pi*j*SM0*SP1"],
            "osc": {},
            "qub": {
                "0": 2,
                "1": 2
            },
            "vars": {
                "v0": 5.00, 
                "v1": 5.1, 
                "j": 0.01, 
                "r": 0.02, 
                "alpha0": -0.33, 
                "alpha1": -0.33
            }
        }
    )";

    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("backend", backendName), std::make_pair("hamiltonian", hamiltonianJson) });    

    auto qubitReg = xacc::qalloc(2);    


    auto randomPulse = std::make_shared<xacc::quantum::Pulse>("myPulse", "d0");
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

    randomPulse->setSamples(samples);
    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    compositeInst->addInstruction(randomPulse);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);

    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}