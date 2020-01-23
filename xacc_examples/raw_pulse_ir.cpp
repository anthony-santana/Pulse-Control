
#include "xacc.hpp"
#include "Pulse.hpp"

// ============================================================================
// This example demonstrate the use of pulse level IR:
// i.e. uses specifically list out what pulses they want to put on each channel.
// ============================================================================

int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    // Use the mock 2Q backend
    const std::string backendName = "Fake2Q";
   
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("sim-mode", "Pulse"), std::make_pair("backend", backendName), std::make_pair("config-json-path", "/home/cades/dev/xacc/quantum/gate/ir/tests/files/test_backends.json") });    


    // Create a simple pulse program using IR
    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    // Create a some pulse IR's: refer to a pulse (by name) from the library and 
    // specify which channel we want to put this pulse on.
    // Note: pulses on the same channel are assumed to be sequential (one after the other)
    auto pulseInst1 = std::make_shared<xacc::quantum::Pulse>("CR90m_u10_3e8d", "d0");
    auto pulseInst2 = std::make_shared<xacc::quantum::Pulse>("CR90m_u11_77cd", "u0");
    auto pulseInst3 = std::make_shared<xacc::quantum::Pulse>("CR90m_u12_3c09", "d1");
    auto pulseInst4 = std::make_shared<xacc::quantum::Pulse>("CR90m_d7_cefe", "u1");
    auto pulseInst5 = std::make_shared<xacc::quantum::Pulse>("CR90m_d6_6b71", "d0");
    auto pulseInst6 = std::make_shared<xacc::quantum::Pulse>("CR90m_u24_1d67", "d1");
    
    // Add those pulse IR's to the composite and execute.
    compositeInst->addInstruction(pulseInst1);
    compositeInst->addInstruction(pulseInst2);
    compositeInst->addInstruction(pulseInst3);
    compositeInst->addInstruction(pulseInst4);
    compositeInst->addInstruction(pulseInst5);
    compositeInst->addInstruction(pulseInst6);
    
    auto qubitReg = xacc::qalloc(2);    
    quaC->execute(qubitReg, compositeInst);
    
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}