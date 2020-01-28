#include "xacc.hpp"
#include "Pulse.hpp"
#include "xacc_service.hpp"


// ============================================================================
// This example demonstrate the use of commnand level IR (in Pulse mode):
// i.e. a command (e.g. X, CX) is translated to a sequence of pulses
// and simulate dynamically. 
// ============================================================================
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    // Use the mock 1Q backend
    const std::string backendName = "Fake1Q";
   
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("backend", backendName), std::make_pair("config-json-path", "/home/cades/dev/xacc/quantum/gate/ir/tests/files/test_backends.json") });    

    // Create a simple pulse program using IR
    auto provider = xacc::getIRProvider("quantum");
    auto program = provider->createComposite("test_program");
    auto x0 = xacc::getContributedService<xacc::Instruction>("pulse::x_0");
    program->addInstructions({ x0 });   

    auto qubitReg = xacc::qalloc(1);    
    quaC->execute(qubitReg, program);
    
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}