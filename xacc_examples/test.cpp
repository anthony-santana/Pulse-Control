// !!! TEMP CODE !!! Thien Nguyen: For testing purposes only

#include "xacc.hpp"

int main (int argc, char** argv) {

	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);

    auto quaC = xacc::getAccelerator("QuaC");
    
    auto qubitReg = xacc::qalloc(2);

	// Create a Program
	auto xasmCompiler = xacc::getCompiler("xasm");
    auto ir = xasmCompiler->compile(R"(__qpu__ void test(qbit q) {
        H(q[0]);
        CX(q[0], q[1]);
    })", quaC);
    
    auto program = ir->getComposite("test");
    quaC->execute(qubitReg, program);
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}