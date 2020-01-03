#pragma once

// Interface to XACC IR: 

// Circuit mode initialization:
// TODO: define more options
__attribute__ ((visibility ("default"))) extern int XACC_QuaC_Initialize(int in_nbQubit);

// Add an IR instruction to the current QuaC circuit
// in_op IR operation (as string) we wish to add.
// param args The operands, if any, for the operation.
// return 0 Success - operatation is added; 
// otherwise (>0) Failure - operatation cannot be added
__attribute__ ((visibility ("default"))) extern int XACC_QuaC_AddInstruction(const char* in_op, const int* in_qbitOperands, int in_qbitOperandCount, int in_argCount, char** in_args);

// Execute the circuit and collect data specified by input params
// Returns: JSON-encoded data of the result.
__attribute__ ((visibility ("default"))) extern const char* XACC_QuaC_ExecuteCircuit(int in_argCount, char** in_args);

// Clean-up any allocated resources
__attribute__ ((visibility ("default"))) extern void XACC_QuaC_Finalize();
