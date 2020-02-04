#pragma once

#define DO_PRAGMA(x) _Pragma(#x)

#ifndef TODO
#define TODO(x) DO_PRAGMA(message("\033[30;43mTODO\033[0m - " #x))
#endif


// Interface to XACC IR: 
// Simulation mode: Circuit (gates) or Pulse (Hamiltonian)
typedef enum {
    CIRCUIT  = 0,
    PULSE = 1
} sim_mode;

typedef enum {
    NONE  = 0,
    MINIMAL = 1, // Important logging
    DEBUG = 2, // More logging
    DEBUG_DIAG = 3 // Very verbose
} log_verbosity;

typedef struct ComplexCoefficient {
    double real;
    double imag;
} ComplexCoefficient;

// Time-stepping data
typedef struct TSData {
    double time;
    int nbChannels;
    double* channelData;
    int nbPops;
    double* populations;
    double* pauliExpectations;
} TSData;

typedef struct XaccPulseChannelProvider PulseChannelProvider;

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

// Pulse simulation initialization:
// Note: we *solve* the master equation using QuaC, not via Monte-Carlo method.
// Hence, we don't need to have the *shots* params.
__attribute__ ((visibility ("default"))) extern int XACC_QuaC_InitializePulseSim(int in_nbQubit, PulseChannelProvider* in_pulseDataProvider, const int* in_qbitDims);


__attribute__ ((visibility ("default"))) extern void XACC_QuaC_SetLogVerbosity(log_verbosity in_verboseConfig);


// Add a decay term (Lindblad)
__attribute__ ((visibility ("default"))) extern void XACC_QuaC_AddQubitDecay(int in_qubitIdx, double in_kappa);

// Run the Pulse simulation and return the expectation values:
// Returns the size of the result array. Caller needs to clean up. 
__attribute__ ((visibility ("default"))) extern int XACC_QuaC_RunPulseSim(double in_dt, double in_stopTime, int in_stepMax, double** out_result, int* out_nbSteps, TSData** out_timeSteppingData);

// ====   Hamiltonian construction API's ====
// Adding a single-operator term to the Hamiltonian:
// (1) Time-independent term: 
// Syntax: coeff * ['X', 'Y', 'Z', 'I', 'SP', 'SM', 'O', 'N']_i
// 'SP' and 'SM' are the sigma plus and sigma minus operators.
// Coefficient is a complex parameter and this term can only act on 1 qubit.
__attribute__ ((visibility ("default"))) extern void XACC_QuaC_AddConstHamiltonianTerm1(const char* in_op, int in_qubitIdx, ComplexCoefficient in_coeff);
// (2) Time-dependent term:
// Similar to (1) but has a time-dependent drive function (double -> double)
// Note: drive signal must have been *mixed* with LO, i.e. it is Re[d(t) * exp(-i * w_LO * t)] = d(t) * cos(w_LO * t)
__attribute__ ((visibility ("default"))) extern void XACC_QuaC_AddTimeDependentHamiltonianTerm1(const char* in_op, int in_qubitIdx, int in_channelId, double in_coefficient);


// Adding a two-operator term to the Hamiltonian:
// (1) Time-independent term: 
__attribute__ ((visibility ("default"))) extern void XACC_QuaC_AddConstHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, ComplexCoefficient in_coeff);
// (2) Time-dependent term:
__attribute__ ((visibility ("default"))) extern void XACC_QuaC_AddTimeDependentHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, int in_channelId, double in_coefficient);

__attribute__ ((visibility ("default"))) extern void XACC_QuaC_SetInitialPopulation(int in_qubitIdx, double in_initialPopulation);

__attribute__ ((visibility ("default"))) extern void XACC_QuaC_DisableAdaptiveTimestepping();


// Some methods to compute useful properties from the density matrix, e.g. to append as extra info to the XACC accelerator buffer
// Calculate the concurrence between 2 qubits
__attribute__ ((visibility ("default"))) extern double XACC_QuaC_CalcConcurrence(int in_qubitIdx1, int in_qubitIdx2);

// ======================================================================