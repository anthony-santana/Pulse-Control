#include "interface_xacc_ir.h"
#include "quac.h"
#include <math.h>
#include "operators.h"
#include "solver.h"
#include "dm_utilities.h"
#include "quantum_gates.h"
#include "petsc.h"
#include <stdbool.h>
#include "macros.h"
#include "PulseControllerHandle.h"

// Global vars 
// Mode of simulation
// Default is circuit mode (quantum gates)
sim_mode g_simulationMode = CIRCUIT;

log_verbosity g_logVerboseLevel = MINIMAL;

void dispatchLog(log_verbosity in_level, const char* in_logFormat, ...) 
{
    // Skip all logs if NONE selected.
    if (g_logVerboseLevel == NONE)
    {
        return;
    }

    const bool shouldLog = (in_level <= g_logVerboseLevel);
    if (shouldLog)
    {
        va_list logArgs;
        va_start(logArgs, in_logFormat);
        vprintf(in_logFormat, logArgs);
        va_end(logArgs);      
    }
}

#define LOG_CRITICAL(format, ...) dispatchLog(MINIMAL, format, ##__VA_ARGS__)
#define LOG_INFO(format, ...) dispatchLog(DEBUG, format, ##__VA_ARGS__)
#define LOG_DEBUG(format, ...) dispatchLog(DEBUG_DIAG, format, ##__VA_ARGS__)

operator* qubits;
Vec psi;
circuit g_circuit;

int g_nbStepCount = 0;
TSData* g_timeSteppingData = NULL;

PulseChannelProvider* g_PulseDataProvider;

// Hard-coded for testing
PetscReal g_timeMax  = 9;
PetscReal g_dt = 0.01;
PetscInt g_stepsMax = 1000;

PetscReal gate_time_step = 1.0;
int nbQubits; 

bool g_wasInitialized = false;

#define ASSERT_QUBIT_INDEX(qubitIdx) \
    if (qubitIdx > nbQubits) \
    { \
        printf("ERROR! Qubit index is out-of-range!\n");\
        exit(1);\
    }

#define ASSERT_PULSE_MODE \
    if (g_simulationMode != PULSE) \
    { \
        printf("ERROR! PULSE simulation mode has not been initialized!\n");\
        exit(1);\
    }

// Generate dummy channel drive functions (signature double(double)) which 
// just refer to g_PulseDataProvider to get the data.
// (each function declared here has an implied channel index)
// Register 2 channels
REGISTER_N_DRIVE_CHANNELS(g_PulseDataProvider, 2);
typedef double channelFunctionType(double time);
// Map from channel name to channel Id
channelFunctionType* GetChannelFunction(const char* in_channelName)
{
    // Simple set-up: 1 qubit with D0 and U0 channels
    if (strcmp(in_channelName, "D0") == 0)
    {
        return GET_DRIVE_CHANNEL_FN(0);
    }
    if (strcmp(in_channelName, "U0") == 0)
    {
        return GET_DRIVE_CHANNEL_FN(1);
    }

    return NULL;
}

// Time-stepping monitor function
PetscErrorCode g_tsDefaultMonitorFunc(TS, PetscInt, PetscReal, Vec, void*);

int XACC_QuaC_Initialize(int in_nbQubit)
{
    // TODO: support passing params (e.g. qubit decay)
    QuaC_initialize(0, NULL);
    nbQubits = in_nbQubit;
    qubits  = malloc(in_nbQubit * sizeof(struct operator));    
    
    for (int i = 0; i < nbQubits; i++)
    {
        create_op(2, &qubits[i]);
    }
    
    // Allocate resources
    create_full_dm(&psi);
    // Initital state
    add_value_to_dm(psi, 0, 0, 1.0);

    // Create circuit
    create_circuit(&g_circuit, -1);
    return 0;
}

int XACC_QuaC_AddInstruction(const char* in_op, const int* in_qbitOperands, int in_qbitOperandCount, int in_argCount, char** in_args)
{
    // Prototype only
    if (strcmp(in_op, "H") == 0)
    {
        if (in_qbitOperandCount == 1)
        {
            add_gate_to_circuit(&g_circuit, 1*gate_time_step, HADAMARD, in_qbitOperands[0]); 
        }
        else
        {
            return -1;
        }        
    }
    if (strcmp(in_op, "CNOT") == 0)
    {
        if (in_qbitOperandCount == 2)
        {
            add_gate_to_circuit(&g_circuit, 1*gate_time_step, CNOT, in_qbitOperands[0], in_qbitOperands[1]); 
        }
        else
        {
            return -1;
        }
    }
    // TODO
    
    // Success
    return 0;
}

const char* XACC_QuaC_ExecuteCircuit(int in_argCount, char** in_args)
{
    start_circuit_at_time(&g_circuit, 0.0);
    time_step(psi, 0.0, g_timeMax, g_dt, g_stepsMax);
    // Debug
    int psi_dims = pow(2, nbQubits);
    print_psi(psi, psi_dims);
    return "";
}

void XACC_QuaC_Finalize()
{
    for (int i=0; i< nbQubits; i++)
    {
        destroy_op(&qubits[i]);
    }
    
    free(qubits);   
    destroy_dm(psi);

    if (g_nbStepCount > 0)
    {
        // Free population data (at each time-step)
        for (int i = 0; i < g_nbStepCount; i++)
        {
            free(g_timeSteppingData[i].populations);
        }
        // Free the array of timestepping data-structures itself.
        free(g_timeSteppingData);
    }

    QuaC_finalize();
}

int XACC_QuaC_InitializePulseSim(int in_nbQubit, double in_dt, double in_stopTime, int in_stepMax, PulseChannelProvider* in_pulseDataProvider)
{
    g_simulationMode = PULSE;
    g_PulseDataProvider = in_pulseDataProvider;
    if (!g_wasInitialized)
    {
        QuaC_initialize(0, NULL);
        g_wasInitialized = true;
    }
   
    nbQubits = in_nbQubit;
    qubits  = malloc(in_nbQubit * sizeof(struct operator));    
    
    for (int i = 0; i < nbQubits; i++)
    {
        create_op(2, &qubits[i]);
        set_initial_pop(qubits[i], 0);
    }
    
    {
        g_dt = in_dt;
        g_timeMax = in_stopTime;
        g_stepsMax = in_stepMax;
    }
   
    return 0;
}

void XACC_QuaC_SetLogVerbosity(log_verbosity in_verboseConfig)
{
    LOG_CRITICAL("Set logging level to %d.\n", in_verboseConfig);
    g_logVerboseLevel = in_verboseConfig;
}

void XACC_QuaC_AddQubitDecay(int in_qubitIdx, double in_kappa)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    add_lin(in_kappa, qubits[in_qubitIdx]);
}

operator* GetQubitOperator(operator qubitOp, const char* in_op)
{
    if (strcmp(in_op, "SM") == 0) {
        // Sigma minus, i.e. itself
        return qubitOp;
    } else if (strcmp(in_op, "SP") == 0) {
        return (qubitOp)->dag;
    } else if (strcmp(in_op, "X") == 0) {
        return (qubitOp)->sig_x;
    } else if (strcmp(in_op, "Y") == 0) {
        return (qubitOp)->sig_y;
    } else if (strcmp(in_op, "Z") == 0) {
        return (qubitOp)->sig_z;
    } else if (strcmp(in_op, "I") == 0) {
        return (qubitOp)->eye;
    } else {
        printf("ERROR! Unknown operator!\n");
        exit(1);
    }
}


void XACC_QuaC_AddConstHamiltonianTerm1(const char* in_op, int in_qubitIdx, ComplexCoefficient in_coeff)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    add_to_ham(in_coeff.real + in_coeff.imag * PETSC_i, GetQubitOperator(qubits[in_qubitIdx], in_op));
}


void XACC_QuaC_AddTimeDependentHamiltonianTerm1(const char* in_op, int in_qubitIdx, const char* in_channelName)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    channelFunctionType* channelFn = GetChannelFunction(in_channelName);
    if (channelFn == NULL)
    {
        printf("ERROR! Unknown drive channel!\n");
        return;
    }
    add_to_ham_time_dep(channelFn, 1, GetQubitOperator(qubits[in_qubitIdx], in_op));
}


void XACC_QuaC_AddConstHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, ComplexCoefficient in_coeff)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx1);
    ASSERT_QUBIT_INDEX(in_qubitIdx2);
    add_to_ham_mult2(in_coeff.real + in_coeff.imag * PETSC_i, GetQubitOperator(qubits[in_qubitIdx1], in_op1), GetQubitOperator(qubits[in_qubitIdx2], in_op2));
}

void XACC_QuaC_AddTimeDependentHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, const char* in_channelName)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx1);
    ASSERT_QUBIT_INDEX(in_qubitIdx2);
    TODO(Implement time-dependent product terms in QuaC)
}

int XACC_QuaC_RunPulseSim(double** out_result, int* out_nbSteps, TSData** out_timeSteppingData)
{
    // Allocate resources
    create_full_dm(&psi);
    set_dm_from_initial_pop(psi);

    set_ts_monitor(g_tsDefaultMonitorFunc);

    g_timeSteppingData =  malloc(g_stepsMax * sizeof(TSData));

    time_step(psi, 0.0, g_timeMax, g_dt, g_stepsMax);
    // Returns the population for each qubit
    int nbResults = get_num_populations();
    *out_result = malloc(nbResults * sizeof(double));
    get_populations(psi, &(*out_result)); 
      
    *out_nbSteps = g_nbStepCount;
    *out_timeSteppingData = g_timeSteppingData;

    return nbResults;
}

PetscErrorCode g_tsDefaultMonitorFunc(TS ts, PetscInt step, PetscReal time, Vec dm, void *ctx)
{
    
    int num_pop = get_num_populations();
    double *populations;
    populations = malloc(num_pop*sizeof(double));
    get_populations(dm, &populations);

    g_timeSteppingData[g_nbStepCount].time = time;
    g_timeSteppingData[g_nbStepCount].nbPops = num_pop;
    g_timeSteppingData[g_nbStepCount].populations = malloc(num_pop * sizeof(double));
    memcpy(g_timeSteppingData[g_nbStepCount].populations, populations, num_pop * sizeof(double));

    if (nid==0)
    {
        // Debug:
        LOG_DEBUG(">> t = %e: ", time);
        for(int i = 0; i < num_pop; i++)
        {
            LOG_DEBUG("%e ", populations[i]);
        }
        LOG_DEBUG("\n");
    }

    free(populations);
    g_nbStepCount++;
    PetscFunctionReturn(0);
}
