#include "interface_xacc_ir.h"
#include "quac.h"
#include "operators.h"
#include "solver.h"
#include "dm_utilities.h"
#include "quantum_gates.h"
#include "petsc.h"

// Global vars 
operator* qubits;
Vec psi;
circuit g_circuit;
// Hard-coded for testing
PetscReal time_max  = 9;
PetscReal dt = 0.01;
PetscInt steps_max = 1000;
PetscReal gate_time_step = 1.0;
int nbQubits; 

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
    time_step(psi, 0.0, time_max, dt, steps_max);
    // Debug
    int psi_dims = pow(2, nbQubits);
    print_dm(psi, psi_dims);
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
    QuaC_finalize();
}
