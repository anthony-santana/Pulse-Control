## IMPORTANT: need to copy the helper functions in examples/aer_vqe_pulse.py
## (to compute the expectation value from the qutrit state vector)
import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc, json, numpy as np
from scipy.linalg import block_diag
# Pauli operators in an extended basis (e.g. qutrit)
def s(n_qubits, dim, site, pauli):
    """
    Parameters:
    -----------
    n_qubits: number of qubits \in (1,2,..)
    dim: dimension 
    site: site for pauli \in (0,1,..)
    pauli: 'X','Y','Z'
    
    Returns:
    --------
    single-site pauli matrix defined over tensor product hilbert space (arbitray dimension)
    """
    oplist = [np.eye(dim)]*n_qubits
    if pauli == 'X': spin_op = block_diag([[0,1],[1,0]],0*np.eye(dim-2))
    elif pauli == 'Y': spin_op = block_diag([[0,-1j],[1j,0]],0*np.eye(dim-2))
    elif pauli == 'Z': spin_op = block_diag([[1,0],[0,-1]],0*np.eye(dim-2))
    else: return TypeError
    oplist[site] = spin_op
    result = oplist[0]
    for i in range(1, n_qubits):
        result = np.kron(oplist[i], result)
    return result

def pauli_op_to_matrix(obs, n_qubits, dim = 3):
    """
    Parameters:
    -----------
    obs: xacc observable (Pauli type)
    n_qubits: number of qubits of the total system. Should be equal to the number of qubits on the device. 
    dim: qubit dimension
    
    Returns:
    --------
    The observable matrix in the extended space (qutrit)
    """
    def getTermCoeff(term):
        if term is None:
            return 0.0
        for op in term:
            coeff = op[1].coeff()
            return coeff
    
    result = getTermCoeff(obs.getIdentitySubTerm()) * np.eye(dim**n_qubits, dtype=complex)
    for term in obs.getNonIdentitySubTerms():
        coeff = getTermCoeff(term)
        termMat = coeff * np.eye(dim**n_qubits)
        for op in term:
            for idx in op[1].ops():
                opName = op[1].ops()[idx]
                termMat = termMat @ s(n_qubits, dim, idx, opName)
            result += termMat
    return result

def calculateExpVal(state_vec, observable):
    """
    Parameters:
    -----------
    state_vec: state vector (1D list/array)
    observable: observable matrix
    
    Returns:
    --------
    expectation value (real)
    """
    state_vec = np.array(state_vec)
    # Check dimension
    if (len(state_vec) != observable.shape[0] or len(state_vec) != observable.shape[1]):
        return ValueError
    temp = observable @ state_vec
    return np.real(state_vec.conj().T @ temp)

# Retrieve the state vector from accelerator buffer result
def getStateVectorFromBuffer(buffer):
    state_vec = buffer["state"]
    result = []
    for pair in state_vec:
        if (len(pair) != 2):
            return TypeError
        else:
            result.append(pair[0] + 1j * pair[1])
    return result

# Validate pulse-level simulation of VQE protocol:
# Note: this is just a way to validate the pulse simulator,
# the pulse program is just a direct lowering of the ansatz circuit.

# Using the pulse simulator (ibmq_armonk)
qpu = xacc.getAccelerator("aer:ibmq_armonk", {"sim-type": "pulse"})
# Construct the Hamiltonian as an XACC-VQE PauliOperator
ham = xacc.getObservable('pauli', '0.70710678118 X0 + 0.70710678118 Z0')

xacc.qasm('''.compiler xasm
.circuit ansatz
.qbit q
.parameters t0, t1, t2
U(q[0], t0, t1, t2);
''')

ansatz = xacc.getCompiled('ansatz')
opt = xacc.getOptimizer('nlopt')
backend_properties = qpu.getProperties()
config = json.loads(backend_properties["config-json"])
# Number of qubits of the backend
nbQubits = config["n_qubits"]

def pulse_opt_func(x):
    buffer = xacc.qalloc(nbQubits)
    program = ansatz.eval(x)
    qpu.execute(buffer, program)
    state_vec = getStateVectorFromBuffer(buffer)
    obs = pauli_op_to_matrix(ham, buffer.size())
    energy = calculateExpVal(state_vec, obs)
    print("E({0}) = {1}".format(x, energy))
    return energy

optimizer = xacc.getOptimizer('nlopt')
result = optimizer.optimize(pulse_opt_func, 3)
print("Optimize:", result)