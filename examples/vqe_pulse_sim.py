## IMPORTANT: need to copy the helper functions in examples/aer_vqe_pulse.py
## (to compute the expectation value from the qutrit state vector)

import xacc, json, numpy as np

# Validate pulse-level simulation of VQE protocol:
# Note: this is just a way to validate the pulse simulator,
# the pulse program is just a direct lowering of the ansatz circuit.

# Using the pulse simulator (ibmq_rome)
qpu = xacc.getAccelerator("aer:ibmq_rome", {"sim-type": "pulse"})
# Construct the Hamiltonian as an XACC-VQE PauliOperator
ham = xacc.getObservable('pauli', '5.907 - 2.1433 X0X1 - 2.1433 Y0Y1 + .21829 Z0 - 6.125 Z1')

xacc.qasm('''.compiler xasm
.circuit ansatz
.qbit q
.parameters t0
X(q[0]);
Ry(q[1], t0);
CNOT(q[1], q[0]);
''')

ansatz = xacc.getCompiled('ansatz')
opt = xacc.getOptimizer('nlopt')
backend_properties = qpu.getProperties()
config = json.loads(backend_properties["config-json"])
# Number of qubits of the backend
nbQubits = config["n_qubits"]

# Just sweep the angles for testing
for a in np.linspace(-np.pi, np.pi, 10):
    buffer = xacc.qalloc(nbQubits)
    program = ansatz.eval([a])
    qpu.execute(buffer, program)
    state_vec = getStateVectorFromBuffer(buffer)
    obs = pauli_op_to_matrix(ham, buffer.size())
    energy = calculateExpVal(state_vec, obs)
    print("E({0}) = {1}".format(a, energy))
