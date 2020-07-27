import sys
sys.path.insert(1, '/home/cades/dev/Pulse_Control/gym_pulsecontrol/')
import xacc_drl
import numpy as np

# Pulse Parameters:
nbQubits = 2
nbSamples = 512
in_bW = 0.02
in_K = 5 # int(2 * env.nbSamples * env.in_bW)
T = 60 + 60 + 500

# Density Matrix for {X[q0], Ry[0.59, q1], CNOT}
expectedDmReal = np.array([
    0, 0, 0, 0,
    0, 0, 0, 0, 
    0, 0, 0.91547034, -0.27818051,
    0, 0, -0.27818051, 0.08452966
], dtype = np.float64)
expectedDmImag = np.zeros(16)
# Used for plot titles only:
gate_operation = 'X[q0], Ry[0.59, q1], CNOT'

# Initialize the Proximal Policy Optimization Module
ppo = xacc_drl.OptimalControl(
    {
    'nbQubits': nbQubits,
    'backend': 'two_qubit',
    'slepian_parameters': [nbSamples, in_bW, in_K, T],
    'expectedDmReal': expectedDmReal,
    'expectedDmImag': expectedDmImag,
    'gate_operation': gate_operation,
    'initial_state': [0, 0]
    }
)

# Execute the algorithm
ppo.execute()







