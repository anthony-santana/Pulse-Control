import sys
sys.path.insert(1, '/home/cades/dev/Pulse_Control/gym_pulsecontrol/')
import xacc_drl
import numpy as np

# Pulse Parameters:
nbQubits = 2
nbSamples = 512
in_bW = 0.02
in_K = 5 #int(2 * nbSamples * in_bW)
T = 50 + 50 + 500 # 130 + 130 + 650 

# Density Matrix for {X[q0], Ry[0.59, q1], CNOT}
expectedDmReal = np.array([
    0, 0, 0, 0,
    0, 0.08572349, 0.27995531, 0, 
    0, 0.27995531, 0.91427651, 0,
    0, 0, 0, 0
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
    'initial_state': [0, 0],
    'channels': ['d0', 'd1'],
    'nbPulses': 2
    }
)

# Execute the algorithm
ppo.execute()







