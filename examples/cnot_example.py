import sys
sys.path.insert(1, '/home/cades/dev/Pulse_Control/gym_pulsecontrol/')
import xacc_drl
import numpy as np

# Pulse Parameters:
nbQubits = 2
nbSamples = 512
in_bW = 0.025
in_K = int(2 * nbSamples * in_bW)
T = 600

# Density Matrix for {CNOT} from |00>:
expectedDmReal = np.array([
    1, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0
], dtype = np.float64)
expectedDmImag = np.zeros(16)
# Used for plot titles only:
gate_operation = 'X[q0], Ry[0.59, q1], CNOT'

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

ppo.execute()