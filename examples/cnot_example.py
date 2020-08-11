import sys
sys.path.insert(1, '/home/cades/dev/Pulse_Control/gym_pulsecontrol/')
import xacc_drl
import numpy as np

# Pulse Parameters:
nbQubits = 2
nbPulses = 1
nbSamples = 512
in_bW = 0.02
in_K = 5 # int(2 * nbSamples * in_bW)
T = 500

# Density Matrix for {CNOT} from |11>:
expectedDmReal = np.array([
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 0
], dtype = np.float64)
expectedDmImag = np.zeros(16)
# Used for plot titles only:
gate_operation = 'CNOT'

ppo = xacc_drl.OptimalControl(
    {
    'nbQubits': nbQubits,
    'nbPulses': nbPulses,
    'backend': 'two_qubit',
    'slepian_parameters': [nbSamples, in_bW, in_K, T],
    'expectedDmReal': expectedDmReal,
    'expectedDmImag': expectedDmImag,
    'gate_operation': gate_operation,
    'initial_state': [1, 1],
    'channels': ['d0'],
    'nbPulses': 1,
    'nsteps': 5 * in_K,
    'nminibatches': in_K
    }
)

ppo.execute()