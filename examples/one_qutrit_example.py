import sys
sys.path.insert(1, '/home/cades/dev/Pulse_Control/gym_pulsecontrol/')
import xacc_drl
import numpy as np

# Pulse Parameters:
nbQubits = 1
nbSamples = 512
in_bW = 0.025
in_K = int(2 * nbSamples * in_bW)
T = 100

# Used for plot titles only:
gate_operation = 'X-Gate'

ppo = xacc_drl.OptimalControl(
    {
    'nbQubits': nbQubits,
    'backend': 'one_qutrit',
    'slepian_parameters': [nbSamples, in_bW, in_K, T],
    'gate_operation': gate_operation,
    }
)

ppo.execute()