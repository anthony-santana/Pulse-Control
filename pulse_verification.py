import sys
from numpy import genfromtxt
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import os
import xacc
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import spectrum

hamiltonianJson = {
    "description": "One-qutrit Hamiltonian.",
    "h_latex": "",
    "h_str": ["(w - 0.5*alpha)*O0", "0.5*alpha*O0*O0", "O*(SM0 + SP0)||D0"],
    "osc": {},
    "qub": {
        "0": 3
    },
    "vars": {
        "w": 31.63772297724,
        "alpha": -1.47969,
        "O": 0.0314
    }
}

# Total time, T, of control pulse
T = 100
# Number of pulse samples
nbSamples = 200

model = xacc.createPulseModel()
# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))
qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024, 'logging-period': 0.1 })
channelConfig = xacc.BackendChannelConfigs()
# Setting resolution of pulse
channelConfig.dt = nbSamples / T 
# Driving on resonance with qubit
channelConfig.loFregs_dChannels = [5.0353]
model.setChannelConfigs(channelConfig)

time_steps = np.arange(nbSamples)
noise_signal = [np.cos(1.47969 * time_steps[i]) for i in range(nbSamples)]
pulseData = genfromtxt('output_files/optimal_pulse549.csv', delimiter=',') + noise_signal
# Add that slepian pulse instruction to XACC
pulseName = 'Slepian' 
print(pulseName)
xacc.addPulse(pulseName, pulseData)   
q = xacc.qalloc(1)
# Create the quantum program that contains the slepian pulse
# and the drive channel (D0) is set on the instruction
provider = xacc.getIRProvider('quantum')
prog = provider.createComposite('pulse')
slepianPulse = provider.createInstruction(pulseName, [0])
slepianPulse.setChannel('d0')
prog.addInstruction(slepianPulse)
# Measure Q0 (using the number of shots that was specified above)
prog.addInstruction(xacc.gate.create("Measure", [0]))

# Run the simulation
qpu.execute(q, prog)
resultProb = q['DensityMatrixDiags'][1]
print(resultProb)

# Do the state population vs. time graphs for a pulse optimized on a noiseless system and then look at how it handles the noisy system