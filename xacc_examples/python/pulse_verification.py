import sys
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

# The Hamiltonian JSON object (OpenPulse format)
# omega0 = 2*pi, rotation speed: 100ns -> pi pulse (assume dt = 1) 
hamiltonianJson = {
        "description": "Hamiltonian of a one-qubit system.\n",
        "h_str": ["-0.5*omega0*Z0", "omegaa*X0||D0"],
        "osc": {},
        "qub": {
            "0": 2
        },
        "vars": {
            "omega0": 6.2831853,
            "omegaa": 0.0314159
        } 
}

# Total time, T, of control pulse
T = 100
# Number of pulse samples
nbSamples = 200
W = 0.05
k = int(2 * nbSamples * W)
n_orders = 15
# Initialize Slepians
Slepians, eigenvalues = spectrum.dpss(nbSamples, (nbSamples*W), k)
Slepians = Slepians[:, 0:n_orders]

model = xacc.createPulseModel()
# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))
qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024 })
channelConfig = xacc.BackendChannelConfigs()
# Setting resolution of pulse
channelConfig.dt = nbSamples / T 
# Driving on resonance with qubit
channelConfig.loFregs_dChannels = [1.0]
model.setChannelConfigs(channelConfig)

weights = np.array([[-2.72625374, 0.17487545, 1.56643592, -3.73380211, 2.52358266, 3.38684644, -4.35706588, 3.63261847, -2.88180233, -3.87335379, -3.04824837,-0.13692838, 2.08190608, 0.30691899, -4.96701513]])

pulseData = (weights * Slepians).sum(axis=1)
# Add that square pulse instruction to XACC
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
qpu.execute(q, prog)

# Run the simulation
qpu.execute(q, prog)
resultProb = q.computeMeasurementProbability('1')
print(resultProb)
