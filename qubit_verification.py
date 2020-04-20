import sys
from numpy import genfromtxt
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import os
import xacc
import json
import spectrum
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
W = 0.02
k = int(2 * nbSamples * W)
n_orders = 4 
# Initialize Slepians
Slepians, eigenvalues = spectrum.dpss(nbSamples, (nbSamples*W), k)
Slepians = Slepians[:, 0:n_orders]

model = xacc.createPulseModel()
# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))
qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024, 'logging-period': 0.1 })
channelConfig = xacc.BackendChannelConfigs()
# Setting resolution of pulse
channelConfig.dt = nbSamples / T 
# Driving on resonance with qubit
channelConfig.loFregs_dChannels = [1.0]
model.setChannelConfigs(channelConfig)

delta = -0.1
# fidelity = np.zeros((20, 50))
# deltas = np.zeros(20)
# for i in range(20):
fidelity = np.zeros((20, 50))
deltas = np.zeros(20)
for i in range(20):
    for j in range(50):
        noise_signal = (1 + delta)
        weights = np.array([-5., 1.65775359, 0.31740352, -3.74176682])
        pulseData = (weights * Slepians).sum(axis=1) * noise_signal
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
        resultProb = q.computeMeasurementProbability('1')
        fidelity[i,j] = resultProb
        print(resultProb)
    delta += 0.01
    deltas[i] = delta
np.savetxt('fidelity_delta.csv', fidelity, delimiter=',')
np.savetxt('deltas.csv', deltas, delimiter=',')