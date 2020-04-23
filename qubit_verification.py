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
T = 16.70261 ##100
# Number of pulse samples
nbSamples = 100
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
fidelity = np.zeros((10, 15))
deltas = np.zeros(10)
fidelity = np.zeros((10, 15))
deltas = np.zeros(10)



#circuit = provider.createComposite('U') 
#hadamard = provider.createInstruction('H', [0])
#circuit.addInstruction(hadamard)

for i in range(10):
    for j in range(15):
        #weights = np.array([-5., 1.65775359, 0.31740352, -3.74176682])
        noise_signal = (1 + delta)
        weights = np.array([0.92441954, 0.31973285, -4.77738279, -2.03994418])
        pulseData = ((weights * Slepians).sum(axis=1)) * noise_signal
        pulseName = 'Slepian' 
        xacc.addPulse(pulseName, pulseData)   
        q = xacc.qalloc(1)
        provider = xacc.getIRProvider('quantum')
        prog = provider.createComposite('pulse')
        slepianPulse = provider.createInstruction(pulseName, [0])
        slepianPulse.setChannel('d0')
        prog.addInstruction(slepianPulse)
        prog.addInstruction(xacc.gate.create("Measure", [0]))
        #qpt = xacc.getAlgorithm('qpt', {'circuit': prog, 'accelerator': qpu, 'optimize-circuit': False})
        #qpt.execute(q)
        # chi_real_vec = [0., 0., 0., 0., 
        #                 0., 2., 0., 0., 
        #                 0., 0., 0., 0.,
        #                 0., 0., 0., 0.]
        #F = qpt.calculate('fidelity', q, {'chi-theoretical-real': chi_real_vec})
        #F = qpt.calculate('fidelity', q, {'chi-theoretical-real':[0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 1., 0., 1.]})
        #print('\nFidelity: ', result)
        qpu.execute(q, prog)
        resultProb = q.computeMeasurementProbability('1')
        fidelity[i,j] = resultProb
        print(resultProb)
    delta += 0.01
    deltas[i] = delta

#np.savetxt('fidelity_delta.csv', fidelity, delimiter=',')
#np.savetxt('deltas.csv', deltas, delimiter=',')