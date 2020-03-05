# Sweeping LO frequency to find qubit resonance freq.
# We need to have the XACC install directory in the Python path.
# Just in case users haven't already done that, set it here.
import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')

# Import xacc and quaC python wrapper
import os
import xacc
import json
# Note to ORNL CADES users:
# Make sure numpy and matplotlib are installed
# If not, these can be installed by:
# sudo apt install python3-numpy
# sudo apt install python3-matplotlib
import numpy as np
import matplotlib
# CADES VM don't have display
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# The Hamiltonian JSON object (OpenPulse format)
hamiltonianJson = {
        "description": "Qubits are modelled as a two level system. System of 2 qubits.\n",
        "h_str": ["wq0/2*Z0", "wq1/2*Z1", "omegad0*X0||D0", "omegad1*X1||D1", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
        "osc": {},
        "qub": {
            "0": 2,
            "1": 2
        },
        "vars": {
            "omegad0": 1.303125,
            "omegad1": 0.97, 
            "wq0": 30.91270129264568,
            "wq1": 30.36010168900955,
            "jq0q1": 0.04
        } 
}

# Create a pulse system model object 
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

if loadResult is True :
    # Simple circuit: X gate (pulse)
    xacc.qasm('''.compiler xasm
    .circuit test
    .qbit q
    X(q[0]);
    Measure(q[0]);
    ''')
    prog = xacc.getCompiled('test')
    qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024 })
    # Load the pulse and cmd-def definitions from the demo JSON file 
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    backendJson = open('../resources/test_backends.json', 'r').read()
    qpu.contributeInstructions(backendJson)  
    
    # Sweep Q1:
    freq1 = np.linspace(4.886, 4.95, 50)
    result1 = np.zeros(freq1.size)
    i = 0
    for freq in freq1:
        channelConfig = xacc.BackendChannelConfigs()
        # dt (time between data samples)
        channelConfig.dt = 3.5555555555555554
        channelConfig.loFregs_dChannels = [freq, 4.83]
        model.setChannelConfigs(channelConfig)
        # Run the simulation
        q = xacc.qalloc(2)
        qpu.execute(q, prog)
        # Store the result (probability of 1)
        result1[i] = q.computeMeasurementProbability('1')
        i = i + 1
    # Plot the result (Q1 only)
    # We can do the same for Q2 (different frequency range)
    plt.scatter(freq1, result1)
    plt.savefig('plotQ1.png')
else :
    print("Failed to load Hamiltonian Json. Please check your input.")
