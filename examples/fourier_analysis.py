import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

import numpy as np
from scipy.fft import fft

# RL Optimized CNOT from |11> with fidelity of 0.9888
alpha =  np.array([ 0.70309901, -1.        ,  0.        ,  0.        ,  0.        ])
pulse = np.array(xacc.SlepianPulse(alpha, 512, 0.02, 5))
fourier = fft(pulse)

plt.plot(fourier)
plt.savefig('/home/cades/dev/Pulse_Control/output_files/cnot_fourier' + '.png')
