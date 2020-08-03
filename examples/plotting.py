import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc
import numpy as np

[ ]

alpha = np.array([-0.16151045, -0.10867336, -0.65445626,  0.24819908,  0.66114694, 1., -0.5727458, -0.6626659,  0.52905434, -0.27316469, 0.47111198, 0.38856477, 0.08256996,  1., -1., -0.47658846, 0.29641816, 1., -0.55198687, -0.52067912])
pulse0 = np.array(xacc.SlepianPulse(alpha, 512, 0.02, 5))

plt.plot(pulse0)
plt.title('{X,Y,CNOT} of T=800, Fidelity = 0.94762')
plt.savefig('/home/cades/dev/Pulse_Control/output_files/test_slepian.png')