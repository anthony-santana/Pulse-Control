import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

import numpy as np 

index = 0
# Density Matrix for {X[q0], Ry[0.59, q1], CNOT}
expectedDmReal = np.array([
    0, 0, 0, 0,
    0, 0.08572349, 0.27995531, 0, 
    0, 0.27995531, 0.91427651, 0,
    0, 0, 0, 0
], dtype = np.float64)
expectedDmImag = np.zeros(16)


def grafs(alpha):

    global index
    global expectedDmReal
    global expectedDmImag

    provider = xacc.getIRProvider('quantum')
    prog = provider.createComposite('pulse_composite')

    model = xacc.createPulseModel()
    model.setQubitInitialPopulation(0, 0.)
    channelConfig = xacc.BackendChannelConfigs()
    channelConfig.dt = 112 / 800
    model.setChannelConfigs(channelConfig)
    qpu = xacc.getAccelerator('QuaC:Default2Q')

    pulseData = np.array(xacc.SlepianPulse(alpha, 112, 0.02, 4))
    pulseName = 'Slepian' + str(index)
    index += 1
    xacc.addPulse(pulseName, pulseData)   
    slepianPulse = provider.createInstruction(pulseName, [0])
    slepianPulse.setChannel('d0')
    prog.addInstruction(slepianPulse)

    q = xacc.qalloc(2)
    q.addExtraInfo("target-dm-real", expectedDmReal)
    q.addExtraInfo("target-dm-imag", expectedDmImag)
    qpu.execute(q, prog)
    fidelityResult = q["fidelity"]
    print(fidelityResult)
    print(alpha)

    return (1. - fidelityResult) #, gradientVector

opt = xacc.getOptimizer("mlpack", {'mlpack-optimizer':'l-bfgs', 'initial-parameters':[.01, .01, .01, .01]})
result = opt.optimize(grafs, 4)