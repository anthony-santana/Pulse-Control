import sys, os, json, gym, numpy as np
import gym_pulsecontrol

from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

from types import MethodType
from stable_baselines.common.policies import MlpPolicy
from stable_baselines import PPO2

class RewardFunctions:

    def one_qubit(self):

        '''
        Single qubit reward function. Just returns the population of the first excited state. 
        Only works for full rotations around X-axis.
        '''

        self.pulseData = np.array(xacc.SlepianPulse(self._state, self.nbSamples, self.in_bW, self.in_K))
        pulseName = 'Slepian' + str(self.index)
        print(pulseName)
        xacc.addPulse(pulseName, self.pulseData)   
        provider = xacc.getIRProvider('quantum')
        prog = provider.createComposite('pulse')
        q = xacc.qalloc(self.nbQubits)
        slepianPulse = provider.createInstruction(pulseName, [0])
        slepianPulse.setChannel('d0')
        prog.addInstruction(slepianPulse)
        prog.addInstruction(xacc.gate.create("Measure", [0]))
        self.qpu.execute(q, prog) 
        fidelityResult = q.computeMeasurementProbability('1')
        print("\nFidelity: {}".format(fidelityResult))
        return  q.computeMeasurementProbability('1')


    def one_qutrit(self):

        '''
        Single qutrit reward function without Quantum Process Tomography. 
        Returns the population of the first excited state, and also prints out the population 
        of the second, leakage state. 
        '''

        self.pulseData = np.array(xacc.SlepianPulse(self._state, self.nbSamples, self.in_bW, self.in_K))
        pulseName = 'Slepian' + str(self.index)
        print(pulseName)
        xacc.addPulse(pulseName, self.pulseData)   
        provider = xacc.getIRProvider('quantum')
        prog = provider.createComposite('pulse')
        q = xacc.qalloc(self.nbQubits)
        slepianPulse = provider.createInstruction(pulseName, [0])
        slepianPulse.setChannel('d0')
        prog.addInstruction(slepianPulse)
        prog.addInstruction(xacc.gate.create("Measure", [0]))
        self.qpu.execute(q, prog)
        fidelityResult = q['DensityMatrixDiags'][1]
        leakage = q['DensityMatrixDiags'][2]
        print("\nFidelity: {}".format(fidelityResult))
        print("\nLeakgge: {}".format(leakage))
        return fidelityResult

    
    def one_qutrit_qpt(self):

        '''
        Single qutrit reward function that uses Quantum Process Tomography to calculate fidelity. 
        '''

        # Create the pulse as weighted sum of Slepian orders
        self.pulseData = np.array(xacc.SlepianPulse(self._state, self.nbSamples, self.in_bW, self.in_K))
        pulseName = 'Slepian' + str(self.index)
        print(pulseName)
        xacc.addPulse(pulseName, self.pulseData) 
        q = xacc.qalloc(self.nbQubits)  
        provider = xacc.getIRProvider('quantum')
        prog = provider.createComposite('pulse_composite')
        slepianPulse = provider.createInstruction(pulseName, [0])
        slepianPulse.setChannel('d0')
        prog.addInstruction(slepianPulse)
        qpt = xacc.getAlgorithm('qpt', {'circuit': prog, 'accelerator': self.qpu, 'optimize-circuit': False})
        qpt.execute(q)
        return qpt.calculate('fidelity', q, {'chi-theoretical-real': self.target_chi})

    
    def two_qubit(self):

        '''
        Two qubit reward function. Returns overlap between target density matrix 
        and the actual calcualted density matrix.
        '''

        # Create the pulse as weighted sum of Slepian orders
        self.pulseData = np.array(xacc.SlepianPulse(self._state, self.nbSamples, self.in_bW, self.in_K))
        pulseName = 'Slepian' + str(self.index)
        print(pulseName)
        xacc.addPulse(pulseName, self.pulseData)   
        provider = xacc.getIRProvider('quantum')
        prog = provider.createComposite('pulse_composite')
        slepianPulse = provider.createInstruction(pulseName, [0])
        slepianPulse.setChannel('d0')
        prog.addInstruction(slepianPulse)
        q = xacc.qalloc(self.nbQubits)
        q.addExtraInfo("target-dm-real", self.expectedDmReal)
        q.addExtraInfo("target-dm-imag", self.expectedDmImag)
        self.qpu.execute(q, prog)
        fidelityResult = q["fidelity"]
        print("\nFidelity: {}".format(fidelityResult))
        return fidelityResult