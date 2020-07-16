import sys, os, json, gym, numpy as np
import gym_pulsecontrol

import spectrum

# Alternative to the following two lines is to run
# from the IDE terminal: export PYTHONPATH=$PYTHONPATH:$HOME/.xacc
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

from types import MethodType
from stable_baselines.common.policies import MlpPolicy
from stable_baselines import PPO2


gym.envs.register(
     id='PulseControl-v0',
     entry_point='gym_pulsecontrol.envs.pulsecontrol_env:PulseEnv',
)

env = gym.make('PulseControl-v0')
env.nbQubits = 2
env.nbSamples = 512
env.T = 2 * np.pi
env.in_bW = 0.025
env.in_K = 4 # int(2 * env.nbSamples * env.in_bW)

# Density Matrix for {X[q0], Ry[0.59, q1], CNOT}
env.expectedDmReal = np.array([
    0.08452966, 0, 0.08452966, 0,
    0, 0, 0, 0,
    0.08452966, 0, 0.08452966, 0.,
    0, 0, 0., 0.
], dtype = np.float64)
env.expectedDmImag = np.zeros(16)
# Used for plot titles only:
env.gate_name = 'X[q0], Ry[0.59, q1], CNOT'

# Create a pulse system model object 
env.model = xacc.createPulseModel()
env.qpu = xacc.getAccelerator('QuaC:Default2Q')
env.channelConfig = xacc.BackendChannelConfigs()
env.channelConfig.dt = env.nbSamples / env.T 
env.model.setChannelConfigs(env.channelConfig)
# Set control and target qubit to 0 -> initial state 00
env.model.setQubitInitialPopulation(0, 0)

def reward_function(self):
    # Create the pulse as weighted sum of Slepian orders
    self.pulseData = np.array(xacc.SlepianPulse(self._state, self.nbSamples, in_bW, in_K))
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
env.reward_function = MethodType(reward_function, env)

drl_model = PPO2('MlpPolicy', env,
            learning_rate=0.0025,
            n_steps=128,
             verbose=0,
             n_cpu_tf_sess=1)
drl_model.learn(total_timesteps=10000)
