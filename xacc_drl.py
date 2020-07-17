import sys, os, json, gym, numpy as np
import gym_pulsecontrol

from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

from types import MethodType
from stable_baselines.common.policies import MlpPolicy
from stable_baselines import PPO2

class OptimalControl:

    def __init__(self, options):
        gym.envs.register(
            id='PulseControl-v0',
            entry_point='gym_pulsecontrol.envs.pulsecontrol_env:PulseEnv',
        )
        self.env = gym.make('PulseControl-v0')
        self.env.nbQubits = options['nbQubits'] if 'nbQubits' in options else print('nbQubits needed')
        self.env.nbSamples = options['slepian_parameters'][0]
        self.env.in_bW = options['slepian_parameters'][1]
        self.env.in_K = options['slepian_parameters'][2]
        self.env.T = options['slepian_parameters'][3]
        self.env.gate_operation = options['gate_operation'] if 'gate_operation' in options else 'Gate'
        self.env.expectedDmReal = options['expectedDmReal'] if 'expectedDmReal' in options else None
        self.env.expectedDmImag = options['expectedDmImag'] if 'expectedDmImag' in options else None
        self.env.initialize()
    
    # Use a hardcoded reward function for now, but the eventual step will be to add a repository of
    # reward functions that users can select from
    # Eventually, make it so that they can create their own custom reward functions
    def reward_function(self):
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

    def execute(self):
        # Create a pulse system model object 
        self.env.model = xacc.createPulseModel()
        self.env.qpu = xacc.getAccelerator('QuaC:Default2Q')
        self.env.channelConfig = xacc.BackendChannelConfigs()
        self.env.channelConfig.dt = self.env.nbSamples / self.env.T 
        self.env.model.setChannelConfigs(self.env.channelConfig)
        # Set control and target qubit to 0 -> initial state 00
        self.env.model.setQubitInitialPopulation(0, 0)

        self.env.reward_function = MethodType(self.reward_function, self.env)

        drl_model = PPO2('MlpPolicy', 
                    self.env,
                    learning_rate=0.0025,
                    n_steps=128,
                    verbose=0,
                    n_cpu_tf_sess=1)
        drl_model.learn(total_timesteps=10000)