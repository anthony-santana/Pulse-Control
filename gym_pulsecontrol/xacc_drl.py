import sys, os, json, gym, numpy as np

sys.path.insert(1, '/home/cades/dev/Pulse_Control/')
import gym_pulsecontrol

sys.path.insert(1, '/home/cades/dev/Pulse_Control/gym_pulsecontrol')
from reward_functions import RewardFunctions 
from quac_backends import Backends

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
        self.env.initial_state = options['initial_state'] if 'initial_state' in options else 0
        self.env.qpt = options['qpt'] if 'qpt' in options else None
        self.env.qutrit = options['qutrit'] if 'qutrit' in options else False
        if self.env.qpt == True:
            self.env.reward_name = options['backend'] + '_qpt'
        else:
            self.env.reward_name = options['backend']
        self.env.initialize()

        # Parameters for PPO. The defaults are just the defaults from the stable_baselines backend.
        self.learning_rate = options['learning_rate'] if 'learning_rate' in options else 0.0025
        self.nsteps = options['nsteps'] if 'nsteps' in options else 128
        self.gamma = options['gamma'] if 'gamma' in options else 0.99
        self.ent_coef = options['ent_coef'] if 'ent_coef' in options else 0.01
        self.vf_coef = options['vf_coef'] if 'vf_coef' in options else 0.5
        self.max_grad_norm = options['max_grad_norm'] if 'max_grad_norm' in options else 0.5
        self.lam = options['lam'] if 'lam' in options else 0.95
        self.nminibatches = options['nminibatches'] if 'nminibatches' in options else 4
        self.noptepochs = options['noptepochs'] if 'noptepochs' in options else 4
        self.cliprange = options['cliprange'] if 'cliprange' in options else 0.2
        self.cliprange_vf = options['cliprange_vf'] if 'cliprange_vf' in options else None


    def execute(self):
        # Create a pulse system model object 
        self.env.model = xacc.createPulseModel()

        if self.env.reward_name == 'one_qubit':
            self.env.qpu = Backends.one_qubit(self)
            self.reward_function = RewardFunctions.one_qubit
        elif self.env.reward_name == 'one_qutrit':
            self.env.qpu = Backends.one_qutrit(self)
            self.reward_function = RewardFunctions.one_qutrit
        elif self.env.reward_name == 'one_qutrit_qpt':
            self.env.qpu = Backends.one_qutrit_qpt(self)
            self.reward_function = RewardFunctions.one_qutrit_qpt
        elif self.env.reward_name == 'two_qubit':
            self.env.qpu = Backends.two_qubit(self)
            self.reward_function = RewardFunctions.two_qubit

        self.env.reward_function = MethodType(self.reward_function, self.env)
    
        drl_model = PPO2('MlpPolicy', 
                    self.env,
                    gamma = self.gamma,
                    learning_rate = self.learning_rate,
                    ent_coef = self.ent_coef,
                    vf_coef = self.vf_coef,
                    n_steps = self.nsteps,
                    max_grad_norm = self.max_grad_norm,
                    lam = self.lam,
                    nminibatches = self.nminibatches,
                    noptepochs = self.noptepochs,
                    cliprange = self.cliprange,
                    cliprange_vf = self.cliprange_vf,
                    verbose=0,
                    n_cpu_tf_sess=1)
        drl_model.learn(total_timesteps=10000)
        drl_model.save("/home/cades/dev/Pulse_Control/output_files/")