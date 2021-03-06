import gym
import xacc
import numpy as np
import scipy as sp
from gym import error, spaces, utils
from gym.utils import seeding

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class PulseEnv(gym.Env):
    metadata = {'render.modes': ['human']}

    def __init__(self):
        self.end_episode_rewards = []
        self.current_reward = []
        self.optimal_pulse = []
        self.optimal_reward = []
        self.optimal_time = []
        self.counter = 0
        self.alpha_len = 0

    def initialize(self):
        self._state = np.zeros(self.in_K * self.nbPulses)

    def reward_function(self):
        pass

    def step(self, action):

        a = action.copy()
        self._state[self.alpha_len] = self._state[self.alpha_len] + a
        self.counter += 1
        reward = self.reward_function()
        print("alpha = ", repr(self._state))

        # a = action.copy()
        # self._state = self._state + a
        # self._state = np.clip(self._state, self.observation_space.low, self.observation_space.high)
        # self.counter += 1
        # reward = self.reward_function()
        # print("REWARD IS ", reward)
        # print("alpha = ", repr(self._state))

        if reward >= 0.9999:
            self.optimal_pulse = self.pulseData.copy()
            self.optimal_reward = reward
            self.optimal_time = self.T
            # print(self._state)
            plt.plot(self.optimal_pulse)
            plt.title(self.gate_operation + ' of Fidelity: ' + str(self.optimal_reward)[0:7] + ' With T = ' + str(self.optimal_time)[0:6])
            plt.ylabel(r'$\Omega(t)$')
            plt.xlabel(' Time Steps ')
            plt.savefig('/home/cades/dev/Pulse_Control/output_files/Optimal_Slepian' + str(self.counter) + '.png')
            np.savetxt('/home/cades/dev/Pulse_Control/output_files/optimal_pulse' + str(self.counter)+ '.csv', self.optimal_pulse, delimiter=',')
            exit()
        # done = bool((np.abs(1.0-reward) < 1e-4))
        done = bool(self.alpha_len == (self.in_K * self.nbPulses)-1)
        self.alpha_len += 1
        next_state = np.copy(self._state)
        return np.array(next_state), reward, done, {}

    def reset(self):
        # self._state = np.zeros(self.in_K * self.nbPulses)
        observation = np.copy(self._state)
        self.alpha_len = 0
        return observation

    def render(self, mode='human'):
        print('Reward=', self.reward_function())

    def close(self):
        pass

    @property
    def action_space(self):
        # return spaces.Box(low = self.action_range[0], high = self.action_range[1], shape=(self.in_K,))
        return spaces.Box(low = -1.0, high = 1.0, shape=(1, )) # shape=(self.in_K * self.nbPulses, )) 

    @property
    def observation_space(self):
        # return spaces.Box(low = self.observation_range[0], high = self.observation_range[1], shape=(self.in_K,))
        return spaces.Box(low = -5.0, high = 5.0, shape=(self.in_K * self.nbPulses, )) 