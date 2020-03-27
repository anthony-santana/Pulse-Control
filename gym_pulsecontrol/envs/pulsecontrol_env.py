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
        self._state = np.zeros(15)
        #self._state = np.zeros(self.slepians_matrix.shape[1]) 
        self.end_episode_rewards = []
        self.current_reward = []
        self.optimal_pulse = []
        self.optimal_reward = []
        self.index = 0

    def reward_function(self):
        self.pulseData = (self._state * self.slepians_matrix).sum(axis=1)
        # Add that square pulse instruction to XACC
        pulseName = 'Slepian' + str(self.index)
        print(pulseName)
        
        xacc.addPulse(pulseName, self.pulseData)   
        # TODO: handle arbitrary number of qubit
        # self.nbQubits contains this information
        q = xacc.qalloc(1)
        # Create the quantum program that contains the slepian pulse
        # and the drive channel (D0) is set on the instruction
        provider = xacc.getIRProvider('quantum')
        prog = provider.createComposite('pulse')
        slepianPulse = provider.createInstruction(pulseName, [0])
        # TODO: need to handle multiple channels
        # we have this information (based on the Hops array)
        slepianPulse.setChannel('d0')
        prog.addInstruction(slepianPulse)
        # TODO: handle multiple qubits
        # Measure Q0 (using the number of shots that was specified above)
        prog.addInstruction(xacc.gate.create("Measure", [0]))
        self.qpu.execute(q, prog)
        
        # TODO: define a generic reward function based on the
        # target unitary matrix (self.targetU) that is sent to here.
        return q.computeMeasurementProbability('1')

    def step(self, action):
        a = action.copy()
        self._state = self._state + a
        self._state = np.clip(self._state, self.observation_space.low, self.observation_space.high)
        self.index += 1
        reward = self.reward_function()
        print("REWARD IS ", reward)
        if reward >= 0.9999:
            self.optimal_pulse = self.pulseData.copy()
            self.optimal_reward = reward
            print(self._state)
            plt.plot(self.optimal_pulse)
            plt.title(' X-Gate of Fidelity: ' + str(self.optimal_reward))
            plt.ylabel(r'$\Omega(t)$')
            plt.xlabel(' Time ')
            plt.savefig('Optimal_Slepian' + str(self.index) + '.png')
            np.savetxt('optimal_pulse' + str(self.index)+ '.csv', self.optimal_pulse, delimiter=',')
            exit()
        done = bool((np.abs(1.0-reward) < 1e-4))
        next_state = np.copy(self._state)
        return np.array(next_state), reward, done, {}

    def reset(self):
        self._state = np.zeros(self.slepians_matrix.shape[1])
        observation = np.copy(self._state)
        return observation

    def render(self, mode='human'):
        print('Reward=', self.reward_function())

    def close(self):
        pass

    @property
    def action_space(self):
        return spaces.Box(low=-0.25, high=0.25, shape=(15,))#shape=(self.n_orders,))

    @property
    def observation_space(self):
        return spaces.Box(low=-5.0, high=5.0, shape=(15,))#shape=(self.n_orders,))