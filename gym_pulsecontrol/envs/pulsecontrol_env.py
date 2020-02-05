import gym
from gym import error, spaces, utils
from gym.utils import seeding

import numpy as np
import scipy as sp

"""
Required Attributes in Main:
1. Qubit frequency
2. Slepian matrix 
3. N for length of Slepian sequence
4. Resolution for resolution of control field
5. Target gate operation
6. Hilbert space dimension
7. "max_steps" -- make sure it's the (max - 1) to prevent indexing errors
"""
class PulseEnv(gym.Env):
    metadata = {'render.modes': ['human']} # could also do 'console'
    def _init_(self):
        super(PulseEnv, self).__init__()
        """
        Initialized the state of the environment
        """
        self.action_space = spaces.Box(low=-10.0, high=10.0, shape=(1,))
        self.observation_space = spaces.Box(low=-5.0, high=5.0, shape=(1,))
        #self.spec = EnvSpec(observation_space=self.observation_space, action_space=self.action_space)
        self.reset()
    
    def Hamiltonian(self, ctrl_amplitude, time):
        """
        Takes the qubit frequency as an attribute and returns the Hamiltonian as a function of time
        """
        #time = raise NotImplementedError
        omega_d = self.qubit_frequency
        B_x = np.real(ctrl_amplitude * np.exp(1j * omega_d * time))
        H_drift = (omega_d / 2) * (operators.I - operators.Z)
        return H_drift + (B_x * operators.X)
    
    def unitary_evolution(self):
        """
        Takes all of the chosen Slepian weights from the episode, calculates the control field, then calculates
        the resulting unitary of the system
        """
        states_vector_np = np.zeros((self.max_steps+1,1))
        states_vector_np[0:len(self.states_vector)] = np.asarray(self.states_vector)
        control_field = (states_vector_np.T * self.slepians_matrix).sum(axis=1)
        U = [1] * self.N
        for i in range(self.N):
            time = i
            ctrl_amplitude = control_field[i]
            H = self.Hamiltonian(ctrl_amplitude, time)
            U[(self.N - 1) - i] = sp.linalg.expm(-1j * self.resolution * H)
        #U = np.asarray(U)
        return self.trotter_multiply(U)
    
    def trotter_multiply(self, U): 
        # Multiply elements one by one 
        arr = np.array([[1,0],[0,1]])
        for x in U: 
             arr = np.matmul(arr, x)
        return arr
        
    def reward_function(self):
        """
        Takes the evolved quantum state and returns its overlap with a target gate operation
        """
        U_targ = self.target_gate
        U_evolved = self.unitary_evolution()
        return (1 / self.d) * np.abs(np.trace(np.dot(np.matrix.getH(U_targ), U_evolved)))
    
    def step(self, action):
        """
        Run one timestep of the environment's dynamics. When end of episode
        is reached, reset() should be called to reset the environment's internal state.
        Input
        -----
        action : an action provided by the environment
        Outputs
        -------
        (observation, reward, done, info)
        observation : agent's observation of the current environment
        reward [Float] : amount of reward due to the previous action
        done : a boolean, indicating whether the episode has ended
        info : a dictionary containing other diagnostic information from the previous action
        """
        a = action.copy()
        a = np.clip(a, self.action_space.low, self.action_space.high)
        self._state = self._state + a
        # May remove this line:
        self._state = np.clip(self._state, self.observation_space.low, self.observation_space.high)
        self.states_vector.append(self._state)
        reward = self.reward_function()
        done = (len(self.states_vector) == self.max_steps)
        next_observation = np.copy(self._state)
        # Optional info call so that we can pass additional info
        info = {}
        return np.array(next_observation), reward, done, info
    
    def reset(self):
        """
        Resets the state of the environment, returning an initial observation.
        Outputs
        -------
        observation : the initial observation of the space. (Initial reward is assumed to be 0.)
        """
        # The observation must be a numpy array
        self._state = (np.asarray(0.0).astype(np.float64)).reshape(1,)
        #self._state = np.random.uniform(self.observation_space.low, self.observation_space.high, size=(1,)).astype(np.float64)
        self.states_vector = []
        observation = np.copy(self._state)
        return observation
    
    def render(self, mode='human'):
        #reward = self.reward_function()
        print('Reward=', self.reward_function())
    
    def close(self):
        pass
    
    @property
    def action_space(self):
        return spaces.Box(low=-10.0, high=10.0, shape=(1,))
    
    @property
    def observation_space(self):
        return spaces.Box(low=-5.0, high=5.0, shape=(1,))
    
    #@property
    #def reward_range(self):
     #   return (0.0, 1.0)
    
    #@property
    #def spec(self):
     #   return EnvSpec(observation_space=self.observation_space, action_space=self.action_space)

class operators:
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    H = np.array([[1, 1], [1, -1]]) /  np.sqrt(2)
    I = np.array([[1, 0], [0, 1]])
    K = np.array([[1, -1j], [1, 1j]]) / np.sqrt(2)