#%%
import gym
import gym_pulsecontrol

import spectrum
import numpy as np

from stable_baselines.common.policies import MlpPolicy
from stable_baselines import TRPO

class operators:
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    H = np.array([[1, 1], [1, -1]]) /  np.sqrt(2)
    I = np.array([[1, 0], [0, 1]])
    K = np.array([[1, -1j], [1, 1j]]) / np.sqrt(2)

tau = 200
N = 1000
resolution = tau/N
W = 0.02
k = int(2*N*W)
n_orders = 36 
# Initialize Slepians
Slepians, eigenvalues = spectrum.dpss(N, (N*W), k)
Slepians = Slepians[:, 0:n_orders]

###############################################################################


env = gym.make('PulseControl-v0')
###############################################################################

env.qubit_frequency = 5.114
env.slepians_matrix = Slepians.copy()
env.N = N
env.resolution = resolution
env.target_gate = operators.H
env.d = 2 # 2-dimensional Hilbert Space
env.max_steps = (n_orders - 1)


model = TRPO('MlpPolicy', env, timesteps_per_batch=n_orders, 
             tensorboard_log='/Users/anthonysantana/gym/gym-PulseControl', 
             full_tensorboard_log=True,verbose=1)
model.learn(total_timesteps=50*n_orders)

#%%
obs = env.reset()
#while True:
for _ in range(n_orders):
    action, _state = model.predict(obs)
    obs, rewards, dones, info = env.step(action)
    env.render()
    
    
    
    
    
    
    
#%%   
# Only run if there are issues with the environment not being registered
#import gym
#env_dict = gym.envs.registration.registry.env_specs.copy()
#for env in env_dict:
 #    if 'PulseControl-v0' in env:
  #        print('Remove {} from registry'.format(env))
   #       del gym.envs.registration.registry.env_specs[env]
#import gym_pulsecontrol'
