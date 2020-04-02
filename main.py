import xacc
import gym
import spectrum
import json
#import gym_pulsecontrol
import numpy as np

from stable_baselines.common.policies import MlpPolicy
from stable_baselines import PPO2

# Total time, T, of control pulse
T = 100
# Number of pulse samples
nbSamples = 200
W = 0.02 #0.05
k = int(2 * nbSamples * W)
n_orders = 4 #15
# Initialize Slepians
Slepians, eigenvalues = spectrum.dpss(nbSamples, (nbSamples*W), k)
Slepians = Slepians[:, 0:n_orders]

gym.envs.register(
     id='PulseControl-v0',
     entry_point='gym_pulsecontrol.envs.pulsecontrol_env:PulseEnv',
)

env = gym.make('PulseControl-v0')
env.slepians_matrix = Slepians.copy()
env.n_orders = n_orders
env.nbQubits = 1
env.nbSamples = nbSamples
env.noise_frequency = 1.47969

hamiltonianJson = {
    "description": "One-qutrit Hamiltonian.",
    "h_latex": "",
    "h_str": ["(w - 0.5*alpha)*O0", "0.5*alpha*O0*O0", "O*(SM0 + SP0)||D0"],
    "osc": {},
    "qub": {
        "0": 3
    },
    "vars": {
        "w": 31.63772297724,
        "alpha": -1.47969,
        "O": 0.0314
    }
}

# Create a pulse system model object 
env.model = xacc.createPulseModel()
# Load the Hamiltonian JSON (string) to the system model
loadResult = env.model.loadHamiltonianJson(json.dumps(hamiltonianJson))
env.qpu = xacc.getAccelerator('QuaC', {'system-model': env.model.name(), 'shots': 1024 })
env.channelConfig = xacc.BackendChannelConfigs()
# Setting resolution of pulse
env.channelConfig.dt = nbSamples / T 
# Drive at resonance: 31.63772297724/(2pi)    
env.channelConfig.loFregs_dChannels = [5.0353]
env.model.setChannelConfigs(env.channelConfig)

drl_model = PPO2('MlpPolicy', env,
            learning_rate=0.0025,
            n_steps=128,
             verbose=0)
drl_model.learn(total_timesteps=10000)
drl_model.save("output_files/Single_Qubit_Model")