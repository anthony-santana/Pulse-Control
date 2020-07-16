import xacc
import gym
import spectrum
import json
#import gym_pulsecontrol
import numpy as np

from types import MethodType
from stable_baselines.common.policies import MlpPolicy
from stable_baselines import PPO2


# Total time, T, of control pulse
T = 100
# Number of pulse samples
nbSamples = 100
W = 0.02
k = int(2 * nbSamples * W)
n_orders = 4
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
env.nbQubits = 2
env.nbSamples = nbSamples
env.T = T
env.T_range = [10.0, 100.0]
#env.reward_bounds = [10.0, 100.0]
env.expectedDmReal = np.array([
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 0
], dtype = np.float64)
env.expectedDmImag = np.zeros(16)

# Create a pulse system model object 
env.model = xacc.createPulseModel()
env.qpu = xacc.getAccelerator('QuaC:Default2Q')
env.channelConfig = xacc.BackendChannelConfigs()
env.channelConfig.dt = nbSamples / env.T 
env.model.setChannelConfigs(env.channelConfig)

def reward_function(self):
    # Running last index of state vector through affine transform to get T
    self.T = self.affine_transform()
    # Changing dt on the backend
    self.channelConfig.dt = nbSamples / self.T 
    self.model.setChannelConfigs(self.channelConfig)
    _state = self._state[:-1]
    # Create the pulse as weighted sum of Slepian orders
    self.pulseData = (_state * self.slepians_matrix).sum(axis=1)
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
drl_model.save("output_files/Double_Qubit_Model")