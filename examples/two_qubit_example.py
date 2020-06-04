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
env.reward_bounds = [10.0, 100.0]
env.expectedDmReal = np.array([
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 0
], dtype = np.float64)
env.expectedDmImag = np.zeros(16)

# Get the two-qubit backend (QuaC default):
# This will load the default two-qubit backend:
'''
{
    "description": "Two-qubit Hamiltonian",
    "h_str": ["_SUM[i,0,1,wq{i}*O{i}]", "_SUM[i,0,1,delta{i}*O{i}*(O{i}-I{i})]", "_SUM[i,0,1,omegad{i}*X{i}||D{i}]", "omegad1*X0||U0", "omegad0*X1||U1", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
    "osc": {},
    "qub": {
        "0": 2,
        "1": 2
    },
    "vars": {
        "wq0": 30.518812656662774, 
        "wq1": 31.238229295532093,
        "delta0": -2.011875935,
        "delta1": -2.008734343,
        "omegad0": -1.703999855,
        "omegad1": -1.703999855,
        "jq0q1": 0.011749557 
    }
}
'''

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