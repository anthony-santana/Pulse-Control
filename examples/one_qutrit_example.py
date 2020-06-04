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
env.target_chi = [0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.] # X-Gate

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
env.channelConfig.loFregs_dChannels = [5.0353]
env.model.setChannelConfigs(env.channelConfig)

def reward_function(self):
    # Create the pulse as weighted sum of Slepian orders
    self.pulseData = (self._state * self.slepians_matrix).sum(axis=1)
    pulseName = 'Slepian' + str(self.index)
    print(pulseName)
    xacc.addPulse(pulseName, self.pulseData)   
    q = xacc.qalloc(self.nbQubits)
    # Create the quantum program that contains the slepian pulse
    # and the drive channel (D0) is set on the instruction
    provider = xacc.getIRProvider('quantum')
    prog = provider.createComposite('pulse')
    slepianPulse = provider.createInstruction(pulseName, [0])
    # TODO: need to handle multiple channels
    # we have this information (based on the Hops array)
    slepianPulse.setChannel('d0')
    prog.addInstruction(slepianPulse)
    qpt = xacc.getAlgorithm('qpt', {'circuit': prog, 'accelerator': self.qpu, 'optimize-circuit': False})
    qpt.execute(q)
    return qpt.calculate('fidelity', q, {'chi-theoretical-real': self.target_chi})
# Passing reward function to the backend
env.reward_function = MethodType(reward_function, env)

drl_model = PPO2('MlpPolicy', env,
            learning_rate=0.0025,
            n_steps=128,
             verbose=0)
drl_model.learn(total_timesteps=10000)
drl_model.save("output_files/Single_Qubit_Model")