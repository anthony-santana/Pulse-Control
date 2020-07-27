import sys, json
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

class Backends:

    def one_qubit(self):
        hamiltonianJson = {
            "description": "Hamiltonian of a one-qubit system.\n",
            "h_str": ["-0.5*omega0*Z0", "omegaa*X0||D0"],
            "osc": {},
            "qub": {
                "0": 2
            },
            "vars": {
                "omega0": 6.2831853,
                "omegaa": 0.0314159
            } 
        }
        # Load the Hamiltonian JSON (string) to the system model
        self.env.loadResult = self.env.model.loadHamiltonianJson(json.dumps(hamiltonianJson))
        self.env.channelConfig = xacc.BackendChannelConfigs()
        self.env.channelConfig.dt = self.env.nbSamples / self.env.T 
        self.env.model.setChannelConfigs(self.env.channelConfig)
        self.env.channelConfig.loFregs_dChannels = [1.0]
        return xacc.getAccelerator('QuaC', {'system-model': self.env.model.name(), 'shots': 1024})
    
    def one_qutrit(self):
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
        # Load the Hamiltonian JSON (string) to the system model
        self.env.loadResult = self.env.model.loadHamiltonianJson(json.dumps(hamiltonianJson))
        self.env.channelConfig = xacc.BackendChannelConfigs()
        self.env.channelConfig.dt = self.env.nbSamples / self.env.T 
        self.env.model.setChannelConfigs(self.env.channelConfig)
        self.env.channelConfig.loFregs_dChannels = [5.0353]
        return xacc.getAccelerator('QuaC', {'system-model': self.env.model.name(), 'shots': 1024})

    def two_qubit(self):
        # Set control and target qubit to 0 -> initial state 00
        self.env.model.setQubitInitialPopulation(self.env.initial_state[0], self.env.initial_state[1])
        self.env.channelConfig = xacc.BackendChannelConfigs()
        self.env.channelConfig.dt = self.env.nbSamples / self.env.T 
        self.env.model.setChannelConfigs(self.env.channelConfig)
        return xacc.getAccelerator('QuaC:Default2Q')

    def two_qutrit(self):
        pass