# XACC Pulse Optimizer plugin implementation that can be used by IR transformation service.
# IMPORTANT: this file must be copied/installed to XACC_INSTALL_DIR/py-plugins
# e.g. $HOME/.xacc/py-plugins

import xacc
from pelix.ipopo.decorators import ComponentFactory, Property, Requires, Provides, \
    Validate, Invalidate, Instantiate

# Import packages required for pulse optimization
import gym
import spectrum
import json
import numpy as np
from stable_baselines.common.policies import MlpPolicy
from stable_baselines import PPO2

# Pulse optimization plugin
@ComponentFactory("py_pulse_optimizer_factory")
@Provides("optimizer")
@Property("_optimizer", "optimizer", "ml_optimizer")
@Property("_name", "name", "ml_optimizer")
@Instantiate("mlopt_instance")
class MlPulseOptimizer(xacc.Optimizer):
    def __init__(self):
        xacc.Optimizer.__init__(self)
        self.options = {}
        self.pulseOpts = None

    def name(self):
        return 'ml_optimizer'

    def setOptions(self, opts):
        self.pulseOpts = opts
        # These params are always present in the option map/dict
        # when using the IR transformation service.
        if 'dimension' in opts:
            self.dimension = opts['dimension']
        # Target unitary matrix
        if 'target-U' in opts:
            self.targetU = opts['target-U']
        # Static Hamiltonian
        if 'static-H' in opts:
            self.H0 = opts['static-H']
        # Control Hamiltonian (list)
        if 'control-H' in opts:
            self.Hops = opts['control-H']
        # Max time horizon
        if 'max-time' in opts:
            self.tMax = opts['max-time']
        # Pulse sample dt (number of samples over the time horizon)
        if 'dt' in opts:
            self.dt = opts['dt']
        if 'hamiltonian-json' in opts:
            self.hamJson = opts['hamiltonian-json']
        # Note: if the method requires specific parameters,
        # we can require those params being specified in the IR transformation options
        # then propagate to here. 

    # This is main entry point that the high-level
    # IR transformation service will call.
    # This needs to return the pair (opt-val, pulses)
    # where opt-val is a floating point number (final value of the cost function);
    # pulses is a single array of all control pulses (one for each control-H term)
    # (appending one pulse array after another).
    def optimize(self):
        # TODO: we can now call any Python lib to
        # perform pulse optimization (marshalling the options/parameters if required)
        # For example, one can use Qutip pulse optimization:
        # Notes about data types: 
        # - targerU: flatten (row-by-row) U matrix into a 1-D array of complex numbers
        # - H0: string-type representation of the static Hamiltonian:
        # e.g.: 0.123 Z0Z1
        # - Hops: array of strings represent terms on the Hamiltonian which can be controlled.
        # Depending on the specific library we use for pulse optimization,
        # we may need to marshal these data types accordingly.
        #print('Target U: ')
        #print(self.targetU)
        #print('Hops: ')
        #print(self.Hops)
        # At the very basic level, we have full access to the Hamiltonian JSON here
        # i.e. we could use any OpenPulse-compatible solver to perform simulation/optimization 
        # if needed for the IR transformation (optimizing for the target unitary matrix).
        #print(self.hamJson)
        
        # Call our Slepian pulse optimization code
        # Total time, T, of control pulse
        T = self.tMax
        # Number of pulse samples
        nbSamples = 200
        W = 0.05
        k = int(2 * nbSamples * W)
        n_orders = 15
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

        # Create a pulse system model object 
        env.model = xacc.createPulseModel()
        # Load the Hamiltonian JSON (string) to the system model
        # Note: This Hamiltonian Json is from the IR transformation.
        loadResult = env.model.loadHamiltonianJson(self.hamJson)
        env.qpu = xacc.getAccelerator('QuaC', {'system-model': env.model.name(), 'shots': 1024 })
        env.channelConfig = xacc.BackendChannelConfigs()
        # Setting resolution of pulse
        env.channelConfig.dt = self.dt
        
        # The target unitary:
        # The optimizer must make use of this Target U in the reward function calculation.
        env.targetU = self.targetU
        env.nbQubits = self.dimension
        # Driving on resonance with qubit
        # TODO: this needs to be passed from the IR transformation as well.
        # Most other optimization methods have their own solver/propagator,
        # this method requires a pulse backend to do optimization.
        # One potential solution is passing the QPU all the way
        # from IR transformation down to here (TODO: tnguyen)
        env.channelConfig.loFregs_dChannels = [1.0]
        env.model.setChannelConfigs(env.channelConfig)

        drl_model = PPO2('MlpPolicy', env,
                    learning_rate=0.0025,
                    verbose=0)
        drl_model.learn(total_timesteps = 400)
        # The optimal pulse sequence
        pulse = env.optimal_pulse
        # Return the final cost functional value and the array of pulse samples.
        return (0.0, pulse)