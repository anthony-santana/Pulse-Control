from ctrl_vqe_utils import *

# Query backend info (dt)
# Aer simulator
qpu = xacc.getAccelerator("aer:ibmq_armonk", {"sim-type": "pulse"})
# IBM backend
#qpu = xacc.getAccelerator("ibm:ibmq_bogota")

backend_properties = qpu.getProperties()
config = json.loads(backend_properties["config-json"])
# Number of qubits of the backend
nbQubits = config["n_qubits"]
ham_vars = config['hamiltonian']['vars']
for i in range(nbQubits):
    if ham_vars["delta" + str(i)] == 0:
        raise Exception('Zero anharmonicity detected.')
# Sampling time (dt)
dt = config["dt"]
#print("line 277:", config)

# Simple test: single qubit
ham = xacc.getObservable('pauli', '0.70710678118 X0 + 0.70710678118 Z0')
# Pulse parameter:
# Simple: single segment, i.e. square pulse.
pulse_opt = PulseOptParams(nb_qubits=1, nb_segments=1, total_length=20.0)
provider = xacc.getIRProvider("quantum")

# Optimization function:
# Note: We need to unpack the flattened x array into {amp, freq, time segment}.
def pulse_opt_func(x):
    # Handle both remote (ibm) and simulator (aer)
    isRemote = (qpu.name() == "ibm")
    
    # This indexing is throwing a bug so  I'm manually 
    # adding all 3 weight terms individually
    # amp = x[0:2] 
    if (pulse_opt.Hanning == True):
        amp = [x[0]]
        amp.append(x[1])
        amp.append(x[2])
        freq = x[3]
        pulse_opt.amplitude = amp
    elif (pulse_opt.Slepian == True):
        amp = [x[0]]
        amp.append(x[1])
        amp.append(x[2])
        # TODO: For loop to automatically append
        # based on nbOrders
        # Index this as x[(nbOrders+1)]
        bW = x[3]
        # Index this as x[(nbOrders+2)]
        phi = x[4]
        freq = 0.0
        pulse_opt.amplitude = amp
        pulse_opt.bandwidth = bW
        pulse_opt.phi = phi
    else:
        amp = x[0]
        freq = x[1]
        pulse_opt.amplitude = [[amp]]
    pulse_opt.freqs = [[freq]]
   
    # Construct the pulse program:
    program = provider.createComposite("vqe_pulse_composite")
    # program = provider.createComposite('pulse')
    # Create the pulse instructions
    # TODO: handle multiple channels....
    if (pulse_opt.Hanning == True):
        pulse_inst = provider.createInstruction(
            "pulse", [0], [], {
                "channel": "d0",
                "samples": pulse_opt.getHanningPulseSamples(0, dt)
            })
    elif (pulse_opt.Slepian == True):
        pulse_inst = provider.createInstruction(
            "pulse", [0], [], {
                "channel": "d0",
                "samples": pulse_opt.getSlepianPulseSamples(0, dt)
            })
    else:
        pulse_inst = provider.createInstruction(
            "pulse", [0], [], {
                "channel": "d0",
                "samples": pulse_opt.getPulseSamples(0, dt)
            })
    program.addInstruction(pulse_inst)
    energy = 0.0

    if not isRemote:
        buffer = xacc.qalloc(nbQubits)
        qpu.execute(buffer, program)
        state_vec = getStateVectorFromBuffer(buffer)
        # print(state_vec)
        obs = pauli_op_to_matrix(ham, buffer.size())
        energy = calculateExpVal(state_vec, obs)
        #print(state_vec, "->", energy)
    else:
        # Observe the pulse program: For use on IBM backend 
        fs_to_exe = ham.observe(program)
        # Execute
        for i in range(len(fs_to_exe)):
            fs = fs_to_exe[i]
            term = ham.getNonIdentitySubTerms()[i]
            coeff = 0.0
            for op in term:
                coeff = op[1].coeff().real
            buffer = xacc.qalloc(1)
            qpu.execute(buffer, fs)
            # print("Exp-Z =", buffer.getExpectationValueZ())
            energy = energy + coeff * buffer.getExpectationValueZ()
    print("Energy(", x, ") =", energy)
   
    # compute probabilities in 2 level approximation with RWA
    #   print(np.abs(expm(-1j * A * omegad0 * X / 2) @ np.array([1., 0]))**2)
    return energy

# pulse_opt.Hanning = True
pulse_opt.Slepian = True 

# Run the optimization loop:
# Optimizer:
# Single segment: 1 amplitude and 1 freq
# Make sure we don't create a pulse with amplitude > 1.0 (error)
optimizer = xacc.getOptimizer('nlopt', {
    "lower-bounds": [-10.0, -10.0, -10.0, 0.025, 0.0],
    "upper-bounds": [10.0, 10.0, 10.0, 0.45, 6.0],
    "initial-parameters": [-10.0,-10.0,-10.0,0.2,0.0],
})
result = optimizer.optimize(pulse_opt_func, 5)
print("Optimization Result:", result)

# optimizer = xacc.getOptimizer('nlopt', {
#     "lower-bounds": [-1.0, -1.0],
#     "upper-bounds": [1.0, 1.0],
#     "maxeval": 100
# })
# result = optimizer.optimize(pulse_opt_func, 2)
# print("Optimization Result:", result)

# Trying with scipy cobyla:
# Initializing each weight term to min. bound (-1.0) and phase
# at zero
# x0 = np.array([0.1,0.1,0.1,0.0])
# opt = optimize.minimize(pulse_opt_func, x0, method='cobyla')
# result = opt.x
# print("Optimization Result:", result)

# Not changing the freq., just varying the amplitude:
# Should see a Rabi oscillation....
# for ampl in np.linspace(0.0, 1.0, 100):
#     val = pulse_opt_func([ampl,0.0])
#     print("E(", ampl, ") =", val)