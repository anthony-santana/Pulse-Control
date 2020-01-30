find_package(XACC REQUIRED)
include_directories(${XACC_INCLUDE_DIRS})

# add_executable(raw_pulse_ir raw_pulse_ir.cpp)
# target_link_libraries(raw_pulse_ir PRIVATE xacc_quac xacc::quantum_gate)

# add_executable(pulse_command_def pulse_command_def.cpp)
# target_link_libraries(pulse_command_def PRIVATE xacc_quac xacc::quantum_gate)


# add_executable(q_ctrl_pulses q_ctrl_pulses.cpp)
# target_link_libraries(q_ctrl_pulses PRIVATE xacc_quac xacc::quantum_gate)

# add_executable(hamiltonian_json_input hamiltonian_json_input.cpp)
# target_link_libraries(hamiltonian_json_input PRIVATE xacc_quac xacc::quantum_gate)

add_executable(one_qubit one_qubit.cpp)
target_link_libraries(one_qubit PRIVATE xacc_quac xacc::quantum_gate)

add_executable(cross_resonant_interaction cross_resonant_interaction.cpp)
target_link_libraries(cross_resonant_interaction PRIVATE xacc_quac xacc::quantum_gate)
