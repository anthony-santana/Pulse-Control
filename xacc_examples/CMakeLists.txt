find_package(XACC REQUIRED)
include_directories(${XACC_INCLUDE_DIRS})
link_libraries(-Wl,-rpath,${XACC_ROOT}/lib -L${XACC_ROOT}/lib)

add_executable(one_qubit one_qubit.cpp)
target_link_libraries(one_qubit PRIVATE xacc_quac xacc::quantum_gate)

add_executable(cross_resonant_interaction cross_resonant_interaction.cpp)
target_link_libraries(cross_resonant_interaction PRIVATE xacc_quac xacc::quantum_gate)

add_executable(pulse_optimization pulse_optimization.cpp)
target_link_libraries(pulse_optimization PRIVATE xacc_quac xacc::quantum_gate)
