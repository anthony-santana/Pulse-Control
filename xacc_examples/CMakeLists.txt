find_package(XACC REQUIRED)
include_directories(${XACC_INCLUDE_DIRS})

add_executable(raw_pulse_ir raw_pulse_ir.cpp)
target_link_libraries(raw_pulse_ir PRIVATE xacc_quac xacc::quantum_gate)

add_executable(pulse_command_def pulse_command_def.cpp)
target_link_libraries(pulse_command_def PRIVATE xacc_quac xacc::quantum_gate)