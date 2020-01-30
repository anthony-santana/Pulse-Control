include_directories(${CMAKE_BINARY_DIR})
add_xacc_test(Simple)
target_link_libraries(SimpleTester xacc xacc_quac xacc::quantum_gate)