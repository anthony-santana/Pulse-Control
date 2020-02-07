include_directories(${CMAKE_BINARY_DIR})
add_xacc_test(Simple)
target_link_libraries(SimpleTester xacc xacc_quac xacc::quantum_gate)

add_xacc_test(OpenPulse)
target_link_libraries(OpenPulseTester xacc xacc_quac xacc::quantum_gate)

add_xacc_test(Circuit)
target_link_libraries(CircuitTester xacc xacc_quac xacc::quantum_gate)

add_xacc_test(DigitalGate)
target_link_libraries(DigitalGateTester xacc xacc_quac xacc::quantum_gate)