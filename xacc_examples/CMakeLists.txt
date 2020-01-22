find_package(XACC REQUIRED)
include_directories(${XACC_INCLUDE_DIRS})

add_executable(test test.cpp)
target_link_libraries(test PRIVATE xacc_quac xacc::quantum_gate)