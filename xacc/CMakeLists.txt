
find_package(XACC REQUIRED)

set (LIBRARY_NAME xacc_quac)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file (GLOB HEADERS *.hpp)
file (GLOB SRC *.cpp Accelerator/*.cpp Utils/*.cpp Executor/*.cpp)

usFunctionGetResourceSource(TARGET ${LIBRARY_NAME} OUT SRC)
usFunctionGenerateBundleInit(TARGET ${LIBRARY_NAME} OUT SRC)

add_library(${LIBRARY_NAME} SHARED ${SRC})
add_dependencies(${LIBRARY_NAME} base)
# Enable more warnings for our C++ code 
target_compile_options(${LIBRARY_NAME} PRIVATE -Wall -Wextra -Wpedantic -Wno-unused-parameter)

set(_bundle_name xacc_quac_accelerator)
set_target_properties(${LIBRARY_NAME} PROPERTIES
    # This is required for every bundle
    COMPILE_DEFINITIONS US_BUNDLE_NAME=${_bundle_name}
    # This is for convenience, used by other CMake functions
    US_BUNDLE_NAME ${_bundle_name}
    )

# Embed meta-data from a manifest.json file
usFunctionEmbedResources(TARGET ${LIBRARY_NAME}
    WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
    FILES
    manifest.json
    )

message ("XACC_QUAC_BASE_INCLUDE_DIR = ${XACC_QUAC_BASE_INCLUDE_DIR}")

target_include_directories(${LIBRARY_NAME} PUBLIC   . 
                                                    ./Accelerator
                                                    ./Utils
                                                    ./Executor
                                                    ${QUAC_INCLUDE_DIRS}
                                                    ${XACC_QUAC_BASE_INCLUDE_DIR}
                                                    ./tpls/cereal/include)


target_link_libraries(${LIBRARY_NAME} PUBLIC xacc::xacc xacc::quantum_gate QuaC base)
set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-shared")

set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")

xacc_configure_plugin_rpath(${LIBRARY_NAME})
set (CMAKE_INSTALL_PREFIX "${XACC_ROOT}" CACHE PATH "default install path" FORCE )

install(TARGETS ${LIBRARY_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/plugins)
add_subdirectory(tests)

