
find_package(XACC REQUIRED)

set (LIBRARY_NAME xacc_quac)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Enable mocking for testing
# We will disable this eventually...
add_compile_definitions(ENABLE_MOCKING)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Mock)

file (GLOB HEADERS *.hpp)
file (GLOB SRC *.cpp Accelerator/*.cpp Utils/*.cpp)

usFunctionGetResourceSource(TARGET ${LIBRARY_NAME} OUT SRC)
usFunctionGenerateBundleInit(TARGET ${LIBRARY_NAME} OUT SRC)

add_library(${LIBRARY_NAME} SHARED ${SRC})
add_dependencies(${LIBRARY_NAME} base)

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
                                                    ${QUAC_INCLUDE_DIRS}
                                                    ${XACC_QUAC_BASE_INCLUDE_DIR})


target_link_libraries(${LIBRARY_NAME} PUBLIC xacc::xacc xacc::quantum_gate ${CMAKE_BINARY_DIR}/libQuaC.so base)
set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-shared")

set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")

xacc_configure_plugin_rpath(${LIBRARY_NAME})
set (CMAKE_INSTALL_PREFIX "${XACC_ROOT}" CACHE PATH "default install path" FORCE )

install(TARGETS ${LIBRARY_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/plugins)
install(FILES ${CMAKE_BINARY_DIR}/libQuaC.so DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)