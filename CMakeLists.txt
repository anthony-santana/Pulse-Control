#***********************************************************************************
# Copyright (c) 2020, UT-Battelle
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the xacc nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#**********************************************************************************/
cmake_minimum_required (VERSION 3.13)

project (QuaC LANGUAGES C CXX)
enable_testing()

set (CMAKE_STANDARD_REQUIRED ON)
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

set (LIBRARY_NAME QuaC)

if (NOT CMAKE_BUILD_TYPE)
  # Default build type is "Release" if not specified 
  message (STATUS "CMAKE_BUILD_TYPE was not specified. Default to \"Release\".")

  set (CMAKE_BUILD_TYPE "Release" CACHE STRING
	   "Choose the type of build, options are: Debug, Release, RelWithDebInfo, MinSizeRel"
	   FORCE
      )
else()
  message (STATUS "Build Type: ${CMAKE_BUILD_TYPE}.")
endif()

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

set(QUAC_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/src)
set (XACC_QUAC_BASE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/base)

add_subdirectory(base)

if (PETSC_DIR)
  set(MPI_C_COMPILER "${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc")
  set(MPI_CXX_COMPILER "${PETSC_DIR}/${PETSC_ARCH}/bin/mpic++")
  set(MPIEXEC_EXECUTABLE "${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec")
  find_package(MPI REQUIRED)
  message(STATUS "Found a suitable MPI compiler ${MPI_CXX_COMPILER}.")
  message(STATUS "Compiler vendor is [${CMAKE_CXX_COMPILER_ID}]")
  message(STATUS "To run MPI programs, use ${MPIEXEC_EXECUTABLE}")
  message(STATUS "Include path: ${MPI_CXX_INCLUDE_DIRS}")
  message(STATUS "Compile flags: ${MPI_CXX_COMPILE_FLAGS}")
  message(STATUS "Link flags: ${MPI_CXX_LINK_FLAGS}")
  message(STATUS "Libraries: ${MPI_CXX_LIBRARIES}")

  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_libraries(${MPI_CXX_LIBRARIES})
  link_libraries(-Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib)
  # TODO: parse these from PETSc Makefile config var 
  link_libraries(-lpetsc -lflapack -lfblas)
  link_libraries(-lstdc++ -ldl -lmpifort -lmpi -lgcc_s -lquadmath)
  link_libraries(-lpthread -lm -lstdc++ -ldl -lmpicxx -lgfortran)

  add_library(QuaC STATIC
    src/dm_utilities.c
    src/kron.c
    src/error_correction.c
    src/interface_xacc_ir.c
    src/operators.c
    src/quac.c
    src/quantum_gates.c
    src/solver.c
  )
  target_include_directories(${LIBRARY_NAME} PUBLIC 
    ${PETSC_DIR}/include 
    ${PETSC_DIR}/${PETSC_ARCH}/include
    ${XACC_QUAC_BASE_INCLUDE_DIR}
  )
  target_compile_options(${LIBRARY_NAME} PRIVATE -fPIC)
  add_dependencies(${LIBRARY_NAME} base)
    
    add_subdirectory(xacc)
    add_subdirectory(xacc_examples)
else()
    message (FATAL_ERROR "Please provide path to PETSc.")
endif()