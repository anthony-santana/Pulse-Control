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

set(QUAC_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/src)
set (XACC_QUAC_BASE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/base)

add_subdirectory(base)

if (PETSC_DIR)
  add_custom_target (libQuaC ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/libQuaC.so)
      add_custom_command (
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libQuaC.so
        COMMAND make PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} clean
        COMMAND make PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} shared_lib
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/libQuaC.so ${CMAKE_CURRENT_BINARY_DIR}/libQuaC.so
        COMMAND rm -f ${CMAKE_SOURCE_DIR}/libQuaC.so
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        DEPENDS base
      )
    
    add_subdirectory(xacc)
    add_subdirectory(xacc_examples)
else()
    message (FATAL_ERROR "Please provide path to PETSc.")
endif()