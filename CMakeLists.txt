# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2013 University of Washington.
#
# Authors:
# Vaughn Iverson
# vsi@uw.edu
# -------------------------------------------------------------------------- #
# This file is part of SEAStAR.
#
# SEAStAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SEAStAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SEAStAR.  If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------- #

# CMake file for full SEASTAR build 
cmake_minimum_required (VERSION 2.8)

# Debug option
option(DEBUG "debug" OFF)

IF(DEBUG)
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING "debug build" FORCE)
ENDIF(DEBUG)

IF(NOT DEBUG)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "default to release build" FORCE)
ENDIF(NOT DEBUG)

MESSAGE("Cmake build type set to ${CMAKE_BUILD_TYPE}")

project (SEASTAR C)

set (BUILD_SHARED_LIBS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
include_directories(${CMAKE_SOURCE_DIR}/src_seastar_shared)

set(UTHASH_INCLUDES ${CMAKE_SOURCE_DIR}/uthash_includes)
include_directories(${UTHASH_INCLUDES})

# FORTIFY_SOURCE is This line is disabled to fix an apparent bug in LVVM / OpenMP support
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_GNU_SOURCE -D_FORTIFY_SOURCE=0 -w -std=c99")

if(NOT CMAKE_GENERATOR STREQUAL "Xcode")

	include(FindOpenMP)

	if(OPENMP_FOUND)

		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

		message(STATUS "System is ${CMAKE_SYSTEM}.")
		if((${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") AND NOT DARWIN_OPENMP_VERIFIED)
			message(STATUS "Testing OpenMP support in compiler for known bugs. This may take a few seconds...")
			try_run(OPENMP_TEST_RESULT OPENMP_TEST_COMPILE_RESULT ${CMAKE_CURRENT_BINARY_DIR}/openmp_test ${CMAKE_SOURCE_DIR}/openmp_test/test.c RUN_OUTPUT_VARIABLE OPENMP_TEST_OUTPUT)
			if (NOT ${OPENMP_TEST_COMPILE_RESULT})
				message(FATAL_ERROR "OpenMP test failed to compile.")
			endif (NOT ${OPENMP_TEST_COMPILE_RESULT})
			if(OPENMP_TEST_RESULT)
				message(WARNING "An error was encountered while testing OpenMP support in the available compiler: ${CMAKE_C_COMPILER} -- ${OPENMP_TEST_OUTPUT}")
				string(REGEX MATCH "[Aa]bort" OPENMP_TEST_ABORT ${OPENMP_TEST_OUTPUT})
				if(OPENMP_TEST_ABORT)
					message(FATAL_ERROR "The default Apple provided compiler on your system has OpenMP bugs.\n You need to install gcc (version 4.4 or later) from MacPorts. See README file for more details.")
				endif(OPENMP_TEST_ABORT)
			else(OPENMP_TEST_RESULT)
				SET(DARWIN_OPENMP_VERIFIED TRUE CACHE BOOL "Darwin OpenMP passed?")
			endif(OPENMP_TEST_RESULT)
		endif((${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") AND NOT DARWIN_OPENMP_VERIFIED)
	else(OPENMP_FOUND)
		message (FATAL_ERROR "OpenMP support is required")
	endif(OPENMP_FOUND)
else(NOT CMAKE_GENERATOR STREQUAL "Xcode")
	if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
		message(WARNING "XCode compilers as of XCode 4.2.1 on OS X Lion 10.7.2 contain known bugs making OpenMP code intermittantly unstable. If \"Abort trap: 6\" errors are encountered, please recompile using the Cmake makefile generator.")
	endif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
	set(CMAKE_XCODE_ATTRIBUTE_GCC_VERSION "com.apple.compilers.llvmgcc42")
	set(CMAKE_XCODE_ATTRIBUTE_ENABLE_OPENMP_SUPPORT YES)
endif(NOT CMAKE_GENERATOR STREQUAL "Xcode")

execute_process(COMMAND bash -c "cd ${CMAKE_SOURCE_DIR}; git describe --match 'v*'" OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE GIT_VERSION RESULT_VARIABLE GIT_RESULT)

# Test if the git command above failed, if so, set a default version string
if (GIT_VERSION)
	add_definitions("-D SS_BUILD_VERSION='\"${GIT_VERSION}\"'")
else (GIT_VERSION)
	message(WARNING "SEAStAR version not detected from git. ${GIT_VERSION}") 
	set(GIT_VERSION "v0.4.9")
	add_definitions("-D SS_BUILD_VERSION='\"${GIT_VERSION}\"'")
endif (GIT_VERSION)

enable_testing()

set(Argtable_DIR "argtable2-13")
include_directories(${CMAKE_SOURCE_DIR}/${Argtable_DIR}/src)
add_subdirectory(${Argtable_DIR})

set(Zlib_DIR "zlib-1.2.8")
include_directories(${CMAKE_SOURCE_DIR}/${Zlib_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/${Zlib_DIR})
add_subdirectory(${Zlib_DIR})

add_subdirectory(src_seastar_shared)

set(COMMON_LIBS "m" "seastar_shared" "argtable2" "z")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_test_data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# Note, it's important to keep these in this order for CTest because 
# the outputs of some tool's tests are used as the inputs of others.
add_subdirectory(src_solid2fastq)
add_subdirectory(src_fastq_nodup)
add_subdirectory(src_trimfastq)
add_subdirectory(src_samplefastq)
add_subdirectory(src_ref_select)
add_subdirectory(src_tetracalc)
add_subdirectory(scripts)

