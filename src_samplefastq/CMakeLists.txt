# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2012 University of Washington.
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

# CMake file for samplefastq
cmake_minimum_required (VERSION 2.6)

project (samplefastq C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(samplefastq samplefastq.c)

target_link_libraries(samplefastq ${COMMON_LIBS})

enable_testing()

add_test(NAME samplefastq_runs COMMAND samplefastq --version)

set_tests_properties(samplefastq_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 1)

add_test(NAME samplefastq_test_data COMMAND samplefastq -f 0.1 -z ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d50 ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_samp)

set_tests_properties(samplefastq_test_data PROPERTIES PASS_REGULAR_EXPRESSION "Matepairs: Before: 61966, After: 6196.*Singlets: Before: 17304 0 After: 1741 0.*Total Reads Processed: 141236, Reads retained: 14133" TIMEOUT 15)

