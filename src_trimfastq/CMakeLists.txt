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

# CMake file for trimfastq
cmake_minimum_required (VERSION 2.6)

project (trimfastq C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(trimfastq trimfastq.c)

target_link_libraries(trimfastq ${COMMON_LIBS})

enable_testing()

add_test(NAME trimfastq_runs COMMAND trimfastq --version)

set_tests_properties(trimfastq_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 1)

add_test(NAME trimfastq_test_data COMMAND trimfastq -z --prefix=lambda_tst.d50 -p 0.5 -l 34 -m 34 -e 3.0 --add_len ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d50)

set_tests_properties(trimfastq_test_data PROPERTIES PASS_REGULAR_EXPRESSION "Low complexity reads discarded: Read1: 157, Read2: 241, Singlets: 0 0.*Matepairs: Before: 79270, After: 61966.*Singlets: Before: 0 0 After: 17304 0.*Read1 singlets: 7569, Read2 singlets: 9735, Original singlets: 0 0.*Total Reads Processed: 158540, Reads retained: 141236" TIMEOUT 15)
