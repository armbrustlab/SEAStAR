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

# CMake file for solid2fastq
cmake_minimum_required (VERSION 2.6)

project (solid2fastq C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(solid2fastq solid2fastq.c)

target_link_libraries(solid2fastq ${COMMON_LIBS})

enable_testing()

add_test(NAME solid2fastq_runs COMMAND solid2fastq --version)

set_tests_properties(solid2fastq_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 1)

add_test(NAME solid2fastq_test_data COMMAND solid2fastq -z --prefix=lambda_tst ${CMAKE_SOURCE_DIR}/test_data/lambda_reads ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads)

set_tests_properties(solid2fastq_test_data PROPERTIES PASS_REGULAR_EXPRESSION "Matepairs: 80374, Matepaired Reads: 160748.*Total singlets: 0.*Total Reads Processed: 160748" TIMEOUT 15)

