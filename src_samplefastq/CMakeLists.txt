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

# CMake file for samplefastq
cmake_minimum_required (VERSION 2.8)

project (samplefastq C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(samplefastq samplefastq.c)

target_link_libraries(samplefastq ${COMMON_LIBS})

if(CMAKE_GENERATOR STREQUAL "Xcode")
    add_custom_command(
        TARGET samplefastq
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/samplefastq
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:samplefastq> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/samplefastq
)
endif(CMAKE_GENERATOR STREQUAL "Xcode")

enable_testing()

add_test(NAME samplefastq_runs COMMAND samplefastq --version)

set_tests_properties(samplefastq_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 3)

add_test(NAME samplefastq_test_data COMMAND bash -c "$<TARGET_FILE:samplefastq> -f 0.1 -z ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90 ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_samp; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_samp.{read1,read2,single}.fastq.gz | openssl md5")

set_tests_properties(samplefastq_test_data PROPERTIES PASS_REGULAR_EXPRESSION "e1a7c29691a8542449dcf06c575c3e02" TIMEOUT 30)

