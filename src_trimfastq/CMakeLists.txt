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

# CMake file for trimfastq
cmake_minimum_required (VERSION 2.8)

project (trimfastq C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(trimfastq trimfastq.c)

target_link_libraries(trimfastq ${COMMON_LIBS})

if(CMAKE_GENERATOR STREQUAL "Xcode")
    add_custom_command(
        TARGET trimfastq
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/trimfastq
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:trimfastq> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/trimfastq
)
endif(CMAKE_GENERATOR STREQUAL "Xcode")

enable_testing()

find_program(GAWK_EXEC "gawk")

add_test(NAME trimfastq_runs COMMAND trimfastq --version)

set_tests_properties(trimfastq_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 3)

add_test(NAME trimfastq_test_data COMMAND bash -c "$<TARGET_FILE:trimfastq> -z --prefix=lambda_tst.d90 -p 0.9 -l 34 -m 34 -e 3.0 --add_len --mates_file ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.{mates,read1,read2,single}.fastq.gz | openssl md5")

set_tests_properties(trimfastq_test_data PROPERTIES PASS_REGULAR_EXPRESSION "61eb67ca4b9458ef260e9fe19b4a09fb" TIMEOUT 30)

add_test(NAME trimfastq_2singlets_test_data COMMAND bash -c "$<TARGET_FILE:trimfastq> -z --prefix=lambda_tst.d50s2 -p 0.5 -s -l 34 -m 34 -e 3.0 --add_len ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d50s2; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d50s2.{read1,read2,single1,single2}.fastq.gz | openssl md5")

set_tests_properties(trimfastq_2singlets_test_data PROPERTIES PASS_REGULAR_EXPRESSION "487481a8b4156040d98295080d5064fa" TIMEOUT 30)

add_test(NAME trimfastq_singlet_input_test_data COMMAND bash -c "$<TARGET_FILE:trimfastq> -z --prefix=lambda_tst.d70 -p 0.7 -l 34 -m 34 -e 3.0 --add_len ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d50s2 ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d70; (gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d70.{read1,read2}.fastq.gz; for f in ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d70.{single1,single2}.fastq.gz; do gunzip -c $f | ${GAWK_EXEC} 'NR %4 == 1 {cur=$0; next} {cur = cur $0} NR % 4 == 0 {x[i++] = cur} END {n = asort(x); for (j=1; j<=n; j++) {print x[j]}}'; done) | openssl md5")

set_tests_properties(trimfastq_singlet_input_test_data PROPERTIES PASS_REGULAR_EXPRESSION "df84824a396d016689ca7a0e35085a15" TIMEOUT 30)
