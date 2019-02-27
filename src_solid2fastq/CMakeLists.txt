# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2014 University of Washington.
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
cmake_minimum_required (VERSION 2.8.5)

project (solid2fastq C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(solid2fastq solid2fastq.c)

target_link_libraries(solid2fastq ${COMMON_LIBS})

if(CMAKE_GENERATOR STREQUAL "Xcode")
    add_custom_command(
        TARGET solid2fastq
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/solid2fastq
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:solid2fastq> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/solid2fastq
)
endif(CMAKE_GENERATOR STREQUAL "Xcode")

enable_testing()

find_program(GAWK_EXEC "gawk")

add_test(NAME solid2fastq_runs COMMAND solid2fastq --version)

set_tests_properties(solid2fastq_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 3)

add_test(NAME solid2fastq_test_data COMMAND bash -c "$<TARGET_FILE:solid2fastq> -z --prefix=lambda_tst ${CMAKE_SOURCE_DIR}/test_data/lambda_reads ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads.{read1,read2}.fastq.gz | openssl md5")

set_tests_properties(solid2fastq_test_data PROPERTIES PASS_REGULAR_EXPRESSION "97925a0daa5e9d31a4b9c746cd0c1ba7" TIMEOUT 30)

# Test suppression of /1 and /2 read ID suffixes
add_test(NAME solid2fastq_test_no_suffix COMMAND bash -c "$<TARGET_FILE:solid2fastq> -z --prefix=lambda_tst --no_suffix ${CMAKE_SOURCE_DIR}/test_data/lambda_reads ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_no_suffix; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_no_suffix.{read1,read2}.fastq.gz | openssl md5")

set_tests_properties(solid2fastq_test_no_suffix PROPERTIES PASS_REGULAR_EXPRESSION "791734e70eb0c93f6db4fa21daf2a263" TIMEOUT 30)

# Create unsynced input fastq to test solid2fastq sorting.  These reads are 
# colorspace, but it's ok to treat them as nucleotide for testing.  To unsync
# R3 (read1) and F3 (read2) skip every fourth R3 read and every third F3 read
add_test(NAME solid2fastq_nt_input_test_data COMMAND bash -c "gzip -dc ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads.read1.fastq.gz | ${GAWK_EXEC} -v EVERY=4 -v HEADSTART=@ -v F3=_F3 -v R3=_R3 'line == 4 {x=0; line=0} NR % 4 == 1 {cnt++;sub(/^[^:]+:/, HEADSTART);sub(/\\/2$/, F3);sub(/\\/1$/, R3)} cnt % EVERY == 0 {x=1} x == 1 {line++; next} {print}' | gzip >${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_nt_unsync_default_R3.fastq.gz; gzip -dc ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads.read2.fastq.gz | ${GAWK_EXEC} -v EVERY=3 -v HEADSTART=@ -v F3=_F3 -v R3=_R3 'line == 4 {x=0; line=0} NR % 4 == 1 {cnt++;sub(/^[^:]+:/, HEADSTART);sub(/\\/1$/, R3);sub(/\\/2$/, F3)} cnt % EVERY == 0 {x=1} x == 1 {line++; next} {print}' | gzip >${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_nt_unsync_default_F3.fastq.gz; $<TARGET_FILE:solid2fastq> -b default -z --prefix=lambda_tst ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_nt_unsync ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_nt; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_nt.{read1,read2,single}.fastq.gz | openssl md5")

set_tests_properties(solid2fastq_nt_input_test_data PROPERTIES PASS_REGULAR_EXPRESSION "5c85d2f526c9cc9fb28068e6d33b9283" TIMEOUT 30)
