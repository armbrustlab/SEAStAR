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

# CMake file for ref_select
cmake_minimum_required (VERSION 2.8)

project (ref_select C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(ref_select ref_select.c)

target_link_libraries(ref_select ${COMMON_LIBS})

if(CMAKE_GENERATOR STREQUAL "Xcode")
    add_custom_command(
        TARGET ref_select
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ref_select
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:ref_select> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ref_select
)
endif(CMAKE_GENERATOR STREQUAL "Xcode")

enable_testing()

add_test(NAME ref_select_runs COMMAND ref_select --version)

set_tests_properties(ref_select_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 3)

add_test(NAME ref_select_test_data COMMAND bash -c "$<TARGET_FILE:ref_select> --per_base --mp_inserts -m --mp_mate_cnt=10 -r ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.single.fastq.gz -r ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read2.fastq.gz -r ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read1.fastq.gz -q -c ${CMAKE_SOURCE_DIR}/test_data/lambda_asm_catalog.txt.gz ${CMAKE_SOURCE_DIR}/test_data/lambda_reads.read2.fastq.gz_lambda_asm.fasta.sam.gz ${CMAKE_SOURCE_DIR}/test_data/lambda_reads.read1.fastq.gz_lambda_asm.fasta.sam.gz ${CMAKE_SOURCE_DIR}/test_data/lambda_reads.single.fastq.gz_lambda_asm.fasta.sam.gz | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm.json | sed 's/\"SEASTAR_tool\":[^}]*},//' | openssl md5")

add_test(NAME ref_select_test_data_no_cat COMMAND bash -c "$<TARGET_FILE:ref_select> --per_base --mp_inserts -m --mp_mate_cnt=10 -r ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.single.fastq.gz -r ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read2.fastq.gz -r ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read1.fastq.gz -q ${CMAKE_SOURCE_DIR}/test_data/lambda_reads.read2.fastq.gz_lambda_asm.fasta.sam.gz ${CMAKE_SOURCE_DIR}/test_data/lambda_reads.read1.fastq.gz_lambda_asm.fasta.sam.gz ${CMAKE_SOURCE_DIR}/test_data/lambda_reads.single.fastq.gz_lambda_asm.fasta.sam.gz | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm_no_cat.json | sed 's/\"SEASTAR_tool\":[^}]*},//' | openssl md5")

set_tests_properties(ref_select_test_data PROPERTIES PASS_REGULAR_EXPRESSION "66fdae801b5415eaf719b66b494505c6" TIMEOUT 30)
set_tests_properties(ref_select_test_data_no_cat PROPERTIES PASS_REGULAR_EXPRESSION "2a19d34c68bbcbde5e9c9622a0e5d895" TIMEOUT 30)

