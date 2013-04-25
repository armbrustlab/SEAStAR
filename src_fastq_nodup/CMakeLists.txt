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

# CMake file for fastq_nodup
cmake_minimum_required (VERSION 2.8)

project (fastq_nodup C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(fastq_nodup fastq_nodup.c)

target_link_libraries(fastq_nodup ${COMMON_LIBS})

if(CMAKE_GENERATOR STREQUAL "Xcode")
    add_custom_command(
        TARGET fastq_nodup
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fastq_nodup
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:fastq_nodup> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fastq_nodup
)
endif(CMAKE_GENERATOR STREQUAL "Xcode")

enable_testing()

add_test(NAME fastq_nodup_runs COMMAND fastq_nodup --version)

set_tests_properties(fastq_nodup_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 3)

add_test(NAME fastq_nodup_test_data COMMAND bash -c "$<TARGET_FILE:fastq_nodup> -z --prefix=lambda_tst.d -l 13 -d 2 -e 3 ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d.{read1,read2}.fastq.gz | openssl md5 ")

set_tests_properties(fastq_nodup_test_data PROPERTIES PASS_REGULAR_EXPRESSION "c4c8331dd1ffbe372ffe326617854ade" TIMEOUT 90)


