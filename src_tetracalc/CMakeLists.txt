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

# CMake file for tetracalc
cmake_minimum_required (VERSION 2.8.5)

project (tetracalc C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_executable(tetracalc tetracalc.c tune_tabs.h)

target_link_libraries(tetracalc ${COMMON_LIBS})

if(CMAKE_GENERATOR STREQUAL "Xcode")
    add_custom_command(
        TARGET tetracalc
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tetracalc
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:tetracalc> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tetracalc
)
endif(CMAKE_GENERATOR STREQUAL "Xcode")

enable_testing()

add_test(NAME tetracalc_runs COMMAND tetracalc --version)

add_test(NAME tetracalc_test_data COMMAND bash -c "$<TARGET_FILE:tetracalc> --merge_tar=0.99 --chunk=50000 ${CMAKE_SOURCE_DIR}/test_data/NC_008209.1.fna ${CMAKE_SOURCE_DIR}/test_data/NC_008044.1.fna ${CMAKE_SOURCE_DIR}/test_data/NC_007205.1.fna | tee ${CMAKE_BINARY_DIR}/output_test_data/3_alphas_clustered.json | openssl md5")

set_tests_properties(tetracalc_runs PROPERTIES PASS_REGULAR_EXPRESSION "version: ${GIT_VERSION}" TIMEOUT 3)

set_tests_properties(tetracalc_test_data PROPERTIES PASS_REGULAR_EXPRESSION "ab56ec052736ca5c2536028c63b50bf6" TIMEOUT 90)
