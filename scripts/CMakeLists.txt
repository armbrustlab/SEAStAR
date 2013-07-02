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

# CMake file for scripts 
cmake_minimum_required (VERSION 2.8)

# Directory for all final scripts
set(SCRIPT_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

find_program(SED_EXEC "sed")
find_program(GAWK_EXEC "gawk")

# N.B. At some point in the future when CMake adds a command to alter file
# permissions the chmod commands used below should be replaced to stay platform
# independent.
# Copy AWK scripts to final script directory
file(GLOB AWK_SOURCES_INPUT RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.awk)
set(AWK_SOURCES_OUTPUT)
foreach(IN_FILE ${AWK_SOURCES_INPUT})
    set(OUT_FILE ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/${IN_FILE})
    list(APPEND AWK_SOURCES_OUTPUT ${OUT_FILE})
    add_custom_command(
        OUTPUT ${OUT_FILE}
        COMMAND bash -c 'echo \"\#!${GAWK_EXEC} -f\" > ${OUT_FILE}'
        COMMAND bash -c '${SED_EXEC} "s/SS_BUILD_VERSION/${GIT_VERSION}/g" ${CMAKE_CURRENT_SOURCE_DIR}/${IN_FILE} >> ${OUT_FILE}'
        COMMAND chmod 555 ${OUT_FILE}
        DEPENDS ${IN_FILE}
    )
endforeach()

# Copy BASH scripts to final script directory
file(GLOB SH_SOURCES_INPUT RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.sh)
set(SH_SOURCES_OUTPUT)
foreach(IN_FILE ${SH_SOURCES_INPUT})
    string(REPLACE ".sh" "" OUT_FILE ${IN_FILE})
    set(OUT_PATH ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/${OUT_FILE})
    list(APPEND SH_SOURCES_OUTPUT ${OUT_PATH})
    add_custom_command(
        OUTPUT ${OUT_PATH}
        COMMAND ${CMAKE_COMMAND} -E remove -f ${OUT_PATH}
        COMMAND bash -c '${SED_EXEC} "s/SS_BUILD_VERSION/${GIT_VERSION}/g" ${CMAKE_CURRENT_SOURCE_DIR}/${IN_FILE} > ${OUT_PATH}'
        COMMAND chmod 555 ${OUT_PATH}
        DEPENDS ${IN_FILE}
    )
endforeach()

# Optionally Transcompile coffeescript to javascript in the scripts directory
find_program(COFFEE_EXEC "coffee")
find_program(NODE_EXEC "node")

file(GLOB COFFEE_SOURCES_INPUT RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.coffee)
set(COFFEE_SOURCES_OUTPUT)
if (COFFEE_EXEC)
    foreach(IN_FILE ${COFFEE_SOURCES_INPUT})
        get_filename_component(BASENAME ${IN_FILE} NAME_WE)
        set(OUT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${BASENAME}.js)
        list(APPEND COFFEE_SOURCES_OUTPUT ${OUT_FILE})
        add_custom_command(
            OUTPUT ${OUT_FILE}
            COMMAND ${CMAKE_COMMAND} -E remove -f ${OUT_FILE}
            COMMAND ${NODE_EXEC} ${COFFEE_EXEC} -c ${CMAKE_CURRENT_SOURCE_DIR}/${IN_FILE}
            COMMAND chmod 444 ${OUT_FILE}
            DEPENDS ${IN_FILE}
        )
    endforeach()
endif (COFFEE_EXEC)

# Copy javascript files to bin directory
set(JS_SOURCES_OUTPUT)
foreach(COFFEE_SRC_FILE ${COFFEE_SOURCES_INPUT})
    get_filename_component(BASENAME ${COFFEE_SRC_FILE} NAME_WE)
    set(IN_FILE ${BASENAME}.js)
    set(OUT_FILE ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/${IN_FILE})
    list(APPEND JS_SOURCES_OUTPUT ${OUT_FILE})
    add_custom_command(
        OUTPUT ${OUT_FILE}
        COMMAND ${CMAKE_COMMAND} -E remove -f ${OUT_FILE}
        COMMAND bash -c '${SED_EXEC} "s/SS_BUILD_VERSION/${GIT_VERSION}/g" ${CMAKE_CURRENT_SOURCE_DIR}/${IN_FILE} > ${OUT_FILE}'
        COMMAND chmod 444 ${OUT_FILE}
        DEPENDS ${IN_FILE}
    )
endforeach()

# Copy contents of node_modules subdirectory to bin/ 
set (NODE_MODS_OUTPUT ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/node_modules)
set (NODE_MODS_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/node_modules)
add_custom_command(
     OUTPUT ${NODE_MODS_OUTPUT}
     COMMAND ${CMAKE_COMMAND} -E copy_directory ${NODE_MODS_INPUT} ${NODE_MODS_OUTPUT}
     DEPENDS ${NODE_MODS_INPUT}
)

# Copy RDP_viz directory to the final script directory
set(RDP_VIZ_FILES ${CMAKE_CURRENT_SOURCE_DIR}/RDP_viz/raphael-2.1.0.js ${CMAKE_CURRENT_SOURCE_DIR}/RDP_viz/tree_bars.html ${CMAKE_CURRENT_SOURCE_DIR}/RDP_viz/tree_bars.js)
set(OUT_DIR ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/RDP_viz)
set(RDP_VIZ_SOURCES_OUTPUT)
foreach(IN_FILE ${RDP_VIZ_FILES})
    get_filename_component(BASENAME ${IN_FILE} NAME)
    set(OUT_FILE ${OUT_DIR}/${BASENAME})
    list(APPEND RDP_VIZ_SOURCES_OUTPUT ${OUT_FILE})
    add_custom_command(
        OUTPUT ${OUT_FILE}
        COMMAND ${CMAKE_COMMAND} -E copy ${IN_FILE} ${OUT_FILE}
        COMMAND chmod 444 ${OUT_FILE}
        DEPENDS ${IN_FILE}
    )
endforeach()

add_custom_target(scripts ALL DEPENDS ${NODE_MODS_OUTPUT} ${JS_SOURCES_OUTPUT} ${COFFEE_SOURCES_OUTPUT} ${AWK_SOURCES_OUTPUT} ${SH_SOURCES_OUTPUT} ${RDP_VIZ_SOURCES_OUTPUT})

enable_testing()

add_test(NAME fastq2csfasta_test_data COMMAND bash -c "gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read2.fastq.gz | ${GAWK_EXEC} -f ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/fastq2csfasta.awk > ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read2.csfasta; gunzip -c ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read1.fastq.gz | ${GAWK_EXEC} -f ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/fastq2csfasta.awk | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read1.csfasta | openssl md5")
set_tests_properties(fastq2csfasta_test_data PROPERTIES PASS_REGULAR_EXPRESSION "9d93279bbd578cbb88d74755d20ee357" TIMEOUT 30)

add_test(NAME csfasta2fasta_test_data COMMAND bash -c "${GAWK_EXEC} -f ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/csfasta2fasta.awk ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read1.csfasta | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read1.fasta | openssl md5")
set_tests_properties(csfasta2fasta_test_data PROPERTIES PASS_REGULAR_EXPRESSION "b3138e9d395a2ac59e8861b617f8f184" TIMEOUT 30)

add_test(NAME csfasta2fasta_paired_test_data COMMAND bash -c "${GAWK_EXEC} -f ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/csfasta2fasta.awk ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read1.csfasta ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.read2.csfasta | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_reads_d90.mates.fasta | openssl md5")
set_tests_properties(csfasta2fasta_paired_test_data PROPERTIES PASS_REGULAR_EXPRESSION "587dcdfa4b69f2422172235fb63faf2a" TIMEOUT 30)

add_test(NAME csfasta2ntfasta_test_data COMMAND bash -c "${NODE_EXEC} ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/graph_ops.js ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm.json FASTA | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm_recon.fna | ${GAWK_EXEC} -v numstarts=4 -f ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/csfasta2ntfasta.awk | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm_recon_4starts.fasta | openssl md5") 
set_tests_properties(csfasta2ntfasta_test_data PROPERTIES PASS_REGULAR_EXPRESSION "744d427f8f9748f990b4b1a87d0fb1f2" TIMEOUT 30)

add_test(NAME graph_ops_scaffolding_test_data COMMAND bash -c "${NODE_EXEC} ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/graph_ops.js ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm.json MST PLUCK PRUNE PUSH INSERT SCAFF DUMP '{\"file\":\"${CMAKE_BINARY_DIR}/output_test_data/lambda_asm_scaffs.json\"}' FASTA '{\"scaff\":true}' | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm_scaffs.fna | openssl md5") 
set_tests_properties(graph_ops_scaffolding_test_data PROPERTIES ENVIRONMENT "PATH=${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}:$ENV{PATH}")
set_tests_properties(graph_ops_scaffolding_test_data PROPERTIES PASS_REGULAR_EXPRESSION "7760c7722f392ffe3946d5022397e062" TIMEOUT 30)

add_test(NAME fix_fasta_line_lengths_test_data COMMAND bash -c "${GAWK_EXEC} --re-interval -v minlen=500 -v maxlen=4000 -v ll=80 -f ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/fix_fasta_line_lengths.awk ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm_recon.fna | tee ${CMAKE_BINARY_DIR}/output_test_data/lambda_asm_recon_fixed_ll.fna | openssl md5")
set_tests_properties(fix_fasta_line_lengths_test_data PROPERTIES PASS_REGULAR_EXPRESSION "69d03ede642d6b5bf005199f4e6ec008" TIMEOUT 30)

add_test(NAME RDP_train_to_heirarchy_test_data COMMAND bash -c "echo -e '0*Root*-1*0*norank\\n1*Bacteria*0*1*domain\\n2*Acidobacteria*1*2*phylum\\n3*Holophagae*2*3*class\\n4*Acanthopleuribacterales*3*4*order\\n5*Acanthopleuribacteraceae*4*5*family\\n6*Acanthopleuribacter*5*6*genus' | ${NODE_EXEC} ${SCRIPT_RUNTIME_OUTPUT_DIRECTORY}/RDP_train_to_tree.js | openssl md5")
set_tests_properties(RDP_train_to_heirarchy_test_data PROPERTIES PASS_REGULAR_EXPRESSION "126f51df16a689518a951e13f3983cf8" TIMEOUT 30)
