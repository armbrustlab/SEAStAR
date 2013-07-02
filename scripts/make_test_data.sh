#!/bin/bash -e
# --------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2013 University of Washington.
#
# Authors:
# Vaughn Iverson
# vsi@uw.edu
# --------------------------------------------------------------------------- #
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
# --------------------------------------------------------------------------- #

# This file is from SEAStAR version: SS_BUILD_VERSION

# Re-create the lambda phage assembled contigs, catalog file and alignment
# from the raw lambda reads as filtered through the CMake test harness.
# To run successfully, this script must be run from within a new subdirectory
# under the SEASTAR source dir, and it must be run after a 'make test' has been 
# run to generate trimmed reads in the output_test_data subdirectory.

# The outputs of this script are new versions of:
#    lambda_asm_catalog.txt.gz
#    lambda_reads.read1.fastq.gz_lambda_asm.fasta.sam.gz
#    lambda_reads.read2.fastq.gz_lambda_asm.fasta.sam.gz
#    lambda_reads.single.fastq.gz_lambda_asm.fasta.sam.gz

# These files are built here and distributed pre-computed with SEAStAR 
# so that an arbitrary build machine does not need to precisely recreate the 
# actions of Velvet and BWA as performed below in order for components 
# that depend on these files to successfully pass their CTest tests.

# Run with colorspace + OpenMP build of Velvet version 1.2.08
velveth_de lambda_asm/ 19 -fastq.gz -shortPaired ../output_test_data/lambda_reads_d90.mates.fastq.gz -short ../output_test_data/lambda_reads_d90.single.fastq.gz > lambda_asm.velveth_de.log 2>&1
velvetg_de lambda_asm/ -scaffolding no -read_trkg no -ins_length auto -ins_length_sd auto -exp_cov 50 -cov_cutoff 5 -min_contig_lgth 75 > lambda_asm.velvetg_de.log 2>&1

csfasta2ntfasta.awk -v catalog_file=lambda_asm_catalog.txt lambda_asm/contigs.fa > lambda_contigs.fna 
gzip -c lambda_asm_catalog.txt > lambda_asm_catalog.txt.gz

# Run with BWA version 0.5.9-r26-dev
bwa index -a is -c lambda_contigs.fna 2> lambda.index.log
bwa aln -c -n 0.001 -l 18 lambda_contigs.fna ../output_test_data/lambda_reads_d90.read1.fastq.gz 2> lambda_trim.read1.aln.log | bwa samse -n 1000000 lambda_contigs.fna - ../output_test_data/lambda_reads_d90.read1.fastq.gz 2> lambda_trim.read1.samse.log | gzip -c > lambda_reads.read1.fastq.gz_lambda_asm.fasta.sam.gz
bwa aln -c -n 0.001 -l 18 lambda_contigs.fna ../output_test_data/lambda_reads_d90.read2.fastq.gz 2> lambda_trim.read2.aln.log | bwa samse -n 1000000 lambda_contigs.fna - ../output_test_data/lambda_reads_d90.read2.fastq.gz 2> lambda_trim.read2.samse.log | gzip -c > lambda_reads.read2.fastq.gz_lambda_asm.fasta.sam.gz
bwa aln -c -n 0.001 -l 18 lambda_contigs.fna ../output_test_data/lambda_reads_d90.mates.fastq.gz 2> lambda_trim.single.aln.log | bwa samse -n 1000000 lambda_contigs.fna - ../output_test_data/lambda_reads_d90.mates.fastq.gz 2> lambda_trim.single.samse.log | gzip -c > lambda_reads.single.fastq.gz_lambda_asm.fasta.sam.gz

# Uncomment line below to copy results of this script to the test_data directory
# replacing the previous version of these files.  Note, this will likely require
# updating the MD5 result hashes in the various downstream CTest cases for 
# 'make test' to run successfully again (e.g. because Velvet doesn't produce
# precisely the same contigs each time it is run...)

# cp lambda_asm_catalog.txt.gz lambda_reads.read1.fastq.gz_lambda_asm.fasta.sam.gz lambda_reads.read2.fastq.gz_lambda_asm.fasta.sam.gz lambda_reads.single.fastq.gz_lambda_asm.fasta.sam.gz ../test_data

