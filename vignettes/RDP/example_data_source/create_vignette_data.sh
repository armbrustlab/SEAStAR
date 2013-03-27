#!/bin/bash -e
# --------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2013 University of Washington.
# 
# Authors:
# Chris Berthiaume
# chrisbee@uw.edu
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
# This script will create the example data set to be used with the RDP vignette.
# It should be run from vignette/RDP (relative to the SEAStAR project source
# root).  Files created by this script will be placed in vignette/RDP.  Files
# used by this script to create the vignette example files can be found in 
# vignette/RDP/prep.
#
# Usage: ./create_vignette_data.sh
#
# Requirements:
# - PWD is in the SEAStAR repository directory vignette/RDP/prep
#
# - art_illumina is in the path (http://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
#   This software must be downloaded separately from this project
#
# - The RDP classifier software zip file is in vignette/RDP/prep (http://sourceforge.net/projects/rdp-classifier/)
#   It should match this shell glob pattern: rdp_classifier*.zip
#   This software must be downloaded separately from this project
#   This script has been tested with version 2.5
#
# - The RDP classifier raw training data zip file is in vignette/RDP/prep (http://sourceforge.net/projects/rdp-classifier/)
#   It should match this shell glob pattern: *rawtrainingdata.zip
#   This software must be downloaded separately from this project
#   This script has been tested with version 9
#
# Simulation reference sequences:
# art_SOLiD will be used to create 49bp SOLiD reads for the following RDP 
# sequences:
#
#            sample1       sample2
# ID         (cov,  %)     (cov,  %)   lineage
# S000000075 (100x, 37.0) (5x,   1.9)  Root; Bacteria; "Deinococcus-Thermus"; Deinococci; Thermales; Thermaceae; Thermus
# S000005187 (50x,  18.5) (10x,  3.7)  Root; Bacteria; "Nitrospira"; "Nitrospira"; "Nitrospirales"; "Nitrospiraceae"; Leptospirillum
# S000498962 (40x,  14.8) (15x,  5.6)  Root; Bacteria; "Proteobacteria"; Gammaproteobacteria; Oceanospirillales; SAR86; "SAR86 clade II"
# S000020432 (25x,  9.3)  (25x,  9.3)  Root; Bacteria; "Fusobacteria"; "Fusobacteria"; "Fusobacteriales"; "Fusobacteriaceae"; Cetobacterium
# S000134123 (25x,  9.3)  (25x,  9.3)  Root; Bacteria; "Deinococcus-Thermus"; Deinococci; Thermales; Thermaceae; Oceanithermus
# S000431293 (15x,  5.6)  (40x,  14.8) Root; Bacteria; "Deinococcus-Thermus"; Deinococci; Deinococcales; Deinococcaceae; Deinobacterium
# S000367885 (10x,  3.7)  (50x,  18.5) Root; Bacteria; "Actinobacteria"; Actinobacteria; Acidimicrobidae; Acidimicrobiales; "Acidimicrobineae"; Acidimicrobiaceae; Acidimicrobium
# S000390161 (5x,   1.9)  (100x, 37.0) Root; Archaea; "Korarchaeota"; "Candidatus Korarchaeum"
#
# Extra reference sequences included in the RDP_slim.fasta database file to
# demonstrate specificity of the selection process.
# S000458920 Root; Bacteria; "Deinococcus-Thermus"; Deinococci; Thermales; Thermaceae; Vulcanithermus
# S000011652 Root; Bacteria; "Nitrospira"; "Nitrospira"; "Nitrospirales"; "Nitrospiraceae"; Thermodesulfovibrio
# S000901535 Root; Bacteria; "Proteobacteria"; Gammaproteobacteria; Oceanospirillales; Oceanospirillaceae; Amphritea
# S000379106 Root; Bacteria; "Fusobacteria"; "Fusobacteria"; "Fusobacteriales"; "Fusobacteriaceae"; Psychrilyobacter
# S000443854 Root; Bacteria; "Actinobacteria"; Actinobacteria; Acidimicrobidae; Acidimicrobiales; "Acidimicrobineae"; Acidimicrobiaceae; Ferrimicrobium
# S000357634 Root; Archaea; "Euryarchaeota"; Halobacteria; Halobacteriales; Halobacteriaceae; Haladaptatus
# S000004750 Root; Bacteria; Firmicutes; Bacilli; Lactobacillales; Aerococcaceae; Facklamia

# Functions ------------------------------------------------------------------ #
# Simulate SOLiD reads with ART
# Arguments:
#   1) input fasta file path
#   2) fold coverage
#   3) aggregated output fastq file prefix
function simulate {
    if [[ $# -lt 3 ]]; then
        echo "Error: Not enough arguments in ${0}:simulate" >&2
        exit 1
    fi
    
    local infile_orig=$1
    local cov=$2
    local outpre=$(basename $3)  # output file prefix, remove dir if present
    
    infile_base=$(basename "$infile_orig")
    infile=""$outdir"/$infile_base"
    
    # art_illumina creates output files in same directory as input files, so 
    # to sandbox results, copy input to output directory
    [[ ! -f "$infile" ]] && cp "$infile_orig" "$infile"
    
    # simulate reads
    art_illumina -i "$infile" -o "${infile%%.fasta}" -l 50 -f "$cov" -q >/dev/null
    
    # Convert art_SOLiD fastq file to SEAStAR style fastq files
    # Also workaround for bug in art_SOLiD 1.0.1 which adds null characters to 
    # some sequences.
    gawk -v output="$outdir/${outpre}_${infile_base%%.fasta}.single.fastq" -f ../../../scripts/fastq_artillumina2seastar.awk "${infile%%.fasta}.fq"
    
    # Append results of this simulation to aggregated file labeled by $outpre
    cat "$outdir/${outpre}_${infile_base%%.fasta}.single.fastq" >>"$outdir/$outpre.single.fastq"
}

# Main ----------------------------------------------------------------------- #
# Check for presence of art_illumina in path
if ! which art_illumina >/dev/null; then
    echo "Error: art_illumina not found in path" >&2
    exit 1
fi

# Output directory for simulated data
outdir=../sim

if [[ -d "$outdir" ]]; then
    echo "$outdir exists, please erase and rerun" >&2
    exit 1
fi
mkdir "$outdir"

# RDP source sequence files used to create simulated read set
files=(S000000075.fasta S000005187.fasta S000498962.fasta S000020432.fasta S000134123.fasta S000431293.fasta S000367885.fasta S000390161.fasta)
# Array of simulated fold coverage values corresponding to elements in "files"
covs=(100 50 40 25 25 15 10 5)

# Simulated reads get placed in fastq subdirectory
# Create sample1 data set
for ((i=0; i < ${#files[@]}; i++ )); do
    simulate "${files[$i]}" "${covs[$i]}" sample1
done
gzip "$outdir"/sample1.single.fastq

# Create sample2 data set with reversed coverage values
for ((i=0; i < ${#files[@]}; i++ )); do
    simulate "${files[$i]}" "${covs[$((${#files[@]} - i - 1))]}" sample2
done
gzip "$outdir"/sample2.single.fastq

# Create a very limited version of RDP sequence database
# This includes the marine uncultured sequence S000498962, but not the version
# with custom lineage
cat S*.fasta >"$outdir"/RDP_slim.fasta

# Clean up intermediate files
# Comment out this command to get a better view of what goes into final files
rm "$outdir"/*.{fq,aln,fastq}
rm "$outdir"/S*.fasta

# Copy the SAR86 sequence file with custom lineage and the SAR86 custom
# hierarchical taxonomy definition lines to the vignette directory.  These
# files will be manually appended to the training set as part of the vignette.
# 
# NB: This step must follow the previous rm, otherwise the fasta file will get
# erased!
cp SAR86_custom_sequence.fasta "$outdir"
cp SAR86_custom_db_taxid.txt "$outdir"
