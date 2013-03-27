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
# RDP pipeline short read alignment wrapper
# This implementation uses BWA as the short read aligner of choice.  Any short
# read aligner can be substituted here as long as it records all places where a
# read aligns in the reference and produces a SAM output file.

# This file is from SEAStAR version: SS_BUILD_VERSION

# Functions ----------------------------------------------------------------- #
function usage {
    echo "Usage: RDP_align [OPTIONS] index fastq_prefix output_directory" >&2
    echo 
    echo "Options: -c         Colorspace data [FALSE]" >&2
    echo "         -h         Print this help text" >&2
    echo
    echo "Arguments: index            BWA reference index prefix" >&2
    echo "           fastq_prefix     Fastq file prefix.  Files can be gzipped with .gz extension" >&2
    echo "           output_directory Output directory" >&2
}

function aln_part {
    local fq;
    local part=$1;
    if [[ -f "${fq_prefix}.${part}.fastq" ]]; then
        fq="${fq_prefix}.${part}.fastq"
    elif [[ -f "${fq_prefix}.${part}.fastq.gz" ]]; then
        fq="${fq_prefix}.${part}.fastq.gz"
    fi
    if [[ -n "$fq" ]]; then
        echo " Aligning fastq file $fq"
        if [[ "$color" -eq 1 ]]; then
            bwa aln -n .001 -l 18 -k 2 -c "$aln_index" "$fq" | bwa samse -n 50000 "$aln_index" - "$fq" | gzip >"${out_prefix}.${part}.sam.gz"
        else
            bwa aln -n .001 -l 18 -k 2 "$aln_index" "$fq" | bwa samse -n 50000 "$aln_index" - "$fq" | gzip >"${out_prefix}.${part}.sam.gz"
        fi
    fi
    
}

# Options and Environment --------------------------------------------------- #
if ! which bwa >/dev/null; then
    echo "Error: bwa not found in path" >&2
    echo >&2
    echo usage
    exit 1
fi

while getopts :ch opt; do 
    case "$opt" in
        c)
            color=1
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Error: invalid option -$OPTARG" >&2
            echo >&2
            usage
            exit 1
            ;;
    esac
done

shift $((OPTIND-1))
if [[ $# -lt 3 ]]; then
    echo "Error: Missing argument(s)" >&2
    echo >&2
    usage
    exit 1
fi

aln_index=$1
fq_prefix=$2
outdir=$3
out_prefix="$outdir"/$(basename "$fq_prefix")

# Main ---------------------------------------------------------------------- #
if [[ ! -d "$outdir" ]]; then
    mkdir -p "$outdir"
fi

for part in single single1 single2 read1 read2; do
    aln_part "$part"
done
