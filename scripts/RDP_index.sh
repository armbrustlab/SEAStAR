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
# RDP pipeline reference indexing wrapper
# This implementation creates a BWA index for the reference fasta file

# This file is from SEAStAR version: SS_BUILD_VERSION

# Functions ----------------------------------------------------------------- #
function usage {
    echo "Usage: RDP_index [OPTIONS] reference_fasta" >&2
    echo 
    echo "Options: -h         Print this help text" >&2
    echo "         -c         Colorspace data [FALSE]" >&2
    echo
    echo "Arguments: reference_fasta  Reference fasta file to index" >&2
}

# Options and Environment --------------------------------------------------- #
if ! which bwa >/dev/null 2>&1; then
    echo "Error: bwa not found in path" >&2
    echo >&2
    usage
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
if [[ $# -lt 1 ]]; then
    echo "Error: Missing argument" >&2
    echo >&2
    usage
    exit 1
fi

ref_fasta=$1

# Main ---------------------------------------------------------------------- #
if [[ "$color" -eq 1 ]]; then
    bwa index -c -a is "$ref_fasta"
else
    bwa index -a is "$ref_fasta"
fi
