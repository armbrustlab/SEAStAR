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
# Run RDP analysis

# This file is from SEAStAR version: SS_BUILD_VERSION

# Functions ----------------------------------------------------------------- #
# Get the real directory path for this script.  If this script is a symbolic
# link, follow the symbolic link to the actual file first.
function real_script_dir {
    link=$(readlink "$0")
    if [[ -z "$link" ]]; then
        # Not a symbolic link
        realpath="$0"
    else
        realpath="$link"
    fi
    thisdir="$(cd "$(dirname "$realpath")" && pwd)"
    echo "$thisdir"
}

function get_vmem_limit {
    local sysname=$(uname)
    local pmem=0
    if [[ "$sysname" == "Linux" ]]; then
        pmem=$(free | awk 'NR == 2 {print $2}')
    elif [[ "$sysname" == "Darwin" ]]; then
        pmem=$(sysctl hw | awk '/hw.memsize:/ {print $2}')
        pmem=$((pmem / 2**10))  # convert from B to KB
    else
        echo "ERROR: System type $sysname not supported" >&2
        exit 1
    fi
    
    local headroom=$((2**10))  # 1 GB of memory head room for other processes
    echo $((pmem - headroom))
}

function get_ssdir {
    local tmpssdir=$(dirname $(which ref_select 2>/dev/null))
    if [[ -z "$tmpssdir" ]]; then
        echo "Error: Could not find SEAStAR bin directory" >&2
        exit 1
    fi
    echo "$tmpssdir"
}

function usage {
    echo "Usage: RDP_go [OPTIONS] <PREFIX>..." >&2
    echo  >&2
    echo "Options: -h         Print this help text" >&2
    echo "         -r DIR     RDP database directory created by RDP_prep" >&2
    echo "         -l INT     Maximum number of reference mappings before a read is rejected from consideration during initial RDP selction step [5000]" >&2
    echo "         -t FLOAT   Minimum bitscore value for a ref sequence to be selected by ref_select during initial RPD selection step [75.0]" >&2
    echo "         -f FLOAT   Minimum bitscore value, as a fraction of the top scoring sequence.  Used by ref_select during initial RDP selection step [0.00125]"  >&2
    echo "         -c         Colorspace data [FALSE]" >&2
    echo "         -d         Debug.  Do not remove intermediate files [FALSE]" >&2
    echo "         -o DIR     Output directory" >&2
    echo >&2
    echo "Arguments: PREFIX   fastq file prefix.  Files may be gzipped with .gz extension." >&2
}

function in2out_prefix {
    local prefix=$(basename "$1")
    echo "$outdir"/${prefix%%/*}
}

function align_prefix {
    local fqprefix=$1
    local prefix=$(in2out_prefix "$fqprefix")
    
    echo "Aligning fastq prefix $fqprefix against full RDP database"
    if [[ "$color" -eq 1 ]]; then
        RDP_align -c "$rdpdir/RDP.fasta.gz" "$fqprefix" "$prefix"_all_taxa_aln
    else
        RDP_align "$rdpdir/RDP.fasta.gz" "$fqprefix" "$prefix"_all_taxa_aln
    fi
    
    echo "Running ref_select"
    ref_select --ref "$rdpdir/RDP.fasta.gz" -a -l "$read_map_limit" -s -t "$bit_thresh" -f "$bit_fraction" "$prefix"_all_taxa_aln/*.sam.gz | gzip >"$prefix"_all_taxa_seastar.json.gz
    
    echo "Creating new RDP database with only selected sequences"
    graph_ops "$prefix"_all_taxa_seastar.json.gz FASTA >"$prefix"_selected.fasta
    
    echo "Creating BWA index for selected RDP sequences"
    if [[ "$color" -eq 1 ]]; then
        bwa index -a is -c "$prefix"_selected.fasta
    else
        bwa index -a is "$prefix"_selected.fasta
    fi
    
    echo "Aligning against selected RDP sequence database"
    if [[ "$color" -eq 1 ]]; then
        RDP_align -c "$prefix"_selected.fasta "$fqprefix" "$prefix"_selected_taxa_aln
    else
        RDP_align "$prefix"_selected.fasta "$fqprefix" "$prefix"_selected_taxa_aln
    fi
    
    echo "Running ref_select"
    ref_select --ref "$prefix"_selected.fasta -t 0.20 -f 0.0085 -s --per_base "$prefix"_selected_taxa_aln/*.sam.gz | gzip >"$prefix"_final_taxa_seastar.json.gz
    
    echo "Extract final selected sequence set"
    graph_ops "$prefix"_final_taxa_seastar.json.gz FASTA '{"abundance":true}' >"$prefix"_final_taxa.fasta
    java -Xmx1g -jar "$rdpdir"/rdp_classifier_2.?/rdp_classifier-2.?.jar -q "$prefix"_final_taxa.fasta -o /dev/stdout -f allrank -t "$rdpdir"/trained/rRNAClassifier.properties >"$prefix"_final_class.txt
    RDP_tree_dev "$rdpdir"/RDP_expand.json 0.1 <"$prefix"_final_class.txt >"$prefix"_final_class.json
    
    # Cleanup some intermediate files
    cleanup "$prefix"_all_taxa_aln
    cleanup "$prefix"_all_taxa_seastar.json.gz
    cleanup "$prefix"_selected.fasta
    cleanup "$prefix"_selected.fasta*
}

function cleanup {
    if [[ "$debug" -eq 0 ]]; then
        if [[ -d "$1" ]]; then
            rm "$1"/*.sam.gz
            rmdir "$1"
        else
            rm "$@"
        fi
    fi
}

# Configuration and Environment --------------------------------------------- #
script_dir=$(real_script_dir)

# This can be set to any alignment wrapper which takes index as arg1, 
# fastq input prefix as arg2, output directory with sam.gz files as arg3.
aligner="RDP_align"

# Virtual memory limit for ulimit -v, in kilobytes.  Protects against
# unintentionally bringing down the server
vmem_limit=$(get_vmem_limit)
ulimit -v "$vmem_limit"

ssdir=$(get_ssdir)

debug=0

# Parse command-line arguments
while getopts :r:l:t:f:o:cdh opt; do 
    case "$opt" in
        r)
            rdpdir=$OPTARG
            ;;
        l)
            read_map_limit=$OPTARG
            ;;
        t)
            bit_thresh=$OPTARG
            ;;
        f)
            bit_fraction=$OPTARG
            ;;
        o)
            outdir=$OPTARG
            ;;
        c)
            color=1
            ;;
        d)
            debug=1
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
    echo "Error: Missing PREFIX argument" >&2
    echo >&2
    usage
    exit 1
fi
prefixes=("$@")

# Simple check to make sure these are prefixes and not full fastq file names
for fqp in "${prefixes[@]}"; do
    if [[ "$fqp" =~ ^.*\.fastq$|^.*\.gz$ ]]; then
        echo "Error: One or more prefix does not appears to be a file name.  Perhaps a file name was entered in error?" >&2
        echo >&2
        usage
        exit 1
    fi
done

if [[ -z "$rdpdir" ]] || [[ -z "$outdir" ]]; then
    echo "Error: missing required option(s)" >&2
    echo >&2
    usage
    exit 1
fi

if [[ -z "$read_map_limit" ]]; then
    read_map_limit="5000"
fi
if [[ -z "$bit_thresh" ]]; then
    bit_thresh="75.0"
fi
if [[ -z "$bit_fraction" ]]; then
    bit_fraction="0.00125"
fi

[[ ! -d "$rdpdir" ]] && echo "Error: directory $rdpdir not found" >&2 && exit 1
if [[ -d "$outdir" ]]; then
    echo "Warning: output directory $outdir exists, erasing contents and continuing" >&2
    rm -rf "$outdir"
fi

# Main ---------------------------------------------------------------------- #
mkdir "$outdir"
final_jsons=()
for fqprefix in "${prefixes[@]}"; do
    align_prefix "$fqprefix"
    final_jsons+=($(basename $(in2out_prefix "$fqprefix")_final_class.json))
done

cd "$outdir"
# Create directory for javascript visualization
mkdir html
cp "$script_dir"/RDP_viz/raphael*.js "$script_dir"/RDP_viz/tree_bars.html "$script_dir"/RDP_viz/tree_bars.js html
RDP_tree_merge_dev "${final_jsons[@]}" | gawk 'BEGIN {print "var tree_struct = "} {print} END {print ";"}' >html/RDP_merged.js
