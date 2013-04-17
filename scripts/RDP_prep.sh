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
# Prepare RDP database

# This file is from SEAStAR version: SS_BUILD_VERSION

# Functions -------------------------------------------------------------------
function get_ssdir {
    local tmpssdir=$(dirname $(which ref_select 2>/dev/null))
    if [[ -z "$tmpssdir" ]]; then
        echo "Error: Could not find SEAStAR bin directory" >&2
        exit 1
    fi
    echo "$tmpssdir"
}

function usage {
    echo "Usage: RDP_prep [OPTIONS]" >&2
    echo  >&2
    echo "Options: -h         Print this help text" >&2
    echo "         -c         Colorspace data (optional) [FALSE]" >&2
    echo "         -d FILE    RDP sequence fasta file.  Can be gzipped with .gz extension." >&2
    echo "         -s DIR     RDP classifier software directory or zip file" >&2
    echo "         -f FILE    RDP classifier training set fasta file.  Can be gzipped with .gz extension." >&2
    echo "         -t FILE    RDP classifier training set taxonomy file.  Can be gzipped with .gz extension." >&2
    echo "         -o DIR     Output directory" >&2
}


# Configuration ---------------------------------------------------------------
ssdir=$(get_ssdir)  # find SEAStAR directory

# fix_fasta_line_lengths.awk options
minlen=1200
maxlen=1550

# Parse command-line arguments
while getopts :o:d:s:f:t:ch opt; do 
    case "$opt" in
        o)
            outdir=$OPTARG
            ;;
        d)
            rdpseqs=$OPTARG
            ;;
        s)
            class=$OPTARG
            ;;
        f)
            trainfasta=$OPTARG
            ;;
        t)
            traintax=$OPTARG
            ;;
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

# Check required options
if [[ -z "$outdir" ]] || [[ -z "$rdpseqs" ]] || [[ -z "$class" ]] || [[ -z "$trainfasta" ]] || [[ -z "$traintax" ]]; then
    echo "Error: missing required option(s)" >&2
    echo >&2
    usage
    exit 1
fi

# Check for required files and directories
[[ ! -f "$rdpseqs" ]] && echo "Error: file $rdpseqs not found" >&2 && exit 1
[[ ! -e "$class" ]] && echo "Error: $class not found" >&2 && exit 1
[[ ! -f "$trainfasta" ]] && echo "Error: file $trainfasta not found" >&2 && exit 1
[[ ! -f "$traintax" ]] && echo "Error: file $traintax not found" >&2 && exit 1
if [[ -d "$outdir" ]]; then
    echo "Warning: output directory $outdir exists, erasing contents and continuing" >&2
    rm -rf "$outdir"
fi

# Main ------------------------------------------------------------------------
mkdir "$outdir"

# Recreate the rdp classifier software directory in output directory.  It's not 
# very big, and this ensures that RDP_go always works with the prepped
# prepped directory, even if the original classifier software moves.
if [[ -d "$class" ]]; then
    cp -r "$class" "$outdir"
else
    unzip -d "$outdir" "$class" -x '__MACOSX/*'
fi
class="$outdir/$(basename "${class%.zip}")"

echo "Filtering RDP database sequences"
if [[ "${rdpseqs: -3}" == ".gz" ]]; then
    gawk -v minlen="$minlen" -v maxlen="$maxlen" -f "$ssdir/fix_fasta_line_lengths.awk" <(gzip -dc "$rdpseqs") 2>"$outdir"/RDP_rejected.txt | gzip >"$outdir"/RDP.fasta.gz
else
    gawk --re-interval -v minlen="$minlen" -v maxlen="$maxlen" -f "$ssdir/fix_fasta_line_lengths.awk" "$rdpseqs" 2>"$outdir"/RDP_rejected.txt | gzip >"$outdir"/RDP.fasta.gz
fi

echo "Training the RDP classifier"
echo "This step may take a few minutes..."
mkdir "$outdir"/trained
if [[ "${traintax: -3}" == ".gz" && "${trainfasta: -3}" == ".gz" ]]; then
    java -Xmx1g -cp "$class"/rdp_classifier-2.?.jar edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker <(gzip -dc "$traintax") <(gzip -dc "$trainfasta") 9 version custom "$outdir"/trained
elif [[ "${traintax: -3}" == ".gz" && "${trainfasta: -3}" != ".gz" ]]; then
    java -Xmx1g -cp "$class"/rdp_classifier-2.?.jar edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker <(gzip -dc "$traintax") "$trainfasta" 9 version custom "$outdir"/trained
elif [[ "${traintax: -3}" != ".gz" && "${trainfasta: -3}" == ".gz" ]]; then
    java -Xmx1g -cp "$class"/rdp_classifier-2.?.jar edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker "$traintax" <(gzip -dc "$trainfasta") 9 version custom "$outdir"/trained
else
    java -Xmx1g -cp "$class"/rdp_classifier-2.?.jar edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker "$traintax" "$trainfasta" 9 version custom "$outdir"/trained
fi
cp "$class"/sampledata/rRNAClassifier.properties "$outdir"/trained

echo "Running RDP_train_to_tree on training set taxonomies"
if [[ "${traintax: -3}" == ".gz" ]]; then
    gzip -dc "$traintax" | RDP_train_to_tree >"$outdir"/RDP_expand.json
else
    cat "$traintax" | RDP_train_to_tree >"$outdir"/RDP_expand.json
fi

echo "Indexing RDP sequence database for alignment"
if [[ "$color" -eq 1 ]]; then
    RDP_index -c "$outdir"/RDP.fasta.gz
else
    RDP_index "$outdir"/RDP.fasta.gz
fi

