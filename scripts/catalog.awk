#!/usr/bin/gawk -f
# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2013 University of Washington.
#
# Authors:
# Chris Berthiaume
# chrisbee@uw.edu
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
# Create a SEAStAR reference sequence catalog file
# Can handle multi-line fasta files with spaces in sequence and empty lines.
#
# usage: catalog.awk ref.fasta >ref.catalog.txt

# This file is from SEAStAR version: SS_BUILD_VERSION

function handle_last_seq(seq, seqid,   fields,   n,   i,   desc) {
    if (! length(seq)) {
        return
    }
    
    gsub(/[ \t]/, "", seq)  # Remove whitespace from sequence
    
    # If no fasta header description, use NA as placeholder string
    n = split(seqid, fields)
    if (n > 1) {
        desc = ""
        i = 2
        while (i < n) {
            desc = desc fields[i] " "
            i++
        }
        desc = desc fields[n]
        printf("%s\t%d\t%s\t%s\n", fields[1], length(seq), fields[1], desc)
    } else {
        printf("%s\t%d\t%s\t%s\n", seqid, length(seq), seqid, "NA")
    }
}

/^[^>]/ || /^ *$/ {
    # Sequence line, blank, or only spaces
    seq = seq $0
    next
}

{
    # ID line
    handle_last_seq(seq, seqid)
    # Reset and get ready for next sequence
    seqid = substr($0, 2)
    seq = ""
}

END {
    handle_last_seq(seq, seqid)
}
