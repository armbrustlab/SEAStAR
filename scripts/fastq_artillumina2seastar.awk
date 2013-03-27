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
# Convert ART simulated illumina fastq file to SEAStAR style fastq format
# 
# Usage: gawk -v output=<out_filename> -f fastq_artillumina2seastar.awk [<in_filename>]
#
#        or when reading fastq from STDIN
#
#        gawk -v output=<out_filename> -f fastq_artillumina2seastar.awk
#
# Example:
# gawk -v output=test.single.fastq -f fastq_artillumina2seastar.awk test.single.fq

# This file is from SEAStAR version: SS_BUILD_VERSION

# Functions ----------------------------------------------------------------- #

function convert_this_read(read,    read_tmp) {
    # read contains the 4 lines of the ART fastq entry, indexed from 1 to 4
    
    # New ID line
    read_tmp = "@" length(read[2]) "|" gensub(/-/, ":", 1, substr(read[1], 2)) "/1"
    print read_tmp > output
    
    # Sequence line
    print read[2] > output
    
    # Separator line
    print "+" > output
    
    # Quality line
    print read[4] > output
}

# Main ---------------------------------------------------------------------- #
BEGIN {
    if (! output) {
        print "Convert ART simulated single-ended illumina fastq file to SEAStAR style fastq format" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "Usage: gawk -v output=<out_filename> -f fastq_artillumina2seastar.awk <in_filename>" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "       or when reading fastq from STDIN" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "       gawk -v output=<out_filename> -f fastq_artillumina2seastar.awk" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "Example:" > "/dev/stderr"
        print "gawk -v output=test.single.fastq -f fastq_artillumina2seastar.awk test.single.fq" > "/dev/stderr"
        exit
    }
}

# Catch last line of read entry
NR % 4 == 0 {
    last_read[4] = $0
    convert_this_read(last_read)
    next
}

# Catch first three lines of read entry
{
    last_read[NR % 4] = $0
}
