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
# Convert ART simulated SOLiD fastq file to SEAStAR style fastq format
# 
# Usage: gawk -v output=<out_filename> [-v prefix=<prefix>] -f fastq_artsolid2seastar.awk [<in_filename>]
#
#        or when reading fastq from STDIN
#
#        gawk -v output=<out_filename> [-v prefix=<prefix>] -f fastq_artsolid2seastar.awk
#
# prefix can be any identifying string that should be placed before the bead
# number of every read, separated by ':'.
# 
# Example:
# gawk -v output=test.single.fastq -v prefix="sample1" -f fastq_artsolid2seastar.awk test.single.fq
#

# This file is from SEAStAR version: SS_BUILD_VERSION

# Functions ----------------------------------------------------------------- #
# From http://www.gnu.org/software/gawk/manual/html_node/Ordinal-Functions.html
function _ord_init(    low, high, i, t)
{
    low = sprintf("%c", 7) # BEL is ascii 7
    if (low == "\a") {    # regular ascii
        low = 0
        high = 127
    } else if (sprintf("%c", 128 + 7) == "\a") {
        # ascii, mark parity
        low = 128
        high = 255
    } else {        # ebcdic(!)
        low = 0
        high = 255
    }
    
    for (i = low; i <= high; i++) {
        t = sprintf("%c", i)
        _ord_[t] = i
   }
}

function ord(str,    c)
{
    # only first character is of interest
    c = substr(str, 1, 1)
    return _ord_[c]
}

function contains_null_char(s) {
    for (i=1; i <= length(s); i++) {
        if (_ord_[substr(s, i, 1)] == 0) {
            return 1
        }
    }
    return 0
}

function convert_this_read(read,    read_tmp) {
    # Read contains the 4 lines of the ART fastq entry, indexed from 1 to 4
    
    # Skip read where sequence contains null characters
    # This is a bug encountered in ART SOLiD read simulator v1.0.1
    if (contains_null_char(read[2])) {
        #print "INFO: Skipping read " read[1] ", sequence contains NULL character" > "/dev/stderr"
        skipped++
        return
    }
    
    # New ID line
    read_tmp = "@"
    read_tmp = read_tmp substr(read[2], 1, 1) num2base[substr(read[2], 2, 1)] # last adaptor and first color call    
    read_tmp = read_tmp "+" (length(read[2]) - 2) # read length
    if (prefix) {
        read_tmp = read_tmp "|" prefix ":" substr(read[1], 2, length(read[1]) - 4) # bead ID
    } else {
        read_tmp = read_tmp "|" substr(read[1], 2, length(read[1]) - 4) # bead ID with tag
    }
    read_tmp = read_tmp mp[substr(read[1], length(read[1]) - 2)] # read position in mate pair
    print read_tmp > output
    
    # New sequence line
    read_tmp = ""
    for (i=3; i <= length(read[2]); i++) {
        read_tmp = read_tmp num2base[substr(read[2], i, 1)]
    }
    print read_tmp > output
    
    # Separator line
    print "+" > output
    
    # Quality line
    print substr(read[4], 2) > output
}

# Main ---------------------------------------------------------------------- #
BEGIN {
    if (! output) {
        print "Convert ART simulated SOLiD fastq file to SEAStAR style fastq format" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "Usage: gawk -v output=<out_filename> [-v prefix=<prefix>] -f fastq_artsolid2seastar.awk <in_filename>" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "       or when reading fastq from STDIN" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "       gawk -v output=<out_filename> [-v prefix=<prefix>] -f fastq_artsolid2seastar.awk" > "/dev/stderr"
        print "" > "/dev/stderr"
        print "prefix can be any identifying string that should be placed before the bead" > "/dev/stderr"
        print "number of every read, separated by ':'." > "/dev/stderr"
        print "Example:" > "/dev/stderr"
        print "gawk -v output=test.single.fastq -v prefix='sample1' -f fastq_artsolid2seastar.awk test.single.fq" > "/dev/stderr"
        exit
    }
    
    num2base[0] = "A"
    num2base[1] = "C"
    num2base[2] = "G"
    num2base[3] = "T"
    mp["_F3"] = "/1"
    mp["_R3"] = "/2"
    
    _ord_init()
    
    skipped = 0
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

END {
    print "INFO: Skipped " skipped " reads in total" > "/dev/stderr"
}
