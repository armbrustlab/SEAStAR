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

# Reads a FASTA file and filters the sequences therein by rejecting those that
# do not fall within a range of lengths [minlen,maxlen], or which contain 
# significant homopolymer runs.  The output FASTA file line length can be  
# selected, allowing this tool (with the filters off) to be used to add or 
# remove sequence linebreaks in a FASTA file.

# Input:   File or stdin - input FASTA file
# Output:  stdout - output FASTA file
#          stderr - Diagnostics about rejected sequences
# Parameters:
# 
# ll = output sequence line length. default = infinte (actually 4 billion)
# maxlen = reject sequences longer than this. default = infinite (actually 4 billion)
# minlen = reject sequences shorter than this. default = 1 
# filter = reject sequences with significant homopolymer runs. Default = 0 (off) 
#          Where significant is a run of nucleotides that a ~25 base colorspace read 
#          of all A's (or more properly all 0s) would align to with 4 or fewer mismatches. 

# This function takes a seqeunce header z, and seqeunce x and processes them
# by filtering on length, optionally rejecting seqeunces with homopolymer runs
# and then writing out sequences with the requested line lengths

# This file is from SEAStAR version: SS_BUILD_VERSION

function process_sequence(s, h,    y) {

        if (h && (length(s) <= maxlen) && (length(s) >= minlen)) {
                if (filter && (s ~ /((ttttt|aaaaa|ccccc|ggggg|nnnnn)[acgtn]{0,2}){4}/)) {
                        printf("Homopoly Rejected: %s %s\n", h, s) > "/dev/stderr";
                } else {
                        print h;
                        for (y = 1; y <= length(s); y += ll) {
                                print substr(s, y, ll);
                        }
                }
        } else if (h) {
                printf("Length Rejected (%d): %s\n", length(s), h) > "/dev/stderr";
        }

	return;
}

BEGIN { 
        # Make damned sure that interval expressions are turned on!
        if (!("aaaa" ~ /a{3,5}/)) {
                print "You MUST use the GAWK --re-interval command-line option for this script to work correctly!!";
                print "You MUST use the GAWK --re-interval command-line option for this script to work correctly!!" > "/dev/stderr";
                exit(1);
        }

        if (ll == 0) { ll = 4000000000; }
        if (maxlen == 0) { maxlen = 4000000000; }
        if (!minlen) { minlen = 1; }
        FS = "_";
} 

# Process new sequence
(substr($1,1,1) == ">") { 

	process_sequence(sequence, head);
	
	# Reset variables for next seqeunce
	sequence= "";
	head = $0;
	next;
}

# Process seqeunce line
{ 
	sequence = sequence $0;  # Append to previous seqeunce lines
}

# Process last sequence
END {

	process_sequence(sequence, head);

}

 
