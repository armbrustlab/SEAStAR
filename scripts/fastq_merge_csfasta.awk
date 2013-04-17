# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2012-2013 University of Washington.
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

# Merge the sequence of a CSFASTA file with a color-space FASTQ file
# The reads in the files must match exactly, except for the actual sequence
# data. This is useful for merging SAET corrected read data in CSFASTA format
# (which lacks quality data) back into the SEAStAR standard FASTQ format.

# Input:   <reads>.fq <reads>.csfasta 
# Output:  stdout

# This file is from SEAStAR version: SS_BUILD_VERSION

BEGIN { 
	# Lookup table
	a["0"] = "A";  
	a["1"] = "C"; 
	a["2"] = "G"; 
	a["3"] = "T"; 
	a["."] = "N"; 

#        DEBUG = 1;

        if (!ARGV[2]) {
           print "Requires two input files: <reads>.fq <reads>.csfasta" > "/dev/stderr"
           exit(1)
        } 
} 

(ARGIND == 1) {
        getline csfasta_head < ARGV[2];
        getline csfasta_seq < ARGV[2];

        read_name = substr($0,5);
        if (read_name != substr(csfasta_head,2)) {
           print "Read header mismatch!" > "/dev/stderr";
           print substr($0,5) " " length(substr($0,5)) > "/dev/stderr";
           print substr(csfasta_head,2) " " length(substr(csfasta_head,2)) > "/dev/stderr";
           exit(1);
        } 
            
        print "@" substr(csfasta_seq,1,1) a[substr(csfasta_seq,2,1)] "+" read_name;
        getline;
        out = "";
#        diff = "";
        for (x = 3; x <= length(csfasta_seq); x++) { 
	     out = out a[substr(csfasta_seq,x,1)]; 
#             if ((DEBUG) && (x > 3)) {
#                new = a[substr(csfasta_seq,x,1)];
#                old = substr($0,x-2,1);
#                diff = diff ((new == old) ? " " : old);
#             }
	};
#	print out "\n+" diff;
	print out "\n+";
        getline;
        getline;
        print $0;
}

