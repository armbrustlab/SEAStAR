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

# Convert single or mated (read1 & read2) FASTA files to a merged FASTA file suitable
# for paired input to Velvet (read2, read1 (reversed))

# Input:   One or two CSFASTA files (if two, they must contain read1 and read2) 
# Output:  stdout

# Variables R1 and R2 control whether read1 and read2 are reversed, respectively.
# By default, read2 is reversed (R2 = 1) and read1 is not (R1 = 0) 
# These may be overridden on the cmd line e.g. -v R2=0

# This file is from SEAStAR version: SS_BUILD_VERSION

BEGIN { 
	# Lookup table
	a["0"] = "A";  
	a["1"] = "C"; 
	a["2"] = "G"; 
	a["3"] = "T"; 
	a["."] = "N"; 

        if ((ARGV[2]) && ((ARGV[1] !~ /read1/) || (ARGV[2] !~ /read2/))) {
           print "Invalid paired input files. Must be <read1> <read2>" > "/dev/stderr"
           exit(1)
        }

        if (ARGV[2]) {
           read2_fn = ARGV[2]
           ARGV[2] = ""
        }

        if (R1 == "") { R1 = 0; }
        if (R2 == "") { R2 = 1; }
} 

{ 
	print ">" substr($0,2); 
	getline; 
	out = ""; 
	# Skip the primer base and first color
        if (!R1) {
            for (x = 3; x <= length($0); x++) { 
	       out = out a[substr($0,x,1)]; 
	    };
        } else {
           for (x = length($0); x >= 3; x--) {
              out = out a[substr($0,x,1)];
           }
        }  
	print out; 

        if (read2_fn) { 
           getline < read2_fn; 
           print ">" substr($0,2); 
           getline < read2_fn; 
           out = ""; 
           if (!R2) {
               for (x = 3; x <= length($0); x++) {
                  out = out a[substr($0,x,1)];
               };
           } else {
              for (x = length($0); x >= 3; x--) {
                 out = out a[substr($0,x,1)];
              }
           }           
           print out;
        }
}

