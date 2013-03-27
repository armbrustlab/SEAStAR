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

# Converts a (colorspace) fastq file to csfasta file (without correspnding _QV.qual)
# This is useful for converting SEAStAR style colorspace FASTQ files back to CSFASTA
# format for use with the SAET tool.

# Input:   File or stdin
# Output:  stdout

# This file is from SEAStAR version: SS_BUILD_VERSION

BEGIN {	
	# Hash to lookup CSFASTA numbers from FASTQ letters 
	a["A"] = "0";  
	a["C"] = "1"; 
	a["G"] = "2"; 
	a["T"] = "3"; 
	a["N"] = "."; 
} 

# For every fourth line
((NR % 4) == 1) { 
	print ">" substr($0,5);	# Header line for the read
	# Prepend the first base and color to the out string
	out = substr($0,2,1) a[substr($0,3,1)]; 
	getline; 
	# Loop over each character and convert 
	for (x = 1; x <= length($0); x++) { 
		out = out a[substr($0,x,1)];
	}; 
	print out; 
}

