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

# This script takes the stdout of ref_select and uses it to select sequences 
# from the reference fasta file used to build the alignment index that ref_select 
# processed. It will optionally include the estimated relative abundance on the 
# header lines of the output fasta file (so that they can be included in the 
# output of the RDP classifier, for example).

# Input file #1:  text stdout output from ref_select
# Input file #2:  fasta file from alignment reference (sequences in ref_select catalog)
# Output:  selected fasta sequences to stdout    

# This file is from SEAStAR version: SS_BUILD_VERSION

BEGIN {
	pop_sum = 0.0;
} 

((ARGIND == 1) && ($5 > 0.0)) {
	b[">" $1] = $5;
	pop_sum += $5;
	c++;
	next
}

(ARGIND != 1) { c2++ }

((FNR % 2 == 1) && ($1 in b)) { 
	if (abundance) {
		print $1 "_" b[$1]*100; 
	} else {
		print $1; 
	}
	getline; 
	print $0
}

END {
	printf("Cluster population: %.4f%% Num taxa: %d %d\n", 100*pop_sum, c, c2) > "/dev/stderr";
}
