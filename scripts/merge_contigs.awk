#!/usr/bin/gawk -f
# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2012-2013 University of Washington.
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


# This script takes a multifasta file with one sequence per contig, 
# where the contigs are in the correct order and orientation for one
# or more scaffolds. Scaffolds are indicated by a white space separated
# sequence number on the FASTA header line:
# >contig_1234 1
# >contig_543 1
# >contig_984 2
# >contig_12 2
# The above FASTA headers describe two scaffolds of two contigs each.
#
# The script outputs a multifasta file with one sequence per scaffold.
# Where the contigs are joined together either by gaps "nnnn..." or
# by overlapping contig end sequence detected in the code below.

# This file is from SEAStAR version: SS_BUILD_VERSION

BEGIN {
	# Setup default values
	if (!n) { n = 6; }	# n = required bases of overlap
	if (!m) { m = 0; }	# m = number of non-ambiguous bases allowed to be trimmed off each contig end
	if (!max) { max = 35; }	# max = maximum overlap to look for
	print "m:", m, "n:", n, "max:", max, "\n" > "/dev/stderr";
} 

# Capture the contig FASTA header
/^>/ {
	prev_head = head; 
	head = $0;

# Detect a new scaffold!	
	if (scaffold && ($2 != scaffold)) {
		# Write the previous scaffold.
		print ">Scaffold_" scaffold;
		print sequence;
		sequence = "";
	}

	# Remember the current scaffold id
	scaffold = $2;
 
	next;
} 

# Process contig sequence
{
	# inc the contig counter
	cont++; 
	# Trim ambiguous bases off the ends of this contig
	match($0, /^[^ACGT]*(.*[ACGT])[^ACGT]*$/, a); 
	contig = a[1]; 
	
	# If a contig has already been added to this scaffold
	if (sequence) {

		l = length(sequence);
		
		# These nested for loops perform the overlap detection
		# Loop for overlaps from max down to n bases.
		for (i = max; i >= n; i--) {
			c = 0;	# Holds the current number of aligning bases
			# Loop through the bases of a potential i base overlap	
			for (j = 1; j <= i; j++) {
				# If this position aligns
				if (substr(sequence, l - i + j, 1) == substr(contig, j, 1)) {
					# increment the alignment counter
					c++;
				} else if ((c < n) || (c < i-(2*m))) {
					# If this position doesn't align, and the previous matches are insufficient
					# Reset the counter...
					c = 0;
				} else {
					# We've just found a match that ends here.
					break;	
				}
			}
			# If a match is found in the i base alignment
			if ((c >= n) && (c >= i-(2*m))){
				# Write diagnostics to stderr	
				print "\n**** Found match! ", cont, i, j, c > "/dev/stderr";
				print substr(sequence, l-max+1) > "/dev/stderr";
				# Determine the overlap bases, set to lowercase
				overlap = tolower(substr(sequence, l-i+(j-c), c));
				# Merge the previous scaffold and current contig around the aligned overlap
				sequence = substr(sequence, 1, l-i+(j-c)-1) overlap substr(contig, c+(j-c));
				# Print the alignment to stderr
				for (x = 1; x < max - i + (j-c); x++) {
					printf(" ") > "/dev/stderr";
				}
				print overlap > "/dev/stderr";
				for (x = 1; x <= max - i; x++) {
					printf(" ") > "/dev/stderr";
				}
				print substr(contig, 1, max) > "/dev/stderr";
				printf("\n") > "/dev/stderr";
				break;
			} else {	# This contig overlap alignment (i bases) did not have a match.
				c = 0;
			}
		}
		
		# If there is no overlapping alignment of these two contigs, just add a gap and the new contig.
		if (c == 0) {
			print "---- no match: ", cont > "/dev/stderr";
			sequence = sequence "nnnnnnnnnnnnnnn" contig;
		}
	} else {  # Handle the first contig per scaffold
		sequence = contig;
	} 
	next; 
}

# Output the final scaffold
END {
	if (scaffold) {
		print ">Scaffold_" scaffold;
		print sequence;
	}
}

