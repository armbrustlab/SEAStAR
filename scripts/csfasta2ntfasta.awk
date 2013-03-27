# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2012 University of Washington.
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

# Convert an input colorspace fasta file to a naive nucleotide fasta file
# This script is useful for converting Velvet output colorspace contigs
# to nucleotide space for use in building a BWA reference index (where the 
# sequences will be immediately translated back to colorspace). By using 
#     -v numstarts=4 
# on the commandline, this script will produce output sequences that can be 
# used to search databases (e.g. using BLAST).  By default sequence headers
# created by Velvet will by shortened to remove coverage information.  This 
# ID conversion can be turned off by using
#     -v keep_header=1

# Input:   File or stdin
# Output:  stdout

# Parameters:
# numstarts = 1 or 4, controls how many nucleotide sequences to produce
#                     for each input colorspace sequence. 1 will produce a
#                     single sequence using the initial nt base "T", which 
#                     has a ~25% chance of producing correct sequence. 4 
#                     will output four nucleotide sequences, one for each
#                     possible starting nt base [ACGT], with one of these
#                     being the correct choice.

# catalog_file = filename  write a ref_select catalog file to this location

# suffix = string to add to the ends of sequence names

# This file is from SEAStAR version: SS_BUILD_VERSION

BEGIN {

	# Colorspace Decoder Ring

        c["AA"] = "A";  c["AC"] = "C";  c["AG"] = "G";  c["AT"] = "T";
        c["CA"] = "C";  c["CC"] = "A";  c["CG"] = "T";  c["CT"] = "G";
        c["GA"] = "G";  c["GC"] = "T";  c["GG"] = "A";  c["GT"] = "C";
        c["TA"] = "T";  c["TC"] = "G";  c["TG"] = "C";  c["TT"] = "A";

        c["NA"] = "T";  c["NC"] = "G";  c["NG"] = "C";  c["NT"] = "A";
        c["AN"] = "N";  c["CN"] = "N";  c["GN"] = "N";  c["TN"] = "N";
        c["NN"] = "N"; 
	
	div = 2; 	
	seqline = 0;
 
	if (!numstarts) { # Default value for numstarts
		numstarts = 1; 
	} else if ((numstarts != 1) && (numstarts != 4)) {
		# Invalid value for numstarts
		print "Invalid value for numstarts parameter, must be 1 or 4." > "/dev/stderr";
		exit(1);	
	}

	# Default case for numstarts = 1
	nt["T"] = "T";

	if (numstarts == 4) {
		nt["A"] = "A";
		nt["C"] = "C";
		nt["G"] = "G";
	}

        if (!suffix) { suffix = ""; }

        FS = "_";
}

# Capture the fasta header line
/^>/ { 
	write_seq(seq);
        if (NF > 2 && !keep_header) {  # Velvet output
           head = $1 "_" $2 suffix;  
           if (catalog_file) {
              id = substr(head, 2);
              desc = substr($0, 2);
           }
        } else {
	   head = $0 suffix;
        }
	seq = "";
	next; 
}

# Accumulate sequence lines

{
	seq = seq $0;
	next;
}

END { 
	write_seq(seq);
}

function write_seq(seq) {
        if (seq) {
                if (catalog_file) {
                   print id "\t" length(seq) "\t" id "\t" desc > catalog_file;
                }

                seq = toupper(seq);  # Make sure the sequence is all upper case

                # Step through all of the starting bases
                for (nn in nt) {
                        last = nn;
                        if (numstarts == 1) {
                             print head;    # Output fasta header, omit the primer base when only using one start 
                        } else {
                             print head "_" last;  # Output fasta header for this case
                        }
                        out = "";

                        # Loop through the colors in the contig and translate
                        for (x = 1; x <= length(seq); x++) {
                                 last = c[last substr(seq,x,1)];
                                 out = out last;
                        }
                        print out;
                }
        }
}

