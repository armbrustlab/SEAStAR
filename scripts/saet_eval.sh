#!/bin/bash -e
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

# Compares csfasta files before and after SAET error correction
# Changes are highlighted with the "|" character
# Stats printed are: read_num, fraction of reads with corrections, corrections per read, corrections per read (only counting those with changes) 
# parms are: original.csfasta corrected.csfasta

# This file is from SEAStAR version: SS_BUILD_VERSION

mkfifo newtmp.fifo
mkfifo oldtmp.fifo
less $1 > oldtmp.fifo &
less $2 > newtmp.fifo &
gawk 'BEGIN { fno = "oldtmp.fifo"; getline o < fno; getline o < fno; } (!(NR % 2) && (o != $0)) { c++; print o; d = ""; for (x = 1; x < length(o); x++) { if (substr($0,x,1) != substr(o,x,1)) { cc++; d = d "|"; } else { d = d " "; } } print d; print $0; print "--", NR/2, c/(NR/2.0), cc/(NR/2.0), cc/c; } (!(NR % 2)) { getline o < fno; getline o < fno  } END { print "----", NR/2, c/(NR/2.0), cc/(NR/2.0), cc/c; }' newtmp.fifo
rm newtmp.fifo
rm oldtmp.fifo
