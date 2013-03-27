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

# Converts an "old style" (colorspace) fastq file to "new style" colorspace 
# fastq file.  For unprocessed reads, it will be faster to simply reconvert
# from CSFASTA and QUAL files to new fastq using the solid2fastq tool.
#
# This script is useful for converting existing fastq files that are the 
# product of more expensive processing (e.g. trimming, deduping, SAET, etc.)  

# Input:   File or stdin
# Output:  stdout

# This file is from SEAStAR version: SS_BUILD_VERSION

{ 
   h = $0; 
   getline s; 
   getline p; 
   getline q; 
   print "@" substr(p,2) "+" substr(h,2); 
   print s; 
   print "+"; 
   print q; 
}
