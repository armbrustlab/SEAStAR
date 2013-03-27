###
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
# along with SEAStAR.  If not, see <http:#www.gnu.org/licenses/>.
# -------------------------------------------------------------------------- #
#
# This script reads in one or more JSON formatted heirarchy files and merges
# them into to a new unified JSON tree formatted output file.
#
# Inputs: JSON formatted heirarchy files in the global variable input_files.
#
# Output: output JSON formatted heirarchy to stdout
#
# Output file format:
# Each taxonomic unit in the RDP input file becomes a branch in the output JSON tree.
# Each branch has the following JSON structure:
#
# name : {
#      "pop" : <fraction of the sample>,
#      "cum" : <cumulative fraction of the sample>,
#      "num" : <number of sequences>,
#      "cnt" : <cumulative number of seqeunces>,
#      "conf" : <mean classifier p-value for sequences in this taxon>,
#      "w_conf" : <population weighted mean classifier p-value for sequences in this taxon>,
#      "level" : <numeric level in heirarchy>,
#      "length" : <branch length of this taxonomic level from parent>
#      "samples" : <sample child nodes containing pop, cum, num, cnt and conf> 
#      "sub" : {<child nodes by name, or empty if a leaf>}
# }
#
###

fs = require('fs')
path = require('path')

# Handle default filename.  

input_files = {}

# Build an object of sample names and corresponding filenames
# Sample name is the string before the first "_" character
for i in process.argv[2..]
   input_files[path.basename(i).split("_")[0]] = i

# Quit without input.  
if Object.keys(input_files).length is 0
   process.stderr.write "No input files found!\n"
   process.exit 1

# Initialize the output object
output_tree = { sample_names : [], sub : {} }

# Walk branch objects and calculates summary stats for cnt and cum
walk = (tree, prev_cnt, prev_cum) ->
   for c of tree.sub
      sub_tree = tree.sub[c]
      sub_tree.cnt = prev_cnt
      sub_tree.cum = prev_cum
      walk(sub_tree, prev_cnt, prev_cum)
      prev_cnt += sub_tree.num
      prev_cum += sub_tree.pop

# Merge a sample into the existing output object
merge = (out_tree, in_tree, sample, prev_cnt, prev_cum) ->
   
   # Walk through all of the subtaxa at this level
   for c of in_tree.sub
      # If the output tree doesn't already have this taxon level
      unless out_tree.sub[c]? 
         # Add it!
         out_tree.sub[c] = in_tree.sub[c]
      
      out_tree.sub[c].samples ?= []
      
      # Make a new object for this sample
      out_tree.sub[c].samples[sample] =   
         num : in_tree.sub[c].num 
         cnt : in_tree.sub[c].cnt
         pop : in_tree.sub[c].pop
         cum : in_tree.sub[c].cum 
         conf : in_tree.sub[c].conf
         w_conf : in_tree.sub[c].w_conf
      
      # Recurse
      merge(out_tree.sub[c], in_tree.sub[c], sample, prev_cnt, prev_cum)
      
      # initialize summary stats properties
      out_tree.sub[c].num = 0
      out_tree.sub[c].cnt = 0
      out_tree.sub[c].pop = 0.0
      out_tree.sub[c].cum = 0.0
      out_tree.sub[c].conf = 0.0
      out_tree.sub[c].w_conf = 0.0		
      
      # Roll-up summary statistics from samples
      for k,s of out_tree.sub[c].samples
         out_tree.sub[c].num += s.num
         out_tree.sub[c].pop += s.pop
         out_tree.sub[c].conf += s.num * s.conf
         out_tree.sub[c].w_conf += s.pop * s.w_conf		
      
      # Round the conf values
      out_tree.sub[c].conf = Math.round(
         100.0 * out_tree.sub[c].conf / out_tree.sub[c].num) / 100.0 
      out_tree.sub[c].w_conf = Math.round(
         100.0 * out_tree.sub[c].w_conf / out_tree.sub[c].pop) / 100.0 

for sample_name,sample_file of input_files
   # Read and parse the JSON RDP taxonomy
   merge(output_tree, JSON.parse(fs.readFileSync(sample_file,'utf8')), output_tree.sample_names.length, 0, 0.0)
   output_tree.sample_names.push(sample_name)

walk(output_tree, 0, 0.0)

# Output JSON to stdout
process.stdout.write JSON.stringify(output_tree, null, 1)

process.exit 0
