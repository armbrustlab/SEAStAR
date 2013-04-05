###
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
# along with SEAStAR.  If not, see <http:#www.gnu.org/licenses/>.
# -------------------------------------------------------------------------- #
###

#
# This script reads a JSON formatted heirarchy file and converts it to NEXUS
# tree format
#
# Inputs: JSON formatted heirarchy file from stdin.
#
# Output: output NEXUS formatted tree to stdout
#
# Method for branch objects that calculates summary stats for conf, cnt and cum
# And prunes branches from the tree with no representative sequences.
#

norm = 1.0;		# This variable is used to normalize all releative populations to 100% total

walk = (tree) ->
   outstr = ''
   
   if tree.sub?	 # If this node in the tree has children
      children = []
      for c,child of tree.sub	# Recursively generate NEXUS strings for the children
         child.label = c
         children.push(walk(child))
      # NEXUS String for this internal node
      outstr = "(#{String(children)})[&pop=#{Math.round(10000.0 * tree.pop/norm)/10000.0},&tax=#{tree.label},&label=#{tree.label}  #{Math.round(10000.0 * tree.pop/norm)/100.0}%]:#{tree.length}"
   else  # NEXUS string for a leaf (taxa) node
      outstr = "'#{tree.conf}   #{tree.label}   #{Math.round(10000.0 * tree.pop/norm)/100.0}%':#{tree.length}"

inputJSON = ''

# Set up reading from STDIN
process.stdin.setEncoding('utf8');
process.stdin.resume();

# Accumulate the entire JSON input
process.stdin.on 'data', (c) ->
    inputJSON = inputJSON.concat(c)

# Parse the JSON, modify the structre a bit to fit the NEXUS format and then output the NEXUS representation
process.stdin.on 'end', () ->     
   input_tree = JSON.parse(inputJSON)
   input_tree.sub.Root.label = "Root";
   norm = input_tree.sub.Root.pop;
   process.stdout.write("#NEXUS\n\n begin trees;\n tree tree_1 = [&R] #{walk(input_tree.sub.Root)};\nend;\n");
   process.exit 0
