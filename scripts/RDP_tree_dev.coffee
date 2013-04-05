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
#
# This script reads in an RDP classifier taxonomic heirarchy file and converts
# it to a JSON tree formatted output file.
#
# Inputs: input modified RDP classifier output on stdin 
#         JSON formatted RDP heirarchy in the file named in RDP_heir_fn (or RDP_expand.json by default)
#         percent abundance below which a genus will be filtered out [0.10]
#
# Output:  output JSON formatted heirarchy to stdout

# Output file format:
# Each taxonomic unit in the RDP input file becomes a branch in the output JSON tree.
# Each branch has the following JSON structure:
#
# name : {
#      "name" : <name of this taxon>,
#      "pop" : <fraction of the sample>,
#      "cum" : <cumulative fraction of the sample>,
#      "num" : <number of sequences>,
#      "cnt" : <cumulative number of seqeunces>,
#      "conf" : <mean classifier p-value for sequences in this taxon>,
#      "level" : <numeric level in heirarchy>,
#      "length" : <branch length of this taxonomic level from parent>
#      "sub" : {<child nodes by name, or empty if a leaf>}
# }
###

fs = require('fs')

# Handle default filename and genus minimum percent abundance
RDP_heir_fn =  process.argv[2] ? 'RDP_expand.json'
minperc = parseFloat(process.argv[3] ? 0.10)

# Read and parse the JSON RDP taxonomy
out_tree = JSON.parse(fs.readFileSync(RDP_heir_fn,'utf8'))

# Default exception handler
process.on 'uncaughtException', (err) -> 
   console.log 'Caught exception: ' + err

# Get stdin ready for reading
process.stdin.setEncoding('utf8')
process.stdin.resume()

# levels present in RDP taxonomic hierarchy
levels = [0, 1, 2, 3, 3.5, 4, 4.5, 5, 5.5, 5.75, 6]

# lookup table of next levels
level_after = {}
i = 0
while i < levels.length - 1
   level_after[levels[i]] = levels[i+1]
   i++

# level name to number lookup
level_name2num = 
   norank     : 0       # trainset 6 hierarchy convention
   rootrank   : 0       # trainset 9+ hierarchy convention
   domain     : 1
   phylum     : 2
   class      : 3 
   subclass   : 3.5
   order      : 4
   suborder   : 4.5
   family     : 5
   subfamily  : 5.5
   supergenus : 5.75 
   genus      : 6

# level number to name lookup
for level_name, level_num in level_name2num
   level_num2name[level_num] = level_name

# Method for branch objects that calculates summary stats for conf, cnt and cum
# And prunes branches from the tree with no representative sequences.
out_tree.walk = (tree, prev_cnt, prev_cum) ->
   if tree.sub? and Object.keys(tree.sub).length > 0   # Make sure there is something to do
      # Shuffle sort child elements by pop (largest, smallest, second largest, second smallest...)
      
      # Create an array and reverse sort it by population.
      tmp = (z for y,z of tree.sub) 
      tmp.sort((b,a) -> a.pop-b.pop)
      
      # Now replace the original sub object
      tree.sub = {}
      # This performs the "shuffle" on the sorted array
      
      while tmp.length
         c = tmp.shift()   # grab first element
         tree.sub[c.name] = c
         c = tmp.pop()     # then last
         tree.sub[c.name] = c if c
         
   for c of tree.sub 
      if tree.sub[c].num != 0 
         sub_tree = tree.sub[c]
         sub_tree.conf = sub_tree.conf / sub_tree.num
         sub_tree.w_conf = sub_tree.w_conf / sub_tree.pop
         sub_tree.cnt = prev_cnt
         sub_tree.cum = prev_cum
         out_tree.walk(sub_tree, prev_cnt, prev_cum)	
         prev_cnt += sub_tree.num
         prev_cum += sub_tree.pop
      else
         delete tree.sub[c]

# Walk all lines in the RDP classifier output, return a list of lines for
# genera which comprise >= minfrac fraction of the 16S abundance
filter = (line_list, minperc) ->
   genera = {}
   for line in line_list
      fields = parse_line line
      genus = fields.tax_names.pop()
      if genus of genera
         genera[genus].percent += fields.percent
         genera[genus].lines.push line
      else
         genera[genus] =
            percent : fields.percent
            lines : [line]
   retlines = []
   for genus, data of genera
      if data.percent >= minperc
         Array::push.apply retlines, data.lines
   retlines

parse_line = (line) ->
   fields = {}
   parts = line.split('\t')                         # Split on semicolons
   [fields.sequence, fields.percent] = parts.shift().split('_')   # Split first field on "_" (remove 1st & 2nd)
   fields.percent = parseFloat(fields.percent)                    # This is its estimated abundance
   
   # Parsed taxonomic names
   parts.shift()
   fields.tax_names = (p.replace(/"/g,'').trim() for p in parts by 3)
   # Parsed taxonomic levels
   parts.shift()
   fields.tax_levels = (p.replace(/"/g,'').trim() for p in parts by 3)
   # Parsed p-value data
   parts.shift()
   fields.tax_pvals = (parseFloat(p) for p in parts by 3)
   fields

# Walk all of the lines in the RDP classifier output 
build = (line) ->
   fields = parse_line line
   
   # Walk down the branches of the RDP taxonomy tree and fill in the 
   # information for this sequence.
   prev = null
   cur = out_tree
   for tax_name, i in fields.tax_names
      # Add info for this sequence to any incertae sedis entries this script
      # already created for previous lines
      prev = cur
      start = prev
      while not (tax_name of cur.sub)
         cur = cur.sub["#{start.name}_incertae_sedis"]
         prev = cur
         cur.num++
         cur.conf += fields.tax_pvals[i]
         cur.w_conf += fields.tax_pvals[i] * fields.percent
         cur.pop += fields.percent
      
      # cur.sub should contain tax_name now
      cur = cur.sub[tax_name]
      
      # Insert incertae sedis entries into lineage for missing levels
      fill_lineage_gap(prev, cur)
      
      # Update this taxon
      cur.num++
      cur.conf += fields.tax_pvals[i]
      cur.w_conf += fields.tax_pvals[i] * fields.percent
      cur.pop += fields.percent
   
   # Now add the sequence itself as a leaf of the lowest taxon level (usually genus, level 6)		
   cur.sub[fields.sequence] =
      pop : fields.percent
      cum : 0.0
      cnt : 0
      num : 1
      conf : fields.tax_pvals[fields.tax_pvals.length-1] 
      w_conf : fields.percent * fields.tax_pvals[fields.tax_pvals.length-1]
      level : 7.0
      length : 1.0
      name : fields.sequence

# Fill in missing taxonomic levels between begin and end with incertae sedis
# copies of begin
fill_lineage_gap = (begin, end) ->
   if not begin?.level? or not end?.level?   # begin or end is not valid taxon
      return
   if level_after[begin.level] == end.level  # no gap to fill
      return
   level = level_after[begin.level]
   cur = begin
   # Now insert taxa to fill gaps between levels
   while level != end.level
      insert_name = "#{begin.name}_incertae_sedis"
      save = cur.sub  # save original sub
      cur.sub = {}
      cur.sub[insert_name] =
            name : insert_name
            pop    : begin.pop
            cum    : begin.cum
            cnt    : begin.cnt
            num    : begin.num
            conf   : begin.conf
            w_conf : begin.w_conf
            level  : level
            length : level - cur.level
      level = level_after[level]
      cur = cur.sub[insert_name]
      cur.sub = save  # put original sub back
   end.length = end.level - cur.level  # recalculate end's length

orig_lines = []
process.stdin.on 'data', do ->
   save = ''
   return (c) ->
      lines = c.split '\n'
      lines[0] = save + lines[0]
      save = lines.pop()
      Array::push.apply orig_lines, lines

process.stdin.on 'end', () ->     # Write the JSON structure to stdout.
   filtered_lines = filter orig_lines, minperc # filter genera by minimum fractional abundance
   build i for i in filtered_lines when i
   # Walk the tree, pruning empty branches and updating cumulative statistics
   out_tree.walk out_tree, 0, 0.0
   # Output JSON to stdout
   process.stdout.write JSON.stringify(out_tree, null, 1)+"\n"
   process.exit 0
