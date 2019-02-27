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

# This script reads in an RDP classifier taxonomic heirarchy file and converts
# it to a JSON tree formatted output file.

# Input: input RDP heirarchy on stdin
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
#      "w_conf" : <population weighted mean classifier p-value for sequences in this taxon>,
#      "level" : <numeric level in heirarchy>,
#      "length" : <branch length of this taxonomic level from parent>
#      "sub" : {<child nodes by name, or empty if a leaf>}
# }

ver = process.version[1..].split('.')
unless ver[1] >= 10 or ver[0] > 0   # Require node version >= 0.10.x
   console.error("ERROR: nodejs version v0.10.0 or greater required.")
   process.exit(1)

taxa = [] # This is a working array of taxa records indexed by taxid

levels =
  rootrank   : 0
  norank     : 0
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

error_codes =
  INVALID_INPUT_LINE : 1

process.on 'uncaughtException', (err) ->
  console.log 'Caught exception: ' + err

process.stdin.resume()
process.stdin.setEncoding('utf8')

# Now loop through lines, building the taxa table and making the heirarchy

build = (line) ->
  return unless (line = line.trim())	# Remove whitespace from ends
  # Split on asterisks
  # process.exit error_codes.INVALID_INPUT_LINE if ([taxid, name, parentid, levnum, level] = line.split("*")).length isnt 5
  return if ([taxid, name, parentid, levnum, level] = line.split("*")).length isnt 5
  taxid = parseInt(taxid)	# taxid and parentid
  parentid = parseInt(parentid)	# are each integers
  name = name.replace(/"/g,"")	# remove quotes from taxa names
  level = level.trim();		# be tolerant of whitespace

  return unless level of levels

  # If this taxon has already been used
  # as a parent (it is out of order)
  # then copy the existing sub taxa

  if taxa[taxid]?
    taxa[taxid].level = levels[level]
    taxa[taxid].name = name
    # Go back and fill in the branch length value for child taxa
    for i, child of taxa[taxid].sub
    # for child in taxa[taxid].children
      child.length = child.level - taxa[taxid].level
  else  # otherwise make a whole new object
    taxa[taxid] =
      name : name
      pop : 0.0
      cum : 0.0
      cnt : 0
      num : 0
      conf : 0.0
      w_conf : 0.0
      level : levels[level]
      length : 0.0
      # children : []
      sub : {}

  # Look to see if this taxon's parent already has a record
  if taxa[parentid]?
    taxa[parentid].sub[name] = taxa[taxid]	# Add this child
    # taxa[parentid].children.push(taxa[taxid])	# Add this child
    taxa[taxid].length = taxa[taxid].level - taxa[parentid].level
  else
    if parentid >= 0 			# If not the root
      taxa[parentid] =
        pop : 0.0
        cum : 0.0
        cnt : 0
        num : 0
        conf : 0.0
        w_conf : 0.0
        level : levels[level]
        length : 0.0
        #  children : [ taxa[taxid] ]
        sub : {}

      taxa[parentid].sub[name] = taxa[taxid]
  return

process.stdin.on 'data', do ->
  save = ''
  return (c) ->
    lines = c.split '\n'
    lines[0] = save + lines[0]
    save = lines.pop()
    build i for i in lines when i

process.stdin.on 'end', () ->     # Write the JSON structure to stdout.
  taxa[0].length = 1
  # print(JSON.stringify(taxa[0]))
  console.log(JSON.stringify({ sub : { Root : taxa[0] } }))
  process.exit 0
