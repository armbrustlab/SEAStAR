###
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
# along with SEAStAR.  If not, see <http:#www.gnu.org/licenses/>.
# -------------------------------------------------------------------------- #
###

#
# This script reads a nucleotide FASTA files with properly oriented contig 
# sequences in scaffold order (as produced by the json_graph_ops "FASTA" 
# command with the "scaff" parameter set to true). 
#
# It also optionally accepts an alternate FASTA files of assembled contigs to
# be used for "healing" gaps between contigs in the first file.  It is possible
# to use more than one alternate assembly for healing by first concatenating
# those assemblies into one file.  IDs do not need to be unique in the healing
# contigs file.
#
# It produces (to STDOUT) an output FASTA file with one sequence per scaffold. 
#
# This script assumes the node.js environment and should be run with additional
# memory for the V8 javascript engine (for example):
#   --max-old-space-size=1900 --max-new-space-size=2048
# 
# Command line: [options] infile.fna [heal_file.fna] 
#

unless process.version.split('.')[1] >= 10   # Require node version v0.x.y to be x >= 10
   console.error("ERROR: nodejs version v0.10.0 or greater required.")
   process.exit(1) 

verbose = false

fs = require('fs')
zlib = require('zlib')

ss_version = "SS_BUILD_VERSION"   # Note, this gets replaced in the build process

scaff_string = ''  # Accumulate decompressed scaff data into this string
heal_string = ''  # Accumulate decompressed heal data into this string

# Lookup-table for reverse complementing sequences
rc_tab = {'A':'T','T':'A','G':'C','C':'G','X':'X','M':'K','K':'M','R':'Y','Y':'R','W':'S','S':'W','V':'B','B':'V','H':'D','D':'H','N':'N','\n':'\n'}

heal_seq = ''
scaffs = []

rev_comp = (st) -> st.split("").reverse().map((b)->rc_tab[b]).join("")

###
# Format the alternate assembly sequences used to "heal" gaps in the main scaffolds
###

process_heal_fasta = (fasta) ->

   # Make a line separated string containing both the forward and reverse complemented 
   # sequences in the fasta string, removing all FASTA formatting, with one sequence per line.
   heal_seq = (s.split("\n")[1..-1].join("").toUpperCase() for s in fasta.split(">")).join("\n")
   # Add the reverse complement of the above string
   heal_seq = rev_comp(heal_seq) + "\n" + heal_seq

###
# Process the scaffolds and write out the results.
###   

process_scaffs = (fasta) ->

   num_contigs = []
   scaff_cnt = 0
   scaff_nums = {}
   scaff_names = []
   
   for s in fasta.split(">")[1..-1]  # Ignore the first empty string
      [name, scaff_string, seqs...] = s.split(/\n| /)
      
      unless scaff_nums[scaff_string]?
         scaff_nums[scaff_string] = scaff_cnt
         scaff_names[scaff_cnt] = scaff_string
         scaff_cnt++

      scaff_num = scaff_nums[scaff_string]
      
      scaffs[scaff_num] ?= []
      num_contigs[scaff_num] ?= 0
      # Trim ambiguity codes from the ends of scaffolds, requiring 10 non-ambiguity reads
      # at the new ends of each scaffold.
      scaffs[scaff_num].push(seqs.join("").toUpperCase().match(/[ACGT]{10}.*[ACGT]{10}/)[0])
      num_contigs[scaff_num]++

   # Try to heal gaps between contigs using the alternate assembly sequence

   if heal_seq
   
      n = heal_n       # n = required bases of overlap
      m = heal_m       # m = number of bases to back away from each contig end
      max = heal_max   # max = maximum overlap to look for
      
      for scaf, scaf_ind in scaffs when scaf
         for contig, y in scaf when y
         
            prev_contig = scaf[y-1]   
            l = prev_contig.length

            # Build a regex string to use to find sequence in the alt assembly string
            r = "#{prev_contig.substr(-(n+m),n)}([^\n]{1,#{max}})#{contig.substr(m,n)}"
            
            if match = heal_seq.match(r)
               # Determine the overlap bases, set to lowercase
               overlap = match[1].toLowerCase()
               # Merge the previous scaffold and current contig around the aligned overlap
               scaf[y] = prev_contig.substr(0,l-m) + overlap + contig.substr(m)
               delete scaf[y-1]

               if verbose
                  console.warn("**** Found match! %d %d", scaf_ind, y)
                  console.warn("%s", prev_contig.substr(-(n+m)))
                  # Print the alignment to stderr
                  pad = ''
                  for x in [1..n]
                     pad = pad + " "
                  console.warn("%s%s", pad, overlap)
                  pad = ''
                  for x in [1..n+overlap.length-m]
                     pad = pad + " "
                  console.warn("%s%s", pad, contig.substr(0, m+n))

   # Remove "empty contig slots" in the scaffold arrays
   for scaf, x in scaffs when scaf
      scaffs[x] = (contig for contig in scaf when contig?)

   # Try to merge ends of remaining adjacent contigs
   n = overlap_n          # n = required bases of overlap
   m = overlap_m          # m = number of non-ambiguous bases allowed to be trimmed off each contig end
   max = overlap_max      # max = maximum overlap to look for
   for scaf, x in scaffs when scaf
      for contig, y in scaf when y 
         
         prev_contig = scaf[y-1]   
         l = prev_contig.length
   
         # These nested for loops perform the overlap detection
         # Loop for overlaps from max down to n bases.
         for i in [max..n]
            c = 0   # Holds the current number of aligning bases
            # Loop through the bases of a potential i base overlap   
            for j in [0...i]
               # If this position aligns
               if (prev_contig.substr(l - i + j, 1) == contig.substr(j, 1)) 
                  c++   # increment the alignment counter
               else if ((c < n) || (c < i-(2*m))) 
                  # If this position doesn't align, and the previous matches are insufficient
                  # Reset the counter...
                  c = 0
               else   # We've just found a match that ends here. 
                  break;   
      
            # If a match is found in the i base alignment
            if ((c >= n) && (c >= i-(2*m)))
               # Determine the overlap bases, set to lowercase
               overlap = prev_contig.substr(l-i+(j-c), c).toLowerCase()
               # Merge the previous scaffold and current contig around the aligned overlap
               scaf[y] = prev_contig.substr(0, l-i+(j-c)) + overlap + contig.substr(c+(j-c))
               delete scaf[y-1]

               if verbose
                  console.warn("**** Found match! %d %d %d %d %d", x, y, i, j, c)
                  console.warn("%s", prev_contig.substr(l-max+1))
                  pad = ''
                  for x in [1...max-i+(j-c)]
                    pad = pad + " "
                  console.warn("%s%s", pad, overlap)
                  pad = ''
                  for x in [1...max - i]
                     pad = pad + " "
                  console.warn("%s%s",pad, contig.substr(0, max))
                  
               break
            else   # This contig overlap alignment (i bases) did not have a match.
               c = 0
   
   # Remove "empty contig slots" in the scaffold arrays
   for scaf, x in scaffs when scaf
      scaffs[x] = (contig for contig in scaf when contig?)
            
   # Optimize ORFs across remaining gaps
   # These regexes match ORFs at the end and beginning of a config, respectively
   rp = /(?:(?:[^T]..)|(?:T(?:(?:[^AG].)|(?:G[^A])|(?:A[^GA]))))+$/
   rc = /^(?:(?:[^T]..)|(?:T(?:(?:[^AG].)|(?:G[^A])|(?:A[^GA]))))+/
   ns = ['','N','NN']  # Used for shifting frames in the loop below

   for scaf, x in scaffs when scaf
      scafout = scaf[0]
      prev_contig = scaf[0]
      rev_prev_contig = rev_comp(prev_contig)
      for contig, y in scaf when y
         rev_contig = rev_comp(contig)
         pfmax = [-1,-1]
         cfmax = [-1,-1]
         prmax = [-1,-1]
         crmax = [-1,-1]

         # Look for maximum ORFs abutting gap in both contigs in all 6 frames 
         for n in [0..2]
            pf = (prev_contig + ns[n]).match(rp)?[0].length ? 0
            pfmax = [n, pf] if pf > pfmax[1]
            cf = (ns[n] + contig).match(rc)?[0].length ? 0
            cfmax = [n, cf] if cf > cfmax[1]
            pr = (ns[n] + rev_prev_contig).match(rc)?[0].length ? 0
            prmax = [n, pr] if pr > prmax[1]
            cr = (rev_contig + ns[n]).match(rp)?[0].length ? 0
            crmax = [n, cr] if cr > crmax[1]
            
         # Which strand is the maximum ORF on?   
         if (pfmax[1]+cfmax[1] > prmax[1]+crmax[1])
            frame = (pfmax[0]+cfmax[0])%3
         else    
            frame = (prmax[0]+crmax[0])%3

         # Add this contig to the scaffold output with a frame adjusted gap
         scafout = scafout + "nnnnnnnnnnnnnnn" + ns[frame] + contig
         
         prev_contig = contig
         rev_prev_contig = rev_contig

      # Write this scaffold 
      process.stdout.write(">#{scaff_names[x]}\n")
      process.stdout.write(scafout)
      process.stdout.write("\n")

###
# Handle command args and process input files
###

output = process.stdout

cmds = []

# Default settings

heal_fn = null
heal_n = 29       # n = required bases of overlap
heal_m = 5        # m = number of non-ambiguous bases allowed to be trimmed off each contig end 
heal_max = 500    # max = maximum overlap to look for 

overlap_n = 6     # n = required bases of overlap
overlap_m = 0     # m = number of non-ambiguous bases allowed to be trimmed off each contig end
overlap_max = 35  # max = maximum overlap to look for

# Handle reading possible modifications to the above defaults

if process.argv.length < 3
   console.error("No inputs provided, use --help for assistance.")
   process.exit(1)

process.argv.shift()   # drop the program and script names
process.argv.shift()   

while process.argv[0][0..1] is '--' or process.argv[0] is '-h'
   
   p = process.argv.shift()
   parm = p.split('=')

   unless p[1]? or p is '-h' or p is '--help'
      console.error("Invalid parameter: #{p}.  No integer value included.")
      process.exit(1)
   
   switch parm[0]
      when '--help', '-h'
         console.warn("\nUsage: seq_scaffold [options] input.fna")
         console.warn('[options] : [--overlap=<n>] [--trim=<n>] [--max=<n>] [--heal=<heafile.fna>] [--heal_overlap=<n>] [--heal_trim=<n>] [--heal_max=<n>]\n')
         console.warn("--overlap=<n> : Minimum number of bases of overlap required to join contigs. Default: #{overlap_n}")
         console.warn("--max=<n> : Maximum number of bases of overlap to look for. Default: #{overlap_max}")
         console.warn("--trim=<n> : Maximum number of non-ambiguous bases to try trimming from contig ends when looking for overlap. Default: #{overlap_m}")
         console.warn("--heal=<healfile.fna> : Use alternatively assembled contigs in the references FASTA file to try to \"heal\" gaps. Default: No Healing")
         console.warn("--heal_overlap=<n> : Minimum number of bases of overlap required to join contigs. Default: #{heal_n}")
         console.warn("--heal_max=<n> : Maximum number of bases of gap to try to heal. Default: #{heal_max}")
         console.warn("--heal_trim=<n> : Maximum number of non-ambiguous bases to try trimming from contig ends during healing. Default: #{heal_m}")
         console.warn("--verbose : Output diagnostic messages to stderr. Default: #{verbose}\n")
         console.warn("SEASTAR Version: #{ss_version}\n")
         
         process.exit(1)
      when '--heal' then heal_fn = parm[1]   
      when '--heal_overlap' then heal_n = Number(parm[1])
      when '--heal_trim' then heal_m = Number(parm[1])
      when '--heal_max' then heal_max = Number(parm[1])
      when '--overlap' then overlap_n = Number(parm[1])
      when '--trim' then overlap_m = Number(parm[1])
      when '--max' then overlap_max = Number(parm[1])
      when '--verbose' then verbose = true
      else
         console.error("Unknown parameter: #{p}")
         process.exit(1)
   
# Read/process scaffold file
inputfn = process.argv[0]

if inputfn == '-'
   input = process.stdin
   input.resume()
else if inputfn.slice(-3) == '.gz'
   input = fs.createReadStream(inputfn)
             .pipe(zlib.createGunzip())
else
   input = fs.createReadStream(inputfn)

input.on('data', (data) -> scaff_string = scaff_string.concat(data))
     .on('end', () ->
        # Read/process heal file first
        if (heal_fn) 
           if heal_fn == '-'
              input = process.stdin
              input.resume()
           else if heal_fn.slice(-3) == '.gz'
              input = fs.createReadStream(heal_fn)
                        .pipe(zlib.createGunzip())
           else
              input = fs.createReadStream(heal_fn)

           input.on('data', (data) -> heal_string = heal_string.concat(data))
                .on('end', () -> 
                   process_heal_fasta(heal_string)
                   process_scaffs(scaff_string)
                   )
        else 
           process_scaffs(scaff_string)
        ) 



