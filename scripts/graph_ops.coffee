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
#
# This script reads in a JSON formatted ref_select output graph and performs
# various transformations on that file.  See the documentation or use the HELP
# command for additional specific information about each transformation.
#
# This script is written in Coffeescript (http://coffeescript.org/) and 
# assumes the node.js (http://nodejs.org/) execution environment. 
#
# It should be run with additional memory allocated for the V8 javascript engine 
# For example:  node --max-old-space-size=8000 --max-new-space-size=8000 ...
# (The included graph_ops (and nodewrap) scripts do this automatically) 
# Also, graph_ops has a highly fragmented memory allocation pattern that benefits
# greatly from the use of node's --always-compact option, which directs the V8
# Garbage Collector to compact allocated memory on every global sweep. This leads
# to improved runtimes and avoidance of GC corner-cases that occasionally may 
# lead to stability issues.
#
# Usage: [-h|--help] [<input.json[.gz]>] [<script.go[.gz]>] [<command> ['{parms}']...]
# 
# Where: <input.json> is an optional datafile to initially LOAD
#        <script.go> is an optional command SCRIPT file
#        <command> is one of the valid commands
#        '{params}' optionally specify parameters for a given command
#        -h for help with running graph_ops
#
# Command list:  (type HELP <command> for more details about a given command):
#
###

unless process.version.split('.')[1] >= 10   # Require node version v0.x.y to be x >= 10
   console.error("ERROR: nodejs version v0.10.0 or greater required.")
   process.exit(1) 

fs = require('fs')       # Node.js built-in filesystem access module
zlib = require('zlib')   # Node.js built-in zlib (de)compression library 
path = require('path')   # Node.js built-in filename path parsing library
repl = require('repl')   # Node.js built-in command line eval loop library
child_process = require('child_process')   # Node.js built-in process spawning library
JSONStream = require('JSONStream')   # npm library for dealing with JSON streams

ss_version = "SS_BUILD_VERSION"   # Note, this gets replaced in the build process

no_data_str = "No data, ending SCRIPT processing."

###
# This array stores the stack of json_graph structures used by the (UN)STASH commands 
###
stash_stack = []

###
# This function handles a file path that contains '~', resolving it to the user's HOME dir
###
resolve_path = (fn) ->
   fn = process.env.HOME + fn[1..] if fn[0] is '~'  # ~ is a bash thing, so needs to be handled separately
   fn = path.normalize(fn)

###
# Helper function to uniformly handle reading from a file, a gzipped file or stdin 
###
read_input_stream = (fn, parse, cb) ->  # fn = Filename, cb = callback passing object or full data buffer
   if fn is '-'                         # '-' indicates STDIN input stream
      input = process.stdin
      input.resume()              # STDIN needs to be "resumed" in node.js (not opened)
   else 
      inputfn = resolve_path(fn)  # Get the path into standard form
      
      # Open stream for reading
      if inputfn.slice(-3) == '.gz'  # If gzipped file, then pipe through gunzip
         input = fs.createReadStream(inputfn)
                   .pipe(zlib.createGunzip())
      else
         input = fs.createReadStream(inputfn)
      
   if parse   # Return stream parsed JSON instead of string buffer
      obj_buffer = null   
      input.pipe(JSONStream.parse())
           .on('root', (data) ->
              cb(null, data))
           .on('error', (err) ->
               console.error("ERROR: JSON read_input_stream could not be opened or parsed '#{fn}' for input.")
               cb(err, null))
   else
      read_buffer = ''   # Initialize read buffer
      # Handlers for events on the the input stream
      input.on('data', (data) -> 
              read_buffer = read_buffer.concat(data))  # Add data to the input buffer
           .on('end', () ->
              cb(null, read_buffer))     # End of data, invoke the callback with the buffer
           .on('error', (err) ->
              console.error("ERROR: read_input_stream could not open file '#{fn}' for input.")
              cb(err, null))

###
# Helper function crates an object to uniformly handle writing to a file, a gzipped file or stdout
# Treat return like a stream, implements on, write and end methods. Smoothes-over the differences
# between the different output possibilities.
###
open_output_stream = (fn, tag = '') ->   # fn = Filename, returns open stream to write to
   
   out = {}
   
   out.on = (sig, cb) ->
      this.o.on(sig, cb)
   
   out.write = (buf, cb) ->
      this.o.write(buf, cb)
      
   out.end = (cb) ->   
      if this.fh?   # Nothing to do for system streams
         this.o.on("error", (err)->cb(err))  
         this.fh.on("close", ()->setImmediate(cb,null))  # Wait for actual I/O device file to close  
         this.o.end()
      else
         cb(null)
   
   if fn  # If there's a filename, then a stream needs to be created
      fn = resolve_path(fn)   # Get the path into standard form

      # The tag functionality substitutes all discrete runs of @ symbols with an optional tag 
      # passed in from outside of a script, etc.
      if tag_rep = fn.match(/(@+)/)
         console.warn("INFO: Substituting '#{tag}' for '#{tag_rep[1]}' in file #{fn}")
         fn = fn.replace(/@+/g,tag)

      # The renumber functionality below triggers special file renaming logic when the 
      # passed name contains one or more '#' characters in a row.  In this case, these 
      # characters are replaced with an incrementing count, zero-padded to the width of 
      # the number of '#' chars in the original string.
      if renumber = fn.match(/([^#]*)(#+)([^#]*)/)
         inc = 0
         renumber[2] = renumber[2].replace(/#/g,'0')
         fn = renumber[1] + renumber[2].slice(inc.toString().length) + inc.toString() + renumber[3]
         while fs.existsSync(fn)
            inc++
            fn = renumber[1] + renumber[2].slice(inc.toString().length) + inc.toString() + renumber[3]
         console.warn("INFO: Opening numbered file #{fn}")
      else if fs.existsSync(fn)
         console.warn("WARNING: Overwriting file #{fn}")
      
      try
         fout = fs.createWriteStream(fn)  # Create the file
      catch error
         console.error("ERROR: Could not open file; #{fn}")
         return null

      if fout  
         if fn.slice(-3) is '.gz'  # If the filename ends in .gz, then pipe output through gzip
            gz = zlib.createGzip()
            gz.pipe(fout)
            out.o = gz
            out.fh = fout
         else          
            out.o = fout
            out.fh = fout
      else
         out = null
   else  # Otherwise STDOUT is ready to go 
      out.o = process.stdout
   
   out   # Return the output object

##
# my_stringify writes JSON output progressively to a stream, avoiding the creation of a 
# giant text buffer (saving memory for large objects)
##

my_stringify = (j, o, l=2, cb) -> 
   trav_cnt = 0
   max_buffer = 32000
   max_recurse = 25

   traverse = (stack, o, cb, buf = '') ->
      trav_cnt++
      if work = stack.pop()
         if (work.lev and (work.obj instanceof Array))  # Handle traversing array objects
            unless work.keys?
               buf = buf + '['
               work.keys = [0...work.obj.length].reverse()
            else if work.keys.length
               buf = buf + ',\n'
            if work.keys.length is 0
               buf = buf + ']'
               traverse(stack, o, cb, buf)
            else
               k = work.keys.pop()
               stack.push({'obj':work.obj,'keys':work.keys,'lev':work.lev},{'obj':work.obj[k],'lev':work.lev-1})
               traverse(stack, o, cb, buf)
         else if (work.lev and (typeof(work.obj) is 'object'))  # Handle traversing normal objects
            unless work.keys?
               buf = buf + '{'
               work.keys = Object.keys(work.obj).reverse()
            else if work.keys.length
               buf = buf + ',\n'
            if work.keys.length is 0
               buf = buf + '}'
               traverse(stack, o, cb, buf)
            else
               k = work.keys.pop()
               stack.push({'obj':work.obj,'keys':work.keys,'lev':work.lev},{'obj':work.obj[k],'lev':work.lev-1})
               buf = buf + "\"#{k}\":"
               traverse(stack, o, cb, buf)
         else  # Handle all other cases
            buf = buf + JSON.stringify(work.obj)    # Buffer for performance, but
            if buf.length > max_buffer or trav_cnt > max_recurse  # guard against excess recursion 
               o.write(buf, (err) ->
                  if err
                     cb(err)
                  else  # Let the I/O get out before next call
                     setImmediate(traverse, stack, o, cb))
            else  # Pass on writing for now
               traverse(stack, o, cb, buf)
      else  # No more work, so clean up and invoke the callback    
         o.write(buf + '\n', (err) ->
            if err
               cb(err)
            else  # Let the I/O get out before next call
               setImmediate(cb, null))
      trav_cnt--

   # Invoke the traverse function with the input it expects
   traverse([{"obj":j,"lev":l}], o, cb)

##
# clone_object makes a deep copy of a JavaScript object (assumes JSON compliant, no cycles, etc)
##

clone_object = (obj) ->
   clone = if (obj instanceof Array) then [] else {}
   for i of obj
      if typeof(obj[i]) is "object"
         clone[i] = clone_object(obj[i])
      else
         clone[i] = obj[i]
   clone

###
# Calc_seq_stats calculates the N50, mean %gc and coverage and contig length statistic 
# for the nodes in the passed-in object. It also returns the total length of all sequences
###
calc_seq_stats = (nodes) ->

   seq_total = 0
   seq_list = []
   cov = 0.0
   gc = 0.0

   for name, n of nodes   # Walk the nodes accumulating stats
      seq_total += n.seq_len
      cov += n.cov*n.seq_len
      gc += (n.pct_gc ? 0)*n.seq_len
      seq_list.push(n)   # Push the node to a list for calculating N50
      
   seq_list.sort((a,b) -> b.seq_len-a.seq_len)   # Descending sort

   running_total = 0

   for l in seq_list   # Walk the sequence list looking for the point where 50% of the 
                       # sequence is accounted for 
       running_total += l.seq_len
       if running_total >= seq_total/2  # Stop when this point is reached
          break
          
   # Half of the sequence is in contigs equal or larger than l  (= N50)     
   [l.seq_len, seq_total, cov/seq_total, gc/seq_total, seq_list[0].name]  

###
# cc_seq_len calculates the total sequence length for all nodes in the provided cc
###
cc_seq_len = (j,cc) -> # j = graph object, cc = single connected component, ie. an array of node ids
   sum_len = 0
   for nid in cc
      sum_len += j.nodes[nid].seq_len
   sum_len   

###
# Bisect returns a function that when called with an accessor function f,
# returns a function that takes an array a and a value x, and returns the
# position in the array where x should be inserted (assuming an ascending 
# sorted array).
###
bisect = (f) ->

   (a, x, lo, hi) ->
   
      lo = 0 if (arguments.length < 3) 
      hi = a.length if (arguments.length < 4) 
      while (lo < hi) 
         mid = lo + hi >> 1
         if (x < f.call(a, a[mid], mid))
            hi = mid 
         else 
            lo = mid + 1
      lo
###
# my_bisect is a function to efficiently find a position in an array of elements sorted by score
###
my_bisect = bisect((d) -> d[0].score)

###
# heuristic scoring function that biases bitscores toward edges that connect
# nodes with similar GC% and Coverage. (Mostly for metagenome assembly) 
###

# Using this the default. It may be (globally) disabled with the "bits" parm of the MST or SST command
use_heuristic = true   

heuristic = (a_gc,b_gc,a_cov,b_cov) -> 
   ((0.08*Math.abs(a_gc - b_gc)) + (0.125*Math.abs(Math.log(a_cov,2) - Math.log(b_cov,2))))

###
# Setup the internal references between nodes and edges in the graph 
# This function is necessary because the JSON representation of the assembly graph 
# cannot contain any references between sub-objects of the datastructure
###
build_graph_refs = (j, nf = null, ef = null) ->
# j = graph object, nf = function to run on all nodes, ef = function to run on all edges

   n_out = (for name, n of j.nodes  # Walk through all of the nodes
      n.id = name
      n.inlinks = []   # These empty arrays are populated below in edge processing 
      n.outlinks = []
      n.links = []
      nf? n
      n)

   e_out = (for e, i in j.edges   # Walk through all of the edges
      # Check to make sure that this edge is consistent with the current nodeset
      # If not, there is a bug somewhere else and this error is fatal.
      unless e.src = j.nodes[e.n1]
         throw "build_graph_refs: Missing node #{e.n1} in edge with #{e.n2}."
      unless e.tar = j.nodes[e.n2]
         throw "build_graph_refs: Missing node #{e.n2} in edge with #{e.n1}."

      # Calculate the edge bitscore, optionally using the %gc and coverage heuristic
      unless e.score? 
         e.score = e.bits
         if use_heuristic
            if e.src.pct_gc? and e.tar.pct_gc?
              s_gc = e.src.pct_gc
              t_gc = e.tar.pct_gc
            else
              s_gc = 0
              t_gc = 0
               
            e.score -= e.bits*heuristic(s_gc, t_gc, e.src.cov, e.tar.cov)
         
      e.index = i
      e)

   # Sort the output edge list by score, several algorithms depend on this (e.g. MST)
   e_out.sort((a,b) -> (b.score - a.score))

   for e in e_out  # Walk the output edges and add edge references to connected nodes
      l_out = [e, e.tar]
      e.src.outlinks.push(l_out)
      e.src.links.push(l_out)
      l_in = [e, e.src]
      e.tar.inlinks.push(l_in)
      e.tar.links.push(l_in)
      ef? e

   [n_out, e_out]   # return the node and edge lists (side effects also accessable in j)

###
# Setup the internal references between removed nodes and edges in the graph 
# Similar to build_graph_refs above, but works on the "removed" pool of nodes and edges
# Including connecting them to "selected" nodes
###
build_removed_graph_refs = (j, nf = null, ef = null) ->
# j = graph object, nf = function to run on all nodes, ef = function to run on all edges

   j.removed_nodes ?= {}   # Create the "removed" datastructures if they are missing.
   j.removed_edges ?= []

   for name, n of j.nodes  # Walk through all of the selected nodes
      n.rem_inlinks = []   # These empty arrays are populated below in edge processing
      n.rem_outlinks = []
      n.rem_links = []
      nf? n

   n_out = (for name, n of j.removed_nodes   # Walk through all of the removeded nodes
      n.id = name
      n.rem_inlinks = []    # These empty arrays are populated below in edge processing
      n.rem_outlinks = []
      n.rem_links = []
      nf? n
      n)

   e_out = (for e, i in j.removed_edges    # Walk through all of the removed edges
      # Check to make sure that this edge is consistent with the current nodeset(s)
      # If not, there is a bug somewhere else and this error is fatal.
      unless e.src = j.nodes[e.n1] or j.removed_nodes[e.n1]
         throw "build_removed_graph_refs: Missing node #{e.n1} in removed edge with #{e.n2}."
      unless e.tar = j.nodes[e.n2] or j.removed_nodes[e.n2]
         throw "build_removed_graph_refs: Missing node #{e.n2} in removed edge with #{e.n1}."
         
      # Calculate the edge bitscore, optionally using the %gc and coverage heuristic   
      unless e.score? 
         e.score = e.bits
         if use_heuristic
            if e.src.pct_gc? and e.tar.pct_gc?
              s_gc = e.src.pct_gc
              t_gc = e.tar.pct_gc
            else
              s_gc = 0
              t_gc = 0
               
            e.score -= e.bits*heuristic(s_gc, t_gc, e.src.cov, e.tar.cov)

      e.index = i
      e)

   # Sort the output edge list by score, several algorithms depend on this (e.g. MST)
   e_out.sort((a,b) -> (b.score - a.score))

   for e in e_out   # Walk the output edges and add edge references to connected nodes
      l_out = [e, e.tar]
      e.src.rem_outlinks.push(l_out)
      e.src.rem_links.push(l_out)
      l_in = [e, e.src]
      e.tar.rem_inlinks.push(l_in)
      e.tar.rem_links.push(l_in)
      ef? e

   # return the removed nodes, and removed / shared edge lists (side effects also accessable in j)
   [n_out, e_out]  

###   
# Remove the internal references between nodes and edges in the graph 
# This undoes all of the linking performed by the "build" functions above, leaving the 
# graph datastructure in a state that can be serialized into JSON
# This function removes references for both selected and removed nodes 
###
remove_graph_refs = (j, nf = null, ef = null) ->
# j = graph object, nf = function to run on all nodes, ef = function to run on all edges
   for name, n of j.nodes
      delete n.id
      delete n.links
      delete n.inlinks
      delete n.outlinks
      delete n.rem_links
      delete n.rem_inlinks
      delete n.rem_outlinks      
      nf? n
      
   j.edges = (for e in j.edges when e?
      delete e.src
      delete e.tar
      delete e.index
      ef? e
      e)

   if j.removed_nodes?
      for name, n of j.removed_nodes
         delete n.id
         delete n.links
         delete n.inlinks
         delete n.outlinks
         delete n.rem_links
         delete n.rem_inlinks
         delete n.rem_outlinks      
         nf? n
   
   if j.removed_edges?
      j.removed_edges = (for e in j.removed_edges when e?
         delete e.src
         delete e.tar
         delete e.index
         ef? e
         e)

###
# Calculate the current set of connected components
###
calc_ccomps = (j, args = { "sortby" : "nodes" }, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done
   if args.help?
      console.warn("
#{args.help} -- Calculate the current graph connected components
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
sortby : <string> -- Specify how to sort the resulting list of CCs.  \"nodes\" or \"sequences\"\n
\n
        Example: #{args.help} {\"sortby\":\"nodes\"} -- Default. Sort CCs in descending\n
        order of number of nodes.\n
\n
        Example: #{args.help} {\"sortby\":\"sequences\"} -- Sort CCs in descending order\n
        of amount of sequence.\n
")
      callback?(null, j)
      return

   # Recalculating ccomps when scaffolds exist is a sign that the graph structure may 
   # have changed, invalidating the existing scaffolds (and clusters based on them...)
   if j.scaffolds?
      console.warn("WARNING: Previous operation changed connected component structure, so existing scaffolds and clusters are being removed.\n")
      delete j.scaffolds
      if j.clusters?
         delete j.clusters

   j.removed_nodes ?= {}

   if args.sortby? and ["nodes","sequences"].indexOf(args.sortby) is -1 
      callback?(new Error('sortby argument must equal either "nodes" or "sequences".'), null)
      return

   try
      [nodes] = build_graph_refs(j) 
   catch err
      callback?(err, null)
      return
            
   j.connected_comps = []

   nodes_avail = {}
   nodes_avail[node.id] = node for node in nodes
   n = nodes[0]
   edge_cand = []
   nodes_added = []

   # Execute a breadth first search of all "non-removed" nodes, to assign nodes to CCOMPS

   while n      
      delete nodes_avail[n.id]
      nodes_added.push(n.id)   

      # Walk through all of the (sorted) links for this node
      # following only those leading to unvisited nodes
      for l in n.links when nodes_avail[l[1].id]?  
         edge_cand.push(l)

      # Take the first link to an unvisited node (perhaps one from a previously visited
      # node, edge_cand is a stack). nodes_avail is checked again below because this node
      # may have been visited via another path in the meantime since l was placed on the 
      # stack
      l = edge_cand.pop()     
      while l and not nodes_avail[l[1].id]?
         l = edge_cand.pop()
               
      # Choose the next node (either via link, or in another unvisited ccomp)         
      if l
         n = l[1]
      else # On to a new component
         # Sort the node ids in the ccomp in descending order by sequence length, then
         # alphabetically by node id with two have the same length
         nodes_added.sort((a,b) -> 
              lendiff = j.nodes[b].seq_len - j.nodes[a].seq_len
              if lendiff is 0
                 if j.nodes[a].id > j.nodes[b].id
                    return 1
                 else
                    return -1 
              else
                 return lendiff
            )
         j.connected_comps.push(nodes_added)
         nodes_added = []
         n = nodes.pop()
         while n and not nodes_avail[n.id]?
            n = nodes.pop()   
   
   remove_graph_refs(j)   
  
   if args.sortby is "sequences"
      # Sort connected components in descending order of amount of sequence in nodes
      j.connected_comps.sort((a,b) -> (cc_seq_len(j,b) - cc_seq_len(j,a)))
   else 
      # Sort connected components in descending order of number of nodes
      j.connected_comps.sort((a,b) -> (b.length - a.length))
 
   callback?(null, j)

###
## Used to sort out the various types of MP edges and how they orient the underlying contigs
###

edge_director = (prev_n, n, e) -> 
   # prev_n is a node, n is another node, e is the edge that connects them
   
   # The logic coded in this function assumes:
   #
   # -- The prev_n node's orientation relative to the reference strand is known,
   #    that is, its "ref_st" property is set to the correct value, true if "on the 
   #    reference strand", false if not.
   # 
   # -- Which node is n1 and which is n2 in the edge "e" is arbitrary at this point
   # 
   # -- The relative orientation of the nodes to each other is known:
   #
   # Meaning of the three possible "e.dir" orientation codes FB, FF and BB 
   # (assigned by refselect):
   # 
   # F and B refer to whether the connecting mate-pairs point off of the 3' (Forward) 
   # or 5' (Backward) end of the node sequence.
   #
   # FB (Forward-Backward) is the correct orientation of two adjacent contigs:
   #
   #     Mates beginning in contig X:          ----F--->  
   #                         Contigs:  5' X======>   Y======> 3'
   #     Mates beginning in contig Y:          <---B---- 
   #
   #  NOTE!  This is equivalent to Backward-Forward so there is no code for "BF" 
   #         since it would describe same relative orientation (just rev comp'ed).
   #         So the below example is also "FB":
   #
   #     Mates beginning in contig X:          ----B--->  
   #                         Contigs:  3' <======X   <======Y 5'
   #     Mates beginning in contig Y:          <---F---- 
   # 
   #
   # FF (Forward-Forward) is the situation where the nodes face toward each other on  
   # opposite strands (invariant of ref strand):
   #
   #     Mates beginning in contig X:          ----F--->  
   #                         Contigs:     X======>   <======Y 
   #     Mates beginning in contig Y:          <---F---- 
   #
   #
   # BB (Backward-Backward) is the situation where the nodes face away from each other on  
   # opposite strands (invariant of ref strand):
   #
   #     Mates beginning in contig X:          ----B--->  
   #                         Contigs:     <======X   Y======> 
   #     Mates beginning in contig Y:          <---B---- 
   #
   #
   # Given all of this information, the task of this function is to change around the
   # parameters of the *edge* only, so that the nodes are properly oriented in the FB 
   # orientation.  Any necessary reverse complementation of actual sequence will be 
   # performed later based on the .ref_str property of each node.

   # In this case, this edge is already properly directed so there's no work to do
   if e.dir is "forward" or e.dir is "pos"
      return e

   e.org_dir = e.dir      # Save the original direction of this edge.
   e.p1 = Math.abs(e.p1)  # Only want absolute values of positions for below calcs
   e.p2 = Math.abs(e.p2)
   
   # This switch statement builds a string out of a list of three parameters because 
   # there is no testable equivalence of arrays in JavaScript, and comparing strings 
   # is more efficient anyway.
   switch [e.dir, prev_n.ref_str, e.n2 is n.id].join("_")

      #  Already properly oriented, both on ref strand 
      #    p==+==> ----> n==+==>   n==+==> ----> p==+==>
      #    (this is what all of the other cases below are working to achieve)
      when    "FB_true_true",       "FB_true_false"   

         if not n.ref_str? or n.ref_str 
            n.ref_str = true   # n is also on the ref_str
            e.dir = "forward"
            e
         else
            null  # Something is wrong!

      #  Already properly oriented, but on wrong strand, so they both need to be reversed 
      #  and swapped in the edge:
      #    p==-==> ----> n==-==>   n==-==> ----> p==-==>
      #      rev     X     rev       rev     X     rev    (rev = reverse node, X = swap positions in edge)
      #    <==+==p <---- <==+==n   <==+==p <---- <==+==n
      
      when   "FB_false_true",        "FB_false_false"   

         if not n.ref_str? or not n.ref_str 
            # Swap the nodes in the edge, and recalculate the mate-pair positions
            # for when the nodes are reverse complemented
            n.ref_str = false
            [e.p1, e.p2] = [e.tar.seq_len-e.p2-1, e.src.seq_len-e.p1-1]
            [e.n1, e.n2] = [e.n2, e.n1]
            [e.src, e.tar] = [e.tar, e.src]
            e.dir = "forward"
            e
         else
            null   # Something is wrong!

      #  Nodes are not properly oriented, but only one of them needs to be reversed
      #  to achieve a proper FB orientation on the correct strand:
      #    p==+==> ----> <==-==n   n==+==> ----> <==-==p    <==-==n ----> p==+==>  <==-==p ----> n==+==>
      #                    rev                     rev        rev                    rev
      #    p==+==> ----> n==+==>   n==+==> ----> p==+==>    n==+==> ----> p==+==>  p==+==> ----> n==+==>
      when    "FF_true_true",       "FF_false_false",       "BB_true_false",      "BB_false_true"

         if not n.ref_str? or n.ref_str is not prev_n.ref_str 
            n.ref_str = not prev_n.ref_str  #  Nodes are on opposite strands
            if prev_n.ref_str ^ (e.n2 is n.id)  # Determine which node to reverse
               e.p1 = e.src.seq_len-e.p1-1
            else
               e.p2 = e.tar.seq_len-e.p2-1
            e.dir = "forward"
            e
         else
            null   # Something is wrong!

      #  Nodes are not properly oriented, but reversing one of them puts them in the 
      #  wrong orientation within the edge, so they also need to be swapped to achieve 
      #  a proper FB orientation on the correct strand:
      #    n==-==> ----> <==+==p   p==-==> ----> <==+==n   <==+==p ----> n==-==>   <==+==n ----> p==-==>
      #      rev     X               rev     X                       X     rev               X     rev
      #    <==+==n <---- <==+==p   <==+=== <---- <==+==n   <==+==p <---- <==+===   <==+==n <---- <==+===
      when   "FF_true_false",       "FF_false_true",        "BB_true_true",       "BB_false_false"

         if not n.ref_str? or n.ref_str is not prev_n.ref_str 
            n.ref_str = not prev_n.ref_str  #  Nodes are on opposite strands
            if prev_n.ref_str ^ (e.n2 is n.id)  # Determine which node to reverse
               e.p1 = e.src.seq_len-e.p1-1
            else
               e.p2 = e.tar.seq_len-e.p2-1   
            # Swap the nodes in the edge   
            [e.p1, e.p2] = [e.p2, e.p1]
            [e.n1, e.n2] = [e.n2, e.n1]
            [e.src, e.tar] = [e.tar, e.src]
            e.dir = "forward"
            e
         else
            null   # Something is wrong!
         
      else
         null   # Something is wrong!

###
# reverse_comp - Reverse complement a sequence 
###

reverse_comp = (seq) ->
   # Look-up table for reverse complement of individual DNA bases
   rc_tab = { 'A':'T', 'T':'A', 'G':'C', 'C':'G', 'X':'X', 'M':'K', 'K':'M', 'R':'Y', 'Y':'R', 'W':'S', 'S':'W', 'V':'B', 'B':'V', 'H':'D', 'D':'H', 'N':'N' }
   seq.toUpperCase().split("").reverse().map((b)->rc_tab[b]).join("")

###
# Make each connected component a directed graph, properly orienting the 
# nodes relative to each other
###

make_directed = (j) ->

   # This function assumes that the graph is acyclic (in a collection of trees) as 
   # produced by the MST or SST algorithms

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

   try
      build_graph_refs(j) 
   catch err
      callback?(err, null)
      return

   # For each connected component
   for cc in j.connected_comps
   
      ccnodes = {}
      edge_cand = []

      # Look for the largest node (by sequence length) 
      max_len = 0
      max_node = null

      for nid in cc
         n = j.nodes[nid]
         ccnodes[nid] = n
         if n.seq_len > max_len
            max_node = n
            max_len = n.seq_len            
   
      # max_node is now deemed to be the "reference strand"
      n = max_node
      n.ref_str = true
      
      # Walk though all of the nodes in this ccomp
      while n

         delete ccnodes[n.id]  # Mark this node as visited (for detecting cycles)
         
         for l in n.links when ccnodes[l[1].id]?    # Walk through all of the links for this node
            l.push(n)
            edge_cand.unshift(l)

         l = edge_cand.pop()  # Next link to traverse (breadth first)    

         n = null
         
         if l
            if not ccnodes[l[1].id]
               throw new Error("make_directed: Cycle detected at node: #{l[1].id}")

            n = l[1]

            # Reorient the edge/nodes to be on the same strand
            unless e = edge_director(l[2], n, l[0])
               throw new Error("make_directed: Improper edge type detected between nodes: #{n.id} and #{l[2].id}")

            # Reverse compliment this node's sequence if it's not on the reference strand
            unless n.ref_str or not n.recon_seq?
               n.recon_seq = reverse_comp(n.recon_seq)
                 
   j.digraph = true
   
   remove_graph_refs(j)

###
## Spanning tree processing using Primm's algorithm
###

maximal_spanning_tree = (j, args = {}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Calculate the Maximal Spanning Tree of all connected components
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
bits : true -- Use raw connection bitscores and not GC% / Coverage adjusted bitscores\n
\n
        Example: #{args.help} {\"bits\":true} -- Use raw connection bitscores from mate-pairing\n
")
      callback?(null, j)
      return

   if args.bits
      use_heuristic = false

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

   j.removed_edges ?= []   # initialize a place to keep removed edges in the output

   try
      [nodes, links] = build_graph_refs(j) 
   catch err
      callback?(err, null)
      return
   
   # Start the MST 

   edges_added = []

   # Walk though each ccomp
   for c in j.connected_comps
      n = j.nodes[c[0]]
      nodes_added = {}
      nodes_added[n.id] = n
      nodes_seen = {}
      edge_cand = []

      # Walk through all of the nodes 
      while n      
         # For each link from this node
         for l in n.links  
            if (not nodes_added[l[1].id]?)
               # If the other node is not reached yet, add this link in sorted order by score
               edge_cand.splice(my_bisect(edge_cand, l[0].score),0,l)
            else unless l[0].done?
               # Otherwise, remove this edge (unless it's already been removed)
               l[0].done = true
               j.removed_edges.push(l[0])

         # Find the unvisited node with the highest scoring edge leading to it 
         next_n = edge_cand.pop()
         while next_n and nodes_added[next_n[1].id]?
            unless next_n[0].done?
               next_n[0].done = true
               j.removed_edges.push(next_n[0])
            next_n = edge_cand.pop()
         
         # if found, keep this edge and setup for the next iteration
         if next_n
            [new_link, n] = next_n
            nodes_added[n.id] = n  
            edges_added.push(new_link)
            new_link.done = true
         else
            n = null
            
   # Keep only the tree edges         
   j.edges = edges_added

   remove_graph_refs(j, null, ((ed) -> delete ed.done))   

   try
      make_directed(j)
   catch err
      callback?(err, null)
      return
   
   callback?(null, j)

###
# Scaffolding spanning tree processing, strives to minimize parallel branches when 
# nodes are short relative to pairing insert size
###
scaffold_spanning_tree = (j, args = {}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Calculate the improved Scaffold Spanning Tree of all connected components
")
      if args.detailed_help?
         console.warn("
\n
        NOTE: This command is generally preferable to the MST command when median contig\n
        length is less than the mean distance between paired reads. That is, when a\n
        relatively large insert size was selected.\n 
\n
Parameters:\n
\n
bits : true -- Use raw connection bitscores and not GC% / Coverage adjusted bitscores\n
\n
        Example: #{args.help} {\"bits\":true} -- Use raw connection bitscores from mate-pairing\n
")
      callback?(null, j)
      return

   if args.bits
      use_heuristic = false

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

   j.removed_edges ?= []   # initialize a place to keep removed edges in the output

   try
      [nodes, links] = build_graph_refs(j, (n) -> n.int_mp_bits = 0.0) 
   catch err
      callback?(err, null)
      return      
   
   # Note the strength of internal edges (indicative of nodes longer than paring insert size)
   for e in j.internal_edges when j.nodes[e.n1]?
      j.nodes[e.n1].int_mp_bits = e.bits

   # Heuristic function to determine if a node needs to be part of the scaffold backbone
   # (to effectively reach nodes further down the line) or if it can safely be relegated 
   # to a (very short) branch
   new_territory = (node) ->
      next_l = null
      for l in node.links
         # See if this node has an "unseen" child or (heuristically) if it's a large 
         # enough contig that it has to be part of the "main" scaffold backbone because 
         # any links managing to hop over it will be statistically weak
         if not nodes_added[l[1].id]? and (not nodes_seen[l[1].id]? or (l[1].int_mp_bits / l[1].links[0][0].bits > 0.02))
            next_l = l
            break
      next_l

   # Start the SST 

   edges_added = []

   # Walk through each ccomp
   for c in j.connected_comps
   
      n = j.nodes[c[0]]
      nodes_added = {}
      nodes_added[n.id] = n
      nodes_seen = {}
      edge_cand = []

      next_n = true

      # For each node in this ccomp
      while next_n      
         next_n = null
         new_link = null
         
         if n  # Skip this step if this node has had its links discarded
            for l in n.links  
               if (not nodes_added[l[1].id]?) and (not (nodes_seen[l[1].id]?) or (l[0].score > nodes_seen[l[1].id][0].score))
                  # If the other node is not reached yet, add this link in sorted order by
                  # score if this is the strongest link to this node that has been seen 
                  nodes_seen[l[1].id] = l  # This is now the best link to this node
                  edge_cand.splice(my_bisect(edge_cand,l[0].score),0,l)
               else unless l[0].done?
                  # Otherwise, remove this edge (unless it's already been removed)
                  l[0].done = true
                  j.removed_edges.push(l[0])

         # Find the unvisited node with the highest scoring edge leading to it
         next_n = edge_cand.pop()       
         while next_n and nodes_added[next_n[1].id]?
            unless next_n[0].done?
               next_n[0].done = true
               j.removed_edges.push(next_n[0])
            next_n = edge_cand.pop()

         # if found, keep this edge and setup for the next iteration
         if next_n
            [new_link, n] = next_n
            nodes_added[n.id] = n  
            edges_added.push(new_link)
            new_link.done = true
            # Unless this node reaches "new territory" (unseen nodes), discard all of its other links
            unless new_territory(n)
               for l in n.links
                  unless l[0].done?
                     l[0].done = true
                     j.removed_edges.push(l[0])
               n = null
   
   j.edges = edges_added

   remove_graph_refs(j, ((nd) -> delete nd.int_mp_bits), ((ed) -> delete ed.done))   

   try
      make_directed(j)
   catch err
      callback?(err, null)
      return
   
   callback?(null, j)

###
# Pluck leaves off the tree, the purpose is to remove short branches before pruning
# to avoid unnecessarily splintering the tree into thousands of branches shorter than
# the mean pairing insert length
###
remove_leaves = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Remove all leaf contig nodes (in or outdegree == 0) from the graph
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
iterate : <int> -- number of iterations of #{args.help} to run. Each is equivalent to\n
       running #{args.help} again as a separate command.  By Default <int> = 3\n
\n
        Example: #{args.help} {\"iterate\":2} -- Equivalent to running:
           #{args.help} {\"iterate\":1}\n
           #{args.help} {\"iterate\":1}\n
\n
min_len : <int> -- Minimum length of isolated nodes to keep.\n
\n
        Example: #{args.help} {\"min_len\":20000} -- Default. Do not \"pluck away\" an\n
        unconnected node containing 20000 or more bases of sequence.
\n")
      callback?(null, j)
      return

   args.iterate ?= 3
   args.min_len ?= 20000

   delete j.connected_comps   # This invalidates the cached ccomps

   j.removed_edges ?= []   # initialize a place to keep removed edges in the output
   j.removed_nodes ?= {}

   j.pluck_iterations ?= 1  # Keep track of home many times PLUCK has been run

   # Number of times to pluck the tree (removing only the ends at each iteration) 
   for i in [j.pluck_iterations..args.iterate+j.pluck_iterations-1]

      try
         [nodes, edges] = build_graph_refs(j) 
      catch err
         callback?(err, null)
         return

      # Remove short single node scaffolds in the first iteration only 
      if i is 1
         # Walk all nodes looking for nodes with zero in or out degree
         for n in nodes when (n.links.length is 0) and (n.seq_len? < args.min_len)
            j.removed_nodes[n.id] = n   # Remove the node itself
            delete j.nodes[n.id]

      # Walk all nodes looking for leaves (degree = 1) with non-removed connected nodes
      for n in nodes when (n.links.length is 1) and (not (n.links[0][1].id in j.removed_nodes))
         j.removed_nodes[n.id] = n   # Remove the node itself
         delete j.nodes[n.id]

      edges_kept = []

      # Keep only edges connecting the remaining nodes 
      for e in edges
         if j.nodes[e.n1]? and j.nodes[e.n2]?
            edges_kept.push(e)
         else
            e.pluck_iteration = i
            j.removed_edges.push(e)
         
      j.edges = edges_kept
            
      remove_graph_refs(j)   
   
   j.pluck_iterations = args.iterate+j.pluck_iterations-1
   callback?(null, j)

###
# Add ends
###
add_ends = (j, args = {}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Extends scaffold ends (reversing the action of PLUCK at scaffold ends)
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
iterate : <int> -- number of iterations of #{args.help} to run. Each is equivalent to\n
       running #{args.help} again as a separate command.  By Default <int> = 3\n
\n
        Example: #{args.help} {\"iterate\":2} -- Equivalent of #{args.help} #{args.help}\n
")
      callback?(null, j)
      return

   args.iterate ?= 3

   unless j.pluck_iterations?
      console.warn("WARNING: PUSH called on unPLUCKed graph.")
      callback?(null, j)
      return
   
   args.iterate = j.pluck_iterations if args.iterate > j.pluck_iterations   

   for i in [j.pluck_iterations..j.pluck_iterations-args.iterate+1]

      calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

      try
         build_graph_refs(j) 
         build_removed_graph_refs(j)
      catch err
         callback?(err, null)
         return

      for c in j.connected_comps
   
         # Find the "head" and "tail" nodes for this scaffold
      
         try 
            [[head_node], [tail_node]] = find_ends(j, c, true)
         catch err
            remove_graph_refs(j, null, ((e) -> delete e.pluck_iteration))
            callback?(err, null)
            return
            
         unless head_node.no_push
            for e in head_node.rem_inlinks when e[0].dir is "forward" and e[0].pluck_iteration is i and not j.nodes[e[0].n1]? 
               j.edges.push(e[0])
               j.nodes[e[1].id] = e[1]
               delete j.removed_nodes[e[1].id]
               delete j.removed_edges[e[0].index]
               break
         
         unless tail_node.no_push
            for e in tail_node.rem_outlinks when e[0].dir is "forward" and e[0].pluck_iteration is i and not j.nodes[e[0].n2]?
               j.edges.push(e[0])
               j.nodes[e[1].id] = e[1]
               delete j.removed_nodes[e[1].id]
               delete j.removed_edges[e[0].index]
               break

      remove_graph_refs(j) unless i is 1

      delete j.connected_comps 
   
   remove_graph_refs(j, null, ((e) -> delete e.pluck_iteration))
   j.pluck_iterations = j.pluck_iterations-args.iterate
   delete j.pluck_iterations if j.pluck_iterations is 0
   callback?(null, j)

###
## prune the tree by cutting off branches
###
cut_branches = (j, args={"verbose":false}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Split the assembly graph at all contig nodes with in or out degree > 1
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
strict : true -- always cut all but the strongest link when there are > 2 in+out links\n
\n
        Example: #{args.help} {\"strict\":false} -- The default case, attempt to determine\n
        which branch(es) to remove to preserve one high weight path across a multiply\n
        linked contig node.\n
\n
        Example: #{args.help} {\"strict\":true} -- Always strictly prune. Equivalent to\n
        {\"ratio\":0.0} below\n
\n
ratio : <float> -- ratio of the strongest to next strongest links to trigger strict pruning\n
\n
        Example: #{args.help} {\"ratio\":0.0} -- Always strictly prune.\n
\n
        Example: #{args.help} {\"ratio\":1.0} -- Never strictly prune.\n
\n
        Example: #{args.help} {\"ratio\":0.2} -- Default. Don't strictly filter when the\n
        score of the highest scoring edge is >= 5x greater than the next highest scoring\n
        link of the same direction (for both in and outlinks)\n
\n
verbose : true -- output diagnostics on STDERR\n
\n
        Example: #{args.help} {\"verbose\":true} -- Generate extra output information\n
")
      callback?(null, j)
      return

   args.strict ?= false  #  If true, always cut all but the strongest link when there are more than 2 in+out links
   args.ratio ?= 0.2     #  Determines the ratio of the strongest to next strongest links to trigger strict pruning 
                         #  e.g. ratio = 0.0 --> always strictly filter.  ratio = 1.0, never strictly filter  
                         #  ratio = 0.2 (default) --> don't strictly filter when the score of the highest scoring edge is >= 5x
                         #  greater than the next highest scoring link of the same direction (for both in and outlinks)

   delete j.connected_comps   # This invalidates the cached ccomps

   j.removed_edges ?= []   # initialize a place to keep removed edges in the output

   try
      [nodes, edges] = build_graph_refs(j) 
   catch err
      callback?(err, null)
      return

   edges_kept = []

   # Walk all nodes looking for nodes with indegree > 1 and/or outdegree > 1
   for n in nodes
    
      if n.inlinks.length > 1
         in_ratio = (n.inlinks[1][0].score / n.inlinks[0][0].score) > args.ratio
         console.warn("Inlinks: #{n.name} #{n.inlinks.length} #{(n.inlinks[1][0].score / n.inlinks[0][0].score)} #{in_ratio}") if args.verbose
         # Mark all but the strongest inlink for removal
         for l in n.inlinks[1..]
            l[0].remove = true
      else
         in_ratio = false  # If there is only one link, the ratio test is automatically met

      if n.outlinks.length > 1
         out_ratio = (n.outlinks[1][0].score / n.outlinks[0][0].score) > args.ratio
         console.warn("Outlinks: #{n.name} #{n.outlinks.length} #{(n.outlinks[1][0].score / n.outlinks[0][0].score)} #{out_ratio}") if args.verbose
         # Mark all but the strongest outlink for removal
         for l in n.outlinks[1..]
            l[0].remove = true 
      else
         out_ratio = false  # If there is only one link, the ratio test is automatically met

      # Strictly filter when it is turned on, and both in and out edges exist, and there were originally more than 2 edges, and the ratio test isn't met
      if args.strict and n.inlinks.length and n.outlinks.length and (n.inlinks.length + n.outlinks.length > 2) and (in_ratio or out_ratio)
         n.no_push = true  # Prevent nodes from being pushed onto this one, since we don't know what is correct.  
         if n.inlinks[0][0].score > n.outlinks[0][0].score
            n.outlinks[0][0].remove = true
            console.warn("Strictly Removing #{n.outlinks[0][0].n1} --> #{n.outlinks[0][0].n2}")
         else    
            n.inlinks[0][0].remove = true
            console.warn("Strictly Removing #{n.inlinks[0][0].n1} --> #{n.inlinks[0][0].n2}")

   # Walk all edges, committing the edge removals marked above
   for e in edges
      if e.remove?
         j.removed_edges.push(e)
         delete e.remove
      else 
         edges_kept.push(e)
            
   j.edges = edges_kept

   remove_graph_refs(j)   
   
   callback?(null, j)

###
## Returns the nodes at the beginning and end of a linear scaffold chain [head, tail] by ccomp
###
find_ends = (j, ccomp, strict = false)->  
      head_nodes = []
      tail_nodes = []
      for cn in ccomp 
         if j.nodes[cn].inlinks.length is 0
            head_nodes.push(j.nodes[cn])
         if j.nodes[cn].outlinks.length is 0
            tail_nodes.push(j.nodes[cn])            

      if strict and (head_nodes.length isnt 1 or tail_nodes.length isnt 1)
         throw Error("find_ends: multiple head/tail nodes found in scaffolded CC.")
         
      [head_nodes, tail_nodes]

###
## Returns [head, tail] lists from find_ends() above for all connected components
## calculates temporary maximal_spanning_tree() to produce this information, but
## restores the global state when done
###
find_all_ends = (j) ->
   
   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps
   
   try
      build_graph_refs(j) 
   catch err
      callback?(err, null)
      return

   output_list = (find_ends(j,c) for c in j.connected_comps)
   
   remove_graph_refs(j)
   
   output_list

###
# heuristically remove edges in a linear ccomp that span coverage and/or GC 
# discontinuities
###
check_connections = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Break linear scaffolds at a GC / coverage discontinuity
")
      if args.detailed_help?
         console.warn("
\n
        NOTE: Each run of #{args.help} will break a given scaffold at no more than one\n
        position. #{args.help} should be run multiple times if multiple misassemblies\n
        are suspected.\n
\n
Parameters:\n
\n
thresh : <float> -- Threshold used to determine whether or not to break a connection\n
\n
        Example: #{args.help} '{\"thresh\":0.5}' -- Default. Medium strength heuristic\n 
        score based on GC% / Coverage statistics of a scaffold on either side of a given\n
        contigs connection edge. The lower the threshold, the more sensitive the filter\n
        is to such discontinuities.\n
")
      callback?(null, j)
      return

   args.thresh ?= 0.5

   split = true
   
   while split

      split = false
      
      calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

      j.removed_edges ?= []   # initialize a place to keep removed edges in the output

      edges_kept = []

      try
         build_graph_refs(j) 
      catch err
         callback?(err, null)
         return

      for c in j.connected_comps
   
         # Find the "head" node
         try
            [[head_node]] = find_ends(j, c, true)
         catch err
            remove_graph_refs(j)
            callback?(err, null)
            return
      
         l = head_node.outlinks[0]
         cnodes = [head_node]
         cedges = []
      
         # Make arrays of nodes and edges in-order
         while l
            cnodes.push(l[1])
            cedges.push(l[0])
            l = l[1].outlinks[0]
      
         # Make arrays of cumulative cov and gc statistics up and down the scaffold

         tot = 0.0      
         len_up = (tot += n.seq_len for n in cnodes)
         tot = 0.0
         gc_up = ((tot += (n.pct_gc ? 0)*n.seq_len)/len_up[i] for n,i in cnodes)
         tot = 0.0
         cov_up = ((tot += n.cov*n.seq_len)/len_up[i] for n,i in cnodes)

         cnodes.reverse()
         tot = 0.0
         len_down = (tot += n.seq_len for n in cnodes)
         tot = 0.0
         gc_down = ((tot += (n.pct_gc ? 0)*n.seq_len)/len_down[i] for n,i in cnodes)
         tot = 0.0
         cov_down = ((tot += n.cov*n.seq_len)/len_down[i] for n,i in cnodes)

         gc_down.reverse()
         cov_down.reverse()
      
         edge_heuristics = (heuristic(gc_up[i],gc_down[i+1],cov_up[i],cov_down[i+1]) for i in [0...cedges.length])

         max_edge_idx = -1
         max_heur = -1
      
         # Find the maximum edge heuristic value and position
         for eh, i in edge_heuristics when eh > max_heur
            max_heur = eh
            max_edge_idx = i
        
         # If this maximum exceeds the given thresh, then remove this edge
         if max_heur > args.thresh
            split = true
            j.removed_edges.push(cedges.splice(max_edge_idx,1)[0])

         edges_kept = edges_kept.concat(cedges)
           
      j.edges = edges_kept

      remove_graph_refs(j)

      delete j.connected_comps   # This invalidates the cached ccomps
   
   callback?(null, j)

###
## filter edges from the graph by bitscore
###
filter_edges = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Remove all edges scoring less than thresh bits
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
thresh : <float> -- Bitscore threshold\n
\n
        Example: #{args.help} {\"thresh\":500.0} -- Default. Remove all edges scoring\n
        less than 500.0 bits.\n
")
      callback?(null, j)
      return

   args.thresh ?= 500.0

   delete j.connected_comps   # This invalidates the cached ccomps

   j.removed_edges ?= []   # initialize a place to keep removed edges in the output

   edges_kept = []

   # Walk all edges
   for e in j.edges 
      # Remove those that don't meet the threshold
      if e.bits < args.thresh
         j.removed_edges.push(e)
      else 
         edges_kept.push(e)
            
   j.edges = edges_kept
   
   callback?(null, j)
   
###
## Output statistics about current graph
###
graph_stats = (j, args = {}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Generate statistics about the assembly graph
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.txt[.gz]\" -- name of txt format file to write statistics to\n
\n
        Example: #{args.help} {\"file\":\"my_graph.txt\"} -- Write stats to the file\n
        my_graph.txt\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Default. Write stats to STDOUT\n
\n
ccdetail : true -- Write connected component details\n
\n
        Example: #{args.help} {\"ccdetail\":true} -- Equivalent to GCC\n
")
      callback?(null, j)
      return

   o = open_output_stream(args.file, args.tag)
   
   # Handle error opening file
   unless o
      callback?(new Error("open_output_stream could not open file '#{args.file}' for output."), null)
      return
   else
      o.on("error", (err)->
         console.error("ERROR: open_output_stream could not write to file '#{args.file}'.")
         callback?(err, null))

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps
   
   o.write("Nodes: #{Object.keys(j.nodes).length}\n")
   o.write("Edges: #{j.edges.length}\n") 
   o.write("Internal Edges: #{j.internal_edges.length}\n")
   o.write("Shared Sequence Edges: #{j.shared_seq_edges.length}\n")
   o.write("Removed Nodes: #{Object.keys(j.removed_nodes).length}\n") if j.removed_nodes?
   o.write("Removed Edges: #{j.removed_edges.length}\n") if j.removed_edges?
   
   [n50, seq_tot, cov, gc] = calc_seq_stats(j.nodes)
   
   o.write("Total sequence length: #{seq_tot}  N50: #{n50}\n")
   o.write("Mean coverage: #{cov.toFixed(1)}  GC content: #{gc.toFixed(1)}%\n")
   
   if j.clusters? and j.scaffolds?
      o.write("\nScaffold Clusters: #{j.clusters.length}\n")
      if args.ccdetail?
         o.write("\nScaffold Cluster details: \nclust\tscaffs\tnodes\tseqlen\tcov\t%GC\tn50\tlongest contig\n")
         for cl, ii in j.clusters
            n_ids = []
            for cc in cl
               n_ids = n_ids.concat(j.scaffolds[cc].nodes)
            [n50, seq_tot, cov, gc, max_node] = calc_seq_stats(j.nodes[nid] for nid in n_ids)
            o.write("#{ii}\t#{cl.length}\t#{n_ids.length}\t#{seq_tot}\t#{cov.toFixed(1)}\t#{gc.toFixed(1)}%\t#{n50}\t#{max_node}\n")
            o.write("\n\tScaffold details:\n\tccnum\tnodes\tseqlen\tcov\t%GC\tn50\tlongest contig\tname\n")
            
            for cc in cl
               [n50, seq_tot, cov, gc, max_node] = calc_seq_stats(j.nodes[nid] for nid in j.scaffolds[cc].nodes)   
               o.write("\t#{j.scaffolds[cc].ccnum}\t#{j.scaffolds[cc].nodes.length}\t#{seq_tot}\t#{cov.toFixed(1)}\t#{gc.toFixed(1)}%\t#{n50}\t#{max_node}\t#{cc}\n")
            o.write("\n")   
   
   if j.connected_comps?
      o.write("\nConnected components: #{j.connected_comps.length}\n")
      if args.ccdetail?
         o.write("\nConnected component details:\nccnum\tnodes\tseqlen\tcov\t%GC\tn50\tlongest contig")
         if j.scaffolds? 
            o.write("\tscaffold name\n")
            scafnames = []
            for sn, s of j.scaffolds 
               scafnames[s.ccnum] = sn
         else 
            o.write("\n")

         for c, i in j.connected_comps
            [n50, seq_tot, cov, gc, max_node] = calc_seq_stats(j.nodes[nid] for nid in c)         
            o.write("#{i}\t#{c.length}\t#{seq_tot}\t#{cov.toFixed(1)}\t#{gc.toFixed(1)}%\t#{n50}\t#{max_node}")
            if scafnames?[i]?
               o.write("\t#{scafnames[i]}\n")
            else   
               o.write("\n")
   
   o.write("\nProcessing steps completed (command history for this graph):\n")
   for p in j.processing?[..-2]
      o.write("#{p[0]}\t#{p[1]}\t#{p[2]}\t#{JSON.stringify(p[3]) if p[3]}\n")

   o.write("\nStash contains #{stash_stack.length} saved graph#{if stash_stack.length isnt 1 then 's' else ''}.\n")
   if args.ccdetail? and stash_stack.length
      o.write("\nCommand histories for all stashed graphs (top to bottom of stack)\n")
      stash_stack.reverse()
      for s, i in stash_stack
         o.write("\nStack position #{stash_stack.length-i} #{if i is 0 then '[top]' else ''}#{if i is stash_stack.length-1 then '[bottom]' else ''}:\n\n")
         for p in s.processing?[..-2]
            o.write("#{p[0]}\t#{p[1]}\t#{p[2]}\t#{JSON.stringify(p[3]) if p[3]}\n")
      stash_stack.reverse() 
   o.end((error) ->
         if error
            callback?(error, null)
         else
            callback?(null, j))
   
###
## Print graph statistics with CC detail and other details
## Shortcut equivalent to:  GC {"ccdetail":true}
###
graph_stats_cc = (j, args ={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Generate statistics about the assembly graph, with details about\n
       each connected component
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.txt[.gz]\" -- name of txt format file to write statistics to\n
\n
        Example: #{args.help} {\"file\":\"my_graph_ccs.txt\"} -- Write stats to the file\n
        my_graph_ccs.txt\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Default. Write stats to STDOUT\n
")
      callback?(null, j)
      return

   args.ccdetail = true
   graph_stats(j, args, callback)

###
# Write DOT output file from active nodes and edges.
#
# There is a fair amount of default "style" applied here with a few options.
#  
# To change styles (or provide alternatives), don't add a lot of code here,
# rather you should look into the "gvpr" language that comes with the graphviz 
# tools... It is purpose built for re-styling DOT files, etc before rendering.
###
export_dot = (j, args = {}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Write the current assembly graph to the graphviz DOT format
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.dot[.gz]\" -- name of dot format file to write graph to\n
\n
        Example: #{args.help} {\"file\":\"my_graph.dot\"} -- Write stats to the file\n
        my_graph.dot\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Default. Write graph to STDOUT\n
\n
detail : true -- Draw labelled contig nodes\n
\n
        Example: #{args.help} {\"detail\":true} -- Draw contigs as ovals containing\n
        the node names of each contig. Note that this disables scaling the node size\n
        by contig length.\n
\n
arrowtype : \"arrowtype_string\" -- Control the type of arrowheads drawn.\n
\n
        Example: #{args.help} {\"arrowtype\":\"normal\"} -- Default. Draw normal arrows\n
        See graphviz documentation at:\n
        http://www.graphviz.org/doc/info/attrs.html#k:arrowType\n
\n
pen_scale : <float> -- Control the relative thickness of the arrow lines.\n
\n
        Example: #{args.help} {\"pen_scale\":0.1} -- Default. Draw normal arrows\n
        See graphviz documentation at:\n
        http://www.graphviz.org/doc/info/attrs.html#d:penwidth\n
\n
const_edge : true -- Draw connection arrow lines at constant width\n
\n
        Example: #{args.help} {\"const_edge\":true} -- Draw edge arrows of constant width\n
        regardless of bitscore.\n
\n
colored_edges : true -- Draw edges colored by the GC% of the connected contigs\n
\n
        Example: #{args.help} {\"colored_edges\":true} -- Draw colored arrows\n
")
      callback?(null, j)
      return

   o = open_output_stream(args.file, args.tag)

   # Handle error opening file
   unless o
      callback?(new Error("open_output_stream could not open file '#{args.file}' for output."), null)
      return
   else
      o.on("error", (err)->
          console.error("ERROR: open_output_stream could not write to file '#{args.file}'.")
          callback?(err, null))
   try
      build_graph_refs(j) 
   catch err
      callback?(err, null) 
      return
   
   args.arrowtype ?= "normal"
   args.pen_scale ?= 0.1

   doublecircle = "shape = \"doublecircle\" style=\"filled\" color = \"red\" penwidth = 3"

   if j.SEASTAR_version?
      o.write "di" if j.digraph
      o.write("graph \"#{args.file ? 'Unnamed'}\" {\n")
      o.write("graph 
[K = 0.5,
 repulsiveforce = 1.5,
 overlap = \"prism\",
 overlap_scaling = 10000.0
];\n")

      unless args.detail
         o.write("node [shape=\"point\"];\n")
      else    
         o.write("node [style=\"filled\"];\n")
         
      o.write("edge [arrowhead=\"#{args.arrowtype}\"];\n")
         
      for name, node of j.nodes
          h = ((node.pct_gc ? 0) - 30.0) / 30.0
          h = 0.0 if h < 0.0
          h = 1.0 if h > 1.0
          o.write("\"#{name}\" 
[bits = #{(node.bits ? 0).toFixed(4)},
 rel_ab = #{(node.rel_ab ? 0).toFixed(15)},
 cov = #{(node.cov ? 0).toFixed(4)},
 cov2 = #{(node.adj_cov ? 0).toFixed(4)},
 uncov = #{((node.pct_uncov ? 0)/100.0).toFixed(4)},
 seq_len = #{(node.seq_len ? 0).toFixed(0)},
 rd_len = #{(node.rd_len ? 0).toFixed(1)},
 gc = #{((node.pct_gc ? 0)/100.0).toFixed(4)},
 n_shr = #{(node.mp_sh ? 0).toFixed(0)},
 n_bwd = #{(node.mp_bwd ? 0).toFixed(0)},
 n_fwd = #{(node.mp_fwd ? 0).toFixed(0)},
 mp_mean = #{(node.mp_ins_mean ? 0).toFixed(4)},
 mp_stdev = #{(node.mp_ins_stdev ? 0).toFixed(4)},
 mp_sl_pairs = #{(node.mp_pairs ? 0).toFixed(0)},
 name = \"#{name}\",
 width = #{0.01*Math.sqrt(node.seq_len)}
 fillcolor = \"#{h},1.0,0.8\"
 color = \"#{h},1.0,0.5\"
 label = \"#{(if ((args.detail is 2) and (node.desc?)) then node.desc else name)}\"
 #{if node.contig_problems? then doublecircle else ""}
];\n") 

      for edge in j.edges
          edge_type = if edge.dir is "pos" or edge.dir is "forward" then "->" else "--"
          if edge.dir is "pos" or edge.interscaffold
             o.write("\"#{edge.n1}\" #{edge_type} \"#{edge.n2}\"
 [type = DEP,
 dir = #{edge.dir},
 penwidth = #{(args.pen_scale*Math.sqrt((if ((not args.const_edge) and (edge.score > 0)) then edge.score else 5)))},
 color = red
];\n")          
          else unless args.colored_edges
             o.write("\"#{edge.n1}\" #{edge_type} \"#{edge.n2}\"
 [type = MP,
 dir = #{edge.dir},
 bits = #{edge.bits.toFixed(3)},
 score = #{edge.score?.toFixed(3)},
 pos1 = #{edge.p1},
 pos2 = #{edge.p2},
 penwidth = #{(args.pen_scale*Math.sqrt((if ((not args.const_edge) and (edge.score > 0)) then edge.score else 1)))},
 arrowsize = #{0.1*(Math.log((if ((not args.const_edge) and (edge.score > 0)) then edge.score else 1), 2))},
 color = black
];\n") 
          else
             h = ((((edge.src.pct_gc ? 0) + (edge.tar.pct_gc ? 0)) / 2) - 30.0) / 30.0
             h = 0.0 if h < 0.0
             h = 1.0 if h > 1.0
             o.write("\"#{edge.n1}\" #{edge_type} \"#{edge.n2}\"
 [type = MP,
 dir = #{edge.dir},
 bits = #{edge.bits.toFixed(3)},
 score = #{edge.score?.toFixed(3)},
 pos1 = #{edge.p1},
 pos2 = #{edge.p2},
 penwidth = #{(args.pen_scale*Math.sqrt((if ((not args.const_edge) and (edge.score > 0)) then edge.score else 1)))},
 arrowsize = #{0.1*(Math.log((if ((not args.const_edge) and (edge.score > 0)) then edge.score else 1), 2))},
 color = \"#{h},1.0,0.5\"
];\n")  

         o.write("}\n")

   remove_graph_refs(j)

   o.end((error) ->
         if error
            callback?(error, null)
         else
            callback?(null, j))

###
# FASTA file output
###
export_fasta = (j, args={"verbose":false}, callback) -> 
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Write sequences contained in the current graph data to a FASTA format file.
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.fasta[.gz]\" -- name of FASTA format file to write sequence to\n
\n
        Example: #{args.help} {\"file\":\"my_seq.fasta\"} -- Write stats to the file\n
        my_seq.fasta\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Default. Write sequence to STDOUT\n
\n
scaff : true -- Output fully scaffolded sequences (using seq_scaffold tool)\n
\n
        Example: #{args.help} {\"scaff\":true} -- Output scaffolded contig sequences\n
\n
        NOTE: Using this option as in the above example will run seq_scaffold with its\n
        default settings.  As an advanced option, this parameter can also accept a string\n
        argument, which will be passed along to the external seq_scaffold tool as its\n
        [options] parameter string. Run with --help for help with the settings offered by\n
        seq_scaffold and the default values.\n
\n
        Example: #{args.help} {\"scaff\":\"--help\"}\n
\n
        Example: #{args.help} {\"scaff\":\"--overlap=7 --heal=othercontigs.fna\"}\n
\n
no_merge_scaffs : true -- Write contig sequences in scaffold order with a scaffold ID in\n
        each header, ready to be provided to the seq_scaffold tool\n
\n
        Example: #{args.help} {\"no_merge_scaffs\":true} -- Output ordered contig sequences\n
\n
abundance : true -- append relative abundance values to the FASTA sequence IDs\n
\n
        Example: #{args.help} {\"abundance\":true} -- Attach abundances\n
")
      callback?(null, j)
      return

   if args.no_merge_scaffs 
      args.scaff = true
      
   if args.scaff and not j.scaffolds?
      callback?(new Error("Can't write scaffolded sequence before the SCAFF command has been run."), null)
      return

   unless args.stream?
      o = open_output_stream(args.file, args.tag)
   else
      o = args.stream   # This is a special internal case, passing in a stream to write to
      
   # Handle error opening file
   unless o
      callback?(new Error("open_output_stream could not open file '#{args.file}' for output."), null)
      return
   else
      o.on("error", (err)->
          console.error("ERROR: open_output_stream could not write to file '#{args.file}'.")
          console.error(err)
          callback?(err, null))

   # Here, just write out the sequences, with optional relative abundance
   unless args.scaff
      for name, node of j.nodes
         if node.recon_seq?
            if args.abundance?
               o.write(">#{name}_#{(node.rel_ab*100).toFixed(6)}\n")
            else
               o.write(">#{name}\n")
            o.write("#{node.recon_seq}\n")
      o.end((err)->
         if err
            callback?(err, null)
         else
            callback?(null, j))

   # Otherwise, we're going to be writing out contigs in scaffold order, with ids to
   # hook them together downstream.
   else
      calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

      try
         build_graph_refs(j) 
      catch err
         callback?(err, null)
         return

      # This function writes the scaffolds to a stream
      write_scaffolds = (out) ->

         for scaf_name, scaf of j.scaffolds
            for nid in scaf.nodes 
               node = j.nodes[nid]
               if node.recon_seq?
                  out.write(">#{node.id} #{scaf_name}\n")
                  out.write("#{node.recon_seq}\n")
      
      # In this case the stream is a file
      if args.no_merge_scaffs
         write_scaffolds(o)
         o.end((err) ->
            if err
               callback?(err, null)
            else
               callback?(null, j))

      # In this case it is the stdin of a script that will actually do the
      # layout, test for overlaps between neighboring contigs, and adjust
      # the number of 'N's in gaps to avoid disrupting the longest ORF that
      # spans each gap.  E.g. seq_scaffold
      # Note, that this functionality, currently written in gawk, should be 
      # ported forward to coffeescript and more tightly integrated here.
      else
         # Decide whether to go with the default options or use the provided arg string
         if typeof(args.scaff) is 'string'
            options = args.scaff.split(/\s+/)
            options.push('-')
         else   
            options = ['-']
            
         scaf = child_process.spawn('seq_scaffold',options,{})
         
         scaf.stdout.setEncoding('ascii')      
         scaf.stderr.setEncoding('ascii')

         scaf.stdout.on('data', (data) ->
            o.write(data)
         ) 
    
         scaf.stderr.on('data', (data) ->
            process.stderr.write(data)
         ) 

         scaf.stdout.on('end', () ->
            o.end((error) ->
               if error
                  callback?(error, null)
               else
                  callback?(null, j)))

         scaf.on('exit', (code)->
            if (code) 
               console.error('ERROR: ' + code)
               callback?(new Error('SEAStAR script "seq_scaffold" could not be executed. Please check that this file is in your PATH.'), null)
         )
         
         write_scaffolds(scaf.stdin)
         scaf.stdin.end()   

      remove_graph_refs(j)
      
      # No callback here because it happens a few lines above in the 'end' event handler
   
###
# old ref_select style table output
#
# This table is used by the 16S abundance pipeline and for doing transcriptomics analysis
###
export_table = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Write ref_select statistics for selected reference sequences to a TSV file
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.tsv[.gz]\" -- Specify a file name to write the TSV format stats table to.\n
        If the filename contains one or more '#' characters in a row, these positions are\n
        replaced with zero-padded digits that will increment each time a file is written to\n
        this filename pattern. If no '#' characters are present, then this command overwrites\n
        any existing file of the same name. (See DUMP command an example using this behavior)\n
\n
        Example: #{args.help} {\"file\":\"my_seq.tsv\"} -- Write stats to the file\n
        my_seq.tsv\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Default. Write stats to STDOUT\n
\n
header : true -- write a header row with labels for each column\n
\n
        Example: #{args.help} {\"header\":true} -- Write a header row. The fields are:\n
\n
        bitscore - Information content of reads aligning with the ref sequence\n
        read_cnt - Number of (possibly fractional) reads aligning with the ref sequence\n
        norm_cnt - Read_cnt normalized to ref sequence length\n
        rel_abun - Relative fractional abundance of (copy number of) the ref sequence to\n
                     all selected sequences (those with bitscores above some thresh)\n
        mean_cov - Mean coverage of the reference sequence by (possibly fractional) reads\n
        read_len - Mean length of reads aligning with this sequence\n
        seq_len  - Length of the ref sequence\n
        pct_gc   - Percent GC content of ref seqeunce (NA if not calculated)\n
        name     - Catalog name of the ref sequence\n
        desc     - Catalog description of the ref sequence\n
")
      callback?(null, j)
      return

   o = open_output_stream(args.file, args.tag)
   
   # Handle error opening file
   unless o
      callback?(new Error("open_output_stream could not open file '#{args.file}' for output."), null)
      return
   else
      o.on("error", (err)->
          console.error("ERROR: open_output_stream could not write to file '#{args.file}'.")
          callback?(err, null))
   
   if args.header
      o.write("id\t
bitscore\t
read_cnt\t
norm_cnt\t
rel_abun\t
mean_cov\t
read_len\t
seq_len\t
pct_gc\t
name\t
desc\n
")
   
   for name, node of j.nodes
      o.write("#{name}\t
#{node.bits}\t
#{node.rd_cnt}\t
#{node.int_cov}\t
#{node.rel_ab}\t
#{node.cov}\t
#{node.rd_len}\t
#{node.seq_len}\t
#{node.pct_gc ? 'NA'}\t
#{node.name}\t
#{node.desc}\n
")

   o.end((error) ->
         if error
            callback?(error, null)
         else
            callback?(null, j))

###
# Reduces the graph to nodes in the object containing node_ids
###
execute_selection = (j, nodes) ->

   edges_kept = []
   removed_edges_kept = []
   removed_nodes_kept = {}
   internal_edges_kept = []
   shared_edges_kept = []

   try
      build_removed_graph_refs(j)
   catch err
      callback?(err, null)
      return

   for id, n of nodes
      for e in n.inlinks
         if nodes[e[1].id]?
            edges_kept.push(e[0])
            
      for e in n.links
         unless nodes[e[1].id]?
            removed_edges_kept.push(e[0])
            removed_nodes_kept[e[1].id] = e[1]
            
      for e in n.rem_inlinks
         if nodes[e[1].id]?
            removed_edges_kept.push(e[0])
         
      for e in n.rem_links   
         unless nodes[e[1].id]?
            removed_edges_kept.push(e[0])
            removed_nodes_kept[e[1].id] = e[1]

   for id, n of removed_nodes_kept
      for e in n.rem_inlinks
         if removed_nodes_kept[e[1].id]?
            removed_edges_kept.push(e[0])

   j.nodes = nodes
   j.edges = edges_kept
   j.removed_edges = removed_edges_kept
   j.removed_nodes = removed_nodes_kept
    
   j.shared_seq_edges = (e for e in j.shared_seq_edges when nodes[e.n1]? or nodes[e.n2]? or removed_nodes_kept[e.n1]? or removed_nodes_kept[e.n2]?)
   j.internal_edges = (e for e in j.internal_edges when nodes[e.n1]? or removed_nodes_kept[e.n1]?)
   
   j

###
# Select ccomps for output or further processing
#
# args (mutually exclusive):
#     ccname = string id of a node to find in a cc
#     ccnames = list of strings of names to find in multiple
#     ccnum = number of a cc to select
#     ccrange = start and end ccnums to select a range e.g. [0,10] or [1,-1]
#     ccnums = list of numbers of ccs to select
#     shift = select all but the first ccomp, like [1,-1] 
#     min_nodes = minimum number of nodes in cc
#     min_seqlen = minimum amount of sequence in nodes in cc
#     sequence = sequence to search for 
#
###
grab_ccomps = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Select specific connected components for further processing
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
ccname : \"contig_name\" -- Select the connected component containing the named contig\n
\n
        Example: #{args.help} {\"ccname\":\"NODE_1234\"} -- Select the connected component\n
        containing the contig named NODE_1234
\n
ccnames : [\"contig_name1\",\"contig_name2\",...] -- Select the connected component(s)\n
        containing the named contigs\n
\n
        Example: #{args.help} {\"ccnames\":[\"NODE_1234\",\"NODE_5678\"]} -- Select the\n
        connected components containing the contigs named NODE_1234 and NODE_5678\n
\n
ccnum : <int> -- Select connected component number <int>\n
\n
        Example: #{args.help} {\"ccnum\":0} -- Default. Select connected component 0.\n
\n
ccnums : [<int>,<int>,...] -- Select the connected components from the list of numbers\n
\n
        Example: #{args.help} {\"ccnums\":[1,2]} -- Select the second and third\n
        connected components (numbering is zero based)\n
\n
ccrange : [<int1>,<int2>] -- Select the connected components numbered in the range\n
        <int1>..<int2> (inclusive). See also: 'shift' option below.\n
\n
        NOTE: <int> may be negative, indicating positions at the end of the list of\n
        connected components.\n 
\n
        Example: #{args.help} {\"ccrange\":[0,5]} -- Select the first 6 connected\n
        components\n
\n
        Example: #{args.help} {\"ccrange\":[-5,-1]} -- Select the last 5 connected\n
        components\n
\n
shift : true -- Select all connected components except the first one.\n
\n
        Example: #{args.help} {\"shift\":true} -- Drop the first connected component.\n
\n
        This is like #{args.help} {\"ccrange\":[1,-1]} except it doesn't generate a\n
        fatal error when there is only one remaining connected component, allowing\n
        processing to potentially continue in any calling SCRIPT commands.\n
\n
NOTE: The following parameters are \"filters\" that are applied to the set of connected\n
components selected by one of the above parameters (or by default, the set of all\n
connected components.) These filters may be used in combination, resulting in a logical\n
\"AND\" relationship (only connected components satisfying all of the filters are\n
selected).\n
\n
min_nodes : <int> -- Select connected components with <int> or more nodes.\n
\n
        Example: #{args.help} {\"min_nodes\":2} -- Select connected components\n
        with 2 or more nodes.\n
\n
min_seqlen : <int> -- Select connected components with <int> or more sequence within nodes.\n
\n
        Example: #{args.help} {\"min_seqlen\":1000} -- Select connected components containing\n
        at least 1000 bases of sequence.\n
\n
sequence : <string> -- Select connected components containing the provided DNA sequence.\n
\n
        NOTE! This isn't BLAST, the sequence must match exactly. Any differences, including\n
        ambiguity codes, etc. will prevent matching. The only extra thing that is done is the\n
        reverse complement of the provided sequence is also searched.\n
\n
        Example: #{args.help} {\"sequence\":\"AGACTAGCAGATATACGATAACGATACGATACGAT\"}\n
        Select connected components containing the provided sequence (or its reverse\n
        complement).\n
\n
sequences : [<string>, ...] -- Like 'sequence' parameter, but takes a list of sequences.\n
\n
        Example: #{args.help} {\"sequences\":[\"AGACTAGCAGATATAC\",\"GATAACGATACGATACGAT\"]}\n
        Select connected components containing any of the provided sequences (or their\n
        reverse complements).\n
")
      callback?(null, j)
      return

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

   try
      build_graph_refs(j) 
   catch err
      callback?(err, null)
      return

   if args.shift
      args.ccrange=[1,-1]
   
   if args.ccrange? # Select a range of CComps by number
      args.ccrange[0] = (j.connected_comps.length + args.ccrange[0]) if args.ccrange[0] < 0
      args.ccrange[1] = (j.connected_comps.length + args.ccrange[1]) if args.ccrange[1] < 0
      args.ccnums = [args.ccrange[0]..args.ccrange[1]]
   else if args.ccname?  # Find the CComp containing node with this name
      for c, i in j.connected_comps when c.indexOf(args.ccname) != -1
         args.ccnums = [i]
         break
   else if args.ccnames? # Find all CComps with nodes with these names 
      args.ccnums = []
      for id in args.ccnames
         for c, i in j.connected_comps when c.indexOf(id) isnt -1 and args.ccnums.indexOf(i) is -1
            args.ccnums.push(i)
            break
   else if args.ccnum? # Select CComp with this number (ordered by number of nodes)
      args.ccnums = [args.ccnum]

   if args.sequence
      args.sequences = [ args.sequence ]
      
   if args.sequences
      ccnums = args.ccnums ? [0...j.connected_comps.length]
      args.ccnums = []
      all_seqs = args.sequences.map(reverse_comp).concat(args.sequences)
      for i in ccnums
         c = j.connected_comps[i]
         for nid in c when (seq = j.nodes[nid].recon_seq) and not (i in args.ccnums)
            for s in all_seqs when seq.indexOf(s) isnt -1
               args.ccnums.push(i)
               break

   if args.min_nodes? # Select all CComps with at least this many nodes
      ccnums = args.ccnums ? [0...j.connected_comps.length]
      args.ccnums = []
      for i in ccnums when j.connected_comps[i].length >= args.min_nodes
         args.ccnums.push(i)

   if args.min_seqlen? # Select all CComps with at least this much sequence with the nodes
      ccnums = args.ccnums ? [0...j.connected_comps.length]
      args.ccnums = []
      for i in ccnums when cc_seq_len(j,j.connected_comps[i]) >= args.min_seqlen
         args.ccnums.push(i)

   unless args.ccnums?  # Default to ccnum 0
      args.ccnums = [0]
   
   unless args.ccnums.length  # Default parameter type is list of ccnums 
      remove_graph_refs(j)
      console.warn "WARN: SELCC: No connected components found matching selection criteria"
      callback?(undefined, null)
      return
   
   nodes = {}
   new_cclist = []
   for cc in args.ccnums
      unless j.connected_comps[cc]?
         remove_graph_refs(j)
         if args.shift
            callback?(null, undefined)
         else
            callback?(new Error("Invalid selected ccnum: #{cc}.  ccnums are zero-based and there are only #{j.connected_comps.length} ccomps in this graph."), null)
         return
      new_cclist.push(j.connected_comps[cc])
      for nid in j.connected_comps[cc]
         nodes[nid] = j.nodes[nid]
         
   execute_selection(j, nodes)
   
   j.connected_comps = new_cclist  # Set new connected components
   
   if j.scaffolds?  # Maintain the current mapping of scaffolds to ccomps 
      j.connected_comps = new_cclist  # Set new connected components
      new_scaffs = {}                 # Build a new scaffold object
      for sn, sc of j.scaffolds       # Walk the current scaffold object and rematch scaffs to ccomps
         nid = sc.nodes[0]            # Choose a node in this scaffold
         # Find a ccomp containing this node
         for cc, ci in j.connected_comps when cc.indexOf(nid) isnt -1  
            sc.ccnum = ci             # The ccnumber has changed, so correct it 
            new_scaffs[sn] = sc       # Copy in the scaffold (keeping the same name)
            break                     # Done processign this scaffold
      j.scaffolds = new_scaffs        # Commit the new scaffold list
      delete j.clusters               # Clusters don't survive this.
   
   remove_graph_refs(j) 
   
   callback?(null, j)

###
# Select nodes and neighbors for output or further processing
###

find_neighbors = (nd, r) ->
   nodeid_list = [nd.id]
   nd.seen = true
   if r > 0
      for [l,n] in nd.links when not n.seen
         nodeid_list = nodeid_list.concat(find_neighbors(n,r-1))
   nodeid_list

###
#   Adds neighboring nodes to the selected list
#
# args:
#     name = string id of a node to find in a cc
#     names = list of strings of names to find in multiple
#     radius = number of hops to take searching for neighoring nodes
#     sequence = DNA sequence to search for.
###
grab_neighbors = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Select contigs from the connected neighborhood(s) of the given contig(s)
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
name : \"contig_name\" -- Use a single contig by name\n
\n
        Example: #{args.help} {\"name\":\"NODE_1234\"} -- Select NODE_1234 (and its\n
        neighbors)\n
\n
names : [\"contig_name1\",\"contig_name2\",...] -- Use multiple contigs by name\n
\n
        Example: #{args.help} {\"names\":[\"NODE_1234\",\"NODE_5678\"]} -- Select these\n
        two contigs (and their neighbors...)\n
\n
radius : <int> -- Size of the neighborhood of contigs to select\n
\n
        Example: #{args.help} {\"radius\":2} -- Select all neighbors and neighbors of\n        
        neighbors\n
\n
        Example: #{args.help} {\"radius\":0} -- Default. Select only the named contig(s)\n
\n
sequence : <string> -- Select contig(s) found to contain the provided DNA sequence.\n
\n
        NOTE! This isn't BLAST, the sequence must match exactly. Any differences, including\n
        ambiguity codes, etc. will prevent matching. The only extra thing that is done is the\n
        reverse complement of the provided sequence is also searched.\n
\n
        Example: #{args.help} {\"sequence\":\"AGACTAGCAGATATACGATAACGATACGATACGAT\"}\n
        Select contig(s) containing the provided sequence (or its reverse complement).\n        
\n
sequences : [<string>, ...] -- Like 'sequence' parameter, but takes a list of sequences.\n
\n
        Example: #{args.help} {\"sequences\":[\"AGACTAGCAGATATAC\",\"GATAACGATACGATACGAT\"]}\n
        Select contigs containing any of the provided sequences (or their reverse\n
        complements).\n
")
      callback?(null, j)
      return

   try
      build_graph_refs(j) 
   catch err
      callback?(err, null)
      return

   args.names = [args.name] if args.name? 

   if args.sequence
      args.sequences = [ args.sequence ]
   if args.sequences
      args.names = []
      all_seqs = args.sequences.map(reverse_comp).concat(args.sequences)
      for nid, n of j.nodes when seq = n.recon_seq 
         for s in all_seqs when seq.indexOf(s) isnt -1
            args.names.push(nid)
            break

   if args.sequence
      args.names = []
      revseq = reverse_comp(args.sequence)
      for nid, n of j.nodes when seq = n.recon_seq
         if seq.indexOf(args.sequence) isnt -1 or seq.indexOf(revseq) isnt -1
            args.names.push(nid)

   args.radius ?= 0

   nodes = {}
   for name in args.names
      for nid in find_neighbors(j.nodes[name], args.radius) 
         nodes[nid] = j.nodes[nid]
        
   execute_selection(j, nodes)
   
   remove_graph_refs(j, (n)->delete n.seen) 

   delete j.connected_comps 
   
   callback?(null, j)

###
# Functions for finding connections between scaffold ends 
###
find_direct_connection = (head, tail, thresh = 0) ->

   try_dir_links = (h,t) -> 
      e = null
      if h and t
         for l in h.rem_links
            if l[1].id is t.id and (e = edge_director(t, h, l[0])) and t.id is e.n1 
               break
            else 
               e = null
      if e
         [e]
      else 
         null

   select_edge = null
   select_head = null
   select_tail = null

   if (edge = try_dir_links(head, tail)) and edge[0]?.bits > thresh
      select_head = head
      select_tail = tail
      select_edge = edge
      thresh = edge[0].bits
      
   if (edge = try_dir_links(head.outlinks[0]?[1], tail)) and edge[0]?.bits > thresh
      select_head = head.outlinks[0][1]
      select_tail = tail
      select_edge = edge
      thresh = edge[0].bits
      
   if (edge = try_dir_links(head, tail.inlinks[0]?[1])) and edge[0]?.bits > thresh
      select_head = head
      select_tail = tail.inlinks[0][1]
      select_edge = edge
      thresh = edge[0].bits
   
   if (edge = try_dir_links(head.outlinks[0]?[1], tail.inlinks[0]?[1])) and edge[0]?.bits > thresh
      select_head = head.outlinks[0][1]         
      select_tail = tail.inlinks[0][1]
      select_edge = edge
      thresh = edge[0].bits
         
   [select_edge, select_head, select_tail, thresh]

###
# Look for indirect connections (via another node not part of either scaffold)
###
find_indirect_connection = (head, tail, thresh = 0) ->
  
   try_indir_links = (h,t) -> 
      et = null
      eh = null
      if h and t
         hnodes = {}
         for l in h.rem_links
            hnodes[l[1].id] = l

         for l in t.rem_links when hnodes[l[1].id]?
            if (et = edge_director(t, l[1], l[0])) and (t.id is et.n1) and
               (eh = edge_director(hnodes[l[1].id][1], h, hnodes[l[1].id][0])) and (h.id is eh.n2)
                  break
            else 
                  et = null
                  eh = null
               
         if (et and eh) then [et, eh] else null
   
   select_result = null
   select_head = null
   select_tail = null
   
   if (result = try_indir_links(head.outlinks[0]?[1], tail.inlinks[0]?[1])) and result?[0]?.bits > thresh and result?[1]?.bits > thresh
      select_head = head.outlinks[0][1]         
      select_tail = tail.inlinks[0][1]
      select_result = result
      thresh = Math.min(result[0].bits, result[1].bits)
   
   if (result = try_indir_links(head, tail.inlinks[0]?[1])) and result?[0]?.bits > thresh and result?[1]?.bits > thresh
      select_head = head
      select_tail = tail.inlinks[0][1]
      select_result = result
      thresh = Math.min(result[0].bits, result[1].bits)
   
   if (result = try_indir_links(head.outlinks[0]?[1], tail)) and result?[0]?.bits > thresh and result?[1]?.bits > thresh
      select_head = head.outlinks[0][1]         
      select_tail = tail
      select_result = result
      thresh = Math.min(result[0].bits, result[1].bits)
   
   if (result = try_indir_links(head, tail)) and result?[0]?.bits > thresh and result?[1]?.bits > thresh
      select_head = head
      select_tail = tail
      select_result = result
      thresh = Math.min(result[0].bits, result[1].bits)

   [select_result, select_head, select_tail, thresh]

###
# Look for both direct and indirect connections
###
find_connection = (head, tail, thresh = 0) ->
   direct = find_direct_connection(head, tail, thresh)
   
   thresh = direct[3] if direct[0]
   
   indirect = find_indirect_connection(head, tail, thresh)

   if indirect[0]
      indirect
   else
      direct

###
# Scaffold end-link
###
scaff_link = (j, args={}, callback) -> 
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Find connections linking scaffold ends
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
thresh : <float> -- Bitscore threshold\n
\n
        Example: #{args.help} {\"thresh\":1000.0}' -- Default. Only reconnect scaffold\n
        ends with connections scoring at least 1000.0 bits\n
")
      callback?(null, j)
      return
   
   args.thresh ?= 1000.0
   
   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps  
   
   try
      build_graph_refs(j) 
      build_removed_graph_refs(j)
   catch err
      callback?(err, null)
      return
   
   # Add ccnums to nodes for later reference
   for c, i in j.connected_comps
      for n in c
         j.nodes[n].ccnum = i

   # Try all pairs of scaffolds
   for c, i in j.connected_comps
      for cc, ii in j.connected_comps when ii >= i
   
         # Find the "head" and "tail" nodes for these scaffold
         
         try
            [[h1], [t1]] = find_ends(j, c, true)
            [[h2], [t2]] = find_ends(j, cc, true)
         catch err
            remove_graph_refs(j, ((n) -> delete n.ccnum)) 
            callback?(err, null)
            return
            
         [e1] = find_connection(h1, t2, args.thresh)

         if e1         
            if e1.length is 1   # Direct connection (one edge)
               unless e1[0].interscaffold
                  e1[0].interscaffold = true
                  j.edges.push(e1[0])
                  delete j.removed_edges[e1[0].index]
                  
               console.warn("Scaffold #{ii} #{e1[0].n1} --> Scaffold #{i} #{e1[0].n2} -- #{e1[0].bits} bits")
               
            else # Indirect connection (two edges and an inbetween node)
               unless e1[0].interscaffold
                  e1[0].interscaffold = true
                  j.edges.push(e1[0])
                  delete j.removed_edges[e1[0].index]
               unless e1[1].interscaffold   
                  e1[1].interscaffold = true
                  j.edges.push(e1[1])
                  delete j.removed_edges[e1[1].index]
                  
               unless j.nodes[e1[1].n1]?   
                  j.nodes[e1[1].n1] = e1[1].src
                  delete j.removed_nodes[e1[1].n1]
                  
               console.warn("Scaffold #{ii} #{e1[0].n1} --> Scaffold #{e1[1].src.ccnum} #{e1[1].n1} --> Scaffold #{i} #{e1[1].n2} -- #{e1[0].bits} + #{e1[1].bits} bits")

         unless i == ii
            [e2] = find_connection(h2, t1, args.thresh)

            if e2
               if e2.length is 1   # Direct connection (one edge)
                  unless e2[0].interscaffold
                     e2[0].interscaffold = true
                     j.edges.push(e2[0])
                     delete j.removed_edges[e2[0].index]
                     
                  console.warn("Scaffold #{i} #{e2[0].n1} --> Scaffold #{ii} #{e2[0].n2} -- #{e2[0].bits} bits")
               else
                  unless e2[0].interscaffold
                     e2[0].interscaffold = true
                     j.edges.push(e2[0])
                     delete j.removed_edges[e2[0].index]
                  unless e2[1].interscaffold   
                     e2[1].interscaffold = true
                     j.edges.push(e2[1])
                     delete j.removed_edges[e2[1].index]

                  unless j.nodes[e2[1].n1]?   
                     j.nodes[e2[1].n1] = e2[1].src
                     delete j.removed_nodes[e2[1].n1]

                  console.warn("Scaffold #{i} #{e2[0].n1} --> Scaffold #{e2[1].src.ccnum} #{e2[1].n1} --> Scaffold #{ii} #{e2[1].n2} -- #{e2[0].bits} + #{e2[1].bits} bits")

   remove_graph_refs(j, ((n) -> delete n.ccnum)) 
   
   callback?(null, j)

###
# Relink
###
relink = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Reconnect previously removed connections between contigs
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
name : \"contig_name\" -- Use a single contig by name\n
\n
        Example: #{args.help} {\"name\":\"NODE_1234\"} -- Restore edges to this contig\n
\n
names : [\"contig_name1\",\"contig_name2\",...] -- Use multiple contigs by name\n
\n
        Example: #{args.help} {\"names\":[\"NODE_1234\",\"NODE_5678\"]} -- Restore edges\n
        to these two contigs (though not necessarily between them...)\n
\n
ccname : \"contig_name\" -- Use all contigs in the connected component containing the\n
        named contig\n
\n
        Example: #{args.help} {\"ccname\":\"NODE_1234\"} -- Restore edges to all contigs\n
        in the connected component containing this contig\n
\n
ccnames : [\"contig_name1\",\"contig_name2\",...] -- Restore edges to all contigs\n
        in the connected components containing these contigs\n
\n
        Example: #{args.help} {\"ccnames\":[\"NODE_1234\",\"NODE_5678\"]} -- Restore edges\n
        to all contigs within the connected compontent(s) containing these two contigs\n
\n
radius : <int> -- Expand the sphere of restored connections to neighbors\n
\n
        Default: #{args.help} {\"radius\":0} -- Restore edges only to neighboring contigs\n
\n
        Example: #{args.help} {\"radius\":2} -- Restore edges to all contigs within two\n
        hops of the selected contigs (including along newly restored paths)\n
\n
complete : true -- Restore connections among all types of contigs\n
\n
        Example: #{args.help} {\"complete\":true} -- Restore connections between contigs\n
        which are both removed from, and part of, currently selected connected components\n
\n
        Default: #{args.help} {\"complete\":false} -- The 'existing' parameter (below)\n
        determines which type of connections are restored.\n
\n
existing : true -- Only restore connections between contigs currently part of selected\n
        connected components. This parameter has no effect when the 'complete' parmeter\n
        (above) is true.\n
\n
        Example: #{args.help} {\"existing\":true} -- Only restore connection between\n
        currently selected contigs\n
\n
        Default: #{args.help} {\"existing\":false} -- Only restore connection between\n
        currently selected and currently unselected contigs; That is, no new connections\n
        within the contigs of currently selected connected components are restored.\n
\n
problems : true -- Restore all connections to contigs that are marked with potential\n
        assembly problems\n
\n
        Example: #{args.help} {\"problems\":true} -- Reconnect all problem contigs to\n
        all of their potential neighbors.  Uses radius:1 and existing:true\n
")
      callback?(null, j)
      return

   args.problems ?= false

   if args.problems
      args.radius = 1
      args.existing = true
      args.names = (nid for nid, n of j.nodes when n.contig_problems?)

   args.complete ?= false
   args.existing ?= false
   args.thresh ?= 0.0
   args.radius ?= 0

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

   try
      build_graph_refs(j) 
      build_removed_graph_refs(j)
   catch err
      callback?(err, null)
      return
      
   node_list = {}

   if args.ccname?
      args.ccnames = [args.ccname]
      
   if args.ccnames?
      args.names ?= []
      for id in args.ccnames
         for c in j.connected_comps when c.indexOf(id) != -1
            args.names = args.names.concat(c)            

   if args.name?
      args.names = [args.name]

   if args.radius and args.names?
      for name in args.names
        for nid in find_neighbors(j.nodes[name], args.radius) 
           node_list[nid] = j.nodes[nid]  
   
   if args.names?
      for id in args.names
         node_list[id] = j.nodes[id]
   else
      node_list = j.nodes

   in_nodes_added = {}
   out_nodes_added = {}
   nodes_added = {}
   edges_added = []

# This loop tries to add removed edges connecting removed nodes to scaffold nodes.
# They are only selected if the mate-pairing gives the proper orientation for the node.
   
   unless args.existing
   
      for id, n of node_list
         for e in n.rem_links when j.removed_nodes[e[1].id]? and edge = edge_director(n, e[1], e[0]) and e[0].score >= args.thresh
            if id is edge.n1
               out_nodes_added[e[1].id] ?= [] 
               out_nodes_added[e[1].id].push(e[0])
            else 
               in_nodes_added[e[1].id] ?= [] 
               in_nodes_added[e[1].id].push(e[0])

      for id, edge_list of in_nodes_added
         edges_added = edges_added.concat(edge_list)
         nodes_added[id] = j.removed_nodes[id]
      for id, edge_list of out_nodes_added
         edges_added = edges_added.concat(edge_list)
         nodes_added[id] = j.removed_nodes[id]

# This code adds properly directed edges between previously removed nodes.
      for e in j.removed_edges when (e? and (nodes_added[e.n1]? and nodes_added[e.n2]?)) and e.score >= args.thresh                        
         if (edge = edge_director(e.src, e.tar, e)) 
            edges_added.push(edge)

# This code optionally adds additional removed edges between existing scaffold nodes 
   if args.complete or args.existing
      for id, n of node_list
         if (not args.ends) or (n.links.length is 1)
            for e in n.rem_inlinks when e[0].score >= args.thresh
               if j.nodes[e[1].id]? and edge = edge_director(n, e[1], e[0])
                   edges_added.push(edge)
                   delete j.removed_edges[edge.index]

# Commit to new nodes.
   for id, n of nodes_added
      j.nodes[id] = n
      delete j.removed_nodes[id]

# Commit to new edges.
   for e in edges_added
      j.edges.push(e)
      delete j.removed_edges[e.index]

   delete j.connected_comps

   remove_graph_refs(j, (n)->delete n.seen) 
   
   callback?(null, j)

###
# Full Ordering -- This is like Relink above, but it inserts extra dependency edges 
# between added nodes based on mate-pair mapping positions in scaffold nodes
###
full_order = (j, args={"verbose":false}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Reconnect contigs with ambiguous placement within a scaffold using mapped\n
          pair position information to resolve ambiguities
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
ccname : \"contig_name\" -- Use the connected component containing the named contig\n
\n
        Example: #{args.help} {\"ccname\":\"NODE_1234\"} -- Insert removed contigs\n
        within the connected component containing this contig\n
\n
ccnames : [\"contig_name1\",\"contig_name2\",...] -- Insert removed contigs within the\n
        connected components containing these contigs\n
\n
thresh : <float> -- Minimum contig connection bitscore to consider\n
\n
        Example: #{args.help} {\"thresh\":250.0} -- Default. Only connections scoring\n
        at least 250.0 bits will count in the calculations\n
\n
min_pos_diff : <int> -- Minimum mapping position difference (within a neighboring contig)\n
        considered reliable enough to use for resolving positional ambiguities.\n
\n
        Example: #{args.help} {\"min_pos_diff\":75} -- Default. Only pairing information\n
        yielding a relative position difference of 75 nucleotides (between candidate\n
        contigs and a neighboring existing contig) will be considered significant enough\n
        to use in determining the relative ordering of the candidate contigs\n
\n
dup_kmer : <int> -- Length of k-mers to use for detecting variant duplicate contigs\n
\n
        Example: #{args.help} {\"dup_kmer\":14} -- Default. Use 14-mers to detect\n
        duplicate variant contigs before inserting them between scaffold backbone contigs.\n

        Example: #{args.help} {\"dup_kmer\":0} -- Disable duplicate variant contig checks.\n
\n
dup_thresh : <float> -- Fraction of kmers from a contig with hits to scaffold backbone.\n
\n
        Example: #{args.help} {\"dup_thresh\":0.15} -- Default. If more than 15% of kmers\n
        for a contig hit kmers from the scaffold backbone contigs, then do not insert\n
        that contig.\n
\n
verbose : true -- output diagnostics on STDERR\n
\n
        Example: #{args.help} {\"verbose\":true} -- Generate extra output information\n
")
      callback?(null, j)
      return

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

   try
      build_graph_refs(j) 
      build_removed_graph_refs(j)
   catch err
      callback?(err, null)
      return

   cc_list = []

   args.thresh ?= 250.0
   args.min_pos_diff ?= 75

   args.dup_kmer ?= 14
   args.dup_thresh ?= 0.15

   if args.ccname?
      args.ccnames = [args.ccname]
      
   if args.ccnames?
      for id in args.ccnames
         for c, i in j.connected_comps when c.indexOf(id) != -1
            cc_list.push(i)            
   else
      cc_list = [0...j.connected_comps.length]

   nodes_added = {}
   edges_added = []
   
   cc_chain_levels = []

# This loop walks through the CCs, and for each CC, walks through scaffolded nodes in
# connection order.  Assumes the CCs are linear scaffolds.

   for idx in cc_list

      in_nodes_added = {}
      out_nodes_added = {}
      shared_seq_nodes = {}
      in_idx = {}
      out_idx = {}
   
      c = j.connected_comps[idx]
      
      try
         [[head], [tail]] = find_ends(j, c, true)
      catch err
         remove_graph_refs(j)
         callback?(err, null)
         return
         
      chain_levels = []
      n = head
      i = 1
      
      while n
         chain_levels[i-1] = n
         n.chain_idx = i
         n.ccnum = idx
         n.cclen = c.length
         i++
         n = n.outlinks[0]?[1]
      
      cc_chain_levels.push(chain_levels)
      
      calc_mers = (k, seq, mers = {}) ->
         for c in [0..seq.length-k]
            mers[seq[c...c+k]] = true
         mers 

      shared_mers = (ref_mers, query_mers) ->
         q_mers = Object.keys(query_mers)
         total = 0
         total++ for q in q_mers when ref_mers[q]?
         total / q_mers.length             

      mers_shared = 0.0
      all_mers_shared = 0.0

      if args.dup_kmer
         all_mers = {}
         for n, l in chain_levels when n.recon_seq?.length >= args.dup_kmer
            all_mers = calc_mers(args.dup_kmer, n.recon_seq, all_mers)
         unless Object.keys(all_mers).length
            console.warn("INSERT parameter dup_kmer = #{args.dup_kmer} but graph has no sequence. Disabling INSERT kmer checking.") if args.verbose
            args.dup_kmer = 0
      
      for n, l in chain_levels
         for e in n.rem_links when j.removed_nodes[e[1].id]? and edge = edge_director(n, e[1], e[0])

            other = if n.id is edge.n1 then edge.n2 else edge.n1

            if args.dup_kmer and e[1].recon_seq?.length >= args.dup_kmer
               contig_mers = calc_mers(args.dup_kmer, e[1].recon_seq)
               all_mers_shared = shared_mers(all_mers, contig_mers)
            
            if all_mers_shared > args.dup_thresh 
               shared_seq_nodes[other] ?= []
               shared_seq_nodes[other].push(n.id)
               console.warn("Shared sequence variant node #{other} rejected from #{n.id} #{args.dup_kmer}_mers shared: #{all_mers_shared}") if args.verbose
            else 
               if n.id is edge.n1 
                  out_nodes_added[e[1].id] ?= [] 
                  out_nodes_added[e[1].id].push(e[0])
                  out_idx[e[1].id] ?= []
                  out_idx[e[1].id].push(n.chain_idx)
               else 
                  in_nodes_added[e[1].id] ?= [] 
                  in_nodes_added[e[1].id].push(e[0])
                  in_idx[e[1].id] ?= []
                  in_idx[e[1].id].push(n.chain_idx)
            
      # Completely remove candidate nodes/edges that are shared sequence variants with 
      # nodes of the existing scaffold backbone, as they cannot be properly inserted in
      # any gap without causing artifactual repeats.
      for id of shared_seq_nodes  
         delete in_nodes_added[id]
         delete in_idx[id]
         delete out_nodes_added[id]
         delete out_idx[id]
         
      # Remove candidate nodes / edges that don't have both in and out edges from the existing scaffold
      # Only consider nodes with both in and out edges to the scaffold nodes
      for id, edge_list of in_nodes_added when not out_nodes_added[id]? 
         delete in_nodes_added[id]
         delete in_idx[id]
      for id, edge_list of out_nodes_added when not in_nodes_added[id]? 
         delete out_nodes_added[id] 
         delete out_idx[id]
   
      chain_minmax = {}   
      for id of in_idx

         # Detect if this node is the smaller contig of a self-linking pair
         # Pass on it if so
         check_list = {}
         for [e,n] in j.removed_nodes[id].rem_links when in_idx[n.id]?
            check_list[n.id] ?= 0
            check_list[n.id]++
         die = false   
         for nid, cnt of check_list when cnt > 1 and j.removed_nodes[id].seq_len < j.removed_nodes[nid].seq_len 
            die = true
            console.warn("Circular variant node #{id} rejected") if args.verbose
            break
         
         break if die
   
         dup = false
   
         if nodes_added[id]?
            j.removed_nodes[id].dup ?= 0
            j.removed_nodes[id].dup++
            node_id = id + "_Dup_#{j.removed_nodes[id].dup}"
            node = {}
            for p of j.removed_nodes[id]
               node[p] = j.removed_nodes[id][p]
            node.name = node_id
            node.id = node_id
            dup = true
         else 
            node_id = id
            node = j.removed_nodes[id]
            
         # Calulate the min and max chain_indices for both the incoming and outgoing edge scaffold nodes
         chain_minmax = 
            ins : 
               min : Math.min.apply(null,out_idx[id])
               max : Math.max.apply(null,out_idx[id])
            outs : 
               min : Math.min.apply(null,in_idx[id])
               max : Math.max.apply(null,in_idx[id])

         if ((chain_minmax.ins.max < chain_minmax.outs.min) and (chain_minmax.outs.max - chain_minmax.ins.min < 7))  # The "normal" case
            edges_added = edges_added.concat(out_nodes_added[id]).concat(in_nodes_added[id])
            nodes_added[node_id] = node
            node.slot = [chain_minmax.ins.min, chain_minmax.outs.max]
            node.ccnum = idx
         else if (not (tail.terminal) and
                  (chain_minmax.ins.min > chain_minmax.outs.max) and 
                  (chain_minmax.outs.max - chain_minmax.outs.min < 4) and 
                  (chain_minmax.ins.max - chain_minmax.ins.min < 4) and 
                  (chain_minmax.ins.min - chain_minmax.outs.max > j.connected_comps[idx].length - 7))    # Nodes entirely "between" the scaffold ends
            edges_added = edges_added.concat(out_nodes_added[id])  # Only add edges from the tail end
            nodes_added[node_id] = node
            node.slot = [chain_minmax.ins.min, Number.MAX_VALUE]
            node.ccnum = idx
         else if ((chain_minmax.outs.min < chain_minmax.ins.min <= chain_minmax.ins.max < chain_minmax.outs.max) and 
                  (chain_minmax.ins.max - chain_minmax.ins.min < 4) and 
                  (chain_minmax.outs.max - chain_minmax.outs.min > j.connected_comps[idx].length - 7))  # Node is within tail end, but reaches across to head
            new_out_edges = (e for e in in_nodes_added[id] when e.tar.chain_idx > chain_minmax.ins.max)
            new_out_idx = (e.tar.chain_idx for e in new_out_edges)
            edges_added = edges_added.concat(out_nodes_added[id]).concat(new_out_edges)  # Only add edges from the tail end
            nodes_added[node_id] = node
            node.slot = [chain_minmax.ins.min, Math.max.apply(null,new_out_idx)]
            node.ccnum = idx
         else if ((chain_minmax.ins.min < chain_minmax.outs.min <= chain_minmax.outs.max < chain_minmax.ins.max) and 
                  (chain_minmax.outs.max - chain_minmax.outs.min < 4) and 
                  (chain_minmax.ins.max - chain_minmax.ins.min > j.connected_comps[idx].length - 7))   # Node is within head end, but reaches across from tail
            new_in_edges = (e for e in out_nodes_added[id] when e.src.chain_idx < chain_minmax.outs.min)
            new_in_idx = (e.src.chain_idx for e in new_in_edges)
            edges_added = edges_added.concat(in_nodes_added[id]).concat(new_in_edges)  # Only add edges from the tail end
            nodes_added[node_id] = node           
            node.slot = [Math.min.apply(null,new_in_idx), chain_minmax.outs.max]
            node.ccnum = idx
         else  # Strange case, reject the node and all edges
            console.warn("Node rejected: #{node_id}") if args.verbose

         # Fix up the a duplicated node name in the edges now that it hasn't been rejected 
         if nodes_added[node_id]? and dup
            for e in out_nodes_added[id]
               e.n2 = node_id
            for e in in_nodes_added[id]
               e.n1 = node_id

# This code adds properly directed edges between previously removed nodes.
   for e in j.removed_edges when (e? and (nodes_added[e.n1]? and nodes_added[e.n2]?))                         
      if ((edge = edge_director(e.src, e.tar, e)) and 
           (edge.src.slot[0] <= edge.tar.slot[0]) and 
           (edge.src.slot[1] <= edge.tar.slot[1]) and
           ((edge.tar.slot[0] - edge.src.slot[1]) < 7) and
           (edge.src.ccnum is edge.tar.ccnum))  
         edges_added.push(edge)

# Commit to new nodes.
   for id, n of nodes_added
      j.nodes[id] = n
      delete j.removed_nodes[id]

# Commit to new edges.
   for e in edges_added
      j.edges.push(e)
      delete j.removed_edges[e.index]

   unless args.no_mp_pos_links?

      remove_graph_refs(j)
      
      try
         build_graph_refs(j, 
            ((n) -> 
               n.out_nodes = {}
               n.in_nodes = {}), 
            ((e) -> 
               e.src.out_nodes[e.n2] = e
               e.tar.in_nodes[e.n1] = e))
      catch err
         callback?(err, null)
         return
         
      dep_edges_added = []
   
      add_dep_edges = (rank_list) ->
         rank_list.sort((a,b)->(a[0]-b[0]))
         new_tail = null
         for [pos, new_head] in rank_list
            if (new_tail and 
                  not new_tail.out_nodes[new_head.id]? and 
                  not new_tail.in_nodes[new_head.id]? and
                  pos - prev_pos > args.min_pos_diff)
               dep_edges_added.push(
                  n1  :  new_tail.id 
                  n2  :  new_head.id
                  dir :  "pos"
                  pos_diff : pos - prev_pos
                  )
            new_tail = new_head
            prev_pos = pos

      for idx in cc_list
         
         chain_levels = cc_chain_levels[idx]
         
         for n, l in chain_levels
   
            base_p = 0
            base_n = 0
            
            cur_lev = {}
            pos_rank = []
    
            prev = chain_levels?[l-1]?.id
            next = chain_levels?[l+1]?.id
   
            if prev
               base_p = n.in_nodes[prev].p1
               for [edge, node] in j.nodes[prev].outlinks
                  cur_lev[edge.n2] = [edge]
   
            if next
               base_n = n.out_nodes[next].p2
               for [edge, node] in j.nodes[next].inlinks
                  cur_lev[edge.n1] ?= []
                  cur_lev[edge.n1][1] = edge
                     
            for id, [e_p, e_n] of cur_lev
               if e_p and e_n and (e_p.bits + e_n.bits > args.thresh)
                  pos = (((e_p.bits/(e_p.bits+e_n.bits))*(e_p.p1-base_p)) + 
                         ((e_n.bits/(e_p.bits+e_n.bits))*(e_n.p2-base_n)))
                  pos_rank.push([pos, e_p.tar]) 
               else if e_p and (e_p.bits > args.thresh) and not (j.nodes[next]?.out_nodes?[e_p.tar.id]?)
                  pos = e_p.p1-base_p
                  pos_rank.push([pos, e_p.tar])   
               else if e_n and (e_n.bits > args.thresh) and not (j.nodes[prev]?.in_nodes?[e_n.src.id]?)
                  pos = e_n.p2-base_n
                  pos_rank.push([pos, e_n.src]) 
   
            add_dep_edges(pos_rank)

      # Check added edges for inconsistencies
      edge_obj = {}
      for e in dep_edges_added
         edge_obj["#{e.n1}_#{e.n2}"] ?= 0
         edge_obj["#{e.n1}_#{e.n2}"]++

      # Commit to new edges.
      for e in dep_edges_added
         unless edge_obj["#{e.n2}_#{e.n1}"] or --edge_obj["#{e.n1}_#{e.n2}"]
            # Check for 3-edge cycles made by a dependency edge
            good = true
         
            for [edge, node] in j.nodes[e.n2].outlinks when node.out_nodes[e.n1]
                  good = false
                  break
         
            j.edges.push(e) if good

   remove_graph_refs(j, (n) -> (
      delete n.in_nodes
      delete n.out_nodes
      delete n.ccnum
      delete n.cclen
      delete n.chain_idx
      delete n.slot)) 

   delete j.connected_comps
   
   callback?(null, j)

###
# Layout scaffolds
###
scaffold = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Lay out a fully linear scaffold from contigs in unambiguously ordered\n
         connected components
")
      if args.detailed_help?
         console.warn("
\n
Parameters: NONE.\n
")
      callback?(null, j)
      return

   if j.scaffolds?
      delete j.scaffolds
      if j.clusters?
         delete j.clusters

   calc_ccomps(j) unless j.connected_comps?  # Make sure there are ccomps

   try
      ccomp_ends = find_all_ends(j)
      build_graph_refs(j, ((n) -> n.in_edges = {}), ((e) -> e.tar.in_edges[e.src.id] = e))
   catch err
      callback?(err, null)
      return
   
   all_tsort = []
   all_scaffs = []
   in_scaffs = {}
   
   for [source_nodes, sink_nodes], cci in ccomp_ends

      # Generate Topological Sort of nodes 
      tsort = []
      for n in source_nodes
         n.in_edges = {}
         
      while n = source_nodes.pop()
         tsort.push(n.id)
         for e in n.outlinks
            delete e[1].in_edges[n.id]
            if Object.keys(e[1].in_edges).length is 0
               source_nodes.push(e[1])
         # Try to rescue a graph with a dependent edge cycle
         if not args.no_rescue and (source_nodes.length is 0) and (tsort.length isnt j.connected_comps[cci].length)
            # Walk all of the nodes in the CCOMP not in the tsort and with one remaining in_edge
            remove_cands = ([j.nodes[id], j.nodes[in_keys[0]]] for id in j.connected_comps[cci] when ((tsort.indexOf(id) is -1) and 
                                                  ((in_keys = Object.keys(j.nodes[id].in_edges)).length is 1) and 
                                                  (j.nodes[id].in_edges[in_keys[0]].pos_diff?)))

            if remove_cands.length
               # If there are any, find the one with the minimum position difference (probably the weakest)
               remove_cands.sort((a,b)->(a[0].in_edges[a[1].id].pos_diff - b[0].in_edges[b[1].id].pos_diff))
               [in_node, out_node] = remove_cands[0]
               console.warn("Rescuing node #{remove_cands[0][0].id}: position dependent edge from #{remove_cands[0][1].id}")
               # remove this edge and proceed
               for edge, ei in out_node.outlinks when edge[0] is in_node.in_edges[out_node.id]
                  out_node.outlinks.splice(ei,1)  # remove the offending out_edge
                  break
               delete in_node.in_edges[out_node.id]
               source_nodes.push(in_node)
         
      all_tsort.push(tsort)
      
      len_to = {}
      
      unless tsort.length
         callback?(new Error("SCAFF: No nodes in topological sort list. Circular scaffold?"), null)
         return
      
      max_id = tsort[0]
      max_len = j.nodes[max_id].seq_len
      
      # Calculate the longest sequence path along the topo sorted nodes
      for id in tsort
         n = j.nodes[id]
         len_to[id] ?= { len : n.seq_len, edge : null }
            
         for [e,n2] in n.outlinks when ((not len_to[n2.id]?) or (len_to[n2.id].len < len_to[id].len + n2.seq_len))
            len_to[n2.id] = { len : len_to[id].len + n2.seq_len, edge : e } 
            if len_to[n2.id].len > max_len
               max_id = n2.id 
               max_len = len_to[n2.id].len

      longest_chain = [j.nodes[max_id]]
      len_to[max_id].edge?.keep = true
      prev = len_to[max_id]
      
      while prev.edge?
         longest_chain.unshift(prev.edge.src)
         prev.edge.keep = true
         prev = len_to[prev.edge.src.id]
 
      for n, i in longest_chain
         all_scaffs.push(n.id)
         in_scaffs[n.id] = i
      
      if (sink_nodes.indexOf(longest_chain[longest_chain.length-1]) is -1)
         end_id = longest_chain[longest_chain.length-1].id
         console.warn("Scaffold Warning: CCOMP #{cci} was not completely traversed. #{end_id} #{(100.0*(tsort.length/j.connected_comps[cci].length)).toFixed(1)}% #{(100.0*(longest_chain.length/j.connected_comps[cci].length)).toFixed(1)}% ")    
      
   for e in j.edges when not e?.keep?
      j.removed_edges.push(e)
      delete j.edges[e.index]

   for id, n of j.nodes when not in_scaffs[id]?
      j.removed_nodes[id] = n
      delete j.nodes[id]

   remove_graph_refs(j, 
      (n) -> (
         delete n.in_edges
         delete n.saved_outlinks),
      (e) -> (
         delete e.keep
      ))
         
   calc_ccomps(j)
   
   try
      build_graph_refs(j) 
   catch err
      callback?(err, null)
      return
   
   # Build the scaffolds data structure
   j.scaffolds = {}
   for c, i in j.connected_comps
      scaf_name = "Scaffold_#{i}"
      j.scaffolds[scaf_name] = {"ccnum":i,"nodes":[]}
      
      try
         [[node]] = find_ends(j,c,true)
      catch err
         remove_graph_refs(j)
         callback?(err, null)
         return
         
      while node
         if node.recon_seq?
            j.scaffolds[scaf_name].nodes.push(node.id)
         node = node.outlinks[0]?[1] 

   remove_graph_refs(j)

   callback?(null, j)

###
# Make explicit edits to the current assembly graph
###
perform_edits = (j, args={}, callback) -> 
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Make manual explicit edits to the selected graph structure
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
rem_nodes : [\"contig1\",\"contig2\",...] -- Remove the specified contig nodes by name\n
\n
        Example: #{args.help} {\"rem_nodes\":[\"NODE_1234\",\"NODE_5678\"]} -- Remove\n
        these two contigs from the graph\n
\n
add_nodes : [\"contig1\",\"contig2\",...] -- Add back the specified contig nodes by name\n
\n
        Example: #{args.help} {\"add_nodes\":[\"NODE_1234\",\"NODE_5678\"]} -- Move\n
        these two contigs from the removed pool to the selected graph\n
\n
rem_edges : [\"contig1\",[\"contig2\",...]] -- Remove the specified connections between\n
        contig1 and any number of others.\n
\n
        Example: #{args.help} {\"rem_edges\":[\"NODE_1234\",[\"NODE_5678\",\"NODE_9\"]]}\n
        -- Remove two connection edges from the graph, both connected to NODE_1234\n
\n
add_edges : [\"contig1\",[\"contig2\",...]] -- Add back the specified connections between\n
        contig1 and any number of others.\n
\n
        Example: #{args.help} {\"add_edges\":[\"NODE_1234\",[\"NODE_5678\",\"NODE_9\"]]}\n
        -- Add two connection edges from the graph, both connected to NODE_1234\n
\n
        NOTE: for rem_edges and add_edges above, if only a single edge is involved, then\n
        the parameter syntax may optionally be flattened. For example:\n
        {\"add_edges\":[\"NODE_1234\",[\"NODE_9\"]]} is equivalent to\n
        {\"add_edges\":[\"NODE_1234\",\"NODE_9\"]}\n
\n
        If more than one node-independent sets of edges are to be added or removed (that\n
        is, those not sharing any contig(s) in common), then multiple calls to #{args.help}\n
        are required to accomplish this task.\n
")
      callback?(null, j)
      return

   try
      build_graph_refs(j) 
      build_removed_graph_refs(j)
   catch err
      callback?(err, null)
      return

   args.rem_nodes ?= []
   args.rem_edges ?= []

   args.add_nodes ?= []
   args.add_edges ?= []

   # Neat little function to return a sensible type
   type = do ->
     classToType = {}
     for name in "Boolean Number String Function Array Date RegExp Undefined Null".split(" ")
       classToType["[object " + name + "]"] = name.toLowerCase()
     (obj) ->
        strType = Object::toString.call(obj)
        classToType[strType] or "object"

   for nid in args.add_nodes when j.removed_nodes[nid]?
      console.warn("EDIT: Adding node #{nid}")
      j.nodes[nid] = j.removed_nodes[nid]
      delete j.removed_nodes[nid]

   for [nid1, n2_list] in args.add_edges when j.nodes[nid1]?
      n2_list = [n2_list] if type(n2_list) is "string"
      for e in j.nodes[nid1].rem_links when (n2_list.indexOf(e[1].id) isnt -1) and j.removed_edges[e[0].index]?
         console.warn("EDIT: Adding edge #{e[0].n1} --> #{e[0].n2} Index: #{e[0].index} #{j.removed_edges[e[0].index].n1} --> #{j.removed_edges[e[0].index].n2}")
         delete j.removed_edges[e[0].index]
         e[0].index = j.edges.push(e[0]) - 1

   for nid in args.rem_nodes when j.nodes[nid]?
      console.warn("EDIT: Deleting node #{nid}")
      for e in j.nodes[nid].links when j.edges[e[0].index]?
         console.warn("EDIT: Removing edge #{e[0].n1} --> #{e[0].n2} Index: #{e[0].index} #{j.edges[e[0].index].n1} --> #{j.edges[e[0].index].n2}")
         delete j.edges[e[0].index]
         e[0].index = j.removed_edges.push(e[0]) - 1
         
      j.removed_nodes[nid] = j.nodes[nid]
      delete j.nodes[nid]

   for [nid1, n2_list] in args.rem_edges when j.nodes[nid1]?
      n2_list = [n2_list] if type(n2_list) is "string"
      for e in j.nodes[nid1].links when (n2_list.indexOf(e[1].id) isnt -1) and j.edges[e[0].index]?
         console.warn("EDIT: Removing edge #{e[0].n1} --> #{e[0].n2} Index: #{e[0].index} #{j.edges[e[0].index].n1} --> #{j.edges[e[0].index].n2}")
         delete j.edges[e[0].index]
         e[0].index = j.removed_edges.push(e[0]) - 1
         

   remove_graph_refs(j) 

   callback?(null, j)

###
# create new nodes from parts of the given nodes
###
cut_node = (j, args={}, callback) -> 
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Create a new node from an existing node using the given coordinates
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
name : \"contig_name\" -- Name of the contig node to be copied\n
\n
        Example: #{args.help} {\"name\":\"NODE_1234\"} -- Use sequence from NODE_1234\n
\n
new_name : \"new_contig_name\" -- Name for the newly created contig node\n
\n
        Example: #{args.help} {\"new_name\":\"NODE_1234a\"} -- New node will be named\n
        NODE_1234a\n
\n
include : [\"contig_name1\", ...] -- Move mate-pair edges from the listed contigs.\n
\n
        Example: #{args.help} {\"name\":\"NODE_1234\",\"include\":[\"NODE_567\",\"NODE_234\"]}\n 
        -- Edges between NODE_1234 and the included contigs will be moved to the new cut\n
        and copied version of NODE_1234.\n
\n
auto_include : true -- Automatically move mate-pair edges falling within the cut sequence.\n
\n
        Example: #{args.help} {\"name\":\"NODE_1234\",\"auto_include\":true}\n 
        -- All selected edges between NODE_1234 and any other contigs will be moved to the\n
         new cut and copied version of NODE_1234. Note that currently removed edges will\n
         not be moved. To accomplish this, the node needs to be RELINKed first.
\n
begin : <int> -- Beginning sequence coordinate for the new node within the original\n
end : <int> -- Ending sequence coordinate for the new node within the original\n
\n
        NOTE: To facilitate trimming sequences to a fixed maximum length, it is\n
        allowable for negative 'start' values and positive 'end' values to be longer than\n
        a sequence. In such cases, the values are set to the start and end of the\n
        sequence, respectively. Positive start and negative end positions must fall\n
        within the sequence, because otherwise the start position will be after the end\n
        position.\n
\n
        Example: #{args.help} {\"begin\":123,\"end\":456} -- The new node will include\n
        sequence from positions 123 to 456\n
\n
        Example: #{args.help} {\"end\":456} -- The new node will include sequence from\n
        position 0 (implied) to 456. end may also be omitted, implying the last position\n
\n
        Example: #{args.help} {\"begin\":-1000} -- New node contains at most the last\n
        1000 bases of any sequence. Sequences under 1000 bases are copied unmodified.\n
")
      callback?(null, j)
      return

   node = j.nodes[args.name]
   new_node = {}

   unless node
      callback?(new Error('CUTND requires a valid node "name" argument.'), null)
      return

   # Maintain compatability with "start" instead of "begin"
   if args.start? and not args.begin?
      console.warn("CUTND: Warning: 'start' parameter is deprecated. Use 'begin' instead.")
      args.begin = args.start 

   args.begin ?= 0
   args.end ?= node.seq_len - 1

   args.include ?= []

   if args.begin < 0
      args.begin = node.seq_len + args.begin

   if args.end < 0
      args.end = node.seq_len + args.end

   if args.begin < 0
      args.begin = 0  

   if args.end >= node.seq_len
      args.end = node.seq_len - 1

   unless args.begin < args.end
      callback?(new Error('CUTND requires valid "begin" and "end" arguments.'), null)
      return

   try
      build_graph_refs(j) 
      build_removed_graph_refs(j)
   catch err
      console.warn "build_graph_refs error in CUTND"
      callback?(err, null) 
      return
   
   k = 1
   args.new_name ?= (while j.nodes[node.name + "_Cut_" + k]?
                        k++
                     node.name + "_Cut_" + k)

   j.nodes[args.new_name] = new_node
   
   console.warn("CUTND: Adding #{args.new_name} from #{args.begin} - #{args.end} of #{node.name}")
   
   for p of node  # Copy parameters
      new_node[p] = node[p]

   # Now fix them up
   new_node.name = args.new_name
   new_node.id = args.new_name
   new_node.seq_len = args.end - args.begin + 1
   
   delete new_node.contig_problems
      
   if new_node.per_nt_cov?
      new_node.per_nt_cov = new_node.per_nt_cov.slice(args.begin,args.end)

   if new_node.per_nt_phys_cov?
      new_node.per_nt_phys_cov = new_node.per_nt_phys_cov.slice(args.begin,args.end)

   if new_node.per_nt_mp_ins?
      new_node.per_nt_mp_ins = new_node.per_nt_mp_ins.slice(args.begin,args.end)

   if new_node.recon_seq?
      new_node.recon_seq = new_node.recon_seq.substring(args.begin,args.end)

   if args.auto_include
      args.include.push(e[1].id) for e in new_node.inlinks when args.begin <= e[0].p2 <= args.end
      args.include.push(e[1].id) for e in new_node.outlinks when args.begin <= e[0].p1 <= args.end

   for e in new_node.inlinks when args.include?.indexOf(e[1].id) isnt -1
      console.warn("CUTND: Moving edge #{e[0].n1} --> #{e[0].n2} to #{args.new_name}")
      e[0].n2 = args.new_name
      if e[0].p2 < args.begin
         e[0].p2 = 0
      else if e[0].p2 > args.end
         e[0].p2 = new_node.seq_len - 1
      else
         e[0].p2 = e[0].p2 - args.begin

   for e in new_node.outlinks when args.include?.indexOf(e[1].id) isnt -1 
      console.warn("CUTND: Moving edge #{e[0].n1} --> #{e[0].n2} to #{args.new_name}")
      e[0].n1 = args.new_name
      if e[0].p1 < args.begin
         e[0].p1 = 0
      else if e[0].p1 > args.end
         e[0].p1 = new_node.seq_len - 1
      else
         e[0].p1 = e[0].p1 - args.begin

   remove_graph_refs(j) 

   callback?(null, j)

###
# Report on problem contigs
###
node_problems = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Generate a report of contigs with likely internal assembly issues
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.txt[.gz]\" -- name of txt format file to write statistics to\n
\n
        Example: #{args.help} {\"file\":\"problems.txt\"} -- Write report to the file\n
        problems.txt\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Default. Write report to STDOUT\n
\n
detail : true -- Provide extra per contig detail\n
\n
        Example: #{args.help} {\"detail\":true}\n
")
      callback?(null, j)
      return

   o = open_output_stream(args.file, args.tag)

   # Handle error opening file
   unless o
      callback?(new Error("open_output_stream could not open file '#{args.file}' for output."), null)
      return
   else
      o.on("error", (err)->
          console.error("ERROR: open_output_stream could not write to file '#{args.file}'.")
          callback?(err, null))

   args.detail ?= false
   node_cnt = 0
   total = 0
   prob_types = {}
   for name, node of j.nodes when node.contig_problems?
      o.write("#{name}\t#{node.contig_problems.length}\n") if args.detail
      node_cnt++
      for problem in node.contig_problems
         total++
         prob_types[problem.type] ?= 0
         prob_types[problem.type]++
         if args.detail
            if problem.type is "Physical coverage break"
               if node.per_nt_phys_cov?
                  min = node.per_nt_phys_cov[problem.start]
                  pos = problem.start
                  for x in [problem.end...problem.start] when node.per_nt_phys_cov[x] < min
                     min = node.per_nt_phys_cov[x]
                     pos = x
               o.write("\t#{problem.type}\t#{problem.start}\t#{problem.end}\t#{pos}\n") 
            else
               o.write("\t#{problem.type}\t#{problem.start}\t#{problem.end}\n") 
   
   o.write("\nProblem nodes: #{node_cnt}\n")
   o.write("Total problems: #{total}\n")
   o.write("Problem breakdown:\n")
   for type, count of prob_types
      o.write("\t#{type}: #{count}\n")
   
   o.end((error) ->
         if error
            callback?(error, null)
         else
            callback?(null, j))
   
###
# Invoke tetracalc on current scaffolds, gather results
###
tetracalc = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done
     
   if args.help?
      console.warn("
#{args.help} -- Cluster scaffolds using the tetracalc tool
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
options : <string> -- Command line options for the tetracalc tool. Otherwise defaults used.\n
\n
        NOTE: Run with --help for help with the settings offered by tetracalc.\n
\n
        Example: #{args.help} {\"options\":\"--merge_tar=0.95 -m 7500\"}\n
\n
        Example: #{args.help} {\"options\":\"--help\"}\n
\n
        Example: #{args.help} {\"options\":\"--fixed -t 0.9 -s 0.8 -r 0.9\"}\n
")
      callback?(null, j)
      return
   
   unless j.scaffolds?
      callback?(new Error("No valid scaffolds found. Run SCAFF first."), null)
      return   

   unless Object.keys(j.scaffolds).length > 1
      callback?(new Error("There must be at least two scaffolds to run CLUST."), null)
      return   

   console.warn("CLUST: invoking tetracalc tool")
      
   cluster_json = ""
   
   if args.options
      options = args.options.split(/\s+/)
      options.push('-')
   else   
      options = ['-']
      
   tet = child_process.spawn('tetracalc',options,{})
   tet.stdout.setEncoding('ascii')      
   tet.stderr.setEncoding('ascii')      

   tet.stdout.on('data', (data) ->
      cluster_json = cluster_json.concat(data) 
   ) 
      
   tet.stderr.on('data', (data) ->
      process.stderr.write(data)
   )   

   tet.stdout.on('end', () ->
      try
         j.clusters = JSON.parse(cluster_json).clusters
      catch err
         callback?(new Error("Failed to parse output of tetracalc program. #{err.message}"), null)
         return
            
      callback?(null, j)
   )

   tet.on('exit', (code) -> 
         if (code) 
            console.error('ERROR: ' + code)
            callback?(new Error('SEAStAR program "tetracalc" could not be executed. Please ensure that it is in your PATH'), null)
   )
   
   export_fasta(j, {"scaff":true,"stream":tet.stdin,"verbose":args.verbose})
   
###
# Select clusters for output or further processing
#
# args (mutually exclusive):
#
# clustnum : Choose a cluster by number
# clustnums : Choose a list of clusters
# clustrange : Choose multiple clusters by a range of numbers 
# shift : Remove the first cluster
#
###
grab_clusts = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Select clusters of scaffolds for further processing
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
clustnum : <int> -- Select a specific cluster\n
\n
        Example: #{args.help} '{\"clustnum\":0}' -- select the first cluster\n
\n
clustnums : [<int1>, <int2>, ...] -- Select specific clusters\n
\n
        Example: #{args.help} '{\"clustnums\":[0,3,5]}' -- select these three clusters\n
\n
clustrange : [<int1>, <int2>] -- Select a range of clusters. See also: 'shift' option below.\n
\n
        Example: #{args.help} '{\"clustrange\":[0,5]}' -- select the first six clusters\n
\n
        Example: #{args.help} '{\"clustrange\":[-5,-1]}' -- select the last five clusters\n
\n
        Example: #{args.help} '{\"clustrange\":[0,-2]}' -- select all clusters except the\n
        last one\n
\n
shift : true -- Select all clusters except the first one.\n
\n
        Example: #{args.help} {\"shift\":true} -- Drop the first cluster.\n
\n
        This is like #{args.help} {\"clustrange\":[1,-1]} except it doesn't generate a\n
        fatal error when there is only one remaining cluster, allowing processing to\n
        potentially continue in any calling SCRIPT commands.\n
\n
exclusive : true -- Remove all scaffolds outside of the selected clusters\n
\n
        Example: #{args.help} '{\"exclusive\":false}' -- Default. Keep all unselected scaffolds.\n
")

      callback?(null, j)
      return

   unless j.clusters?  # Make sure there are clusters
      callback?(new Error("No valid scaffold clusters found. Run CLUST first."), null)
      return

   try
      build_graph_refs(j) 
   catch err
      callback?(err, null)
      return

   if args.shift
      args.clustrange=[1,-1]

   if args.clustnum?
        args.clustnums = [args.clustnum]
   else if args.clustrange?
        args.clustrange[0] = (j.clusters.length + args.clustrange[0]) if args.clustrange[0] < 0
        args.clustrange[1] = (j.clusters.length + args.clustrange[1]) if args.clustrange[1] < 0
        args.clustnums = [args.clustrange[0]..args.clustrange[1]]
   else
      args.clustnums ?= [0]
   
   new_cclist = []
   new_clusts = []
   new_scaffs = {}
   
   for cl in args.clustnums
      unless j.clusters[cl]?
         remove_graph_refs(j)
         if args.shift
            callback?(null, undefined)
         else
            callback?(new Error("Invalid selected clustnum: #{cl}. clustnums are zero-based and there are only #{j.clusters.length} clusters in this graph."), null)
         return
      new_clusts.push(j.clusters[cl])   
      for scaf_id in j.clusters[cl]
         if args.exclusive 
            new_cclist.push(j.scaffolds[scaf_id].ccnum)
         new_scaffs[scaf_id] = j.scaffolds[scaf_id]
         
   if args.exclusive      
      grab_ccomps(j,{"ccnums" : new_cclist})

   j.clusters = new_clusts
   j.scaffolds = new_scaffs
 
   remove_graph_refs(j) 
   
   callback?(null, j)
   
###
# DUMP -- Output current JSON structure to a file
#
# args:
#
# file : Filename to write
#
###
write_json = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Output current data structure to a file, or by default to STDOUT
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.json[.gz]\" -- Specify a file name to write the JSON format sequence graph to.\n
\n
        If the filename contains one or more `#` characters in a row, these positions are\n
        replaced with zero-padded digits that will increment each time a file is written to\n
        this filename pattern. If no `#` characters are present, then this command overwrites\n
        any existing file of the same name.\n
\n
        Example: #{args.help} {\"file\":\"my_assembly.json\"} -- Write data to the file\n
        my_assembly.json\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Write data piped to STDOUT (default)\n
\n
        Example: #{args.help} {\"file\":\"my_assembly_###.json\"} -- Write data to the file:\n
        my_assembly_000.json (first time this is run)\n
\n
        Example: #{args.help} {\"file\":\"my_assembly_###.json\"} -- Run again, write data\n
        to the file: my_assembly_001.json\n 
")
      callback?(null, j)
      return

   o = open_output_stream(args.file, args.tag)
   
   # Handle error opening file
   unless o
      callback?(new Error("open_output_stream could not open file '#{args.file}' for output."), null)
      return
   else
      o.on("error", (err)->
          console.error("ERROR: open_output_stream could not write to file '#{args.file}'.")
          callback?(err, null))

   my_stringify(j, o, 2, (err) ->
      if err
         callback?(err, null)
      else       
         o.end((error) ->
            if error
               callback?(error, null)
            else
               callback?(null, j)))

###
# read_json -- 
#
# args:
#
# file : Filename to read
#
###
read_json = (j, args={'file':'-'}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Input JSON sequence graph data from a file, or by default from STDIN.
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename.json[.gz]\" -- Specify the name of a JSON format sequence graph file\n
\n
        Example: #{args.help} {\"file\":\"my_assembly.json\"} -- Read data from the file\n
        my_assembly.json\n
\n
        Example: #{args.help} {\"file\":\"-\"} -- Read data piped from STDIN (default)\n
\n
        NOTE: On the command line, if the first parameter isn't a valid command string,\n
        and it ends in `.json[.gz]`, then it is assmued to be the name of a JSON file,\n
        and an implicit `LOAD` command will be run using that filename.\n
")
      callback?(null, j)
      return

   read_input_stream(args.file, true, (error, read_buffer) -> 
           if error or not read_buffer
                 console.error("ERROR: Empty buffer returned for file: #{args.file}")
                 callback?(error or new Error("Empty buffer returned for file: #{args.file}"), null)
                 return
           else
              j = read_buffer
              j.processing ?= []
              j.processing.push(["$",ss_version,'LOAD',clone_object(args)])
              callback?(null, j)) 

###
# command_file -- Implementation of the SCRIPT command 
###
command_file = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Use contents of file as a series of commands, or read from the console if\n
          no file provided
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
file : \"filename[.gz]\" -- Read commands from the named file\n
\n
        Example: #{args.help} {\"file\":\"my_script.go\"} -- Read and run commands from\n
        the file my_script.go\n
\n
        Example: #{args.help} -- Enter an interactive (command line) session at this point\n
        typing commands one at a time\n
\n
tag : \"filename_part\" -- Part of a filename to include in files written from this script\n
          NOTE: tag renaming is not used within interactive SCRIPT command sessions.\n
\n
        Example: #{args.help} {\"file\":\"my_script.go\", \"tag\":\"Run1\"} -- Add the \n
        string \"Run1\" in place of the '@' character in the filenames of any output files\n
        written by commands in the script my_script.go (e.g. \"@_output.json\" would\n
        become \"Run1_output.json\").\n
\n
        Example: #{args.help} {\"file\":\"my_script.go\", \"tag\":\"Sample_A\"} -- Add the \n
        string \"Sample_A\" in place of runs of '@' characters anywhere in the file path\n
        of files written by commands in the script my_script.go, for example:\n
        \"~/data/samples/@@@/@@@_output.json\" would become:\n
        \"~/data/samples/Sample_A/Sample_A_output.json\"\n
")
      callback?(null, j)
      return

   queue_commands = (error, read_buffer) ->
   
      if error or not read_buffer
         callback?(error or new Error("empty read_buffer after reading SCRIPT file."), null)
         return
      
      lines = read_buffer.split("\n")
      script_cmds = []
      for l in lines
         [cmd, argstr] = l.trim().replace(/\t/," ").split(" ",2)   
         cmd = cmd.toUpperCase()
         if typeof(commands[cmd]) is 'function'
            if argstr        
               try
                  cmd_args = JSON.parse(argstr)
               catch err
                  console.error("ERROR: Could not parse JSON arguments '#{argstr}' to '#{cmd}' command")
                  callback?(err, null)
                  return
            else
               cmd_args = null               
            
            if cmd_args?.file? and args.tag? and not (cmd is 'LOAD')
               cmd_args.tag ?= args.tag 
               
            if (cmd is 'SCRIPT') and not cmd_args?.file? 
               callback?(new Error("Invalid command: You may not launch an interactive SCRIPT session from within a SCRIPT."), null)
               return
            else 
               script_cmds.push([commands[cmd], cmd_args, cmd])
         else if cmd and cmd[0] isnt '#'
            callback?(new Error("Invalid command: '#{cmd}'"), null)
            return 

      read_buffer = null

      script_next_cmd = (err, return_j) ->
         if err
            err_str = "SCRIPT aborted due to error in command execution."
            unless err.message is no_data_str  # Ignore, it will be regenerated below if necessary 
               unless err.message is err_str
                  console.error("ERROR terminating SCRIPT processing: #{err.message}\n")
                  err = Error(err_str)
               callback?(err, null)
               return
   
         if (script_cmds.length)  # If there are commands remaining to process
            [cmd_func, cmd_args, cmd] = script_cmds.shift()    # Prepare the next command

            unless return_j?.nodes? or cmd is 'LOAD' or cmd is 'SCRIPT' or cmd is 'UNSTASH' or cmd is 'HELP'
               callback?(new Error(no_data_str), {"processing" : (j?.processing or [])})  # Empty JSON object, except preserve command history
               return   

            # Record history of this command in the graph object   
            return_j?.processing?.push([args.file, ss_version, cmd, clone_object(cmd_args)]) unless cmd is 'HELP'  
            console.warn("Executing (#{args.file}) #{cmd} #{JSON.stringify(cmd_args)}")  # Output progress to stderr
            setImmediate(cmd_func, return_j, cmd_args, script_next_cmd)  # Invoke the command with the graph object, args, and callback function
         else
            callback?(null, return_j)

      script_next_cmd(null, j) 

   if args.file?
      read_input_stream(args.file, false, queue_commands)
   else
      unless process.stderr._type? and process.stderr._type is 'tty'  
         console.log("ERROR: Interactive use of the SCRIPT command requires that STDERR not be redirected.\n")
         callback?(new Error("Interactive use of the SCRIPT command requires that STDERR not be redirected."), null)
         return

      console.warn("\nEntering interactive SCRIPT mode, SEAStAR Version: #{ss_version}\nFor general help, quit and rerun with the -h option. Type 'HELP' for help with commands.\nType '.save <filename>' to save commands from this session to an output file.\nPress <Ctrl-D> to quit.\n")
      if process.stdout._type? and process.stdout._type isnt 'tty'
         console.warn("NOTE: You appear to be redirecting the output of this session, no interactive command results will be visible.\n")
      
      repl.start({prompt:">> ", output:process.stderr, input:process.stdin, terminal:true, ignoreUndefined:true, eval:(l,cx,fn,cb) ->
         [cmd, argstr...] = l[1..-2].trim().replace(/\t|\n/," ").split(" ")
         argstr = argstr.join(' ').replace(/^['"]|["']$/g,"") 
         cmd = cmd.toUpperCase()
         if cmd is 'HELP'
            cmd_func = commands[cmd]
            if argstr?.match(/{.*}/)  # Does it look like JSON?
               try
                  arg_json = JSON.parse(argstr)
               catch err  # JSON didn't parse into a valid object, so die
                  console.error("\nERROR: Could not parse JSON arguments '#{argstr}' to '#{cmd}' command\n")
                  cb(null, undefined)
                  return
               setImmediate(cmd_func, j, arg_json, (j) -> cb(null, undefined))
            else
               setImmediate(cmd_func, j, {'topic':argstr}, (j) -> cb(null, undefined))
            
         else if typeof(commands[cmd]) is 'function'

            if argstr        
               try
                  args = JSON.parse(argstr)
               catch err
                  console.warn("\nInvalid arguments '#{argstr}' to '#{cmd}'. Could not parse JSON string.\n")
                  cb(null, undefined)
                  return
            else
               args = null
            
            if (cmd is 'SCRIPT') and not args?.file?
               console.warn("\nInvalid command: You may not launch nested interactive SCRIPT sessions.\n")
               cb(null, undefined)
            else unless (j?.nodes? or cmd is 'LOAD' or cmd is 'SCRIPT' or cmd is 'UNSTASH') or cmd is 'HELP'
               console.warn("\nERROR: No data is currently loaded. Use the LOAD command to read a datafile\n")
               cb(null, undefined)
            else
               cmd_func = commands[cmd]
               setImmediate(cmd_func, j, args, (err, return_j) ->
                  if err 
                     console.error("ERROR: #{err.message}\n")
                     console.warn("\nCommand failed: Please examine the above error messages and try again.\n")
                  else unless return_j
                     console.warn("\nWARNING: Previous command returned no output.\n")
                     j = {"processing" : (j?.processing or [])}  # Maintain the command history
                  else
                     j = return_j
                     j?.processing?.push([">>", ss_version, cmd, clone_object(args)])
                  cb(null, undefined)
               )
         else 
            console.warn("\nInvalid command: '#{cmd}'\n") if cmd
            cb(null, undefined)
         
      }).on('exit',() -> 
               console.warn("\n")
               callback?(null, j))

###
# output_help -- Implementation of the HELP command
###
output_help = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- List all valid commands, or provide detailed help for a specific command\n
        with #{args.help} <COMMAND>
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
topic : \"command\" -- Provide detailed help for a specific command.\n
\n
        Example: #{args.help} {\"topic\":\"#{args.help}\"} -- Provide help about the\n
        #{args.help} command itself\n
\n
        Using the special word 'ALL' (which is not a valid command) in place of the\n
        command string will cause detailed help to be printed for all commands.\n
\n
        Example: #{args.help} {\"topic\":\"ALL\"} -- Provide detail help about commands\n
\n
        NOTE: The #{args.help} command when used as the first and only command on the\n
        command-line (or within an interactive SCRIPT session) may include the topic\n
        command directly as the second command line parameter\n
\n
        Example: #{args.help} #{args.help} is equivalent to #{args.help} {\"topic\":\"#{args.help}\"} in this situation\n
\n
cmd_line : true -- Provide detailed help for UNIX shell command line use.\n
\n
        Example: #{args.help} {\"cmd_line\":\"true\"} -- UNIX command line help.\n
\n
        This help can also be obtained from the shell by running with the -h option.\n
")
      callback?(null, j)
      return

   cmd_list = () ->
         for c of commands
            if typeof(commands[c]) is 'function'
               commands[c](j,{'help':c})
            else
               console.warn(commands[c])
         console.warn('')

   if args.topic?
      args.topic = args.topic.toUpperCase()
      if typeof(commands[args.topic]) is 'function'
         console.warn('')
         commands[args.topic](j,{'help':args.topic,'detailed_help':true})
      else if args.topic.toUpperCase() is 'ALL'
         cmd_list()
         for c of commands
            if typeof(commands[c]) is 'function'
               console.warn('==========================================================================================\n')        
               commands[c](j,{'help':c,'detailed_help':true})
      else 
         console.warn("\nCommand '#{args.topic}' not found, please select a valid command from below:\n")
         cmd_list()
   else if args.cmd_line?  
         console.warn("
#{commands['prelude1']}
\n
Usage: graph_ops [<input.json[.gz]>] [<script.go[.gz]>] [<command> ['{params}']]...\n
\n
where: <input.json[.gz]> is an optional datafile to initially LOAD\n
       <script.go[.gz]> is an optional SCRIPT to initially execute\n
       <command> is a valid command (see below)\n
       '{params}' optionally specify parameters for a given command\n
\n
Multiple commands with optional parameters may be provided in succession for execution.\n
\n
Running this program without any options will launch an interactive mode,\n
equivalent to:  graph_ops SCRIPT\n  
\n
For a list of valid commands: use: graph_ops HELP\n
\n
For detailed help for a specific command, use: graph_ops HELP <command>\n
\n
SEAStAR Version: #{ss_version}\n
")
      else
         cmd_list()
         
   callback?(null, j)

###
# push_stash -- 
###
push_stash = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Copy (push) the current graph to the top of a stack in memory
")
      if args.detailed_help?
         console.warn("
\n
Parameters:  None.\n
")
      callback?(null, j)
      return

   stash_stack.push(clone_object(j))

   console.warn("#{stash_stack.length} graph#{if stash_stack.length isnt 1 then 's' else ''} in stash.")

   callback?(null, j) 

###
# pop_stash -- 
###
pop_stash = (j, args={}, callback) ->
# j = json graph object, args = arguments from command, callback = function to call when done

   if args.help?
      console.warn("
#{args.help} -- Restore (pop) the current graph from the top of the stack (undo changes) 
")
      if args.detailed_help?
         console.warn("
\n
Parameters:\n
\n
free : true -- Free the graph at the top of the stack without restoring its state.\n
\n
        Example: #{args.help} {\"free\":true} -- Discard most recently STASHed graph\n
\n
")
      callback?(null, j)
      return

   unless stash_stack.length > 0
      callback?(new Error("UNSTASH failed, stash is empty."), null)
      return
   else

      if args.free
         stash_stack.pop()  # Don't restore state
      else 
         hist = j.processing or []  # Preserve history across UNSTASH
         j = undefined
         j = stash_stack.pop()
         j.processing = hist  # Restore history
      
      console.warn("#{stash_stack.length} graph#{if stash_stack.length isnt 1 then 's' else ''} remain#{if stash_stack.length is 1 then 's' else ''} in stash.")
      callback?(null, j) 

###
# ************************************************************************
# Mapping of command strings to implememnting functions defined above
#
# All functions below implement the same calling interface:
#
# function(j, args, callback), where:
# 
# j        :  the JSON graph object to work on
# args     :  an object of arguments to this command  
# callback :  a function to call when the operation is complete
#
###

# Lower case attributes below can't be executed
commands = 
   'prelude1'  : '\n
This program reads json format data files produced by the SEAStAR ref_select program.\n
It implements a variety of commands for manipulating this data for assembly,\n
visualization or output to a variety of file formats.\n'
   'prelude2' : '
For help with UNIX command line options, run the program with -h or --help\n\n
All commands below accept parameters in the form: {"parm1":value1,"parm2":value2...}\n
For examples and/or detailed help with specific commands, type:  HELP <command>\n
'
   'spacer1'  : '==============================='
   'text1'    : 'File I/O commands\n'
   'LOAD'     : read_json
   'DUMP'     : write_json
   'TABLE'    : export_table
   'FASTA'    : export_fasta
   'DOT'      : export_dot
   'spacer2'  : '\n==============================='
   'text2'    : 'Assembly pipeline commands (in this order)\n'
   'MST'      : maximal_spanning_tree
   'SST'      : scaffold_spanning_tree
   'text3'    : '   *** SST is an optional replacement for MST, useful for metagenomes'
   'text3a'   : '   *** or relatively short contigs'
   'PLUCK'    : remove_leaves
   'PRUNE'    : cut_branches
   'SLICE'    : check_connections
   'text4'    : '   *** SLICE is optional, but recommended for metagenomes'
   'PUSH'     : add_ends
   'INSERT'   : full_order
   'SCAFF'    : scaffold
   'CLUST'    : tetracalc
   'SELCLUST' : grab_clusts
   'spacer3'  : '\n==============================='
   'text6'    : 'Assembly graph filter/edit utilities\n'
   'CCOMPS'   : calc_ccomps
   'SELCC'    : grab_ccomps
   'EDGFLT'   : filter_edges
   'RELINK'   : relink
   'SELND'    : grab_neighbors
   'SCAFLNK'  : scaff_link
   'PROBS'    : node_problems
   'EDIT'     : perform_edits
   'CUTND'    : cut_node
   'spacer4'  : '\n==============================='
   'text7'    : 'Information and control\n'
   'GC'       : graph_stats
   'GCC'      : graph_stats_cc
   'STASH'    : push_stash
   'UNSTASH'  : pop_stash
   'SCRIPT'   : command_file
   'HELP'     : output_help

###   
# This function gets called once the entire JSON data structure is available
# It parses the commandline options and builds a callback_list that is used
# to chain command calls together in an 'event driven' friendly manner.
# It ends by calling the first command on the list via the call_next_cmd function.
###
process_commands = (callback_list) ->
   json = {}

   # This function is used to pass control to the next command as a callback 
   call_next_cmd = (err, j) ->
      
      if err
         console.error("ERROR: #{err.message}\n")
         process.exit(1)
      
      if (callback_list.length)  # If there are commands remaining to process
         [cmd_func, cmd_args, cmd] = callback_list.shift()  # Prepare the next command
         unless j?.nodes? or cmd is 'LOAD' or cmd is 'SCRIPT' or cmd is 'UNSTASH' or cmd is 'HELP' 
            console.error("ERROR: No data, ending command line processing.\n")
            process.exit(1)
            
         # Record history of this command in the graph object
         j?.processing?.push(["$", ss_version, cmd, clone_object(cmd_args)]) unless cmd is 'HELP'  
         console.warn("Executing #{cmd} #{JSON.stringify(cmd_args)}")  # Output progress to stderr
         setImmediate(cmd_func, j, cmd_args, call_next_cmd)  # Invoke the command with the graph object, args, and callback function

   # Make the initial call with an empty graph object
   # The callback list is available via closure in the call_next_cmd function above
   call_next_cmd(null, json)
      
###
# Handle commandline args
###

cmds = []  # Initialize the command list
callback_list = [] # Initialize the callback_list

# If no parameters, drop into interactive SCRIPT mode
if process.argv.length < 3  
      callback_list.push([commands['SCRIPT'], {}, 'SCRIPT'])
      
# Special case for HELP (or some common variant thereof) is first parameter on command line
else if (process.argv[2].toUpperCase() is 'HELP')
   if process.argv[3]?
      if process.argv[3].match(/{.*}/)
         callback_list.push([commands['HELP'], JSON.parse(process.argv[3]), 'HELP'])
      else
         callback_list.push([commands['HELP'], {'topic':process.argv[3]}, 'HELP'])
   else
      callback_list.push([commands['HELP'], {}, 'HELP'])
else if (process.argv[2].match(/^-?-h(elp)?$/i))
      callback_list.push([commands['HELP'], {'cmd_line':true}, 'HELP'])
else   # In this case, the first argument is a filename, which triggers an implicit LOAD command
   start_arg = 2  # First arg is assumed to be a command, but maybe its a file to LOAD or SCRIPT
   unless (process.argv[2].toUpperCase() is 'LOAD') or (process.argv[2].toUpperCase() is 'SCRIPT')
      if process.argv[2].match(/^(-|(.+\.json(\.gz)?))$/)   # Does this look like a data file?
         callback_list.push([commands['LOAD'], {'file':process.argv[2]}, 'LOAD'])
         start_arg++  # First arg is a JSON filename, so skip it below
         unless process.argv[start_arg]?  # If there are no explicit commands, go interactive
            process.argv[start_arg] = 'SCRIPT'
      if process.argv[start_arg].match(/.+\.go(\.gz)?$/) # Does this look like a SCRIPT file?
         callback_list.push([commands['SCRIPT'], {'file':process.argv[start_arg]}, 'SCRIPT'])
         start_arg++  # Next (maybe still first) arg is a SCRIPT filename, so skip it below
      
   prev_cmd = false   
   # Loop through the command line arguments, building up a command string to pass to the control loop   
   for cmd in process.argv[start_arg..-1]
      if prev_cmd and cmd.match(/{.*}/)    # If this looks like a JSON string, it's a parameter, so add it to the previous command
         try
            callback_list[-1..][0][1] = JSON.parse(cmd)
         catch err  # JSON didn't parse into a valid object, so die
            console.error("ERROR: Could not parse JSON arguments '#{cmd}' to #{callback_list[-1..][0][2]} command")
            process.exit(1)
         prev_cmd = false
      else if typeof(commands[cmd.toUpperCase()]) is 'function' 
         # If this is a valid command string, then add 
         cmd = cmd.toUpperCase() 
         callback_list.push([commands[cmd], null, cmd])
         prev_cmd = true
      else # Don't know what this is, so die.
         console.error("ERROR: Invalid command line command: '#{cmd}'\n")
         process.exit(1)
   
# Start running commands
process_commands(callback_list)
