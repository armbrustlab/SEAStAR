// -------------------------------------------------------------------------- //
// Center for Environmental Genomics
// Copyright (C) 2009-2012 University of Washington.
//
// Authors:
// Vaughn Iverson
// vsi@uw.edu
// -------------------------------------------------------------------------- //
// This file is part of SEAStAR.
//
// SEAStAR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SEAStAR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SEAStAR.  If not, see <http://www.gnu.org/licenses/>.
// -------------------------------------------------------------------------- //

// This script reads in an RDP classifier taxonomic heirarchy file and converts
// it to a JSON tree formatted output file.

// Inputs: input modified RDP classifier output on stdin 
//         JSON formatted RDP heirarchy in the file named in RDP_heir_fn (or RDP_expand.json by default)

// Output: output JSON formatted heirarchy to stdout

// Output file format:
// Each taxonomic unit in the RDP input file becomes a branch in the output JSON tree.
// Each branch has the following JSON structure:
//
// name : {
//      "pop" : <fraction of the sample>,
//      "cum" : <cumulative fraction of the sample>,
//      "num" : <number of sequences>,
//      "cnt" : <cumulative number of seqeunces>,
//      "conf" : <mean classifier p-value for sequences in this taxon>,
//      "w_conf" : <population weighted mean classifier p-value for sequences in this taxon>,
//      "level" : <numeric level in heirarchy>,
//      "length" : <branch length of this taxonomic level from parent>
//      "sub" : {<child nodes by name, or empty if a leaf>}
// }

// Handle default filename.  
if (typeof RDP_heir_fn == 'undefined') {
	RDP_heir_fn = 'RDP_expand.json';
}

// Read and parse the JSON RDP taxonomy
var out_tree = JSON.parse(read(RDP_heir_fn));

// Method for branch objects that calculates summary stats for conf, cnt and cum
// And prunes branches from the tree with no representative sequences.
out_tree.walk = function (tree, prev_cnt, prev_cum) {
 
	for (var c in tree.sub) {

	        // Shuffle sort child elements by pop (largest, smallest, second largest, second smallest, etc.)
               	var x = 0;
               	var tmp = [];
               	// Create an array so we can reverse sort it.
               	for (var y in tree.sub) {
                       	tmp[x] = tree.sub[y];
//                       	tmp[x].name = y;
                       	x++;
               	}
               	tmp.sort(function (b,a) { return a.pop-b.pop; } );
       		// Now replace the original sub object
       		tree.sub = {};
               	// This performs the "shuffle" on the sorted array
               	for (x = 0; x < Math.round(tmp.length / 2); x++) {
               		tree.sub[tmp[x].name] = tmp[x];
//                       	delete tmp[x].name;
                       	var z = tmp.length - x - 1;
                       	if (z != x) {
                       		tree.sub[tmp[z].name] = tmp[z];
//                      	        	delete tmp[z].name;
               	        }
                }
		break;	// No need to loop, this was just to test that sub had something in it.
	}

	for (var c in tree.sub) {
		if (tree.sub[c].num != 0) {
			var sub_tree = tree.sub[c];
			sub_tree.conf = sub_tree.conf / sub_tree.num;
			sub_tree.w_conf = sub_tree.w_conf / sub_tree.pop;
			sub_tree.cnt = prev_cnt;
			sub_tree.cum = prev_cum;
			out_tree.walk(sub_tree, prev_cnt, prev_cum);	
			prev_cnt += sub_tree.num;
			prev_cum += sub_tree.pop;
		} else {
			delete tree.sub[c];
		}
	}
}
	
var line = "";

// Walk all of the lines in the RDP classifier output 

//while (line = readln()) {
while (line = readline()) {

        line = line.trim();     // Remove whitespace from ends

        if (line) {     // If not blank
		var parts = line.split(";");    // Split on semicolons
		var firstpart = parts[0].split("_"); // Split first field on "_"
		var sequence = firstpart[0];		// This is the sequence name
		var percent = parseFloat(firstpart[1]);	// This is its estimated abundance
		var tax_levels = [];	// These arrays will accept parsed taxonomy 
		var tax_pvals = [];	// and p-value data
		
		// Loop over the taxon levels and initialize the above arrays
		for (var x = 2; x < parts.length; x+=2) {
			tax_levels[x/2-1] = parts[x].replace(/"/g,"").trim();
			tax_pvals[x/2-1] = parseFloat(parts[x+1]);
		}

		var cur = out_tree;
		var lev = 0;

		// Walk down the branches of the RDP taxonomy tree and fill in the 
		// information for this sequence.

		for (lev = 0; lev < tax_levels.length; lev++) {
			cur = cur.sub[tax_levels[lev]];
			cur.num++;
			cur.conf += tax_pvals[lev];
			cur.w_conf += tax_pvals[lev] * percent;
			cur.pop += percent;
		}

		// Now add the sequence itself as a leaf of the lowest taxon level (usually genus, level 6)		
		cur.sub[sequence] = { 
					pop : percent,
                                        cum : 0.0,
                                        cnt : 0,
                                        num : 1,
                                        conf : tax_pvals[tax_pvals.length-1], 
					w_conf : percent * tax_pvals[tax_pvals.length-1],
                                        level : 7.0,
                                        length : 1.0,
                                        name : sequence
//                                        sub : {},
				    };
	}
}

// Walk the tree, pruning empty branches and updating cumulative statistics

out_tree.walk(out_tree, 0, 0.0);

// Output JSON to stdout

print(JSON.stringify(out_tree, null, 1));

quit(0);

