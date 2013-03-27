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

// Input: input RDP heirarchy on stdin 
// Output:  output JSON formatted heirarchy to stdout

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
//      "level" : <numeric level in heirarchy>,
//      "length" : <branch length of this taxonomic level from parent>
//      "sub" : {<child nodes by name, or empty if a leaf>}
// }

var taxa = []; // This is a working array of taxa records indexed by taxid

var levels = { 	"norank" : 0, "domain" : 1, "phylum" : 2, "class" : 3, 
           	"subclass" : 3.5, "order" : 4, "suborder" : 4.5, 
           	"family" : 5, "subfamily" : 5.5, "supergenus" : 5.75, 
           	"genus" : 6 };

var error_codes = { "INVALID_INPUT_LINE" : 1 };

// Now loop through lines, building the taxa table and making the heirarchy 

var line = "";

//while (line = readln()) {
while (line = readline()) {
	
	line = line.trim();	// Remove whitespace from ends
	
	if (line) {	// If not blank
		var parts = line.split("*");	// Split on asterisks
		if (parts.length != 5) {	// There must be 5 columns
			quit(error_codes["INVALID_INPUT_LINE"]);
		}
		var taxid = parseInt(parts[0]);		// taxid and parentid 
		var parentid = parseInt(parts[2]);	// are each integers
		var name = parts[1].replace(/"/g,"");	// remove quotes from taxa names
		parts[4] = parts[4].trim();		// be tolerant of whitespace

		if (parts[4] in levels) {	// If this taxon has already been used 
						// as a parent (it is out of order)
			if (taxid in taxa) {	// then copy the existing sub taxa 
				taxa[taxid].level = levels[parts[4]];  
				// Go back and fill in the branch length value for child taxa	
				for (var i in taxa[taxid].sub) {
					taxa[taxid].sub[i].length = taxa[taxid].sub[i].level - taxa[taxid].level;
				}

			} else {		// otherwise make a whole new object
				taxa[taxid] = {	pop : 0.0,
						cum : 0.0, 
						cnt : 0,
						num : 0,
						conf : 0.0, 
						w_conf : 0.0, 
						level : levels[parts[4]],
						length : 0.0,
						sub : {} 
						};
			}
			// Look to see if this taxon's parent already has a record
			if (parentid in taxa) {
				taxa[parentid].sub[name] = taxa[taxid];	// Add this child
				taxa[taxid].length = taxa[taxid].level - taxa[parentid].level;
			} else if (parentid >= 0) {			// If not the root
				taxa[parentid] = { 
						pop : 0.0,
                                                cum : 0.0, 
                                                cnt : 0,
                                                num : 0,
                                                conf : 0.0,
                                                w_conf : 0.0,
						level : 0.0,
						length : 0.0,
						sub : {} 
						};	// Add a placeholder
				taxa[parentid].sub[name] = taxa[taxid];
			}
		}
	}
}

// Write the JSON structure to stdout.
taxa[0].length = 1; 
print(JSON.stringify({ sub : { Root : taxa[0] } }));

quit(0);

