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

// This script reads a JSON formatted heirarchy file and converts it to NEXUS
// tree format

// Inputs: JSON formatted heirarchy file from stdin.

// Output: output NEXUS formatted tree to stdout

// Method for branch objects that calculates summary stats for conf, cnt and cum
// And prunes branches from the tree with no representative sequences.
walk = function (tree) {
	var outstr = "";

	if (tree.hasOwnProperty("sub")) {
		var children = [];
	        for (var c in tree.sub) {
			tree.sub[c].label = c;
                	children.push(walk(tree.sub[c]));
        	}
		outstr = "(" + String.concat(children) + ")" + "[&pop=" + Math.round(10000.0 * tree.pop/norm)/10000.0 + ",&tax=" + tree.label + ",&label=" + tree.label + "  " + Math.round(10000.0 * tree.pop/norm)/100.0 + "%]:" + tree.length;
	} else {
		outstr = "'" + tree.conf + "   " + tree.label + "   " + Math.round(10000.0 * tree.pop/norm)/100.0 + "%':" + tree.length;
	}

	return(outstr);
}

var input_JSON = "";
var line = "";
while (line = readline()) {
	input_JSON = input_JSON + line;
}

var input_tree = JSON.parse(input_JSON);

// Output NEXUS to stdout

input_tree.sub.Root.label = "Root";
var norm = input_tree.sub.Root.pop;
putstr("#NEXUS\n\n begin trees;\n tree tree_1 = [&R] " + walk(input_tree.sub.Root) + ";\nend;\n");

quit(0);

