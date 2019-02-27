// Generated by CoffeeScript 1.12.7

/*
 * -------------------------------------------------------------------------- #
 * Center for Environmental Genomics
 * Copyright (C) 2009-2013 University of Washington.
 *
 * Authors:
 * Vaughn Iverson
 * vsi@uw.edu
 * -------------------------------------------------------------------------- #
 * This file is part of SEAStAR.
 *
 * SEAStAR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SEAStAR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SEAStAR.  If not, see <http:#www.gnu.org/licenses/>.
 * -------------------------------------------------------------------------- #
 */

(function() {
  var build, error_codes, levels, taxa, ver;

  ver = process.version.slice(1).split('.');

  if (!(ver[1] >= 10 || ver[0] > 0)) {
    console.error("ERROR: nodejs version v0.10.0 or greater required.");
    process.exit(1);
  }

  taxa = [];

  levels = {
    rootrank: 0,
    norank: 0,
    domain: 1,
    phylum: 2,
    "class": 3,
    subclass: 3.5,
    order: 4,
    suborder: 4.5,
    family: 5,
    subfamily: 5.5,
    supergenus: 5.75,
    genus: 6
  };

  error_codes = {
    INVALID_INPUT_LINE: 1
  };

  process.on('uncaughtException', function(err) {
    return console.log('Caught exception: ' + err);
  });

  process.stdin.resume();

  process.stdin.setEncoding('utf8');

  build = function(line) {
    var child, i, level, levnum, name, parentid, ref, ref1, taxid;
    if (!(line = line.trim())) {
      return;
    }
    if ((ref = line.split("*"), taxid = ref[0], name = ref[1], parentid = ref[2], levnum = ref[3], level = ref[4], ref).length !== 5) {
      return;
    }
    taxid = parseInt(taxid);
    parentid = parseInt(parentid);
    name = name.replace(/"/g, "");
    level = level.trim();
    if (!(level in levels)) {
      return;
    }
    if (taxa[taxid] != null) {
      taxa[taxid].level = levels[level];
      taxa[taxid].name = name;
      ref1 = taxa[taxid].sub;
      for (i in ref1) {
        child = ref1[i];
        child.length = child.level - taxa[taxid].level;
      }
    } else {
      taxa[taxid] = {
        name: name,
        pop: 0.0,
        cum: 0.0,
        cnt: 0,
        num: 0,
        conf: 0.0,
        w_conf: 0.0,
        level: levels[level],
        length: 0.0,
        sub: {}
      };
    }
    if (taxa[parentid] != null) {
      taxa[parentid].sub[name] = taxa[taxid];
      taxa[taxid].length = taxa[taxid].level - taxa[parentid].level;
    } else {
      if (parentid >= 0) {
        taxa[parentid] = {
          pop: 0.0,
          cum: 0.0,
          cnt: 0,
          num: 0,
          conf: 0.0,
          w_conf: 0.0,
          level: levels[level],
          length: 0.0,
          sub: {}
        };
        taxa[parentid].sub[name] = taxa[taxid];
      }
    }
  };

  process.stdin.on('data', (function() {
    var save;
    save = '';
    return function(c) {
      var i, j, len, lines, results;
      lines = c.split('\n');
      lines[0] = save + lines[0];
      save = lines.pop();
      results = [];
      for (j = 0, len = lines.length; j < len; j++) {
        i = lines[j];
        if (i) {
          results.push(build(i));
        }
      }
      return results;
    };
  })());

  process.stdin.on('end', function() {
    taxa[0].length = 1;
    console.log(JSON.stringify({
      sub: {
        Root: taxa[0]
      }
    }));
    return process.exit(0);
  });

}).call(this);
