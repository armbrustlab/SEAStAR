/*
# --------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2013 University of Washington.
#
# Authors:
# Vaughn Iverson
# vsi@uw.edu
# --------------------------------------------------------------------------- #
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
# along with SEAStAR.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------------- #
*/
window.onload = function () {

// From geometry.js -- Example 14-2 from JavaScript: The Definitive Guide, Fifth Edition
var Geometry = {};
if (window.innerWidth) { // All browsers but IE
    Geometry.getViewportWidth = function() { return window.innerWidth; };
    Geometry.getViewportHeight = function() { return window.innerHeight; };
    Geometry.getHorizontalScroll = function() { return window.pageXOffset; };
    Geometry.getVerticalScroll = function() { return window.pageYOffset; };
}
else if (document.documentElement && document.documentElement.clientWidth) {
    // These functions are for IE6 when there is a DOCTYPE
    Geometry.getViewportWidth =
        function() { return document.documentElement.clientWidth; };
    Geometry.getViewportHeight = 
        function() { return document.documentElement.clientHeight; };
    Geometry.getHorizontalScroll = 
        function() { return document.documentElement.scrollLeft; };
    Geometry.getVerticalScroll = 
        function() { return document.documentElement.scrollTop; };
}
else if (document.body.clientWidth) {
    // These are for IE4, IE5, and IE6 without a DOCTYPE
    Geometry.getViewportWidth =
        function() { return document.body.clientWidth; };
    Geometry.getViewportHeight =
        function() { return document.body.clientHeight; };
    Geometry.getHorizontalScroll =
        function() { return document.body.scrollLeft; };
    Geometry.getVerticalScroll = 
        function() { return document.body.scrollTop; };
};

var	r = Raphael(0, 95, Geometry.getViewportWidth(), Geometry.getViewportHeight()-95),
	w_inc = 120*(Geometry.getViewportWidth()/1800.0),
	major_gap = 35*(Geometry.getViewportWidth()/1800.0),
	minor_gap = 7*(Geometry.getViewportWidth()/1800.0),
	w_gap = -1*(Geometry.getViewportWidth()/1800.0),
	x_offset = 30*(Geometry.getViewportWidth()/1800.0),
	y_offset = 75,
	font_param = 9.0*(Geometry.getViewportWidth()/1800.0),
	hide_level = 1.0,
	fig_h = (Geometry.getViewportHeight()-220.0),
	code = document.getElementById("code"), 
	download = document.getElementById("download_link"), 
	checkbox = document.getElementById("percent"),
	level_f = document.getElementById("level_form"),
	samp_f = document.getElementById("sample_form"),
	show_levels = [], 
	show_samples = [],
	num_samples = 0,
	num_levels = 0,
	num_samp_in_data = tree_struct.sample_names.length,
	serializer = new XMLSerializer();
	
function redraw()
{
 r.remove();
 window.onload();
}	
	
function init_levels(f,show_l) {
	for (var l = 0; l < f.elements.length; l++) {
		var c = f.elements[l];
		c.onclick = redraw;
		show_l[l] = c.checked;
		if (c.checked) { num_levels++; }
	}
}

function init_samples(f,show_s) {
	for (var l = 0; l < f.elements.length; l++) {
		var c = f.elements[l];
		if (l >= num_samp_in_data) {
			c.checked = false;
			show_s[l] = false;
		} else {
			c.onclick = redraw;
			show_s[l] = c.checked;
			if (c.checked) { num_samples++; }
		}
	}
}
	
init_levels(level_f, show_levels);	
init_samples(samp_f, show_samples);	

if (num_samples == 1) {
	major_gap = major_gap / 3.0; 
	minor_gap = 0;
}

window.onresize = redraw;	
	
var s_tree = null;

show_percent = checkbox.checked;

/**
 * Draw a bar chart into an <svg> element.
 * Arguments:
 */
function barWedge(R, s, startval, val, total, x1, y, h, w, n_samp, lev, color, border, label, node, rot) {
	    
    var start = startval/total;
	var thick = val/total;   
    var end = start + thick;

//    var x1 = x + Math.ceil(lev)*n_samp*w;

    var y1 = y + h*start;
    var w1 = w;
    var h1 = h*thick;

	if (Math.floor(lev) != lev) {
		w1 = w/2;
	}
    
    var tx = x1 + w1/2.0;
    var ty = y1 + h1/2.0;

   var n = val * 100.0 / total;
   
   //var n = val;

	if (label != "") {
		var wedge = R.rect(x1, y1, w1, h1).attr({ fill : color, stroke : border, "stroke-width" : "1"}).toBack();
		node.wedge = wedge;
		s.push(wedge);	

		var txt_col = "black";
		if (node.w_conf < 0.50) {
			txt_col = "white";
		}

//  Shows percentages too		
		var txt;
		if (show_percent) {
			txt = R.text(tx, ty, label+" ("+n.toFixed(1)+"%)").attr({font: font_param+'px "Helvetica"', stroke : "none", fill : txt_col}).toFront();
		} else {
			txt = R.text(tx, ty, label).attr({font: font_param+'px "Helvetica"', stroke : "none", fill : txt_col}).toFront();
		}
		
		if (Math.floor(lev) != lev) {
			txt.hide();
		}
		
		if (rot) {
			txt.rotate(-90);
		}
		
//   		if (n < 1.0) { 
   		if (0) { 
   			txt.hide();
			txt.mouseover((function () {var t = txt; return function () { t.show(); }; })()).mouseout((function () {var t = txt; return function () { t.hide(); }; })());
			wedge.mouseover((function () {var t = txt; return function () { t.show(); }; })()).mouseout((function () {var t = txt; return function () { t.hide(); }; })());
   		}
	
		s.push(txt);
		node.txt = txt;
		
		if (node.sub.length == 0) {
			wedge.dblclick((function () {var l = node.label; return function () { window.open("http://greengenes.lbl.gov/cgi-bin/show_one_record_v2.pl?prokMSA_id=" + l); }; })());
			txt.dblclick((function () {var l = node.label; return function () { window.open("http://greengenes.lbl.gov/cgi-bin/show_one_record_v2.pl?prokMSA_id=" + l); }; })());
		} else {
			// txt.click((function () { var n = node; return function () { sideTree(node); }; })());
			// wedge.click((function () { var n = node; return function () { sideTree(node); }; })());
		}

	} else {
		var wedge = R.rect(x1, y1, w1, h1).attr({stroke : border, "stroke-width" : "3"}).toFront();
 		s.push(wedge);
	}
	return wedge;
};

var tot = tree_struct.sub.Root.pop;
var cur_s_tree = null;

var main_set = r.set();
var side_set = null;
var bottom_legend = ["domain", "phylum", "class", "order", "family", "genus / environmental clade"];

var greek = { "α-" : /Alpha/, "β-" : /Beta/ , "γ-" : /Gamma/, "δ-" : /Delta/, "ε-" : /Epsilon/ };
var flip = false;
var last = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
var prev_lev = [x_offset,x_offset,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
var lev_width = [0.0,0.25,1.0,1.0,1.0,1.0,1.0,1.0,0.0];
var l_cnt = [0,0,0,0,0,0,0,0,0];
var b_cnt = [0,0,0,0,0,0,0,0,0];
var inc = 0;
var s_inc = 0;
var l_inc = 0;

function walk_tree(node, s_inc, y, s, t, sample, rev, lev_inc, name) {
	if (node.level < 7) {
		var b = "black";
		// var color = Raphael.hsl2rgb(1.0-((node.cum+node.pop/2.0)/tot), 0.2+0.8*(node.level)/7.0, 0.5+0.2*node.w_conf).hex;
		var color = Raphael.hsl2rgb(1.0-((node.cum+node.pop/2.0)/tot), 0.1+0.9*(node.level)/7.0, 0.75).hex;
		node.label = name;
		for (var letter in greek) {
			node.label = node.label.replace(greek[letter],letter);
		}

		var node_pop = node.samples[sample].pop;

		if (node.hasOwnProperty("sub")) {
			for (var z in node.sub) { if (node.sub[z].samples[sample]) { walk_tree(node.sub[z], s_inc, y, s, t, sample, rev, (show_levels[node.level-1] ? lev_inc+1 : lev_inc), z); } }
			
			if ((node.level > 0) && (node_pop >= hide_level)) {
				if (show_levels[node.level-1]) {
					barWedge(r, s, last[node.level], node_pop, t, s_inc*(minor_gap+w_inc*lev_width[node.level]) + prev_lev[lev_inc], y, fig_h, w_inc*lev_width[node.level], num_samples, lev_inc, color, b, node.label, node, (node.level == 1)); 

				}
				last[node.level] += node_pop;

			} else {
				l_cnt[node.level] += 1;
			};
			
			if (node.level == 3) {
				if (node_pop >= hide_level) {
					inc = show_levels[node.level - 1] ? lev_inc + 1 : lev_inc;				
					// inc = lev_inc;				
					for (var l = node.level+1; l < 7; l++) { 
						if (show_levels[l-1]) {
							if (l_cnt[l] != 0) {
								barWedge(r, s, last[l], last[node.level] - last[l], t, s_inc*(minor_gap+w_inc*lev_width[l]) + prev_lev[inc], y, fig_h, w_inc*lev_width[l], num_samples, inc, "#c0c0c0", b,  l_cnt[l] +" more "+node.label.replace(/bacteria/,"-"), node, 0);
								last[l] = last[node.level]; 
								l_cnt[l] = 0;
							}
							inc++;
						}
					}
				} else {
					for (var l = 2; l < 7; l++) { 
						b_cnt[l] += l_cnt[l];
						l_cnt[l] = 0;
					}
				}
			}
			
			if (node.level == 1) {
				inc = show_levels[node.level - 1] ? lev_inc + 1 : lev_inc;
				for (var l = node.level+1; l < 7; l++) { 
					if (show_levels[l-1]) {
						if (b_cnt[l] != 0) {
							barWedge(r, s, last[l], last[node.level] - last[l], t, s_inc*(minor_gap+w_inc*lev_width[l]) + prev_lev[inc], y, fig_h, w_inc*lev_width[l], num_samples, inc, "#f0f0f0", b,"other "+node.label, node, 0);
							last[l] = last[node.level]; 
							b_cnt[l] = 0;
						}
						inc++;
					}
				}
			}
		}
	}
};

var prev_l = 0;

for (var lev = 1; lev <= bottom_legend.length; lev++) {
		
		if (show_levels[lev-1]) {
			l_inc++;
			if (l_inc >= 2) {
				prev_lev[l_inc] = prev_lev[l_inc-1] + major_gap + (num_samples-1)*minor_gap + num_samples*w_inc*lev_width[prev_l];
			}
			prev_l = lev;
		}
}

l_inc++;
prev_lev[l_inc] = prev_lev[l_inc-1] + major_gap + (num_samples-1)*minor_gap + num_samples*w_inc*lev_width[prev_l];


s_inc = 0;
for (var samp = 0; samp < num_samp_in_data; samp++) {
	if (show_samples[num_samp_in_data - samp - 1]) {
		var inv_samp = num_samples-s_inc-1
		walk_tree(tree_struct.sub.Root, inv_samp, y_offset, main_set, tree_struct.sub.Root.samples[num_samp_in_data-samp-1].pop, num_samp_in_data-samp-1, flip, 1, "Root");
		last = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
		l_cnt = [0,0,0,0,0,0,0,0,0];
		s_inc++;
	}
}

// var bar_pos_x = x_offset+num_levels*num_samples*w_inc+(num_levels-1)*major_gap+(num_samples-1)*minor_gap + w_inc/5.0;
var bar_pos_x = prev_lev[l_inc] - major_gap + w_inc/10.0;


if (!show_percent) {
	r.path("M{0},{1}L{2},{3}", bar_pos_x, y_offset+fig_h*0.1, bar_pos_x, y_offset+fig_h*0.15).attr({ stroke : "black", "stroke-width" : 0.3*font_param+"px" });

	r.text(bar_pos_x+w_inc/7.5, y_offset+fig_h*0.125, "5%").attr({font: 1.2*font_param+"px 'Helvetica'", stroke : "none", fill : "#000000"});
}

if (num_samples) {
	inc = 1;
	l_inc = 0;
	for (var x in bottom_legend) {
		l_inc++;
		if (show_levels[x]) {
			r.text((prev_lev[inc] + prev_lev[inc+1] - major_gap)/2, y_offset*0.75, bottom_legend[x]).attr({font: 1.4*font_param+"px 'Helvetica'", stroke : "none", fill : "#000000"}).toFront();

			s_inc = 0;
			for (var samp = 0; samp < num_samp_in_data; samp++) {
				if (show_samples[samp]) {
					r.text(lev_width[l_inc]*w_inc/2 + s_inc*(minor_gap + lev_width[l_inc]*w_inc) + prev_lev[inc], fig_h+1.3*y_offset, tree_struct.sample_names[samp]).attr({font: 1.4*font_param+"px 'Helvetica'", stroke : "none", fill : "#000000"}).toFront();
					s_inc++;
				}
			}
			inc++;
		}
	}
}

download.onmouseover = function () {
	var svg_xml = serializer.serializeToString(r.canvas);
	console.log(svg_xml);
	this.setAttribute('href-lang', "image/svg+xml");
	// btoa with utf-8 strings can cause Character Out Of Range exception
	// See https://developer.mozilla.org/en-US/docs/DOM/window.btoa#Unicode_Strings for source of this solution
	this.setAttribute('href', 'data:image/svg+xml;base64,\n' + btoa(unescape(encodeURIComponent(svg_xml))));
};

checkbox.onclick = redraw;

};
