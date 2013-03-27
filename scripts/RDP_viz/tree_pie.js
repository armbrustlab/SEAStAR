cur_s_tree = null;
cur_s_tree_name = "";

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

	var dim = ((Geometry.getViewportHeight()-50 < (Geometry.getViewportWidth()/1.5)) ? Geometry.getViewportHeight()-50 : (Geometry.getViewportWidth()/1.5)),
	cent = dim/2.0,
	R = 200.0*(dim/1000.0),
	rad_inc = 50.0*(dim/1000.0),
	w = 30.0*(dim/1000.0), 	
	r = Raphael(0, 50, 1.5*dim, dim),
	font_param = 8.0*(dim/1000.0), 
	code = document.getElementById("code"), 
	button = document.getElementById("SVGbutton"), 
	serializer = new XMLSerializer();
	
function resize()
{
 r.remove();
 window.onload();
}

window.onresize = resize;	
	
var s_tree = null;

/**
 * Draw a pie chart into an <svg> element.
 * Arguments:
 *   canvas: the SVG element (or the id of that element) to draw into.
 *   data: an array of numbers to chart, one for each wedge of the pie.
 *   cx, cy, r: the center and radius of the pie
 *   colors: an array of HTML color strings, one for each wedge
 *   labels: an array of labels to appear in the legend, one for each wedge
 *   lx, ly: the upper-left corner of the chart legend
 */
 
function pieWedge(R, s, startval, val, total, cx, cy, r, r2, color, border, label, node, side) {
    
    // Now figure out how big this slice of pie is.  Angles in radians.
    var angle = Math.PI; 

	if (startval == 0.0) {
		// This works around a drawing bug when startval == 0 
		startval = 0.0000001;
	}
    
    if (val != total) {
	    angle = val/total*Math.PI*2.0;
	}

    var start = startval/total*Math.PI*2.0;    
    
    var end = start + angle;

    // Compute the two points where our wedge intersects the circle
    // These formulas are chosen so that an angle of 0 is at 12 o'clock
    // and positive angles increase clockwise.
    var x1 = cx + r * Math.sin(start);
    var y1 = cy - r * Math.cos(start);
    var x2 = cx + r * Math.sin(end);
    var y2 = cy - r * Math.cos(end);

    // Compute the two points where our wedge intersects the circle
    // These formulas are chosen so that an angle of 0 is at 12 o'clock
    // and positive angles increase clockwise.
    var x1_2 = cx + r2 * Math.sin(start);
    var y1_2 = cy - r2 * Math.cos(start);
    var x2_2 = cx + r2 * Math.sin(end);
    var y2_2 = cy - r2 * Math.cos(end);

    // Compute the two points where our wedge intersects the outer circle
    // These formulas are chosen so that an angle of 0 is at 12 o'clock
    // and positive angles increase clockwise.
    var x1_3 = cx + rad_inc*7.0 * Math.sin(start);
    var y1_3 = cy - rad_inc*7.0 * Math.cos(start);
    var x2_3 = cx + rad_inc*7.0 * Math.sin(end);
    var y2_3 = cy - rad_inc*7.0 * Math.cos(end);

	// Compute text label location
    var tx = cx + (r+r2)/2.0 * Math.sin((start + end)/2.0);
	var ty = cy - (r+r2)/2.0 * Math.cos((start + end)/2.0);

    // This is a flag for angles larger than than a half circle
    var big = 0;
    if (end - start > Math.PI) big = 1;
    
    if (val != total) {
       var d = "M " + x1_2 + "," + y1_2 +  // Start at circle center
            " L " + x1 + "," + y1 +     // Draw line to (x1,y1)
            " A " + r + "," + r +       // Draw an arc of radius r
            " 0 " + big + " 1 " +       // Arc details...
            x2 + "," + y2 +             // Arc goes to to (x2,y2)
            " L " + x2_2 + "," + y2_2 + // Draw line to (x1,y1)
            " A " + r2 + "," + r2 +     // Draw an arc of radius r2
            " 1 " + big + " 0 " +       // Arc details...
            x1_2 + "," + y1_2;           // Arc goes to to (x2,y2)            
      } else {
       var d = "M " + x1_2 + "," + y1_2 +  // Start at circle center
            " A " + r2 + "," + r2 +       // Draw an arc of radius r
            " 0 " + big + " 1 " +       // Arc details...
            x2_2 + "," + y2_2 +         // Arc goes to to (x2,y2)
            " A " + r2 + "," + r2 +       // Draw an arc of radius r
            " 0 " + big + " 1 " +       // Arc details...
            x1_2 + "," + y1_2 +         // Arc goes to to (x2,y2)
            " M " + x1 + "," + y1 +     // Draw line to (x1,y1)
            " A " + r + "," + r +     // Draw an arc of radius r2
            " 1 " + big + " 0 " +       // Arc details...
            x2 + "," + y2 +           // Arc goes to to (x2,y2)    
            " A " + r + "," + r +     // Draw an arc of radius r2
            " 1 " + big + " 0 " +       // Arc details...
            x1 + "," + y1;           // Arc goes to to (x2,y2)    
      }

	if (label != "") {
		var wedge = R.path(d).attr({ fill : color, stroke : border, "stroke-width" : 0.1*font_param+"px"}).toBack();
 		node.wedge = wedge;
 		s.push(wedge);	
		
		var n = node.pop/tot * 100.0;
		var txt_angle = 0; 
		if (val != total) { 
			txt_angle = ((start+end)/2.0)*(360/(Math.PI*2));
			if ((txt_angle < 0.0) || (txt_angle >= 360.0)) {  txt_angle -= 360*Math.floor(txt_angle/360.0)  }
			if (node.level == 7) { txt_angle += ((txt_angle > 180) ? 90 : -90) }
			else { txt_angle += (((txt_angle > 90) && (txt_angle < 270)) ? -180 : 0) }
		}
		var txt_col = "black";
		if (node.w_conf < 0.50) {
			txt_col = "white";
		}
		var attr2 = {font: 1.25*font_param+"px 'Helvetica'", stroke : "none", fill : txt_col};
    		var txt = R.text(tx, ty, label).attr(attr2).rotate(txt_angle,true).toFront().hide();
    		txt.mouseover((function () {var t = txt; return function () { t.show(); }; })()).mouseout((function () {var t = txt; return function () { t.hide(); }; })());
 		s.push(txt);
	    	node.txt = txt;
		wedge.mouseover((function () {var t = txt; return function () { t.show(); }; })()).mouseout((function () {var t = txt; return function () { t.hide(); }; })());
		if (node.level == 7.0) {
			txt.dblclick((function () {var l = label; return function () { if (l.match(/^S[0-9]+$/)) { window.open("http://rdp.cme.msu.edu/hierarchy/detail.jsp?format=genbank&seqid=" + l); } }; })());
		}
		if (side != node.level) {
	    		txt.click((function () { var n = node; return function () { sideTree(node, label); }; })());
			wedge.click((function () { var n = node; return function () { sideTree(node, label); }; })());
		}
	} else {
		var wedge = R.path(d).attr({ stroke : border, "stroke-width" : +0.3*font_param+"px"}).toFront();
 		s.push(wedge);
	}

	return wedge;
};

var ring_gap = -1;

var tot = tree_struct.sub.Root.pop;

var main_set = r.set();
var side_set = null;

function sideTree(node, name) {

	if ((cur_s_tree) && (cur_s_tree.outline)) {
		cur_s_tree.outline.remove();
		cur_s_tree.outline = null;
		if (side_set) { side_set.remove(); }
	}

	cur_s_tree = node;
	cur_s_tree_name = name;
	cur_s_tree.outline = pieWedge(r, main_set, node.cum, node.pop, tot, R*1.85, cent, (node.level-1)*rad_inc, 7*rad_inc, "", "white", "", node, 0); 
	side_set = r.set();	
	walk_tree(node, R*5.60, cent, side_set, node.pop, name, node.level)
    	var attr2 = {font: 2.0*font_param+"px 'Helvetica'", stroke : "none", fill : "#FFFFFF"}; 
	var n = node.pop/tot * 100.0;
   	var txt = r.text(R*5.60, cent, name + "\n" + n.toFixed(2) + "%").attr(attr2);
	side_set.push(txt);
}

function walk_tree(node, cx, cy, s, t, name, side) {
	var b = "gray";
	var color = Raphael.hsb2rgb((node.cum+node.pop/2.0)/tot, 0.1+0.9*(node.level)/7.0, 0.1+0.8*node.w_conf).hex;
		if (node.level > 0) { 
			pieWedge(r, s, node.cum, node.pop, t, cx, cy, (node.level-1)*rad_inc, (node.level)*rad_inc-ring_gap, color, b, name, node, side);
		};
		for (var x in node.sub) { 
			walk_tree(node.sub[x], cx, cy, s, t, x, side);
		}
};

walk_tree(tree_struct.sub.Root, R*1.85, cent, main_set, tot, "Root", 0);

if (cur_s_tree) {
	sideTree(cur_s_tree, cur_s_tree_name);
}

button.onclick = function () {
	var str = serializer.serializeToString(r.canvas);
	code.value = str; 
};

};
