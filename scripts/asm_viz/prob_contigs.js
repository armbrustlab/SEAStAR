var i_contig = 0;
var nodes = new Array();
var nodes_name2index = new Array();
var width = 1300;
var height = 500;
var xpad = 50;
var ypad = 20;
var axisPad = 3;
var min_px_per_datum = 2;
var step = 1;
var bit_thresh = 10.0;
var button_w = 120;
var button_h = Math.floor(button_w/1.618);
var button_pad = 4;
var c;  // contig
var edges;  // edges for current contig
var i_switch = 0;
var switch_choices = ["per_nt_cov", "per_nt_phys_cov", "per_nt_mp_ins"];
var prob_colors = {
    "Physical coverage break": "red",
    "Collapsed duplication": "yellow",
    "Small insert": "blue"
};
var global_mean_cov;
var global_dup_thresh;
var prob_mark_data;
var data;


var info = d3.select("#info");

// Control area (buttons, sliders, etc)
var control = d3.select("#control").append("svg")
    .attr("width", width)
    .attr("height", button_h + 2*ypad);

// Main plot area
var svg_viz = d3.select("#viz").append("svg")
    .attr("width", width)
    .attr("height", height);

// Scales
var xScale = d3.scale.linear()
    .domain([0, 1000])
    .range([0, width-xpad*2]);
xScale.axis = d3.svg.axis().scale(xScale);

var yScale = d3.scale.linear()
    .domain([0, 1000])
    .range([height-ypad*3, 0]);
yScale.axis = d3.svg.axis().scale(yScale).orient("left");

// SVG groups
// y axis
var yAxis = svg_viz.append("g")
    .attr("class", "axis yAxis")
    .attr("transform", "translate(" + (xpad-axisPad) + "," + ypad + ")")
    .call(yScale.axis);

// x axis
var xAxis = svg_viz.append("g")
    .attr("class", "axis xAxis")
    .attr("transform", "translate(" + xpad + "," + (height-ypad+axisPad) + ")")
    .call(xScale.axis);

// line plot with filled area, [coverage, physical coverage, mp insert size]
var path_group = svg_viz.append("g")
    .attr("class", "contig_paths")
    .attr("transform", "translate(" + xpad + ", " + ypad + ")")
    .attr("fill", "steelblue");

// edge lines
var epath_group = svg_viz.append("g")
    .attr("class", "epaths")
    .attr("transform", "translate(" + xpad + ", " + ypad + ")")
    .attr("opacity", .5);

// problem markers
var prob_group = svg_viz.append("g")
    .attr("class", "problems")
    .attr("transform", "translate(" + xpad + ", " + ypad + ")")
    .attr("fill-opacity", .3)
    .attr("stroke", "black");

// global mean coverage line
var global_cov_group = svg_viz.append("g")
    .attr("class", "global_cov")
    .attr("transform", "translate(" + xpad + ", " + ypad + ")")
    .attr("opacity", .5)
    .attr("stroke-width", 2)
    .attr("stroke-dasharray", 5)
    .attr("stroke", "red");

// global duplicate threshold line
var global_dup_group = svg_viz.append("g")
    .attr("class", "global_dup")
    .attr("transform", "translate(" + xpad + ", " + ypad + ")")
    .attr("opacity", .5)
    .attr("stroke-width", 2)
    .attr("stroke-dasharray", 5)
    .attr("stroke", "green");

// Buttons
var button_group = control.append("g")
    .attr("class", "buttons")
    .attr("transform", "translate(" + xpad + ", " + ypad + ")");

button_group.append("rect")  // background for buttons, good when overlaid
    .attr("width", button_w*2 + button_pad)
    .attr("height", button_h)
    .attr("fill", "white");

button_group.append("text")
    .attr("class", "back")
    .attr("x", Math.floor(button_w/2))
    .attr("y", Math.floor(button_h/2))
    .text("back");
var switch_text = button_group.append("text")
    .attr("class", "switch")
    .attr("x", Math.floor(button_w + button_pad + button_w/2))
    .attr("y", Math.floor(button_h/2))
    .text(switch_choices[i_switch]);
button_group.append("text")
    .attr("class", "next")
    .attr("x", Math.floor(2 * (button_w + button_pad) + button_w/2))
    .attr("y", Math.floor(button_h/2))
    .text("next");
var back_button = button_group.append("rect")
    .attr("class", "back")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", button_w)
    .attr("height", button_h)
    .attr("fill", "white")
    .attr("opacity", .3)
    .attr("stroke", "black");
var switch_button = button_group.append("rect")
    .attr("class", "switch")
    .attr("x", button_w + button_pad)
    .attr("y", 0)
    .attr("width", button_w)
    .attr("height", button_h)
    .attr("fill", "white")
    .attr("opacity", .3)
    .attr("stroke", "black");
var next_button = button_group.append("rect")
    .attr("class", "next")
    .attr("x", 2 * (button_w + button_pad))
    .attr("y", 0)
    .attr("width", button_w)
    .attr("height", button_h)
    .attr("fill", "white")
    .attr("opacity", .3)
    .attr("stroke", "black");

// Area path generator function
var area = d3.svg.area()
    .x(function(d, index) { return xScale(index); })
    .y0(height-3*ypad)
    .y1(function(d) { return yScale(d); })
    .interpolate("monotone");

// line for edges
var eline = d3.svg.line()
    .x(function(d) { return xScale(Math.floor(Math.abs(d.x)/step)); })
    .y(function(d) { return yScale(d.y); })
    .interpolate("basis");

// lines for global mean coverage and global duplicate threshold
var gline = d3.svg.line()  // global coverage lines
    .x(function(d) { return xScale(d.x); })
    .y(function(d) { return yScale(d.y); })
    .interpolate("basis");


function viz() {
    update_data();
    update_scales();
    update_axes();
    update_header_text();
    update_main_graph();
    update_global_lines();
    update_edges();
    update_problem_markers();
    update_buttons();
}

function update_data() {
    update_current_contig_and_anchor();
    c = nodes[i_contig].node;
    edges = nodes[i_contig].edges;
    data = subsample(c[switch_choices[i_switch]]);
}

function update_current_contig_and_anchor() {
    // If anchor is set and real contig name use that to update current contig
    var anchor_contig_name = unescape(self.document.location.hash.substring(1));
    if (anchor_contig_name) {
        // Anchor has been set so go to the correct contig
        if (nodes_name2index.hasOwnProperty(anchor_contig_name)) {
            i_contig = nodes_name2index[anchor_contig_name];
        }
    }
    // Make sure anchor reflects current contig
    update_anchor();
}

function update_anchor() {
    self.document.location.hash = "#" + nodes[i_contig].node.chunk_name;
}

function update_scales() {
    // Recalculate scales
    xScale.domain([0, data.length]).range([0, width-xpad*2]);
    yScale.domain([0, d3.max([d3.max(data), global_dup_thresh])]).range([height-ypad*3, 0]);
    // Update the x axis tick format function to account for subsampling
    xScale.axis.tickFormat(function(t) { return t * step; });
}

function update_header_text() {
    var header = info.selectAll(".header").data([c]);
    header
        .enter()
        .append("pre")
            .attr("class", "header");
    header
        .exit()
        .remove();
    header
        .transition()
        .text(function(d) {
            var s = "Problem contig\n";
            s += "  name == " + d.chunk_name + "\n";
            s += "  length == " + d.seq_len + "\n";
            for (var cp in d.contig_problems) {
                var prob = d.contig_problems[cp];
                s += "    " + prob.type + "\n";
                s += "      " + prob.start + " - " + prob.end + "\n";
            }
            return s;
        });
}

function update_axes() {
    // Redraw axes
    xAxis
        .transition()
        .duration(500)
        .ease("exp-in-out")
        .call(xScale.axis);
    
    yAxis
        .transition()
        .duration(500)
        .ease("exp-in-out")
        .call(yScale.axis);
}

function update_main_graph() {
    var paths = path_group.selectAll("path").data([data]);
    paths
        .enter()
        .append("path");
    paths
        .exit()
        .remove();
    paths
        .transition()
        .duration(500)
        .ease("exp-in-out")
        .attr("d", area);
}

function update_global_lines() {
    var gcov, gdup;
    if (switch_choices[i_switch] == "per_nt_mp_ins") {
        gcov = global_cov_group.selectAll("path").data([]);
        gdup = global_dup_group.selectAll("path").data([]);
    } else {
        var global_cov_data = [{"x": 0, "y": global_mean_cov},
                           {"x": data.length-1, "y": global_mean_cov}];
        var global_dup_data = [{"x": 0, "y": global_dup_thresh},
                           {"x": data.length-1, "y": global_dup_thresh}];
        gcov = global_cov_group.selectAll("path").data([global_cov_data]);
        gdup = global_dup_group.selectAll("path").data([global_dup_data]);
    }
    
    gcov
        .enter()
        .append("path");
    gcov
        .exit()
        .remove();
    gcov
        .transition()
        .duration(500)
        .ease("exp-in-out")
        .attr("d", gline);
    
    gdup
        .enter()
        .append("path");
    gdup
        .exit()
        .remove();
    gdup
        .transition()
        .duration(500)
        .ease("exp-in-out")
        .attr("d", gline);
}

function update_problem_markers() {
    if (! c.hasOwnProperty("contig_problems")) {
        var probs = prob_group.selectAll(".problem").data(new Array());
    } else {
        var probs = prob_group.selectAll(".problem").data(c.contig_problems);
    }
    probs
        .enter()
        .append("rect");
    probs
        .exit()
        .remove();
    probs
        .transition()
        .duration(500)
        .ease("exp-in-out")
        .attr("class", "problem")
        .attr("x", function(d) { return xScale(Math.floor(d.start/step)); })
        .attr("y", height-ypad*2.8)
        .attr("width", function(d) {
            return xScale(Math.floor((d.end/step) - (d.start/step)));
        })
        .attr("height", 15)
        .attr("fill", function(d) {
            return prob_colors.hasOwnProperty(d.type) ? prob_colors[d.type] : "black";
        });
}

function update_edges() {
    // Draw edges to other contigs and write edge summary text
    var e, other_name, other_pos, edge_pos;
    var edge_path_data = new Array();
    var edge_text_data = new Array();
    
    // Filter edges by bit score and build path and text data objects
    for (var e_index in edges) {
        e = edges[e_index];
        if (e.bits > bit_thresh) {
            if (e.n1 != c.name) {
                other_name = e.n1;
                other_pos = e.p1;
                edge_pos = e.p2;
            } else {
                other_name = e.n2;
                other_pos = e.p2;
                edge_pos = e.p1;
            }
            edge_path_data.push([{"x": edge_pos, "y": 0},
                                 {"x": edge_pos, "y": d3.max([d3.max(data), global_dup_thresh])}]);
            edge_text_data.push(e);
        }
    }
    
    // Draw edges
    var epaths = epath_group.selectAll(".epath").data(edge_path_data);
    epaths
        .enter()
        .append("path")
            .attr("class", "epath")
            .attr("stroke", "blue")
            .attr("stroke-width", 2);
    epaths
        .exit()
        .remove();
    epaths
        .transition()
        .duration(500)
        .ease("exp-in-out")
        .attr("d", eline);
    
    // Write edge JSON objects
    var etexts = info.selectAll(".etext").data(edge_text_data);
    etexts
        .enter()
        .append("pre")
            .attr("class", "etext");
    etexts
        .exit()
        .remove();
    etexts
        .transition()
        .text(function(d) { return JSON.stringify(d, null, "  "); });
}

function update_buttons() {
    // Click back button to return to previous contig
    back_button.on("click", function() {
        i_contig = (i_contig + nodes.length - 1) % nodes.length;
        update_anchor();
        viz();
    });
    
    // Click next button to advance to next contig
    next_button.on("click", function() {
        i_contig = (i_contig + 1) % nodes.length;
        update_anchor();
        viz();
    });
    
    // Switch between different data sets for this contig
    switch_button.on("click", function() {
        i_switch = (i_switch + 1) % switch_choices.length;
        switch_text
            .transition()
            .text(switch_choices[i_switch]);
        viz();
    });
}

function subsample(d) {
    if (d.length <= (width / min_px_per_datum)) {
        // no need to subsample
        step = 1;
        return d;
    }
    step = Math.ceil(d.length / (width / min_px_per_datum));
    var i, j, sum;
    var new_data = Array();
    for (i=0; i < d.length - step; i+=step) {
        sum = 0;
        for (j=0; j < step; j++) {
            sum += d[i+j];
        }
        new_data.push(sum/step);
    }
    sum = 0;
    for (j=i; j < d.length; j++) {
        sum += d[j];
    }
    new_data.push(sum/(d.length-i));
    return new_data;
}

// json_file should be set in HTML document
d3.json(json_file,
    function(j) {
        // Grab global data stats
        global_mean_cov = j.run_stats.mean_coverage;
        global_dup_thresh = j.run_stats.coverage_duplication_threshold;
        
        // Turn node objects into array
        var e;
        var i = 0;
        for (var n in j.nodes) {
            if (j.nodes[n].hasOwnProperty("contig_problems")) {
                for (var cp in j.nodes[n].contig_problems) {
                    //if (j.nodes[n].contig_problems[cp].type == "Physical coverage break") {
                    if (1) {
                        j.nodes[n].chunk_name = n;
                        nodes[i] = {"node": j.nodes[n], 
                                    "edges": new Array()};
                        nodes_name2index[n] = i;
                        for (var e_index in j.edges) {
                            e = j.edges[e_index];
                            if (e.n1 == nodes[i].node.name || e.n2 == nodes[i].node.name) {
                                nodes[i].edges.push(e);
                            }
                        }
                        i++;
                        break;
                    }
                }
            }
        }
        viz();
    }
);
