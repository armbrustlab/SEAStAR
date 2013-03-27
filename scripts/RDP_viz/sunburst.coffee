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
#
# Renders an RDP JSON tree as an in-browser interactive "sunburst diagram" 
# (essentially a zoomable heirarchical pie chart).
#
# This program uses the D3.js data visualization library: http://mbostock.github.com/d3/
#
# Input: This script reads a JSON formatted tree file, such as those created by 
# the RDP_train_to_tree or RDP_tree scripts.  
#
# Output: On-screen and Downloadable SVG file for the visualized tree structure
#
# Embedding: A minimal .html file is require to make this work in an (HTML5 compliant) browser:
#
#<!DOCTYPE html>
#<html>
#  <head>
#    <script type="text/javascript" src="http://mbostock.github.com/d3/d3.min.js"></script>
#  </head>
#  <body>
#     <a href="" id="download_link" >"Right click, save/download" to download SVG file</a>
#     <script src="./RDP_merged.js" type="text/javascript" charset="utf-8"></script>
#     <script type="text/javascript" src="sunburst.js"></script>
#  </body>
#</html>
#
###

#        font_param = 8.0*(dim/1000.0),

# Change this to determine which dataset to load (via XMLHTTPRequest, inside d3.js)
json_filename = "./GG2_RDP.json"

# Automatically scale to the size of the browser window
w = window.innerWidth
h = window.innerHeight
r = Math.min(w, h) / 2
r -= 0.05*r   # leave a little margin

max_level = 7.0  #  Levels in the RDP taxonomy (sub/super are fractional)

# Set up d3 scaling functions for polor coordinates, mapping 0..1 to:
# x : angle in radians
# y : radius from the center in pixels

xscale = d3.scale.linear().range([0, 2 * Math.PI]).clamp(true)
yscale = d3.scale.pow().exponent(1.2).domain([0, 1]).range([0, r])

duration = 1000   # how long animated transitions should take (ms)

# Create a new DIV element for rendering into
chart = d3.select("body").append("div")
          .attr("class", "chart")

# Create the parent SVG object
svg = d3.select(".chart").append("svg")
    .attr("width", w)
    .attr("height", h)
    .style("background-color", "Black")

# Create a top level group for global transformations, place in center of the canvas
vis = svg.append("svg:g")
    .attr("transform", "translate(#{w/2},#{h/2})scale(1)")

# Hook an "export svg" function to an anchor with a specific id
link = d3.select("#download_link")
    .on("mouseover", () -> 
       html = d3.select("svg")    # get the SVG document
          .attr("version", 1.1)
          .attr("xmlns", "http://www.w3.org/2000/svg")
          .node().parentNode.innerHTML;
       # serialize SVG into the href of the provided anchor with base64 encoding
       # btoa with utf-8 strings can cause Character Out Of Range exception
       # See https://developer.mozilla.org/en-US/docs/DOM/window.btoa#Unicode_Strings for source of this solution
       d3.select(this)
          .attr("href-lang", "image/svg+xml")
          .attr("href", "data:image/svg+xml;base64,\n" + btoa(unescape(encodeURIComponent(html)))))

# Function called to rescale/reposition the rendering within the current window dimensions
window.onresize = () -> 
    new_w = window.innerWidth
    new_h = window.innerHeight
    new_r = Math.min(new_w, new_h) / 2 - 50
    svg.attr("width", new_w).attr("height", new_h)
    vis.attr("transform", "translate(#{new_w / 2},#{new_h / 2})scale(#{new_r/r})")

# Function used to parse and interpret the RDP tree datastructure 
partition = d3.layout.partition()
    .children((d) ->               # Child accessor function
       (for name, child of d.sub   # Walk the sub object and produce an array of children
          child.name ?= name       # Turn the propertry name into a name property of the object 
          child))
    .sort(null)
    .value((d) -> d.pop)    # The size of each node scale with population. For count, set to 1 

# Function remapping the parition coordinates using RDP levels (not all branches are equal length)
arc = d3.svg.arc()   # Note that the functions below "clamp" for the transition animations
    .startAngle((d) -> Math.max(0, Math.min(2 * Math.PI, xscale(d.x))))
    .endAngle((d) -> Math.max(0, Math.min(2 * Math.PI, xscale(d.x + d.dx))))
    .innerRadius((d) -> 
       Math.max(0, 
          yscale(d.y = (if d.parent? then ((0.5+d.parent.level)/(max_level+0.5)) else 0.0))))
    .outerRadius((d) -> 
       d.dy = (0.5+d.level)/(max_level+0.5) - d.y
       Math.max(0, yscale(d.y + d.dy)))

# Function to determine the proper radius for the center of node's text label 
text_x = (d) -> yscale((if d.y then (d.y + d.dy/2) else 0.0))

# Function to determine the SVG transform string to position a text label
text_xform = (d) -> 
    angle = ((xscale(d.x + d.dx / 2) - Math.PI / 2) / Math.PI * 180)
    twist = (if d.dx / (xscale.domain()[1]-xscale.domain()[0]) < 0.04 then 0.0 else 90.00)
    twist -= 180.0 if 90.0 < angle + twist < 270.0
    "rotate(#{angle})translate(#{text_x(d)},0)rotate(#{twist})"

# Function to calculate the percieved brightness of an RGB color
brightness = (rgb) -> (rgb.r * .299 + rgb.g * .587 + rgb.b * .114)

# Function to maximum radius reached by the children of a node
maxY = (d) -> (if d.children then Math.max.apply(Math, d.children.map(maxY)) else d.y + d.dy)

# Returns a d3 transition "tween" function for arc paths
arcTween = (d) ->
   my = maxY(d)  # this works by dynamically changing the scaling functions 
   xd = d3.interpolate(xscale.domain(), [d.x, d.x + d.dx])
   yd = d3.interpolate(yscale.domain(), [d.y, my])
   yr = d3.interpolate(yscale.range(), [(if d.y then 20 else 0), r])
   (d) -> ((t) -> ( 
      xscale.domain(xd(t))
      yscale.domain(yd(t)).range(yr(t))
      arc(d)))

# Function to determine whether one node is a descendant of another
isDescendantOf = (p, c) -> (
   if (p is c) then true else (if c.parent? then isDescendantOf(p, c.parent) else false))

# Load the JSON input data and set the visualization in motion!
#d3.json(json_filename, (json) -> 
do () ->
  json = tree_struct
  json.sub.Root.name = "Root";	# Name the root node
   
  nodes = partition.nodes(json.sub.Root)    # Generate a d3 parition on the input tree

  # Function called when arc wedge (d) is clicked 
  zoom = (d) -> 
     path.attr("display","")   # Enable all paths to be visible 
     
     path.transition().duration(duration)   #  Execute the path transition
         .attrTween("d", arcTween(d))       #  Animate the arcs
         .each("end", (e) ->                                  #  When the animation is over
            if not (isDescendantOf(d,e) or (e is d.parent))   #  Hide the paths that aren't
               d3.select(this).attr("display","none"))        #  decendents of the current root 
 
     text.transition().duration(duration)   #  Execute the text transition
         .attrTween("transform", (d) ->     #  Let d3 transition the text SVG transforms
            () -> text_xform(d))
         .each("end", (e) ->                       # When done animating, hide text for nodes that
            if (isDescendantOf(d,e))               # aren't the current root or a descendant
               d3.select(this).attr("display","")
            else 
               d3.select(this).attr("display","none"))

  path = vis.selectAll("path")     #  Create paths (the arc / wedges) from the data in the 
      .data(nodes)                 #  partitioned nodes
      .enter().append("svg:path")
      .attr("d", arc)
      .attr("fill-rule", "evenodd")
      .style("stroke", "#fff")      
      .style("fill", (d) ->     #  Color from the HSV "color wheel" space
         d.color = d3.hsl(360.0 * (d.x + d.dx / 2), 0.1+0.9*d.y, (0.1+0.8*d.w_conf)/2.0))
      #  Corresponding text labels are hidden until moused over
      .on("mouseover", (d,i) -> d3.select(text[0][i]).attr("visibility", "visible"))
      .on("mouseout", (d,i) -> d3.select(text[0][i]).attr("visibility", "hidden"))
      .on("click", zoom)
      
  text = vis.selectAll("text")    # Create and preposition hidden text labels for all nodes
      .data(nodes)
      .enter().append("svg:text")
      .style("opacity", 1)        # Change color from black to white if the wedge is too dim
      .style("fill", (d) -> (if brightness(d3.rgb(d.color)) < 90 then "#eee" else "#000"))
      .attr("dy", "0.3em")
      .attr("text-anchor","middle")
      .attr("transform", text_xform) 
      .attr("visibility", "hidden")
      .text((d) -> d.name)
      #  Corresponding text labels are hidden until moused over
      .on("mouseover", (d,i) -> d3.select(text[0][i]).attr("visibility", "visible"))
      .on("mouseout", (d,i) -> d3.select(text[0][i]).attr("visibility", "hidden"))
      .on("click", zoom)
          
  text.filter((d) -> d.level is max_level)   # Add an even to leaves only to open a new page
      .on("dblclick", (d) ->                 # for the corresponding RDP sequence when the name
         if d.name.match(/^S[0-9]+$/)        # is of the form e.g. S833457938 
            window.open("http://rdp.cme.msu.edu/hierarchy/detail.jsp?format=genbank&seqid=#{d.name}"))
