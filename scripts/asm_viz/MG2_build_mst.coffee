# Determine window size
w = window.innerWidth
h = window.innerHeight
r = Math.min(w, h) / 2

gc_hue = d3.scale.linear().domain([30.0,60.0]).range([0.0,360.0]).clamp(true)
link_scale = d3.scale.log()

force = d3.layout.force()
   .charge(-40)
   .linkDistance(1)
   .linkStrength((d)-> link_scale(d.bits))
   .gravity(0.05)
   .size([w, h])

svg = d3.select("#chart").append("svg")
   .attr("width", w)
   .attr("height", h)
 
bb_scale = 1 
 
# Create a top level group for global transformations, place in center of the canvas
vis = svg.append("g")
    .attr("transform", "translate(#{w/2},#{h/2})scale(1)")
 
vis.append("g").attr("id","links") 
vis.append("g").attr("id","nodes")

# Function called to rescale/reposition the rendering within the current window dimensions
window.onresize = () -> 
    new_w = window.innerWidth
    new_h = window.innerHeight
    new_r = Math.min(new_w, new_h) / 2
    svg.attr("width", new_w).attr("height", new_h)
    vis.attr("transform", "translate(#{new_w/2},#{new_h/2})scale(#{bb_scale*new_r/r})")
 
my_bisect = d3.bisector((d) -> d.bits).right

d3.json("MG2.json", (json) ->

   nodes = (for name,n of json.nodes
      n.id = name
      n.links = []
      n)

   # link_pool = (e for e in json.edges when e.bits > 100.0)

   link_pool = (e for e in json.edges)

   link_pool.sort((a,b) -> d3.descending(a.bits,b.bits))
 
   for l in link_pool
      l.source = json.nodes[l.node1]
      l.target = json.nodes[l.node2]
      l.source.links.push([l,l.target])
      l.target.links.push([l,l.source])
   
   link_scale.domain(d3.extent(link_pool,(d) -> d.bits)).range([1,5])
 
   links = []
 
   force.nodes(nodes)
      .links(links)
      .start()

   link = svg.select("#links").selectAll("line.link")
      .data(links)
      .enter().append("line")
      .attr("class", "link")
      .style("stroke-width", (d) -> Math.log(d.bits,2)-4)
 
   first_link = d3.select(link)
 
   node = svg.select("#nodes").selectAll("circle.node")
      .data(nodes)
      .enter().append("circle")
      .attr("class", "node")
      .attr("r", (d) -> Math.log(d.sequence_length,2)-3)
      .style("fill", (d) -> d3.hsl(gc_hue(d.percent_gc),1,0.45))
      .call(force.drag)
 
   node.append("title")
      .text((d) -> d.id)
 
   force.on("tick", () ->
      bb = d3.max([d3.max(nodes, (d) -> Math.abs(d.x-0.5*w)), d3.max(nodes, (d) -> Math.abs(d.y-0.5*h))])  
      
      bb_scale = r/(1.05*bb) 
      window.onresize()
    
      link.attr("x1", (d) -> d.source.x - 0.5*w)
          .attr("y1", (d) -> d.source.y - 0.5*h)
          .attr("x2", (d) -> d.target.x - 0.5*w)
          .attr("y2", (d) -> d.target.y - 0.5*h)
 
      node.attr("cx", (d) -> d.x - 0.5*w)
          .attr("cy", (d) -> d.y - 0.5*h)
   )

   timer = window.setInterval(do () ->
      n = nodes[0]
      nodes_added = {}
      nodes_added[n.id] = n
      nodes_seen = {}
      edge_cand = []
      
      () ->
               
         # Add edges into their properly sorted positions unless the other node is already added
         for l in n.links 
            unless (nodes_added[l[0].source.id]? && nodes_added[l[0].target.id]?)
               edge_cand.splice(my_bisect(edge_cand,l[0].bits),0,l[0])  # MST  
               # edge_cand.push(l)   # Depth first
               # edge_cand.unshift(l)  # Breadth first 
               
         # Go through sorted candidate edges looking for nodes to insert or edges to reject
         while (e = edge_cand.pop())
            if (nodes_added[e.source.id]? && not nodes_added[e.target.id]?)
               n = e.target
               break
            else if (not nodes_added[e.source.id]? && nodes_added[e.target.id]?)
               n = e.source
               break
             
         new_link = e
         nodes_added[n.id] = n

         if new_link
            
            links.push(new_link)
            force.links(links).start()
     
            svg.select("#links").selectAll("line.link").data(links)
               .enter().append("line")
               .attr("class", "link")
               .style("stroke-width", (d) -> Math.log(d.bits,2)-4)
               .attr("x1", (d) -> d.source.x - 0.5*w)
               .attr("y1", (d) -> d.source.y - 0.5*h)
               .attr("x2", (d) -> d.target.x - 0.5*w)
               .attr("y2", (d) -> d.target.y - 0.5*h)
            link = svg.select("#links").selectAll("line.link")
         else 
            window.clearInterval(timer)
   ,50) 
#
#   window.setTimeout(() -> 
#      time_count = 0
#      timer = window.setInterval(() -> 
#         force.resume()
#         time_count++
#         force.gravity(time_count*.01)
#         if time_count > 10
#            window.clearInterval(timer)
#            link_scale.range([0,10])
#            force.start()
#      ,3000)
#   ,5000)
#   
)
