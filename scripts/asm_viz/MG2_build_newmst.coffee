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
   .gravity(0.1)
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
 
my_bisect = d3.bisector((d) -> d[0].bits).right

d3.json("MG2.json", (json) ->

   nodes = (for name,n of json.nodes
      n.id = name
      n.links = []
      n.internal_mp_bits = 0.0
      n)

   for e in json.internal_edges
      json.nodes[e.node1].internal_mp_bits = e.bits

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
      dead_ends = []
      edge_cand = []
      
      () ->
         next_n = null
         new_link = null
         
         new_territory = (node) ->
            next_l = null
            for l in node.links                    # See if this node has an "unseen" child
               if not nodes_added[l[1].id]? and (not nodes_seen[l[1].id]? or (l[1].internal_mp_bits / l[1].links[0][0].bits > 0.02))
                  next_l = l
                  break
            next_l

         if n                           # If this node has an unseen node
            for l in n.links                    # Walk through all of the (sorted) links for this node
               if nodes_added[l[1].id]? 
                  continue
                  
               if not nodes_seen[l[1].id]? or (l[0].bits > nodes_seen[l[1].id][0].bits)
                  nodes_seen[l[1].id] = l
                  edge_cand.splice(my_bisect(edge_cand,l[0].bits),0,l)
               
         next_n = edge_cand.pop()         # select a seen node to try
         while next_n and nodes_added[next_n[1].id]?
            next_n = edge_cand.pop()
               
         if next_n
            [new_link, n] = next_n
            nodes_added[n.id] = n  # Otherwise add this node to the selected pool
            unless new_territory(next_n[1])
               n = null
               
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
   ,10) 
)
