width = 3000
height = 3000
gc_hue = d3.scale.linear().domain([30.0,60.0]).range([0.0,360.0]).clamp(true)
link_scale = d3.scale.log()

force = d3.layout.force()
   .charge(-40)
   .linkDistance(1)
   .linkStrength((d)-> link_scale(d.bits))
   .gravity(0.01)
   .size([width, height])
 
svg = d3.select("#chart").append("svg")
   .attr("width", width)
   .attr("height", height)
 
svg.append("g").attr("id","links") 
svg.append("g").attr("id","nodes")
 
d3.json("MG2.json", (json) ->

   for l in json.edges
      l.source = json.nodes[l.node1]
      l.target = json.nodes[l.node2]
      
   nodes = (for name,n of json.nodes
      n.id = name
      n)
      
   link_pool = (e for e in json.edges when e.bits > 250.0)
   
   link_pool.sort((a,b) -> d3.descending(a.bits,b.bits))
   
   link_scale.domain(d3.extent(link_pool,(d) -> d.bits)).range([0,16])
 
   links = link_pool
 
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
      link.attr("x1", (d) -> d.source.x)
          .attr("y1", (d) -> d.source.y)
          .attr("x2", (d) -> d.target.x)
          .attr("y2", (d) -> d.target.y)
 
      node.attr("cx", (d) -> d.x)
          .attr("cy", (d) -> d.y)
   )


   window.setTimeout(() -> 
      time_count = 0
      timer = window.setInterval(() -> 
         force.resume()
         time_count++
         force.gravity(time_count*.01)
         if time_count > 10
            window.clearInterval(timer)
            link_scale.range([0,10])
            force.start()
      ,3000)
   ,5000)
   
)
