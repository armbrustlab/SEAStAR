# fastq_sub.coffee -- Coffeescript fastq subdivider
# 
# Arguments: <input filename> <number of parts> <part number to produce>
#
#

fs = require('fs')
util = require('util')

# Catch exceptions
process.on 'uncaughtException', (err) -> 
  process.stderr.write "Caught exception: #{err}\n"
  process.exit 1

# grab command line args:  
[fn,num_parts,part_num] = process.argv[2..4]  

# parse integers from cmd line strings
num_parts = parseInt(num_parts)
part_num = parseInt(part_num)

# validate inputs
if num_parts < 1 or not (0 < part_num <= num_parts)
  process.stderr.write "Error: number of parts must be > 0 and part number(s) must be [1...num_parts]\n"
  process.exit 1

len = 0
start = 0

# read file statistics
stats = fs.statSync(fn)

# validate stream type
if not stats.isFile()
   process.stderr.write "Error: Input must be a seekable file.\n"
   process.exit 1
else  # calculate the starting and ending positions
   len = stats.size
   start = Math.floor(len*(part_num-1)/num_parts)
   end = Math.floor(len*(part_num)/num_parts)
   process.stderr.write "Part #{part_num} of #{num_parts} of File: #{fn}\n"
   process.stderr.write "File is #{len} bytes long, starting at byte #{start}\n"

# Opening the stream and parsing it for lines, etc....
instream = fs.createReadStream(fn,{flags:'r',encoding:'ascii',bufferSize:Math.pow(2,16),start:start})
             .on('open',(fd) -> 
                process.stderr.write "Stream opened on file descriptor #{fd}\n")
             .on('data', do () ->
                save = ''         # Store unused characters from last line of previous buffer 
                save_lines = []   # Store unused lines (parts of a read) from previous buffer
                data_used = 0     # How much data has been consumed so far
                return (c) ->     # Function to call when a new data buffer is available
                   c = save + c   # prepend previous string
                   lines = c.split '\n'   # make an array of lines
                   lines = save_lines.concat(lines)  # prepend previous lines
                   save = lines.pop()                # save any partial line
                   # If there isn't a whole read at the beginning, drop it.
                   while lines.length > 4 and (lines[0][0] != '@' or lines[2][0] != '+')
                       data_used += lines[0].length + 1
                       lines.shift()
                   # Handle one read
                   while ((lines.length >= 4) and (start+data_used < end))
                       data_used += l.length+1 for l in lines[0..3]
                       # data_used += lines[0..3].reduce(((x,y)->(x+y.length)),4)
                       per_read(lines.splice(0,4))
                   # decide when to end    
                   if start+data_used >= end
                      process.stderr.write ("#{data_used} bytes consumed. Should stop at #{end}\n")   
                      process.exit 1 
                   
                   save_lines = lines   # Save the remaining lines for next time
                )

# Called per read for output etc.
per_read = (ra) ->
   process.stdout.write "#{i}\n" for i in ra
