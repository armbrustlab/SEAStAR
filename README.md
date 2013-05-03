<link href="style.css" media="screen" rel="stylesheet" type="text/css" />

SEAStAR - A framework for the analysis of next-generation metagenomes (and more)
==============================

The Basics
------------------------------

####SEAStAR is a package of tools supporting the construction of complete analysis pipelines for next-generation (Illumina&reg;, SOLiD&trade;) sequencing data generated from environmental samples.  
##### It includes high-performance tools for dealing with:

* Converting between file formats (CSFASTA -> FASTQ)
* Trimming raw reads for quality (with tuning support)
* PCR de-duplication of paired reads (without reference sequences)
* Selecting and estimating the relative abundance of sequences from large reference databases (e.g. 16S rDNA)
* Sub-sampling paired FASTQ files randomly, or based on reads included in (or excluded from) reference alignments
* Converting assembled color-space (SOLiD) contigs to nucleotide-space
* Connecting assembled contigs together via paired reads (constructing an assembly graph)
* Splitting complicated metagenomic assembly graphs into well-supported scaffolds
* Binning scaffolds by organism using tetra-nucleotide statistics

#####SEAStAR works with, but does not supply:

* Short-read sequence aligners (e.g. [BWA](http://bio-bwa.sourceforge.net), [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml))
* De novo contig assemblers (e.g. [Velvet](http://www.ebi.ac.uk/~zerbino/velvet/))
* Tools for visualizing assembly graphs (e.g. [GraphViz](http://www.graphviz.org/), [ZGRViewer](http://zvtm.sourceforge.net/zgrviewer.html))
* 16S Taxonomic classifiers (e.g. [RDP Classifier](http://sourceforge.net/projects/rdp-classifier/)) 

You can find out more about SEAStAR on its [Armbrust Lab Homepage](http://armbrustlab.ocean.washington.edu/SEAStAR) and obtain news regarding updates and related info by following [@SEAStAR_meta](https://twitter.com/SEAStAR_meta) on Twitter. 

This file contains information on how to build and install the SEAStAR tools. For information on using the tools themselves, please see the included SEAStAR User Guide file.

License
------------------------------
SEAStAR is released under the GPLv3 license, a copy of which is provided in the included file "COPYING". By using this software, you are agreeing to be bound by the terms of this license.

Installation
------------------------------
The instructions that follow are for building the SEAStAR tools from source code. However, if you'd initially like to try out SEAStAR by working through our "Quick Start" tutorial, we provide a "ready-to-go" [SEAStAR Virtual Machine appliance](http://armbrustlab.ocean.washington.edu/node/305) image that includes all of the necessary tools and a working sample dataset that can be run within VirtualBox without the need to compile any of the tools on your computer. This pre-built VM is intended as an aid to learning and we strongly advise against (and will not provide any help for) trying to use it for analysis of real datasets.

SEAStAR is designed to build and run on any 64-bit Unix-like system, including Linux and Mac OS X version 10.7 or later. Many components of SEAStAR are optimized for multiple CPU cores and require substantial memory. We recommend a machine with a minimum of 4 CPU cores and 32 GB of RAM to run these tools.  Depending on your datasets and what you are trying to do (e.g. de novo assembly) you may require a substantially more powerful machine than this minimum recommendation. 

The SEAStAR package has dependencies on a small number of software packages that you (or your system administrator) may need to install or update. The process described in the next section will notify you if your system is missing any of these components.

####Required Tools:
 
* [gcc](http://gcc.gnu.org) -- version 4.4 or newer, supporting [OpenMP](http://openmp.org). 
* [cmake](http://www.cmake.org) -- version 2.8 or newer
* [node.js](http://nodejs.org) -- version 0.10 or newer
* [gawk](http://www.gnu.org/software/gawk/) -- version 3.1.5 or newer (version 4.0 or newer recommended))

Additional instructions are available below for fulfilling these requirements for Mac OS X, and for programmers wishing to make modifications to the included source code.

Once you have the above packages: To build SEAStAR using Unix style command line tools, run the following commands from the directory where all files generated in the build process should be placed (including executables). This is your "destination tree". 

    cmake [dir] 
    make

Where [dir] is the path to the root of the SEASTAR source tree (where this README file is found). 

If the path "." is used for [dir] above (run from the "source tree"), then the binary and source tree will be the same (an "in-source build"). After a successful make, executables will be found in the bin/ subdirectory. 

This directory (the bin subdirectory of the destination tree) should be added to your PATH environment variable, so that newly built tools can be found from your data analysis directories:

    export PATH=$PATH:[dest_dir]/bin   # Where [dest_dir] is the fully qualified path to your destination tree directory 

To test the newly built components:

    make test

If any tests fail, do not use the executables!

To clean all files generated in the source directory for an in-source build (this will only work for git checked-out repositories):

    git clean -fxd

For an out-of-source build you can simply delete the destination tree directory and start again.

Additional installation details for Mac OS X
------------------------------

For Mac OS X users: To fulfill the above requirements, you will first need to download and install Apple's XCode developer package (using the [App store](https://developer.apple.com/xcode/index.php)), and then we recommend installing the other required packages using [MacPorts](http://www.macports.org/).

Visit the link below to download and install MacPorts.
> http://www.macports.org/install.php

Then run the following commands to install the required packages:

    sudo port selfupdate
    sudo port install cmake
    sudo port install node
    sudo port install gawk
    sudo port install gcc44

Note, it is possible to use a newer version of gcc (such as gcc46) if you already have it installed or for some other reason. However, most of our OS X testing has been performed with gcc version 4.4.

It may be possible to use [HomeBrew](http://mxcl.github.com/homebrew/) instead of Macports to install these packages, but we haven't tested it.

You will need to define an environment variable to explicitly tell cmake which compiler to use. Note that this must be done each time you start a command line session where you wish to run cmake again (or add it to your shell startup file, e.g. .bashrc in your home directory).  For example:
 
    export CC=/opt/local/bin/gcc-mp-4.4  # change this if you are using a different version!

####An important note about compilers on Mac OS X (starting with Xcode 4.2 on OS X 10.7 and later versions):

There are known bugs in OpenMP (multi-core processor support) in the gcc compiler Apple supplies with Xcode on OS X Lion (or later).  This is why we specify above that you install gcc separately with MacPorts. The cmake script provided checks OS X systems to see if the OpenMP support is working correctly with the default (or specified) C compiler. If you receive an error when trying to build that says "You need to install gcc (version 4.4 or later) from MacPorts" it is because our build system is attempting to use the XCode gcc compiler, and not the one you installed from MacPorts.

For Developers
------------------------------

Some of the included JavaScript (.js) files are automatically generated from [CoffeeScript](http://coffeescript.org) source files (CoffeeScript is a [transcompiled](http://en.wikipedia.org/wiki/Source-to-source_compiler) dialect of JavaScript with Python-like syntax.) If you wish to modify these components, please edit the .coffee files in the scripts/ subdirectory of the source tree. The make system will automatically regenerate the .js files in the bin/ subdirectory of the destination tree. To successfully transcompile these files, you will need the to install the CoffeeScript package for node.js:

    sudo npm install -g coffee-script

It is sometimes useful to build with GCC debug flags turned on.  To achieve this follow the normal cmake build procedure with one additional user defined cmake cache entry:

    cmake -D DEBUG=ON [dir]

For Mac OS X Developers
------------------------------

The following section covers using cmake to build XCode project files. However, we do not recommend using executables built by XCode for anything other than development and debugging purposes due to the aforementioned bugs in the XCode compilers. This may change in the future, but as of now beware of these known issues with OpenMP in XCode.

To make XCode project files (for Mac OS X only):

    cmake -G Xcode [dir] 

Where [dir] is the path to the root of the destination binary tree. If the path '.' is used, then the binary and source tree will be the same (i.e. an in-source build). You may then load the SEAStAR project file into XCode to build, debug, etc.

Alternatively, an xcode project may be built on the command line as (choosing Debug or Release as appropriate):

    xcodebuild -alltargets -configuration [Debug|Release] 

A word of Warning: once the project is imported into XCode, the destination tree will not be backed-up by Time Machine on OS X. For in-source builds, the binary and source trees are the same directory, so Time Machine will not back up your source code changes if you develop and build within the source tree. For this reason, it is highly advisable to do out-of-source builds when developing in XCode, unless you back up your local git repository via a mechanism other than Time Machine.




