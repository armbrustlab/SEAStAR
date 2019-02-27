<link href="style.css" media="screen" rel="stylesheet" type="text/css" />

[![Build Status](https://travis-ci.org/armbrustlab/SEAStAR.svg?branch=master)](https://travis-ci.org/armbrustlab/SEAStAR)

SEAStAR - A framework for the analysis of next-generation metagenomes
------------------------------

### The Basics

#### SEAStAR is a package of tools supporting the construction of complete analysis pipelines for next-generation (Illumina&reg;, SOLiD&trade;) sequencing data generated from environmental samples.
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
* Identifying small circular scaffolds that are likely virus or plasmid genomes

##### SEAStAR works with, but does not supply:

* Short-read sequence aligners (e.g. [BWA](http://bio-bwa.sourceforge.net), [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml))
* De novo contig assemblers (e.g. [Velvet](http://www.ebi.ac.uk/~zerbino/velvet/))
* Tools for visualizing assembly graphs (e.g. [GraphViz](http://www.graphviz.org/), [ZGRViewer](http://zvtm.sourceforge.net/zgrviewer.html))
* 16S Taxonomic classifiers (e.g. [RDP Classifier](http://sourceforge.net/projects/rdp-classifier/))

You can find out more about SEAStAR on its [Armbrust Lab Homepage](https://armbrustlab.ocean.washington.edu/tools/seastar/).

This file contains information on how to build and install the SEAStAR tools. For information on using the tools themselves, please see the included SEAStAR User Guide file.

### License
SEAStAR is released under the GPLv3 license, a copy of which is provided in the included file "COPYING". By using this software, you are agreeing to be bound by the terms of this license.

### Installation
The instructions that follow are for building the SEAStAR tools from source code.

SEAStAR is designed to build and run on any 64-bit Unix-like system, including Linux and Mac OS X version 10.7 or later. Many components of SEAStAR are optimized for multiple CPU cores and require substantial memory. We recommend a machine with a minimum of 4 CPU cores and 32 GB of RAM to run these tools.  Depending on your datasets and what you are trying to do (e.g. de novo assembly) you may require a substantially more powerful machine than this minimum recommendation.

The SEAStAR package has dependencies on a small number of software packages that you (or your system administrator) may need to install or update. The process described in the next section will notify you if your system is missing any of these components.

#### Required Tools:

* [gcc](http://gcc.gnu.org) -- version 4.2 or newer, supporting [OpenMP](http://openmp.org) (version 4.7 recommended)
* [cmake](http://www.cmake.org) -- version 2.8.5 or newer
* [node.js](http://nodejs.org) -- version 0.10 or newer
* [gawk](http://www.gnu.org/software/gawk/) -- version 3.1.5 or newer (version 4.0.2 recommended)

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

### Additional installation details for Mac OS X

For Mac OS X users: To fulfill the above requirements, you will first need to download and install Apple's "Command Line Developer Tools".

    xcode-select --install

And then we recommend installing the other required packages using [HomeBrew](https://brew.sh/) or [MacPorts](http://www.macports.org/).

#### HomeBrew (preferred)

Visit the link below to download and install HomeBrew.
> https://brew.sh/

Then run the following commands to install the required packages:

    brew update
    brew install cmake
    brew install node   # Node may also optionally be installed using nvm
    brew install gawk
    brew install gcc@8

#### MacPorts

Visit the link below to download and install MacPorts.
> http://www.macports.org/install.php

Then run the following commands to install the required packages:

    sudo port selfupdate
    sudo port install cmake
    sudo port install node  # Node may also optionally be installed using nvm
    sudo port install gawk
    sudo port install gcc8  # or whatever version you may prefer

### An important note about compilers on Mac OS X :

Xcode's default Clang-based compiler does not support OpenMP (a standard for writing efficiently parallelized C code); this is why we specify above that you must install the gcc compiler. The cmake script provided checks OS X systems to see if the OpenMP support is working correctly with the default (or specified) C compiler. If you receive an error when trying to build that says "You need to install gcc (version 4.4 or later)" it is because our build system is attempting to use the Xcode compiler, and not the one you installed using HomeBrew or MacPorts.

You will need to define an environment variable to explicitly tell cmake which compiler to use. Note that this must be done each time you start a command line session where you wish to run cmake again (or add it to your shell startup file, e.g. .bash_profile in your home directory).  For example:

    export CC=/usr/local/bin/gcc-8  # Homebrew

or

    export CC=/opt/local/bin/gcc-8  # MacPorts

Change the numbers above if you are using a different version!

Note: It may also be possible, with some more work, to use a newer version of clang than provided by Apple to compile with OpenMP support, but we mave not tested this. Both HomeBrew and MacPorts enable installation of LLVM 7 (which includes the clang C compiler). Have fun with that!

### For Developers

Some of the included JavaScript (.js) files are automatically generated from [CoffeeScript](http://coffeescript.org) source files (CoffeeScript is a [transcompiled](http://en.wikipedia.org/wiki/Source-to-source_compiler) dialect of JavaScript with Python-like syntax.) If you wish to modify these components, please edit the .coffee files in the scripts/ subdirectory of the source tree. The make system will automatically regenerate the .js files in the bin/ subdirectory of the destination tree. To successfully transcompile these files, you will need the to install the CoffeeScript package for node.js:

    sudo npm install -g coffee-script

It is sometimes useful to build with GCC debug flags turned on.  To achieve this follow the normal cmake build procedure with one additional user defined cmake cache entry:

    cmake -D DEBUG=ON [dir]
