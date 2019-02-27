#!/bin/bash
# --------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2015 University of Washington.
#
# Authors:
# Chris Berthiaume
# chrisbee@uw.edu
#
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
#
# Wrapper shell script to make any node JavaScript file executable on the PATH.
# This script sets memory limits for the V8 JavaScript engine based on physical
# memory size, leaving some amount of memory for other processes.
#
# usage: nodewrap JSFILE [ARGS]

# This file is from SEAStAR version: SS_BUILD_VERSION

# Get the real directory path for this script.  If this script is a symbolic
# link, follow the symbolic link to the actual file first.
function real_script_dir {
    link=$(readlink "$0")
    if [[ -z "$link" ]]; then
        # Not a symbolic link
        realpath="$0"
    else
        realpath="$link"
    fi
    thisdir="$(cd "$(dirname "$realpath")" && pwd)"
    echo "$thisdir"
}

if [[ $# -eq 0 ]]; then
    echo "ERROR: No arguments provided" >&2
    echo "usage: nodewrap JSFILE [ARGS]" >&2
    echo
    echo "JSFILE must be in the same directory as this script" >&2
    exit 1
fi

# Find directory of this script, use it to create full path to JavaScript file
thisdir=$(real_script_dir)
jsfile="$thisdir"/"$1"
shift  # shift JavaScript file from argument list

# Get system physical memory size in megabytes.
sysname=$(uname)
if [[ "$sysname" == "Linux" ]]; then
    pmem=$(free | awk 'NR == 2 {print $2}')
    pmem=$((pmem / 2**10))  # convert from KB to MB
elif [[ "$sysname" == "Darwin" ]]; then
    pmem=$(sysctl hw | awk '/hw.memsize:/ {print $2}')
    pmem=$((pmem / 2**20))  # convert from B to MB
else
    echo "ERROR: System type $sysname not supported" >&2
    exit 1
fi

if node --version | grep -q 'v0.10'
then
    newspace='--max-new-space-size=16384'
else
    newspace='--max_semi_space_size=16'
fi

headroom=$((2**10))  # 1 GB of physical memory head room for other processes
v8default=$((2**10))  # The size of the 64bit V8 engine default memory limit
maxold=$((pmem - headroom))

if [[ $maxold -gt $v8default ]]; then
    node --always_compact "$newspace" --max-old-space-size="$maxold" "$jsfile" "$@"
else
    node --always_compact "$newspace" "$jsfile" "$@"
fi
