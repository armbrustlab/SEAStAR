#!/bin/bash
# --------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2013 University of Washington.
# 
# Authors:
# Chris Berthiaume
# chrisbee@uw.edu
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
# Wrapper Bash shell script for RDP_tree_dev.js
# 
# Usage: RDP_tree_dev.js RDP_trainset_tree.json <sample_class.txt >sample_class.json

# This file is from SEAStAR version: SS_BUILD_VERSION

# JavaScript filename to run with node
js=RDP_tree_dev.js

# Run node script with custom memory settings
nodewrap "$js" "$@"
