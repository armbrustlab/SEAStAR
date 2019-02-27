# -------------------------------------------------------------------------- #
# Center for Environmental Genomics
# Copyright (C) 2009-2014 University of Washington.
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
# along with SEAStAR.  If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------- #

# CMake file for seastar_shared
cmake_minimum_required (VERSION 2.8.5)

project (seastar_shared C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake_modules/)

add_library(seastar_shared seastar_shared.c)

target_link_libraries(seastar_shared)

