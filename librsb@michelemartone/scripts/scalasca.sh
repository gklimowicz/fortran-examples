#!/bin/sh
#
# Copyright (C) 2008-2015 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

# These suggestions are intended for the librsb developer usage.
# E.g.: with scalasca-2.1, scorep-1.3, cube-4.2.3 .
# (dependencies are: papi-5.3, libqt, ... )
#
# sh autogen.sh
# ./configure CC=gcc        FC=gfortran CFLAGS='-fopenmp' --disable-dependency-tracking
# make        clean         # not to interfere with code generation commands
# make        CC='scorep gcc' FC='scorep gfortran' LD='scorep ld' -j 16
#
# Profiling
# scan ./rsbench ... 
# square ... # <give as argument the newly created analysis directory containing *cube* files>
#
# To collect PAPI events:
# export SCOREP_METRIC_PAPI=PAPI_L1_TCM,PAPI_L2_TCM SCOREP_ENABLE_TRACING=false
# ./rsbench ... 
# square ... # <give as argument the newly created analysis directory containing *cube* files>
# 
