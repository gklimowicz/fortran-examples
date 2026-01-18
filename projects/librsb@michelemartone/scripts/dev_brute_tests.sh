#!/bin/bash
#
# Copyright (C) 2008-2019 Michele Martone
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

# This script is intended for the librsb developer usage.

set -e
set -x

/sbin/ldconfig -p | grep libefence.so
# in the following line, libefence may fail when librsb allocates e.g. 0 bytes (it happens)
LD_PRELOAD=libefence.so.0.0 rsbench --help
LD_PRELOAD=libefence.so.0.0 rsbench --version
LD_PRELOAD=libefence.so.0.0 rsbench --memory-benchmark
LD_PRELOAD=libefence.so.0.0 rsbench -Q1 # || exit 255 # efence does not stand zero sized reallocs
LD_PRELOAD=libefence.so.0.0 rsbench -B # note: efence does not stand zero sized reallocs
