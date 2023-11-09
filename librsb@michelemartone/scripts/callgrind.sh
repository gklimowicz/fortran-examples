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

# This script is intended for the librsb developer usage.
CMD=" ./rsbench  -oa -Ob --dense 100 --compare-competitors --verbose -R -qH --no-want-ancillary-execs -n 2"
valgrind  --tool=callgrind --cache-sim=yes --dump-instr=yes  --separate-threads=yes $CMD
# now one can use tools as e.g.: qcachegrind  to analyze 
