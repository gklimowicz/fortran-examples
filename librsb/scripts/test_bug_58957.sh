#!/bin/bash
#
# Copyright (C) 2008-2020 Michele Martone
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

# https://savannah.gnu.org/bugs/index.php?58957
set -e
set -x
test -x rsbench
cat > crash.mtx << EOF
%%MatrixMarket matrix coordinate real symmetric
10 10 0
EOF
./rsbench --plot-matrix -f crash.mtx # if buggy, crashes
