#!/bin/bash
#
# Copyright (C) 2008-2021 Michele Martone
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
which octave
which gcc
which m4
which nproc
./autogen.sh
./configure CFLAGS=-O0 FC=' '
make distclean
export CC='gcc' FC=' ' CFLAGS='-O0 -pipe '
./configure --enable-allocator-wrapper --enable-doc-build --disable-debug --prefix=`pwd`/local/librsb-nofortran --enable-matrix-types=blas --with-zlib
make clean
make cleanall
make -j `nproc`
make dist
make qqtests
make install
make itests
