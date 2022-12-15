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
which nproc
find -iname '*.gcda' -or -iname '*.info' -exec rm '{}' ';'

if true; then
#if ! test -f Makefile || ! grep CFLAGS.*coverage Makefile 2>&1 > /dev/null  ; then
	sh autogen.sh
	export CFLAGS=${CFLAGS:='--coverage -O0 -ggdb -pipe'}
	export CXXFLAGS=${CXXFLAGS:='--coverage -O0 -ggdb -pipe'}
	export FCFLAGS=${FCFLAGS:='--coverage -O0 -ggdb -pipe'}
	export LIBS='-lgcov'
	#./configure --enable-c-examples --enable-fortran-examples --prefix=`pwd`/local/librsb-coverage --enable-matrix-types=blas --with-zlib --disable-shared
	./configure --enable-c-examples --enable-fortran-examples --prefix=`pwd`/local/librsb-coverage --enable-matrix-types=all --with-zlib --disable-shared --enable-octave-testing --enable-debug --enable-debug-getenvs "$@"
	make clean
fi

#if ! grep CFLAGS.*coverage Makefile 2>&1 > /dev/null  ; then
#	echo "[!] Re-Make'ing with coverage flags !"
#	make clean
#	make -j `nproc` \
#		"`grep  ^CFLAGS Makefile | sed  's/^CFLAGS *=/CFLAGS= --coverage /g'`" \
#		"`grep ^FCFLAGS Makefile | sed 's/^FCFLAGS *=/FCFLAGS=--coverage /g'`" \
#		"`grep ^FCFLAGS Makefile | sed 's/^CXXFLAGS *=/CXXFLAGS=--coverage /g'`" \
#		"`grep ^'LIBS\>' Makefile | sed 's/^LIBS *=/LIBS=-lgcov /g'`" \
#	rsbench
#else
#	true;
#fi
cd examples 
cd -
if test -z "${RSB_WANT_CONFIGURE_ONLY}" ; then
if test -z "${RSB_WANT_LOW_MEMORY_BUILD}" ; then
	make -j `nproc` rsbench sbtc sbtf
	make -j `nproc`
else
	make -j `nproc` all
	make -j `nproc` sbtc
fi
fi

#rm -f *.gcda        *.gcov
