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
if ! test -f Makefile ; then
	sh autogen.sh
	export CFLAGS='-O2 -pg -pipe '
	export FCFLAGS='-O2 -pg -pipe '
	./configure --enable-c-examples --enable-fortran-examples --prefix=`pwd`/local/librsb-profile --enable-matrix-types=blas --with-zlib \
		'--disable-shared' \
	"$@"
fi
if test -n "${RSB_WANT_CONFIGURE_ONLY}" ; then
	exit
fi
if ! grep '^CFLAGS =.*-pg' Makefile 2>&1 > /dev/null ; then
	make clean
	
	make -j `nproc`
	make -j `nproc` \
		"`grep ^CFLAGS  Makefile | sed  's/^CFLAGS *=/CFLAGS= -pg /g;s/ --coverage//g;s/-O0/-O2/g;s/ -ggdb / /g;s/ -g //g;s/ -lgcov//g;'`" \
		"`grep ^FCFLAGS Makefile | sed 's/^FCFLAGS *=/FCFLAGS=-pg /g;s/ --coverage//g;s/-O0/-O2/g;s/ -ggdb / /g;s/ -g //g;s/ -lgcov//g;'`" \
		"`grep ^LDFLAGS Makefile | sed 's/^LDFLAGS *=/LDFLAGS=-pg -static /g;s/ --coverage//g;s/-O0/-O2/g;s/ -ggdb / /g;s/ -g //g;s/ -lgcov//g;'`" \
		rsbench
else
	make -j `nproc`
	make -j `nproc` \
		"`grep ^LDFLAGS Makefile | sed 's/^LDFLAGS *=/LDFLAGS=-pg -static /g;s/ --coverage//g;s/-O0/-O2/g;s/ -ggdb / /g;s/ -g //g;s/ -lgcov//g;'`" \
		rsbench
fi
make install
make itests
rm -f gmon.out
./rsbench -Q 10 || exit
gprof rsbench | head -20
if test -f ~/src/scripts-ext/gprof2dot.py && which dot ; then
	gprof rsbench | ~/src/scripts-ext/gprof2dot.py | dot -Tps > profile.eps
fi
