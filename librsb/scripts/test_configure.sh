#!/bin/bash
#
# Copyright (C) 2008-2022 Michele Martone
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

set -e
set -x
XFAIL=0
function xfail() { echo XFAIL=$((++XFAIL)); }
test -f ./configure
rm -fR   /dev/shm/librsb_configure_tests
mkdir -p /dev/shm/librsb_configure_tests
cd       /dev/shm/librsb_configure_tests
make dist VERSION=tests -C ${OLDPWD} && tar xzf ${OLDPWD}/librsb-tests.tar.gz --transform 's/librsb-tests.//g'
if test -x librsbpp/configure -a -x rsblib/configure -a -x rsbtest/configure ; then HAVESUBDIRS=1; else unset HAVESUBDIRS; fi

export CXX=' ' # same as always adding "--without-librsbpp --without-rsblib --without-rsbtest" -- here purpose is speed up tests

test `./configure --help           | wc -l` -gt 230
test `./configure --help=recursive | wc -l` -gt 150
test `./configure --help=short     | wc -l` -gt 150

export CFLAGS=-O0\ -pipe
export CXXFLAGS=-O0\ -pipe

./configure
grep 'CFLAGS.*:.*-O0'   config.log
grep 'CXXFLAGS.*:.*-O0' config.log

if test x${HAVESUBDIRS} != x ; then
	./configure CXX=' '
	grep 'Build librsbpp.*no.$' config.log
	grep 'Build rsblib.*no.$' config.log
	grep 'Build rsbtest.*no.$' config.log
fi

./configure CFLAGS=-O1 CXXFLAGS=-O1
grep 'CFLAGS.*:.*-O1'   config.log
grep 'CXXFLAGS.*:.*-O1' config.log

./configure --enable-openmp 'ac_cv_prog_c_openmp=-DMY_OMP_FLAGS'
grep 'OPENMP_CFLAGS.*:.*-DMY_OMP_FLAGS' config.log

./configure OPENMP_CFLAGS=-DMY_OMP_FLAGS
grep 'OPENMP_CFLAGS.*:.*-DMY_OMP_FLAGS' config.log

./configure FC=' '     | grep 'FC.*disable.*: *$'
./configure FC=' '     | grep 'FC.*disable.*: *$'
./configure FC='false' | grep 'FC.*disable.*: false*$'

unset RSB_USER_SET_MEM_HIERARCHY_INFO
./configure --prefix=`pwd`/cfgt
grep -F "prefix='`pwd`/cfgt'" config.log 
grep -E "^RSB_USER_SET_MEM_HIERARCHY_INFO=''" config.log

./configure --with-c99-flag
grep -E 'CFLAGS=.*-std=c99'   config.log

./configure --without-c99-flag
grep -E 'CFLAGS=.*-std=c99' config.log && xfail

./configure --enable-openmp
grep -E 'CFLAGS.*-fopenmp' config.log
grep -E 'FCFLAGS.*-fopenmp' config.log

./configure --disable-openmp
( ! grep -E 'CFLAGS.*-fopenmp' config.log ; )
( ! grep -E 'FCFLAGS.*-fopenmp' config.log; )

( ! ./configure                           --non-existing-option ; )
( ! ./configure  --enable-option-checking --non-existing-option ; )
( ! ./configure --disable-option-checking --non-existing-option ; )
./configure                          --with-non-existing-option
./configure --enable-option-checking --with-non-existing-option
./configure --enable-option-checking --enable-non-existing-option
./configure                          --enable-non-existing-option
./configure  --help | grep -n oski   && xfail # removed in 1.3
./configure  --help | grep -n papi   && xfail # removed in 1.3
./configure  --help | grep -n likwid && xfail # removed in 1.3
./configure  --help | grep -- --enable-long-indices
./configure  --help | grep -- --enable-debug-getenvs
./configure  --help | grep -- --enable-fortran-module-install
./configure  --help | grep -- --disable-sparse-blas-interface
./configure  --help | grep -- --disable-fortran-examples
./configure  --help | grep -- --with-librsbpp
./configure  --help | grep -- --with-rsblib
./configure  --help | grep -- --with-rsbtest
./configure  --help | grep -- --enable-matrix-types
./configure  --help | grep -- --enable-allocator-wrapper
./configure  --help | grep -- --enable-doc-build
./configure  --help | grep -- --enable-debug
./configure  --help | grep -- --enable-debug-getenvs

./configure --with-mkl=-L`pwd`
grep -E "^RSB_RSBENCH_LIBS.*-L" config.log | grep -F "`pwd`"

./configure --with-mkl=-lmy_m
grep -E "^RSB_RSBENCH_LIBS.*-lmy_m" config.log 

./configure --with-mkl="-L`pwd` -lmy_m"
grep -E "^RSB_RSBENCH_LIBS.*-lmy_m" config.log | grep -F "`pwd`"

if which g++ && test x${HAVESUBDIRS} == x ; then
	( ! ./configure              --with-mkl-include=`pwd`/my_mkl         | grep -F "`pwd`"; )
else
	(   ./configure              --with-mkl-include=`pwd`/my_mkl CXX=g++ | grep -F "`pwd`"; )
fi

./configure --with-mkl=-lmy_m --with-mkl-include=`pwd`/my_mkl | grep -F "`pwd`"/my_mkl | grep -- -I;

./configure --with-memhinfo=bbb --with-memhinfo=aaa
grep -E "^RSB_USER_SET_MEM_HIERARCHY_INFO='aaa'" config.log

./configure --with-memhinfo=aaa RSB_USER_SET_MEM_HIERARCHY_INFO=bbb
grep -E "^RSB_USER_SET_MEM_HIERARCHY_INFO='aaa'" config.log

RSB_USER_SET_MEM_HIERARCHY_INFO=bbb ./configure --with-memhinfo=aaa 
grep -E "^RSB_USER_SET_MEM_HIERARCHY_INFO='aaa'" config.log

RSB_USER_SET_MEM_HIERARCHY_INFO=bbb ./configure
grep -E "^RSB_USER_SET_MEM_HIERARCHY_INFO=''" config.log

./configure RSB_USER_SET_MEM_HIERARCHY_INFO='L2:4/64/512K,L1:8/64/24K'
grep -E "^RSB_USER_SET_MEM_HIERARCHY_INFO=''" config.log

if test "x$MKLROOT" != x ; then
    if test -n "`which icc`" ; then
        ./configure --with-mkl --enable-openmp CC=gcc;
        grep 'RSB_RSBENCH_LIBS.*=.*mkl_gnu_thread'              Makefile;
        ( ! grep 'RSB_RSBENCH_LIBS.*=.*mkl_intel_thread'        Makefile )
        ( ! grep 'RSB_RSBENCH_LIBS.*=.*mkl_sequential'          Makefile )
  
        ./configure --with-mkl --disable-openmp CC=gcc;
        grep    'RSB_RSBENCH_LIBS.*=.*mkl_sequential'           Makefile;
        ( ! grep 'RSB_RSBENCH_LIBS.*=.*mkl_gnu_thread'          Makefile )
        ( ! grep 'RSB_RSBENCH_LIBS.*=.*mkl_intel_thread'        Makefile )
    fi
  
    if test -n "`which icc`" ; then
        ./configure --with-mkl --enable-openmp CC=icc;
        grep 'RSB_RSBENCH_LIBS.*=.*mkl_intel_thread'              Makefile;
        ( ! grep 'RSB_RSBENCH_LIBS.*=.*mkl_gnu_thread'        Makefile )
        ( ! grep 'RSB_RSBENCH_LIBS.*=.*mkl_sequential'          Makefile )
    fi
fi

echo XFAIL:$((XFAIL))
true

#cleanup can wait (temp files are good for inspection)
#cd -
#rm -fR       /dev/shm/librsb_configure_tests
