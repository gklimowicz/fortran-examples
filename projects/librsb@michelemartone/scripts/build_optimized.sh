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
which m4
which nproc
if test -z "${RSB_WANT_CONFIGURE_ONLY}" ; then
	./autogen.sh
fi
./configure FCFLAGS=-O0 CFLAGS=-O0
make distclean
export WANT_GCC_MKL='' WANT_ICC_MKL=''
export LIBS
export CFLAGS
export CC=${CC:=`which icc   || which gcc      || which clang || false`}
export FC=${FC:=`which ifort || which gfortran || which flang || true `}
CFLAGS=-Ofast
CFLAGS+=' -march=native'
CFLAGS+=' -mtune=native'
export FCFLAGS="${CFLAGS}"
export FCFLAGS="${FCFLAGS}" CFLAGS="${CFLAGS}"
if [[ "`basename ${CC}`" =~ icc ]] ; then
	test -n "${MKLROOT}" && export WANT_ICC_MKL=1
else
	CFLAGS+=' -pipe'
	#${CC} -c ${CFLAGS} --help=optimizers > O3-opts.txt
	test -n "${MKLROOT}" && export WANT_GCC_MKL='1'
fi
test -z "${CXXFLAGS}" && export CXXFLAGS="${CFLAGS}"
./configure --prefix=`pwd`/local/librsb-optimized --enable-matrix-types=blas --with-zlib --disable-c-examples --disable-fortran-examples \
	${WANT_GCC_MKL:+ --with-mkl="-static -rpath ${MKLROOT}/lib/intel64 -L${MKLROOT}/lib/intel64 -fopenmp -lpthread -Wl,--start-group,-lmkl_intel_lp64,-lmkl_gnu_thread,-lmkl_core,--end-group" --with-mkl-include=${MKLROOT}/include/ }	\
	${WANT_ICC_MKL:+ --with-mkl="-static -rpath ${MKLROOT}/lib/intel64 -L${MKLROOT}/lib/intel64 -fopenmp -lpthread -Wl,--start-group,-lmkl_intel_lp64,-lmkl_intel_thread,-lmkl_core,--end-group" --with-mkl-include=${MKLROOT}/include/ }	\
	"$@"
if test -z "${RSB_WANT_CONFIGURE_ONLY}" ; then
	make clean
	make cleanall
	make -j `nproc`
	make qqtests
	make install
	make itests
	make tests # needed for dist
	make dist
fi
