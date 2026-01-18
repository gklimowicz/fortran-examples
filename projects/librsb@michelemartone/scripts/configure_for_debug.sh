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

# This script is an example of configuring the library for debug purposes.
set -e
set -x

export CC=${CC:=`which icc   || which gcc      || which clang || false`}
export FC=${FC:=`which ifort || which gfortran || which flang || true `}

GCC_WERR_CFLAGS="-Werror=discarded-qualifiers -Werror=incompatible-pointer-types -Werror=unused-value -Werror=parentheses -pedantic -Werror=implicit-function-declaration -Werror=strict-prototypes"
GCC_WERR_CFLAGS+=" -pedantic-errors"
GCC_WERR_CFLAGS+=" -Werror=declaration-after-statement"
GCC_WERR_CFLAGS+=" -Werror=unused-function"
GCC_WERR_CFLAGS+=" -Werror=uninitialized"
GCC_WERR_CFLAGS+=" -Werror=maybe-uninitialized"

for w in $GCC_WERR_CFLAGS; do # remove unsupported flags
	if ${CC} --help=warning | grep -- "${w/-*=/}" ; then
		echo "good, has ${w}"
	else
		GCC_WERR_CFLAGS="${GCC_WERR_CFLAGS/$w/}"
	fi
done

CFLAGS_FOR_CONFIGURE="-Wno-error=strict-prototypes -Wno-error=builtin-declaration-mismatch" # these flags seem to strict for configure: strcpy may fail from being detected
for w in $CFLAGS_FOR_CONFIGURE; do # remove unsupported flags
	if ${CC} --help=warning | grep -- "${w/-Wno-error=/}" ; then
		echo "good, has ${w}"
	else
		CFLAGS_FOR_CONFIGURE="${CFLAGS_FOR_CONFIGURE/$w/}"
	fi
done

if test -z "${RSB_WANT_CONFIGURE_ONLY}" ; then
	./autogen.sh
fi
./configure '--enable-allocator-wrapper'  '--enable-debug' \
       	FC=${FC} \
       	CC=${CC} \
	'CFLAGS=-O0 -fPIC -ggdb -pipe -Wall -Wredundant-decls -Wno-switch -Wdisabled-optimization -Wdeclaration-after-statement   -Wpointer-arith -Wstrict-prototypes '"${GCC_WERR_CFLAGS} ${CFLAGS_FOR_CONFIGURE} ${CFLAGS}" \
	"FCFLAGS=-O0 -fPIC -ggdb ${FCFLAGS}"	\
	"CXXFLAGS=-O0 -fPIC -ggdb ${CXXFLAGS}"	\
	'--enable-librsb-stats' \
	'--enable-rsb-num-threads' \
	'--enable-zero-division-checks-on-solve' \
	'--disable-shared' \
	'--enable-internal-headers-install' \
	'--enable-internals-error-verbosity=1'  \
	'--enable-interface-error-verbosity=2'  \
	'--enable-io-level=7'			\
	'--enable-debug-getenvs'		\
	${FC:+--enable-fortran-module-install}	\
	'--prefix'=`pwd`/local/librsb-debug	\
	'--with-baselib-cflags= -Wall -Werror=all -Werror=pedantic -Werror=unused-variable    -Werror=unused-but-set-variable ' \
	'--with-kernels-cflags= -Wall -Werror=all -Werror=pedantic -Wno-error=unused-variable -Wno-error=unused-but-set-variable ' \
	\
	'--with-zlib' $(which doxygen >/dev/null && echo --enable-doc-build;)	\
	"$@"

sed -i "s/$CFLAGS_FOR_CONFIGURE/${CFLAGS_FOR_CONFIGURE//-Wno-/   -W}/g" Makefile */Makefile config.log rsb-config.h # restore stricter flags which could otherwise have broken ./configure
