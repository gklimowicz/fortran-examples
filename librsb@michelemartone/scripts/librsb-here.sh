#!/bin/bash
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

RSBB_URL=svn+ssh://user@host/repository/trunk
RSBB_SVN=${RSBB_SVN:-svn}
RSBB_SVN_OPTIONS=${RSBB_SVN_OPTIONS:-}
RSBB_RUN_AUTOGEN=1
RSBB_RUN_CONFIGURE=1
RSBB_CONFIGURE_ARGS=${RSBB_CONFIGURE_ARGS:-}
RSBB_CONFIGURE_ADD=${RSBB_CONFIGURE_ADD:-}
RSBB_RUN_MAKE_CLEAN=1
RSBB_RUN_MAKE=1
RSBB_RUN_MAKE_QTESTS=1
RSBB_RUN_MAKE_INSTALL=1
#
env | grep ^RSBB
#
RSD=`pwd`/librsb-src
RBD=`pwd`/librsb-usr
#
$RSBB_SVN ${RSBB_SVN_OPTIONS} --force co ${RSBB_URL} $RSD || exit -1
cd $RSD || exit -1
if test configure.ac -nt configure ; then sh autogen.sh || exit -1 ; fi
	if test x$RSBB_RUN_AUTOGEN = x1 ; then
		sh autogen.sh || exit -1
	fi
	if test x$RSBB_RUN_CONFIGURE = x1 ; then
		$ECHO ./configure CC=${CC} CFLAGS="${CFLAGS} -fPIC" ${RSBB_CONFIGURE_ARGS} ${RSBB_CONFIGURE_ADD} --prefix=${RBD}
		#touch config.h -r is.h
	fi
	if test x$RSBB_RUN_MAKE_CLEAN = x1 ; then
		$ECHO make clean || exit -1
	fi
	if test x$RSBB_RUN_MAKE = x1 ; then
		$ECHO make || exit -1
	fi
	if test x$RSBB_RUN_MAKE_QTESTS = x1 ; then
		$ECHO make qtests || exit -1
	fi
	if test x$RSBB_RUN_MAKE_INSTALL = x1 ; then
		$ECHO make install || exit -1
	fi
cd -
