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
# It extracts headers out of MatrixMarket matrix files

if=$1
od=$2

fail() { echo "[!] $@" ; exit ; }
if   test x = x"$if" ; then fail "" ; fi
if   test x = x"$od" ; then fail "" ; fi
if ! test -f   "$if" ; then fail "" ; fi
if ! test -d   "$od" ; then fail "" ; fi


bf=`basename $if`
of=$od/$bf.head
echo "$bf > $of"  
grep ^% -C 1 $if > $of


