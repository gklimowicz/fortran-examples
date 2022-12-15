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
set -o pipefail
# This script is intended for the librsb developer usage.
if test x"${srcdir}" = x ; then srcdir=. ; fi
if cat ${srcdir}/examples/*.c | grep '^.\{71,\}' ; then echo 'Some source code exceeds 71 chars!'; grep -n '^.\{71,\}' ${srcdir}/examples/*.c ; exit 255 ; else true ; fi
if cat ${srcdir}/README       | grep '^[^	].\{80,\}' ; then echo 'Some source code exceeds 81 chars!'; grep -n '^[^	].\{80,\}' ${srcdir}/README       ; exit 255 ; else true ; fi
test `${srcdir}/rsbench -h | wc -l` -ge 61
test `${srcdir}/rsbench -h | wc -c` -ge 1966
test `${srcdir}/rsbench -oa -Ob -h | wc -l` -ge 157
test `${srcdir}/rsbench -oa -Ob -h | wc -c` -ge 4600
exit 0;
