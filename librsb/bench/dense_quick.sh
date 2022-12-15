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

export RSBDB_MAXMEM=${RSBDB_MAXMEM:-}
if test -z "${RSBDB_MAXMEM}" ; then RSBDB_MAXMEM=`cat /proc/meminfo  | grep MemTotal  | sed 's/MemTotal: *\([0-9]*\)\s*kB/\1/g'` ; RSBDB_MAXMEM=$((${RSBDB_MAXMEM}/1024)) ;  fi
if test -z "${RSBDB_MAXMEM}" ; then RSBDB_MAXMEM=`cat /proc/meminfo  | grep MemTotal  | sed 's/MemTotal: *\([0-9]*\)\s*MB/\1/g'` ; RSBDB_MAXMEM=${RSBDB_MAXMEM} ;  fi
export RSBDB_MF='8'
export RSBDB_BPNZ=${RSBDB_BPNZ:-32}
#echo ${RSBDB_MAXMEM}
maxdim=`perl -e "print int(sqrt(((${RSBDB_MAXMEM}*1024*1024)/${RSBDB_BPNZ})/${RSBDB_MF}))"`
dimv=4;
dims=`while test $dimv -lt $maxdim ; do echo $dimv $(( ($dimv*99)/70 )) ;dimv=$((dimv*2)); done` 
#echo $dims
up=12
export RSBDB_DIMS="$dims"
export RSBDB_SPACINGS='1'
export RSBDB_INCXS='1'
RSBDB_TIMES=${RSBDB_TIMES:-30}
#export| grep RSBDB_
bench/dense.sh \ 
