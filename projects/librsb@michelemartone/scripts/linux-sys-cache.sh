#!/bin/sh
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

# This script tries to produce information about caches.
# It may fail easily, because the Linux sys interface is volatile.
# It is used as an alternative way to get cache info from the system on pre-production systems, where the sysconf interface does not work effectively.

ncpu=`ls -d /sys/devices/system/cpu/cpu* | wc -l`

if test "x$ncpu" = x ; then
	echo "" 1>&2
       	exit
fi

if test ! -d "/sys/devices/system/cpu/cpu0/cache" ; then exit ; fi 
ncache=`ls -d /sys/devices/system/cpu/cpu0/cache/index* | wc -l`

if test "x$ncache" = x ; then
	if test x"`which lstopo 2>&1 | wc -l`" != x ; then
		# lstopo calling convention as of version 1.11.12
		CS=''
		CL=4/64/
		L1=`lstopo --of console | grep L1 | head -1 | sed  's/^\(.*L1d\? L.0 (\)\([0-9]\+.B\).*$/\2/g'`
		if test -n "$L1" ; then CS="L1:$CL$L1"; fi
		L2=`lstopo --of console | grep L2 | head -1 | sed  's/^\(.*L2 L.0 (\)\([0-9]\+.B\).*$/\2/g'`
		if test -n "$L2" ; then CS="L2:$CL$L2,$CS"; fi
		L3=`lstopo --of console | grep L3 | head -1 | sed  's/^\(.*L3 L.0 (\)\([0-9]\+.B\).*$/\2/g'`
		if test -n "$L3" ; then CS="L3:$CL$L3,$CS"; fi
		echo "$CS"
	else
		echo "" 1>&2
	fi
       	exit
fi

cacheinfo=""

#for n in `seq 1 $ncache`
for cf in /sys/devices/system/cpu/cpu0/cache/index*
do
	if test ! -d "/sys/devices/system/cpu/cpu0/cache" ; then continue ; fi 
	#echo $cf
	tp=`cat $cf/type`
	if test "x$tp" = x"Data" || test "x$tp" = x"Unified"  ; then
		sz=`cat $cf/size`
		as=`cat $cf/ways_of_associativity`
		ls=`if test -f $cf/coherency_line_size ; then cat $cf/coherency_line_size; else echo 64; fi` # bogus also ok
		lv=`cat $cf/level`
		if test "x$lv" = x"1" ; then
			cacheinfo="L$lv:$as/$ls/$sz" # a,b,c parameters
		else
			cacheinfo="L$lv:$as/$ls/$sz,$cacheinfo" # a,b,c parameters
		fi
	fi
done

#echo "ncpu:$ncpu ncache:$ncache cacheinfo:$cacheinfo"
echo "$cacheinfo"

# Examples:
# L3:16/64/12288K;L2:8/64/256K;L1:8/64/32K;
# L2:4/64/512K;L1:8/64/32K;

