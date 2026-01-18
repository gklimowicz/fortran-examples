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

r=$1
c=$2
rs=$3
cs=$4

help()
{
	echo "#usage examples:" 
	echo "#dense , r=3 c=3 rs=1 cs=1 ./$0"
	echo "./$0 3 3"
	echo "#strided, r=3 c=3 rs=1 cs=1 ./$0"
	echo "./$0 3 3 1 1 "
	exit -1
}


[ -z "$r" ] && help
[ -z "$c" ] && help

echo "%%MatrixMarket matrix coordinate real general"
echo "$r $c $(((r/rs)*(c/cs)))"
for((i=1;i<=$r;i+=rs))
do
for((j=1;j<=$c;j+=cs))
do
	#echo $i $j 1
	echo $i $j $i.$j
done
done

