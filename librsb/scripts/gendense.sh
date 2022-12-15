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
b=$3

help()
{
	echo "#usage examples:" 
	echo "#dense, r=3 c=3 ./$0"
	echo "./$0 3 3"
	echo "#banded, r=3 c=3 b=2 ./$0"
	echo "./$0 3 3 2"
	exit -1
}

[ -z "$r" ] && help
[ -z "$c" ] && help

# b should be less than r


if test -z "$b" ;
then
	echo "%%MatrixMarket matrix coordinate real general"
	echo "$r $c $((r*c))"
	for((i=1;i<=$r;++i))
	do
	for((j=1;j<=$c;++j))
	do
		echo $i $j 1
	done
	done
else
	if [ $((r!=c)) == 1 ]
	then
		echo "banded matrices should be square!"
		exit -1
	fi

	if [ $((b>=r)) == 1 ]
	then
		echo "band width cannot exceed matrix size-1!"
		exit -1
	fi

	echo "%%MatrixMarket matrix coordinate real general"

	echo "$r $c $(((2*b+1)*r-b*(b+1)))"
	for((i=1;i<=$r;++i))
	do
	l=$((i-b<1?1:i-b))
	u=$((i+b>c?c:i+b))
	for((j=l;j<=u;++j))
	do
		echo $i $j 1
	done
	done
fi


