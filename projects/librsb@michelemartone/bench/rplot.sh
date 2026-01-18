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
# Recursion plot script.

e=./rsbench

make $e || exit

md=~/matrices
echo $md

for f in $md/*.mtx ; do
	#echo $f
	hn=`hostname`
	d=`date +%Y%m%d`
	# FIXME : this substitution will not work on SP5's bash
	m=${f/.mtx/}
	m="`basename $m`"
	for p in "-E" "-C" " " ; do
	for h in "-H" " " ; do
	for l in "-L" " " ; do
	for D in "-D" " " ; do
		# -T will force lower triangular
		flags="$p $h $l $D"
		flagsn=${flags// }
		ofn="recplot-$hn-$d-$m-$flagsn.eps"
		echo "$e --plot-matrix -f $f -aRzd -r1 -c1 -Fbr -T $flags  > $ofn"
		$e --plot-matrix -f $f -aRzd -r1 -c1 -Fbr -T $flags  >  $ofn
	done
	done
	done
	done
done

