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
# It converts complex/integer/pattern Matrix Market matrices in stdin to real Matrix Market matrices on stdout.
# needs sed, cat, bash
shopt -s nocasematch
p2r()
{
read line
if [[ $line =~ pattern ]]
then
	printf "%s\n" "${line/ pattern/ real}"
	while read line && [[ $line =~ ^% ]] ; do printf "%s\n" "$line" ; done
	printf "%s\n" "% matrix adapted from pattern to real by the $0 script, `date`" 
	printf "%s\n" "$line"
	while read line ; do printf "%s 1\n" "$line" ; done
else
	printf "%s\n" "$line"
	cat
fi
}

i2r()
{
read line
if [[ $line =~ integer ]]
then
	printf "%s\n" "${line/ integer/ real}"
	while read line && [[ $line =~ ^% ]] ; do printf "%s\n" "$line" ; done
	printf "%s\n" "% matrix adapted from integer to real by the $0 script, `date`" 
	printf "%s\n" "$line"
	while read line ; do printf "%s\n" "$line" ; done
else
	printf "%s\n" "$line"
	cat
fi
}

c2r()
{
read line
if [[ $line =~ complex ]]
then
	printf "%s\n" "${line/ complex/ real}"
	while read line && [[ $line =~ ^% ]] ; do printf "%s\n" "$line" ; done
	printf "%s\n" "% matrix adapted from complex to real by the $0 script, `date`" 
	printf "%s\n" "$line"
	cat | sed 's/\s[^ 	]*$//g'
else
	printf "%s\n" "$line"
	cat
fi
}

i2r | p2r | c2r
