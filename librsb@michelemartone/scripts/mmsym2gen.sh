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
sed 's/^\([0-9]\+\)\(\s\+\)\([0-9]\+\)\(\s\+\)\(.*\)/\1\2\3\4\5\n\3\2\1\4\5/g'  | uniq
