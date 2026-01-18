# Copyright (C) 2008-2021 Michele Martone
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

#
# Parallel Sparse Blas Tester Code.
#

source("./sbtg.m")

res=[rsb_octave_license("f")];
printf("%s",res);


res=[rsb_octave_license("f")];
#res=[res,"!\n"];
res=[res,"! Parallel Sparse BLAS fortran interface testing code\n"];
res=[res,"!\n"];
res=[res,findent,"module psb_mvsv_tester\n"];
res=[res,findent,"contains\n"];
res=[res,"!\n"];
printf("%s",res);
all_test("p","decl");
res="";
res=[res,"!\n"];
res=[res,findent,"end module psb_mvsv_tester\n"];
res=[res,"!\n"];
printf("%s",res);
