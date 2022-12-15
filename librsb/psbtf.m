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

#
# Parallel Sparse Blas Tester Code.
#

source("./sbtg.m")

res=rsb_octave_license("f");
res=[res,"!\n"];
res=[res,"! Parallel Sparse BLAS fortran interface testing code\n"];
res=[res,"!\n"];
res=[res,"!\n"];
res=[res,"!> @cond INNERDOC\n"];
res=[res,findent,"PROGRAM main\n\n"];
#res=[res,findent,"PROGRAM main\n\n",findent,"INTEGER :: res,passed=0,failed=0;\n"];
res=[res,findent,"USE psb_base_mod\n"];
res=[res,findent,"USE psb_mvsv_tester\n"];
res=[res,findent,"IMPLICIT NONE\n"];
res=[res,findent,"INTEGER :: res,passed=0,failed=0,fi=0\n"];
res=[res,findent,"INTEGER            :: ictxt, iam=-1, np=-1\n"];
res=[res,findent,"CHARACTER(LEN=psb_fidasize_) :: afmt\n"];
#res=sprintf("%s%s",res,findent,"INTEGER            :: ictxt, iam, np\n");
res=sprintf("%s%s",res,findent,"CALL psb_init(ictxt)\n");
res=sprintf("%s%s",res,findent,"CALL psb_info(ictxt,iam,np)\n");
res=[res,findent,"IF(iam<0)THEN\n"];
res=[res,findent,findent,"GOTO 9999\n"];
res=[res,findent,"ENDIF\n"];
res=[res,findent,"DO fi=1,2\n"];
res=[res,findent,"IF(fi.EQ.1)afmt=psb_csr_afmt_\n"];
res=[res,findent,"IF(fi.EQ.2)afmt=psb_coo_afmt_\n"];
printf("%s",res);
all_test("p","CALL");
res=["" ,findent,"ENDDO\n"];
res=[res,"9999",findent,"CONTINUE\n"];
res=[res,findent,"PRINT *,\"PASSED:\",passed\n"];
res=[res,findent,"PRINT *,\"FAILED:\",failed\n"];
res=sprintf("%s%s",res,findent,"CALL psb_exit(ictxt)\n");
res=[res,findent,"END PROGRAM\n"];
res=[res,"!> @endcond\n"];
res=sprintf("%s%s",res,"\n");
res=sprintf("%s%s",res,"\n\n");
res=sprintf("%s%s",res,"");
printf("%s",res);
