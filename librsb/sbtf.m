# Copyright (C) 2008-2019 Michele Martone
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
# Sparse Blas Tester Code.
#

source("./sbtg.m")

res=rsb_octave_license("f");
res=[res,"!\n"];
res=[res,rsb_octave_doc_f_header];
res=[res,"!\n"];
res=[res,"! Sparse BLAS fortran interface testing code\n"];
res=[res,"!\n"];
res=[res,"! FIXME: missing library initialization!\n"];
res=[res,"! FIXME: using zero based indices is only partially supprted!\n"];
res=[res,"!\n"];
printf("%s",res);
all_test("f","decl");
res=sprintf("%s%s","","!\n");
res=[res,findent,"PROGRAM main\n\n",findent,"USE rsb\n",findent,"INTEGER :: passed=0,failed=0,errval=0\n"];
res=sprintf("%s%s",res,findent,"info = rsb_lib_init(RSB_NULL_INIT_OPTIONS)\n");
res=[res,findent,"IF(info.NE.0)THEN\n"];
res=[res,findent,"STOP 1\n"];
res=[res,findent,"ENDIF\n"];
printf("%s",res);
all_test("f","CALL");
res=["" ,findent,"PRINT *,\"PASSED:\",passed\n"];
res=[res,findent,"PRINT *,\"FAILED:\",failed\n"];
res=[res,findent,"info = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS)\n"];
#res=[res,findent,"IF(info.NE.0)THEN\n"];
#res=[res,findent,"STOP 1\n"];
#res=[res,findent,"ENDIF\n"];
#res=[res,findent,"IF(failed.NE.0)THEN\n"];
#res=[res,findent,"STOP 1\n"];
#res=[res,findent,"ENDIF\n"];
res=[res,findent,"IF(failed.NE.0.OR.info.NE.0)THEN\n"];
res=[res,findent,"STOP 1\n"];
res=[res,findent,"ENDIF\n"];
res=[res,findent,"END PROGRAM main\n"];
res=sprintf("%s%s",res,"\n");
res=sprintf("%s%s",res,"\n\n");
res=sprintf("%s%s",res,rsb_octave_doc_f_footer);
printf("%s",res);

