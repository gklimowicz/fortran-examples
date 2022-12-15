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

res=rsb_octave_license("c");
res=sprintf("%s%s",res,"/* @cond INNERDOC */\n");
res=sprintf("%s%s",res,rsb_octave_doc_c_header);
res=sprintf("%s%s",res,"\
#include <stdio.h>\n\
#ifdef SBTC_USE_RSB_H\n\
#include <rsb.h>\n\
#endif /* SBTC_USE_RSB_H */\n\
#include <complex.h>\n\
#ifdef RSB_RSB_H_INCLUDED\n\
#include \"rsb_internals.h\"\n\
#define RSB_BLAS_SUPPORT_EMPTY 1\n\
#define RSB_BLAS_SUPPORTED_TYPE(T) ((errval = rsb__BLAS_is_type_supported(T)) != RSB_ERR_UNSUPPORTED_TYPE) \n\
#endif /* RSB_RSB_H_INCLUDED */\n\
#ifndef RSB_RSB_H_INCLUDED\n\
#include <blas_sparse.h>\n\
#define RSB_PROGRAM_SUCCESS 0\n\
#define RSB_PROGRAM_ERROR (-1)\n\
#define RSB_ERR_NO_ERROR 0\n\
#define RSB_ERROR printf\n\
#define RSB_WITH_SPARSE_BLAS_INTERFACE 1\n\
#define RSB_BLAS_SUPPORTED_TYPE(T) 1\n\
/*#define rsb_err_t int*/\n\
#define RSB_ERR_UNSUPPORTED_TYPE 0x004\n\
/*#define rsb_nnz_idx_t int*/\n\
/*#define rsb_coo_idx_t int*/\n\
#define RSB_BLAS_ERROR -1\n\
#define RSB_BLAS_NO_ERROR 0\n\
#define RSB_BLAS_SUPPORT_EMPTY 0\n\
#define rsb__debug_print_vectors_diff(A1,A2,A3,A4,A5,A6,A7) RSB_ERR_NO_ERROR\n\
int rsb__do_are_same(void*v1_,void*v2_,int n,int typecode,int s1,int s2){ char*v1=(char*)v1_,*v2=(char*)v2_; int vi,bi,bs; switch(typecode){case('S'):bs=4;break;case('C'): case('D'):bs=8;break;case('Z'):bs=16;break;default: return RSB_ERR_NO_ERROR; } for(vi=0;vi< n;++vi) for(bi=0;bi<bs;++bi) if(v1[vi*bs*s1+bi] != v2[vi*bs*s2+bi]) return RSB_BLAS_ERROR; return RSB_ERR_NO_ERROR;}\n\
#endif /* RSB_RSB_H_INCLUDED */\n\
const int rsb_sbtc_print_vec(const void*v,int n,int typecode){ float*fv=(float*)v; double*dv=(double*)v; int vi,fl=1,fi; if(typecode=='C' || typecode=='Z')fl=2; if(typecode=='S' || typecode=='C')for(vi=0;vi<n;++vi){for(fi=0;fi<fl;++fi)printf(\"%f \" ,fv[vi*fl+fi]);printf(\"\\n\");} if(typecode=='D' || typecode=='Z')for(vi=0;vi<n;++vi){for(fi=0;fi<fl;++fi)printf(\" %lf\",dv[vi*fl+fi]);printf(\"\\n\");} ; return RSB_ERR_NO_ERROR;}\n\
");
printf("%s",res);
#quit
res=sprintf("%s%s","" ,"#if RSB_WITH_SPARSE_BLAS_INTERFACE\n");
res=sprintf("%s%s",res,all_test("c","decl"));
res=sprintf("%s%s",res,"#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */\n");
res=sprintf("%s%s",res,"\nint main(void)\n {\n int errval;int passed=0,failed=0,skipped=0;\n");
res=sprintf("%s%s",res,"#ifdef RSB_RSB_H_INCLUDED\n");
res=sprintf("%s%s",res,"\tif( rsb_lib_init(RSB_NULL_INIT_OPTIONS) != RSB_ERR_NO_ERROR) goto err;\n");
#res=sprintf("%s%s",res,"\tif( rsb_blas_mini_tester() != RSB_ERR_NO_ERROR) goto err;\n");
res=sprintf("%s%s",res,"#endif /* RSB_RSB_H_INCLUDED */\n");
res=sprintf("%s%s",res,"#if RSB_WITH_SPARSE_BLAS_INTERFACE\n");
printf("%s",res);
all_test("c","CALL");
res=sprintf("%s","printf(\"	PASSED:%d\\n	SKIPPED:%d (tests for BLAS types/matrix types excluded at configure/make time are skipped)\\n	FAILED:%d (if any check failed, this may indicate a bug)\\n\",passed,skipped,failed);\n");
res=sprintf("%s%s",res,"#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */\n");
#res=sprintf("%s%s",res,"\n	if(failed) goto err;\n");
#res=sprintf("%s%s",res,"\n	return RSB_PROGRAM_SUCCESS;\n err: return RSB_PROGRAM_ERROR;\n");
res=sprintf("%s%s",res,"\n	if (rsb_perror(NULL,rsb_lib_exit( RSB_NULL_EXIT_OPTIONS )) != RSB_ERR_NO_ERROR) { failed=-1; }\n");
res=sprintf("%s%s%s",res,"\n	if (failed)\n","err:		return RSB_PROGRAM_ERROR;\n	else\n		return RSB_PROGRAM_SUCCESS;\n");
res=sprintf("%s%s",res,"\n}\n");
res=sprintf("%s%s",res,"");
res=sprintf("%s%s",res,"/* @endcond */\n");
printf("%s",res);

