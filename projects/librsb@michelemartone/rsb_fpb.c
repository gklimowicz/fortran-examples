/*                                                                                                                            

Copyright (C) 2008-2019 Michele Martone

This file is part of librsb.

librsb is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

librsb is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with librsb; see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

*/
/* @cond INNERDOC  */
/**
 * @file
 * @author Michele Martone
 * @brief Floating point microbenchmarks.
 */

#include "rsb_common.h"
#include "rsb.h"
#include <limits.h>

#define RSB_DECLARE_FPB_F(FN) rsb_err_t FN(rsb__times_t times, size_t bs, rsb_bool_t aloud)	/* times to perform, buffer size */
#define RSB_M 1000000.0
typedef rsb_err_t (*rsb_fpb_fp_t)(rsb__times_t,size_t,rsb_bool_t);	/* floating point benchmark function pointer type */

#define RSB_DEFINE_FPB_F(FNAME,FPBEXPR,FPBNAME) \
static RSB_DECLARE_FPB_F(FNAME) \
{ \
	/* \ingroup gr_internals */ \
	rsb_nnz_idx_t N; rsb__times_t t,T=times; rsb_time_t dt; \
	rsb_type_t mtca[]=RSB_MATRIX_TYPE_CODES_ARRAY; rsb_char_t*mtcn[]=RSB_MATRIX_TYPES_ARRAY; \
	rsb_int ti=0; void *p=NULL; \
	p = rsb__calloc(bs); \
	if(!p) return RSB_ERR_ENOMEM; \
	for(ti=0;ti<RSB_IMPLEMENTED_TYPES;++ti) \
	{ \
		rsb_char_t *tn=mtcn[ti]; rsb_type_t typecode=mtca[ti]; \
		rsb_aligned_t alpha[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE]; \
 \
		/* RSB_INFO("bs:%zd, T:%zd\n",bs,T); */  \
 \
		N=bs/RSB_SIZEOF(typecode); \
		rsb__util_set_area_to_converted_integer(&alpha,typecode,1); \
		rsb__util_set_array_to_converted_integer(p,typecode,N,1,1); \
		dt = - rsb_time(); \
		for(t=0;t<T;++t) \
			FPBEXPR; \
		dt += rsb_time(); \
		if(aloud) \
			RSB_INFO("#op\ttype\tbs\tpasses\telements\tMOPS\n"), \
			RSB_INFO("%s\t%s\t%zd\t%zd\t%zd\t%f\n",FPBNAME,tn,bs,(size_t)times,(size_t)N,(((1.0/dt)*N)*T)/RSB_M); \
	} \
	RSB_CONDITIONAL_FREE(p); \
	return RSB_ERR_NO_ERROR; \
}

RSB_INTERNALS_COMMON_HEAD_DECLS

/* This is horrible and sad, sad, I know. */
RSB_DEFINE_FPB_F(rsb_fpb_add,rsb__util_vector_add(p,&alpha,typecode,N),"ADD")
RSB_DEFINE_FPB_F(rsb_fpb_sum,rsb__util_vector_sum(&alpha,p,typecode,N),"SUM")
RSB_DEFINE_FPB_F(rsb_fpb_mul,rsb__cblas_Xscal(typecode,N,&alpha,p,1),"MUL")
RSB_DEFINE_FPB_F(rsb_fpb_neg,rsb__util_do_negate(p,typecode,N),"NEG")
RSB_DEFINE_FPB_F(rsb_fpb_inc,rsb__vector_increase_by_one(p,typecode,N),"INC")
RSB_DEFINE_FPB_F(rsb_fpb_sqr,rsb__util_vector_sqrt(p,typecode,N),"SQR")
RSB_DEFINE_FPB_F(rsb_fpb_div,rsb__util_vector_div(p,&alpha,typecode,N),"DIV")

rsb_err_t rsb__fp_benchmark(void)
{
	/**
	 * Benchmark the floating point units.
	 * You should call rsb_lib_init(RSB_NULL_INIT_OPTIONS) before.
	 *
	 * may add: pow, log
	 * benchmark for ops in the L1
	 * strided ops
	 *
	 */
	rsb_fpb_fp_t fpba[]={rsb_fpb_add,rsb_fpb_sum,rsb_fpb_mul,rsb_fpb_neg,rsb_fpb_inc,rsb_fpb_sqr,rsb_fpb_div};
	rsb_int i;
//	size_t bs = rsb__get_lastlevel_c_size();
//	size_t bs = rsb__get_first_level_c_size();
	rsb_int_t cln = rsb__get_cache_levels_num(),cli;
	rsb__times_t times,mtimes = RSB_MEMSCAN_MIN_TIMES,Mtimes = RSB_MEMSCAN_MAX_TIMES;
	rsb_time_t mt = rsb__getenv_real_t("RSB_FPBENCH_MULTITYPE_TIME", RSB_BENCHMARK_MIN_SECONDS);
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	for(cli=1;cli<=cln;++cli)
	for(i=0;i<sizeof(fpba)/sizeof(rsb_fpb_fp_t);++i)
	{
		size_t bs = rsb__get_lnc_size(cli);

		if(!bs)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
		RSB_INFO("#probing for an iterations count (to a total of %f s) .. \n",mt);
		for(times=mtimes;times<Mtimes;times*=2)
		{
			rsb_bool_t aloud = RSB_BOOL_FALSE;
			rsb_time_t dt;

			dt = - rsb_time();
			fpba[i](times,bs,aloud);
			dt += rsb_time();
			if(dt>mt)
			{
				aloud = RSB_BOOL_TRUE,
				fpba[i](times,bs,aloud);
				break;	/* break the inner loop, go for another benchmark */
			}
		}
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
