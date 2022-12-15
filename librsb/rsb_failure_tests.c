/*

Copyright (C) 2008-2021 Michele Martone

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
 * @brief
 *
 * failure testing code.
 * \internal
 *
 * */
#include "rsb_common.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

#define RSB_WANT_FAILURE_TESTER 1
#define RSB_FAILURE_NOTICE(STMT) STMT;RSB_INFO("Injecting failure:\n%s\n",#STMT)
#define RSB_VARIATION_NOTICE(STMT) STMT;RSB_INFO("Injecting variation:\n%s\n",#STMT)
#define RSB_FT_FFL_PRINTF printf("In %s located in %20s:%d :\n",__func__,__FILE__,__LINE__)
#if 0
static rsb_err_t rsb_do_meminfo(void)
{
	/*!
	   \ingroup rsb_doc_library
	  
	   Write to the info stream (see #RSB_IO_WANT_OUTPUT_STREAM) some memory allocation info.

	   \warning \rsb_warn_soon_to_be_deprecated_msg 
	   \return \rsberrcodemsg
	 */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	RSB_INFO("allocated %zd bytes of memory in %zu chunks\n",rsb__get_g_rsb_memory_count(),rsb__get_g_rsb_allocations_count());
#endif
	return RSB_ERR_NO_ERROR;
}
#define RSB_ALLOC_INFO() RSB_FT_FFL_PRINTF ;rsb_do_meminfo()
#else
#define RSB_ALLOC_INFO() RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
#endif
#define RSB_MTX_FREE_PARANOID(MTXAP) if(MTXAP) { RSB_ALLOC_INFO(); RSB_MTX_FREE(MTXAP); }
#define RSB_FREE_PARANOID(PTR) if(PTR) { RSB_ALLOC_INFO(); RSB_CONDITIONAL_FREE(PTR); }

static rsb_bool_t rsb_random_event(rsb_int out_of) 
{
	rsb_bool_t res = RSB_BOOL_FALSE;

	RSB_ALLOC_INFO();
	if(!out_of)
		goto rf;
	if( rsb__rand_coo_index(out_of) + 1 == out_of)
		res = RSB_BOOL_TRUE;
rf:
	RSB_ALLOC_INFO();
	return res;
}

rsb_err_t rsb_blas_failure_tester(const rsb_char_t*tds)
{
	/**
	 * \ingroup gr_internals
	 * FIXME: still unfinished. this tester shall perform randomized failure-prone operations.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_WANT_FAILURE_TESTER
	struct rsb_limiter ls;
	struct rsb_mtx_t * mtxAp=NULL, * mtxBp=NULL, * mtxCp=NULL;
	rsb_coo_idx_t*IA = NULL;
	rsb_coo_idx_t*JA = NULL;
	rsb_type_t typecodea[] = RSB_MATRIX_TYPE_CODES_ARRAY;
	const int typecodei = 0;
	const int af = 1;/* allow failures */
	rsb_type_t typecode = typecodea[typecodei];
	void*VA = NULL;
	rsb_coo_idx_t nr = 16, nc = nr , nm = RSB_MAX(nr,nc);
	const rsb_coo_idx_t bs = 1, maxdim = 16, mindim = 4;
	const rsb_coo_idx_t lbw = mindim, ubw = mindim;
	rsb_nnz_idx_t nnz = 0;
	//const int fipi=0; /* if -1: probability 0; else, probability 1/fipi */ 
	const rsb_flags_t flags = RSB_FLAG_DEFAULT_MATRIX_FLAGS;
	void * xp = NULL, * yp = NULL;
	rsb_aligned_t alpha[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE], beta[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	const void *alphap = &alpha,* betap = &beta;
	const int /*sp=20,*/hp = 10,lp =100, mp = 1000,ip = 10000;
	const int maxasym = 5;

	RSB_INFO("BASIC FAILURE BASED TEST: BEGIN\n");
	RSB_ALLOC_INFO();
	for(errval = rsb__limiter_init_from_str(&ls,tds);rsb__limiter_continue(&ls);rsb__limiter_step(&ls))
	{
		RSB_ALLOC_INFO();
		typecode = typecodea[typecodei];
	/* FIXME: please write here */
	if(!mtxAp)
	{
		/* random size */
		nr = mindim+rsb__rand_coo_index(maxdim)+lbw+ubw;
		if(rsb_random_event(lp*af))
			{RSB_VARIATION_NOTICE(nc*=1+rsb__rand_coo_index(maxasym));}
		else
			nc = nr;/* square matrix */
		nm = RSB_MAX(nr,nc);
		RSB_INFO("Create a %zd x %zd matrix...\n",(rsb_printf_int_t)nr,(rsb_printf_int_t)nc);

		if(RSB_SOME_ERROR(errval = rsb__generate_blocked_banded_coo(nr,1,lbw,ubw,&IA,&JA,&VA,&nnz,typecode)))
		{
		       	RSB_PERR_GOTO(dca,RSB_ERRM_ES); 
	       	}
		/* try to instantiate one, with a probability of failure injection */
		/* if(...) ... */
		/* change typecode, possibly to an invalid one */
		/* zero a parameter */
		/* dealloc an array and nullify its pointer */
		/* zeros in numbers */
		/* too big indices */
		/* negative indices */
		/* inject a NaN */
		//if(rsb_random_event(lp*af)) {nc=nr = RSB_MAX_MATRIX_DIM+1;}
		//if(rsb_random_event(lp*af)) {nnz=0;}
		if(rsb_random_event(lp*af)) {RSB_FAILURE_NOTICE(typecode = RSB_NUMERICAL_TYPE_INVALID_TYPE);}
		/**/
		/**/
		if(!IA||!JA||!VA)
		{
			RSB_PERR_GOTO(dca,RSB_ERRM_ES);
		}
		RSB_ALLOC_INFO();
		mtxAp = rsb__do_mtx_alloc_from_coo_inplace(VA,IA,JA,nnz,typecode,nr,nc,bs,bs,flags,&errval);
		RSB_ALLOC_INFO();

		if(!mtxAp || RSB_SOME_ERROR(errval))
		{
			/**/
			rsb__do_perror(NULL,errval);
			RSB_PERR_GOTO(dca,RSB_ERRM_ES);
		}
		else
		{
			xp = rsb__calloc_vector(nm,typecode);
			yp = rsb__calloc_vector(nm,typecode);
			if(!xp || !yp)
			{
				RSB_FREE_PARANOID(xp);
				RSB_FREE_PARANOID(yp);
				RSB_PERR_GOTO(dca,RSB_ERRM_ES);
			}
		}
		RSB_ALLOC_INFO();
		goto ndca;
dca:
		RSB_INFO("At:\n");
		rsb__limiter_info(&ls);
		RSB_INFO("Freeing matrix due to error\n");
		RSB_FREE_PARANOID(xp); RSB_FREE_PARANOID(yp);
		RSB_FREE_PARANOID(IA); RSB_FREE_PARANOID(JA); RSB_FREE_PARANOID(VA);
ndca:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
	}
	else
	{
		rsb_trans_t transA = RSB_TRANSPOSITION_N, transB = transA;
		rsb_coo_idx_t incx=1,incy=1;
		RSB_INFO("Use the matrix...\n");
		if(rsb_random_event(ip*af)) {RSB_FAILURE_NOTICE(RSB_MTX_FREE_PARANOID(mtxAp));}
		if(rsb_random_event(lp*af)) {RSB_FAILURE_NOTICE(RSB_FREE_PARANOID(xp));}
		if(rsb_random_event(lp*af)) {RSB_FAILURE_NOTICE(RSB_FREE_PARANOID(yp));}
		if(rsb_random_event(ip*af)) {RSB_FAILURE_NOTICE(transA = RSB_INVALID_TRANS);}
		if(rsb_random_event(mp*af)) {RSB_FAILURE_NOTICE(incx=-1);}
		if(rsb_random_event(mp*af)) {RSB_FAILURE_NOTICE(incy=-1);}
		if(RSB_SOME_ERROR(errval = rsb_do_spmv(transA,alphap,mtxAp,xp,incx,betap,yp,incy)))
		{
			rsb__do_perror(NULL,errval);
			RSB_PERR_GOTO(mdca,RSB_ERRM_ES);
		}
		if(rsb_random_event(2*af))
		{
			RSB_MTX_FREE_PARANOID(mtxBp);
		RSB_ALLOC_INFO();
			/* if( ( mtxBp = rsb__mtx_clone_simple(mtxAp) ) != NULL ) */
			if( rsb__mtx_clone(&mtxBp, mtxAp->typecode,RSB_TRANSPOSITION_N,NULL,mtxAp,flags) == RSB_ERR_NO_ERROR )
			{
		RSB_ALLOC_INFO();
				rsb__do_perror(NULL,errval);
				RSB_PERR_GOTO(mdca,RSB_ERRM_ES);
			}
		}
		if(rsb_random_event(2*af))
		{
			RSB_MTX_FREE_PARANOID(mtxCp);
			mtxCp = rsb__do_matrix_mul(typecode,transA,alphap,mtxAp,transB,betap,mtxBp,&errval);
			if(RSB_SOME_ERROR(errval))
			{
				rsb__do_perror(NULL,errval);
				RSB_PERR_GOTO(mdca,RSB_ERRM_ES);
			}
		}
		if(rsb_random_event(hp*af))
		{
		       	RSB_FAILURE_NOTICE(goto mdca);
		}
		goto nmdca;
mdca:
		rsb__limiter_info(&ls);
		RSB_MTX_FREE_PARANOID(mtxAp);
		RSB_MTX_FREE_PARANOID(mtxBp);
		RSB_FREE_PARANOID(xp);
		RSB_FREE_PARANOID(yp);
		RSB_FREE_PARANOID(IA);
		RSB_FREE_PARANOID(JA);
		RSB_FREE_PARANOID(VA);
nmdca:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
	}
	}
			RSB_MTX_FREE_PARANOID(mtxAp);
			RSB_MTX_FREE_PARANOID(mtxBp);
			RSB_MTX_FREE_PARANOID(mtxCp);
			RSB_FREE_PARANOID(xp);
		       	RSB_FREE_PARANOID(yp);
			RSB_FREE_PARANOID(IA);
		       	RSB_FREE_PARANOID(JA);
		       	RSB_FREE_PARANOID(VA);
	rsb__limiter_info(&ls);
	RSB_ALLOC_INFO();
	errval = rsb__do_check_leak();
	RSB_INFO("BASIC FAILURE BASED TEST: END\n");
#endif /* RSB_WANT_FAILURE_TESTER */
	return errval;
}

/* @endcond */
