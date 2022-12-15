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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains some functions testing accuracy.
 * */

#include "rsb_internals.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

#define RSB_NORMWISE_BACKWARD_ERROR_TOLERANCE 1.e-6
#define RSB_WANT_VERBOSE_ACCURACY_TESTS 0

rsb_err_t rsb__vectors_reinit(void *rhs, void *out, rsb_type_t typecode, rsb_nnz_idx_t rn, rsb_nnz_idx_t on, size_t incr, size_t inco) 
{
	if(RSB_SOME_ERROR(rsb__cblas_Xscal(typecode,on,NULL,out,inco))){ return RSB_ERR_INTERNAL_ERROR;  }
	if(RSB_SOME_ERROR(rsb__fill_with_ones( rhs,typecode,rn,incr))){ return RSB_ERR_INTERNAL_ERROR; }
	return RSB_ERR_NO_ERROR;
}

void * rsb__calloc_vector(rsb_nnz_idx_t n, rsb_type_t typecode)
{
	const size_t so = RSB_SIZEOF(typecode);
	return rsb__calloc(so*n);
}

void * rsb__malloc_vector(rsb_nnz_idx_t n, rsb_type_t typecode)
{
	const size_t so = RSB_SIZEOF(typecode);
	return rsb__malloc(so*n);
}

void * rsb__realloc_vector(void* p, rsb_nnz_idx_t n, rsb_type_t typecode)
{
	const size_t so = RSB_SIZEOF(typecode);
	return rsb__realloc(p,so*n);
}

rsb_err_t rsb__init_rsb_struct_from_coo(struct rsb_mtx_t *mtxAp, const struct rsb_coo_mtx_t *coop)
{
	const rsb_err_t errval = RSB_ERR_NO_ERROR;

	rsb__init_struct(mtxAp);
	mtxAp->bpntr=coop->IA;
	mtxAp->bindx=coop->JA;
	mtxAp->VA=coop->VA;
	mtxAp->typecode=coop->typecode;
	mtxAp->nr=coop->nr;
	mtxAp->nc=coop->nc;
	mtxAp->nnz=coop->nnz;
	mtxAp->br=1;
	mtxAp->bc=1;
	mtxAp->cpntr=0;
	mtxAp->rpntr=0;
	mtxAp->roff=0;
	mtxAp->coff=0;
	mtxAp->broff=0;
	mtxAp->bcoff=0;

	return errval;
}

rsb_err_t rsb__do_spmv_fullword_coo(const struct rsb_coo_mtx_t*coop, rsb_flags_t flags, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA)
{
	/* FIXME: UNFINISHED  */
	struct rsb_mtx_t mtxA;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb__init_rsb_struct_from_coo(&mtxA,coop);
	if(RSB_SOME_ERROR(errval)) goto err;
	flags = RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS|RSB_DO_FLAG_FILTEROUT((flags),RSB_DO_FLAGS_EXTRACT_STORAGE(flags));
	if((errval = rsb__do_set_init_storage_flags(&mtxA,flags))!=RSB_ERR_NO_ERROR)
		goto err;
	RSB_DO_FLAG_ADD(mtxA.flags,RSB_DO_FLAG_FILTERONLY(flags,RSB_FLAGS_RSB_AGNOSTIC));
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_non_recursive(&mtxA,x,y,alphap,betap,incx,incy,transA RSB_DEFAULT_INNER_NRHS_SPMV_ARGS	) );
err:
	return errval;
}

static rsb_err_t rsb_do_check_normwise_backward_error(const struct rsb_mtx_t*mtxAp, const void *X, const void *AX, const void * B, rsb_trans_t transA)
{
	/* normwise backward error in the infinity norm */
	/* Note: this won't work for integer, obviously.  */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_aligned_t Xinorm[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t Binorm[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t Ainorm[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t denominator[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t err[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t eps[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_type_t typecode = mtxAp->typecode;
	rsb_coo_idx_t tm = RSB_MTX_TRANSPOSED_ROWS(mtxAp,transA);
	rsb_coo_idx_t tk = RSB_MTX_TRANSPOSED_COLS(mtxAp,transA);

	RSB_NUMERICAL_TYPE_SET_ELEMENT_FROM_DOUBLE(eps,RSB_NORMWISE_BACKWARD_ERROR_TOLERANCE,typecode);
	if(RSB_SOME_ERROR(errval = rsb__do_matrix_norm(mtxAp,Ainorm,RSB_EXTF_NORM_INF))){RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(RSB_SOME_ERROR(errval = rsb__vector_sum_of_abs(Xinorm,X,typecode,tk))){RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(RSB_SOME_ERROR(errval = rsb__vector_sum_of_abs(Binorm,B,typecode,tm))){RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(RSB_SOME_ERROR(errval = rsb__vector_mult(Ainorm,Xinorm,denominator,typecode,1))){RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(RSB_SOME_ERROR(errval = rsb__util_vector_add(denominator,Binorm,typecode,1))){RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(RSB_IS_ELEMENT_ZERO(denominator,typecode))
	{
		errval = RSB__ERR_ZERO_INF_NORM;
#if RSB_OUT_ERR_VERBOSITY>=3
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#else /* RSB_OUT_ERR_VERBOSITY */
		goto err;
#endif /* RSB_OUT_ERR_VERBOSITY */
	}
	if(RSB_SOME_ERROR(errval = rsb__vector_sum_of_abs_diffs(err,AX,B,typecode,tk))){RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(RSB_SOME_ERROR(errval = rsb__util_vector_div(err,denominator,typecode,1))){RSB_PERR_GOTO(err,RSB_ERRM_ES);}

	if(!RSB_IS_ELEMENT_LESS_THAN(err,eps,typecode))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		if(RSB_WANT_VERBOSE_ACCURACY_TESTS)
			RSB_ERROR("error is %lg, more than %lg!\n",*(double*)(&err[0]),*(double*)(&eps[0])); /* consider using rsb__debug_print_value */
		goto err;
	}
	else
	{
		if(RSB_WANT_VERBOSE_ACCURACY_TESTS)
			RSB_STDOUT("error is %lg, less than %lg.\n",*(double*)(&err[0]),*(double*)(&eps[0]));
	}
err:
	return errval;
}

rsb_err_t rsb__do_spmv_accuracy_test(const struct rsb_coo_mtx_t*coop, rsb_thread_t * ca, rsb_thread_t cn, rsb_flags_t flags)
{
	/* FIXME: UNFINISHED  */
	/* this is mostly a debug function */
	/* FIXME: what about incx/incy/trans ? */
	/* FIXME: shal support ca=NULL, cn=0  */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_thread_t ci=0;
	void*X=NULL,*Y=NULL,*Z=NULL;
	rsb_coo_idx_t vd;
	rsb_coo_idx_t incx=1,incy=1;
	rsb_trans_t transA = RSB_DEFAULT_TRANSPOSITION;
	rsb_aligned_t alpha[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t twoalpha[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t zsum[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_thread_t default_ca[2],default_cn=2;
	void*alphap=alpha,*betap=NULL/*,*twoalphap=twoalpha*/;
	if(RSB_SOME_ERROR(errval = rsb__util_is_valid_coo_struct(coop)))
	{
		goto err;
	}
	rsb__util_set_area_to_converted_integer(alpha,coop->typecode,1);
	rsb__util_set_area_to_converted_integer(twoalpha,coop->typecode,2);
	vd = RSB_MAX(coop->nr,coop->nc);
	X = rsb__malloc_vector(vd,coop->typecode);
	Y = rsb__malloc_vector(vd,coop->typecode);
	Z = rsb__calloc_vector(vd,coop->typecode);
	if(!X||!Y||!Z){RSB_ERROR("vector allocation problems!\n");RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_ENOMEM);goto err;}
	rsb__fill_with_ones(X,coop->typecode,vd,1);
	/* we compute spmv on this coo instance */
	errval = rsb__do_spmv_fullword_coo(coop,flags,X,Z,alphap,betap,incx,incy,transA);
	if(RSB_SOME_ERROR(errval)){RSB_ERROR("!\n");goto err;}
	errval = rsb__util_vector_sum(zsum,Z,coop->typecode,vd);
	if(RSB_SOME_ERROR(errval)){RSB_ERROR("!\n");goto err;}
	rsb__vector_to_abs(zsum,coop->typecode,1);

	/* FIXME: need clone + cleanup here! */
	if(cn==0 || ca==NULL)
	{
		ca = default_ca;
		cn = default_cn;
		ca[0] = 1;
		ca[1] = rsb_get_num_threads();
		if(ca[1] <= ca[0]) cn = 1;
	}
	for(ci=0;ci<cn;++ci)
	{
		struct rsb_mtx_t*mtxAp=NULL;
		rsb_aligned_t ysum[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
		rsb__set_num_threads(ca[ci]);
		if(RSB_WANT_VERBOSE_ACCURACY_TESTS)
		RSB_STDOUT("for %d threads:\n",ca[ci]);
		mtxAp = rsb__do_mtx_alloc_from_coo_const(coop->VA,coop->IA,coop->JA,coop->nnz,coop->typecode,coop->nr,coop->nc,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags,&errval);
		if(RSB_SOME_ERROR(errval)){RSB_ERROR("!\n");goto err;}
		if(!mtxAp){errval = RSB_ERR_INTERNAL_ERROR;RSB_ERROR("!\n");goto err;}
		rsb__cblas_Xscal(coop->typecode,vd,NULL,Y,1);
		/* TODO: may iterate the following in order to test indeterminism */
		if(mtxAp->nnz != coop->nnz)
		{
			/* input not cleaned up? */
			if(RSB_WANT_VERBOSE_ACCURACY_TESTS)
				RSB_STDOUT("%ld vs %ld nnz ? input not cleaned up?\n",(long int)mtxAp->nnz,(long int)coop->nnz);
			errval = RSB_ERR_BADARGS;
			goto ierr;
		}
		errval = rsb_do_spmv(transA,alphap,mtxAp,X,incx,betap,Y,incy);
		if(RSB_SOME_ERROR(errval)){RSB_ERROR("!\n");goto ierr;}
		errval = rsb__util_vector_sum(ysum,Y,coop->typecode,vd);
		if(RSB_SOME_ERROR(errval)){RSB_ERROR("!\n");goto ierr;}
		rsb__vector_to_abs(ysum,coop->typecode,1);
		//rsb__debug_print_vector(ysum,1,coop->typecode,1);
		if(!RSB_SOME_ERROR(rsb__do_are_same(ysum,zsum,1,coop->typecode,1,1)))
		{
			/* same result (no numerical error at all). no further check really needed. */
			if(RSB_WANT_VERBOSE_ACCURACY_TESTS)
				RSB_STDOUT("Bitwise identical values.\n");
			errval = rsb_do_check_normwise_backward_error(mtxAp,X,Y,Z,transA);
			if(errval==RSB__ERR_ZERO_INF_NORM)
				errval = RSB_ERR_NO_ERROR; /* ok, happens */
			if(RSB_SOME_ERROR(errval))
				RSB_PERR_GOTO(ierr,RSB_ERRM_ES);
		}
		else
		{
			if(RSB_WANT_VERBOSE_ACCURACY_TESTS)
				RSB_STDOUT("Non-identical values. Checking backward error.\n");
			errval = rsb_do_check_normwise_backward_error(mtxAp,X,Y,Z,transA);
			if(errval==RSB__ERR_ZERO_INF_NORM)
				RSB_PERR_GOTO(ierr,"Error: Null infinite norm and different spmv results !?\n");
			if(RSB_SOME_ERROR(errval)){RSB_PERR_GOTO(ierr,RSB_ERRM_ES);}
		}
ierr:
		RSB_MTX_FREE(mtxAp);
	}
err:
	RSB_CONDITIONAL_FREE(X);
	RSB_CONDITIONAL_FREE(Y);
	RSB_CONDITIONAL_FREE(Z);
	return errval;
}

/* @endcond */

