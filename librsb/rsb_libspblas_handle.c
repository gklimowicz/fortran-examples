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
/* @cond INNERDOC */
/**
 * \internal
 * @file
 * @author Michele Martone
 * @brief
 * 
 * Sparse BLAS extras (and internals).
 * \internal
 * 
 * */
/*
 * TODO: support for blas_field_type, blas_base_type, blas_sort_type, ...
 * 	 error reporting, input sanitizing, error handling, ...
 * */
/*
*/
/*   #include "blas_sparse/blas_sparse.h"*/
//#include "blas_sparse/blas_enum.h"
#include "rsb.h"
#include "rsb_libspblas.h"
#include "rsb_internals.h"
/* #include "libspblas_handle.h" */
#include "rsb_psblas.h"
#include "rsb_do.h"
/*  #include "blas_sparse/blas_sparse_proto.h"*/

RSB_INTERNALS_COMMON_HEAD_DECLS

#define RSB_WANT_SPARSE_BLAS_EXTENSIONS 1

/* #define RSB_BLAS_INVALID_VAL blas_invalid_handle */

#define RSB_ATPNAME_ANY 0
#define RSB_SPBLAS_DEF_TUNING_ROUNDS RSB_CONST_MAX_TUNING_ROUNDS
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING 
#define RSB_SPB_AT_OP(MTXAP,RNT,HINT,NRHS,ORDER,ALPHAP,BETAP,LHS,RHS,LDC,LDB,OPTYPE) 	\
	if((MTXAP) && (HINT) == RSB_SPB_THR_STR_AUTO_NEXTOP /* ... next operation */ ) \
	{		\
		/* errval = */ rsb__tune_spxx(&(MTXAP), NULL, &(RNT), RSB_SPBLAS_DEF_TUNING_ROUNDS, RSB_CONST_DEF_MS_AT_AUTO_STEPS, RSB_CONST_DEF_MS_AT_AUTO_STEPS,RSB_CONST_AT_OP_SAMPLES_MIN, RSB_CONST_AT_OP_SAMPLES_MAX, 0, trans, ALPHAP, NULL, NRHS, ORDER, NULL, LDB, BETAP, NULL, LDC, OPTYPE, NULL, NULL, NULL, RSB_AUT0_TUNING_SILENT, NULL, NULL, NULL, NULL, NULL); \
		if(RHS == NULL || LHS == NULL) { brv = RSB_BLAS_NO_ERROR; goto err; /* wanted just tuning */ } \
		(HINT) = RSB_SPB_THREADS_DEFAULT; \
	}
#endif /* RSB_BLAS_WANT_EXPERIMENTAL_TUNING */

#define RSB_BLAS_IS_ATPNAME_OFF(PNAME) ( (PNAME) == blas_rsb_spmv_n_autotuning_off || (PNAME) == blas_rsb_spmv_t_autotuning_off || (PNAME) == blas_rsb_spmv_autotuning_off )
#define RSB_BLAS_IS_ATPNAME_ON(PNAME)  ( (PNAME) == blas_rsb_spmv_n_autotuning_on  || (PNAME) == blas_rsb_spmv_t_autotuning_on  || (PNAME) == blas_rsb_spmv_autotuning_on || (PNAME) == blas_rsb_autotune_next_operation )
#define RSB_BLAS_IS_ATPNAME_ANY(PNAME)  ( RSB_ATPNAME_ANY == RSB_ATPNAME_ANY )
#define RSB_BLAS_IS_ATPNAME(PNAME) ( RSB_BLAS_IS_ATPNAME_OFF(PNAME) || RSB_BLAS_IS_ATPNAME_ON(PNAME) )
#define RSB_BLAS_ALLOW_MTX_UPD 1 /* this allows updates in either blas_rsb_duplicates_sum or blas_rsb_duplicates_ovw style; after usds(), blas_rsb_duplicates_ovw is restored and if desired, blas_rsb_duplicates_sum should be set again. */
#define RSB_TUNING_NEW_STYLE 1 /* FIXME: not yet active */

static struct 
{
	struct rsb_blas_sparse_matrix_t * bsms;
	size_t n;
	size_t next_handle;
} rsb_blas_handles;

               /* Service Routines */

static int rsb_compar_vbr_blas_sparse_matrix_t(const void * ap, const void * bp)
{
	/**
	 \ingroup gr_internals
	 */
	const blas_sparse_matrix ha = ((struct rsb_blas_sparse_matrix_t*)ap)->handle;
	const blas_sparse_matrix hb = ((struct rsb_blas_sparse_matrix_t*)bp)->handle;

        return
                 ( ha >  hb ) ? 1 :
                 (( ha == hb ) ? 0 : -1);
}

static struct rsb_blas_sparse_matrix_t * rsb__BLAS_matrix_retrieve(blas_sparse_matrix handle)
{
	/**
	 \ingroup gr_internals
	 */
	/*
	 shall retrieve the internal data structure associated to the handle
	 */
	struct rsb_blas_sparse_matrix_t key;

	if ( handle == blas_invalid_handle )
		return NULL;

	key.handle = handle;
	return bsearch(&key,rsb_blas_handles.bsms,rsb_blas_handles.n,sizeof(struct rsb_blas_sparse_matrix_t),rsb_compar_vbr_blas_sparse_matrix_t);
}

static size_t rsb__BLAS_matrix_retrieve_index(blas_sparse_matrix handle)
{
	/**
	 \ingroup gr_internals
	 */
	/*
	 shall retrieve the internal data structure associated to the handle
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;

	bsm = rsb__BLAS_matrix_retrieve(handle);
	if(!bsm)
		return RSB_BLAS_HANDLE_INVALID;
	return bsm-rsb_blas_handles.bsms;
}

rsb_err_t rsb__BLAS_is_type_supported(rsb_char_t c)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_err_t errval = RSB_ERR_UNSUPPORTED_TYPE;

	switch(c)
	{
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
		case('c'): case('C'):
#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
		case('d'): case('D'):
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
		case('z'): case('Z'):
#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
#ifdef RSB_NUMERICAL_TYPE_FLOAT
		case('s'): case('S'):
#endif /* RSB_NUMERICAL_TYPE_FLOAT */
		errval = RSB_ERR_NO_ERROR;
		break;
#ifdef RSB_NUMERICAL_TYPE_INT
		case('i'): case('I'): /* NOTE: BLAS-unsupported, but for e.g. rsb_mtx_alloc_from_coo_begin is OK */
#endif /* RSB_NUMERICAL_TYPE_INT */
		default:
		break;
	};
	return errval;
}

struct rsb_mtx_t * rsb__BLAS_inner_matrix_retrieve(blas_sparse_matrix handle)
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_mtx_t * mtxAp = NULL;

	if(rsb_blas_handles.n==1)
	{
		/* this is a trick */
	       	mtxAp = rsb_blas_handles.bsms[0].mtxAp;
	}
	else
	{
		struct rsb_blas_sparse_matrix_t * bsk = NULL;
	       	bsk = rsb__BLAS_matrix_retrieve(handle);
		if(bsk)
			mtxAp = bsk->mtxAp;
	}
	return mtxAp;
}

rsb_err_t rsb__BLAS_handles_free(void)
{
	/**
	 \ingroup gr_internals
	 This function shall be called as a finalizer.
	 It can be dangerous, if called before initialization or before other sparse BLAS-related operations.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(rsb_blas_handles.bsms==NULL)
		;
	else
	{
		size_t i;
		const size_t n = rsb_blas_handles.n;

		for(i=0;i<n;++i)
		{
			const blas_sparse_matrix handle = rsb_blas_handles.bsms[n-i-1].handle;
			rsb__BLAS_Xusds( handle );
		}
		RSB_CONDITIONAL_FREE(rsb_blas_handles.bsms);
	}
	RSB_DO_ERR_RETURN(errval)
}

static blas_sparse_matrix rsb__BLAS_handle_alloc(void)
{
	/**
	 \ingroup gr_internals
	 */
	/*
	 shall allocate a new matrix handle
	 NOTE: this allocates an initial pool of descriptors which shall be freed at a point.
	 */
	struct rsb_blas_sparse_matrix_t * nbsms = NULL;
	size_t n;
	blas_sparse_matrix handle = RSB_BLAS_INVALID_VAL;

	RSB_DEBUG_ASSERT(rsb_blas_handles.n <= RSB_BLAS_MATRICES_MAX);
	RSB_DEBUG_ASSERT(rsb_blas_handles.n>=0);

	if(rsb_blas_handles.n >= RSB_BLAS_MATRICES_MAX)
		RSB_PERR_GOTO(err,"matrix limit reached")

	if(rsb_blas_handles.bsms==NULL)
		nbsms = rsb__calloc(sizeof(struct rsb_blas_sparse_matrix_t)*(rsb_blas_handles.n+1)),
		rsb_blas_handles.next_handle = RSB_BLAS_FIRST_HANDLE;/* handles reset */
	else
		nbsms = rsb__realloc(rsb_blas_handles.bsms,sizeof(struct rsb_blas_sparse_matrix_t)*(rsb_blas_handles.n+1));

	if(!nbsms)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	else
		rsb_blas_handles.bsms=nbsms;
	/* we have one extra struct in, now. note that we should blank it before use  */
	handle = rsb_blas_handles.next_handle<RSB_BLAS_FIRST_HANDLE?RSB_BLAS_FIRST_HANDLE:rsb_blas_handles.next_handle;
	for(n=0;n<rsb_blas_handles.n && rsb_blas_handles.bsms[n].handle<handle;++n)
		;/* TODO: this is inefficient. should fix this */

	/*  inserting at n will keep the vector sorted ascendent */
	if(n != rsb_blas_handles.n)	
	{
		RSB_MEMMOVE(rsb_blas_handles.bsms+n+1,rsb_blas_handles.bsms+n,sizeof(struct rsb_blas_sparse_matrix_t)*(rsb_blas_handles.n-n));
	}

	rsb_blas_handles.next_handle++;
	if(rsb_blas_handles.next_handle>RSB_BLAS_LAST_HANDLE)
		rsb_blas_handles.next_handle = RSB_BLAS_FIRST_HANDLE;

	RSB_BZERO(rsb_blas_handles.bsms+n,sizeof(rsb_blas_handles.bsms[n]));/* we blank the new handle */
	rsb_blas_handles.bsms[n].handle=handle;
	rsb_blas_handles.n++;
	return handle;
err:
	return RSB_BLAS_INVALID_VAL;
}

blas_sparse_matrix rsb__BLAS_handle_free(blas_sparse_matrix handle)
{
	/**
	 \ingroup gr_internals
	 */
	/*
	 shall free a new handle
	 */
	size_t n;
	struct rsb_blas_sparse_matrix_t * nbsms = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if(rsb_blas_handles.n<1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	n = rsb__BLAS_matrix_retrieve_index(handle);
	if(n>RSB_BLAS_MATRICES_MAX)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if(n<rsb_blas_handles.n-1)
		RSB_MEMMOVE(rsb_blas_handles.bsms+n,rsb_blas_handles.bsms+n+1,sizeof(struct rsb_blas_sparse_matrix_t)*(rsb_blas_handles.n-(n+1)));

	rsb_blas_handles.n--;

	nbsms = rsb__realloc(rsb_blas_handles.bsms,sizeof(struct rsb_blas_sparse_matrix_t)*(rsb_blas_handles.n));
	if(rsb_blas_handles.n==0 || (nbsms!=NULL))
		rsb_blas_handles.bsms=nbsms;

	retval = handle;
err:
	return retval;
}

               /* Matrix assembly routines */

blas_sparse_matrix rsb__BLAS_new_matrix_begin(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnzest, rsb_type_t typecode, rsb_coo_idx_t br, rsb_coo_idx_t bc, const rsb_coo_idx_t*rbp, const rsb_coo_idx_t*cbp)
{
	/**
	 \ingroup gr_internals
	 */
	/* 
	 * shall allocate a new matrix descriptor and handle
	 * */
	blas_sparse_matrix handle = RSB_BLAS_INVALID_VAL;
	struct rsb_blas_sparse_matrix_t * bsm = NULL;

	if(rbp && cbp)
	{
		m=k=0;
		RSB_FCOO_ASUM(m,rbp,0,br);
		RSB_FCOO_ASUM(k,cbp,0,bc);
		nnzest=1+RSB_MAX(m,k);
	}

	if(!RSB_ARE_VALID_MATRIX_INIT_PARS(m,k,nnzest?nnzest:1,typecode))
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	RSB_DEBUG_ASSERT(nnzest>0);

	if(br<0 || bc<0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)


	if( (rbp && !cbp) || ( cbp && !rbp))
		RSB_PERR_GOTO(err,RSB_ERRM_ES)


	if( (handle = rsb__BLAS_handle_alloc()) == RSB_BLAS_INVALID_VAL)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	if( (bsm = rsb__BLAS_matrix_retrieve(handle) ) == NULL )
		RSB_PERR_GOTO(errh,RSB_ERRM_ES)

	RSB_BZERO_P(bsm); // should be already blanked, unless cycling in the descriptors array

	if( rbp && cbp && br && bc )
	{
		bsm->rbp = rsb__clone_area_with_extra(rbp,sizeof(rsb_coo_idx_t)*(br),sizeof(rsb_coo_idx_t),0);
		bsm->cbp = rsb__clone_area_with_extra(cbp,sizeof(rsb_coo_idx_t)*(bc),sizeof(rsb_coo_idx_t),0);
		if( (!bsm->rbp) || (!bsm->cbp) )
			RSB_PERR_GOTO(errr,RSB_ERRM_ES)
		bsm->rbp[0]=bsm->cbp[0]=0;
		rsb__do_prefix_sum_coo_idx_t(bsm->rbp,br+1);
		rsb__do_prefix_sum_coo_idx_t(bsm->cbp,bc+1);
	}

	if(br==0)
		br = 1;

	if(bc==0)
		bc = 1;

	bsm->symmetry=blas_general;
	bsm->diag_type=blas_non_unit_diag;
	bsm->mtxAp=NULL;
	bsm->coomatrix.nnz=nnzest;
	bsm->coomatrix.nr=m;
	bsm->coomatrix.nc=k;
	bsm->coomatrix.typecode=typecode;
	bsm->nnzin=0;
	bsm->k=br;
	bsm->l=bc;
	bsm->handle=handle;
	bsm->type=blas_new_handle;
	bsm->dupstra = blas_rsb_duplicates_sum;
	bsm->fmt_hint = blas_rsb_rep_rsb;
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING 
	bsm->opt_mvn_hint = bsm->opt_mvt_hint = RSB_SPB_THREADS_DEFAULT;
#endif /* RSB_BLAS_WANT_EXPERIMENTAL_TUNING */
	switch(typecode)
	{
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
		case('d'): case('D'):
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
		case('z'): case('Z'):
#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
		bsm->fprecision = blas_double_precision;break;
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
		case('c'): case('C'):
#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
#ifdef RSB_NUMERICAL_TYPE_FLOAT
		case('s'): case('S'):
#endif /* RSB_NUMERICAL_TYPE_FLOAT */
		bsm->fprecision = blas_single_precision;break;
#ifdef RSB_NUMERICAL_TYPE_INT
		case('i'): case('I'): /* for e.g. rsb_mtx_alloc_from_coo_begin */
#endif /* RSB_NUMERICAL_TYPE_INT */
		bsm->fprecision = 'I';break; /* FIXME; may foresee 'blas_int_precision' extension */
		default:
		RSB_PERR_GOTO(errr,RSB_ERRM_ES)
	       	break;
	};
	/* allocate temporary resources */
	if(rsb__allocate_coo_matrix_t(&bsm->coomatrix)==NULL)
		RSB_PERR_GOTO(errr,RSB_ERRM_ES)

	bsm->type=blas_open_handle;

	return handle;
errr:
	RSB_CONDITIONAL_FREE(bsm->rbp);
	RSB_CONDITIONAL_FREE(bsm->cbp);
errh:
	rsb__BLAS_handle_free(handle);
err:
	return RSB_BLAS_INVALID_VAL;
}

static blas_sparse_matrix rsb__BLAS_new_matrix_expand_store(struct rsb_blas_sparse_matrix_t * bsm, rsb_nnz_idx_t to_nnz)
{
	/**
	 \ingroup gr_internals
	 */
	/* expand the temporary store for inserting data */
	RSB_DEBUG_ASSERT(bsm);

	if( bsm->type!=blas_open_handle )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	RSB_DEBUG_ASSERT(bsm->coomatrix.IA);
	RSB_DEBUG_ASSERT(bsm->coomatrix.JA);
	RSB_DEBUG_ASSERT(bsm->coomatrix.VA);

	if( bsm->nnzin == 0 && to_nnz == 0 && bsm->coomatrix.nnz != 0 )
	{
		/* when more nonzoeroes were estimated to come than they did */
		rsb__destroy_coo_matrix_t(&bsm->coomatrix);
		bsm->coomatrix.nnz = 0;
	}
	else
		if( /*  ( to_nnz != bsm->coomatrix.nnz ) && */rsb__reallocate_coo_matrix_t(&bsm->coomatrix, to_nnz)==NULL && to_nnz>0)
		{
			RSB_PERR_GOTO(err,"Failed reallocation from %zd to %zd nonzeroes\n", (size_t) bsm->coomatrix.nnz, (size_t) to_nnz)
		}

	return bsm->handle;
err:
	return RSB_BLAS_INVALID_VAL;
}

static blas_sparse_matrix rsb__BLAS_new_matrix_expand_store_try(struct rsb_blas_sparse_matrix_t * bsm, rsb_nnz_idx_t to_nnz_min, rsb_nnz_idx_t to_nnz_max)
{
	/* try opportunistically an upper limit first; if no success, try the lower one */
	/* TODO: on the long run, shall substitute rsb__BLAS_new_matrix_expand_store */
	blas_sparse_matrix handle = RSB_BLAS_INVALID_VAL;

	handle = rsb__BLAS_new_matrix_expand_store(bsm, to_nnz_max);
	if( handle == RSB_BLAS_INVALID_VAL && to_nnz_min < to_nnz_max )
		handle = rsb__BLAS_new_matrix_expand_store(bsm, to_nnz_min);
	return handle;
}

               /* Nonzeroes insertion routines */

static blas_sparse_matrix rsb__BLAS_new_matrix_insert_block(struct rsb_blas_sparse_matrix_t * bsm, const void * val, rsb_blas_int_t row_stride, rsb_blas_int_t col_stride, rsb_blas_int_t i, rsb_blas_int_t j)
{
	/**
	 \ingroup gr_internals
	 No check is performed on the block size arrays.
	 */
	rsb_coo_idx_t ii = 0, jj = 0, ob = 0;
	int rb = 0, cb = 0, roff = 0, coff = 0;
	size_t nnz = 0,es = 0;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	RSB_DEBUG_ASSERT(bsm);
        es = RSB_NUMERICAL_TYPE_SIZE(bsm->coomatrix.typecode);

	ob = bsm->off;
	i -= ob;
	j -= ob;
	if(bsm->rbp)
	{
		rb=bsm->rbp[i+1]-bsm->rbp[i];
		roff=bsm->rbp[i];
	}
	else
	{
		rb=bsm->k;
		roff=rb*i;
	}

	if(bsm->cbp)
	{
		cb=bsm->cbp[j+1]-bsm->cbp[j];
		coff=bsm->cbp[j];
	}
	else
	{
		cb=bsm->l;
		coff=cb*j;
	}

	nnz=rb*cb;

	if(nnz==0)
	{
		/* in different cases; e.g.: not a blocked mtxAp, zero sized block, ... */
		RSB_ERROR(RSB_ERRM_BNCS);
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(bsm->nnzin+nnz > bsm->coomatrix.nnz)
		if(rsb__BLAS_new_matrix_expand_store(bsm,RSB_MAX(bsm->nnzin+nnz,2*bsm->coomatrix.nnz))==RSB_BLAS_INVALID_VAL)
			RSB_PERR_GOTO(err,RSB_ERRM_ES)

	/* please note that structural zeroes are not supported by the sparse blas interface,
	 * so we don't use 
	 *  if(!RSB_DO_FLAG_HAS(bsm->flags, RSB_FLAG_DISCARD_ZEROS) ... )
	 */
	for (ii=0; ii<rb; ii++)
	{
		for (jj=0; jj<cb; jj++)
		{
			const rsb_nnz_idx_t nzoff=bsm->nnzin+ii*cb+jj;
			const rsb_byte_t*eval=((const rsb_byte_t*)val)+es*(ii*row_stride+jj*col_stride);
			bsm->coomatrix.IA[nzoff]=roff+ii;
			bsm->coomatrix.JA[nzoff]=coff+jj;
			if(!RSB_IS_ELEMENT_ZERO(eval,bsm->coomatrix.typecode))
 		 		rsb__memcpy(((rsb_byte_t*)bsm->coomatrix.VA)+es*nzoff,eval,es);
		}
	}
	bsm->nnzin += nnz;

	retval = RSB_BLAS_NO_ERROR;
err:
	return retval;
}

static blas_sparse_matrix rsb__BLAS_new_matrix_insert_row(struct rsb_blas_sparse_matrix_t * bsm, rsb_blas_int_t i, rsb_blas_int_t nnz, const void * val, const rsb_blas_int_t *jndx )
{
	/**
	 \ingroup gr_internals
	 */
	rsb_blas_int_t k;
	int ob = 0;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	RSB_DEBUG_ASSERT(bsm);

	ob = bsm->off;

	if(bsm->nnzin+nnz > bsm->coomatrix.nnz)
		if(rsb__BLAS_new_matrix_expand_store(bsm,RSB_MAX(bsm->nnzin+nnz,2*bsm->coomatrix.nnz))==RSB_BLAS_INVALID_VAL)
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
  	rsb__numerical_memcpy(bsm->coomatrix.typecode,bsm->coomatrix.VA,bsm->nnzin,val,0,nnz);
	for(k=0;k<nnz;++k)
		bsm->coomatrix.IA[bsm->nnzin+k] = i-ob,
		bsm->coomatrix.JA[bsm->nnzin+k] = jndx[k]-ob;
	bsm->nnzin += nnz;

	retval = RSB_BLAS_NO_ERROR;
err:
	return retval;
}

static blas_sparse_matrix rsb__BLAS_new_matrix_insert_col(struct rsb_blas_sparse_matrix_t * bsm, rsb_blas_int_t j, rsb_blas_int_t nnz, const void * val, const rsb_blas_int_t *indx )
{
	/**
	 \ingroup gr_internals
	 */
	/* append coo data */
	/* ok */
	int ob = 0;
	rsb_blas_int_t k;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	RSB_DEBUG_ASSERT(bsm);

	if( bsm->nnzin+nnz > bsm->coomatrix.nnz )
		if( rsb__BLAS_new_matrix_expand_store(bsm,RSB_MAX(bsm->nnzin+nnz,2*bsm->coomatrix.nnz)) == RSB_BLAS_INVALID_VAL )
			RSB_PERR_GOTO(err,RSB_ERRM_ES)

	ob=bsm->off;

  	rsb__numerical_memcpy(bsm->coomatrix.typecode,bsm->coomatrix.VA,bsm->nnzin,val,0,nnz);
	for(k=0;k<nnz;++k)
		bsm->coomatrix.IA[bsm->nnzin+k]=indx[k]-ob,
		bsm->coomatrix.JA[bsm->nnzin+k]=j-ob;
	bsm->nnzin += nnz;

	retval = RSB_BLAS_NO_ERROR;
err:
	return retval;
}

static blas_sparse_matrix rsb__BLAS_new_matrix_insert_clique(struct rsb_blas_sparse_matrix_t * bsm, const rsb_blas_int_t k, const rsb_blas_int_t l, const void * val, const rsb_blas_int_t row_stride, const rsb_blas_int_t col_stride, const rsb_blas_int_t *indx, const rsb_blas_int_t *jndx )
{
	/**
	 \ingroup gr_internals
	 append coo data
	*/
	const size_t nnz = k*l;
	int i = 0, j = 0, ob = 0;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	RSB_DEBUG_ASSERT(bsm);

#if 0
	RSB_ERROR(RSB_ERRM_CIIAUF);
	RSB_PERR_GOTO(err,RSB_ERRM_ES)
	/* FIXME: unfinished! */
#endif

	if( bsm->nnzin+nnz > bsm->coomatrix.nnz )
		if( rsb__BLAS_new_matrix_expand_store(bsm,RSB_MAX(bsm->nnzin+nnz,2*bsm->coomatrix.nnz)) == RSB_BLAS_INVALID_VAL )
			RSB_PERR_GOTO(err,RSB_ERRM_ES)

	ob=bsm->off;

	for (i=0; i<k; i++)
	{
		for (j=0; j<l; j++)
		{
			bsm->coomatrix.IA[bsm->nnzin+i*l+j] = indx[i]-ob,
			bsm->coomatrix.JA[bsm->nnzin+i*l+j] = jndx[j]-ob;
 		 	rsb__numerical_memcpy(bsm->coomatrix.typecode,bsm->coomatrix.VA,(bsm->nnzin+(i*l+j)),val,(i*row_stride+j*col_stride),1); /* FIXME: this is just one element */
		}
	}
	bsm->nnzin += nnz;
	retval = RSB_BLAS_NO_ERROR;
err:
	return retval;
}

static rsb_blas_int_t rsb__BLAS_new_matrix_insert_entries( struct rsb_blas_sparse_matrix_t * bsm, const rsb_blas_int_t nnz, const char *val, const rsb_blas_int_t * indx, const rsb_blas_int_t * jndx)
{
	/**
	 \ingroup gr_internals
	 */
	int ob = 0;
	rsb_blas_int_t k;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	RSB_DEBUG_ASSERT(bsm);

#if RSB_BLAS_ALLOW_MTX_UPD
	if( bsm->type == blas_valid_handle )
	{
		rsb_err_t errval = RSB_ERR_NO_ERROR;

		errval = rsb__do_set_elements(bsm->mtxAp, val, indx, jndx, nnz,  
			bsm->base == blas_one_base ? RSB_FLAG_FORTRAN_INDICES_INTERFACE : RSB_FLAG_NOFLAGS);
		if( RSB_SOME_ERROR( errval ))
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		goto ok;
	}
#endif /* RSB_BLAS_ALLOW_MTX_UPD */

	if(bsm->nnzin+nnz > bsm->coomatrix.nnz)
		if( rsb__BLAS_new_matrix_expand_store_try(bsm, (bsm->nnzin+nnz), RSB_MAX(bsm->nnzin+nnz,2*bsm->coomatrix.nnz)) == RSB_BLAS_INVALID_VAL )
		{
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}

	ob = bsm->off;
  	rsb__numerical_memcpy(bsm->coomatrix.typecode,bsm->coomatrix.VA,bsm->nnzin,val,0,nnz);
	for(k=0;k<nnz;++k)
		bsm->coomatrix.IA[bsm->nnzin+k]=indx[k]-ob,
		bsm->coomatrix.JA[bsm->nnzin+k]=jndx[k]-ob;

	bsm->nnzin += nnz;
ok:
	retval = RSB_BLAS_NO_ERROR;
err:
	return retval;
}

static rsb_flags_t rsb__BLAS_new_matrix_finish_flags(struct rsb_blas_sparse_matrix_t * bsm)
{
	/* sets rsb flags based on user-set properties */
	rsb_flags_t flags = RSB_FLAG_NOFLAGS | RSB_FLAG_SORT_INPUT | RSB_FLAG_OWN_PARTITIONING_ARRAYS /* | RSB_FLAG_WANT_BCSS_STORAGE */ ;

	switch(bsm->symmetry)
	{
		case(blas_lower_symmetric):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER_SYMMETRIC);
		break;
		case(blas_lower_hermitian):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER_HERMITIAN);
		break;
		case(blas_upper_symmetric):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER_SYMMETRIC);
		break;
		case(blas_upper_hermitian):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER_HERMITIAN);
		break;
		case(blas_lower_triangular):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER_TRIANGULAR);
		break;
		case(blas_upper_triangular):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER_TRIANGULAR);
		break;
		case(blas_triangular):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
		break;
		case(blas_general):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_NOFLAGS);
		break;
#if 0
		case(blas_symmetric):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_SYMMETRIC);
		break;
		case(blas_hermitian):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_HERMITIAN);
		break;
		case(blas_triangular):
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
		break;
#else
		default:
		return RSB_INVALID_FLAGS;
#endif
	}

	if( bsm->dupstra == blas_rsb_duplicates_ovw )
	{
		RSB_DO_FLAG_SUBST(flags,RSB_FLAG_ALL_DUPLICATE_FLAGS,RSB_FLAG_DUPLICATES_KEEP_LAST);
	}
	else
	if( bsm->dupstra == blas_rsb_duplicates_sum )
	{
		RSB_DO_FLAG_SUBST(flags,RSB_FLAG_ALL_DUPLICATE_FLAGS,RSB_FLAG_DUPLICATES_SUM);
	}

	if(bsm->diag_type==blas_unit_diag)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT);

//	if(bsm->diag_type==blas_one_base) // this is complicated: should we keep these flags in both inner and outer interface ?
//		RSB_DO_FLAG_ADD(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
	switch(bsm->fmt_hint)
	{
#if RSB_WANT_SPARSE_BLAS_EXTENSIONS
		case(blas_rsb_rep_coo ): RSB_DO_FLAG_ADD(flags,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS); 
 break;
		case(blas_rsb_rep_csr ): 
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS); 
#if !RSB_OLD_COO_CRITERIA
			
			if (bsm->nnzin > bsm->coomatrix.nr)
			{
				/* Need enough nnz (nnz>nr) for csr to fit in coo arrays. 
				   FIXME: if duplicates within these nnzin (which may be <nnz) may break this condition.
 				*/
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS); 
				RSB_DO_FLAG_DEL(flags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS); 
			}
#endif
		break;
		case(blas_rsb_rep_rec ): RSB_DO_FLAG_ADD(flags,RSB_FLAG_QUAD_PARTITIONING); break;
		case(blas_rsb_rep_hwi ): RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES); break;
#endif /* RSB_WANT_SPARSE_BLAS_EXTENSIONS */
		case(blas_rsb_rep_rsb ):
		default:
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_DEFAULT_MATRIX_FLAGS);
	}
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS);
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS);

	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_QUAD_PARTITIONING);
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES);
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_COO_STORAGE);
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS);

	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES);		// NEW
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS);	// NEW
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES_COO);
	//RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_HALF_DETECTED_CACHE);	// EXPERIMENTAL
	return flags;
}

static rsb_blas_int_t rsb__BLAS_autotune( struct rsb_blas_sparse_matrix_t * bsm, rsb_blas_int_t pname )
{
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING 
	if( bsm == NULL )
	       goto err;

	if( bsm->type != blas_open_handle && bsm->type != blas_valid_handle )
	       goto err;

	if(RSB_BLAS_IS_ATPNAME_OFF(pname))
	{
		if( pname == blas_rsb_spmv_autotuning_off )
			bsm->opt_mvn_hint = bsm->opt_mvt_hint = RSB_SPB_THREADS_DEFAULT;
		if( pname == blas_rsb_spmv_n_autotuning_off )
			bsm->opt_mvn_hint = RSB_SPB_THREADS_DEFAULT;
		if( pname == blas_rsb_spmv_t_autotuning_off )
			bsm->opt_mvt_hint = RSB_SPB_THREADS_DEFAULT;
		goto done;
	}

	if(RSB_BLAS_IS_ATPNAME_ON(pname))
	{
#if RSB_TUNING_NEW_STYLE
		if( pname == blas_rsb_spmv_autotuning_on )
			bsm->opt_mvn_hint = /* bsm->opt_mvt_hint = */ RSB_SPB_THR_STR_AUTO;
		if( pname == blas_rsb_spmv_n_autotuning_on )
			bsm->opt_mvn_hint = RSB_SPB_THR_STR_AUTO;
		if( pname == blas_rsb_spmv_t_autotuning_on )
			bsm->opt_mvt_hint = RSB_SPB_THR_STR_AUTO;
		if( pname == blas_rsb_autotune_next_operation )
			bsm->opt_mvn_hint = bsm->opt_mvt_hint = RSB_SPB_THR_STR_AUTO_NEXTOP ;
#else /* RSB_TUNING_NEW_STYLE */
		if( pname == blas_rsb_spmv_autotuning_on )
			bsm->opt_mvn_hint = /* bsm->opt_mvt_hint = */ RSB_SPB_THREADS_AUTO;
		if( pname == blas_rsb_spmv_n_autotuning_on )
			bsm->opt_mvn_hint = RSB_SPB_THREADS_AUTO;
		if( pname == blas_rsb_spmv_t_autotuning_on )
			bsm->opt_mvt_hint = RSB_SPB_THREADS_AUTO;
		if( pname == blas_rsb_autotune_next_operation )
			bsm->opt_mvn_hint = bsm->opt_mvt_hint = RSB_SPB_THR_STR_AUTO_NEXTOP ;
#endif /* RSB_TUNING_NEW_STYLE */
		
		if( bsm->type == blas_open_handle )
			goto done;
		else
			;/* continue */
	}

	if( bsm->type == blas_valid_handle )
	{
		if( RSB_BLAS_IS_ATPNAME_ANY(pname) )
		{
			rsb_err_t errval = RSB_ERR_NO_ERROR;
			int nont = 0, tont = 0;
			const rsb_int_t mnt = rsb__set_num_threads(RSB_THREADS_GET_MAX);
			struct rsb_mtx_t * mtxOp = bsm->mtxAp;

			if( bsm->opt_mvt_hint == RSB_SPB_THREADS_AUTO )
			{
				errval = rsb__do_tune_spmm(NULL,NULL,&tont,2*mnt,10.0/mnt,RSB_TRANSPOSITION_T,NULL,mtxOp,1,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,NULL,0,NULL,NULL,0);
				if( RSB_SOME_ERROR(errval))
					goto err;
				bsm->opt_mvt_hint = tont;
			}

			if( bsm->opt_mvn_hint == RSB_SPB_THREADS_AUTO )
			{
				errval = rsb__do_tune_spmm(NULL,NULL,&nont,2*mnt,10.0/mnt,RSB_TRANSPOSITION_N,NULL,mtxOp,1,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,NULL,0,NULL,NULL,0);
				if( RSB_SOME_ERROR(errval)) 
					goto err;
				bsm->opt_mvn_hint = nont;
			}

#if RSB_TUNING_NEW_STYLE
			if( bsm->opt_mvt_hint == RSB_SPB_THR_STR_AUTO )
			{
				errval = rsb__tune_spxx(&mtxOp, NULL, &tont, RSB_SPBLAS_DEF_TUNING_ROUNDS, RSB_CONST_DEF_MS_AT_AUTO_STEPS, RSB_CONST_DEF_MS_AT_AUTO_STEPS, RSB_CONST_AT_OP_SAMPLES_MIN, RSB_CONST_AT_OP_SAMPLES_MAX, 0, RSB_TRANSPOSITION_T, NULL, NULL, 1, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, NULL, 0, NULL, NULL, 0, rsb_op_spmv, NULL, NULL, NULL, RSB_AUT0_TUNING_SILENT, NULL, NULL, NULL, NULL, NULL);
				bsm->opt_mvt_hint = tont;
			}

			if( bsm->opt_mvn_hint == RSB_SPB_THR_STR_AUTO )
			{
				errval = rsb__tune_spxx(&mtxOp, NULL, &nont, RSB_SPBLAS_DEF_TUNING_ROUNDS, RSB_CONST_DEF_MS_AT_AUTO_STEPS, RSB_CONST_DEF_MS_AT_AUTO_STEPS, RSB_CONST_AT_OP_SAMPLES_MIN, RSB_CONST_AT_OP_SAMPLES_MAX, 0, RSB_TRANSPOSITION_N, NULL, NULL, 1, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, NULL, 0, NULL, NULL, 0, rsb_op_spmv, NULL, NULL, NULL, RSB_AUT0_TUNING_SILENT, NULL, NULL, NULL, NULL, NULL);
				bsm->opt_mvn_hint = nont;
			}

			if( bsm->mtxAp != mtxOp && mtxOp )
		       	{
#if !RSB_AT_DESTROYS_MTX
			       	RSB_CONDITIONAL_FREE(bsm->mtxAp);
#endif /* RSB_AT_DESTROYS_MTX */
				bsm->mtxAp = mtxOp;
		       	}

			if( bsm->opt_mvn_hint == RSB_SPB_THR_STR_AUTO_NEXTOP )
				; /* will optimize at operation time */

			if( RSB_SOME_ERROR(errval))
				goto err;
#endif /* RSB_TUNING_NEW_STYLE */
			return RSB_BLAS_NO_ERROR;
		}
	}
done:
	return RSB_BLAS_NO_ERROR;
err:
	return RSB_BLAS_ERROR;
#else
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	return RSB_BLAS_ERROR_UNSUPPORTED;
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	return RSB_BLAS_NO_ERROR; /* bogus success */
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#endif
}

static blas_sparse_matrix rsb__BLAS_new_matrix_finish(struct rsb_blas_sparse_matrix_t * bsm, const rsb_flags_t * flagsp)
{
	/**
	 \ingroup gr_internals
	 shall finish the build status of a new matrix  
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;
	rsb_flags_t flags;

	if(!bsm || bsm->type != blas_open_handle)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	/* trim the extra storage */
	if(rsb__BLAS_new_matrix_expand_store(bsm,bsm->nnzin)==RSB_BLAS_INVALID_VAL)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	flags = flagsp ? *flagsp:rsb__BLAS_new_matrix_finish_flags(bsm);
	bsm->mtxAp = rsb__mtx_alloc_inner(bsm->coomatrix.VA,bsm->coomatrix.IA,bsm->coomatrix.JA,bsm->coomatrix.nnz,0,0,bsm->coomatrix.typecode,bsm->coomatrix.nr,bsm->coomatrix.nc,bsm->k,bsm->l,flags,&errval);
	if(!bsm->mtxAp)
	{
		rsb__destroy_coo_matrix_t(&bsm->coomatrix);
		//bsm->coomatrix.IA=bsm->coomatrix.JA=bsm->coomatrix.VA=NULL;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
#if RSB_ALLOW_STDOUT
	if(0)
		RSB_STDOUT("sparse blas allocated (%zd x %zd) @ %p with flags 0x%x (coo:%d, csr:%d), storage: %x\n",
			(rsb_printf_int_t)bsm->mtxAp->nr, (rsb_printf_int_t)bsm->mtxAp->nc, (void*)bsm->mtxAp, bsm->mtxAp->flags,
			RSB_DO_FLAG_HAS(bsm->mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE),
			RSB_DO_FLAG_HAS(bsm->mtxAp->flags,RSB_FLAG_WANT_BCSS_STORAGE),
			bsm->mtxAp->matrix_storage
			);
#endif /* RSB_ALLOW_STDOUT */

/*	rsb__do_print_matrix_stats(bsm->mtxAp,RSB_CONST_DUMP_DEFAULT); */
/*	rsb__do_print_matrix_stats(bsm->mtxAp,RSB_CONST_DUMP_TIMES); */
/*	rsb__do_print_matrix_stats(bsm->mtxAp,RSB_CONST_DUMP_RECURSION); */

	bsm->type = blas_valid_handle;

#if RSB_BLAS_ALLOW_MTX_UPD
	RSB_DO_FLAG_SUBST(bsm->mtxAp->flags,RSB_FLAG_ALL_DUPLICATE_FLAGS,RSB_FLAG_DUPLICATES_KEEP_LAST);
#endif /* RSB_BLAS_ALLOW_MTX_UPD */

#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING 
	if( bsm->opt_mvn_hint != RSB_SPB_THREADS_DEFAULT || bsm->opt_mvt_hint != RSB_SPB_THREADS_DEFAULT )
		rsb__BLAS_autotune( bsm, RSB_ATPNAME_ANY ); /* FIXME: for now, no error reporting here */
#endif

	retval = bsm->handle;
err:
	return retval;
}

static blas_sparse_matrix rsb__BLAS_matrix_destroy(blas_sparse_matrix handle)
{
	/**
	 \ingroup gr_internals
	 shall destroy a matrix  
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if( (bsm = rsb__BLAS_matrix_retrieve(handle) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	RSB_DEBUG_ASSERT(bsm);
	switch(bsm->type)
	{
		case blas_new_handle:
		/* nothing to do  */
		break;
		case blas_valid_handle:
			rsb__do_mtx_free(bsm->mtxAp);
		break;
		case blas_open_handle:
			rsb__destroy_coo_matrix_t(&bsm->coomatrix);
		break;
		default:
		/* nothing to do  */
		break;
	}

	RSB_CONDITIONAL_FREE(bsm->rbp);
	RSB_CONDITIONAL_FREE(bsm->cbp);

	retval = rsb__BLAS_handle_free(handle);
err:
	return retval;
}

blas_sparse_matrix rsb__BLAS_Xuscr_begin( rsb_blas_int_t m, rsb_blas_int_t n, rsb_type_t typecode)
{
	/**
	 \ingroup gr_internals
	 \rsb_spblasl2_cr_begin_msg
	 */
	blas_sparse_matrix retval = rsb__BLAS_new_matrix_begin(m,n,1+RSB_MAX(m,n), typecode,1,1,NULL,NULL);
	return retval;
}

blas_sparse_matrix rsb__BLAS_Xuscr_block_begin( rsb_blas_int_t Mb, rsb_blas_int_t Nb, rsb_blas_int_t k, rsb_blas_int_t l, rsb_type_t typecode)
{
	/**
	 \ingroup gr_internals
	 */
	blas_sparse_matrix retval = rsb__BLAS_new_matrix_begin(Mb*k,Nb*l,1+RSB_MAX(Mb*k,Nb*l), typecode,k,l,NULL,NULL);
	return retval;
}

blas_sparse_matrix rsb__BLAS_Xuscr_variable_block_begin( rsb_blas_int_t Mb, rsb_blas_int_t Nb, const rsb_blas_int_t *k, const rsb_blas_int_t *l, rsb_type_t typecode)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if(!k || !l)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	retval = rsb__BLAS_new_matrix_begin(0,0,0,typecode,Mb,Nb,k,l);
err:
	return retval;
}

rsb_blas_int_t rsb__BLAS_Xuscr_insert_entry( blas_sparse_matrix A, const void * valp, rsb_blas_int_t i, rsb_blas_int_t j )
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( bsm->type != blas_open_handle)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	retval = rsb__BLAS_new_matrix_insert_entries(bsm,1,valp,&i,&j);
err:
	return retval;
}

rsb_blas_int_t rsb__BLAS_Xuscr_insert_entries( blas_sparse_matrix A, rsb_blas_int_t nz, const void * val, const rsb_blas_int_t *indx, const rsb_blas_int_t *jndx )
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	if( ( bsm->type != blas_open_handle )
#if RSB_BLAS_ALLOW_MTX_UPD
			&& ( bsm->type != blas_valid_handle )
#endif /* RSB_BLAS_ALLOW_MTX_UPD */
		)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	retval = rsb__BLAS_new_matrix_insert_entries( bsm, nz, val, indx, jndx );
err:
	return retval;
}

rsb_blas_int_t rsb__BLAS_Xuscr_insert_col( blas_sparse_matrix A, rsb_blas_int_t j, rsb_blas_int_t nz, const void * val, const rsb_blas_int_t *indx )
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( bsm->type != blas_open_handle)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( RSB_INVALID_COO_INDEX(j) || RSB_INVALID_NNZ_INDEX(nz) )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	retval = rsb__BLAS_new_matrix_insert_col(bsm, j, nz, val, indx );
err:
	return retval;
}

rsb_blas_int_t rsb__BLAS_Xuscr_insert_row( blas_sparse_matrix A, rsb_blas_int_t i, rsb_blas_int_t nz, const void * val, const rsb_blas_int_t *jndx )
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( bsm->type != blas_open_handle)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( RSB_INVALID_COO_INDEX(i) || RSB_INVALID_NNZ_INDEX(nz) )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	retval = rsb__BLAS_new_matrix_insert_row(bsm, i, nz, val, jndx );
err:
	return retval;
}

rsb_blas_int_t rsb__BLAS_Xuscr_insert_clique( blas_sparse_matrix A, const rsb_blas_int_t k, const rsb_blas_int_t l, const void * val, const rsb_blas_int_t row_stride, const rsb_blas_int_t col_stride, const rsb_blas_int_t *indx, const rsb_blas_int_t *jndx )
{
	/**
	 \ingroup gr_internals
	FIXME: this prototype does not respect the standard !
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( bsm->type != blas_open_handle)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	retval = rsb__BLAS_new_matrix_insert_clique(bsm, k, l, val, row_stride, col_stride, indx, jndx );
err:
	return retval;
}

rsb_blas_int_t rsb__BLAS_Xuscr_insert_block( blas_sparse_matrix A, const void * val, rsb_blas_int_t row_stride, rsb_blas_int_t col_stride, rsb_blas_int_t i, rsb_blas_int_t j )
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retval = RSB_BLAS_INVALID_VAL;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( bsm->type != blas_open_handle)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	retval = rsb__BLAS_new_matrix_insert_block(bsm, val, row_stride, col_stride, i, j);
err:
	return retval;
}

               /* Completion of Construction Routines */

rsb_blas_int_t rsb__BLAS_Xuscr_end_flagged( blas_sparse_matrix A, const rsb_flags_t*flagsp)
{
	/**
	 \ingroup gr_internals
	 matrix finish contruction phase
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t errval = RSB_BLAS_ERROR; 

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( bsm->type != blas_open_handle )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	errval = rsb__BLAS_new_matrix_finish(bsm,flagsp) == A ? RSB_BLAS_NO_ERROR : RSB_BLAS_ERROR;
err:
	return errval;

}

rsb_blas_int_t rsb__BLAS_Xuscr_end( blas_sparse_matrix A )
{
	rsb_blas_int_t retval = rsb__BLAS_Xuscr_end_flagged(A,NULL);
	return retval;
}

               /* Matrix Property Routines */

rsb_blas_int_t rsb__BLAS_usgp( blas_sparse_matrix A, rsb_blas_int_t pname )
{
	/**
	 \ingroup gr_internals
	 matrix property get
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_blas_int_t retcode = RSB_BLAS_ERROR;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		goto ret;

	switch (pname)
	{
		/*  */
		case (blas_num_rows) : retcode=bsm->coomatrix.nr; break;
		case (blas_num_cols) : retcode=bsm->coomatrix.nc; break;
		case (blas_num_nonzeros) : if (bsm->mtxAp != NULL) rsb__BLAS_Xusget_matrix_nnz(A,&retcode) ; else retcode = bsm->coomatrix.nnz; break;
		//case (blas_num_nonzeros) : retcode=bsm->coomatrix.nnz; break;
		/*  */
		case (blas_complex) : retcode=bsm->field == blas_complex; break;
		case (blas_real) : retcode=bsm->field == blas_real; break;
		case (blas_single_precision) : retcode=bsm->fprecision == blas_single_precision ? 1:0; break;
		case (blas_double_precision) : retcode=bsm->fprecision == blas_double_precision ? 1:0; break;
		/*  */
		case (blas_triangular) : retcode=( /*bsm->symmetry == blas_triangular || */bsm->symmetry == blas_lower_triangular || bsm->symmetry == blas_upper_triangular ); break;
		case (blas_lower_triangular) : retcode=bsm->symmetry == blas_lower_triangular; break;
		case (blas_upper_triangular) : retcode=bsm->symmetry == blas_upper_triangular; break;
		/*  */
		case (blas_general) : retcode=bsm->symmetry == blas_general; break;
		case (blas_lower_symmetric) : retcode=( bsm->symmetry == blas_lower_symmetric ); break;
		case (blas_upper_symmetric) : retcode=( bsm->symmetry == blas_upper_symmetric ); break;
		case (blas_symmetric) : retcode=( /*bsm->symmetry == blas_symmetric || */bsm->symmetry == blas_lower_symmetric || bsm->symmetry == blas_upper_symmetric ); break;
		/* case (blas_hermitian) : retcode=bsm->symmetry == blas_hermitian; break; */
		case (blas_hermitian) : retcode=( /*bsm->symmetry == blas_hermitian || */bsm->symmetry == blas_lower_hermitian || bsm->symmetry == blas_upper_hermitian ); break;
		/*  */
		case (blas_zero_base) : retcode=bsm->base == blas_zero_base; break;
		case (blas_one_base) : retcode=bsm->base == blas_one_base; break;
		/*  */
		case (blas_rowmajor) : retcode=bsm->order == blas_rowmajor; break;
		case (blas_colmajor) : retcode=bsm->order == blas_colmajor; break;
		/*  */
		case (blas_new_handle) : retcode=bsm->type == blas_new_handle; break;
		case (blas_open_handle) : retcode=bsm->type == blas_open_handle; break;
		case (blas_valid_handle) : retcode=bsm->type == blas_valid_handle;
		case (blas_invalid_handle) : retcode=bsm->type != blas_valid_handle; /* FIXME */
		/* the following occur in the NIST 1.02 version */
		/* case (blas_unassembled) : retcode = (bsm->mtxAp == NULL); */
		case (blas_unassembled) : retcode=bsm->type != blas_valid_handle;
		case (blas_regular) : retcode = 0;
		case (blas_irregular) : retcode = 1;
		case (blas_block) : retcode = 0;

		break;
#if RSB_WANT_SPARSE_BLAS_EXTENSIONS
		case (blas_rsb_duplicates_ovw) :
		case (blas_rsb_duplicates_sum) :
		retcode = bsm->dupstra ;
		break;
		case (blas_rsb_rep_csr) :
		case (blas_rsb_rep_coo) :
		case (blas_rsb_rep_rsb) :
		case (blas_rsb_rep_hwi) :
		case (blas_rsb_rep_rec) :
		retcode = bsm->fmt_hint;
		break;
#endif /* RSB_WANT_SPARSE_BLAS_EXTENSIONS */
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING 
		case (blas_rsb_spmv_autotuning_on ):
		case (blas_rsb_spmv_autotuning_off ):
		/* NOTE: we return values for untransposed here. */
		case (blas_rsb_spmv_n_autotuning_off ):
		case (blas_rsb_spmv_n_autotuning_on ):
		case (blas_rsb_autotune_next_operation):
			retcode = bsm->opt_mvn_hint > 0 ?  bsm->opt_mvn_hint : 0;
		break;
		case (blas_rsb_spmv_t_autotuning_off ):
		case (blas_rsb_spmv_t_autotuning_on ):
			retcode = bsm->opt_mvt_hint > 0 ?  bsm->opt_mvt_hint : 0;
		break;
#else
		case (blas_rsb_spmv_autotuning_on ):
		case (blas_rsb_spmv_autotuning_off ):
		case (blas_rsb_spmv_n_autotuning_off ):
		case (blas_rsb_spmv_n_autotuning_on ):
		case (blas_rsb_spmv_t_autotuning_off ):
		case (blas_rsb_spmv_t_autotuning_on ):
		case (blas_rsb_autotune_next_operation):
#if RSB_WANT_OMP_RECURSIVE_KERNELS
			retcode = RSB_BLAS_ERROR_UNSUPPORTED;
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
			retcode = RSB_CONST_MIN_SUPPORTED_CORES;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		break;
#endif /* RSB_BLAS_WANT_EXPERIMENTAL_TUNING */
	}
ret:
	return retcode;
}

rsb_blas_int_t rsb__BLAS_Xusget_rows_nnz( blas_sparse_matrix A, rsb_blas_int_t fr, rsb_blas_int_t lr, rsb_blas_int_t * nnzp)
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( (bsm = rsb__BLAS_matrix_retrieve(A) ) == NULL )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	if( bsm->type != blas_valid_handle )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	if( bsm->base == blas_one_base )
		--fr,--lr;

	*nnzp = rsb__dodo_get_rows_nnz(bsm->mtxAp,fr,lr,RSB_FLAG_C_INDICES_INTERFACE,&errval);
	return RSB_ERROR_TO_BLAS_ERROR(errval);
err:
	return RSB_BLAS_INVALID_VAL;
}

rsb_blas_int_t rsb__BLAS_ussp( blas_sparse_matrix A, rsb_blas_int_t pname )
{
	/**
	 \ingroup gr_internals
	 */
	/* TODO: table 3.4 in the Sparse BLAS standard from 2001 has values not present in the reference implementation ! */
	struct rsb_blas_sparse_matrix_t * bsm = NULL;

	bsm = rsb__BLAS_matrix_retrieve(A);

	if( RSB_BLAS_IS_ATPNAME(pname) )
	{
		rsb_blas_int_t ret = rsb__BLAS_autotune( bsm, pname ); /* This property can be set at any time */
		return ret;
	}

#if RSB_BLAS_ALLOW_MTX_UPD
	if( ( bsm != NULL ) && ( bsm->type == blas_valid_handle ) && ( bsm->mtxAp != NULL ) )
	{
		switch (pname)
		{
			case (blas_rsb_duplicates_ovw) :
			RSB_DO_FLAG_SUBST(bsm->mtxAp->flags,RSB_FLAG_ALL_DUPLICATE_FLAGS,RSB_FLAG_DUPLICATES_KEEP_LAST);
			bsm->dupstra = pname;
			break;
			case (blas_rsb_duplicates_sum) :
			RSB_DO_FLAG_SUBST(bsm->mtxAp->flags,RSB_FLAG_ALL_DUPLICATE_FLAGS,RSB_FLAG_DUPLICATES_SUM);
			bsm->dupstra = pname;
			break;
			default:
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
		goto ok;
	}
#endif /* RSB_BLAS_ALLOW_MTX_UPD */

	if( ( bsm == NULL ) || bsm->type != blas_open_handle )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	switch (pname)
	{	
		/* Extension properties are ones influencing assembly. */
#if RSB_WANT_SPARSE_BLAS_EXTENSIONS
		case (blas_rsb_duplicates_ovw) :
		case (blas_rsb_duplicates_sum) :
			bsm->dupstra = pname;
			goto ok;
		break;
		case (blas_rsb_rep_csr) :
		case (blas_rsb_rep_coo) :
		case (blas_rsb_rep_rsb) :
		case (blas_rsb_rep_hwi) :
		case (blas_rsb_rep_rec) :
			bsm->fmt_hint = pname;
			goto ok;
		break;
		default: RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
#endif /* RSB_WANT_SPARSE_BLAS_EXTENSIONS */
	}

	if( bsm->nnzin != 0 )
	{
		/*
		 According to [dv_2002]:
		  "Calls to USSP should be made after a call to the BEGIN routine but before
		   the first call to an INSERT routine for the same handle."
		   */
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	switch (pname)
	{
		case (blas_general)       : bsm->symmetry = blas_general; break;
		case (blas_one_base)       : bsm->off=1; bsm->base=blas_one_base; break;
		case (blas_zero_base)       : bsm->base=blas_zero_base; break;

		case (blas_non_unit_diag) : bsm->diag_type=blas_non_unit_diag; break;
		case (blas_unit_diag) : bsm->diag_type=blas_unit_diag; break;
#if 1
		case (blas_complex) :
		case (blas_real)	:
		case (blas_double_precision) :
		case (blas_single_precision) :
		/* FIXME : should return error only on type different from the already set one */
		/* FIXME : unless what the standard mandates is really conversion ?!  */
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		break;
#endif
#if 0
		/* FIXME: the following do not exist, but may be useful  */
		case (blas_no_repeated_indices ) :
			return RSB_BLAS_ERROR_UNIMPLEMENTED;
		case (blas_repeated_indices ) :
			return RSB_BLAS_ERROR_UNIMPLEMENTED;
		break;
#endif
		case (blas_triangular) : bsm->symmetry=blas_triangular; break;
		case (blas_lower_triangular) : bsm->symmetry=blas_lower_triangular; break;
		case (blas_upper_triangular) : bsm->symmetry=blas_upper_triangular; break;
		case (blas_symmetric)       : return RSB_BLAS_ERROR_WRONG_USGP_ARG;	/* TODO */ break;
		case (blas_lower_symmetric) : bsm->symmetry=blas_lower_symmetric; break;
		case (blas_upper_symmetric) : bsm->symmetry=blas_upper_symmetric; break;
		case (blas_hermitian) : return RSB_BLAS_ERROR_WRONG_USGP_ARG;	/* TODO */ break;
		case (blas_lower_hermitian) : bsm->symmetry=blas_lower_hermitian; break;
		case (blas_upper_hermitian) : bsm->symmetry=blas_upper_hermitian; break;
		case (blas_rowmajor) :
			bsm->order=blas_rowmajor;
		break;
		case (blas_colmajor) :
			return RSB_BLAS_ERROR_UNIMPLEMENTED;	/* TODO */
		break;
		case (blas_regular) : bsm->sparsity_optimization_type=blas_regular; break;
#if 1
		/* FIXME: we interpret the following as hints, but in future we may use them */
		case (blas_block) :
		case (blas_irregular) :
		case (blas_unassembled) :
		break;
#endif
		/* ... */
		default:
		return RSB_BLAS_ERROR;
	}
ok:
	return RSB_BLAS_NO_ERROR;
err:
	return RSB_BLAS_ERROR;
}

               /* Destruction Routine */

rsb_blas_int_t rsb__BLAS_Xusds( blas_sparse_matrix A )
{
	/**
	 */
	/*
	 * Destroys the given matrix
	 * */
	const rsb_blas_int_t res = ( rsb__BLAS_matrix_destroy(A) == RSB_BLAS_INVALID_VAL ) ? RSB_BLAS_ERROR : RSB_BLAS_NO_ERROR;

	return res;
}

rsb_trans_t rsb__do_psblas_trans_to_rsb_trans(const char trans)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_trans_t rtrans = RSB_INVALID_FLAGS;

	switch(trans)
	{
		case(RSB_PSBLAS_TRANS_N):
		rtrans = RSB_TRANSPOSITION_N;
		break;
		case(RSB_PSBLAS_TRANS_T):
		rtrans = RSB_TRANSPOSITION_T;
		break;
		case(RSB_PSBLAS_TRANS_C):
		rtrans = RSB_TRANSPOSITION_C;
		break;
	}
	return rtrans;
}

rsb_trans_t rsb__blas_trans_to_rsb_trans(enum blas_trans_type trans)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_trans_t rtrans = RSB_INVALID_FLAGS;

	switch(trans)
	{
		case(blas_no_trans):
		rtrans = RSB_TRANSPOSITION_N;
		break;
		case(blas_trans):
		rtrans = RSB_TRANSPOSITION_T;
		break;
		case(blas_conj_trans):
		rtrans = RSB_TRANSPOSITION_C;
		break;
	}
	return rtrans;
}

rsb_order_t rsb__blas_order_to_rsb_order(enum blas_order_type order)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_order_t rorder = RSB_FLAG_WANT_ROW_MAJOR_ORDER;

	switch(order)
	{
		case(blas_colmajor):
		rorder = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
		break;
		case(blas_rowmajor):
		rorder = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
		break;
	}
	return rorder;
}

int rsb__BLAS_Xusrows_scale(blas_sparse_matrix A,const void * d,enum blas_trans_type trans)
{
	/**
	 \ingroup gr_internals
	 \rsb_spblasl2e_usrows_scale_msg
	 */
	struct rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	rsb_err_t errval = rsb__do_scal(mtxAp,d,rsb__blas_trans_to_rsb_trans(trans));

	return RSB_ERROR_TO_BLAS_ERROR(errval);
}

int rsb__BLAS_Xusget_diag(blas_sparse_matrix A,void * d)
{
	/**
	 \ingroup gr_internals
	 \rsb_spblasl2e_usget_diag_msg
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const struct rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);

	errval = rsb__do_matrix_compute(mtxAp, d, RSB_EXTF_DIAG);
	return RSB_ERROR_TO_BLAS_ERROR(errval);
}

int rsb__BLAS_Xusget_rows_sparse(blas_sparse_matrix A,void *  VA, rsb_blas_int_t * IA, rsb_blas_int_t * JA, rsb_blas_int_t * nnz, rsb_blas_int_t fr, rsb_blas_int_t lr)
{
	/**
	 \ingroup gr_internals
	 \rsb_spblasl2e_usget_rows_sparse_msg
	 */
        rsb_err_t errval = RSB_ERR_NO_ERROR;
	const struct rsb_blas_sparse_matrix_t * bsm = rsb__BLAS_matrix_retrieve(A);

	if(!bsm || !bsm->mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	errval = rsb__do_get_rows_sparse(RSB_TRANSPOSITION_N,NULL,bsm->mtxAp,VA,IA,JA,fr,lr,nnz,
			RSB_FLAG_SORT_INPUT|(
			bsm->base == blas_one_base ? RSB_FLAG_FORTRAN_INDICES_INTERFACE : RSB_FLAG_NOFLAGS));
err:
	return RSB_ERROR_TO_BLAS_ERROR(errval);
}

int rsb__BLAS_Xusget_matrix_nnz(blas_sparse_matrix A, rsb_blas_int_t * nnzAp)
{
	/**
	 \ingroup gr_internals
	 \rsb_spblasl2e_usget_matrix_nnz_msg
	 */
	const struct rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);

	if(!mtxAp)
		return RSB_ERROR_TO_BLAS_ERROR(RSB_ERR_BADARGS);
	*nnzAp = mtxAp->nnz;
	return RSB_BLAS_NO_ERROR;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__BLAS_Xusget_rows_sums(blas_sparse_matrix A, void * rs, enum blas_trans_type trans)
{
	/**
	 \ingroup gr_internals
	 */
	struct rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	rsb_err_t errval = rsb__do_rowssums(mtxAp,rsb__blas_trans_to_rsb_trans(trans),rs);

	return RSB_ERROR_TO_BLAS_ERROR(errval);
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

int rsb__BLAS_Xusget_infinity_norm(blas_sparse_matrix A, void * in, enum blas_trans_type trans)
{
	/**
	 \ingroup gr_internals
	 */
	const struct rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	/* rsb_err_t errval = rsb__do_matrix_norm(mtxAp,in,RSB_EXTF_NORM_INF); */

	rsb_err_t errval = rsb__do_matrix_norm(mtxAp,in,RSB_EXTF_NORM_INF);
	return RSB_ERROR_TO_BLAS_ERROR(errval);
}

int rsb__BLAS_Xusset_elements(blas_sparse_matrix A, const rsb_blas_int_t * ia, const rsb_blas_int_t *ja, const void *  va, rsb_blas_int_t nnz)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;

	mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	if(!mtxAp)
	{
		errval = RSB_ERROR_TO_BLAS_ERROR(RSB_ERR_BADARGS);
		goto err;
	}

	errval = rsb__do_set_coo_elements(mtxAp,va,ia,ja,nnz);
err:
	return RSB_ERROR_TO_BLAS_ERROR(errval);
}

int rsb__BLAS_Xusset_element(blas_sparse_matrix A,rsb_blas_int_t i, rsb_blas_int_t j, const void * v)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);

	errval = rsb__do_set_coo_element(mtxAp,v,i,j);
	return RSB_ERROR_TO_BLAS_ERROR(errval);
}

int rsb__BLAS_Xusget_element(blas_sparse_matrix A,rsb_blas_int_t i, rsb_blas_int_t j, void * v)
{
	/**
	 \ingroup gr_internals
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if 1
	struct rsb_blas_sparse_matrix_t * bsm = NULL;
	if( ( bsm = rsb__BLAS_matrix_retrieve(A)) != NULL && bsm->base == blas_one_base ) 
		i-=1, j-=1;
	errval = rsb__do_get_coo_element(bsm->mtxAp,v,i,j);
#else
	struct rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	errval = rsb__do_get_coo_element(mtxAp,v,i,j);
#endif
	return RSB_ERROR_TO_BLAS_ERROR(errval);
}

int rsb__BLAS_Xusmm(enum blas_trans_type transA, const void * alphap, blas_sparse_matrix A, const void * b, rsb_blas_int_t ldb, const void * betap, void * c, rsb_blas_int_t ldc, rsb_blas_int_t nrhs, enum blas_order_type order)
{
	/**
	 	Multiplies by multivector, accumulating in a multivector and scaling it.
	*/
	//const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	const rsb_trans_t trans = rsb__blas_trans_to_rsb_trans(transA);
	int brv = RSB_BLAS_ERROR;
	struct rsb_blas_sparse_matrix_t * bsm = rsb__BLAS_matrix_retrieve(A);
	rsb_order_t rorder = rsb__blas_order_to_rsb_order(order);

	if(bsm)
	{
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING
		rsb_thread_t ornt = rsb_get_num_threads();
		if((bsm->opt_mvn_hint) == RSB_SPB_THR_STR_AUTO_NEXTOP )
			ornt = -ornt; /* want threads tuning */
		RSB_SPB_AT_OP(bsm->mtxAp,ornt,bsm->opt_mvn_hint,nrhs,rorder,alphap,betap,c,b,ldc,ldb,rsb_op_spmv)
		RSB_SPB_THREADS_PUSH
#endif /* RSB_BLAS_WANT_EXPERIMENTAL_TUNING */
		if(!bsm->mtxAp)
			goto err;
		brv = RSB_ERROR_TO_BLAS_ERROR(rsb__do_spmm_general(bsm->mtxAp,b,c,alphap,betap,1,1,trans,RSB_OP_FLAG_DEFAULT,rorder,nrhs,ldb,ldc));
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING
		RSB_SPB_THREADS_POP
#endif /* RSB_BLAS_WANT_EXPERIMENTAL_TUNING */
	}
err:
	return brv;
}

int rsb__BLAS_Xusmv(enum blas_trans_type transA, const void * alphap, blas_sparse_matrix A, const void * Xp, rsb_blas_int_t incX, const void * betap, void * Yp, rsb_blas_int_t incY)
{
	/**
	 	Multiplies by a vector, accumulating in a vector and scaling it.
		\f$y \leftarrow \alpha A   x + \beta y  ,\f$
		\f$y \leftarrow \alpha A^T x + \beta y,\f$
		\f$y \leftarrow \alpha A^H x + \beta y\f$
	*/
	//const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	const rsb_trans_t trans = rsb__blas_trans_to_rsb_trans(transA);
	int brv = RSB_BLAS_ERROR;
	struct rsb_blas_sparse_matrix_t * bsm = rsb__BLAS_matrix_retrieve(A);

	if(bsm)
	{
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING
		rsb_thread_t ornt = rsb_get_num_threads();
		if((bsm->opt_mvn_hint) == RSB_SPB_THR_STR_AUTO_NEXTOP )
			ornt = -ornt; /* want threads tuning */
		RSB_SPB_AT_OP(bsm->mtxAp,ornt,bsm->opt_mvn_hint,1,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,alphap,betap,Yp,Xp,0,0,rsb_op_spmv)
		RSB_SPB_THREADS_PUSH
#endif /* RSB_BLAS_WANT_EXPERIMENTAL_TUNING */
		if(!bsm->mtxAp)
			goto err;
		brv = RSB_ERROR_TO_BLAS_ERROR(rsb__do_spmv_general(trans,alphap,bsm->mtxAp,Xp,incX,betap,Yp,incY,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS));
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING
		RSB_SPB_THREADS_POP
#endif /* RSB_BLAS_WANT_EXPERIMENTAL_TUNING */
	}
err:
	return brv;
}

int rsb__BLAS_Xussv(enum blas_trans_type transT, void * alpha, blas_sparse_matrix T, void * x, rsb_blas_int_t incx)
{
	/**
	 	Solves triangular system by a vector, scaling the result.
		 \f$x \leftarrow \alpha T^{-1}x,\f$
		 \f$x \leftarrow \alpha T^{-T}x,\f$
		 \f$x \leftarrow \alpha T^{-H}x\f$
	*/
	const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(T);
	const rsb_trans_t trans = rsb__blas_trans_to_rsb_trans(transT);

	return RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsv(trans,alpha,mtxAp,x,incx,x,incx));
}

static rsb_flags_t rsb__flags_from_props(rsb_bool_t is_hermitian, rsb_bool_t is_symmetric, rsb_bool_t is_lower, rsb_bool_t is_upper, rsb_type_t typecode)
{
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

	if(!RSB_IS_MATRIX_TYPE_COMPLEX(typecode) && ( is_hermitian == RSB_BOOL_TRUE ) )
	{
		is_hermitian = RSB_BOOL_FALSE;
		is_symmetric = RSB_BOOL_TRUE;
	}
	if(is_hermitian == RSB_BOOL_TRUE && !RSB_EXPERIMENTAL_EXPAND_SYMMETRIC_MATRICES_BY_DEFAULT)
	{
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_HERMITIAN);
	}
	if(is_symmetric == RSB_BOOL_TRUE && !RSB_EXPERIMENTAL_EXPAND_SYMMETRIC_MATRICES_BY_DEFAULT)
	{
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_SYMMETRIC);
	}

	if( (!RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER)) && (!RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER)) )
	{
		if(is_upper)
 			RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
		if(is_lower)
 			RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
	}
	if( RSB_NAND( RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER), RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER)) )
 			RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
	return flags;
}

static rsb_blas_int_t rsb__mtx_flags_usgp( blas_sparse_matrix A, rsb_flags_t flags )
{
	/**
	 \ingroup gr_internals
	 matrix property set
	TODO: missing checks.
	 */
	rsb_blas_int_t retcode = RSB_BLAS_ERROR;

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_TRIANGULAR))
		retcode = rsb__BLAS_ussp( A, blas_upper_triangular);

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_TRIANGULAR))
		retcode = rsb__BLAS_ussp( A, blas_lower_triangular);

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_HERMITIAN))
		retcode = rsb__BLAS_ussp( A, blas_upper_hermitian);

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_HERMITIAN))
		retcode = rsb__BLAS_ussp( A, blas_lower_hermitian);

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_SYMMETRIC))
		retcode = rsb__BLAS_ussp( A, blas_upper_symmetric);

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_SYMMETRIC))
		retcode = rsb__BLAS_ussp( A, blas_lower_symmetric);

	/* blas_triangular, blas_general, blas_symmetric, blas_hermitian, blas_zero_base, blas_one_base, blas_rowmajor, blas_colmajor */
	return retcode;
}

/* The following shall be documented. */
blas_sparse_matrix rsb__load_spblas_matrix_file_as_matrix_market(const rsb_char_t * filename, rsb_type_t typecode )
{
	/**
	 	\ingroup gr_internals
	 	Loads a BLAS Sparse matrix from a Matrix Market file.
		This is a \librsb extension.

		Sets either blas_upper_triangular, blas_lower_triangular, blas_upper_hermitian, blas_lower_hermitian, blas_upper_symmetric or blas_lower_symmetric property according to the loaded file.
	 */
	struct rsb_coo_mtx_t coo;
	blas_sparse_matrix A = blas_invalid_handle /*RSB_BLAS_INVALID_VAL*/;
	rsb_bool_t is_symmetric = /*RSB_BOOL_MAYBE*/RSB_BOOL_FALSE;
	rsb_bool_t is_hermitian = /*RSB_BOOL_MAYBE*/RSB_BOOL_FALSE;
	rsb_bool_t is_pattern = RSB_BOOL_MAYBE;
	rsb_bool_t is_lower = RSB_BOOL_MAYBE, is_upper = RSB_BOOL_MAYBE;
	rsb_bool_t is_vector = RSB_BOOL_FALSE;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_MATRIX_FLAGS | RSB_FLAG_NOFLAGS;

	RSB_BZERO_P(&coo);
	coo.typecode = typecode;
	/* coo.typecode = RSB_NUMERICAL_TYPE_INVALID_TYPE; */

	if( rsb__util_mm_info_matrix_f(filename,&coo.nr,&coo.nc,&coo.nnz,&coo.typecode,&is_symmetric,&is_hermitian,&is_pattern,&is_lower,&is_upper,&is_vector) )
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	A = rsb__BLAS_Xuscr_begin(coo.nr,coo.nc,coo.typecode);
	if( A == RSB_BLAS_INVALID_VAL )
		RSB_PERR_GOTO(derr,RSB_ERRM_ES) /* TODO: this particular case should be made illegal. */

	if( A == blas_invalid_handle )
		RSB_PERR_GOTO(derr,RSB_ERRM_ES) /* TODO: this particular case should be made official. */

	if( RSB_SOME_ERROR(rsb__util_mm_load_matrix_f(filename,&coo.IA,&coo.JA,&coo.VA,&coo.nr,&coo.nc,&coo.nnz,coo.typecode,RSB_FLAG_NOFLAGS,&is_lower,&is_upper)) )
		RSB_PERR_GOTO(derr,RSB_ERRM_ES)

	flags = rsb__flags_from_props(is_hermitian, is_symmetric, is_lower, is_upper, typecode); /* redundancy: rsb__BLAS_Xuscr_end_flagged does not use these */
	rsb__mtx_flags_usgp(A, flags);
	
	if( rsb__BLAS_Xuscr_insert_entries(A,coo.nnz,coo.VA,coo.IA,coo.JA) == RSB_BLAS_INVALID_VAL )
		RSB_PERR_GOTO(derr,RSB_ERRM_ES)

	if( rsb__BLAS_Xuscr_end_flagged(A,&flags) /*rsb__BLAS_Xuscr_end(A)*/ == RSB_BLAS_INVALID_VAL )
		RSB_PERR_GOTO(derr,RSB_ERRM_ES)

	goto ok;
derr:
	/* FIXME: missing proper deallocation and program consistence in case of error */
	/* better error reporting is needed */
	RSB_ERROR(RSB_ERRM_ES);
	rsb__destroy_coo_matrix_t(&coo);
err:
	return blas_invalid_handle /*RSB_BLAS_INVALID_VAL*/;
ok:
	rsb__destroy_coo_matrix_t(&coo);
	return A; 
}

/* @endcond */
