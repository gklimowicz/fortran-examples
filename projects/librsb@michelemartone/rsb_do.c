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
 * Implementation of the interface functions.
 * \internal
 *
 * */
/* TODO: introduce RSB_MSG_BADARGS_ERROR(ERRVAL,MSG,BAN) */
#include "rsb_common.h"
#include "rsb_util.h"
#include "rsb.h"
#include "rsb_unroll.h"
#ifdef RSB_HAVE_SYS_UTSNAME_H 
#include <sys/utsname.h>	/* uname */
#endif /* RSB_HAVE_SYS_UTSNAME_H */
#include "rsb_do.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

static rsb_err_t rsb__do_prec_build(struct rsb_mtx_t ** mtxLpp, struct rsb_mtx_t ** mtxUpp, const struct rsb_mtx_t * mtxAp)
{
	/* 
	 * FIXME: UNFINISHED, UNTESTED
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *L=NULL,*U=NULL;
	struct rsb_coo_mtx_t csr,lcoo,ucoo;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_MATRIX_FLAGS;

	csr.VA=NULL;
	csr.IA=NULL;
	csr.JA=NULL;

	RSB_DEBUG_ASSERT(mtxAp  );
	RSB_DEBUG_ASSERT(mtxLpp );
	RSB_DEBUG_ASSERT(mtxUpp);

	if(mtxAp->nr != mtxAp->nc)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	csr.nr=mtxAp->nr;
	csr.nc=mtxAp->nc;
	csr.nnz = RSB_MAX(mtxAp->nnz,RSB_MAX(csr.nr,csr.nc)+1);
	csr.typecode=mtxAp->typecode;
	ucoo=csr;
	lcoo=csr;
	if(rsb__allocate_coo_matrix_t(&csr)!=&csr)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	csr.nnz=mtxAp->nnz;
	/* a more efficient routine would: */
	/* build a fullword CSR clone (using a RSB-to-CSR constructor) */
	/* perform ILU on the CSR struct */
	/* build RSB clones for the L and U CSR parts (using a modified CSR-to-RSB constructor; possibly building them together, in one shot) */
	/* FIXME: TODO */
	/* a better 'quick hack' solution would: */
	/* build a fullword CSR clone (using a RSB-to-CSR constructor) */
	if((errval = rsb__do_get_csr(mtxAp->typecode,mtxAp,csr.VA,csr.IA,csr.JA,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS))!=RSB_ERR_NO_ERROR)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	/* perform ILU on the CSR struct */
	if((errval = rsb__prec_csr_ilu0(&csr))!=RSB_ERR_NO_ERROR)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	/* perform CSR to COO conversion */
	rsb__do_count_tri_in_csr(&csr,&lcoo.nnz,&ucoo.nnz);

	ucoo.nnz=RSB_MAX(ucoo.nnz,ucoo.nr+1);
	if(rsb__allocate_coo_matrix_t(&ucoo)!=&ucoo) { errval = RSB_ERR_ENOMEM; goto err; }

	lcoo.nnz=RSB_MAX(lcoo.nnz,lcoo.nr+1);

	if(rsb__allocate_coo_matrix_t(&lcoo)!=&lcoo) { errval = RSB_ERR_ENOMEM; goto err; }
	rsb__do_count_tri_in_csr(&csr,&lcoo.nnz,&ucoo.nnz);

	rsb__do_copy_tri_from_csr_to_coo(&csr,&lcoo,&ucoo);

	/* allocating L and U */
	L = rsb__do_mtx_alloc_from_coo_const(lcoo.VA,lcoo.IA,lcoo.JA,lcoo.nnz,lcoo.typecode,lcoo.nr,lcoo.nc,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags|RSB_FLAG_LOWER_TRIANGULAR,&errval);
	if(!L)
	{
		rsb__destroy_coo_matrix_t(&lcoo);
		RSB_BZERO_P(&lcoo);
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	U = rsb__do_mtx_alloc_from_coo_const(ucoo.VA,ucoo.IA,ucoo.JA,ucoo.nnz,ucoo.typecode,ucoo.nr,ucoo.nc,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags|RSB_FLAG_UPPER_TRIANGULAR,&errval);
	if(!U)
	{
		rsb__destroy_coo_matrix_t(&ucoo);
		RSB_BZERO_P(&ucoo);
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	/*RSB_ERROR(RSB_ERRM_WUF);*/

	*mtxLpp=L;
	*mtxUpp=U;
	goto ret;
err:
	RSB_MTX_FREE(L);
	RSB_MTX_FREE(U);
	rsb__do_perror(NULL,errval);
	errval = RSB_ERR_BADARGS;
ret:
	rsb__destroy_coo_matrix_t(&lcoo);
	rsb__destroy_coo_matrix_t(&ucoo); 
	rsb__destroy_coo_matrix_t(&csr);
	return errval;
}

rsb_err_t rsb__do_get_preconditioner(void *opd, const struct rsb_mtx_t * mtxAp, rsb_precf_t prec_flags, const void *ipd)/* FIXME: temporary interface */
{
	// FIXME: UNFINISHED, UNTESTED
	rsb_err_t errval = RSB_ERR_GENERIC_ERROR;
	struct rsb_mtx_t * LU[2]={NULL,NULL};

	if(!opd || !mtxAp)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	errval = rsb__do_prec_build(&LU[0],&LU[1],mtxAp);
	rsb__memcpy(opd,LU,sizeof(LU));
err:
	return errval;
}
#if 0
rsb_err_t rsb_do_prec_apply(const struct rsb_prec_t * precp, const void *r)
{
	/*
	 * given the preconditioner P, computes \f$ r \leftarrow {P}^{-1} r.\f$
	 * \f$ r' = {P}^{-1} r \f$
	 * \f$ r' = {L \cdot U}^{-1} r \f$
	 * \f$ r' = {U}^{-1} {L}^{-1} r \f$
	 * \f$ r' = ({U}^{-1} ({L}^{-1} r)) \f$
	 * */
	//rsb__debug_print_vector(r,precp->L->nr,precp->L->typecode,1);
	double one=1.0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_spsv(RSB_TRANSPOSITION_N,&one,precp->L,r,1,r,1));
	if(RSB_SOME_ERROR(errval))
		goto err;
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_spsv(RSB_TRANSPOSITION_N,&one,precp->U,r,1,r,1));
	if(RSB_SOME_ERROR(errval))
		goto err;

	//rsb__debug_print_vector(r,precp->L->nr,precp->L->typecode,1);
err:
	if(RSB_SOME_ERROR(errval))
		rsb__do_perror(NULL,errval);
	return errval;
}
#endif

rsb_err_t rsb__do_get_rows_sparse(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags)
{
	/* TODO: having a return scaled rows would be an efficient feature. */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t off = 0;

#if RSB_ALLOW_ZERO_DIM 
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
	{
		goto err; /* FIXME: skipping checks on ldB, ldC, op_flags*/
	}
#endif

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}

	if(VA == NULL && JA == NULL && IA == NULL && rnzp != NULL)
	{
		*rnzp = rsb__dodo_get_rows_nnz(mtxAp, frA, lrA,flags,&errval);
		goto err;
	}

	if(!rnzp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"user did not supply a results nonzeroes pointer\n");
	}

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
	{
		off=1;
		lrA--;
		frA--;
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
	}

	if(frA<0 || lrA>mtxAp->nr)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}

	if(!VA)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_NULL_VA);
	}

	if(frA>lrA)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}

	if(RSB_DOES_TRANSPOSE(transA))
		RSB_SWAP(rsb_coo_idx_t*,IA,JA);

	*rnzp=0;
#if 0
        errval = rsb__do_get_rows_sparse_rec(mtxAp,VA,frA,lrA,IA,JA,rnzp,off,off);
	if(flags & RSB_FLAG_SORT_INPUT)
		errval = rsb__util_sort_row_major_inner(VA,IA,JA,*rnzp,mtxAp->nr+off,mtxAp->nc+off,mtxAp->typecode,flags|RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR);
#if 0
	if( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
		rsb__util_nnz_array_to_fortran_indices(IA,*rnzp),
		rsb__util_nnz_array_to_fortran_indices(JA,*rnzp);
#endif
#else
#if 1
	{
		rsb_coo_idx_t i;
		for(i=frA;i<=lrA;++i)
		        RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_rows_sparse_rec(mtxAp,VA,i,i,IA,JA,rnzp,off,off));
	}
#else
	{
		/* this is too slow for many leaf matrices */
		rsb_submatrix_idx_t si;
		rsb_coo_idx_t i;
		for(i=frA;i<=lrA;++i)
		for(si=0;si<mtxAp->all_leaf_matrices_n;++si)
		        RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_rows_sparse_rec(mtxAp->all_leaf_matrices[si].mtxlp,VA,i,i,IA,JA,rnzp,off,off));
	}
#endif
#endif
	RSB_DEBUG_ASSERT(rsb__dodo_get_rows_nnz(mtxAp,frA,lrA,flags,NULL)==*rnzp);
	if(RSB_DOES_TRANSPOSE(transA))
	{
		RSB_SWAP(rsb_coo_idx_t*,IA,JA);
		/* swapping back and ready for sorting. if now, would get column major. */
		if(RSB_SOME_ERROR(errval = rsb__util_sort_row_major_inner(VA,IA,JA,*rnzp,mtxAp->nr,mtxAp->nc,mtxAp->typecode,flags)))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_EM);
		}
	}

	if(RSB_DOES_CONJUGATE(transA))
	{
		RSB_DO_ERROR_CUMULATE(errval,rsb__util_do_conjugate(VA,mtxAp->typecode,*rnzp));
	}

	if(alphap)
	{
		RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xscal(mtxAp->typecode,*rnzp,alphap,VA,1));
	}
err:
	if(RSB_SOME_ERROR(errval))
	RSB_ERROR(RSB_ERRM_NL);
	return errval;
}

rsb_err_t rsb__do_scal(struct rsb_mtx_t * mtxAp, const void * d, rsb_trans_t trans)
{
	/* TODO : what should be the semantics of scaling a symmetric matrix ? */
	/* FIXME : and error handling ? **/
	rsb_err_t errval = RSB_ERR_UNSUPPORTED_OPERATION;

#ifdef RSB_HAVE_OPTYPE_SCALE

#if RSB_ALLOW_ZERO_DIM 
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
	{
		errval = RSB_ERR_NO_ERROR;
		goto err; /* FIXME: skipping checks on ldB, ldC, op_flags*/
	}
#endif

	if(!mtxAp)
	//	return RSB_ERR_NO_ERROR;
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}

	if( rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix=NULL;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			rsb_coo_idx_t off;
			if(RSB_DOES_NOT_TRANSPOSE(trans))
				off=submatrix->roff-mtxAp->roff;
			else
				off=submatrix->coff-mtxAp->coff;
	
			rsb__do_scal(submatrix,((rsb_byte_t*)d)+mtxAp->el_size*off,trans);
		}
		//return RSB_ERR_NO_ERROR;
		{
			errval = RSB_ERR_NO_ERROR;
			goto err;
		}
	}
	else
		errval = rsb__do_scale(mtxAp,trans,d);
#else /* RSB_HAVE_OPTYPE_SCALE */
#endif /* RSB_HAVE_OPTYPE_SCALE */
err:
	return errval;
}

rsb_err_t rsb__dodo_getdiag( const struct rsb_mtx_t * mtxAp, void * diagonal )
{
  	// FIXME : missing documentation and error checks!
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}
	if(!RSB_BLOCK_CROSSED_BY_DIAGONAL(0,0,mtxAp->nr,mtxAp->nc))
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	if(1)
	//if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_COO )
	{
		// FIXME: THIS IS SLOW, TEMPORARY
		rsb_coo_idx_t i;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		long nt = rsb_get_num_threads();
		const int gdc = RSB_DIVIDE_IN_CHUNKS(mtxAp->nr,nt);
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		#pragma omp parallel for schedule(static,gdc) reduction(|:errval)  RSB_NTC
		for(i=0;i<mtxAp->nr;++i)
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_coo_element(mtxAp,((rsb_char_t*)diagonal)+mtxAp->el_size*(i),i,i));
		errval = RSB_ERR_NO_ERROR;
	}
#if 0
	else
	{
		RSB_BZERO(diagonal,mtxAp->el_size*RSB_MTX_DIAG_SIZE(mtxAp));
		errval = rsb_do_getdiag(mtxAp,diagonal);
	}
#endif
err:
	return errval;
}

static rsb_err_t rsb_do_elemental_scale(struct rsb_mtx_t * mtxAp, const void * alphap)
{
	/*!
	   \ingroup rsb_doc_matrix_handling

	   Computes \f$ A \leftarrow \alpha A \f$.

	   \param \rsb_mtxt_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \return \rsberrcodemsg
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * submatrix=NULL;

	if(!mtxAp)
       	{
		errval = RSB_ERR_BADARGS;
	       	goto err;
	}
#if RSB_WANT_PARALLEL_ELEMENTAL_OPS
	if(1)
		errval = rsb__do_elemental_scale_parallel(mtxAp,alphap);
#else /* RSB_WANT_PARALLEL_ELEMENTAL_OPS */
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t smi;
		//#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC
		RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
			RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xscal(submatrix->typecode,submatrix->nnz,alphap,submatrix->VA,1));
	}
#endif /* RSB_WANT_PARALLEL_ELEMENTAL_OPS */
	else
		RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xscal(mtxAp->typecode,mtxAp->nnz,alphap,mtxAp->VA,1));
err:
	return errval;
}

static rsb_err_t rsb_do_elemental_scale_inv(struct rsb_mtx_t * mtxAp, const void * alphap)
{
	/*!
	   \ingroup rsb_doc_matrix_handling

	   Computes \f$ A \leftarrow \frac{1}{\alpha} A \f$.

	   \param \rsb_mtxt_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \return \rsberrcodemsg
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * submatrix=NULL;

	if(!mtxAp)
       	{
		errval = RSB_ERR_BADARGS;
	       	goto err;
	}

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t smi;
		//#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC
		RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
			RSB_DO_ERROR_CUMULATE(errval,rsb__vector_scale_inv(submatrix->VA,alphap,submatrix->typecode,submatrix->nnz));
	}
	else
		RSB_DO_ERROR_CUMULATE(errval,rsb__vector_scale_inv(mtxAp->VA,alphap,mtxAp->typecode,mtxAp->nnz));
err:
	return errval;
}

static rsb_err_t rsb_do_elemental_pow(struct rsb_mtx_t * mtxAp, const void * alphap)
{
	/*!
	   \ingroup rsb_doc_matrix_handling

	   Raises each matrix element to the given power.

	   \param \rsb_mtxt_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \return \rsberrcodemsg
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * submatrix=NULL;

	if(!mtxAp)
       	{
		errval = RSB_ERR_BADARGS;
	       	goto err;
	}
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t smi;
		//#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC
		RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
			RSB_DO_ERROR_CUMULATE(errval,rsb__util_vector_pow(submatrix->VA,submatrix->typecode,alphap,submatrix->nnz));
	}
	else
		RSB_DO_ERROR_CUMULATE(errval,rsb__util_vector_pow(mtxAp->VA,mtxAp->typecode,alphap,mtxAp->nnz));
err:
	return errval;
}

static rsb_err_t rsb_dodo_negation(struct rsb_mtx_t * mtxAp)
{
#ifdef RSB_HAVE_OPTYPE_NEGATION
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp)
		{errval = RSB_ERR_BADARGS;goto err;}

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix=NULL;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb_dodo_negation(submatrix));
	}
	else
		errval = rsb_do_negation(mtxAp,0xf1c57415,RSB_DEFAULT_TRANSPOSITION);
#else /* RSB_HAVE_OPTYPE_NEGATION */
	/* FIXME : eliminate negation as mop ! */
//	return RSB_ERR_UNSUPPORTED_OPERATION;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if(!mtxAp)
		{errval = RSB_ERR_BADARGS;goto err;}

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix=NULL;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb_dodo_negation(submatrix));
	}
	else
		/* FIXME : assuming elements are contiguous ! */
		errval = rsb__util_do_negate(mtxAp->VA,mtxAp->typecode,mtxAp->element_count);
#endif /* RSB_HAVE_OPTYPE_NEGATION */
err:
	return errval;
}

rsb_err_t rsb__do_elemental_unop(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags)
{
	// FIXME: untested
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp) {errval = RSB_ERR_BADARGS; goto err;}
	switch(elop_flags)
	{
		case(RSB_ELOPF_NEG):
		errval = rsb_dodo_negation(mtxAp);
		break;
		//#define RSB_ELOPF_SQRT		0x00000010		/*!< Elemental square root (usable with rsb_mtx_elemental_unop). */
		/*
		case(RSB_ELOPF_SQRT):
		errval=....(mtxAp);
		break;
		*/
	/*	case(RSB_ELOPF_TRANS):
		errval = rsb_transpose(&mtxAp);
		break;
		case(RSB_ELOPF_HTRANS):
		errval = rsb_htranspose(&mtxAp);
		break;*/
		default:
		{errval = RSB_ERR_BADARGS; goto err;}
	}
err:
	return errval;
}

rsb_err_t rsb__do_elemental_binop(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags, const void * opp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_trans_t transA = RSB_TRANSPOSITION_N; 
	void * topp=NULL;

	if(!mtxAp) {errval = RSB_ERR_BADARGS; goto err;}
	if(!opp) {errval = RSB_ERR_BADARGS; goto err;}

	switch(elop_flags)
	{
		case(RSB_ELOPF_SCALE_COLS_REAL):
		case(RSB_ELOPF_SCALE_COLS):
		transA = RSB_TRANSPOSITION_T; 
		break;
		default:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
	}

	switch(elop_flags)
	{
		case(RSB_ELOPF_SCALE_COLS_REAL):
		case(RSB_ELOPF_SCALE_ROWS_REAL):
		if( RSB_IS_MATRIX_TYPE_COMPLEX(mtxAp->typecode) )
		{
			/* Note: this is inefficient */
			const rsb_type_t typecode = RSB_NUMERICAL_TYPE_REAL_TYPE(mtxAp->typecode);
			if( NULL == (topp = rsb__calloc_vector(mtxAp->nr,mtxAp->typecode)) )
			{
				errval = RSB_ERR_ENOMEM;
				goto err;
			}
			errval = rsb__cblas_Xcopy(typecode,mtxAp->nr,opp,1,topp,2);
		}		
		case(RSB_ELOPF_SCALE_COLS):
		case(RSB_ELOPF_SCALE_ROWS):
		errval = rsb__do_scal(mtxAp,topp?topp:opp,transA);
		break;
		case(RSB_ELOPF_MUL):
		errval = rsb_do_elemental_scale(mtxAp,opp);
		break;
		case(RSB_ELOPF_DIV):
		errval = rsb_do_elemental_scale_inv(mtxAp,opp);
		break;
		case(RSB_ELOPF_POW):
		errval = rsb_do_elemental_pow(mtxAp,opp);
		break;
		default:
		{
			errval = RSB_ERR_BADARGS;
		       	goto err;
		}
	}
err:
	RSB_CONDITIONAL_FREE(topp);
	return errval;
}

rsb_nnz_idx_t rsb__dodo_get_rows_nnz(const struct rsb_mtx_t *mtxAp, rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_flags_t flags, rsb_err_t * errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t rnz = 0;

	RSB_DEBUG_ASSERT(fr <= lr);

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		lr--,fr--;
	errval = rsb__do_get_rows_nnz(mtxAp,fr,lr,&rnz);
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	if(RSB_SOME_ERROR(errval))
		rnz=0;
	return rnz;
}

#if RSB_WANT_PARALLEL_ELEMENTAL_OPS
rsb_err_t rsb__do_elemental_scale_parallel(struct rsb_mtx_t * mtxAp, const void * alphap);
{
	/**
		\ingroup gr_internals
		TODO: move to somewhere else
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_thread_t wet = rsb_get_num_threads();

	if(RSB_UNLIKELY(mtxAp->nnz<wet*RSB_MIN_THREAD_MEMCPY_NNZ))
	{
		RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xscal(mtxAp->typecode,mtxAp->nnz,alphap,((rsb_byte_t*)mtxAp->VA),1));
	}
	else
	{
		rsb_nnz_idx_t wi;
		size_t cnz=(mtxAp->nnz+wet-1)/wet;	/* chunk size */
		#pragma omp parallel for schedule(static,1) RSB_NTC
		for(wi=0;wi<wet;++wi)
		{
			size_t coff=wi*cnz;
			size_t cnnz=(wi<wet-1)?cnz:mtxAp->nnz-((wet-1)*cnz);
			printf("%d nz on %d\n",cnnz,wi);
			RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xscal(mtxAp->typecode,cnnz,alphap,((rsb_byte_t*)mtxAp->VA)+RSB_SIZEOF(mtxAp->typecode)*mtxAp->coff,1));
		}
	}
}
#endif /* RSB_WANT_PARALLEL_ELEMENTAL_OPS */

rsb_err_t rsb__do_matrix_add_to_dense(const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_nnz_idx_t ldB, rsb_nnz_idx_t nr, rsb_nnz_idx_t nc, rsb_bool_t rowmajor, void * Bp)
{
	//  FIXME: could this be documented in two groups (mops and unfinished) at the same time ?
	//  TODO: what about supporting transA ?
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_aligned_t pone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];

	if(! (nr>0 && nc>0 && mtxAp) )
	{
	       	errval = RSB_ERR_BADARGS;
	       	goto err;
	}

	if(!alphap)
	{
		rsb__util_set_area_to_converted_integer(&pone[0],mtxAp->typecode,+1);
		alphap = &pone[0];
	}

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{

		rsb_submatrix_idx_t smi;
		#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC if(mtxAp->nnz > RSB_MIN_NNZ_FOR_PARALLEL_ADD_TO_DENSE)
		RSB_SUBMATRIX_FOREACH_LEAF_IDX(mtxAp,smi)
		{
			const struct rsb_mtx_t * submatrix = (mtxAp)->all_leaf_matrices[smi].mtxlp;
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_add_submatrix_to_dense(submatrix,alphap,Bp,ldB,nr,nc,rowmajor));
		}
	}
	else
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_add_submatrix_to_dense(mtxAp,alphap,Bp,ldB,nr,nc,rowmajor));
err:
	return errval;
}

rsb_err_t rsb__do_switch_rsb_mtx_to_csr_sorted(struct rsb_mtx_t * mtxAp, void ** VAP, rsb_coo_idx_t ** IAP, rsb_coo_idx_t ** JAP, rsb_flags_t flags)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t coo;
	const rsb_nnz_idx_t nnz=mtxAp?mtxAp->nnz:0;
	const rsb_coo_idx_t m=mtxAp?mtxAp->nr:0;
	//const rsb_coo_idx_t k=mtxAp?mtxAp->nc:0;
	const rsb_flags_t mflags=mtxAp?mtxAp->flags:RSB_FLAG_NOFLAGS;

	if(!mtxAp)
       	{
	       	errval = RSB_ERR_BADARGS;
	       	goto err;
       	}

	if(!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{
	       	errval = RSB_ERR_BADARGS;
	       	goto err;
       	}

	if(!IAP || !JAP || !VAP)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_VIJP);
	}
	errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_sorted(mtxAp,&coo);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb__util_compress_to_row_pointers_array(NULL,nnz,m,mflags,flags,coo.IA);
	if(RSB_SOME_ERROR(errval))
		goto err;
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		rsb__util_nnz_array_to_fortran_indices(coo.JA,nnz);
	*JAP=coo.JA;
	*IAP=coo.IA;
	*VAP=coo.VA;
err:
	return errval;
}

rsb_err_t rsb__do_get_csr(rsb_type_t typecode, const struct rsb_mtx_t *mtxAp, rsb_byte_t * VA, rsb_nnz_idx_t * RP, rsb_coo_idx_t * JA, rsb_flags_t flags)
{
	/* NOTE this writes more than mtxAp->nnz elements! */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp)
	{
	       	errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
       	}

	if((mtxAp->typecode != typecode && typecode != RSB_NUMERICAL_TYPE_SAME_TYPE ) || (mtxAp->flags != flags))/* FIXME: condition on cflags is unnecessarily restrictive */
	{
		const rsb_flags_t flagsC = flags | (mtxAp->flags & RSB_FLAG_ALL_STRUCTURAL_FLAGS);
		struct rsb_mtx_t * mtxCp = NULL;
		errval = rsb__mtx_clone(&mtxCp, typecode, RSB_TRANSPOSITION_N, NULL, mtxAp, flagsC);
		if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		errval = rsb__dodo_get_csr(mtxCp,&VA,&RP,&JA);
		RSB_MTX_FREE(mtxCp);
	}
	else
		errval = rsb__dodo_get_csr(mtxAp,&VA,&RP,&JA);

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
       	}
	/* Note: would be nice to move C -> Fortran indices semantics to rsb__dodo_get_csr */
	if(flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE)
		rsb__util_nnz_array_to_fortran_indices(RP,mtxAp->nr+1),
		rsb__util_nnz_array_to_fortran_indices(JA,mtxAp->nnz);
err:
	return errval;
}

rsb_err_t rsb__do_get_matrix_info(const struct rsb_mtx_t *mtxAp, enum rsb_mif_t miflags, void* info, size_t buflen)
{
	/*!
	   \ingroup FIXME 
	   \warning \rsb_warn_unfinished_msg
		FIXME: UNFINISHED, UNTESTED
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_real_t rrv=0;
	size_t szv=0;
	rsb_coo_idx_t civ=0;
	rsb_nnz_idx_t niv=0;
	rsb_flags_t fiv=0;
	rsb_type_t tiv=0;
	rsb_blk_idx_t biv=0;
	char*cis=(char*)info;

	if(!mtxAp || !info)
	{
	       	errval = RSB_ERR_BADARGS;
	       	goto err;
	}
	switch(miflags)
	{
		case RSB_MIF_INDEX_STORAGE_IN_BYTES__TO__SIZE_T:
		{
		       	szv = rsb__get_index_storage_amount(mtxAp);
		       	if(buflen<=0) *(size_t*)info = szv;
		        else snprintf(cis,buflen,"%zd",szv);
	       	}
		break;
		case RSB_MIF_INDEX_STORAGE_IN_BYTES_PER_NNZ__TO__RSB_REAL_T:
		{
			const size_t isa = rsb__get_index_storage_amount(mtxAp);
			if ( mtxAp->nnz > 0 )
				rrv = ((rsb_real_t)isa) / ((rsb_real_t)mtxAp->nnz);
			else
				rrv = 2 * sizeof(rsb_coo_idx_t);
			if(buflen<=0)
				*(rsb_real_t*) info=rrv;
			else
				snprintf(cis,buflen,"%lg",rrv);
		}
		break;
		case RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T:
		{
		       	civ = (mtxAp->nr);
	                if(buflen<=0) *(rsb_coo_idx_t*)info = civ;
		        else snprintf(cis,buflen,"%ld",(long int)civ);
	       	}
		break;
		case RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T:
		{
		       	civ = (mtxAp->nc);
	                if(buflen<=0) *(rsb_coo_idx_t*)info = civ;
		        else snprintf(cis,buflen,"%ld",(long int)civ);
	       	}
		break;
		case RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T:
		{
		       	niv = (mtxAp->nnz);
	                if(buflen<=0) *(rsb_nnz_idx_t*)info = niv;
		        else snprintf(cis,buflen,"%ld",(long int)niv);
	       	}
		break;
		case RSB_MIF_TOTAL_SIZE__TO__SIZE_T:
		{
		       	szv = rsb__get_sizeof(mtxAp);
		       	if(buflen<=0) *(size_t*)info = szv;
		        else snprintf(cis,buflen,"%zd",szv);
	       	}
		break;
		case RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T:
		{
		       	fiv = (mtxAp->flags);
		       	if(buflen<=0) *(rsb_flags_t*)info = fiv;
		        else snprintf(cis,buflen,"%ld",(long int)fiv);
	       	}
		break;
		case RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T:
		{
		       	tiv = (mtxAp->typecode);
		       	if(buflen<=0) *(rsb_type_t*)info = tiv;
		        else snprintf(cis,buflen,"%ld",(long int)tiv);
	       	}
		break;
		case RSB_MIF_MATRIX_INFO__TO__CHAR_P:				
		{
		       	tiv = (mtxAp->typecode);
		       	if(buflen<=0) { errval = RSB_ERR_BADARGS; goto err; }
		        else snprintf(cis,buflen,RSB_PRINTF_MTX_SUMMARY_ARGS(mtxAp));

#if 0
			if(cis && *cis) /* anonymize pointer address (no address dump) */
			{
				char * lbp = strstr(cis,"[0x");
				if(lbp)
					for(lbp+=3;*lbp&&*lbp!=']';++lbp)
						*lbp='F';
			}
#endif
	       	}
		break;
		case RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T:				
		{
		       	biv = (mtxAp->all_leaf_matrices_n);
		       	if(buflen<=0) *(rsb_blk_idx_t*)info = biv;
		        else snprintf(cis,buflen,"%d",biv);
	       	}
		break;
		default:
		errval = RSB_ERR_GENERIC_ERROR;
	}
err:
	return errval;
}

rsb_err_t rsb__do_check_leak(void)
{
	/*!
	   \ingroup rsb_doc_library
	  
	   Called after \ref rsb_lib_exit(), will report on the standard output stream
	   (see #RSB_IO_WANT_OUTPUT_STREAM) whether some previously allocated
	   memory area was not freed by librsb.
	   \n
	   Will report leak information only if built with the #RSB_DISABLE_ALLOCATOR_WRAPPER symbol undefined.
	   \n
	   Will return #RSB_ERR_NO_ERROR on no leak; an error otherwise.
	   \n
	  
	   \warning \rsb_warn_soon_to_be_deprecated_msg 
	   \return \rsberrcodemsg
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(rsb__get_g_rsb_memory_count())
	{
		RSB_INFO("WARNING: allocated memory  : %zu : POSSIBLE MEMORY LEAK\n",rsb__get_g_rsb_memory_count());
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
	}

	if(rsb__get_g_rsb_allocations_count())
	{
		RSB_INFO("WARNING: allocations count : %zu : POSSIBLE MEMORY LEAK\n",rsb__get_g_rsb_allocations_count());
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
	}
	return errval;
}

rsb_err_t rsb__do_matrix_norm(const struct rsb_mtx_t * mtxAp , void * np, enum rsb_extff_t flags)
{
	rsb_err_t errval = RSB_ERR_BADARGS;

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(ret,RSB_ERRM_ES);
	}

	switch(flags)
	{
		case RSB_EXTF_NORM_ONE:
		errval = rsb__do_infinity_norm(mtxAp,np,RSB_BOOL_FALSE,RSB_TRANSPOSITION_T);
		break;
		case RSB_EXTF_NORM_TWO:
		errval = rsb__cblas_Xnrm2(mtxAp->typecode,mtxAp->nnz,rsb__do_get_first_submatrix(mtxAp)->VA,1,np);
		break;
		case RSB_EXTF_NORM_INF:
		errval = rsb__do_infinity_norm(mtxAp,np,RSB_BOOL_FALSE,RSB_TRANSPOSITION_N);
		break;
		default:
		RSB_PERR_GOTO(ret,RSB_ERRM_ES);
		break;
	}
ret:
	return errval;
}

rsb_err_t rsb__do_matrix_compute(const struct rsb_mtx_t * mtxAp , void * dp, enum rsb_extff_t flags)
{
	rsb_err_t errval = RSB_ERR_BADARGS;

#if RSB_ALLOW_ZERO_DIM 
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
	{
		errval = RSB_ERR_NO_ERROR;
		goto ret; /* FIXME: skipping further error checks */
	}
#endif

	if( mtxAp == NULL || ( mtxAp->nnz > 0 && dp == NULL ) )
		goto ret;

	switch(flags)
	{
		case RSB_EXTF_SUMS_ROW:
		errval = rsb__do_rows_sums_inner(mtxAp,dp,RSB_BOOL_FALSE,RSB_TRANSPOSITION_N);
		break;
		case RSB_EXTF_SUMS_COL:
		errval = rsb__do_rows_sums_inner(mtxAp,dp,RSB_BOOL_FALSE,RSB_TRANSPOSITION_T);
		break;
		case RSB_EXTF_ASUMS_ROW:
		errval = rsb__do_absolute_rows_sums(mtxAp,dp);
		break;
		case RSB_EXTF_ASUMS_COL:
		errval = rsb__do_absolute_columns_sums(mtxAp,dp);
		break;
		case RSB_EXTF_DIAG:
		errval = rsb__dodo_getdiag(mtxAp,dp);
		break;
		default:
		RSB_PERR_GOTO(ret,RSB_ERRM_ES);
		break;
	}
ret:
	return errval;
}

rsb_err_t rsb__do_load_vector_file_as_matrix_market(const rsb_char_t * filename, rsb_type_t typecode, void * yp, rsb_coo_idx_t *yvlp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!filename || ((!yvlp) && (!yp)))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}
	if(yvlp)
	{
		/* FIXME: temporarily ignoring second dimension! */
		rsb_coo_idx_t yvk=0, yvm=0;
		rsb_bool_t is_vector = RSB_BOOL_FALSE;

		if(RSB_SOME_ERROR(errval = rsb__util_mm_info_matrix_f(filename,&yvm,&yvk,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&is_vector)) )
		{
			RSB_PERR_GOTO(err,RSB_ERRM_EM);
		}
		*yvlp=yvm;
	}
	if(yp)
	{
		/* printf("stub: reading in %s...\n",filename); */
		rsb_nnz_idx_t vnz=0;

		errval = rsb__util_mm_load_vector_f(filename,&yp,&vnz,typecode);
	}
err:
	return errval;
}

struct rsb_mtx_t * rsb__dodo_load_matrix_file_as_matrix_market(const rsb_char_t * filename, rsb_flags_t flags, rsb_type_t typecode, rsb_err_t *errvalp)
{
	// FIXME
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

       	errval = rsb__do_load_matrix_file_as_matrix_market(&mtxAp,filename,flags,typecode);
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

rsb_bool_t rsb__do_was_initialized(void)
{
	/*!
	   \ingroup rsb_doc_library
	  
	   Call this function to know whether the library had already been initialized or not.
	   \n
	   This function is mainly intended to be used in between \ref rsb_lib_exit() and \ref rsb_lib_init() calls,
	   or generally after one or more calls to \ref rsb_lib_init() were already been done.
	   \n
	   It is not meant to be called before the 'first' initialization ever, unless 
	   the user is sure this library was built on a system which supports default
	   initialization to zero of static variables (which indeed is supported by most standards;
	   e.g.: ANSI C: http://flash-gordon.me.uk/ansi.c.txt ).  
	  
	   \return #RSB_BOOL_TRUE if it was initialized, #RSB_BOOL_FALSE otherwise.
	 */
	/* TODO: redocument elsewhere! redundant function! */
	return (rsb_global_session_handle.rsb_g_initialized == RSB_BOOL_TRUE) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;
}

static rsb_err_t rsb__do_switch_rsb_mtx_to_coo_unsorted(struct rsb_mtx_t * mtxAp, void ** VAP, rsb_coo_idx_t ** IAP, rsb_coo_idx_t ** JAP, rsb_flags_t flags)
{
	struct rsb_coo_mtx_t coo;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_ASSERT( RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS) );

	errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_unsorted(mtxAp,&coo);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		rsb__util_nnz_array_to_fortran_indices(coo.IA,coo.nnz),
		rsb__util_nnz_array_to_fortran_indices(coo.JA,coo.nnz);
	*JAP = coo.JA;
	*IAP = coo.IA;
	*VAP = coo.VA;

	return errval;
}

static rsb_err_t rsb__do_switch_rsb_mtx_to_coo_sorted(struct rsb_mtx_t * mtxAp, void ** VAP, rsb_coo_idx_t ** IAP, rsb_coo_idx_t ** JAP, rsb_flags_t flags)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t coo;
	const rsb_nnz_idx_t nnz = mtxAp ? mtxAp-> nnz:0;

	RSB_ASSERT( RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS) );

	RSB_BZERO_P(&coo);
	errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_sorted(mtxAp, &coo);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err, RSB_ERRM_EM);

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		rsb__util_nnz_array_to_fortran_indices(coo.IA, nnz),
		rsb__util_nnz_array_to_fortran_indices(coo.JA, nnz);

	*JAP = coo.JA;
	*IAP = coo.IA;
	*VAP = coo.VA;
err:
	return errval;
}

rsb_err_t rsb__do_switch_rsb_mtx_to_coo(struct rsb_mtx_t * mtxAp, void ** VAP, rsb_coo_idx_t ** IAP, rsb_coo_idx_t ** JAP, rsb_flags_t flags)
{
	rsb_err_t errval = RSB_ERR_BADARGS;

	if(!mtxAp)
		RSB_PERR_GOTO(err, RSB_ERRM_E_MTXAP"\n");

	/* Purpose of the following is avoidance of internally allocated memory leakage. */
	/* TODO: As an improvement, one may relax this constraint when the allocation wrapper is off. */
	if(!RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
		RSB_PERR_GOTO(err, RSB_ERRM_IMNIP);

	if(!IAP || !JAP || !VAP)
		RSB_PERR_GOTO(err, RSB_ERRM_EM);

	if(RSB_DO_FLAG_HAS(flags, RSB_FLAG_SORTED_INPUT))
		errval = rsb__do_switch_rsb_mtx_to_coo_sorted(mtxAp, VAP, IAP, JAP, flags);
	else
		errval = rsb__do_switch_rsb_mtx_to_coo_unsorted(mtxAp, VAP, IAP, JAP, flags);
err:
	return errval;
}

#if RSB_WANT_COO_BEGIN 
struct rsb_mtx_t * rsb__do_mtx_alloc_from_coo_begin(rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_flags_t flags, rsb_err_t * errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	blas_sparse_matrix bmtxA = RSB_BLAS_INVALID_VAL;

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC) && RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_BFSAH);
	}

	rsb__init_struct(mtxAp = rsb__calloc(sizeof(struct rsb_mtx_t)));
	if(!mtxAp)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP"\n");
	}
	RSB_MTX_SET_HBDF(mtxAp);
	bmtxA = mtxAp->RSB_MTX_BDF = rsb__BLAS_Xuscr_begin(nrA,ncA,typecode);
	if( mtxAp->RSB_MTX_BDF == RSB_BLAS_INVALID_VAL )
	{
		errval = RSB_ERR_GENERIC_ERROR;
		RSB_CONDITIONAL_FREE(mtxAp);
		RSB_PERR_GOTO(err,RSB_ERRM_IPEWIEM);
	}

	/* FIXME : the following need an improvement  */
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE)) rsb__BLAS_ussp( bmtxA, blas_one_base);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT)) rsb__BLAS_ussp( bmtxA, blas_unit_diag );
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_TRIANGULAR)) rsb__BLAS_ussp( bmtxA, blas_lower_triangular);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_TRIANGULAR)) rsb__BLAS_ussp( bmtxA, blas_upper_triangular);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_TRIANGULAR)) rsb__BLAS_ussp( bmtxA, blas_triangular); /* ask for detection */
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_SYMMETRIC)) rsb__BLAS_ussp( bmtxA, blas_lower_symmetric);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_SYMMETRIC)) rsb__BLAS_ussp( bmtxA, blas_upper_symmetric);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_HERMITIAN)) rsb__BLAS_ussp( bmtxA, blas_lower_hermitian);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_HERMITIAN)) rsb__BLAS_ussp( bmtxA, blas_upper_hermitian);
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

rsb_err_t rsb__do_mtx_alloc_from_coo_end(struct rsb_mtx_t ** mtxApp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	blas_sparse_matrix bmtxA = RSB_BLAS_INVALID_VAL;
	struct rsb_mtx_t * mtxBp = NULL;
	struct rsb_mtx_t * mtxAp = NULL;

	if(!mtxApp || !*mtxApp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAPP);
	}

	mtxAp = *mtxApp ;

	if( !RSB_MTX_HBDF( mtxAp ) )
	{
		/* errval = RSB_ERR_NO_ERROR; */
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_DNSAMIWAFCB);
	}

	bmtxA = RSB_MTX_HBDFH(mtxAp);
	/* FIXME: missing serious check on mtxAp->flags ! */
	if( rsb__BLAS_Xuscr_end_flagged(bmtxA,NULL) == RSB_BLAS_INVALID_VAL )
	{
	       	errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_PFTM);
		/* FIXME: insufficient cleanup */
	}
	mtxBp = rsb__BLAS_inner_matrix_retrieve(bmtxA);
	*mtxApp = mtxBp;
	rsb__free(mtxAp);
	rsb__BLAS_handle_free(bmtxA); /* ignoring return value ... */
err:
	return errval;
}
#endif /* RSB_WANT_COO_BEGIN */

rsb_err_t rsb__do_upd_vals(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags, const void * omegap)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	switch(elop_flags)
	{
		case(RSB_ELOPF_MUL):
		case(RSB_ELOPF_DIV):
		case(RSB_ELOPF_POW):
		case(RSB_ELOPF_SCALE_ROWS):
		case(RSB_ELOPF_SCALE_COLS):
		case(RSB_ELOPF_SCALE_ROWS_REAL):
		case(RSB_ELOPF_SCALE_COLS_REAL):
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_elemental_binop(mtxAp,elop_flags,omegap));
		break;
		case(RSB_ELOPF_NEG):
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_elemental_unop(mtxAp,elop_flags));
		break;
		default: {errval = RSB_ERR_BADARGS; goto err;}
	}
err:
	return errval;
}

rsb_err_t rsb__do_mtx_get_info(const struct rsb_mtx_t *mtxAp, enum rsb_mif_t miflags, void* minfop)
{
	rsb_err_t errval = RSB_ERR_UNIMPLEMENTED_YET;
	errval = rsb__do_get_matrix_info(mtxAp,miflags,minfop,0);
	return errval;
}

rsb_err_t rsb__do_file_mtx_save(const struct rsb_mtx_t * mtxAp, const rsb_char_t * filename)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	errval = rsb__do_print_matrix_stats(mtxAp,RSB_CONST_DUMP_MATRIX_MARKET,filename);
	return errval;
}

rsb_err_t rsb__do_vec_save(const rsb_char_t * filename, rsb_type_t typecode, const void * Yp, rsb_coo_idx_t yvl)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	FILE*stream = NULL;
	const int incX = 1;
	// TODO: improve error handling here.

       	if(filename == NULL)
		stream = stdout;
	else
		stream = fopen(filename,"w");

	errval = rsb__debug_print_vector_extra(Yp,yvl,typecode,incX,0x1,stream);

	if(filename != NULL && stream)
		fclose(stream);

	return errval;
}

struct rsb_mtx_t * rsb__do_mtx_alloc_from_csr_inplace (void *VA, rsb_coo_idx_t * RP, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp )
{
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if RSB_ALLOW_EMPTY_MATRICES
	if( nnzA > 0 )
#endif /* RSB_ALLOW_EMPTY_MATRICES */
		errval = rsb__util_uncompress_row_pointers_array(RP,nrA,flagsA,flagsA,RP);
	/* now RP is effectively IA, nnzA elements long */
	if(RSB_SOME_ERROR(errval))
		RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	if( RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		rsb_util_coo_arrays_sub(RP,JA,1,1,nnzA),
		RSB_DO_FLAG_DEL(flagsA,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
	RSB_DO_FLAG_ADD(flagsA,RSB_FLAG_SORTED_INPUT);
	if(errval == RSB_ERR_NO_ERROR)
		mtxAp = rsb__do_mtx_alloc_from_coo_inplace(VA,RP,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,&errval);
	return mtxAp;
}

rsb_err_t rsb__do_file_mtx_rndr(void * pmp, const char * filename, rsb_coo_idx_t pmlWidth, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( pmlWidth != pmHeight )
	{
		/* enforce square rendering limitation of rsb_file_mtx_rndr (FIXME) */
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	switch(rflags)
	{
		case(RSB_MARF_RGB):
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_pixmap_RGB_from_matrix(filename,pmp,pmWidth,pmHeight));
		break;
		case(RSB_MARF_EPS):
		// rsb_dump_postscript_from_mtx_t(fd,mtxAp,1,1,pmWidth,pmHeight,0);
		// RSB_DO_ERROR_CUMULATE(errval,rsb__dump_postscript_from_matrix(filename,1,1,pmWidth,pmHeight,0));
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_UNIMPLEMENTED_YET);
		//RSB_DO_ERROR_CUMULATE(errval,rsb__dump_postscript_recursion_from_matrix(filename,1,1,pmWidth,pmHeight,RSB_FLAG_NOFLAGS,1,1,0,RSB_NUMERICAL_TYPE_DEFAULT));
		break;
		default: {errval = RSB_ERR_UNIMPLEMENTED_YET; goto err;}
	}
err:
	return errval;
}


/* @endcond */
