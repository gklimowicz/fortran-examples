
/* @cond INNERDOC */
/**
 * @file
 * @brief
 * Auxiliary functions.
 */

/*

Copyright (C) 2008-2022 Michele Martone

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
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

 */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb_common.h"
static rsb_err_t rsb_do_csr_ilu0_DOUBLE(struct rsb_coo_mtx_t * coop){
	/**
	 * \ingroup gr_internals
		FIXME: INCOMPLETE, EXPERIMENTAL, TEMPORARILY HERE
		On exit, the matrix will contain the L and U factors of a pattern preserving incomplete LU factorization (ILU 0).
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t i;

{
	double *VA = coop->VA;
	const rsb_coo_idx_t *PA = coop->IA;
	const rsb_coo_idx_t *JA = coop->JA;
	for(i=1;i<coop->nr;++i)
	{
		const rsb_nnz_idx_t ifp = PA[i],ilp = PA[i+1],irnz = ilp-ifp;
		rsb_nnz_idx_t idp = RSB_MARKER_NNZ_VALUE,ikp = RSB_MARKER_NNZ_VALUE;
		if(irnz)
		{

			idp = rsb__nnz_split_coo_bsearch(JA+ifp,i,irnz)+ifp;
			assert(idp<=ilp);
			assert(idp>=ifp);
			for(ikp=ifp;ikp<idp;++ikp)// k = 1...i-1
			{
				/* FIXME: write a sparse vectors dot product macro and apply it here */
				const rsb_nnz_idx_t k = JA[ikp],kfp = PA[k],klp = PA[k+1],krnz = klp-kfp;
				const int kdp = rsb__nnz_split_coo_bsearch(JA+kfp,k,krnz)+kfp;
				rsb_nnz_idx_t kjp = kfp,ijp = ikp+1;
				VA[ikp]/=VA[kdp];
				/* FIXME: to optimize this phase, we should loop on the shorter row */
				for(;ijp<ilp;++ijp)// j = k+1...n
				{
					for(;JA[kjp]<JA[ijp] && kjp<klp;++kjp)
						;
					if(kjp==klp)
						goto out;
					/* JA[kjp]>=JA[ijp] */
					for(;JA[kjp]>JA[ijp] && ijp<ilp;++ijp)
						;
					if(ijp==ilp)
						goto out;
					/* JA[kjp]==JA[ijp] */
					VA[ijp]-=VA[ikp]*VA[kjp];
				}
out:
				RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
				
			}
		}
	}
}
	RSB_DO_ERR_RETURN(errval)
}
static rsb_err_t rsb_do_csr_ilu0_FLOAT(struct rsb_coo_mtx_t * coop){
	/**
	 * \ingroup gr_internals
		FIXME: INCOMPLETE, EXPERIMENTAL, TEMPORARILY HERE
		On exit, the matrix will contain the L and U factors of a pattern preserving incomplete LU factorization (ILU 0).
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t i;

{
	float *VA = coop->VA;
	const rsb_coo_idx_t *PA = coop->IA;
	const rsb_coo_idx_t *JA = coop->JA;
	for(i=1;i<coop->nr;++i)
	{
		const rsb_nnz_idx_t ifp = PA[i],ilp = PA[i+1],irnz = ilp-ifp;
		rsb_nnz_idx_t idp = RSB_MARKER_NNZ_VALUE,ikp = RSB_MARKER_NNZ_VALUE;
		if(irnz)
		{

			idp = rsb__nnz_split_coo_bsearch(JA+ifp,i,irnz)+ifp;
			assert(idp<=ilp);
			assert(idp>=ifp);
			for(ikp=ifp;ikp<idp;++ikp)// k = 1...i-1
			{
				/* FIXME: write a sparse vectors dot product macro and apply it here */
				const rsb_nnz_idx_t k = JA[ikp],kfp = PA[k],klp = PA[k+1],krnz = klp-kfp;
				const int kdp = rsb__nnz_split_coo_bsearch(JA+kfp,k,krnz)+kfp;
				rsb_nnz_idx_t kjp = kfp,ijp = ikp+1;
				VA[ikp]/=VA[kdp];
				/* FIXME: to optimize this phase, we should loop on the shorter row */
				for(;ijp<ilp;++ijp)// j = k+1...n
				{
					for(;JA[kjp]<JA[ijp] && kjp<klp;++kjp)
						;
					if(kjp==klp)
						goto out;
					/* JA[kjp]>=JA[ijp] */
					for(;JA[kjp]>JA[ijp] && ijp<ilp;++ijp)
						;
					if(ijp==ilp)
						goto out;
					/* JA[kjp]==JA[ijp] */
					VA[ijp]-=VA[ikp]*VA[kjp];
				}
out:
				RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
				
			}
		}
	}
}
	RSB_DO_ERR_RETURN(errval)
}
static rsb_err_t rsb_do_csr_ilu0_FLOAT_COMPLEX(struct rsb_coo_mtx_t * coop){
	/**
	 * \ingroup gr_internals
		FIXME: INCOMPLETE, EXPERIMENTAL, TEMPORARILY HERE
		On exit, the matrix will contain the L and U factors of a pattern preserving incomplete LU factorization (ILU 0).
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t i;

{
	float complex *VA = coop->VA;
	const rsb_coo_idx_t *PA = coop->IA;
	const rsb_coo_idx_t *JA = coop->JA;
	for(i=1;i<coop->nr;++i)
	{
		const rsb_nnz_idx_t ifp = PA[i],ilp = PA[i+1],irnz = ilp-ifp;
		rsb_nnz_idx_t idp = RSB_MARKER_NNZ_VALUE,ikp = RSB_MARKER_NNZ_VALUE;
		if(irnz)
		{

			idp = rsb__nnz_split_coo_bsearch(JA+ifp,i,irnz)+ifp;
			assert(idp<=ilp);
			assert(idp>=ifp);
			for(ikp=ifp;ikp<idp;++ikp)// k = 1...i-1
			{
				/* FIXME: write a sparse vectors dot product macro and apply it here */
				const rsb_nnz_idx_t k = JA[ikp],kfp = PA[k],klp = PA[k+1],krnz = klp-kfp;
				const int kdp = rsb__nnz_split_coo_bsearch(JA+kfp,k,krnz)+kfp;
				rsb_nnz_idx_t kjp = kfp,ijp = ikp+1;
				VA[ikp]/=VA[kdp];
				/* FIXME: to optimize this phase, we should loop on the shorter row */
				for(;ijp<ilp;++ijp)// j = k+1...n
				{
					for(;JA[kjp]<JA[ijp] && kjp<klp;++kjp)
						;
					if(kjp==klp)
						goto out;
					/* JA[kjp]>=JA[ijp] */
					for(;JA[kjp]>JA[ijp] && ijp<ilp;++ijp)
						;
					if(ijp==ilp)
						goto out;
					/* JA[kjp]==JA[ijp] */
					VA[ijp]-=VA[ikp]*VA[kjp];
				}
out:
				RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
				
			}
		}
	}
}
	RSB_DO_ERR_RETURN(errval)
}
static rsb_err_t rsb_do_csr_ilu0_DOUBLE_COMPLEX(struct rsb_coo_mtx_t * coop){
	/**
	 * \ingroup gr_internals
		FIXME: INCOMPLETE, EXPERIMENTAL, TEMPORARILY HERE
		On exit, the matrix will contain the L and U factors of a pattern preserving incomplete LU factorization (ILU 0).
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t i;

{
	double complex *VA = coop->VA;
	const rsb_coo_idx_t *PA = coop->IA;
	const rsb_coo_idx_t *JA = coop->JA;
	for(i=1;i<coop->nr;++i)
	{
		const rsb_nnz_idx_t ifp = PA[i],ilp = PA[i+1],irnz = ilp-ifp;
		rsb_nnz_idx_t idp = RSB_MARKER_NNZ_VALUE,ikp = RSB_MARKER_NNZ_VALUE;
		if(irnz)
		{

			idp = rsb__nnz_split_coo_bsearch(JA+ifp,i,irnz)+ifp;
			assert(idp<=ilp);
			assert(idp>=ifp);
			for(ikp=ifp;ikp<idp;++ikp)// k = 1...i-1
			{
				/* FIXME: write a sparse vectors dot product macro and apply it here */
				const rsb_nnz_idx_t k = JA[ikp],kfp = PA[k],klp = PA[k+1],krnz = klp-kfp;
				const int kdp = rsb__nnz_split_coo_bsearch(JA+kfp,k,krnz)+kfp;
				rsb_nnz_idx_t kjp = kfp,ijp = ikp+1;
				VA[ikp]/=VA[kdp];
				/* FIXME: to optimize this phase, we should loop on the shorter row */
				for(;ijp<ilp;++ijp)// j = k+1...n
				{
					for(;JA[kjp]<JA[ijp] && kjp<klp;++kjp)
						;
					if(kjp==klp)
						goto out;
					/* JA[kjp]>=JA[ijp] */
					for(;JA[kjp]>JA[ijp] && ijp<ilp;++ijp)
						;
					if(ijp==ilp)
						goto out;
					/* JA[kjp]==JA[ijp] */
					VA[ijp]-=VA[ikp]*VA[kjp];
				}
out:
				RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
				
			}
		}
	}
}
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__prec_ilu0(struct rsb_mtx_t * mtxAp){
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t coo;

	if(!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_ES);
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	if(     !rsb__is_terminal_recursive_matrix(mtxAp) ||
		!rsb__is_css_matrix(mtxAp) ||
		(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES) ||
		rsb__is_symmetric(mtxAp) ||
		!rsb__is_square(mtxAp) ||
		rsb__submatrices(mtxAp)!=1 ||
 		RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT)
		)
	{
		if(!rsb__is_terminal_recursive_matrix(mtxAp))
			RSB_ERROR(RSB_ERRM_NT);
		if(!rsb__is_css_matrix(mtxAp))
			RSB_ERROR("not a CSR matrix!\n");
		/*if(rsb__is_root_matrix(mtxAp)!=RSB_BOOL_TRUE)
			RSB_ERROR("non-root matrix!\n");*/
		if( rsb__submatrices(mtxAp)!=1)
			RSB_ERROR("non-single submatrix (%d)!\n",rsb__submatrices(mtxAp));
		if( rsb__is_symmetric(mtxAp))
			RSB_ERROR("matrix is symmetric!\n");
		if(!rsb__is_square(mtxAp))
			RSB_ERROR(RSB_ERRM_NON_SQUARE);
		RSB_ERROR(RSB_ERRM_ES);
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	if(mtxAp->nr==1)
		goto err;
	if((errval = rsb__project_rsb_to_coo(mtxAp,&coo))!=RSB_ERR_NO_ERROR)
		goto err;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
		return rsb_do_csr_ilu0_DOUBLE(&coo);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT  )
		return rsb_do_csr_ilu0_FLOAT(&coo);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
		return rsb_do_csr_ilu0_FLOAT_COMPLEX(&coo);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
		return rsb_do_csr_ilu0_DOUBLE_COMPLEX(&coo);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	errval = RSB_ERR_INTERNAL_ERROR;
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__prec_csr_ilu0(struct rsb_coo_mtx_t * coop){
	// FIXME: temporary
	if(coop->nr==1)
		goto err;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( coop->typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
				return rsb_do_csr_ilu0_DOUBLE(coop);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( coop->typecode == RSB_NUMERICAL_TYPE_FLOAT  )
				return rsb_do_csr_ilu0_FLOAT(coop);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( coop->typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
				return rsb_do_csr_ilu0_FLOAT_COMPLEX(coop);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( coop->typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
				return rsb_do_csr_ilu0_DOUBLE_COMPLEX(coop);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
err:
	return RSB_ERR_INTERNAL_ERROR;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */


/* @endcond */
