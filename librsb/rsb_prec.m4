dnl
dnl
dnl	@author: Michele Martone
dnl
ifelse(LIBMMVBR_INCLUDED_PREC_M4,1,`',`
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
/* @cond INNERDOC */
dnl
/**
 * @file
 * @brief
 * Auxiliary functions.
 */
RSB_M4_HEADER_MESSAGE()dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_PREC_H_INCLUDED
#define RSB_PREC_H_INCLUDED
')
dnl
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

dnl
#include "rsb_common.h"
dnl #include "rsb_internals.h"
dnl #include "rsb_types.h"
dnl 
dnl
dnl
dnl	FIXME : COMMENT THIS FILE
dnl	-------------------------
dnl
dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
',`dnl
dnl
dnl `rsb_err_t rsb_do_csr_ilu0_'touppercase(RSB_M4_CHOPSPACES(mtype))`(struct rsb_mtx_t * mtxAp)'dnl
`static rsb_err_t rsb_do_csr_ilu0_'touppercase(RSB_M4_CHOPSPACES(mtype))`(struct rsb_coo_mtx_t * coop)'dnl
dnl `rsb_err_t rsb_do_csr_ilu0_'touppercase(RSB_M4_CHOPSPACES(mtype))(mtype `*VA, const rsb_coo_idx_t *PA, const rsb_coo_idx_t *JA)'dnl
{
	/**
	 * \ingroup gr_internals
		FIXME: INCOMPLETE, EXPERIMENTAL, TEMPORARILY HERE
		On exit, the matrix will contain the L and U factors of a pattern preserving incomplete LU factorization (ILU 0).
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t i;

{
	mtype *VA = coop->VA;
	const rsb_coo_idx_t *PA = coop->IA;
	const rsb_coo_idx_t *JA = coop->JA;
dnl	const rsb_coo_idx_t *PA = mtxAp->bpntr;
dnl	const rsb_coo_idx_t *JA = mtxAp->bindx;
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
dnl err:
	RSB_DO_ERR_RETURN(errval)
}
dnl
dnl
')dnl
')dnl
dnl

dnl const void * rsb__prec_ilu0(struct rsb_mtx_t * mtxAp)`'dnl
rsb_err_t rsb__prec_ilu0(struct rsb_mtx_t * mtxAp)`'dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`dnl
{
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
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( mtxAp->typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
		return rsb_do_csr_ilu0_`'touppercase(RSB_M4_CHOPSPACES(mtype))(&coo);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	errval = RSB_ERR_INTERNAL_ERROR;
err:
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl

rsb_err_t rsb__prec_csr_ilu0(struct rsb_coo_mtx_t * coop)`'dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`dnl
{
	// FIXME: temporary
	if(coop->nr==1)
		goto err;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( coop->typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
		dnl return rsb_do_csr_ilu0_`'touppercase(RSB_M4_CHOPSPACES(mtype))(coo->VA,coo->IA,coo->JA);
		return rsb_do_csr_ilu0_`'touppercase(RSB_M4_CHOPSPACES(mtype))(coop);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
err:
	return RSB_ERR_INTERNAL_ERROR;
}
')dnl
dnl

dnl
#ifdef __cplusplus
}
#endif /* __cplusplus */
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_PREC_H_INCLUDED */
')
')
dnl
/* @endcond */
dnl
