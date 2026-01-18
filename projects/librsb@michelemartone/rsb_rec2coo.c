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
 * @brief Code for matrix format conversion. 
 * */
#include "rsb_common.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

static rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo_leaf(struct rsb_mtx_t * mtxAp, rsb_bool_t do_shift)
{
	/**
		\ingroup gr_internals
		TODO: move somewhere else
		TODO: need temporary memory to pass to e.g. rsb__do_switch_compressed_array_to_fullword_coo() and thus avoid allocations.

	// to free the unnecessary data:
	// RSB_CONDITIONAL_FREE(mtxAp
	// RSB_CONDITIONAL_FREE(mtxAp->all_leaf_matrices)
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(rsb__is_coo_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
			rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*)mtxAp->bpntr,mtxAp->nnz,do_shift?mtxAp->roff:0),
			rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*)mtxAp->bindx,mtxAp->nnz,do_shift?mtxAp->coff:0);
		goto err;
	}
	errval = rsb__do_switch_compressed_array_to_fullword_coo(mtxAp->bpntr,mtxAp->Mdim,do_shift?mtxAp->roff:0,NULL);
	mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCOR;
	RSB_DO_FLAG_DEL(mtxAp->flags,(RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS));
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
		rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*)mtxAp->bindx,mtxAp->nnz,do_shift?mtxAp->coff:0);
err:
	RSB_DO_FLAG_SUBST(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES),(RSB_FLAG_WANT_COO_STORAGE));
	// rsb__do_print_matrix_stats(mtxAp, RSB_CONST_DUMP_MATRIX_MARKET RSB_CONST_DUMP_RECURSION_BRIEF, NULL);
	// rsb__do_print_matrix_stats(mtxAp, RSB_CONST_DUMP_MATRIX_MARKET RSB_CONST_DUMP_RECURSION_BRIEF, NULL);
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo(struct rsb_mtx_t * mtxAp, rsb_bool_t do_shift)
{
	/**
		\ingroup gr_internals
		TODO: move somewhere else
		TODO: flags checks
		FIXME: UNTESTED

	// to free the unnecessary data:
	// RSB_CONDITIONAL_FREE(mtxAp)
	// RSB_CONDITIONAL_FREE(mtxAp->all_leaf_matrices)
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(RSB_UNLIKELY(!mtxAp))
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_BADARGS;
	}

	if(RSB_UNLIKELY(rsb__is_root_matrix(mtxAp)))
	{
		if(!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS))
			errval = RSB_ERR_BADARGS;
		else
			errval = rsb__do_switch_recursive_matrix_to_fullword_storage(mtxAp);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
	}

	/*deleted by rsb__do_switch_recursive_matrix_to_fullword_storage*/
	//RSB_DO_FLAG_ADD(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
			if(submatrix)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo(submatrix,do_shift));
	}
	else
		errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo_leaf(mtxAp,do_shift);
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo_parallel(struct rsb_mtx_t * mtxAp, rsb_bool_t do_shift)
{
	/**
		\ingroup gr_internals
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_submatrix_idx_t all_leaf_matrices_n = mtxAp->all_leaf_matrices_n;
	rsb_submatrix_idx_t n;

	//rsb__do_print_matrix_stats(mtxAp, RSB_CONST_DUMP_RECURSION_BRIEF, NULL);
	#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC
	for(n=0;n<all_leaf_matrices_n;++n)
	{
		struct rsb_mtx_t *submatrix = mtxAp->all_leaf_matrices[n].mtxlp;
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo_leaf(submatrix,do_shift));
	}
	//rsb__do_print_matrix_stats(mtxAp, RSB_CONST_DUMP_MATRIX_MARKET , NULL);
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_sorted(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop)
{
	/**
		\ingroup gr_internals
		 Makes sense only for in place allocated.
		On exit, the pointer matrix is deallocated.
		FIXME: Here it would make sense to use a recursive merge algorithm.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *fsm = NULL;
	//rsb_flags_t flags;
	struct rsb_coo_mtx_t coo;
	int wmb = 1; /* want merge based (new: 20140727) */

	RSB_BZERO_P(&coo);

	if(RSB_UNLIKELY(!mtxAp))
	{
		RSB_ERROR(RSB_ERRM_E_MTXAP);
		return RSB_ERR_BADARGS;
	}
	
	if(mtxAp->all_leaf_matrices_n == 1)
		wmb = 0; /* merge routine will not convert a single leaf's format */

#if 0
	fsm = rsb__do_get_first_submatrix(mtxAp);
	if(!fsm)
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	flags = mtxAp->flags;
	errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo(mtxAp,RSB_BOOL_TRUE);
	RSB_CONDITIONAL_FREE(mtxAp->all_leaf_matrices);
	RSB_INIT_COO_FROM_MTX(coop,mtxAp);
	RSB_BIND_COO_TO_MTX(coop,fsm);
	RSB_CONDITIONAL_FREE(mtxAp);
	//if((errval = rsb__util_sort_row_major_parallel(coop->VA,coop->IA,coop->JA,coop->nnz,coop->nr,coop->nc,coop->typecode,flags))!=RSB_ERR_NO_ERROR)
	if((errval = rsb__util_sort_row_major_bucket_based_parallel(coop->VA,coop->IA,coop->JA,coop->nnz,coop->nr,coop->nc,coop->typecode,flags))!=RSB_ERR_NO_ERROR)
		goto err;
#else
	if(wmb)
	{
		errval = rsb__leaves_merge_multiple(mtxAp, NULL, NULL, NULL, 0, 1);

		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err, RSB_ERRM_ES);
		}

		fsm = rsb__do_get_first_submatrix(mtxAp);
		if(!fsm)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err, RSB_ERRM_ES);
		}
	
		RSB_INIT_COO_FROM_MTX(coop, mtxAp);
		RSB_BIND_COO_TO_MTX(coop, fsm);
		RSB_ASSERT(coop->VA || coop->nnz == 0);
		RSB_ASSERT(coop->IA || coop->nnz == 0);
		RSB_ASSERT(coop->JA || coop->nnz == 0);
	}
	else
	{
		fsm = rsb__do_get_first_submatrix(mtxAp);
		if(!fsm)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err, RSB_ERRM_ES);
		}
	
		RSB_INIT_CXX_FROM_MTX(&coo, mtxAp);
		coo.nr = coo.nc = 0; // to have rsb__allocate_coo_matrix_t allocate nnz and not more
		if(rsb__allocate_coo_matrix_t(&coo) != &coo)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err, RSB_ERRM_ES);
		}
		errval = rsb__do_get_rows_sparse(RSB_TRANSPOSITION_N, NULL, mtxAp, coo.VA, coo.IA, coo.JA, 0, mtxAp->nr-1, &coo.nnz, RSB_FLAG_NOFLAGS);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err, RSB_ERRM_ES);
		}
		RSB_INIT_COO_FROM_MTX(coop, mtxAp);
		RSB_BIND_COO_TO_MTX(coop, fsm);
		RSB_COO_MEMCPY_parallel(coop->VA, coop->IA, coop->JA, coo.VA, coo.IA, coo.JA, 0, 0, coo.nnz, mtxAp->el_size);
		rsb__destroy_coo_matrix_t(&coo);
	}
	fsm->VA = NULL;
	fsm->bpntr = NULL;
	fsm->bindx = NULL;
	rsb__destroy_inner(mtxAp);
#endif
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_unsorted(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop)
{
	/**
		\ingroup gr_internals
		TODO: move somewhere else
		FIXME: UNTESTED,TEMPORARY, makes sense only for in place allocated
		this conversion gives you sorted coordinates.
		on exit, the pointer matrix is deallocated
		FIXME: error behaviour is undefined
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//struct rsb_coo_mtx_t coo;
	struct rsb_mtx_t *fsm = NULL;

	if(RSB_UNLIKELY(!mtxAp))
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_BADARGS;
	}
#if 0
	RSB_INIT_CXX_FROM_MTX(&coo,mtxAp);
	if(rsb__allocate_coo_matrix_t(&coo)!=&coo)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	errval = rsb__do_get_rows_sparse(RSB_TRANSPOSITION_N,NULL,mtxAp,coo.VA,coo.IA,coo.JA,0,mtxAp->nr-1,&coo.nnz,RSB_FLAG_NOFLAGS);
	if(RSB_SOME_ERROR(errval)) goto err;
	//rsb__destroy_inner(mtxAp);
	rsb__do_mtx_free(mtxAp);
	coop->VA = coo.VA;
	coop->IA = coo.IA;
	coop->JA = coo.JA;
	RSB_INIT_COO_FROM_MTX(coop,&coo);
//	mtxAp->VA = coo.VA;
//	mtxAp->bpntr = coo.IA;
//	mtxAp->bindx = coo.JA;
#else
	fsm = rsb__do_get_first_submatrix(mtxAp);
	if(!fsm)
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo_parallel(mtxAp,RSB_BOOL_TRUE);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_BIND_COO_TO_MTX(coop,fsm);
	RSB_INIT_COO_FROM_MTX(coop,mtxAp);
	fsm->VA = NULL;
	fsm->bpntr = NULL;
	fsm->bindx = NULL;
	rsb__destroy_inner(mtxAp);
#endif
err:
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
