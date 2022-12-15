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
 * This source file contains functions for sparse matrices sum.
 */
/* FIXME : UNFINISHED, UNCHECKED, UNSECURED, preliminar code
 * TODO : spscatter
 * */

#include "rsb_internals.h"

#define RSB_SPSUM_VERBOSITY 0

struct rsb_mtx_t * rsb__do_matrix_sum(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_err_t * errvalp)
{
	/*!
	 * \todo: unfinished
	 * TODO: overflows are possible; need checks.
	 * TODO: need a specialized approach for symmetric matrices (e.g.: sum of two symmetric ones is symmetric)!
	 * TODO: what about RSB_NUMERICAL_TYPE_SAME_TYPE ?
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_SPSUM_VERBOSITY
	rsb_nnz_idx_t rnz=0,an,bn,cn;
#endif /* RSB_SPSUM_VERBOSITY */
	struct rsb_coo_mtx_t cooa,coob,cooc;
	struct rsb_mtx_t * mtxCp = NULL;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS|RSB_FLAG_DISCARD_ZEROS|RSB_FLAG_SORTED_INPUT;
	rsb_coo_idx_t tam,tak,tbm,tbk;

	RSB_BZERO_P(&cooa);
	RSB_BZERO_P(&coob);
	RSB_BZERO_P(&cooc);

	if( !mtxAp /*|| !alphap || !betap*/ || !mtxBp )
	{
		/* Note: alphap==NULL and betap==NULL are allowed */
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	tam = RSB_MTX_TRANSPOSED_ROWS(mtxAp,transA);
	tak = RSB_MTX_TRANSPOSED_COLS(mtxAp,transA);
	tbm = RSB_MTX_TRANSPOSED_ROWS(mtxBp,transB);
	tbk = RSB_MTX_TRANSPOSED_COLS(mtxBp,transB);

	if( ( tam != tbm ) || ( tak != tbk ) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(!mtxAp) { errval = RSB_ERR_GENERIC_ERROR; RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(!mtxBp) { errval = RSB_ERR_GENERIC_ERROR; RSB_PERR_GOTO(err,RSB_ERRM_ES);}

	if(!RSB_IS_VALID_NNZ_SUM(mtxAp->nnz,mtxBp->nnz))
	{
		errval = RSB_ERR_LIMITS;
		RSB_PERR_GOTO(err,"number of matrices sum nnz may exceed maximum allowed.\n");
	}

	/*
	 * TODO: if same type, same transposition, same matrices, one may simply clone and scale.
	 * */
	cooc.nnz = 2*RSB_MAX(mtxAp->nnz+mtxBp->nnz,tam+1)+2*(tam+1); /* FIXME: this is excess allocation for symmetry handling */
	cooc.typecode = typecode;
	if(rsb__callocate_coo_matrix_t(&cooc)!=&cooc) { RSB_PERR_GOTO(err,RSB_ERRM_ES); }

	/* ...->flags&RSB_FLAG_ANY_SYMMETRY is added filtering to avoid expansion */
	RSB_DO_ERROR_CUMULATE(errval,rsb__clone_coo(mtxAp,transA,alphap,typecode,&cooa,flags|(mtxAp->flags&(RSB_FLAG_ANY_SYMMETRY))));
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES); }

	RSB_DO_ERROR_CUMULATE(errval,rsb__clone_coo(mtxBp,transB, betap,typecode,&coob,flags|(mtxBp->flags&(RSB_FLAG_ANY_SYMMETRY))));
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES); }

	if(RSB_DOES_TRANSPOSE(transA))
	{
		errval = rsb__util_sort_row_major_inner(cooa.VA,cooa.IA,cooa.JA,cooa.nnz,cooa.nr,cooa.nc,typecode,RSB_FLAG_WANT_ROW_MAJOR_ORDER);
		if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES); }
	}

	if(RSB_DOES_TRANSPOSE(transB))
	{
		errval = rsb__util_sort_row_major_inner(coob.VA,coob.IA,coob.JA,coob.nnz,coob.nr,coob.nc,typecode,RSB_FLAG_WANT_ROW_MAJOR_ORDER);
		if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES); }
	}

	RSB_COO_MEMCPY(cooc.VA,cooc.IA,cooc.JA,cooa.VA,cooa.IA,cooa.JA,0       ,0,cooa.nnz,RSB_SIZEOF(typecode));
	RSB_COO_MEMCPY(cooc.VA,cooc.IA,cooc.JA,coob.VA,coob.IA,coob.JA,cooa.nnz,0,coob.nnz,RSB_SIZEOF(typecode));

	cooc.nnz=cooa.nnz+coob.nnz;
	RSB_DO_FLAG_DEL(flags,RSB_FLAG_SORTED_INPUT);

	{
		rsb_nnz_idx_t dnz = 0;
		errval = rsb__cor_merge_dups(typecode, cooc.VA, cooc.IA, cooc.JA, 0, cooa.nnz, coob.nnz, 0, 1, &dnz, NULL);
		cooc.nnz-=dnz;
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT);
	}

#if RSB_SPSUM_VERBOSITY
	RSB_STDOUT("sum output will have %d nnz\n",rnz);
#endif /* RSB_SPSUM_VERBOSITY */
	cooc.nr = mtxAp->nr;
	cooc.nc = mtxAp->nc;
	mtxCp = rsb__do_mtx_alloc_from_coo_inplace(cooc.VA,cooc.IA,cooc.JA,cooc.nnz,cooc.typecode,cooc.nr,cooc.nc,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags|RSB_FLAG_DUPLICATES_SUM,&errval);
	if(!mtxCp||RSB_SOME_ERROR(errval))
	{
	       	RSB_PERR_GOTO(err,RSB_ERRM_ES);
       	}
	RSB_DO_FLAG_DEL(mtxCp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
	RSB_DO_FLAG_DEL(mtxCp->flags,RSB_FLAG_DUPLICATES_SUM);

#if RSB_SPSUM_VERBOSITY
	RSB_STDOUT("sum output will have %d nnz\n",rnz);
#endif /* RSB_SPSUM_VERBOSITY */
	goto ok;
err:
	rsb__do_perror(NULL,errval);
	RSB_ERROR("!\n");
	rsb__destroy_coo_matrix_t(&cooc);
ok:
	rsb__destroy_coo_matrix_t(&cooa);
	rsb__destroy_coo_matrix_t(&coob);
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	RSB_DO_MTX_RETURN(mtxCp,errval);
}
/* @endcond */
