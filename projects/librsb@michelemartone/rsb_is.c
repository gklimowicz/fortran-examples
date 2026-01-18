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
 * This source file contains matrix info getter functions.
 * */

#include "rsb_internals.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

rsb_bool_t rsb__is_coo_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */
	rsb_bool_t is;
	RSB_DEBUG_ASSERT(mtxAp);

	is = (
#ifdef RSB_MATRIX_STORAGE_BCOR
	 (mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR)||
#endif /* RSB_MATRIX_STORAGE_BCOR */
#ifdef RSB_MATRIX_STORAGE_BCOC
	 (mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOC)||
#endif /* RSB_MATRIX_STORAGE_BCOC */
	  0
	)?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
	return is;
}

rsb_bool_t rsb__is_square(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */
	RSB_DEBUG_ASSERT(mtxAp);

	return (mtxAp->nr == mtxAp->nc)?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
}

rsb_bool_t rsb__is_hermitian(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */

	return (rsb__get_hermitian_flag(mtxAp))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__is_triangle(rsb_flags_t flags)
{
	return (rsb__is_lower_triangle(flags) | rsb__is_upper_triangle(flags));
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/

rsb_bool_t rsb__is_lower_triangle(rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * */
	return (RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER|RSB_FLAG_TRIANGULAR));
}

rsb_bool_t rsb__is_upper_triangle(rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * */
	return (RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER|RSB_FLAG_TRIANGULAR));
}

rsb_bool_t rsb__is_symmetric(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */
	return (rsb__get_symmetry_flag(mtxAp))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
}

rsb_bool_t rsb__is_not_unsymmetric(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */
	if(rsb__get_hermitian_flag(mtxAp) || rsb__get_symmetry_flag(mtxAp))
		return RSB_BOOL_TRUE;
	else
		return RSB_BOOL_FALSE;
}

rsb_bool_t rsb__is_csr_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 */
	const rsb_bool_t is_csr = RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_BCSS_STORAGE);

	return is_csr;
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__is_bcss_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given matrix is row or column block major
	 *
	 * */
	rsb_bool_t ret = 0;

	if(!mtxAp)
		return ret;
	ret = 
#ifdef RSB_MATRIX_STORAGE_BCSR
		mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSR ||
#endif /* RSB_MATRIX_STORAGE_BCSR */
#ifdef RSB_MATRIX_STORAGE_BCSC
		mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSC ||
#endif /* RSB_MATRIX_STORAGE_BCSC */
		 0;

#if RSB_EXPERIMENTAL_USE_PURE_BCSS
	if(ret)
	{
		RSB_ASSERT(mtxAp->br>0);
		RSB_ASSERT(mtxAp->bc>0);
	}
#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
	return ret;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_bool_t rsb__is_css_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given matrix is not CSR or CSC.
	 *
	 * */
	//rsb_bool_t ret = RSB_BOOL_FALSE;
	rsb_blk_idx_t br, bc;

	rsb__get_blocking_size(mtxAp, &br, &bc);

	return ( br==1 && bc==1 ) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;
}

rsb_bool_t rsb__is_bcsr_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given matrix is row block major.
	 *
	 * */
	rsb_bool_t ret = RSB_BOOL_FALSE;

#ifdef RSB_MATRIX_STORAGE_BCSR
	if( ( mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSR ) != 0 ) ret = RSB_BOOL_TRUE;
#endif /* RSB_MATRIX_STORAGE_BCSR */

#if RSB_EXPERIMENTAL_USE_PURE_BCSS
	if(ret)
	{
		RSB_ASSERT(mtxAp->br>0);
		RSB_ASSERT(mtxAp->bc>0);
	}
#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
	return ret;
}

#if RSB_WANT_BCSC_LEAVES  
rsb_bool_t rsb__is_bcsc_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given matrix is column block major.
	 *
	 * */
	rsb_bool_t ret = 0;

	if(!mtxAp)
		return ret;
	
	ret = 
#ifdef RSB_MATRIX_STORAGE_BCSC
		mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSC ||
#endif /* RSB_MATRIX_STORAGE_BCSC */
	0;

#if RSB_EXPERIMENTAL_USE_PURE_BCSS
	if(ret)
	{
		RSB_ASSERT(mtxAp->br>0);
		RSB_ASSERT(mtxAp->bc>0);
	}
#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
	return ret;
}
#endif /* RSB_WANT_BCSC_LEAVES */

rsb_bool_t rsb__have_fixed_blocks_matrix_flags(rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given flags are for a fixed block partitioning.
	 * */
	return RSB_DO_FLAG_HAS_INTERSECTION(flags,( RSB_FLAG_WANT_FIXED_BLOCKING_VBR | RSB_FLAG_WANT_BCSS_STORAGE | RSB_FLAG_WANT_COO_STORAGE ));
}

#ifdef RSB_FLAG_WANT_LINKED_STORAGE
rsb_bool_t rsb__have_linked_storage(const rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given flags are for a linked lists storage.
	 */
#ifdef RSB_FLAG_WANT_LINKED_STORAGE
	return RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_LINKED_STORAGE);
#else /* RSB_FLAG_WANT_LINKED_STORAGE */
	return RSB_BOOL_FALSE;
#endif /* RSB_FLAG_WANT_LINKED_STORAGE */
}
#endif

rsb_bool_t rsb__is_terminal_recursive_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given matrix is terminal
	 * FIXME : is this function really needed ?
	 * FIXME : should return one for terminal of non recursive ?
	 * TODO rsb__is_terminal_recursive_matrix -> rsb__is_terminal_matrix or rsb__is_leaf_matrix
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	int smc = 0;

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			++smc;

	return (smc==0);
}

rsb_bool_t rsb__is_recursive_matrix(rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given flags are for a recursive storage.
	 */
	return (RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING));
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__is_fixed_block_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return nonzero if the given matrix was partitioned with a fixed blocking,
	 *         thus enabling optimized operations on it.
	 * */
	if(!mtxAp)
		return 0;
	
	if( rsb__have_fixed_blocks_matrix_flags(mtxAp->flags) )
		return 1;

	/* FIXME : is this ok ? */
	if(
#ifdef RSB_MATRIX_STORAGE_VBR
		mtxAp->matrix_storage & RSB_MATRIX_STORAGE_VBR ||
#endif /* RSB_MATRIX_STORAGE_VBR */
#ifdef RSB_MATRIX_STORAGE_VBC
		mtxAp->matrix_storage & RSB_MATRIX_STORAGE_VBC ||
#endif/* RSB_MATRIX_STORAGE_VBC */
		0 )
		return RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_FIXED_BLOCKING_VBR);

	else
	return 
#ifdef RSB_MATRIX_STORAGE_BCSR
		mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSR ||
#endif /* RSB_MATRIX_STORAGE_BCSR */
#ifdef RSB_MATRIX_STORAGE_BCSC
		mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSC ||
#endif /* RSB_MATRIX_STORAGE_BCSC */
		0;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_bool_t rsb__util_are_flags_suitable_for_optimized_1x1_constructor(rsb_flags_t flags)
{
	/*!
	 	\ingroup gr_internals
		FIXME : temporary
	*/
	return	(RSB_DO_FLAG_HAS(  flags,RSB_FLAG_WANT_BCSS_STORAGE)  &&
		 (!RSB_DO_FLAG_HAS(flags,RSB_FLAG_AUTO_BLOCKING))    &&
		 (!RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING )) );
}

rsb_bool_t rsb__is_root_matrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 	\ingroup gr_internals
	*/
	return (!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_NON_ROOT_MATRIX))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
}

static rsb_bool_t rsb__mtx_chk_intrnl(const struct rsb_mtx_t *mtxAp, const struct rsb_mtx_t *mtxRp)
{
	/*!
	 	\ingroup gr_internals

		Check function, useful when debugging or re-developing core functionality.

		FIXME: will die in the presence of the RSB_FLAG_FORTRAN_INDICES_INTERFACE flag
		TODO: move to rsb__mtx_check.c or rsb_chk.c
		TODO: invoke rsb__check_bounds.
	*/
	rsb_bool_t return_ok = RSB_BOOL_FALSE;

	if( RSB_INVALID_COO_INDEX(mtxAp->Mdim)  || RSB_INVALID_COO_INDEX(mtxAp->nr) || RSB_INVALID_NNZ_INDEX(mtxAp->nnz) )
	{
		RSB_PERR_GOTO(err,RSB_ERRM_BADDIM);
	}
	
	if ( rsb__get_hermitian_flag(mtxAp) && rsb__get_symmetry_flag(mtxAp) )
	{
		RSB_PERR_GOTO(err,"bad flags: matrix at once hermitian and symmetric!\n");
	}

	if( (mtxAp)->all_leaf_matrices_n )
	if( rsb__submatrices(mtxAp) < (mtxAp)->all_leaf_matrices_n )
	{
	       	RSB_ERROR("more leaf submatrices %ld than overall %ld is inconsistent!", (mtxAp)->all_leaf_matrices_n, rsb__submatrices(mtxAp) );
		RSB_ERROR(RSB_ERRM_NL);
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

#if RSB_MERCY_FOR_LEGACY_INTERFACE
	if(
			rsb_get_matrix_n_rows(mtxAp)!=mtxAp->nr || 
			rsb_get_matrix_n_columns(mtxAp)!=mtxAp->nc || 
			rsb_get_matrix_n_rows(NULL)!=RSB_DEFAULT_UNDEFINED_COO_VALUE  || 
			rsb_get_matrix_n_columns(NULL)!=RSB_DEFAULT_UNDEFINED_COO_VALUE 
			)
#else /* RSB_MERCY_FOR_LEGACY_INTERFACE */
	if(0)
#endif /* RSB_MERCY_FOR_LEGACY_INTERFACE */
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	/* if(rsb__is_recursive_matrix(mtxAp->flags) && ! RSB_DO_FLAG_HAS_INTERSECTION(mtxAp->flags,RSB_FLAG_NON_ROOT_MATRIX ) ) */
	if( rsb__is_recursive_matrix(mtxAp->flags) )
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix;

		/* RSB_ASSERT( mtxAp->VA == NULL && mtxAp->bpntr == NULL && mtxAp->bindx == NULL ); */

		if( rsb__is_root_matrix(mtxAp) ) /* 20140921 */
		{
			rsb_submatrix_idx_t smi;

			RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
			{
				if(!submatrix)
				{
			       		RSB_PERR_GOTO(err,"leaf node %d (of %d) is NULL !?\n",smi,mtxAp->all_leaf_matrices_n);
				}

				if( rsb__is_root_matrix(submatrix) )
				{
			       		RSB_PERR_GOTO(err,"leaf node %d (of %d) has root flag in flags 0x%x=%d !?\n",smi,mtxAp->all_leaf_matrices_n,submatrix->flags,submatrix->flags);
				}

				if( submatrix->nnz > 0 && submatrix->VA == NULL )
				{
			       		RSB_PERR_GOTO(err,"leaf node %d (of %d) has %d nonzeroes and VA=NULL !?\n",smi,mtxAp->all_leaf_matrices_n,submatrix->nnz);
				}

				if( rsb__is_recursive_matrix(submatrix->flags) )
				{
			       		RSB_PERR_GOTO(err,"leaf node %d (of %d) has quad-partitioning flag in flags 0x%x=%d !?\n",smi,mtxAp->all_leaf_matrices_n,submatrix->flags,submatrix->flags);
				}

				if( submatrix->all_leaf_matrices != NULL )
				{
			       		RSB_PERR_GOTO(err,"leaf node %d (of %d) has a non-NULL submatrices pointer !?\n",smi,mtxAp->all_leaf_matrices_n);
				}

				if( submatrix->all_leaf_matrices_n > 0 )
				{
			       		RSB_PERR_GOTO(err,"leaf node %d (of %d) has a submatrices count of %d !?\n",smi,mtxAp->all_leaf_matrices_n,submatrix->all_leaf_matrices_n);
				}
			}
		}

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			if(RSB__IS_SUBM_MISPLACED(submatrix,mtxAp) )
			{
				RSB_STDOUT(RSB_PRINTF_MTX_SUMMARIZED_ARGS("Maybe misplaced submatrix: ",submatrix,RSB_ERRM_NL));
			       	RSB_PERR_GOTO(err,"submatrix %p at %d %d (root+%zd) seems misplaced\n",submatrix,submatrix->roff,submatrix->coff,submatrix-mtxRp);
			}

			if(!rsb__mtx_chk_intrnl(submatrix,mtxRp))
			{
			       	RSB_PERR_GOTO(err,"submatrix %p at %d %d (root+%zd) seems corrupted\n",submatrix,submatrix->roff,submatrix->coff,submatrix-mtxRp);
			}
		}
	}
	else
	{
		rsb_nnz_idx_t n;

		if(!RSB_IS_MATRIX_STORAGE_ALLOWED_FOR_LEAF(mtxAp->matrix_storage))
		{
		       	RSB_PERR_GOTO(err,RSB_ERRM_ES)
	       	}

		if(
					RSB_INVALID_COO_INDEX(mtxAp->roff) || 
					RSB_INVALID_COO_INDEX(mtxAp->coff) || 
					RSB_INVALID_COO_INDEX(mtxAp->nr) || 
					RSB_INVALID_COO_INDEX(mtxAp->nc) || 
					(mtxAp->roff>mtxAp->broff) || 
					(mtxAp->coff>mtxAp->bcoff) || 
					(mtxAp->nr<mtxAp->bm) || 
					(mtxAp->nc<mtxAp->bk) || 
					0
					)
		{
			RSB_PERR_GOTO(err,"submatrix bounds ([%d .. %d ... %d .. %d, %d .. %d ... %d .. %d]) are wrong!\n",
					mtxAp->roff,
					mtxAp->broff,
					mtxAp->bm,
					mtxAp->nr,
					mtxAp->coff,
					mtxAp->bcoff,
					mtxAp->bk,
					mtxAp->nc
					);
		}

		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
		{
			if(rsb__is_coo_matrix(mtxAp))
			if( (!RSB_INDICES_FIT_IN_HALFWORD(mtxAp->nr,mtxAp->nc)) || 0)
			{
				RSB_PERR_GOTO(err,"coo submatrix bounds are wrong, given the halfword indices!\n");
			}

			if(rsb__is_csr_matrix(mtxAp))
			if( (!RSB_INDEX_FIT_IN_HALFWORD(mtxAp->nc)) || 0)
			{
				RSB_PERR_GOTO(err,"csr submatrix bounds are wrong, given the halfword indices!\n");
			}
		}

		if(rsb__is_coo_matrix(mtxAp))
		{
			const rsb_coo_idx_t mai = 0; /* minimal allowed index */ /* FIXME: if one-based, this shall be 1 ! */

			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
			{
				RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)

				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					if( IA[n] < mai || JA[n] < mai )
					{
						RSB_PERR_GOTO(err,"negative halfword COO indices @%d: %d<%d || %d,%d!\n", n,IA[n],mai,JA[n],mai);
					}

					if( IA[n]>=mtxAp->Mdim || JA[n]>=mtxAp->mdim )
					{
						RSB_PERR_GOTO(err,"bad halfword COO indices @%d: %d>=%d || %d>=%d!\n", n,IA[n],mtxAp->Mdim,JA[n],mtxAp->mdim); 
					}
				}

				if( rsb__util_is_halfword_coo_array_sorted_up_partial_order(IA,mtxAp->nnz) != RSB_BOOL_TRUE )
				{
					RSB_PERR_GOTO(err,"halfword COO input is not sorted! \n");
				}
				goto ok;
			}
			else
			{
				RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)

				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					if( IA[n] < mai || JA[n] < mai )
					{
						RSB_PERR_GOTO(err,"negative fullword COO indices @%d: %d<%d || %d,%d!\n", n,IA[n],mai,JA[n],mai);
					}

					if( IA[n]>=mtxAp->Mdim || JA[n]>=mtxAp->mdim )
					{
						RSB_PERR_GOTO(err,"bad fullword COO indices @%d: %d>=%d || %d>=%d!\n",
							n,IA[n],mtxAp->Mdim,JA[n],mtxAp->mdim);
					}
				}

				if( rsb__util_is_nnz_array_sorted_up_partial_order(IA,mtxAp->nnz) != RSB_BOOL_TRUE )
				{
					RSB_PERR_GOTO(err,"fullword COO input is not sorted! \n");
				}

				goto ok;
			}
		}
	
		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_WANT_COO_STORAGE))
		    && !RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
			;
		//	goto ok;//{RSB_PERR_GOTO(err,"full word COO is not allowed on a leaf matrix!\n");}
		//
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE))
		{
			if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES_COO)))
			{
			}
			else
			{
				// FIXME: I am not sure whether this code is ever executed.
				RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					if( IA[n]>=mtxAp->Mdim || JA[n]>=mtxAp->mdim )
					{
						RSB_PERR_GOTO(err,"bad fullword COO indices @%d: %d>=%d || %d>=%d!\n",
							n,IA[n],mtxAp->Mdim,JA[n],mtxAp->mdim);
					}
				}
			}
			goto ok;
		}
		else
				;/* ok */
#if 0
//		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES_CSR))
//		    && RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES_COO)))
//				{RSB_PERR_GOTO(err,"both halfword COO and halfword CSR is not allowed on a leaf matrix!\n");}
//		else
//				;/* ok */
#endif

		if(!rsb__is_csr_matrix(mtxAp))
		{
			RSB_PERR_GOTO(err,"not a csr matrix ?\n");
		}
		if(mtxAp->element_count != mtxAp->nnz)
		{
			RSB_PERR_GOTO(err,RSB_ERRM_EM);
		}
#if RSB_WANT_DBC
		if(mtxAp->element_count != mtxAp->block_count)
		{
			RSB_PERR_GOTO(err,RSB_ERRM_EM);
		}
#endif
#if RSB_OLD_COO_CRITERIA
		if(!mtxAp->bpntr)
		{
			RSB_PERR_GOTO(err,"!bpntr!\n");
		}
		if(!mtxAp->bindx)
		{
			RSB_PERR_GOTO(err,"!bindx!\n");
		}
		if(mtxAp->bpntr[0]!=0)
		{
			RSB_PERR_GOTO(err,"bpntr[0]!=0!\n");
		}
#else
		if(mtxAp->nnz != 0 && !mtxAp->bpntr)
		{
			RSB_PERR_GOTO(err,"!bpntr!\n");
		}
		if(mtxAp->nnz != 0 && !mtxAp->bindx)
		{
			RSB_PERR_GOTO(err,"!bindx!\n");
		}
		if(mtxAp->bpntr && mtxAp->nnz > 0 && mtxAp->bpntr[0]!=0)
		{
			RSB_PERR_GOTO(err,"bpntr[0]!=0!\n");
		}
		if(mtxAp->nnz==0)
			goto ok;
#endif
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE))
		{
			if(!RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
			{
				if(RSB_SOME_ERROR( rsb__util_is_sorted_coo_as_row_major(mtxAp->bpntr,mtxAp->bindx,mtxAp->nnz,mtxAp->typecode,NULL,mtxAp->flags)) )
				{
					RSB_PERR_GOTO(err,"COO matrix seems unsorted!\n");
				}
			}
			else
			{
				/* FIXME: missing */		
			}
		}
		else
		{
			if(mtxAp->bpntr[mtxAp->Mdim]!=mtxAp->nnz)
			{
				RSB_PERR_GOTO(err,"%d=bpntr[Mdim]!=nnz=%d\n",(int)mtxAp->bpntr[mtxAp->Mdim],(int)mtxAp->nnz);
			}
			if(!rsb__util_is_nnz_array_sorted_up_partial_order(mtxAp->bpntr,mtxAp->Mdim+1))
			{
				RSB_PERR_GOTO(err,"bpntr seems unsorted!\n");
			}
		}

		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_CSR_RESERVED))
			for(n=0;RSB_LIKELY(n<mtxAp->Mdim);++n)
			if( RSB_UNLIKELY( mtxAp->bpntr[n+1] - mtxAp->bpntr[n] > mtxAp->mdim ) )
			{
				RSB_ERROR("invalid CSR pointer:  mtxAp->bpntr[%d+1] - mtxAp->bpntr[%d] > mtxAp->mdim: %d - %d > %d !\n",n,n,mtxAp->bpntr[n+1],mtxAp->bpntr[n],mtxAp->mdim);
				RSB_PERR_GOTO(err,"!\n");
			}


		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES_COO)))
		{
			/* FIXME: write me */
		}
		else
		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES_CSR)))
		{
			for(n=0;RSB_LIKELY(n<mtxAp->Mdim);++n)
			if(!rsb__util_is_halfword_coo_array_sorted_up(
						((rsb_half_idx_t*)mtxAp->bindx)+mtxAp->bpntr[n],mtxAp->bpntr[n+1]-mtxAp->bpntr[n]))
			{
				RSB_PERR_GOTO(err,"(halfword) bindx seems unsorted!\n");
			}
		}
		else
		{
			if(RSB_SOME_ERROR(rsb__csr_chk(mtxAp->bpntr,mtxAp->bindx,mtxAp->Mdim,mtxAp->mdim,mtxAp->nnz,0)))
			{
				RSB_PERR_GOTO(err,"CSR submatrix seems corrupt!\n");
			}
		}
	}
ok:
	return_ok = RSB_BOOL_TRUE;
err:
	if( return_ok == RSB_BOOL_FALSE)
	{
		RSB_ERROR(RSB_PRINTF_MTX_SUMMARIZED_ARGS(RSB_ERRM_BS,mtxAp,RSB_ERRM_NL));
	}
	return return_ok;
}

rsb_bool_t rsb__mtx_chk(const struct rsb_mtx_t *mtxAp)
{
	return rsb__mtx_chk_intrnl(mtxAp, mtxAp);
}

rsb_bool_t rsb__do_is_matrix_binary_loaded(const struct rsb_mtx_t * mtxAp)
{
	rsb_bool_t is_bio; // binary I/O matrix
#if 0
	struct rsb_mtx_t *fsm = rsb__do_get_first_submatrix(mtxAp);
	is_bio=!((long)mtxAp<((long)fsm->bpntr) || (long)(mtxAp)>=((long)fsm->bpntr+mtxAp->nnz));
#else
	is_bio = RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_FIX_FOR_BINARY_LOADED_MATRIX);
#endif
	return is_bio;
}


/* @endcond */
