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
 * @brief Sparse matrices assembling code.
 * @author Michele Martone
 * */

#include "rsb_common.h"
#define RSB_WANT_PRINT_WARNING_ON_DISCARDED_NNZ 0
RSB_INTERNALS_COMMON_HEAD_DECLS

struct rsb_mtx_t * rsb__mtx_alloc_inner(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_type_t typecode, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_blk_idx_t Mb, rsb_blk_idx_t Kb, rsb_flags_t flags, rsb_err_t * errvalp)
{
	/*!

	   Allocates a blocked sparse matrix being recursively partitioned, 
	   but in a data structure which is specified by flags, and 
	   thus not necessarily with exact BCSR/BCSC leaves.

	   \return a valid matrix pointer or NULL
	*/

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;

	RSB_DEBUG_ASSERT(roff>=-1 && coff>=-1); /* for Fortran */

	if((m==0 || k==0) && nnz>0)
	{
		/* as a special case, we detect the m and k boundaries, if nnz>0 and m or k are zero */
		/* TODO: shall use rsb__util_coo_alloc_copy_and_stats instead */
		if(m==0 && IA) {m = rsb__util_find_coo_max_index_val(IA,nnz)+roff+1;}
		if(k==0 && JA) {k = rsb__util_find_coo_max_index_val(JA,nnz)+coff+1;}
		//printf("rc %d %d %d \n",m,k,nnz);
	}

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
	{
		RSB_PERR_GOTO(err,"!\n");
	}

	if(RSB__FLAG_HAS_UNSPECIFIED_TRIANGLE(flags))
		RSB_DO_FLAG_ADD(flags,rsb__do_detect_and_add_triangular_flags(IA,JA,nnz,flags));
	
	if(roff && IA)
		rsb__util_coo_array_add(IA,nnz,roff);
	if(coff && JA)
		rsb__util_coo_array_add(JA,nnz,coff);

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORT_INPUT);
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_OWN_PARTITIONING_ARRAYS);	/* this is in order to free p_r and p_c with the matrix itself, and ignore original flag on this topic */
	if(
			(m<RSB_MIN_MATRIX_DIM||k<RSB_MIN_MATRIX_DIM||nnz<RSB_MIN_MATRIX_NNZ) ||
			(m>RSB_MAX_MATRIX_DIM||k>RSB_MAX_MATRIX_DIM||nnz>RSB_MAX_MATRIX_NNZ) 
	)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
#if RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT
	if(RSB_DO_TOOFEWNNZFORRCSR(nnz,RSB_MIN(m,k)))
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_COO_STORAGE);
#endif /* RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT */
	if( RSB_HAVE_GOOD_PARMS_FOR_IN_PLACE_RCSR(m,k,nnz,flags)
#if RSB_ALLOW_EMPTY_MATRICES
	|| (nnz==0)
#endif /* RSB_ALLOW_EMPTY_MATRICES */
			)
	{
		if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_NON_ROOT_MATRIX))
		{
			mtxAp = rsb__allocate_recursive_sparse_matrix_from_row_major_coo(VA,IA,JA,m,k,nnz,typecode,NULL,flags,errvalp);
			if(errvalp && RSB_SOME_ERROR(*errvalp))
				RSB_ERROR("%s\n", rsb__get_errstr_ptr(*errvalp));
			return mtxAp;
		}
	}
#if RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT
	errval = RSB_ERR_INTERNAL_ERROR;
	RSB_ERROR("Trying to call unsupported parameters combination (nr:%ld nc:%ld nnz:%ld)!\n",(long int)m,(long int)k,(long int)nnz);
	rsb__debug_print_flags(flags);
	RSB_PERR_GOTO(err,RSB_ERRM_INTERNAL_ERROR );
#endif /* RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT */

	if(mtxAp)
		return mtxAp;
	else
		goto err;
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return NULL;
}

rsb_err_t rsb__do_cleanup_nnz(void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t *onnzp, rsb_type_t typecode, rsb_flags_t flags)
{
	/* 
	 * onnzp can be NULL
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(nnz==0) /* diagonal implicit, for example */
		goto ok;

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
	{
		rsb_nnz_idx_t discarded = 0, gap = 0;
		errval = rsb__weed_out_diagonal(VA,IA,JA,nnz,typecode,&gap,&discarded);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,"failed compacting non diagonal elements !\n");
		}
		RSB_DEBUG_ASSERT(discarded>=0);
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz));
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz-discarded));
		if(RSB_WANT_PRINT_WARNING_ON_DISCARDED_NNZ && discarded>0)
			;//RSB_INFO("#RSB_FLAG_UNIT_DIAG_IMPLICIT (EXPERIMENTAL) : discarded %zd diagonal elements\n",(rsb_printf_int_t)discarded);
		nnz -= discarded;
	}

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_TRIANGULAR) && roff==coff && roff==0)
	{
		rsb_nnz_idx_t discarded = 0, gap = 0;
		errval = rsb__weed_out_non_lowtri(VA,IA,JA,nnz,typecode,&gap,&discarded);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,"failed compacting non lower triangular elements !\n");
		}
		RSB_DEBUG_ASSERT(discarded>=0);
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz));
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz-discarded));
		if(RSB_WANT_PRINT_WARNING_ON_DISCARDED_NNZ && discarded>0)
			RSB_INFO("#RSB_FLAG_LOWER_TRIANGULAR (EXPERIMENTAL) : discarded %zd non lower triangular\n",(rsb_printf_int_t)discarded);
		nnz -= discarded;
	}

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_TRIANGULAR) && roff==coff && roff==0)
	{
		rsb_nnz_idx_t discarded = 0, gap = 0;
		errval = rsb__weed_out_non_upptri(VA,IA,JA,nnz,typecode,&gap,&discarded);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,"failed compacting non upper triangular elements !\n");
		}
		RSB_DEBUG_ASSERT(discarded>=0);
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz));
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz-discarded));
		if(RSB_WANT_PRINT_WARNING_ON_DISCARDED_NNZ && discarded>0)
			RSB_INFO("#RSB_FLAG_UPPER_TRIANGULAR (EXPERIMENTAL) : discarded %zd non upper triangular\n",(rsb_printf_int_t)discarded);
		nnz -= discarded;
	}

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_DISCARD_ZEROS))
	{
		rsb_nnz_idx_t discarded = 0, gap = 0;
		errval = rsb__util_compact_nonzeros(VA,IA,JA,nnz,typecode,&gap,&discarded,RSB_FLAG_NOFLAGS);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,"failed compacting nonzeros!\n");
		}
		RSB_DEBUG_ASSERT(discarded>=0);
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz));
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz-discarded));
		if(RSB_WANT_PRINT_WARNING_ON_DISCARDED_NNZ && discarded>0)
			RSB_INFO("#RSB_FLAG_DISCARD_ZEROS (EXPERIMENTAL) : discarded %zd zeros\n",(rsb_printf_int_t)discarded);
		nnz -= discarded;
	}

	if(1)
	{
		rsb_nnz_idx_t discarded = 0, gap = 0;
		errval = rsb__do_util_compact_out_of_range(VA,IA,JA,nnz,roff,coff,m,k,typecode,&gap,&discarded);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,"failed compacting out of range!\n");
		}
		RSB_DEBUG_ASSERT(discarded>=0);
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz));
		RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz-discarded));
		if(RSB_WANT_PRINT_WARNING_ON_DISCARDED_NNZ && discarded>0)
			RSB_INFO("#: discarded %zd nonzeroes with out of range coordinates\n",(rsb_printf_int_t)discarded);
		nnz -= discarded;
	}
ok:
	RSB_SET_IF_NOT_NULL(onnzp,nnz);
err:
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
