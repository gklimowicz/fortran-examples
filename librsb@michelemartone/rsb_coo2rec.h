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
/*!
 * @file
 * @brief  Recursive Sparse matrices assembling code.
 * @author Michele Martone
 * */
#ifndef RSB_COO2RCSR_H_INCLUDED
#define RSB_COO2RCSR_H_INCLUDED
#include "rsb_common.h"
#define RSB_HAVE_GOOD_PARMS_FOR_RCSR(R,C,NNZ,FLAGS) \
	((!(NNZ==0)) || (RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_UNIT_DIAG_IMPLICIT))) 

#define RSB_HAVE_GOOD_PARMS_FOR_EMPTY(R,C,NNZ,FLAGS) ( \
		((NNZ==0) && RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_WANT_COO_STORAGE)) \
		)
#define RSB_HAVE_GOOD_PARMS_FOR_IN_PLACE_RCSR(R,C,NNZ,FLAGS) ( \
	!RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER) && \
	((!RSB_DO_TOOFEWNNZFORRCSR(NNZ,R)) || \
	( \
		(RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_WANT_COO_STORAGE)) || \
		(RSB_DO_FLAGS_EXTRACT_STORAGE(FLAGS)==RSB_FLAG_NOFLAGS) \
	)) && \
	RSB_HAVE_GOOD_PARMS_FOR_RCSR(R,C,NNZ,FLAGS) \
	)
#define RSB_DO_TOOFEWNNZFORRCSR(NNZ,M) ((NNZ)<(2*(M+1)))	/*  */
#define RSB_DO_TOOFEWNNZFORCSR(NNZ,M)  ((NNZ)<(1*(M+1)))	/*  */
struct rsb_mtx_t * rsb__allocate_recursive_sparse_matrix_from_row_major_coo(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_err_t *errvalp);
void rsb__do_set_in_place_submatrices_offsets(struct rsb_mtx_t *RSB_RESTRICT submatrices, rsb_submatrix_idx_t cmc, rsb_char_t *RSB_RESTRICT  VA, rsb_coo_idx_t *RSB_RESTRICT  IA, rsb_coo_idx_t *RSB_RESTRICT JA, size_t el_size);
rsb_err_t rsb__do_switch_recursive_matrix_to_fullword_storage(struct rsb_mtx_t * mtxAp);
rsb_err_t rsb__project_rsb_to_coo(const struct rsb_mtx_t *mtxAp, struct rsb_coo_mtx_t *coop);
rsb_err_t rsb__compute_bounded_box(struct rsb_mtx_t * mtxAp);
#define RSB_STDOUT_MATRIX_SUMMARY_ARGS(M) RSB_PRINTF_MTX_SUMMARY_ARGS(M)
#define RSB_STDOUT_MATRIX_SUMMARY(M)  RSB_STDOUT(RSB_STDOUT_MATRIX_SUMMARY_ARGS(M))

#define RSB_FPRINTF_MATRIX_SUMMARY(FP,M)  RSB_FPRINTF(FP,RSB_STDOUT_MATRIX_SUMMARY_ARGS(M))

#define RSB_PRINTF_COO_MATRIX_SUMMARY_ARGS(CM)  \
			"(%ld x %ld)[%p] @ (? , ?) (%ld nnz, %.2lg nnz/r) flags 0x??, typecode: %x:",		\
				(long int)(CM)->nr, (long int)(CM)->nc, (const void*)(CM),								\
			       	(long int)(CM)->nnz,									\
			       	((double)(CM)->nnz)/(CM)->nr,							\
				CM->typecode

#define RSB_STDOUT_COO_MATRIX_SUMMARY_ARGS(M) RSB_PRINTF_COO_MATRIX_SUMMARY_ARGS(M)

#if RSB_ALLOW_STDOUT
#define RSB_ERROR_MATRIX_SUMMARY(M)  RSB_STDOUT(RSB_STDOUT_MATRIX_SUMMARY_ARGS(M))
#define RSB_STDOUT_COO_MATRIX_SUMMARY(CM)  RSB_STDOUT(RSB_STDOUT_COO_MATRIX_SUMMARY_ARGS(CM))
#define RSB_INFO_MATRIX_SUMMARY  RSB_STDOUT_MATRIX_SUMMARY
#else /* RSB_ALLOW_STDOUT */
#define RSB_INFO_MATRIX_SUMMARY(M)  RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS  
#define RSB_ERROR_MATRIX_SUMMARY(M) RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS  
#endif /* RSB_ALLOW_STDOUT */
#define RSB_OLD_COO_CRITERIA (RSB_LIBRSB_VER < 10300)

#endif /* RSB_COO2RCSR_H_INCLUDED */
/* @endcond */
