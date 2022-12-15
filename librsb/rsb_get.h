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
 * @author Michele Martone
 * @brief
 * This source file contains matrix getter functions.
 * */
#ifndef RSB_GET_H_INCLUDED
#define RSB_GET_H_INCLUDED

#include "rsb_internals.h"

rsb_err_t rsb__dodo_get_csr(const struct rsb_mtx_t *mtxAp, rsb_byte_t ** VA, rsb_nnz_idx_t ** RP, rsb_coo_idx_t ** JA);
rsb_err_t rsb__do_get_csc(const struct rsb_mtx_t *mtxAp, rsb_byte_t ** VA, rsb_nnz_idx_t ** RP, rsb_coo_idx_t ** IA);
rsb_err_t rsb__do_get_coo(const struct rsb_mtx_t *mtxAp, rsb_byte_t ** VA, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, rsb_flags_t flags);
rsb_err_t rsb__do_get_coo_noalloc(const struct rsb_mtx_t *mtxAp, rsb_byte_t * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * nnzp, rsb_flags_t flags );
rsb_err_t rsb__do_get_row_dense(const struct rsb_mtx_t *mtxAp, void * row, rsb_blk_idx_t rowindex);
rsb_err_t rsb__do_get_rows_dense(const struct rsb_mtx_t *mtxAp, void * row, rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t  *rnz, rsb_flags_t flags);
rsb_err_t rsb__do_get_rows_sparse_rec(const struct rsb_mtx_t *mtxAp , void * RSB_RESTRICT VA , rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT rnz, rsb_coo_idx_t ioff, rsb_coo_idx_t joff);
rsb_err_t rsb__do_get_columns_sparse(const struct rsb_mtx_t *mtxAp , void * RSB_RESTRICT VA , rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT CP, rsb_coo_idx_t ioff, rsb_coo_idx_t joff);
rsb_err_t rsb__do_get_rows_nnz(const struct rsb_mtx_t *mtxAp, rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_nnz_idx_t  *rnz);
rsb_coo_idx_t rsb__mtx_strict_diagonal_size(const struct rsb_mtx_t *mtxAp);
rsb_err_t rsb__do_infinity_norm(const struct rsb_mtx_t *mtxAp , void * infinity_norm, const rsb_bool_t do_testing, rsb_trans_t transA);
double rsb__do_get_matrix_fillin(const struct rsb_mtx_t *mtxAp);
rsb_nnz_idx_t rsb__do_get_matrix_nnz(const struct rsb_mtx_t *mtxAp);
rsb_long_t rsb__submatrices(const struct rsb_mtx_t * mtxAp);
void rsb__get_blocking_size(const struct rsb_mtx_t * mtxAp, rsb_blk_idx_t *brp, rsb_blk_idx_t *bcp);
rsb_err_t rsb__get_physical_blocking_size(const struct rsb_mtx_t * mtxAp, rsb_blk_idx_t *brp, rsb_blk_idx_t *bcp);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__get_blocking_size_as_row_major(const struct rsb_mtx_t * mtxAp, rsb_blk_idx_t *bMp, rsb_blk_idx_t *bmp);
size_t rsb__do_get_max_blocksize(const struct rsb_mtx_t * mtxAp);
/* void rsb__do_extract_nnz_from_block(void * blockpointer, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t br, rsb_coo_idx_t bc, rsb_coo_idx_t baserow, rsb_coo_idx_t basecolumn, rsb_type_t typecode, size_t el_size, rsb_nnz_idx_t *nnzp ); */
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_submatrix_idx_t rsb__get_recursive_matrix_depth(const struct rsb_mtx_t *mtxAp);
rsb_flags_t rsb__get_symmetry_flag(const struct rsb_mtx_t *mtxAp);
rsb_flags_t rsb__get_symmetry_type_flag(const struct rsb_mtx_t *mtxAp);
rsb_flags_t rsb__get_hermitian_flag(const struct rsb_mtx_t *mtxAp);
RSB_INLINE rsb_coo_idx_t rsb__do_get_rows_of(const struct rsb_mtx_t *mtxAp, rsb_trans_t transA);
RSB_INLINE rsb_coo_idx_t rsb__do_get_columns_of(const struct rsb_mtx_t *mtxAp, rsb_trans_t transA);
rsb_flags_t rsb__get_diagonal_type_flag(const struct rsb_mtx_t *mtxAp);
rsb_err_t rsb__do_absolute_rows_sums( const struct rsb_mtx_t * mtxAp , void * d);
rsb_err_t rsb__do_rows_sums_inner(const struct rsb_mtx_t *mtxAp , rsb_byte_t * row_sums, rsb_bool_t do_testing, rsb_trans_t transA);
rsb_err_t rsb__do_absolute_columns_sums( const struct rsb_mtx_t * mtxAp , void * d);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_get_block_sparse_pattern(const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t fr, rsb_coo_idx_t lr, rsb_coo_idx_t fc, rsb_coo_idx_t lc, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, const rsb_coo_idx_t * IREN, const rsb_coo_idx_t * JREN, rsb_nnz_idx_t *rnz, rsb_flags_t flags );
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_err_t rsb__do_get_block_sparse(const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t fr, rsb_coo_idx_t lr, rsb_coo_idx_t fc, rsb_coo_idx_t lc, const rsb_coo_idx_t * IREN, const rsb_coo_idx_t * JREN, rsb_nnz_idx_t *rnz, rsb_flags_t flags );
rsb_nnz_idx_t rsb__do_get_block_nnz(const struct rsb_mtx_t *mtxAp, rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_blk_idx_t fc, rsb_blk_idx_t lc, rsb_flags_t flags, rsb_err_t * errvalp);

#endif /* RSB_GET_H_INCLUDED */
/* @endcond */
