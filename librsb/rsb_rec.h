/*                                                                                                                            

Copyright (C) 2008-2019 Michele Martone

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
 * @file
 * @brief Recursion handling code
 * @author Michele Martone
 * */

#ifndef RSB_REC_H_INCLUDED
#define RSB_REC_H_INCLUDED

#include "rsb_internals.h"

enum rsb_op_t{	/* FIXME : temporary, experimental  */
            rsb_op_spmv = 1,
            rsb_op_spsvl = 2,
            rsb_op_spsvlt = 3,
            rsb_op_spsvu = 4,
            rsb_op_spsvut = 5,
            rsb_op_get_csr = 6,
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
            rsb_op_ata = 7,
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
            rsb_op_nop = 0
};
#define RSB_PTR_SHIFT(DPTR,MPTR,PPTR,TC)	/* FIXME: temporarily here */ \
		if( (PPTR) > (MPTR) )	\
		{ ( DPTR) = TC (( (const rsb_byte_t *) (DPTR) ) + ( ( (const rsb_byte_t *) (PPTR) ) - ( (const rsb_byte_t *) (MPTR) )) ); }	\
		else	\
		{ ( DPTR) = TC (( (const rsb_byte_t *) (DPTR) ) - ( ( (const rsb_byte_t *) (MPTR) ) - ( (const rsb_byte_t *) (PPTR) )) ); }

#define RSB_TMP_OVERALLOC_MTX 4 /* 1 < RSB_TMP_OVERALLOC_MTX < 4; a temporary measure */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__init_set_quad_submatrices_info(const struct rsb_mtx_partitioning_info_t * pinfop, struct rsb_mtx_t ** matrices, rsb_nnz_idx_t uuk, rsb_nnz_idx_t mk, rsb_nnz_idx_t uk, rsb_nnz_idx_t lk, rsb_nnz_idx_t llk, rsb_coo_idx_t mB, rsb_coo_idx_t kB, rsb_coo_idx_t roff, rsb_coo_idx_t coff);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_err_t rsb__get_array_of_leaf_matrices(struct rsb_mtx_t *mtxAp, struct rsb_translated_matrix_t ** tmatricesp, rsb_submatrix_idx_t *countp);
rsb_err_t rsb__fill_array_of_leaf_matrices(const struct rsb_translated_matrix_t *tmatrix, struct rsb_translated_matrix_t *matrices, rsb_submatrix_idx_t * n);
rsb_err_t rsb__sort_array_of_leaf_matrices(const struct rsb_translated_matrix_t *rmatrix,struct rsb_translated_matrix_t *matrices, rsb_submatrix_idx_t n, enum rsb_op_t op);
int rsb__compar_rcsr_matrix_for_spsvl(const void * ap, const void * bp);
size_t rsb__get_index_storage_amount(const struct rsb_mtx_t *mtxAp);
rsb_submatrix_idx_t rsb__get_diagonal_elements_count(const struct rsb_mtx_t *mtxAp);
rsb_submatrix_idx_t rsb__get_diagonal_submatrices_count(const struct rsb_mtx_t *mtxAp);
rsb_err_t rsb__sort_array_of_leaf_matrices_for_ussv(const struct rsb_mtx_t * mtxAp, struct rsb_translated_matrix_t *leaf_matrices, rsb_submatrix_idx_t n, rsb_trans_t transl);
rsb_err_t rsb__leaves_merge(struct rsb_mtx_t * RSB_RESTRICT mtxAp, rsb_submatrix_idx_t manp, rsb_time_t * RSB_RESTRICT stp, rsb_time_t *RSB_RESTRICT atp, rsb_time_t *RSB_RESTRICT ltp, const int wv, int kc);
rsb_err_t rsb__leaves_merge_multiple(struct rsb_mtx_t *mtxAp, rsb_time_t *stp, rsb_time_t *atp, rsb_time_t *ltp, const int wv, int kc);
rsb_err_t rsb__mtx_split(struct rsb_mtx_t * RSB_RESTRICT mtxAp, rsb_submatrix_idx_t manp, rsb_time_t * RSB_RESTRICT stp, rsb_time_t * RSB_RESTRICT atp, rsb_time_t * RSB_RESTRICT ltp, const int wv, int kc);
rsb_err_t rsb__mtx_realloc_with_spare_leaves(struct rsb_mtx_t **mtxApp, rsb_submatrix_idx_t slc);

#endif /* RSB_REC_H_INCLUDED */
/* @endcond */
