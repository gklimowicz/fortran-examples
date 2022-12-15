/*                                                                                                                            

Copyright (C) 2008-2015 Michele Martone

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
 * This source file contains functions for sparse triangular solve.
 * */

#ifndef RSB_SPTRSV_H_INCLUDED
#define RSB_SPTRSV_H_INCLUDED
#include "rsb_internals.h"		/* */

/* #define RSB_ENABLE_INNER_NRHS_SPSV 1 */
#define RSB_ENABLE_INNER_NRHS_SPSV RSB_ENABLE_INNER_NRHS_SPMV 
#if RSB_ENABLE_INNER_NRHS_SPSV
#define RSB_INNER_NRHS_SPSV_ARGS	,const rsb_int_t nrhs, /*const size_t outtot, const size_t rhstot,*/ const size_t outnri, const size_t rhsnri
#define RSB_INNER_NRHS_SPSV_ARGS_IDS	,nrhs/*,outtot,rhstot*/,outnri,rhsnri
#define RSB_OUTER_NRHS_SPSV_ARGS	,const rsb_int_t nrhs, const size_t outnri, const size_t rhsnri
#define RSB_OUTER_NRHS_SPSV_ARGS_IDS	,nrhs,outnri,rhsnri
#else
#define RSB_INNER_NRHS_SPSV_ARGS	
#define RSB_INNER_NRHS_SPSV_ARGS_IDS
#define RSB_OUTER_NRHS_SPSV_ARGS_IDS
#endif /* RSB_ENABLE_INNER_NRHS_SPSV */

#define RSB_DEFAULT_INNER_NRHS_SPSV_ARGS	,1,/*0,0,*/0,0
#define RSB_DEFAULT_OUTER_NRHS_SPSV_ARGS	,1,0,0

rsb_err_t rsb__do_get_submatrices_for_ussv( const struct rsb_mtx_t * mtxAp, struct rsb_translated_matrix_t ** all_leaf_matricesp, rsb_submatrix_idx_t * all_leaf_matrices_np, rsb_trans_t transT);
void rsb__submatrices_exclude_nontriangular(struct rsb_translated_matrix_t * all_leaf_matrices, rsb_submatrix_idx_t * all_leaf_matrices_np, const struct rsb_mtx_t * mtxAp);
rsb_err_t rsb__do_spsv(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, const void * Xp, rsb_coo_idx_t incX, void * Yp, rsb_coo_idx_t incY);
rsb_err_t rsb__do_spsv_general(rsb_trans_t transl, const void * alphap, const struct rsb_mtx_t * mtxAp, const void * x, rsb_coo_idx_t incx, void * y, rsb_coo_idx_t incy, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPSV_ARGS);
rsb_err_t rsb__do_get_submatrices_block_for_get_csr(const struct rsb_mtx_t * mtxAp, struct rsb_translated_matrix_t ** all_leaf_matricesp, rsb_submatrix_idx_t * all_leaf_matrices_np);
#endif /* RSB_SPTRSV_H_INCLUDED */
/* @endcond */
