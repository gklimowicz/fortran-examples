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
 * This source file contains matrix generating functions.
 * */

#ifndef RSB_GEN_H_INCLUDED
#define RSB_GEN_H_INCLUDED

#include "rsb_internals.h"	/* rsb_coo_mtx_t */
/* #define RSB_HSQUARE(X) (((X)*(X+1))/2) */
#define RSB_HSQUARE(X) ((X)<2?(((X)*(X+1))/2):( ((X)*(X-1))/2 + X ) )
#define RSB_NNZ_OF_BANDED(DIM,L,U) (((DIM)*(1+(L)+(U)))-(RSB_HSQUARE(L))-(RSB_HSQUARE(U)))
/* struct rsb_mtx_t * rsb_generate_diagonal(const rsb_coo_idx_t rows, double * timep, rsb_type_t typecode); */
struct rsb_mtx_t * rsb__generate_banded(const rsb_blk_idx_t br, const rsb_blk_idx_t bc, const rsb_coo_idx_t rows, const rsb_coo_idx_t cols, rsb_coo_idx_t bw, double * timep, rsb_type_t typecode);
#if RSB_WANT_DO_LOCK_TEST
struct rsb_mtx_t * rsb__generate_dense_lower_triangular(const rsb_coo_idx_t dim, double * timep, rsb_type_t typecode);
#endif /* RSB_WANT_DO_LOCK_TEST */
rsb_err_t rsb__generate_dense_lower_triangular_coo(rsb_nnz_idx_t dim, rsb_nnz_idx_t spacing, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_nnz_idx_t *nnzp, rsb_type_t typecode);
rsb_err_t rsb__generate_dense_full(rsb_nnz_idx_t dim_r, rsb_nnz_idx_t dim_c, rsb_nnz_idx_t spacing, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_nnz_idx_t *nnzp, rsb_type_t typecode);
rsb_err_t rsb__generate_blocked_banded_coo(rsb_nnz_idx_t dim, rsb_nnz_idx_t spacing, rsb_nnz_idx_t lbw, rsb_nnz_idx_t ubw, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_nnz_idx_t *nnzp, rsb_type_t typecode);
struct rsb_mtx_t * rsb__generate_blocked_banded(const rsb_blk_idx_t br, const rsb_blk_idx_t bc, const rsb_coo_idx_t rows, const rsb_coo_idx_t cols, const rsb_coo_idx_t bw, double * timep, rsb_type_t typecode,rsb_bool_t want_lowtri);
void rsb__do_fill_with_diag(void *VA, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, rsb_coo_idx_t ioff, rsb_coo_idx_t joff, rsb_nnz_idx_t nzoff, rsb_type_t typecode, rsb_nnz_idx_t nnz);
rsb_err_t rsb__generate_blocked_banded_mtx(rsb_nnz_idx_t dim, rsb_nnz_idx_t spacing, rsb_nnz_idx_t lbw, rsb_nnz_idx_t ubw, struct rsb_mtx_t ** mtxApp, rsb_type_t typecode);
#endif /* RSB_GEN_H_INCLUDED */
/* @endcond */
