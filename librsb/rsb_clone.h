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
/*
 * @author Michele Martone
 */
#ifndef RSB_CLONE_H_INCLUDED
#define RSB_CLONE_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb_common.h"

void * rsb__clone_area(const void *src, size_t size);
rsb_err_t rsb__util_coo_alloc(void **RSB_RESTRICT VAp, rsb_coo_idx_t ** RSB_RESTRICT IAp, rsb_coo_idx_t ** RSB_RESTRICT JAp, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_bool_t do_calloc);
rsb_err_t rsb__util_coo_alloc_copy_and_stats(void **RSB_RESTRICT VAp, rsb_coo_idx_t ** RSB_RESTRICT IAp, rsb_coo_idx_t ** RSB_RESTRICT JAp, const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t*RSB_RESTRICT mp, rsb_coo_idx_t*RSB_RESTRICT kp, rsb_nnz_idx_t nnz, rsb_nnz_idx_t ennz, rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t iflags, rsb_flags_t*RSB_RESTRICT flagsp);
void * rsb__clone_area_with_extra(const void *src, size_t csize, size_t bsize, size_t esize);
/* void * rsb__clone_area_parallel(const void *src, size_t size); */
struct rsb_mtx_t *rsb__mtx_clone_simple(const struct rsb_mtx_t *mtxAp);
rsb_err_t rsb__mtx_clone(struct rsb_mtx_t ** mtxBpp, rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_flags_t flags);
/*void * rsb__clone_inner(const struct rsb_mtx_t *mtxAp, struct rsb_mtx_t *new_matrix); */
rsb_err_t rsb__clone_coo(const struct rsb_mtx_t * mtxAp, rsb_trans_t transA, const void *alphap, rsb_type_t typecode, struct rsb_coo_mtx_t*dcoop, rsb_flags_t flags/*, rsb_extff_t eflags*/);
#if !RSB_AT_DESTROYS_MTX
rsb_err_t rsb__mtx_transplant_from_clone(struct rsb_mtx_t ** mtxDpp, struct rsb_mtx_t * mtxSp);
#endif  /* RSB_AT_DESTROYS_MTX */
size_t rsb__submatrices_max_ptr_diff(const struct rsb_mtx_t * mtxAp);
rsb_err_t rsb__mtx_shift_leaf_ptrs(struct rsb_mtx_t *RSB_RESTRICT  mtxCp, const struct rsb_mtx_t *RSB_RESTRICT  mtxAp, long smc);
void * rsb__clone_area_guided(void * RSB_RESTRICT dst, const void *RSB_RESTRICT src, size_t size, size_t nmemb, const struct rsb_mtx_t *RSB_RESTRICT mtxAp, const rsb_thread_t * RSB_RESTRICT cta, const rsb_thread_t nct, rsb_err_t * errvalp);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_CLONE_H_INCLUDED */
/* @endcond */
