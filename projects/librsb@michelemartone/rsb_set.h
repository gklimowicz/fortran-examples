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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * Matrix setter/getter functions.
 * */

#ifndef RSB_SET_H_INCLUDED
#define RSB_SET_H_INCLUDED

#include "rsb_internals.h"

rsb_err_t rsb__do_set_coo_element(struct rsb_mtx_t * mtxAp, const void * vp, const rsb_coo_idx_t i, const rsb_coo_idx_t j);
rsb_err_t rsb__do_upd_coo_element(struct rsb_mtx_t * mtxAp, const void * vp, const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_flags_t flags);
rsb_err_t rsb__do_get_coo_element(const struct rsb_mtx_t * mtxAp, void * vp, rsb_coo_idx_t i, rsb_coo_idx_t j);
rsb_err_t rsb__do_set_coo_elements(struct rsb_mtx_t * mtxAp, const void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz);
const void * rsb__do_coo_element_inner_address(const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t i, rsb_coo_idx_t j);
rsb_err_t rsb__do_reverse_odd_rows(struct rsb_mtx_t * mtxAp);
rsb_err_t rsb__do_zsort_coo_submatrices(struct rsb_mtx_t * mtxAp);
#endif /* RSB_SET_H_INCLUDED */
rsb_err_t rsb__do_get_nnz_element(const struct rsb_mtx_t * mtxAp, void * vp, rsb_coo_idx_t*ip, rsb_coo_idx_t*jp, rsb_nnz_idx_t nzi);
/* @endcond */
