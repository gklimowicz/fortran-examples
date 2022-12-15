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
 * This source file contains searching functions.
 * */

#ifndef RSB_SEARCH_H_INCLUDED
#define RSB_SEARCH_H_INCLUDED

#include "rsb_swt.h"	/* rsb_half_idx_t */
#include "rsb_common.h"	/* rsb_coo_mtx_t */

rsb_nnz_idx_t rsb__nnz_split_nnz_bsearch(const rsb_nnz_idx_t*A,const rsb_nnz_idx_t S,const rsb_nnz_idx_t n);
rsb_nnz_idx_t rsb__nnz_split_coo_bsearch(const rsb_coo_idx_t*A,const rsb_coo_idx_t S,const rsb_nnz_idx_t n);
rsb_nnz_idx_t rsb__nnz_split_hcoo_bsearch(const rsb_half_idx_t *A, const rsb_half_idx_t S,const rsb_nnz_idx_t n);
rsb_nnz_idx_t rsb__seek_nnz_idx_t(const rsb_nnz_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n);
rsb_nnz_idx_t rsb__seek_coo_idx_t(const rsb_coo_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n);
rsb_nnz_idx_t rsb__seek_nnz_idx_t_linear(const rsb_nnz_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n);
rsb_nnz_idx_t rsb__seek_half_idx_t(const rsb_half_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n);
#endif /* RSB_SEARCH_H_INCLUDED */
/* @endcond */
