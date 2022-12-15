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
 * This source file contains functions for COO handling.
 * */

#ifndef RSB_COO_H_INCLUDED
#define RSB_COO_H_INCLUDED

#include "rsb_internals.h"

rsb_nnz_idx_t rsb__weed_out_duplicates(rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, void *VA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags);
rsb_nnz_idx_t rsb__check_for_zeros(const void * RSB_RESTRICT VA, rsb_nnz_idx_t nnz, rsb_type_t typecode);
rsb_nnz_idx_t rsb__check_for_nonzeros(const void * RSB_RESTRICT VA, rsb_nnz_idx_t nnz, rsb_type_t typecode);
rsb_err_t rsb__util_compact_nonzeros(void *RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp, rsb_flags_t flags);
rsb_err_t rsb__util_compact_marked_coo_array( rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT VA, rsb_nnz_idx_t nnz, size_t el_size, rsb_coo_idx_t fd, rsb_nnz_idx_t * movedp, rsb_nnz_idx_t * movesp);
rsb_err_t rsb__weed_out_non_lowtri(void *RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp );
rsb_err_t rsb__weed_out_non_upptri(void *RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp );
rsb_err_t rsb__weed_out_diagonal(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp );
rsb_err_t rsb__do_util_compact_out_of_range(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t roff, rsb_coo_idx_t  coff, rsb_coo_idx_t Mdim, rsb_coo_idx_t mdim, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp );

#endif /* RSB_COO_H_INCLUDED */
/* @endcond */
