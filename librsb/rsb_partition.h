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
/**
 * @file
 * @author Michele Martone
 * @brief
 * Auxiliary functionalities header file.
 */
#ifndef RSB_PARTITION_H_INCLUDED
#define RSB_PARTITION_H_INCLUDED

#include "rsb_internals.h"

int rsb__util_nnz2vbr(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const rsb_blk_idx_t rows, const rsb_blk_idx_t columns, rsb_blk_idx_t **p_rp, rsb_blk_idx_t **p_cp, rsb_blk_idx_t *M_b, rsb_blk_idx_t *K_b, int blockrowsize, int blockcolumnsize);
rsb_bool_t rsb__should_rejoin_small_leaf(rsb_nnz_idx_t nnz, rsb_nnz_idx_t mk, rsb_nnz_idx_t uk, rsb_nnz_idx_t lk, rsb_type_t typecode);

#endif /* RSB_PARTITION_H_INCLUDED */
/* @endcond */
