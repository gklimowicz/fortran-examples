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
 * @brief Sparse matrices assembling code.
 * */

#ifndef RSB_ASM_H_INCLUDED
#define RSB_ASM_H_INCLUDED
#include "rsb_common.h"
struct rsb_mtx_t * rsb__mtx_alloc_inner(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_type_t typecode, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_blk_idx_t Mb, rsb_blk_idx_t Kb, rsb_flags_t flags, rsb_err_t * errvalp);
#endif /* RSB_ASM_H_INCLUDED */
/* @endcond */
