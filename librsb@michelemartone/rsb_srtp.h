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
 * This source file contains parallel sorting functions.
 * */

#ifndef RSB_PSORT_H_INCLUDED
#define RSB_PSORT_H_INCLUDED

#include "rsb_common.h"

rsb_err_t rsb__util_sort_row_major_inner(void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const  rsb_type_t typecode , const rsb_flags_t flags /*, void * WA, size_t wb */);
rsb_err_t rsb__util_sort_row_major_parallel(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t nr, rsb_coo_idx_t nc,  rsb_type_t typecode, rsb_flags_t flags);
rsb_err_t rsb__util_sort_row_major_bucket_based_parallel(void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const  rsb_type_t typecode , const rsb_flags_t flags /*, void * WA, size_t wb */);

#endif /* RSB_PSORT_H_INCLUDED */
/* @endcond */
