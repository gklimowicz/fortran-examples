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
 * This source file contains some functions testing accuracy.
 * */

#ifndef RSB_TEST_ACCURACY_H_INCLUDED
#define RSB_TEST_ACCURACY_H_INCLUDED

#include "rsb_internals.h"
rsb_err_t rsb__vectors_reinit(void *rhs, void *out, rsb_type_t typecode, rsb_nnz_idx_t rn, rsb_nnz_idx_t on, size_t incr, size_t inco); 
void * rsb__calloc_vector(rsb_nnz_idx_t n, rsb_type_t typecode);
void * rsb__malloc_vector(rsb_nnz_idx_t n, rsb_type_t typecode);
void * rsb__realloc_vector(void* p, rsb_nnz_idx_t n, rsb_type_t typecode);
rsb_err_t rsb__do_spmv_accuracy_test(const struct rsb_coo_mtx_t*coop, rsb_thread_t * ca, rsb_thread_t cn, rsb_flags_t flags);
rsb_err_t rsb__do_spmv_fullword_coo(const struct rsb_coo_mtx_t*coop, rsb_flags_t flags, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA);
rsb_err_t rsb__init_rsb_struct_from_coo(struct rsb_mtx_t *mtxAp, const struct rsb_coo_mtx_t *coop);
#endif /* RSB_TEST_ACCURACY_H_INCLUDED */
/* @endcond */
