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
 * This source file contains functions for sparse matrix multiply.
 * FIXME : unfinished code.
 * */

#ifndef RSB_SPGEMM_H_INCLUDED
#define RSB_SPGEMM_H_INCLUDED
#include "rsb_internals.h"	/* rsb_coo_mtx_t */
rsb_err_t rsb__do_spgemm_test_code(const int argc, char * const argv[]);
struct rsb_mtx_t * rsb__do_matrix_mul(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_err_t * errvalp);
rsb_err_t rsb__do_spgemm_to_dense(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_coo_idx_t ldc, rsb_coo_idx_t nr, rsb_coo_idx_t nc, rsb_bool_t isccolmajor, void *cVA, rsb_time_t * dtp, size_t * opsp);
#endif /* RSB_SPGEMM_H_INCLUDED */
/* @endcond */
