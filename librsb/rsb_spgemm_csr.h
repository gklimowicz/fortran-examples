/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains some CSR sparse matrices multiplication code.
 * */

/*

Copyright (C) 2008-2022 Michele Martone

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
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

 */

#ifndef RSB_SPGEMM_COO_H_INCLUDED
#define RSB_SPGEMM_COO_H_INCLUDED
#include "rsb_internals.h"


rsb_err_t rsb__do_util_csr_csr_sparse_mul_serial(rsb_nnz_idx_t * PA, rsb_coo_idx_t * JA, void *VA_, const rsb_nnz_idx_t *ARP, const rsb_nnz_idx_t *BRP, const rsb_coo_idx_t *AJA, const rsb_coo_idx_t *BJA, const void * aVA_, const void * bVA_, const rsb_coo_idx_t cm, const rsb_coo_idx_t ck, rsb_nnz_idx_t * p, void * acc_, rsb_nnz_idx_t * opsp , rsb_type_t typecode, const rsb_coo_idx_t afr, const rsb_coo_idx_t ars)
;

rsb_err_t rsb__do_util_csr_csr_dense_mul_serial(rsb_coo_idx_t ldc, rsb_coo_idx_t nr, rsb_coo_idx_t nc, rsb_bool_t isccolmajor, void *cVA_, const rsb_nnz_idx_t *ARP, const rsb_nnz_idx_t *BRP, const rsb_coo_idx_t *AJA, const rsb_coo_idx_t *BJA, const void * aVA_, const void * bVA_, const rsb_coo_idx_t cm, const rsb_coo_idx_t ck, rsb_nnz_idx_t * opsp , rsb_type_t typecode, const rsb_coo_idx_t afr, const rsb_coo_idx_t ars)
;

#endif /* RSB_SPGEMM_COO_H_INCLUDED */
/* @endcond */
