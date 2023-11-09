/* @cond INNERDOC */
/**
 * @file
 * @brief
 * Permutation functions.
 */

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



#ifndef RSB_PERMUTE_H_INCLUDED
#define RSB_PERMUTE_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "rsb_common.h"

rsb_err_t rsb__do_permute_values_in_place_with_coo_index(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type);
rsb_err_t rsb__do_permute_values_in_place_with_nnz_index(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type);
rsb_err_t rsb__do_permute_values_with_coo_index( void * rVA, const void *VA, rsb_coo_idx_t * rIA, const rsb_coo_idx_t * IA, rsb_coo_idx_t * rJA, const rsb_coo_idx_t * JA, const rsb_coo_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type);

rsb_err_t rsb__do_permute_rows_with_coo_index( rsb_coo_idx_t * IA, const rsb_coo_idx_t * K, rsb_nnz_idx_t nnz);

rsb_err_t rsb__do_permute_values_with_nnz_index( void * rVA, const void *VA, rsb_coo_idx_t * rIA, const rsb_coo_idx_t * IA, rsb_coo_idx_t * rJA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t typecode);

void rsb__ip_reord(rsb_nnz_idx_t n, void * VAp, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * P, rsb_type_t typecode);
void rsb__util_do_scatter_rows(void * RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, const void * RSB_RESTRICT iVA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT PA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode);
rsb_err_t rsb__chk_permute(void);
#if 0
void rsb__util_do_scatter_rows(void * RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT PA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode);
#endif /* 0 */



#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_PERMUTE_H_INCLUDED */

/* @endcond */
