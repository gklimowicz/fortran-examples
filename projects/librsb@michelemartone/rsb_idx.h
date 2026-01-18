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
 * This source file contains functions for COO handling.
 * */

#ifndef RSB_IDX_H_INCLUDED
#define RSB_IDX_H_INCLUDED

#include "rsb_internals.h"

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_nnz_array_set_sequence(rsb_nnz_idx_t * p, rsb_nnz_idx_t n, rsb_nnz_idx_t o, rsb_nnz_idx_t i);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void rsb__util_coo_array_set_sequence(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t o, rsb_coo_idx_t i);
void rsb__util_coo_array_set(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a);
void rsb__util_nnz_array_set(rsb_nnz_idx_t * p, rsb_nnz_idx_t n, rsb_nnz_idx_t a);
//void rsb__util_coo_array_mul(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a);
void rsb__util_coo_arrays_mul(rsb_coo_idx_t * p, rsb_coo_idx_t * q, rsb_coo_idx_t a, rsb_coo_idx_t b, rsb_nnz_idx_t n);
void rsb__util_coo_array_add(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a);
void rsb__util_hcoo_array_add(rsb_half_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_coo_arrays_add(rsb_coo_idx_t * p, rsb_coo_idx_t * q, rsb_coo_idx_t a, rsb_coo_idx_t b, rsb_nnz_idx_t n);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void rsb_util_coo_arrays_sub(rsb_coo_idx_t * p, rsb_coo_idx_t * q, rsb_coo_idx_t a, rsb_coo_idx_t b, rsb_nnz_idx_t n);
#if RSB_WANT_BCSC_LEAVES  
void rsb__util_nnz_array_add_array(rsb_nnz_idx_t * p, const rsb_nnz_idx_t * q, rsb_nnz_idx_t n);
#endif /* RSB_WANT_BCSC_LEAVES */
void rsb__util_coo_array_sub(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t s);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_coo_array_to_fortran_indices(rsb_coo_idx_t * p, rsb_nnz_idx_t n);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void rsb__util_coo_array_to_fortran_indices_parallel(rsb_coo_idx_t * p, rsb_nnz_idx_t n);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_coo_array_from_fortran_indices(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_bool_t want_parallel);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void rsb__util_coo_upper_to_lower_symmetric(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_coo_lower_to_upper_symmetric(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz);
void rsb__util_nnz_array_from_fortran_indices(rsb_coo_idx_t * p, rsb_nnz_idx_t n);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void rsb__util_nnz_array_to_fortran_indices(rsb_coo_idx_t * p, rsb_nnz_idx_t n);
rsb_bool_t rsb__util_coo_check_if_triangle_non_empty(const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_flags_t flags);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_flags_t rsb__util_coo_determine_uplo_flags(const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t nnz);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_err_t rsb__util_coo_check_if_has_diagonal_elements(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_bool_t *has_diagonal_elements, rsb_bool_t wv);
void rsb__util_coo_array_copy_trans_add(rsb_coo_idx_t * d, const rsb_coo_idx_t * s, rsb_nnz_idx_t n, rsb_coo_idx_t a);
void rsb__util_hcoo_array_copy_trans_add(rsb_coo_idx_t * d, const rsb_half_idx_t * s, rsb_nnz_idx_t n, rsb_coo_idx_t a);
rsb_bool_t rsb__util_reverse_halfword_coo_array(rsb_half_idx_t* p, rsb_nnz_idx_t n);
rsb_bool_t rsb__util_reverse_fullword_coo_array(rsb_coo_idx_t* p, rsb_nnz_idx_t n);
rsb_bool_t rsb__util_is_coo_array_sorted_up_partial_order(const rsb_coo_idx_t * p, const rsb_nnz_idx_t n);
rsb_bool_t rsb__util_is_halfword_coo_array_sorted_up_partial_order(const rsb_half_idx_t * p, const rsb_nnz_idx_t n);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__util_is_coo_array_sorted_up(const rsb_coo_idx_t * p, const rsb_nnz_idx_t n);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/
rsb_bool_t rsb__util_is_halfword_coo_array_sorted_up(const rsb_half_idx_t* p, const rsb_nnz_idx_t n);
rsb_bool_t rsb__util_is_nnz_array_sorted_up(const rsb_nnz_idx_t * p, const rsb_nnz_idx_t n);
rsb_int_t rsb__util_find_max_index_val(const rsb_int_t * p, rsb_int_t n);
rsb_int_t rsb__util_find_min_index_val(const rsb_int_t * p, rsb_int_t n);
rsb_nnz_idx_t rsb__util_find_coo_max_index_val(const rsb_nnz_idx_t * p, rsb_nnz_idx_t n);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_nnz_idx_t rsb__util_find_coo_min_index_val(const rsb_nnz_idx_t * p, rsb_nnz_idx_t n);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_coo_idx_t rsb__util_find_coo_val_idx(const rsb_coo_idx_t * p, const rsb_nnz_idx_t n, const rsb_coo_idx_t m);
rsb_err_t rsb__util_find_extremal_half_index_val(const rsb_half_idx_t * RSB_RESTRICT p, rsb_nnz_idx_t n, rsb_coo_idx_t lb, rsb_coo_idx_t ub, rsb_half_idx_t *lf, rsb_half_idx_t * RSB_RESTRICT uf);
rsb_err_t rsb__util_find_extremal_full_index_val(const rsb_coo_idx_t * RSB_RESTRICT p, rsb_nnz_idx_t n, rsb_coo_idx_t lb, rsb_coo_idx_t ub, rsb_coo_idx_t * RSB_RESTRICT lf, rsb_coo_idx_t * RSB_RESTRICT uf);
void rsb__util_coo_array_renumber(rsb_coo_idx_t * a, const rsb_coo_idx_t * iren, rsb_nnz_idx_t n, rsb_flags_t aflags, rsb_flags_t pflags, rsb_flags_t oflags);
rsb_err_t rsb__util_uncompress_row_pointers_array(const rsb_coo_idx_t * pa, rsb_nnz_idx_t n, rsb_flags_t iflags, rsb_flags_t oflags, rsb_coo_idx_t * ta);
rsb_err_t rsb__util_compress_to_row_pointers_array(rsb_coo_idx_t * RSB_RESTRICT pa, rsb_nnz_idx_t nz, rsb_coo_idx_t m, rsb_flags_t iflags, rsb_flags_t oflags, rsb_coo_idx_t * ta);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__debug_print_index_vector(const rsb_coo_idx_t * v1, size_t n);
rsb_err_t rsb__debug_print_index_vectors_diff(const rsb_coo_idx_t * v1, const rsb_coo_idx_t * v2, size_t n, int onlyfirst);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_bool_t rsb__util_is_nnz_array_sorted_up_partial_order(const rsb_nnz_idx_t * p, const rsb_nnz_idx_t n);
rsb_int_t rsb__util_find_max_index(const rsb_int_t * p, rsb_int_t n);
rsb_int_t rsb__util_find_min_index(const rsb_int_t * p, rsb_int_t n);
#define RSB_IA_MEMCPY(ID,IS,DOFF,SOFF,NNZ,I0) \
	rsb__util_coo_array_copy_trans_add(((rsb_coo_idx_t*)(ID))+(DOFF),((rsb_coo_idx_t*)(IS))+(SOFF),NNZ,I0)
#define RSB_IA_MEMCPY_H(ID,IS,DOFF,SOFF,NNZ,I0) \
	rsb__util_hcoo_array_copy_trans_add(((rsb_coo_idx_t*)(ID))+(DOFF),((rsb_half_idx_t*)(IS))+(SOFF),NNZ,I0)

#endif /* RSB_IDX_H_INCLUDED */
/* @endcond */
