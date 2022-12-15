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
 * This source file contains sorting functions.
 * */

#ifndef RSB_SRT_H_INCLUDED
#define RSB_SRT_H_INCLUDED

#include "rsb_internals.h"	/* rsb_coo_mtx_t */

rsb_err_t rsb__do_util_sortcoo(
	void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode,
	const struct rsb_mtx_partitioning_info_t * pinfop , rsb_flags_t flags, void * WA, size_t wb);

rsb_nnz_idx_t rsb__asymmetric_z_index( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl);
#if RSB_OBSOLETE_QUARANTINE
void rsb__asymmetric_z_nnz_indices( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl, rsb_nnz_idx_t * a , rsb_nnz_idx_t * b );
#endif /* RSB_OBSOLETE_QUARANTINE */

#if 0
rsb_err_t rsb__do_index_based_bcsr_msort( 
	rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, void * VA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	,void * WA, size_t wb
	);

rsb_err_t rsb__do_index_based_recursive_bcsr_sort( 
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	);

rsb_err_t rsb__do_nnz_index_based_sort_and_permute( 
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	);
#endif

void rsb__do_util_compact_permutation_nnz_idx_t_array(rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz);
#if 0
void rsb__do_util_compact_permutation_coo_idx_t_array(rsb_coo_idx_t * K, rsb_nnz_idx_t nnz);
rsb_err_t rsb__do_coo_index_sort_on_rows_array_make( 
	rsb_coo_idx_t * K, const rsb_coo_idx_t * IA,
	const rsb_coo_idx_t nr, const rsb_coo_idx_t br,
	const rsb_nnz_idx_t nnz, const rsb_type_t typecode);

rsb_err_t rsb__do_nnz_index_based_bcsr_msort( 
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	,void * WA, size_t wb);

rsb_err_t rsb__do_double_pass_nnz_index_based_bcsr_msort( 
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags);/* FIXME */

rsb_err_t rsb__do_nnz_index_sort_array_make( 
	rsb_nnz_idx_t * K, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t roffset,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags,
	int want_recursive_sort
	,enum rsb_op_flags_t op_flags
	/*, int want_rows_sort */);
#endif

#if 0
rsb_err_t rsb__do_double_coo_index_sort_array_make( 
	rsb_coo_idx_t * K, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t roffset,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags,
	int want_recursive_sort
	,enum rsb_op_flags_t op_flags
	/*, int want_rows_sort */);
#endif

rsb_nnz_idx_t rsb__nearest_power_of_two( const rsb_nnz_idx_t n );

#if RSB_OBSOLETE_QUARANTINE
void rsb__asymmetric_z_indices_encode( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl , rsb_coo_idx_t *h, rsb_coo_idx_t *l);
#endif /* RSB_OBSOLETE_QUARANTINE */
rsb_err_t rsb__do_index_based_z_morton_sort( 
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode
	,enum rsb_op_flags_t op_flags
	);

rsb_err_t rsb__do_index_based_bcsr_sort( 
	rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	,void * WA, size_t wb
	);

#define RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_ONE_PASS(NNZ,M,K,BR,BC) (((NNZ)+1) * sizeof(rsb_nnz_idx_t)  * 2)
#define RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_TWO_PASS(NNZ,M,K,BR,BC)  \
	RSB_MAX((((NNZ)+1) * sizeof(rsb_coo_idx_t)  * 2),(K)*(BR) * sizeof(rsb_nnz_idx_t) * 2)


#define RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT(NNZ,M,K,BR,BC) \
	RSB_MAX( \
		RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_ONE_PASS(NNZ,M,K,BR,BC), \
		RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_TWO_PASS(NNZ,M,K,BR,BC))
#endif /* RSB_SRT_H_INCLUDED */
/* @endcond */
