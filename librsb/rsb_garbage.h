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
 * This source file contains functions for benchmarking, integration testing, and miscellaneous.
 * */

#ifndef RSB_GARBAGE_H_INCLUDED
#define RSB_GARBAGE_H_INCLUDED

#include "rsb_internals.h"	/* rsb_coo_mtx_t */
#include "rsb_common.h"

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__test_bitmap_driver(rsb_coo_idx_t r, rsb_coo_idx_t c);
int rsb__test_fill_matrix_nnz(rsb_type_t typecode, rsb_nnz_idx_t nnz, void *VA );
rsb_err_t rsb__test_fill_matrix_coords(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t allow_duplicates);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
int rsb__test_dump_main(const int argc,rsb_char_t *const argv[]);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__test_gen_and_print_matrix(rsb_type_t typecode, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
int rsb_test_main_block_partitioned_matrix_stats(int argc,rsb_char_t *argv[]);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_blk_idx_t rsb__rand_blk_index(rsb_blk_idx_t max_plus_one);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_coo_idx_t rsb__rand_coo_index(rsb_coo_idx_t max_plus_one);
rsb_flags_t rsb__sample_program_options_get_flags(int c, const rsb_char_t * optarg);
int rsb__dump_postscript(const int argc, rsb_char_t * const argv[]);
rsb_err_t rsb__oski_estimate_bcsr_fill_from_coo(/*  const*/ rsb_coo_idx_t * IA, /*const*/ rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_fillin_t * efillinmap );
rsb_err_t rsb__do_column_expand(rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t * kp, rsb_int factor);
rsb_err_t rsb__do_print_some_vector_stats(const void * p, rsb_type_t typecode, rsb_nnz_idx_t m, rsb_nnz_idx_t inc);

#define RSB_SAMPLE_PROGRAM_OPTIONS_GET_FLAGS "q:QLECHDVARisF:PT:"
#define RSB_FLAG_DEFAULT_STORAGE RSB_FLAG_WANT_BCSS_STORAGE

#endif /* RSB_GARBAGE_H_INCLUDED */
/* @endcond */
