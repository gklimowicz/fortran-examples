/* @cond INNERDOC */


/**
 * @file
 * @brief
 * Auxiliary functions.
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


#ifndef RSB_UTIL_H_INCLUDED
#define RSB_UTIL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb_common.h"
/* non blas-like functions */

rsb_err_t rsb__util_m4_sanity_check(void);
const void * rsb__util_increase_by_one(void *p, rsb_nnz_idx_t n, rsb_type_t typecode);
void rsb__util_set_area_to_fraction_of_integer(void *p, const int alphai, rsb_type_t typecode);
void rsb__util_set_area_to_negated_fraction(void *p, const void *alpha, rsb_type_t typecode);
void rsb__util_set_area_to_converted_integer(void *p, rsb_type_t typecode, const rsb_int n);
rsb_coo_idx_t * rsb__util_get_partitioning_array( size_t bs, size_t X , rsb_blk_idx_t * X_b, rsb_flags_t flags);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__vector_diff(void * c, const void * a, const void * b, rsb_type_t type, size_t n);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */




rsb_err_t rsb__vector_norm_strided(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
;
rsb_err_t rsb__util_vector_sum_strided(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
;
rsb_err_t rsb__util_vector_sum(void * c, const void * a, rsb_type_t type, size_t n)
;
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */





rsb_err_t rsb__util_vector_add(void * a, const void * alphap, rsb_type_t type, size_t n)
;
rsb_err_t rsb__util_vector_div(void * a, const void * alphap, rsb_type_t type, size_t n)
;
rsb_err_t rsb__vector_increase_by_one(void * a, rsb_type_t type, size_t n)
;
rsb_err_t rsb__util_vector_pow(void * a, rsb_type_t type, const void *y, size_t n)
;
rsb_err_t rsb__util_vector_sqrt(void * a, rsb_type_t type, size_t n)
;
rsb_err_t rsb__vector_scale_inv(void * a, const void * alphap, rsb_type_t type, size_t n)
;
rsb_err_t rsb__vector_sum_of_abs_diffs(void * c, const void * a, const void * b, rsb_type_t type, size_t n)
;
rsb_err_t rsb__vector_sum_of_abs(void * c, const void * a, rsb_type_t type, size_t n)
;
rsb_err_t rsb__vector_to_abs(void * a, rsb_type_t type, size_t n)
;


rsb_err_t rsb__util_set_array_to_converted_integer(void *p, rsb_type_t typecode, const rsb_nnz_idx_t n, const rsb_nnz_idx_t incp, const rsb_int v)
;
rsb_err_t rsb__vectors_left_sum_reduce_and_zero(void * d, void * s, const rsb_type_t typecode, const size_t n, const size_t incd, const size_t off)
;


rsb_err_t rsb__cblas_Xaxpy(rsb_type_t type, size_t n, const void * alphap, const void * x, const int incx, void * y, const int incy)
;
rsb_err_t rsb__vector_mult(const void * a, const void * b, void * c, rsb_type_t type, size_t n)
;
rsb_err_t rsb__xcopy(void * a, const void * b, rsb_nnz_idx_t toi, rsb_nnz_idx_t foi, rsb_nnz_idx_t n,size_t el_size)
;
rsb_err_t rsb__do_are_similar_parametric(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy, int extra_decimals)
;
rsb_err_t rsb__do_are_similar(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
;
rsb_err_t rsb__do_are_same(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
;

rsb_err_t rsb__sqrt_of_sum_of_fabs_diffs(const void * a, const void * b, void *err, rsb_type_t type, size_t n)
;
rsb_err_t rsb__fill_with_increasing_values(void * array, rsb_type_t type, size_t n)
;
rsb_err_t rsb__util_do_conjugate(void * array, rsb_type_t type, size_t n)
;
rsb_err_t rsb__util_do_negate(void * array, rsb_type_t type, size_t n)
;
rsb_err_t rsb__util_find_min(void * minp, const void * array, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
;
rsb_err_t rsb__util_find_max(void * maxp, const void * array, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
;
rsb_err_t rsb__util_drop_to_zero_if_above_threshold(void * array, rsb_type_t type, size_t n, const void * threshold)
;
rsb_nnz_idx_t rsb__util_count_positive(void * array, rsb_type_t type, size_t n)
;
rsb_nnz_idx_t rsb__util_count_negative(void * array, rsb_type_t type, size_t n)
;
rsb_err_t rsb__util_drop_to_zero_if_under_threshold(void * array, rsb_type_t type, size_t n, const void * threshold)
;
rsb_err_t rsb__fill_with_ones(void * array, rsb_type_t type, size_t n, size_t incx);
rsb_err_t rsb__debug_print_vectors_diff_fd(const void * v1, const void * v2, size_t n, rsb_type_t type, size_t incx, size_t incy, int onlyfirst, FILE*fd);
rsb_err_t rsb__debug_print_vectors_diff(const void * v1, const void * v2, size_t n, rsb_type_t type, size_t incx, size_t incy, int onlyfirst);
rsb_err_t rsb__debug_print_value(const void * v, rsb_type_t type);
rsb_err_t rsb__debug_print_vector_extra(const void * v1, size_t n, rsb_type_t type, size_t inc, int style, FILE*stream);
rsb_err_t rsb__debug_print_vector(const void * v1, size_t n, rsb_type_t type, size_t inc);
rsb_err_t rsb__debug_print_vectors(const void * v1, const void * v2, size_t n, size_t incx, size_t incy, rsb_type_t type);

#ifdef RSB_WANT_OSKI_BENCHMARKING 
rsb_err_t rsb__do_account_sorted_optimized_css(
	 const rsb_coo_idx_t * MIndx, const rsb_coo_idx_t * mIndx,
	 const rsb_coo_idx_t Mdim, const rsb_coo_idx_t mdim,
	 const rsb_nnz_idx_t nnz, rsb_nnz_idx_t * elements_per_block_row, rsb_nnz_idx_t * blocks_per_block_row
)
;
#endif /* RSB_WANT_OSKI_BENCHMARKING */

#if RSB_OBSOLETE_QUARANTINE_UNUSED

rsb_err_t rsb__do_account_sorted_optimized(
	 struct rsb_mtx_t * mtxAp,
	 const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	 const rsb_coo_idx_t Idim, const rsb_coo_idx_t Jdim,
	 const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop,
rsb_nnz_idx_t * elements_per_block_row, 
rsb_nnz_idx_t * blocks_per_block_row
)
;
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_err_t rsb__do_insert_sorted_optimized_css( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * MIndx, const rsb_coo_idx_t * mIndx, const rsb_nnz_idx_t nnz)
;
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_insert_sorted_optimized( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop)
;
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_block(rsb_type_t type, const void * VA, rsb_blk_idx_t roff, rsb_blk_idx_t coff, rsb_blk_idx_t rows, rsb_blk_idx_t cols )
;
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_blocks(const struct rsb_mtx_t *mtxAp)
;
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__test_print_csr(rsb_type_t type, rsb_flags_t flags, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t want_header, FILE*stream)
;
rsb_err_t rsb__test_print_coo_mm(rsb_type_t type, rsb_flags_t flags, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t want_header, FILE*stream)
;
/*static*/ /*inline*/ size_t rsb__do_sizeof(rsb_type_t type);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_coo_sum( struct rsb_coo_mtx_t*coocp, const void *alphap, const struct rsb_coo_mtx_t*cooap, const void *betap,  const struct rsb_coo_mtx_t*coobp)
;
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__cor_merge_dups(rsb_type_t typecode, void* RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t offB, rsb_nnz_idx_t nnzB, rsb_nnz_idx_t nnzC, const int wv, int wp, rsb_nnz_idx_t *onzp, struct rsb_coo_mtx_t*RSB_RESTRICT coop)
;
rsb_err_t rsb__do_copy_converted_scaled(const void *RSB_RESTRICT  src, void *RSB_RESTRICT dst, const void *RSB_RESTRICT  alphap, rsb_type_t stype,rsb_type_t dtype, size_t nnz, rsb_trans_t transA)
;
rsb_err_t rsb__util_csc2csr(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t*flagsp)
;
rsb_err_t rsb__util_coo_copy_and_stats(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, rsb_coo_idx_t*m, rsb_coo_idx_t*k, const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t iflags, rsb_flags_t*flagsp)
;
rsb_err_t rsb__util_coo_copy(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo)
;
/* sparse blas level 1 equivalent functions */

int rsb__BLAS_Xusdot(const rsb_type_t typecode, const enum blas_conj_type conj_arg, const rsb_blas_int_t nz, const void*x, const rsb_blas_int_t*indx, const void*y, const rsb_blas_int_t incy, void*r, const enum blas_base_type index_base)
;
int rsb__BLAS_Xusaxpy(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*alpha, const void*x, const rsb_blas_int_t*indx, const void*y, const rsb_blas_int_t incy, const enum blas_base_type index_base)
;
int rsb__BLAS_Xusga(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*y, const rsb_blas_int_t incy, void*x, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
;
int rsb__BLAS_Xusgz(const rsb_type_t typecode, const rsb_blas_int_t nz, void*y, const rsb_blas_int_t incy, void*x, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
;
int rsb__BLAS_Xussc(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*x, void*y, const rsb_blas_int_t incy, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
;
/* blas level 1 equivalent functions */

rsb_err_t rsb__cblas_Xcopy(rsb_type_t typecode, rsb_nnz_idx_t n, const void * x, rsb_nnz_idx_t incx, void * y, rsb_nnz_idx_t incy)
;
rsb_err_t rsb__cblas_Xnrm2(rsb_type_t type, size_t n, const void * a, rsb_nnz_idx_t incA, void * c);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__cblas_Xdotu_sub(rsb_type_t type, size_t n, const void * x, rsb_nnz_idx_t incx, const void * y, rsb_nnz_idx_t incy, void *dotu);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__cblas_Xscal(rsb_type_t type, size_t n, const void * alphap, void * a, size_t stride);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__coo_insertion_sort(rsb_type_t typecode, void* VB, rsb_coo_idx_t * IB, rsb_coo_idx_t * JB, rsb_nnz_idx_t offA, rsb_nnz_idx_t nnzA)
;
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void rsb__coo_to_lr( void * RSB_RESTRICT VBu, rsb_coo_idx_t*RSB_RESTRICT IB, rsb_coo_idx_t*RSB_RESTRICT JB, void * RSB_RESTRICT VAu, rsb_coo_idx_t*RSB_RESTRICT IA, rsb_coo_idx_t*RSB_RESTRICT JA, rsb_coo_idx_t mj, rsb_nnz_idx_t nnzA, rsb_nnz_idx_t nzoffB, rsb_nnz_idx_t nzoffA, rsb_nnz_idx_t*RSB_RESTRICT nzlp, rsb_nnz_idx_t*RSB_RESTRICT nzrp, rsb_coo_idx_t iadd, rsb_coo_idx_t jadd, rsb_type_t typecode)
;
rsb_err_t rsb__util_testall(void)
;
#define RSB__MAX_SIZEOF	( \
	RSB_MAX(sizeof(float),\
	RSB_MAX(sizeof(double),\
	RSB_MAX(sizeof(float complex),\
	RSB_MAX(sizeof(double complex),\
	0)))))

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_UTIL_H_INCLUDED */

/* @endcond */
