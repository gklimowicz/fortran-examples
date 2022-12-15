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
/*
 * @author Michele Martone
 */
#ifndef RSB_INTERNALS_H_INCLUDED
#define RSB_INTERNALS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb_common.h"

/**
 * @file
 * @brief
 * Low level routines and tools for our sparse matrix formats implementations.
 *
 */
rsb_bool_t rsb__are_coo_matrices_equal(const struct rsb_coo_mtx_t *cm1, const struct rsb_coo_mtx_t *cm2);
rsb_bool_t rsb__are_coo_matrices_both_empty(const struct rsb_coo_mtx_t *cm1, rsb_flags_t flags1, const struct rsb_coo_mtx_t *cm2, rsb_flags_t flags2);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__are_coo_matrices_equal_or_both_empty(const struct rsb_coo_mtx_t *cm1, rsb_flags_t flags1, const struct rsb_coo_mtx_t *cm2, rsb_flags_t flags2);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void * rsb__destroy_coo_matrix_t(struct rsb_coo_mtx_t *cmp);
void * rsb__allocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp);
void * rsb__xallocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp, rsb_bool_t want_calloc, rsb_flags_t flags);
void * rsb__callocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp);
void * rsb__reallocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp, rsb_nnz_idx_t nnnz);


/*
 * Please note that the sole use of this function is the major bottleneck during matrix creation.
 * When thinking about optimizing matrix creation, come back here: this routine eats up to 90% 
 * of the time required for matrix creation.
 * */
int rsb__nnz_coord_compar(const void *key, const void *am);


/* initialization, destroying */
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void * rsb__init_options_t(struct rsb_options_t *o);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void * rsb__init_struct(struct rsb_mtx_t *mtxAp);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__fill_struct(struct rsb_mtx_t *mtxAp, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_type_t typecode, rsb_flags_t flags);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
void * rsb__fill_coo_struct(struct rsb_coo_mtx_t *mtxAp, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode);
void * rsb__init_blank_pointers(struct rsb_mtx_t *mtxAp);
void * rsb__transpose_coo_matrix_t(struct rsb_coo_mtx_t *cmp);
void * rsb__do_mtx_free(struct rsb_mtx_t *mtxAp);
size_t rsb__get_sizeof(const struct rsb_mtx_t *mtxAp );
void * rsb__destroy_inner(struct rsb_mtx_t *mtxAp);
#if RSB_WANT_BITMAP
void * rsb__destroy_options_t(struct rsb_options_t *o);
#endif /* RSB_WANT_BITMAP */
rsb_err_t rsb__set_init_flags_and_stuff( struct rsb_mtx_t *mtxAp, struct rsb_options_t * o, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_nnz_idx_t block_count, rsb_nnz_idx_t element_count, rsb_type_t typecode, rsb_flags_t flags );
rsb_err_t rsb__do_set_init_storage_flags(struct rsb_mtx_t *mtxAp, rsb_flags_t flags);

/* allocation */
rsb_bitmap_data_t * rsb__allocate_bitmap(rsb_blk_idx_t rows, rsb_blk_idx_t cols);
rsb_bitmap_data_t * rsb__allocate_bitvector(rsb_blk_idx_t bits);
struct rsb_mtx_t * rsb__allocate_from_coo_sorted(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t m, rsb_coo_idx_t k, struct rsb_options_t * o, rsb_type_t typecode, rsb_flags_t flags, rsb_err_t *errvalp);
struct rsb_mtx_t * rsb__allocate_css_from_coo_sorted(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t m, rsb_coo_idx_t k, struct rsb_options_t * o, rsb_type_t typecode, rsb_flags_t flags, rsb_err_t *errvalp);
#ifdef RSB_WANT_OSKI_BENCHMARKING 
rsb_err_t rsb__allocate_csr_arrays_from_coo_sorted(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_type_t typecode, void **VAp, rsb_coo_idx_t ** indxp, rsb_nnz_idx_t ** indptrp);
rsb_err_t rsb__allocate_csc_arrays_from_coo_sorted(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_type_t typecode, void **VAp, rsb_coo_idx_t ** indxp, rsb_nnz_idx_t ** indptrp);
#endif /* RSB_WANT_OSKI_BENCHMARKING */

/* bit handling */
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_blk_idx_t rsb__bitmap_bit_count(const rsb_bitmap_data_t *bitmap, const rsb_blk_idx_t rows, const rsb_blk_idx_t cols);

/* check */
rsb_err_t rsb__recheck_insertion(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, const struct rsb_mtx_t *mtxAp, const struct rsb_options_t *o);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
const void * rsb__is_valid_options_t(const struct rsb_options_t *o, rsb_coo_idx_t m, rsb_coo_idx_t k);

/* misc */
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void* rsb__get_block_address( rsb_blk_idx_t blockrow, rsb_blk_idx_t blockcolumn, const struct rsb_mtx_t *mtxAp);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

size_t rsb__get_g_rsb_memory_count(void);/* rsb_sys.c */

rsb_err_t rsb__compute_partial_fillin_for_nnz_fractions(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,const  rsb_nnz_idx_t * nnz, const rsb_nnz_idx_t nnzn, struct rsb_mtx_partitioning_info_t * pinfop, size_t * element_countp, size_t * block_countp);
#if 0
rsb_err_t rsb__compute_partial_fillin_for_nnz_fraction(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,const  rsb_nnz_idx_t nnz, struct rsb_mtx_partitioning_info_t * pinfop, size_t * element_countp, size_t * block_countp);
#endif

#define RSB__USE_MTX_PARTITIONING_INFO_T 0 /* rsb_mtx_partitioning_info_t is obsolete */

rsb_err_t rsb__fprint_matrix_implementation_code(const struct rsb_mtx_t *mtxAp, const rsb_char_t * op, rsb_flags_t inflags, FILE*fd);
const rsb_char_t * rsb__sprint_matrix_implementation_code2(const struct rsb_mtx_t *mtxAp, rsb_char_t * buf, rsb_flags_t inflags);
rsb_err_t rsb__util_get_tn_array(const rsb_char_t* optarg, int* bxlp, rsb_blk_idx_t **bxvp);
rsb_err_t rsb__util_get_bx_array_or_default(const rsb_char_t defsym, const rsb_char_t* defarg, const rsb_char_t* optarg, int* bxlp, rsb_blk_idx_t **bxvp);
rsb_err_t rsb__util_get_bx_array(const rsb_char_t* optarg, int* bxlp, rsb_blk_idx_t **bxvp);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_nnz_idx_t rsb__util_atonnz(const rsb_char_t * optarg);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
//rsb_long_t rsb__util_atol(const rsb_char_t *nptr);
rsb_real_t rsb__util_atof(const rsb_char_t *nptr);
int  rsb__util_atoi(const rsb_char_t *nptr);
const rsb_char_t *rsb__basename(const rsb_char_t *path);
size_t rsb__util_strlen(const rsb_char_t *s);
#if RSB__USE_MTX_PARTITIONING_INFO_T
rsb_err_t rsb__do_is_valid_pinfo_t(const struct rsb_mtx_partitioning_info_t * pinfop);
#endif /* RSB__USE_MTX_PARTITIONING_INFO_T */
rsb_err_t rsb__print_configuration_string(const char *pn, rsb_char_t * cs, rsb_bool_t wci);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_blk_idx_t rsb__recursive_middle_block_index(rsb_blk_idx_t i);
rsb_err_t rsb__recursive_split_point_parms_get(const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t * moff, rsb_coo_idx_t * koff);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
#if RSB__USE_MTX_PARTITIONING_INFO_T
rsb_err_t rsb__do_get_blocking_from_pinfo(const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_blk_idx_t *mbp, rsb_blk_idx_t *kbp);
#endif

/* fill */
rsb_err_t rsb__do_insert_sorted( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop);
rsb_err_t rsb__do_account_sorted( struct rsb_mtx_t * mtxAp, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_nnz_idx_t * elements_per_block_row, rsb_nnz_idx_t * blocks_per_block_row);
rsb_long_t rsb__terminal_recursive_matrix_count(const struct rsb_mtx_t *mtxAp);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__copy_css_arrays(const void *iVA, const rsb_coo_idx_t * iINDX, const rsb_coo_idx_t * iXA, const rsb_nnz_idx_t nnz, rsb_coo_idx_t X, rsb_type_t typecode, void *oVA, rsb_coo_idx_t * oINDX, rsb_nnz_idx_t * oXA);
rsb_long_t rsb__terminal_recursive_matrix_count_with_flags(const struct rsb_mtx_t *mtxAp, rsb_flags_t flags);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
// rsb_long_t rsb__terminal_recursive_matrix_count_with_flags_but(const struct rsb_mtx_t *mtxAp, rsb_flags_t flags, rsb_flags_t nflags);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__recursive_middle_index(const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t * M_bp, rsb_coo_idx_t * K_bp );
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_trans_t rsb__do_transpose_transposition(rsb_trans_t transA);
struct rsb_mtx_t * rsb__do_get_first_submatrix(const struct rsb_mtx_t *mtxAp);

#if 0
rsb_err_t rsb_spmm_inner(const struct rsb_mtx_t * mtxAp, const void * mrhs, void *mout, rsb_int_t bstride, rsb_int_t cstride, rsb_int_t nrhs, rsb_trans_t transA);
#endif
rsb_long_t rsb__terminal_recursive_matrix_count_with_storage_and_flags(const struct rsb_mtx_t *mtxAp, rsb_fmt_t matrix_storage, rsb_flags_t flags);
rsb_long_t rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(const struct rsb_mtx_t *mtxAp, rsb_fmt_t matrix_storage, rsb_flags_t flags);
rsb_err_t rsb__do_compute_terminal_nnz_min_max_avg_count(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * minnz, rsb_nnz_idx_t * maxnz, rsb_nnz_idx_t * avgnz);
rsb_err_t rsb__do_compute_terminal_nnz_min_max_count(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * minnz, rsb_nnz_idx_t * maxnz);
rsb_char_t rsb__do_get_symmetry_char(const struct rsb_mtx_t *mtxAp);
rsb_flags_t rsb__do_flip_uplo_flags(rsb_flags_t flags);
rsb_flags_t rsb__do_detect_and_add_triangular_flags(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_flags_t flags);
rsb_trans_t rsb__do_transposition_from_char(rsb_char_t tc);
rsb_err_t rsb__do_load_matrix_file_as_matrix_market(struct rsb_mtx_t ** mtxApp, const rsb_char_t * filename, rsb_flags_t flags, rsb_type_t typecode);
rsb_err_t rsb__get_row_dense(const struct rsb_mtx_t * mtxAp, void* row, rsb_coo_idx_t i );
rsb_err_t rsb__do_cleanup_nnz(void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t *onnzp, rsb_type_t typecode, rsb_flags_t flags);
#if 0
rsb_err_t rsb_spmv_aa  (const struct rsb_mtx_t * mtxAp, const void * x, void * y, rsb_trans_t transA);
rsb_err_t rsb_spmv_sa  (const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, rsb_trans_t transA);
rsb_err_t rsb_spmv_unua  (const struct rsb_mtx_t * mtxAp, const void * x, void * y, rsb_trans_t transA);
rsb_err_t rsb_spmv_az  (const struct rsb_mtx_t * mtxAp, const void * x, void * y, rsb_trans_t transA);
rsb_err_t rsb_spmv_uxux  (const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA);
#endif
/*rsb_err_t rsb_cssm(struct rsb_mtx_t * mtxAp, void * x, const void * y, void * alpha, void * beta, rsb_trans_t transT);*/
#if 0
rsb_err_t rsb_spmm_az (const struct rsb_mtx_t * mtxAp, const void * mrhs, void *mout, rsb_int_t bstride, rsb_int_t cstride, rsb_int_t nrhs, rsb_trans_t transA);

rsb_err_t rsb_spmm_sxsx(const struct rsb_mtx_t * mtxAp, const void * b, void * c, rsb_nnz_idx_t ldb, rsb_nnz_idx_t ldc, rsb_coo_idx_t nrhs, rsb_trans_t transA, const void * alphap, const void * betap, rsb_flags_t order);
#endif
/* rsb_err_t rsb_spmv_sxsx(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_coo_idx_t incx, rsb_coo_idx_t incy); */

/*rsb_err_t rsb__get_row_dense(const struct rsb_mtx_t * mtxAp, void* row, rsb_coo_idx_t i );*/
/*rsb_err_t rsb_get_rows_dense(const struct rsb_mtx_t * mtxAp, void* row, rsb_coo_idx_t fr, rsb_coo_idx_t lr, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t *rnz, rsb_flags_t flags );*/
rsb_err_t rsb__do_set_elements(struct rsb_mtx_t * mtxAp, const void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags);
struct rsb_mtx_t * rsb__load_matrix_file_as_binary(const rsb_char_t * filename, rsb_err_t *errvalp);
rsb_err_t rsb_spmv_uaua(const struct rsb_mtx_t * mtxAp, const void * rhs, void * out, rsb_trans_t transA);
rsb_err_t rsb_spmv_uauz(const struct rsb_mtx_t * mtxAp, const void * rhs, void * out, rsb_trans_t transA);
rsb_err_t rsb__do_spsm(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * betap, const void * b, rsb_nnz_idx_t ldb, void * c, rsb_nnz_idx_t ldc);
rsb_err_t rsb__print_matrix_unsorted_coo(const struct rsb_mtx_t *mtxAp);
rsb_err_t rsb__util_sort_row_major_buffered(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k,  rsb_type_t typecode , rsb_flags_t flags, void * WA, size_t wb );
rsb_err_t rsb__set_ldX_for_spmm(rsb_trans_t transA, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, rsb_nnz_idx_t * ldBp, rsb_nnz_idx_t * ldCp);
rsb_err_t rsb__do_spmm(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * b, rsb_nnz_idx_t ldb, const void * betap, void * c, rsb_nnz_idx_t ldc, enum rsb_op_flags_t op_flags);
rsb_err_t rsb__do_spmm_general(const struct rsb_mtx_t * mtxAp, const void * b, void * c, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA, enum rsb_op_flags_t op_flags, rsb_flags_t order,const rsb_int_t nrhs, const size_t outnri, const size_t rhsnri);
rsb_err_t rsb__do_transpose(struct rsb_mtx_t ** mtxApp, rsb_bool_t want_conj);
rsb_err_t rsb__do_get_elements(const struct rsb_mtx_t * mtxAp, void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags);
rsb_err_t rsb__stropts_set(const rsb_char_t *opn, const rsb_char_t *arg);/* FIXME: in stropts.c */
rsb_err_t rsb__do_set_initopt_as_string(const rsb_char_t *opn, const rsb_char_t *arg);
rsb_err_t rsb__do_lib_get_info_str(int what, rsb_char_t* sbuf, size_t buflen);
struct rsb_mtx_t * rsb__do_mtx_alloc_from_coo_inplace(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb__do_mtx_alloc_from_coo_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb__do_mtx_alloc_from_csr_const(const void *VA, const rsb_coo_idx_t * RP, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb__do_mtx_alloc_from_csc_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * CP, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp);
int rsb__util_atoi_km10(const rsb_char_t *nptr);
int rsb__util_atoi_km2(const rsb_char_t *nptr);
void rsb__cat_compver(rsb_char_t * buf);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_spata(const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_char_t *rsb__util_strcat(rsb_char_t *dest, const rsb_char_t *src);

#define RSB__FLAGS_PRINTF_ARGS(FLAGS)	/* print matrix flags in human readable form  */ \
	"( 0x%x = { rec:%d coo:%d css:%d hw:%d ic:%d fi:%d symflags:" RSB_PRINTF_FLAGS_FMT_STR " } )\n", \
		(int)(FLAGS), \
		(int)(RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_QUAD_PARTITIONING)), \
		(int)(RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_WANT_COO_STORAGE)), \
		(int)(RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_WANT_BCSS_STORAGE)), \
		(int)(RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_USE_HALFWORD_INDICES)), \
		(int)(RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS)), \
		(int)(RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_FORTRAN_INDICES_INTERFACE)), \
		RSB_PRINTF_FLAGS_ARGS((FLAGS))  
void rsb__debug_print_flags(rsb_flags_t flags);
rsb_err_t rsb__print_matrix_stats(const struct rsb_mtx_t * mtxAp);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_INTERNALS_H_INCLUDED */
/* @endcond */
