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
/**
 * @file
 * @author Michele Martone
 * @brief  Sparse BLAS interface internals
 * */
#ifndef LIBSPBLAS_HANDLE_H_INCLUDED
#define LIBSPBLAS_HANDLE_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb.h"
#include "blas_sparse/blas_enum.h"
#include "rsb_libspblas.h"
#include "rsb_internals.h"

#define RSB_BLAS_FIRST_HANDLE 1024
#define RSB_BLAS_LAST_HANDLE (RSB_MAX_VALUE_FOR_TYPE(blas_sparse_matrix)-1)
#define RSB_BLAS_MATRICES_MAX (RSB_BLAS_LAST_HANDLE-RSB_BLAS_FIRST_HANDLE+1)
#define RSB_BLAS_HANDLE_INVALID (RSB_BLAS_LAST_HANDLE+1)
#define RSB_BLAS_INT_MAX RSB_MAX_VALUE_FOR_TYPE(rsb_blas_int_t)

#define RSB_BLAS_HANDLE_TO_RSB_ERROR(HANDLE) (HANDLE==(RSB_BLAS_HANDLE_INVALID)?RSB_ERR_GENERIC_ERROR:RSB_ERR_NO_ERROR)

#define RSB_ERROR_TO_BLAS_ERROR(E) (((E)==RSB_ERR_NO_ERROR)?(RSB_BLAS_NO_ERROR):(RSB_BLAS_ERROR))
#define RSB_BLAS_ERROR_TO_RSB_ERROR(E) (((E)==RSB_BLAS_NO_ERROR)?(RSB_ERR_NO_ERROR):(RSB_ERR_GENERIC_ERROR))
#define RSB_BLAS_WANT_EXPERIMENTAL_TUNING RSB_WANT_OMP_RECURSIVE_KERNELS	/* FIXME: this is experimental */

#define RSB_BLAS_UPLO_CHAR(UPLO) ((UPLO)==blas_lower_triangular?'L':'U')
#define RSB_BLAS_DIAG_CHAR(DIAG) ((DIAG)==blas_non_unit_diag?'E':'I')
#define RSB_BLAS_INVALID_VAL -1

#define RSB_WMSG_INIT	"Are you sure to have initialized the library? Seems not! Expect a crash..." /* TODO: move out of here */
#define RSB_CHECK_LIB_INIT_CHECK	if(!rsb__do_was_initialized()) RSB_WARN(RSB_WMSG_INIT "\n"); /* TODO: this shall become an error */
#define RSB_SPB_INTERFACE_PREAMBLE RSB_INTERFACE_PREAMBLE RSB_CHECK_LIB_INIT_CHECK
#define RSB_SPB_INTERFACE_RETURN(EXP) { int istat = EXP; RSB_INTERFACE_ENDCMD RSB_DO_ERR_MANIFEST_INTERFACE(RSB_BLAS_ERROR_TO_RSB_ERROR(istat)) return istat; }
#define RSB_SPB_INTERFACE_RETURN_HDL(EXP) { int handle = EXP; RSB_INTERFACE_ENDCMD RSB_DO_ERR_MANIFEST_INTERFACE(RSB_BLAS_HANDLE_TO_RSB_ERROR(handle)) return handle; }
#define RSB_SPB_INTERFACE_RETURN_VOID() { RSB_INTERFACE_ENDCMD return; }
#define RSB_SPB_INTERFACE_RETURN_EXP(EXP) { RSB_INTERFACE_ENDCMD return (EXP); }

/*typedef rsb_blas_sparse_matrix_handle_t blas_sparse_matrix;*/
typedef int rsb_blas_pname_t; /* here and in generated code */

/*!
 * \ingroup rsb_doc_sparse_blas
 * \brief An internal, helper structure.
 * \internal
 */
struct rsb_blas_sparse_matrix_t
{
	struct rsb_mtx_t * mtxAp;
	struct rsb_coo_mtx_t coomatrix;
	rsb_nnz_idx_t nnzin;
	blas_sparse_matrix handle;
	/* rsb_blas_int_t prop ;*/
	int k, l, off;
	rsb_blas_int_t *rbp,*cbp;
	enum blas_handle_type   type;
	enum blas_diag_type diag_type;
	enum blas_symmetry_type symmetry;
	enum blas_base_type base;
	enum blas_order_type order;
	enum blas_field_type fprecision;
	enum blas_field_type field;
	enum blas_sparsity_optimization_type sparsity_optimization_type;
	enum blas_rsb_ext_type dupstra;
	enum blas_rsb_ext_type fmt_hint;
#if RSB_BLAS_WANT_EXPERIMENTAL_TUNING 
	rsb_thread_t opt_mvn_hint, opt_mvt_hint;
#endif
/*  	we should also deal with :
blas_order_type 
blas_trans_type 
blas_uplo_type  
blas_diag_type 
blas_side_type 
blas_cmach_type 
blas_norm_type 
blas_sort_type 
blas_conj_type 
blas_jrot_type 
blas_prec_type 
blas_base_type 
blas_symmetry_type 
blas_field_type 
blas_size_type 
blas_handle_type
blas_sparsity_optimization_type 
 */
};



rsb_err_t rsb__BLAS_is_type_supported(rsb_char_t c);
struct rsb_mtx_t * rsb__BLAS_inner_matrix_retrieve(blas_sparse_matrix handle);
blas_sparse_matrix rsb__BLAS_Xuscr_begin( rsb_blas_int_t m, rsb_blas_int_t n, rsb_type_t typecode);
blas_sparse_matrix rsb__BLAS_Xuscr_block_begin( rsb_blas_int_t Mb, rsb_blas_int_t Nb, rsb_blas_int_t k, rsb_blas_int_t l, rsb_type_t typecode);
blas_sparse_matrix rsb__BLAS_Xuscr_variable_block_begin( rsb_blas_int_t Mb, rsb_blas_int_t Nb, const rsb_blas_int_t *k, const rsb_blas_int_t *l, rsb_type_t typecode);
rsb_blas_int_t rsb__BLAS_Xuscr_insert_entry( blas_sparse_matrix A, const void * valp, rsb_blas_int_t i, rsb_blas_int_t j );
rsb_blas_int_t rsb__BLAS_Xuscr_insert_entries( blas_sparse_matrix A, rsb_blas_int_t nz, const void * val, const rsb_blas_int_t *indx, const rsb_blas_int_t *jndx );
rsb_blas_int_t rsb__BLAS_Xuscr_insert_col( blas_sparse_matrix A, rsb_blas_int_t j, rsb_blas_int_t nz, const void * val, const rsb_blas_int_t *indx );
rsb_blas_int_t rsb__BLAS_Xuscr_insert_row( blas_sparse_matrix A, rsb_blas_int_t i, rsb_blas_int_t nz, const void * val, const rsb_blas_int_t *jndx );
rsb_blas_int_t rsb__BLAS_Xuscr_insert_clique( blas_sparse_matrix A, const rsb_blas_int_t k, const rsb_blas_int_t l, const void * val, const rsb_blas_int_t row_stride, const rsb_blas_int_t col_stride, const rsb_blas_int_t *indx, const rsb_blas_int_t *jndx );
rsb_blas_int_t rsb__BLAS_Xuscr_insert_block( blas_sparse_matrix A, const void * val, rsb_blas_int_t row_stride, rsb_blas_int_t col_stride, rsb_blas_int_t i, rsb_blas_int_t j);
rsb_blas_int_t rsb__BLAS_Xuscr_end( blas_sparse_matrix A);
rsb_blas_int_t rsb__BLAS_Xuscr_end_flagged( blas_sparse_matrix A, const rsb_flags_t*flagsp);

rsb_blas_int_t rsb__BLAS_usgp( blas_sparse_matrix A, rsb_blas_int_t pname );
rsb_blas_int_t rsb__BLAS_ussp( blas_sparse_matrix A, rsb_blas_int_t pname );
rsb_blas_int_t rsb__BLAS_Xusds( blas_sparse_matrix A );
rsb_trans_t rsb__blas_trans_to_rsb_trans(enum blas_trans_type trans);
rsb_trans_t rsb__do_psblas_trans_to_rsb_trans(const char trans);
rsb_order_t rsb__blas_order_to_rsb_order(enum blas_order_type order);
blas_sparse_matrix rsb__BLAS_new_matrix_begin(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnzest, rsb_type_t typecode, rsb_coo_idx_t br, rsb_coo_idx_t bc, const rsb_coo_idx_t*rbp, const rsb_coo_idx_t*cbp);
rsb_err_t rsb__BLAS_handles_free(void);

int rsb__BLAS_Xusrows_scale(blas_sparse_matrix A,const void * d,enum blas_trans_type trans);
int rsb__BLAS_Xusget_diag(blas_sparse_matrix A,void * d);
int rsb__BLAS_Xusget_rows_sparse(blas_sparse_matrix A,void *  VA, rsb_blas_int_t * IA, rsb_blas_int_t * JA, rsb_blas_int_t * nnz, rsb_blas_int_t fr, rsb_blas_int_t lr);
int rsb__BLAS_Xusget_matrix_nnz(blas_sparse_matrix A, rsb_blas_int_t * nnz);
int rsb__BLAS_Xusget_infinity_norm(blas_sparse_matrix A, void * in, enum blas_trans_type trans);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__BLAS_Xusget_rows_sums(blas_sparse_matrix A, void * rs, enum blas_trans_type trans);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
int rsb__BLAS_Xusset_elements(blas_sparse_matrix A,const rsb_blas_int_t * ia, const rsb_blas_int_t *ja, const void *  va, rsb_blas_int_t nnz);
int rsb__BLAS_Xusset_element(blas_sparse_matrix A,rsb_blas_int_t i, rsb_blas_int_t j, const void * v);
int rsb__BLAS_Xusget_element(blas_sparse_matrix A,rsb_blas_int_t i, rsb_blas_int_t j, void * v);
int rsb__BLAS_Xusmv(enum blas_trans_type transA, const void * alphap, blas_sparse_matrix A, const void * Xp, rsb_blas_int_t incX, const void * betap, void * Yp, rsb_blas_int_t incY);
int rsb__BLAS_Xusmm(enum blas_trans_type transA, const void * alphap, blas_sparse_matrix A, const void * b, rsb_blas_int_t ldb, const void * betap, void * c, rsb_blas_int_t ldc, rsb_blas_int_t nrhs, enum blas_order_type order);
int rsb__BLAS_Xussv(enum blas_trans_type transT, void * alpha, blas_sparse_matrix T, void * x, rsb_blas_int_t incx);
void blas_usgp_f_( blas_sparse_matrix*A, rsb_blas_int_t * pname, rsb_blas_int_t * istat );
void blas_ussp_f_( blas_sparse_matrix*A, rsb_blas_int_t * pname, rsb_blas_int_t * istat );
blas_sparse_matrix rsb__load_spblas_matrix_file_as_matrix_market(const rsb_char_t * filename, rsb_type_t typecode );
rsb_blas_int_t rsb__BLAS_Xusget_rows_nnz( blas_sparse_matrix A, rsb_blas_int_t fr, rsb_blas_int_t lr, rsb_blas_int_t * nnzp);
blas_sparse_matrix rsb__BLAS_handle_free(blas_sparse_matrix handle);

#define RSB_BLAS_STDOUT_MATRIX_SUMMARY(A)						\
	{										\
		struct rsb_mtx_t *RSB_DUMMY_ID=rsb__BLAS_inner_matrix_retrieve(A);	\
		if(RSB_DUMMY_ID)RSB_STDOUT_MATRIX_SUMMARY(RSB_DUMMY_ID );		\
		/*else            RSB_STDOUT_MATRIX_SUMMARY(RSB_DUMMY_MTX);*/		\
	}

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* LIBSPBLAS_HANDLE_H_INCLUDED */
/* @endcond */
