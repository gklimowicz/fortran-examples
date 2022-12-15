

/*!
        @file
        @author Michele Martone

	@brief  This file implements Sparse BLAS for librsb.
	Supported types  :(float,double,float complex,double complex) .
	Unsupported types:() .
	Level 1 ops      :(dot,axpy,ga,gz,sc) .
	Level 2 ops      :(mv,sv) .
	Level 3 ops      :(mm,sm) .
*/

#ifndef BLAS_ENUM_H
#define BLAS_ENUM_H

  /* Enumerated types */

/*! Used to specify a dense array's elements layout. */
enum blas_order_type {
            blas_rowmajor = 101, /*!< Row major. */
            blas_colmajor = 102  /*!< Column major. */ };

/*! Used to specify a transposition operator to a matrix operand. */
enum blas_trans_type {
            blas_no_trans   = 111, /*!< No transposition. */
            blas_trans      = 112, /*!< Transposition. */
            blas_conj_trans = 113  /*!< Transposition and conjugation. */ };

/*! Specifies (#BLAS_ussp) or inquiries (#BLAS_usgp) upper or lower triangularity of a matrix. */
enum blas_uplo_type  {
            blas_upper = 121, /*!< Upper triangular matrix. */
            blas_lower = 122  /*!< Lower triangular matrix. */ };

/*! Specifies (#BLAS_ussp) or inquiries (#BLAS_usgp) whether the diagonal of a matrix is (implicitly) unitary or not. */
enum blas_diag_type {
            blas_non_unit_diag = 131,  /*!< Unit diagional matrix. */
            blas_unit_diag     = 132   /*!< Non unit diagional matrix (the default). */ };

/*! Unused/Unsupported. */
enum blas_side_type {
            blas_left_side  = 141, /*!< Unsupported. */ 
            blas_right_side = 142  /*!< Unsupported. */  };

/*! Unused/Unsupported. */
enum blas_cmach_type {
            blas_base      = 151, /*!< Unsupported. */ 
            blas_t         = 152, /*!< Unsupported. */ 
            blas_rnd       = 153, /*!< Unsupported. */ 
            blas_ieee      = 154, /*!< Unsupported. */ 
            blas_emin      = 155, /*!< Unsupported. */ 
            blas_emax      = 156, /*!< Unsupported. */ 
            blas_eps       = 157, /*!< Unsupported. */ 
            blas_prec      = 158, /*!< Unsupported. */ 
            blas_underflow = 159, /*!< Unsupported. */ 
            blas_overflow  = 160, /*!< Unsupported. */ 
            blas_sfmin     = 161  /*!< Unsupported. */ };

/*! Unused/Unsupported. */
enum blas_norm_type {
            blas_one_norm       = 171, /*!< Unsupported. */ 
            blas_real_one_norm  = 172, /*!< Unsupported. */ 
            blas_two_norm       = 173, /*!< Unsupported. */ 
            blas_frobenius_norm = 174, /*!< Unsupported. */ 
            blas_inf_norm       = 175, /*!< Unsupported. */ 
            blas_real_inf_norm  = 176, /*!< Unsupported. */ 
            blas_max_norm       = 177, /*!< Unsupported. */ 
            blas_real_max_norm  = 178  /*!< Unsupported. */ };

/*! Unused/Unsupported. */
enum blas_sort_type {
            blas_increasing_order = 181,  /*!< Unsupported. */ 
            blas_decreasing_order = 182   /*!< Unsupported. */  };

/*! Unused/Unsupported. */
enum blas_conj_type {
            blas_conj    = 191, /*!< Unsupported. */
            blas_no_conj = 192  /*!< Unsupported. */ };

/*! Unused/Unsupported. */
enum blas_jrot_type {
            blas_jrot_inner  = 201, /*!< Unsupported. */
            blas_jrot_outer  = 202, /*!< Unsupported. */
            blas_jrot_sorted = 203  /*!< Unsupported. */ };

/*! Unused/Unsupported. */
enum blas_prec_type {
            blas_prec_single     = 211, /*!< Unsupported. */
            blas_prec_double     = 212, /*!< Unsupported. */
            blas_prec_indigenous = 213, /*!< Unsupported. */
            blas_prec_extra      = 214  /*!< Unsupported. */ };

/*! Index base (valid at matrix build/modify time). */
enum blas_base_type {
            blas_zero_base = 221, /*!< Zero based indices (default when matrix created using the C interface). */
            blas_one_base  = 222  /*!< Zero based indices (default when matrix created using the Fortran interface). */ };

/*! Symmetry properties. If not specified otherwise, valid for the both of #BLAS_ussp and #BLAS_usgp.
 */
enum blas_symmetry_type {
            blas_general          = 231, /*!< General unsymmetric matrix (default). For #BLAS_usgp only. */
            blas_symmetric        = 232, /*!< Symmetric matrix (either #blas_lower_symmetric or #blas_upper_symmetric). For #BLAS_usgp only. */
            blas_hermitian        = 233, /*!< Hermitian matrix (either #blas_lower_hermitian or #blas_upper_hermitian). For #BLAS_usgp only. */
            blas_triangular       = 234, /*!< Triangular matrix (either #blas_lower_triangular or #blas_upper_triangular). For #BLAS_usgp only. */
            blas_lower_triangular = 235, /*!< Lower triangular matrix. */
            blas_upper_triangular = 236, /*!< Upper triangular matrix. */
            blas_lower_symmetric  = 237, /*!< Lower symmetric matrix. */
            blas_upper_symmetric  = 238, /*!< Upper symmetric matrix. */
            blas_lower_hermitian  = 239, /*!< Lower hermitian matrix. */
            blas_upper_hermitian  = 240  /*!< Upper hermitian matrix. */ };

/*! Numerical field type; can be used with #BLAS_usgp to inquiry about a matrix numerical type (1 will be returned in case of success, 0 in case of failure). */
enum blas_field_type {
            blas_complex          = 241, /*!< Will succeed if matrix is of 'C' or 'Z' type. */
            blas_real             = 242, /*!< Will succeed if matrix is of 'S' or 'D' type. */
            blas_double_precision = 243, /*!< Will succeed if matrix is of 'D' or 'Z' type. */
            blas_single_precision = 244  /*!< Will succeed if matrix is of 'S' or 'C' type. */ };

/*! Quantities that can be obtained via #BLAS_usgp. */
enum blas_size_type {
            blas_num_rows      = 251, /*!< Get the matrix rows count. */
            blas_num_cols      = 252, /*!< Get the matrix columns count. */
            blas_num_nonzeros  = 253  /*!< Get the matrix nonzeros count. */ };

/*! The following are not fully implemented. Usable with #BLAS_usgp. */
enum blas_handle_type{
            blas_invalid_handle = 261, /*!< Used to check whether the handle is invalid. */
			blas_new_handle     = 262, /*!< Will give 1 if the handle is new. */
			blas_open_handle    = 263, /*!< will give 1 if the handle is open. */
			blas_valid_handle   = 264  /*!< Will give 1 if the handle is valid (that is, after #BLAS_duscr_end/#BLAS_zuscr_end/#BLAS_cuscr_end/#BLAS_zuscr_end). */ };

/*! The following are usable with #BLAS_usgp only. */
enum blas_sparsity_optimization_type {
            blas_regular       = 271, /*!< Will give 0. */
            blas_irregular     = 272, /*!< Will give 1. */
            blas_block         = 273, /*!< Will give 0. */
            blas_unassembled   = 274  /*!< Complementary to #blas_valid_handle. */ };

/*! Properties suitable to be used with #BLAS_ussp/#BLAS_usgp. All of these are not in the Sparse BLAS standard. */
enum blas_rsb_ext_type {
            blas_rsb_spmv_autotuning_on   = 6660,	/*!< Turn on executing threads autotuning for #BLAS_dusmv, #BLAS_zusmv, #BLAS_susmv, #BLAS_cusmv. As an extension to the standard, the autotuning properties can be turned on/off at any time; if the autotuning feature has not been enabled at build time, using these properties will make the call fail. For more information, see #rsb_tune_spmm. (EXPERIMENTAL) */
            blas_rsb_spmv_autotuning_off  = 6661,	/*!< Turn off executing threads autotuning for #BLAS_dusmv, #BLAS_zusmv, #BLAS_susmv, #BLAS_cusmv. See #blas_rsb_spmv_autotuning_on. (EXPERIMENTAL) */
            blas_rsb_spmv_n_autotuning_on   = 6662,	/*!< Turn on executing threads autotuning for untransposed #BLAS_dusmv, #BLAS_zusmv, #BLAS_susmv, #BLAS_cusmv. See #blas_rsb_spmv_autotuning_on. (EXPERIMENTAL) */
            blas_rsb_spmv_n_autotuning_off  = 6663,	/*!< Turn on executing threads autotuning for untransposed #BLAS_dusmv, #BLAS_zusmv, #BLAS_susmv, #BLAS_cusmv. See #blas_rsb_spmv_autotuning_on. (EXPERIMENTAL) */
            blas_rsb_spmv_t_autotuning_on   = 6664,	/*!< Turn on executing threads autotuning for transposed #BLAS_dusmv, #BLAS_zusmv, #BLAS_susmv, #BLAS_cusmv. See #blas_rsb_spmv_autotuning_on. (EXPERIMENTAL) */
            blas_rsb_spmv_t_autotuning_off  = 6665,	/*!< Turn on executing threads autotuning for transposed #BLAS_dusmv, #BLAS_zusmv, #BLAS_susmv, #BLAS_cusmv. See #blas_rsb_spmv_autotuning_on. (EXPERIMENTAL) */
            blas_rsb_autotune_next_operation= 6666,	/*!< Turn on executing threads autotuning for the next operation among #BLAS_dusmv, #BLAS_zusmv, #BLAS_susmv, #BLAS_cusmv). See #blas_rsb_spmv_autotuning_on. (EXPERIMENTAL) */
            blas_rsb_rep_rec         = 9993,	/*!< Request/check for recursive representation. */
            blas_rsb_rep_hwi         = 9994,	/*!< Request/check for half-word indices. */
            blas_rsb_rep_rsb         = 9995,	/*!< Request/check for RSB representation. */
            blas_rsb_rep_csr         = 9996,	/*!< Request/check for CSR representation. */
            blas_rsb_rep_coo         = 9997,	/*!< Request/check for COO representation. */
            blas_rsb_duplicates_ovw   = 9998,	/*!< Request/check for duplicate nonzeroes overwriting policy. */
            blas_rsb_duplicates_sum   = 9999 	/*!< Request/check for duplicate nonzeroes summation policy. */
};

#endif
   /* BLAS_ENUM_H */
#include "rsb_libspblas_handle.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* @cond INNERDOC  */
RSB_INTERNALS_COMMON_HEAD_DECLS
/* @endcond */
               /* Level 1 Computational Routines */
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susdot rsb__wp_BLAS_susdot
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susdot(const enum blas_conj_type conj, const int nnz, const float * x,
		const int *indx, const float * y, const int incy, float * r,
		const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusdot(RSB_NUMERICAL_TYPE_FLOAT ,conj,nnz,x,indx,y,incy,r,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susdot_ rsb__wp_blas_susdot_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susdot_(const enum blas_conj_type*conj,const int*nnz,const float *x,const int *indx,const float *y,const int*incy,float *r,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_susdot(*conj,*nnz,x,indx,y,*incy,r,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusdot rsb__wp_BLAS_dusdot
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusdot(const enum blas_conj_type conj, const int nnz, const double * x,
		const int *indx, const double * y, const int incy, double * r,
		const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusdot(RSB_NUMERICAL_TYPE_DOUBLE ,conj,nnz,x,indx,y,incy,r,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusdot_ rsb__wp_blas_dusdot_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusdot_(const enum blas_conj_type*conj,const int*nnz,const double *x,const int *indx,const double *y,const int*incy,double *r,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_dusdot(*conj,*nnz,x,indx,y,*incy,r,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusdot rsb__wp_BLAS_cusdot
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusdot(const enum blas_conj_type conj, const int nnz, const void *x,
		const int *indx, const void *y, const int incy, void *r,
		const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusdot(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ,conj,nnz,x,indx,y,incy,r,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusdot_ rsb__wp_blas_cusdot_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusdot_(const enum blas_conj_type*conj,const int*nnz,const void *x,const int *indx,const void *y,const int*incy,void *r,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_cusdot(*conj,*nnz,x,indx,y,*incy,r,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusdot rsb__wp_BLAS_zusdot
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusdot(const enum blas_conj_type conj, const int nnz, const void *x,
		const int *indx, const void *y, const int incy, void *r,
		const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusdot(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ,conj,nnz,x,indx,y,incy,r,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusdot_ rsb__wp_blas_zusdot_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusdot_(const enum blas_conj_type*conj,const int*nnz,const void *x,const int *indx,const void *y,const int*incy,void *r,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_dot_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_zusdot(*conj,*nnz,x,indx,y,*incy,r,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susaxpy rsb__wp_BLAS_susaxpy
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susaxpy(const int nnz, float  alpha, const float * x, const int *indx,
                 float * y, const int incy, const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusaxpy(RSB_NUMERICAL_TYPE_FLOAT ,nnz,&alpha,x,indx,y,incy,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susaxpy_ rsb__wp_blas_susaxpy_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susaxpy_(const int*nnz,float*alpha,const float *x,const int *indx,float *y,const int*incy,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_susaxpy(*nnz,*alpha,x,indx,y,*incy,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusaxpy rsb__wp_BLAS_dusaxpy
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusaxpy(const int nnz, double  alpha, const double * x, const int *indx,
                 double * y, const int incy, const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusaxpy(RSB_NUMERICAL_TYPE_DOUBLE ,nnz,&alpha,x,indx,y,incy,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusaxpy_ rsb__wp_blas_dusaxpy_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusaxpy_(const int*nnz,double*alpha,const double *x,const int *indx,double *y,const int*incy,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_dusaxpy(*nnz,*alpha,x,indx,y,*incy,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusaxpy rsb__wp_BLAS_cusaxpy
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusaxpy(const int nnz, const void * alpha, const void *x, const int *indx,
                 void *y, const int incy, const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusaxpy(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ,/* FIXME: this is an exception; shall use a formal substitution technique, rather */nnz,alpha,x,indx,y,incy,index_base ))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusaxpy_ rsb__wp_blas_cusaxpy_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusaxpy_(const int*nnz,const void *alpha,const void *x,const int *indx,void *y,const int*incy,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_cusaxpy(*nnz,alpha,x,indx,y,*incy,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusaxpy rsb__wp_BLAS_zusaxpy
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusaxpy(const int nnz, const void * alpha, const void *x, const int *indx,
                 void *y, const int incy, const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusaxpy(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ,/* FIXME: this is an exception; shall use a formal substitution technique, rather */nnz,alpha,x,indx,y,incy,index_base ))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusaxpy_ rsb__wp_blas_zusaxpy_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusaxpy_(const int*nnz,const void *alpha,const void *x,const int *indx,void *y,const int*incy,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_axpy_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_zusaxpy(*nnz,alpha,x,indx,y,*incy,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susga rsb__wp_BLAS_susga
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susga(const int nnz, const float * y, const int incy, float * x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusga(RSB_NUMERICAL_TYPE_FLOAT ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susga_ rsb__wp_blas_susga_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susga_(const int*nnz,const float *y,const int*incy,float *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_susga(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusga rsb__wp_BLAS_dusga
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusga(const int nnz, const double * y, const int incy, double * x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusga(RSB_NUMERICAL_TYPE_DOUBLE ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusga_ rsb__wp_blas_dusga_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusga_(const int*nnz,const double *y,const int*incy,double *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_dusga(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusga rsb__wp_BLAS_cusga
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusga(const int nnz, const void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusga(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusga_ rsb__wp_blas_cusga_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusga_(const int*nnz,const void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_cusga(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusga rsb__wp_BLAS_zusga
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusga(const int nnz, const void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusga(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusga_ rsb__wp_blas_zusga_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusga_(const int*nnz,const void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_ga_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_zusga(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susgz rsb__wp_BLAS_susgz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susgz(const int nnz, float * y, const int incy, float * x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusgz(RSB_NUMERICAL_TYPE_FLOAT ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susgz_ rsb__wp_blas_susgz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susgz_(const int*nnz,float *y,const int*incy,float *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_susgz(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusgz rsb__wp_BLAS_dusgz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusgz(const int nnz, double * y, const int incy, double * x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusgz(RSB_NUMERICAL_TYPE_DOUBLE ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusgz_ rsb__wp_blas_dusgz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusgz_(const int*nnz,double *y,const int*incy,double *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_dusgz(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusgz rsb__wp_BLAS_cusgz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusgz(const int nnz, void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusgz(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusgz_ rsb__wp_blas_cusgz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusgz_(const int*nnz,void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_cusgz(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusgz rsb__wp_BLAS_zusgz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusgz(const int nnz, void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusgz(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ,nnz,y,incy,x,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusgz_ rsb__wp_blas_zusgz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusgz_(const int*nnz,void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_gz_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_zusgz(*nnz,y,*incy,x,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_sussc rsb__wp_BLAS_sussc
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_sussc(const int nnz, const float * x, float * y, const int incy, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xussc(RSB_NUMERICAL_TYPE_FLOAT ,nnz,x,y,incy,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_sussc_ rsb__wp_blas_sussc_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_sussc_(const int*nnz,const float *x,float *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_sussc(*nnz,x,y,*incy,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dussc rsb__wp_BLAS_dussc
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dussc(const int nnz, const double * x, double * y, const int incy, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xussc(RSB_NUMERICAL_TYPE_DOUBLE ,nnz,x,y,incy,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dussc_ rsb__wp_blas_dussc_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dussc_(const int*nnz,const double *x,double *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_dussc(*nnz,x,y,*incy,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cussc rsb__wp_BLAS_cussc
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cussc(const int nnz, const void *x, void *y, const int incy, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xussc(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ,nnz,x,y,incy,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cussc_ rsb__wp_blas_cussc_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cussc_(const int*nnz,const void *x,void *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_cussc(*nnz,x,y,*incy,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zussc rsb__wp_BLAS_zussc
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zussc(const int nnz, const void *x, void *y, const int incy, const int *indx,
              const enum blas_base_type index_base)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_return_msg
	  \warning \rsb_spblasl1_msg
	*/
	RSB_SPB_INTERFACE_PREAMBLE

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xussc(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ,nnz,x,y,incy,indx,index_base))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zussc_ rsb__wp_blas_zussc_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zussc_(const int*nnz,const void *x,void *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat)
{

	/*!
	  \ingroup rsb_doc_sparse_blas
	  \rsb_spblasl1_sc_msg\rsb_spblas_istat_msg
	  \warning \rsb_spblasl1_msg
	*/

	int istatv = BLAS_zussc(*nnz,x,y,*incy,indx,*index_base );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
}

               /* Level 2 Computational Routines */

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susmv rsb__wp_BLAS_susmv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susmv(const enum blas_trans_type transA, float alpha,
    const blas_sparse_matrix A, const float * x, const int incx, float * y, const int incy)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const float beta = ((float)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmv(transA,&alpha,A,x,incx,&beta,y,incy))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susmv_ rsb__wp_blas_susmv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susmv_(const enum blas_trans_type*transA,float*alpha,const blas_sparse_matrix*A,const float *x,const int*incx,float *y,const int*incy,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susmv(*transA,*alpha,*A,x,*incx,y,*incy);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusmv rsb__wp_BLAS_dusmv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusmv(const enum blas_trans_type transA, double alpha,
    const blas_sparse_matrix A, const double * x, const int incx, double * y, const int incy)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const double beta = ((double)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmv(transA,&alpha,A,x,incx,&beta,y,incy))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusmv_ rsb__wp_blas_dusmv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusmv_(const enum blas_trans_type*transA,double*alpha,const blas_sparse_matrix*A,const double *x,const int*incx,double *y,const int*incy,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusmv(*transA,*alpha,*A,x,*incx,y,*incy);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusmv rsb__wp_BLAS_cusmv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusmv(const enum blas_trans_type transA, const void *alpha,
    const blas_sparse_matrix A, const void *x, const int incx, void *y, const int incy)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const float complex beta = ((float complex)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmv(transA,alpha,A,x,incx,&beta,y,incy))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusmv_ rsb__wp_blas_cusmv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusmv_(const enum blas_trans_type*transA,const void *alpha,const blas_sparse_matrix*A,const void *x,const int*incx,void *y,const int*incy,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusmv(*transA,alpha,*A,x,*incx,y,*incy);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusmv rsb__wp_BLAS_zusmv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusmv(const enum blas_trans_type transA, const void *alpha,
    const blas_sparse_matrix A, const void *x, const int incx, void *y, const int incy)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const double complex beta = ((double complex)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmv(transA,alpha,A,x,incx,&beta,y,incy))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusmv_ rsb__wp_blas_zusmv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusmv_(const enum blas_trans_type*transA,const void *alpha,const blas_sparse_matrix*A,const void *x,const int*incx,void *y,const int*incy,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusmv(*transA,alpha,*A,x,*incx,y,*incy);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}


#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_sussv rsb__wp_BLAS_sussv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_sussv(enum blas_trans_type transT, float alpha,
    const blas_sparse_matrix T, float * x, const int incx)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(T);
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsv(rsb__blas_trans_to_rsb_trans(transT),&alpha,mtxAp,x,incx,x,incx)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_sussv_ rsb__wp_blas_sussv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_sussv_(enum blas_trans_type*transT,float*alpha,const blas_sparse_matrix*T,float *x,const int*incx,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_sussv(*transT,*alpha,*T,x,*incx);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dussv rsb__wp_BLAS_dussv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dussv(enum blas_trans_type transT, double alpha,
    const blas_sparse_matrix T, double * x, const int incx)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(T);
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsv(rsb__blas_trans_to_rsb_trans(transT),&alpha,mtxAp,x,incx,x,incx)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dussv_ rsb__wp_blas_dussv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dussv_(enum blas_trans_type*transT,double*alpha,const blas_sparse_matrix*T,double *x,const int*incx,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dussv(*transT,*alpha,*T,x,*incx);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cussv rsb__wp_BLAS_cussv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cussv(enum blas_trans_type transT, const void *alpha,
    const blas_sparse_matrix T, void *x, const int incx)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(T);
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsv(rsb__blas_trans_to_rsb_trans(transT),alpha,mtxAp,x,incx,x,incx)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cussv_ rsb__wp_blas_cussv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cussv_(enum blas_trans_type*transT,const void *alpha,const blas_sparse_matrix*T,void *x,const int*incx,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cussv(*transT,alpha,*T,x,*incx);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zussv rsb__wp_BLAS_zussv
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zussv(enum blas_trans_type transT, const void *alpha,
    const blas_sparse_matrix T, void *x, const int incx)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(T);
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsv(rsb__blas_trans_to_rsb_trans(transT),alpha,mtxAp,x,incx,x,incx)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zussv_ rsb__wp_blas_zussv_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zussv_(enum blas_trans_type*transT,const void *alpha,const blas_sparse_matrix*T,void *x,const int*incx,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sv_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zussv(*transT,alpha,*T,x,*incx);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}


               /* Level 3 Computational Routines */

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susmm rsb__wp_BLAS_susmm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, float alpha, const blas_sparse_matrix A, const float * b, const int ldb,
       float *  c, const int ldc)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const float beta = ((float)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmm(transA,&alpha,A,b,ldb,&beta,c,ldc,nrhs,order))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susmm_ rsb__wp_blas_susmm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,float*alpha,const blas_sparse_matrix*A,const float *b,const int*ldb,float *c,const int*ldc,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susmm(*order,*transA,*nrhs,*alpha,*A,b,*ldb,c,*ldc);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusmm rsb__wp_BLAS_dusmm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, double alpha, const blas_sparse_matrix A, const double * b, const int ldb,
       double *  c, const int ldc)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const double beta = ((double)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmm(transA,&alpha,A,b,ldb,&beta,c,ldc,nrhs,order))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusmm_ rsb__wp_blas_dusmm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,double*alpha,const blas_sparse_matrix*A,const double *b,const int*ldb,double *c,const int*ldc,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusmm(*order,*transA,*nrhs,*alpha,*A,b,*ldb,c,*ldc);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusmm rsb__wp_BLAS_cusmm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, const void *alpha, const blas_sparse_matrix A, const void *b, const int ldb,
       void * c, const int ldc)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const float complex beta = ((float complex)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmm(transA,alpha,A,b,ldb,&beta,c,ldc,nrhs,order))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusmm_ rsb__wp_blas_cusmm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,const void *alpha,const blas_sparse_matrix*A,const void *b,const int*ldb,void *c,const int*ldc,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusmm(*order,*transA,*nrhs,alpha,*A,b,*ldb,c,*ldc);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusmm rsb__wp_BLAS_zusmm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, const void *alpha, const blas_sparse_matrix A, const void *b, const int ldb,
       void * c, const int ldc)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const double complex beta = ((double complex)(1.0));
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmm(transA,alpha,A,b,ldb,&beta,c,ldc,nrhs,order))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusmm_ rsb__wp_blas_zusmm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,const void *alpha,const blas_sparse_matrix*A,const void *b,const int*ldb,void *c,const int*ldc,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_mm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusmm(*order,*transA,*nrhs,alpha,*A,b,*ldb,c,*ldc);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}


#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_sussm rsb__wp_BLAS_sussm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_sussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, float alpha, const blas_sparse_matrix T, float * b, const int ldb)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const float beta = ((float)(0));
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsm(rsb__blas_trans_to_rsb_trans(transT),&alpha,rsb__BLAS_inner_matrix_retrieve(T),nrhs,rsb__blas_order_to_rsb_order(order),&beta,b,ldb,b,ldb)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_sussm_ rsb__wp_blas_sussm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_sussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,float*alpha,const blas_sparse_matrix*T,float *b,const int*ldb,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_sussm(*order,*transT,*nrhs,*alpha,*T,b,*ldb);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dussm rsb__wp_BLAS_dussm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, double alpha, const blas_sparse_matrix T, double * b, const int ldb)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const double beta = ((double)(0));
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsm(rsb__blas_trans_to_rsb_trans(transT),&alpha,rsb__BLAS_inner_matrix_retrieve(T),nrhs,rsb__blas_order_to_rsb_order(order),&beta,b,ldb,b,ldb)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dussm_ rsb__wp_blas_dussm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,double*alpha,const blas_sparse_matrix*T,double *b,const int*ldb,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dussm(*order,*transT,*nrhs,*alpha,*T,b,*ldb);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cussm rsb__wp_BLAS_cussm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, const void *alpha, const blas_sparse_matrix T, void *b, const int ldb)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const float complex beta = ((float complex)(0));
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsm(rsb__blas_trans_to_rsb_trans(transT),alpha,rsb__BLAS_inner_matrix_retrieve(T),nrhs,rsb__blas_order_to_rsb_order(order),&beta,b,ldb,b,ldb)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cussm_ rsb__wp_blas_cussm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,const void *alpha,const blas_sparse_matrix*T,void *b,const int*ldb,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cussm(*order,*transT,*nrhs,alpha,*T,b,*ldb);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zussm rsb__wp_BLAS_zussm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, const void *alpha, const blas_sparse_matrix T, void *b, const int ldb)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

{
	const double complex beta = ((double complex)(0));
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsm(rsb__blas_trans_to_rsb_trans(transT),alpha,rsb__BLAS_inner_matrix_retrieve(T),nrhs,rsb__blas_order_to_rsb_order(order),&beta,b,ldb,b,ldb)))
}
	}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zussm_ rsb__wp_blas_zussm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,const void *alpha,const blas_sparse_matrix*T,void *b,const int*ldb,int*istat)
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_sm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zussm(*order,*transT,*nrhs,alpha,*T,b,*ldb);
	RSB_SET_IF_NOT_NULL(istat,istatv);
}


               /* Handle Management Routines */
               /*             +              */
               /* Creation Routines */
               /*             +              */
               /* Insertion Routines */
               /*             +              */
               /* Completion of Construction Routines */
               /*             +              */
               /* Matrix Property Routines */
               /*             +              */
               /* Destruction Routine */

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_begin rsb__wp_BLAS_suscr_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_suscr_begin( int m, int n )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_begin(m,n,RSB_NUMERICAL_TYPE_FLOAT ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_begin_ rsb__wp_blas_suscr_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_begin_( int*m,int*n,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_suscr_begin(*m,*n );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_begin rsb__wp_BLAS_duscr_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_duscr_begin( int m, int n )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_begin(m,n,RSB_NUMERICAL_TYPE_DOUBLE ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_begin_ rsb__wp_blas_duscr_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_begin_( int*m,int*n,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_duscr_begin(*m,*n );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_begin rsb__wp_BLAS_cuscr_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_cuscr_begin( int m, int n )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_begin(m,n,RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_begin_ rsb__wp_blas_cuscr_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_begin_( int*m,int*n,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_cuscr_begin(*m,*n );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_begin rsb__wp_BLAS_zuscr_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_zuscr_begin( int m, int n )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_begin(m,n,RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_begin_ rsb__wp_blas_zuscr_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_begin_( int*m,int*n,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_begin_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_zuscr_begin(*m,*n );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_block_begin rsb__wp_BLAS_suscr_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_suscr_block_begin( int Mb, int Nb, int k, int l )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_block_begin(Mb,Nb,k,l,RSB_NUMERICAL_TYPE_FLOAT ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_block_begin_ rsb__wp_blas_suscr_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_block_begin_( int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_suscr_block_begin(*Mb,*Nb,*k,*l );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_block_begin rsb__wp_BLAS_duscr_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_duscr_block_begin( int Mb, int Nb, int k, int l )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_block_begin(Mb,Nb,k,l,RSB_NUMERICAL_TYPE_DOUBLE ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_block_begin_ rsb__wp_blas_duscr_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_block_begin_( int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_duscr_block_begin(*Mb,*Nb,*k,*l );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_block_begin rsb__wp_BLAS_cuscr_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_cuscr_block_begin( int Mb, int Nb, int k, int l )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_block_begin(Mb,Nb,k,l,RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_block_begin_ rsb__wp_blas_cuscr_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_block_begin_( int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_cuscr_block_begin(*Mb,*Nb,*k,*l );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_block_begin rsb__wp_BLAS_zuscr_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_zuscr_block_begin( int Mb, int Nb, int k, int l )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_block_begin(Mb,Nb,k,l,RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_block_begin_ rsb__wp_blas_zuscr_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_block_begin_( int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_block_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_zuscr_block_begin(*Mb,*Nb,*k,*l );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_variable_block_begin rsb__wp_BLAS_suscr_variable_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_suscr_variable_block_begin( int Mb, int Nb,
		const int *K, const int *L )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_variable_block_begin(Mb,Nb,K,L,RSB_NUMERICAL_TYPE_FLOAT ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_variable_block_begin_ rsb__wp_blas_suscr_variable_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_variable_block_begin_( int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_suscr_variable_block_begin(*Mb,*Nb,K,L );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_variable_block_begin rsb__wp_BLAS_duscr_variable_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_duscr_variable_block_begin( int Mb, int Nb,
		const int *K, const int *L )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_variable_block_begin(Mb,Nb,K,L,RSB_NUMERICAL_TYPE_DOUBLE ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_variable_block_begin_ rsb__wp_blas_duscr_variable_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_variable_block_begin_( int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_duscr_variable_block_begin(*Mb,*Nb,K,L );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_variable_block_begin rsb__wp_BLAS_cuscr_variable_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_cuscr_variable_block_begin( int Mb, int Nb,
		const int *K, const int *L )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_variable_block_begin(Mb,Nb,K,L,RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_variable_block_begin_ rsb__wp_blas_cuscr_variable_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_variable_block_begin_( int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_cuscr_variable_block_begin(*Mb,*Nb,K,L );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_variable_block_begin rsb__wp_BLAS_zuscr_variable_block_begin
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
blas_sparse_matrix BLAS_zuscr_variable_block_begin( int Mb, int Nb,
		const int *K, const int *L )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN_HDL(rsb__BLAS_Xuscr_variable_block_begin(Mb,Nb,K,L,RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_variable_block_begin_ rsb__wp_blas_zuscr_variable_block_begin_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_variable_block_begin_( int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_vbr_msg\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg
         */

	int istatv = BLAS_zuscr_variable_block_begin(*Mb,*Nb,K,L );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_end rsb__wp_BLAS_suscr_end
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_suscr_end( blas_sparse_matrix A )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_end(A))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_end_ rsb__wp_blas_suscr_end_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_end_( blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_suscr_end(*A );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_end rsb__wp_BLAS_duscr_end
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_duscr_end( blas_sparse_matrix A )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_end(A))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_end_ rsb__wp_blas_duscr_end_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_end_( blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_duscr_end(*A );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_end rsb__wp_BLAS_cuscr_end
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cuscr_end( blas_sparse_matrix A )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_end(A))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_end_ rsb__wp_blas_cuscr_end_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_end_( blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cuscr_end(*A );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_end rsb__wp_BLAS_zuscr_end
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zuscr_end( blas_sparse_matrix A )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_end(A))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_end_ rsb__wp_blas_zuscr_end_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_end_( blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zuscr_end(*A );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_insert_entry rsb__wp_BLAS_suscr_insert_entry
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_suscr_insert_entry( blas_sparse_matrix A, float  val, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entry(A,&val,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_insert_entry_ rsb__wp_blas_suscr_insert_entry_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_insert_entry_( blas_sparse_matrix*A,float*val,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_suscr_insert_entry(*A,*val,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_insert_entry rsb__wp_BLAS_duscr_insert_entry
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_duscr_insert_entry( blas_sparse_matrix A, double  val, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entry(A,&val,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_insert_entry_ rsb__wp_blas_duscr_insert_entry_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_insert_entry_( blas_sparse_matrix*A,double*val,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_duscr_insert_entry(*A,*val,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_insert_entry rsb__wp_BLAS_cuscr_insert_entry
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cuscr_insert_entry( blas_sparse_matrix A, const void * val, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entry(A,val,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_insert_entry_ rsb__wp_blas_cuscr_insert_entry_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_insert_entry_( blas_sparse_matrix*A,const void *val,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cuscr_insert_entry(*A,val,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_insert_entry rsb__wp_BLAS_zuscr_insert_entry
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zuscr_insert_entry( blas_sparse_matrix A, const void * val, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entry(A,val,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_insert_entry_ rsb__wp_blas_zuscr_insert_entry_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_insert_entry_( blas_sparse_matrix*A,const void *val,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zuscr_insert_entry(*A,val,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_insert_entries rsb__wp_BLAS_suscr_insert_entries
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_suscr_insert_entries( blas_sparse_matrix A, int nnz, const float * val,
                            const int *indx, const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entries(A,nnz,val,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_insert_entries_ rsb__wp_blas_suscr_insert_entries_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_insert_entries_( blas_sparse_matrix*A,int*nnz,const float *val,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_suscr_insert_entries(*A,*nnz,val,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_insert_entries rsb__wp_BLAS_duscr_insert_entries
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_duscr_insert_entries( blas_sparse_matrix A, int nnz, const double * val,
                            const int *indx, const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entries(A,nnz,val,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_insert_entries_ rsb__wp_blas_duscr_insert_entries_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_insert_entries_( blas_sparse_matrix*A,int*nnz,const double *val,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_duscr_insert_entries(*A,*nnz,val,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_insert_entries rsb__wp_BLAS_cuscr_insert_entries
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cuscr_insert_entries( blas_sparse_matrix A, int nnz, const void *val,
                            const int *indx, const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entries(A,nnz,val,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_insert_entries_ rsb__wp_blas_cuscr_insert_entries_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_insert_entries_( blas_sparse_matrix*A,int*nnz,const void *val,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cuscr_insert_entries(*A,*nnz,val,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_insert_entries rsb__wp_BLAS_zuscr_insert_entries
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zuscr_insert_entries( blas_sparse_matrix A, int nnz, const void *val,
                            const int *indx, const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_entries(A,nnz,val,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_insert_entries_ rsb__wp_blas_zuscr_insert_entries_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_insert_entries_( blas_sparse_matrix*A,int*nnz,const void *val,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zuscr_insert_entries(*A,*nnz,val,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_insert_col rsb__wp_BLAS_suscr_insert_col
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_suscr_insert_col( blas_sparse_matrix A, int j, int nnz,
                           const float * val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_col(A,j,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_insert_col_ rsb__wp_blas_suscr_insert_col_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_insert_col_( blas_sparse_matrix*A,int*j,int*nnz,const float *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_suscr_insert_col(*A,*j,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_insert_col rsb__wp_BLAS_duscr_insert_col
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_duscr_insert_col( blas_sparse_matrix A, int j, int nnz,
                           const double * val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_col(A,j,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_insert_col_ rsb__wp_blas_duscr_insert_col_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_insert_col_( blas_sparse_matrix*A,int*j,int*nnz,const double *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_duscr_insert_col(*A,*j,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_insert_col rsb__wp_BLAS_cuscr_insert_col
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cuscr_insert_col( blas_sparse_matrix A, int j, int nnz,
                           const void *val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_col(A,j,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_insert_col_ rsb__wp_blas_cuscr_insert_col_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_insert_col_( blas_sparse_matrix*A,int*j,int*nnz,const void *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cuscr_insert_col(*A,*j,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_insert_col rsb__wp_BLAS_zuscr_insert_col
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zuscr_insert_col( blas_sparse_matrix A, int j, int nnz,
                           const void *val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_col(A,j,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_insert_col_ rsb__wp_blas_zuscr_insert_col_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_insert_col_( blas_sparse_matrix*A,int*j,int*nnz,const void *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zuscr_insert_col(*A,*j,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_insert_row rsb__wp_BLAS_suscr_insert_row
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_suscr_insert_row( blas_sparse_matrix A, int i, int nnz,
                           const float * val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_row(A,i,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_insert_row_ rsb__wp_blas_suscr_insert_row_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_insert_row_( blas_sparse_matrix*A,int*i,int*nnz,const float *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_suscr_insert_row(*A,*i,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_insert_row rsb__wp_BLAS_duscr_insert_row
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_duscr_insert_row( blas_sparse_matrix A, int i, int nnz,
                           const double * val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_row(A,i,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_insert_row_ rsb__wp_blas_duscr_insert_row_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_insert_row_( blas_sparse_matrix*A,int*i,int*nnz,const double *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_duscr_insert_row(*A,*i,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_insert_row rsb__wp_BLAS_cuscr_insert_row
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cuscr_insert_row( blas_sparse_matrix A, int i, int nnz,
                           const void *val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_row(A,i,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_insert_row_ rsb__wp_blas_cuscr_insert_row_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_insert_row_( blas_sparse_matrix*A,int*i,int*nnz,const void *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cuscr_insert_row(*A,*i,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_insert_row rsb__wp_BLAS_zuscr_insert_row
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zuscr_insert_row( blas_sparse_matrix A, int i, int nnz,
                           const void *val, const int *indx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_row(A,i,nnz,val,indx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_insert_row_ rsb__wp_blas_zuscr_insert_row_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_insert_row_( blas_sparse_matrix*A,int*i,int*nnz,const void *val,const int *indx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zuscr_insert_row(*A,*i,*nnz,val,indx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_insert_clique rsb__wp_BLAS_suscr_insert_clique
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_suscr_insert_clique( blas_sparse_matrix A, const int k, const int l,
                       const float * val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_clique(A,k,l,val,row_stride,col_stride,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_insert_clique_ rsb__wp_blas_suscr_insert_clique_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_insert_clique_( blas_sparse_matrix*A,const int*k,const int*l,const float *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_suscr_insert_clique(*A,*k,*l,val,*row_stride,*col_stride,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_insert_clique rsb__wp_BLAS_duscr_insert_clique
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_duscr_insert_clique( blas_sparse_matrix A, const int k, const int l,
                       const double * val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_clique(A,k,l,val,row_stride,col_stride,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_insert_clique_ rsb__wp_blas_duscr_insert_clique_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_insert_clique_( blas_sparse_matrix*A,const int*k,const int*l,const double *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_duscr_insert_clique(*A,*k,*l,val,*row_stride,*col_stride,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_insert_clique rsb__wp_BLAS_cuscr_insert_clique
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cuscr_insert_clique( blas_sparse_matrix A, const int k, const int l,
                       const void *val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_clique(A,k,l,val,row_stride,col_stride,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_insert_clique_ rsb__wp_blas_cuscr_insert_clique_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_insert_clique_( blas_sparse_matrix*A,const int*k,const int*l,const void *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cuscr_insert_clique(*A,*k,*l,val,*row_stride,*col_stride,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_insert_clique rsb__wp_BLAS_zuscr_insert_clique
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zuscr_insert_clique( blas_sparse_matrix A, const int k, const int l,
                       const void *val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_clique(A,k,l,val,row_stride,col_stride,indx,jndx))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_insert_clique_ rsb__wp_blas_zuscr_insert_clique_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_insert_clique_( blas_sparse_matrix*A,const int*k,const int*l,const void *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zuscr_insert_clique(*A,*k,*l,val,*row_stride,*col_stride,indx,jndx );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_suscr_insert_block rsb__wp_BLAS_suscr_insert_block
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_suscr_insert_block( blas_sparse_matrix A, const float * val,
                        int row_stride, int col_stride, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_block(A,val,row_stride,col_stride,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_suscr_insert_block_ rsb__wp_blas_suscr_insert_block_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_suscr_insert_block_( blas_sparse_matrix*A,const float *val,int*row_stride,int*col_stride,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_suscr_insert_block(*A,val,*row_stride,*col_stride,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_duscr_insert_block rsb__wp_BLAS_duscr_insert_block
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_duscr_insert_block( blas_sparse_matrix A, const double * val,
                        int row_stride, int col_stride, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_block(A,val,row_stride,col_stride,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_duscr_insert_block_ rsb__wp_blas_duscr_insert_block_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_duscr_insert_block_( blas_sparse_matrix*A,const double *val,int*row_stride,int*col_stride,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_duscr_insert_block(*A,val,*row_stride,*col_stride,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cuscr_insert_block rsb__wp_BLAS_cuscr_insert_block
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cuscr_insert_block( blas_sparse_matrix A, const void *val,
                        int row_stride, int col_stride, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_block(A,val,row_stride,col_stride,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cuscr_insert_block_ rsb__wp_blas_cuscr_insert_block_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cuscr_insert_block_( blas_sparse_matrix*A,const void *val,int*row_stride,int*col_stride,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cuscr_insert_block(*A,val,*row_stride,*col_stride,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zuscr_insert_block rsb__wp_BLAS_zuscr_insert_block
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zuscr_insert_block( blas_sparse_matrix A, const void *val,
                        int row_stride, int col_stride, int i, int j )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_insert_block(A,val,row_stride,col_stride,i,j))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zuscr_insert_block_ rsb__wp_blas_zuscr_insert_block_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zuscr_insert_block_( blas_sparse_matrix*A,const void *val,int*row_stride,int*col_stride,int*i,int*j,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zuscr_insert_block(*A,val,*row_stride,*col_stride,*i,*j );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}



#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_uscr_end rsb__wp_BLAS_uscr_end
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_uscr_end( blas_sparse_matrix A )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xuscr_end(A))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_uscr_end_ rsb__wp_blas_uscr_end_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_uscr_end_( blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_cr_end_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_uscr_end(*A );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_usds rsb__wp_BLAS_usds
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_usds( blas_sparse_matrix A )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_ds_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusds(A))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_usds_ rsb__wp_blas_usds_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_usds_( blas_sparse_matrix*A,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2_ds_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_usds(*A );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susrows_scale rsb__wp_BLAS_susrows_scale
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susrows_scale( blas_sparse_matrix A,const float *  d, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusrows_scale(A,d,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susrows_scale_ rsb__wp_blas_susrows_scale_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susrows_scale_( blas_sparse_matrix*A,const float *d,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susrows_scale(*A,d,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusrows_scale rsb__wp_BLAS_dusrows_scale
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusrows_scale( blas_sparse_matrix A,const double *  d, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusrows_scale(A,d,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusrows_scale_ rsb__wp_blas_dusrows_scale_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusrows_scale_( blas_sparse_matrix*A,const double *d,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusrows_scale(*A,d,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusrows_scale rsb__wp_BLAS_cusrows_scale
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusrows_scale( blas_sparse_matrix A,const void * d, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusrows_scale(A,d,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusrows_scale_ rsb__wp_blas_cusrows_scale_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusrows_scale_( blas_sparse_matrix*A,const void *d,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusrows_scale(*A,d,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusrows_scale rsb__wp_BLAS_zusrows_scale
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusrows_scale( blas_sparse_matrix A,const void * d, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusrows_scale(A,d,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusrows_scale_ rsb__wp_blas_zusrows_scale_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusrows_scale_( blas_sparse_matrix*A,const void *d,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usrows_scale_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusrows_scale(*A,d,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susget_diag rsb__wp_BLAS_susget_diag
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susget_diag( blas_sparse_matrix A,float *  d )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_diag(A,d))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susget_diag_ rsb__wp_blas_susget_diag_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susget_diag_( blas_sparse_matrix*A,float *d,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susget_diag(*A,d );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusget_diag rsb__wp_BLAS_dusget_diag
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusget_diag( blas_sparse_matrix A,double *  d )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_diag(A,d))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusget_diag_ rsb__wp_blas_dusget_diag_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusget_diag_( blas_sparse_matrix*A,double *d,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusget_diag(*A,d );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusget_diag rsb__wp_BLAS_cusget_diag
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusget_diag( blas_sparse_matrix A,void * d )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_diag(A,d))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusget_diag_ rsb__wp_blas_cusget_diag_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusget_diag_( blas_sparse_matrix*A,void *d,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusget_diag(*A,d );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusget_diag rsb__wp_BLAS_zusget_diag
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusget_diag( blas_sparse_matrix A,void * d )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_diag(A,d))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusget_diag_ rsb__wp_blas_zusget_diag_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusget_diag_( blas_sparse_matrix*A,void *d,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_diag_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusget_diag(*A,d );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susget_rows_nnz rsb__wp_BLAS_susget_rows_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susget_rows_nnz( blas_sparse_matrix A, const int fr, const int lr, int * nnzp )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_nnz(A,fr,lr,nnzp))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susget_rows_nnz_ rsb__wp_blas_susget_rows_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susget_rows_nnz_( blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susget_rows_nnz(*A,*fr,*lr,nnzp );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusget_rows_nnz rsb__wp_BLAS_dusget_rows_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusget_rows_nnz( blas_sparse_matrix A, const int fr, const int lr, int * nnzp )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_nnz(A,fr,lr,nnzp))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusget_rows_nnz_ rsb__wp_blas_dusget_rows_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusget_rows_nnz_( blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusget_rows_nnz(*A,*fr,*lr,nnzp );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusget_rows_nnz rsb__wp_BLAS_cusget_rows_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusget_rows_nnz( blas_sparse_matrix A, const int fr, const int lr, int * nnzp )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_nnz(A,fr,lr,nnzp))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusget_rows_nnz_ rsb__wp_blas_cusget_rows_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusget_rows_nnz_( blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusget_rows_nnz(*A,*fr,*lr,nnzp );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusget_rows_nnz rsb__wp_BLAS_zusget_rows_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusget_rows_nnz( blas_sparse_matrix A, const int fr, const int lr, int * nnzp )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_nnz(A,fr,lr,nnzp))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusget_rows_nnz_ rsb__wp_blas_zusget_rows_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusget_rows_nnz_( blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusget_rows_nnz(*A,*fr,*lr,nnzp );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susget_rows_sparse rsb__wp_BLAS_susget_rows_sparse
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susget_rows_sparse( blas_sparse_matrix A, float *  VA, int * IA, int * JA, int * nnz, const int fr, const int lr )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_sparse(A,VA,IA,JA,nnz,fr,lr))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susget_rows_sparse_ rsb__wp_blas_susget_rows_sparse_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susget_rows_sparse_( blas_sparse_matrix*A,float *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susget_rows_sparse(*A,VA,IA,JA,nnz,*fr,*lr );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusget_rows_sparse rsb__wp_BLAS_dusget_rows_sparse
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusget_rows_sparse( blas_sparse_matrix A, double *  VA, int * IA, int * JA, int * nnz, const int fr, const int lr )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_sparse(A,VA,IA,JA,nnz,fr,lr))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusget_rows_sparse_ rsb__wp_blas_dusget_rows_sparse_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusget_rows_sparse_( blas_sparse_matrix*A,double *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusget_rows_sparse(*A,VA,IA,JA,nnz,*fr,*lr );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusget_rows_sparse rsb__wp_BLAS_cusget_rows_sparse
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusget_rows_sparse( blas_sparse_matrix A, void * VA, int * IA, int * JA, int * nnz, const int fr, const int lr )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_sparse(A,VA,IA,JA,nnz,fr,lr))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusget_rows_sparse_ rsb__wp_blas_cusget_rows_sparse_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusget_rows_sparse_( blas_sparse_matrix*A,void *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusget_rows_sparse(*A,VA,IA,JA,nnz,*fr,*lr );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusget_rows_sparse rsb__wp_BLAS_zusget_rows_sparse
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusget_rows_sparse( blas_sparse_matrix A, void * VA, int * IA, int * JA, int * nnz, const int fr, const int lr )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_rows_sparse(A,VA,IA,JA,nnz,fr,lr))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusget_rows_sparse_ rsb__wp_blas_zusget_rows_sparse_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusget_rows_sparse_( blas_sparse_matrix*A,void *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_rows_sparse_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusget_rows_sparse(*A,VA,IA,JA,nnz,*fr,*lr );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susget_matrix_nnz rsb__wp_BLAS_susget_matrix_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susget_matrix_nnz( blas_sparse_matrix A,int * nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_matrix_nnz(A,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susget_matrix_nnz_ rsb__wp_blas_susget_matrix_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susget_matrix_nnz_( blas_sparse_matrix*A,int *nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susget_matrix_nnz(*A,nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusget_matrix_nnz rsb__wp_BLAS_dusget_matrix_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusget_matrix_nnz( blas_sparse_matrix A,int * nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_matrix_nnz(A,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusget_matrix_nnz_ rsb__wp_blas_dusget_matrix_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusget_matrix_nnz_( blas_sparse_matrix*A,int *nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusget_matrix_nnz(*A,nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusget_matrix_nnz rsb__wp_BLAS_cusget_matrix_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusget_matrix_nnz( blas_sparse_matrix A,int * nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_matrix_nnz(A,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusget_matrix_nnz_ rsb__wp_blas_cusget_matrix_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusget_matrix_nnz_( blas_sparse_matrix*A,int *nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusget_matrix_nnz(*A,nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusget_matrix_nnz rsb__wp_BLAS_zusget_matrix_nnz
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusget_matrix_nnz( blas_sparse_matrix A,int * nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_matrix_nnz(A,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusget_matrix_nnz_ rsb__wp_blas_zusget_matrix_nnz_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusget_matrix_nnz_( blas_sparse_matrix*A,int *nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_matrix_nnz_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusget_matrix_nnz(*A,nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susget_infinity_norm rsb__wp_BLAS_susget_infinity_norm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susget_infinity_norm( blas_sparse_matrix A,float * in, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_infinity_norm(A,in,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susget_infinity_norm_ rsb__wp_blas_susget_infinity_norm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susget_infinity_norm_( blas_sparse_matrix*A,float *in,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susget_infinity_norm(*A,in,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusget_infinity_norm rsb__wp_BLAS_dusget_infinity_norm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusget_infinity_norm( blas_sparse_matrix A,double * in, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_infinity_norm(A,in,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusget_infinity_norm_ rsb__wp_blas_dusget_infinity_norm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusget_infinity_norm_( blas_sparse_matrix*A,double *in,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusget_infinity_norm(*A,in,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusget_infinity_norm rsb__wp_BLAS_cusget_infinity_norm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusget_infinity_norm( blas_sparse_matrix A,void *in, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_infinity_norm(A,in,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusget_infinity_norm_ rsb__wp_blas_cusget_infinity_norm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusget_infinity_norm_( blas_sparse_matrix*A,void *in,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusget_infinity_norm(*A,in,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusget_infinity_norm rsb__wp_BLAS_zusget_infinity_norm
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusget_infinity_norm( blas_sparse_matrix A,void *in, const enum blas_trans_type trans )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_infinity_norm(A,in,trans))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusget_infinity_norm_ rsb__wp_blas_zusget_infinity_norm_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusget_infinity_norm_( blas_sparse_matrix*A,void *in,const enum blas_trans_type*trans,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_infinity_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusget_infinity_norm(*A,in,*trans );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susset_elements rsb__wp_BLAS_susset_elements
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susset_elements( blas_sparse_matrix A,const int * ia, const int *ja, const float *  va, const int nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_elements(A,ia,ja,va,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susset_elements_ rsb__wp_blas_susset_elements_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susset_elements_( blas_sparse_matrix*A,const int *ia,const int *ja,const float *va,const int*nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susset_elements(*A,ia,ja,va,*nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusset_elements rsb__wp_BLAS_dusset_elements
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusset_elements( blas_sparse_matrix A,const int * ia, const int *ja, const double *  va, const int nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_elements(A,ia,ja,va,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusset_elements_ rsb__wp_blas_dusset_elements_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusset_elements_( blas_sparse_matrix*A,const int *ia,const int *ja,const double *va,const int*nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusset_elements(*A,ia,ja,va,*nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusset_elements rsb__wp_BLAS_cusset_elements
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusset_elements( blas_sparse_matrix A,const int * ia, const int *ja, const void * va, const int nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_elements(A,ia,ja,va,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusset_elements_ rsb__wp_blas_cusset_elements_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusset_elements_( blas_sparse_matrix*A,const int *ia,const int *ja,const void *va,const int*nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusset_elements(*A,ia,ja,va,*nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusset_elements rsb__wp_BLAS_zusset_elements
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusset_elements( blas_sparse_matrix A,const int * ia, const int *ja, const void * va, const int nnz )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_elements(A,ia,ja,va,nnz))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusset_elements_ rsb__wp_blas_zusset_elements_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusset_elements_( blas_sparse_matrix*A,const int *ia,const int *ja,const void *va,const int*nnz,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_elements_norm_msg.\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusset_elements(*A,ia,ja,va,*nnz );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susset_element rsb__wp_BLAS_susset_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susset_element( blas_sparse_matrix A,const int i, const int j, float *  v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susset_element_ rsb__wp_blas_susset_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susset_element_( blas_sparse_matrix*A,const int*i,const int*j,float *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susset_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusset_element rsb__wp_BLAS_dusset_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusset_element( blas_sparse_matrix A,const int i, const int j, double *  v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusset_element_ rsb__wp_blas_dusset_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusset_element_( blas_sparse_matrix*A,const int*i,const int*j,double *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusset_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusset_element rsb__wp_BLAS_cusset_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusset_element( blas_sparse_matrix A,const int i, const int j, void * v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusset_element_ rsb__wp_blas_cusset_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusset_element_( blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusset_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusset_element rsb__wp_BLAS_zusset_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusset_element( blas_sparse_matrix A,const int i, const int j, void * v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusset_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusset_element_ rsb__wp_blas_zusset_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusset_element_( blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usset_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusset_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_susget_element rsb__wp_BLAS_susget_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_susget_element( blas_sparse_matrix A,const int i, const int j, float *  v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_susget_element_ rsb__wp_blas_susget_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_susget_element_( blas_sparse_matrix*A,const int*i,const int*j,float *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_susget_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_dusget_element rsb__wp_BLAS_dusget_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_dusget_element( blas_sparse_matrix A,const int i, const int j, double *  v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_dusget_element_ rsb__wp_blas_dusget_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_dusget_element_( blas_sparse_matrix*A,const int*i,const int*j,double *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_dusget_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_cusget_element rsb__wp_BLAS_cusget_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_cusget_element( blas_sparse_matrix A,const int i, const int j, void * v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_cusget_element_ rsb__wp_blas_cusget_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_cusget_element_( blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_cusget_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define BLAS_zusget_element rsb__wp_BLAS_zusget_element
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
int BLAS_zusget_element( blas_sparse_matrix A,const int i, const int j, void * v )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_return_msg
         */
	RSB_SPB_INTERFACE_PREAMBLE

	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusget_element(A,i,j,v))
}
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define blas_zusget_element_ rsb__wp_blas_zusget_element_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_zusget_element_( blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat )
{
         /*!
           \ingroup rsb_doc_sparse_blas
           \rsb_spblasl2e_usget_element_norm_msg\rsb_spblas_istat_msg
         */

	int istatv = BLAS_zusget_element(*A,*i,*j,v );
	RSB_SET_IF_NOT_NULL(istat,istatv);
}





int BLAS_usgp( blas_sparse_matrix A, rsb_blas_pname_t pname ) /*  FIXME: temporarily here */
{
	/**
	 \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_gp_msg
	 \rsb_spblas_return_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	RSB_SPB_INTERFACE_RETURN_EXP(rsb__BLAS_usgp(A,pname))
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define  blas_usgp_	rsb__wp_blas_usgp_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_usgp_( blas_sparse_matrix*A, rsb_blas_pname_t *pname, int * istat ) /*  FIXME: temporarily here */
{
	/** \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_gp_msg
	 \rsb_spblas_istat_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	*istat=BLAS_usgp(*A,*pname); /*  FIXME: temporarily here */
	RSB_SPB_INTERFACE_RETURN_VOID()
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
#define  blas_ussp_	rsb__wp_blas_ussp_
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
void blas_ussp_( blas_sparse_matrix*A, rsb_blas_pname_t *pname, int * istat ) /*  FIXME: temporarily here */
{
	/**
	 \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_sp_msg
	 \rsb_spblas_istat_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	*istat=BLAS_ussp(*A,*pname); /*  FIXME: temporarily here */
	RSB_SPB_INTERFACE_RETURN_VOID()
}

int BLAS_ussp( blas_sparse_matrix A, rsb_blas_pname_t pname ) /*  FIXME: temporarily here */
{
	/**
	 \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_sp_msg
	 \rsb_spblas_return_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_ussp(A,pname))
}



blas_sparse_matrix rsb_blas_file_mtx_load(const rsb_char_t * filename, rsb_type_t typecode )
{
	/*!
 	\ingroup rsb_doc_sparse_blas
	Load Matrix Market matrix file of specified type to a matrix, and return Sparse BLAS handler.

	\param \rsb_filename_inp_param_msg
	\param \rsb_type_o_param_msg
	\return \rsb_spblasl2_A_msg_ftn_rt

	\n
	 */
	RSB_SPB_INTERFACE_PREAMBLE
{
	RSB_SPB_INTERFACE_RETURN_EXP( rsb__load_spblas_matrix_file_as_matrix_market(filename,typecode) )
}
}

struct rsb_mtx_t * rsb_blas_get_mtx(blas_sparse_matrix A)
{
	/*!
 	\ingroup rsb_doc_sparse_blas
	\rsb_BLAS_get_mtx_msg

	\rsb_spblasl2_A_msg
	\return \rsbmtxpmessage_bg

	\n
	
	\rsb_BLAS_get_mtx_example
	\see_rsb_BLAS_get_mtx_msg
	\rsb_BLAS_get_mtx_msg_todo
	\rsb_BLAS_get_mtx_msg_note
	\rsb_BLAS_get_mtx_msg_warn
	 */
	RSB_SPB_INTERFACE_PREAMBLE
{
	struct rsb_mtx_t * mtxAp = NULL;
	mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	RSB_SPB_INTERFACE_RETURN_EXP( mtxAp )
}
}

#ifdef __cplusplus
}
#endif  /* __cplusplus */



