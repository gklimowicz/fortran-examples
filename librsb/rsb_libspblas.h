

/*!
        @file
        @author Michele Martone

	@brief  This file specifies the Sparse BLAS interface to librsb.
	Supported types  :(float,double,float complex,double complex) .
	Unsupported types:() .
	Level 1 ops      :(dot,axpy,ga,gz,sc) .
	Level 2 ops      :(mv,sv) .
	Level 3 ops      :(mm,sm) .
*/

#ifndef RSB_LIBSPBLAS_H_INCLUDED
#define RSB_LIBSPBLAS_H_INCLUDED
#ifndef RSB_RSB_H_INCLUDED
#error "You are using Sparse BLAS headers from librsb -- You should include <rsb.h> first!"
#endif /* RSB_RSB_H_INCLUDED */
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

/** the sparse matrix descriptor type */
typedef int blas_sparse_matrix;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

               /* Level 1 Computational Routines */
int BLAS_susdot(const enum blas_conj_type conj, const int nnz, const float * x,
		const int *indx, const float * y, const int incy, float * r,
		const enum blas_base_type index_base);
void blas_susdot_(const enum blas_conj_type*conj,const int*nnz,const float *x,const int *indx,const float *y,const int*incy,float *r,const enum blas_base_type*index_base,int*istat);
int BLAS_dusdot(const enum blas_conj_type conj, const int nnz, const double * x,
		const int *indx, const double * y, const int incy, double * r,
		const enum blas_base_type index_base);
void blas_dusdot_(const enum blas_conj_type*conj,const int*nnz,const double *x,const int *indx,const double *y,const int*incy,double *r,const enum blas_base_type*index_base,int*istat);
int BLAS_cusdot(const enum blas_conj_type conj, const int nnz, const void *x,
		const int *indx, const void *y, const int incy, void *r,
		const enum blas_base_type index_base);
void blas_cusdot_(const enum blas_conj_type*conj,const int*nnz,const void *x,const int *indx,const void *y,const int*incy,void *r,const enum blas_base_type*index_base,int*istat);
int BLAS_zusdot(const enum blas_conj_type conj, const int nnz, const void *x,
		const int *indx, const void *y, const int incy, void *r,
		const enum blas_base_type index_base);
void blas_zusdot_(const enum blas_conj_type*conj,const int*nnz,const void *x,const int *indx,const void *y,const int*incy,void *r,const enum blas_base_type*index_base,int*istat);

int BLAS_susaxpy(const int nnz, float  alpha, const float * x, const int *indx,
                 float * y, const int incy, const enum blas_base_type index_base);
void blas_susaxpy_(const int*nnz,float*alpha,const float *x,const int *indx,float *y,const int*incy,const enum blas_base_type*index_base,int*istat);
int BLAS_dusaxpy(const int nnz, double  alpha, const double * x, const int *indx,
                 double * y, const int incy, const enum blas_base_type index_base);
void blas_dusaxpy_(const int*nnz,double*alpha,const double *x,const int *indx,double *y,const int*incy,const enum blas_base_type*index_base,int*istat);
int BLAS_cusaxpy(const int nnz, const void * alpha, const void *x, const int *indx,
                 void *y, const int incy, const enum blas_base_type index_base);
void blas_cusaxpy_(const int*nnz,const void *alpha,const void *x,const int *indx,void *y,const int*incy,const enum blas_base_type*index_base,int*istat);
int BLAS_zusaxpy(const int nnz, const void * alpha, const void *x, const int *indx,
                 void *y, const int incy, const enum blas_base_type index_base);
void blas_zusaxpy_(const int*nnz,const void *alpha,const void *x,const int *indx,void *y,const int*incy,const enum blas_base_type*index_base,int*istat);

int BLAS_susga(const int nnz, const float * y, const int incy, float * x, const int *indx,
              const enum blas_base_type index_base);
void blas_susga_(const int*nnz,const float *y,const int*incy,float *x,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_dusga(const int nnz, const double * y, const int incy, double * x, const int *indx,
              const enum blas_base_type index_base);
void blas_dusga_(const int*nnz,const double *y,const int*incy,double *x,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_cusga(const int nnz, const void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base);
void blas_cusga_(const int*nnz,const void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_zusga(const int nnz, const void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base);
void blas_zusga_(const int*nnz,const void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat);

int BLAS_susgz(const int nnz, float * y, const int incy, float * x, const int *indx,
              const enum blas_base_type index_base);
void blas_susgz_(const int*nnz,float *y,const int*incy,float *x,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_dusgz(const int nnz, double * y, const int incy, double * x, const int *indx,
              const enum blas_base_type index_base);
void blas_dusgz_(const int*nnz,double *y,const int*incy,double *x,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_cusgz(const int nnz, void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base);
void blas_cusgz_(const int*nnz,void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_zusgz(const int nnz, void *y, const int incy, void *x, const int *indx,
              const enum blas_base_type index_base);
void blas_zusgz_(const int*nnz,void *y,const int*incy,void *x,const int *indx,const enum blas_base_type*index_base,int*istat);

int BLAS_sussc(const int nnz, const float * x, float * y, const int incy, const int *indx,
              const enum blas_base_type index_base);
void blas_sussc_(const int*nnz,const float *x,float *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_dussc(const int nnz, const double * x, double * y, const int incy, const int *indx,
              const enum blas_base_type index_base);
void blas_dussc_(const int*nnz,const double *x,double *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_cussc(const int nnz, const void *x, void *y, const int incy, const int *indx,
              const enum blas_base_type index_base);
void blas_cussc_(const int*nnz,const void *x,void *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat);
int BLAS_zussc(const int nnz, const void *x, void *y, const int incy, const int *indx,
              const enum blas_base_type index_base);
void blas_zussc_(const int*nnz,const void *x,void *y,const int*incy,const int *indx,const enum blas_base_type*index_base,int*istat);

               /* Level 2 Computational Routines */

int BLAS_susmv(const enum blas_trans_type transA, float alpha,
    const blas_sparse_matrix A, const float * x, const int incx, float * y, const int incy);

void blas_susmv_(const enum blas_trans_type*transA,float*alpha,const blas_sparse_matrix*A,const float *x,const int*incx,float *y,const int*incy,int*istat);

int BLAS_dusmv(const enum blas_trans_type transA, double alpha,
    const blas_sparse_matrix A, const double * x, const int incx, double * y, const int incy);

void blas_dusmv_(const enum blas_trans_type*transA,double*alpha,const blas_sparse_matrix*A,const double *x,const int*incx,double *y,const int*incy,int*istat);

int BLAS_cusmv(const enum blas_trans_type transA, const void *alpha,
    const blas_sparse_matrix A, const void *x, const int incx, void *y, const int incy);

void blas_cusmv_(const enum blas_trans_type*transA,const void *alpha,const blas_sparse_matrix*A,const void *x,const int*incx,void *y,const int*incy,int*istat);

int BLAS_zusmv(const enum blas_trans_type transA, const void *alpha,
    const blas_sparse_matrix A, const void *x, const int incx, void *y, const int incy);

void blas_zusmv_(const enum blas_trans_type*transA,const void *alpha,const blas_sparse_matrix*A,const void *x,const int*incx,void *y,const int*incy,int*istat);


int BLAS_sussv(enum blas_trans_type transT, float alpha,
    const blas_sparse_matrix T, float * x, const int incx);

void blas_sussv_(enum blas_trans_type*transT,float*alpha,const blas_sparse_matrix*T,float *x,const int*incx,int*istat);

int BLAS_dussv(enum blas_trans_type transT, double alpha,
    const blas_sparse_matrix T, double * x, const int incx);

void blas_dussv_(enum blas_trans_type*transT,double*alpha,const blas_sparse_matrix*T,double *x,const int*incx,int*istat);

int BLAS_cussv(enum blas_trans_type transT, const void *alpha,
    const blas_sparse_matrix T, void *x, const int incx);

void blas_cussv_(enum blas_trans_type*transT,const void *alpha,const blas_sparse_matrix*T,void *x,const int*incx,int*istat);

int BLAS_zussv(enum blas_trans_type transT, const void *alpha,
    const blas_sparse_matrix T, void *x, const int incx);

void blas_zussv_(enum blas_trans_type*transT,const void *alpha,const blas_sparse_matrix*T,void *x,const int*incx,int*istat);


               /* Level 3 Computational Routines */

int BLAS_susmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, float alpha, const blas_sparse_matrix A, const float * b, const int ldb,
       float *  c, const int ldc);

void blas_susmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,float*alpha,const blas_sparse_matrix*A,const float *b,const int*ldb,float *c,const int*ldc,int*istat);

int BLAS_dusmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, double alpha, const blas_sparse_matrix A, const double * b, const int ldb,
       double *  c, const int ldc);

void blas_dusmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,double*alpha,const blas_sparse_matrix*A,const double *b,const int*ldb,double *c,const int*ldc,int*istat);

int BLAS_cusmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, const void *alpha, const blas_sparse_matrix A, const void *b, const int ldb,
       void * c, const int ldc);

void blas_cusmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,const void *alpha,const blas_sparse_matrix*A,const void *b,const int*ldb,void *c,const int*ldc,int*istat);

int BLAS_zusmm(const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, const void *alpha, const blas_sparse_matrix A, const void *b, const int ldb,
       void * c, const int ldc);

void blas_zusmm_(const enum blas_order_type*order,const enum blas_trans_type*transA,const int*nrhs,const void *alpha,const blas_sparse_matrix*A,const void *b,const int*ldb,void *c,const int*ldc,int*istat);


int BLAS_sussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, float alpha, const blas_sparse_matrix T, float * b, const int ldb);

void blas_sussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,float*alpha,const blas_sparse_matrix*T,float *b,const int*ldb,int*istat);

int BLAS_dussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, double alpha, const blas_sparse_matrix T, double * b, const int ldb);

void blas_dussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,double*alpha,const blas_sparse_matrix*T,double *b,const int*ldb,int*istat);

int BLAS_cussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, const void *alpha, const blas_sparse_matrix T, void *b, const int ldb);

void blas_cussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,const void *alpha,const blas_sparse_matrix*T,void *b,const int*ldb,int*istat);

int BLAS_zussm(const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, const void *alpha, const blas_sparse_matrix T, void *b, const int ldb);

void blas_zussm_(const enum blas_order_type*order,const enum blas_trans_type*transT,const int*nrhs,const void *alpha,const blas_sparse_matrix*T,void *b,const int*ldb,int*istat);


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

blas_sparse_matrix BLAS_suscr_begin(int m, int n);
void blas_suscr_begin_(int*m,int*n,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_duscr_begin(int m, int n);
void blas_duscr_begin_(int*m,int*n,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_cuscr_begin(int m, int n);
void blas_cuscr_begin_(int*m,int*n,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_zuscr_begin(int m, int n);
void blas_zuscr_begin_(int*m,int*n,blas_sparse_matrix*A,int*istat);

blas_sparse_matrix BLAS_suscr_block_begin(int Mb, int Nb, int k, int l);
void blas_suscr_block_begin_(int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_duscr_block_begin(int Mb, int Nb, int k, int l);
void blas_duscr_block_begin_(int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_cuscr_block_begin(int Mb, int Nb, int k, int l);
void blas_cuscr_block_begin_(int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_zuscr_block_begin(int Mb, int Nb, int k, int l);
void blas_zuscr_block_begin_(int*Mb,int*Nb,int*k,int*l,blas_sparse_matrix*A,int*istat);

blas_sparse_matrix BLAS_suscr_variable_block_begin(int Mb, int Nb,
		const int *K, const int *L);
void blas_suscr_variable_block_begin_(int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_duscr_variable_block_begin(int Mb, int Nb,
		const int *K, const int *L);
void blas_duscr_variable_block_begin_(int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_cuscr_variable_block_begin(int Mb, int Nb,
		const int *K, const int *L);
void blas_cuscr_variable_block_begin_(int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat);
blas_sparse_matrix BLAS_zuscr_variable_block_begin(int Mb, int Nb,
		const int *K, const int *L);
void blas_zuscr_variable_block_begin_(int*Mb,int*Nb,const int *K,const int *L,blas_sparse_matrix*A,int*istat);

int BLAS_suscr_end(blas_sparse_matrix A);
void blas_suscr_end_(blas_sparse_matrix*A,int*istat);
int BLAS_duscr_end(blas_sparse_matrix A);
void blas_duscr_end_(blas_sparse_matrix*A,int*istat);
int BLAS_cuscr_end(blas_sparse_matrix A);
void blas_cuscr_end_(blas_sparse_matrix*A,int*istat);
int BLAS_zuscr_end(blas_sparse_matrix A);
void blas_zuscr_end_(blas_sparse_matrix*A,int*istat);

int BLAS_suscr_insert_entry(blas_sparse_matrix A, float  val, int i, int j);
void blas_suscr_insert_entry_(blas_sparse_matrix*A,float*val,int*i,int*j,int*istat);
int BLAS_duscr_insert_entry(blas_sparse_matrix A, double  val, int i, int j);
void blas_duscr_insert_entry_(blas_sparse_matrix*A,double*val,int*i,int*j,int*istat);
int BLAS_cuscr_insert_entry(blas_sparse_matrix A, const void * val, int i, int j);
void blas_cuscr_insert_entry_(blas_sparse_matrix*A,const void *val,int*i,int*j,int*istat);
int BLAS_zuscr_insert_entry(blas_sparse_matrix A, const void * val, int i, int j);
void blas_zuscr_insert_entry_(blas_sparse_matrix*A,const void *val,int*i,int*j,int*istat);

int BLAS_suscr_insert_entries(blas_sparse_matrix A, int nnz, const float * val,
                            const int *indx, const int *jndx);
void blas_suscr_insert_entries_(blas_sparse_matrix*A,int*nnz,const float *val,const int *indx,const int *jndx,int*istat);
int BLAS_duscr_insert_entries(blas_sparse_matrix A, int nnz, const double * val,
                            const int *indx, const int *jndx);
void blas_duscr_insert_entries_(blas_sparse_matrix*A,int*nnz,const double *val,const int *indx,const int *jndx,int*istat);
int BLAS_cuscr_insert_entries(blas_sparse_matrix A, int nnz, const void *val,
                            const int *indx, const int *jndx);
void blas_cuscr_insert_entries_(blas_sparse_matrix*A,int*nnz,const void *val,const int *indx,const int *jndx,int*istat);
int BLAS_zuscr_insert_entries(blas_sparse_matrix A, int nnz, const void *val,
                            const int *indx, const int *jndx);
void blas_zuscr_insert_entries_(blas_sparse_matrix*A,int*nnz,const void *val,const int *indx,const int *jndx,int*istat);

int BLAS_suscr_insert_col(blas_sparse_matrix A, int j, int nnz,
                           const float * val, const int *indx);
void blas_suscr_insert_col_(blas_sparse_matrix*A,int*j,int*nnz,const float *val,const int *indx,int*istat);
int BLAS_duscr_insert_col(blas_sparse_matrix A, int j, int nnz,
                           const double * val, const int *indx);
void blas_duscr_insert_col_(blas_sparse_matrix*A,int*j,int*nnz,const double *val,const int *indx,int*istat);
int BLAS_cuscr_insert_col(blas_sparse_matrix A, int j, int nnz,
                           const void *val, const int *indx);
void blas_cuscr_insert_col_(blas_sparse_matrix*A,int*j,int*nnz,const void *val,const int *indx,int*istat);
int BLAS_zuscr_insert_col(blas_sparse_matrix A, int j, int nnz,
                           const void *val, const int *indx);
void blas_zuscr_insert_col_(blas_sparse_matrix*A,int*j,int*nnz,const void *val,const int *indx,int*istat);

int BLAS_suscr_insert_row(blas_sparse_matrix A, int i, int nnz,
                           const float * val, const int *indx);
void blas_suscr_insert_row_(blas_sparse_matrix*A,int*i,int*nnz,const float *val,const int *indx,int*istat);
int BLAS_duscr_insert_row(blas_sparse_matrix A, int i, int nnz,
                           const double * val, const int *indx);
void blas_duscr_insert_row_(blas_sparse_matrix*A,int*i,int*nnz,const double *val,const int *indx,int*istat);
int BLAS_cuscr_insert_row(blas_sparse_matrix A, int i, int nnz,
                           const void *val, const int *indx);
void blas_cuscr_insert_row_(blas_sparse_matrix*A,int*i,int*nnz,const void *val,const int *indx,int*istat);
int BLAS_zuscr_insert_row(blas_sparse_matrix A, int i, int nnz,
                           const void *val, const int *indx);
void blas_zuscr_insert_row_(blas_sparse_matrix*A,int*i,int*nnz,const void *val,const int *indx,int*istat);

int BLAS_suscr_insert_clique(blas_sparse_matrix A, const int k, const int l,
                       const float * val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx);
void blas_suscr_insert_clique_(blas_sparse_matrix*A,const int*k,const int*l,const float *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat);
int BLAS_duscr_insert_clique(blas_sparse_matrix A, const int k, const int l,
                       const double * val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx);
void blas_duscr_insert_clique_(blas_sparse_matrix*A,const int*k,const int*l,const double *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat);
int BLAS_cuscr_insert_clique(blas_sparse_matrix A, const int k, const int l,
                       const void *val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx);
void blas_cuscr_insert_clique_(blas_sparse_matrix*A,const int*k,const int*l,const void *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat);
int BLAS_zuscr_insert_clique(blas_sparse_matrix A, const int k, const int l,
                       const void *val, const int row_stride,
                       const int col_stride, const int *indx,
                       const int *jndx);
void blas_zuscr_insert_clique_(blas_sparse_matrix*A,const int*k,const int*l,const void *val,const int*row_stride,const int*col_stride,const int *indx,const int *jndx,int*istat);

int BLAS_suscr_insert_block(blas_sparse_matrix A, const float * val,
                        int row_stride, int col_stride, int i, int j);
void blas_suscr_insert_block_(blas_sparse_matrix*A,const float *val,int*row_stride,int*col_stride,int*i,int*j,int*istat);
int BLAS_duscr_insert_block(blas_sparse_matrix A, const double * val,
                        int row_stride, int col_stride, int i, int j);
void blas_duscr_insert_block_(blas_sparse_matrix*A,const double *val,int*row_stride,int*col_stride,int*i,int*j,int*istat);
int BLAS_cuscr_insert_block(blas_sparse_matrix A, const void *val,
                        int row_stride, int col_stride, int i, int j);
void blas_cuscr_insert_block_(blas_sparse_matrix*A,const void *val,int*row_stride,int*col_stride,int*i,int*j,int*istat);
int BLAS_zuscr_insert_block(blas_sparse_matrix A, const void *val,
                        int row_stride, int col_stride, int i, int j);
void blas_zuscr_insert_block_(blas_sparse_matrix*A,const void *val,int*row_stride,int*col_stride,int*i,int*j,int*istat);



int BLAS_uscr_end(blas_sparse_matrix A);
void blas_uscr_end_(blas_sparse_matrix*A,int*istat);
int BLAS_usds(blas_sparse_matrix A);
void blas_usds_(blas_sparse_matrix*A,int*istat);

int BLAS_susrows_scale(blas_sparse_matrix A,const float *  d, const enum blas_trans_type trans);
void blas_susrows_scale_(blas_sparse_matrix*A,const float *d,const enum blas_trans_type*trans,int*istat);
int BLAS_dusrows_scale(blas_sparse_matrix A,const double *  d, const enum blas_trans_type trans);
void blas_dusrows_scale_(blas_sparse_matrix*A,const double *d,const enum blas_trans_type*trans,int*istat);
int BLAS_cusrows_scale(blas_sparse_matrix A,const void * d, const enum blas_trans_type trans);
void blas_cusrows_scale_(blas_sparse_matrix*A,const void *d,const enum blas_trans_type*trans,int*istat);
int BLAS_zusrows_scale(blas_sparse_matrix A,const void * d, const enum blas_trans_type trans);
void blas_zusrows_scale_(blas_sparse_matrix*A,const void *d,const enum blas_trans_type*trans,int*istat);

int BLAS_susget_diag(blas_sparse_matrix A,float *  d);
void blas_susget_diag_(blas_sparse_matrix*A,float *d,int*istat);
int BLAS_dusget_diag(blas_sparse_matrix A,double *  d);
void blas_dusget_diag_(blas_sparse_matrix*A,double *d,int*istat);
int BLAS_cusget_diag(blas_sparse_matrix A,void * d);
void blas_cusget_diag_(blas_sparse_matrix*A,void *d,int*istat);
int BLAS_zusget_diag(blas_sparse_matrix A,void * d);
void blas_zusget_diag_(blas_sparse_matrix*A,void *d,int*istat);

int BLAS_susget_rows_nnz(blas_sparse_matrix A, const int fr, const int lr, int * nnzp);
void blas_susget_rows_nnz_(blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat);
int BLAS_dusget_rows_nnz(blas_sparse_matrix A, const int fr, const int lr, int * nnzp);
void blas_dusget_rows_nnz_(blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat);
int BLAS_cusget_rows_nnz(blas_sparse_matrix A, const int fr, const int lr, int * nnzp);
void blas_cusget_rows_nnz_(blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat);
int BLAS_zusget_rows_nnz(blas_sparse_matrix A, const int fr, const int lr, int * nnzp);
void blas_zusget_rows_nnz_(blas_sparse_matrix*A,const int*fr,const int*lr,int *nnzp,int*istat);

int BLAS_susget_rows_sparse(blas_sparse_matrix A, float *  VA, int * IA, int * JA, int * nnz, const int fr, const int lr);
void blas_susget_rows_sparse_(blas_sparse_matrix*A,float *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat);
int BLAS_dusget_rows_sparse(blas_sparse_matrix A, double *  VA, int * IA, int * JA, int * nnz, const int fr, const int lr);
void blas_dusget_rows_sparse_(blas_sparse_matrix*A,double *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat);
int BLAS_cusget_rows_sparse(blas_sparse_matrix A, void * VA, int * IA, int * JA, int * nnz, const int fr, const int lr);
void blas_cusget_rows_sparse_(blas_sparse_matrix*A,void *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat);
int BLAS_zusget_rows_sparse(blas_sparse_matrix A, void * VA, int * IA, int * JA, int * nnz, const int fr, const int lr);
void blas_zusget_rows_sparse_(blas_sparse_matrix*A,void *VA,int *IA,int *JA,int *nnz,const int*fr,const int*lr,int*istat);

int BLAS_susget_matrix_nnz(blas_sparse_matrix A,int * nnz);
void blas_susget_matrix_nnz_(blas_sparse_matrix*A,int *nnz,int*istat);
int BLAS_dusget_matrix_nnz(blas_sparse_matrix A,int * nnz);
void blas_dusget_matrix_nnz_(blas_sparse_matrix*A,int *nnz,int*istat);
int BLAS_cusget_matrix_nnz(blas_sparse_matrix A,int * nnz);
void blas_cusget_matrix_nnz_(blas_sparse_matrix*A,int *nnz,int*istat);
int BLAS_zusget_matrix_nnz(blas_sparse_matrix A,int * nnz);
void blas_zusget_matrix_nnz_(blas_sparse_matrix*A,int *nnz,int*istat);

int BLAS_susget_infinity_norm(blas_sparse_matrix A,float * in, const enum blas_trans_type trans);
void blas_susget_infinity_norm_(blas_sparse_matrix*A,float *in,const enum blas_trans_type*trans,int*istat);
int BLAS_dusget_infinity_norm(blas_sparse_matrix A,double * in, const enum blas_trans_type trans);
void blas_dusget_infinity_norm_(blas_sparse_matrix*A,double *in,const enum blas_trans_type*trans,int*istat);
int BLAS_cusget_infinity_norm(blas_sparse_matrix A,void *in, const enum blas_trans_type trans);
void blas_cusget_infinity_norm_(blas_sparse_matrix*A,void *in,const enum blas_trans_type*trans,int*istat);
int BLAS_zusget_infinity_norm(blas_sparse_matrix A,void *in, const enum blas_trans_type trans);
void blas_zusget_infinity_norm_(blas_sparse_matrix*A,void *in,const enum blas_trans_type*trans,int*istat);

int BLAS_susset_elements(blas_sparse_matrix A,const int * ia, const int *ja, const float *  va, const int nnz);
void blas_susset_elements_(blas_sparse_matrix*A,const int *ia,const int *ja,const float *va,const int*nnz,int*istat);
int BLAS_dusset_elements(blas_sparse_matrix A,const int * ia, const int *ja, const double *  va, const int nnz);
void blas_dusset_elements_(blas_sparse_matrix*A,const int *ia,const int *ja,const double *va,const int*nnz,int*istat);
int BLAS_cusset_elements(blas_sparse_matrix A,const int * ia, const int *ja, const void * va, const int nnz);
void blas_cusset_elements_(blas_sparse_matrix*A,const int *ia,const int *ja,const void *va,const int*nnz,int*istat);
int BLAS_zusset_elements(blas_sparse_matrix A,const int * ia, const int *ja, const void * va, const int nnz);
void blas_zusset_elements_(blas_sparse_matrix*A,const int *ia,const int *ja,const void *va,const int*nnz,int*istat);

int BLAS_susset_element(blas_sparse_matrix A,const int i, const int j, float *  v);
void blas_susset_element_(blas_sparse_matrix*A,const int*i,const int*j,float *v,int*istat);
int BLAS_dusset_element(blas_sparse_matrix A,const int i, const int j, double *  v);
void blas_dusset_element_(blas_sparse_matrix*A,const int*i,const int*j,double *v,int*istat);
int BLAS_cusset_element(blas_sparse_matrix A,const int i, const int j, void * v);
void blas_cusset_element_(blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat);
int BLAS_zusset_element(blas_sparse_matrix A,const int i, const int j, void * v);
void blas_zusset_element_(blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat);

int BLAS_susget_element(blas_sparse_matrix A,const int i, const int j, float *  v);
void blas_susget_element_(blas_sparse_matrix*A,const int*i,const int*j,float *v,int*istat);
int BLAS_dusget_element(blas_sparse_matrix A,const int i, const int j, double *  v);
void blas_dusget_element_(blas_sparse_matrix*A,const int*i,const int*j,double *v,int*istat);
int BLAS_cusget_element(blas_sparse_matrix A,const int i, const int j, void * v);
void blas_cusget_element_(blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat);
int BLAS_zusget_element(blas_sparse_matrix A,const int i, const int j, void * v);
void blas_zusget_element_(blas_sparse_matrix*A,const int*i,const int*j,void *v,int*istat);






#define BLAS_ussp rsb__wp__BLAS_ussp
#define BLAS_usgp rsb__wp__BLAS_usgp
int BLAS_ussp( blas_sparse_matrix A, int pname );
int BLAS_usgp( blas_sparse_matrix A, int pname );



blas_sparse_matrix rsb_blas_file_mtx_load(const rsb_char_t * filename, rsb_type_t typecode ); /* This is a librsb extension. */


struct rsb_mtx_t * rsb_blas_get_mtx(blas_sparse_matrix A);

#ifdef __cplusplus
}
#endif  /* __cplusplus */


#endif /* RSB_LIBSPBLAS_H_INCLUDED */


