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
