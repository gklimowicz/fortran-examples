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
/*! @file
 *  @author Michele Martone
 *  @brief Implementation of the library user interface specified in rsb.h and rsb_types.h.
 */
/*
 *  The user interface functions and data structures for this library implementations.
 *  (functions declared as static are not intended to be part of the user interface)
 *  Internals should not be present in this file.
 *  As a good rule, each function defined in this file (as well as any internal) shall NOT call another similar, but only its internal wrapper.
 * */
#include "rsb_internals.h"
#include <stdio.h>
#include "rsb_do.h"

struct rsb_session_handle_t rsb_global_session_handle;

#define RSB_INTERFACE_RETURN_MTX_ERRP(MTXAP,ERRVAL,ERRVALP) \
	                                 RSB_INTERFACE_ENDCMD \
	RSB_CONDITIONAL_ERRPSET(ERRVALP,ERRVAL) RSB_DO_MTX_RETURN_INTERFACE(MTXAP,ERRVAL);
#define RSB_INTERFACE_RETURN_MTX(MTXAP)  RSB_INTERFACE_ENDCMD return MTXAP;
#define RSB_INTERFACE_RETURN_ERR(ERRVAL) 	RSB_INTERFACE_ENDCMD RSB_DO_ERR_RETURN_INTERFACE(ERRVAL)
/* #define RSB_INTERFACE_RETURN_ERR_SILENT(ERRVAL) RSB_INTERFACE_ENDCMD return (ERRVAL); */
#define RSB_INTERFACE_RETURN_VAL(VAL)    RSB_INTERFACE_ENDCMD {return (VAL);}

#define RSB_INITIALIZE_CHECK_MTX_ERRP(ERRVALP)	if( ! rsb__do_was_initialized() ) { RSB_ERROR(RSB_ERRM_UL); RSB_INTERFACE_RETURN_MTX_ERRP(NULL,RSB_ERR_UNSUPPORTED_OPERATION,ERRVALP) }

#define RSB_NULL_IF_NUL(P) (((P)&&!*(P))?NULL:(P)) /* convert '' to NULL */

/*!
 * \internal
 * This library, currently, can be used by only one master thread.
 * Therefore it uses no handle for library execution instances.
 * */

rsb_err_t rsb_lib_init(struct rsb_initopts * iop)
{
	/*!
	   \ingroup rsb_doc_library rsb_doc_rsb 
 	   \brief
	   This is the library initialization function.
	   \n
	   It must be called only once before using any other library function.
	   \n
	   It is allowed to call it again after \ref rsb_lib_exit().
	   \n
	   To fine-tune the library behaviour, one may specify a number of options via the \c iop parameter.
	   \n
	   Options may be specified also after \ref rsb_lib_init() by calling \ref rsb_lib_reinit().
	   \n
	   One may call #RSB_REINIT_SINGLE_VALUE_GET  with flag  #RSB_IO_WANT_IS_INITIALIZED_MARKER  to verify whether the library has been initialized or not.
	   \n
	   If the \c RSB_NUM_THREADS environment variable is set, \ref rsb_lib_init() uses it and sets the number of active threads, thus overriding what detected by the OpenMP runtime (e.g. \c OMP_NUM_THREADS).
	  
	   \param \rsb_io_str_msg
	   \return \rsberrcodemsg

	   An example snippet declaring an error variable accumulator at program's beginning:
           \snippet examples/snippets.c Declare error codes variable
           and initializing the library soon thereafter:
           \snippet examples/snippets.c Initialize the library
	   \see_lib_init
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_init(iop);
	RSB_INTERFACE_RETURN_ERR(errval)
}

/* @cond INNERDOC  */
/* TODO: this is a "in development" function, not yet declared in rsb.h ; shall make it official when complete */
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__lib_get_info_str(int what, rsb_char_t* sbuf, size_t buflen)
{
	/* \see_lib_init */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_lib_get_info_str(what,sbuf,buflen);
	RSB_INTERFACE_RETURN_ERR(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
/* @endcond */

rsb_err_t rsb_lib_set_opt(enum rsb_opt_t iof, const void*iop)
{
	/*!
	   \ingroup rsb_doc_library rsb_doc_rsb

	 Sets value of a library option.
	 \n
 	 A value specified by the request flag \c iof  will be fetched from \c *iop and will be used to update the selected option in the library internal state.

	 \rsb_iof_param_msg
	 \rsb_iop_out_param_msg

	 Example snip:
	 \snippet examples/snippets.c Setting a single optional library parameter

	 \see \rsb_iof_macros
	 \see_lib_init
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	RSB_DO_REINIT_SINGLE_VALUE_C_IOP(iof,iop,RSB_IO_SPECIFIER_SET,errval);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_lib_get_opt(enum rsb_opt_t iof, void*iop)
{
	/*!
	   \ingroup rsb_doc_library rsb_doc_rsb

	 Gets value of a library option.
	 \n
 	 A value specified by the request flag \c iof  will be fetched from the library internal state and \c *iop will be updated accordingly.

	 \rsb_iof_param_msg
	 \rsb_iop_out_param_msg
	 \see \rsb_iof_macros
	 \see_lib_init
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	RSB_DO_REINIT_SINGLE_VALUE_GET(iof,iop,errval);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_lib_set_opt_str(const rsb_char_t* opnp, const rsb_char_t* opvp)
{
	/*!
	   \ingroup rsb_doc_library rsb_doc_rsb

	   Specifies individual library options in order to fine-tune the library behaviour.
	   \n
	   Both the option name and the value shall be expressed as strings, identical to their preprocessor identifiers (see #rsb_opt_t ).
	   The \c opnp string will be translated internally to the corresponding request flag values, and the passed value will be parsed out of the \c opvp string.
	   \n
	  
	   \param \rsb_io_str_msg_opnp
	   \param \rsb_io_str_msg_opvp
	   \return \rsberrcodemsg

           \snippet examples/snippets.c rsb_lib_reinit__rsb_lib_set_opt_str_snip

	   \see_lib_init
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_set_initopt_as_string(opnp,opvp);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_lib_reinit(struct rsb_initopts * iop)
{
	/*!
	   \ingroup rsb_doc_library rsb_doc_rsb
	  
	   Changes the library operation options which were set at initialization time either by a user or as defaults.
	   \n
	   Not all options may be supported, depending on build time library settings. 
	   \n
	   If an unsupported option was specified, an appropriate error (e.g.: #RSB_ERR_UNSUPPORTED_OPERATION) will be returned.  
	   \n
	   On the first error, option processing is interrupted and the remaining options (if any) are not processed.
	   \n
	   Program execution may continue safely even if an error code is returned (that is, library status should be consistent).
	   \n
	   
	   \param \rsb_io_str_msg
	   \return \rsberrcodemsg

           \snippet examples/snippets.c rsb_lib_reinit__rsb_lib_set_opt_str_snip

	   \see_lib_init
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_reinit(iop);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_lib_exit(struct rsb_initopts * iop)
{
	/*!
	   \ingroup rsb_doc_library rsb_doc_rsb
	  
	   Finalize \librsb.
	   \n
	   #rsb_lib_exit should be called after having freed all matrices.
	   \n
	   If not all of the data structures were properly deallocated before, this function may still attempt finalizing the library and return the #RSB_ERR_MEMORY_LEAK error code (this depends on the \c --enable-allocator-wrapper configure time option).
	   Any allocated memory will be lost (\librsb does not keep track of allocated matrices).
	   \n
	   Internal library state will be cleared.
	   After this call, it is legal to initialize the library again, by calling \ref rsb_lib_init().
	   \n
	   On an error, the library state may be inconsistent, so it is advisable to either
	   terminate program execution (rather than forcing a new initialization with \ref rsb_lib_init()).
	   \n
	   Parameter  \c iop  is reserved for future use; for now it is safe to pass #RSB_NULL_EXIT_OPTIONS.
	   \n
	   It should be safe to call \ref rsb_lib_exit() more than once.
	   \n

	   \param \rsb_io_str_msg
	   \return \rsberrcodemsg

	   An example snippet declaring an error variable accumulator at program's beginning:
           \snippet examples/snippets.c Declare error codes variable
           and finalizing the library at program's end:
           \snippet examples/snippets.c Finalize the library
	   \see_lib_init
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_exit();
	RSB_INTERFACE_RETURN_ERR(errval)
}

struct rsb_mtx_t * rsb_mtx_alloc_from_coo_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp)
{
	/*!
 	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   Given as input COO arrays \c VA,IA,JA, allocates and assembles an RSB matrix using separate arrays.
	  
	   \param \rsb_ro_va_ia_ja_desc_msg
	   \param \rsb_nnzA_inp_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg If any of \rsb_nrA or \rsb_ncA is zero, it will be detected on the basis of the \c IA and \c JA arrays and \rsb_flagsA.
	   \param \rsb_nrbows_A_sparse_inp_param_msg
	   \param \rsb_flagsa_coc_param_msg
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage

	   Example snip:
	   \snippet examples/snippets.c Allocate matrix without error flags check
	   And another, with duplicate sum flags:
	   \snippet examples/snippets.c Allocate matrix with error flags check
	   And yet another, allocating a triangular matrix:
	   \snippet examples/snippets.c Allocate a matrix with triangular flags

           \see_lib_alloc
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,&errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp,errval,errvalp);
}

struct rsb_mtx_t * rsb_mtx_alloc_from_coo_inplace(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp)
{
	/*!
	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb
	   
	   \rsb_mtx_alloc_coo_inplace_msg
	   \n
	   \rsb_note_assume_nnz_sized

	   \param \rsb_rw_va_ia_ja_desc_msg
	   \param \rsb_nnzA_inp_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg
	   \param \rsb_nrbows_A_sparse_inp_param_msg
	   \param \rsb_flagsa_coi_param_msg
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage

	   \see_lib_alloc
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxAp = rsb__do_mtx_alloc_from_coo_inplace(VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,&errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp,errval,errvalp);
}

struct rsb_mtx_t * rsb_mtx_free(struct rsb_mtx_t * mtxAp)
{
	/*!
	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   Frees a previously allocated sparse matrix structure.
	   \n
	   In the case the matrix has the #RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS flag, the main three data arrays \rsb_va_ia_ja_decl will not be freed by #rsb_mtx_free (see \rsb_lib_alloc_in_place).

	   \param \rsb_mtxt_inp_param_msg_a
	   \return \rsb_ret_null

           Example freeing a sparse matrix:
           \snippet examples/snippets.c Free a sparse matrix
	   \see_lib_alloc
	 */
	struct rsb_mtx_t * mtxBp = NULL;
	RSB_INTERFACE_PREAMBLE
	mtxBp = rsb__do_mtx_free(mtxAp);
	RSB_INTERFACE_RETURN_MTX(mtxBp);
}

rsb_err_t rsb_mtx_clone(struct rsb_mtx_t ** mtxBpp, rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_flags_t flags)
{
	/*!
	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   This function clones a given matrix, allocating a fresh data structure or overwriting an existing one.
	   \n
	   Target type (specified by \c typecode) can be different from that in the matrix.
	   \c
	   If \c alphap=NULL, the cloned matrix will not be scaled.
	   \n
	   This new structure will be completely separated and independent from the original one.
	   \n
	   Examples:
	 */
	 /**
\code{.c}
// will clone the matrix exactly
errval = rsb_mtx_clone(&mtxBp,RSB_NUMERICAL_TYPE_SAME_TYPE,RSB_TRANSPOSITION_N,NULL,mtxAp,RSB_FLAG_IDENTICAL_FLAGS);
// will clone the transpose of the matrix
errval = rsb_mtx_clone(&mtxBp,RSB_NUMERICAL_TYPE_SAME_TYPE,RSB_TRANSPOSITION_T,NULL,mtxAp,RSB_FLAG_IDENTICAL_FLAGS);
// will clone the lower triangle of the matrix
errval = rsb_mtx_clone(&mtxBp,RSB_NUMERICAL_TYPE_SAME_TYPE,RSB_TRANSPOSITION_N,NULL,mtxAp,RSB_FLAG_TRIANGULAR|RSB_FLAG_LOWER);
\endcode
	 */
	 /**
	   \param \rsb_mtxtpp_inp_param_msg_b If \c *mtxBpp==NULL, a fresh clone will be assigned there; if not, the existing matrix structure will be freed and allocated to host the new one. The case \c *mtxBpp==mtxAp is supported.
	   \param \rsb_type_o_param_msg
	   \param \rsb_transa_inp_param_msg
	   \param \rsb_alpha_s_inp_param_msg Of the type code of \c mtxAp.
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_flags_stru_fla_msg
	   \return \rsberrcodemsg

	   Example snip:
	   \snippet examples/snippets.c Clone and transpose a sparse matrix

	   \see_lib_alloc
	 */
	/* FIXME: what if RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS ? */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__mtx_clone(mtxBpp,typecode,transA,alphap,mtxAp,flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

#if 0
rsb_err_t rsb_get_rows_dense(const struct rsb_mtx_t * mtxAp, void* row, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags )
{
        /*!
	 * \ingroup rsb_doc_matrix_operations rsb_doc_rsb
	 * \return \rsberrcodemsg
         * FIXME : missing test case, document
         * */
        return rsb__do_get_rows_dense(mtxAp,row,frA,lrA,IA,JA,rnzp,flags);
}
#endif



#define RSB_EXPOSE_NEW_GENERAL_INTERFACE 1	/* temporary (internals) to delimit the new interface which supersedes the deprecated one */
#if RSB_EXPOSE_NEW_GENERAL_INTERFACE
#if 0
/* #define RSB_EXTF_NONE		0x00000000*/			/* */
#define RSB_EXTF_SLOWTRI	0x00000001			/*!< Flag values for extracting the strictly lower submatrix*/
#define RSB_EXTF_SUPPTRI	0x00000002			/*!< Flag values for extracting the strictly upper submatrix .*/
#define RSB_EXTF_DIAG		0x00000004			/*!< Flag values for extracting the diagonal submatrix.*/
#define RSB_EXTF_LOWTRI		(RSB_EXTF_SLOWTRI|RSB_EXTF_DIAG)/*!< Flag values for extracting the lower submatrix.*/
#define RSB_EXTF_UPPTRI		(RSB_EXTF_SUPPTRI|RSB_EXTF_DIAG)/*!< Flag values for extracting the upper submatrix.*/
#define RSB_EXTF_OFFDIAG	(RSB_EXTF_SUPPTRI|RSB_EXTF_SLOWTRI)/*!< Flag values for extracting the whole matrix.*/
#define RSB_EXTF_EXPSYMM	0x00000008			/*!< */
#define RSB_EXTF_EXPDIAG	0x00000010			/*!< */
#define RSB_EXTF_EXP		(RSB_EXTF_EXPDIAG|RSB_EXTF_EXPSYMM)	/*!< */
#define RSB_EXTF_ALL		(RSB_EXTF_OFFDIAG|RSB_EXTF_DIAG)/*!< Flag values for extracting the whole matrix.*/
#define RSB_EXTF_EXPALL		(RSB_EXTF_ALL|RSB_EXTF_EXP)	/*!< */
#define RSB_EXTF_DEFAULT	RSB_EXTF_ALL			/*!< Flag values for extracting the whole matrix. */
/* #define RSB_EXTF_SYMMEXP	0x00000020*/
rsb_err_t rsb_get_submatrix_as_coo(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t *mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags);/* NEW, unfinished */

rsb_err_t rsb_get_submatrix_as_coo(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t *mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags/* , rsb_extff_t eflags*/)/* NEW, unfinished */
{
	/*!
	   \ingroup rsb_doc_matrix_conversion rsb_doc_rsb

	   Extracts a submatrix.
	   Call this function with VA,IA,JA NULL in order to get nonzeroes count.

	   \param \rsb_type_param_msg
	   \param \rsb_transa_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_wr_va_ia_ja_desc_msg
	   \param \rsb_inp_rnz_msg
	   \param \rsb_flags_idc_param_msg
	   \return \rsberrcodemsg

	   \warning \rsb_warn_unfinished_msg 
	   \warning \rsb_warn_unfinished_flags_doc_msg
	 */
	 /*
	   \todo: Shall document eflags.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb_do_get_submatrix_as_coo(typecode, transA, alphap, mtxAp, VA, IA, JA, rnzp, flags/* , eflags*/);
	RSB_INTERFACE_RETURN_ERR(errval)
}
#endif
#endif /* RSB_EXPOSE_NEW_GENERAL_INTERFACE */

#if 0
rsb_err_t rsb_spmv_nt(const void *alphap, const struct rsb_mtx_t * mtxAp, const void * x1p, const void * x2p, rsb_coo_idx_t incX, const void * betap, void * y1p, void * y2p, rsb_coo_idx_t incY);
rsb_err_t rsb_spmv_ata(const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY);
rsb_err_t rsb_spmv_power(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp,  rsb_int_t exp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY);

rsb_err_t rsb_spmv_nt(const void *alphap, const struct rsb_mtx_t * mtxAp, const void * x1p, const void * x2p, rsb_coo_idx_t incX, const void * betap, void * y1p, void * y2p, rsb_coo_idx_t incY)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Computes \f$Y_1 \leftarrow \beta Y_1 + \alpha {A}     \cdot X_1 \f$
	   and      \f$Y_2 \leftarrow \beta Y_2 + \alpha {A}^{T} \cdot X_2 \f$.

	   \param \rsb_beta_inp_param_msg
	   \param \rsb_transa_inp_param_msg
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_incx_inp_param_msg
	   \param \rsb_incy_inp_param_msg
	   \param \rsb_y1y2_inp_param_msg
	   \param \rsb_x1x2_inp_param_msg
	   \return \rsberrcodemsg

	   \warning \rsb_warn_untested_msg
	 */

	// FIXME: this is only a placeholder, waiting for a combined implementation.
	// once done, should speedup methods like Biconjugate Gradient (BiCG).
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb_spmv(RSB_TRANSPOSITION_N,alphap,mtxAp,x1p,incX,betap,y1p,incY)|
		rsb_spmv(RSB_TRANSPOSITION_T,alphap,mtxAp,x2p,incX,betap,y2p,incY);
	RSB_INTERFACE_RETURN_ERR(errval)
}
#endif

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
/* Postponed... (from rsb.h) */
rsb_err_t rsb_spata(const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
/* Postponed... */
rsb_err_t rsb_spata(const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Computes \f$Y \leftarrow \beta Y + \alpha {A}^{T} {A} \cdot X \f$.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_x_inp_param_msg
	   \param \rsb_y_out_param_msg
	   \param \rsb_incx_inp_param_msg
	   \param \rsb_incy_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \return \rsberrcodemsg

	   \warning \rsb_warn_unimplemented_msg
	   \warning \rsb_warn_untested_msg
	 */
	/* FIXME: untested; details to be finished ! */
	rsb_err_t errval = RSB_ERR_UNIMPLEMENTED_YET;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_spata(alphap, mtxAp, Xp, incX, betap, Yp, incY);
	RSB_INTERFACE_RETURN_ERR(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if 0
rsb_err_t rsb_spmv_power(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp,  rsb_int_t exp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Y, rsb_coo_idx_t incY)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Computes \f$Y \leftarrow \beta Y + \alpha ({A}^{T})^{exp} {A} \cdot X \f$.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_x_inp_param_msg
	   \param \rsb_y_out_param_msg
	   \param \rsb_incx_inp_param_msg
	   \param \rsb_incy_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \param \rsb_exp_inp_param_msg
	   \return \rsberrcodemsg

	   \warning \rsb_warn_unimplemented_msg
	   \warning \rsb_warn_untested_msg
	 */

	rsb_err_t errval = RSB_ERR_UNIMPLEMENTED_YET;
	RSB_INTERFACE_PREAMBLE
	RSB_INTERFACE_RETURN_ERR(errval)
	// FIXME: this is only a placeholder, waiting for a combined implementation.
}
#endif

rsb_err_t rsb_spmv(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Multiplies a sparse matrix \f$opa(A)\f$ by a vector \f$X\f$, updating vector \f$Y\f$.
	   \n
	   Computes \f$Y \leftarrow \beta Y + \alpha \cdot opa(A) \cdot X \f$.
	   \n
	   It is not allowed to supply same \c Xp and \c Yp  (that is, \c Xp==Yp).
	   \n

	   \rsb_transa_mtx_msg
	   \rsb_num_threads

	   \param \rsb_transa_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_x_inp_param_msg
	   \param \rsb_incx_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \param \rsb_y_out_param_msg
	   \param \rsb_incy_inp_param_msg
	   \return \rsberrcodemsg

	   Example snip:
	   \snippet examples/snippets.c Multiply a sparse matrix by a dense vector

	   \rsb_librsbpp_env
	   \see_lib_spmx
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
       	errval = rsb_do_spmv(transA, alphap, mtxAp, Xp, incX, betap, Yp, incY);
	RSB_INTERFACE_RETURN_ERR(errval)
}

#if 0
rsb_err_t rsb_spmv_sa(const struct rsb_mtx_t * mtxAp, const void * Xp, void * Yp, const void *alphap, rsb_trans_t transA)
{
	/*!
	 * \ingroup rsb_doc_matrix_operations rsb_doc_rsb
	 * computes \f$Y \leftarrow Y + \alpha op(A) \cdot X \f$
	 * \return \rsberrcodemsg
	 * 
	 * */
	if(!alphap || !mtxAp)
		return RSB_ERR_BADARGS;
	return rsb__do_spmv_general(transA,alphap,mtxTp,Xp,1,NULL,Yp,1,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
}
#endif

rsb_err_t rsb_spsv(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, const void * Xp, rsb_coo_idx_t incX, void * Yp, rsb_coo_idx_t incY)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Computes \f$Y \leftarrow \alpha \cdot opt( T )^{-1} \cdot X \f$, with upper or lower triangular \f$T\f$.
	   It is allowed to supply same \c Xp and \c Yp  (that is, \c Xp==Yp).

	   \rsb_transt_mtx_msg

	   \param \rsb_transt_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_inp_param_msg_t
	   \param \rsb_x_inp_param_msg
	   \param \rsb_incx_inp_param_msg
	   \param \rsb_y_out_param_msg
	   \param \rsb_incy_inp_param_msg
	   \return \rsberrcodemsg
	   \rsb_spsv_no_zero

	   Example backsolving a triangular system:
	   \snippet examples/snippets.c Backsolve a triangular system

	   \see_lib_spsx
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_spsv(transT, alphap, mtxTp, Xp, incX, Yp, incY);
	RSB_INTERFACE_RETURN_ERR(errval)
}

#if 0
static rsb_err_t rsb__do_spsv_sxsx(const struct rsb_mtx_t * mtxAp, void * Yp, const void * alphap, rsb_coo_idx_t incX, rsb_trans_t transl)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb
	   computes \f$Y \leftarrow \alpha op(A)^{-1} \cdot Y \f$.
	   \return \rsberrcodemsg
	  
	   It is allowed to use rhs == out, but in this case beta should be set to 1 and incX=incY, or the result will be undefined.
	 */
	return rsb__do_spsv_general(transl,alphap,mtxAp,Yp,1,Yp,1,RSB_OP_FLAG_DEFAULT RSB_INNER_NRHS_SPSV_ARGS_IDS);
}
#endif

#if 0
rsb_err_t rsb_spmv_uxux(const struct rsb_mtx_t * mtxAp, const void * Xp, void * Yp, const void *alphap, const void * betap, rsb_trans_t transA)
{
	/*!
	 * \ingroup rsb_doc_matrix_operations rsb_doc_rsb
	 * computes \f$Y \leftarrow \beta \cdot Y + \alpha\cdot A\cdot X\f$
	 * \return \rsberrcodemsg
	 * */
	if(!alphap || !betap)
		return RSB_ERR_BADARGS;
	return rsb__do_spmv_general(transA,alphap,mtxAp,Xp,1,betap,Yp,1,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
}
#endif

#if 0
rsb_err_t rsb_spmm_az(const struct rsb_mtx_t * mtxAp, const void * mrhs, void *mout, rsb_int_t bstride, rsb_int_t cstride, rsb_int_t nrhs, rsb_trans_t transA)
{
	/*!
	 * \ingroup rsb_doc_matrix_operations rsb_doc_rsb
	 * computes \f$Y \leftarrow op(A) \cdot X \f$
	 * when X is a multi-vector with nrhs elements, mrhs elements having stride bstride and mout elements having stride cstride
	 * \return \rsberrcodemsg
	 * */
	 /* FIXME : and error detection ? **/
#ifdef RSB_HAVE_OPTYPE_SPMM_AZ
	if(!mtxAp || !mout)
		return -1;

	rsb__cblas_Xscal(mtxAp->typecode,nrhs*mtxAp->nr,NULL,mout,1);	/*FIXME:temporary*/

	return rsb_spmm_inner(mtxAp,mrhs,mout,bstride,cstride,nrhs,transA);
#else
	return RSB_ERR_UNSUPPORTED_OPERATION;
#endif
}

rsb_err_t rsb_spmm_sxsx(const struct rsb_mtx_t * mtxAp, const void * Bp, void * Cp, rsb_nnz_idx_t ldB, rsb_nnz_idx_t ldC, rsb_coo_idx_t nrhs, rsb_trans_t transA, const void * alphap, const void * betap, rsb_flags_t order)
{
	/*!
	   \return \rsberrcodemsg
	 */
	return rsb__do_spmm(transA,alphap,mtxAp,nrhs,order,Bp,ldB,betap,Cp,ldC,RSB_OP_FLAG_DEFAULT);
}
#endif

// rsb_err_t rsb_spsm_sxsx(const struct rsb_mtx_t * mtxAp, void * Bp, rsb_nnz_idx_t ldB, rsb_coo_idx_t nrhs, rsb_trans_t transT, const void * alphap, const void * betap, rsb_flags_t order)

rsb_err_t rsb_spsm(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * betap, const void * Bp, rsb_nnz_idx_t ldB, void * Cp, rsb_nnz_idx_t ldC)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Computes \f$Y \leftarrow \alpha \cdot opt( T )^{-1} \cdot B \f$, with upper or lower triangular \f$T\f$.

	   \rsb_transt_mtx_msg

	   \param \rsb_transt_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_inp_param_msg_t
	   \param \rsb_nrhs_inp_param_msg
	   \param \rsb_order_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \param \rsb_b_inp_param_msg
	   \param \rsb_ldb_inp_param_msg
	   \param \rsb_c_inp_param_msg
	   \param \rsb_ldc_inp_param_msg
	   \return \rsberrcodemsg
	   \see_lib_spsx
	 */
	   // \param \rsb_incx_inp_param_msg \param \rsb_incy_inp_param_msg
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_spsm(transT,alphap,mtxTp,nrhs,order,betap,Bp,ldB,Cp,ldC);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_coo_sort(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA,  rsb_type_t typecode, rsb_flags_t flagsA )
{
	/*!
	   \ingroup gr_util rsb_doc_rsb

	   Sorts row-major the given COO input arrays representing a sparse matrix \f$A\f$.

	   \param \rsb_wr_va_ia_ja_desc_msg
	   \param \rsb_nnzA_inp_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg If any of \rsb_nrA or \rsb_ncA is zero, no sort occurs.
	   \param \rsb_type_param_msg
	   \param \rsb_flagsa_inp_param_msg If unsure, use #RSB_FLAG_NOFLAGS.
	   \return \rsberrcodemsg
	   \see_lib_util

	   \note By invoking with swapped \c IA and \c JA (and swapping \c nrA and \c ncA as well) one can obtain column major order.
	 */
	/* \warning \rsb_warn_unfinished_flags_doc_msg */
	/* Does it support Fortran flags ? */
	/* In the future, one may reuse this interface for:
	 * - cleaning up nonzeroes
	 * - sorting in different ways
	 * - compacting
	 * - checking only if input is sorted; e.g.: using the RSB_FLAG_SORTED_INPUT flag
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
#if 0
	/* This is not the default and shall be rechecked. */
	errval = rsb__util_sort_row_major_buffered(VA,IA,JA,nnzA,nrA,ncA,typecode,flags,NULL,0);
#else
	/* This is the default, well tested. */
	errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnzA,nrA,ncA,typecode,flagsA);
#endif
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_coo_cleanup(rsb_coo_idx_t* nnzp, void* VA, rsb_coo_idx_t* IA, rsb_coo_idx_t* JA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_type_t typecode, rsb_flags_t flagsA )
{
	/*!
	   \ingroup gr_util rsb_doc_rsb

	   Compacts the given COO input arrays representing a sparse matrix \f$A\f$.
	   Will either sum together duplicates or use the last one, depending on whether #RSB_FLAG_DUPLICATES_KEEP_LAST or #RSB_FLAG_DUPLICATES_SUM is present in flagsA.

	   It is important that the input is sorted and flagsA shall contain RSB_FLAG_SORTED_INPUT, otherwise the algorithm's complexity will be quadratic.

	   \param nnzp  Pointer to the number of nonzeroes after the cleanup. 
	   \param \rsb_wr_va_ia_ja_desc_msg
	   \param \rsb_nnzA_inp_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_flagsa_inp_param_msg If unsure, use #RSB_FLAG_NOFLAGS.
	   \return \rsberrcodemsg
	   \see_lib_util
	   \see rsb_coo_sort

	   \warning	This is an experimental librsb-1.3 function.

	   \note By invoking with swapped \c IA and \c JA (and swapping \c nrA and \c ncA as well) one can obtain column major order.

	 */
	 /**
	   Examples:

	   \snippet examples/snippets.c COO cleanup 1

	   \snippet examples/snippets.c COO cleanup 2
	 */

	/* TODO: "By invoking" .. in doc/Doxyfile. */
	/* Cleans up and compacts the given COO input arrays representing a sparse matrix \f$A\f$.  Will delete any nonzero non compliant with the given flags. */
	/* Relevant flags are: RSB_FLAG_UNIT_DIAG_IMPLICIT, RSB_FLAG_LOWER_TRIANGULAR, RSB_FLAG_UPPER_TRIANGULAR, RSB_FLAG_DISCARD_ZEROS. */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t nnz = 0;

	RSB_INTERFACE_PREAMBLE
	/* rsb_coo_idx roff=0,coff=0;
	errval = rsb__do_cleanup_nnz(VA,IA,JA,nnzA,roff,coff,nrA,ncA,nnzp,typecode,flagsA); */
	/* TODO: if sorted shall look for duplicates */
	if( nnzp )
		nnz = rsb__weed_out_duplicates(IA,JA,VA,nnzA,typecode,flagsA);
	RSB_SET_IF_NOT_NULL(nnzp,nnz);
	/* TODO: shall accumulate duplicates */
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_file_mtx_get_dims(const char * filename, rsb_coo_idx_t* nrp, rsb_coo_idx_t *ncp, rsb_coo_idx_t *nzp, rsb_flags_t*flagsp)
{
	/*!
	   Reads structural information (dimensions, structural flags) for a matrix file into user specified (and optionally \c NULL) variables.

	   \ingroup rsb_doc_misc rsb_doc_rsb
	   \param \rsb_filename_inp_param_msg
	   \param \rsb_nrcowsp_inp_param_msg
	   \param \rsb_nnzp_inp_param_msg
	   \param \rsb_flagsp_inp_param_msg
	   \return \rsberrcodemsg If read dimensions are illegal (see #rsb_coo_idx_t,#rsb_nnz_idx_t), #RSB_ERR_LIMITS will be returned.

	   Example getting dimensions of a sparse matrix stored in a Matrix Market file:
	   \snippet examples/snippets.c Get dimensions of sparse matrix stored in Matrix Market file

	   \rsb_matrixmarketonlynote_m
	   \note Upper/lower flags will not be reported; hermitiannes do.
	   \see_lib_get
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_file_mtx_get_dims(filename, nrp, ncp, nzp, flagsp);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_perror(void *stream, rsb_err_t errval)
{
	/*!
	   \ingroup rsb_doc_error_handling rsb_doc_rsb
	  
	   Prints out to the specified \c stream  a string corresponding to the error code (using \c <stdio.h>'s \c fprintf).
	   If \c stream==NULL, will print out to the default output stream; see #RSB_IO_WANT_OUTPUT_STREAM .
	   
	   \param stream A \c (FILE*) pointer, as declared in \c <stdio.h>; can be \c NULL.
	   \param \rsb_errval_inp_param_msg
	   \return \rsberrcodemsg
	   \see_lib_error
	 */
	// \warning \rsb_warn_soon_to_be_updated_msg.
	//   \todo : Should use all bits of the errval variable.
	//   \todo : Should rename the function or make a new one matching perror().
	//   \todo : Could invoke this function from rsb_strerror_r(*,NULL,*)
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_perror(stream,errval);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_strerror_r(rsb_err_t errval, rsb_char_t * buf, size_t buflen)
{
	/*!
	   \ingroup rsb_doc_error_handling rsb_doc_rsb

	   Writes a textual description of an error code in the specified string buffer.
	   No more than buflen characters will be written (comprehensive of the terminating \c NUL character).
	   No action will be performed on \c #RSB_ERR_NO_ERROR.
	   Notice too that error flags cannot be added in the way flags are (e.g. \c (#RSB_ERR_GENERIC_ERROR|#RSB_ERR_BADARGS) evaluates to \c #RSB_ERR_GENERIC_ERROR).
	  
	   \param \rsb_errval_inp_param_msg
	   \param \rsb_buf_inp_param_msg
	   \param \rsb_buflen_inp_param_msg

	   \return \rsberrcodemsg

           Examples:
           \snippet examples/snippets.c Copy error message to string
	or
           \snippet examples/snippets.c snip__rsb_strerror_r

	   \see_lib_error
	 */
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_strerror_r(errval,buf,buflen);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_upd_vals(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags, const void * omegap)
{
	/*!
	   \ingroup rsb_doc_matrix_handling rsb_doc_rsb

	   \f$ A \leftarrow op (A,\Omega) \f$
	   Updates the matrix \f$A\f$ by applying either a row-wise or an elemental operation \f$op\f$, which is determined by \c elop_flags.
	   If an unary operation is selected, \c omegap can be \c NULL.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_flags_elop_param_msg
	   \param \rsb_omega_inp_param_msg
	   \return \rsberrcodemsg

           Example snip:
	   \snippet examples/snippets.c snip__rsb_mtx_upd_vals

	   \see_lib_set
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_upd_vals(mtxAp, elop_flags, omegap);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_set_vals(struct rsb_mtx_t * mtxAp, const void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	/*!
	   \ingroup rsb_doc_matrix_handling rsb_doc_rsb

	   Updates the specified matrix elements, if found in the nonzero pattern.

	   In the special case of a matrix in assembly state (that is, one that has been created as empty with #rsb_mtx_alloc_from_coo_begin() and not yet assembled with #rsb_mtx_alloc_from_coo_end() ) all the supplied matrix elements will be accepted: whether already present or not.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_ro_va_ia_ja_desc_msg
	   \param \rsb_nnz_inp_param_msg
	   \param \rsb_flags_setv_inp_param_msg
	   \return \rsberrcodemsg

	   \see_lib_set
	 */

	/* FIXME: new, UNFINISHED */
	/* FIXME: shall document what will do on out-of-pattern elements */
	/* should support sum, max, etc .. */
//	RSB_ERROR("!!\n");
//	if(flags == RSB_FLAG_DUPLICATES_SUM)
//		return RSB_ERR_UNIMPLEMENTED_YET;
//	else
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_set_elements(mtxAp,VA,IA,JA,nnz,flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_vals(const struct rsb_mtx_t * mtxAp, void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	/*!
	   \ingroup rsb_doc_matrix_handling rsb_doc_rsb

	   Gets the specified matrix elements, if found.
	   Please note that unlike #rsb_mtx_set_vals, the matrix has to be fully assembled here.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_wr_va_rd_ia_ja_desc_msg
	   \param \rsb_nnz_inp_param_msg
	   \param \rsb_flags_getv_inp_param_msg
	   \return \rsberrcodemsg

	   Example snip:
           \snippet examples/snippets.c snip__rsb_mtx_get_vals

           \see_lib_get
	 */

	/* may return an ...UNFINALIZED... error here ... */
	/* TODO: could document better error behaviour (e.g.: what if all updated except one ? ) */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_get_elements(mtxAp,VA,IA,JA,nnz,flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_file_mtx_save(const struct rsb_mtx_t * mtxAp, const rsb_char_t * filename)
{
	/*!
	   \ingroup rsb_doc_input_output rsb_doc_rsb

	   Saves the given matrix to the specified matrix file.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_filename_out_param_msg
	   \return \rsberrcodemsg

	   \warning \rsb_warn_flags_not_complete_msg

	   \rsb_matrixmarketonlynote_m

           Example, printing a matrix to standard output:
           \snippet examples/snippets.c Print a matrix to standard output
           \see_lib_info
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_file_mtx_save(mtxAp,RSB_NULL_IF_NUL(filename));
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_file_vec_save(const rsb_char_t * filename, rsb_type_t typecode, const void * Yp, rsb_coo_idx_t yvl)
{
	/*!
	   \ingroup rsb_doc_input_output rsb_doc_rsb

	   Saves a dense vector to the specified file, using the numerical type representation as specified by the user.
	   This function assumes \c Yp!=NULL and \c yvl>0.

	   \param \rsb_filename_inv_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_y_out_param_msg
	   \param \rsb_yvl_param_msg
	   \return \rsberrcodemsg

	   \rsb_matrixmarketonlynote_v

	   Example printing to standard output:
	   \snippet examples/snippets.c Print to stdout an nrA-long numerical vector

           \see_lib_info
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_vec_save(filename, typecode, Yp, yvl);
	RSB_INTERFACE_RETURN_ERR(errval);
}

rsb_err_t rsb_file_vec_load(const rsb_char_t * filename, rsb_type_t typecode, void * Yp, rsb_coo_idx_t *yvlp)
{
	/*!
	   \ingroup rsb_doc_input_output rsb_doc_rsb

	   Loads a dense vector from the specified file, using the numerical type representation as specified by the user.
	   This function is intended to be called in two steps: first with \c Yp=NULL, in order to write the vector length to \c *yvlp ; then, with \c yvlp=NULL, to get \c Yp written.

	   \param \rsb_filename_inv_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_y_inp_param_msg
	   \param \rsb_yvlp_param_msg
	   \return \rsberrcodemsg

	   Example loading vector matrix from file
	   \snippet examples/snippets.c Load vector matrix from file

	   \rsb_matrixmarketonlynote_v
           \see_lib_info
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_load_vector_file_as_matrix_market(filename,typecode,Yp,yvlp);
	RSB_INTERFACE_RETURN_ERR(errval);
}

struct rsb_mtx_t * rsb_file_mtx_load(const rsb_char_t * filename, rsb_flags_t flagsA, rsb_type_t typecode, rsb_err_t *errvalp)
{
	/*!
	   \ingroup rsb_doc_input_output rsb_doc_rsb

	   Loads a sparse matrix from the specified matrix file, assembling it in the format specified by \rsb_flags, using the numerical type representation as specified by the user.
	   \rsb_extra_internal_verbosity

	   \param \rsb_filename_inp_param_msg
	   \param \rsb_flagsa_inp_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage

	   \rsb_matrixmarketonlynote_m

           Example loading a matrix from a Matrix Market file:
           \snippet examples/snippets.c Load a matrix from Matrix Market file
	   \see_lib_info
	 */
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxAp = rsb__dodo_load_matrix_file_as_matrix_market(filename, flagsA, typecode, &errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp,errval,errvalp);
}

struct rsb_mtx_t * rsb_sppsp(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_err_t * errvalp)
{
	/*!
	   \ingroup rsb_doc_matrix_handling rsb_doc_rsb

	   Computes the weighted sum of two sparse matrices, returning a new matrix:
	   \f$C \leftarrow \alpha\cdot transA(A) + \beta\cdot transB{B} \f$
	   Symmetry flags are ignored in this operation.

	   \rsb_transa_mtx_msg
	   \rsb_transb_mtx_msg

	   \param \rsb_type_param_msg
	   \param \rsb_transa_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_abi_param_msg_a
	   \param \rsb_transb_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \param \rsb_mtxt_abi_param_msg_b
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage

	   Example snip:
	   \snippet examples/snippets.c snip__rsb_sppsp

	   \see_lib_gemm

	   \warning \rsb_warn_not_th_tested_msg
	   \warning \rsb_warn_unoptimized_msg 
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxCp = NULL;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxCp = rsb__do_matrix_sum(typecode,transA,alphap,mtxAp,transB,betap,mtxBp,&errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxCp,errval,errvalp);
}

struct rsb_mtx_t * rsb_spmsp(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_err_t * errvalp)
{
	/*!
	   \ingroup rsb_doc_matrix_handling rsb_doc_rsb

	   Computes the weighted product of two sparse matrices in a new sparse matrix (also known as SpGEMM operation):
	   \f$C \leftarrow \alpha \cdot opa(A) \cdot \beta \cdot opb(B) \f$
	   Symmetry/Hermitian flags are ignored by this operation.

	   \rsb_transa_mtx_msg
	   \rsb_transb_mtx_msg

	   \param \rsb_type_param_msg
	   \param \rsb_transa_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_abi_param_msg_a
	   \param \rsb_transb_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \param \rsb_mtxt_abi_param_msg_b
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage

	   Example snip:
	   \snippet examples/snippets.c snip__rsb_spmsp

	   \warning Parameters \c alphap,betap,transA,transB  are not yet taken in consideration. The following defaults are valid: \f$\alpha=1.0\f$ and \f$\beta=1.0\f$, and \c transA=transB=#RSB_TRANSPOSITION_N.

	   \see_lib_gemm
	 */
	/* FIXME: NEW, UNFINISHED, UNTESTED, UNSECURED */
	/* \warning \rsb_warn_not_th_tested_msg \warning \rsb_warn_unoptimized_msg  */
	struct rsb_mtx_t * mtxCp = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxCp = rsb__do_matrix_mul(typecode,transA,alphap,mtxAp,transB,betap,mtxBp,&errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxCp,errval,errvalp);
}

rsb_err_t rsb_mtx_add_to_dense(const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_nnz_idx_t ldB, rsb_nnz_idx_t nrB, rsb_nnz_idx_t ncB, rsb_bool_t rowmajorB, void * Bp)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb
	   
	   Dense matrix B is updated by adding scaled sparse matrix \f${A}\f$ to it:
	   \f$B \leftarrow B + \alpha {A} \f$

	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_abi_param_msg_a
	   \param \rsb_ldb_inp_param_msg
	   \param \rsb_nrcows_B_dense_inp_param_msg
	   \param \rsb_rowmajor_B_inp_param_msg
	   \param \rsb_dmtx_abi_param_msg_b
	   \return \rsberrcodemsg

           Example snip:
	   \snippet examples/snippets.c snip__rsb_mtx_add_to_dense

	   \note Please note that it suffices to 'transpose' \c Bp's description parameters to get \f$A\f$ transposed summed in.
	   \note Symmetry is currently not expanded.
	   \note Threaded, for large enough matrices.
	   \see_lib_gemm
	 */
	/* TODO: add transA */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_matrix_add_to_dense(alphap, mtxAp, ldB, nrB, ncB, rowmajorB, Bp);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_trans_t rsb_psblas_trans_to_rsb_trans(const char psbtrans)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb
	
	    Translate a PSBLAS transposition value character to a \librsb one.
	    \n
	    See the PSBLAS library website/documentation for valid input values.

	   \param \rsb_psb_trans_inp_param_msg 
	   \return A valid transposition code; that is #RSB_TRANSPOSITION_N for 'N', #RSB_TRANSPOSITION_T for 'T', RSB_TRANSPOSITION_C for 'C',  (See \ref matrix_transposition_flags_section).

           Example snip:
	   \snippet examples/snippets.c main_rsb_psblas_trans_to_rsb_trans

	   \see_lib_psblas
	 */
	rsb_trans_t rsbtrans = RSB_TRANSPOSITION_INVALID;
	RSB_INTERFACE_PREAMBLE
	rsbtrans = rsb__do_psblas_trans_to_rsb_trans(psbtrans);
	RSB_INTERFACE_RETURN_VAL(rsbtrans)
}

struct rsb_mtx_t * rsb_mtx_alloc_from_csr_const(const void *VA, const rsb_coo_idx_t * RP, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp)
{
	/*!
 	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   Given input read only CSR format arrays, allocates and assembles an RSB matrix (stored in separate arrays).
	  
	   \param \rsb_ro_va_rp_ja_desc_msg
	   \param \rsb_nnzA_inp_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg
	   \param \rsb_nrbows_A_sparse_inp_param_msg
	   \param \rsb_flagsa_csr_param_msg
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage
	   \see_lib_alloc
	 */
	// FIXME: flags and index and alloc mangling, here 
	// FIXME: UNTESTED, AND NNZ<M ?
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxAp = rsb__do_mtx_alloc_from_csr_const(VA,RP,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,&errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp,errval,errvalp);
}

struct rsb_mtx_t * rsb_mtx_alloc_from_csc_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * CP, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp)
{
	/*!
 	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   Given input read only CSC format arrays, allocates and assembles an RSB matrix (stored in separate arrays).
	  
	   \param \rsb_ro_va_ia_cp_desc_msg
	   \param \rsb_nnzA_inp_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg
	   \param \rsb_nrbows_A_sparse_inp_param_msg
	   \param \rsb_flagsa_csc_param_msg
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage
	   \see_lib_alloc

           Example:
           \snippet examples/snippets.c snip__rsb_mtx_alloc_from_csc_const
	 */
	// FIXME: flags and index and alloc mangling, here 
	// FIXME: UNTESTED, AND NNZ<M ?
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
       	mtxAp = rsb__do_mtx_alloc_from_csc_const(VA,IA,CP,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,&errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp,errval,errvalp);
}

struct rsb_mtx_t * rsb_mtx_alloc_from_csr_inplace (void *VA, rsb_nnz_idx_t * RP, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp )
{
	/*!
 	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   \rsb_mtx_alloc_csr_inplace_msg
	   \n
	   \rsb_note_assume_nnz_sized
	  
	   \param \rsb_wr_va_rp_ja_desc_msg
	   \param \rsb_nnzA_inp_param_msg
	   \param \rsb_type_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg
	   \param \rsb_nrbows_A_sparse_inp_param_msg
	   \param \rsb_flagsa_csr_param_msg
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxpmessage
	   \see_lib_alloc
	 */
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxAp = rsb__do_mtx_alloc_from_csr_inplace (VA, RP, JA, nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval );
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp,errval,errvalp);
}

rsb_err_t rsb_mtx_switch_to_csr(struct rsb_mtx_t * mtxAp, void ** VAp, rsb_coo_idx_t ** IAp, rsb_coo_idx_t ** JAp, rsb_flags_t flags)
{
	/*!
 	   \ingroup rsb_doc_matrix_conversion rsb_doc_rsb
	   
	   Switches the matrix to the CSR format, in-place. 

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_wr_va_ia_ja_p_desc_msg
	   \param \rsb_flags_idc_param_msg Flags #RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS are forbidden.
	   \return \rsberrcodemsg

	   \note \rsb_note_switch_in_place
	   \warning \rsb_warn_not_th_tested_msg

           Example:
           \snippet examples/snippets.c snip__rsb_mtx_switch_to_csr

	   \see_lib_conv
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_switch_rsb_mtx_to_csr_sorted(mtxAp, VAp, IAp, JAp, flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_coo(const struct rsb_mtx_t * mtxAp, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_flags_t flags )
{
	rsb_nnz_idx_t nnz = 0;
	/*! 
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Returns the matrix converted in a coordinate storage format.
	   \n
	   Elements will be stored in no particular order.
	   \n
	   If there are structural or fill-in zero elements, these will be skipped.
	   \n
	   Writes as many entries as there are nonzeroes (use #rsb_mtx_get_info(mtxAp,#RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T,&nnz)) to find out how many in order to allocate the arrays correctly.
	  
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_wr_va_ia_ja_desc_msg
	   \param \rsb_flags_getco_inp_param_msg
	   \return \rsberrcodemsg

	   \see_lib_get
	   */
	   /*
	    No more than mtxAp->nnz elements will be written.
	   \todo Allow optional VA,IA,JA, for pattern matrices, or other purposes.
	   */
	 // FIXME: does not support misc flags !
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_get_coo_noalloc(mtxAp,VA,IA,JA,&nnz,flags);
//err:
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_csr(rsb_type_t typecode, const struct rsb_mtx_t *mtxAp, void * VA, rsb_nnz_idx_t * RP, rsb_coo_idx_t * JA, rsb_flags_t flags )
{
	/*!
 	   \ingroup rsb_doc_matrix_conversion rsb_doc_rsb

	   Fills the given arrays with the matrix expressed in the CSR format. 

	   \param \rsb_type_param_msg
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_wo_va_rp_ja_desc_msg
	   \param \rsb_flags_getcs_inp_param_msg
	   \return \rsberrcodemsg

	   \see_lib_get
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_get_csr(typecode,mtxAp,VA,RP,JA,flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_rows_sparse(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags)
{
        /*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Writes to the given COO arrays the specified submatrix.

	   Invoke with \c VA,IA,JA  set to \c NULL  in order to get the nonzeroes count written to \c *rnzp, and know how large the arrays should be.

	   \rsb_IA_can_null_msg (in this case it will be ignored).
	   The written rows are ordered.
	  
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_wr_va_rd_ia_ja_desc_msg
	   \param \rsb_inp_frlr_msg
	   \param \rsb_inp_rnz_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_transa_inp_param_msg
	   \param \rsb_flags_getrs_inp_param_msg
	   \return \rsberrcodemsg
	   
           Example snip:
	   \snippet examples/snippets.c snip__rsb_mtx_get_rows_sparse

	   \see_lib_get
         */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_get_rows_sparse(transA, alphap, mtxAp, VA, IA, JA, frA, lrA, rnzp, flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_coo_block(const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_coo_idx_t fcA, rsb_coo_idx_t lcA, const rsb_coo_idx_t * IREN, const rsb_coo_idx_t * JREN, rsb_nnz_idx_t *rnzp, rsb_flags_t flags )
{
	/*!
	   \ingroup rsb_doc_matrix_conversion rsb_doc_rsb
	   Writes in COO format the specified submatrix.
	   \n
	   Works in two stages: first the user invokes it with \c VA,IA,JA set to \c NULL  to get \c *rnzp.
	   Then the \c VA,IA,JA arrays can be allocated, and the function called again, this time with \c rnzp=NULL but the \c VA,IA,JA arrays pointers non \c NULL (or at least, one of them).

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_wr_va_ia_ja_desc_msg
	   \param \rsb_inp_frlr_msg
	   \param \rsb_inp_fclc_msg
	   \param \rsb_xren_inp_param_msg
	   \param \rsb_inp_rnz_msg
	   \param \rsb_flags_getcb_inp_param_msg
	   \return \rsberrcodemsg

	   Examples:
           \snippet examples/snippets.c Extract one sparse matrix block
*/
/**
	   And other examples:
\code{.c}
// get nnz count first
errval=rsb_mtx_get_coo_block(mtxAp,NULL,NULL,NULL,frA,lrA,fcA,lcA,NULL,NULL,&rnz,flags )
// allocate VA, IA, JA to rnz elements
...
// get the  rnz  values then
errval=rsb_mtx_get_coo_block(mtxAp,  VA,  IA,  JA,frA,lrA,fcA,lcA,NULL,NULL,NULL,flags )
\endcode
*/
/**
	   \warning Expect this function to change soon (e.g.: have scaling parameters, etc.). Contact the author if you intend to use it.
	   \see_lib_get
	 */
	/* \rsb_VA_can_null_msg (in such case, only the pattern information is extracted). */
	/* FIXME: shall test rsb_mtx_get_coo_block with VA=NULL */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_get_block_sparse(mtxAp,VA,IA,JA,frA,lrA,fcA,lcA,IREN,JREN,rnzp,flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_spmm(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Updates a dense matrix with the product of sparse matrix by dense matrix;
	   that is, computes \f$ C \leftarrow \beta\cdot C + \alpha\cdot opa(A) \cdot B \f$.

	   \rsb_transa_mtx_msg
	   \rsb_num_threads
	   \rsb_spmm_compact_nrhs_msg

	   \param \rsb_transa_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_nrhs_inp_param_msg
	   \param \rsb_order_inp_param_msg
	   \param \rsb_b_inp_param_msg
	   \param \rsb_ldb_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \param \rsb_c_inp_param_msg
	   \param \rsb_ldc_inp_param_msg
 	   \return \rsberrcodemsg
	   \rsb_librsbpp_env
	   \see_lib_spmx
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_spmm(transA,alphap,mtxAp,nrhs,order,Bp,ldB,betap,Cp,ldC,RSB_OP_FLAG_DEFAULT);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_spmsp_to_dense(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp , rsb_nnz_idx_t ldC, rsb_nnz_idx_t nrC, rsb_nnz_idx_t ncC, rsb_bool_t rowmajorC, void *Cp)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb

	   Computes the product of sparse matrices and adds it to a dense matrix:
	   \f$C \leftarrow \alpha opa(A) \cdot \beta \cdot opb(B) \f$.

	   \rsb_transa_mtx_msg
	   \rsb_transb_mtx_msg

	   \param \rsb_type_param_msg
	   \param \rsb_transa_inp_param_msg
	   \param \rsb_alpha_inp_param_msg
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_transb_inp_param_msg
	   \param \rsb_beta_inp_param_msg
	   \param \rsb_mtxt_inp_param_msg_b
	   \param \rsb_ldc_inp_param_msg
	   \param \rsb_nrcows_C_dense_inp_param_msg
	   \param \rsb_rowmajor_C_inp_param_msg
	   \param \rsb_dmtx_abi_param_msg_c
 	   \return \rsberrcodemsg

	   \warning Parameters \c alphap,betap,transA,transB  are not yet taken in consideration. The following defaults are valid: \f$\alpha=1.0\f$ and \f$\beta=1.0\f$, and \c transA=transB=#RSB_TRANSPOSITION_N.

           Example snip:
	   \snippet examples/snippets.c snip__rsb_spmsp_to_dense

	   \see_lib_gemm
	 */
	/* \todo \rsb_todo_unfinished_inc_msg */
	/* \warning \rsb_warn_unfinished_msg \warning \rsb_warn_unfinished_noerr_msg */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_spgemm_to_dense(typecode,transA,alphap,mtxAp,transB,betap,mtxBp,ldC,nrC,ncC,!rowmajorC,Cp,NULL,NULL);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_rndr(const char * filename, const struct rsb_mtx_t*mtxAp, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags)
{
	/*!
	   \ingroup rsb_doc_matrix_operations rsb_doc_rsb
	   Renders a matrix to a file.
	   Currently, only Encapsulated Postscript (EPS) is supported.
	
	   \param \rsb_filename_out_param_msg
	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_render_pmwidth_inp_param_msg
	   \param \rsb_render_pmheight_inp_param_msg
	   \param \rsb_render_rflags_inp_param_msg

           Example rendering a sparse matrix to Postscript:
	   \snippet examples/snippets.c Render a Sparse matrix to Postscript

	   Setting environment variable \c RSB_USE_HOSTNAME=0 prevents hostname being in the EPS plot internal comments.

	   \see_lib_rndr
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_mtx_render(RSB_NULL_IF_NUL(filename), mtxAp, pmWidth, pmHeight, rflags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_file_mtx_rndr(void * pmp, const char * filename, rsb_coo_idx_t pmlWidth, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb

	   Renders as pixel map the matrix contained in a matrix file.

	   \param \rsb_render_pmp_inp_param_msg
	   \param \rsb_filename_inp_param_msg
	   \param \rsb_render_pmlwidth_inp_param_msg
	   \param \rsb_render_pmwidth_inp_param_msg
	   \param \rsb_render_pmheight_inp_param_msg
	   \param \rsb_render_rflags_inp_param_msg
	   \return \rsberrcodemsg
	   
	   \note At the time being, \c pmlWidth is required to be equal to \c pmWidth.

	   Example rendering a matrix from a Matrix Market file to a pixelmap in memory:
	   \snippet examples/snippets.c Render matrix from Matrix Market pixelmap in memory
		
	   \see_lib_rndr
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_file_mtx_rndr(pmp, filename, pmlWidth, pmWidth, pmHeight, rflags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_switch_to_coo(struct rsb_mtx_t * mtxAp, void ** VAp, rsb_coo_idx_t ** IAp, rsb_coo_idx_t ** JAp, rsb_flags_t flags)
{
	/*!
 	   \ingroup rsb_doc_matrix_conversion rsb_doc_rsb

	   Switches a matrix to COO arrays in place.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_wr_va_ia_ja_p_desc_msg
	   \param \rsb_flags_swcoo_inp_param_msg
	   \return \rsberrcodemsg

	   \note \rsb_note_switch_in_place
	   \warning \rsb_warn_not_th_tested_msg

           Example:
           \snippet examples/snippets.c snip__rsb_mtx_switch_to_coo

	   \see_lib_conv
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_switch_rsb_mtx_to_coo(mtxAp, VAp, IAp, JAp, flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_prec(void *opdp, const struct rsb_mtx_t * mtxAp, rsb_precf_t prec_flags, const void *ipdp)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb

	   A function computing a simple preconditioner out of \c mtxAp.

	   \param opdp Preconditioner data pointer (output).
	   \param \rsb_mtxt_inp_param_msg_a
	   \param prec_flags Valid preconditioner request flags (currently, only #RSB_PRECF_ILU0 is supported; for it, \c *opdp will be overwritten with two \c rsb_mtx_t pointers, respectively a lower and an upper matrix.).
	   \param ipdp  Preconditioner data pointer (input) (ignored at the moment).

	   \return \rsberrcodemsg

           Example:
           \snippet examples/snippets.c snip__rsb_mtx_get_prec

	   \note Matrix should be square, have at least two rows, and have at least one nonzero.
	   \see_lib_get
	*/
	/*
	   \warning \rsb_warn_not_th_tested_msg
	*/
	/* FIXME: temporary interface */
	rsb_err_t errval = RSB_ERR_UNIMPLEMENTED_YET;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_get_preconditioner(opdp,mtxAp,prec_flags,ipdp);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_info(const struct rsb_mtx_t *mtxAp, enum rsb_mif_t miflags, void* minfop)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb

	   \rsb_mtx_getinfo_msg.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_miflags_inp_param_msg
	   \param \rsb_minfop_inp_param_msg

	   \return \rsberrcodemsg

	   Example snip:
	   \snippet examples/snippets.c snip__rsb_mtx_get_info

	   \see_lib_info
	*/
	rsb_err_t errval = RSB_ERR_UNIMPLEMENTED_YET;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_mtx_get_info(mtxAp, miflags, minfop);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_info_str(const struct rsb_mtx_t *mtxAp, const rsb_char_t *mis, void* minfop, size_t buflen)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb

	   \rsb_mtx_getinfo_msg, via a string form query.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param mis A string specifying any identifier among the matrix info ones. See #rsb_mif_t for a list of valid identifiers that can be supplied in string form.
	   \param \rsb_minfop_inp_param_msg
	   \param buflen If greater than 0, \c minfop will be treated as a string of length \c buflen and filled with the desired value via the standard \c snprintf() function.

	   \return \rsberrcodemsg

	   Example snip:
	   \snippet examples/snippets.c Get an info string for the matrix

	   \see_lib_info
	*/
	/* \warning \rsb_warn_not_th_tested_msg */
	rsb_err_t errval = RSB_ERR_UNIMPLEMENTED_YET;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_get_matrix_info_from_string(mtxAp,mis,minfop,buflen);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_nrm(const struct rsb_mtx_t * mtxAp , void * Np, enum rsb_extff_t flags)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb

	   Computes a matrix norm (either infinite-norm or or 2-norm or 1-norm).

	   \param \rsb_mtxt_inp_param_msg_a
	   \param Np  Points to a scalar value which will be overwritten with the selected norm.
	   \param flags Either #RSB_EXTF_NORM_ONE or #RSB_EXTF_NORM_TWO or #RSB_EXTF_NORM_INF.

	   In case of a complex type, only the real part will be written to \c Np. 

	   \return \rsberrcodemsg
	   \see_lib_get
	*/
	rsb_err_t errval = RSB_ERR_BADARGS;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_matrix_norm(mtxAp, Np, flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_mtx_get_vec(const struct rsb_mtx_t * mtxAp , void * Dp, enum rsb_extff_t flags)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb

	   Will overwrite a supplied array with a specific vector quantity.

	   \param \rsb_mtxt_inp_param_msg_a
	   \param \rsb_d_inp_param_msg
	   \param flags Either one of the different extraction filter flags (e.g.: #RSB_EXTF_DIAG, #RSB_EXTF_SUMS_ROW, ...) .
	   \return \rsberrcodemsg
	   \see_lib_get
	*/
	rsb_err_t errval = RSB_ERR_BADARGS;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_matrix_compute(mtxAp,Dp,flags);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_time_t rsb_time(void)
{
	/*!
	   \ingroup rsb_doc_misc rsb_doc_rsb

	   Returns the current time in seconds.
	   This function is meant to be used for computing wall clock time intervals (e.g.: for benchmarking purposes). 
	   The user should not rely on this function for absolute time computations.

	   \return A value for the current time, in seconds.
	   \see_lib_util
	 */
	return rsb__do_time();
}

#if RSB_WANT_COO_BEGIN 
struct rsb_mtx_t * rsb_mtx_alloc_from_coo_begin(rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_flags_t flagsA, rsb_err_t * errvalp)
{
	/*!
 	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   Creates an empty matrix structure in assembly state.
	   The user then populates it using #rsb_mtx_set_vals() repeatedly; then assembles it with #rsb_mtx_alloc_from_coo_end().
	  
	   \param \rsb_nnzA_inp_param_msg_i
	   \param \rsb_type_param_msg
	   \param \rsb_nrcows_A_sparse_inp_param_msg
	   \param \rsb_flagsa_coc_param_msg
	   \param \rsb_errvp_inp_param_msg
	   \return \rsbmtxapmessage
	   \warning \rsb_warn_not_th_tested_msg
	   \see_lib_alloc
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxAp = rsb__do_mtx_alloc_from_coo_begin(nnzA,typecode,nrA,ncA,flagsA,&errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp,errval,errvalp);
}

rsb_err_t rsb_mtx_alloc_from_coo_end(struct rsb_mtx_t ** mtxApp)
{
	/*!
 	   \ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	   Assembles RSB arrays for a matrix in build state created with #rsb_mtx_alloc_from_coo_begin() and populated with #rsb_mtx_set_vals().
	   \n
	   After assembly, any operation on the matrix is allowed.
	  
	   \param \rsb_mtxt_inp_param_msg_i
	   \return \rsberrcodemsg
	   \warning \rsb_warn_not_th_tested_msg
	   \note Note that the memory location of the matrix will be changed by this call, and the (old) \c *mtxApp  address value will be not valid anymore.
	   \see_lib_alloc
	 */
	rsb_err_t errval = RSB_ERR_BADARGS;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_mtx_alloc_from_coo_end(mtxApp);
	RSB_INTERFACE_RETURN_ERR(errval)
}
#endif

#if 0
rsb_err_t rsb_tune_wrt(struct rsb_mtx_t ** mtxOpp, rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, const struct rsb_mtx_t * mtxAp)
{
	/*!
 	\ingroup rsb_doc_matrix_assembly rsb_doc_rsb

	Tunes matrix with respect to a user specified "benchmark" or "performance oracle" function.
	...
	\rsb_version_12
	*/

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* rsb__tune_spxx_bos (...) */
	return errval;
}
#endif

rsb_err_t rsb_tune_spmm(struct rsb_mtx_t ** mtxOpp, rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC)
{
	/*!
 	\ingroup rsb_doc_matrix_assembly rsb_doc_rsb

       	An auto-tuner: optimizes either the matrix instance, the thread count or both for the #rsb_spmm operation.

	\rsb_tune__doc_msg
	\param \rsb_tune_mtxOpp_iou_param_msg
	\param \rsb_tune_sfp_iou_param_msg
	\param \rsb_tune_tnp_iou_param_msg
	\param \rsb_tune_maxr_iou_param_msg
	\param \rsb_tune_maxt_iou_param_msg
	\param \rsb_transa_inp_param_msg
	\param \rsb_alpha_inp_param_msg
	\param \rsb_mtxt_inp_param_msg_a
	\param \rsb_nrhs_inp_param_msg
	\param \rsb_order_inp_param_msg
	\param \rsb_b_tune_inp_param_msg
	\param \rsb_ldb_inp_param_msg
	\param \rsb_beta_inp_param_msg
	\param \rsb_c_tune_inp_param_msg
	\param \rsb_ldc_inp_param_msg
	\return \rsberrcodemsg

	 */
	 /**
	   Examples:
\code{.c}
// obtain best thread count for mtxAp:
errval = rsb_tune_spmm(NULL  ,&sf,&tn ,maxr,maxt,transA,&alpha,mtxAp,nrhs,order,Bp,ldB,&beta,Cp,ldC);

// obtain best thread count for mtxAp; Bp and Cp will be allocated by the tuner:
errval = rsb_tune_spmm(NULL  ,&sf,&tn ,maxr,maxt,transA,&alpha,mtxAp,nrhs,order,NULL,0,&beta,NULL,0);

// obtain best clone of mtxAp (for current thread count):
assert(mtxOp == NULL && mtxAp != NULL);
errval = rsb_tune_spmm(&mtxOp,&sf,NULL,maxr,maxt,transA,&alpha,mtxAp,nrhs,order,Bp,ldB,&beta,Cp,ldC);

// obtain best clone of mtxAp and best thread count:
assert(mtxOp == NULL && mtxAp != NULL);
errval = rsb_tune_spmm(&mtxOp,&sf,&tn ,maxr,maxt,transA,&alpha,mtxAp,nrhs,order,Bp,ldB,&beta,Cp,ldC);

// replace mtxAp with best clone (if any):
errval = rsb_tune_spmm(&mtxAp,&sf,NULL,maxr,maxt,transA,&alpha,NULL ,nrhs,order,Bp,ldB,&beta,Cp,ldC);

// replace mtxAp with best clone (if any) and obtain best thread count:
errval = rsb_tune_spmm(&mtxAp,&sf,&tn ,maxr,maxt,transA,&alpha,NULL ,nrhs,order,Bp,ldB,&beta,Cp,ldC);

// illegal call:
assert(mtxOp != NULL && mtxAp != NULL);
errval = rsb_tune_spmm(&mtxOp,&sf,&tn ,maxr,maxt,transA,&alpha,mtxAp,nrhs,order,Bp,ldB,&beta,Cp,ldC);
\endcode
	 */
	 /**
	\warning 
	\rsb_tune_warning_doc_msg
	\todo
	\rsb_tune_todo_doc_msg
	\rsb_tune__doc_images
	\see_lib_spmx
	*/
	rsb_err_t errval = RSB_ERR_BADARGS;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_tune_spmm( mtxOpp, sfp, tnp, maxr, maxt, transA, alphap, mtxAp, nrhs, order, Bp, ldB, betap, Cp, ldC);
	RSB_INTERFACE_RETURN_ERR(errval)
}

rsb_err_t rsb_tune_spsm(struct rsb_mtx_t ** mtxOpp, rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC)
{
	/*!
 	\ingroup rsb_doc_matrix_assembly rsb_doc_rsb

       	An auto-tuner: optimizes either the matrix instance, the thread count or both for the #rsb_spsm operation.

	\rsb_tune__doc_msg
	\param \rsb_tune_mtxOpp_iou_param_msg
	\param \rsb_tune_sfp_iou_param_msg
	\param \rsb_tune_tnp_iou_param_msg
	\param \rsb_tune_maxr_iou_param_msg
	\param \rsb_tune_maxt_iou_param_msg
	\param \rsb_transa_inp_param_msg
	\param \rsb_alpha_inp_param_msg
	\param \rsb_mtxt_inp_param_msg_a
	\param \rsb_nrhs_inp_param_msg
	\param \rsb_order_inp_param_msg
	\param \rsb_b_tune_inp_param_msg
	\param \rsb_ldb_inp_param_msg
	\param \rsb_beta_inp_param_msg
	\param \rsb_c_tune_inp_param_msg
	\param \rsb_ldc_inp_param_msg
	\return \rsberrcodemsg

	\rsb_spsv_no_zero
	\warning 
	\rsb_tune_warning_doc_msg
	\todo
	\rsb_tune_todo_doc_msg
	\rsb_tune__doc_images
	\see_lib_spsx
	\see rsb_tune_spmm
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_INTERFACE_PREAMBLE
	errval = rsb__do_tune_spsm( mtxOpp, sfp, tnp, maxr, maxt, transA, alphap, mtxAp, nrhs, order, Bp, ldB, betap, Cp, ldC);
	RSB_INTERFACE_RETURN_ERR(errval)
}

/*
struct rsb_mtx_t * rsb_BLAS_get_mtx(blas_sparse_matrix handle)
{
	struct rsb_mtx_t * mtxAp = NULL;
	RSB_INTERFACE_PREAMBLE
	mtxAp = rsb_do_BLAS_get_mtx(handle);
	RSB_INTERFACE_RETURN_MTX(mtxAp);
}
*/

