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
#ifndef RSB_RSB_DO_H_INCLUDED
#define RSB_RSB_DO_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb_common.h"

/**
 * @file
 * @brief
 * Implementation of the interface functions.
 *
 */

rsb_err_t rsb__do_get_rows_sparse(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags);
rsb_err_t rsb__do_scal(struct rsb_mtx_t * mtxAp, const void * d, rsb_trans_t trans);
rsb_err_t rsb__dodo_getdiag( const struct rsb_mtx_t * mtxAp, void * diagonal );
rsb_err_t rsb__do_elemental_binop(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags, const void * alphap);
rsb_err_t rsb__do_elemental_unop(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags);
rsb_nnz_idx_t rsb__dodo_get_rows_nnz(const struct rsb_mtx_t *mtxAp, rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_flags_t flags, rsb_err_t * errvalp);
#define RSB_WANT_PARALLEL_ELEMENTAL_OPS 0 /* FIXME: temporary ! */
#if RSB_WANT_PARALLEL_ELEMENTAL_OPS
rsb_err_t rsb__do_elemental_scale_parallel(struct rsb_mtx_t * mtxAp, const void * alphap);
#endif /* RSB_WANT_PARALLEL_ELEMENTAL_OPS */
rsb_err_t rsb__do_matrix_add_to_dense(const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_nnz_idx_t ldb, rsb_nnz_idx_t nr, rsb_nnz_idx_t nc, rsb_bool_t rowmajor, void * Bp);
rsb_err_t rsb__do_switch_rsb_mtx_to_csr_sorted(struct rsb_mtx_t * mtxAp, void ** VAP, rsb_coo_idx_t ** IAP, rsb_coo_idx_t ** JAP, rsb_flags_t flags);
rsb_err_t rsb__do_get_preconditioner(void *opd, const struct rsb_mtx_t * mtxAp, rsb_precf_t prec_flags, const void *ipd);/* FIXME: temporary interface */
rsb_err_t rsb__do_get_csr(rsb_type_t typecode, const struct rsb_mtx_t *mtxAp, rsb_byte_t * VA, rsb_nnz_idx_t * RP, rsb_coo_idx_t * JA, rsb_flags_t flags);
rsb_err_t rsb__do_get_matrix_info(const struct rsb_mtx_t *mtxAp, enum rsb_mif_t miflags, void* info, size_t buflen);
rsb_err_t rsb__do_check_leak(void);
rsb_err_t rsb__do_matrix_norm(const struct rsb_mtx_t * mtxAp , void * np, enum rsb_extff_t flags);
rsb_err_t rsb__do_load_vector_file_as_matrix_market(const rsb_char_t * filename, rsb_type_t typecode, void * yp, rsb_coo_idx_t *yvlp);
struct rsb_mtx_t * rsb__dodo_load_matrix_file_as_matrix_market(const rsb_char_t * filename, rsb_flags_t flags, rsb_type_t typecode, rsb_err_t *errvalp);
rsb_bool_t rsb__do_was_initialized(void);
rsb_err_t rsb__do_matrix_compute(const struct rsb_mtx_t * mtxAp , void * dp, enum rsb_extff_t flags);
rsb_err_t rsb__do_switch_rsb_mtx_to_coo(struct rsb_mtx_t * mtxAp, void ** VAP, rsb_coo_idx_t ** IAP, rsb_coo_idx_t ** JAP, rsb_flags_t flags);
struct rsb_mtx_t * rsb__do_mtx_alloc_from_coo_begin(rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_flags_t flags, rsb_err_t * errvalp);
rsb_err_t rsb__do_mtx_alloc_from_coo_end(struct rsb_mtx_t ** mtxAp);

/* TODO: this is a "secret" function, not declared in rsb.h ; shall make it official some day */
rsb_err_t rsb__lib_get_info_str(int what, rsb_char_t* sbuf, size_t buflen);
rsb_err_t rsb__do_upd_vals(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags, const void * omegap);
rsb_err_t rsb__do_mtx_get_info(const struct rsb_mtx_t *mtxAp, enum rsb_mif_t miflags, void* minfop);
rsb_err_t rsb__do_file_mtx_save(const struct rsb_mtx_t * mtxAp, const rsb_char_t * filename);
struct rsb_mtx_t * rsb__do_mtx_alloc_from_csr_inplace (void *VA, rsb_coo_idx_t * RP, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp );
rsb_err_t rsb__do_file_mtx_rndr(void * pmp, const char * filename, rsb_coo_idx_t pmlWidth, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags);
rsb_err_t rsb__do_vec_save(const rsb_char_t * filename, rsb_type_t typecode, const void * Yp, rsb_coo_idx_t yvl);

#define RSB_ERR_DEFAULT_INTERFACE_ERROR RSB_ERR_GENERIC_ERROR
#if RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE
/* please note that the code is likely to fail self-consistency tests, if writing to stderr */
#define RSB_DEBUG_VERBOSE_INTERFACE_NOTICE	{ if(rsb_global_session_handle.rsb_g_verbose_interface)RSB_STDERR("In file %20s (in %s) at line %10d:\n",__FILE__,__func__,__LINE__);}
#else
#define RSB_DEBUG_VERBOSE_INTERFACE_NOTICE	{}
#endif /* RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE */

#if RSB_WANT_LIBRSB_TIMER
#define RSB_INTERFACE_TIMER_DCLS rsb_time_t etime, tetime = rsb_global_session_handle.etime;
#define RSB_INTERFACE_TIMER_CMDS { etime = -rsb__do_time(); }
#define RSB_INTERFACE_TIMER_ENDC { etime += rsb__do_time(); rsb_global_session_handle.etime = etime + tetime; }
#else
#define RSB_INTERFACE_TIMER_DCLS
#define RSB_INTERFACE_TIMER_CMDS
#define RSB_INTERFACE_TIMER_ENDC
#endif

#define RSB_INTERFACE_PREAMBLE_DCLS RSB_INTERFACE_TIMER_DCLS
#define RSB_INTERFACE_PREAMBLE_CMDS RSB_INTERFACE_TIMER_CMDS RSB_DEBUG_VERBOSE_INTERFACE_NOTICE
#define RSB_INTERFACE_PREAMBLE RSB_INTERFACE_PREAMBLE_DCLS RSB_INTERFACE_PREAMBLE_CMDS
#define RSB_INTERFACE_ENDCMD RSB_INTERFACE_TIMER_ENDC

#define RSB_REPORTABLE_ERROR(ERRVAL) (RSB_SOME_ERROR(ERRVAL) && (ERRVAL)!=RSB_ERR_ELEMENT_NOT_FOUND) /* Real error, i.e. won't treat 'not found' as an error.. */
#if RSB_OUT_ERR_VERBOSITY
# if RSB_OUT_ERR_VERBOSITY==1
# define RSB_DO_ERR_MANIFEST_INTERFACE(ERRVAL) if(rsb_global_session_handle.error_stream!=NULL)if(RSB_REPORTABLE_ERROR(ERRVAL)){RSB_ERROR(RSB_ERRM_NL);rsb__do_perror(NULL,ERRVAL);/* don't put rsb_perror here or infinite recursion will arise :-) */}
# endif /* RSB_OUT_ERR_VERBOSITY */
# if RSB_OUT_ERR_VERBOSITY==2
# define RSB_DO_ERR_MANIFEST_INTERFACE(ERRVAL) if(RSB_REPORTABLE_ERROR(ERRVAL)){RSB_ERROR(RSB_ERRM_NL);rsb__do_perror(NULL,ERRVAL);rsb__print_trace();}
# endif /* RSB_OUT_ERR_VERBOSITY */
# if RSB_OUT_ERR_VERBOSITY>=3 && RSB_OUT_ERR_VERBOSITY<=98
/* it would be better to put the following error in the configure script */
# error Error verbosity (set at configure time with --enable-interface-error-verbosity) shall be either 0,1,2,99 !
# endif /* RSB_OUT_ERR_VERBOSITY */
# if RSB_OUT_ERR_VERBOSITY==99
# define RSB_DO_ERR_MANIFEST_INTERFACE(ERRVAL) {if(RSB_REPORTABLE_ERROR(ERRVAL)){RSB_ERROR(RSB_ERRM_NL);rsb__do_perror(NULL,ERRVAL);RSB_STDOUT("Terminating program now.\n");rsb__print_trace();RSB_EXIT(RSB_ERR_TO_PROGRAM_ERROR(ERRVAL));}}
# endif /* RSB_OUT_ERR_VERBOSITY */
#else /* RSB_OUT_ERR_VERBOSITY */
# define RSB_DO_ERR_MANIFEST_INTERFACE(ERRVAL) 
#endif /* RSB_OUT_ERR_VERBOSITY */
#define RSB_DO_ERR_RETURN_INTERFACE(ERRVAL) {RSB_DO_ERR_MANIFEST_INTERFACE(ERRVAL) return (ERRVAL);}
#define RSB_DO_MTX_RETURN_INTERFACE(MATRIX,ERRVAL) {RSB_DO_ERR_MANIFEST_INTERFACE(ERRVAL) return (MATRIX);}

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_RSB_DO_H_INCLUDED */
/* @endcond */
