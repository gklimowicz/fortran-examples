/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains MKL interfacing functions.
 * */


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

#ifndef RSB_RSB_MKL_H_INCLUDED
#define RSB_RSB_MKL_H_INCLUDED
#include "rsb_internals.h"
#if RSB_WANT_MKL
#ifdef RSB_INCLUDE_MKL_MKL_H
#include <mkl/mkl.h>
#include <mkl/mkl_blas.h>	/* dgemm, ... */
#include <mkl/mkl_spblas.h>
#else
#include <mkl.h>
#include <mkl_blas.h>	/* dgemm, ... */
#include <mkl_spblas.h>
#endif
/* #include <mkl_types.h> */
/* #include <mkl_service.h> */ /* mkl_get_version */
#if defined(RSB_WANT_MKL) && defined(INTEL_MKL_VERSION) && ( (INTEL_MKL_VERSION) >= 110302 )
#define RSB_WANT_MKL_INSPECTOR 1
#endif /* RSB_WANT_MKL */


#include "rsb_tune.h" /* rsb_tattr_t */



rsb_err_t rsb__mkl_gemv(rsb_type_t typecode, const void * Mp, const void*Bp, void*Xp, rsb_nnz_idx_t mdim, rsb_coo_idx_t vdim, rsb_coo_idx_t*udimp);



rsb_err_t rsb__mkl_coo_spmv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags);
rsb_err_t rsb__mkl_coo_spmm(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nrhs, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags);
rsb_err_t rsb__mkl_coo_spsv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags);
rsb_err_t rsb__do_mkl_csr_spmv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags);

rsb_err_t rsb__mkl_csr_spmv_bench(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp);
rsb_err_t rsb__do_mkl_csr_spmm(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags);
rsb_err_t rsb__mkl_csr_spmm_bench(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp);
rsb_err_t rsb__do_mkl_csr_spsv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags);
rsb_err_t rsb__do_mkl_csr_spsm(const void *VA, const MKL_INT m, const MKL_INT nrhs, const MKL_INT * IP, const MKL_INT *JA, const void * b, void * c, const void *alphap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, const MKL_INT ldb, const MKL_INT ldc);
rsb_err_t rsb__mkl_csr_spsv_bench(const void *VA, const MKL_INT m, const MKL_INT k/*, const MKL_INT n*/, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, /*const MKL_INT ldb,*/ void * c,/* const MKL_INT ldc,*/ const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp);
rsb_err_t rsb__mkl_csr_spsm_bench(const void *VA, const MKL_INT m, const MKL_INT nrhs, const MKL_INT * IP, const MKL_INT *JA, const void * b, void * c, const void *alphap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, const MKL_INT ldb, const MKL_INT ldc, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp);
rsb_err_t rsb__mkl_coo2csr(const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const void *IVA, const MKL_INT * IIA, const MKL_INT *IJA, const void *OVA, const MKL_INT * OIA, const MKL_INT *OJA, rsb_type_t typecode, const MKL_INT mib);
rsb_int_t rsb__print_mkl_version(char *pn, rsb_bool_t los)
;
#if RSB_WANT_MKL_INSPECTOR

rsb_err_t rsb__mkl_inspector_csr_spmv_bench(void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, MKL_INT * IP, MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp, rsb_int_t hintlvl)
;
rsb_err_t rsb__mkl_inspector_csr_spmm_bench(void *VA, const MKL_INT m, const MKL_INT k, const sparse_layout_t layout, const MKL_INT nrhs, const MKL_INT nnz, MKL_INT * IP, MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp, rsb_int_t hintlvl)
;
#endif /* RSB_WANT_MKL_INSPECTOR */

#endif /* RSB_WANT_MKL */
#endif  /* RSB_RSB_MKL_H_INCLUDED */
/* @endcond */

