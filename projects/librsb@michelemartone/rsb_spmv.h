/* @cond INNERDOC */
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block coordinates format.
 Kernels unrolled, with no loops, for only user-specified blockings.
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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains functions for sparse recursive multicore matrix vector multiplication.
 */
/*
 * FIXME: many beta-related ops are NOT parallel, and this is BAD.
 *
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



#ifndef RSB_SPMV_H_INCLUDED
#define RSB_SPMV_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "rsb_internals.h"
#include "rsb_lock.h"
#if RSB_USE_LIBRSBPP
#include "librsbpp.h"
#endif /* RSB_USE_LIBRSBPP */
#if RSB_USE_MKL
#include "rsb_mkl.h"
#endif /* RSB_USE_MKL */
#define RSB_ENABLE_INNER_NRHS_SPMV 1
#if RSB_ENABLE_INNER_NRHS_SPMV
#define RSB_INNER_NRHS_SPMV_ARGS	,const rsb_int_t nrhs, /*const size_t outtot, const size_t rhstot,*/ const size_t outnri, const size_t rhsnri
#define RSB_INNER_NRHS_SPMV_ARGS_IDS	,nrhs/*,outtot,rhstot*/,outnri,rhsnri
#define RSB_INNER_NRHS_SPMV_YSCA_IDS	,nrhs,outnri
#define RSB_OUTER_NRHS_SPMV_ARGS	,const rsb_int_t nrhs, const size_t outnri, const size_t rhsnri
#define RSB_OUTER_NRHS_SPMV_ARGS_IDS	,nrhs,outnri,rhsnri
#else /* RSB_ENABLE_INNER_NRHS_SPMV */
#define RSB_INNER_NRHS_SPMV_ARGS	
#define RSB_INNER_NRHS_SPMV_YSCA_IDS		/* */
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
#define RSB_DEFAULT_INNER_NRHS_SPMV_ARGS	,1,/*0,0,*/0,0
#define RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS	,1,0,0

rsb_err_t rsb__do_spmv_non_recursive(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA RSB_INNER_NRHS_SPMV_ARGS)
;
rsb_err_t rsb__do_spmv_recursive_parallel(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPMV_ARGS)
;
rsb_err_t rsb__do_spmv_recursive_serial(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA RSB_INNER_NRHS_SPMV_ARGS)
;
rsb_err_t rsb__do_spmv_general(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, const void * x, rsb_coo_idx_t incx, const void * betap, void * y, rsb_coo_idx_t incy, enum rsb_op_flags_t op_flags RSB_OUTER_NRHS_SPMV_ARGS)
;

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_SPMV_H_INCLUDED */

/* @endcond */
