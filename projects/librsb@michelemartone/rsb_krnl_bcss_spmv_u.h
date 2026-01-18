/* @cond INNERDOC */
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block compressed sparse stripes format.
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
#ifndef RSB_BCSS_SPMV_U_H_INCLUDED
#define RSB_BCSS_SPMV_U_H_INCLUDED
#include "rsb.h"
#include "rsb_common.h"
#include "rsb_internals.h"



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_C__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_H__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_C__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_H__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_C__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_H__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_C__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_H__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_C__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_H__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tN_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tT_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tC_r1_c1_uu_sU_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tN_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tT_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tC_r1_c1_uu_sS_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tN_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tT_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tC_r1_c1_uu_sH_dE_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tN_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tT_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tC_r1_c1_uu_sU_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tN_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tT_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tC_r1_c1_uu_sS_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tN_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tT_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_C__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_H__tC_r1_c1_uu_sH_dI_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_C__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_H__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_C__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_H__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_C__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_H__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_C__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_H__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_C__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_H__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tN_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tT_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tC_r1_c1_uu_sU_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tN_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tT_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tC_r1_c1_uu_sS_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tN_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tT_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tC_r1_c1_uu_sH_dE_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tN_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tT_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tC_r1_c1_uu_sU_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tN_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tT_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tC_r1_c1_uu_sS_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tN_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tT_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_C__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_H__tC_r1_c1_uu_sH_dI_uG(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uaua_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uauz_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_uxua_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_unua_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sasa_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);



rsb_err_t rsb__BCSR_spmv_sxsa_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy);


#endif /* RSB_BCSS_SPMV_U_H_INCLUDED */
/* @endcond */
