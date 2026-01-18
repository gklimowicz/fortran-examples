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
#ifndef RSB_BCSS_MISC_U_H_INCLUDED
#define RSB_BCSS_MISC_U_H_INCLUDED
#include "rsb_internals.h"
#include "rsb.h"


rsb_err_t rsb__BCSR_infty_norm_double_C__tN_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tN_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tN_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tN_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tT_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tT_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tT_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tT_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tC_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tC_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tC_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tC_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tN_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tN_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tN_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tN_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tT_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tT_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tT_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tT_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tC_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tC_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tC_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tC_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tN_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tN_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tN_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tN_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tT_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tT_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tT_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tT_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tC_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C__tC_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tC_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H__tC_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tN_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tN_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tN_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tN_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tT_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tT_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tT_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tT_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tC_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tC_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tC_r1_c1_uu_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tC_r1_c1_uu_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tN_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tN_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tN_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tN_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tT_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tT_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tT_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tT_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tC_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tC_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tC_r1_c1_uu_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tC_r1_c1_uu_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tN_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tN_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tN_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tN_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tT_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tT_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tT_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tT_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tC_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C__tC_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tC_r1_c1_uu_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H__tC_r1_c1_uu_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_double_C__tN_r1_c1_uu_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tN_r1_c1_uu_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tN_r1_c1_uu_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tN_r1_c1_uu_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tT_r1_c1_uu_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tT_r1_c1_uu_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tT_r1_c1_uu_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tT_r1_c1_uu_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tC_r1_c1_uu_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tC_r1_c1_uu_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tC_r1_c1_uu_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tC_r1_c1_uu_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tN_r1_c1_uu_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tN_r1_c1_uu_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tN_r1_c1_uu_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tN_r1_c1_uu_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tT_r1_c1_uu_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tT_r1_c1_uu_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tT_r1_c1_uu_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tT_r1_c1_uu_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tC_r1_c1_uu_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tC_r1_c1_uu_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tC_r1_c1_uu_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tC_r1_c1_uu_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tN_r1_c1_uu_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tN_r1_c1_uu_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tN_r1_c1_uu_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tN_r1_c1_uu_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tT_r1_c1_uu_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tT_r1_c1_uu_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tT_r1_c1_uu_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tT_r1_c1_uu_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tC_r1_c1_uu_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C__tC_r1_c1_uu_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tC_r1_c1_uu_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H__tC_r1_c1_uu_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_infty_norm_float_C__tN_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tN_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tN_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tN_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tT_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tT_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tT_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tT_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tC_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tC_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tC_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tC_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tN_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tN_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tN_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tN_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tT_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tT_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tT_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tT_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tC_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tC_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tC_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tC_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tN_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tN_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tN_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tN_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tT_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tT_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tT_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tT_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tC_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C__tC_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tC_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H__tC_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tN_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tN_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tN_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tN_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tT_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tT_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tT_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tT_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tC_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tC_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tC_r1_c1_uu_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tC_r1_c1_uu_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tN_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tN_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tN_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tN_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tT_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tT_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tT_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tT_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tC_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tC_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tC_r1_c1_uu_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tC_r1_c1_uu_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tN_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tN_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tN_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tN_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tT_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tT_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tT_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tT_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tC_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C__tC_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tC_r1_c1_uu_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H__tC_r1_c1_uu_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_float_C__tN_r1_c1_uu_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tN_r1_c1_uu_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tN_r1_c1_uu_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tN_r1_c1_uu_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tT_r1_c1_uu_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tT_r1_c1_uu_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tT_r1_c1_uu_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tT_r1_c1_uu_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tC_r1_c1_uu_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tC_r1_c1_uu_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tC_r1_c1_uu_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tC_r1_c1_uu_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tN_r1_c1_uu_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tN_r1_c1_uu_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tN_r1_c1_uu_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tN_r1_c1_uu_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tT_r1_c1_uu_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tT_r1_c1_uu_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tT_r1_c1_uu_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tT_r1_c1_uu_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tC_r1_c1_uu_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tC_r1_c1_uu_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tC_r1_c1_uu_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tC_r1_c1_uu_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tN_r1_c1_uu_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tN_r1_c1_uu_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tN_r1_c1_uu_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tN_r1_c1_uu_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tT_r1_c1_uu_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tT_r1_c1_uu_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tT_r1_c1_uu_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tT_r1_c1_uu_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tC_r1_c1_uu_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C__tC_r1_c1_uu_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tC_r1_c1_uu_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H__tC_r1_c1_uu_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tN_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tN_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tN_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tN_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tT_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tT_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tT_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tT_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tC_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tC_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tC_r1_c1_uu_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tC_r1_c1_uu_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tN_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tN_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tN_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tN_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tT_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tT_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tT_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tT_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tC_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tC_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tC_r1_c1_uu_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tC_r1_c1_uu_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tN_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tN_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tN_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tN_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tT_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tT_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tT_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tT_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tC_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C__tC_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tC_r1_c1_uu_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H__tC_r1_c1_uu_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_float_complex_C__tN_r1_c1_uu_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tN_r1_c1_uu_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tN_r1_c1_uu_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tN_r1_c1_uu_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tT_r1_c1_uu_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tT_r1_c1_uu_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tT_r1_c1_uu_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tT_r1_c1_uu_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tC_r1_c1_uu_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tC_r1_c1_uu_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tC_r1_c1_uu_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tC_r1_c1_uu_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tN_r1_c1_uu_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tN_r1_c1_uu_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tN_r1_c1_uu_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tN_r1_c1_uu_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tT_r1_c1_uu_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tT_r1_c1_uu_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tT_r1_c1_uu_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tT_r1_c1_uu_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tC_r1_c1_uu_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tC_r1_c1_uu_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tC_r1_c1_uu_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tC_r1_c1_uu_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tN_r1_c1_uu_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tN_r1_c1_uu_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tN_r1_c1_uu_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tN_r1_c1_uu_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tT_r1_c1_uu_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tT_r1_c1_uu_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tT_r1_c1_uu_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tT_r1_c1_uu_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tC_r1_c1_uu_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C__tC_r1_c1_uu_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tC_r1_c1_uu_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H__tC_r1_c1_uu_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tN_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tN_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tN_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tN_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tT_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tT_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tT_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tT_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tC_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tC_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tC_r1_c1_uu_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tC_r1_c1_uu_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tN_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tN_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tN_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tN_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tT_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tT_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tT_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tT_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tC_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tC_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tC_r1_c1_uu_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tC_r1_c1_uu_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tN_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tN_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tN_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tN_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tT_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tT_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tT_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tT_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tC_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C__tC_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tC_r1_c1_uu_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H__tC_r1_c1_uu_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_double_complex_C__tN_r1_c1_uu_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tN_r1_c1_uu_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tN_r1_c1_uu_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tN_r1_c1_uu_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tT_r1_c1_uu_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tT_r1_c1_uu_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tT_r1_c1_uu_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tT_r1_c1_uu_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tC_r1_c1_uu_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tC_r1_c1_uu_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tC_r1_c1_uu_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tC_r1_c1_uu_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tN_r1_c1_uu_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tN_r1_c1_uu_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tN_r1_c1_uu_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tN_r1_c1_uu_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tT_r1_c1_uu_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tT_r1_c1_uu_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tT_r1_c1_uu_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tT_r1_c1_uu_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tC_r1_c1_uu_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tC_r1_c1_uu_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tC_r1_c1_uu_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tC_r1_c1_uu_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tN_r1_c1_uu_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tN_r1_c1_uu_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tN_r1_c1_uu_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tN_r1_c1_uu_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tT_r1_c1_uu_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tT_r1_c1_uu_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tT_r1_c1_uu_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tT_r1_c1_uu_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tC_r1_c1_uu_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C__tC_r1_c1_uu_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tC_r1_c1_uu_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H__tC_r1_c1_uu_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tN_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tN_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tN_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tN_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tT_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tT_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tT_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tT_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tC_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tC_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tC_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tC_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tN_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tN_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tN_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tN_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tT_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tT_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tT_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tT_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tC_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tC_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tC_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tC_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tN_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tN_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tN_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tN_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tT_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tT_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tT_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tT_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tC_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_C_u_tC_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tC_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_H_u_tC_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tN_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tN_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tN_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tN_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tT_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tT_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tT_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tT_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tC_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tC_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tC_sU_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tC_sU_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tN_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tN_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tN_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tN_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tT_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tT_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tT_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tT_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tC_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tC_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tC_sS_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tC_sS_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tN_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tN_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tN_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tN_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tT_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tT_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tT_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tT_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tC_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_C_u_tC_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tC_sH_dE_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_H_u_tC_sH_dI_uG(const double * VA, double * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_double_C_u_tN_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tN_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tN_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tN_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tT_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tT_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tT_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tT_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tC_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tC_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tC_sU_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tC_sU_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tN_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tN_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tN_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tN_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tT_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tT_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tT_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tT_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tC_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tC_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tC_sS_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tC_sS_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tN_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tN_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tN_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tN_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tT_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tT_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tT_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tT_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tC_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_C_u_tC_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tC_sH_dE_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_scale_double_H_u_tC_sH_dI_uG(double * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double *scale_factors);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tN_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tN_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tN_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tN_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tT_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tT_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tT_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tT_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tC_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tC_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tC_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tC_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tN_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tN_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tN_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tN_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tT_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tT_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tT_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tT_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tC_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tC_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tC_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tC_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tN_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tN_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tN_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tN_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tT_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tT_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tT_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tT_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tC_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_C_u_tC_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tC_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_H_u_tC_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tN_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tN_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tN_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tN_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tT_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tT_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tT_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tT_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tC_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tC_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tC_sU_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tC_sU_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tN_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tN_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tN_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tN_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tT_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tT_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tT_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tT_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tC_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tC_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tC_sS_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tC_sS_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tN_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tN_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tN_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tN_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tT_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tT_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tT_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tT_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tC_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_C_u_tC_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tC_sH_dE_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_H_u_tC_sH_dI_uG(const float * VA, float * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_float_C_u_tN_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tN_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tN_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tN_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tT_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tT_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tT_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tT_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tC_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tC_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tC_sU_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tC_sU_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tN_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tN_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tN_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tN_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tT_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tT_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tT_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tT_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tC_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tC_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tC_sS_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tC_sS_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tN_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tN_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tN_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tN_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tT_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tT_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tT_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tT_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tC_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_C_u_tC_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tC_sH_dE_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_scale_float_H_u_tC_sH_dI_uG(float * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float *scale_factors);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tN_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tN_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tN_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tN_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tT_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tT_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tT_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tT_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tC_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tC_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tC_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tC_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tN_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tN_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tN_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tN_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tT_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tT_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tT_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tT_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tC_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tC_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tC_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tC_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tN_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tN_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tN_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tN_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tT_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tT_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tT_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tT_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tC_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_C_u_tC_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tC_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_float_complex_H_u_tC_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tN_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tN_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tN_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tN_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tT_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tT_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tT_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tT_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tC_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tC_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tC_sU_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tC_sU_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tN_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tN_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tN_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tN_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tT_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tT_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tT_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tT_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tC_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tC_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tC_sS_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tC_sS_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tN_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tN_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tN_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tN_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tT_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tT_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tT_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tT_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tC_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_C_u_tC_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tC_sH_dE_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_float_complex_H_u_tC_sH_dI_uG(const float complex * VA, float complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tN_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tN_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tN_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tN_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tT_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tT_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tT_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tT_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tC_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tC_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tC_sU_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tC_sU_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tN_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tN_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tN_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tN_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tT_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tT_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tT_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tT_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tC_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tC_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tC_sS_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tC_sS_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tN_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tN_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tN_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tN_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tT_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tT_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tT_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tT_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tC_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_C_u_tC_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tC_sH_dE_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_scale_float_complex_H_u_tC_sH_dI_uG(float complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const float complex *scale_factors);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tN_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tN_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tN_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tN_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tT_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tT_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tT_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tT_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tC_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tC_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tC_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tC_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tN_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tN_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tN_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tN_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tT_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tT_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tT_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tT_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tC_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tC_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tC_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tC_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tN_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tN_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tN_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tN_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tT_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tT_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tT_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tT_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tC_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_C_u_tC_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tC_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_infty_norm_double_complex_H_u_tC_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tN_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tN_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tN_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tN_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tT_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tT_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tT_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tT_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tC_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tC_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tC_sU_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tC_sU_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tN_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tN_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tN_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tN_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tT_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tT_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tT_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tT_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tC_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tC_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tC_sS_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tC_sS_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tN_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tN_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tN_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tN_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tT_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tT_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tT_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tT_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tC_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_C_u_tC_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tC_sH_dE_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_rowssums_double_complex_H_u_tC_sH_dI_uG(const double complex * VA, double complex * row_sums, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tN_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tN_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tN_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tN_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tT_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tT_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tT_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tT_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tC_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tC_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tC_sU_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tC_sU_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tN_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tN_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tN_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tN_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tT_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tT_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tT_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tT_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tC_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tC_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tC_sS_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tC_sS_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tN_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tN_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tN_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tN_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tT_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tT_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tT_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tT_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tC_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_C_u_tC_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tC_sH_dE_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);



rsb_err_t rsb__BCSR_scale_double_complex_H_u_tC_sH_dI_uG(double complex * VA, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags, const double complex *scale_factors);


#endif /* RSB_BCSS_MISC_U_H_INCLUDED */
/* @endcond */
