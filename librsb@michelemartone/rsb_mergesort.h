

/* @cond INNERDOC */
/**
 * @file
 * @brief
 * Sorting functions.
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


#ifndef RSB_MERGESORT_H_INCLUDED
#define RSB_MERGESORT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#include "rsb.h"
#include "rsb_common.h"
#include "rsb_internals.h"


#if RSB_OBSOLETE_QUARANTINE


rsb_err_t rsb__do_mergesort_CSR(
	rsb_coo_idx_t *iarray,
	rsb_coo_idx_t *jarray,
	void *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *iresult,
	rsb_coo_idx_t *jresult,
	void *result,
	rsb_type_t type);
rsb_err_t rsb__do_mergesort_BCSR(
	rsb_coo_idx_t *iarray,
	rsb_coo_idx_t *jarray,
	void *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t mb,
	rsb_coo_idx_t kb,
	rsb_coo_idx_t *iresult,
	rsb_coo_idx_t *jresult,
	void *result,
	rsb_type_t type);
rsb_err_t rsb__do_mergesort_VBR(
	rsb_coo_idx_t *iarray,
	rsb_coo_idx_t *jarray,
	rsb_coo_idx_t *biarray,
	rsb_coo_idx_t *bjarray,
	void *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *iresult,
	rsb_coo_idx_t *jresult,
	rsb_coo_idx_t *biresult,
	rsb_coo_idx_t *bjresult,
	void *result,
	rsb_type_t type);
rsb_err_t rsb_do_mergesort_double_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	double *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double *restrict result)
;

void rsb_do_merge_double_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const double* left, const double* restrict right,  double* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_float_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	float *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float *restrict result)
;

void rsb_do_merge_float_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const float* left, const float* restrict right,  float* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_float_complex_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	float complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float complex *restrict result)
;

void rsb_do_merge_float_complex_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const float complex* left, const float complex* restrict right,  float complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_double_complex_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	double complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double complex *restrict result)
;

void rsb_do_merge_double_complex_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const double complex* left, const double complex* restrict right,  double complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_double_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	double *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double *restrict result)
;

void rsb_do_merge_double_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const double* left, const double* restrict right,  double* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_float_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	float *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float *restrict result)
;

void rsb_do_merge_float_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const float* left, const float* restrict right,  float* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_float_complex_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	float complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float complex *restrict result)
;

void rsb_do_merge_float_complex_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const float complex* left, const float complex* restrict right,  float complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_double_complex_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	double complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double complex *restrict result)
;

void rsb_do_merge_double_complex_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const double complex* left, const double complex* restrict right,  double complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_double_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	double *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	double *restrict result)
;

void rsb_do_merge_double_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const double* left, const double* restrict right,  double* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_float_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	float *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	float *restrict result)
;

void rsb_do_merge_float_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const float* left, const float* restrict right,  float* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_float_complex_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	float complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	float complex *restrict result)
;

void rsb_do_merge_float_complex_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const float complex* left, const float complex* restrict right,  float complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;
rsb_err_t rsb_do_mergesort_double_complex_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	double complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	double complex *restrict result)
;

void rsb_do_merge_double_complex_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const double complex* left, const double complex* restrict right,  double complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )

;

#endif /* RSB_OBSOLETE_QUARANTINE */

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_MERGESORT_H_INCLUDED */

/* @endcond */
