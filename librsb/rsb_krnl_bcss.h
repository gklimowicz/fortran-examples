/* @cond INNERDOC */
/*!
 @file
 @brief

 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block compressed sparse stripes format.
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
#ifndef RSB_BCSS_H_INCLUDED
#define RSB_BCSS_H_INCLUDED
#include "rsb_internals.h"
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uaua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uaua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uauz_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uauz_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_uxua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_unua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_unua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sasa_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sasa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_uxua_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_uxua_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spmv_sxsa_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spmv_sxsa_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_spsv_sxsx_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_spsv_sxsx_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_infty_norm_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_infty_norm_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_infty_norm_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_rowssums_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_rowssums_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_rowssums_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tN_r1_c1_ul_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tT_r1_c1_ul_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_l_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tC_r1_c1_ul_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_l_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_l_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_N(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tN_r1_c1_uu_sN_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_T(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tT_r1_c1_uu_sT_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
/* a macro is faster than a switch construct */
#define RSB_double_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_float_complex_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_float_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_double_complex_kernel_dispatcher_BCSR_scale_u_C(R,C) ( (R)==(1) && (C)==(1)?  rsb__BCSR_scale_double_complex___tC_r1_c1_uu_sC_d_u \
 :  (\
 )) 

/* a macro is faster than a switch construct */
#define RSB_type_kernel_dispatcher_BCSR_scale_u_BCSR(TYPE,R,C) \
(  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE  ? (void*)RSB_double_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT  ? (void*)RSB_float_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ? (void*)RSB_float_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
  (TYPE)==RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ? (void*)RSB_double_complex_kernel_dispatcher_BCSR_scale_u_BCSR(R,C) : \
NULL ) 
#endif /* RSB_BCSS_H_INCLUDED */
/* @endcond */
