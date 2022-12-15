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
 * This source file is an adaptation of Gilles Gouaillardet OpenMP+fgets_unlocked suggested implementation
 * This code is still experimental and untested.
 */
#include "rsb_internals.h"
#include "rsb_lock.h"

#ifndef RSB_OMPIO_H_INCLUDED
#define RSB_OMPIO_H_INCLUDED

#define RSB_MAX_FGETS_LINES 1536
#if RSB_WANT_OMPIO_SUPPORT
void rsb_ompio_DOUBLE (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,double**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re);void rsb_ompio_FLOAT (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,float**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re);void rsb_ompio_FLOAT_COMPLEX (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,float complex**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re);void rsb_ompio_DOUBLE_COMPLEX (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,double complex**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re);void rsb_ompio_PATTERN (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re);
#endif	/* RSB_WANT_OMPIO_SUPPORT */

#endif	/* RSB_OMPIO_H_INCLUDED */

/* @endcond */
