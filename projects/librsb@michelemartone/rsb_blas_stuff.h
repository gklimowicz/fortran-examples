/*                                                                                                                            

Copyright (C) 2008-2015 Michele Martone

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
/**
 * @file
 * @author Michele Martone
 * @brief
 *
 * BLAS-like stuff
 * */
#ifndef RSB_BLAS_STUFF_H_INCLUDED
#define RSB_BLAS_STUFF_H_INCLUDED

#include "rsb_internals.h"

void rsb__BLAS_Xaxpy_parallel(rsb_int_t n, const void *alphap, void * a, rsb_int_t inca, const void * b, rsb_int_t incb, rsb_type_t typecode);
void rsb__cblas_Xscal_parallel(rsb_type_t typecode, size_t n, const void * alphap, void * a, size_t stride);

#endif /* RSB_BLAS_STUFF_H_INCLUDED */
/* @endcond */
