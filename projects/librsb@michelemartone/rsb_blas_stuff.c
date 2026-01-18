/*                                                                                                                            

Copyright (C) 2008-2019 Michele Martone

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
/* @cond INNERDOC  */
/**
 * @file
 * @author Michele Martone
 * @brief
 *
 * BLAS-like stuff
 * */
#include "rsb_blas_stuff.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

void rsb__BLAS_Xaxpy_parallel(rsb_int_t n, const void *alphap, void * a, rsb_int_t inca, const void * b, rsb_int_t incb, rsb_type_t typecode)
{
	/**
		\ingroup gr_internals
		alphap can be NULL
	 	a <- a + alpha * b
	 */
	const rsb_nnz_idx_t wet = rsb_get_num_threads(); /* want executing threads */

	if(RSB_UNLIKELY(n<wet*RSB_MIN_THREAD_XAXPY_NNZ))/* at least RSB_MIN_THREAD_MEMCPY_NNZ nnz to trigger memcpy */
	{
		rsb__cblas_Xaxpy(typecode,n,alphap,b,incb,a,inca);
	}
	else
	{
		rsb_nnz_idx_t wi, cnz = (wet+n-1)/wet;	/* chunk size */
		size_t es = RSB_SIZEOF(typecode);

		#pragma omp parallel for schedule(static,1) RSB_NTC 
		for(wi=0;wi<wet;++wi)
		{
			rsb_nnz_idx_t coff = wi*cnz;
			rsb_nnz_idx_t cnnz = (wi<wet-1)?cnz:n-((wet-1)*cnz);
			rsb__cblas_Xaxpy(typecode,cnnz,alphap,((rsb_byte_t*)b)+es*coff*incb,incb,((rsb_byte_t*)a)+es*coff*inca,inca);
		}
	}
}

void rsb__cblas_Xscal_parallel(rsb_type_t typecode, size_t n, const void * alphap, void * a, size_t stride)
{
	/**
		\ingroup gr_internals
		alphap can be NULL
	 	a <- alpha * a
	 */
	const rsb_nnz_idx_t wet = rsb_get_num_threads(); /* want executing threads */

	if(RSB_UNLIKELY(n<wet*RSB_MIN_THREAD_XAXPY_NNZ))/* at least RSB_MIN_THREAD_MEMCPY_NNZ nnz to trigger memcpy */
	{
		rsb__cblas_Xscal(typecode,n,alphap,a,stride);
	}
	else
	{
		rsb_nnz_idx_t wi,cnz = (wet+n-1)/wet;	/* chunk size */
		size_t es = RSB_SIZEOF(typecode);

		#pragma omp parallel for schedule(static,1) RSB_NTC 
		for(wi=0;wi<wet;++wi)
		{
			rsb_nnz_idx_t coff = wi*cnz;
			rsb_nnz_idx_t cnnz = (wi<wet-1)?cnz:n-((wet-1)*cnz);
			rsb__cblas_Xscal(typecode,cnnz,alphap,((rsb_byte_t*)a)+es*coff*stride,stride);
		}
	}
}

/* @endcond */
