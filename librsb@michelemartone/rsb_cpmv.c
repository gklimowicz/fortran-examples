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
/* @cond INNERDOC  */
 /**
 * @file
 * @brief Copy/Move primitives 
 * @author Michele Martone
 * */
#include "rsb_common.h"

/*!
 An integer type for thread indices (internal).
 */
typedef rsb_thread_t rsb_tc_t;

RSB_INTERNALS_COMMON_HEAD_DECLS

void RSB_BZERO_parallel(void * p, size_t n)
{
	/**
		\ingroup gr_internals
		TODO: move to somewhere else
	 */
	rsb_char_t * cp = p;
	const rsb_tc_t wet = rsb_get_num_threads(); /* want executing threads */

	if(RSB_UNLIKELY(n<wet*RSB_MIN_THREAD_BZERO_BYTES))
	{
		RSB_BZERO(cp,n);
	}
	else
	{
		rsb_nnz_idx_t wi;
		size_t cn = (n+wet-1)/wet;	/* chunk size */
		#pragma omp parallel for schedule(static,1) RSB_NTC
		for(wi=0;wi<wet;++wi)
		{
			size_t coff = wi*cn;
			size_t cnn = (wi<wet-1)?cn:n-((wet-1)*cn);
			RSB_BZERO(cp+coff,cnn);
		}
	}
}

void RSB_A_MEMCPY_parallel(void * RSB_RESTRICT ID, const void * RSB_RESTRICT IS, size_t DOFF, size_t SOFF, size_t NNZ, size_t ES)
{
	/**
		\ingroup gr_internals
		TODO: move to somewhere else
	 */
	const rsb_tc_t wet = rsb_get_num_threads(); /* want executing threads */

	RSB_DEBUG_ASSERT(RSB_MIN_THREAD_MEMCPY_NNZ);

	if(RSB_UNLIKELY(NNZ<wet*RSB_MIN_THREAD_MEMCPY_NNZ))/* at least RSB_MIN_THREAD_MEMCPY_NNZ nnz to trigger memcpy */
	{
		RSB_A_MEMCPY(ID,IS,DOFF,SOFF,NNZ,ES);
	}
	else
	{
		rsb_nnz_idx_t wi;
		size_t cnz = (NNZ+wet-1)/wet;	/* chunk size */
		#pragma omp parallel for schedule(static,1) RSB_NTC
		for(wi=0;wi<wet;++wi)
		{
			size_t coff = wi*cnz;
			size_t cnnz = (wi<wet-1)?cnz:NNZ-((wet-1)*cnz);
			RSB_A_MEMCPY(ID,IS,DOFF+coff,SOFF+coff,cnnz,ES);
		}
	}
}

void RSB_COA_MEMCPY_parallel(void * ID, const void * IS, size_t DOFF, size_t SOFF, size_t NNZ)
{
	/**
		\ingroup gr_internals
	 */
	RSB_A_MEMCPY_parallel(ID,IS,DOFF,SOFF,NNZ,sizeof(rsb_coo_idx_t));
}

/* @endcond */
