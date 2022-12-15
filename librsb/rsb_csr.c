/*

Copyright (C) 2008-2021 Michele Martone

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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains functions for CSR handling.
 * */
#include "rsb_internals.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

static rsb_err_t rsb_is_correctly_built_csr_matrix(const rsb_nnz_idx_t * PA, const rsb_coo_idx_t * JA, const rsb_coo_idx_t nrA, const rsb_coo_idx_t ncA, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t ib)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t ri;
	const rsb_coo_idx_t pb = ib;
	const rsb_coo_idx_t jb = ib; /* ib = index base */

	if(!PA ||!JA || RSB_INVALID_COO_INDEX(nrA)|| RSB_INVALID_COO_INDEX(ncA)|| RSB_INVALID_NNZ_INDEX(nnz))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"PA:%p JA:%p nrA:%zd ncA:%zd nnzA:%zd\n",PA,JA,(rsb_printf_int_t)nrA,(rsb_printf_int_t)ncA,(rsb_printf_int_t)nnz);
	}

	if(PA[nrA]!=nnz+pb)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"PA[nrA]=%zd vs nnzA=%zd (pb=%zd)\n",(rsb_printf_int_t)PA[nrA],(rsb_printf_int_t)nnz,(rsb_printf_int_t)pb);
	}

	for(ri=0;ri<nrA;++ri)
	{
		rsb_nnz_idx_t ni;

		if(PA[ri]>PA[ri+1])
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err,"PA[%zd]>PA[%zd]: %zd>%zd (row off its bounds)\n",(rsb_printf_int_t)ri,(rsb_printf_int_t)(ri+1),(rsb_printf_int_t)(PA[ri]),(rsb_printf_int_t)(PA[ri+1]));
		}
		if(PA[ri+1]-PA[ri] > ncA)
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
		for(ni=PA[ri]-pb;ni<PA[ri+1]-pb;++ni)
		{
			if(ni+1<PA[ri+1]-pb)
				if(JA[ni]>=JA[ni+1])
				{
					errval = RSB_ERR_BADARGS;
					RSB_PERR_GOTO(err,"i=%zd JA[%zd]>=JA[%zd]: %zd>=%zd (adjacent duplicates)\n",(rsb_printf_int_t)ri,(rsb_printf_int_t)ni,(rsb_printf_int_t)(ni+1),(rsb_printf_int_t)JA[ni],(rsb_printf_int_t)JA[ni+1]);
				}
			if(JA[ni]-jb>=ncA)
		       	{
				errval = RSB_ERR_BADARGS;
				RSB_PERR_GOTO(err,"i=%zd  JA[%zd]-%zd>=ncA: %zd >= %zd (column exceeding matrix)\n",(rsb_printf_int_t)ri,(rsb_printf_int_t)ni,(rsb_printf_int_t)jb,(rsb_printf_int_t)JA[ni],(rsb_printf_int_t)ncA);
			}
		}
	}
err:
        RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__csr_chk(const rsb_nnz_idx_t * RSB_RESTRICT IP, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t ib)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	errval = rsb_is_correctly_built_csr_matrix(IP, JA, nrA, ncA, nnzA, ib);
	return errval;
}

rsb_err_t rsb__csc_chk(const rsb_nnz_idx_t * RSB_RESTRICT IP, const rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t ib)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	errval = rsb__csr_chk(IP,IA,nrA,ncA,nnzA,ib);
	return errval;
}

/* @endcond */
