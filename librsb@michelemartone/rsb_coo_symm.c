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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains functions for COO symmetry handling.
 * */

#include "rsb_internals.h"

rsb_err_t rsb__reallocate_with_symmetry( rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void **coo, rsb_nnz_idx_t * nnz, rsb_type_t typecode )
{
	/*!
	 * \ingroup gr_internals
	 * Assuming that for a symmetric matrix we are given 
	 * these arrays containing (i,j) pairs, with no (j,i) pairs at all (except j=i),
	 * we reallocate arrays (if possible) and fill them with the symmetric elements
	 * with no duplicate.
	 *
	 * note : this is a slow service/debug function, not a high performance one.
	 * */
	rsb_coo_idx_t * new_IA, *new_JA;
	void * ncoo;
	rsb_nnz_idx_t nnnz;
	size_t i,odel = 0, el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);/* off diagonal elements */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!IA || !JA || !*IA || !*JA || !nnz || !*nnz || !el_size)
		return RSB_ERR_BADARGS;
	
	nnnz = *nnz;
	for(i=0;i<nnnz;++i)
	{
		if((*IA)[i]!=(*JA)[i])++odel;
	}
	nnnz += odel;

	if(!odel)
	{
		return RSB_ERR_NO_ERROR;/* a diagonal matrix */
	}

	errval = rsb__util_coo_alloc(&ncoo,&new_IA,&new_JA,nnnz,typecode,RSB_BOOL_TRUE);
	if(RSB_SOME_ERROR(errval))
		return RSB_ERR_ENOMEM;
		
	RSB_COO_MEMCPY(ncoo,new_IA,new_JA,*coo,*IA,*JA,0,0,*nnz,el_size);

	odel = 0;
	for(i=0;i<*nnz;++i)
	{
		if((*IA)[i]!=(*JA)[i])
		{
			new_IA[*nnz+odel] = new_JA[i];
			new_JA[*nnz+odel] = new_IA[i];
/*			new_IA[*nnz+odel] = 1;
			new_JA[*nnz+odel] = 1;*/
			rsb__memcpy(((char*)ncoo) + (*nnz+odel) * el_size, ((const char*)ncoo) + i * el_size ,el_size);
			++odel;
			//RSB_ERROR("%d %d  ..\n",new_IA[i],new_JA[i]);
		}
	}

	rsb__free(*IA);
	rsb__free(*JA);
	rsb__free(*coo);

	*IA = new_IA;
	*JA = new_JA;
	*coo = ncoo;
	*nnz += odel;

	return errval;
}


/* @endcond */
