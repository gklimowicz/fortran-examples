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
 /**
 * @file
 * @brief Code for matrix format conversion. 
 * @author Michele Martone
 * */
#include "rsb_common.h"

rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_csr(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop)
{
	/**
		\ingroup gr_internals
		Makes sense only for in place allocated.
		On exit, pointed matrix is deallocated.
 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_ASSERT(mtxAp!=NULL);

	if(RSB_DO_TOOFEWNNZFORCSR(mtxAp->nnz,mtxAp->nr))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_sorted(mtxAp,coop);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	errval = rsb__do_switch_fullword_array_to_compressed(coop->IA,coop->nnz,coop->nr);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_csc(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop)
{
	/**
		\ingroup gr_internals
		Makes sense only for in place allocated.
		On exit, pointed matrix is deallocated.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t coo;

	RSB_ASSERT(mtxAp!=NULL);

	if(RSB_DO_TOOFEWNNZFORCSR(mtxAp->nnz,mtxAp->nc))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	RSB_INIT_CXX_FROM_MTX(&coo,mtxAp);
	coo.nr=coo.nc=0; // to have rsb__allocate_coo_matrix_t allocate nnz and not more
	if(rsb__allocate_coo_matrix_t(&coo)!=&coo)
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	rsb__util_coo_array_set(coo.IA,coo.nnz,0);
	errval = rsb__do_get_csc(mtxAp,(rsb_byte_t**)(&coo.VA),&coo.JA,&coo.IA);
	coo.nr=mtxAp->nr;
	coo.nc=mtxAp->nc;
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	coop->typecode=mtxAp->typecode;
	rsb__do_mtx_free(mtxAp);
	coop->nnz=coo.nnz;
	coop->VA=coo.VA;
	coop->IA=coo.IA;
	coop->JA=coo.JA;
	coop->nr=coo.nr;
	coop->nc=coo.nc;
err:
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
