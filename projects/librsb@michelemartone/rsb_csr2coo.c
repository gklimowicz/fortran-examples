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
 * @brief CSR to COO conversion code
 * @author Michele Martone
 * */
#include "rsb_common.h"

void rsb__do_prefix_sum_coo_idx_t(rsb_nnz_idx_t *IA, rsb_nnz_idx_t nnz)
{
	/* FIXME: shall optimize */
	rsb_nnz_idx_t i;
	for(i=1;RSB_LIKELY(i<nnz);++i)
		IA[i] += IA[i-1];
}

rsb_err_t rsb__do_switch_fullword_array_to_compressed(rsb_nnz_idx_t *IA, rsb_nnz_idx_t nnz, rsb_nnz_idx_t m)
{
		/**
			see rsb__do_switch_compressed_array_to_fullword_coo
 			FIXME: need a no-calloc version
	 		TODO: rsb__do_switch_fullword_array_to_compressed -> rsb__idx_fia2fpa
  		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		rsb_nnz_idx_t i;
		rsb_coo_idx_t * IP = NULL;
		IP = rsb__calloc(sizeof(rsb_coo_idx_t)*(m+1));
		if(!IP)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
#if 0
		for(i=0;RSB_LIKELY(i<nnz);++i)
			if(IA[i]>=m || IA[i]<0)
			{
				errval = RSB_ERR_BADARGS;
				RSB_PERR_GOTO(err,"0 <= IA[%d]=%d < m=%d  ?\n",i,IA[i],m);
			}
#endif
		for(i=0;RSB_LIKELY(i<nnz);++i)
			IP[IA[i]+1]++;
		for(i=0;RSB_LIKELY(i<m);++i)
			IP[i+1] += IP[i];
		RSB_COA_MEMCPY(IA,IP,0,0,m+1);
err:
		RSB_CONDITIONAL_FREE(IP);
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_switch_compressed_array_to_fullword_coo(rsb_nnz_idx_t *RSB_RESTRICT IP, rsb_nnz_idx_t m, rsb_coo_idx_t off, rsb_coo_idx_t *RSB_RESTRICT TA)
{
		/**
	 		Requires m+1 temporary space.
			see rsb__do_switch_fullword_array_to_compressed
	 		TODO: rsb__do_switch_compressed_array_to_fullword_coo -> rsb__idx_fpa2fia
  		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		rsb_nnz_idx_t /*k,*/li,ri;
		//rsb_nnz_idx_t nnz = IP[m+1];
		rsb_coo_idx_t i;
		rsb_coo_idx_t * RSB_RESTRICT IA = TA;

		if(!IA)
		{
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
			#pragma omp critical (rsb__idx_fpa2fia)
#endif
			IA = rsb__malloc(sizeof(rsb_coo_idx_t)*(m+1));
		}
		if(!IA)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
		RSB_COA_MEMCPY(IA,IP,0,0,m+1);
		for(i=0;RSB_LIKELY(i<m);++i)
		{
			ri = IA[i+1];
			li = IA[i];
			rsb__util_coo_array_set(IP+li,ri-li,i+off);
		}
err:
		if(IA!=TA)
		{
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
			#pragma omp critical (rsb__idx_fpa2fia)
#endif
			RSB_CONDITIONAL_FREE(IA);
		}
	RSB_DO_ERR_RETURN(errval)
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_switch_in_place_csr_to_in_place_coo(struct rsb_mtx_t * mtxAp, rsb_bool_t do_shift)
{
	/**
		\ingroup gr_internals
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t li,ri;
	rsb_coo_idx_t i;
	// IA needs expansion
	rsb_coo_idx_t * IA = NULL;

	if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
	{
		rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*)(mtxAp->bindx),mtxAp->nnz,0);
	}
	else
	{
	}

	if(rsb__is_coo_matrix(mtxAp))
	{
		// FIXME: TODO (nothing todo)
		goto err;
	}
	IA = rsb__malloc(sizeof(rsb_coo_idx_t)*(mtxAp->Mdim+1));
	if(!IA)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
		errval = RSB_ERR_ENOMEM;
	}
	RSB_COA_MEMCPY(IA,mtxAp->bpntr,0,0,mtxAp->Mdim+1);
	for(i=0;RSB_LIKELY(i<mtxAp->Mdim);++i)
	{
		ri = IA[i+1];
		li = IA[i];
		rsb__util_coo_array_set(mtxAp->bpntr+li,ri-li,i);
	}
	if(do_shift)
	{
		// IA (mtxAp->bpntr) needs displacement of mtxAp->roff
		// JA (mtxAp->bindx) needs displacement of mtxAp->coff
		rsb__util_coo_arrays_add(mtxAp->bpntr, mtxAp->bindx, mtxAp->roff, mtxAp->coff, mtxAp->nnz);
	}
	// VA is opaque to us: no processing is needed
	RSB_CONDITIONAL_FREE(IA);
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_nnz_idx_t rsb__do_count_lowtri_in_csr(const struct rsb_coo_mtx_t *csrp)
{
	register rsb_coo_idx_t i;
	register rsb_nnz_idx_t lnz = 0;
	const rsb_coo_idx_t *IA = csrp->IA;
	const rsb_coo_idx_t *JA = csrp->JA;
	for(i=0;i<csrp->nr;++i)
	{
		register rsb_nnz_idx_t nnz0 = IA[i+0];
		register rsb_nnz_idx_t nnz1 = IA[i+1];
		lnz += rsb__nnz_split_coo_bsearch(JA+nnz0,i+1,nnz1-nnz0);
	}
	return lnz;
}

rsb_nnz_idx_t rsb__do_count_upptri_in_csr(const struct rsb_coo_mtx_t *csrp)
{
	register rsb_coo_idx_t i;
	register rsb_nnz_idx_t unz = 0;
	const rsb_coo_idx_t *IA = csrp->IA;
	const rsb_coo_idx_t *JA = csrp->JA;
	for(i=0;i<csrp->nr;++i)
	{
		register rsb_nnz_idx_t nnz0 = IA[i+0];
		register rsb_nnz_idx_t nnz1 = IA[i+1];
		unz += nnz1-nnz0-rsb__nnz_split_coo_bsearch(JA+nnz0,i,nnz1-nnz0);
	}
	return unz;
}

rsb_nnz_idx_t rsb__do_copy_lowtri_from_csr_to_coo(const struct rsb_coo_mtx_t *csrp, struct rsb_coo_mtx_t *coop)
{
	register rsb_coo_idx_t i;
	register rsb_nnz_idx_t lnz = 0;
	const rsb_coo_idx_t *IA = csrp->IA;
	const rsb_coo_idx_t *JA = csrp->JA;
	const rsb_coo_idx_t *VA = csrp->VA;
	size_t el_size = RSB_SIZEOF(csrp->typecode);
	for(i=0;i<csrp->nr;++i)
	{
		register rsb_nnz_idx_t nnz0 = IA[i+0];
		register rsb_nnz_idx_t nnz1 = IA[i+1];
		nnz1 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,i+1,nnz1-nnz0);
		RSB_CSR2COO_MEMCPY(coop->VA,coop->IA,coop->JA,VA,i,JA,lnz,nnz0,nnz1-nnz0,el_size);
		lnz += nnz1-nnz0;
	}
	return lnz;
}

rsb_nnz_idx_t rsb__do_copy_upptri_from_csr_to_coo(const struct rsb_coo_mtx_t *csrp, struct rsb_coo_mtx_t *coop)
{
	register rsb_coo_idx_t i;
	register rsb_nnz_idx_t unz = 0;
	const rsb_coo_idx_t *IA = csrp->IA;
	const rsb_coo_idx_t *JA = csrp->JA;
	const rsb_coo_idx_t *VA = csrp->VA;
	size_t el_size = RSB_SIZEOF(csrp->typecode);
	for(i=0;i<csrp->nr;++i)
	{
		register rsb_nnz_idx_t nnz0 = IA[i+0];
		register rsb_nnz_idx_t nnz1 = IA[i+1];
		nnz0 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,i,nnz1-nnz0);
		RSB_CSR2COO_MEMCPY(coop->VA,coop->IA,coop->JA,VA,i,JA,unz,nnz0,nnz1-nnz0,el_size);
		unz += nnz1-nnz0;
	}
	return unz;
}

rsb_nnz_idx_t rsb__do_count_tri_in_csr(const struct rsb_coo_mtx_t *csrp, rsb_nnz_idx_t *lnzp, rsb_nnz_idx_t *unzp)
{
	/* FIXME: should optimize */
	if(lnzp)
		*lnzp = rsb__do_count_lowtri_in_csr(csrp);
	if(unzp)
		*unzp = rsb__do_count_upptri_in_csr(csrp);
	return (lnzp?*lnzp:0)+(unzp?*unzp:0);
}
rsb_nnz_idx_t rsb__do_copy_tri_from_csr_to_coo(const struct rsb_coo_mtx_t *csrp, struct rsb_coo_mtx_t *lcoop, struct rsb_coo_mtx_t *ucoop)
{
	/* FIXME: should optimize */
	return rsb__do_copy_lowtri_from_csr_to_coo(csrp,lcoop)+rsb__do_copy_upptri_from_csr_to_coo(csrp,ucoop);
}

/* @endcond */
