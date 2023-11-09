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
 * @brief Functions for sparse format switch.
 * */
#include "rsb_internals.h"		/* */
#include "rsb_swt.h"		/* */

rsb_err_t rsb__do_switch_leaf(struct rsb_mtx_t * mtxAp, rsb_fmt_t matrix_storage, rsb_flags_t flags, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t *TA)
{
	/*
	 * In place switch of rows ordered COO to COO or CSR, either halfword-compressed or not.
	 * Does not require submatrix bounds to be computed.
	 * If *TA, no allocations shall originate from here.
	 * If a reasonable conversion (e.g. no to-CSR conversion with nnzA<nrA+1) is being requested, (sizeof(rsb_coo_idx_t) * RSB_MIN(mtxMp->nnz,1+mtxMp->nr) ) should suffice for TA.
	 * TODO: this function calls OpenMP-enabled functions (e.g.: rsb__util_compress_to_row_pointers_array); fix this in an appropriate way.
	 * It's assumed that roff==coff==0 on splitting, otherwise merging is ongoing, and roff/coff are being subtracted from mtxAp->roff,mtxAp->coff.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *IA = mtxAp->bpntr, *JA = mtxAp->bindx;
	const rsb_nnz_idx_t nnzA = mtxAp->nnz;
	const rsb_coo_idx_t nrA = mtxAp->nr, ncA = mtxAp->nc;

	/* RSB_STDOUT("switch with off %d/%d, flag %d, ms %d.\n", roff, coff, flags & RSB_FLAG_USE_HALFWORD_INDICES, matrix_storage); */

	if(matrix_storage == RSB_MATRIX_STORAGE_AUTO )
	{
		matrix_storage = RSB_MATRIX_STORAGE_BCOR;
		if( nnzA >= nrA+1 && nnzA >= ncA+1 )
			matrix_storage = RSB_MATRIX_STORAGE_BCSR;
		/*
		 * May consider enabling:
		if( RSB_INDICES_FIT_IN_HALFWORD(nrA, ncA))
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES);
		*/
	}

	if(!RSB_INDICES_FIT_IN_HALFWORD(nrA, ncA))
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_USE_HALFWORD_INDICES);

	switch(matrix_storage)
	{
		case( RSB_MATRIX_STORAGE_BCSR ):	/* ... -> CSR */
			if( roff != 0 || coff != 0 )
			{
				errval = RSB_ERR_BADARGS;
				RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}
			if( nnzA < nrA+1 )
			{
				errval = RSB_ERR_BADARGS;
				RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}
			switch(flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				case(RSB_FLAG_USE_HALFWORD_INDICES): /* ... -> HCSR */

					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCSR) /* CSR -> HCSR */
					{
						/* row pointers are ok */
						if(!(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
							rsb__do_switch_array_to_halfword_coo(JA,nnzA,0);
						if( (mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
							; /* columns indices are ok */
					}
					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR) /* COO -> HCSR */
					{
						if( (mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
						{
							rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*) IA,nnzA,0);
						}
						errval = rsb__util_compress_to_row_pointers_array(TA,nnzA,nrA,RSB_FLAG_C_INDICES_INTERFACE,RSB_FLAG_C_INDICES_INTERFACE,IA);
						if(!(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
						{
							rsb__do_switch_array_to_halfword_coo(JA,nnzA,0);
						}
					}
					mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCSR;
					RSB_DO_FLAG_SUBST(mtxAp->flags,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES);
					RSB_DO_FLAG_ADD(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES_CSR));
				break;

				case(RSB_FLAG_USE_FULLWORD_INDICES):	/* -> FCSR */

					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCSR) /* CSR -> FCSR */
					{
						/* row pointers are ok */
						if(!(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
							; /* all done: CSR -> FCSR */
						if( (mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
							rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*) JA,nnzA,0); /* HCSR -> FCSR */
					}

					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR) /* COO -> FCSR */
					{
						if( (mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)) /* HCOO -> FCSR */
						{
							rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*) IA,nnzA,0); /* HCSR -> FCSR */
							rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*) JA,nnzA,0); /* HCSR -> FCSR */
						}
						/* FCOO -> FCSR */
						errval = rsb__util_compress_to_row_pointers_array(TA,nnzA,nrA,RSB_FLAG_C_INDICES_INTERFACE,RSB_FLAG_C_INDICES_INTERFACE,IA);
					}
					mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCSR;
					RSB_DO_FLAG_SUBST(mtxAp->flags,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS);
					RSB_DO_FLAG_ADD(mtxAp->flags,(RSB_FLAG_USE_CSR_RESERVED)); /* ! */
				break;
				default:
					errval = RSB_ERR_UNIMPLEMENTED_YET;
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
				break;
			}
			if(nnzA)
			{
				RSB_ASSERT( mtxAp->bpntr[0] == 0 );
				RSB_ASSERT( mtxAp->bpntr[nrA] == nnzA );
			}
		break;

		case( RSB_MATRIX_STORAGE_BCOR ): /* COO -> ... */
			switch(flags & RSB_FLAG_USE_HALFWORD_INDICES) /* COO -> H... */
			{
				case(RSB_FLAG_USE_HALFWORD_INDICES):	/* -> HCOO */
					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCSR)
					{
						errval = rsb__do_switch_compressed_array_to_fullword_coo(IA,nrA,roff,TA);
						rsb__do_switch_array_to_halfword_coo(IA,nnzA,0);
					}
					else
					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR)
					{
						if( (mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
							rsb__util_hcoo_array_add((rsb_half_idx_t*) IA,nnzA,roff);
						else
							rsb__do_switch_array_to_halfword_coo(IA,nnzA,roff);
					}

					if(!(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
					{
						rsb__do_switch_array_to_halfword_coo(JA,nnzA,coff);
					}
					else
					if( (mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES))
					{
						rsb__util_hcoo_array_add((rsb_half_idx_t*) JA,nnzA,coff);
					}
					mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCOR;
					RSB_DO_FLAG_SUBST(mtxAp->flags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES_COO);
				break;

				case(RSB_FLAG_USE_FULLWORD_INDICES):	/* -> FCOO */
					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCSR)
					{
						errval = rsb__do_switch_compressed_array_to_fullword_coo(IA,nrA,roff,TA);
					}
					if(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR)
					{
						if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
							rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*) IA,nnzA,roff);
						else
							rsb__util_coo_array_add(IA,nnzA,roff);
					}
					if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
						rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*) JA,nnzA,coff);
					else
						rsb__util_coo_array_add(JA,nnzA,coff);
					mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCOR;
					RSB_DO_FLAG_SUBST(mtxAp->flags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS);
				break;
				default:
					errval = RSB_ERR_UNIMPLEMENTED_YET;
					RSB_PERR_GOTO(ret,RSB_ERRM_ES);
				break;
			}
		break;
	}
ret:
	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_QUAD_PARTITIONING);
err:
	if(mtxAp && ( errval != RSB_ERR_NO_ERROR ) )
	{
		mtxAp->roff -= roff;
		mtxAp->coff -= coff;
		RSB_DEBUG_ASSERT(  (!mtxAp) || rsb__mtx_chk(mtxAp)==RSB_BOOL_TRUE );
	}
	return errval;
}

rsb_bool_t rsb__do_is_candidate_size_for_halfword_coo(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: get rid of this.
	 * */
#if 0
{
	rsb_coo_idx_t i,j,ij;
	i=m;
	j=k;
	ij = RSB_COO_HALFWORDS_VALUES_PACK(i,j);
	i = RSB_COO_HALFWORDS_VALUES_UNPACK_LI(ij);
	j = RSB_COO_HALFWORDS_VALUES_UNPACK_UJ(ij);
	RSB_INFO("(%d %d) -> (%d %d) (%d)\n",m,k,i,j,ij);
}
#endif
	rsb_bool_t is = RSB_BOOL_FALSE;

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES_COO))
		is=(!RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(m) && !RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(k));
	else
		is = RSB_BOOL_FALSE;
	return is;
}

rsb_bool_t rsb__do_is_candidate_size_for_halfword_csr(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: get rid of this.
	 * */
	rsb_bool_t is = RSB_BOOL_FALSE;

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES))
		is=(/*!RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(m) && */!RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(k));
	else
		is = RSB_BOOL_FALSE;
	return is;
}

rsb_bool_t rsb__do_is_candidate_size_for_halfword(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: get rid of this.
	 * */
	rsb_bool_t is = RSB_BOOL_FALSE;

	is = rsb__do_is_candidate_size_for_halfword_csr(m,k,nnz,flags) || rsb__do_is_candidate_size_for_halfword_coo(m,k,flags);
	return is;
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_is_candidate_for_fullword_coo(const struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: get rid of this.
	 * */
	rsb_bool_t is = RSB_BOOL_FALSE;

	if(!mtxAp || !rsb__is_terminal_recursive_matrix(mtxAp) || !rsb__is_css_matrix(mtxAp) /* || rsb__is_not_unsymmetric(mtxAp)*/)
		return RSB_BOOL_FALSE;

	if( RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE))
		is = RSB_BOOL_TRUE;
	else
		is = RSB_BOOL_FALSE;
	return is;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/

rsb_bool_t rsb__do_is_candidate_for_halfword_coo(const struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: get rid of this.
	 * */
	rsb_bool_t is = RSB_BOOL_FALSE;

	if(!mtxAp || !rsb__is_terminal_recursive_matrix(mtxAp) || !rsb__is_coo_matrix(mtxAp) )
		goto ret;

#if RSB_OLD_COO_CRITERIA
	if((mtxAp->nnz/mtxAp->Mdim) > RSB_CONST_MIN_NNZ_PER_ROW_FOR_COO_SWITCH)
		goto ret;
#endif

	is = rsb__do_is_candidate_size_for_halfword_coo(mtxAp->nr,mtxAp->nc,mtxAp->flags);
ret:
	return is;
}

rsb_err_t rsb__do_is_candidate_for_halfword_csr(const struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: get rid of this.
	 * */
	if(!mtxAp || !rsb__is_terminal_recursive_matrix(mtxAp) || (!rsb__is_css_matrix(mtxAp)/* || rsb__is_not_unsymmetric(mtxAp)*/
			/*&& !rsb__is_bcss_matrix(mtxAp)*/ ))
		return RSB_BOOL_FALSE;

	return rsb__do_is_candidate_size_for_halfword_csr(mtxAp->nr,mtxAp->nc,mtxAp->nnz,mtxAp->flags);
}

void rsb__do_switch_array_to_fullword_coo(rsb_half_idx_t *hp, rsb_nnz_idx_t n, const rsb_coo_idx_t off)
{
        /*! 
         * \ingroup gr_experimental
         * */
#if 0
        /* FIXME: with icc -fast, this produce bad results (on an array of length 2 with [0,1], produces zeros)! */
        rsb_coo_idx_t *p=(rsb_coo_idx_t*)hp;
        register rsb_nnz_idx_t k;
        for(k=n;k>0;--k)
                p[k-1]=(rsb_coo_idx_t) hp[k-1];
#else
#if !defined(__INTEL_COMPILER)
	register	/* with debug compile mode on, icc -O0 had problems here, too */ 
#endif /* __INTEL_COMPILER */
        rsb_nnz_idx_t k;
	
        if(n<1)
		return;
	if(off==0)
#if defined(__INTEL_COMPILER)
	/* using Intel(R) C Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 12.0.0.084 Build 20101006 we noticed a wrong operation (zeroes and/or junk ones were computed), if not using the 'novector' pragma. */
	#pragma novector
#endif /* __INTEL_COMPILER */
        for(k=n;RSB_LIKELY(k>1);--k)
        {   
                ((rsb_coo_idx_t*)hp)[k-1]=hp[k-1];
        }   
	else
#if defined(__INTEL_COMPILER)
	/* using Intel(R) C Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 12.0.0.084 Build 20101006 we noticed a wrong operation (zeroes and/or junk ones were computed), if not using the 'novector' pragma. */
	#pragma novector
#endif /* __INTEL_COMPILER */
        for(k=n;RSB_LIKELY(k>1);--k)
        {   
                ((rsb_coo_idx_t*)hp)[k-1]=off+hp[k-1];
        }   
        ((rsb_coo_idx_t*)hp)[0]=off+hp[0];
#endif
}

void rsb__do_switch_array_to_halfword_coo(rsb_coo_idx_t *p, rsb_nnz_idx_t n, const rsb_half_idx_t off)
{
	/*!
	 * \ingroup gr_internals
	 * */
	rsb_half_idx_t *hp=(rsb_half_idx_t*)p;
	register rsb_nnz_idx_t k;

	if(off)
		for(k=0;RSB_LIKELY(k<n);++k)
			hp[k]=((rsb_half_idx_t)p[k])+off;
	else
		for(k=0;RSB_LIKELY(k<n);++k)
			hp[k]=((rsb_half_idx_t)p[k]);
}

rsb_err_t rsb__do_switch_to_halfword_csr(struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp || !rsb__do_is_candidate_for_halfword_csr(mtxAp))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
/*	RSB_INFO("HCSR for %d %d\n",mtxAp->roff,mtxAp->coff); */
	rsb__do_switch_array_to_halfword_coo(mtxAp->bindx,mtxAp->nnz,0);
	RSB_DO_FLAG_SUBST(mtxAp->flags,RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_switch_to_halfword_coo(struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp || !rsb__do_is_candidate_for_halfword_coo(mtxAp))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
#if 0
  	/* RSB_INFO("HCOO for %d %d\n",mtxAp->roff,mtxAp->coff); */
	for(i=0;i<mtxAp->Mdim;++i)
	{
		for(k=mtxAp->bpntr[i];k<mtxAp->bpntr[i+1]  ;++k)
		{
		       	j=mtxAp->bindx[k];
			ij = RSB_COO_HALFWORDS_VALUES_PACK(i,j);
			mtxAp->bindx[k]=ij;
#if 0
			RSB_ASSERT(RSB_COO_HALFWORDS_VALUES_UNPACK_LI(ij)==i);
			RSB_ASSERT(RSB_COO_HALFWORDS_VALUES_UNPACK_UJ(ij)==j);
#endif
		}
	}
#else
	rsb__do_switch_array_to_halfword_coo(mtxAp->bindx,mtxAp->nnz,0);
	rsb__do_switch_array_to_halfword_coo(mtxAp->bpntr,mtxAp->nnz,0);
#endif
	RSB_DO_FLAG_SUBST(mtxAp->flags,RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES,RSB_FLAG_USE_HALFWORD_INDICES_COO);
err:
	RSB_DO_ERR_RETURN(errval)
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_switch_to_fullword_csr(struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: get rid of this.
	 * TODO:RENAME: rsb__do_switch_to_fullword_csr -> rsb__mtx_rsb2csr
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_half_idx_t *hbindx;

	if(!mtxAp || !rsb__do_is_candidate_for_halfword_csr(mtxAp))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	hbindx=(rsb_half_idx_t*)mtxAp->bindx;

	rsb__do_switch_array_to_fullword_coo(hbindx,mtxAp->nnz,0);
	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES);
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/

rsb_err_t rsb__do_switch_to_fullword_zcoo(struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * TODO: obsolete this.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	register rsb_nnz_idx_t k;

	if(!mtxAp || !rsb__do_is_candidate_for_halfword_coo(mtxAp))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	for(k=0;k<mtxAp->nnz;++k)
	{
		mtxAp->bindx[k]=RSB_COO_HALFWORDS_VALUES_UNPACK_UJ(mtxAp->bindx[k]);
	}
	mtxAp->bindx[mtxAp->nnz]=0;
	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES);
err:
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
