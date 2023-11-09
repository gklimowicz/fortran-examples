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
 * This source file contains functions for COO handling and check.
 * */
#include "rsb_internals.h"

rsb_err_t rsb__util_is_valid_coo_array(const rsb_coo_idx_t * p, rsb_nnz_idx_t n)
{
	register rsb_nnz_idx_t k;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_DEBUG_ASSERT(p);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	for(k=0;RSB_LIKELY(k<n);++k)
		if(!RSB_IS_VALID_COO_INDEX(p[k]))
		{
			errval = RSB_ERR_GENERIC_ERROR;
			RSB_PERR_GOTO(err,"%zd : %zd\n",(rsb_printf_int_t)k,(rsb_printf_int_t)p[k]);
		}
err:
		return errval;
}

static rsb_err_t rsb__util_are_valid_coo_arrays(const rsb_coo_idx_t * p, const rsb_coo_idx_t * q, rsb_nnz_idx_t n)
{
	return
		(rsb__util_is_valid_coo_array(p,n)==RSB_ERR_NO_ERROR && 
		 rsb__util_is_valid_coo_array(q,n)==RSB_ERR_NO_ERROR ) ?
		RSB_ERR_NO_ERROR : RSB_ERR_GENERIC_ERROR;
}

#define RSB_WANT_PINFO 0

static rsb_err_t rsb__util_is_sorted_coo(const rsb_coo_idx_t *MIndx, const rsb_coo_idx_t * mIndx, rsb_nnz_idx_t nnz, rsb_type_t typecode/*, const struct rsb_mtx_partitioning_info_t * pinfop**/, rsb_flags_t flags )
{
	/*!
	 * \ingroup gr_internals
	 * A function to check if a nonzeros array is block-sorted.
	 *
	 *  FIXME : does it work with recursive ordering ?
	 * \return 0 if sorted, an error code in the other cases.
	 */
	/* Note : are you sure this is the only check for all setups ? */
	/* Note : this algorithm can be improved in plenty of ways */
	rsb_coo_idx_t i = 0,j = 0;
	rsb_nnz_idx_t k = 0;
//	const rsb_coo_idx_t *IA = NULL,*JA = NULL;
	const rsb_coo_idx_t *Mbndx = NULL,*mbndx = NULL;
#if RSB_WANT_PINFO
	rsb_coo_idx_t Mdim = 0,mdim = 0;
#endif /* RSB_WANT_PINFO */
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	cosnt rsb_bool_t want_recursive_sort = flags & RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING;
	rsb_coo_idx_t Mb = 1;
	rsb_coo_idx_t Kb = 1;
#endif

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_UNSUPPORTED_FEATURE;
	}
	
	if( !mIndx || !MIndx || nnz < 0 )
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_BADARGS;
	}

	if( 0 == (RSB_NUMERICAL_TYPE_SIZE(typecode)) )
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_UNSUPPORTED_TYPE;
	}

#if !RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS
	if(!pinfop)
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_BADARGS;
	}
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
	if(nnz<2)
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_NO_ERROR;
	}

#if RSB_WANT_PINFO
	if(pinfop)
	{
		Mbndx = pinfop->rpntr;
		mbndx = pinfop->cpntr;
		Mdim = pinfop->M_b;
		mdim = pinfop->K_b;
	}
	else
#endif /* RSB_WANT_PINFO */
	{
		Mbndx = MIndx;
		mbndx = mIndx;
	}

#if RSB_WANT_PINFO
        if(pinfop && ( !pinfop->rpntr || !pinfop->cpntr ) )
        {
                //errval = RSB_ERR_INTERNAL_ERROR;
                goto oops;
        }
#endif /* RSB_WANT_PINFO */
	
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	if(rsb__have_fixed_blocks_matrix_flags(flags) && mbndx && Mbndx)
	{
			/* FIXME */
			Kb = mbndx[1]-mbndx[0];
			Mb = Mbndx[1]-Mbndx[0];
	}

	if( want_recursive_sort && !rsb__have_fixed_blocks_matrix_flags(flags) )
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_UNIMPLEMENTED_YET;
	}
	
	if( want_recursive_sort )
	{
		/* FIXME : does not handle column transposition */
		int ml = 0, kl = 0;

		rsb_coo_idx_t Idim = (pinfop->nr+(Mb-1))/Mb;
		rsb_coo_idx_t Jdim = (pinfop->nc+(Kb-1))/Kb;

		while( (1<<ml) < Idim ) ml++;
		while( (1<<kl) < Jdim ) kl++;

		for( k=0;k<nnz-1;++k)
#if 0
		/* this is not the same! */
		rsb__asymmetric_z_index( 
			RSB_GET_BLOCK_ROW_FOR_NZ_(MIndx+k+0,Mbndx,Mdim),
			RSB_GET_BLOCK_COL_FOR_NZ_(mIndx+k+0,mbndx,mdim), pinfop->nr, pinfop->nc, ml, kl )
			>
		rsb__asymmetric_z_index( 
			RSB_GET_BLOCK_ROW_FOR_NZ_(MIndx+k+1,Mbndx,Mdim),
			RSB_GET_BLOCK_COL_FOR_NZ_(mIndx+k+1,mbndx,mdim), pinfop->nr, pinfop->nc, ml, kl ))
#else
		if(
			rsb__asymmetric_z_index((MIndx[k+0]/Mb),(mIndx[k+0]/Kb),Idim,Jdim,ml,kl)>
			rsb__asymmetric_z_index((MIndx[k+1]/Mb),(mIndx[k+1]/Kb),Idim,Jdim,ml,kl))
#endif
			goto oops;
		goto ok;
	}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

#if 0
	{
		rsb_nnz_idx_t i;
		for(i=0;i<nnz-1;++i)
		{
			if(mIndx[i]>mIndx[i+1]||(mIndx[i]==mIndx[i+1]&&MIndx[i]>MIndx[i+1]))
			{
				RSB_INFO("nnz %d : (%d,%d)\n",i  ,mIndx[i  ],MIndx[i  ]);
				RSB_INFO("nnz %d : (%d,%d)\n",i+1,mIndx[i+1],MIndx[i+1]);
				return RSB_ERR_GENERIC_ERROR;
			}
		}
	}
	else
#endif
	{
#if 1
		/* NEW : NEED COMMENTS : FIXME */
		if(nnz<1)
			goto ok;
		k = 0;
		i = 0;j = 0;

#if RSB_WANT_PINFO
		if(!pinfop)/* 1x1 */
#endif /* RSB_WANT_PINFO */
		for( k=1;k<nnz;++k )
		{
/*			RSB_DEBUG_ASSERT( MIndx[k-1] >= 0 );
			RSB_DEBUG_ASSERT( MIndx[k-1] <= MIndx[k] );
			RSB_DEBUG_ASSERT(!( mIndx[k-1] > mIndx[k] && MIndx[k-1] >= MIndx[k] ));*/

			if( MIndx[k-1] < 0 )
			{
				
				RSB_STDERR("for k=%zd\n",(rsb_printf_int_t)(k-1));
				RSB_STDERR("row index (%zd) is smaller than any one of ours\n",(size_t)MIndx[k-1]);
				goto oops1;
			}

			if( MIndx[k-1] > MIndx[k] )
			{
				RSB_STDERR("for k=%zd\n",(rsb_printf_int_t)(k-1));
				RSB_STDERR("row index (%zd) is bigger than any one of ours\n",(size_t)MIndx[k-1]);
				goto oops1;
			}

			if( mIndx[k-1] > mIndx[k] && MIndx[k-1] >= MIndx[k] )
			{
				RSB_STDERR("for k=%zd\n",(rsb_printf_int_t)(k-1));
				RSB_STDERR("col index (%zd) is bigger than any one of ours\n",(size_t)mIndx[k-1]);
				goto oops1;
			}
		}
#if RSB_WANT_PINFO
		else
		for( k=0;k<nnz;++k )
		{
			rsb_blk_idx_t li = i;

			if( MIndx[k] < Mbndx[0] )
			{
				RSB_STDERR("row index (%zd) is smaller than any one of ours (%zd)\n",(size_t)MIndx[k],(size_t)Mbndx[0]);
				goto oops;
			}

			while( i<Mdim && MIndx[k] > Mbndx[i+1] )
				++i;

			if( MIndx[k] > Mbndx[i+1] )
			{
				RSB_STDERR("row index (%zd) is bigger than any one of ours (%zd)\n",(size_t)MIndx[k],(size_t)Mbndx[i+1]);
				goto oops;
			}

			/* next block row index is ok */

			if(i>li)
				j = 0;	/* new block row */

			if( mIndx[k] < mbndx[0] )
			{
				RSB_STDERR("col index (%zd) is smaller than any one of ours (%zd)\n",(size_t)mIndx[k],(size_t)mbndx[0]);
				goto oops;
			}

			while( j<mdim && mIndx[k] > mbndx[j+1] )
				++j;

			if( mIndx[k] > mbndx[j+1] )
			{
				RSB_STDERR("col index (%zd) is bigger than any one of ours (%zd)\n",(size_t)mIndx[k],(size_t)mbndx[j+1]);
				goto oops;
			}
		}
#endif /* RSB_WANT_PINFO */
#else
		/* quite slow */
		JA = mIndx;
		IA = MIndx;
		for(i=0;i<pinfop->M_b;++i)
		for(j=0;j<pinfop->K_b;++j)
		{
			if(k>=nnz)
				goto k_nnz;/* 'premature' exit : empty last block */

			if(IA[k]<pinfop->rpntr[i] /* || ( IA[k]>=pinfop->rpntr[i] && JA[k]<pinfop->cpntr[j] )*/ )
			{
				RSB_ERROR("nnz %d : %d < %d (block row %d)\n",k,IA[k],pinfop->rpntr[i],i);
				RSB_ERROR("nnz %d : %d <?%d (block col %d)\n",k, JA[k],pinfop->cpntr[i],j);
				goto oops;/* this block should have been seen before */
			}

			/* if any, scan nnz's in this block */
			while(	k<nnz &&
				JA[k]>=pinfop->cpntr[j] && JA[k]< pinfop->cpntr[j+1] &&
				IA[k]>=pinfop->rpntr[i] &&  IA[k]< pinfop->rpntr[i+1] ) ++k;
			/* to the next block, even if this did not match */
		}
#endif
	}
	goto k_nnz;

k_nnz:
	if(k!=nnz)
	{
		RSB_STDERR("block sorting does not seem to be ok:\n");
		RSB_STDERR("empty last block ?\n");
		RSB_STDERR("element %zd %zd encountered at %zd'th (out of %zd) nnz's block (%zd %zd) (%zd - %zd , %zd - %zd)\n",
		(size_t)MIndx[k],(size_t)mIndx[k],(size_t)k,(size_t)nnz, (size_t)i,(size_t)j,(size_t)Mbndx[i],(size_t)Mbndx[i+1]-1,(size_t)mbndx[j],(size_t)mbndx[j+1]-1);
		goto err;
	}
	else
	{
		if(RSB_WANT_VERBOSE_MESSAGES)
			RSB_STDERR("block sorting seems ok\n");
			/* all ok */
	}
ok:
	return RSB_ERR_NO_ERROR;
#if RSB_WANT_PINFO
oops:
#endif /* RSB_WANT_PINFO */
	RSB_ERROR(	"block sorting does not seem to be ok:\n"
		 	"resurgent block ?\n"
			"element %zd %zd encountered at %zd'th (out of %zd) nnz's block (%zd %zd) (%zd - %zd , %zd - %zd)\n",
		(size_t)MIndx[k],(size_t)mIndx[k],(size_t)k,(size_t)nnz, (size_t)i,(size_t)j,(size_t)Mbndx[i],(size_t)Mbndx[i+1]-1,(size_t)mbndx[j],(size_t)mbndx[j+1]-1);
	goto err;
oops1:
	RSB_ERROR("block sorting does not seem to be ok..\n");
err:
	return RSB_ERR_GENERIC_ERROR;
}

rsb_err_t rsb__util_is_sorted_coo_as_row_major(const rsb_coo_idx_t *IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags )
{
	// Note: TOOD: swap names of these two functions.
	RSB_DO_FLAG_DEL(flags,RSB_INTERNAL_FLAG_CSR_SORTING_MASK);	/* NEW */

	if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
		return rsb__util_is_sorted_coo(IA,JA,nnz,typecode,flags);
	else
	{
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER);
		return rsb__util_is_sorted_coo(JA,IA,nnz,typecode,flags);
	}
}

rsb_err_t rsb__util_is_valid_coo_struct(const struct rsb_coo_mtx_t*coop)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if((!coop) || (!RSB_IS_VALID_NNZ_INDEX(coop->nnz)) || (!RSB_IS_VALID_COO_INDEX(coop->nr)) || (!RSB_IS_VALID_COO_INDEX(coop->nc)))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	if(RSB_MATRIX_UNSUPPORTED_TYPE(coop->typecode))
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	errval = rsb__util_are_valid_coo_arrays(coop->IA,coop->JA,coop->nnz);
err:
	return errval;
}

/* @endcond */
