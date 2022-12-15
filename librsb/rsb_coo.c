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
 * This source file contains functions for COO handling.
 * */
#include "rsb_internals.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

#define RSB_MEMMOVE_BASED_DUPS_HANDLING 0

static rsb_nnz_idx_t rsb_weed_out_duplicates_unsorted(rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, void *RSB_RESTRICT VA, rsb_nnz_idx_t nnz, rsb_type_t typecode)
{
	/*!
	 * \ingroup gr_internals
	 * Weeds out duplicate coordinate elements.
	 * returns the true nnz after weeding out duplicates
	 *
	 * \note : basic, unoptimized implementation.
	 * \note : if needed, could enhance this routine by restructuring and using rsb__util_compact_marked_coo_array
	 */
	rsb_nnz_idx_t i,k,dups = 0;
	size_t el_size = 0;

	if(!IA || !JA || RSB_INVALID_NNZ_INDEX(nnz) )
		return 0;

	el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);
	RSB_DEBUG_ASSERT(el_size);

	for(k=0  ;k<nnz;++k)
	for(i=k+1;i<nnz;++i)
	{
		if( IA[k]==IA[i] && JA[k]==JA[i] )
		{
			/* this is a debug method, therefore it is stupid */
			RSB_MEMMOVE(IA+i,IA+i+1,sizeof(rsb_coo_idx_t)*(nnz-i-1));
			RSB_MEMMOVE(JA+i,JA+i+1,sizeof(rsb_coo_idx_t)*(nnz-i-1));
			/* note that it is legal and ok to move 0 for the next operation */
			RSB_MEMMOVE(
				((rsb_byte_t*)(VA))+i*el_size,
				((rsb_byte_t*)(VA))+(i+1)*el_size,el_size*(nnz-i-1));
			++dups;
			--nnz;
		}
	}
	return nnz;
}

static rsb_nnz_idx_t rsb_weed_out_duplicates_from_sorted(rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, void *RSB_RESTRICT VA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * Weeds out duplicate coordinate elements.
	 *
	 * \return the true nnz after weeding out duplicates
	 * \note : assumes input is sorted.
	 * \note : if input is not really sorted, will remove only contiguous duplicates.
	 *
	 * TODO : could suffer of overflow, because in principle rsb_coo_idx_t != rsb_nnz_idx_t .
	 * Only works for total orders (thus, no blocked orderings!).
	 * TODO : make this parallel.
	 */
	size_t el_size = 0;
#if (!RSB_MEMMOVE_BASED_DUPS_HANDLING)
	const rsb_coo_idx_t marker = RSB_MARKER_COO_VALUE; 	
	rsb_coo_idx_t fd = marker,ld = marker;    /* first duplicate sequence, last duplicate sequence */
	rsb_nnz_idx_t k = 0, dups = 0, moved = 0, moves = 0;
#endif /* RSB_MEMMOVE_BASED_DUPS_HANDLING */
	if(!IA || !JA || RSB_INVALID_NNZ_INDEX(nnz) || nnz < 2)
	{
		goto ret;
	}

	el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);

#if (!RSB_MEMMOVE_BASED_DUPS_HANDLING)

	for(k=0  ;k<nnz-1;  )
	if(IA[k]==IA[k+1] && JA[k]==JA[k+1] )
	{
		/* we found a duplicate pair */
		rsb_coo_idx_t ldc; /* local duplicates count */
		rsb_byte_t*lp = ((rsb_byte_t*)(VA))+k*el_size,*rp;
		ldc = 1;

		while( k+ldc<nnz-1 && IA[k+ldc] == IA[k+ldc+1] && JA[k+ldc] == JA[k+ldc+1] )
		{
			/* we look for more dups */
			++ldc;
		}
		rp = lp+ldc*el_size;
	//	RSB_INFO("dup: %d: %d %d (%d x)\n",k,IA[k],JA[k],ldc);
#ifdef RSB_FLAG_DUPLICATES_SUM
		//	RSB_ERROR("%d..%d..%d\n",k,ldc,nnz);
		if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_DUPLICATES_SUM))
		{
			// FIXME: missing error handling
			/*errval|=*/rsb__util_vector_sum_strided(lp,lp,typecode,ldc+1,1);
		}
		else
#endif /* RSB_FLAG_DUPLICATES_SUM */
#ifdef RSB_FLAG_DUPLICATES_KEEP_LAST
		if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_DUPLICATES_KEEP_LAST))
		{
			RSB_MEMCPY(lp,rp,el_size);
		}
		else
#endif /* RSB_FLAG_DUPLICATES_KEEP_LAST */
		{
			/* the first element is the one remaining */
		}

		if(fd==marker)
		{
			/* if there (at k) we have the first duplicates sequence, we store its index just after */
			fd = k+1;
		}
		else
		{
			/* if this is not the first one, we advertise this sequence index in JA[ld] */
			JA[ld] = k+1;
		}
		/* we write the current sequence length in I[k+1] */
		IA[k+1] = ldc;

		/* we advance */
		ld = k+1;
		k += ldc+1;
		dups += ldc;
	}
	else
		++k;

	//RSB_ERROR("! %d dups\n",dups);
	/* no dups ? nothing to do. */
	if(!dups)
		goto ret;

	/* we mark the last duplicate sequence as such */
	JA[ld] = marker;

	/* ok, we are ready for compacting the sequence */
	rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,el_size,fd,&moved,&moves);
	//RSB_INFO("%d nnz - %d dups\n",nnz,dups);
	nnz -= dups;
	goto ret;
#else /* RSB_MEMMOVE_BASED_DUPS_HANDLING */
	el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);
	/* very, very slow */
	for(k=0;RSB_LIKELY(k<nnz-1);++k)
	if(RSB_UNLIKELY(IA[k]==IA[k+1] && JA[k]==JA[k+1] ))
	{
		/* this is a debug method, therefore it is stupid */
		RSB_MEMMOVE(IA+k,IA+k+1,sizeof(rsb_coo_idx_t)*(nnz-k-1));
		RSB_MEMMOVE(JA+k,JA+k+1,sizeof(rsb_coo_idx_t)*(nnz-k-1));
		/* note that it is legal and ok to move 0 for the next operation */
		RSB_MEMMOVE(
			((rsb_byte_t*)(VA))+ k   *el_size,
			((rsb_byte_t*)(VA))+(k+1)*el_size,el_size*(nnz-k-1));
		++dups;
		--nnz;
		RSB_ERROR("dup: %d: %d %d\n",k,IA[k],JA[k]);
	}
#endif /* RSB_MEMMOVE_BASED_DUPS_HANDLING */
ret:
	return nnz;
}

rsb_nnz_idx_t rsb__weed_out_duplicates(rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, void *RSB_RESTRICT VA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * Weeds out duplicate coordinate elements.
	 *
	 * \return the true nnz after weeding out duplicates
	 */
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SORTED_INPUT))
		return rsb_weed_out_duplicates_from_sorted(IA,JA,VA,nnz,typecode,flags);
	else
		return rsb_weed_out_duplicates_unsorted(IA,JA,VA,nnz,typecode);
}

rsb_nnz_idx_t rsb__check_for_zeros(const void * VA, rsb_nnz_idx_t nnz, rsb_type_t typecode)
{
	/*!
	 * \ingroup gr_internals
	 * Checks for zero elements.
	 *
	 * Note : basic, unoptimized implementation.
	 */
	rsb_nnz_idx_t k,zeros = 0;
	size_t el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);

	if(!VA || RSB_INVALID_NNZ_INDEX(nnz) || RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		return RSB_ERR_BADARGS;
	}

	for(k=0;RSB_LIKELY(k<nnz);++k)
		if( RSB_IS_ELEMENT_ZERO(((rsb_byte_t*)VA) + k * el_size , typecode ))
		{
			++zeros;
			/* could be improved, of course */
		}
	return zeros;
}

rsb_nnz_idx_t rsb__check_for_nonzeros(const void * VA, rsb_nnz_idx_t nnz, rsb_type_t typecode)
{
	/*!
	 * \ingroup gr_internals
	 * Checks for non zero elements.
	 *
	 * Note : basic, unoptimized implementation.
	 */
	return nnz-rsb__check_for_zeros(VA,nnz,typecode);
}

rsb_err_t rsb__util_compact_marked_coo_array( rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, void *RSB_RESTRICT VA, rsb_nnz_idx_t nnz, size_t el_size, rsb_coo_idx_t fd, rsb_nnz_idx_t * movedp, rsb_nnz_idx_t * movesp)
{
	/*!
		\ingroup gr_internals
		\return the number of moved elements
		The same technique could be used for in-place BCSR element displacement.
		Compacts IA,JA,VA, deleting [x] and expecting:
		
		    ...@fd|<- 1.....D ->|   @N |< 2>|           @NNZ
		IA: ...[D][x][x]...[x][x][*][2][x][x][*][*][*]...[*]
		JA: ...[N][x][x]...[x][x][*][M][x][x][*][*][*]...[*]
		VA: ...[*][x][x]...[x][x][*][*][x][x][*][*][*]...[*]
		
		M is a marker value.
		Note: an improvement would be to make jumps in JA to be relative instead of absolute.
	*/
	rsb_nnz_idx_t k = 0,moved = 0,moves = 0;
	const rsb_coo_idx_t marker = RSB_MARKER_COO_VALUE; 	
	rsb_coo_idx_t nld = 0,ld = 0;
	rsb_byte_t* vp = VA;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_DEBUG_ASSERT(IA);
	RSB_DEBUG_ASSERT(JA );
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(fd));
	
	if(!IA || !JA || RSB_INVALID_NNZ_INDEX(nnz) )
		return RSB_ERR_BADARGS;

	for( ld=fd,k=fd ; RSB_LIKELY(JA[ld] != marker ); ld=nld )
	{
		const rsb_coo_idx_t ldc = IA[ld];       /* local marked count */
		const rsb_coo_idx_t lnd = JA[ld]-(ld+ldc);      /* local non marked */
		nld = JA[ld];

		RSB_DEBUG_ASSERT(fd >=0);
		RSB_DEBUG_ASSERT(ld >=0);
		RSB_DEBUG_ASSERT(ld <nnz);
		RSB_DEBUG_ASSERT(ldc>0);
		RSB_DEBUG_ASSERT(lnd>=0);
		
		if(0)
		RSB_INFO("k : %ld  ld : %ld  lnd : %ld  ldc : %ld  JA[ld] : %ld   nld : %ld\n", (long)k,(long)ld,(long)lnd,(long)ldc,(long)JA[ld],(long)nld);
		if(0)
		RSB_INFO("(%zd .. %zd) <- (%zd .. %zd)\n", 
			(rsb_printf_int_t)(k),
			(rsb_printf_int_t)(k+(lnd-1)),
			(rsb_printf_int_t)(ldc+ld),
			(rsb_printf_int_t)(ldc+ld+(lnd-1))
			);

		RSB_MEMMOVE(IA+k,IA+ld+ldc,(lnd) * sizeof(rsb_coo_idx_t));
		RSB_MEMMOVE(JA+k,JA+ld+ldc,(lnd) * sizeof(rsb_coo_idx_t));
		RSB_MEMMOVE(vp+(el_size*k),vp+el_size*(ld+ldc), (lnd) * el_size);
		k += (lnd);
		moved += (lnd);
		moves++;
	}
	
	if( RSB_LIKELY(JA[ld] == marker ) )
	{
		/* JA[ld]==marker (last marked sequence ) */
		const rsb_coo_idx_t ldc = IA[ld];
		const rsb_coo_idx_t lnd = (nnz-(ld+ldc));	/* local marked count, local non marked */

		if(0)
		RSB_INFO("k : %ld  ld : %ld  lnd : %ld  ldc : %ld  JA[ld] : %ld   nld : %ld\n", (long)k,(long)ld,(long)lnd,(long)ldc,(long)JA[ld],(long)nld);
	
		RSB_DEBUG_ASSERT(lnd>=0);
		if(lnd)
		{
			RSB_MEMMOVE(IA+k,IA+ld+ldc,lnd * sizeof(rsb_coo_idx_t));
			RSB_MEMMOVE(JA+k,JA+ld+ldc,lnd * sizeof(rsb_coo_idx_t));
			RSB_MEMMOVE(vp+(el_size*k),vp+el_size*(ld+ldc), lnd * el_size);
			moved += (lnd);
			moves++;

			if(0)
			RSB_INFO("(%zd .. %zd) <- (%zd .. %zd)\n", 
				(rsb_printf_int_t)(k),
				(rsb_printf_int_t)(k+(lnd-1)),
				(rsb_printf_int_t)(ldc+ld),
				(rsb_printf_int_t)(ldc+ld+(lnd-1))
				);
		}
	}
	if(movesp)
		*movesp = moves;
	if(movedp)
		*movedp = moved;

	if(0)
		RSB_STDERR("performed %zd moves, moved %zd elements out of %zd\n",(rsb_printf_int_t)moves,(rsb_printf_int_t)moved,(rsb_printf_int_t)nnz);

	RSB_DEBUG_ASSERT(moved>=0 );
	RSB_DEBUG_ASSERT(moved<=nnz);
	RSB_DEBUG_ASSERT(moves>=0 );
	RSB_DEBUG_ASSERT(moves<=nnz);
	RSB_DEBUG_ASSERT(moves<=moved);
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_util_compact_nonzeros(void *RSB_RESTRICT VA, rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp, rsb_flags_t flags  )
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Will compact the non zero input coefficients of type typecode.
	 *
	 * \param VA	a pointer to a valid coefficients array
	 * \param IA	a pointer to a valid rows coefficients array
	 * \param JA	a pointer to a valid columns coefficients array
	 * \param nnz	the coefficients count
	 * \param typecode	the coefficients typecode
	 * \param gapp	a pointer where the cut off elements number will be stored
	 * \return the number of discarded elements (0 or more) or an error code otherwise
	 *
	 * Note: this documentation is obsolete.
	 * Note : this code is slow, both algorithmically and not (it is debug stuff).
	 * TODO: shall use flags for verbosity
	 * */
	size_t el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);	/* missing unsupported typecode check */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	const rsb_coo_idx_t marker = RSB_MARKER_COO_VALUE; 	
	rsb_coo_idx_t fz = marker,lz = marker;    /* first zero sequence, last zero sequence */
	rsb_nnz_idx_t k = 0,zeros = 0,holes = 0,dzeros = 0;
	const int verbose = 0 * (flags != RSB_FLAG_NOFLAGS); /* FIXME */

	if(!VA || !IA || !JA || !gapp || RSB_INVALID_NNZ_INDEX(nnz) )
	{
		errval = RSB_ERR_BADARGS;
		{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}

	for(k=0  ;RSB_LIKELY(k<nnz);  )
	if( RSB_UNLIKELY(RSB_IS_ELEMENT_ZERO(((rsb_byte_t*)VA)+el_size*k,typecode)) )
	{
		/* we found a zero */
		rsb_coo_idx_t lzc = 0; /* local zeros count */
		int iod = (IA[k+lzc]==JA[k+lzc])?1:0;
		lzc = 1;
		dzeros += iod;

		if(verbose)
			RSB_STDOUT("zero: %ld:  r: %ld  c: %ld (%ld x, diag=%c)\n",(long int)(k+lzc),(long int)(IA[k+lzc]),(long int)(JA[k+lzc]),(long int)(lzc),iod?'y':'n');
		while( k+lzc<nnz && RSB_IS_ELEMENT_ZERO(((rsb_byte_t*)VA)+el_size*(k+lzc),typecode) )
		{
			iod = (IA[k+lzc]==JA[k+lzc])?1:0;
			/* we look for more zeros */
			if(verbose)
				RSB_STDOUT("zero: %ld:  r: %ld  c: %ld (%ld x, diag=%c)\n",(long int)(k+lzc),(long int)IA[k+lzc],(long int)JA[k+lzc],(long int)lzc,iod?'y':'n');
			++lzc;
		}
		holes += (k+1+lzc!=nnz);	/* we do not count bottom duplicates as a hole */
		
		if( RSB_UNLIKELY( fz == marker ) )
		{
			/* if this (at k) we have the first duplicate sequence, we keep its index */
			fz = k;
		}
		else
		{
			/* if this is not the first one, we advertise this sequence index in JA[lz] */
			JA[lz] = k;
		}
		/* we write the current one length in I[k] */
		IA[k] = lzc;

		/* we advance */
		lz = k;
		k += lzc;
		zeros += lzc;
	}
	else
		++k;

	/* no zeros ? nothing to do. */
	if(zeros<=0)
	{
		/* scrap scrap */
	}
	else
	{
		rsb_nnz_idx_t moved = 0,moves = 0;
		/* we mark the last zero sequence as such */
		JA[lz] = marker;

		/* ok, we are ready for compacting the sequence */
		errval = rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,el_size,fz,&moved,&moves);
/*		if(moves!=holes) // will trigger false positive error in cases like (0,0)  <- (1,1)
		{
			RSB_ERROR("%zd != %zd\n",(rsb_printf_int_t)moves,(rsb_printf_int_t)holes);
                	return RSB_ERR_INTERNAL_ERROR;
		}*/
		if(RSB_SOME_ERROR(errval))
		{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}

	if(discardedp)
		*discardedp = zeros;
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_util_compact_out_of_range(void *RSB_RESTRICT VA, rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t roff, rsb_coo_idx_t  coff, rsb_coo_idx_t Mdim, rsb_coo_idx_t mdim, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp )
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Will compact the non out of range input coefficients of type typecode
	 *
	 * \param VA	a pointer to a valid coefficients array
	 * \param IA	a pointer to a valid rows coefficients array
	 * \param JA	a pointer to a valid columns coefficients array
	 * \param nnz	the coefficients count
	 * \param typecode	the coefficients typecode
	 * \param gapp	a pointer where the cut off elements number will be stored
	 * \return the number of discarded elements (0 or more) or an error code otherwise
	 *
	 * Note that documentation is out of date. 
	 * Note : this code is slow, both algorithmically and not (it is debug stuff).
	 * */
	size_t el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);	/* missing unsupported typecode check */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	const rsb_coo_idx_t marker = RSB_MARKER_COO_VALUE; 	
	rsb_coo_idx_t fz = marker,lz = marker;    /* first zero sequence, last zero sequence */
	rsb_nnz_idx_t k = 0,zeros = 0,holes = 0;


	if(!VA || !IA || !JA || !gapp || RSB_INVALID_NNZ_INDEX(nnz) )
	{
		errval = RSB_ERR_BADARGS;
	       	RSB_PERR_GOTO(err,"Bad argument(s): VA=%p IA=%p JA=%p gapp=%p nnz=%d\n",VA,IA,JA,gapp,(int)nnz);
	}

	for(k=0  ;RSB_LIKELY(k<nnz);  )
	if(RSB_UNLIKELY( IA[k]<roff || IA[k]>=roff+Mdim || JA[k]<coff || JA[k]>=coff+mdim ))
	{
		/* we found an out of range element */
		rsb_coo_idx_t lzc; /* local out of range count */
		lzc = 1;

		while( k+lzc<nnz && ( IA[k+lzc]<roff || IA[k+lzc]>=roff+Mdim || JA[k+lzc]<coff || JA[k+lzc]>=coff+mdim ) )
		{
			/* we look for more */
			++lzc;
		}
		holes += (k+1+lzc!=nnz);	/* we do not count bottom duplicates as a hole */

	//	RSB_INFO("zero: %d: %d %d (%d x)\n",k,IA[k],JA[k],ldc);
		
		if(RSB_UNLIKELY(fz==marker))
		{
			/* if this (at k) we have the first duplicate sequence, we keep its index */
			fz = k;
		}
		else
		{
			/* if this is not the first one, we advertise this sequence index in JA[lz] */
			JA[lz] = k;
		}
		/* we write the current one length in I[k] */
		IA[k] = lzc;

		/* we advance */
		lz = k;
		k += lzc;
		zeros += lzc;
	}
	else
		++k;

	/* no zeros ? nothing to do. */
	if(zeros<=0)
	{
		/* scrap scrap */
	}
	else
	{
		rsb_nnz_idx_t moved = 0,moves = 0;
		/* we mark the last zero sequence as such */
		JA[lz] = marker;

		/* ok, we are ready for compacting the sequence */
		errval = rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,el_size,fz,&moved,&moves);
/*		if(moves!=holes) // will trigger false positive error in cases like (0,0)  <- (1,1)
		{
			RSB_ERROR("%zd != %zd\n",(rsb_printf_int_t)moves,(rsb_printf_int_t)holes);
                	return RSB_ERR_INTERNAL_ERROR;
		}*/
		if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(discardedp)
		*discardedp = zeros;
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__util_compact_nonzeros(void *RSB_RESTRICT VA, rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp, rsb_flags_t flags )
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Will sort the input coefficients of type typecode without
	 * changing their relative order, except for coefficients which are zero.
	 *
	 * \param VA	a pointer to a valid coefficients array
	 * \param IA	a pointer to a valid rows coefficients array
	 * \param JA	a pointer to a valid columns coefficients array
	 * \param nnz	the coefficients count
	 * \param typecode	the coefficients typecode
	 * \return the number of discarded elements (0 or more) or an error code otherwise
	 * 
	 * Note that documentation is out of date. 
	 * */
	if(RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		return RSB_ERR_UNSUPPORTED_TYPE;
	}
	return rsb_do_util_compact_nonzeros(VA, IA,  JA, nnz, typecode, gapp, discardedp, flags);
}

rsb_err_t rsb__weed_out_non_upptri(void *RSB_RESTRICT VA, rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp )
{
	/*!
	 * \ingroup gr_internals
	 *
	 * */
	return rsb__weed_out_non_lowtri(VA,JA,IA,nnz,typecode,gapp,discardedp);
}

rsb_err_t rsb__weed_out_non_lowtri(void *RSB_RESTRICT VA, rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp )
{
	/*!
	 * \ingroup gr_internals
	 *
	 * FIXME: remove gapp argument.
	 * */
	size_t el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);	/* missing unsupported typecode check */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	const rsb_coo_idx_t marker = RSB_MARKER_COO_VALUE; 	
	rsb_coo_idx_t fz = marker,lz = marker;    /* first zero sequence, last zero sequence */
	rsb_nnz_idx_t k = 0,zeros = 0,holes = 0;


	if(!VA || !IA || !JA /*|| !gapp*/ || RSB_INVALID_NNZ_INDEX(nnz) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	for(k=0  ;RSB_LIKELY(k<nnz);  )
	if(RSB_UNLIKELY(IA[k]<JA[k]))
	{
		/* we found a zero */
		rsb_coo_idx_t lzc; /* local zeros count */
		lzc = 1;

		while( RSB_UNLIKELY( k+lzc<nnz && IA[k+lzc] < JA[k+lzc] ) )
		{
			/* we look for more zeros */
			++lzc;
		}
		holes += (k+1+lzc!=nnz);	/* we do not count bottom duplicates as a hole */

		//RSB_INFO("up tri: %d: %d %d (%d x)\n",k,IA[k],JA[k],lzc);
		
		if(RSB_UNLIKELY(fz==marker))
		{
			/* if this (at k) we have the first duplicate sequence, we keep its index */
			fz = k;
		}
		else
		{
			/* if this is not the first one, we advertise this sequence index in JA[lz] */
			JA[lz] = k;
		}
		/* we write the current one length in I[k] */
		IA[k] = lzc;

		/* we advance */
		lz = k;
		k += lzc;
		zeros += lzc;
	}
	else
		++k;

	/* no zeros ? nothing to do. */
	if(zeros<=0)
	{
		/* scrap scrap */
	}
	else
	{
		rsb_nnz_idx_t moved = 0,moves = 0;
		/* we mark the last zero sequence as such */
		JA[lz] = marker;

		/* ok, we are ready for compacting the sequence */
		errval = rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,el_size,fz,&moved,&moves);
/*		if(moves!=holes) // will trigger false positive error in cases like (0,0)  <- (1,1)
		{
			RSB_ERROR("%zd != %zd\n",(rsb_printf_int_t)moves,(rsb_printf_int_t)holes);
                	return RSB_ERR_INTERNAL_ERROR;
		}*/
		if(RSB_SOME_ERROR(errval))
		{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}

	if(discardedp)
		*discardedp = zeros;
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__weed_out_diagonal(void *RSB_RESTRICT VA, rsb_coo_idx_t *RSB_RESTRICT IA, rsb_coo_idx_t *RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_nnz_idx_t *RSB_RESTRICT gapp, rsb_nnz_idx_t * RSB_RESTRICT discardedp )
{
	/*!
	 * \ingroup gr_internals
	 *
	 * */
	size_t el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);	/* missing unsupported typecode check */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	const rsb_coo_idx_t marker = RSB_MARKER_COO_VALUE; 	
	rsb_coo_idx_t fz = marker,lz = marker;    /* first zero sequence, last zero sequence */
	rsb_nnz_idx_t k = 0,zeros = 0,holes = 0;


	if(!VA || !IA || !JA || !gapp || (nnz!=0 && RSB_INVALID_NNZ_INDEX(nnz)) )
	{
		errval = RSB_ERR_BADARGS;
		{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}
	if(nnz==0)
	{goto err;};/* nothing to do */

	for(k=0  ;RSB_LIKELY(k<nnz);  )
	if(RSB_UNLIKELY(IA[k]==JA[k]))
	{
		/* we found a zero */
		rsb_coo_idx_t lzc; /* local zeros count */
		lzc = 1;

		while( k+lzc<nnz && IA[k+lzc] == JA[k+lzc] )
		{
			/* we look for more zeros */
			++lzc;
		}
		holes += (k+1+lzc!=nnz);	/* we do not count bottom duplicates as a hole */

		//RSB_INFO("up tri: %d: %d %d (%d x)\n",k,IA[k],JA[k],lzc);
		
		if(RSB_UNLIKELY(fz==marker))
		{
			/* if this (at k) we have the first duplicate sequence, we keep its index */
			fz = k;
		}
		else
		{
			/* if this is not the first one, we advertise this sequence index in JA[lz] */
			JA[lz] = k;
		}
		/* we write the current one length in I[k] */
		IA[k] = lzc;

		/* we advance */
		lz = k;
		k += lzc;
		zeros += lzc;
	}
	else
		++k;

	/* no zeros ? nothing to do. */
	if(zeros<=0)
	{
		/* scrap scrap */
	}
	else
	{
		rsb_nnz_idx_t moved = 0,moves = 0;
		/* we mark the last zero sequence as such */
		JA[lz] = marker;

		/* ok, we are ready for compacting the sequence */
		errval = rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,el_size,fz,&moved,&moves);
/*		if(moves!=holes) // will trigger false positive error in cases like (0,0)  <- (1,1)
		{
			RSB_ERROR("%zd != %zd\n",(rsb_printf_int_t)moves,(rsb_printf_int_t)holes);
                	return RSB_ERR_INTERNAL_ERROR;
		}*/
		if(RSB_SOME_ERROR(errval))
		{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}

	if(discardedp)
		*discardedp = zeros;
err:
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
