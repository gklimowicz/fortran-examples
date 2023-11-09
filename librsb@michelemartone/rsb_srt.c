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
 * This source file contains many sorting functions, with variations on the index type and algorithm.
 **/
/*
 * 
 * FIXME : There are idiosyncracies in these sorting routines.
 *         Certain ones work with an in place algorithm, certain
 *         ones not, and certain ones need a bigger index type.
 * */
#include "rsb_internals.h"	/* rsb_coo_mtx_t	*/
#include "rsb_msort_up.h"	/* msort_up		*/
#ifdef RSB_HAVE_GSL
#include <gsl/gsl_sort.h>
#endif /* RSB_HAVE_GSL */

#define RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT 0 /* TODO: this is transition path out */
RSB_INTERNALS_COMMON_HEAD_DECLS

#if RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT
static rsb_err_t rsb__do_index_based_bcsr_msort( 
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	,void * WA, size_t wb
	);
#endif /* RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT */
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
static rsb_err_t rsb__do_index_based_recursive_bcsr_sort( 
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	);
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */
#if RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT
static rsb_err_t rsb__do_nnz_index_based_sort_and_permute(
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	);
#endif /* RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT */

rsb_err_t rsb__do_util_sortcoo(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t nr, rsb_coo_idx_t nc, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop , rsb_flags_t flags, void * WA, size_t wb)
{
	/**
	 * \ingroup gr_internals
	 * Will sort the input coefficients of type typecode.
	 *
	 * \param VA	a pointer to a valid coefficients array
	 * \param IA	a pointer to a valid rows coefficients array
	 * \param JA	a pointer to a valid columns coefficients array
	 * \param nnz	the coefficients count
	 * \param typecode	the coefficients typecode
	 * \param pinfop	valid partitioning info structure pointer or NULL
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * If pinfop is not NULL, its partitioning information will be used to
	 * sort the arrays in a blockwise fashion (i.e.: inside a block,
	 * coefficients relative order doesn't matter).
	 * 
	 * Note that the RSB_SORT_IN_PLACE flag isn' supported, as it is not suitable for mergesort
	 *
	 * */
	void *rVA=NULL;
	rsb_coo_idx_t * rIA=NULL,*rJA=NULL;
	rsb_coo_idx_t * bIA=NULL, *bJA=NULL;	/* block coordinates for each nonzero */
#if RSB__USE_MTX_PARTITIONING_INFO_T
	rsb_coo_idx_t * brIA=NULL,*brJA=NULL;
#endif
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t t = RSB_TIME_ZERO;
#if RSB__USE_MTX_PARTITIONING_INFO_T
	rsb_blk_idx_t br = 1, bc = 1;
#else
	const rsb_blk_idx_t br = 1, bc = 1;	/* default, if ! pinfop*/
#endif
	/* aliases */
#if RSB_OBSOLETE_QUARANTINE
	rsb_coo_idx_t * bMindx=NULL, *bmindx=NULL;
	rsb_coo_idx_t * brMindx=NULL,*brmindx=NULL;
#endif
//	enum rsb_op_flags_t op_flags = RSB_OP_FLAG_WANT_PARALLEL_SORT;
//	enum rsb_op_flags_t op_flags = RSB_OP_FLAG_WANT_SERIAL_SORT;
	enum rsb_op_flags_t op_flags = RSB_WANT_OMP_RECURSIVE_KERNELS?RSB_OP_FLAG_WANT_PARALLEL_SORT:RSB_OP_FLAG_WANT_SERIAL_SORT;
	
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	if( rsb_global_session_handle.asm_sort_method == 2 ) /* EXPERIMENTAL, undocumented. */
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING);
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT))
	{
		// Unsupported. This limitation is here since certain functions are in RSB_OBSOLETE_QUARANTINE .
		errval = RSB_ERR_UNIMPLEMENTED_YET;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	{
		errval = RSB_ERR_UNSUPPORTED_FEATURE;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	RSB_DEBUG_ASSERT(RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_BCSS_STORAGE));
	RSB_DEBUG_ASSERT(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT));

	if(nnz==0)
		goto err;/* a special case */

	if(!VA || !IA || !JA || RSB_INVALID_NNZ_INDEX(nnz) )
		return RSB_ERR_BADARGS;

#if !RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS 
//	if(!pinfop)return RSB_ERR_BADARGS;// pinfop is optional (as used now in the lib)
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS  */

	t = - rsb_time();
#if RSB__USE_MTX_PARTITIONING_INFO_T
	if( pinfop )
		rsb__do_get_blocking_from_pinfo(pinfop, flags, &br, &bc);
#endif

	RSB_DEBUG_ASSERT(!pinfop);
#if RSB__USE_MTX_PARTITIONING_INFO_T
	if( pinfop && ( flags & RSB_FLAG_SHOULD_DEBUG ) )
	{
		if((errval = rsb__do_is_valid_pinfo_t(pinfop))!=RSB_ERR_NO_ERROR)
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ESIIB);
		}
		else
		{
			if(RSB_WANT_VERBOSE_MESSAGES)
				RSB_STDERR("sorting input seems ok \n");
		}
	}	
#endif /* RSB__USE_MTX_PARTITIONING_INFO_T */

#if RSB__USE_MTX_PARTITIONING_INFO_T
	if( flags & RSB_FLAG_WANT_BCSS_STORAGE )
		/* FIXME : NEW . IT BREAKS NON BCSS */
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT);
#endif /* RSB__USE_MTX_PARTITIONING_INFO_T */

#if RSB__USE_MTX_PARTITIONING_INFO_T
	if( ! ( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT ) )
	if( pinfop && (( flags & RSB_FLAG_WANT_BCSS_STORAGE ) == 0) )
	{
		rsb_time_t p;
		rsb_nnz_idx_t nzi=0;
		p = - rsb_time();

		bIA    = rsb__malloc(sizeof(rsb_coo_idx_t)*(nnz));
		bJA    = rsb__malloc(sizeof(rsb_coo_idx_t)*(nnz));
		brIA   = rsb__malloc(sizeof(rsb_coo_idx_t)*(nnz));
		brJA   = rsb__malloc(sizeof(rsb_coo_idx_t)*(nnz));

		if(!brIA || !brJA ){RSB_PERR_GOTO(err,RSB_ERRM_ES)}

		for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
		{
			bJA[nzi]=RSB_GET_BLOCK_COL_FOR_NZ_(JA+nzi,pinfop->cpntr,pinfop->K_b);
			bIA[nzi]=RSB_GET_BLOCK_ROW_FOR_NZ_(IA+nzi,pinfop->rpntr,pinfop->M_b);

	                RSB_DEBUG_ASSERT(bIA[nzi]>=0);
	                RSB_DEBUG_ASSERT(bJA[nzi]>=0);
		}
		p += rsb_time();
		if(RSB_WANT_VERBOSE_MESSAGES)
			RSB_STDERR("slow pre-sorting took %lg seconds\n",p);
	}
#endif /* RSB__USE_MTX_PARTITIONING_INFO_T */

	if( ! ( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT ) )
	{
		#pragma omp critical (rsb_coo_sort)
		errval = rsb__util_coo_alloc(&rVA,&rIA,&rJA,nnz,typecode,RSB_BOOL_FALSE);
	}
#if RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT
	else
	{
		rIA = IA;
		rJA = JA;
		rVA = VA;
	}
#endif /* */

	if(!rIA || !rJA || !rVA)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

#if RSB_WANT_INDEX_BASED_SORT && defined(RSB_MATRIX_STORAGE_BCSR)
	if( ( flags & RSB_FLAG_WANT_BCSS_STORAGE ) || ( flags & RSB_FLAG_WANT_FIXED_BLOCKING_VBR ) )
	{

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
		if( flags & RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING )
		{
			errval = rsb__do_index_based_recursive_bcsr_sort(IA,JA,VA,rIA,rJA,rVA,nr,nc,br,bc,nnz,typecode,flags,op_flags);
		}
		else
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */
		{
#if RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT
			if( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT )
			{
				errval = rsb__do_index_based_bcsr_msort(rIA,rJA,rVA,nr,nc,br,bc,nnz,typecode,flags,op_flags,WA,wb);
			}
			else
#endif
			{
				errval = rsb__do_index_based_bcsr_sort(IA,JA,VA,rIA,rJA,rVA,nr,nc,br,bc,nnz,typecode,flags,op_flags,WA,wb);
			}
		}

		if(errval == RSB_ERR_NO_ERROR)
		{	
			if( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT )
				/* FIXME : in this case we'll have a segfault, because we did not allocate copies ! */
				goto sorted_in_place;
			goto sorted;
		}
		else
		if(errval == RSB_ERR_LIMITS && !( flags & RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING ) )
		{
			/* WARNING : switching back to our traditional, slower sorting */
			errval = RSB_ERR_NO_ERROR;
		}
		else
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
	}
#endif /* RSB_WANT_INDEX_BASED_SORT && defined(RSB_MATRIX_STORAGE_BCSR) */
	
#if RSB_OBSOLETE_QUARANTINE
#ifdef RSB_MATRIX_STORAGE_BCSR
	if( ( flags & RSB_FLAG_WANT_BCSS_STORAGE ) != 0)
	{
		/* FIXME : temporary fix (need ad-hoc variables)  */

		if(bc==1 && br==1)
			rsb__do_mergesort_CSR( IA, JA, VA, nnz, rIA, rJA, rVA, typecode);
		else
			rsb__do_mergesort_BCSR( IA, JA, VA, nnz, br,bc, rIA, rJA, rVA, typecode);
	}
	else
#endif /* RSB_MATRIX_STORAGE_BCSR */
	if( bIA && bJA )
		rsb__do_mergesort_VBR( IA, JA, bIA, bJA, VA, nnz, rIA, rJA, brIA, brJA, rVA, typecode);
	else
		rsb__do_mergesort_CSR( IA, JA, VA, nnz, rIA, rJA, rVA, typecode);
#else /* RSB_OBSOLETE_QUARANTINE */
	errval = RSB_ERR_INTERNAL_ERROR;
	RSB_PERR_GOTO(err,RSB_ERRM_ES)
#endif /* RSB_OBSOLETE_QUARANTINE */
sorted:

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if( ! ( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT ) )
	{
		/* we copy back the sorted arrays to the input arrays */
		RSB_COO_MEMCPY(VA,IA,JA,rVA,rIA,rJA,0,0,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode));
	}

sorted_in_place:

	RSB_DEBUG_ASSERT(!pinfop);
#if RSB__USE_MTX_PARTITIONING_INFO_T
	if( pinfop && ( flags & RSB_FLAG_SHOULD_DEBUG ) )
	{
		if((errval = rsb__do_is_valid_pinfo_t(pinfop))!=RSB_ERR_NO_ERROR)
		{
			RSB_PERR_GOTO(err,RSB_ERRM_SLSIB);
		}
		else
		{
			if(RSB_WANT_VERBOSE_MESSAGES)
				RSB_STDERR("sorting seems ok (1/2)\n");
		}
	}	
#endif /* RSB__USE_MTX_PARTITIONING_INFO_T */
	
/*	if( ( flags & RSB_FLAG_SHOULD_DEBUG ) && ! ( flags & RSB_FLAG_QUAD_PARTITIONING ) )
	{
		RSB_WARN("skipping recursive sort check!\n");
	}*/

/*	if( ( flags & RSB_FLAG_SHOULD_DEBUG ) && ! ( flags & RSB_FLAG_QUAD_PARTITIONING ) )*/
	if( flags & RSB_FLAG_SHOULD_DEBUG )
	{
		if(RSB_WANT_VERBOSE_MESSAGES)
			RSB_STDERR("just sorted. let's check\n");

		errval = rsb__util_is_sorted_coo_as_row_major(IA,  JA, nnz, typecode, pinfop, flags);

		if(RSB_SOME_ERROR(errval))
		{
			RSB_STDERR(RSB_ERRM_SLIINS);
			goto err;
		}
		else
		{
			if(RSB_WANT_VERBOSE_MESSAGES)
				RSB_STDERR("sorting seems ok (2/2)\n");
		}
	}

err:
	t += rsb_time();

	if(RSB_SOME_ERROR(errval))
		RSB_ERROR("!\n");

//	if(RSB_WANT_VERBOSE_MESSAGES)
//		RSB_STDERR("matrix sorted in %lf seconds \n", t);
	if( ! ( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT ) )
	#pragma omp critical (rsb_coo_sort)
	{
		RSB_CONDITIONAL_FREE(bIA);
		RSB_CONDITIONAL_FREE(bJA);
#if RSB__USE_MTX_PARTITIONING_INFO_T
		RSB_CONDITIONAL_FREE(brIA);
		RSB_CONDITIONAL_FREE(brJA);
#endif
		RSB_CONDITIONAL_FREE(rIA);
		RSB_CONDITIONAL_FREE(rJA);
		RSB_CONDITIONAL_FREE(rVA);
	}
	RSB_DO_ERR_RETURN(errval)
}

static int rsb_compar_coo_idx_t(const void * ap, const void * bp)
{
	/**
		\ingroup gr_internals
	*/
	const rsb_coo_idx_t a=*(rsb_coo_idx_t*)ap;
	const rsb_coo_idx_t b=*(rsb_coo_idx_t*)bp;

        return
                 ( a >  b ) ? 1 :
                 (( a == b ) ? 0 : -1);
}

static int rsb_compar_nnz_idx_t(const void * ap, const void * bp)
{
	/**
		\ingroup gr_internals
	*/
	const rsb_nnz_idx_t a=*(rsb_nnz_idx_t*)ap;
	const rsb_nnz_idx_t b=*(rsb_nnz_idx_t*)bp;

        return
                 ( a >  b ) ? 1 :
                 (( a == b ) ? 0 : -1);
}

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
static inline rsb_nnz_idx_t rsb_coo_index_bit_interleave(rsb_coo_idx_t o, rsb_coo_idx_t e)
{
	/**
		\ingroup gr_internals
		
		Interleaves two index values.
		Could be performed in one assembly instruction, if available.

		assumes sizeof(rsb_nnz_idx_t) >= sizeof(rsb_coo_idx_t)
		assumes no branches will occur, as this can be optimized at compile time.
	*/
	rsb_nnz_idx_t i = 0, O=o, E=e;

	//RSB_ASSERT(sizeof(rsb_nnz_idx_t) >= sizeof(rsb_coo_idx_t));

	RSB_DEBUG_ASSERT(O>=0);
	RSB_DEBUG_ASSERT(E>=0);

	if (sizeof(rsb_nnz_idx_t)==1)
	{
		E = (E | (E << 2)) & 0x33;
		E = (E | (E << 1)) & 0x55;
		O = (O | (O << 2)) & 0x33;
		O = (O | (O << 1)) & 0x55;
	}
	else
	if (sizeof(rsb_nnz_idx_t)==2)
	{
		E = (E | (E << 4)) & 0x0F0F;
		E = (E | (E << 2)) & 0x3333;
		E = (E | (E << 1)) & 0x5555;
		O = (O | (O << 4)) & 0x0F0F;
		O = (O | (O << 2)) & 0x3333;
		O = (O | (O << 1)) & 0x5555;
	}
	else
	if (sizeof(rsb_nnz_idx_t)==4)
	{
		E = (E | (E << 8)) & 0x00FF00FF;
		E = (E | (E << 4)) & 0x0F0F0F0F;
		E = (E | (E << 2)) & 0x33333333;
		E = (E | (E << 1)) & 0x55555555;
		O = (O | (O << 8)) & 0x00FF00FF;
		O = (O | (O << 4)) & 0x0F0F0F0F;
		O = (O | (O << 2)) & 0x33333333;
		O = (O | (O << 1)) & 0x55555555;
	}
	else
	if (sizeof(rsb_nnz_idx_t)==8)
	{
		E = (E | (E <<16)) & 0x0000FFFF0000FFFF;
		E = (E | (E << 8)) & 0x00FF00FF00FF00FF;
		E = (E | (E << 4)) & 0x0F0F0F0F0F0F0F0F;
		E = (E | (E << 2)) & 0x3333333333333333;
		E = (E | (E << 1)) & 0x5555555555555555;
		O = (O | (O <<16)) & 0x0000FFFF0000FFFF;
		O = (O | (O << 8)) & 0x00FF00FF00FF00FF;
		O = (O | (O << 4)) & 0x0F0F0F0F0F0F0F0F;
		O = (O | (O << 2)) & 0x3333333333333333;
		O = (O | (O << 1)) & 0x5555555555555555;
	}
	else
	{
		RSB_ERROR(RSB_ERRM_FYCITINS);
		/* FIXME : fatal! */
	}

	i = (E | (O << 1));
/*	if(i<0)
	{
		printf("OVERFLOW %d %d %d %d %d\n",i,e,o,E,O);
	}*/
	RSB_DEBUG_ASSERT((i & ~-1)>=0);
	return i;
}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
static void rsb_expand_coo_indices( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl, rsb_coo_idx_t * mzp, rsb_coo_idx_t * kzp)
{
	rsb_coo_idx_t mz=0,kz=0;
	rsb_coo_idx_t mh=nr,kh=nc,mc=i,kc=j,khb,mhb;
	register int im=0,ik=0;
	const register int lm=ml,lk=kl;

	RSB_DEBUG_ASSERT(mzp);
	RSB_DEBUG_ASSERT(kzp);

	while( mh >= 2 )
	{
		mhb=mh;
		mh=(mh+1)/2;
		if(mc >= mh)
		{
			mz = (mz<<1) | 1;
			mc-= mh;
			mh = mhb-mh;
		}
		else
		{
			mz = (mz<<1);
		}
		++im;
	}

	while( kh >= 2 )
	{
		khb=kh;
		kh=(kh+1)/2;
		if(kc >= kh)
		{
			kz = (kz<<1) | 1;
			kc-= kh;
			kh = khb-kh;
		}
		else
		{
			kz = (kz<<1);
		}
		++ik;
	}

//			RSB_STDERR("shifts : %d %d  %d %d\n",lm,lk,im,ik);
//			RSB_STDERR("shifts : %d %d\n",lm,lk);
//			RSB_STDERR("shifts : %d %d\n",im,ik);
//			RSB_STDERR("Z : %d %d -> %d %d : %d\n",i,j,mz,kz,-1);
#if 0
			/* FIXME : seems like REMOVING these assertions slows down the code a lot ! */
			RSB_ASSERT(lm>=im);
			RSB_ASSERT(lk>=ik);
#endif
	mz<<=(lm-im);
	kz<<=(lk-ik);
//			if(lm-im<lk-ik)mz<<=(lk-ik); else kz<<=(lm-im);
//			mz<<=(lm); kz<<=(lk);

	if(lm<lk)mz<<=(lk-lm); else kz<<=(lm-lk);

	*kzp=kz;
	*mzp=mz;
}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
static inline void rsb_asymmetric_z_indices( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl , rsb_coo_idx_t *h, rsb_coo_idx_t *l)
{
	/**
		Interleaves two rsb_coo_idx_t words, bitwise, regardless sizeof(rsb_coo_idx_t).
	in:
		e:	   o:
		+-+-+-+-+  +-+-+-+-+
		|0|2|4|6|  |1|3|5|7|
		+-+-+-+-+  +-+-+-+-+
	out:
		+-+-+-+-+  +-+-+-+-+
		|0|1|2|3|  |4|5|6|7|
		+-+-+-+-+  +-+-+-+-+

		FIXME : UNFINISHED
	*/

#if 0
	RSB_DEBUG_ASSERT(h);
	RSB_DEBUG_ASSERT(l);
	
	unsigned a;
	int b;
	a = ( (~0) << 16 ) >>16;
	b = ( (~0) << 16 ) >>16;
	printf("%x %x\n",a,b);
	a = ( (~0) << 16 ) ;
	b = ( (~0) << 16 ) ;
	printf("%x %x\n",a,b);
	a = ~ a ;
	b = ~ b ;
	printf("%x %x\n",a,b);
	a = ( (~0) >> 16 ) ;
	b = ( (~0) >> 16 ) ;
	printf("%x %x\n",a,b);

#else
	/* the compiler should be smart enough here */
	rsb_coo_idx_t mz=0,kz=0;
	const int hcb=(sizeof(rsb_coo_idx_t)*RSB_CHAR_BIT)/2;	/* half coo bytes */
	const rsb_coo_idx_t fm =~((rsb_coo_idx_t)0);	/* full bits mask */
	const rsb_coo_idx_t lm = ~(fm<<(hcb-1));
	const rsb_coo_idx_t hm = ~lm;
	const int hs = hcb-1;
	RSB_DEBUG_ASSERT(h);
	RSB_DEBUG_ASSERT(l);
	

//	printf("%x %x %x\n",fm,lm,hs);


//	printf("halfword : %x %x %x: \n",lm,~lm,1<<30);

	rsb_expand_coo_indices( i, j, nr, nc, ml, kl, &mz, &kz );
	
	/* to avoid trouble, we move the highest bit from l to h */
	/* FIXME : we deliberately ignore the highest two bits of h */
	*l = rsb_coo_index_bit_interleave( mz&lm     , kz&lm     );
	*h = rsb_coo_index_bit_interleave((mz&hm)>>hs,(kz&hm)>>hs);
	RSB_DEBUG_ASSERT(*h>=0);
	RSB_DEBUG_ASSERT(*l>=0);
#endif

}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

#if RSB_OBSOLETE_QUARANTINE
void rsb__asymmetric_z_indices_encode( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl , rsb_coo_idx_t *h, rsb_coo_idx_t *l)
{
	rsb_asymmetric_z_indices(i,j,nr,nc,ml,kl,h,l);
}
#endif /* RSB_OBSOLETE_QUARANTINE */

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
rsb_nnz_idx_t rsb__nearest_power_of_two( const rsb_nnz_idx_t n )
{
	/**
	 	\ingroup gr_internals
		\param nr a positive (nr>0) index
		\returns the nearest power of two not less than nr, if possible (higher bit unset), otherwise the highest bit
	*/
	register int bits=1;
//	RSB_DEBUG_ASSERT(i>0);

	while(n>>(bits+1))
		++bits;
	if((1<<bits) < n)
		return 1 << (bits+1);
	else
		return n;/* n was a power of two */
}
#endif

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
rsb_nnz_idx_t rsb__asymmetric_z_index( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl)
{
	/**
	 	\ingroup gr_internals

		Given (block) coordinate indices, computes a suitable asymmetric Z index for the coordinate.
		The order defined by this index is suitable to recursively partition matrices.

		FIXME : there are probably still limits and overflow related problems !
		\todo : should be inline ?
		
		\todo this routine is all about bit mangling, and it is performance critical, although non obvious.
   		

                \code 
		If sorting 2D coordinates with this index, the elements will get sorted as follows:

 			->->->->->->->->	--->--->--->--->	--------------->
 			<-<-<-<-<-<-<-<-	 /    /   /    /	            /  
 			->->->->->->->->	/    /   /    /                   /  
 			<-<-<-<-<-<-<-<-	<---<---<---<---	        /       
 			->->->->->->->->	--->--->--->--->	      /        
 			<-<-<-<-<-<-<-<-	 /    /   /    /	    /          
 			->->->->->->->->	/    /   /    /           /          
 			<-<-<-<-<-<-<-<-	<---<---<---<---	<---------------
  
                \endcode 
	*/
	rsb_coo_idx_t mz=0,kz=0;
#if 1
	rsb_expand_coo_indices( i, j, nr, nc, ml, kl, &mz, &kz );
	return rsb_coo_index_bit_interleave(mz,kz);
#else
	rsb_coo_idx_t mh=nr,kh=nc,mc=i,kc=j,khb,mhb;
	register int im=0,ik=0,lm=ml,lk=kl;

	while( mh >= 2 )
	{
		mhb=mh;
		mh=(mh+1)/2;
		if(mc >= mh)
		{
			mz = (mz<<1) | 1;
			mc-= mh;
			mh = mhb-mh;
		}
		else
		{
			mz = (mz<<1);
		}
		++im;
	}

	while( kh >= 2 )
	{
		khb=kh;
		kh=(kh+1)/2;
		if(kc >= kh)
		{
			kz = (kz<<1) | 1;
			kc-= kh;
			kh = khb-kh;
		}
		else
		{
			kz = (kz<<1);
		}
		++ik;
	}

//			RSB_STDERR("shifts : %d %d  %d %d\n",lm,lk,im,ik);
//			RSB_STDERR("shifts : %d %d\n",lm,lk);
//			RSB_STDERR("shifts : %d %d\n",im,ik);
//			RSB_STDERR("Z : %d %d -> %d %d : %d\n",i,j,mz,kz,-1);
#if 0
			/* FIXME : seems like REMOVING these assertions slows down the code a lot ! */
			RSB_ASSERT(lm>=im);
			RSB_ASSERT(lk>=ik);
#endif
	mz<<=(lm-im);
	kz<<=(lk-ik);
//			if(lm-im<lk-ik)mz<<=(lk-ik); else kz<<=(lm-im);
//			mz<<=(lm); kz<<=(lk);

	if(lm<lk)mz<<=(lk-lm); else kz<<=(lm-lk);
//			if(lm<lk)kz>>=(lk-lm); else mz>>=(lm-lk);
//			if(im<ik)kz>>=(ik-im); else mz>>=(im-ik);
//			if(im<ik)mz<<=(ik-im); else kz<<=(im-ik);
//			if(lm<lk)kz<<=(lk-ik); else mz<<=(lm-ik);
//			if(lm>lk)mz<<=(lk); else kz<<=(lm);
//			if(im>ik)kz<<=(im-ik); else mz<<=(ik-im);
			//RSB_STDERR("Z : %d %d -> %d %d \n",IA[n]/br,JA[n]/bc,mz,kz);

//			RSB_STDERR("Z : %d %d -> %d %d : %d\n",i,j,mz,kz,rsb_coo_index_bit_interleave(mz,kz));
	return rsb_coo_index_bit_interleave(mz,kz);
#endif
}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */


#if RSB_OBSOLETE_QUARANTINE
void rsb__asymmetric_z_nnz_indices( const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_coo_idx_t nr, rsb_coo_idx_t nc	, int ml, int kl, rsb_nnz_idx_t * a , rsb_nnz_idx_t * b )
{
	/**
		\note : this function was written for cases in which sizeof(rsb_coo_idx_t)==sizeof(rsb_nnz_idx_t)
	*/
	rsb_coo_idx_t mz=0,kz=0;

	RSB_DEBUG_ASSERT(sizeof(rsb_coo_idx_t)==sizeof(rsb_nnz_idx_t));
	rsb_expand_coo_indices( i, j, nr, nc, ml, kl, &mz, &kz );
}
#endif /* RSB_OBSOLETE_QUARANTINE */

static rsb_err_t rsb__do_coo_index_sort_on_rows_array_make( 
	rsb_coo_idx_t * K, const rsb_coo_idx_t * IA,
	const rsb_coo_idx_t nr, const rsb_coo_idx_t br,
	const rsb_nnz_idx_t nnz, const rsb_type_t typecode)
{
	/**
		Fill an array with block indices corresponding to given row indices.

		Block row size is fixed.
		\param br block row size
		\param nr row count
		\param nnz length of IA
		\param IA the row index array
		\param typecode
	*/
	const rsb_coo_idx_t Idim = (nr+(br-1))/br;
	rsb_nnz_idx_t n;

	if(!IA || RSB_INVALID_BLK_INDEX(br) || RSB_INVALID_COO_INDEX(nr) || RSB_INVALID_COO_INDEX(Idim) || RSB_INVALID_NNZ_INDEX(nnz))
		return RSB_ERR_BADARGS;

#if RSB__USE_MTX_PARTITIONING_INFO_T
	if(br!=1)
		for(n=0;RSB_LIKELY(n<nnz);++n)
		{
			K[2*n+0] =((IA[n]+0)/br);
			K[2*n+1] =n;

			RSB_DEBUG_ASSERT(K[2*n+0]>=0);
		}
	else
#endif
		for(n=0;RSB_LIKELY(n<nnz);++n)
		{
			K[2*n+0] =(IA[n]);
			K[2*n+1] =n;

			RSB_DEBUG_ASSERT(IA[n]<Idim);
			RSB_DEBUG_ASSERT(K[2*n+0]>=0);
		}

	return RSB_ERR_NO_ERROR;
}

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
static rsb_err_t rsb__do_double_coo_index_sort_array_make( 
	rsb_coo_idx_t * K, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t roffset,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags
	,int want_recursive_sort
	,enum rsb_op_flags_t op_flags
	/*, int want_rows_sort */)
{
	/**
		Fill an array with block indices corresponding to given row and column indices.

		TODO: clean up and document.
	*/
	const rsb_coo_idx_t Idim = (nr+(br-1))/br;
	const rsb_coo_idx_t Jdim = (nc+(bc-1))/bc;
	rsb_nnz_idx_t n;

	if( want_recursive_sort )
	{
		int ml=0, kl=0;

		while( (1<<ml) < Idim ) ml++;
		while( (1<<kl) < Jdim ) kl++;

		RSB_DEBUG_ASSERT(ml>=0);
		RSB_DEBUG_ASSERT(kl>=0);

		//op_flags = RSB_OP_FLAG_WANT_PARALLEL_SORT;
		//op_flags = RSB_OP_FLAG_WANT_SERIAL_SORT;
		/* note : integer division is quite fast .. */
		if(op_flags == RSB_OP_FLAG_WANT_PARALLEL_SORT)
		{
		if(br==1 && bc==1)
		{
#pragma omp parallel for schedule(static) RSB_NTC
			for(n=0;n<nnz;++n)
		//	for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				rsb_asymmetric_z_indices(IA[n],JA[n],Idim,Jdim,ml,kl,K+2*n,K+2*n+1);
//				printf("%d %d -> %d %d\n",mz,kz,*h,*l);
				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
				RSB_DEBUG_ASSERT(K[2*n+1]>=0);
			}
		}
		else
		{
#pragma omp parallel for schedule(static) RSB_NTC
			for(n=0;n<nnz;++n)
		//	for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				rsb_asymmetric_z_indices((IA[n])/br,(JA[n])/bc,Idim,Jdim,ml,kl,K+2*n,K+2*n+1);

				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
				RSB_DEBUG_ASSERT(K[2*n+1]>=0);
			}
		}
		}
		else
		{
		if(br==1 && bc==1)
		{
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				rsb_asymmetric_z_indices(IA[n],JA[n],Idim,Jdim,ml,kl,K+2*n,K+2*n+1);
//				printf("%d %d -> %d %d\n",mz,kz,*h,*l);
				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
				RSB_DEBUG_ASSERT(K[2*n+1]>=0);
			}
		}
		else
		{
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				rsb_asymmetric_z_indices((IA[n])/br,(JA[n])/bc,Idim,Jdim,ml,kl,K+2*n,K+2*n+1);

				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
				RSB_DEBUG_ASSERT(K[2*n+1]>=0);
			}
		}
		}
	}
	else
	{
		/* FIXME : UNFINISHED */
		return RSB_ERR_UNIMPLEMENTED_YET;
	}
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

static rsb_err_t rsb__do_nnz_index_sort_array_make( 
	rsb_nnz_idx_t * K, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t roffset,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags,
	int want_recursive_sort,
	/*, int want_rows_sort */
	enum rsb_op_flags_t op_flags
	)
{
	/**
		Fill an array with block indices corresponding to given row and column indices.
		Note: If want_recursive_sort==0 one must ensure that nr*nc fits in rsb_nnz_idx_t to avoid overflow.

		TODO: clean up and document.
	*/
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	const rsb_coo_idx_t Idim = (nr+(br-1))/br;
#endif
	const rsb_coo_idx_t Jdim = (nc+(bc-1))/bc;
	rsb_nnz_idx_t n;

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	if( want_recursive_sort )
		return RSB_ERR_INTERNAL_ERROR;
#endif
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	if( want_recursive_sort )
	{
		int ml=0, kl=0;

		while( (1<<ml) < Idim ) ml++;
		while( (1<<kl) < Jdim ) kl++;

		RSB_DEBUG_ASSERT(ml>=0);
		RSB_DEBUG_ASSERT(kl>=0);
		if(op_flags == RSB_OP_FLAG_WANT_PARALLEL_SORT)
		{
		/* note : integer division is quite fast .. */
		if(br==1 && bc==1)
		{
			#pragma omp parallel for schedule(static) RSB_NTC
		//	for(n=0;RSB_LIKELY(n<nnz);++n)
			for(n=0;n<nnz;++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				K[2*n+0]=rsb__asymmetric_z_index((IA[n]),(JA[n]),Idim,Jdim,ml,kl);
				K[2*n+1]=n;
				RSB_DEBUG_ASSERT(K[2*n+0]>=0); //if(!RSB_IS_SIGNED(rsb_coo_idx_t)) { RSB_DEBUG_ASSERT(K[2*n+0]>=0); }
			}
		}
		else
		{
			#pragma omp parallel for schedule(static) RSB_NTC
		//	for(n=0;RSB_LIKELY(n<nnz);++n)
			for(n=0;n<nnz;++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				K[2*n+0]=rsb__asymmetric_z_index((IA[n])/br,(JA[n])/bc,Idim,Jdim,ml,kl);
				K[2*n+1]=n;
				RSB_DEBUG_ASSERT(K[2*n+0]>=0); //if(!RSB_IS_SIGNED(rsb_coo_idx_t)) { RSB_DEBUG_ASSERT(K[2*n+0]>=0); }
			}
		}
		}
		else
		{
		/* note : integer division is quite fast .. */
		if(br==1 && bc==1)
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				K[2*n+0]=rsb__asymmetric_z_index((IA[n]),(JA[n]),Idim,Jdim,ml,kl);
				K[2*n+1]=n;
				RSB_DEBUG_ASSERT(K[2*n+0]>=0); //if(!RSB_IS_SIGNED(rsb_coo_idx_t)) { RSB_DEBUG_ASSERT(K[2*n+0]>=0); }
			}
		else
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				/* this is Z block sort balanced, unlike rsb_coo_index_bit_interleave  */
				K[2*n+0]=rsb__asymmetric_z_index((IA[n])/br,(JA[n])/bc,Idim,Jdim,ml,kl);
				K[2*n+1]=n;
				RSB_DEBUG_ASSERT(K[2*n+0]>=0); //if(!RSB_IS_SIGNED(rsb_coo_idx_t)) { RSB_DEBUG_ASSERT(K[2*n+0]>=0); }
			}
		}
#if 0
	else
	if( want_rows_sort )
	{
		/* FIXME : still unused */
		if(br==1)
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				K[2*n+0] =(IA[n]);
				K[2*n+1] =n;

                                RSB_DEBUG_ASSERT(IA[n]<Idim);
				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
			}
		else
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				K[2*n+0] =((IA[n]+0)/br);
				K[2*n+1] =n;

				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
			}
	}
#endif
	}
	else
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */
	{
		/* Note: the following may easily overflow if too large dimensions given. */
		RSB_DEBUG_ASSERT(roffset >= 0);
		if(br==1 && bc==1)
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				K[2*n+0] =(IA[n]-roffset);
				K[2*n+0]*= Jdim;
				K[2*n+0]+=(JA[n]);
				K[2*n+1] =n;

                                RSB_DEBUG_ASSERT(JA[n]<Jdim);

				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
			}
		else
			for(n=0;RSB_LIKELY(n<nnz);++n)
			{
				K[2*n+0] =((IA[n]-roffset)/br);
				K[2*n+0]*=Jdim;
				K[2*n+0]+=((JA[n]+0)/bc);
				K[2*n+1] =n;

				RSB_DEBUG_ASSERT(K[2*n+0]>=0);
			}
	}

	return RSB_ERR_NO_ERROR;
}

static void rsb__do_util_compact_permutation_coo_idx_t_array(rsb_coo_idx_t * K, rsb_nnz_idx_t nnz)
{
	/**
	 	\ingroup gr_internals

		Compacts a permutation vector (overwriting every first of two element with the one following).
		Requires an array of 2*nnz elements;
	*/
	rsb_nnz_idx_t n;
	/* compacting K into its first half */
	for(n=0;RSB_LIKELY(n<nnz);++n)
		K[n]=K[2*n+1];
}

void rsb__do_util_compact_permutation_nnz_idx_t_array(rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz)
{
	/**
	 	\ingroup gr_internals

		Compacts a permutation vector (overwriting every first of two element with the one following).
		Requires an array of 2*nnz elements;
	*/
	rsb_nnz_idx_t n;
	/* compacting K into its first half */
	for(n=0;RSB_LIKELY(n<nnz);++n)
		K[n]=K[2*n+1];
}

#if RSB_OBSOLETE_QUARANTINE
static rsb_err_t rsb__do_double_pass_coo_index_based_bcsr_msort( 
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags)
{
	/* FIXME */
	/**
		\ingroup gr_internals

		FIXME : FINISH ME ME
	*/
#if 0
	rsb_nnz_idx_t n1=0,n2=0;
	rsb_nnz_idx_t *K;
	size_t el_size = RSB_SIZEOF(typecode);
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	K = rsb__malloc( (nnz+1) * sizeof(rsb_nnz_idx_t)  * 2 );
	if(!K)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	errval = rsb_do_nnz_index_sort_on_rows_array_make(K,rIA,nr,br,nnz,typecode);

	errval = rsb__do_nnz_index_based_sort_and_permute(rIA,rJA,rVA,rIA,rJA,rVA,K,nnz,typecode,flags);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	RSB_CONDITIONAL_FREE(K);
	K = rsb__malloc( nc*br * sizeof(rsb_nnz_idx_t)  * 2 );
	/* FIXME ! overflows could still happen here! */

	if(!K)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	while(n1!=nnz)
	{
		/* FIXME : need specialized code here */
		while( n2+1<nnz && IA[n1]/br==IA[n2+1]/br )
			++n2;

		errval = rsb__do_nnz_index_sort_array_make(K,IA+n1,JA+n1,nr,nc,IA[n1],br,bc,(n2+1)-n1,typecode,flags,0,op_flags);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
		errval = rsb__do_nnz_index_based_sort_and_permute(IA+n1,JA+n1,RSB_VA_OFFSET_POINTER(VA,mtxAp->el_size,n1),rIA+n1,rJA+n1,RSB_VA_OFFSET_POINTER(rVA,mtxAp->el_size,n1),K,(n2+1)-n1,typecode,flags);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}

		++n2;n1=n2;
	}
err:
	RSB_CONDITIONAL_FREE(K);

	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_GENERIC_ERROR;
#endif
}
#endif /* RSB_OBSOLETE_QUARANTINE */

#if RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT
static rsb_err_t rsb__do_nnz_index_based_bcsr_msort( 
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	,void * WA, size_t wb)
{
	/**
		\ingroup gr_internals

		FIXME : could optimize a bit more!
		FIXME : should implement double pass msort!
	*/
	/* nothing to do for RSB_FLAG_WANT_COLUMN_MAJOR_ORDER :
	 the calling routine should already have swapped input arguments accordingly */
	rsb_nnz_idx_t * K=NULL;
	rsb_coo_idx_t Idim;
	rsb_coo_idx_t Jdim;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t st=0,it=0;

	RSB_DEBUG_ASSERT(nr);
	RSB_DEBUG_ASSERT(nc);
	RSB_DEBUG_ASSERT(br);
	RSB_DEBUG_ASSERT(bc);

	Idim=(nr+(br-1))/br;
	Jdim=(nc+(bc-1))/bc;

	RSB_DEBUG_ASSERT(Idim>0);
	RSB_DEBUG_ASSERT(Jdim>0);

	if(nnz<2)
		goto err;

	if(br<1 || bc<1 || nr<1 || nc<1)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	
	/* this check implies a cast to rsb_nnz_idx_t */
	if(
	 RSB_NNZ_MUL_OVERFLOW(Idim,Jdim) != 0  ||
	 RSB_NNZ_MUL_OVERFLOW(Idim,Idim) != 0  ||
	 RSB_NNZ_MUL_OVERFLOW(Jdim,Jdim) != 0 
	)

	{
		/* should resort to a double pass algorithm */
		return RSB_ERR_LIMITS;
	}

	if(WA && wb >= RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT(nnz,nr,nc,br,bc))
		K=WA;
	else
		K = rsb__malloc(RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT(nnz,nr,nc,br,bc));

	if(!K)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	it = - rsb_time();
	errval = rsb__do_nnz_index_sort_array_make(K,rIA,rJA,nr,nc,0,br,bc,nnz,typecode,flags,0,op_flags);
	it += rsb_time();
	
	if(RSB_SOME_ERROR(errval))
	{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	st = - rsb_time();
	errval = rsb__do_nnz_index_based_sort_and_permute(rIA,rJA,rVA,rIA,rJA,rVA,K,nnz,typecode,flags,op_flags);
	st += rsb_time();

	if(RSB_SOME_ERROR(errval))
	{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

err:
	if(WA && wb >= RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT(nnz,nr,nc,br,bc))
		;
	else
		RSB_CONDITIONAL_FREE(K);

	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT */

#if RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT
static rsb_err_t rsb__do_index_based_bcsr_msort( 
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	,void * WA, size_t wb
	)
{
	/**
		\ingroup gr_internals
	*/
	/* nothing to do for RSB_FLAG_WANT_COLUMN_MAJOR_ORDER :
	 the calling routine should already have swapped input arguments accordingly */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( rsb_global_session_handle.asm_sort_method != 3 ) // internal, for test purposes only (FIXME)
		errval = rsb__do_nnz_index_based_bcsr_msort(rIA,rJA,rVA,nr,nc,br,bc,nnz,typecode,flags,op_flags,WA,wb);

	if(errval == RSB_ERR_LIMITS || rsb_global_session_handle.asm_sort_method == 3)
	{
		void *tVA = NULL;
		rsb_coo_idx_t * tIA = NULL, *tJA = NULL;

		errval = rsb__util_coo_alloc(&rVA,&rIA,&rJA,nnz,typecode,RSB_BOOL_FALSE);
		if(!RSB_SOME_ERROR(errval))
		{
			errval = rsb__do_index_based_bcsr_sort(rIA,rJA,rVA,tIA,tJA,tVA,nr,nc,br,bc,nnz,typecode,flags,op_flags,WA,wb);
			RSB_COO_MEMCPY(rVA,rIA,rJA,tVA,tIA,tJA,0,0,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode));
			//return rsb__do_double_pass_coo_index_based_bcsr_msort( rIA, rJA, rVA, rIA, rJA, rVA, nr, nc, br, bc, nnz, typecode, flags);
		}
		RSB_CONDITIONAL_FREE(tIA);
		RSB_CONDITIONAL_FREE(tJA);
		RSB_CONDITIONAL_FREE(tVA);
	}

	if(RSB_SOME_ERROR(errval))
		RSB_ERROR(RSB_ERRM_ES);
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT */

#if RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT
static rsb_err_t rsb__do_nnz_index_based_sort_and_permute(
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	)
{
	/**
	 	\ingroup gr_internals

		Sort and permute with a rsb_nnz_idx_t index-based permutation array.
	*/

	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	const rsb_bool_t want_msort=1;	/* FIXME : should choose it in some other way */

	if(!want_msort)
	{
#ifdef RSB_HAVE_GSL
		gsl_heapsort( K , (size_t) nnz, 2*sizeof(rsb_nnz_idx_t), &rsb_compar_nnz_idx_t );
#else /* RSB_HAVE_GSL */
		qsort( K , (size_t) nnz, 2*sizeof(rsb_nnz_idx_t), &rsb_compar_nnz_idx_t );
#endif /* RSB_HAVE_GSL */
		rsb__do_util_compact_permutation_nnz_idx_t_array(K, nnz);

		if(flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT)
			errval = rsb__do_permute_values_in_place_with_nnz_index(
				rVA,rIA,rJA,K,nnz,typecode);
		else
			errval = rsb__do_permute_values_with_nnz_index(
				rVA, VA, rIA, IA, rJA, JA, K, nnz, typecode);
	}
	else
#endif
	{
		rsb_bool_t was_already_sorted=0;
		rsb__do_util_compact_permutation_nnz_idx_t_array(K+1, nnz-1);	/* FIXME : a hack ! */

		was_already_sorted = RSB_SOME_ERROR(rsb__do_msort_up(nnz,K,K+nnz))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;

		if(flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT)
			;/* if in place, no data to copy */
		else
		{
			/* if not in place, we copy first */
			RSB_COO_MEMCPY(rVA,rIA,rJA,VA,IA,JA,0,0,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode));
		}

		if(!was_already_sorted)
			rsb__ip_reord(nnz, rVA, rIA, rJA, K+nnz, typecode);
	}
	
	if(RSB_SOME_ERROR(errval))
	{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB__TOLERATE__RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT */

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
static rsb_err_t rsb__do_index_based_recursive_bcsr_sort( 
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags
	,enum rsb_op_flags_t op_flags
	)
{
	/**
		\ingroup gr_internals
		
		An index based coordinate sorting routine.
		Usually faster than merge sort.
		Will allocate 2 * nnz * sizeof(rsb_nnz_idx_t) bytes for a permutation vector.
		FIXME : in some cases will alocate 3 * nnz * sizeof(rsb_nnz_idx_t) bytes for a permutation vector.

		\attention : limited to smaller matrices (will exit, in case) due to overflow problems.
		\todo : it could be modified to work around the potential overflow problem, but
			should need an estimate of the maximum nnz per row amount.
		FIXME : needs more error checks (e.g: overflow of 2*nnz index .. )
		FIXME : it will believe no overflow is possible if last positive bit is found set
	*/
	/* nothing to do for RSB_FLAG_WANT_COLUMN_MAJOR_ORDER :
	 the calling routine should already have swapped input arguments accordingly */
	rsb_nnz_idx_t * K=NULL;
	rsb_coo_idx_t Idim;
	rsb_coo_idx_t Jdim;
	//rsb_nnz_idx_t n;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t it=0;
#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
	rsb_nnz_idx_t nIdim,nJdim;
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

	RSB_DEBUG_ASSERT(nr);
	RSB_DEBUG_ASSERT(nc);
	RSB_DEBUG_ASSERT(br);
	RSB_DEBUG_ASSERT(bc);

	Idim=(nr+(br-1))/br;
	Jdim=(nc+(bc-1))/bc;

	RSB_DEBUG_ASSERT(Idim>0);
	RSB_DEBUG_ASSERT(Jdim>0);

	if( ! ( flags & RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING ) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(nnz<2)
	{
		if(nnz==1 && rIA && rJA && rVA)
		{
			/* nothing to sort; only one copy is needed. */
			rIA[0]=IA[0];
			rJA[0]=JA[0];
			rsb__memcpy(rVA,VA,RSB_NUMERICAL_TYPE_SIZE(typecode));
		}
		/* now it is ok to return */
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(br<1 || bc<1 || nr<1 || nc<1)
	{
		/* FIXME : use macros for this check */
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	
	/*
		FIXME : coo overflow is NOT a menace.
		nnz_index is THE problem.
	*/
	if(0)
	{
	size_t a=Idim,b=Jdim,c=a*b;
	printf("OFLOW ? %ld %ld -> %ld  %d %d %d  %zd %zd %zd\n",(long)(Idim),(long)Jdim,(long)(Idim*Jdim),
//		(size_t)Idim*(size_t)(Jdim),
//		(size_t)Idim*(size_t)(Jdim)== (size_t)(RSB_MAX_MATRIX_DIM),
		(Idim*Jdim)== (size_t)(RSB_MAX_MATRIX_DIM),
		(Idim*Jdim) < (size_t)(RSB_MAX_MATRIX_DIM),
		(size_t)((size_t)Idim)*((size_t)Jdim) > (size_t)(Idim*Jdim),
		(size_t)(((size_t)Idim)*((size_t)Jdim)) , (size_t)(Idim*Jdim),
		c
		);
	}

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING || RSB_WANT_DO_LOCK_TEST /* like RSB_OBSOLETE_QUARANTINE */
	nIdim = rsb__nearest_power_of_two(Idim);
	nJdim = rsb__nearest_power_of_two(Jdim);

	if(
	 rsb_global_session_handle.asm_sort_method == 2 ||
	 nIdim<Idim || nJdim<Jdim || /* FIXME: this is the check that should not be */
	 RSB_NNZ_MUL_OVERFLOW(nIdim,nJdim) != 0  ||
	 RSB_NNZ_MUL_OVERFLOW(nIdim,nIdim) != 0  ||
	 RSB_NNZ_MUL_OVERFLOW(nJdim,nJdim) != 0 
	)
	{
		/* FIXME : NEW */
		goto double_coo_index;

/*
		errval = RSB_ERR_LIMITS;
		RSB_PERR_GOTO(err,"ERROR : index overflow\n");*/
	}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

	if( RSB_NNZ_MUL_OVERFLOW(Idim,Jdim) != 0 )/* this check implies a cast to rsb_nnz_idx_t */
	{
		/* overflow. should work around this. */
/*		RSB_INFO("NO OVERFLOW ? : %ld * %ld = %ld (%ld)  (%zd), Idim=%ld Jdim=%ld\n",
			(long)nr,(long)nc,(long)(nr*nc),(long)(nr)*(long)(nc),(size_t)((size_t)nr)*((size_t)nc),(long)Idim,(long)Jdim);
*/
		errval = RSB_ERR_LIMITS;
/*
		RSB_ERROR(RSB_ERRM_WOPSTASA);
*/
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	goto single_nnz_index;

single_nnz_index:

	K = rsb__malloc( ((nnz+2)+nnz) * sizeof(rsb_nnz_idx_t) );

	if(!K)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	it = - rsb_time();
	errval = rsb__do_nnz_index_sort_array_make(K,IA,JA,nr,nc,0,br,bc,nnz,typecode,flags,1,op_flags);
	it += rsb_time();

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	errval = rsb__do_nnz_index_based_sort_and_permute(IA,JA,VA,rIA,rJA,rVA,K,nnz,typecode,flags,op_flags);

	goto done;

#if RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING /* like RSB_OBSOLETE_QUARANTINE */
double_coo_index:
	/* NEW : for msort_up2 only */
	K = rsb__malloc( (nnz+2) * sizeof(rsb_nnz_idx_t) + nnz * ( 2 * sizeof(rsb_coo_idx_t) ) );

	if(!K)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	it = - rsb_time();
	errval = rsb__do_double_coo_index_sort_array_make((rsb_coo_idx_t*)K,IA,JA,nr,nc,0,br,bc,nnz,typecode,flags,1,op_flags);
	it += rsb_time();

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	{
		rsb_bool_t was_already_sorted=0;

		was_already_sorted = rsb__do_msort_up2coo(nnz,K,((rsb_coo_idx_t*)K)+2*nnz);

		if(flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT)
			;/* if in place, no data to copy */
		else
		{
			/* if not in place, we copy first */
			RSB_COO_MEMCPY(rVA,rIA,rJA,VA,IA,JA,0,0,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode));
		}

		if(!was_already_sorted)
			rsb__ip_reord(nnz, rVA, rIA, rJA, ((rsb_coo_idx_t*)K)+2*nnz, typecode);
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */
	/* FIXME : UNFINISHED */

	goto done;
done:
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
err:
	RSB_CONDITIONAL_FREE(K);

	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING */

rsb_err_t rsb__do_index_based_bcsr_sort( 
	rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_coo_idx_t br, rsb_coo_idx_t bc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode,
	rsb_flags_t flags,
	enum rsb_op_flags_t op_flags,
	void * WA, size_t wb
	)
{
	/**
		\ingroup gr_internals
		
		An index based coordinate sorting routine.
		Usually faster than merge sort.
		Will allocate 2 * nnz * sizeof(rsb_nnz_idx_t) bytes for a permutation vector.
		Can use optionally specified working area WA fitting wb bytes (see code below for wb requirements).


		Input shall be in  IA, JA, VA
		Ouput shall be in rIA,rJA,rVA
		Note that IA, JA, VA is used as swap area and therefore overwritten.
		Aliasing forbidden.

		\attention : limited to smaller matrices (will exit, in case) due to overflow problems.
		\todo : it could be modified to work around the potential overflow problem, but
			should need an estimate of the maximum nnz per row amount.
		FIXME : needs more error checks (e.g: overflow of 2*nnz index .. )
		NOTE: Only reason this is not static is the testing driver.
		TODO: it's ugly and unclear; shall remedy.
		NOTE: nothing to do for RSB_FLAG_WANT_COLUMN_MAJOR_ORDER: the calling routine should already have swapped input arguments accordingly
	*/
	rsb_nnz_idx_t * KP=NULL;
	rsb_coo_idx_t * CP=NULL;
	rsb_coo_idx_t Idim;
	rsb_coo_idx_t Jdim;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t st = RSB_TIME_ZERO, pt = RSB_TIME_ZERO, it = RSB_TIME_ZERO;
	rsb_bool_t want_two_pass_sort = 0;	/* cannot be combined with recursive sort for now */

	if( flags & RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	RSB_DEBUG_ASSERT(rVA!=VA);
	RSB_DEBUG_ASSERT(rIA!=IA);
	RSB_DEBUG_ASSERT(rJA!=JA);


	if(nnz<2)
	{
		if( nnz==1 && rIA && rJA && rVA && rIA!=IA && rJA!=JA && rVA!=VA )
		{
			/* nothing to sort; only one copy is needed. */
			rIA[0]=IA[0];
			rJA[0]=JA[0];
			rsb__memcpy(rVA,VA,RSB_NUMERICAL_TYPE_SIZE(typecode));
		}
		/* now it is ok to return */
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(br<1 || bc<1 || nr<1 || nc<1)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	
	RSB_DEBUG_ASSERT(nr);
	RSB_DEBUG_ASSERT(nc);
	RSB_DEBUG_ASSERT(br);
	RSB_DEBUG_ASSERT(bc);

	Idim=(nr+(br-1))/br;
	Jdim=(nc+(bc-1))/bc;

	RSB_DEBUG_ASSERT(Idim>0);
	RSB_DEBUG_ASSERT(Jdim>0);

	if(
	 RSB_NNZ_MUL_OVERFLOW(Idim,Jdim) != 0  ||
	 RSB_NNZ_MUL_OVERFLOW(Idim,Idim) != 0  ||
	 RSB_NNZ_MUL_OVERFLOW(Jdim,Jdim) != 0 
	)
		want_two_pass_sort = 1;

	if(want_two_pass_sort)
	{
		const size_t mwb = RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_TWO_PASS(nnz,nr,nc,br,bc);

		if(WA && wb >= mwb)
			CP = WA;
		else
			#pragma omp critical (rsb_coo_sort)
			CP = rsb__malloc(mwb);
		if(!CP)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
		it = - rsb_time();
		errval = rsb__do_coo_index_sort_on_rows_array_make(CP,IA,nr,br,nnz,typecode);
		it += rsb_time();
	}
	else
	{
		const size_t mwb = RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_ONE_PASS(nnz,nr,nc,br,bc);

		if(WA && wb >= mwb)
			KP = WA;
		else
			#pragma omp critical (rsb_coo_sort)
			KP = rsb__malloc(mwb);
		if(!KP)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
		it = - rsb_time();
		errval = rsb__do_nnz_index_sort_array_make(KP,IA,JA,nr,nc,0,br,bc,nnz,typecode,flags,0,op_flags);
		it += rsb_time();
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	st = - rsb_time();
#ifdef RSB_HAVE_GSL
	/* uhm, slow */
	if(want_two_pass_sort)
		gsl_heapsort( CP , (size_t) nnz, 2*sizeof(rsb_coo_idx_t), &rsb_compar_coo_idx_t );
	else
		gsl_heapsort( KP , (size_t) nnz, 2*sizeof(rsb_nnz_idx_t), &rsb_compar_nnz_idx_t );
#else /* RSB_HAVE_GSL */
	if(want_two_pass_sort)
		qsort( CP , (size_t) nnz, 2*sizeof(rsb_coo_idx_t), &rsb_compar_coo_idx_t );
	else
		qsort( KP , (size_t) nnz, 2*sizeof(rsb_nnz_idx_t), &rsb_compar_nnz_idx_t );
#endif /* RSB_HAVE_GSL */
	st += rsb_time();

	pt = - rsb_time();
	/* compacting CP or KP into its first half (discard comparison index: only ordering remains) */
	if(want_two_pass_sort)
		rsb__do_util_compact_permutation_coo_idx_t_array(CP, nnz);
	else
		rsb__do_util_compact_permutation_nnz_idx_t_array(KP , nnz);

	if(want_two_pass_sort)
	{
		if(flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT)
		{
			errval = rsb__do_permute_values_in_place_with_coo_index(rVA, rIA, rJA, CP, nnz, typecode);
			RSB_COO_MEMCPY(VA,IA,JA,rVA,rIA,rJA,0,0,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode));
		}
		else
		{
			errval = rsb__do_permute_values_with_coo_index(rVA, VA, rIA, IA, rJA, JA, CP, nnz, typecode);
			RSB_COO_MEMCPY(VA,IA,JA,rVA,rIA,rJA,0,0,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode));
		}
	}
	else
	{
		if(flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT)
			errval = rsb__do_permute_values_in_place_with_nnz_index(rVA, rIA, rJA, KP, nnz, typecode);
		else
			errval = rsb__do_permute_values_with_nnz_index(rVA, VA, rIA, IA, rJA, JA, KP, nnz, typecode);
	}
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(want_two_pass_sort)
	{
		/* second pass: by row */
		rsb_nnz_idx_t n1=0,n2=0;
		const size_t el_size = RSB_SIZEOF(typecode);

		if(!CP)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
		while(n1!=nnz)
		{
#if 0
			while( n2+1<nnz && IA[n1]/br==IA[n2+1]/br )
				++n2;
#else
			/* EXPERIMENTAL: determine end of row (TODO: shall better use function here) */
			/*
				we don't know in advance how many elements belong to this block row.
				we first go fast forward, then slow down :)
			 */
			rsb_nnz_idx_t delta=1;

			while( n2+delta<nnz && IA[n1]/br==IA[n2+delta]/br )
				n2+=delta, delta*=2;

			/* now, n2+delta>=nnz  ||  IA[n1]/br!=IA[n2+delta]/br */
	                RSB_DEBUG_ASSERT(n2+delta>=nnz  ||  IA[n1]/br!=IA[n2+delta]/br);
			/* if delta == 0, we are done. */
			while( delta>0 )
			{
				if( n2>=nnz || IA[n1]/br!=IA[n2]/br )
					delta/=2, n2-=delta;
				else
				if( n2+delta<nnz && IA[n1]/br==IA[n2+delta]/br )
					n2+=delta, delta/=2;
				else
					delta/=2;
			}
			RSB_DEBUG_ASSERT( n2  < nnz && IA[n1]/br==IA[n2  ]/br );
			RSB_DEBUG_ASSERT( n2+1>=nnz || IA[n1]/br!=IA[n2+1]/br );
	                RSB_DEBUG_ASSERT(n2<nnz);
	                RSB_DEBUG_ASSERT(IA[n1]/br==IA[n2]/br);
#endif
			/* NOTE: need a specialized code, once rsb_nnz_idx_t != rsb_coo_idx_t */
	                RSB_DEBUG_ASSERT(sizeof(rsb_nnz_idx_t) == sizeof(rsb_coo_idx_t));

			if(0)RSB_INFO("sorting : %zd .. %zd\n",(rsb_printf_int_t)n1,(rsb_printf_int_t)n2);
			//consider alternative:
			//errval = rsb__do_coo_index_sort_on_rows_array_make(CP,JA+n1,nc,bc,(n2+1)-n1,typecode);
			errval = rsb__do_nnz_index_sort_array_make(CP,IA+n1,JA+n1,nr,nc,IA[n1],br,bc,(n2+1)-n1,typecode,flags,0,op_flags);
			if(RSB_SOME_ERROR(errval))
			{
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
			}
			qsort( CP , (size_t) (n2+1)-n1, 2*sizeof(rsb_coo_idx_t), &rsb_compar_coo_idx_t );
			rsb__do_util_compact_permutation_coo_idx_t_array(CP, (n2+1)-n1);
			if(flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT)
				errval = rsb__do_permute_values_in_place_with_coo_index(RSB_VA_OFFSET_POINTER(rVA,el_size,n1), rIA+n1, rJA+n1, CP, (n2+1)-n1, typecode);
			else
				errval = rsb__do_permute_values_with_coo_index( RSB_VA_OFFSET_POINTER(rVA,el_size,n1), RSB_VA_OFFSET_POINTER(VA,el_size,n1), rIA+n1, IA+n1, rJA+n1, JA+n1, CP, (n2+1)-n1, typecode);

			if(RSB_SOME_ERROR(errval))
			{
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
			}

			++n2;
			n1=n2;
		}
	}
	pt += rsb_time();
	RSB_DEBUG_ASSERT(pt >= RSB_TIME_ZERO);
	RSB_DEBUG_ASSERT(st >= RSB_TIME_ZERO);
	RSB_DEBUG_ASSERT(it >= RSB_TIME_ZERO);

//	RSB_INFO("#sorting : nnz/s : %lg\n",((double)nnz)/(st));
//	RSB_INFO("#sorting : nnz*log(nnz)/s : %lg\n",((double)nnz)*log((double)nnz)/(st));

	// sorting/permutation time = 4 ~ 10 
	if( RSB_WANT_VERBOSE_MESSAGES )
	RSB_INFO(	"# sorting times:\n"
			"#index init 		: %lgs\n"
			"#index sorting (qsort)	: %lgs\n"
			"#data permutation	: %lgs\n"
			"#sorting/permutation	: %lgs\n",
			it,st,pt,st/pt
		);
	//RSB_INFO("#allocation time   : %lg\n",at);

err:
	if(CP && want_two_pass_sort)
	{
		if( wb >= RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_TWO_PASS(nnz,nr,nc,br,bc))
			;
		else
			#pragma omp critical (rsb_coo_sort)
			RSB_CONDITIONAL_FREE(CP);
	} 
	if(KP  && !want_two_pass_sort)
	{
		if( wb >= RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT_ONE_PASS(nnz,nr,nc,br,bc))
			;
		else
			#pragma omp critical (rsb_coo_sort)
			RSB_CONDITIONAL_FREE(KP);
	} 

	RSB_DO_ERR_RETURN(errval)
} /* rsb__do_index_based_bcsr_sort */

static RSB_INLINE rsb_coo_idx_t rsb_do_undilate_coo_odd(rsb_coo_idx_t w)
{
	rsb_coo_idx_t E=w;
	if (sizeof(rsb_coo_idx_t)==1)
	{
		E = (E & 0x11) | ((E & 0x44)>>1);
		E = (E & 0x0F) | ((E & 0x30)>>2);
	}
	else
	if (sizeof(rsb_coo_idx_t)==2)
	{
		E = (E & 0x1111) | ((E & 0x4444)>>1);
		E = (E & 0x0303) | ((E & 0x3030)>>2);
		E = (E & 0x000F) | ((E & 0x0F00)>>4);
	}
	else
	if (sizeof(rsb_coo_idx_t)==4)
	{
		E = (E & 0x11111111) | ((E & 0x44444444)>>1);
		E = (E & 0x03030303) | ((E & 0x30303030)>>2);
		E = (E & 0x000F000F) | ((E & 0x0F000F00)>>4);
		E = (E & 0x000000FF) | ((E & 0x00FF0000)>>8);
	}
	else
	if (sizeof(rsb_coo_idx_t)==8)
	{
		E = (E & 0x1111111111111111) | ((E & 0x4444444444444444)>>1 );
		E = (E & 0x0303030303030303) | ((E & 0x3030303030303030)>>2 );
		E = (E & 0x000F000F000F000F) | ((E & 0x0F000F000F000F00)>>4 );
		E = (E & 0x000000FF000000FF) | ((E & 0x00FF000000FF0000)>>8 );
		E = (E & 0x000000000000FFFF) | ((E & 0x000000FFFF000000)>>16);
		RSB_ERROR(RSB_ERRM_FYCITINS);
	}
	else
	{
		RSB_ERROR(RSB_ERRM_FYCITINS);
		/* FIXME : fatal! */
	}
	return E;
}

static RSB_INLINE rsb_coo_idx_t rsb_do_undilate_coo_even(rsb_coo_idx_t w)
{
	return rsb_do_undilate_coo_odd(w>>1);
}

static RSB_INLINE rsb_coo_idx_t rsb_do_dilate_coo(rsb_coo_idx_t w)
{
	rsb_coo_idx_t E=w;
	if (sizeof(rsb_coo_idx_t)==1)
	{
		E = (E | (E << 2)) & 0x33;
		E = (E | (E << 1)) & 0x55;
	}
	else
	if (sizeof(rsb_coo_idx_t)==2)
	{
		E = (E | (E << 4)) & 0x0F0F;
		E = (E | (E << 2)) & 0x3333;
		E = (E | (E << 1)) & 0x5555;
	}
	else
	if (sizeof(rsb_coo_idx_t)==4)
	{
		E = (E | (E << 8)) & 0x00FF00FF;
		E = (E | (E << 4)) & 0x0F0F0F0F;
		E = (E | (E << 2)) & 0x33333333;
		E = (E | (E << 1)) & 0x55555555;
	}
	else
	if (sizeof(rsb_coo_idx_t)==8)
	{
		E = (E | (E <<16)) & 0x0000FFFF0000FFFF;
		E = (E | (E << 8)) & 0x00FF00FF00FF00FF;
		E = (E | (E << 4)) & 0x0F0F0F0F0F0F0F0F;
		E = (E | (E << 2)) & 0x3333333333333333;
		E = (E | (E << 1)) & 0x5555555555555555;
	}
	else
	{
		RSB_ERROR(RSB_ERRM_FYCITINS);
		/* FIXME : fatal! */
	}
	return E;
}

#define RSB_COO_INDEX_LO_MASK	0x0000FFFF
#define RSB_COO_INDEX_HI_MASK	0xFFFF0000
#define RSB_COO_INDEX_HBITSOF    ((RSB_CHAR_BIT*(sizeof(rsb_coo_idx_t))/2))
#define RSB_COO_INDEX_HI_SHIFTED(X) (((X)&RSB_COO_INDEX_HI_MASK)>>RSB_COO_INDEX_HBITSOF)
#define RSB_COO_INDEX_EVEN_MASK	0xAAAAAAAA
#define RSB_COO_INDEX_ODD_MASK	0x55555555

static RSB_INLINE rsb_coo_idx_t RSB_Z_2_COO_HI_WORD(rsb_coo_idx_t i, rsb_coo_idx_t j)
{
	return rsb_do_undilate_coo_even(j) |(rsb_do_undilate_coo_even(i)<<RSB_COO_INDEX_HBITSOF);
}
static RSB_INLINE rsb_coo_idx_t RSB_Z_2_COO_LO_WORD(rsb_coo_idx_t i, rsb_coo_idx_t j)
{
	return rsb_do_undilate_coo_odd(j) |(rsb_do_undilate_coo_odd(i)<<RSB_COO_INDEX_HBITSOF);
}
static RSB_INLINE rsb_coo_idx_t RSB_COO_2_Z_HI_WORD(rsb_coo_idx_t i, rsb_coo_idx_t j)
{
	return ( (rsb_do_dilate_coo(RSB_COO_INDEX_HI_SHIFTED(i))<<1)| (rsb_do_dilate_coo(RSB_COO_INDEX_HI_SHIFTED(j))));
}
static RSB_INLINE rsb_coo_idx_t RSB_COO_2_Z_LO_WORD(rsb_coo_idx_t i, rsb_coo_idx_t j)
{
	return ( (rsb_do_dilate_coo(i&RSB_COO_INDEX_LO_MASK)<<1)| (rsb_do_dilate_coo(j&RSB_COO_INDEX_LO_MASK) ));
}	

rsb_err_t rsb__do_index_based_z_morton_sort( 
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA,
	rsb_coo_idx_t * rIA, rsb_coo_idx_t * rJA, void * rVA,
	rsb_coo_idx_t nr, rsb_coo_idx_t nc,
	rsb_nnz_idx_t nnz,
	rsb_type_t typecode
	,enum rsb_op_flags_t op_flags
	)
{
	/**
		\ingroup gr_internals
	*/
	//rsb_nnz_idx_t * K=NULL;
	//rsb_coo_idx_t Idim;
	//rsb_coo_idx_t Jdim;
	//const rsb_coo_idx_t br=1,bc=1;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	
	//rsb_nnz_idx_t n;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//rsb_time_t it=0;
	//rsb_nnz_idx_t nIdim,nJdim;

	RSB_DEBUG_ASSERT(nr);
	RSB_DEBUG_ASSERT(nc);
	//RSB_DEBUG_ASSERT(br);
	//RSB_DEBUG_ASSERT(bc);

	{
		rsb_nnz_idx_t n=0;
		rsb_coo_idx_t t;
		for(n=0;n<nnz;++n)
		{
//			printf("%d: %0x %0x -> ",n,rIA[n],rJA[n]);
			t = RSB_COO_2_Z_HI_WORD(rIA[n],rJA[n]);
			rJA[n]=RSB_COO_2_Z_LO_WORD(rIA[n],rJA[n]);
			rIA[n]=t;
//			printf(" %0x %0x\n",rIA[n],rJA[n]) ;
		}
	//	for(n=0;n<nnz;++n) printf("%d: %d %d\n",n,rIA[n],rJA[n]);
		errval = rsb__util_sort_row_major_inner(rVA,rIA,rJA,nnz,nr,nc,typecode,flags);
		for(n=0;n<nnz;++n)
		{
//			printf("%x: %0x %0x -> ",n,rIA[n],rJA[n]);
			t = RSB_Z_2_COO_HI_WORD(rIA[n],rJA[n]);
			rJA[n]=RSB_Z_2_COO_LO_WORD(rIA[n],rJA[n]);
			rIA[n]=t;
//			printf(" %0x %0x\n",rIA[n],rJA[n]) ;
		}
//		for(n=0;n<nnz;++n) printf("%d: %d %d\n",n,rIA[n],rJA[n]);
	}
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
