/*                                                                                                                            

Copyright (C) 2008-2020 Michele Martone

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
 * This source file contains searching functions.
 * */
#include "rsb_internals.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

rsb_nnz_idx_t rsb__nnz_split_hcoo_bsearch(const rsb_half_idx_t *A, const rsb_half_idx_t S,const rsb_nnz_idx_t n)
{
	/*!
	 * \ingroup gr_internals
	 *
 	 * \return the found index, or 0
	 *
	 * Performs a binary search in the given sorted array to find the first
  	 * element which is >= S, and returns its index.
	 * \note : n=0 is not allowed
	*/
	register rsb_nnz_idx_t l=0,h=n-1,mid;

	if(n<1)
		return 0;

	if( S > A[h]  )
		return n;/* no such element */

	if( A[l]>=S )/* the point we look for could be before */
		return l;

	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(S));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));
	RSB_DEBUG_ASSERT(n>1);

	do
	{
		mid=l + ((h+1)-l)/2;

		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(A[mid]));
		RSB_DEBUG_ASSERT( A[mid] >= A[0] );
		RSB_DEBUG_ASSERT( A[mid] <= A[n-1] );

//		RSB_INFO("hop %d, at %d\n",hop,base+hop);
//		RSB_INFO("h %d l %d m %d\n",h,l,mid);
		if( A[mid]<S )/* the point we search is after this */
		{
			l=mid;
//		RSB_INFO("*+\n");
		}
		else
		if( A[mid]>=S )/* the point we look for could be before */
		{
			if(h==mid)
				goto ok;
			h=mid;
//			RSB_INFO("*-\n");
		}
	}
	while(RSB_LIKELY(l!=h));
	
	ok:

//	RSB_INFO(" \n");

	return mid;
}

rsb_nnz_idx_t rsb__nnz_split_nnz_bsearch(const rsb_nnz_idx_t*A,const rsb_nnz_idx_t S,const rsb_nnz_idx_t n)
{
	/*!
	 * \ingroup gr_internals
	 *
 	 * \return the found index, or 0
	 *
	 * Performs a binary search in the given sorted array to find the first
  	 * element which is >= S, and returns its index.
	 * \note : n<=0 is not allowed
	*/
	register rsb_nnz_idx_t l=0,h=n-1,mid;

	if(n<1)
		return 0;

	if( S > A[h]  )
		return n;/* no such element */

	if( A[l]>=S )/* the point we look for could be before */
		return l;

	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(S));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));
	RSB_DEBUG_ASSERT(n>0);

	do
	{
		mid=l + ((h+1)-l)/2;

		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(A[mid]));
		RSB_DEBUG_ASSERT( A[mid] >= A[0] );
		RSB_DEBUG_ASSERT( A[mid] <= A[n-1] );

//		RSB_INFO("hop %d, at %d\n",hop,base+hop);
//		RSB_INFO("h %d l %d m %d\n",h,l,mid);
		if( A[mid]<S )/* the point we search is after this */
		{
			l=mid;
//		RSB_INFO("*+\n");
		}
		else
		if( A[mid]>=S )/* the point we look for could be before */
		{
			if(h==mid)
				goto ok;
			h=mid;
//			RSB_INFO("*-\n");
		}
	}
	while(RSB_LIKELY(l!=h));
	
	ok:

//	RSB_INFO(" \n");

	return mid;
}

rsb_nnz_idx_t rsb__nnz_split_coo_bsearch(const rsb_coo_idx_t*A,const rsb_coo_idx_t S,const rsb_nnz_idx_t n)
{
	/*!
	 * \ingroup gr_internals
	 *
 	 * \return the found index, or 0
	 *
	 * Performs a binary search in the given sorted array to find the first
  	 * element which is >= S, and returns its index.
	 * \note : n<=0 is not allowed
	*/
	register rsb_nnz_idx_t l=0,h=n-1,mid;

	if(n<1)
		return 0;

	if( S > A[h]  )
		return n;/* no such element */

	if( A[l]>=S )/* the point we look for could be before */
		return l;

	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(S));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));
	RSB_DEBUG_ASSERT(n>0);

	do
	{
		mid=l + ((h+1)-l)/2;

		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(A[mid]));
		RSB_DEBUG_ASSERT( A[mid] >= A[0] );
		RSB_DEBUG_ASSERT( A[mid] <= A[n-1] );

//		RSB_INFO("hop %d, at %d\n",hop,base+hop);
//		RSB_INFO("h %d l %d m %d\n",h,l,mid);
		if( A[mid]<S )/* the point we search is after this */
		{
			l=mid;
//		RSB_INFO("*+\n");
		}
		else
		if( A[mid]>=S )/* the point we look for could be before */
		{
			if(h==mid)
				goto ok;
			h=mid;
//			RSB_INFO("*-\n");
		}
	}
	while(RSB_LIKELY(l!=h));
	
	ok:

//	RSB_INFO(" \n");

	return mid;
}

rsb_nnz_idx_t rsb__seek_coo_idx_t(const rsb_coo_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n)
{
	/**
	 	\ingroup gr_internals

		Assuming p an array monotonically, increasingly sorted, will return the first index in it containing 
		a value equal to v among the first n.
		Returns RSB_MARKER_NNZ_VALUE if no such value exists.
		Ready to handle unsigned indices.
	*/
	rsb_nnz_idx_t d=n;
	rsb_nnz_idx_t e=1;

	RSB_DEBUG_ASSERT(v>=0);
	RSB_DEBUG_ASSERT(n>=0);
	RSB_DEBUG_ASSERT(p);

	if(n<1)
		goto head_not_found;
	if(n==1)
		{d=0;goto maybe;}
	/* n>1 */

	while(e<((n-1)/2)+1)
		e*=2;
	d=e;
	e/=2;

	/* now e < n */
	RSB_DEBUG_ASSERT(e<n);
	RSB_DEBUG_ASSERT(d>0);
	RSB_DEBUG_ASSERT(d>=e);

	do
	{
		if(v>=p[d])
		{
			if(v==p[d])
				goto head_found;
			else
				if(d+e<n)
					d+=e;
		}
		else
			d-=e;
		e/=2;
	}
	while(RSB_LIKELY(e>0));
maybe:
	if(v==p[d])
		goto head_found;

	if(d==1 && p[d=0]==v)
		goto head_found;

head_not_found:
		return RSB_MARKER_NNZ_VALUE;
head_found:
		return d;
}

rsb_nnz_idx_t rsb__seek_half_idx_t(const rsb_half_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n)
{
	/**
	 	\ingroup gr_internals

		Assuming p an array monotonically, increasingly sorted, will return the first index in it containing 
		a value equal to v among the first n.
		Returns RSB_MARKER_NNZ_VALUE if no such value exists.
		Ready to handle unsigned indices.
	*/
	rsb_nnz_idx_t d=n;
	rsb_nnz_idx_t e=1;

	RSB_DEBUG_ASSERT(v>=0);
	RSB_DEBUG_ASSERT(n>=0);
	RSB_DEBUG_ASSERT(p);

	if(n<1)
		goto head_not_found;
	if(n==1)
		{d=0;goto maybe;}
	/* n>1 */

	while(e<((n-1)/2)+1)
		e*=2;
	d=e;
	e/=2;

	/* now e < n */
	RSB_DEBUG_ASSERT(e<n);
	RSB_DEBUG_ASSERT(d>0);
	RSB_DEBUG_ASSERT(d>=e);

	do
	{
		if(v>=p[d])
		{
			if(v==p[d])
				goto head_found;
			else
				if(d+e<n)
					d+=e;
		}
		else
			d-=e;
		e/=2;
	}
	while(RSB_LIKELY(e>0));
maybe:
	if(v==p[d])
		goto head_found;

	if(d==1 && p[d=0]==v)
		goto head_found;

head_not_found:
		return RSB_MARKER_NNZ_VALUE;
head_found:
		return d;
}

rsb_nnz_idx_t rsb__seek_nnz_idx_t(const rsb_nnz_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n)
{
	/**
	 	\ingroup gr_internals

		Assuming p an array monotonically, increasingly sorted, will return the first index in it containing 
		a value equal to v among the first n.
		Returns RSB_MARKER_NNZ_VALUE if no such value exists.
		Ready to handle unsigned indices.
	*/
#if 1
	rsb_nnz_idx_t d=n;
	rsb_nnz_idx_t e=1;

	RSB_DEBUG_ASSERT(v>=0);
	RSB_DEBUG_ASSERT(n>=0);
	RSB_DEBUG_ASSERT(p);

	if(n<1)
		goto head_not_found;
	if(n==1)
		{d=0;goto maybe;}
	/* n>1 */

	while(e<((n-1)/2)+1)
		e*=2;
	d=e;
	e/=2;

	/* now e < n */
	RSB_DEBUG_ASSERT(e<n);
	RSB_DEBUG_ASSERT(d>0);
	RSB_DEBUG_ASSERT(d>=e);

	do
	{
		if(v>=p[d])
		{
			if(v==p[d])
				goto head_found;
			else
				if(d+e<n)
					d+=e;
		}
		else
			d-=e;
		e/=2;
	}
	while(RSB_LIKELY(e>0));
maybe:
	if(v==p[d])
		goto head_found;

	if(d==1 && p[d=0]==v)
		goto head_found;

head_not_found:
		return RSB_MARKER_NNZ_VALUE;
head_found:
		return d;
#else
	RSB_DEBUG_ASSERT(v>=0);
	RSB_DEBUG_ASSERT(n>=0);
	RSB_DEBUG_ASSERT(p);
	rsb_nnz_idx_t k=0;
	/* fallback, slow */
	for(k=0;RSB_LIKELY(k<n);++k)
		if(p[k]==v)
			return k;
#endif
	return RSB_MARKER_NNZ_VALUE;
}

rsb_nnz_idx_t rsb__seek_nnz_idx_t_linear(const rsb_nnz_idx_t *p, rsb_nnz_idx_t v, rsb_nnz_idx_t n)
{
	/**
	 	\ingroup gr_internals
	*/
	rsb_nnz_idx_t k=0;
	RSB_DEBUG_ASSERT(v>=0);
	RSB_DEBUG_ASSERT(n>=0);
	RSB_DEBUG_ASSERT(p);
	/* fallback, slow */
	for(k=0;RSB_LIKELY(k<n);++k)
		if(p[k]==v)
			return k;
	return RSB_MARKER_NNZ_VALUE;
}

/* @endcond */

