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
 * @author Michele Martone
 * @brief Sorting routines.
 */
#include "rsb_msort_up.h"
#include "rsb_internals.h"

static rsb_nnz_idx_t iabs(rsb_nnz_idx_t i)
{
	/**
	 * \ingroup gr_internals
	 * Fortran style iabs().
  	 */
	return i<0?-i:i;
}

static rsb_nnz_idx_t isign(rsb_nnz_idx_t i,rsb_nnz_idx_t j)
{
	/**
	 * \ingroup gr_internals
	 * Fortran style lsign():
	 * "Returns the absolute value of A times the sign of B."
  	 */
	register rsb_nnz_idx_t ia=i<0?-i:i;
	return j>=0?ia:-ia;
}

rsb_err_t rsb__do_msort_up(rsb_nnz_idx_t n, const rsb_nnz_idx_t * RSB_RESTRICT k, rsb_nnz_idx_t * RSB_RESTRICT l)
{
	/**
	 	\ingroup gr_internals
		Adapted C code from PSBLAS Fortran msort_up routine.
		\param n
		\param k an n   sized array, for input indices
		\param l an n+2 sized array, for output permutation links
	*/
	/* integer k(n),l(0:n+1) */
	rsb_nnz_idx_t p,q,s,t;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* intrinsic iabs,isign */
	/*
		first step: we are preparing ordered sublists, exploiting
		what order was already in the input data; negative links
	  	mark the end of the sublists
	*/

	--k;	/* pointer fix for C */

	l[0] = 1;
	t = n + 1;
	for(p = 1; RSB_LIKELY(p<=n - 1);++p)
	{
		if (k[p] <= k[p+1])
		{
			l[p] = p + 1;
		}
		else
		{
			l[t] = - (p+1);
			t = p;
		}
	}
	l[t] = 0;
	l[n] = 0;
	// see if the input was already sorted
	if (l[n+1] == 0)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	else
	{
		l[n+1] = iabs(l[n+1]);
	}

mergepass:
	/*
		 otherwise, begin a pass through the list.
		 throughout all the subroutine we have:
		  p, q: pointing to the sublists being merged
		  s: pointing to the most recently processed record
		  t: pointing to the end of previously completed sublist
	*/
	s = 0;
	t = n + 1;
	p = l[s];
	q = l[t];
	if (RSB_UNLIKELY(q == 0)) goto mergepass_exit;

outer:

	if (k[p] > k[q])
	{
		l[s] = isign(q,l[s]);
		s = q;
		q = l[q];
		if (q > 0)
		{
			while(1)
			{
				if (k[p]<= k[q]) goto outer;
				s = q;
				q = l[q];
				if (q <= 0) break;
			}
		}
		l[s] = p;
		s = t;
		while(1)
		{
			t = p;
			p = l[p];
			if (p <= 0) break;
		}
	}
	else
	{
		l[s] = isign(p,l[s]);
		s = p;
		p = l[p];
		if (p>0)
		{
			while(1)
			{
				if (k[p] > k[q]) goto outer;
				s = p;
				p = l[p];
				if (p <= 0) break;
			}
		}
		//  otherwise, one sublist ended, and we append to it the rest
		// of the other one.
		l[s] = q;
		s = t;
		while(1)
		{
			t = q;
			q = l[q];
			if (q <= 0) break;
		}
	}

	p = -p;
	q = -q;
	if (q == 0)
	{

		l[s] = isign(p,l[s]);
		l[t] = 0;
		goto outer_out;
	}

	goto outer;
outer_out:

	goto mergepass;
mergepass_exit:
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_msort_up2coo(rsb_nnz_idx_t n, const rsb_coo_idx_t * k, rsb_nnz_idx_t * l)
{
	/**
	 	\ingroup gr_internals
		Adapted C code from PSBLAS Fortran msort_up routine.
		Modified to handle lexicographical order.
		The only difference with rsb__do_msort_up is the comparison on k.

		\param n
		\param k an 2*n   sized array, for input indices
		\param l an n+2 sized array, for output permutation links
	*/
	/* integer k(n),l(0:n+1) */
	rsb_nnz_idx_t p,q,s,t;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* intrinsic iabs,isign */
	/*
		first step: we are preparing ordered sublists, exploiting
		what order was already in the input data; negative links
	  	mark the end of the sublists
	*/

	--k;	/* pointer fix for C */
	--k;	/* pointer fix for C */

	l[0] = 1;
	t = n + 1;
	for(p = 1; p<=n - 1;++p)
	{
/* lexicographical order, less than or equal */
#define RSB_LO_LTOE(H1,L1,H2,L2)	(H1<H2 || (H1==H2 && L1<=L2))
/* lexicographical order, greater than */
#define RSB_LO_GT(H1,L1,H2,L2)		(H1>H2 || (H1==H2 && L1> L2))

		if( RSB_LO_LTOE(k[2*p],k[2*p+1],k[2*p+2],k[2*p+2+1]) )
		{
			l[p] = p + 1;
		}
		else
		{
			l[t] = - (p+1);
			t = p;
		}
	}
	l[t] = 0;
	l[n] = 0;
	// see if the input was already sorted
	if (l[n+1] == 0)
	{
		errval = RSB_ERR_BADARGS;
		goto err ;
	}
	else
	{
		l[n+1] = iabs(l[n+1]);
	}

mergepass:
	/*
		 otherwise, begin a pass through the list.
		 throughout all the subroutine we have:
		  p, q: pointing to the sublists being merged
		  s: pointing to the most recently processed record
		  t: pointing to the end of previously completed sublist
	*/
	s = 0;
	t = n + 1;
	p = l[s];
	q = l[t];
	if (RSB_UNLIKELY(q == 0)) goto mergepass_exit;

outer:

	if( RSB_LO_GT(k[2*p],k[2*p+1],k[2*q],k[2*q+1]) )
	{
		l[s] = isign(q,l[s]);
		s = q;
		q = l[q];
		if (q > 0)
		{
			while(1)
			{
				if( RSB_LO_LTOE(k[2*p],k[2*p+1],k[2*q],k[2*q+1]) ) goto outer;
				s = q;
				q = l[q];
				if (q <= 0) break;
			}
		}
		l[s] = p;
		s = t;
		while(1)
		{
			t = p;
			p = l[p];
			if (p <= 0) break;
		}
	}
	else
	{
		l[s] = isign(p,l[s]);
		s = p;
		p = l[p];
		if (p>0)
		{
			while(1)
			{
				if( RSB_LO_GT(k[2*p],k[2*p+1],k[2*q],k[2*q+1]) ) goto outer;
				s = p;
				p = l[p];
				if (p <= 0) break;
			}
		}
		//  otherwise, one sublist ended, and we append to it the rest
		// of the other one.
		l[s] = q;
		s = t;
		while(1)
		{
			t = q;
			q = l[q];
			if (q <= 0) break;
		}
	}

	p = -p;
	q = -q;
	if (q == 0)
	{

		l[s] = isign(p,l[s]);
		l[t] = 0;
		goto outer_out;
	}

	goto outer;
outer_out:

	goto mergepass;
mergepass_exit:
err:
	RSB_DO_ERR_RETURN(errval)
}

/* @endcond */
