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

#define RSB_LIKELY_OMP(EXP) EXP

RSB_INTERNALS_COMMON_HEAD_DECLS

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_nnz_array_set_sequence(rsb_nnz_idx_t * p, rsb_nnz_idx_t n, rsb_nnz_idx_t o, rsb_nnz_idx_t i)
{
	/*!
		\ingroup gr_internals
		TODO: rsb__util_nnz_array_set_sequence -> rsb__nnz_idx_set_sequence
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(o>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(o));

	for(k=0;RSB_LIKELY(k<n);++k)
	{
		p[k] = o+k*i;
	}
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void rsb__util_coo_array_set_sequence(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t o, rsb_coo_idx_t i)
{
	/*!
		\ingroup gr_internals
		TODO: rsb__util_coo_array_set_sequence -> rsb__coo_idx_set_sequence
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(o>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(o));

	for(k=0;RSB_LIKELY(k<n);++k)
	{
		p[k] = o+k*i;
	}
}

void rsb__util_nnz_array_set(rsb_nnz_idx_t * p, rsb_nnz_idx_t n, rsb_nnz_idx_t a)
{
	/*!
		\ingroup gr_internals
		Adds s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(a>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(a));

	for(k=0;RSB_LIKELY(k<n);++k)
	{
		p[k] = a;
	}
}

void rsb__util_coo_array_set(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a)
{
	/*!
		\ingroup gr_internals
		Adds s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(a>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(a) || a==RSB_MARKER_COO_VALUE);

	for(k=0;RSB_LIKELY(k<n);++k)
	{
		p[k] = a;
	}
}

void rsb__util_hcoo_array_copy_trans_add(rsb_coo_idx_t * d, const rsb_half_idx_t * s, rsb_nnz_idx_t n, rsb_coo_idx_t a)
{
	/*!
		\ingroup gr_internals
		Adds s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(s);
	RSB_DEBUG_ASSERT(d);
	RSB_DEBUG_ASSERT(a>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	for(k=0;RSB_LIKELY(k+3<n);k+=4)
	{
		/* RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(d[k])); */
		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(a+s[k]));
		d[k+0] = a+s[k+0];
		d[k+1] = a+s[k+1];
		d[k+2] = a+s[k+2];
		d[k+3] = a+s[k+3];
		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(d[k]));
	}
	for(;(k<n);k+=1)
		d[k] = s[k]+a;
}

void rsb__util_coo_array_copy_trans_add(rsb_coo_idx_t * d, const rsb_coo_idx_t * s, rsb_nnz_idx_t n, rsb_coo_idx_t a)
{
	/*!
		\ingroup gr_internals
		Adds s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(s);
	RSB_DEBUG_ASSERT(d);
	RSB_DEBUG_ASSERT(a>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	for(k=0;RSB_LIKELY(k+3<n);k+=4)
	{
		/* RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(d[k])); */
		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(a+s[k]));
		d[k+0] = a+s[k+0];
		d[k+1] = a+s[k+1];
		d[k+2] = a+s[k+2];
		d[k+3] = a+s[k+3];
		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(d[k]));
	}
	for(;(k<n);k+=1)
		d[k] = s[k]+a;
}

static void rsb__util_coo_array_mul(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a)
{
	/*!
		\ingroup gr_internals
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(a>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(a)
		for(k=0;RSB_LIKELY(k<n);++k)
		{
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
			p[k] *= a;
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
		}
}

void rsb__util_coo_arrays_mul(rsb_coo_idx_t * p, rsb_coo_idx_t * q, rsb_coo_idx_t a, rsb_coo_idx_t b, rsb_nnz_idx_t n)
{
	/*!
		multiply COO arrays p[n] and q[n] respectively by a and b.
	*/
	rsb__util_coo_array_mul(p,n,a);
	rsb__util_coo_array_mul(q,n,b);
}

void rsb__util_coo_array_add(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a)
{
	/*!
		\ingroup gr_internals
		Adds s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	/* RSB_DEBUG_ASSERT(a>=0); */
	RSB_DEBUG_ASSERT(a>=-1); /* FIXME: -1 is sometimes necessary for Fortran cases ... */
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(a)
		for(k=0;RSB_LIKELY(k<n);++k)
		{
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
			p[k] += a;
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
		}
}

void rsb__util_hcoo_array_add(rsb_half_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t a)
{
	/*!
		\ingroup gr_internals
		Adds s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(a>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(a)
		for(k=0;RSB_LIKELY(k<n);++k)
		{
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
			p[k] += a;
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
		}
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_coo_arrays_add(rsb_coo_idx_t * p, rsb_coo_idx_t * q, rsb_coo_idx_t a, rsb_coo_idx_t b, rsb_nnz_idx_t n)
{
	/*!
		TODO : document
	*/
	rsb__util_coo_array_add(p,n,a);
	rsb__util_coo_array_add(q,n,b);
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void rsb_util_coo_arrays_sub(rsb_coo_idx_t * p, rsb_coo_idx_t * q, rsb_coo_idx_t a, rsb_coo_idx_t b, rsb_nnz_idx_t n)
{
	/*!
		TODO : document
	*/
	rsb__util_coo_array_sub(p,n,a);
	rsb__util_coo_array_sub(q,n,b);
}

void rsb__util_nnz_array_add(rsb_nnz_idx_t * p, rsb_nnz_idx_t n, rsb_nnz_idx_t a)
{
	/*!
		\ingroup gr_internals
		Subtracts s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(a>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(a)
		for(k=0;RSB_LIKELY(k<n);++k)
		{
			p[k] += a;
			RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(p[k]));
		}
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
static void rsb__util_nnz_array_sub(rsb_nnz_idx_t * p, rsb_nnz_idx_t n, rsb_nnz_idx_t s)
{
	/*!
		\ingroup gr_internals
		Subtracts s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(s>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(s)
		for(k=0;RSB_LIKELY(k<n);++k)
		{
			p[k] -= s;
			RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(p[k]));
		}
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void rsb__util_coo_array_sub(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_coo_idx_t s)
{
	/*!
		\ingroup gr_internals
		Subtracts s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(s>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(s)
		for(k=0;RSB_LIKELY(k<n);++k)
		{
			p[k] -= s;
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
		}
}

void rsb__util_nnz_array_to_fortran_indices(rsb_nnz_idx_t * p, rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
		Adds 1 to the given matrix coordinate indices array.
	*/
	rsb__util_nnz_array_add(p, n, 1);
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_coo_array_to_fortran_indices(rsb_coo_idx_t * p, rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
		Adds 1 to the given matrix coordinate indices array.
	*/
	rsb__util_coo_array_add(p, n, 1);
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void rsb__util_coo_array_to_fortran_indices_parallel(rsb_coo_idx_t * p, rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
		Adds 1 to the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	const rsb_nnz_idx_t mcs = RSB_MINIMUM_VECOP_OMP_CHUNK; 
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(n>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(p)
	{
		#pragma omp parallel for schedule(guided,mcs) shared(p)   RSB_NTC
		for(k=0;RSB_LIKELY_OMP(k<n);++k)
		{
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
			++p[k];
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
		}
	}
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_nnz_array_from_fortran_indices(rsb_nnz_idx_t * p, rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
		Subtracts 1 from the given matrix coordinate indices array.
	*/
	rsb__util_nnz_array_sub(p, n, 1);
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if RSB_OBSOLETE_QUARANTINE_UNUSED
// OK, but unused 
void rsb__util_coo_array_from_fortran_indices(rsb_coo_idx_t * p, rsb_nnz_idx_t n, rsb_bool_t want_parallel)
{
	/*!
		\ingroup gr_internals
		Subtracts 1 from the given matrix coordinate indices array.
		TODO: parallelize
	*/
	register rsb_nnz_idx_t k;
	const rsb_nnz_idx_t mcs = RSB_MINIMUM_VECOP_OMP_CHUNK; 

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(p>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	if(!p)
		return;

	if(want_parallel)
	{
		#pragma omp parallel for schedule(guided,mcs) shared(p)   RSB_NTC
		for(k=0;RSB_LIKELY_OMP(k<n);++k)
		{
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
			--p[k];
			RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(p[k]));
		}
	}
	else
		rsb__util_coo_array_sub(p, n, 1);
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_flags_t rsb__util_coo_determine_uplo_flags(const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t nnz)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t n;
	const rsb_flags_t tflags = RSB_FLAG_UPPER_TRIANGULAR|RSB_FLAG_LOWER_TRIANGULAR;
	rsb_flags_t flags = tflags;
	
	for(n=0;n<nnz ;++n)
	if(IA[n]<JA[n])
	{
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_LOWER);
		for(;RSB_LIKELY(n<nnz) ;++n)
		if(IA[n]>JA[n])
		{
			RSB_DO_FLAG_DEL(flags,RSB_FLAG_UPPER);
			goto done;
		}
	}
	else
	if(IA[n]>JA[n])
	{
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_UPPER);
		for(;RSB_LIKELY(n<nnz) ;++n)
		if(IA[n]<JA[n])
		{
			RSB_DO_FLAG_DEL(flags,RSB_FLAG_LOWER);
			goto done;
		}
	}
done:
	if((flags&tflags)==RSB_FLAG_TRIANGULAR )
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_TRIANGULAR);
	return flags;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_bool_t rsb__util_coo_check_if_triangle_non_empty(const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t n;
	
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER))
		for(n=0;RSB_LIKELY(n<nnz);++n)
			if(IA[n]<JA[n])
				return RSB_BOOL_TRUE;

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER))
		for(n=0;RSB_LIKELY(n<nnz);++n)
			if(IA[n]>JA[n])
				return RSB_BOOL_TRUE;

	return RSB_BOOL_FALSE;
}

void rsb__util_coo_upper_to_lower_symmetric(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz)
{
	/*!
		\ingroup gr_internals
		A transpose function.
	*/
	rsb_nnz_idx_t n;
	
	for(n=0;n<nnz;++n)
		if(IA[n]<JA[n])
			RSB_SWAP(rsb_coo_idx_t,IA[n],JA[n]);
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__util_coo_lower_to_upper_symmetric(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz)
{
	/*!
		\ingroup gr_internals
	*/
	rsb__util_coo_upper_to_lower_symmetric(JA, IA, nnz);
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__util_coo_check_if_has_diagonal_elements(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_bool_t *has_diagonal_elements, rsb_bool_t wv)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_bitmap_data_t * bmap = NULL;
	rsb_nnz_idx_t n,cnnz = 0;

	if(RSB_INVALID_NNZ_INDEX(nnz) || RSB_INVALID_COO_INDEX(m) || !IA || !JA || !has_diagonal_elements)
		return RSB_ERR_BADARGS;

	if(m>nnz)
		goto missing_diagonal_element; /* trivially */
	/* we allow for duplicates, though (instead a count would be enough) */

	bmap = rsb__allocate_bitmap(1,m);

	if(!bmap)
		return RSB_ERR_ENOMEM;

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		if(RSB_UNLIKELY(IA[n]==JA[n]))
		{
			if(	IA[n]<0 || IA[n]>m || JA[n]<0 || JA[n]>m)
			{
				RSB_ERROR(RSB_ERRMSG_BADCOO"\n");
				goto badinput;
			}
			RSB_BITMAP_SET(bmap,1,nnz,0,IA[n]);
			++cnnz;
		}
	}

	if(cnnz<m)
	if(wv)
	{
		RSB_STDERR("Missing %zd diagonal elements.\n",(size_t)(m-cnnz));
	}

	for(n=0;RSB_LIKELY(n<m);++n)
		if(!RSB_BITMAP_GET(bmap,1,nnz,0,n))
		{
			if(wv)
				RSB_STDERR("Missing element %zd.\n",(size_t)n);
			goto missing_diagonal_element;
		}
goto ok;

ok:
	RSB_CONDITIONAL_FREE(bmap);
	*has_diagonal_elements = RSB_BOOL_TRUE;
	return RSB_ERR_NO_ERROR;

missing_diagonal_element:
	RSB_CONDITIONAL_FREE(bmap);
	*has_diagonal_elements = RSB_BOOL_FALSE;
	return RSB_ERR_NO_ERROR;

badinput:
	RSB_CONDITIONAL_FREE(bmap);
	return RSB_ERR_BADARGS;
}

rsb_bool_t rsb__util_is_halfword_coo_array_sorted_up(const rsb_half_idx_t* p, const rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i;
	if(n<2)
		return RSB_BOOL_TRUE;
	for(i=1;RSB_LIKELY(i<n);++i)
		if(RSB_UNLIKELY(p[i-1]>=p[i]))
			return RSB_BOOL_FALSE;
	return RSB_BOOL_TRUE;
}

rsb_bool_t rsb__util_is_halfword_coo_array_sorted_up_partial_order(const rsb_half_idx_t * p, const rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i;
	if(n<2)
		return RSB_BOOL_TRUE;
	for(i=1;RSB_LIKELY(i<n);++i)
		if(RSB_UNLIKELY(p[i-1]>p[i]))
		{
			/* RSB_STDOUT("n=%d, p[%d-1]>p[%d] -- %d > %d\n",n,i,i,p[i-1],p[i]); */
			return RSB_BOOL_FALSE;
		}
	return RSB_BOOL_TRUE;
}

rsb_bool_t rsb__util_is_nnz_array_sorted_up_partial_order(const rsb_nnz_idx_t * p, const rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i;
	if(n<2)
		goto yes;
	for(i=1;RSB_LIKELY(i<n);++i)
		if(RSB_UNLIKELY(p[i-1]>p[i]))
		{
			/* RSB_STDOUT("n=%d, p[%d-1]>p[%d] -- %d > %d\n",n,i,i,p[i-1],p[i]); */
			goto no;
		}
yes:
	return RSB_BOOL_TRUE;
no:
	/* RSB_STDOUT("%d=p[%d]>=p[%d]=%d\n",p[i-1],i-1,i,p[i]); */
	return RSB_BOOL_FALSE;
}

rsb_bool_t rsb__util_is_coo_array_sorted_up_partial_order(const rsb_nnz_idx_t * p, const rsb_nnz_idx_t n)
{
	return rsb__util_is_nnz_array_sorted_up_partial_order(p, n);
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__util_is_coo_array_sorted_up(const rsb_coo_idx_t * p, const rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i;
	if(n<2)
		goto yes;
	for(i=1;RSB_LIKELY(i<n);++i)
		if(RSB_UNLIKELY(p[i-1]>=p[i]))
			goto no;
yes:
	return RSB_BOOL_TRUE;
no:
	/* TODO: if this becomes a debugging function, one can use a RSB_ERROR printout instead */
	/* RSB_STDOUT("%d=p[%d]>=p[%d]=%d\n",p[i-1],i-1,i,p[i]); */
	return RSB_BOOL_FALSE;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/

rsb_bool_t rsb__util_is_nnz_array_sorted_up(const rsb_nnz_idx_t * p, const rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i;
	if(n<2)
		goto yes;
	for(i=1;RSB_LIKELY(i<n);++i)
		if(RSB_UNLIKELY(p[i-1]>=p[i]))
			goto no;
yes:
	return RSB_BOOL_TRUE;
no:
	return RSB_BOOL_FALSE;
}

#if RSB_WANT_BCSC_LEAVES  
void rsb__util_nnz_array_add_array(rsb_nnz_idx_t * p, const rsb_nnz_idx_t * q, rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
		Adds s from the given matrix coordinate indices array.
	*/
	register rsb_nnz_idx_t k;

	RSB_DEBUG_ASSERT(p || n==0);
	RSB_DEBUG_ASSERT(q || n==0);
	RSB_DEBUG_ASSERT(n>=0);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(n));

	for(k=0;RSB_LIKELY(k<n);++k)
	{
		p[k] += q[k];
	}
}
#endif /* RSB_WANT_BCSC_LEAVES */

rsb_int_t rsb__util_find_max_index(const rsb_int_t * p, rsb_int_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_int_t i,l=0;
	rsb_int_t m;
	if(n<1)
		goto ret;
	m = p[l];
	for(i=1;RSB_LIKELY(i<n);++i)
		if(p[i]>m)
			m = p[i], l = i;
ret:
	return l;
}

rsb_int_t rsb__util_find_min_index(const rsb_int_t * p, rsb_int_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_int_t i,l=0;
	rsb_int_t m;
	if(n<1)
		goto ret;
	m = p[l];
	for(i=1;RSB_LIKELY(i<n);++i)
		if(p[i]<m)
			m = p[i], l = i;
ret:
	return l;
}

rsb_nnz_idx_t rsb__util_find_coo_max_index_val(const rsb_nnz_idx_t * p, rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i;
	rsb_nnz_idx_t m;
	if(n<1)
		return RSB_MARKER_COO_VALUE;
	m = p[0];
	for(i=1;RSB_LIKELY(i<n);++i)
		if(p[i]>m)
			m = p[i];
	return m;
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_nnz_idx_t rsb__util_find_coo_min_index_val(const rsb_nnz_idx_t * p, rsb_nnz_idx_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i;
	rsb_nnz_idx_t m;
	if(n<1)
		return RSB_MARKER_COO_VALUE;
	m = p[0];
	for(i=1;RSB_LIKELY(i<n);++i)
		if(p[i]<m)
			m = p[i];
	return m;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_int_t rsb__util_find_max_index_val(const rsb_int_t * p, rsb_int_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_int_t i;
	rsb_int_t m;
	if(n<1)
		return RSB_MARKER_INT_VALUE;
	m = p[0];
	for(i=1;RSB_LIKELY(i<n);++i)
		if(p[i]>m)
			m = p[i];
	return m;
}

rsb_int_t rsb__util_find_min_index_val(const rsb_int_t * p, rsb_int_t n)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_int_t i;
	rsb_int_t m;
	if(n<1)
		return RSB_MARKER_INT_VALUE;
	m = p[0];
	for(i=1;RSB_LIKELY(i<n);++i)
		if(p[i]<m)
			m = p[i];
	return m;
}

rsb_coo_idx_t rsb__util_find_coo_val_idx(const rsb_coo_idx_t * p, const rsb_nnz_idx_t n, const rsb_coo_idx_t m)
{
	/*!
		\ingroup gr_internals
		Return index of m if found, n otherwise.
	*/
	rsb_nnz_idx_t i;

	for ( i=0; RSB_LIKELY(i<n); ++i )
		if( p[i] == m )
			return i;
	return n;
}

void rsb__util_coo_array_renumber(rsb_coo_idx_t * a, const rsb_coo_idx_t * iren, rsb_nnz_idx_t n, rsb_flags_t aflags, rsb_flags_t pflags, rsb_flags_t oflags)
{
	/*!
		\ingroup gr_internals
	*/
	rsb_nnz_idx_t i = 0;
	rsb_coo_idx_t ioff = 0,ooff = 0,oooff = 0;

	if( aflags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
		ioff = 1;
	if( pflags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
		ooff = 1;
	if( oflags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
		oooff = 1;

	for(i=0;RSB_LIKELY(i<n);++i)
		a[i] = iren[a[i]-ioff]-ooff+oooff;
}

rsb_err_t rsb__util_compress_to_row_pointers_array(rsb_coo_idx_t * RSB_RESTRICT pa, rsb_nnz_idx_t nz, rsb_coo_idx_t m, rsb_flags_t iflags, rsb_flags_t oflags, rsb_coo_idx_t * ta)
{
	/* 
	 * Note that this routine invokes OpenMP.
	 * Requires m+1 temporary space.
	 * TODO: rsb__util_compress_to_row_pointers_array -> rsb__idx_fia2fpa
	 * */
	rsb_nnz_idx_t i;
	rsb_bool_t sa = RSB_BOOL_TRUE;
	rsb_nnz_idx_t ifo = ( iflags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )?1:0;
	rsb_nnz_idx_t ofo = ( oflags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )?1:0;

	if(!ta)
		return RSB_ERR_BADARGS;
	if(!pa)
		pa = rsb__calloc((m+1)*sizeof(rsb_coo_idx_t));
	else 
		sa = RSB_BOOL_FALSE,
		RSB_BZERO(pa,(m+1)*sizeof(rsb_coo_idx_t));
	if(!pa)
		goto err;
	for(i=0;RSB_LIKELY(i<nz);++i)
		pa[1+ta[i]-ifo]++;
	for(i=0;RSB_LIKELY(i<m );++i)
		pa[i+1] += pa[i]; /* TODO: need a prefix sum routine */
	if(ofo)
		for(i=0;RSB_LIKELY(i<m+1);++i)
			pa[i]++;

	RSB_COA_MEMCPY_parallel(ta,pa,0,0,m+1);
	if(sa)
		RSB_CONDITIONAL_FREE(pa);
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_ENOMEM; 
}

rsb_err_t rsb__util_uncompress_row_pointers_array(const rsb_coo_idx_t * pa, rsb_nnz_idx_t n, rsb_flags_t iflags, rsb_flags_t oflags, rsb_coo_idx_t * ta)
{
	/*
	 TODO: write a version to exploit no-pointer-aliasing (if available)
	 */
	rsb_nnz_idx_t i,nz;
	rsb_bool_t sa = RSB_BOOL_TRUE; /* same array */
	rsb_nnz_idx_t ifo = ( iflags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )?1:0;
	rsb_nnz_idx_t ofo = ( oflags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )?1:0;
	rsb_coo_idx_t * ota = ta, *tmp = NULL;

	if(!pa || !ta)
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_BADARGS;
	}
	nz = pa[n]-ifo;
	if(nz==0)
		goto ret;
	if(ta==pa)
	{
		if(RSB_LIKELY(n+1<nz))
			pa = tmp = rsb__clone_area(pa,sizeof(rsb_coo_idx_t)*(n+1));
		else
			tmp = ta = rsb__clone_area(pa,sizeof(rsb_coo_idx_t)* nz  );
	}
	else 
		sa = RSB_BOOL_FALSE;
	if((!ta) || (!pa))
		goto err;
	/* TODO: this shall be parallel! */
	for(i=0;RSB_LIKELY(i<n);++i)
		rsb__util_coo_array_set(ta+(pa[i]-ifo),(pa[i+1]-pa[i]),i+ofo);
	if(sa && !(n+1<nz))
		RSB_COA_MEMCPY_parallel(ota,ta,0,0,nz);
	RSB_CONDITIONAL_FREE(tmp);
ret:
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_ENOMEM; 
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__debug_print_index_vector(const rsb_coo_idx_t * v1, size_t n){
	/*! 
	 **/
#if RSB_ALLOW_STDOUT
	size_t i = 0;
	int ioff = 1,voff = 1;
	if(!v1)
		return RSB_ERR_BADARGS;

	for(i=0;i<n ;++i) 
		RSB_STDOUT("%zd : %ld \n",(rsb_printf_int_t)(i+ioff),(long int)(v1[i]+voff));
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__debug_print_index_vectors_diff(const rsb_coo_idx_t * v1, const rsb_coo_idx_t * v2, size_t n, int onlyfirst){
	/*! 
	 **/
#if RSB_ALLOW_STDOUT
	size_t i,differing = 0;
	if(!v1 || !v2)return RSB_ERR_BADARGS;

	RSB_STDERR("\t indices vectors diff :\n");
	
		for(i=0;i<n ;++i) 
			if(v1[i]!=v2[i]){
		differing++;
		if((onlyfirst==0)||(onlyfirst>differing))
		RSB_STDOUT("%zd : %ld %ld \n",(rsb_printf_int_t)i,						(long int)v1[i],(long int)v2[i]		);
			}
		if(differing>onlyfirst)RSB_STDOUT("...and %zd more ...\n",(rsb_printf_int_t)(differing-onlyfirst));
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__util_find_extremal_half_index_val(const rsb_half_idx_t * RSB_RESTRICT p, rsb_nnz_idx_t n, rsb_coo_idx_t lb, rsb_coo_idx_t ub, rsb_half_idx_t *lf, rsb_half_idx_t * RSB_RESTRICT uf)
{
	/* TODO: remove the useless 'ub' argument */
	/* TODO: this is a naive implementation; need a better one */
	rsb_half_idx_t vm = RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t),vM = 0;
	rsb_nnz_idx_t i = 0;

	for(i=0;i<n;++i)
		vm = RSB_MIN(vm,p[i]), vM = RSB_MAX(vM,p[i]);
	if(lf)
		*lf = vm;
	if(uf)
		*uf = vM;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_find_extremal_full_index_val(const rsb_coo_idx_t * RSB_RESTRICT p, rsb_nnz_idx_t n, rsb_coo_idx_t lb, rsb_coo_idx_t ub, rsb_coo_idx_t * RSB_RESTRICT lf, rsb_coo_idx_t * RSB_RESTRICT uf)
{
	/* TODO: remove the useless 'ub' argument */
	/* TODO: this is a naive implementation; need a better one */
	rsb_coo_idx_t vm = RSB_MAX_VALUE_FOR_TYPE(rsb_coo_idx_t),vM = 0;
	rsb_nnz_idx_t i = 0;

	for(i=0;i<n;++i)
		vm = RSB_MIN(vm,p[i]), vM = RSB_MAX(vM,p[i]);
	RSB_ASSIGN_IF_DP(lf,vm);
	RSB_ASSIGN_IF_DP(uf,vM);
	return RSB_ERR_NO_ERROR;

}

rsb_bool_t rsb__util_reverse_halfword_coo_array(rsb_half_idx_t* p, rsb_nnz_idx_t n)
{
	rsb_nnz_idx_t nzi;
	--n;
	for(nzi=0;nzi<(n+1)/2;++nzi)
		RSB_SWAP(rsb_half_idx_t,p[n-nzi],p[nzi]);
	return RSB_BOOL_TRUE;
}

rsb_bool_t rsb__util_reverse_fullword_coo_array(rsb_coo_idx_t* p, rsb_nnz_idx_t n)
{
	rsb_nnz_idx_t nzi;
	--n;
	for(nzi=0;nzi<(n+1)/2;++nzi)
		RSB_SWAP(rsb_coo_idx_t,p[n-nzi],p[nzi]);
	return RSB_BOOL_TRUE;
}

/* @endcond */
