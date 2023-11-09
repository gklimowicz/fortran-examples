/* @cond INNERDOC */
/**
 * @file
 * @brief
 * Permutation functions.
 */

/*

Copyright (C) 2008-2022 Michele Martone

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
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

 */



#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "rsb_common.h"

rsb_err_t rsb__do_permute_values_in_place_with_coo_index(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type){
	/**
		This is in place swapping, and it is much slower than not in place.
		
		FIXME : document 

		Sadly, this is O(nnz^2).
		... Or qsorts order ?
	 */
	rsb_coo_idx_t n;/* this is the case where coo cannot overflow */

	switch(type)
	{
		/* supported (double,float,float complex,double complex) */
	case RSB_NUMERICAL_TYPE_DOUBLE 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_coo_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(double,((double*)VA)[n],((double*)VA)[t]);
	}
			break;
	case RSB_NUMERICAL_TYPE_FLOAT 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_coo_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(float,((float*)VA)[n],((float*)VA)[t]);
	}
			break;
	case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_coo_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(float complex,((float complex*)VA)[n],((float complex*)VA)[t]);
	}
			break;
	case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_coo_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(double complex,((double complex*)VA)[n],((double complex*)VA)[t]);
	}
			break;
			/* unsupported type */
		default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__do_permute_values_in_place_with_nnz_index(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type){
		/*	
			This is in place swapping, and it is much slower than not in place.
			
			FIXME : document and finish (s/double/ * /).

			Sadly, this is O(nnz^2).
			... Or qsorts order ?
		 */
	rsb_nnz_idx_t n;

	switch(type)
	{
		/* supported (double,float,float complex,double complex) */
	case RSB_NUMERICAL_TYPE_DOUBLE 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_nnz_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(double,((double*)VA)[n],((double*)VA)[t]);
	}
			break;
	case RSB_NUMERICAL_TYPE_FLOAT 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_nnz_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(float,((float*)VA)[n],((float*)VA)[t]);
	}
			break;
	case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_nnz_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(float complex,((float complex*)VA)[n],((float complex*)VA)[t]);
	}
			break;
	case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:

	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		register rsb_nnz_idx_t t,m=n;
		RSB_DEBUG_ASSERT(K[n]>=0);
		if(K[n]==n)
			continue;
		if(K[n]>n)
			t=K[n];
		else
		{
			/* follow swap chain */
			while(K[K[m]]<n)
				m=K[m];

			t=K[K[m]];
			RSB_DEBUG_ASSERT(t>=0);

#if RSB_DEBUG_SORT_STUFF 
			K[K[m]]=-1;	// just a debug measure
#endif /* RSB_DEBUG_SORT_STUFF */
			K[m]=t;
		}
		/* perform the swap */
		RSB_SWAP(rsb_coo_idx_t,IA[n],IA[t]);
		RSB_SWAP(rsb_coo_idx_t,JA[n],JA[t]);
		RSB_SWAP(double complex,((double complex*)VA)[n],((double complex*)VA)[t]);
	}
			break;
			/* unsupported type */
		default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__do_permute_values_with_coo_index( void * rVA, const void *VA, rsb_coo_idx_t * rIA, const rsb_coo_idx_t * IA, rsb_coo_idx_t * rJA, const rsb_coo_idx_t * JA, const rsb_coo_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type)	{
		/*
		 * FIXME : UNOPTIMIZED !
		 */
		rsb_coo_idx_t i;/* in this algorithm, coo cannot overflow */

		/* should permute here */
		for(i=0;RSB_LIKELY(i<nnz);++i)
		{
			RSB_DEBUG_ASSERT(K[i]>=0);
			RSB_DEBUG_ASSERT(K[i]<nnz);

			rIA [i]=IA [K[i]];
			rJA [i]=JA [K[i]];
		}

		switch(type)
		{
			/* supported (double,float,float complex,double complex) */
			case RSB_NUMERICAL_TYPE_DOUBLE 	:
		{
for(i=0;i+15<nnz;i+=16){
((double*)rVA)[i+0 ]=((double*)VA)[K[(i+0 )]];
	((double*)rVA)[i+1 ]=((double*)VA)[K[(i+1 )]];
	((double*)rVA)[i+2 ]=((double*)VA)[K[(i+2 )]];
	((double*)rVA)[i+3 ]=((double*)VA)[K[(i+3 )]];
	((double*)rVA)[i+4 ]=((double*)VA)[K[(i+4 )]];
	((double*)rVA)[i+5 ]=((double*)VA)[K[(i+5 )]];
	((double*)rVA)[i+6 ]=((double*)VA)[K[(i+6 )]];
	((double*)rVA)[i+7 ]=((double*)VA)[K[(i+7 )]];
	((double*)rVA)[i+8 ]=((double*)VA)[K[(i+8 )]];
	((double*)rVA)[i+9 ]=((double*)VA)[K[(i+9 )]];
	((double*)rVA)[i+10 ]=((double*)VA)[K[(i+10 )]];
	((double*)rVA)[i+11 ]=((double*)VA)[K[(i+11 )]];
	((double*)rVA)[i+12 ]=((double*)VA)[K[(i+12 )]];
	((double*)rVA)[i+13 ]=((double*)VA)[K[(i+13 )]];
	((double*)rVA)[i+14 ]=((double*)VA)[K[(i+14 )]];
	((double*)rVA)[i+15 ]=((double*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((double*)rVA)[i+0 ]=((double*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			case RSB_NUMERICAL_TYPE_FLOAT 	:
		{
for(i=0;i+15<nnz;i+=16){
((float*)rVA)[i+0 ]=((float*)VA)[K[(i+0 )]];
	((float*)rVA)[i+1 ]=((float*)VA)[K[(i+1 )]];
	((float*)rVA)[i+2 ]=((float*)VA)[K[(i+2 )]];
	((float*)rVA)[i+3 ]=((float*)VA)[K[(i+3 )]];
	((float*)rVA)[i+4 ]=((float*)VA)[K[(i+4 )]];
	((float*)rVA)[i+5 ]=((float*)VA)[K[(i+5 )]];
	((float*)rVA)[i+6 ]=((float*)VA)[K[(i+6 )]];
	((float*)rVA)[i+7 ]=((float*)VA)[K[(i+7 )]];
	((float*)rVA)[i+8 ]=((float*)VA)[K[(i+8 )]];
	((float*)rVA)[i+9 ]=((float*)VA)[K[(i+9 )]];
	((float*)rVA)[i+10 ]=((float*)VA)[K[(i+10 )]];
	((float*)rVA)[i+11 ]=((float*)VA)[K[(i+11 )]];
	((float*)rVA)[i+12 ]=((float*)VA)[K[(i+12 )]];
	((float*)rVA)[i+13 ]=((float*)VA)[K[(i+13 )]];
	((float*)rVA)[i+14 ]=((float*)VA)[K[(i+14 )]];
	((float*)rVA)[i+15 ]=((float*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((float*)rVA)[i+0 ]=((float*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:
		{
for(i=0;i+15<nnz;i+=16){
((float complex*)rVA)[i+0 ]=((float complex*)VA)[K[(i+0 )]];
	((float complex*)rVA)[i+1 ]=((float complex*)VA)[K[(i+1 )]];
	((float complex*)rVA)[i+2 ]=((float complex*)VA)[K[(i+2 )]];
	((float complex*)rVA)[i+3 ]=((float complex*)VA)[K[(i+3 )]];
	((float complex*)rVA)[i+4 ]=((float complex*)VA)[K[(i+4 )]];
	((float complex*)rVA)[i+5 ]=((float complex*)VA)[K[(i+5 )]];
	((float complex*)rVA)[i+6 ]=((float complex*)VA)[K[(i+6 )]];
	((float complex*)rVA)[i+7 ]=((float complex*)VA)[K[(i+7 )]];
	((float complex*)rVA)[i+8 ]=((float complex*)VA)[K[(i+8 )]];
	((float complex*)rVA)[i+9 ]=((float complex*)VA)[K[(i+9 )]];
	((float complex*)rVA)[i+10 ]=((float complex*)VA)[K[(i+10 )]];
	((float complex*)rVA)[i+11 ]=((float complex*)VA)[K[(i+11 )]];
	((float complex*)rVA)[i+12 ]=((float complex*)VA)[K[(i+12 )]];
	((float complex*)rVA)[i+13 ]=((float complex*)VA)[K[(i+13 )]];
	((float complex*)rVA)[i+14 ]=((float complex*)VA)[K[(i+14 )]];
	((float complex*)rVA)[i+15 ]=((float complex*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((float complex*)rVA)[i+0 ]=((float complex*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:
		{
for(i=0;i+15<nnz;i+=16){
((double complex*)rVA)[i+0 ]=((double complex*)VA)[K[(i+0 )]];
	((double complex*)rVA)[i+1 ]=((double complex*)VA)[K[(i+1 )]];
	((double complex*)rVA)[i+2 ]=((double complex*)VA)[K[(i+2 )]];
	((double complex*)rVA)[i+3 ]=((double complex*)VA)[K[(i+3 )]];
	((double complex*)rVA)[i+4 ]=((double complex*)VA)[K[(i+4 )]];
	((double complex*)rVA)[i+5 ]=((double complex*)VA)[K[(i+5 )]];
	((double complex*)rVA)[i+6 ]=((double complex*)VA)[K[(i+6 )]];
	((double complex*)rVA)[i+7 ]=((double complex*)VA)[K[(i+7 )]];
	((double complex*)rVA)[i+8 ]=((double complex*)VA)[K[(i+8 )]];
	((double complex*)rVA)[i+9 ]=((double complex*)VA)[K[(i+9 )]];
	((double complex*)rVA)[i+10 ]=((double complex*)VA)[K[(i+10 )]];
	((double complex*)rVA)[i+11 ]=((double complex*)VA)[K[(i+11 )]];
	((double complex*)rVA)[i+12 ]=((double complex*)VA)[K[(i+12 )]];
	((double complex*)rVA)[i+13 ]=((double complex*)VA)[K[(i+13 )]];
	((double complex*)rVA)[i+14 ]=((double complex*)VA)[K[(i+14 )]];
	((double complex*)rVA)[i+15 ]=((double complex*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((double complex*)rVA)[i+0 ]=((double complex*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			/* unsupported type */
			default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
		}
		return RSB_ERR_NO_ERROR;
	}


rsb_err_t rsb__do_permute_rows_with_coo_index( rsb_coo_idx_t * IA, const rsb_coo_idx_t * K, rsb_nnz_idx_t nnz)	{
		/*
		 * FIXME : UNOPTIMIZED !
		 */
		rsb_coo_idx_t i;/* in this algorithm, coo cannot overflow */

		/* should permute here */
		for(i=0;RSB_LIKELY(i<nnz);++i)
		{
			RSB_DEBUG_ASSERT(K[i]>=0);
			RSB_DEBUG_ASSERT(K[i]<nnz);

			IA [i]=K[IA[i]];
		}
		return RSB_ERR_NO_ERROR;
	}


rsb_err_t rsb__do_permute_values_with_nnz_index( void * rVA, const void *VA, rsb_coo_idx_t * rIA, const rsb_coo_idx_t * IA, rsb_coo_idx_t * rJA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t typecode)	{
		/*
		 * FIXME : UNOPTIMIZED !
		 */
		rsb_nnz_idx_t i;

		/* should permute here */
		for(i=0;RSB_LIKELY(i<nnz);++i)
		{
			RSB_DEBUG_ASSERT(K[i]>=0);
			RSB_DEBUG_ASSERT(K[i]<nnz);

			rIA [i]=IA [K[i]];
			rJA [i]=JA [K[i]];
		}

		switch(typecode)
		{
			/* supported (double,float,float complex,double complex) */
			case RSB_NUMERICAL_TYPE_DOUBLE 	:
		{
for(i=0;i+15<nnz;i+=16){
((double*)rVA)[i+0 ]=((double*)VA)[K[(i+0 )]];
	((double*)rVA)[i+1 ]=((double*)VA)[K[(i+1 )]];
	((double*)rVA)[i+2 ]=((double*)VA)[K[(i+2 )]];
	((double*)rVA)[i+3 ]=((double*)VA)[K[(i+3 )]];
	((double*)rVA)[i+4 ]=((double*)VA)[K[(i+4 )]];
	((double*)rVA)[i+5 ]=((double*)VA)[K[(i+5 )]];
	((double*)rVA)[i+6 ]=((double*)VA)[K[(i+6 )]];
	((double*)rVA)[i+7 ]=((double*)VA)[K[(i+7 )]];
	((double*)rVA)[i+8 ]=((double*)VA)[K[(i+8 )]];
	((double*)rVA)[i+9 ]=((double*)VA)[K[(i+9 )]];
	((double*)rVA)[i+10 ]=((double*)VA)[K[(i+10 )]];
	((double*)rVA)[i+11 ]=((double*)VA)[K[(i+11 )]];
	((double*)rVA)[i+12 ]=((double*)VA)[K[(i+12 )]];
	((double*)rVA)[i+13 ]=((double*)VA)[K[(i+13 )]];
	((double*)rVA)[i+14 ]=((double*)VA)[K[(i+14 )]];
	((double*)rVA)[i+15 ]=((double*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((double*)rVA)[i+0 ]=((double*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			case RSB_NUMERICAL_TYPE_FLOAT 	:
		{
for(i=0;i+15<nnz;i+=16){
((float*)rVA)[i+0 ]=((float*)VA)[K[(i+0 )]];
	((float*)rVA)[i+1 ]=((float*)VA)[K[(i+1 )]];
	((float*)rVA)[i+2 ]=((float*)VA)[K[(i+2 )]];
	((float*)rVA)[i+3 ]=((float*)VA)[K[(i+3 )]];
	((float*)rVA)[i+4 ]=((float*)VA)[K[(i+4 )]];
	((float*)rVA)[i+5 ]=((float*)VA)[K[(i+5 )]];
	((float*)rVA)[i+6 ]=((float*)VA)[K[(i+6 )]];
	((float*)rVA)[i+7 ]=((float*)VA)[K[(i+7 )]];
	((float*)rVA)[i+8 ]=((float*)VA)[K[(i+8 )]];
	((float*)rVA)[i+9 ]=((float*)VA)[K[(i+9 )]];
	((float*)rVA)[i+10 ]=((float*)VA)[K[(i+10 )]];
	((float*)rVA)[i+11 ]=((float*)VA)[K[(i+11 )]];
	((float*)rVA)[i+12 ]=((float*)VA)[K[(i+12 )]];
	((float*)rVA)[i+13 ]=((float*)VA)[K[(i+13 )]];
	((float*)rVA)[i+14 ]=((float*)VA)[K[(i+14 )]];
	((float*)rVA)[i+15 ]=((float*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((float*)rVA)[i+0 ]=((float*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:
		{
for(i=0;i+15<nnz;i+=16){
((float complex*)rVA)[i+0 ]=((float complex*)VA)[K[(i+0 )]];
	((float complex*)rVA)[i+1 ]=((float complex*)VA)[K[(i+1 )]];
	((float complex*)rVA)[i+2 ]=((float complex*)VA)[K[(i+2 )]];
	((float complex*)rVA)[i+3 ]=((float complex*)VA)[K[(i+3 )]];
	((float complex*)rVA)[i+4 ]=((float complex*)VA)[K[(i+4 )]];
	((float complex*)rVA)[i+5 ]=((float complex*)VA)[K[(i+5 )]];
	((float complex*)rVA)[i+6 ]=((float complex*)VA)[K[(i+6 )]];
	((float complex*)rVA)[i+7 ]=((float complex*)VA)[K[(i+7 )]];
	((float complex*)rVA)[i+8 ]=((float complex*)VA)[K[(i+8 )]];
	((float complex*)rVA)[i+9 ]=((float complex*)VA)[K[(i+9 )]];
	((float complex*)rVA)[i+10 ]=((float complex*)VA)[K[(i+10 )]];
	((float complex*)rVA)[i+11 ]=((float complex*)VA)[K[(i+11 )]];
	((float complex*)rVA)[i+12 ]=((float complex*)VA)[K[(i+12 )]];
	((float complex*)rVA)[i+13 ]=((float complex*)VA)[K[(i+13 )]];
	((float complex*)rVA)[i+14 ]=((float complex*)VA)[K[(i+14 )]];
	((float complex*)rVA)[i+15 ]=((float complex*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((float complex*)rVA)[i+0 ]=((float complex*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:
		{
for(i=0;i+15<nnz;i+=16){
((double complex*)rVA)[i+0 ]=((double complex*)VA)[K[(i+0 )]];
	((double complex*)rVA)[i+1 ]=((double complex*)VA)[K[(i+1 )]];
	((double complex*)rVA)[i+2 ]=((double complex*)VA)[K[(i+2 )]];
	((double complex*)rVA)[i+3 ]=((double complex*)VA)[K[(i+3 )]];
	((double complex*)rVA)[i+4 ]=((double complex*)VA)[K[(i+4 )]];
	((double complex*)rVA)[i+5 ]=((double complex*)VA)[K[(i+5 )]];
	((double complex*)rVA)[i+6 ]=((double complex*)VA)[K[(i+6 )]];
	((double complex*)rVA)[i+7 ]=((double complex*)VA)[K[(i+7 )]];
	((double complex*)rVA)[i+8 ]=((double complex*)VA)[K[(i+8 )]];
	((double complex*)rVA)[i+9 ]=((double complex*)VA)[K[(i+9 )]];
	((double complex*)rVA)[i+10 ]=((double complex*)VA)[K[(i+10 )]];
	((double complex*)rVA)[i+11 ]=((double complex*)VA)[K[(i+11 )]];
	((double complex*)rVA)[i+12 ]=((double complex*)VA)[K[(i+12 )]];
	((double complex*)rVA)[i+13 ]=((double complex*)VA)[K[(i+13 )]];
	((double complex*)rVA)[i+14 ]=((double complex*)VA)[K[(i+14 )]];
	((double complex*)rVA)[i+15 ]=((double complex*)VA)[K[(i+15 )]];
	}
for(     ;i<nnz;++i){ ((double complex*)rVA)[i+0 ]=((double complex*)VA)[K[(i+0 )]];
	 }
}

			
			break;
			/* unsupported type */
			default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
		}
		return RSB_ERR_NO_ERROR;
	}


void rsb__ip_reord(rsb_nnz_idx_t n, void * VAp, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * P, rsb_type_t typecode){
	/**
		This is an adapted PSBLAS psb_ip_reord_d1i2 routine.
	*/
	
	switch(typecode)
	{
			/* supported (double,float,float complex,double complex) */
	case RSB_NUMERICAL_TYPE_DOUBLE 	:
	{
		rsb_coo_idx_t isw1, isw2;
		rsb_nnz_idx_t lswap, lp, k;
		double swap;
		double * VA=VAp;

		lp = P[0];
		k  = 1;
		while(1)
		{
			if (RSB_UNLIKELY((lp==0) || (k>n))) break;
			while(1)
			{
				if (lp >= k) break;
				lp = P[lp];
			}
			lswap    = P[lp];
			P[lp]  = P[k];
			P[k]   = lp;
			--lp;
			--k;
			swap   = VA[lp];
			VA[lp] = VA[k];
			VA[k]  = swap;
			isw1   = IA[lp];
			IA[lp] = IA[k];
			IA[k]  = isw1;
			isw2   = JA[lp];
			JA[lp] = JA[k];
			JA[k]  = isw2;
			++k;
			lp = lswap ;
			k  = k + 1;
		}
	}
		break;
		case RSB_NUMERICAL_TYPE_FLOAT 	:
	{
		rsb_coo_idx_t isw1, isw2;
		rsb_nnz_idx_t lswap, lp, k;
		float swap;
		float * VA=VAp;

		lp = P[0];
		k  = 1;
		while(1)
		{
			if (RSB_UNLIKELY((lp==0) || (k>n))) break;
			while(1)
			{
				if (lp >= k) break;
				lp = P[lp];
			}
			lswap    = P[lp];
			P[lp]  = P[k];
			P[k]   = lp;
			--lp;
			--k;
			swap   = VA[lp];
			VA[lp] = VA[k];
			VA[k]  = swap;
			isw1   = IA[lp];
			IA[lp] = IA[k];
			IA[k]  = isw1;
			isw2   = JA[lp];
			JA[lp] = JA[k];
			JA[k]  = isw2;
			++k;
			lp = lswap ;
			k  = k + 1;
		}
	}
		break;
		case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:
	{
		rsb_coo_idx_t isw1, isw2;
		rsb_nnz_idx_t lswap, lp, k;
		float complex swap;
		float complex * VA=VAp;

		lp = P[0];
		k  = 1;
		while(1)
		{
			if (RSB_UNLIKELY((lp==0) || (k>n))) break;
			while(1)
			{
				if (lp >= k) break;
				lp = P[lp];
			}
			lswap    = P[lp];
			P[lp]  = P[k];
			P[k]   = lp;
			--lp;
			--k;
			swap   = VA[lp];
			VA[lp] = VA[k];
			VA[k]  = swap;
			isw1   = IA[lp];
			IA[lp] = IA[k];
			IA[k]  = isw1;
			isw2   = JA[lp];
			JA[lp] = JA[k];
			JA[k]  = isw2;
			++k;
			lp = lswap ;
			k  = k + 1;
		}
	}
		break;
		case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:
	{
		rsb_coo_idx_t isw1, isw2;
		rsb_nnz_idx_t lswap, lp, k;
		double complex swap;
		double complex * VA=VAp;

		lp = P[0];
		k  = 1;
		while(1)
		{
			if (RSB_UNLIKELY((lp==0) || (k>n))) break;
			while(1)
			{
				if (lp >= k) break;
				lp = P[lp];
			}
			lswap    = P[lp];
			P[lp]  = P[k];
			P[k]   = lp;
			--lp;
			--k;
			swap   = VA[lp];
			VA[lp] = VA[k];
			VA[k]  = swap;
			isw1   = IA[lp];
			IA[lp] = IA[k];
			IA[k]  = isw1;
			isw2   = JA[lp];
			JA[lp] = JA[k];
			JA[k]  = isw2;
			++k;
			lp = lswap ;
			k  = k + 1;
		}
	}
		break;
	
		/* unsupported type */
		default :
			return;
	}


}

void rsb__util_do_scatter_rows(void * RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, const void * RSB_RESTRICT iVA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT PA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode){
	/**
		This is an adapted PSBLAS psb_ip_reord_d1i2 routine.
	*/
	
	switch(typecode)
	{
			/* supported (double,float,float complex,double complex) */
	case RSB_NUMERICAL_TYPE_DOUBLE 	:
	{
		rsb_nnz_idx_t nzi;
		double*VA=(double*)oVA;
		{
for(nzi=0;nzi+15<nnz;nzi+=16){
VA[PA[IA[nzi+0 ]]]=((double*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		VA[PA[IA[nzi+1 ]]]=((double*)iVA)[nzi+1 ];
			oIA[PA[IA[nzi+1 ]]]=IA[nzi+1 ];
			oJA[PA[IA[nzi+1 ]]]=JA[nzi+1 ];
			PA[IA[nzi+1 ]]++;
		VA[PA[IA[nzi+2 ]]]=((double*)iVA)[nzi+2 ];
			oIA[PA[IA[nzi+2 ]]]=IA[nzi+2 ];
			oJA[PA[IA[nzi+2 ]]]=JA[nzi+2 ];
			PA[IA[nzi+2 ]]++;
		VA[PA[IA[nzi+3 ]]]=((double*)iVA)[nzi+3 ];
			oIA[PA[IA[nzi+3 ]]]=IA[nzi+3 ];
			oJA[PA[IA[nzi+3 ]]]=JA[nzi+3 ];
			PA[IA[nzi+3 ]]++;
		VA[PA[IA[nzi+4 ]]]=((double*)iVA)[nzi+4 ];
			oIA[PA[IA[nzi+4 ]]]=IA[nzi+4 ];
			oJA[PA[IA[nzi+4 ]]]=JA[nzi+4 ];
			PA[IA[nzi+4 ]]++;
		VA[PA[IA[nzi+5 ]]]=((double*)iVA)[nzi+5 ];
			oIA[PA[IA[nzi+5 ]]]=IA[nzi+5 ];
			oJA[PA[IA[nzi+5 ]]]=JA[nzi+5 ];
			PA[IA[nzi+5 ]]++;
		VA[PA[IA[nzi+6 ]]]=((double*)iVA)[nzi+6 ];
			oIA[PA[IA[nzi+6 ]]]=IA[nzi+6 ];
			oJA[PA[IA[nzi+6 ]]]=JA[nzi+6 ];
			PA[IA[nzi+6 ]]++;
		VA[PA[IA[nzi+7 ]]]=((double*)iVA)[nzi+7 ];
			oIA[PA[IA[nzi+7 ]]]=IA[nzi+7 ];
			oJA[PA[IA[nzi+7 ]]]=JA[nzi+7 ];
			PA[IA[nzi+7 ]]++;
		VA[PA[IA[nzi+8 ]]]=((double*)iVA)[nzi+8 ];
			oIA[PA[IA[nzi+8 ]]]=IA[nzi+8 ];
			oJA[PA[IA[nzi+8 ]]]=JA[nzi+8 ];
			PA[IA[nzi+8 ]]++;
		VA[PA[IA[nzi+9 ]]]=((double*)iVA)[nzi+9 ];
			oIA[PA[IA[nzi+9 ]]]=IA[nzi+9 ];
			oJA[PA[IA[nzi+9 ]]]=JA[nzi+9 ];
			PA[IA[nzi+9 ]]++;
		VA[PA[IA[nzi+10 ]]]=((double*)iVA)[nzi+10 ];
			oIA[PA[IA[nzi+10 ]]]=IA[nzi+10 ];
			oJA[PA[IA[nzi+10 ]]]=JA[nzi+10 ];
			PA[IA[nzi+10 ]]++;
		VA[PA[IA[nzi+11 ]]]=((double*)iVA)[nzi+11 ];
			oIA[PA[IA[nzi+11 ]]]=IA[nzi+11 ];
			oJA[PA[IA[nzi+11 ]]]=JA[nzi+11 ];
			PA[IA[nzi+11 ]]++;
		VA[PA[IA[nzi+12 ]]]=((double*)iVA)[nzi+12 ];
			oIA[PA[IA[nzi+12 ]]]=IA[nzi+12 ];
			oJA[PA[IA[nzi+12 ]]]=JA[nzi+12 ];
			PA[IA[nzi+12 ]]++;
		VA[PA[IA[nzi+13 ]]]=((double*)iVA)[nzi+13 ];
			oIA[PA[IA[nzi+13 ]]]=IA[nzi+13 ];
			oJA[PA[IA[nzi+13 ]]]=JA[nzi+13 ];
			PA[IA[nzi+13 ]]++;
		VA[PA[IA[nzi+14 ]]]=((double*)iVA)[nzi+14 ];
			oIA[PA[IA[nzi+14 ]]]=IA[nzi+14 ];
			oJA[PA[IA[nzi+14 ]]]=JA[nzi+14 ];
			PA[IA[nzi+14 ]]++;
		VA[PA[IA[nzi+15 ]]]=((double*)iVA)[nzi+15 ];
			oIA[PA[IA[nzi+15 ]]]=IA[nzi+15 ];
			oJA[PA[IA[nzi+15 ]]]=JA[nzi+15 ];
			PA[IA[nzi+15 ]]++;
		}
for(     ;nzi<nnz;++nzi){ VA[PA[IA[nzi+0 ]]]=((double*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		 }
}

	}
		break;
		case RSB_NUMERICAL_TYPE_FLOAT 	:
	{
		rsb_nnz_idx_t nzi;
		float*VA=(float*)oVA;
		{
for(nzi=0;nzi+15<nnz;nzi+=16){
VA[PA[IA[nzi+0 ]]]=((float*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		VA[PA[IA[nzi+1 ]]]=((float*)iVA)[nzi+1 ];
			oIA[PA[IA[nzi+1 ]]]=IA[nzi+1 ];
			oJA[PA[IA[nzi+1 ]]]=JA[nzi+1 ];
			PA[IA[nzi+1 ]]++;
		VA[PA[IA[nzi+2 ]]]=((float*)iVA)[nzi+2 ];
			oIA[PA[IA[nzi+2 ]]]=IA[nzi+2 ];
			oJA[PA[IA[nzi+2 ]]]=JA[nzi+2 ];
			PA[IA[nzi+2 ]]++;
		VA[PA[IA[nzi+3 ]]]=((float*)iVA)[nzi+3 ];
			oIA[PA[IA[nzi+3 ]]]=IA[nzi+3 ];
			oJA[PA[IA[nzi+3 ]]]=JA[nzi+3 ];
			PA[IA[nzi+3 ]]++;
		VA[PA[IA[nzi+4 ]]]=((float*)iVA)[nzi+4 ];
			oIA[PA[IA[nzi+4 ]]]=IA[nzi+4 ];
			oJA[PA[IA[nzi+4 ]]]=JA[nzi+4 ];
			PA[IA[nzi+4 ]]++;
		VA[PA[IA[nzi+5 ]]]=((float*)iVA)[nzi+5 ];
			oIA[PA[IA[nzi+5 ]]]=IA[nzi+5 ];
			oJA[PA[IA[nzi+5 ]]]=JA[nzi+5 ];
			PA[IA[nzi+5 ]]++;
		VA[PA[IA[nzi+6 ]]]=((float*)iVA)[nzi+6 ];
			oIA[PA[IA[nzi+6 ]]]=IA[nzi+6 ];
			oJA[PA[IA[nzi+6 ]]]=JA[nzi+6 ];
			PA[IA[nzi+6 ]]++;
		VA[PA[IA[nzi+7 ]]]=((float*)iVA)[nzi+7 ];
			oIA[PA[IA[nzi+7 ]]]=IA[nzi+7 ];
			oJA[PA[IA[nzi+7 ]]]=JA[nzi+7 ];
			PA[IA[nzi+7 ]]++;
		VA[PA[IA[nzi+8 ]]]=((float*)iVA)[nzi+8 ];
			oIA[PA[IA[nzi+8 ]]]=IA[nzi+8 ];
			oJA[PA[IA[nzi+8 ]]]=JA[nzi+8 ];
			PA[IA[nzi+8 ]]++;
		VA[PA[IA[nzi+9 ]]]=((float*)iVA)[nzi+9 ];
			oIA[PA[IA[nzi+9 ]]]=IA[nzi+9 ];
			oJA[PA[IA[nzi+9 ]]]=JA[nzi+9 ];
			PA[IA[nzi+9 ]]++;
		VA[PA[IA[nzi+10 ]]]=((float*)iVA)[nzi+10 ];
			oIA[PA[IA[nzi+10 ]]]=IA[nzi+10 ];
			oJA[PA[IA[nzi+10 ]]]=JA[nzi+10 ];
			PA[IA[nzi+10 ]]++;
		VA[PA[IA[nzi+11 ]]]=((float*)iVA)[nzi+11 ];
			oIA[PA[IA[nzi+11 ]]]=IA[nzi+11 ];
			oJA[PA[IA[nzi+11 ]]]=JA[nzi+11 ];
			PA[IA[nzi+11 ]]++;
		VA[PA[IA[nzi+12 ]]]=((float*)iVA)[nzi+12 ];
			oIA[PA[IA[nzi+12 ]]]=IA[nzi+12 ];
			oJA[PA[IA[nzi+12 ]]]=JA[nzi+12 ];
			PA[IA[nzi+12 ]]++;
		VA[PA[IA[nzi+13 ]]]=((float*)iVA)[nzi+13 ];
			oIA[PA[IA[nzi+13 ]]]=IA[nzi+13 ];
			oJA[PA[IA[nzi+13 ]]]=JA[nzi+13 ];
			PA[IA[nzi+13 ]]++;
		VA[PA[IA[nzi+14 ]]]=((float*)iVA)[nzi+14 ];
			oIA[PA[IA[nzi+14 ]]]=IA[nzi+14 ];
			oJA[PA[IA[nzi+14 ]]]=JA[nzi+14 ];
			PA[IA[nzi+14 ]]++;
		VA[PA[IA[nzi+15 ]]]=((float*)iVA)[nzi+15 ];
			oIA[PA[IA[nzi+15 ]]]=IA[nzi+15 ];
			oJA[PA[IA[nzi+15 ]]]=JA[nzi+15 ];
			PA[IA[nzi+15 ]]++;
		}
for(     ;nzi<nnz;++nzi){ VA[PA[IA[nzi+0 ]]]=((float*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		 }
}

	}
		break;
		case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:
	{
		rsb_nnz_idx_t nzi;
		float complex*VA=(float complex*)oVA;
		{
for(nzi=0;nzi+15<nnz;nzi+=16){
VA[PA[IA[nzi+0 ]]]=((float complex*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		VA[PA[IA[nzi+1 ]]]=((float complex*)iVA)[nzi+1 ];
			oIA[PA[IA[nzi+1 ]]]=IA[nzi+1 ];
			oJA[PA[IA[nzi+1 ]]]=JA[nzi+1 ];
			PA[IA[nzi+1 ]]++;
		VA[PA[IA[nzi+2 ]]]=((float complex*)iVA)[nzi+2 ];
			oIA[PA[IA[nzi+2 ]]]=IA[nzi+2 ];
			oJA[PA[IA[nzi+2 ]]]=JA[nzi+2 ];
			PA[IA[nzi+2 ]]++;
		VA[PA[IA[nzi+3 ]]]=((float complex*)iVA)[nzi+3 ];
			oIA[PA[IA[nzi+3 ]]]=IA[nzi+3 ];
			oJA[PA[IA[nzi+3 ]]]=JA[nzi+3 ];
			PA[IA[nzi+3 ]]++;
		VA[PA[IA[nzi+4 ]]]=((float complex*)iVA)[nzi+4 ];
			oIA[PA[IA[nzi+4 ]]]=IA[nzi+4 ];
			oJA[PA[IA[nzi+4 ]]]=JA[nzi+4 ];
			PA[IA[nzi+4 ]]++;
		VA[PA[IA[nzi+5 ]]]=((float complex*)iVA)[nzi+5 ];
			oIA[PA[IA[nzi+5 ]]]=IA[nzi+5 ];
			oJA[PA[IA[nzi+5 ]]]=JA[nzi+5 ];
			PA[IA[nzi+5 ]]++;
		VA[PA[IA[nzi+6 ]]]=((float complex*)iVA)[nzi+6 ];
			oIA[PA[IA[nzi+6 ]]]=IA[nzi+6 ];
			oJA[PA[IA[nzi+6 ]]]=JA[nzi+6 ];
			PA[IA[nzi+6 ]]++;
		VA[PA[IA[nzi+7 ]]]=((float complex*)iVA)[nzi+7 ];
			oIA[PA[IA[nzi+7 ]]]=IA[nzi+7 ];
			oJA[PA[IA[nzi+7 ]]]=JA[nzi+7 ];
			PA[IA[nzi+7 ]]++;
		VA[PA[IA[nzi+8 ]]]=((float complex*)iVA)[nzi+8 ];
			oIA[PA[IA[nzi+8 ]]]=IA[nzi+8 ];
			oJA[PA[IA[nzi+8 ]]]=JA[nzi+8 ];
			PA[IA[nzi+8 ]]++;
		VA[PA[IA[nzi+9 ]]]=((float complex*)iVA)[nzi+9 ];
			oIA[PA[IA[nzi+9 ]]]=IA[nzi+9 ];
			oJA[PA[IA[nzi+9 ]]]=JA[nzi+9 ];
			PA[IA[nzi+9 ]]++;
		VA[PA[IA[nzi+10 ]]]=((float complex*)iVA)[nzi+10 ];
			oIA[PA[IA[nzi+10 ]]]=IA[nzi+10 ];
			oJA[PA[IA[nzi+10 ]]]=JA[nzi+10 ];
			PA[IA[nzi+10 ]]++;
		VA[PA[IA[nzi+11 ]]]=((float complex*)iVA)[nzi+11 ];
			oIA[PA[IA[nzi+11 ]]]=IA[nzi+11 ];
			oJA[PA[IA[nzi+11 ]]]=JA[nzi+11 ];
			PA[IA[nzi+11 ]]++;
		VA[PA[IA[nzi+12 ]]]=((float complex*)iVA)[nzi+12 ];
			oIA[PA[IA[nzi+12 ]]]=IA[nzi+12 ];
			oJA[PA[IA[nzi+12 ]]]=JA[nzi+12 ];
			PA[IA[nzi+12 ]]++;
		VA[PA[IA[nzi+13 ]]]=((float complex*)iVA)[nzi+13 ];
			oIA[PA[IA[nzi+13 ]]]=IA[nzi+13 ];
			oJA[PA[IA[nzi+13 ]]]=JA[nzi+13 ];
			PA[IA[nzi+13 ]]++;
		VA[PA[IA[nzi+14 ]]]=((float complex*)iVA)[nzi+14 ];
			oIA[PA[IA[nzi+14 ]]]=IA[nzi+14 ];
			oJA[PA[IA[nzi+14 ]]]=JA[nzi+14 ];
			PA[IA[nzi+14 ]]++;
		VA[PA[IA[nzi+15 ]]]=((float complex*)iVA)[nzi+15 ];
			oIA[PA[IA[nzi+15 ]]]=IA[nzi+15 ];
			oJA[PA[IA[nzi+15 ]]]=JA[nzi+15 ];
			PA[IA[nzi+15 ]]++;
		}
for(     ;nzi<nnz;++nzi){ VA[PA[IA[nzi+0 ]]]=((float complex*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		 }
}

	}
		break;
		case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:
	{
		rsb_nnz_idx_t nzi;
		double complex*VA=(double complex*)oVA;
		{
for(nzi=0;nzi+15<nnz;nzi+=16){
VA[PA[IA[nzi+0 ]]]=((double complex*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		VA[PA[IA[nzi+1 ]]]=((double complex*)iVA)[nzi+1 ];
			oIA[PA[IA[nzi+1 ]]]=IA[nzi+1 ];
			oJA[PA[IA[nzi+1 ]]]=JA[nzi+1 ];
			PA[IA[nzi+1 ]]++;
		VA[PA[IA[nzi+2 ]]]=((double complex*)iVA)[nzi+2 ];
			oIA[PA[IA[nzi+2 ]]]=IA[nzi+2 ];
			oJA[PA[IA[nzi+2 ]]]=JA[nzi+2 ];
			PA[IA[nzi+2 ]]++;
		VA[PA[IA[nzi+3 ]]]=((double complex*)iVA)[nzi+3 ];
			oIA[PA[IA[nzi+3 ]]]=IA[nzi+3 ];
			oJA[PA[IA[nzi+3 ]]]=JA[nzi+3 ];
			PA[IA[nzi+3 ]]++;
		VA[PA[IA[nzi+4 ]]]=((double complex*)iVA)[nzi+4 ];
			oIA[PA[IA[nzi+4 ]]]=IA[nzi+4 ];
			oJA[PA[IA[nzi+4 ]]]=JA[nzi+4 ];
			PA[IA[nzi+4 ]]++;
		VA[PA[IA[nzi+5 ]]]=((double complex*)iVA)[nzi+5 ];
			oIA[PA[IA[nzi+5 ]]]=IA[nzi+5 ];
			oJA[PA[IA[nzi+5 ]]]=JA[nzi+5 ];
			PA[IA[nzi+5 ]]++;
		VA[PA[IA[nzi+6 ]]]=((double complex*)iVA)[nzi+6 ];
			oIA[PA[IA[nzi+6 ]]]=IA[nzi+6 ];
			oJA[PA[IA[nzi+6 ]]]=JA[nzi+6 ];
			PA[IA[nzi+6 ]]++;
		VA[PA[IA[nzi+7 ]]]=((double complex*)iVA)[nzi+7 ];
			oIA[PA[IA[nzi+7 ]]]=IA[nzi+7 ];
			oJA[PA[IA[nzi+7 ]]]=JA[nzi+7 ];
			PA[IA[nzi+7 ]]++;
		VA[PA[IA[nzi+8 ]]]=((double complex*)iVA)[nzi+8 ];
			oIA[PA[IA[nzi+8 ]]]=IA[nzi+8 ];
			oJA[PA[IA[nzi+8 ]]]=JA[nzi+8 ];
			PA[IA[nzi+8 ]]++;
		VA[PA[IA[nzi+9 ]]]=((double complex*)iVA)[nzi+9 ];
			oIA[PA[IA[nzi+9 ]]]=IA[nzi+9 ];
			oJA[PA[IA[nzi+9 ]]]=JA[nzi+9 ];
			PA[IA[nzi+9 ]]++;
		VA[PA[IA[nzi+10 ]]]=((double complex*)iVA)[nzi+10 ];
			oIA[PA[IA[nzi+10 ]]]=IA[nzi+10 ];
			oJA[PA[IA[nzi+10 ]]]=JA[nzi+10 ];
			PA[IA[nzi+10 ]]++;
		VA[PA[IA[nzi+11 ]]]=((double complex*)iVA)[nzi+11 ];
			oIA[PA[IA[nzi+11 ]]]=IA[nzi+11 ];
			oJA[PA[IA[nzi+11 ]]]=JA[nzi+11 ];
			PA[IA[nzi+11 ]]++;
		VA[PA[IA[nzi+12 ]]]=((double complex*)iVA)[nzi+12 ];
			oIA[PA[IA[nzi+12 ]]]=IA[nzi+12 ];
			oJA[PA[IA[nzi+12 ]]]=JA[nzi+12 ];
			PA[IA[nzi+12 ]]++;
		VA[PA[IA[nzi+13 ]]]=((double complex*)iVA)[nzi+13 ];
			oIA[PA[IA[nzi+13 ]]]=IA[nzi+13 ];
			oJA[PA[IA[nzi+13 ]]]=JA[nzi+13 ];
			PA[IA[nzi+13 ]]++;
		VA[PA[IA[nzi+14 ]]]=((double complex*)iVA)[nzi+14 ];
			oIA[PA[IA[nzi+14 ]]]=IA[nzi+14 ];
			oJA[PA[IA[nzi+14 ]]]=JA[nzi+14 ];
			PA[IA[nzi+14 ]]++;
		VA[PA[IA[nzi+15 ]]]=((double complex*)iVA)[nzi+15 ];
			oIA[PA[IA[nzi+15 ]]]=IA[nzi+15 ];
			oJA[PA[IA[nzi+15 ]]]=JA[nzi+15 ];
			PA[IA[nzi+15 ]]++;
		}
for(     ;nzi<nnz;++nzi){ VA[PA[IA[nzi+0 ]]]=((double complex*)iVA)[nzi+0 ];
			oIA[PA[IA[nzi+0 ]]]=IA[nzi+0 ];
			oJA[PA[IA[nzi+0 ]]]=JA[nzi+0 ];
			PA[IA[nzi+0 ]]++;
		 }
}

	}
		break;
	
		/* unsupported type */
		default :
			return;
	}
}

rsb_err_t rsb__chk_permute(void){
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
	
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	double VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_coo_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	float VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_coo_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	float complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_coo_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	double complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_coo_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	double VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_coo_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	float VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_coo_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	float complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_coo_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	double complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_coo_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_coo_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}


{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	double VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_nnz_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	float VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_nnz_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	float complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_nnz_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	double complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {0,1,2};
	rsb_coo_idx_t JA[] = {0,1,2};
	rsb_nnz_idx_t  K[] = {0,1,2};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	double VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	float VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	float complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	double complex VA[] = {1,1,1};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_in_place_with_nnz_index(VA, IA, JA, K, nnz, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}


{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	double rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const float VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	float rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	float complex rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	double complex rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}


{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	double rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const float VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	float rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	float complex rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	double complex rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_coo_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}


{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_nnz_idx_t  K[] = {2,1,0};
	double rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_nnz_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const float VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_nnz_idx_t  K[] = {2,1,0};
	float rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_nnz_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_nnz_idx_t  K[] = {2,1,0};
	float complex rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_nnz_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_nnz_idx_t  K[] = {2,1,0};
	double complex rVA[] = {1,1,1};
	rsb_coo_idx_t rIA[] = {0,0,0};
	rsb_coo_idx_t rJA[] = {0,0,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_values_with_nnz_index(rVA,VA, rIA,IA, rJA,JA, K, nnz, typecode);

	if( rIA[0]!=0 || rIA[1]!=1 || rIA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}


{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	double VA[] = {2,1,0};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  P[] = {3,2,1,0,0};
	const rsb_nnz_idx_t nnz = 3;

	rsb__ip_reord(nnz, VA, IA, JA, P, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	float VA[] = {2,1,0};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  P[] = {3,2,1,0,0};
	const rsb_nnz_idx_t nnz = 3;

	rsb__ip_reord(nnz, VA, IA, JA, P, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	float complex VA[] = {2,1,0};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  P[] = {3,2,1,0,0};
	const rsb_nnz_idx_t nnz = 3;

	rsb__ip_reord(nnz, VA, IA, JA, P, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	double complex VA[] = {2,1,0};
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_coo_idx_t JA[] = {2,1,0};
	rsb_nnz_idx_t  P[] = {3,2,1,0,0};
	const rsb_nnz_idx_t nnz = 3;

	rsb__ip_reord(nnz, VA, IA, JA, P, typecode);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}


{
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_rows_with_coo_index(IA, K, nnz);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_rows_with_coo_index(IA, K, nnz);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_rows_with_coo_index(IA, K, nnz);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}
{
	rsb_coo_idx_t IA[] = {2,1,0};
	rsb_nnz_idx_t  K[] = {2,1,0};
	const rsb_nnz_idx_t nnz = 3;

	errval = rsb__do_permute_rows_with_coo_index(IA, K, nnz);

	if( IA[0]!=0 || IA[1]!=1 || IA[2]!=2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
}


err:
	return errval;
}

#if 0
void rsb__util_do_scatter_rows(void * RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT PA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode){
	/**
		This is an adapted PSBLAS psb_ip_reord_d1i2 routine.
	*/
	
	switch(typecode)
	{
			/* supported (double,float,float complex,double complex) */
	case RSB_NUMERICAL_TYPE_DOUBLE 	:
	{
		rsb_nnz_idx_t n;
		for(n=0;RSB_LIKELY(n<nnz);++n)
			((double*)oVA)[PA[IA[n]]]=((double*)VA)[n],
			oIA[PA[IA[n]]]=IA[n],
			oJA[PA[IA[n]]]=JA[n],
			PA[IA[n]]++;
		}
		break;
		case RSB_NUMERICAL_TYPE_FLOAT 	:
	{
		rsb_nnz_idx_t n;
		for(n=0;RSB_LIKELY(n<nnz);++n)
			((float*)oVA)[PA[IA[n]]]=((float*)VA)[n],
			oIA[PA[IA[n]]]=IA[n],
			oJA[PA[IA[n]]]=JA[n],
			PA[IA[n]]++;
		}
		break;
		case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:
	{
		rsb_nnz_idx_t n;
		for(n=0;RSB_LIKELY(n<nnz);++n)
			((float complex*)oVA)[PA[IA[n]]]=((float complex*)VA)[n],
			oIA[PA[IA[n]]]=IA[n],
			oJA[PA[IA[n]]]=JA[n],
			PA[IA[n]]++;
		}
		break;
		case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:
	{
		rsb_nnz_idx_t n;
		for(n=0;RSB_LIKELY(n<nnz);++n)
			((double complex*)oVA)[PA[IA[n]]]=((double complex*)VA)[n],
			oIA[PA[IA[n]]]=IA[n],
			oJA[PA[IA[n]]]=JA[n],
			PA[IA[n]]++;
		}
		break;
	
		/* unsupported type */
		default :
			return;
	}
}

#endif /* 0 */



#ifdef __cplusplus
}
#endif  /* __cplusplus */

/* @endcond */
