dnl
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
dnl
/* @cond INNERDOC */
dnl
/**
 * @file
 * @brief
 * Permutation functions.
 */
RSB_M4_HEADER_MESSAGE()dnl


dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_PERMUTE_H_INCLUDED
#define RSB_PERMUTE_H_INCLUDED
')
dnl
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
dnl
#include "rsb_common.h"
dnl 

rsb_err_t rsb__do_permute_values_in_place_with_coo_index(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
		This is in place swapping, and it is much slower than not in place.
		
		FIXME : document 

		Sadly, this is O(nnz^2).
		... Or qsorts order ?
	 */
	rsb_coo_idx_t n;/* this is the case where coo cannot overflow */

	switch(type)
	{
		/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
	case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:

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
		RSB_SWAP(type,((type*)VA)[n],((type*)VA)[t]);
	}
			break;
')dnl
			/* unsupported type */
		default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

rsb_err_t rsb__do_permute_values_in_place_with_nnz_index(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
		/*	
			This is in place swapping, and it is much slower than not in place.
			
			FIXME : document and finish (s/double/ * /).

			Sadly, this is O(nnz^2).
			... Or qsorts order ?
		 */
	rsb_nnz_idx_t n;

	switch(type)
	{
		/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
	case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:

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
		RSB_SWAP(type,((type*)VA)[n],((type*)VA)[t]);
	}
			break;
')dnl
			/* unsupported type */
		default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

rsb_err_t rsb__do_permute_values_with_coo_index( void * rVA, const void *VA, rsb_coo_idx_t * rIA, const rsb_coo_idx_t * IA, rsb_coo_idx_t * rJA, const rsb_coo_idx_t * JA, const rsb_coo_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t type)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
	{
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
			/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`nnz',` ((type*)rVA)[i+LI]=((type*)VA)[K[(i+LI)]];
	')
			
			break;
')dnl
			/* unsupported type */
			default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
		}
		return RSB_ERR_NO_ERROR;
	}
')dnl


rsb_err_t rsb__do_permute_rows_with_coo_index( rsb_coo_idx_t * IA, const rsb_coo_idx_t * K, rsb_nnz_idx_t nnz)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
	{
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
')dnl


rsb_err_t rsb__do_permute_values_with_nnz_index( void * rVA, const void *VA, rsb_coo_idx_t * rIA, const rsb_coo_idx_t * IA, rsb_coo_idx_t * rJA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t * K, rsb_nnz_idx_t nnz, rsb_type_t typecode)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
	{
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
			/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`nnz',` ((type*)rVA)[i+LI]=((type*)VA)[K[(i+LI)]];
	')
			
			break;
')dnl
			/* unsupported type */
			default :
				return RSB_ERR_UNSUPPORTED_TYPE	;
		}
		return RSB_ERR_NO_ERROR;
	}
')dnl


void rsb__ip_reord(rsb_nnz_idx_t n, void * VAp, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * P, rsb_type_t typecode)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
		This is an adapted PSBLAS psb_ip_reord_d1i2 routine.
	*/
	
	switch(typecode)
	{
			/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
	case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:
	{
		rsb_coo_idx_t isw1, isw2;
		rsb_nnz_idx_t lswap, lp, k;
		type swap;
		type * VA=VAp;

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
	')
		/* unsupported type */
		default :
			return;
	}


}
')dnl

void rsb__util_do_scatter_rows(void * RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, const void * RSB_RESTRICT iVA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT PA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
		This is an adapted PSBLAS psb_ip_reord_d1i2 routine.
	*/
	
	switch(typecode)
	{
			/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
	case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:
	{
		rsb_nnz_idx_t nzi;
		type*VA=(type*)oVA;
		RSB_M4_SIMPLE_LOOP_UNROLL(`nzi',`LI',`0',`nnz',` 
			VA[PA[IA[nzi+LI]]]=((type*)iVA)[nzi+LI];
			oIA[PA[IA[nzi+LI]]]=IA[nzi+LI];
			oJA[PA[IA[nzi+LI]]]=JA[nzi+LI];
			PA[IA[nzi+LI]]++;
		')
	}
		break;
	')
		/* unsupported type */
		default :
			return;
	}
}
')dnl

rsb_err_t rsb__chk_permute(void)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
	
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	type VA[] = {1,1,1};
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
')dnl

foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	type VA[] = {1,1,1};
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
')dnl


foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	type VA[] = {1,1,1};
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
')dnl

foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	type VA[] = {1,1,1};
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
')dnl


foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	const type VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	type rVA[] = {1,1,1};
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
')dnl


foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	const type VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_coo_idx_t  K[] = {2,1,0};
	type rVA[] = {1,1,1};
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
')dnl


foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	const type VA[] = {2,1,0};
	const rsb_coo_idx_t IA[] = {2,1,0};
	const rsb_coo_idx_t JA[] = {2,1,0};
	const rsb_nnz_idx_t  K[] = {2,1,0};
	type rVA[] = {1,1,1};
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
')dnl


foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type);
	type VA[] = {2,1,0};
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
')dnl


foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
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
')dnl


err:
	return errval;
}
')dnl

#if 0
void rsb__util_do_scatter_rows(void * RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT PA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
		This is an adapted PSBLAS psb_ip_reord_d1i2 routine.
	*/
	
	switch(typecode)
	{
			/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
	case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:
	{
		rsb_nnz_idx_t n;
		for(n=0;RSB_LIKELY(n<nnz);++n)
			((type*)oVA)[PA[IA[n]]]=((type*)VA)[n],
			oIA[PA[IA[n]]]=IA[n],
			oJA[PA[IA[n]]]=JA[n],
			PA[IA[n]]++;
		}
		break;
	')
		/* unsupported type */
		default :
			return;
	}
}
')dnl

#endif /* 0 */



dnl
#ifdef __cplusplus
}
#endif  /* __cplusplus */
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_PERMUTE_H_INCLUDED */
')
dnl
/* @endcond */
dnl
