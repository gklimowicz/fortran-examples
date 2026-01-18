/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains some CSR sparse matrices multiplication code.
 * */

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
#include "rsb_internals.h"


rsb_err_t rsb__do_util_csr_csr_sparse_mul_serial(rsb_nnz_idx_t * PA, rsb_coo_idx_t * JA, void *VA_, const rsb_nnz_idx_t *ARP, const rsb_nnz_idx_t *BRP, const rsb_coo_idx_t *AJA, const rsb_coo_idx_t *BJA, const void * aVA_, const void * bVA_, const rsb_coo_idx_t cm, const rsb_coo_idx_t ck, rsb_nnz_idx_t * p, void * acc_, rsb_nnz_idx_t * opsp , rsb_type_t typecode, const rsb_coo_idx_t afr, const rsb_coo_idx_t ars)
{
	rsb_nnz_idx_t cblocks=0; 
	rsb_nnz_idx_t ops=0; 
	rsb_coo_idx_t ai,aj;
	rsb_coo_idx_t al,bl,cl;
	rsb_coo_idx_t bj;
	//rsb_coo_idx_t bi;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	double *VA=VA_,*acc=acc_;
	const double * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;
		rsb_nnz_idx_t marker;

		//assert(cblocks==PA[ai]);	// this is true on the serial execution of this loop
		cblocks=PA[ai];		// this shall work even in a parallel execution of this loop (with differing acc/p arrays)
		marker=cblocks+1;
		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];
				if(p[bj]<marker)
					p[bj]=marker,
					(JA)[cblocks++]=bj,
					acc[bj] =aVA[al]*bVA[bl];
				else
					acc[bj]+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}

		for(cl=(PA)[ai];cl<(PA)[ai+1];++cl)
		{
			((double*)(VA))[cl]=acc[(JA)[cl]];	/* FIXME */
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	float *VA=VA_,*acc=acc_;
	const float * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;
		rsb_nnz_idx_t marker;

		//assert(cblocks==PA[ai]);	// this is true on the serial execution of this loop
		cblocks=PA[ai];		// this shall work even in a parallel execution of this loop (with differing acc/p arrays)
		marker=cblocks+1;
		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];
				if(p[bj]<marker)
					p[bj]=marker,
					(JA)[cblocks++]=bj,
					acc[bj] =aVA[al]*bVA[bl];
				else
					acc[bj]+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}

		for(cl=(PA)[ai];cl<(PA)[ai+1];++cl)
		{
			((float*)(VA))[cl]=acc[(JA)[cl]];	/* FIXME */
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	float complex *VA=VA_,*acc=acc_;
	const float complex * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;
		rsb_nnz_idx_t marker;

		//assert(cblocks==PA[ai]);	// this is true on the serial execution of this loop
		cblocks=PA[ai];		// this shall work even in a parallel execution of this loop (with differing acc/p arrays)
		marker=cblocks+1;
		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];
				if(p[bj]<marker)
					p[bj]=marker,
					(JA)[cblocks++]=bj,
					acc[bj] =aVA[al]*bVA[bl];
				else
					acc[bj]+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}

		for(cl=(PA)[ai];cl<(PA)[ai+1];++cl)
		{
			((float complex*)(VA))[cl]=acc[(JA)[cl]];	/* FIXME */
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	double complex *VA=VA_,*acc=acc_;
	const double complex * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;
		rsb_nnz_idx_t marker;

		//assert(cblocks==PA[ai]);	// this is true on the serial execution of this loop
		cblocks=PA[ai];		// this shall work even in a parallel execution of this loop (with differing acc/p arrays)
		marker=cblocks+1;
		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];
				if(p[bj]<marker)
					p[bj]=marker,
					(JA)[cblocks++]=bj,
					acc[bj] =aVA[al]*bVA[bl];
				else
					acc[bj]+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}

		for(cl=(PA)[ai];cl<(PA)[ai+1];++cl)
		{
			((double complex*)(VA))[cl]=acc[(JA)[cl]];	/* FIXME */
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__do_util_csr_csr_dense_mul_serial(rsb_coo_idx_t ldc, rsb_coo_idx_t nr, rsb_coo_idx_t nc, rsb_bool_t isccolmajor, void *cVA_, const rsb_nnz_idx_t *ARP, const rsb_nnz_idx_t *BRP, const rsb_coo_idx_t *AJA, const rsb_coo_idx_t *BJA, const void * aVA_, const void * bVA_, const rsb_coo_idx_t cm, const rsb_coo_idx_t ck, rsb_nnz_idx_t * opsp , rsb_type_t typecode, const rsb_coo_idx_t afr, const rsb_coo_idx_t ars)
{
	rsb_nnz_idx_t ops=0; 
	rsb_coo_idx_t ai,aj;
	rsb_coo_idx_t al,bl;
	//rsb_coo_idx_t bi;
	rsb_coo_idx_t bj;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	double *cVA=cVA_;
	const double * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;

		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];

				RSB_BLOCK_X_MAJOR_REFERENCE(cVA,ldc,ai,bj,isccolmajor)+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	float *cVA=cVA_;
	const float * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;

		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];

				RSB_BLOCK_X_MAJOR_REFERENCE(cVA,ldc,ai,bj,isccolmajor)+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	float complex *cVA=cVA_;
	const float complex * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;

		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];

				RSB_BLOCK_X_MAJOR_REFERENCE(cVA,ldc,ai,bj,isccolmajor)+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	double complex *cVA=cVA_;
	const double complex * aVA=aVA_,*bVA=bVA_; ;
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;

		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];

				RSB_BLOCK_X_MAJOR_REFERENCE(cVA,ldc,ai,bj,isccolmajor)+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}


/* @endcond */
