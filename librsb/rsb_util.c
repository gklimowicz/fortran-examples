/* @cond INNERDOC */


/**
 * @file
 * @brief
 * Auxiliary functions.
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
/* non blas-like functions */

rsb_err_t rsb__util_m4_sanity_check(void){
	/**
		There are bugs in the m4 macros or a bad m4 implementation which will trigger this test to fail.
		We are interested in catching them, as we should rely on a sane m4 environment.
	*/
	
	if(
		0!=0 ||
		1!=1 || 
		1!=1 || 
		0!=0 ||
		0!=0 ||
		0!=0 ||
		0!=0 ||
		1!=1 ||
		0!=0 ||
		1!=1 ||
		1!=1 ||
		1!=1 ||
		0
		)
		goto err;
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_INTERNAL_ERROR;
}

const void * rsb__util_increase_by_one(void *p, rsb_nnz_idx_t n, rsb_type_t typecode){
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  ) {(((double*)p)[n])+=1;return p;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  ) {(((float*)p)[n])+=1;return p;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ) {(((float complex*)p)[n])+=1;return p;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ) {(((double complex*)p)[n])+=1;return p;}
	else 
#endif
	return NULL;
}

void rsb__util_set_area_to_fraction_of_integer(void *p, const int alphai, rsb_type_t typecode){
	/*
		alpha NULL will imply 1
	*/
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  ) {*(double*)p = 1;*(double*)p/=alphai;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  ) {*(float*)p = 1;*(float*)p/=alphai;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ) {*(float complex*)p = 1;*(float complex*)p/=alphai;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ) {*(double complex*)p = 1;*(double complex*)p/=alphai;}
	else 
#endif
	return;
}

void rsb__util_set_area_to_negated_fraction(void *p, const void *alpha, rsb_type_t typecode){
	/*
		alpha NULL will imply 1
	*/
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  ) {*(double*)p = -1;if(alpha)*(double*)p/=(*(double*)alpha);}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  ) {*(float*)p = -1;if(alpha)*(float*)p/=(*(float*)alpha);}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ) {*(float complex*)p = -1;if(alpha)*(float complex*)p/=(*(float complex*)alpha);}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ) {*(double complex*)p = -1;if(alpha)*(double complex*)p/=(*(double complex*)alpha);}
	else 
#endif
	return;
}

void rsb__util_set_area_to_converted_integer(void *p, rsb_type_t typecode, const rsb_int n){
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  ) {*(double*)p = (double)n;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  ) {*(float*)p = (float)n;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ) {*(float complex*)p = (float complex)n;}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ) {*(double complex*)p = (double complex)n;}
	else 
#endif
	return;
}

rsb_coo_idx_t * rsb__util_get_partitioning_array( size_t bs, size_t X , rsb_blk_idx_t * X_b, rsb_flags_t flags){
	/*!
	 * Given a block size (be it rows or columns), an element size X in bytes,
	 * and a dimension (rows or columns), returns an array containing the 
	 * indices of the elements in each block.
	 *
	 * Therefore, the allocated arrays 
	 *
	 * \param bs	the block size
	 * \param X	the rows or columns count
	 * \param X_b	on output, the allocated array elements count : (X+bs-1)/bs
	 * \return NULL on error;  a valid array pointer on success
	 *
	 * FIXME : why not size_t ? or maybe rsb_size_t ?
	 * */
	size_t i;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t * p_x = NULL;

	*X_b = (X+bs-1)/bs;

	/* WARNING : 1 is the extreme limit before overflow :) */
	if( ( ((size_t)(*X_b)) < ((size_t)((X+bs-1)/bs))) || (RSB_BLK_ADD_OVERFLOW(*X_b,1)) )
	{
		/* overflow. should print some message. */
		errval = RSB_ERR_LIMITS;goto err;
	}

	p_x = rsb__malloc(sizeof(rsb_coo_idx_t)*(*X_b+1));
	if(! p_x) goto err;
	/* note: should use some perrno some day */

	/* note the last block size : it is the same, regardless congruences */
	{
for(i=0;i+15<*X_b;i+=16){
p_x[i+0 ] = (i+0 )*bs;
	p_x[i+1 ] = (i+1 )*bs;
	p_x[i+2 ] = (i+2 )*bs;
	p_x[i+3 ] = (i+3 )*bs;
	p_x[i+4 ] = (i+4 )*bs;
	p_x[i+5 ] = (i+5 )*bs;
	p_x[i+6 ] = (i+6 )*bs;
	p_x[i+7 ] = (i+7 )*bs;
	p_x[i+8 ] = (i+8 )*bs;
	p_x[i+9 ] = (i+9 )*bs;
	p_x[i+10 ] = (i+10 )*bs;
	p_x[i+11 ] = (i+11 )*bs;
	p_x[i+12 ] = (i+12 )*bs;
	p_x[i+13 ] = (i+13 )*bs;
	p_x[i+14 ] = (i+14 )*bs;
	p_x[i+15 ] = (i+15 )*bs;
	}
for(     ;i<*X_b;++i){ p_x[i+0 ] = (i+0 )*bs;
	 }
}


	/* FIXME : this point should be remarked and documented way better ! */
	if(flags&(RSB_FLAG_WANT_BCSS_STORAGE|RSB_FLAG_WANT_FIXED_BLOCKING_VBR))
		p_x[*X_b] = *X_b*bs;	/* the last element of p_x is the index of the last matrix row/column    + 1  */
	else
		p_x[*X_b] = X;	/* the last element of p_x is the index of the last matrix row/column    + 1  */
	
	return p_x;
err:
	RSB_CONDITIONAL_FREE(p_x);
	rsb__do_perror(NULL,errval);
	return NULL;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__vector_diff(void * c, const void * a, const void * b, rsb_type_t type, size_t n){
	/*!
	 * c <- a-b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see daxpy,dcopy in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double*ta = a,*tb = b;double *tc = c;
		{
for(i=0;i+15<n;i+=16){
		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
			tc[i+1 ] = ta[i+1 ]-tb[i+1 ];
			tc[i+2 ] = ta[i+2 ]-tb[i+2 ];
			tc[i+3 ] = ta[i+3 ]-tb[i+3 ];
			tc[i+4 ] = ta[i+4 ]-tb[i+4 ];
			tc[i+5 ] = ta[i+5 ]-tb[i+5 ];
			tc[i+6 ] = ta[i+6 ]-tb[i+6 ];
			tc[i+7 ] = ta[i+7 ]-tb[i+7 ];
			tc[i+8 ] = ta[i+8 ]-tb[i+8 ];
			tc[i+9 ] = ta[i+9 ]-tb[i+9 ];
			tc[i+10 ] = ta[i+10 ]-tb[i+10 ];
			tc[i+11 ] = ta[i+11 ]-tb[i+11 ];
			tc[i+12 ] = ta[i+12 ]-tb[i+12 ];
			tc[i+13 ] = ta[i+13 ]-tb[i+13 ];
			tc[i+14 ] = ta[i+14 ]-tb[i+14 ];
			tc[i+15 ] = ta[i+15 ]-tb[i+15 ];
	}
for(     ;i<n;++i){ 		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float*ta = a,*tb = b;float *tc = c;
		{
for(i=0;i+15<n;i+=16){
		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
			tc[i+1 ] = ta[i+1 ]-tb[i+1 ];
			tc[i+2 ] = ta[i+2 ]-tb[i+2 ];
			tc[i+3 ] = ta[i+3 ]-tb[i+3 ];
			tc[i+4 ] = ta[i+4 ]-tb[i+4 ];
			tc[i+5 ] = ta[i+5 ]-tb[i+5 ];
			tc[i+6 ] = ta[i+6 ]-tb[i+6 ];
			tc[i+7 ] = ta[i+7 ]-tb[i+7 ];
			tc[i+8 ] = ta[i+8 ]-tb[i+8 ];
			tc[i+9 ] = ta[i+9 ]-tb[i+9 ];
			tc[i+10 ] = ta[i+10 ]-tb[i+10 ];
			tc[i+11 ] = ta[i+11 ]-tb[i+11 ];
			tc[i+12 ] = ta[i+12 ]-tb[i+12 ];
			tc[i+13 ] = ta[i+13 ]-tb[i+13 ];
			tc[i+14 ] = ta[i+14 ]-tb[i+14 ];
			tc[i+15 ] = ta[i+15 ]-tb[i+15 ];
	}
for(     ;i<n;++i){ 		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex*ta = a,*tb = b;float complex *tc = c;
		{
for(i=0;i+15<n;i+=16){
		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
			tc[i+1 ] = ta[i+1 ]-tb[i+1 ];
			tc[i+2 ] = ta[i+2 ]-tb[i+2 ];
			tc[i+3 ] = ta[i+3 ]-tb[i+3 ];
			tc[i+4 ] = ta[i+4 ]-tb[i+4 ];
			tc[i+5 ] = ta[i+5 ]-tb[i+5 ];
			tc[i+6 ] = ta[i+6 ]-tb[i+6 ];
			tc[i+7 ] = ta[i+7 ]-tb[i+7 ];
			tc[i+8 ] = ta[i+8 ]-tb[i+8 ];
			tc[i+9 ] = ta[i+9 ]-tb[i+9 ];
			tc[i+10 ] = ta[i+10 ]-tb[i+10 ];
			tc[i+11 ] = ta[i+11 ]-tb[i+11 ];
			tc[i+12 ] = ta[i+12 ]-tb[i+12 ];
			tc[i+13 ] = ta[i+13 ]-tb[i+13 ];
			tc[i+14 ] = ta[i+14 ]-tb[i+14 ];
			tc[i+15 ] = ta[i+15 ]-tb[i+15 ];
	}
for(     ;i<n;++i){ 		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex*ta = a,*tb = b;double complex *tc = c;
		{
for(i=0;i+15<n;i+=16){
		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
			tc[i+1 ] = ta[i+1 ]-tb[i+1 ];
			tc[i+2 ] = ta[i+2 ]-tb[i+2 ];
			tc[i+3 ] = ta[i+3 ]-tb[i+3 ];
			tc[i+4 ] = ta[i+4 ]-tb[i+4 ];
			tc[i+5 ] = ta[i+5 ]-tb[i+5 ];
			tc[i+6 ] = ta[i+6 ]-tb[i+6 ];
			tc[i+7 ] = ta[i+7 ]-tb[i+7 ];
			tc[i+8 ] = ta[i+8 ]-tb[i+8 ];
			tc[i+9 ] = ta[i+9 ]-tb[i+9 ];
			tc[i+10 ] = ta[i+10 ]-tb[i+10 ];
			tc[i+11 ] = ta[i+11 ]-tb[i+11 ];
			tc[i+12 ] = ta[i+12 ]-tb[i+12 ];
			tc[i+13 ] = ta[i+13 ]-tb[i+13 ];
			tc[i+14 ] = ta[i+14 ]-tb[i+14 ];
			tc[i+15 ] = ta[i+15 ]-tb[i+15 ];
	}
for(     ;i<n;++i){ 		tc[i+0 ] = ta[i+0 ]-tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

static rsb_err_t rsb_vector_norm_square(void * c, const void * a, rsb_type_t type, size_t n)
{
	/*!
	 * c <- a^T*a
         *
	 * \param a	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double*ta = a;double *tc = c;
		tc[0] = ((double)(0));
		{
for(i=0;i+15<n;i+=16){
		tc[0]+=ta[i+0 ]*ta[i+0 ];
			tc[0]+=ta[i+1 ]*ta[i+1 ];
			tc[0]+=ta[i+2 ]*ta[i+2 ];
			tc[0]+=ta[i+3 ]*ta[i+3 ];
			tc[0]+=ta[i+4 ]*ta[i+4 ];
			tc[0]+=ta[i+5 ]*ta[i+5 ];
			tc[0]+=ta[i+6 ]*ta[i+6 ];
			tc[0]+=ta[i+7 ]*ta[i+7 ];
			tc[0]+=ta[i+8 ]*ta[i+8 ];
			tc[0]+=ta[i+9 ]*ta[i+9 ];
			tc[0]+=ta[i+10 ]*ta[i+10 ];
			tc[0]+=ta[i+11 ]*ta[i+11 ];
			tc[0]+=ta[i+12 ]*ta[i+12 ];
			tc[0]+=ta[i+13 ]*ta[i+13 ];
			tc[0]+=ta[i+14 ]*ta[i+14 ];
			tc[0]+=ta[i+15 ]*ta[i+15 ];
	}
for(     ;i<n;++i){ 		tc[0]+=ta[i+0 ]*ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float*ta = a;float *tc = c;
		tc[0] = ((float)(0));
		{
for(i=0;i+15<n;i+=16){
		tc[0]+=ta[i+0 ]*ta[i+0 ];
			tc[0]+=ta[i+1 ]*ta[i+1 ];
			tc[0]+=ta[i+2 ]*ta[i+2 ];
			tc[0]+=ta[i+3 ]*ta[i+3 ];
			tc[0]+=ta[i+4 ]*ta[i+4 ];
			tc[0]+=ta[i+5 ]*ta[i+5 ];
			tc[0]+=ta[i+6 ]*ta[i+6 ];
			tc[0]+=ta[i+7 ]*ta[i+7 ];
			tc[0]+=ta[i+8 ]*ta[i+8 ];
			tc[0]+=ta[i+9 ]*ta[i+9 ];
			tc[0]+=ta[i+10 ]*ta[i+10 ];
			tc[0]+=ta[i+11 ]*ta[i+11 ];
			tc[0]+=ta[i+12 ]*ta[i+12 ];
			tc[0]+=ta[i+13 ]*ta[i+13 ];
			tc[0]+=ta[i+14 ]*ta[i+14 ];
			tc[0]+=ta[i+15 ]*ta[i+15 ];
	}
for(     ;i<n;++i){ 		tc[0]+=ta[i+0 ]*ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex*ta = a;float complex *tc = c;
		tc[0] = ((float complex)(0));
		{
for(i=0;i+15<n;i+=16){
		tc[0]+=ta[i+0 ]*ta[i+0 ];
			tc[0]+=ta[i+1 ]*ta[i+1 ];
			tc[0]+=ta[i+2 ]*ta[i+2 ];
			tc[0]+=ta[i+3 ]*ta[i+3 ];
			tc[0]+=ta[i+4 ]*ta[i+4 ];
			tc[0]+=ta[i+5 ]*ta[i+5 ];
			tc[0]+=ta[i+6 ]*ta[i+6 ];
			tc[0]+=ta[i+7 ]*ta[i+7 ];
			tc[0]+=ta[i+8 ]*ta[i+8 ];
			tc[0]+=ta[i+9 ]*ta[i+9 ];
			tc[0]+=ta[i+10 ]*ta[i+10 ];
			tc[0]+=ta[i+11 ]*ta[i+11 ];
			tc[0]+=ta[i+12 ]*ta[i+12 ];
			tc[0]+=ta[i+13 ]*ta[i+13 ];
			tc[0]+=ta[i+14 ]*ta[i+14 ];
			tc[0]+=ta[i+15 ]*ta[i+15 ];
	}
for(     ;i<n;++i){ 		tc[0]+=ta[i+0 ]*ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex*ta = a;double complex *tc = c;
		tc[0] = ((double complex)(0));
		{
for(i=0;i+15<n;i+=16){
		tc[0]+=ta[i+0 ]*ta[i+0 ];
			tc[0]+=ta[i+1 ]*ta[i+1 ];
			tc[0]+=ta[i+2 ]*ta[i+2 ];
			tc[0]+=ta[i+3 ]*ta[i+3 ];
			tc[0]+=ta[i+4 ]*ta[i+4 ];
			tc[0]+=ta[i+5 ]*ta[i+5 ];
			tc[0]+=ta[i+6 ]*ta[i+6 ];
			tc[0]+=ta[i+7 ]*ta[i+7 ];
			tc[0]+=ta[i+8 ]*ta[i+8 ];
			tc[0]+=ta[i+9 ]*ta[i+9 ];
			tc[0]+=ta[i+10 ]*ta[i+10 ];
			tc[0]+=ta[i+11 ]*ta[i+11 ];
			tc[0]+=ta[i+12 ]*ta[i+12 ];
			tc[0]+=ta[i+13 ]*ta[i+13 ];
			tc[0]+=ta[i+14 ]*ta[i+14 ];
			tc[0]+=ta[i+15 ]*ta[i+15 ];
	}
for(     ;i<n;++i){ 		tc[0]+=ta[i+0 ]*ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

static rsb_err_t rsb_vector_norm(void * c, const void * a, rsb_type_t type, size_t n)
{
	/*!
	 * c <- sqrt(a^T*a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	rsb_err_t errval;
	if(!c)
		return RSB_ERR_BADARGS;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		double*cp = (double*)c;
		errval = rsb_vector_norm_square(cp,a,type,n);
		*cp = sqrt(*cp);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		float*cp = (float*)c;
		errval = rsb_vector_norm_square(cp,a,type,n);
		*cp = sqrtf(*cp);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		float complex*cp = (float complex*)c;
		errval = rsb_vector_norm_square(cp,a,type,n);
		*cp = csqrtf(*cp);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		double complex*cp = (double complex*)c;
		errval = rsb_vector_norm_square(cp,a,type,n);
		*cp = csqrt(*cp);
	}
	else 
#endif
		errval = RSB_ERR_UNSUPPORTED_TYPE;
	RSB_DO_ERR_RETURN(errval)
}


rsb_err_t rsb__vector_norm_strided(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
{
	/*!
	 * c <- sqrt(a^T*a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	rsb_err_t errval;
	if(!c)
		return RSB_ERR_BADARGS;
	if(inc==1)
		return rsb_vector_norm(c,a,type,n);

	errval = RSB_ERR_INTERNAL_ERROR;
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__util_vector_sum_strided(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
{
	/*!
	 * c <- sum(a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	/* See also rsb__cblas_Xnrm2 */

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		register double acc = ((double)(0)); const double*ta = a; double*tc = c;
	{
for(i=0;i+15<n;i+=16){
		acc+=ta[(i+0 )*inc];
			acc+=ta[(i+1 )*inc];
			acc+=ta[(i+2 )*inc];
			acc+=ta[(i+3 )*inc];
			acc+=ta[(i+4 )*inc];
			acc+=ta[(i+5 )*inc];
			acc+=ta[(i+6 )*inc];
			acc+=ta[(i+7 )*inc];
			acc+=ta[(i+8 )*inc];
			acc+=ta[(i+9 )*inc];
			acc+=ta[(i+10 )*inc];
			acc+=ta[(i+11 )*inc];
			acc+=ta[(i+12 )*inc];
			acc+=ta[(i+13 )*inc];
			acc+=ta[(i+14 )*inc];
			acc+=ta[(i+15 )*inc];
	}
for(     ;i<n;++i){ 		acc+=ta[(i+0 )*inc];
	 }
}
; 
		tc[0] = acc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		register float acc = ((float)(0)); const float*ta = a; float*tc = c;
	{
for(i=0;i+15<n;i+=16){
		acc+=ta[(i+0 )*inc];
			acc+=ta[(i+1 )*inc];
			acc+=ta[(i+2 )*inc];
			acc+=ta[(i+3 )*inc];
			acc+=ta[(i+4 )*inc];
			acc+=ta[(i+5 )*inc];
			acc+=ta[(i+6 )*inc];
			acc+=ta[(i+7 )*inc];
			acc+=ta[(i+8 )*inc];
			acc+=ta[(i+9 )*inc];
			acc+=ta[(i+10 )*inc];
			acc+=ta[(i+11 )*inc];
			acc+=ta[(i+12 )*inc];
			acc+=ta[(i+13 )*inc];
			acc+=ta[(i+14 )*inc];
			acc+=ta[(i+15 )*inc];
	}
for(     ;i<n;++i){ 		acc+=ta[(i+0 )*inc];
	 }
}
; 
		tc[0] = acc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		register float complex acc = ((float complex)(0)); const float complex*ta = a; float complex*tc = c;
	{
for(i=0;i+15<n;i+=16){
		acc+=ta[(i+0 )*inc];
			acc+=ta[(i+1 )*inc];
			acc+=ta[(i+2 )*inc];
			acc+=ta[(i+3 )*inc];
			acc+=ta[(i+4 )*inc];
			acc+=ta[(i+5 )*inc];
			acc+=ta[(i+6 )*inc];
			acc+=ta[(i+7 )*inc];
			acc+=ta[(i+8 )*inc];
			acc+=ta[(i+9 )*inc];
			acc+=ta[(i+10 )*inc];
			acc+=ta[(i+11 )*inc];
			acc+=ta[(i+12 )*inc];
			acc+=ta[(i+13 )*inc];
			acc+=ta[(i+14 )*inc];
			acc+=ta[(i+15 )*inc];
	}
for(     ;i<n;++i){ 		acc+=ta[(i+0 )*inc];
	 }
}
; 
		tc[0] = acc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		register double complex acc = ((double complex)(0)); const double complex*ta = a; double complex*tc = c;
	{
for(i=0;i+15<n;i+=16){
		acc+=ta[(i+0 )*inc];
			acc+=ta[(i+1 )*inc];
			acc+=ta[(i+2 )*inc];
			acc+=ta[(i+3 )*inc];
			acc+=ta[(i+4 )*inc];
			acc+=ta[(i+5 )*inc];
			acc+=ta[(i+6 )*inc];
			acc+=ta[(i+7 )*inc];
			acc+=ta[(i+8 )*inc];
			acc+=ta[(i+9 )*inc];
			acc+=ta[(i+10 )*inc];
			acc+=ta[(i+11 )*inc];
			acc+=ta[(i+12 )*inc];
			acc+=ta[(i+13 )*inc];
			acc+=ta[(i+14 )*inc];
			acc+=ta[(i+15 )*inc];
	}
for(     ;i<n;++i){ 		acc+=ta[(i+0 )*inc];
	 }
}
; 
		tc[0] = acc;
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_vector_sum(void * c, const void * a, rsb_type_t type, size_t n)
{
	/*!
	 * c <- sum(a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double*ta = a; double*tc = c; tc[0] = ((double)(0));
	{
for(i=0;i+15<n;i+=16){
	tc[0]+=ta[i+0 ];
		tc[0]+=ta[i+1 ];
		tc[0]+=ta[i+2 ];
		tc[0]+=ta[i+3 ];
		tc[0]+=ta[i+4 ];
		tc[0]+=ta[i+5 ];
		tc[0]+=ta[i+6 ];
		tc[0]+=ta[i+7 ];
		tc[0]+=ta[i+8 ];
		tc[0]+=ta[i+9 ];
		tc[0]+=ta[i+10 ];
		tc[0]+=ta[i+11 ];
		tc[0]+=ta[i+12 ];
		tc[0]+=ta[i+13 ];
		tc[0]+=ta[i+14 ];
		tc[0]+=ta[i+15 ];
	}
for(     ;i<n;++i){ 	tc[0]+=ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float*ta = a; float*tc = c; tc[0] = ((float)(0));
	{
for(i=0;i+15<n;i+=16){
	tc[0]+=ta[i+0 ];
		tc[0]+=ta[i+1 ];
		tc[0]+=ta[i+2 ];
		tc[0]+=ta[i+3 ];
		tc[0]+=ta[i+4 ];
		tc[0]+=ta[i+5 ];
		tc[0]+=ta[i+6 ];
		tc[0]+=ta[i+7 ];
		tc[0]+=ta[i+8 ];
		tc[0]+=ta[i+9 ];
		tc[0]+=ta[i+10 ];
		tc[0]+=ta[i+11 ];
		tc[0]+=ta[i+12 ];
		tc[0]+=ta[i+13 ];
		tc[0]+=ta[i+14 ];
		tc[0]+=ta[i+15 ];
	}
for(     ;i<n;++i){ 	tc[0]+=ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex*ta = a; float complex*tc = c; tc[0] = ((float complex)(0));
	{
for(i=0;i+15<n;i+=16){
	tc[0]+=ta[i+0 ];
		tc[0]+=ta[i+1 ];
		tc[0]+=ta[i+2 ];
		tc[0]+=ta[i+3 ];
		tc[0]+=ta[i+4 ];
		tc[0]+=ta[i+5 ];
		tc[0]+=ta[i+6 ];
		tc[0]+=ta[i+7 ];
		tc[0]+=ta[i+8 ];
		tc[0]+=ta[i+9 ];
		tc[0]+=ta[i+10 ];
		tc[0]+=ta[i+11 ];
		tc[0]+=ta[i+12 ];
		tc[0]+=ta[i+13 ];
		tc[0]+=ta[i+14 ];
		tc[0]+=ta[i+15 ];
	}
for(     ;i<n;++i){ 	tc[0]+=ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex*ta = a; double complex*tc = c; tc[0] = ((double complex)(0));
	{
for(i=0;i+15<n;i+=16){
	tc[0]+=ta[i+0 ];
		tc[0]+=ta[i+1 ];
		tc[0]+=ta[i+2 ];
		tc[0]+=ta[i+3 ];
		tc[0]+=ta[i+4 ];
		tc[0]+=ta[i+5 ];
		tc[0]+=ta[i+6 ];
		tc[0]+=ta[i+7 ];
		tc[0]+=ta[i+8 ];
		tc[0]+=ta[i+9 ];
		tc[0]+=ta[i+10 ];
		tc[0]+=ta[i+11 ];
		tc[0]+=ta[i+12 ];
		tc[0]+=ta[i+13 ];
		tc[0]+=ta[i+14 ];
		tc[0]+=ta[i+15 ];
	}
for(     ;i<n;++i){ 	tc[0]+=ta[i+0 ];
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_err_t rsb__vector_mult_sum(const void * a, const void * b, void * c, rsb_type_t type, size_t n, const int inca, const int incb)
{
	/*!
	 * c <- sum(a*b)
	 * It is allowed to give c == a or c == b or a==b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * 
	 * p.s.: this routine is, numerically speaking, a crime!
	 * 
	 * */
	size_t i;
	if(a==b && inca==incb)
		return rsb_vector_norm_square_strided(c,a,type,n,inca);
	if(inca == 1 && incb == 1)
	{
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double*tb = b; const double*ta = a; double*tc = c,cacc = ((double)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[i+0 ]*tb[i+0 ];
			cacc+=ta[i+1 ]*tb[i+1 ];
			cacc+=ta[i+2 ]*tb[i+2 ];
			cacc+=ta[i+3 ]*tb[i+3 ];
			cacc+=ta[i+4 ]*tb[i+4 ];
			cacc+=ta[i+5 ]*tb[i+5 ];
			cacc+=ta[i+6 ]*tb[i+6 ];
			cacc+=ta[i+7 ]*tb[i+7 ];
			cacc+=ta[i+8 ]*tb[i+8 ];
			cacc+=ta[i+9 ]*tb[i+9 ];
			cacc+=ta[i+10 ]*tb[i+10 ];
			cacc+=ta[i+11 ]*tb[i+11 ];
			cacc+=ta[i+12 ]*tb[i+12 ];
			cacc+=ta[i+13 ]*tb[i+13 ];
			cacc+=ta[i+14 ]*tb[i+14 ];
			cacc+=ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 		cacc+=ta[i+0 ]*tb[i+0 ];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float*tb = b; const float*ta = a; float*tc = c,cacc = ((float)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[i+0 ]*tb[i+0 ];
			cacc+=ta[i+1 ]*tb[i+1 ];
			cacc+=ta[i+2 ]*tb[i+2 ];
			cacc+=ta[i+3 ]*tb[i+3 ];
			cacc+=ta[i+4 ]*tb[i+4 ];
			cacc+=ta[i+5 ]*tb[i+5 ];
			cacc+=ta[i+6 ]*tb[i+6 ];
			cacc+=ta[i+7 ]*tb[i+7 ];
			cacc+=ta[i+8 ]*tb[i+8 ];
			cacc+=ta[i+9 ]*tb[i+9 ];
			cacc+=ta[i+10 ]*tb[i+10 ];
			cacc+=ta[i+11 ]*tb[i+11 ];
			cacc+=ta[i+12 ]*tb[i+12 ];
			cacc+=ta[i+13 ]*tb[i+13 ];
			cacc+=ta[i+14 ]*tb[i+14 ];
			cacc+=ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 		cacc+=ta[i+0 ]*tb[i+0 ];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex*tb = b; const float complex*ta = a; float complex*tc = c,cacc = ((float complex)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[i+0 ]*tb[i+0 ];
			cacc+=ta[i+1 ]*tb[i+1 ];
			cacc+=ta[i+2 ]*tb[i+2 ];
			cacc+=ta[i+3 ]*tb[i+3 ];
			cacc+=ta[i+4 ]*tb[i+4 ];
			cacc+=ta[i+5 ]*tb[i+5 ];
			cacc+=ta[i+6 ]*tb[i+6 ];
			cacc+=ta[i+7 ]*tb[i+7 ];
			cacc+=ta[i+8 ]*tb[i+8 ];
			cacc+=ta[i+9 ]*tb[i+9 ];
			cacc+=ta[i+10 ]*tb[i+10 ];
			cacc+=ta[i+11 ]*tb[i+11 ];
			cacc+=ta[i+12 ]*tb[i+12 ];
			cacc+=ta[i+13 ]*tb[i+13 ];
			cacc+=ta[i+14 ]*tb[i+14 ];
			cacc+=ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 		cacc+=ta[i+0 ]*tb[i+0 ];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex*tb = b; const double complex*ta = a; double complex*tc = c,cacc = ((double complex)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[i+0 ]*tb[i+0 ];
			cacc+=ta[i+1 ]*tb[i+1 ];
			cacc+=ta[i+2 ]*tb[i+2 ];
			cacc+=ta[i+3 ]*tb[i+3 ];
			cacc+=ta[i+4 ]*tb[i+4 ];
			cacc+=ta[i+5 ]*tb[i+5 ];
			cacc+=ta[i+6 ]*tb[i+6 ];
			cacc+=ta[i+7 ]*tb[i+7 ];
			cacc+=ta[i+8 ]*tb[i+8 ];
			cacc+=ta[i+9 ]*tb[i+9 ];
			cacc+=ta[i+10 ]*tb[i+10 ];
			cacc+=ta[i+11 ]*tb[i+11 ];
			cacc+=ta[i+12 ]*tb[i+12 ];
			cacc+=ta[i+13 ]*tb[i+13 ];
			cacc+=ta[i+14 ]*tb[i+14 ];
			cacc+=ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 		cacc+=ta[i+0 ]*tb[i+0 ];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	else
	{
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double*tb = b; const double*ta = a; double*tc = c,cacc = ((double)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
			cacc+=ta[inca*(i+1 )]*tb[incb*(i+1 )];
			cacc+=ta[inca*(i+2 )]*tb[incb*(i+2 )];
			cacc+=ta[inca*(i+3 )]*tb[incb*(i+3 )];
			cacc+=ta[inca*(i+4 )]*tb[incb*(i+4 )];
			cacc+=ta[inca*(i+5 )]*tb[incb*(i+5 )];
			cacc+=ta[inca*(i+6 )]*tb[incb*(i+6 )];
			cacc+=ta[inca*(i+7 )]*tb[incb*(i+7 )];
			cacc+=ta[inca*(i+8 )]*tb[incb*(i+8 )];
			cacc+=ta[inca*(i+9 )]*tb[incb*(i+9 )];
			cacc+=ta[inca*(i+10 )]*tb[incb*(i+10 )];
			cacc+=ta[inca*(i+11 )]*tb[incb*(i+11 )];
			cacc+=ta[inca*(i+12 )]*tb[incb*(i+12 )];
			cacc+=ta[inca*(i+13 )]*tb[incb*(i+13 )];
			cacc+=ta[inca*(i+14 )]*tb[incb*(i+14 )];
			cacc+=ta[inca*(i+15 )]*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float*tb = b; const float*ta = a; float*tc = c,cacc = ((float)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
			cacc+=ta[inca*(i+1 )]*tb[incb*(i+1 )];
			cacc+=ta[inca*(i+2 )]*tb[incb*(i+2 )];
			cacc+=ta[inca*(i+3 )]*tb[incb*(i+3 )];
			cacc+=ta[inca*(i+4 )]*tb[incb*(i+4 )];
			cacc+=ta[inca*(i+5 )]*tb[incb*(i+5 )];
			cacc+=ta[inca*(i+6 )]*tb[incb*(i+6 )];
			cacc+=ta[inca*(i+7 )]*tb[incb*(i+7 )];
			cacc+=ta[inca*(i+8 )]*tb[incb*(i+8 )];
			cacc+=ta[inca*(i+9 )]*tb[incb*(i+9 )];
			cacc+=ta[inca*(i+10 )]*tb[incb*(i+10 )];
			cacc+=ta[inca*(i+11 )]*tb[incb*(i+11 )];
			cacc+=ta[inca*(i+12 )]*tb[incb*(i+12 )];
			cacc+=ta[inca*(i+13 )]*tb[incb*(i+13 )];
			cacc+=ta[inca*(i+14 )]*tb[incb*(i+14 )];
			cacc+=ta[inca*(i+15 )]*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex*tb = b; const float complex*ta = a; float complex*tc = c,cacc = ((float complex)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
			cacc+=ta[inca*(i+1 )]*tb[incb*(i+1 )];
			cacc+=ta[inca*(i+2 )]*tb[incb*(i+2 )];
			cacc+=ta[inca*(i+3 )]*tb[incb*(i+3 )];
			cacc+=ta[inca*(i+4 )]*tb[incb*(i+4 )];
			cacc+=ta[inca*(i+5 )]*tb[incb*(i+5 )];
			cacc+=ta[inca*(i+6 )]*tb[incb*(i+6 )];
			cacc+=ta[inca*(i+7 )]*tb[incb*(i+7 )];
			cacc+=ta[inca*(i+8 )]*tb[incb*(i+8 )];
			cacc+=ta[inca*(i+9 )]*tb[incb*(i+9 )];
			cacc+=ta[inca*(i+10 )]*tb[incb*(i+10 )];
			cacc+=ta[inca*(i+11 )]*tb[incb*(i+11 )];
			cacc+=ta[inca*(i+12 )]*tb[incb*(i+12 )];
			cacc+=ta[inca*(i+13 )]*tb[incb*(i+13 )];
			cacc+=ta[inca*(i+14 )]*tb[incb*(i+14 )];
			cacc+=ta[inca*(i+15 )]*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex*tb = b; const double complex*ta = a; double complex*tc = c,cacc = ((double complex)(0));
		{
for(i=0;i+15<n;i+=16){
		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
			cacc+=ta[inca*(i+1 )]*tb[incb*(i+1 )];
			cacc+=ta[inca*(i+2 )]*tb[incb*(i+2 )];
			cacc+=ta[inca*(i+3 )]*tb[incb*(i+3 )];
			cacc+=ta[inca*(i+4 )]*tb[incb*(i+4 )];
			cacc+=ta[inca*(i+5 )]*tb[incb*(i+5 )];
			cacc+=ta[inca*(i+6 )]*tb[incb*(i+6 )];
			cacc+=ta[inca*(i+7 )]*tb[incb*(i+7 )];
			cacc+=ta[inca*(i+8 )]*tb[incb*(i+8 )];
			cacc+=ta[inca*(i+9 )]*tb[incb*(i+9 )];
			cacc+=ta[inca*(i+10 )]*tb[incb*(i+10 )];
			cacc+=ta[inca*(i+11 )]*tb[incb*(i+11 )];
			cacc+=ta[inca*(i+12 )]*tb[incb*(i+12 )];
			cacc+=ta[inca*(i+13 )]*tb[incb*(i+13 )];
			cacc+=ta[inca*(i+14 )]*tb[incb*(i+14 )];
			cacc+=ta[inca*(i+15 )]*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 		cacc+=ta[inca*(i+0 )]*tb[incb*(i+0 )];
	 }
}
; 
		*tc = cacc;
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	return RSB_ERR_NO_ERROR;
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

static rsb_err_t rsb_fill_with_zeros_nostride(void * array, rsb_type_t type, size_t n)
{
	/*!
	 * \ingroup gr_vec
	 * Will zero the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  ){
	double*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = ((double)(0));ta[i+1 ] = ((double)(0));ta[i+2 ] = ((double)(0));ta[i+3 ] = ((double)(0));ta[i+4 ] = ((double)(0));ta[i+5 ] = ((double)(0));ta[i+6 ] = ((double)(0));ta[i+7 ] = ((double)(0));ta[i+8 ] = ((double)(0));ta[i+9 ] = ((double)(0));ta[i+10 ] = ((double)(0));ta[i+11 ] = ((double)(0));ta[i+12 ] = ((double)(0));ta[i+13 ] = ((double)(0));ta[i+14 ] = ((double)(0));ta[i+15 ] = ((double)(0));}
for(     ;i<n;++i){ ta[i+0 ] = ((double)(0)); }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  ){
	float*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = ((float)(0));ta[i+1 ] = ((float)(0));ta[i+2 ] = ((float)(0));ta[i+3 ] = ((float)(0));ta[i+4 ] = ((float)(0));ta[i+5 ] = ((float)(0));ta[i+6 ] = ((float)(0));ta[i+7 ] = ((float)(0));ta[i+8 ] = ((float)(0));ta[i+9 ] = ((float)(0));ta[i+10 ] = ((float)(0));ta[i+11 ] = ((float)(0));ta[i+12 ] = ((float)(0));ta[i+13 ] = ((float)(0));ta[i+14 ] = ((float)(0));ta[i+15 ] = ((float)(0));}
for(     ;i<n;++i){ ta[i+0 ] = ((float)(0)); }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ){
	float complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = ((float complex)(0));ta[i+1 ] = ((float complex)(0));ta[i+2 ] = ((float complex)(0));ta[i+3 ] = ((float complex)(0));ta[i+4 ] = ((float complex)(0));ta[i+5 ] = ((float complex)(0));ta[i+6 ] = ((float complex)(0));ta[i+7 ] = ((float complex)(0));ta[i+8 ] = ((float complex)(0));ta[i+9 ] = ((float complex)(0));ta[i+10 ] = ((float complex)(0));ta[i+11 ] = ((float complex)(0));ta[i+12 ] = ((float complex)(0));ta[i+13 ] = ((float complex)(0));ta[i+14 ] = ((float complex)(0));ta[i+15 ] = ((float complex)(0));}
for(     ;i<n;++i){ ta[i+0 ] = ((float complex)(0)); }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ){
	double complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = ((double complex)(0));ta[i+1 ] = ((double complex)(0));ta[i+2 ] = ((double complex)(0));ta[i+3 ] = ((double complex)(0));ta[i+4 ] = ((double complex)(0));ta[i+5 ] = ((double complex)(0));ta[i+6 ] = ((double complex)(0));ta[i+7 ] = ((double complex)(0));ta[i+8 ] = ((double complex)(0));ta[i+9 ] = ((double complex)(0));ta[i+10 ] = ((double complex)(0));ta[i+11 ] = ((double complex)(0));ta[i+12 ] = ((double complex)(0));ta[i+13 ] = ((double complex)(0));ta[i+14 ] = ((double complex)(0));ta[i+15 ] = ((double complex)(0));}
for(     ;i<n;++i){ ta[i+0 ] = ((double complex)(0)); }
}
}
	else 
#endif
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

static rsb_err_t rsb_fill_with_zeros(void * array, rsb_type_t type, size_t n, size_t incx)
{
	/*!
	 * \ingroup gr_vec
	 * Will zero the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(incx==1)
		return rsb_fill_with_zeros_nostride(array,type,n);

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  ){
	double*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[(i+0 )*incx] = ((double)(0));ta[(i+1 )*incx] = ((double)(0));ta[(i+2 )*incx] = ((double)(0));ta[(i+3 )*incx] = ((double)(0));ta[(i+4 )*incx] = ((double)(0));ta[(i+5 )*incx] = ((double)(0));ta[(i+6 )*incx] = ((double)(0));ta[(i+7 )*incx] = ((double)(0));ta[(i+8 )*incx] = ((double)(0));ta[(i+9 )*incx] = ((double)(0));ta[(i+10 )*incx] = ((double)(0));ta[(i+11 )*incx] = ((double)(0));ta[(i+12 )*incx] = ((double)(0));ta[(i+13 )*incx] = ((double)(0));ta[(i+14 )*incx] = ((double)(0));ta[(i+15 )*incx] = ((double)(0));}
for(     ;i<n;++i){ ta[(i+0 )*incx] = ((double)(0)); }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  ){
	float*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[(i+0 )*incx] = ((float)(0));ta[(i+1 )*incx] = ((float)(0));ta[(i+2 )*incx] = ((float)(0));ta[(i+3 )*incx] = ((float)(0));ta[(i+4 )*incx] = ((float)(0));ta[(i+5 )*incx] = ((float)(0));ta[(i+6 )*incx] = ((float)(0));ta[(i+7 )*incx] = ((float)(0));ta[(i+8 )*incx] = ((float)(0));ta[(i+9 )*incx] = ((float)(0));ta[(i+10 )*incx] = ((float)(0));ta[(i+11 )*incx] = ((float)(0));ta[(i+12 )*incx] = ((float)(0));ta[(i+13 )*incx] = ((float)(0));ta[(i+14 )*incx] = ((float)(0));ta[(i+15 )*incx] = ((float)(0));}
for(     ;i<n;++i){ ta[(i+0 )*incx] = ((float)(0)); }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ){
	float complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[(i+0 )*incx] = ((float complex)(0));ta[(i+1 )*incx] = ((float complex)(0));ta[(i+2 )*incx] = ((float complex)(0));ta[(i+3 )*incx] = ((float complex)(0));ta[(i+4 )*incx] = ((float complex)(0));ta[(i+5 )*incx] = ((float complex)(0));ta[(i+6 )*incx] = ((float complex)(0));ta[(i+7 )*incx] = ((float complex)(0));ta[(i+8 )*incx] = ((float complex)(0));ta[(i+9 )*incx] = ((float complex)(0));ta[(i+10 )*incx] = ((float complex)(0));ta[(i+11 )*incx] = ((float complex)(0));ta[(i+12 )*incx] = ((float complex)(0));ta[(i+13 )*incx] = ((float complex)(0));ta[(i+14 )*incx] = ((float complex)(0));ta[(i+15 )*incx] = ((float complex)(0));}
for(     ;i<n;++i){ ta[(i+0 )*incx] = ((float complex)(0)); }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ){
	double complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[(i+0 )*incx] = ((double complex)(0));ta[(i+1 )*incx] = ((double complex)(0));ta[(i+2 )*incx] = ((double complex)(0));ta[(i+3 )*incx] = ((double complex)(0));ta[(i+4 )*incx] = ((double complex)(0));ta[(i+5 )*incx] = ((double complex)(0));ta[(i+6 )*incx] = ((double complex)(0));ta[(i+7 )*incx] = ((double complex)(0));ta[(i+8 )*incx] = ((double complex)(0));ta[(i+9 )*incx] = ((double complex)(0));ta[(i+10 )*incx] = ((double complex)(0));ta[(i+11 )*incx] = ((double complex)(0));ta[(i+12 )*incx] = ((double complex)(0));ta[(i+13 )*incx] = ((double complex)(0));ta[(i+14 )*incx] = ((double complex)(0));ta[(i+15 )*incx] = ((double complex)(0));}
for(     ;i<n;++i){ ta[(i+0 )*incx] = ((double complex)(0)); }
}
}
	else 
#endif
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

static rsb_err_t rsb_vector_scale(void * a, const void * alphap, rsb_type_t type, size_t n)
{
	/*!
	 * a <- a * alpha
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param alphap scaling value (if NULL assumed to be zero)
	 * \param n	the input array length
	 * \note see dscal in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(alphap==NULL || RSB_IS_ELEMENT_ZERO(alphap,type))
		return rsb_fill_with_zeros(a,type,n,1);
		
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double alpha = *(double*)alphap; double*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]*=alpha;
		ta[i+1 ]*=alpha;
		ta[i+2 ]*=alpha;
		ta[i+3 ]*=alpha;
		ta[i+4 ]*=alpha;
		ta[i+5 ]*=alpha;
		ta[i+6 ]*=alpha;
		ta[i+7 ]*=alpha;
		ta[i+8 ]*=alpha;
		ta[i+9 ]*=alpha;
		ta[i+10 ]*=alpha;
		ta[i+11 ]*=alpha;
		ta[i+12 ]*=alpha;
		ta[i+13 ]*=alpha;
		ta[i+14 ]*=alpha;
		ta[i+15 ]*=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]*=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float alpha = *(float*)alphap; float*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]*=alpha;
		ta[i+1 ]*=alpha;
		ta[i+2 ]*=alpha;
		ta[i+3 ]*=alpha;
		ta[i+4 ]*=alpha;
		ta[i+5 ]*=alpha;
		ta[i+6 ]*=alpha;
		ta[i+7 ]*=alpha;
		ta[i+8 ]*=alpha;
		ta[i+9 ]*=alpha;
		ta[i+10 ]*=alpha;
		ta[i+11 ]*=alpha;
		ta[i+12 ]*=alpha;
		ta[i+13 ]*=alpha;
		ta[i+14 ]*=alpha;
		ta[i+15 ]*=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]*=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex alpha = *(float complex*)alphap; float complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]*=alpha;
		ta[i+1 ]*=alpha;
		ta[i+2 ]*=alpha;
		ta[i+3 ]*=alpha;
		ta[i+4 ]*=alpha;
		ta[i+5 ]*=alpha;
		ta[i+6 ]*=alpha;
		ta[i+7 ]*=alpha;
		ta[i+8 ]*=alpha;
		ta[i+9 ]*=alpha;
		ta[i+10 ]*=alpha;
		ta[i+11 ]*=alpha;
		ta[i+12 ]*=alpha;
		ta[i+13 ]*=alpha;
		ta[i+14 ]*=alpha;
		ta[i+15 ]*=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]*=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex alpha = *(double complex*)alphap; double complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]*=alpha;
		ta[i+1 ]*=alpha;
		ta[i+2 ]*=alpha;
		ta[i+3 ]*=alpha;
		ta[i+4 ]*=alpha;
		ta[i+5 ]*=alpha;
		ta[i+6 ]*=alpha;
		ta[i+7 ]*=alpha;
		ta[i+8 ]*=alpha;
		ta[i+9 ]*=alpha;
		ta[i+10 ]*=alpha;
		ta[i+11 ]*=alpha;
		ta[i+12 ]*=alpha;
		ta[i+13 ]*=alpha;
		ta[i+14 ]*=alpha;
		ta[i+15 ]*=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]*=alpha;
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

static rsb_err_t rsb_strided_vector_scale(void * a, const void * alphap, rsb_type_t type, size_t n, size_t stride)
{
	/*!
	 * a <- a * alpha
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see dscal in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(stride==1)
		return rsb_vector_scale(a,alphap,type,n);
	if(alphap==NULL || RSB_IS_ELEMENT_ZERO(alphap,type))
		return rsb_fill_with_zeros(a,type,n,stride);

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double alpha = *(double*)alphap; double*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[stride*(i+0 )]*=alpha;
			ta[stride*(i+1 )]*=alpha;
			ta[stride*(i+2 )]*=alpha;
			ta[stride*(i+3 )]*=alpha;
			ta[stride*(i+4 )]*=alpha;
			ta[stride*(i+5 )]*=alpha;
			ta[stride*(i+6 )]*=alpha;
			ta[stride*(i+7 )]*=alpha;
			ta[stride*(i+8 )]*=alpha;
			ta[stride*(i+9 )]*=alpha;
			ta[stride*(i+10 )]*=alpha;
			ta[stride*(i+11 )]*=alpha;
			ta[stride*(i+12 )]*=alpha;
			ta[stride*(i+13 )]*=alpha;
			ta[stride*(i+14 )]*=alpha;
			ta[stride*(i+15 )]*=alpha;
	}
for(     ;i<n;++i){ 		ta[stride*(i+0 )]*=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float alpha = *(float*)alphap; float*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[stride*(i+0 )]*=alpha;
			ta[stride*(i+1 )]*=alpha;
			ta[stride*(i+2 )]*=alpha;
			ta[stride*(i+3 )]*=alpha;
			ta[stride*(i+4 )]*=alpha;
			ta[stride*(i+5 )]*=alpha;
			ta[stride*(i+6 )]*=alpha;
			ta[stride*(i+7 )]*=alpha;
			ta[stride*(i+8 )]*=alpha;
			ta[stride*(i+9 )]*=alpha;
			ta[stride*(i+10 )]*=alpha;
			ta[stride*(i+11 )]*=alpha;
			ta[stride*(i+12 )]*=alpha;
			ta[stride*(i+13 )]*=alpha;
			ta[stride*(i+14 )]*=alpha;
			ta[stride*(i+15 )]*=alpha;
	}
for(     ;i<n;++i){ 		ta[stride*(i+0 )]*=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex alpha = *(float complex*)alphap; float complex*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[stride*(i+0 )]*=alpha;
			ta[stride*(i+1 )]*=alpha;
			ta[stride*(i+2 )]*=alpha;
			ta[stride*(i+3 )]*=alpha;
			ta[stride*(i+4 )]*=alpha;
			ta[stride*(i+5 )]*=alpha;
			ta[stride*(i+6 )]*=alpha;
			ta[stride*(i+7 )]*=alpha;
			ta[stride*(i+8 )]*=alpha;
			ta[stride*(i+9 )]*=alpha;
			ta[stride*(i+10 )]*=alpha;
			ta[stride*(i+11 )]*=alpha;
			ta[stride*(i+12 )]*=alpha;
			ta[stride*(i+13 )]*=alpha;
			ta[stride*(i+14 )]*=alpha;
			ta[stride*(i+15 )]*=alpha;
	}
for(     ;i<n;++i){ 		ta[stride*(i+0 )]*=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex alpha = *(double complex*)alphap; double complex*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[stride*(i+0 )]*=alpha;
			ta[stride*(i+1 )]*=alpha;
			ta[stride*(i+2 )]*=alpha;
			ta[stride*(i+3 )]*=alpha;
			ta[stride*(i+4 )]*=alpha;
			ta[stride*(i+5 )]*=alpha;
			ta[stride*(i+6 )]*=alpha;
			ta[stride*(i+7 )]*=alpha;
			ta[stride*(i+8 )]*=alpha;
			ta[stride*(i+9 )]*=alpha;
			ta[stride*(i+10 )]*=alpha;
			ta[stride*(i+11 )]*=alpha;
			ta[stride*(i+12 )]*=alpha;
			ta[stride*(i+13 )]*=alpha;
			ta[stride*(i+14 )]*=alpha;
			ta[stride*(i+15 )]*=alpha;
	}
for(     ;i<n;++i){ 		ta[stride*(i+0 )]*=alpha;
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_vector_add(void * a, const void * alphap, rsb_type_t type, size_t n)
{
	/*!
	 * a <- a + alpha
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
		
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double alpha = *(double*)alphap; double*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[i+0 ]+=alpha;
			ta[i+1 ]+=alpha;
			ta[i+2 ]+=alpha;
			ta[i+3 ]+=alpha;
			ta[i+4 ]+=alpha;
			ta[i+5 ]+=alpha;
			ta[i+6 ]+=alpha;
			ta[i+7 ]+=alpha;
			ta[i+8 ]+=alpha;
			ta[i+9 ]+=alpha;
			ta[i+10 ]+=alpha;
			ta[i+11 ]+=alpha;
			ta[i+12 ]+=alpha;
			ta[i+13 ]+=alpha;
			ta[i+14 ]+=alpha;
			ta[i+15 ]+=alpha;
	}
for(     ;i<n;++i){ 		ta[i+0 ]+=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float alpha = *(float*)alphap; float*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[i+0 ]+=alpha;
			ta[i+1 ]+=alpha;
			ta[i+2 ]+=alpha;
			ta[i+3 ]+=alpha;
			ta[i+4 ]+=alpha;
			ta[i+5 ]+=alpha;
			ta[i+6 ]+=alpha;
			ta[i+7 ]+=alpha;
			ta[i+8 ]+=alpha;
			ta[i+9 ]+=alpha;
			ta[i+10 ]+=alpha;
			ta[i+11 ]+=alpha;
			ta[i+12 ]+=alpha;
			ta[i+13 ]+=alpha;
			ta[i+14 ]+=alpha;
			ta[i+15 ]+=alpha;
	}
for(     ;i<n;++i){ 		ta[i+0 ]+=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex alpha = *(float complex*)alphap; float complex*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[i+0 ]+=alpha;
			ta[i+1 ]+=alpha;
			ta[i+2 ]+=alpha;
			ta[i+3 ]+=alpha;
			ta[i+4 ]+=alpha;
			ta[i+5 ]+=alpha;
			ta[i+6 ]+=alpha;
			ta[i+7 ]+=alpha;
			ta[i+8 ]+=alpha;
			ta[i+9 ]+=alpha;
			ta[i+10 ]+=alpha;
			ta[i+11 ]+=alpha;
			ta[i+12 ]+=alpha;
			ta[i+13 ]+=alpha;
			ta[i+14 ]+=alpha;
			ta[i+15 ]+=alpha;
	}
for(     ;i<n;++i){ 		ta[i+0 ]+=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex alpha = *(double complex*)alphap; double complex*ta = a;
		{
for(i=0;i+15<n;i+=16){
		ta[i+0 ]+=alpha;
			ta[i+1 ]+=alpha;
			ta[i+2 ]+=alpha;
			ta[i+3 ]+=alpha;
			ta[i+4 ]+=alpha;
			ta[i+5 ]+=alpha;
			ta[i+6 ]+=alpha;
			ta[i+7 ]+=alpha;
			ta[i+8 ]+=alpha;
			ta[i+9 ]+=alpha;
			ta[i+10 ]+=alpha;
			ta[i+11 ]+=alpha;
			ta[i+12 ]+=alpha;
			ta[i+13 ]+=alpha;
			ta[i+14 ]+=alpha;
			ta[i+15 ]+=alpha;
	}
for(     ;i<n;++i){ 		ta[i+0 ]+=alpha;
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_vector_div(void * a, const void * alphap, rsb_type_t type, size_t n)
{
	/*!
	 * this is a benchmark-oriented function only..
	 * \return \rsberrcodemsg
	 * */
	size_t i;
		
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double alpha = *(double*)alphap; double*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]/=alpha;
		ta[i+1 ]/=alpha;
		ta[i+2 ]/=alpha;
		ta[i+3 ]/=alpha;
		ta[i+4 ]/=alpha;
		ta[i+5 ]/=alpha;
		ta[i+6 ]/=alpha;
		ta[i+7 ]/=alpha;
		ta[i+8 ]/=alpha;
		ta[i+9 ]/=alpha;
		ta[i+10 ]/=alpha;
		ta[i+11 ]/=alpha;
		ta[i+12 ]/=alpha;
		ta[i+13 ]/=alpha;
		ta[i+14 ]/=alpha;
		ta[i+15 ]/=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]/=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float alpha = *(float*)alphap; float*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]/=alpha;
		ta[i+1 ]/=alpha;
		ta[i+2 ]/=alpha;
		ta[i+3 ]/=alpha;
		ta[i+4 ]/=alpha;
		ta[i+5 ]/=alpha;
		ta[i+6 ]/=alpha;
		ta[i+7 ]/=alpha;
		ta[i+8 ]/=alpha;
		ta[i+9 ]/=alpha;
		ta[i+10 ]/=alpha;
		ta[i+11 ]/=alpha;
		ta[i+12 ]/=alpha;
		ta[i+13 ]/=alpha;
		ta[i+14 ]/=alpha;
		ta[i+15 ]/=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]/=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex alpha = *(float complex*)alphap; float complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]/=alpha;
		ta[i+1 ]/=alpha;
		ta[i+2 ]/=alpha;
		ta[i+3 ]/=alpha;
		ta[i+4 ]/=alpha;
		ta[i+5 ]/=alpha;
		ta[i+6 ]/=alpha;
		ta[i+7 ]/=alpha;
		ta[i+8 ]/=alpha;
		ta[i+9 ]/=alpha;
		ta[i+10 ]/=alpha;
		ta[i+11 ]/=alpha;
		ta[i+12 ]/=alpha;
		ta[i+13 ]/=alpha;
		ta[i+14 ]/=alpha;
		ta[i+15 ]/=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]/=alpha;
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex alpha = *(double complex*)alphap; double complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]/=alpha;
		ta[i+1 ]/=alpha;
		ta[i+2 ]/=alpha;
		ta[i+3 ]/=alpha;
		ta[i+4 ]/=alpha;
		ta[i+5 ]/=alpha;
		ta[i+6 ]/=alpha;
		ta[i+7 ]/=alpha;
		ta[i+8 ]/=alpha;
		ta[i+9 ]/=alpha;
		ta[i+10 ]/=alpha;
		ta[i+11 ]/=alpha;
		ta[i+12 ]/=alpha;
		ta[i+13 ]/=alpha;
		ta[i+14 ]/=alpha;
		ta[i+15 ]/=alpha;
	}
for(     ;i<n;++i){ 	ta[i+0 ]/=alpha;
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__vector_increase_by_one(void * a, rsb_type_t type, size_t n)
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
		
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{ double*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=((double)(1.0));
		ta[i+1 ]+=((double)(1.0));
		ta[i+2 ]+=((double)(1.0));
		ta[i+3 ]+=((double)(1.0));
		ta[i+4 ]+=((double)(1.0));
		ta[i+5 ]+=((double)(1.0));
		ta[i+6 ]+=((double)(1.0));
		ta[i+7 ]+=((double)(1.0));
		ta[i+8 ]+=((double)(1.0));
		ta[i+9 ]+=((double)(1.0));
		ta[i+10 ]+=((double)(1.0));
		ta[i+11 ]+=((double)(1.0));
		ta[i+12 ]+=((double)(1.0));
		ta[i+13 ]+=((double)(1.0));
		ta[i+14 ]+=((double)(1.0));
		ta[i+15 ]+=((double)(1.0));
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=((double)(1.0));
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{ float*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=((float)(1.0));
		ta[i+1 ]+=((float)(1.0));
		ta[i+2 ]+=((float)(1.0));
		ta[i+3 ]+=((float)(1.0));
		ta[i+4 ]+=((float)(1.0));
		ta[i+5 ]+=((float)(1.0));
		ta[i+6 ]+=((float)(1.0));
		ta[i+7 ]+=((float)(1.0));
		ta[i+8 ]+=((float)(1.0));
		ta[i+9 ]+=((float)(1.0));
		ta[i+10 ]+=((float)(1.0));
		ta[i+11 ]+=((float)(1.0));
		ta[i+12 ]+=((float)(1.0));
		ta[i+13 ]+=((float)(1.0));
		ta[i+14 ]+=((float)(1.0));
		ta[i+15 ]+=((float)(1.0));
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=((float)(1.0));
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{ float complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=((float complex)(1.0));
		ta[i+1 ]+=((float complex)(1.0));
		ta[i+2 ]+=((float complex)(1.0));
		ta[i+3 ]+=((float complex)(1.0));
		ta[i+4 ]+=((float complex)(1.0));
		ta[i+5 ]+=((float complex)(1.0));
		ta[i+6 ]+=((float complex)(1.0));
		ta[i+7 ]+=((float complex)(1.0));
		ta[i+8 ]+=((float complex)(1.0));
		ta[i+9 ]+=((float complex)(1.0));
		ta[i+10 ]+=((float complex)(1.0));
		ta[i+11 ]+=((float complex)(1.0));
		ta[i+12 ]+=((float complex)(1.0));
		ta[i+13 ]+=((float complex)(1.0));
		ta[i+14 ]+=((float complex)(1.0));
		ta[i+15 ]+=((float complex)(1.0));
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=((float complex)(1.0));
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{ double complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=((double complex)(1.0));
		ta[i+1 ]+=((double complex)(1.0));
		ta[i+2 ]+=((double complex)(1.0));
		ta[i+3 ]+=((double complex)(1.0));
		ta[i+4 ]+=((double complex)(1.0));
		ta[i+5 ]+=((double complex)(1.0));
		ta[i+6 ]+=((double complex)(1.0));
		ta[i+7 ]+=((double complex)(1.0));
		ta[i+8 ]+=((double complex)(1.0));
		ta[i+9 ]+=((double complex)(1.0));
		ta[i+10 ]+=((double complex)(1.0));
		ta[i+11 ]+=((double complex)(1.0));
		ta[i+12 ]+=((double complex)(1.0));
		ta[i+13 ]+=((double complex)(1.0));
		ta[i+14 ]+=((double complex)(1.0));
		ta[i+15 ]+=((double complex)(1.0));
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=((double complex)(1.0));
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_vector_pow(void * a, rsb_type_t type, const void *y, size_t n)
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(!a || !y)
		goto err;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		double ty = *(double*)y,*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = pow(ta[i+0 ],ty);
		ta[i+1 ] = pow(ta[i+1 ],ty);
		ta[i+2 ] = pow(ta[i+2 ],ty);
		ta[i+3 ] = pow(ta[i+3 ],ty);
		ta[i+4 ] = pow(ta[i+4 ],ty);
		ta[i+5 ] = pow(ta[i+5 ],ty);
		ta[i+6 ] = pow(ta[i+6 ],ty);
		ta[i+7 ] = pow(ta[i+7 ],ty);
		ta[i+8 ] = pow(ta[i+8 ],ty);
		ta[i+9 ] = pow(ta[i+9 ],ty);
		ta[i+10 ] = pow(ta[i+10 ],ty);
		ta[i+11 ] = pow(ta[i+11 ],ty);
		ta[i+12 ] = pow(ta[i+12 ],ty);
		ta[i+13 ] = pow(ta[i+13 ],ty);
		ta[i+14 ] = pow(ta[i+14 ],ty);
		ta[i+15 ] = pow(ta[i+15 ],ty);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = pow(ta[i+0 ],ty);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		float ty = *(float*)y,*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = powf(ta[i+0 ],ty);
		ta[i+1 ] = powf(ta[i+1 ],ty);
		ta[i+2 ] = powf(ta[i+2 ],ty);
		ta[i+3 ] = powf(ta[i+3 ],ty);
		ta[i+4 ] = powf(ta[i+4 ],ty);
		ta[i+5 ] = powf(ta[i+5 ],ty);
		ta[i+6 ] = powf(ta[i+6 ],ty);
		ta[i+7 ] = powf(ta[i+7 ],ty);
		ta[i+8 ] = powf(ta[i+8 ],ty);
		ta[i+9 ] = powf(ta[i+9 ],ty);
		ta[i+10 ] = powf(ta[i+10 ],ty);
		ta[i+11 ] = powf(ta[i+11 ],ty);
		ta[i+12 ] = powf(ta[i+12 ],ty);
		ta[i+13 ] = powf(ta[i+13 ],ty);
		ta[i+14 ] = powf(ta[i+14 ],ty);
		ta[i+15 ] = powf(ta[i+15 ],ty);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = powf(ta[i+0 ],ty);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		float complex ty = *(float complex*)y,*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = cpowf(ta[i+0 ],ty);
		ta[i+1 ] = cpowf(ta[i+1 ],ty);
		ta[i+2 ] = cpowf(ta[i+2 ],ty);
		ta[i+3 ] = cpowf(ta[i+3 ],ty);
		ta[i+4 ] = cpowf(ta[i+4 ],ty);
		ta[i+5 ] = cpowf(ta[i+5 ],ty);
		ta[i+6 ] = cpowf(ta[i+6 ],ty);
		ta[i+7 ] = cpowf(ta[i+7 ],ty);
		ta[i+8 ] = cpowf(ta[i+8 ],ty);
		ta[i+9 ] = cpowf(ta[i+9 ],ty);
		ta[i+10 ] = cpowf(ta[i+10 ],ty);
		ta[i+11 ] = cpowf(ta[i+11 ],ty);
		ta[i+12 ] = cpowf(ta[i+12 ],ty);
		ta[i+13 ] = cpowf(ta[i+13 ],ty);
		ta[i+14 ] = cpowf(ta[i+14 ],ty);
		ta[i+15 ] = cpowf(ta[i+15 ],ty);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = cpowf(ta[i+0 ],ty);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		double complex ty = *(double complex*)y,*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = cpow(ta[i+0 ],ty);
		ta[i+1 ] = cpow(ta[i+1 ],ty);
		ta[i+2 ] = cpow(ta[i+2 ],ty);
		ta[i+3 ] = cpow(ta[i+3 ],ty);
		ta[i+4 ] = cpow(ta[i+4 ],ty);
		ta[i+5 ] = cpow(ta[i+5 ],ty);
		ta[i+6 ] = cpow(ta[i+6 ],ty);
		ta[i+7 ] = cpow(ta[i+7 ],ty);
		ta[i+8 ] = cpow(ta[i+8 ],ty);
		ta[i+9 ] = cpow(ta[i+9 ],ty);
		ta[i+10 ] = cpow(ta[i+10 ],ty);
		ta[i+11 ] = cpow(ta[i+11 ],ty);
		ta[i+12 ] = cpow(ta[i+12 ],ty);
		ta[i+13 ] = cpow(ta[i+13 ],ty);
		ta[i+14 ] = cpow(ta[i+14 ],ty);
		ta[i+15 ] = cpow(ta[i+15 ],ty);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = cpow(ta[i+0 ],ty);
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
err:
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_vector_sqrt(void * a, rsb_type_t type, size_t n)
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(!a)goto err;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{double*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = sqrt(ta[i+0 ]);
		ta[i+1 ] = sqrt(ta[i+1 ]);
		ta[i+2 ] = sqrt(ta[i+2 ]);
		ta[i+3 ] = sqrt(ta[i+3 ]);
		ta[i+4 ] = sqrt(ta[i+4 ]);
		ta[i+5 ] = sqrt(ta[i+5 ]);
		ta[i+6 ] = sqrt(ta[i+6 ]);
		ta[i+7 ] = sqrt(ta[i+7 ]);
		ta[i+8 ] = sqrt(ta[i+8 ]);
		ta[i+9 ] = sqrt(ta[i+9 ]);
		ta[i+10 ] = sqrt(ta[i+10 ]);
		ta[i+11 ] = sqrt(ta[i+11 ]);
		ta[i+12 ] = sqrt(ta[i+12 ]);
		ta[i+13 ] = sqrt(ta[i+13 ]);
		ta[i+14 ] = sqrt(ta[i+14 ]);
		ta[i+15 ] = sqrt(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = sqrt(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{float*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = sqrtf(ta[i+0 ]);
		ta[i+1 ] = sqrtf(ta[i+1 ]);
		ta[i+2 ] = sqrtf(ta[i+2 ]);
		ta[i+3 ] = sqrtf(ta[i+3 ]);
		ta[i+4 ] = sqrtf(ta[i+4 ]);
		ta[i+5 ] = sqrtf(ta[i+5 ]);
		ta[i+6 ] = sqrtf(ta[i+6 ]);
		ta[i+7 ] = sqrtf(ta[i+7 ]);
		ta[i+8 ] = sqrtf(ta[i+8 ]);
		ta[i+9 ] = sqrtf(ta[i+9 ]);
		ta[i+10 ] = sqrtf(ta[i+10 ]);
		ta[i+11 ] = sqrtf(ta[i+11 ]);
		ta[i+12 ] = sqrtf(ta[i+12 ]);
		ta[i+13 ] = sqrtf(ta[i+13 ]);
		ta[i+14 ] = sqrtf(ta[i+14 ]);
		ta[i+15 ] = sqrtf(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = sqrtf(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{float complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = csqrtf(ta[i+0 ]);
		ta[i+1 ] = csqrtf(ta[i+1 ]);
		ta[i+2 ] = csqrtf(ta[i+2 ]);
		ta[i+3 ] = csqrtf(ta[i+3 ]);
		ta[i+4 ] = csqrtf(ta[i+4 ]);
		ta[i+5 ] = csqrtf(ta[i+5 ]);
		ta[i+6 ] = csqrtf(ta[i+6 ]);
		ta[i+7 ] = csqrtf(ta[i+7 ]);
		ta[i+8 ] = csqrtf(ta[i+8 ]);
		ta[i+9 ] = csqrtf(ta[i+9 ]);
		ta[i+10 ] = csqrtf(ta[i+10 ]);
		ta[i+11 ] = csqrtf(ta[i+11 ]);
		ta[i+12 ] = csqrtf(ta[i+12 ]);
		ta[i+13 ] = csqrtf(ta[i+13 ]);
		ta[i+14 ] = csqrtf(ta[i+14 ]);
		ta[i+15 ] = csqrtf(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = csqrtf(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{double complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = csqrt(ta[i+0 ]);
		ta[i+1 ] = csqrt(ta[i+1 ]);
		ta[i+2 ] = csqrt(ta[i+2 ]);
		ta[i+3 ] = csqrt(ta[i+3 ]);
		ta[i+4 ] = csqrt(ta[i+4 ]);
		ta[i+5 ] = csqrt(ta[i+5 ]);
		ta[i+6 ] = csqrt(ta[i+6 ]);
		ta[i+7 ] = csqrt(ta[i+7 ]);
		ta[i+8 ] = csqrt(ta[i+8 ]);
		ta[i+9 ] = csqrt(ta[i+9 ]);
		ta[i+10 ] = csqrt(ta[i+10 ]);
		ta[i+11 ] = csqrt(ta[i+11 ]);
		ta[i+12 ] = csqrt(ta[i+12 ]);
		ta[i+13 ] = csqrt(ta[i+13 ]);
		ta[i+14 ] = csqrt(ta[i+14 ]);
		ta[i+15 ] = csqrt(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = csqrt(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
err:
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__vector_scale_inv(void * a, const void * alphap, rsb_type_t type, size_t n)
{
	/*!
	 * a <- 1/a * alpha
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see dscal in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	if(!alphap)
		return RSB_ERR_BADARGS;
		
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		double alphai = ((double)(1.0))/(*(double*)alphap);
		return rsb_vector_scale(a,&alphai,type,n);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		float alphai = ((float)(1.0))/(*(float*)alphap);
		return rsb_vector_scale(a,&alphai,type,n);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		float complex alphai = ((float complex)(1.0))/(*(float complex*)alphap);
		return rsb_vector_scale(a,&alphai,type,n);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		double complex alphai = ((double complex)(1.0))/(*(double complex*)alphap);
		return rsb_vector_scale(a,&alphai,type,n);
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__vector_sum_of_abs_diffs(void * c, const void * a, const void * b, rsb_type_t type, size_t n)
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double*ap = a,*bp = b;
		double ac = ((double)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += fabs(ap[i+0 ]-bp[i+0 ]);
				ac += fabs(ap[i+1 ]-bp[i+1 ]);
				ac += fabs(ap[i+2 ]-bp[i+2 ]);
				ac += fabs(ap[i+3 ]-bp[i+3 ]);
				ac += fabs(ap[i+4 ]-bp[i+4 ]);
				ac += fabs(ap[i+5 ]-bp[i+5 ]);
				ac += fabs(ap[i+6 ]-bp[i+6 ]);
				ac += fabs(ap[i+7 ]-bp[i+7 ]);
				ac += fabs(ap[i+8 ]-bp[i+8 ]);
				ac += fabs(ap[i+9 ]-bp[i+9 ]);
				ac += fabs(ap[i+10 ]-bp[i+10 ]);
				ac += fabs(ap[i+11 ]-bp[i+11 ]);
				ac += fabs(ap[i+12 ]-bp[i+12 ]);
				ac += fabs(ap[i+13 ]-bp[i+13 ]);
				ac += fabs(ap[i+14 ]-bp[i+14 ]);
				ac += fabs(ap[i+15 ]-bp[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += fabs(ap[i+0 ]-bp[i+0 ]);
		 }
}
; 
		*((double*)(c)) = ac;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float*ap = a,*bp = b;
		float ac = ((float)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += fabsf(ap[i+0 ]-bp[i+0 ]);
				ac += fabsf(ap[i+1 ]-bp[i+1 ]);
				ac += fabsf(ap[i+2 ]-bp[i+2 ]);
				ac += fabsf(ap[i+3 ]-bp[i+3 ]);
				ac += fabsf(ap[i+4 ]-bp[i+4 ]);
				ac += fabsf(ap[i+5 ]-bp[i+5 ]);
				ac += fabsf(ap[i+6 ]-bp[i+6 ]);
				ac += fabsf(ap[i+7 ]-bp[i+7 ]);
				ac += fabsf(ap[i+8 ]-bp[i+8 ]);
				ac += fabsf(ap[i+9 ]-bp[i+9 ]);
				ac += fabsf(ap[i+10 ]-bp[i+10 ]);
				ac += fabsf(ap[i+11 ]-bp[i+11 ]);
				ac += fabsf(ap[i+12 ]-bp[i+12 ]);
				ac += fabsf(ap[i+13 ]-bp[i+13 ]);
				ac += fabsf(ap[i+14 ]-bp[i+14 ]);
				ac += fabsf(ap[i+15 ]-bp[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += fabsf(ap[i+0 ]-bp[i+0 ]);
		 }
}
; 
		*((float*)(c)) = ac;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex*ap = a,*bp = b;
		float complex ac = ((float complex)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += cabsf(ap[i+0 ]-bp[i+0 ]);
				ac += cabsf(ap[i+1 ]-bp[i+1 ]);
				ac += cabsf(ap[i+2 ]-bp[i+2 ]);
				ac += cabsf(ap[i+3 ]-bp[i+3 ]);
				ac += cabsf(ap[i+4 ]-bp[i+4 ]);
				ac += cabsf(ap[i+5 ]-bp[i+5 ]);
				ac += cabsf(ap[i+6 ]-bp[i+6 ]);
				ac += cabsf(ap[i+7 ]-bp[i+7 ]);
				ac += cabsf(ap[i+8 ]-bp[i+8 ]);
				ac += cabsf(ap[i+9 ]-bp[i+9 ]);
				ac += cabsf(ap[i+10 ]-bp[i+10 ]);
				ac += cabsf(ap[i+11 ]-bp[i+11 ]);
				ac += cabsf(ap[i+12 ]-bp[i+12 ]);
				ac += cabsf(ap[i+13 ]-bp[i+13 ]);
				ac += cabsf(ap[i+14 ]-bp[i+14 ]);
				ac += cabsf(ap[i+15 ]-bp[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += cabsf(ap[i+0 ]-bp[i+0 ]);
		 }
}
; 
		*((float complex*)(c)) = ac;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex*ap = a,*bp = b;
		double complex ac = ((double complex)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += cabs(ap[i+0 ]-bp[i+0 ]);
				ac += cabs(ap[i+1 ]-bp[i+1 ]);
				ac += cabs(ap[i+2 ]-bp[i+2 ]);
				ac += cabs(ap[i+3 ]-bp[i+3 ]);
				ac += cabs(ap[i+4 ]-bp[i+4 ]);
				ac += cabs(ap[i+5 ]-bp[i+5 ]);
				ac += cabs(ap[i+6 ]-bp[i+6 ]);
				ac += cabs(ap[i+7 ]-bp[i+7 ]);
				ac += cabs(ap[i+8 ]-bp[i+8 ]);
				ac += cabs(ap[i+9 ]-bp[i+9 ]);
				ac += cabs(ap[i+10 ]-bp[i+10 ]);
				ac += cabs(ap[i+11 ]-bp[i+11 ]);
				ac += cabs(ap[i+12 ]-bp[i+12 ]);
				ac += cabs(ap[i+13 ]-bp[i+13 ]);
				ac += cabs(ap[i+14 ]-bp[i+14 ]);
				ac += cabs(ap[i+15 ]-bp[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += cabs(ap[i+0 ]-bp[i+0 ]);
		 }
}
; 
		*((double complex*)(c)) = ac;
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__vector_sum_of_abs(void * c, const void * a, rsb_type_t type, size_t n)
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double*ap = a;
		double ac = ((double)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += fabs(ap[i+0 ]);
				ac += fabs(ap[i+1 ]);
				ac += fabs(ap[i+2 ]);
				ac += fabs(ap[i+3 ]);
				ac += fabs(ap[i+4 ]);
				ac += fabs(ap[i+5 ]);
				ac += fabs(ap[i+6 ]);
				ac += fabs(ap[i+7 ]);
				ac += fabs(ap[i+8 ]);
				ac += fabs(ap[i+9 ]);
				ac += fabs(ap[i+10 ]);
				ac += fabs(ap[i+11 ]);
				ac += fabs(ap[i+12 ]);
				ac += fabs(ap[i+13 ]);
				ac += fabs(ap[i+14 ]);
				ac += fabs(ap[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += fabs(ap[i+0 ]);
		 }
}
; 
		*((double*)(c)) = ac;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float*ap = a;
		float ac = ((float)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += fabsf(ap[i+0 ]);
				ac += fabsf(ap[i+1 ]);
				ac += fabsf(ap[i+2 ]);
				ac += fabsf(ap[i+3 ]);
				ac += fabsf(ap[i+4 ]);
				ac += fabsf(ap[i+5 ]);
				ac += fabsf(ap[i+6 ]);
				ac += fabsf(ap[i+7 ]);
				ac += fabsf(ap[i+8 ]);
				ac += fabsf(ap[i+9 ]);
				ac += fabsf(ap[i+10 ]);
				ac += fabsf(ap[i+11 ]);
				ac += fabsf(ap[i+12 ]);
				ac += fabsf(ap[i+13 ]);
				ac += fabsf(ap[i+14 ]);
				ac += fabsf(ap[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += fabsf(ap[i+0 ]);
		 }
}
; 
		*((float*)(c)) = ac;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex*ap = a;
		float complex ac = ((float complex)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += cabsf(ap[i+0 ]);
				ac += cabsf(ap[i+1 ]);
				ac += cabsf(ap[i+2 ]);
				ac += cabsf(ap[i+3 ]);
				ac += cabsf(ap[i+4 ]);
				ac += cabsf(ap[i+5 ]);
				ac += cabsf(ap[i+6 ]);
				ac += cabsf(ap[i+7 ]);
				ac += cabsf(ap[i+8 ]);
				ac += cabsf(ap[i+9 ]);
				ac += cabsf(ap[i+10 ]);
				ac += cabsf(ap[i+11 ]);
				ac += cabsf(ap[i+12 ]);
				ac += cabsf(ap[i+13 ]);
				ac += cabsf(ap[i+14 ]);
				ac += cabsf(ap[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += cabsf(ap[i+0 ]);
		 }
}
; 
		*((float complex*)(c)) = ac;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex*ap = a;
		double complex ac = ((double complex)(0));
		{
for(i=0;i+15<n;i+=16){
		ac += cabs(ap[i+0 ]);
				ac += cabs(ap[i+1 ]);
				ac += cabs(ap[i+2 ]);
				ac += cabs(ap[i+3 ]);
				ac += cabs(ap[i+4 ]);
				ac += cabs(ap[i+5 ]);
				ac += cabs(ap[i+6 ]);
				ac += cabs(ap[i+7 ]);
				ac += cabs(ap[i+8 ]);
				ac += cabs(ap[i+9 ]);
				ac += cabs(ap[i+10 ]);
				ac += cabs(ap[i+11 ]);
				ac += cabs(ap[i+12 ]);
				ac += cabs(ap[i+13 ]);
				ac += cabs(ap[i+14 ]);
				ac += cabs(ap[i+15 ]);
		}
for(     ;i<n;++i){ 		ac += cabs(ap[i+0 ]);
		 }
}
; 
		*((double complex*)(c)) = ac;
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__vector_to_abs(void * a, rsb_type_t type, size_t n)
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{double*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = fabs(ta[i+0 ]);
		ta[i+1 ] = fabs(ta[i+1 ]);
		ta[i+2 ] = fabs(ta[i+2 ]);
		ta[i+3 ] = fabs(ta[i+3 ]);
		ta[i+4 ] = fabs(ta[i+4 ]);
		ta[i+5 ] = fabs(ta[i+5 ]);
		ta[i+6 ] = fabs(ta[i+6 ]);
		ta[i+7 ] = fabs(ta[i+7 ]);
		ta[i+8 ] = fabs(ta[i+8 ]);
		ta[i+9 ] = fabs(ta[i+9 ]);
		ta[i+10 ] = fabs(ta[i+10 ]);
		ta[i+11 ] = fabs(ta[i+11 ]);
		ta[i+12 ] = fabs(ta[i+12 ]);
		ta[i+13 ] = fabs(ta[i+13 ]);
		ta[i+14 ] = fabs(ta[i+14 ]);
		ta[i+15 ] = fabs(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = fabs(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{float*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = fabsf(ta[i+0 ]);
		ta[i+1 ] = fabsf(ta[i+1 ]);
		ta[i+2 ] = fabsf(ta[i+2 ]);
		ta[i+3 ] = fabsf(ta[i+3 ]);
		ta[i+4 ] = fabsf(ta[i+4 ]);
		ta[i+5 ] = fabsf(ta[i+5 ]);
		ta[i+6 ] = fabsf(ta[i+6 ]);
		ta[i+7 ] = fabsf(ta[i+7 ]);
		ta[i+8 ] = fabsf(ta[i+8 ]);
		ta[i+9 ] = fabsf(ta[i+9 ]);
		ta[i+10 ] = fabsf(ta[i+10 ]);
		ta[i+11 ] = fabsf(ta[i+11 ]);
		ta[i+12 ] = fabsf(ta[i+12 ]);
		ta[i+13 ] = fabsf(ta[i+13 ]);
		ta[i+14 ] = fabsf(ta[i+14 ]);
		ta[i+15 ] = fabsf(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = fabsf(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{float complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = cabsf(ta[i+0 ]);
		ta[i+1 ] = cabsf(ta[i+1 ]);
		ta[i+2 ] = cabsf(ta[i+2 ]);
		ta[i+3 ] = cabsf(ta[i+3 ]);
		ta[i+4 ] = cabsf(ta[i+4 ]);
		ta[i+5 ] = cabsf(ta[i+5 ]);
		ta[i+6 ] = cabsf(ta[i+6 ]);
		ta[i+7 ] = cabsf(ta[i+7 ]);
		ta[i+8 ] = cabsf(ta[i+8 ]);
		ta[i+9 ] = cabsf(ta[i+9 ]);
		ta[i+10 ] = cabsf(ta[i+10 ]);
		ta[i+11 ] = cabsf(ta[i+11 ]);
		ta[i+12 ] = cabsf(ta[i+12 ]);
		ta[i+13 ] = cabsf(ta[i+13 ]);
		ta[i+14 ] = cabsf(ta[i+14 ]);
		ta[i+15 ] = cabsf(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = cabsf(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{double complex*ta = a;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ] = cabs(ta[i+0 ]);
		ta[i+1 ] = cabs(ta[i+1 ]);
		ta[i+2 ] = cabs(ta[i+2 ]);
		ta[i+3 ] = cabs(ta[i+3 ]);
		ta[i+4 ] = cabs(ta[i+4 ]);
		ta[i+5 ] = cabs(ta[i+5 ]);
		ta[i+6 ] = cabs(ta[i+6 ]);
		ta[i+7 ] = cabs(ta[i+7 ]);
		ta[i+8 ] = cabs(ta[i+8 ]);
		ta[i+9 ] = cabs(ta[i+9 ]);
		ta[i+10 ] = cabs(ta[i+10 ]);
		ta[i+11 ] = cabs(ta[i+11 ]);
		ta[i+12 ] = cabs(ta[i+12 ]);
		ta[i+13 ] = cabs(ta[i+13 ]);
		ta[i+14 ] = cabs(ta[i+14 ]);
		ta[i+15 ] = cabs(ta[i+15 ]);
	}
for(     ;i<n;++i){ 	ta[i+0 ] = cabs(ta[i+0 ]);
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

static rsb_err_t rsb_alpha_sum(void * a, const void * b, const void * alphap, rsb_type_t type, size_t n)
{
	/*!
	 * a <- a + alpha * b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see daxpy in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double alpha = alphap ? *(double*)alphap : ((double)(1.0));
	double*ta = a; const double*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=alpha*tb[i+0 ];
		ta[i+1 ]+=alpha*tb[i+1 ];
		ta[i+2 ]+=alpha*tb[i+2 ];
		ta[i+3 ]+=alpha*tb[i+3 ];
		ta[i+4 ]+=alpha*tb[i+4 ];
		ta[i+5 ]+=alpha*tb[i+5 ];
		ta[i+6 ]+=alpha*tb[i+6 ];
		ta[i+7 ]+=alpha*tb[i+7 ];
		ta[i+8 ]+=alpha*tb[i+8 ];
		ta[i+9 ]+=alpha*tb[i+9 ];
		ta[i+10 ]+=alpha*tb[i+10 ];
		ta[i+11 ]+=alpha*tb[i+11 ];
		ta[i+12 ]+=alpha*tb[i+12 ];
		ta[i+13 ]+=alpha*tb[i+13 ];
		ta[i+14 ]+=alpha*tb[i+14 ];
		ta[i+15 ]+=alpha*tb[i+15 ];
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=alpha*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float alpha = alphap ? *(float*)alphap : ((float)(1.0));
	float*ta = a; const float*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=alpha*tb[i+0 ];
		ta[i+1 ]+=alpha*tb[i+1 ];
		ta[i+2 ]+=alpha*tb[i+2 ];
		ta[i+3 ]+=alpha*tb[i+3 ];
		ta[i+4 ]+=alpha*tb[i+4 ];
		ta[i+5 ]+=alpha*tb[i+5 ];
		ta[i+6 ]+=alpha*tb[i+6 ];
		ta[i+7 ]+=alpha*tb[i+7 ];
		ta[i+8 ]+=alpha*tb[i+8 ];
		ta[i+9 ]+=alpha*tb[i+9 ];
		ta[i+10 ]+=alpha*tb[i+10 ];
		ta[i+11 ]+=alpha*tb[i+11 ];
		ta[i+12 ]+=alpha*tb[i+12 ];
		ta[i+13 ]+=alpha*tb[i+13 ];
		ta[i+14 ]+=alpha*tb[i+14 ];
		ta[i+15 ]+=alpha*tb[i+15 ];
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=alpha*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex alpha = alphap ? *(float complex*)alphap : ((float complex)(1.0));
	float complex*ta = a; const float complex*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=alpha*tb[i+0 ];
		ta[i+1 ]+=alpha*tb[i+1 ];
		ta[i+2 ]+=alpha*tb[i+2 ];
		ta[i+3 ]+=alpha*tb[i+3 ];
		ta[i+4 ]+=alpha*tb[i+4 ];
		ta[i+5 ]+=alpha*tb[i+5 ];
		ta[i+6 ]+=alpha*tb[i+6 ];
		ta[i+7 ]+=alpha*tb[i+7 ];
		ta[i+8 ]+=alpha*tb[i+8 ];
		ta[i+9 ]+=alpha*tb[i+9 ];
		ta[i+10 ]+=alpha*tb[i+10 ];
		ta[i+11 ]+=alpha*tb[i+11 ];
		ta[i+12 ]+=alpha*tb[i+12 ];
		ta[i+13 ]+=alpha*tb[i+13 ];
		ta[i+14 ]+=alpha*tb[i+14 ];
		ta[i+15 ]+=alpha*tb[i+15 ];
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=alpha*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex alpha = alphap ? *(double complex*)alphap : ((double complex)(1.0));
	double complex*ta = a; const double complex*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[i+0 ]+=alpha*tb[i+0 ];
		ta[i+1 ]+=alpha*tb[i+1 ];
		ta[i+2 ]+=alpha*tb[i+2 ];
		ta[i+3 ]+=alpha*tb[i+3 ];
		ta[i+4 ]+=alpha*tb[i+4 ];
		ta[i+5 ]+=alpha*tb[i+5 ];
		ta[i+6 ]+=alpha*tb[i+6 ];
		ta[i+7 ]+=alpha*tb[i+7 ];
		ta[i+8 ]+=alpha*tb[i+8 ];
		ta[i+9 ]+=alpha*tb[i+9 ];
		ta[i+10 ]+=alpha*tb[i+10 ];
		ta[i+11 ]+=alpha*tb[i+11 ];
		ta[i+12 ]+=alpha*tb[i+12 ];
		ta[i+13 ]+=alpha*tb[i+13 ];
		ta[i+14 ]+=alpha*tb[i+14 ];
		ta[i+15 ]+=alpha*tb[i+15 ];
	}
for(     ;i<n;++i){ 	ta[i+0 ]+=alpha*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__util_set_array_to_converted_integer(void *p, rsb_type_t typecode, const rsb_nnz_idx_t n, const rsb_nnz_idx_t incp, const rsb_int v)
{
	/*!
	 * */
	size_t i;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	double*tp = p; const double tv = (double)v;
	{
for(i=0;i+15<n;i+=16){
	tp[((i+0 )*incp)] = tv;
		tp[((i+1 )*incp)] = tv;
		tp[((i+2 )*incp)] = tv;
		tp[((i+3 )*incp)] = tv;
		tp[((i+4 )*incp)] = tv;
		tp[((i+5 )*incp)] = tv;
		tp[((i+6 )*incp)] = tv;
		tp[((i+7 )*incp)] = tv;
		tp[((i+8 )*incp)] = tv;
		tp[((i+9 )*incp)] = tv;
		tp[((i+10 )*incp)] = tv;
		tp[((i+11 )*incp)] = tv;
		tp[((i+12 )*incp)] = tv;
		tp[((i+13 )*incp)] = tv;
		tp[((i+14 )*incp)] = tv;
		tp[((i+15 )*incp)] = tv;
	}
for(     ;i<n;++i){ 	tp[((i+0 )*incp)] = tv;
	 }
}
; 
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	float*tp = p; const float tv = (float)v;
	{
for(i=0;i+15<n;i+=16){
	tp[((i+0 )*incp)] = tv;
		tp[((i+1 )*incp)] = tv;
		tp[((i+2 )*incp)] = tv;
		tp[((i+3 )*incp)] = tv;
		tp[((i+4 )*incp)] = tv;
		tp[((i+5 )*incp)] = tv;
		tp[((i+6 )*incp)] = tv;
		tp[((i+7 )*incp)] = tv;
		tp[((i+8 )*incp)] = tv;
		tp[((i+9 )*incp)] = tv;
		tp[((i+10 )*incp)] = tv;
		tp[((i+11 )*incp)] = tv;
		tp[((i+12 )*incp)] = tv;
		tp[((i+13 )*incp)] = tv;
		tp[((i+14 )*incp)] = tv;
		tp[((i+15 )*incp)] = tv;
	}
for(     ;i<n;++i){ 	tp[((i+0 )*incp)] = tv;
	 }
}
; 
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	float complex*tp = p; const float complex tv = (float complex)v;
	{
for(i=0;i+15<n;i+=16){
	tp[((i+0 )*incp)] = tv;
		tp[((i+1 )*incp)] = tv;
		tp[((i+2 )*incp)] = tv;
		tp[((i+3 )*incp)] = tv;
		tp[((i+4 )*incp)] = tv;
		tp[((i+5 )*incp)] = tv;
		tp[((i+6 )*incp)] = tv;
		tp[((i+7 )*incp)] = tv;
		tp[((i+8 )*incp)] = tv;
		tp[((i+9 )*incp)] = tv;
		tp[((i+10 )*incp)] = tv;
		tp[((i+11 )*incp)] = tv;
		tp[((i+12 )*incp)] = tv;
		tp[((i+13 )*incp)] = tv;
		tp[((i+14 )*incp)] = tv;
		tp[((i+15 )*incp)] = tv;
	}
for(     ;i<n;++i){ 	tp[((i+0 )*incp)] = tv;
	 }
}
; 
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	double complex*tp = p; const double complex tv = (double complex)v;
	{
for(i=0;i+15<n;i+=16){
	tp[((i+0 )*incp)] = tv;
		tp[((i+1 )*incp)] = tv;
		tp[((i+2 )*incp)] = tv;
		tp[((i+3 )*incp)] = tv;
		tp[((i+4 )*incp)] = tv;
		tp[((i+5 )*incp)] = tv;
		tp[((i+6 )*incp)] = tv;
		tp[((i+7 )*incp)] = tv;
		tp[((i+8 )*incp)] = tv;
		tp[((i+9 )*incp)] = tv;
		tp[((i+10 )*incp)] = tv;
		tp[((i+11 )*incp)] = tv;
		tp[((i+12 )*incp)] = tv;
		tp[((i+13 )*incp)] = tv;
		tp[((i+14 )*incp)] = tv;
		tp[((i+15 )*incp)] = tv;
	}
for(     ;i<n;++i){ 	tp[((i+0 )*incp)] = tv;
	 }
}
; 
	}
	else
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__vectors_left_sum_reduce_and_zero(void * d, void * s, const rsb_type_t typecode, const size_t n, const size_t incd, const size_t off)
{
	/*!
	 * d[off:off+n-1] <- d[off:off+n-1] + s[off:off+n-1] 
	 * s[off:off+n-1] <- 0
         *
	 * \param array	an array pointer
	 * \param typecode	a valid type code
	 * \param incd	the stride of d
	 * \param off offset in the vectors
	 * \return \rsberrcodemsg
	 * */
	size_t i;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	double*td = d,*ts = s;
	{
for(i=0;i+15<n;i+=16){
	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((double)(0));
		td[(off+i+1 )*incd]+=ts[(off+i+1 )];
	ts[(off+i+1 )] = ((double)(0));
		td[(off+i+2 )*incd]+=ts[(off+i+2 )];
	ts[(off+i+2 )] = ((double)(0));
		td[(off+i+3 )*incd]+=ts[(off+i+3 )];
	ts[(off+i+3 )] = ((double)(0));
		td[(off+i+4 )*incd]+=ts[(off+i+4 )];
	ts[(off+i+4 )] = ((double)(0));
		td[(off+i+5 )*incd]+=ts[(off+i+5 )];
	ts[(off+i+5 )] = ((double)(0));
		td[(off+i+6 )*incd]+=ts[(off+i+6 )];
	ts[(off+i+6 )] = ((double)(0));
		td[(off+i+7 )*incd]+=ts[(off+i+7 )];
	ts[(off+i+7 )] = ((double)(0));
		td[(off+i+8 )*incd]+=ts[(off+i+8 )];
	ts[(off+i+8 )] = ((double)(0));
		td[(off+i+9 )*incd]+=ts[(off+i+9 )];
	ts[(off+i+9 )] = ((double)(0));
		td[(off+i+10 )*incd]+=ts[(off+i+10 )];
	ts[(off+i+10 )] = ((double)(0));
		td[(off+i+11 )*incd]+=ts[(off+i+11 )];
	ts[(off+i+11 )] = ((double)(0));
		td[(off+i+12 )*incd]+=ts[(off+i+12 )];
	ts[(off+i+12 )] = ((double)(0));
		td[(off+i+13 )*incd]+=ts[(off+i+13 )];
	ts[(off+i+13 )] = ((double)(0));
		td[(off+i+14 )*incd]+=ts[(off+i+14 )];
	ts[(off+i+14 )] = ((double)(0));
		td[(off+i+15 )*incd]+=ts[(off+i+15 )];
	ts[(off+i+15 )] = ((double)(0));
	}
for(     ;i<n;++i){ 	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((double)(0));
	 }
}
; 
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	float*td = d,*ts = s;
	{
for(i=0;i+15<n;i+=16){
	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((float)(0));
		td[(off+i+1 )*incd]+=ts[(off+i+1 )];
	ts[(off+i+1 )] = ((float)(0));
		td[(off+i+2 )*incd]+=ts[(off+i+2 )];
	ts[(off+i+2 )] = ((float)(0));
		td[(off+i+3 )*incd]+=ts[(off+i+3 )];
	ts[(off+i+3 )] = ((float)(0));
		td[(off+i+4 )*incd]+=ts[(off+i+4 )];
	ts[(off+i+4 )] = ((float)(0));
		td[(off+i+5 )*incd]+=ts[(off+i+5 )];
	ts[(off+i+5 )] = ((float)(0));
		td[(off+i+6 )*incd]+=ts[(off+i+6 )];
	ts[(off+i+6 )] = ((float)(0));
		td[(off+i+7 )*incd]+=ts[(off+i+7 )];
	ts[(off+i+7 )] = ((float)(0));
		td[(off+i+8 )*incd]+=ts[(off+i+8 )];
	ts[(off+i+8 )] = ((float)(0));
		td[(off+i+9 )*incd]+=ts[(off+i+9 )];
	ts[(off+i+9 )] = ((float)(0));
		td[(off+i+10 )*incd]+=ts[(off+i+10 )];
	ts[(off+i+10 )] = ((float)(0));
		td[(off+i+11 )*incd]+=ts[(off+i+11 )];
	ts[(off+i+11 )] = ((float)(0));
		td[(off+i+12 )*incd]+=ts[(off+i+12 )];
	ts[(off+i+12 )] = ((float)(0));
		td[(off+i+13 )*incd]+=ts[(off+i+13 )];
	ts[(off+i+13 )] = ((float)(0));
		td[(off+i+14 )*incd]+=ts[(off+i+14 )];
	ts[(off+i+14 )] = ((float)(0));
		td[(off+i+15 )*incd]+=ts[(off+i+15 )];
	ts[(off+i+15 )] = ((float)(0));
	}
for(     ;i<n;++i){ 	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((float)(0));
	 }
}
; 
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	float complex*td = d,*ts = s;
	{
for(i=0;i+15<n;i+=16){
	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((float complex)(0));
		td[(off+i+1 )*incd]+=ts[(off+i+1 )];
	ts[(off+i+1 )] = ((float complex)(0));
		td[(off+i+2 )*incd]+=ts[(off+i+2 )];
	ts[(off+i+2 )] = ((float complex)(0));
		td[(off+i+3 )*incd]+=ts[(off+i+3 )];
	ts[(off+i+3 )] = ((float complex)(0));
		td[(off+i+4 )*incd]+=ts[(off+i+4 )];
	ts[(off+i+4 )] = ((float complex)(0));
		td[(off+i+5 )*incd]+=ts[(off+i+5 )];
	ts[(off+i+5 )] = ((float complex)(0));
		td[(off+i+6 )*incd]+=ts[(off+i+6 )];
	ts[(off+i+6 )] = ((float complex)(0));
		td[(off+i+7 )*incd]+=ts[(off+i+7 )];
	ts[(off+i+7 )] = ((float complex)(0));
		td[(off+i+8 )*incd]+=ts[(off+i+8 )];
	ts[(off+i+8 )] = ((float complex)(0));
		td[(off+i+9 )*incd]+=ts[(off+i+9 )];
	ts[(off+i+9 )] = ((float complex)(0));
		td[(off+i+10 )*incd]+=ts[(off+i+10 )];
	ts[(off+i+10 )] = ((float complex)(0));
		td[(off+i+11 )*incd]+=ts[(off+i+11 )];
	ts[(off+i+11 )] = ((float complex)(0));
		td[(off+i+12 )*incd]+=ts[(off+i+12 )];
	ts[(off+i+12 )] = ((float complex)(0));
		td[(off+i+13 )*incd]+=ts[(off+i+13 )];
	ts[(off+i+13 )] = ((float complex)(0));
		td[(off+i+14 )*incd]+=ts[(off+i+14 )];
	ts[(off+i+14 )] = ((float complex)(0));
		td[(off+i+15 )*incd]+=ts[(off+i+15 )];
	ts[(off+i+15 )] = ((float complex)(0));
	}
for(     ;i<n;++i){ 	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((float complex)(0));
	 }
}
; 
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	double complex*td = d,*ts = s;
	{
for(i=0;i+15<n;i+=16){
	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((double complex)(0));
		td[(off+i+1 )*incd]+=ts[(off+i+1 )];
	ts[(off+i+1 )] = ((double complex)(0));
		td[(off+i+2 )*incd]+=ts[(off+i+2 )];
	ts[(off+i+2 )] = ((double complex)(0));
		td[(off+i+3 )*incd]+=ts[(off+i+3 )];
	ts[(off+i+3 )] = ((double complex)(0));
		td[(off+i+4 )*incd]+=ts[(off+i+4 )];
	ts[(off+i+4 )] = ((double complex)(0));
		td[(off+i+5 )*incd]+=ts[(off+i+5 )];
	ts[(off+i+5 )] = ((double complex)(0));
		td[(off+i+6 )*incd]+=ts[(off+i+6 )];
	ts[(off+i+6 )] = ((double complex)(0));
		td[(off+i+7 )*incd]+=ts[(off+i+7 )];
	ts[(off+i+7 )] = ((double complex)(0));
		td[(off+i+8 )*incd]+=ts[(off+i+8 )];
	ts[(off+i+8 )] = ((double complex)(0));
		td[(off+i+9 )*incd]+=ts[(off+i+9 )];
	ts[(off+i+9 )] = ((double complex)(0));
		td[(off+i+10 )*incd]+=ts[(off+i+10 )];
	ts[(off+i+10 )] = ((double complex)(0));
		td[(off+i+11 )*incd]+=ts[(off+i+11 )];
	ts[(off+i+11 )] = ((double complex)(0));
		td[(off+i+12 )*incd]+=ts[(off+i+12 )];
	ts[(off+i+12 )] = ((double complex)(0));
		td[(off+i+13 )*incd]+=ts[(off+i+13 )];
	ts[(off+i+13 )] = ((double complex)(0));
		td[(off+i+14 )*incd]+=ts[(off+i+14 )];
	ts[(off+i+14 )] = ((double complex)(0));
		td[(off+i+15 )*incd]+=ts[(off+i+15 )];
	ts[(off+i+15 )] = ((double complex)(0));
	}
for(     ;i<n;++i){ 	td[(off+i+0 )*incd]+=ts[(off+i+0 )];
	ts[(off+i+0 )] = ((double complex)(0));
	 }
}
; 
	}
	else
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}


static rsb_err_t rsb_alpha_sum_strided(void * a, const void * b, const void * alphap, rsb_type_t type, size_t n, int inca, int incb)
{
	/*!
	 * a <- a + alpha * b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see daxpy in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(inca == 1 && incb == 1)
		return rsb_alpha_sum(a,b,alphap,type,n);
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double alpha = alphap ? *(double*)alphap : ((double)(1.0));
	double*ta = a; const double*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
		ta[inca*(i+1 )]+=alpha*tb[incb*(i+1 )];
		ta[inca*(i+2 )]+=alpha*tb[incb*(i+2 )];
		ta[inca*(i+3 )]+=alpha*tb[incb*(i+3 )];
		ta[inca*(i+4 )]+=alpha*tb[incb*(i+4 )];
		ta[inca*(i+5 )]+=alpha*tb[incb*(i+5 )];
		ta[inca*(i+6 )]+=alpha*tb[incb*(i+6 )];
		ta[inca*(i+7 )]+=alpha*tb[incb*(i+7 )];
		ta[inca*(i+8 )]+=alpha*tb[incb*(i+8 )];
		ta[inca*(i+9 )]+=alpha*tb[incb*(i+9 )];
		ta[inca*(i+10 )]+=alpha*tb[incb*(i+10 )];
		ta[inca*(i+11 )]+=alpha*tb[incb*(i+11 )];
		ta[inca*(i+12 )]+=alpha*tb[incb*(i+12 )];
		ta[inca*(i+13 )]+=alpha*tb[incb*(i+13 )];
		ta[inca*(i+14 )]+=alpha*tb[incb*(i+14 )];
		ta[inca*(i+15 )]+=alpha*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float alpha = alphap ? *(float*)alphap : ((float)(1.0));
	float*ta = a; const float*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
		ta[inca*(i+1 )]+=alpha*tb[incb*(i+1 )];
		ta[inca*(i+2 )]+=alpha*tb[incb*(i+2 )];
		ta[inca*(i+3 )]+=alpha*tb[incb*(i+3 )];
		ta[inca*(i+4 )]+=alpha*tb[incb*(i+4 )];
		ta[inca*(i+5 )]+=alpha*tb[incb*(i+5 )];
		ta[inca*(i+6 )]+=alpha*tb[incb*(i+6 )];
		ta[inca*(i+7 )]+=alpha*tb[incb*(i+7 )];
		ta[inca*(i+8 )]+=alpha*tb[incb*(i+8 )];
		ta[inca*(i+9 )]+=alpha*tb[incb*(i+9 )];
		ta[inca*(i+10 )]+=alpha*tb[incb*(i+10 )];
		ta[inca*(i+11 )]+=alpha*tb[incb*(i+11 )];
		ta[inca*(i+12 )]+=alpha*tb[incb*(i+12 )];
		ta[inca*(i+13 )]+=alpha*tb[incb*(i+13 )];
		ta[inca*(i+14 )]+=alpha*tb[incb*(i+14 )];
		ta[inca*(i+15 )]+=alpha*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex alpha = alphap ? *(float complex*)alphap : ((float complex)(1.0));
	float complex*ta = a; const float complex*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
		ta[inca*(i+1 )]+=alpha*tb[incb*(i+1 )];
		ta[inca*(i+2 )]+=alpha*tb[incb*(i+2 )];
		ta[inca*(i+3 )]+=alpha*tb[incb*(i+3 )];
		ta[inca*(i+4 )]+=alpha*tb[incb*(i+4 )];
		ta[inca*(i+5 )]+=alpha*tb[incb*(i+5 )];
		ta[inca*(i+6 )]+=alpha*tb[incb*(i+6 )];
		ta[inca*(i+7 )]+=alpha*tb[incb*(i+7 )];
		ta[inca*(i+8 )]+=alpha*tb[incb*(i+8 )];
		ta[inca*(i+9 )]+=alpha*tb[incb*(i+9 )];
		ta[inca*(i+10 )]+=alpha*tb[incb*(i+10 )];
		ta[inca*(i+11 )]+=alpha*tb[incb*(i+11 )];
		ta[inca*(i+12 )]+=alpha*tb[incb*(i+12 )];
		ta[inca*(i+13 )]+=alpha*tb[incb*(i+13 )];
		ta[inca*(i+14 )]+=alpha*tb[incb*(i+14 )];
		ta[inca*(i+15 )]+=alpha*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex alpha = alphap ? *(double complex*)alphap : ((double complex)(1.0));
	double complex*ta = a; const double complex*tb = b;
	{
for(i=0;i+15<n;i+=16){
	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
		ta[inca*(i+1 )]+=alpha*tb[incb*(i+1 )];
		ta[inca*(i+2 )]+=alpha*tb[incb*(i+2 )];
		ta[inca*(i+3 )]+=alpha*tb[incb*(i+3 )];
		ta[inca*(i+4 )]+=alpha*tb[incb*(i+4 )];
		ta[inca*(i+5 )]+=alpha*tb[incb*(i+5 )];
		ta[inca*(i+6 )]+=alpha*tb[incb*(i+6 )];
		ta[inca*(i+7 )]+=alpha*tb[incb*(i+7 )];
		ta[inca*(i+8 )]+=alpha*tb[incb*(i+8 )];
		ta[inca*(i+9 )]+=alpha*tb[incb*(i+9 )];
		ta[inca*(i+10 )]+=alpha*tb[incb*(i+10 )];
		ta[inca*(i+11 )]+=alpha*tb[incb*(i+11 )];
		ta[inca*(i+12 )]+=alpha*tb[incb*(i+12 )];
		ta[inca*(i+13 )]+=alpha*tb[incb*(i+13 )];
		ta[inca*(i+14 )]+=alpha*tb[incb*(i+14 )];
		ta[inca*(i+15 )]+=alpha*tb[incb*(i+15 )];
	}
for(     ;i<n;++i){ 	ta[inca*(i+0 )]+=alpha*tb[incb*(i+0 )];
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__cblas_Xaxpy(rsb_type_t type, size_t n, const void * alphap, const void * x, const int incx, void * y, const int incy)
{
	/*!
	 * y <- y + alpha * x
         */
	return rsb_alpha_sum_strided(y,x,alphap,type,n,incy,incx);
}

rsb_err_t rsb__vector_mult(const void * a, const void * b, void * c, rsb_type_t type, size_t n)
{
	/*!
	 * c <- a*b
	 * It is allowed to give c == a or c == b or a == b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * 
	 * FIXME : useless ?
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double*ta = a; const double*tb = b; double*tc = c;
	{
for(i=0;i+15<n;i+=16){
	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
		tc[i+1 ] = ta[i+1 ]*tb[i+1 ];
		tc[i+2 ] = ta[i+2 ]*tb[i+2 ];
		tc[i+3 ] = ta[i+3 ]*tb[i+3 ];
		tc[i+4 ] = ta[i+4 ]*tb[i+4 ];
		tc[i+5 ] = ta[i+5 ]*tb[i+5 ];
		tc[i+6 ] = ta[i+6 ]*tb[i+6 ];
		tc[i+7 ] = ta[i+7 ]*tb[i+7 ];
		tc[i+8 ] = ta[i+8 ]*tb[i+8 ];
		tc[i+9 ] = ta[i+9 ]*tb[i+9 ];
		tc[i+10 ] = ta[i+10 ]*tb[i+10 ];
		tc[i+11 ] = ta[i+11 ]*tb[i+11 ];
		tc[i+12 ] = ta[i+12 ]*tb[i+12 ];
		tc[i+13 ] = ta[i+13 ]*tb[i+13 ];
		tc[i+14 ] = ta[i+14 ]*tb[i+14 ];
		tc[i+15 ] = ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float*ta = a; const float*tb = b; float*tc = c;
	{
for(i=0;i+15<n;i+=16){
	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
		tc[i+1 ] = ta[i+1 ]*tb[i+1 ];
		tc[i+2 ] = ta[i+2 ]*tb[i+2 ];
		tc[i+3 ] = ta[i+3 ]*tb[i+3 ];
		tc[i+4 ] = ta[i+4 ]*tb[i+4 ];
		tc[i+5 ] = ta[i+5 ]*tb[i+5 ];
		tc[i+6 ] = ta[i+6 ]*tb[i+6 ];
		tc[i+7 ] = ta[i+7 ]*tb[i+7 ];
		tc[i+8 ] = ta[i+8 ]*tb[i+8 ];
		tc[i+9 ] = ta[i+9 ]*tb[i+9 ];
		tc[i+10 ] = ta[i+10 ]*tb[i+10 ];
		tc[i+11 ] = ta[i+11 ]*tb[i+11 ];
		tc[i+12 ] = ta[i+12 ]*tb[i+12 ];
		tc[i+13 ] = ta[i+13 ]*tb[i+13 ];
		tc[i+14 ] = ta[i+14 ]*tb[i+14 ];
		tc[i+15 ] = ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex*ta = a; const float complex*tb = b; float complex*tc = c;
	{
for(i=0;i+15<n;i+=16){
	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
		tc[i+1 ] = ta[i+1 ]*tb[i+1 ];
		tc[i+2 ] = ta[i+2 ]*tb[i+2 ];
		tc[i+3 ] = ta[i+3 ]*tb[i+3 ];
		tc[i+4 ] = ta[i+4 ]*tb[i+4 ];
		tc[i+5 ] = ta[i+5 ]*tb[i+5 ];
		tc[i+6 ] = ta[i+6 ]*tb[i+6 ];
		tc[i+7 ] = ta[i+7 ]*tb[i+7 ];
		tc[i+8 ] = ta[i+8 ]*tb[i+8 ];
		tc[i+9 ] = ta[i+9 ]*tb[i+9 ];
		tc[i+10 ] = ta[i+10 ]*tb[i+10 ];
		tc[i+11 ] = ta[i+11 ]*tb[i+11 ];
		tc[i+12 ] = ta[i+12 ]*tb[i+12 ];
		tc[i+13 ] = ta[i+13 ]*tb[i+13 ];
		tc[i+14 ] = ta[i+14 ]*tb[i+14 ];
		tc[i+15 ] = ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex*ta = a; const double complex*tb = b; double complex*tc = c;
	{
for(i=0;i+15<n;i+=16){
	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
		tc[i+1 ] = ta[i+1 ]*tb[i+1 ];
		tc[i+2 ] = ta[i+2 ]*tb[i+2 ];
		tc[i+3 ] = ta[i+3 ]*tb[i+3 ];
		tc[i+4 ] = ta[i+4 ]*tb[i+4 ];
		tc[i+5 ] = ta[i+5 ]*tb[i+5 ];
		tc[i+6 ] = ta[i+6 ]*tb[i+6 ];
		tc[i+7 ] = ta[i+7 ]*tb[i+7 ];
		tc[i+8 ] = ta[i+8 ]*tb[i+8 ];
		tc[i+9 ] = ta[i+9 ]*tb[i+9 ];
		tc[i+10 ] = ta[i+10 ]*tb[i+10 ];
		tc[i+11 ] = ta[i+11 ]*tb[i+11 ];
		tc[i+12 ] = ta[i+12 ]*tb[i+12 ];
		tc[i+13 ] = ta[i+13 ]*tb[i+13 ];
		tc[i+14 ] = ta[i+14 ]*tb[i+14 ];
		tc[i+15 ] = ta[i+15 ]*tb[i+15 ];
	}
for(     ;i<n;++i){ 	tc[i+0 ] = ta[i+0 ]*tb[i+0 ];
	 }
}
; 
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__xcopy(void * a, const void * b, rsb_nnz_idx_t toi, rsb_nnz_idx_t foi, rsb_nnz_idx_t n,size_t el_size)
{
	/*!
	 * a[toi:toi+n] <- b[foi:foi+n] 
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 * */
	rsb__memcpy(((rsb_byte_t*)a)+el_size*toi,((const rsb_byte_t*)b)+el_size*foi,el_size*n);
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__do_are_similar_parametric(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy, int extra_decimals)
{
	/*!
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 *
	 * For cases like 1+0I differing from 1-0I .
	 * */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		rsb_nnz_idx_t i;
		const double *a = ap;
		const double *b = bp;
		const double threshold = 1e-6 * pow(10*((double)(1.0)),(double)extra_decimals);

		for(i=0;i<n;++i)
		{
			const double av = a[incx*(i)];
			const double bv = b[incy*(i)];
			if( av - bv )
			{
				const double aav = fabs(av);
				const double abv = fabs(bv);
				if( av && fabs((aav - abv)) / fabs(av) > threshold )
					goto differing;
				if( bv && fabs((aav - abv)) / fabs(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		rsb_nnz_idx_t i;
		const float *a = ap;
		const float *b = bp;
		const float threshold = 1e-5 * powf(10*((float)(1.0)),(float)extra_decimals);

		for(i=0;i<n;++i)
		{
			const float av = a[incx*(i)];
			const float bv = b[incy*(i)];
			if( av - bv )
			{
				const float aav = fabsf(av);
				const float abv = fabsf(bv);
				if( av && fabsf((aav - abv)) / fabsf(av) > threshold )
					goto differing;
				if( bv && fabsf((aav - abv)) / fabsf(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		rsb_nnz_idx_t i;
		const float complex *a = ap;
		const float complex *b = bp;
		const float threshold = 1e-4 * cpowf(10*((float complex)(1.0)),(float complex)extra_decimals);

		for(i=0;i<n;++i)
		{
			const float complex av = a[incx*(i)];
			const float complex bv = b[incy*(i)];
			if( av - bv )
			{
				const float complex aav = cabsf(av);
				const float complex abv = cabsf(bv);
				if( av && cabsf((aav - abv)) / cabsf(av) > threshold )
					goto differing;
				if( bv && cabsf((aav - abv)) / cabsf(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		rsb_nnz_idx_t i;
		const double complex *a = ap;
		const double complex *b = bp;
		const double threshold = 1e-6 * cpow(10*((double complex)(1.0)),(double complex)extra_decimals);

		for(i=0;i<n;++i)
		{
			const double complex av = a[incx*(i)];
			const double complex bv = b[incy*(i)];
			if( av - bv )
			{
				const double complex aav = cabs(av);
				const double complex abv = cabs(bv);
				if( av && cabs((aav - abv)) / cabs(av) > threshold )
					goto differing;
				if( bv && cabs((aav - abv)) / cabs(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
	return RSB_ERR_UNSUPPORTED_TYPE;
differing:
	return RSB_ERR_GENERIC_ERROR;
}

rsb_err_t rsb__do_are_similar(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
{
	/*!
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 *
	 * For cases like 1+0I differing from 1-0I .
	 * */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		rsb_nnz_idx_t i;
		const double *a = ap;
		const double *b = bp;
		const double threshold = 1e-6;

		for(i=0;i<n;++i)
		{
			const double av = a[incx*(i)];
			const double bv = b[incy*(i)];
			if( av - bv )
			{
				const double aav = fabs(av);
				const double abv = fabs(bv);
				if( av && fabs((aav - abv)) / fabs(av) > threshold )
					goto differing;
				if( bv && fabs((aav - abv)) / fabs(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		rsb_nnz_idx_t i;
		const float *a = ap;
		const float *b = bp;
		const float threshold = 1e-5;

		for(i=0;i<n;++i)
		{
			const float av = a[incx*(i)];
			const float bv = b[incy*(i)];
			if( av - bv )
			{
				const float aav = fabsf(av);
				const float abv = fabsf(bv);
				if( av && fabsf((aav - abv)) / fabsf(av) > threshold )
					goto differing;
				if( bv && fabsf((aav - abv)) / fabsf(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		rsb_nnz_idx_t i;
		const float complex *a = ap;
		const float complex *b = bp;
		const float threshold = 1e-4;

		for(i=0;i<n;++i)
		{
			const float complex av = a[incx*(i)];
			const float complex bv = b[incy*(i)];
			if( av - bv )
			{
				const float complex aav = cabsf(av);
				const float complex abv = cabsf(bv);
				if( av && cabsf((aav - abv)) / cabsf(av) > threshold )
					goto differing;
				if( bv && cabsf((aav - abv)) / cabsf(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		rsb_nnz_idx_t i;
		const double complex *a = ap;
		const double complex *b = bp;
		const double threshold = 1e-6;

		for(i=0;i<n;++i)
		{
			const double complex av = a[incx*(i)];
			const double complex bv = b[incy*(i)];
			if( av - bv )
			{
				const double complex aav = cabs(av);
				const double complex abv = cabs(bv);
				if( av && cabs((aav - abv)) / cabs(av) > threshold )
					goto differing;
				if( bv && cabs((aav - abv)) / cabs(bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
	return RSB_ERR_UNSUPPORTED_TYPE;
differing:
	return RSB_ERR_GENERIC_ERROR;
}

rsb_err_t rsb__do_are_same(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
{
	return rsb__do_are_similar_parametric(ap, bp, n, typecode, incx, incy, 0);
}

static rsb_err_t rsb__xcopy_strided_typed(void * a, const void * b, rsb_nnz_idx_t toi, rsb_nnz_idx_t foi, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
{
	/*!
	 * a[toi:toi+n] <- b[foi:foi+n] 
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 * */
	if(incx==1 && incy==1)
		return rsb__xcopy(a,b,toi,foi,n,RSB_SIZEOF(typecode));
	/* else */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	rsb_nnz_idx_t i;
	double *ap = a; const double *bp = b;
	ap+=toi;
	bp+=foi;
	{
for(i=0;i+15<n;i+=16){
	ap[(i+0 )*incx] = bp[(i+0 )*incy];
		ap[(i+1 )*incx] = bp[(i+1 )*incy];
		ap[(i+2 )*incx] = bp[(i+2 )*incy];
		ap[(i+3 )*incx] = bp[(i+3 )*incy];
		ap[(i+4 )*incx] = bp[(i+4 )*incy];
		ap[(i+5 )*incx] = bp[(i+5 )*incy];
		ap[(i+6 )*incx] = bp[(i+6 )*incy];
		ap[(i+7 )*incx] = bp[(i+7 )*incy];
		ap[(i+8 )*incx] = bp[(i+8 )*incy];
		ap[(i+9 )*incx] = bp[(i+9 )*incy];
		ap[(i+10 )*incx] = bp[(i+10 )*incy];
		ap[(i+11 )*incx] = bp[(i+11 )*incy];
		ap[(i+12 )*incx] = bp[(i+12 )*incy];
		ap[(i+13 )*incx] = bp[(i+13 )*incy];
		ap[(i+14 )*incx] = bp[(i+14 )*incy];
		ap[(i+15 )*incx] = bp[(i+15 )*incy];
	}
for(     ;i<n;++i){ 	ap[(i+0 )*incx] = bp[(i+0 )*incy];
	 }
}
; 
		return RSB_ERR_NO_ERROR;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	rsb_nnz_idx_t i;
	float *ap = a; const float *bp = b;
	ap+=toi;
	bp+=foi;
	{
for(i=0;i+15<n;i+=16){
	ap[(i+0 )*incx] = bp[(i+0 )*incy];
		ap[(i+1 )*incx] = bp[(i+1 )*incy];
		ap[(i+2 )*incx] = bp[(i+2 )*incy];
		ap[(i+3 )*incx] = bp[(i+3 )*incy];
		ap[(i+4 )*incx] = bp[(i+4 )*incy];
		ap[(i+5 )*incx] = bp[(i+5 )*incy];
		ap[(i+6 )*incx] = bp[(i+6 )*incy];
		ap[(i+7 )*incx] = bp[(i+7 )*incy];
		ap[(i+8 )*incx] = bp[(i+8 )*incy];
		ap[(i+9 )*incx] = bp[(i+9 )*incy];
		ap[(i+10 )*incx] = bp[(i+10 )*incy];
		ap[(i+11 )*incx] = bp[(i+11 )*incy];
		ap[(i+12 )*incx] = bp[(i+12 )*incy];
		ap[(i+13 )*incx] = bp[(i+13 )*incy];
		ap[(i+14 )*incx] = bp[(i+14 )*incy];
		ap[(i+15 )*incx] = bp[(i+15 )*incy];
	}
for(     ;i<n;++i){ 	ap[(i+0 )*incx] = bp[(i+0 )*incy];
	 }
}
; 
		return RSB_ERR_NO_ERROR;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	rsb_nnz_idx_t i;
	float complex *ap = a; const float complex *bp = b;
	ap+=toi;
	bp+=foi;
	{
for(i=0;i+15<n;i+=16){
	ap[(i+0 )*incx] = bp[(i+0 )*incy];
		ap[(i+1 )*incx] = bp[(i+1 )*incy];
		ap[(i+2 )*incx] = bp[(i+2 )*incy];
		ap[(i+3 )*incx] = bp[(i+3 )*incy];
		ap[(i+4 )*incx] = bp[(i+4 )*incy];
		ap[(i+5 )*incx] = bp[(i+5 )*incy];
		ap[(i+6 )*incx] = bp[(i+6 )*incy];
		ap[(i+7 )*incx] = bp[(i+7 )*incy];
		ap[(i+8 )*incx] = bp[(i+8 )*incy];
		ap[(i+9 )*incx] = bp[(i+9 )*incy];
		ap[(i+10 )*incx] = bp[(i+10 )*incy];
		ap[(i+11 )*incx] = bp[(i+11 )*incy];
		ap[(i+12 )*incx] = bp[(i+12 )*incy];
		ap[(i+13 )*incx] = bp[(i+13 )*incy];
		ap[(i+14 )*incx] = bp[(i+14 )*incy];
		ap[(i+15 )*incx] = bp[(i+15 )*incy];
	}
for(     ;i<n;++i){ 	ap[(i+0 )*incx] = bp[(i+0 )*incy];
	 }
}
; 
		return RSB_ERR_NO_ERROR;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	rsb_nnz_idx_t i;
	double complex *ap = a; const double complex *bp = b;
	ap+=toi;
	bp+=foi;
	{
for(i=0;i+15<n;i+=16){
	ap[(i+0 )*incx] = bp[(i+0 )*incy];
		ap[(i+1 )*incx] = bp[(i+1 )*incy];
		ap[(i+2 )*incx] = bp[(i+2 )*incy];
		ap[(i+3 )*incx] = bp[(i+3 )*incy];
		ap[(i+4 )*incx] = bp[(i+4 )*incy];
		ap[(i+5 )*incx] = bp[(i+5 )*incy];
		ap[(i+6 )*incx] = bp[(i+6 )*incy];
		ap[(i+7 )*incx] = bp[(i+7 )*incy];
		ap[(i+8 )*incx] = bp[(i+8 )*incy];
		ap[(i+9 )*incx] = bp[(i+9 )*incy];
		ap[(i+10 )*incx] = bp[(i+10 )*incy];
		ap[(i+11 )*incx] = bp[(i+11 )*incy];
		ap[(i+12 )*incx] = bp[(i+12 )*incy];
		ap[(i+13 )*incx] = bp[(i+13 )*incy];
		ap[(i+14 )*incx] = bp[(i+14 )*incy];
		ap[(i+15 )*incx] = bp[(i+15 )*incy];
	}
for(     ;i<n;++i){ 	ap[(i+0 )*incx] = bp[(i+0 )*incy];
	 }
}
; 
		return RSB_ERR_NO_ERROR;
	}
	else 
#endif
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__sqrt_of_sum_of_fabs_diffs(const void * a, const void * b, void *err, rsb_type_t type, size_t n)
{
	/*!
	 * Will compute the square root of the sum of the squares of the vectors elements differences.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * 
	 * FIXME
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double*ta = a; const double*tb = b;
	*((double*)err) = ((double)(0));
	{
for(i=0;i+15<n;i+=16){
	*((double*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
		*((double*)(err))+=(ta[i+1 ]-tb[i+1 ])*(ta[i+1 ]-tb[i+1 ]);
		*((double*)(err))+=(ta[i+2 ]-tb[i+2 ])*(ta[i+2 ]-tb[i+2 ]);
		*((double*)(err))+=(ta[i+3 ]-tb[i+3 ])*(ta[i+3 ]-tb[i+3 ]);
		*((double*)(err))+=(ta[i+4 ]-tb[i+4 ])*(ta[i+4 ]-tb[i+4 ]);
		*((double*)(err))+=(ta[i+5 ]-tb[i+5 ])*(ta[i+5 ]-tb[i+5 ]);
		*((double*)(err))+=(ta[i+6 ]-tb[i+6 ])*(ta[i+6 ]-tb[i+6 ]);
		*((double*)(err))+=(ta[i+7 ]-tb[i+7 ])*(ta[i+7 ]-tb[i+7 ]);
		*((double*)(err))+=(ta[i+8 ]-tb[i+8 ])*(ta[i+8 ]-tb[i+8 ]);
		*((double*)(err))+=(ta[i+9 ]-tb[i+9 ])*(ta[i+9 ]-tb[i+9 ]);
		*((double*)(err))+=(ta[i+10 ]-tb[i+10 ])*(ta[i+10 ]-tb[i+10 ]);
		*((double*)(err))+=(ta[i+11 ]-tb[i+11 ])*(ta[i+11 ]-tb[i+11 ]);
		*((double*)(err))+=(ta[i+12 ]-tb[i+12 ])*(ta[i+12 ]-tb[i+12 ]);
		*((double*)(err))+=(ta[i+13 ]-tb[i+13 ])*(ta[i+13 ]-tb[i+13 ]);
		*((double*)(err))+=(ta[i+14 ]-tb[i+14 ])*(ta[i+14 ]-tb[i+14 ]);
		*((double*)(err))+=(ta[i+15 ]-tb[i+15 ])*(ta[i+15 ]-tb[i+15 ]);
	}
for(     ;i<n;++i){ 	*((double*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
	 }
}
; 
	*((double*)err) = sqrt((*((double*)err)));
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float*ta = a; const float*tb = b;
	*((float*)err) = ((float)(0));
	{
for(i=0;i+15<n;i+=16){
	*((float*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
		*((float*)(err))+=(ta[i+1 ]-tb[i+1 ])*(ta[i+1 ]-tb[i+1 ]);
		*((float*)(err))+=(ta[i+2 ]-tb[i+2 ])*(ta[i+2 ]-tb[i+2 ]);
		*((float*)(err))+=(ta[i+3 ]-tb[i+3 ])*(ta[i+3 ]-tb[i+3 ]);
		*((float*)(err))+=(ta[i+4 ]-tb[i+4 ])*(ta[i+4 ]-tb[i+4 ]);
		*((float*)(err))+=(ta[i+5 ]-tb[i+5 ])*(ta[i+5 ]-tb[i+5 ]);
		*((float*)(err))+=(ta[i+6 ]-tb[i+6 ])*(ta[i+6 ]-tb[i+6 ]);
		*((float*)(err))+=(ta[i+7 ]-tb[i+7 ])*(ta[i+7 ]-tb[i+7 ]);
		*((float*)(err))+=(ta[i+8 ]-tb[i+8 ])*(ta[i+8 ]-tb[i+8 ]);
		*((float*)(err))+=(ta[i+9 ]-tb[i+9 ])*(ta[i+9 ]-tb[i+9 ]);
		*((float*)(err))+=(ta[i+10 ]-tb[i+10 ])*(ta[i+10 ]-tb[i+10 ]);
		*((float*)(err))+=(ta[i+11 ]-tb[i+11 ])*(ta[i+11 ]-tb[i+11 ]);
		*((float*)(err))+=(ta[i+12 ]-tb[i+12 ])*(ta[i+12 ]-tb[i+12 ]);
		*((float*)(err))+=(ta[i+13 ]-tb[i+13 ])*(ta[i+13 ]-tb[i+13 ]);
		*((float*)(err))+=(ta[i+14 ]-tb[i+14 ])*(ta[i+14 ]-tb[i+14 ]);
		*((float*)(err))+=(ta[i+15 ]-tb[i+15 ])*(ta[i+15 ]-tb[i+15 ]);
	}
for(     ;i<n;++i){ 	*((float*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
	 }
}
; 
	*((float*)err) = sqrtf((*((float*)err)));
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex*ta = a; const float complex*tb = b;
	*((float complex*)err) = ((float complex)(0));
	{
for(i=0;i+15<n;i+=16){
	*((float complex*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
		*((float complex*)(err))+=(ta[i+1 ]-tb[i+1 ])*(ta[i+1 ]-tb[i+1 ]);
		*((float complex*)(err))+=(ta[i+2 ]-tb[i+2 ])*(ta[i+2 ]-tb[i+2 ]);
		*((float complex*)(err))+=(ta[i+3 ]-tb[i+3 ])*(ta[i+3 ]-tb[i+3 ]);
		*((float complex*)(err))+=(ta[i+4 ]-tb[i+4 ])*(ta[i+4 ]-tb[i+4 ]);
		*((float complex*)(err))+=(ta[i+5 ]-tb[i+5 ])*(ta[i+5 ]-tb[i+5 ]);
		*((float complex*)(err))+=(ta[i+6 ]-tb[i+6 ])*(ta[i+6 ]-tb[i+6 ]);
		*((float complex*)(err))+=(ta[i+7 ]-tb[i+7 ])*(ta[i+7 ]-tb[i+7 ]);
		*((float complex*)(err))+=(ta[i+8 ]-tb[i+8 ])*(ta[i+8 ]-tb[i+8 ]);
		*((float complex*)(err))+=(ta[i+9 ]-tb[i+9 ])*(ta[i+9 ]-tb[i+9 ]);
		*((float complex*)(err))+=(ta[i+10 ]-tb[i+10 ])*(ta[i+10 ]-tb[i+10 ]);
		*((float complex*)(err))+=(ta[i+11 ]-tb[i+11 ])*(ta[i+11 ]-tb[i+11 ]);
		*((float complex*)(err))+=(ta[i+12 ]-tb[i+12 ])*(ta[i+12 ]-tb[i+12 ]);
		*((float complex*)(err))+=(ta[i+13 ]-tb[i+13 ])*(ta[i+13 ]-tb[i+13 ]);
		*((float complex*)(err))+=(ta[i+14 ]-tb[i+14 ])*(ta[i+14 ]-tb[i+14 ]);
		*((float complex*)(err))+=(ta[i+15 ]-tb[i+15 ])*(ta[i+15 ]-tb[i+15 ]);
	}
for(     ;i<n;++i){ 	*((float complex*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
	 }
}
; 
	*((float complex*)err) = csqrtf((*((float complex*)err)));
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex*ta = a; const double complex*tb = b;
	*((double complex*)err) = ((double complex)(0));
	{
for(i=0;i+15<n;i+=16){
	*((double complex*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
		*((double complex*)(err))+=(ta[i+1 ]-tb[i+1 ])*(ta[i+1 ]-tb[i+1 ]);
		*((double complex*)(err))+=(ta[i+2 ]-tb[i+2 ])*(ta[i+2 ]-tb[i+2 ]);
		*((double complex*)(err))+=(ta[i+3 ]-tb[i+3 ])*(ta[i+3 ]-tb[i+3 ]);
		*((double complex*)(err))+=(ta[i+4 ]-tb[i+4 ])*(ta[i+4 ]-tb[i+4 ]);
		*((double complex*)(err))+=(ta[i+5 ]-tb[i+5 ])*(ta[i+5 ]-tb[i+5 ]);
		*((double complex*)(err))+=(ta[i+6 ]-tb[i+6 ])*(ta[i+6 ]-tb[i+6 ]);
		*((double complex*)(err))+=(ta[i+7 ]-tb[i+7 ])*(ta[i+7 ]-tb[i+7 ]);
		*((double complex*)(err))+=(ta[i+8 ]-tb[i+8 ])*(ta[i+8 ]-tb[i+8 ]);
		*((double complex*)(err))+=(ta[i+9 ]-tb[i+9 ])*(ta[i+9 ]-tb[i+9 ]);
		*((double complex*)(err))+=(ta[i+10 ]-tb[i+10 ])*(ta[i+10 ]-tb[i+10 ]);
		*((double complex*)(err))+=(ta[i+11 ]-tb[i+11 ])*(ta[i+11 ]-tb[i+11 ]);
		*((double complex*)(err))+=(ta[i+12 ]-tb[i+12 ])*(ta[i+12 ]-tb[i+12 ]);
		*((double complex*)(err))+=(ta[i+13 ]-tb[i+13 ])*(ta[i+13 ]-tb[i+13 ]);
		*((double complex*)(err))+=(ta[i+14 ]-tb[i+14 ])*(ta[i+14 ]-tb[i+14 ]);
		*((double complex*)(err))+=(ta[i+15 ]-tb[i+15 ])*(ta[i+15 ]-tb[i+15 ]);
	}
for(     ;i<n;++i){ 	*((double complex*)(err))+=(ta[i+0 ]-tb[i+0 ])*(ta[i+0 ]-tb[i+0 ]);
	 }
}
; 
	*((double complex*)err) = csqrt((*((double complex*)err)));
	}
	else 
#endif
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__fill_with_increasing_values(void * array, rsb_type_t type, size_t n)
{
	/*!
	 * \ingroup gr_vec
	 * FIXME : document me
	 * starts with one.
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{ 
	double*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = (const double)(i+0 +1);ta[i+1 ] = (const double)(i+1 +1);ta[i+2 ] = (const double)(i+2 +1);ta[i+3 ] = (const double)(i+3 +1);ta[i+4 ] = (const double)(i+4 +1);ta[i+5 ] = (const double)(i+5 +1);ta[i+6 ] = (const double)(i+6 +1);ta[i+7 ] = (const double)(i+7 +1);ta[i+8 ] = (const double)(i+8 +1);ta[i+9 ] = (const double)(i+9 +1);ta[i+10 ] = (const double)(i+10 +1);ta[i+11 ] = (const double)(i+11 +1);ta[i+12 ] = (const double)(i+12 +1);ta[i+13 ] = (const double)(i+13 +1);ta[i+14 ] = (const double)(i+14 +1);ta[i+15 ] = (const double)(i+15 +1);}
for(     ;i<n;++i){ ta[i+0 ] = (const double)(i+0 +1); }
}

	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{ 
	float*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = (const float)(i+0 +1);ta[i+1 ] = (const float)(i+1 +1);ta[i+2 ] = (const float)(i+2 +1);ta[i+3 ] = (const float)(i+3 +1);ta[i+4 ] = (const float)(i+4 +1);ta[i+5 ] = (const float)(i+5 +1);ta[i+6 ] = (const float)(i+6 +1);ta[i+7 ] = (const float)(i+7 +1);ta[i+8 ] = (const float)(i+8 +1);ta[i+9 ] = (const float)(i+9 +1);ta[i+10 ] = (const float)(i+10 +1);ta[i+11 ] = (const float)(i+11 +1);ta[i+12 ] = (const float)(i+12 +1);ta[i+13 ] = (const float)(i+13 +1);ta[i+14 ] = (const float)(i+14 +1);ta[i+15 ] = (const float)(i+15 +1);}
for(     ;i<n;++i){ ta[i+0 ] = (const float)(i+0 +1); }
}

	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{ 
	float complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = (const float complex)(i+0 +1);ta[i+1 ] = (const float complex)(i+1 +1);ta[i+2 ] = (const float complex)(i+2 +1);ta[i+3 ] = (const float complex)(i+3 +1);ta[i+4 ] = (const float complex)(i+4 +1);ta[i+5 ] = (const float complex)(i+5 +1);ta[i+6 ] = (const float complex)(i+6 +1);ta[i+7 ] = (const float complex)(i+7 +1);ta[i+8 ] = (const float complex)(i+8 +1);ta[i+9 ] = (const float complex)(i+9 +1);ta[i+10 ] = (const float complex)(i+10 +1);ta[i+11 ] = (const float complex)(i+11 +1);ta[i+12 ] = (const float complex)(i+12 +1);ta[i+13 ] = (const float complex)(i+13 +1);ta[i+14 ] = (const float complex)(i+14 +1);ta[i+15 ] = (const float complex)(i+15 +1);}
for(     ;i<n;++i){ ta[i+0 ] = (const float complex)(i+0 +1); }
}

	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{ 
	double complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = (const double complex)(i+0 +1);ta[i+1 ] = (const double complex)(i+1 +1);ta[i+2 ] = (const double complex)(i+2 +1);ta[i+3 ] = (const double complex)(i+3 +1);ta[i+4 ] = (const double complex)(i+4 +1);ta[i+5 ] = (const double complex)(i+5 +1);ta[i+6 ] = (const double complex)(i+6 +1);ta[i+7 ] = (const double complex)(i+7 +1);ta[i+8 ] = (const double complex)(i+8 +1);ta[i+9 ] = (const double complex)(i+9 +1);ta[i+10 ] = (const double complex)(i+10 +1);ta[i+11 ] = (const double complex)(i+11 +1);ta[i+12 ] = (const double complex)(i+12 +1);ta[i+13 ] = (const double complex)(i+13 +1);ta[i+14 ] = (const double complex)(i+14 +1);ta[i+15 ] = (const double complex)(i+15 +1);}
for(     ;i<n;++i){ ta[i+0 ] = (const double complex)(i+0 +1); }
}

	}
	else 
#endif
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_do_conjugate(void * array, rsb_type_t type, size_t n)
{
	/*!
	 * \ingroup gr_vec
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
		return RSB_ERR_NO_ERROR;
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
		return RSB_ERR_NO_ERROR;
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		float complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = conjf(ta[i+0 ]);ta[i+1 ] = conjf(ta[i+1 ]);ta[i+2 ] = conjf(ta[i+2 ]);ta[i+3 ] = conjf(ta[i+3 ]);ta[i+4 ] = conjf(ta[i+4 ]);ta[i+5 ] = conjf(ta[i+5 ]);ta[i+6 ] = conjf(ta[i+6 ]);ta[i+7 ] = conjf(ta[i+7 ]);ta[i+8 ] = conjf(ta[i+8 ]);ta[i+9 ] = conjf(ta[i+9 ]);ta[i+10 ] = conjf(ta[i+10 ]);ta[i+11 ] = conjf(ta[i+11 ]);ta[i+12 ] = conjf(ta[i+12 ]);ta[i+13 ] = conjf(ta[i+13 ]);ta[i+14 ] = conjf(ta[i+14 ]);ta[i+15 ] = conjf(ta[i+15 ]);}
for(     ;i<n;++i){ ta[i+0 ] = conjf(ta[i+0 ]); }
}

	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		double complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = conj(ta[i+0 ]);ta[i+1 ] = conj(ta[i+1 ]);ta[i+2 ] = conj(ta[i+2 ]);ta[i+3 ] = conj(ta[i+3 ]);ta[i+4 ] = conj(ta[i+4 ]);ta[i+5 ] = conj(ta[i+5 ]);ta[i+6 ] = conj(ta[i+6 ]);ta[i+7 ] = conj(ta[i+7 ]);ta[i+8 ] = conj(ta[i+8 ]);ta[i+9 ] = conj(ta[i+9 ]);ta[i+10 ] = conj(ta[i+10 ]);ta[i+11 ] = conj(ta[i+11 ]);ta[i+12 ] = conj(ta[i+12 ]);ta[i+13 ] = conj(ta[i+13 ]);ta[i+14 ] = conj(ta[i+14 ]);ta[i+15 ] = conj(ta[i+15 ]);}
for(     ;i<n;++i){ ta[i+0 ] = conj(ta[i+0 ]); }
}

	}
	else 
#endif
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_do_negate(void * array, rsb_type_t type, size_t n)
{
	/*!
	 * \ingroup gr_vec
	 * Will negate the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
{ 
	double*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = -ta[i+0 ];ta[i+1 ] = -ta[i+1 ];ta[i+2 ] = -ta[i+2 ];ta[i+3 ] = -ta[i+3 ];ta[i+4 ] = -ta[i+4 ];ta[i+5 ] = -ta[i+5 ];ta[i+6 ] = -ta[i+6 ];ta[i+7 ] = -ta[i+7 ];ta[i+8 ] = -ta[i+8 ];ta[i+9 ] = -ta[i+9 ];ta[i+10 ] = -ta[i+10 ];ta[i+11 ] = -ta[i+11 ];ta[i+12 ] = -ta[i+12 ];ta[i+13 ] = -ta[i+13 ];ta[i+14 ] = -ta[i+14 ];ta[i+15 ] = -ta[i+15 ];}
for(     ;i<n;++i){ ta[i+0 ] = -ta[i+0 ]; }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
{ 
	float*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = -ta[i+0 ];ta[i+1 ] = -ta[i+1 ];ta[i+2 ] = -ta[i+2 ];ta[i+3 ] = -ta[i+3 ];ta[i+4 ] = -ta[i+4 ];ta[i+5 ] = -ta[i+5 ];ta[i+6 ] = -ta[i+6 ];ta[i+7 ] = -ta[i+7 ];ta[i+8 ] = -ta[i+8 ];ta[i+9 ] = -ta[i+9 ];ta[i+10 ] = -ta[i+10 ];ta[i+11 ] = -ta[i+11 ];ta[i+12 ] = -ta[i+12 ];ta[i+13 ] = -ta[i+13 ];ta[i+14 ] = -ta[i+14 ];ta[i+15 ] = -ta[i+15 ];}
for(     ;i<n;++i){ ta[i+0 ] = -ta[i+0 ]; }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{ 
	float complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = -ta[i+0 ];ta[i+1 ] = -ta[i+1 ];ta[i+2 ] = -ta[i+2 ];ta[i+3 ] = -ta[i+3 ];ta[i+4 ] = -ta[i+4 ];ta[i+5 ] = -ta[i+5 ];ta[i+6 ] = -ta[i+6 ];ta[i+7 ] = -ta[i+7 ];ta[i+8 ] = -ta[i+8 ];ta[i+9 ] = -ta[i+9 ];ta[i+10 ] = -ta[i+10 ];ta[i+11 ] = -ta[i+11 ];ta[i+12 ] = -ta[i+12 ];ta[i+13 ] = -ta[i+13 ];ta[i+14 ] = -ta[i+14 ];ta[i+15 ] = -ta[i+15 ];}
for(     ;i<n;++i){ ta[i+0 ] = -ta[i+0 ]; }
}
}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{ 
	double complex*ta = array;
{
for(i=0;i+15<n;i+=16){
ta[i+0 ] = -ta[i+0 ];ta[i+1 ] = -ta[i+1 ];ta[i+2 ] = -ta[i+2 ];ta[i+3 ] = -ta[i+3 ];ta[i+4 ] = -ta[i+4 ];ta[i+5 ] = -ta[i+5 ];ta[i+6 ] = -ta[i+6 ];ta[i+7 ] = -ta[i+7 ];ta[i+8 ] = -ta[i+8 ];ta[i+9 ] = -ta[i+9 ];ta[i+10 ] = -ta[i+10 ];ta[i+11 ] = -ta[i+11 ];ta[i+12 ] = -ta[i+12 ];ta[i+13 ] = -ta[i+13 ];ta[i+14 ] = -ta[i+14 ];ta[i+15 ] = -ta[i+15 ];}
for(     ;i<n;++i){ ta[i+0 ] = -ta[i+0 ]; }
}
}
	else 
#endif
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_find_min(void * minp, const void * array, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(n<1)return RSB_ERR_BADARGS;
	if(inc<1)return RSB_ERR_BADARGS;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{const double * ap = array;double *mp = minp;
	*mp = *ap;for(i = 1;i<n;++i){if(fabs(ap[i*inc])<fabs(*mp) )*mp = ap[i*inc];
	}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{const float * ap = array;float *mp = minp;
	*mp = *ap;for(i = 1;i<n;++i){if(fabsf(ap[i*inc])<fabsf(*mp) )*mp = ap[i*inc];
	}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{const float complex * ap = array;float complex *mp = minp;
	*mp = *ap;for(i = 1;i<n;++i){if(cabsf(ap[i*inc])<cabsf(*mp) )*mp = ap[i*inc];
	}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{const double complex * ap = array;double complex *mp = minp;
	*mp = *ap;for(i = 1;i<n;++i){if(cabs(ap[i*inc])<cabs(*mp) )*mp = ap[i*inc];
	}}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_find_max(void * maxp, const void * array, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(n<1)return RSB_ERR_BADARGS;
	if(inc<1)return RSB_ERR_BADARGS;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{const double * ap = array;double *mp = maxp;
	*mp = *ap;for(i=1;i<n;++i){if(fabs(ap[i*inc])>fabs(*mp))*mp = ap[i*inc];
	}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{const float * ap = array;float *mp = maxp;
	*mp = *ap;for(i=1;i<n;++i){if(fabsf(ap[i*inc])>fabsf(*mp))*mp = ap[i*inc];
	}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{const float complex * ap = array;float complex *mp = maxp;
	*mp = *ap;for(i=1;i<n;++i){if(cabsf(ap[i*inc])>cabsf(*mp))*mp = ap[i*inc];
	}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{const double complex * ap = array;double complex *mp = maxp;
	*mp = *ap;for(i=1;i<n;++i){if(cabs(ap[i*inc])>cabs(*mp))*mp = ap[i*inc];
	}}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_drop_to_zero_if_above_threshold(void * array, rsb_type_t type, size_t n, const void * threshold)
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{const double th = (*(const double*)(threshold)); double*ta = array;
	for(i = 0;i<n;++i)
	{if(fabs(th)<fabs(ta[i]))ta[i] = ((double)(0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{const float th = (*(const float*)(threshold)); float*ta = array;
	for(i = 0;i<n;++i)
	{if(fabsf(th)<fabsf(ta[i]))ta[i] = ((float)(0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{const float complex th = (*(const float complex*)(threshold)); float complex*ta = array;
	for(i = 0;i<n;++i)
	{if(cabsf(th)<cabsf(ta[i]))ta[i] = ((float complex)(0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{const double complex th = (*(const double complex*)(threshold)); double complex*ta = array;
	for(i = 0;i<n;++i)
	{if(cabs(th)<cabs(ta[i]))ta[i] = ((double complex)(0));}}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_nnz_idx_t rsb__util_count_positive(void * array, rsb_type_t type, size_t n)
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i, c = 0;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{	double*ta = array;
		 for(i=0;i<n;++i)
			c+=((ta[i])>(double)0);
	}else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{	float*ta = array;
		 for(i=0;i<n;++i)
			c+=((ta[i])>(float)0);
	}else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{	float complex*ta = array;
		 for(i=0;i<n;++i)
			c+=(crealf(ta[i])>(float)0);
	}else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{	double complex*ta = array;
		 for(i=0;i<n;++i)
			c+=(creal(ta[i])>(double)0);
	}else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return c;
}

rsb_nnz_idx_t rsb__util_count_negative(void * array, rsb_type_t type, size_t n)
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i, c = 0;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{	double*ta = array;
		 for(i=0;i<n;++i)
			c+=((ta[i])<(double)0);
	}else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{	float*ta = array;
		 for(i=0;i<n;++i)
			c+=((ta[i])<(float)0);
	}else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{	float complex*ta = array;
		 for(i=0;i<n;++i)
			c+=(crealf(ta[i])<(float)0);
	}else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{	double complex*ta = array;
		 for(i=0;i<n;++i)
			c+=(creal(ta[i])<(double)0);
	}else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return c;
}

rsb_err_t rsb__util_drop_to_zero_if_under_threshold(void * array, rsb_type_t type, size_t n, const void * threshold)
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  ) {
	const double th = (*(double*)(threshold)); double*ta = ((double*)(array));
	for(i=0;i<n;++i){if(fabs(th)>fabs(ta[i]))ta[i] = ((double)(0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  ) {
	const float th = (*(float*)(threshold)); float*ta = ((float*)(array));
	for(i=0;i<n;++i){if(fabsf(th)>fabsf(ta[i]))ta[i] = ((float)(0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ) {
	const float complex th = (*(float complex*)(threshold)); float complex*ta = ((float complex*)(array));
	for(i=0;i<n;++i){if(cabsf(th)>cabsf(ta[i]))ta[i] = ((float complex)(0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ) {
	const double complex th = (*(double complex*)(threshold)); double complex*ta = ((double complex*)(array));
	for(i=0;i<n;++i){if(cabs(th)>cabs(ta[i]))ta[i] = ((double complex)(0));}}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__fill_with_ones(void * array, rsb_type_t type, size_t n, size_t incx){
	/*!
	 * \ingroup gr_vec
	 * Will set to one the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * TODO:RENAME: rsb__fill_with_ones -> rsb__val_fill_with_ones.
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  ){
	double*ta = ((double*)(array));
 for(i=0;i<n;++i) {ta[i*incx] = ((double)(1.0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  ){
	float*ta = ((float*)(array));
 for(i=0;i<n;++i) {ta[i*incx] = ((float)(1.0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ){
	float complex*ta = ((float complex*)(array));
 for(i=0;i<n;++i) {ta[i*incx] = ((float complex)(1.0));}}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ){
	double complex*ta = ((double complex*)(array));
 for(i=0;i<n;++i) {ta[i*incx] = ((double complex)(1.0));}}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__debug_print_vectors_diff_fd(const void * v1, const void * v2, size_t n, rsb_type_t type, size_t incx, size_t incy, int onlyfirst, FILE*fd){
	/*! 
	 * A debug function for printing the difference of two vectors of a specified type, in parallel.
	 * FIXME : It should take into account thresholds specific to each numerical type.
	 **/
#if RSB_ALLOW_STDOUT
	size_t i, differing = 0;
	if(!v1 || !v2)return RSB_ERR_BADARGS;

	/*RSB_STDERR("\t vectors diff :\n"); */
	RSB_FPRINTF(fd,"\t vectors diff :\n");
	
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double *v1p = v1,*v2p = v2;
		const double th = 1e-6;
		for(i=0;i<n ;++i) 
												if(fabs((double)(v1p[i*incx]-v2p[i*incy]))>th)/*FIXME : incomplete check*/
{		differing++;
		if((onlyfirst==0)||(onlyfirst>differing))
		RSB_FPRINTF(fd,"%zd : "RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING" "RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING"\n",(rsb_printf_int_t)i,								v1p[i*incx],v2p[i*incy]		);
}
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float *v1p = v1,*v2p = v2;
		const float th = 1e-5;
		for(i=0;i<n ;++i) 
												if(fabs((double)(v1p[i*incx]-v2p[i*incy]))>th)/*FIXME : incomplete check*/
{		differing++;
		if((onlyfirst==0)||(onlyfirst>differing))
		RSB_FPRINTF(fd,"%zd : "RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING" "RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING"\n",(rsb_printf_int_t)i,								v1p[i*incx],v2p[i*incy]		);
}
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex *v1p = v1,*v2p = v2;
		const float th = 1e-4;
		for(i=0;i<n ;++i) 
						if(crealf(v1p[i*incx])-crealf(v2p[i*incy])>th)/*FIXME : incomplete check*/{		differing++;
		if((onlyfirst==0)||(onlyfirst>differing))
		RSB_FPRINTF(fd,"%zd : "RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING" "RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING"\n",(rsb_printf_int_t)i,						crealf(v1p[i*incx]),cimagf(v1p[i*incx]),crealf(v2p[i*incy]),cimagf(v2p[i*incy])		);
}
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex *v1p = v1,*v2p = v2;
		const double th = 1e-6;
		for(i=0;i<n ;++i) 
				if(creal(v1p[i*incx])-creal(v2p[i*incy])>th)/*FIXME : incomplete check*/{		differing++;
		if((onlyfirst==0)||(onlyfirst>differing))
		RSB_FPRINTF(fd,"%zd : "RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING" "RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING"\n",(rsb_printf_int_t)i,				creal(v1p[i*incx]),cimag(v1p[i*incx]),creal(v2p[i*incy]),cimag(v2p[i*incy])		);
}
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	if(differing>onlyfirst)RSB_FPRINTF(fd,"...(for a total of %zd differing entries)...\n",(rsb_printf_int_t)(differing-onlyfirst));
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}

rsb_err_t rsb__debug_print_vectors_diff(const void * v1, const void * v2, size_t n, rsb_type_t type, size_t incx, size_t incy, int onlyfirst){
	return rsb__debug_print_vectors_diff_fd(v1, v2, n, type, incx, incy, onlyfirst, RSB_DEFAULT_STREAM);
}

rsb_err_t rsb__debug_print_value(const void * v, rsb_type_t type){
	/*! 
	 **/
#if RSB_ALLOW_STDOUT
	if(!v)return RSB_ERR_BADARGS;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double *v1p = v;
		RSB_STDOUT(RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING,								v1p[0]		);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float *v1p = v;
		RSB_STDOUT(RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING,								v1p[0]		);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex *v1p = v;
		RSB_STDOUT(RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING,						crealf(v1p[0]),cimagf(v1p[0])		);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex *v1p = v;
		RSB_STDOUT(RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING,				creal(v1p[0]),cimag(v1p[0])		);
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}

rsb_err_t rsb__debug_print_vector_extra(const void * v1, size_t n, rsb_type_t type, size_t inc, int style, FILE*stream){
	/*! 
	 * A debug function for printing two vectors of a specified type, in parallel.
	 **/
#if RSB_ALLOW_STDOUT
	rsb_nnz_idx_t i;
	int want_header = ( style == 0x1 );
	const char * ts = RSB_IS_MATRIX_TYPE_COMPLEX(type)?"complex":"real";
	const char * ss = RSB_SYMMETRY_STRING(RSB_FLAG_NOFLAGS);
	
	if( n <= 0 )
		goto errb;

	if(!v1 || !stream)
		goto errb;

	/*if(!want_header)
		RSB_STDERR("\t vectors  :\n");*/
	
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double *v1p = v1;
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix array %s %s\n%zd %zd\n",ts,ss,(rsb_printf_int_t)n,(rsb_printf_int_t)1);
		for(i=0;i<n;++i) 
		RSB_FPRINTF(stream,RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING "\n",								v1p[i*inc]		);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float *v1p = v1;
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix array %s %s\n%zd %zd\n",ts,ss,(rsb_printf_int_t)n,(rsb_printf_int_t)1);
		for(i=0;i<n;++i) 
		RSB_FPRINTF(stream,RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING "\n",								v1p[i*inc]		);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex *v1p = v1;
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix array %s %s\n%zd %zd\n",ts,ss,(rsb_printf_int_t)n,(rsb_printf_int_t)1);
		for(i=0;i<n;++i) 
		RSB_FPRINTF(stream,RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING "\n",						crealf(v1p[i*inc]),cimagf(v1p[i*inc])		);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex *v1p = v1;
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix array %s %s\n%zd %zd\n",ts,ss,(rsb_printf_int_t)n,(rsb_printf_int_t)1);
		for(i=0;i<n;++i) 
		RSB_FPRINTF(stream,RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING "\n",				creal(v1p[i*inc]),cimag(v1p[i*inc])		);
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
errb:
	return RSB_ERR_BADARGS;
}

rsb_err_t rsb__debug_print_vector(const void * v1, size_t n, rsb_type_t type, size_t inc){
	return rsb__debug_print_vector_extra(v1, n, type, inc, 0x0, stdout);
}

rsb_err_t rsb__debug_print_vectors(const void * v1, const void * v2, size_t n, size_t incx, size_t incy, rsb_type_t type){
	/*! 
	 * A debug function for printing two vectors of a specified type, in parallel.
	 **/
#if RSB_ALLOW_STDOUT
	size_t i;
	if(!v1 || !v2)return RSB_ERR_BADARGS;

	RSB_STDERR("\t vectors  :\n");
	
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double *v1p = v1,*v2p = v2;
		for(i=0;i<n;++i) 
		RSB_STDOUT(RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING" "RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING"\n",v1p[(i)*incx],v2p[(i)*incy]);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float *v1p = v1,*v2p = v2;
		for(i=0;i<n;++i) 
		RSB_STDOUT(RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING" "RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING"\n",v1p[(i)*incx],v2p[(i)*incy]);
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex *v1p = v1,*v2p = v2;
		for(i=0;i<n;++i) 
		RSB_STDOUT(RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING" "RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING"\n",crealf(v1p[(i)*incx]),cimagf(v1p[(i)*incx]),crealf(v2p[(i)*incy]),cimagf(v2p[(i)*incy]));
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex *v1p = v1,*v2p = v2;
		for(i=0;i<n;++i) 
		RSB_STDOUT(RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING" "RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING"\n",creal(v1p[(i)*incx]),cimag(v1p[(i)*incx]),creal(v2p[(i)*incy]),cimag(v2p[(i)*incy]));
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}


#ifdef RSB_WANT_OSKI_BENCHMARKING 
rsb_err_t rsb__do_account_sorted_optimized_css(
	 const rsb_coo_idx_t * MIndx, const rsb_coo_idx_t * mIndx,
	 const rsb_coo_idx_t Mdim, const rsb_coo_idx_t mdim,
	 const rsb_nnz_idx_t nnz, rsb_nnz_idx_t * elements_per_block_row, rsb_nnz_idx_t * blocks_per_block_row
)
{
	/**
	 	\ingroup gr_internals

		elements_per_block_row and blocks_per_block_row arrays should be blanked.
		FIXME : missing error handling.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t n = 0;

	if(blocks_per_block_row)
	for(n=0;n<nnz;++n)
	{
		RSB_DEBUG_ASSERT(MIndx[n]<Mdim);
		RSB_DEBUG_ASSERT(MIndx[n]>=0);
		elements_per_block_row[MIndx[n]]++;
		blocks_per_block_row  [MIndx[n]]++;
	}
	else
	for(n=0;n<nnz;++n)
	{
		RSB_DEBUG_ASSERT(MIndx[n]<Mdim);
		RSB_DEBUG_ASSERT(MIndx[n]>=0);
		elements_per_block_row[MIndx[n]]++;
	}
	RSB_DO_ERR_RETURN(errval)
}

#endif /* RSB_WANT_OSKI_BENCHMARKING */

#if RSB_OBSOLETE_QUARANTINE_UNUSED

rsb_err_t rsb__do_account_sorted_optimized(
	 struct rsb_mtx_t * mtxAp,
	 const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	 const rsb_coo_idx_t Idim, const rsb_coo_idx_t Jdim,
	 const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop,
rsb_nnz_idx_t * elements_per_block_row, 
rsb_nnz_idx_t * blocks_per_block_row
)
{
	/**
	 *	\ingroup gr_internals
	 * 	FIXME : document this
	 */
	rsb_coo_idx_t blockrows = 0;
	rsb_coo_idx_t blockcolumns = 0;
	rsb_coo_idx_t baserow = 0;
	rsb_coo_idx_t basecolumn = 0;
	const rsb_coo_idx_t *Mpntr = NULL;
	const rsb_coo_idx_t *mpntr = NULL;
	const rsb_coo_idx_t *MIndx = NULL;
	const rsb_coo_idx_t *mIndx = NULL;
	rsb_blk_idx_t mI = 0, MI = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;	/* will index a nnz sized array */
	int K = 0;
	

	if(nnz==0)
	{
		/* FIXME: new case, incomplete (useful for diagonal implicit matrices) */
		return RSB_ERR_NO_ERROR;
	}

#if RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS
	if(!pinfop)
	{
		/* a performance fix */
		if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			return rsb__do_account_sorted_optimized_css(JA,IA,Jdim,Idim,nnz,elements_per_block_row,blocks_per_block_row);
		else
			return rsb__do_account_sorted_optimized_css(IA,JA,Idim,Jdim,nnz,elements_per_block_row,blocks_per_block_row);
	}
#endif
	
	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	{
		mpntr = pinfop->rpntr;
		Mpntr = pinfop->cpntr;
		mIndx = IA;
		MIndx = JA;
	}
	else
	{
		Mpntr = pinfop->rpntr;
		mpntr = pinfop->cpntr;
		MIndx = IA;
		mIndx = JA;
	}

	/*	storage BCOR	*/
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCOR )
{
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif


	k = mI = MI = K=0;
	while( MIndx[k] >= Mpntr[MI+1] )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= mpntr[mI+1] )++mI;	/* skipping preceding block columns .. */
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
	baserow = Mpntr[MI];
	basecolumn = mpntr[mI];
	*elements_per_block_row = 0;
	*blocks_per_block_row   = 0;	
	elements_per_block_row[MI*0] += blockrows * blockcolumns;
	blocks_per_block_row[MI]   +=1;

	while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k,(rsb_printf_int_t) (MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */

			while( mIndx[k] >= mpntr[mI+1] )++mI;
			blockcolumns = mpntr[mI+1] - mpntr[mI];
			basecolumn = mpntr[mI];

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];
			}
			else
			{
				/* same block row  */
			}
			*elements_per_block_row += blockrows * blockcolumns;
			blocks_per_block_row[MI]   +=1;
			++K;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */

			while( MIndx[k] >= Mpntr[MI+1] )++MI;
			blockrows    = Mpntr[MI+1] - Mpntr[MI];
			baserow = Mpntr[MI];

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */

				mI = 0;
				while( mIndx[k] >= mpntr[mI+1] )++mI;
				blockcolumns = mpntr[mI+1] - mpntr[mI];
				basecolumn = mpntr[mI];
			}
			else
			{
				/* new row block, same column  */
			}
			/* get rid of this var : elements_per_block_row */
			*elements_per_block_row += blockrows * blockcolumns;
			blocks_per_block_row[MI]   +=1;
			++K;
		}
		else
		{
			/* same block row for sure */
		}
		++k;
	}
	errval = RSB_ERR_NO_ERROR;goto ret;
	}
	/*	storage BCSR	*/
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCSR )
{
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif


	k = mI = MI = K=0;
	while( MIndx[k] >= (blockrows   *(MI+1)) )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= (blockcolumns*(mI+1)) )++mI;	/* skipping preceding block columns .. */
	blockrows    = (blockrows   *(MI+1)) - (blockrows   *(MI));
	blockcolumns = (blockcolumns*(mI+1)) - (blockcolumns*(mI));
	baserow = (blockrows   *(MI));
	basecolumn = (blockcolumns*(mI));
	*elements_per_block_row = 0;
	*blocks_per_block_row   = 0;	
	elements_per_block_row[MI*0] += blockrows * blockcolumns;
	blocks_per_block_row[MI]   +=1;

	while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k,(rsb_printf_int_t) (MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
			mI = mIndx[k]/blockcolumns;
			basecolumn = (blockcolumns*(mI));

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));
			}
			else
			{
				/* same block row  */
			}
			*elements_per_block_row += blockrows * blockcolumns;
			blocks_per_block_row[MI]   +=1;
			++K;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
			MI = MIndx[k]/blockrows;
			baserow = (blockrows   *(MI));

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = mIndx[k]/blockcolumns;
				basecolumn = (blockcolumns*(mI));
			}
			else
			{
				/* new row block, same column  */
			}
			/* get rid of this var : elements_per_block_row */
			*elements_per_block_row += blockrows * blockcolumns;
			blocks_per_block_row[MI]   +=1;
			++K;
		}
		else
		{
			/* same block row for sure */
		}
		++k;
	}
	errval = RSB_ERR_NO_ERROR;goto ret;
	}
	errval = RSB_ERR_INTERNAL_ERROR;
ret:	return errval;
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_err_t rsb__do_insert_sorted_optimized_css( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * MIndx, const rsb_coo_idx_t * mIndx, const rsb_nnz_idx_t nnz)
{
	/**
	 	\ingroup gr_internals

		elements_per_block_row and blocks_per_block_row arrays should be blanked.
		FIXME : missing error handling.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t n = 0;

	/* in case of RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR, they are equal */
	if(mtxAp->VA != VA)
		rsb__memcpy(mtxAp->VA  ,VA  ,mtxAp->el_size*nnz);

	for(n=0;n<nnz+1;++n)
		mtxAp->indptr[n] = n;

	for(n=0;n<mtxAp->nnz;++n)
		mtxAp->bindx [n] = mIndx[n];
	mtxAp->bindx [nnz] = 0;

	// should also set bindx, indptr, 
	RSB_DO_ERR_RETURN(errval)
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_insert_sorted_optimized( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop)
{
	/*
	 *	FIXME ! UNFINISHED 
	 * 	and please note that linked format is incomplete, so it does not support well block column major
	 */
	rsb_coo_idx_t blockrows = 0;
	rsb_coo_idx_t blockcolumns = 0;
	rsb_coo_idx_t baserow = 0;
	rsb_coo_idx_t basecolumn = 0;
	rsb_nnz_idx_t *indptr = mtxAp->indptr;
	rsb_nnz_idx_t *bindx = mtxAp->bindx;
	const rsb_coo_idx_t *Mpntr = NULL;
	const rsb_coo_idx_t *mpntr = NULL;
	const rsb_coo_idx_t *MIndx = NULL;
	const rsb_coo_idx_t *mIndx = NULL;
	rsb_blk_idx_t mI = 0, MI = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;	/* will index a nnz sized array */
	rsb_nnz_idx_t K = 0;

	if(nnz==0)
	{
		/* FIXME: new case, incomplete (useful for diagonal implicit matrices) */
		K = 0;		/* if nnz == 0 then K == 0 */
		bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
		return RSB_ERR_NO_ERROR;
	}


#if RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS
	if(!pinfop)
	{
		/* a performance fix */
		if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			return rsb__do_insert_sorted_optimized_css( mtxAp, VA, JA, IA, nnz );
		else
			return rsb__do_insert_sorted_optimized_css( mtxAp, VA, IA, JA, nnz );
	}
#endif

	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	{
		mpntr = pinfop->rpntr;
		Mpntr = pinfop->cpntr;
		mIndx = IA;
		MIndx = JA;
	}
	else
	{
		Mpntr = pinfop->rpntr;
		mpntr = pinfop->cpntr;
		MIndx = IA;
		mIndx = JA;
	}


	/*	type double, storage BCOR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCOR )
{
	double * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= Mpntr[MI+1] )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= mpntr[mI+1] )++mI;	/* skipping preceding block columns .. */
	baserow = Mpntr[MI];
	basecolumn = mpntr[mI];
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */

			while( mIndx[k] >= mpntr[mI+1] )++mI;
			blockcolumns = mpntr[mI+1] - mpntr[mI];
			basecolumn = mpntr[mI];

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;

				while( mIndx[k] >= mpntr[mI+1] )++mI;
				blockcolumns = mpntr[mI+1] - mpntr[mI];
				basecolumn = mpntr[mI];
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	indptr[ K ]
	/*K * blockrows * blockcolumns*/
	/*RSB_BLOCK_OFFSET(mtxAp,K)/mtxAp->el_size*/ /* FIXME : unfinished ! */ 
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const double*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	/*	type float, storage BCOR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCOR )
{
	float * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= Mpntr[MI+1] )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= mpntr[mI+1] )++mI;	/* skipping preceding block columns .. */
	baserow = Mpntr[MI];
	basecolumn = mpntr[mI];
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */

			while( mIndx[k] >= mpntr[mI+1] )++mI;
			blockcolumns = mpntr[mI+1] - mpntr[mI];
			basecolumn = mpntr[mI];

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;

				while( mIndx[k] >= mpntr[mI+1] )++mI;
				blockcolumns = mpntr[mI+1] - mpntr[mI];
				basecolumn = mpntr[mI];
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	indptr[ K ]
	/*K * blockrows * blockcolumns*/
	/*RSB_BLOCK_OFFSET(mtxAp,K)/mtxAp->el_size*/ /* FIXME : unfinished ! */ 
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const float*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	/*	type float complex, storage BCOR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCOR )
{
	float complex * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= Mpntr[MI+1] )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= mpntr[mI+1] )++mI;	/* skipping preceding block columns .. */
	baserow = Mpntr[MI];
	basecolumn = mpntr[mI];
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */

			while( mIndx[k] >= mpntr[mI+1] )++mI;
			blockcolumns = mpntr[mI+1] - mpntr[mI];
			basecolumn = mpntr[mI];

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;

				while( mIndx[k] >= mpntr[mI+1] )++mI;
				blockcolumns = mpntr[mI+1] - mpntr[mI];
				basecolumn = mpntr[mI];
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	indptr[ K ]
	/*K * blockrows * blockcolumns*/
	/*RSB_BLOCK_OFFSET(mtxAp,K)/mtxAp->el_size*/ /* FIXME : unfinished ! */ 
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const float complex*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	/*	type double complex, storage BCOR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCOR )
{
	double complex * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= Mpntr[MI+1] )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= mpntr[mI+1] )++mI;	/* skipping preceding block columns .. */
	baserow = Mpntr[MI];
	basecolumn = mpntr[mI];
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */

			while( mIndx[k] >= mpntr[mI+1] )++mI;
			blockcolumns = mpntr[mI+1] - mpntr[mI];
			basecolumn = mpntr[mI];

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */

				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows    = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;

				while( mIndx[k] >= mpntr[mI+1] )++mI;
				blockcolumns = mpntr[mI+1] - mpntr[mI];
				basecolumn = mpntr[mI];
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	indptr[ K ]
	/*K * blockrows * blockcolumns*/
	/*RSB_BLOCK_OFFSET(mtxAp,K)/mtxAp->el_size*/ /* FIXME : unfinished ! */ 
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const double complex*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	/*	type double, storage BCSR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCSR )
{
	double * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= (blockrows   *(MI+1)) )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= (blockcolumns*(mI+1)) )++mI;	/* skipping preceding block columns .. */
	baserow = (blockrows   *(MI));
	basecolumn = (blockcolumns*(mI));
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
			mI = mIndx[k]/blockcolumns;
			basecolumn = (blockcolumns*(mI));

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;
				mI = mIndx[k]/blockcolumns;
				basecolumn = (blockcolumns*(mI));
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	K * blockrows * blockcolumns
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const double*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	/*	type float, storage BCSR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCSR )
{
	float * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= (blockrows   *(MI+1)) )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= (blockcolumns*(mI+1)) )++mI;	/* skipping preceding block columns .. */
	baserow = (blockrows   *(MI));
	basecolumn = (blockcolumns*(mI));
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
			mI = mIndx[k]/blockcolumns;
			basecolumn = (blockcolumns*(mI));

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;
				mI = mIndx[k]/blockcolumns;
				basecolumn = (blockcolumns*(mI));
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	K * blockrows * blockcolumns
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const float*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	/*	type float complex, storage BCSR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCSR )
{
	float complex * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= (blockrows   *(MI+1)) )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= (blockcolumns*(mI+1)) )++mI;	/* skipping preceding block columns .. */
	baserow = (blockrows   *(MI));
	basecolumn = (blockcolumns*(mI));
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
			mI = mIndx[k]/blockcolumns;
			basecolumn = (blockcolumns*(mI));

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;
				mI = mIndx[k]/blockcolumns;
				basecolumn = (blockcolumns*(mI));
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	K * blockrows * blockcolumns
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const float complex*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	/*	type double complex, storage BCSR	*/
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	if( mtxAp->matrix_storage==RSB_MATRIX_STORAGE_BCSR )
{
	double complex * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif

	while( MIndx[k] >= (blockrows   *(MI+1)) )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= (blockcolumns*(mI+1)) )++mI;	/* skipping preceding block columns .. */
	baserow = (blockrows   *(MI));
	basecolumn = (blockcolumns*(mI));
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */



	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
			mI = mIndx[k]/blockcolumns;
			basecolumn = (blockcolumns*(mI));

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
				MI = MIndx[k]/blockrows;
				baserow = (blockrows   *(MI));

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;
				mI = mIndx[k]/blockcolumns;
				basecolumn = (blockcolumns*(mI));
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);


		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += 
	K * blockrows * blockcolumns
;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += (MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const double complex*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
	errval = RSB_ERR_INTERNAL_ERROR;
	return errval;
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_block(rsb_type_t type, const void * VA, rsb_blk_idx_t roff, rsb_blk_idx_t coff, rsb_blk_idx_t rows, rsb_blk_idx_t cols )
{
	/*!
	 * Will dump to stdout a dense matrix.
	 * Used for debugging purposes.
	 *
	 * FIXME : should be integrated with the macro subsystem in util.m4, and support column major order, and debugged.
	 */
#if RSB_ALLOW_STDOUT
	register rsb_coo_idx_t i, j;

	if(RSB_BLK_MUL_OVERFLOW(rows,cols))
		return RSB_ERR_LIMITS;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if(type == RSB_NUMERICAL_TYPE_DOUBLE )
	{
		for(i=0;i<rows;++i)for(j=0;j<cols;++j)
		if(((double*)VA)[cols*i+j]!=((double)(0)) )
		{ RSB_STDOUT(""
		"%zd"/* FIXME : this could be any index type! */
		"\t"
		"%zd"
		"\t"
		RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING
		"\n",(rsb_printf_int_t)(roff+i+1),(rsb_printf_int_t)(coff+j+1),
((double*)VA)[cols*i+j]);
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if(type == RSB_NUMERICAL_TYPE_FLOAT )
	{
		for(i=0;i<rows;++i)for(j=0;j<cols;++j)
		if(((float*)VA)[cols*i+j]!=((float)(0)) )
		{ RSB_STDOUT(""
		"%zd"/* FIXME : this could be any index type! */
		"\t"
		"%zd"
		"\t"
		RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING
		"\n",(rsb_printf_int_t)(roff+i+1),(rsb_printf_int_t)(coff+j+1),
((float*)VA)[cols*i+j]);
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if(type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	{
		for(i=0;i<rows;++i)for(j=0;j<cols;++j)
		if(((float complex*)VA)[cols*i+j]!=((float complex)(0)) )
		{ RSB_STDOUT(""
		"%zd"/* FIXME : this could be any index type! */
		"\t"
		"%zd"
		"\t"
		RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING
		"\n",(rsb_printf_int_t)(roff+i+1),(rsb_printf_int_t)(coff+j+1),
crealf(((float complex*)VA)[cols*i+j]),cimagf(((float complex*)VA)[cols*i+j]));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if(type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	{
		for(i=0;i<rows;++i)for(j=0;j<cols;++j)
		if(((double complex*)VA)[cols*i+j]!=((double complex)(0)) )
		{ RSB_STDOUT(""
		"%zd"/* FIXME : this could be any index type! */
		"\t"
		"%zd"
		"\t"
		RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING
		"\n",(rsb_printf_int_t)(roff+i+1),(rsb_printf_int_t)(coff+j+1),
creal(((double complex*)VA)[cols*i+j]),cimag(((double complex*)VA)[cols*i+j]));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_blocks(const struct rsb_mtx_t *mtxAp)
{
	return RSB_ERR_UNIMPLEMENTED_YET;
#if 0
	/*! 
	 * \ingroup gr_internals
	 * A debug function for printing out the matrix structure.
	 *
	 * FIXME : UNFINISHED
	 * Note : it is extremely slow.
	 **/
	rsb_blk_idx_t i,j;
	if(!mtxAp)return RSB_ERR_BADARGS;
	if(!mtxAp->options)return RSB_ERR_BADARGS;

	RSB_STDERR("\t block structure :\n");
	
	/* this prints out the matrix blocks nnz structure */
	for(i=0;i<mtxAp->M_b;++i)
	{
		for(j=0;j<mtxAp->K_b;++j)
		if((RSB_BITMAP_GET(mtxAp->options->bitmap,mtxAp->M_b,mtxAp->K_b,i,j)))
		{
			RSB_STDERR("1");
		}
		else
		{
			RSB_STDERR("0");
		}
		RSB_STDERR("\n");
	}
	return RSB_ERR_NO_ERROR;
#endif
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__test_print_csr(rsb_type_t type, rsb_flags_t flags, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t want_header, FILE*stream)
{
	/**
	 * \ingroup gr_internals
	 * Dumps out a whole matrix, from its CSR representation.
	 * 
	 * Warning : the nonzeros should be sorted on input.
	 */
#if RSB_ALLOW_STDOUT
	rsb_coo_idx_t k;
	if( !stream )goto err;
	if( !IA )goto err;
	if( ( !JA || !VA ) && nnz>0  )goto err;

	RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t)rows);
	/* RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t) nnz); */
	for(k=0;k<rows+1;++k) { RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t)(IA[k]+1)); }
	for(k=0;k<nnz   ;++k) { RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t)(JA[k]+1)); }
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if(type == RSB_NUMERICAL_TYPE_DOUBLE )
	{
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING
				"\n"
				,((double*)VA)[k]);
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if(type == RSB_NUMERICAL_TYPE_FLOAT )
	{
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING
				"\n"
				,((float*)VA)[k]);
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if(type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	{
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING
				"\n"
				,crealf(((float complex*)VA)[k]),cimagf(((float complex*)VA)[k]));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if(type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	{
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING
				"\n"
				,creal(((double complex*)VA)[k]),cimag(((double complex*)VA)[k]));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
err:
	return RSB_ERR_GENERIC_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}

rsb_err_t rsb__test_print_coo_mm(rsb_type_t type, rsb_flags_t flags, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t want_header, FILE*stream)
{
	/**
	 * \ingroup gr_internals
	 * Dumps out a whole matrix, from its coordinates, in matrix market format.
	 * 
	 * Warning : the nonzeros should be sorted on input.
	 */
#if RSB_ALLOW_STDOUT
	rsb_coo_idx_t k;
	const char * ts = RSB_IS_MATRIX_TYPE_COMPLEX(type)?"complex":"real";
	const char * ss = RSB_SYMMETRY_STRING(flags);
	
	if( !stream )
	{
		goto err;
	}

	if( ( !IA || !JA || !VA ) && nnz > 0 )
		goto err;
	if( rows < 0 || cols < 0 || nnz < 0 )
		goto err;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if(type == RSB_NUMERICAL_TYPE_DOUBLE )
	{
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix coordinate %s %s\n%zd %zd %zd\n",ts,ss,(rsb_printf_int_t)rows,(rsb_printf_int_t)cols,(rsb_printf_int_t)nnz);
/*		for(k=0;k<nnz;++k) { RSB_FPRINTF(stream,"%6zd %6zd %20g\n",(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),((float*)VA)[k]); }*/
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				"%zd"
				"\t"
				"%zd"
				"\t"
				RSB_MATRIX_STORAGE_DOUBLE_PRINTF_STRING
				"\n"
				,(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),((double*)VA)[k]);
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if(type == RSB_NUMERICAL_TYPE_FLOAT )
	{
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix coordinate %s %s\n%zd %zd %zd\n",ts,ss,(rsb_printf_int_t)rows,(rsb_printf_int_t)cols,(rsb_printf_int_t)nnz);
/*		for(k=0;k<nnz;++k) { RSB_FPRINTF(stream,"%6zd %6zd %20g\n",(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),((float*)VA)[k]); }*/
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				"%zd"
				"\t"
				"%zd"
				"\t"
				RSB_MATRIX_STORAGE_FLOAT_PRINTF_STRING
				"\n"
				,(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),((float*)VA)[k]);
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if(type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	{
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix coordinate %s %s\n%zd %zd %zd\n",ts,ss,(rsb_printf_int_t)rows,(rsb_printf_int_t)cols,(rsb_printf_int_t)nnz);
/*		for(k=0;k<nnz;++k) { RSB_FPRINTF(stream,"%6zd %6zd %20g\n",(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),((float*)VA)[k]); }*/
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				"%zd"
				"\t"
				"%zd"
				"\t"
				RSB_MATRIX_STORAGE_FLOAT_COMPLEX_PRINTF_STRING
				"\n"
				,(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),crealf(((float complex*)VA)[k]),cimagf(((float complex*)VA)[k]));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if(type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	{
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix coordinate %s %s\n%zd %zd %zd\n",ts,ss,(rsb_printf_int_t)rows,(rsb_printf_int_t)cols,(rsb_printf_int_t)nnz);
/*		for(k=0;k<nnz;++k) { RSB_FPRINTF(stream,"%6zd %6zd %20g\n",(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),((float*)VA)[k]); }*/
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				"%zd"
				"\t"
				"%zd"
				"\t"
				RSB_MATRIX_STORAGE_DOUBLE_COMPLEX_PRINTF_STRING
				"\n"
				,(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),creal(((double complex*)VA)[k]),cimag(((double complex*)VA)[k]));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
err:
	return RSB_ERR_GENERIC_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}

/*static*/ /*inline*/ size_t rsb__do_sizeof(rsb_type_t type)	{
		/*
		 * FIXME : UNUSED ?
		 */
		size_t so = 0;
		switch(type)
		{
			/* supported (double,float,float complex,double complex) */
			case RSB_NUMERICAL_TYPE_DOUBLE 	:
				so = sizeof(double);
			break;
			case RSB_NUMERICAL_TYPE_FLOAT 	:
				so = sizeof(float);
			break;
			case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:
				so = sizeof(float complex);
			break;
			case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:
				so = sizeof(double complex);
			break;
			/* unsupported type */
			default :
			RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS 
		}
		return so;
	}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_coo_sum( struct rsb_coo_mtx_t*coocp, const void *alphap, const struct rsb_coo_mtx_t*cooap, const void *betap,  const struct rsb_coo_mtx_t*coobp)
{
	struct rsb_coo_mtx_t cooa = *cooap, coob = *coobp, cooc = *coocp;
	rsb_nnz_idx_t /*rnz = 0,*/an, bn, cn;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if(cooa.typecode == RSB_NUMERICAL_TYPE_DOUBLE )
	{
	double alpha = alphap?*(double*)alphap:((double)(1.0));
	double beta  = betap ?*(double*)betap :((double)(1.0));
	for(cn = 0, an = 0, bn = 0;an<cooa.nnz || bn<coob.nnz;)
	{
		rsb_nnz_idx_t ap = an, bp = bn;
		if(cooa.IA[an]==coob.IA[bn] && cooa.JA[an]==coob.JA[bn])
			cooc.IA[cn] = cooa.IA[an],cooc.JA[cn] = cooa.JA[an],
			((double*)cooc.VA)[cn] = alpha * ((double*)cooa.VA)[an] + beta * ((double*)coob.VA)[bn],
			ap = an, bp = bn, ++cn, ++an, ++bn;

		for(;an<cooa.nnz && cooa.IA[an]==cooa.IA[ap] && cooa.JA[an]==cooa.JA[ap] ;++an)
			//RSB_STDOUT("x> %d %d\n",cooa.IA[an],cooa.JA[an])
			((double*)cooc.VA)[cn] += alpha * ((double*)cooa.VA)[an];

		for(;bn<coob.nnz && coob.IA[bn]==coob.IA[bp] && coob.JA[bn]==coob.JA[bp] ;++bn)
			//RSB_STDOUT("x> %d %d\n",coob.IA[bn],coob.JA[bn])
			((double*)cooc.VA)[cn] += beta  * ((double*)coob.VA)[bn];

		if( bn<coob.nnz )
		for(;an<cooa.nnz && (cooa.IA[an]<coob.IA[bn] ||
			       	(cooa.IA[an] <= coob.IA[bn] && cooa.JA[an]<coob.JA[bn]))
			       	;++an)
				//RSB_STDOUT("-> %d %d\n",cooa.IA[an],cooa.JA[an]),
			cooc.IA[cn] = cooa.IA[an], cooc.JA[cn] = cooa.JA[an],
			((double*)cooc.VA)[cn] = alpha * ((double*)cooa.VA)[an],
			++cn;

		if( an<cooa.nnz )
		for(;bn<coob.nnz && (cooa.IA[an]>coob.IA[bn] ||
			       	(cooa.IA[an]>=coob.IA[bn] && cooa.JA[an]>coob.JA[bn]))
			       	;++bn)
			//	RSB_STDOUT("-> %d %d\n",coob.IA[bn],coob.JA[bn]),
			cooc.IA[cn] = coob.IA[bn],cooc.JA[cn] = coob.JA[bn],
			((double*)cooc.VA)[cn] = beta * ((double*)coob.VA)[bn],
			++cn;
		//RSB_STDOUT("? %d %d\n",an,bn);
	}
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if(cooa.typecode == RSB_NUMERICAL_TYPE_FLOAT )
	{
	float alpha = alphap?*(float*)alphap:((float)(1.0));
	float beta  = betap ?*(float*)betap :((float)(1.0));
	for(cn = 0, an = 0, bn = 0;an<cooa.nnz || bn<coob.nnz;)
	{
		rsb_nnz_idx_t ap = an, bp = bn;
		if(cooa.IA[an]==coob.IA[bn] && cooa.JA[an]==coob.JA[bn])
			cooc.IA[cn] = cooa.IA[an],cooc.JA[cn] = cooa.JA[an],
			((float*)cooc.VA)[cn] = alpha * ((float*)cooa.VA)[an] + beta * ((float*)coob.VA)[bn],
			ap = an, bp = bn, ++cn, ++an, ++bn;

		for(;an<cooa.nnz && cooa.IA[an]==cooa.IA[ap] && cooa.JA[an]==cooa.JA[ap] ;++an)
			//RSB_STDOUT("x> %d %d\n",cooa.IA[an],cooa.JA[an])
			((float*)cooc.VA)[cn] += alpha * ((float*)cooa.VA)[an];

		for(;bn<coob.nnz && coob.IA[bn]==coob.IA[bp] && coob.JA[bn]==coob.JA[bp] ;++bn)
			//RSB_STDOUT("x> %d %d\n",coob.IA[bn],coob.JA[bn])
			((float*)cooc.VA)[cn] += beta  * ((float*)coob.VA)[bn];

		if( bn<coob.nnz )
		for(;an<cooa.nnz && (cooa.IA[an]<coob.IA[bn] ||
			       	(cooa.IA[an] <= coob.IA[bn] && cooa.JA[an]<coob.JA[bn]))
			       	;++an)
				//RSB_STDOUT("-> %d %d\n",cooa.IA[an],cooa.JA[an]),
			cooc.IA[cn] = cooa.IA[an], cooc.JA[cn] = cooa.JA[an],
			((float*)cooc.VA)[cn] = alpha * ((float*)cooa.VA)[an],
			++cn;

		if( an<cooa.nnz )
		for(;bn<coob.nnz && (cooa.IA[an]>coob.IA[bn] ||
			       	(cooa.IA[an]>=coob.IA[bn] && cooa.JA[an]>coob.JA[bn]))
			       	;++bn)
			//	RSB_STDOUT("-> %d %d\n",coob.IA[bn],coob.JA[bn]),
			cooc.IA[cn] = coob.IA[bn],cooc.JA[cn] = coob.JA[bn],
			((float*)cooc.VA)[cn] = beta * ((float*)coob.VA)[bn],
			++cn;
		//RSB_STDOUT("? %d %d\n",an,bn);
	}
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if(cooa.typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	{
	float complex alpha = alphap?*(float complex*)alphap:((float complex)(1.0));
	float complex beta  = betap ?*(float complex*)betap :((float complex)(1.0));
	for(cn = 0, an = 0, bn = 0;an<cooa.nnz || bn<coob.nnz;)
	{
		rsb_nnz_idx_t ap = an, bp = bn;
		if(cooa.IA[an]==coob.IA[bn] && cooa.JA[an]==coob.JA[bn])
			cooc.IA[cn] = cooa.IA[an],cooc.JA[cn] = cooa.JA[an],
			((float complex*)cooc.VA)[cn] = alpha * ((float complex*)cooa.VA)[an] + beta * ((float complex*)coob.VA)[bn],
			ap = an, bp = bn, ++cn, ++an, ++bn;

		for(;an<cooa.nnz && cooa.IA[an]==cooa.IA[ap] && cooa.JA[an]==cooa.JA[ap] ;++an)
			//RSB_STDOUT("x> %d %d\n",cooa.IA[an],cooa.JA[an])
			((float complex*)cooc.VA)[cn] += alpha * ((float complex*)cooa.VA)[an];

		for(;bn<coob.nnz && coob.IA[bn]==coob.IA[bp] && coob.JA[bn]==coob.JA[bp] ;++bn)
			//RSB_STDOUT("x> %d %d\n",coob.IA[bn],coob.JA[bn])
			((float complex*)cooc.VA)[cn] += beta  * ((float complex*)coob.VA)[bn];

		if( bn<coob.nnz )
		for(;an<cooa.nnz && (cooa.IA[an]<coob.IA[bn] ||
			       	(cooa.IA[an] <= coob.IA[bn] && cooa.JA[an]<coob.JA[bn]))
			       	;++an)
				//RSB_STDOUT("-> %d %d\n",cooa.IA[an],cooa.JA[an]),
			cooc.IA[cn] = cooa.IA[an], cooc.JA[cn] = cooa.JA[an],
			((float complex*)cooc.VA)[cn] = alpha * ((float complex*)cooa.VA)[an],
			++cn;

		if( an<cooa.nnz )
		for(;bn<coob.nnz && (cooa.IA[an]>coob.IA[bn] ||
			       	(cooa.IA[an]>=coob.IA[bn] && cooa.JA[an]>coob.JA[bn]))
			       	;++bn)
			//	RSB_STDOUT("-> %d %d\n",coob.IA[bn],coob.JA[bn]),
			cooc.IA[cn] = coob.IA[bn],cooc.JA[cn] = coob.JA[bn],
			((float complex*)cooc.VA)[cn] = beta * ((float complex*)coob.VA)[bn],
			++cn;
		//RSB_STDOUT("? %d %d\n",an,bn);
	}
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if(cooa.typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	{
	double complex alpha = alphap?*(double complex*)alphap:((double complex)(1.0));
	double complex beta  = betap ?*(double complex*)betap :((double complex)(1.0));
	for(cn = 0, an = 0, bn = 0;an<cooa.nnz || bn<coob.nnz;)
	{
		rsb_nnz_idx_t ap = an, bp = bn;
		if(cooa.IA[an]==coob.IA[bn] && cooa.JA[an]==coob.JA[bn])
			cooc.IA[cn] = cooa.IA[an],cooc.JA[cn] = cooa.JA[an],
			((double complex*)cooc.VA)[cn] = alpha * ((double complex*)cooa.VA)[an] + beta * ((double complex*)coob.VA)[bn],
			ap = an, bp = bn, ++cn, ++an, ++bn;

		for(;an<cooa.nnz && cooa.IA[an]==cooa.IA[ap] && cooa.JA[an]==cooa.JA[ap] ;++an)
			//RSB_STDOUT("x> %d %d\n",cooa.IA[an],cooa.JA[an])
			((double complex*)cooc.VA)[cn] += alpha * ((double complex*)cooa.VA)[an];

		for(;bn<coob.nnz && coob.IA[bn]==coob.IA[bp] && coob.JA[bn]==coob.JA[bp] ;++bn)
			//RSB_STDOUT("x> %d %d\n",coob.IA[bn],coob.JA[bn])
			((double complex*)cooc.VA)[cn] += beta  * ((double complex*)coob.VA)[bn];

		if( bn<coob.nnz )
		for(;an<cooa.nnz && (cooa.IA[an]<coob.IA[bn] ||
			       	(cooa.IA[an] <= coob.IA[bn] && cooa.JA[an]<coob.JA[bn]))
			       	;++an)
				//RSB_STDOUT("-> %d %d\n",cooa.IA[an],cooa.JA[an]),
			cooc.IA[cn] = cooa.IA[an], cooc.JA[cn] = cooa.JA[an],
			((double complex*)cooc.VA)[cn] = alpha * ((double complex*)cooa.VA)[an],
			++cn;

		if( an<cooa.nnz )
		for(;bn<coob.nnz && (cooa.IA[an]>coob.IA[bn] ||
			       	(cooa.IA[an]>=coob.IA[bn] && cooa.JA[an]>coob.JA[bn]))
			       	;++bn)
			//	RSB_STDOUT("-> %d %d\n",coob.IA[bn],coob.JA[bn]),
			cooc.IA[cn] = coob.IA[bn],cooc.JA[cn] = coob.JA[bn],
			((double complex*)cooc.VA)[cn] = beta * ((double complex*)coob.VA)[bn],
			++cn;
		//RSB_STDOUT("? %d %d\n",an,bn);
	}
	}
	else 
#endif
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__cor_merge_dups(rsb_type_t typecode, void* RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t offB, rsb_nnz_idx_t nnzB, rsb_nnz_idx_t nnzC, const int wv, int wp, rsb_nnz_idx_t *onzp, struct rsb_coo_mtx_t*RSB_RESTRICT coop)
{
	/**
		See rsb__cor_merge.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VB = NULL, *VC = NULL, *VT = NULL;
	rsb_coo_idx_t * IB = NULL, *JB = NULL;
	rsb_coo_idx_t * IC = NULL, *JC = NULL;
	rsb_coo_idx_t * IT = NULL, *JT = NULL;
	rsb_nnz_idx_t bi = 0, ci = 0, ti = 0;
	rsb_nnz_idx_t b0 = 0, c0 = 0;
	rsb_nnz_idx_t onz = 0;
	struct rsb_coo_mtx_t coo;
	size_t es = RSB_SIZEOF(typecode);

	if( nnzB == 0 || nnzC == 0 )
	{
		goto ret;
	}

	b0 = offB;
	c0 = offB + nnzB;
	VB = RSB_TYPED_OFF_PTR(typecode,VA,b0);
	VC = RSB_TYPED_OFF_PTR(typecode,VA,c0);
	IB = IA + b0;
	IC = IA + c0;
	JB = JA + b0;
	JC = JA + c0;

	RSB_BZERO_P(&coo);
	coo.nnz = nnzB + nnzC;
	coo.typecode = typecode;

	if( coop && coop->nnz)
	{
		coo = *coop;
		coo.nnz = nnzB + nnzC; /* necessary */
	}
	else
	{
		if( NULL == rsb__allocate_coo_matrix_t(&coo) )
			goto err;
	}

	IT = coo.IA;
	JT = coo.JA;
	VT = coo.VA;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if(typecode == RSB_NUMERICAL_TYPE_DOUBLE )
	{
	double * vT = VT;
	double * vB = VB;
	double * vC = VC;

again_double:

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

	/* FIXME: this works as RSB_FLAG_DUPLICATES_SUM but should support either merge, last, first, ...  */
       	if   ( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi] + vC[ci];
		++bi,++ci,++ti;
		++onz;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi] + vC[ci];
		++bi;
		++ci;
		++ti;
		++onz;
	}

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC && bi < nnzB )
		goto again_double;

again_once_double:

       	if   ( bi<nnzB && ci==nnzC )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci==nnzC && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

       	if   ( ci<nnzC && bi==nnzB )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( ci<nnzC && bi==nnzB && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti]+= vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC || bi < nnzB )
		goto again_once_double;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if(typecode == RSB_NUMERICAL_TYPE_FLOAT )
	{
	float * vT = VT;
	float * vB = VB;
	float * vC = VC;

again_float:

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

	/* FIXME: this works as RSB_FLAG_DUPLICATES_SUM but should support either merge, last, first, ...  */
       	if   ( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi] + vC[ci];
		++bi,++ci,++ti;
		++onz;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi] + vC[ci];
		++bi;
		++ci;
		++ti;
		++onz;
	}

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC && bi < nnzB )
		goto again_float;

again_once_float:

       	if   ( bi<nnzB && ci==nnzC )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci==nnzC && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

       	if   ( ci<nnzC && bi==nnzB )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( ci<nnzC && bi==nnzB && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti]+= vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC || bi < nnzB )
		goto again_once_float;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if(typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	{
	float complex * vT = VT;
	float complex * vB = VB;
	float complex * vC = VC;

again_float_complex:

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

	/* FIXME: this works as RSB_FLAG_DUPLICATES_SUM but should support either merge, last, first, ...  */
       	if   ( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi] + vC[ci];
		++bi,++ci,++ti;
		++onz;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi] + vC[ci];
		++bi;
		++ci;
		++ti;
		++onz;
	}

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC && bi < nnzB )
		goto again_float_complex;

again_once_float_complex:

       	if   ( bi<nnzB && ci==nnzC )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci==nnzC && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

       	if   ( ci<nnzC && bi==nnzB )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( ci<nnzC && bi==nnzB && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti]+= vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC || bi < nnzB )
		goto again_once_float_complex;
	}
	else 
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if(typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	{
	double complex * vT = VT;
	double complex * vB = VB;
	double complex * vC = VC;

again_double_complex:

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

	/* FIXME: this works as RSB_FLAG_DUPLICATES_SUM but should support either merge, last, first, ...  */
       	if   ( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi] + vC[ci];
		++bi,++ci,++ti;
		++onz;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi] + vC[ci];
		++bi;
		++ci;
		++ti;
		++onz;
	}

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC && bi < nnzB )
		goto again_double_complex;

again_once_double_complex:

       	if   ( bi<nnzB && ci==nnzC )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci==nnzC && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

       	if   ( ci<nnzC && bi==nnzB )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( ci<nnzC && bi==nnzB && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti]+= vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC || bi < nnzB )
		goto again_once_double_complex;
	}
	else 
#endif
		errval = RSB_ERR_INTERNAL_ERROR;

	coo.nnz -= onz;
	RSB_COA_MEMCPY(IA,IT,offB,0,(coo.nnz));
	RSB_COA_MEMCPY(JA,JT,offB,0,(coo.nnz));
	if(wp)
	{
		RSB_A_MEMCPY_parallel(  VA,VT,offB,0,(coo.nnz),es);
	}
	else
	{
		RSB_A_MEMCPY(  VA,VT,offB,0,(coo.nnz),es);
	}
	RSB_ASSERT(rsb__util_is_coo_array_sorted_up_partial_order(IA,coo.nnz));
	goto done;
err:
	errval = RSB_ERR_ENOMEM;
done:
	if( coop && coop->nnz)
		;
	else
		rsb__destroy_coo_matrix_t(&coo);
	RSB_ASSIGN_IF(onzp,onz);
ret:
	return errval;
}

rsb_err_t rsb__do_copy_converted_scaled(const void *RSB_RESTRICT  src, void *RSB_RESTRICT dst, const void *RSB_RESTRICT  alphap, rsb_type_t stype,rsb_type_t dtype, size_t nnz, rsb_trans_t transA)
{
	/*!
	 * Copies scaled and conj-transposed.
	 * alpha according to src code type.
	 * \return \rsberrcodemsg
	 * */
	rsb_nnz_idx_t nzi;

	if((!dst) || (!src))
		return RSB_ERR_BADARGS;

	if( stype == RSB_NUMERICAL_TYPE_DOUBLE  && dtype == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double alpha = alphap?*(double*)alphap:((double)(1.0));
		const double*tsrc = src;
		double*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double)(alpha*tsrc[nzi]);
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_DOUBLE  && dtype == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const double alpha = alphap?*(double*)alphap:((double)(1.0));
		const double*tsrc = src;
		float*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float)(alpha*tsrc[nzi]);
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_DOUBLE  && dtype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const double alpha = alphap?*(double*)alphap:((double)(1.0));
		const double*tsrc = src;
		float complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*conjf(tsrc[nzi])) + 0*I;
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*tsrc[nzi]) + 0*I;
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_DOUBLE  && dtype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double alpha = alphap?*(double*)alphap:((double)(1.0));
		const double*tsrc = src;
		double complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*conj(tsrc[nzi])) + 0*I;
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*tsrc[nzi]) + 0*I;
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT  && dtype == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const float alpha = alphap?*(float*)alphap:((float)(1.0));
		const float*tsrc = src;
		double*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double)(alpha*tsrc[nzi]);
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT  && dtype == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float alpha = alphap?*(float*)alphap:((float)(1.0));
		const float*tsrc = src;
		float*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float)(alpha*tsrc[nzi]);
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT  && dtype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float alpha = alphap?*(float*)alphap:((float)(1.0));
		const float*tsrc = src;
		float complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*conjf(tsrc[nzi])) + 0*I;
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*tsrc[nzi]) + 0*I;
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT  && dtype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const float alpha = alphap?*(float*)alphap:((float)(1.0));
		const float*tsrc = src;
		double complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*conj(tsrc[nzi])) + 0*I;
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*tsrc[nzi]) + 0*I;
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const float complex alpha = alphap?*(float complex*)alphap:((float complex)(1.0));
		const float complex*tsrc = src;
		double*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = crealf((double)(alpha*tsrc[nzi]));
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float complex alpha = alphap?*(float complex*)alphap:((float complex)(1.0));
		const float complex*tsrc = src;
		float*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = crealf((float)(alpha*tsrc[nzi]));
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex alpha = alphap?*(float complex*)alphap:((float complex)(1.0));
		const float complex*tsrc = src;
		float complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*conjf(tsrc[nzi]));
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*tsrc[nzi]);
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const float complex alpha = alphap?*(float complex*)alphap:((float complex)(1.0));
		const float complex*tsrc = src;
		double complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*conj(tsrc[nzi]));
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*tsrc[nzi]);
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double complex alpha = alphap?*(double complex*)alphap:((double complex)(1.0));
		const double complex*tsrc = src;
		double*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = creal((double)(alpha*tsrc[nzi]));
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const double complex alpha = alphap?*(double complex*)alphap:((double complex)(1.0));
		const double complex*tsrc = src;
		float*tdst = dst;
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = creal((float)(alpha*tsrc[nzi]));
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const double complex alpha = alphap?*(double complex*)alphap:((double complex)(1.0));
		const double complex*tsrc = src;
		float complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*conjf(tsrc[nzi]));
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (float complex)(alpha*tsrc[nzi]);
	}
	else 
	if( stype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  && dtype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex alpha = alphap?*(double complex*)alphap:((double complex)(1.0));
		const double complex*tsrc = src;
		double complex*tdst = dst;
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*conj(tsrc[nzi]));
		else
			for(nzi=0;nzi<nnz;++nzi) tdst[nzi] = (double complex)(alpha*tsrc[nzi]);
	}
	else 
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_csc2csr(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t*flagsp)
{
	/*!
	 * */
	rsb_nnz_idx_t nzi = 0, nzo;
	rsb_coo_idx_t nr, nc;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_bool_t islowtri = RSB_BOOL_TRUE, isupptri = RSB_BOOL_TRUE;
	rsb_nnz_idx_t lowtrin = 0, upptrin = 0;

	RSB_BZERO(oIA, sizeof(*oIA)*(m+1));
	oIA[0] = offo;
	for(nzi=0;nzi<nnz;++nzi)
		oIA[IA[nzi]-offi+1]++;
	for(nr=0;nr<m;++nr)
		oIA[nr+1]+=oIA[nr];
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	for(nc=0;nc<k;++nc)
	for(nzi = JA[nc]-offi;nzi<JA[nc+1]-offi;++nzi)
	{
		nzo = oIA[IA[nzi]-offi]++;
		oJA[nzo] = nc+offo;
		((double*)oVA)[nzo] = ((const double*)VA)[nzi];
	}
	else 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	for(nc=0;nc<k;++nc)
	for(nzi = JA[nc]-offi;nzi<JA[nc+1]-offi;++nzi)
	{
		nzo = oIA[IA[nzi]-offi]++;
		oJA[nzo] = nc+offo;
		((float*)oVA)[nzo] = ((const float*)VA)[nzi];
	}
	else 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	for(nc=0;nc<k;++nc)
	for(nzi = JA[nc]-offi;nzi<JA[nc+1]-offi;++nzi)
	{
		nzo = oIA[IA[nzi]-offi]++;
		oJA[nzo] = nc+offo;
		((float complex*)oVA)[nzo] = ((const float complex*)VA)[nzi];
	}
	else 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	for(nc=0;nc<k;++nc)
	for(nzi = JA[nc]-offi;nzi<JA[nc+1]-offi;++nzi)
	{
		nzo = oIA[IA[nzi]-offi]++;
		oJA[nzo] = nc+offo;
		((double complex*)oVA)[nzo] = ((const double complex*)VA)[nzi];
	}
	else 
		return RSB_ERR_UNSUPPORTED_TYPE;
	for(nc=0;nc<k;++nc)
	for(nzi=JA[nc]-offi;nzi<JA[nc+1]-offi;++nzi)
	{
		oIA[IA[nzi]-offi]--;
		if(IA[nzi]-offi>nc)
			lowtrin++;
		else
			if(IA[nzi]-offi<nc)
				upptrin++;
	}
	if(upptrin)
		islowtri = RSB_BOOL_FALSE;
	if(lowtrin)
		isupptri = RSB_BOOL_FALSE;
	if(isupptri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
	if(islowtri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
	if( RSB_XOR(upptrin,lowtrin) )
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
	if(*flagsp) RSB_DO_FLAG_ADD(*flagsp,flags);
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_coo_copy_and_stats(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, rsb_coo_idx_t*m, rsb_coo_idx_t*k, const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t iflags, rsb_flags_t*flagsp)
{
	/*!
         * FIXME: unfinished! shall support also typecode-based zeros removal
	 * */
	rsb_nnz_idx_t nzi = 0;
	rsb_coo_idx_t maxi = 0,maxj = 0;
	rsb_bool_t islowtri = RSB_BOOL_TRUE,isupptri = RSB_BOOL_TRUE;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_nnz_idx_t lowtrin = 0,upptrin = 0;

	if(nnz<1)
		goto done;

	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
{
	rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
	maxi = i, maxj = j;
	((double*)oVA)[nzi] = ((double*)VA)[nzi];
	oIA[nzi] = i-offi+offo;
	oJA[nzi] = j-offi+offo;
	lowtrin |= (i>j), upptrin |= (i<j);
	for(nzi=1;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi],j = JA[nzi];
		maxi = RSB_MAX(maxi, i);
		maxj = RSB_MAX(maxj, j);
		((double*)oVA)[nzi] = ((double*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
		lowtrin |= (i>j);
		upptrin |= (i<j);
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
{
	rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
	maxi = i, maxj = j;
	((float*)oVA)[nzi] = ((float*)VA)[nzi];
	oIA[nzi] = i-offi+offo;
	oJA[nzi] = j-offi+offo;
	lowtrin |= (i>j), upptrin |= (i<j);
	for(nzi=1;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi],j = JA[nzi];
		maxi = RSB_MAX(maxi, i);
		maxj = RSB_MAX(maxj, j);
		((float*)oVA)[nzi] = ((float*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
		lowtrin |= (i>j);
		upptrin |= (i<j);
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{
	rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
	maxi = i, maxj = j;
	((float complex*)oVA)[nzi] = ((float complex*)VA)[nzi];
	oIA[nzi] = i-offi+offo;
	oJA[nzi] = j-offi+offo;
	lowtrin |= (i>j), upptrin |= (i<j);
	for(nzi=1;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi],j = JA[nzi];
		maxi = RSB_MAX(maxi, i);
		maxj = RSB_MAX(maxj, j);
		((float complex*)oVA)[nzi] = ((float complex*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
		lowtrin |= (i>j);
		upptrin |= (i<j);
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{
	rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
	maxi = i, maxj = j;
	((double complex*)oVA)[nzi] = ((double complex*)VA)[nzi];
	oIA[nzi] = i-offi+offo;
	oJA[nzi] = j-offi+offo;
	lowtrin |= (i>j), upptrin |= (i<j);
	for(nzi=1;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi],j = JA[nzi];
		maxi = RSB_MAX(maxi, i);
		maxj = RSB_MAX(maxj, j);
		((double complex*)oVA)[nzi] = ((double complex*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
		lowtrin |= (i>j);
		upptrin |= (i<j);
	}
}
	else
		return RSB_ERR_UNSUPPORTED_TYPE;
	if(upptrin)
		islowtri = RSB_BOOL_FALSE;
	if(lowtrin)
		isupptri = RSB_BOOL_FALSE;
	if(isupptri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
	if(islowtri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
	if( RSB_XOR(upptrin,lowtrin) )
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
	if(flagsp) RSB_DO_FLAG_ADD(*flagsp,flags);
	if(flagsp)
		if(RSB_DO_FLAG_HAS_INTERSECTION(*flagsp,RSB_FLAG_ANY_SYMMETRY))
			maxi = maxj = RSB_MAX(maxi,maxj);
	if(m) *m = maxi+1;
	if(k) *k = maxj+1;
done:
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__util_coo_copy(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo)
{
	/*!
         * FIXME: unfinished! shall support also typecode-based zeros removal
	 * */
	rsb_nnz_idx_t nzi = 0;

	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
{
	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
		((double*)oVA)[nzi] = ((double*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
{
	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
		((float*)oVA)[nzi] = ((float*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{
	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
		((float complex*)oVA)[nzi] = ((float complex*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{
	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
		((double complex*)oVA)[nzi] = ((double complex*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
	}
}
	else
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

/* sparse blas level 1 equivalent functions */

int rsb__BLAS_Xusdot(const rsb_type_t typecode, const enum blas_conj_type conj_arg, const rsb_blas_int_t nz, const void*x, const rsb_blas_int_t*indx, const void*y, const rsb_blas_int_t incy, void*r, const enum blas_base_type index_base)
{
	/*!
		\rsb_spblasl1_dot_msg
		\rsb_warn_untested_msg
	*/
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
{
	double*xa = (double*)x;
	double*ya = (double*)y;
	double*rp = (double*)r;
	double ac = ((double)(0));
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += xa[nzi] * ya[xi*incy];
	}
	RSB_SET_IF_NOT_NULL(rp,ac);
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
{
	float*xa = (float*)x;
	float*ya = (float*)y;
	float*rp = (float*)r;
	float ac = ((float)(0));
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += xa[nzi] * ya[xi*incy];
	}
	RSB_SET_IF_NOT_NULL(rp,ac);
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{
	float complex*xa = (float complex*)x;
	float complex*ya = (float complex*)y;
	float complex*rp = (float complex*)r;
	float complex ac = ((float complex)(0));
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	if( conj_arg == blas_conj )
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += conjf(xa[nzi]) * ya[xi*incy];
	}
	else
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += xa[nzi] * ya[xi*incy];
	}
	RSB_SET_IF_NOT_NULL(rp,ac);
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{
	double complex*xa = (double complex*)x;
	double complex*ya = (double complex*)y;
	double complex*rp = (double complex*)r;
	double complex ac = ((double complex)(0));
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	if( conj_arg == blas_conj )
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += conj(xa[nzi]) * ya[xi*incy];
	}
	else
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += xa[nzi] * ya[xi*incy];
	}
	RSB_SET_IF_NOT_NULL(rp,ac);
}
	else
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

int rsb__BLAS_Xusaxpy(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*alpha, const void*x, const rsb_blas_int_t*indx, const void*y, const rsb_blas_int_t incy, const enum blas_base_type index_base)
{
	/*!
		\rsb_spblasl1_axpy_msg
		\rsb_warn_untested_msg
	*/
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
{
	const double*xa = (const double*)x;
	double*ya = (double*)y;
	const double alphav = *(double*)alpha;
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[nzi*incy] += alphav*xa[xi];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
{
	const float*xa = (const float*)x;
	float*ya = (float*)y;
	const float alphav = *(float*)alpha;
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[nzi*incy] += alphav*xa[xi];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{
	const float complex*xa = (const float complex*)x;
	float complex*ya = (float complex*)y;
	const float complex alphav = *(float complex*)alpha;
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[nzi*incy] += alphav*xa[xi];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{
	const double complex*xa = (const double complex*)x;
	double complex*ya = (double complex*)y;
	const double complex alphav = *(double complex*)alpha;
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[nzi*incy] += alphav*xa[xi];
	}
}
	else
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

int rsb__BLAS_Xusga(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*y, const rsb_blas_int_t incy, void*x, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
{
	/*!
		\rsb_spblasl1_ga_msg
		\rsb_warn_untested_msg
	*/
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
{
	double*xa = (double*)x;
	const double*ya = (const double*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
{
	float*xa = (float*)x;
	const float*ya = (const float*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{
	float complex*xa = (float complex*)x;
	const float complex*ya = (const float complex*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{
	double complex*xa = (double complex*)x;
	const double complex*ya = (const double complex*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
	}
}
	else
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

int rsb__BLAS_Xusgz(const rsb_type_t typecode, const rsb_blas_int_t nz, void*y, const rsb_blas_int_t incy, void*x, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
{
	/*!
		\rsb_spblasl1_gz_msg
		\rsb_warn_untested_msg
	*/
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
{
	double*xa = (double*)x;
	double*ya = (double*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
		ya[xi*incy] = ((double)(0));
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
{
	float*xa = (float*)x;
	float*ya = (float*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
		ya[xi*incy] = ((float)(0));
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{
	float complex*xa = (float complex*)x;
	float complex*ya = (float complex*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
		ya[xi*incy] = ((float complex)(0));
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{
	double complex*xa = (double complex*)x;
	double complex*ya = (double complex*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
		ya[xi*incy] = ((double complex)(0));
	}
}
	else
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

int rsb__BLAS_Xussc(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*x, void*y, const rsb_blas_int_t incy, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
{
	/*!
		\rsb_spblasl1_sc_msg
		\rsb_warn_untested_msg
	*/
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
{
	const double*xa = (const double*)x;
	double*ya = (double*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[xi*incy] = xa[nzi];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
{
	const float*xa = (const float*)x;
	float*ya = (float*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[xi*incy] = xa[nzi];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
{
	const float complex*xa = (const float complex*)x;
	float complex*ya = (float complex*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[xi*incy] = xa[nzi];
	}
}
	else
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
{
	const double complex*xa = (const double complex*)x;
	double complex*ya = (double complex*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[xi*incy] = xa[nzi];
	}
}
	else
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

/* blas level 1 equivalent functions */

rsb_err_t rsb__cblas_Xcopy(rsb_type_t typecode, rsb_nnz_idx_t n, const void * x, rsb_nnz_idx_t incx, void * y, rsb_nnz_idx_t incy)
{
	return rsb__xcopy_strided_typed(y,x,0,0,n,typecode,incy,incx);
}

rsb_err_t rsb__cblas_Xnrm2(rsb_type_t type, size_t n, const void * a, rsb_nnz_idx_t incA, void * c){
	/*!
	 * c <- sqrt(sum(|a_i|^2))
         *
	 * \param a	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see dznrm2 in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
	const double*ta = a;double *tc = c,acc = ((double)(0)),tmp = ((double)(0));
	{
for(i=0;i+15<n;i+=16){
	acc = fabs(ta[(i+0 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+1 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+2 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+3 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+4 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+5 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+6 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+7 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+8 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+9 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+10 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+11 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+12 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+13 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+14 )*incA]);tmp += acc*acc;
		acc = fabs(ta[(i+15 )*incA]);tmp += acc*acc;
	}
for(     ;i<n;++i){ 	acc = fabs(ta[(i+0 )*incA]);tmp += acc*acc;
	 }
}
; 
	tc[0] = (sqrt(tmp));
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( type == RSB_NUMERICAL_TYPE_FLOAT  )
	{
	const float*ta = a;float *tc = c,acc = ((float)(0)),tmp = ((float)(0));
	{
for(i=0;i+15<n;i+=16){
	acc = fabsf(ta[(i+0 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+1 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+2 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+3 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+4 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+5 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+6 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+7 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+8 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+9 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+10 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+11 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+12 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+13 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+14 )*incA]);tmp += acc*acc;
		acc = fabsf(ta[(i+15 )*incA]);tmp += acc*acc;
	}
for(     ;i<n;++i){ 	acc = fabsf(ta[(i+0 )*incA]);tmp += acc*acc;
	 }
}
; 
	tc[0] = (sqrtf(tmp));
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
	const float complex*ta = a;float *tc = c,acc = ((float complex)(0)),tmp = ((float complex)(0));
	{
for(i=0;i+15<n;i+=16){
	acc = cabsf(ta[(i+0 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+1 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+2 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+3 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+4 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+5 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+6 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+7 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+8 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+9 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+10 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+11 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+12 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+13 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+14 )*incA]);tmp += acc*acc;
		acc = cabsf(ta[(i+15 )*incA]);tmp += acc*acc;
	}
for(     ;i<n;++i){ 	acc = cabsf(ta[(i+0 )*incA]);tmp += acc*acc;
	 }
}
; 
	tc[0] = crealf(csqrtf(tmp));
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
	const double complex*ta = a;double *tc = c,acc = ((double complex)(0)),tmp = ((double complex)(0));
	{
for(i=0;i+15<n;i+=16){
	acc = cabs(ta[(i+0 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+1 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+2 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+3 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+4 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+5 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+6 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+7 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+8 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+9 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+10 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+11 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+12 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+13 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+14 )*incA]);tmp += acc*acc;
		acc = cabs(ta[(i+15 )*incA]);tmp += acc*acc;
	}
for(     ;i<n;++i){ 	acc = cabs(ta[(i+0 )*incA]);tmp += acc*acc;
	 }
}
; 
	tc[0] = creal(csqrt(tmp));
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__cblas_Xdotu_sub(rsb_type_t type, size_t n, const void * x, rsb_nnz_idx_t incx, const void * y, rsb_nnz_idx_t incy, void *dotu){
	/*!
	 * */
	return rsb__vector_mult_sum(x,y,dotu,type,n,incx,incy);
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__cblas_Xscal(rsb_type_t type, size_t n, const void * alphap, void * a, size_t stride){
	/*!
	 * a <- a * alpha
	 * */
	return rsb_strided_vector_scale(a,alphap,type,n,stride);
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__coo_insertion_sort(rsb_type_t typecode, void* VB, rsb_coo_idx_t * IB, rsb_coo_idx_t * JB, rsb_nnz_idx_t offA, rsb_nnz_idx_t nnzA)
{
	/* only for *small* arrays, where allocation of a temporary array is not justified */
	rsb_coo_idx_t * IA = NULL, *JA = NULL;
	rsb_nnz_idx_t i, j;

	IA = IB + offA;
	JA = JB + offA;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		double * VA = (double*) RSB_TYPED_OFF_PTR(typecode,VB,offA);
		for(i=1;i<nnzA;++i)
		for(j=i;j>0 && RSB_COO_LT(IA[j],JA[j],IA[j-1],JA[j-1]);--j)
		{
			RSB_SWAP(rsb_coo_idx_t,IA[j],IA[j-1]);
			RSB_SWAP(rsb_coo_idx_t,JA[j],JA[j-1]);
			RSB_SWAP(double        ,VA[j],VA[j-1]);
		}
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		float * VA = (float*) RSB_TYPED_OFF_PTR(typecode,VB,offA);
		for(i=1;i<nnzA;++i)
		for(j=i;j>0 && RSB_COO_LT(IA[j],JA[j],IA[j-1],JA[j-1]);--j)
		{
			RSB_SWAP(rsb_coo_idx_t,IA[j],IA[j-1]);
			RSB_SWAP(rsb_coo_idx_t,JA[j],JA[j-1]);
			RSB_SWAP(float        ,VA[j],VA[j-1]);
		}
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		float complex * VA = (float complex*) RSB_TYPED_OFF_PTR(typecode,VB,offA);
		for(i=1;i<nnzA;++i)
		for(j=i;j>0 && RSB_COO_LT(IA[j],JA[j],IA[j-1],JA[j-1]);--j)
		{
			RSB_SWAP(rsb_coo_idx_t,IA[j],IA[j-1]);
			RSB_SWAP(rsb_coo_idx_t,JA[j],JA[j-1]);
			RSB_SWAP(float complex        ,VA[j],VA[j-1]);
		}
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		double complex * VA = (double complex*) RSB_TYPED_OFF_PTR(typecode,VB,offA);
		for(i=1;i<nnzA;++i)
		for(j=i;j>0 && RSB_COO_LT(IA[j],JA[j],IA[j-1],JA[j-1]);--j)
		{
			RSB_SWAP(rsb_coo_idx_t,IA[j],IA[j-1]);
			RSB_SWAP(rsb_coo_idx_t,JA[j],JA[j-1]);
			RSB_SWAP(double complex        ,VA[j],VA[j-1]);
		}
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void rsb__coo_to_lr( void * RSB_RESTRICT VBu, rsb_coo_idx_t*RSB_RESTRICT IB, rsb_coo_idx_t*RSB_RESTRICT JB, void * RSB_RESTRICT VAu, rsb_coo_idx_t*RSB_RESTRICT IA, rsb_coo_idx_t*RSB_RESTRICT JA, rsb_coo_idx_t mj, rsb_nnz_idx_t nnzA, rsb_nnz_idx_t nzoffB, rsb_nnz_idx_t nzoffA, rsb_nnz_idx_t*RSB_RESTRICT nzlp, rsb_nnz_idx_t*RSB_RESTRICT nzrp, rsb_coo_idx_t iadd, rsb_coo_idx_t jadd, rsb_type_t typecode)
{
	/*
	 * Given COO arrays matrices A and (temporary) B, store the coefficients left of the mj-th column before the one coming after it, respecting their original ordering.
	 * A serial function.
	 * Each of arrays VAu, IA, JA fit nzoffA+nnzA elements.
	 * Each of arrays VBu, IB, JB fit nzoffB+nnzA elements.
	 * *nzlp will be set to the count of the elements left  of mj.
	 * *nzrp will be set to the count of the elements right of mj.
	 * iadd will be added to the each of the nnzA IA element.
	 * jadd will be added to the each of the nnzA JA element which is >= mj.
	 * */
	rsb_nnz_idx_t nzl = 0, nzr = 0, nzi = 0;

	RSB_DEBUG_ASSERT(IA!=IB);
	RSB_DEBUG_ASSERT(JA!=JB);

	IA += nzoffA;
	JA += nzoffA;

	IB += nzoffB;
	JB += nzoffB;
	
switch(typecode)
{
			/* supported (double,float,float complex,double complex) */
case RSB_NUMERICAL_TYPE_DOUBLE 	:
{
	double * VA = ((double *)VAu) + nzoffA; 
	double * VB = ((double *)VBu) + nzoffB; 

	RSB_DEBUG_ASSERT(VA!=VB);

	/* visit IA, JA, VA */
	for(nzi=0;nzi<nnzA;++nzi)
	{
		if( JA[nzi] < mj )
		{
			/* elements left of mj go to the beginning of IB,JB,VB */
			IB[nzl] = IA[nzi] + iadd;
			JB[nzl] = JA[nzi];
			VB[nzl] = VA[nzi];
			nzl++;
		}
		else
		{
			/* elements right of mj go to the end of IB,JB,VB, reversed */
			nzr++;
			IB[nnzA-nzr] = IA[nzi] + iadd;
			JB[nnzA-nzr] = JA[nzi] + jadd;
			VB[nnzA-nzr] = VA[nzi];
		}
	}

	/* first nzl elements lay left of mj; the last nzr are right of mj, reversed */
	RSB_DEBUG_ASSERT( nzl+nzr == nnzA );

	/* copy left elements back to A */
	for(nzi=0;nzi<nzl ;++nzi)
	{
		IA[nzi] = IB[nzi];
		JA[nzi] = JB[nzi];
		VA[nzi] = VB[nzi];
	}
	
	/* copy right elements back to A, reversed to original relative order */
	for(     ;nzi<nnzA;++nzi)
	{
		IA[nzi] = IB[nnzA-(1+nzi-nzl)];
		JA[nzi] = JB[nnzA-(1+nzi-nzl)];
		VA[nzi] = VB[nnzA-(1+nzi-nzl)];
	}
}
	break;
case RSB_NUMERICAL_TYPE_FLOAT 	:
{
	float * VA = ((float *)VAu) + nzoffA; 
	float * VB = ((float *)VBu) + nzoffB; 

	RSB_DEBUG_ASSERT(VA!=VB);

	/* visit IA, JA, VA */
	for(nzi=0;nzi<nnzA;++nzi)
	{
		if( JA[nzi] < mj )
		{
			/* elements left of mj go to the beginning of IB,JB,VB */
			IB[nzl] = IA[nzi] + iadd;
			JB[nzl] = JA[nzi];
			VB[nzl] = VA[nzi];
			nzl++;
		}
		else
		{
			/* elements right of mj go to the end of IB,JB,VB, reversed */
			nzr++;
			IB[nnzA-nzr] = IA[nzi] + iadd;
			JB[nnzA-nzr] = JA[nzi] + jadd;
			VB[nnzA-nzr] = VA[nzi];
		}
	}

	/* first nzl elements lay left of mj; the last nzr are right of mj, reversed */
	RSB_DEBUG_ASSERT( nzl+nzr == nnzA );

	/* copy left elements back to A */
	for(nzi=0;nzi<nzl ;++nzi)
	{
		IA[nzi] = IB[nzi];
		JA[nzi] = JB[nzi];
		VA[nzi] = VB[nzi];
	}
	
	/* copy right elements back to A, reversed to original relative order */
	for(     ;nzi<nnzA;++nzi)
	{
		IA[nzi] = IB[nnzA-(1+nzi-nzl)];
		JA[nzi] = JB[nnzA-(1+nzi-nzl)];
		VA[nzi] = VB[nnzA-(1+nzi-nzl)];
	}
}
	break;
case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 	:
{
	float complex * VA = ((float complex *)VAu) + nzoffA; 
	float complex * VB = ((float complex *)VBu) + nzoffB; 

	RSB_DEBUG_ASSERT(VA!=VB);

	/* visit IA, JA, VA */
	for(nzi=0;nzi<nnzA;++nzi)
	{
		if( JA[nzi] < mj )
		{
			/* elements left of mj go to the beginning of IB,JB,VB */
			IB[nzl] = IA[nzi] + iadd;
			JB[nzl] = JA[nzi];
			VB[nzl] = VA[nzi];
			nzl++;
		}
		else
		{
			/* elements right of mj go to the end of IB,JB,VB, reversed */
			nzr++;
			IB[nnzA-nzr] = IA[nzi] + iadd;
			JB[nnzA-nzr] = JA[nzi] + jadd;
			VB[nnzA-nzr] = VA[nzi];
		}
	}

	/* first nzl elements lay left of mj; the last nzr are right of mj, reversed */
	RSB_DEBUG_ASSERT( nzl+nzr == nnzA );

	/* copy left elements back to A */
	for(nzi=0;nzi<nzl ;++nzi)
	{
		IA[nzi] = IB[nzi];
		JA[nzi] = JB[nzi];
		VA[nzi] = VB[nzi];
	}
	
	/* copy right elements back to A, reversed to original relative order */
	for(     ;nzi<nnzA;++nzi)
	{
		IA[nzi] = IB[nnzA-(1+nzi-nzl)];
		JA[nzi] = JB[nnzA-(1+nzi-nzl)];
		VA[nzi] = VB[nnzA-(1+nzi-nzl)];
	}
}
	break;
case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 	:
{
	double complex * VA = ((double complex *)VAu) + nzoffA; 
	double complex * VB = ((double complex *)VBu) + nzoffB; 

	RSB_DEBUG_ASSERT(VA!=VB);

	/* visit IA, JA, VA */
	for(nzi=0;nzi<nnzA;++nzi)
	{
		if( JA[nzi] < mj )
		{
			/* elements left of mj go to the beginning of IB,JB,VB */
			IB[nzl] = IA[nzi] + iadd;
			JB[nzl] = JA[nzi];
			VB[nzl] = VA[nzi];
			nzl++;
		}
		else
		{
			/* elements right of mj go to the end of IB,JB,VB, reversed */
			nzr++;
			IB[nnzA-nzr] = IA[nzi] + iadd;
			JB[nnzA-nzr] = JA[nzi] + jadd;
			VB[nnzA-nzr] = VA[nzi];
		}
	}

	/* first nzl elements lay left of mj; the last nzr are right of mj, reversed */
	RSB_DEBUG_ASSERT( nzl+nzr == nnzA );

	/* copy left elements back to A */
	for(nzi=0;nzi<nzl ;++nzi)
	{
		IA[nzi] = IB[nzi];
		JA[nzi] = JB[nzi];
		VA[nzi] = VB[nzi];
	}
	
	/* copy right elements back to A, reversed to original relative order */
	for(     ;nzi<nnzA;++nzi)
	{
		IA[nzi] = IB[nnzA-(1+nzi-nzl)];
		JA[nzi] = JB[nnzA-(1+nzi-nzl)];
		VA[nzi] = VB[nnzA-(1+nzi-nzl)];
	}
}
	break;
	/* unsupported type */
	default :
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS 
}

	*nzlp = nzl;
	*nzrp = nzr;
}
rsb_err_t rsb__util_testall(void)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_blk_idx_t inc;
	RSB_STDOUT("MTX PRINT TEST BEGIN\n");
{
	const double iVA [] = {1,2};
	const rsb_coo_idx_t iIA[] = {0,1};
	const rsb_coo_idx_t iJA[] = {0,1};
	double oVA [] = {-1,-2};
	rsb_coo_idx_t oJA[] = {-1,-2};
	rsb_coo_idx_t oIA[] = {-1,-2};
	rsb_coo_idx_t m=0,k=0;
	const rsb_nnz_idx_t nnz=2;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const rsb_coo_idx_t offi=0, offo=0;
	double exp = ((double)(0));
	const double alpha = 2;
	double res = 1;
 	rsb_flags_t iflags=RSB_FLAG_NOFLAGS, *flagsp=&iflags;

	inc = 1;

	// {-1,-2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, &m, &k, nnz, typecode, offi, offo, iflags, flagsp);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, NULL, NULL, nnz, typecode, offi, offo, iflags, NULL);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+2,+4}
	errval = rsb__cblas_Xscal(typecode, nnz, &alpha, oVA, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=6) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+2,+4} -> {+4,+8}
	errval = rsb_strided_vector_scale(oVA,&alpha,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=12) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+8} -> {+2,+1}
	oVA[0] /= 2; // avoid pow roundoff errors
	oVA[1] /= 8; // avoid pow roundoff errors
	errval = rsb__util_vector_pow(oVA, typecode, &alpha, nnz);
	// {+2,+1} -> {+4,+1}
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=5) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+1} -> {+1,+1}
	errval = rsb__util_vector_pow(oVA, typecode, &exp, nnz);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+1} -> {+2,+1}
	errval = rsb_strided_vector_scale(oVA, &alpha, typecode, nnz*0+1, inc*2); // test for inc==2
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__test_print_coo_mm(typecode, iflags, oIA, oJA, oVA, m, k, nnz, RSB_BOOL_TRUE, RSB_DEFAULT_STREAM);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const float iVA [] = {1,2};
	const rsb_coo_idx_t iIA[] = {0,1};
	const rsb_coo_idx_t iJA[] = {0,1};
	float oVA [] = {-1,-2};
	rsb_coo_idx_t oJA[] = {-1,-2};
	rsb_coo_idx_t oIA[] = {-1,-2};
	rsb_coo_idx_t m=0,k=0;
	const rsb_nnz_idx_t nnz=2;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const rsb_coo_idx_t offi=0, offo=0;
	float exp = ((float)(0));
	const float alpha = 2;
	float res = 1;
 	rsb_flags_t iflags=RSB_FLAG_NOFLAGS, *flagsp=&iflags;

	inc = 1;

	// {-1,-2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, &m, &k, nnz, typecode, offi, offo, iflags, flagsp);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, NULL, NULL, nnz, typecode, offi, offo, iflags, NULL);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+2,+4}
	errval = rsb__cblas_Xscal(typecode, nnz, &alpha, oVA, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=6) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+2,+4} -> {+4,+8}
	errval = rsb_strided_vector_scale(oVA,&alpha,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=12) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+8} -> {+2,+1}
	oVA[0] /= 2; // avoid pow roundoff errors
	oVA[1] /= 8; // avoid pow roundoff errors
	errval = rsb__util_vector_pow(oVA, typecode, &alpha, nnz);
	// {+2,+1} -> {+4,+1}
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=5) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+1} -> {+1,+1}
	errval = rsb__util_vector_pow(oVA, typecode, &exp, nnz);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+1} -> {+2,+1}
	errval = rsb_strided_vector_scale(oVA, &alpha, typecode, nnz*0+1, inc*2); // test for inc==2
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__test_print_coo_mm(typecode, iflags, oIA, oJA, oVA, m, k, nnz, RSB_BOOL_TRUE, RSB_DEFAULT_STREAM);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const float complex iVA [] = {1,2};
	const rsb_coo_idx_t iIA[] = {0,1};
	const rsb_coo_idx_t iJA[] = {0,1};
	float complex oVA [] = {-1,-2};
	rsb_coo_idx_t oJA[] = {-1,-2};
	rsb_coo_idx_t oIA[] = {-1,-2};
	rsb_coo_idx_t m=0,k=0;
	const rsb_nnz_idx_t nnz=2;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const rsb_coo_idx_t offi=0, offo=0;
	float complex exp = ((float complex)(0));
	const float complex alpha = 2;
	float complex res = 1;
 	rsb_flags_t iflags=RSB_FLAG_NOFLAGS, *flagsp=&iflags;

	inc = 1;

	// {-1,-2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, &m, &k, nnz, typecode, offi, offo, iflags, flagsp);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, NULL, NULL, nnz, typecode, offi, offo, iflags, NULL);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+2,+4}
	errval = rsb__cblas_Xscal(typecode, nnz, &alpha, oVA, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=6) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+2,+4} -> {+4,+8}
	errval = rsb_strided_vector_scale(oVA,&alpha,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=12) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+8} -> {+2,+1}
	oVA[0] /= 2; // avoid pow roundoff errors
	oVA[1] /= 8; // avoid pow roundoff errors
	errval = rsb__util_vector_pow(oVA, typecode, &alpha, nnz);
	// {+2,+1} -> {+4,+1}
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=5) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+1} -> {+1,+1}
	errval = rsb__util_vector_pow(oVA, typecode, &exp, nnz);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+1} -> {+2,+1}
	errval = rsb_strided_vector_scale(oVA, &alpha, typecode, nnz*0+1, inc*2); // test for inc==2
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__test_print_coo_mm(typecode, iflags, oIA, oJA, oVA, m, k, nnz, RSB_BOOL_TRUE, RSB_DEFAULT_STREAM);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const double complex iVA [] = {1,2};
	const rsb_coo_idx_t iIA[] = {0,1};
	const rsb_coo_idx_t iJA[] = {0,1};
	double complex oVA [] = {-1,-2};
	rsb_coo_idx_t oJA[] = {-1,-2};
	rsb_coo_idx_t oIA[] = {-1,-2};
	rsb_coo_idx_t m=0,k=0;
	const rsb_nnz_idx_t nnz=2;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const rsb_coo_idx_t offi=0, offo=0;
	double complex exp = ((double complex)(0));
	const double complex alpha = 2;
	double complex res = 1;
 	rsb_flags_t iflags=RSB_FLAG_NOFLAGS, *flagsp=&iflags;

	inc = 1;

	// {-1,-2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, &m, &k, nnz, typecode, offi, offo, iflags, flagsp);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, NULL, NULL, nnz, typecode, offi, offo, iflags, NULL);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+2,+4}
	errval = rsb__cblas_Xscal(typecode, nnz, &alpha, oVA, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=6) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+2,+4} -> {+4,+8}
	errval = rsb_strided_vector_scale(oVA,&alpha,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=12) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+8} -> {+2,+1}
	oVA[0] /= 2; // avoid pow roundoff errors
	oVA[1] /= 8; // avoid pow roundoff errors
	errval = rsb__util_vector_pow(oVA, typecode, &alpha, nnz);
	// {+2,+1} -> {+4,+1}
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=5) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+1} -> {+1,+1}
	errval = rsb__util_vector_pow(oVA, typecode, &exp, nnz);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+1} -> {+2,+1}
	errval = rsb_strided_vector_scale(oVA, &alpha, typecode, nnz*0+1, inc*2); // test for inc==2
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__test_print_coo_mm(typecode, iflags, oIA, oJA, oVA, m, k, nnz, RSB_BOOL_TRUE, RSB_DEFAULT_STREAM);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

	RSB_STDOUT("MTX PRINT TEST END\n");

for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double VA [] = {1,0,2,0,3,0};
	const double VAL[17] = {1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1};
	double nrm2[2] = {0,-1};
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);

	// TODO: maybe get rid of redundance here.
	errval = rsb__util_vector_sum_strided(&nrm2[0], VA, type, n, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, n, VA, inc, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(nrm2[0]!=nrm2[0]) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, 17, VAL, 1, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(roundf(nrm2[1]) != 4 ) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE ;
	double VA [] = {+1,-1,+3,-3,+4,-4};
	const double VAU[] = {+0,+0,+3,-3,+4,-4};
	const double VAA[] = {+0,+0,+3,-3,+0,+0};
	double soad = ((double)(0));
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);
	#if defined(__INTEL_COMPILER) /* meant to circumvent e.g. icc -O2 -fp-model fast=2 (observed on 19.1.2.254 or 2021.3.0) */
	        const double lthr = 3 * 0.95, hthr = 3 * 1.05;
	#else
	        const double lthr = 3, hthr = 3;
	#endif

	errval = rsb__util_drop_to_zero_if_under_threshold(VA,type,n,&lthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAU,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_drop_to_zero_if_above_threshold(VA,type,n,&hthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAA,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT ;
	const float VA [] = {1,0,2,0,3,0};
	const float VAL[17] = {1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1};
	float nrm2[2] = {0,-1};
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);

	// TODO: maybe get rid of redundance here.
	errval = rsb__util_vector_sum_strided(&nrm2[0], VA, type, n, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, n, VA, inc, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(nrm2[0]!=nrm2[0]) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, 17, VAL, 1, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(roundf(nrm2[1]) != 4 ) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT ;
	float VA [] = {+1,-1,+3,-3,+4,-4};
	const float VAU[] = {+0,+0,+3,-3,+4,-4};
	const float VAA[] = {+0,+0,+3,-3,+0,+0};
	float soad = ((float)(0));
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);
	#if defined(__INTEL_COMPILER) /* meant to circumvent e.g. icc -O2 -fp-model fast=2 (observed on 19.1.2.254 or 2021.3.0) */
	        const float lthr = 3 * 0.95, hthr = 3 * 1.05;
	#else
	        const float lthr = 3, hthr = 3;
	#endif

	errval = rsb__util_drop_to_zero_if_under_threshold(VA,type,n,&lthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAU,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_drop_to_zero_if_above_threshold(VA,type,n,&hthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAA,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex VA [] = {1,0,2,0,3,0};
	const float complex VAL[17] = {1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1};
	float complex nrm2[2] = {0,-1};
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);

	// TODO: maybe get rid of redundance here.
	errval = rsb__util_vector_sum_strided(&nrm2[0], VA, type, n, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, n, VA, inc, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(nrm2[0]!=nrm2[0]) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, 17, VAL, 1, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(roundf(nrm2[1]) != 4 ) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	float complex VA [] = {+1,-1,+3,-3,+4,-4};
	const float complex VAU[] = {+0,+0,+3,-3,+4,-4};
	const float complex VAA[] = {+0,+0,+3,-3,+0,+0};
	float complex soad = ((float complex)(0));
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);
	#if defined(__INTEL_COMPILER) /* meant to circumvent e.g. icc -O2 -fp-model fast=2 (observed on 19.1.2.254 or 2021.3.0) */
	        const float complex lthr = 3 * 0.95, hthr = 3 * 1.05;
	#else
	        const float complex lthr = 3, hthr = 3;
	#endif

	errval = rsb__util_drop_to_zero_if_under_threshold(VA,type,n,&lthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAU,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_drop_to_zero_if_above_threshold(VA,type,n,&hthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAA,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex VA [] = {1,0,2,0,3,0};
	const double complex VAL[17] = {1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1};
	double complex nrm2[2] = {0,-1};
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);

	// TODO: maybe get rid of redundance here.
	errval = rsb__util_vector_sum_strided(&nrm2[0], VA, type, n, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, n, VA, inc, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(nrm2[0]!=nrm2[0]) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, 17, VAL, 1, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(roundf(nrm2[1]) != 4 ) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	double complex VA [] = {+1,-1,+3,-3,+4,-4};
	const double complex VAU[] = {+0,+0,+3,-3,+4,-4};
	const double complex VAA[] = {+0,+0,+3,-3,+0,+0};
	double complex soad = ((double complex)(0));
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);
	#if defined(__INTEL_COMPILER) /* meant to circumvent e.g. icc -O2 -fp-model fast=2 (observed on 19.1.2.254 or 2021.3.0) */
	        const double complex lthr = 3 * 0.95, hthr = 3 * 1.05;
	#else
	        const double complex lthr = 3, hthr = 3;
	#endif

	errval = rsb__util_drop_to_zero_if_under_threshold(VA,type,n,&lthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAU,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_drop_to_zero_if_above_threshold(VA,type,n,&hthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAA,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

RSB_STDOUT("DIFF PRINT TEST BEGIN\n");
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double VAU[] = {+0,+0,+3,-3,+4,-4};
	const double VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0]));

	// TODO: beef this test up; add e.g. quiet option to rsb__debug_print_vectors_diff
	errval = rsb__debug_print_vectors_diff(VAU,VAA,n,type,1,1,0);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT ;
	const float VAU[] = {+0,+0,+3,-3,+4,-4};
	const float VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0]));

	// TODO: beef this test up; add e.g. quiet option to rsb__debug_print_vectors_diff
	errval = rsb__debug_print_vectors_diff(VAU,VAA,n,type,1,1,0);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex VAU[] = {+0,+0,+3,-3,+4,-4};
	const float complex VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0]));

	// TODO: beef this test up; add e.g. quiet option to rsb__debug_print_vectors_diff
	errval = rsb__debug_print_vectors_diff(VAU,VAA,n,type,1,1,0);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex VAU[] = {+0,+0,+3,-3,+4,-4};
	const double complex VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0]));

	// TODO: beef this test up; add e.g. quiet option to rsb__debug_print_vectors_diff
	errval = rsb__debug_print_vectors_diff(VAU,VAA,n,type,1,1,0);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double VAU[] = {+0,+0,+3,-3,+4,-4};
	const double VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0])*inc);

	errval = rsb__debug_print_vectors(VAU,VAA,n,1,inc,type);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__debug_print_vector(VAU,n,type,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT ;
	const float VAU[] = {+0,+0,+3,-3,+4,-4};
	const float VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0])*inc);

	errval = rsb__debug_print_vectors(VAU,VAA,n,1,inc,type);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__debug_print_vector(VAU,n,type,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex VAU[] = {+0,+0,+3,-3,+4,-4};
	const float complex VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0])*inc);

	errval = rsb__debug_print_vectors(VAU,VAA,n,1,inc,type);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__debug_print_vector(VAU,n,type,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex VAU[] = {+0,+0,+3,-3,+4,-4};
	const double complex VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0])*inc);

	errval = rsb__debug_print_vectors(VAU,VAA,n,1,inc,type);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__debug_print_vector(VAU,n,type,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
RSB_STDOUT("DIFF PRINT TEST END\n");

{
	const int inc = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE ;
	double v1[] = {+1,+0,+5,-2,+4,-6};
	double v2[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t off = 1;
	const rsb_nnz_idx_t n = sizeof(v1)/(sizeof(v1[0])*inc)-off;
	const double zero = ((double)(0)), sum = 2;
	double res = zero;
	rsb_nnz_idx_t cnt;

	cnt = rsb__util_count_positive(v1,type,n+off);
	if(cnt!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	cnt = rsb__util_count_negative(v1,type,n+off);
	if(cnt!=2) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	errval = rsb__vectors_left_sum_reduce_and_zero(v1,v2,type,n,inc,off);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v1+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=sum) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v2+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=zero) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const int inc = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT ;
	float v1[] = {+1,+0,+5,-2,+4,-6};
	float v2[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t off = 1;
	const rsb_nnz_idx_t n = sizeof(v1)/(sizeof(v1[0])*inc)-off;
	const float zero = ((float)(0)), sum = 2;
	float res = zero;
	rsb_nnz_idx_t cnt;

	cnt = rsb__util_count_positive(v1,type,n+off);
	if(cnt!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	cnt = rsb__util_count_negative(v1,type,n+off);
	if(cnt!=2) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	errval = rsb__vectors_left_sum_reduce_and_zero(v1,v2,type,n,inc,off);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v1+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=sum) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v2+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=zero) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const int inc = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	float complex v1[] = {+1,+0,+5,-2,+4,-6};
	float complex v2[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t off = 1;
	const rsb_nnz_idx_t n = sizeof(v1)/(sizeof(v1[0])*inc)-off;
	const float complex zero = ((float complex)(0)), sum = 2;
	float complex res = zero;
	rsb_nnz_idx_t cnt;

	cnt = rsb__util_count_positive(v1,type,n+off);
	if(cnt!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	cnt = rsb__util_count_negative(v1,type,n+off);
	if(cnt!=2) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	errval = rsb__vectors_left_sum_reduce_and_zero(v1,v2,type,n,inc,off);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v1+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=sum) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v2+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=zero) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const int inc = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	double complex v1[] = {+1,+0,+5,-2,+4,-6};
	double complex v2[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t off = 1;
	const rsb_nnz_idx_t n = sizeof(v1)/(sizeof(v1[0])*inc)-off;
	const double complex zero = ((double complex)(0)), sum = 2;
	double complex res = zero;
	rsb_nnz_idx_t cnt;

	cnt = rsb__util_count_positive(v1,type,n+off);
	if(cnt!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	cnt = rsb__util_count_negative(v1,type,n+off);
	if(cnt!=2) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	errval = rsb__vectors_left_sum_reduce_and_zero(v1,v2,type,n,inc,off);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v1+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=sum) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v2+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=zero) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

#if RSB_WITH_SPARSE_BLAS_INTERFACE
{
	/* BLAS types and interfaces only here  */
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT ;
	const float x[] = {+1,+0,+5,-2,+4,-6};
	float y[] = {+2,+0,+3,-2,+0,+0};
	float w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const float zero = ((float)(0));
	float alpha = 2;
	float dotr = zero;
	float res = zero;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	int istat;
#endif
	const enum blas_base_type bzb = blas_zero_base;
	const enum blas_conj_type bc = blas_conj;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_susdot_(&bc, &nz, x, &(indx[0]), y, &incy, &dotr, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusdot(type, bj, nz, x, indx, y, incy, &dotr, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_susaxpy_(&nz,&alpha, x, &(indx[0]), y, &incy, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_susga_(&nz, x, &incy, y, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, bzb);
#endif
	errval=RSB_BLAS_ERROR_TO_RSB_ERROR(istat);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_susgz_(&nz, y, &incy, w, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_sussc_(&nz,x,y,&incy,&(indx[0]),&bzb,&istat);
#else
	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,bzb); 
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	/* BLAS types and interfaces only here  */
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double x[] = {+1,+0,+5,-2,+4,-6};
	double y[] = {+2,+0,+3,-2,+0,+0};
	double w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const double zero = ((double)(0));
	double alpha = 2;
	double dotr = zero;
	double res = zero;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	int istat;
#endif
	const enum blas_base_type bzb = blas_zero_base;
	const enum blas_conj_type bc = blas_conj;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_dusdot_(&bc, &nz, x, &(indx[0]), y, &incy, &dotr, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusdot(type, bj, nz, x, indx, y, incy, &dotr, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_dusaxpy_(&nz,&alpha, x, &(indx[0]), y, &incy, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_dusga_(&nz, x, &incy, y, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, bzb);
#endif
	errval=RSB_BLAS_ERROR_TO_RSB_ERROR(istat);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_dusgz_(&nz, y, &incy, w, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_dussc_(&nz,x,y,&incy,&(indx[0]),&bzb,&istat);
#else
	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,bzb); 
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	/* BLAS types and interfaces only here  */
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex x[] = {+1,+0,+5,-2,+4,-6};
	float complex y[] = {+2,+0,+3,-2,+0,+0};
	float complex w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const float complex zero = ((float complex)(0));
	float complex alpha = 2;
	float complex dotr = zero;
	float complex res = zero;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	int istat;
#endif
	const enum blas_base_type bzb = blas_zero_base;
	const enum blas_conj_type bc = blas_conj;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_cusdot_(&bc, &nz, x, &(indx[0]), y, &incy, &dotr, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusdot(type, bj, nz, x, indx, y, incy, &dotr, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_cusaxpy_(&nz,&alpha, x, &(indx[0]), y, &incy, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_cusga_(&nz, x, &incy, y, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, bzb);
#endif
	errval=RSB_BLAS_ERROR_TO_RSB_ERROR(istat);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_cusgz_(&nz, y, &incy, w, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_cussc_(&nz,x,y,&incy,&(indx[0]),&bzb,&istat);
#else
	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,bzb); 
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	/* BLAS types and interfaces only here  */
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex x[] = {+1,+0,+5,-2,+4,-6};
	double complex y[] = {+2,+0,+3,-2,+0,+0};
	double complex w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const double complex zero = ((double complex)(0));
	double complex alpha = 2;
	double complex dotr = zero;
	double complex res = zero;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	int istat;
#endif
	const enum blas_base_type bzb = blas_zero_base;
	const enum blas_conj_type bc = blas_conj;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_zusdot_(&bc, &nz, x, &(indx[0]), y, &incy, &dotr, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusdot(type, bj, nz, x, indx, y, incy, &dotr, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_zusaxpy_(&nz,&alpha, x, &(indx[0]), y, &incy, &bzb, &istat);
#else
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_zusga_(&nz, x, &incy, y, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, bzb);
#endif
	errval=RSB_BLAS_ERROR_TO_RSB_ERROR(istat);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_zusgz_(&nz, y, &incy, w, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	blas_zussc_(&nz,x,y,&incy,&(indx[0]),&bzb,&istat);
#else
	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,bzb); 
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */

{
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double x[] = {+1,+0,+5,-2,+4,-6};
	double y[] = {+2,+0,+3,-2,+0,+0};
	double w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const double zero = ((double)(0));
	double alpha = 2;
	double dotr = zero;
	double res = zero;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

	errval = rsb__BLAS_Xusdot(type, blas_conj, nz, x, indx, y, incy, &dotr, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,blas_zero_base); 
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT ;
	const float x[] = {+1,+0,+5,-2,+4,-6};
	float y[] = {+2,+0,+3,-2,+0,+0};
	float w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const float zero = ((float)(0));
	float alpha = 2;
	float dotr = zero;
	float res = zero;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

	errval = rsb__BLAS_Xusdot(type, blas_conj, nz, x, indx, y, incy, &dotr, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,blas_zero_base); 
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex x[] = {+1,+0,+5,-2,+4,-6};
	float complex y[] = {+2,+0,+3,-2,+0,+0};
	float complex w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const float complex zero = ((float complex)(0));
	float complex alpha = 2;
	float complex dotr = zero;
	float complex res = zero;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

	errval = rsb__BLAS_Xusdot(type, blas_conj, nz, x, indx, y, incy, &dotr, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,blas_zero_base); 
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const int incy = 1;
	const rsb_type_t type = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex x[] = {+1,+0,+5,-2,+4,-6};
	double complex y[] = {+2,+0,+3,-2,+0,+0};
	double complex w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const double complex zero = ((double complex)(0));
	double complex alpha = 2;
	double complex dotr = zero;
	double complex res = zero;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

	errval = rsb__BLAS_Xusdot(type, blas_conj, nz, x, indx, y, incy, &dotr, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,blas_zero_base); 
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double x[] = {+1,+0,+5,-2,+4,-6};
	double y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double alpha = 2;
	double res = ((double)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const double x[] = {+1,+0,+5,-2,+4,-6};
	float y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double alpha = 2;
	float res = ((float)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const double x[] = {+1,+0,+5,-2,+4,-6};
	float complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double alpha = 2;
	float complex res = ((float complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double x[] = {+1,+0,+5,-2,+4,-6};
	double complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double alpha = 2;
	double complex res = ((double complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const float x[] = {+1,+0,+5,-2,+4,-6};
	double y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float alpha = 2;
	double res = ((double)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const float x[] = {+1,+0,+5,-2,+4,-6};
	float y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float alpha = 2;
	float res = ((float)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float x[] = {+1,+0,+5,-2,+4,-6};
	float complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float alpha = 2;
	float complex res = ((float complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const float x[] = {+1,+0,+5,-2,+4,-6};
	double complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float alpha = 2;
	double complex res = ((double complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const float complex x[] = {+1,+0,+5,-2,+4,-6};
	double y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float complex alpha = 2;
	double res = ((double)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const float complex x[] = {+1,+0,+5,-2,+4,-6};
	float y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float complex alpha = 2;
	float res = ((float)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const float complex x[] = {+1,+0,+5,-2,+4,-6};
	float complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float complex alpha = 2;
	float complex res = ((float complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const float complex x[] = {+1,+0,+5,-2,+4,-6};
	double complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	float complex alpha = 2;
	double complex res = ((double complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE ;
	const double complex x[] = {+1,+0,+5,-2,+4,-6};
	double y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double complex alpha = 2;
	double res = ((double)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT ;
	const double complex x[] = {+1,+0,+5,-2,+4,-6};
	float y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double complex alpha = 2;
	float res = ((float)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_FLOAT_COMPLEX ;
	const double complex x[] = {+1,+0,+5,-2,+4,-6};
	float complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double complex alpha = 2;
	float complex res = ((float complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
{
	const rsb_type_t stypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const rsb_type_t dtypecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX ;
	const double complex x[] = {+1,+0,+5,-2,+4,-6};
	double complex y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	double complex alpha = 2;
	double complex res = ((double complex)(0));

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}



err:
	return errval;;
}


#ifdef __cplusplus
}
#endif  /* __cplusplus */

/* @endcond */
