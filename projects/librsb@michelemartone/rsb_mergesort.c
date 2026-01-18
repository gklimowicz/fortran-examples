

/* @cond INNERDOC */
/**
 * @file
 * @brief
 * Sorting functions.
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


#include "rsb.h"
#include "rsb_common.h"
#include "rsb_internals.h"


#if RSB_OBSOLETE_QUARANTINE

RSB_INTERNALS_COMMON_HEAD_DECLS

rsb_err_t rsb__do_mergesort_CSR(
	rsb_coo_idx_t *iarray,
	rsb_coo_idx_t *jarray,
	void *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *iresult,
	rsb_coo_idx_t *jresult,
	void *result,
	rsb_type_t type)
{
	/*!
	 * \ingroup gr_util
	 *	This function will sort the non zero elements of a sparse matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 * 	\param length  the input  arrays length
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  mtype array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	FIXME : UNDOCUMENTED
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */


	if(type == RSB_NUMERICAL_TYPE_DOUBLE )
	return rsb_do_mergesort_double_CSR(iarray, jarray,
array, length,
iresult, jresult,
result);
	else
	if(type == RSB_NUMERICAL_TYPE_FLOAT )
	return rsb_do_mergesort_float_CSR(iarray, jarray,
array, length,
iresult, jresult,
result);
	else
	if(type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	return rsb_do_mergesort_float_complex_CSR(iarray, jarray,
array, length,
iresult, jresult,
result);
	else
	if(type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	return rsb_do_mergesort_double_complex_CSR(iarray, jarray,
array, length,
iresult, jresult,
result);
	else
	return RSB_ERR_UNSUPPORTED_TYPE;
}

rsb_err_t rsb__do_mergesort_BCSR(
	rsb_coo_idx_t *iarray,
	rsb_coo_idx_t *jarray,
	void *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t mb,
	rsb_coo_idx_t kb,
	rsb_coo_idx_t *iresult,
	rsb_coo_idx_t *jresult,
	void *result,
	rsb_type_t type)
{
	/*!
	 * \ingroup gr_util
	 *	This function will sort the non zero elements of a sparse matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 * 	\param length  the input  arrays length
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  mtype array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array

	 *	FIXME : UNDOCUMENTED
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */


	if(type == RSB_NUMERICAL_TYPE_DOUBLE )
	return rsb_do_mergesort_double_BCSR(iarray, jarray,
	mb,kb,array, length,
iresult, jresult,
result);
	else
	if(type == RSB_NUMERICAL_TYPE_FLOAT )
	return rsb_do_mergesort_float_BCSR(iarray, jarray,
	mb,kb,array, length,
iresult, jresult,
result);
	else
	if(type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	return rsb_do_mergesort_float_complex_BCSR(iarray, jarray,
	mb,kb,array, length,
iresult, jresult,
result);
	else
	if(type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	return rsb_do_mergesort_double_complex_BCSR(iarray, jarray,
	mb,kb,array, length,
iresult, jresult,
result);
	else
	return RSB_ERR_UNSUPPORTED_TYPE;
}

rsb_err_t rsb__do_mergesort_VBR(
	rsb_coo_idx_t *iarray,
	rsb_coo_idx_t *jarray,
	rsb_coo_idx_t *biarray,
	rsb_coo_idx_t *bjarray,
	void *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *iresult,
	rsb_coo_idx_t *jresult,
	rsb_coo_idx_t *biresult,
	rsb_coo_idx_t *bjresult,
	void *result,
	rsb_type_t type)
{
	/*!
	 * \ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 * 	\param length  the input  arrays length
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  mtype array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 * 	\param biarray  the input  block row    indices array
	 * 	\param bjarray  the input  block column indices array
	 * 	\param biresult the output block row    indices array
	 * 	\param bjresult the output block column indices array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */


	if(type == RSB_NUMERICAL_TYPE_DOUBLE )
	return rsb_do_mergesort_double_VBR(iarray, jarray,
	biarray,bjarray,array, length,
iresult, jresult,
	biresult,bjresult,result);
	else
	if(type == RSB_NUMERICAL_TYPE_FLOAT )
	return rsb_do_mergesort_float_VBR(iarray, jarray,
	biarray,bjarray,array, length,
iresult, jresult,
	biresult,bjresult,result);
	else
	if(type == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
	return rsb_do_mergesort_float_complex_VBR(iarray, jarray,
	biarray,bjarray,array, length,
iresult, jresult,
	biresult,bjresult,result);
	else
	if(type == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
	return rsb_do_mergesort_double_complex_VBR(iarray, jarray,
	biarray,bjarray,array, length,
iresult, jresult,
	biresult,bjresult,result);
	else
	return RSB_ERR_UNSUPPORTED_TYPE;
}

rsb_err_t rsb_do_mergesort_double_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	double *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      double matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  double array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	double * left  ;
	double * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(double*)result = *(double*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_double_CSR
	( ileft, jleft,
		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_double_CSR
	(iright, jright,
		right, length-middle,  iresult+middle  ,jresult+middle,
	((double*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(double)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((double*)result)+middle ,sizeof(double)*(length-middle));

	rsb_do_merge_double_CSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_double_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const double* left, const double* restrict right,  double* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our CSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
		
		ileft[left_index] < iright[right_index] ||
		(	ileft[left_index] == iright[right_index]	&&
			jleft[left_index] <= jright[right_index]	)
		)
				{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_float_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	float *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      float matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  float array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	float * left  ;
	float * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(float*)result = *(float*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_float_CSR
	( ileft, jleft,
		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_float_CSR
	(iright, jright,
		right, length-middle,  iresult+middle  ,jresult+middle,
	((float*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(float)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((float*)result)+middle ,sizeof(float)*(length-middle));

	rsb_do_merge_float_CSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_float_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const float* left, const float* restrict right,  float* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our CSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
		
		ileft[left_index] < iright[right_index] ||
		(	ileft[left_index] == iright[right_index]	&&
			jleft[left_index] <= jright[right_index]	)
		)
				{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_float_complex_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	float complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float complex *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      float complex matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  float complex array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	float complex * left  ;
	float complex * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(float complex*)result = *(float complex*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_float_complex_CSR
	( ileft, jleft,
		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_float_complex_CSR
	(iright, jright,
		right, length-middle,  iresult+middle  ,jresult+middle,
	((float complex*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(float complex)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((float complex*)result)+middle ,sizeof(float complex)*(length-middle));

	rsb_do_merge_float_complex_CSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_float_complex_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const float complex* left, const float complex* restrict right,  float complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our CSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
		
		ileft[left_index] < iright[right_index] ||
		(	ileft[left_index] == iright[right_index]	&&
			jleft[left_index] <= jright[right_index]	)
		)
				{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_double_complex_CSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	double complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double complex *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      double complex matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  double complex array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	double complex * left  ;
	double complex * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(double complex*)result = *(double complex*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_double_complex_CSR
	( ileft, jleft,
		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_double_complex_CSR
	(iright, jright,
		right, length-middle,  iresult+middle  ,jresult+middle,
	((double complex*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(double complex)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((double complex*)result)+middle ,sizeof(double complex)*(length-middle));

	rsb_do_merge_double_complex_CSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_double_complex_CSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
		const double complex* left, const double complex* restrict right,  double complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our CSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
		
		ileft[left_index] < iright[right_index] ||
		(	ileft[left_index] == iright[right_index]	&&
			jleft[left_index] <= jright[right_index]	)
		)
				{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_double_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	double *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      double matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  double array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	double * left  ;
	double * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(double*)result = *(double*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_double_BCSR
	( ileft, jleft,
 mb, kb, 		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_double_BCSR
	(iright, jright,
 mb, kb, 		right, length-middle,  iresult+middle  ,jresult+middle,
	((double*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(double)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((double*)result)+middle ,sizeof(double)*(length-middle));

	rsb_do_merge_double_BCSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,
	mb,kb,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_double_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const double* left, const double* restrict right,  double* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our BCSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
			
		ileft[left_index]/mb < iright[right_index]/mb ||
		(	ileft[left_index]/mb == iright[right_index]/mb	&&
			jleft[left_index]/kb <= jright[right_index]/kb	)
		)
			{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_float_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	float *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      float matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  float array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	float * left  ;
	float * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(float*)result = *(float*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_float_BCSR
	( ileft, jleft,
 mb, kb, 		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_float_BCSR
	(iright, jright,
 mb, kb, 		right, length-middle,  iresult+middle  ,jresult+middle,
	((float*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(float)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((float*)result)+middle ,sizeof(float)*(length-middle));

	rsb_do_merge_float_BCSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,
	mb,kb,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_float_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const float* left, const float* restrict right,  float* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our BCSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
			
		ileft[left_index]/mb < iright[right_index]/mb ||
		(	ileft[left_index]/mb == iright[right_index]/mb	&&
			jleft[left_index]/kb <= jright[right_index]/kb	)
		)
			{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_float_complex_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	float complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	float complex *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      float complex matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  float complex array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	float complex * left  ;
	float complex * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(float complex*)result = *(float complex*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_float_complex_BCSR
	( ileft, jleft,
 mb, kb, 		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_float_complex_BCSR
	(iright, jright,
 mb, kb, 		right, length-middle,  iresult+middle  ,jresult+middle,
	((float complex*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(float complex)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((float complex*)result)+middle ,sizeof(float complex)*(length-middle));

	rsb_do_merge_float_complex_BCSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,
	mb,kb,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_float_complex_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const float complex* left, const float complex* restrict right,  float complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our BCSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
			
		ileft[left_index]/mb < iright[right_index]/mb ||
		(	ileft[left_index]/mb == iright[right_index]/mb	&&
			jleft[left_index]/kb <= jright[right_index]/kb	)
		)
			{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_double_complex_BCSR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
	double complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	double complex *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      double complex matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  double complex array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
	size_t tn=0;
	double complex * left  ;
	double complex * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
		*(double complex*)result = *(double complex*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_double_complex_BCSR
	( ileft, jleft,
 mb, kb, 		left,   middle,
	        iresult  ,       jresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_double_complex_BCSR
	(iright, jright,
 mb, kb, 		right, length-middle,  iresult+middle  ,jresult+middle,
	((double complex*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(double complex)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
	RSB_MEMCPY( right, ((double complex*)result)+middle ,sizeof(double complex)*(length-middle));

	rsb_do_merge_double_complex_BCSR		(
			ileft,iright,iresult,
			jleft,jright,jresult,
	mb,kb,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_double_complex_BCSR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,		const double complex* left, const double complex* restrict right,  double complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our BCSR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
			
		ileft[left_index]/mb < iright[right_index]/mb ||
		(	ileft[left_index]/mb == iright[right_index]/mb	&&
			jleft[left_index]/kb <= jright[right_index]/kb	)
		)
			{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_double_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	double *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	double *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      double matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  double array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 * 	\param biarray  the input  block row    indices array
	 * 	\param bjarray  the input  block column indices array
	 * 	\param biresult the output block row    indices array
	 * 	\param bjresult the output block column indices array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;

	rsb_coo_idx_t * bileft  ;
	rsb_coo_idx_t * biright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;

	rsb_coo_idx_t * bjleft  ;
	rsb_coo_idx_t * bjright ;
	size_t tn=0;
	double * left  ;
	double * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;

		*biresult = *biarray;
		*bjresult = *bjarray;
		*(double*)result = *(double*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;


	bileft  = biarray;
	bjleft  = bjarray;
	biright = biarray+middle;
	bjright = bjarray+middle;
	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_double_VBR
	( ileft, jleft,
		bileft, bjleft,
		left,   middle,
	        iresult  ,       jresult,
		biresult, bjresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_double_VBR
	(iright, jright,
		biright, bjright,
		right, length-middle,  iresult+middle  ,jresult+middle,
		biresult+middle, bjresult+middle,
	((double*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);

	RSB_MEMCPY(bileft ,biresult       ,so*middle);
	RSB_MEMCPY(bjleft ,bjresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(double)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));

	RSB_MEMCPY(biright ,biresult+middle       ,so*(length-middle));
	RSB_MEMCPY(bjright ,bjresult+middle       ,so*(length-middle));
	RSB_MEMCPY( right, ((double*)result)+middle ,sizeof(double)*(length-middle));

	rsb_do_merge_double_VBR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			bileft, biright, biresult,
			bjleft, bjright, bjresult,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_double_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const double* left, const double* restrict right,  double* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our VBR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,BIEL,BJEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	biresult[result_index]=(BIEL);  \
	bjresult[result_index]=(BJEL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	biresult[result_index]=bileft[left_index];  \
	bjresult[result_index]=bjleft[left_index];  \
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	biresult[result_index]=biright[right_index];  \
	bjresult[result_index]=bjright[right_index];  \
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
	
		bileft[left_index] < biright[right_index] ||
		(	bileft[left_index] == biright[right_index]	&&
			bjleft[left_index] <= bjright[right_index]	)
		)
					{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_float_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	float *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	float *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      float matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  float array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 * 	\param biarray  the input  block row    indices array
	 * 	\param bjarray  the input  block column indices array
	 * 	\param biresult the output block row    indices array
	 * 	\param bjresult the output block column indices array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;

	rsb_coo_idx_t * bileft  ;
	rsb_coo_idx_t * biright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;

	rsb_coo_idx_t * bjleft  ;
	rsb_coo_idx_t * bjright ;
	size_t tn=0;
	float * left  ;
	float * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;

		*biresult = *biarray;
		*bjresult = *bjarray;
		*(float*)result = *(float*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;


	bileft  = biarray;
	bjleft  = bjarray;
	biright = biarray+middle;
	bjright = bjarray+middle;
	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_float_VBR
	( ileft, jleft,
		bileft, bjleft,
		left,   middle,
	        iresult  ,       jresult,
		biresult, bjresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_float_VBR
	(iright, jright,
		biright, bjright,
		right, length-middle,  iresult+middle  ,jresult+middle,
		biresult+middle, bjresult+middle,
	((float*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);

	RSB_MEMCPY(bileft ,biresult       ,so*middle);
	RSB_MEMCPY(bjleft ,bjresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(float)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));

	RSB_MEMCPY(biright ,biresult+middle       ,so*(length-middle));
	RSB_MEMCPY(bjright ,bjresult+middle       ,so*(length-middle));
	RSB_MEMCPY( right, ((float*)result)+middle ,sizeof(float)*(length-middle));

	rsb_do_merge_float_VBR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			bileft, biright, biresult,
			bjleft, bjright, bjresult,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_float_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const float* left, const float* restrict right,  float* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our VBR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,BIEL,BJEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	biresult[result_index]=(BIEL);  \
	bjresult[result_index]=(BJEL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	biresult[result_index]=bileft[left_index];  \
	bjresult[result_index]=bjleft[left_index];  \
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	biresult[result_index]=biright[right_index];  \
	bjresult[result_index]=bjright[right_index];  \
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
	
		bileft[left_index] < biright[right_index] ||
		(	bileft[left_index] == biright[right_index]	&&
			bjleft[left_index] <= bjright[right_index]	)
		)
					{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_float_complex_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	float complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	float complex *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      float complex matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  float complex array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 * 	\param biarray  the input  block row    indices array
	 * 	\param bjarray  the input  block column indices array
	 * 	\param biresult the output block row    indices array
	 * 	\param bjresult the output block column indices array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;

	rsb_coo_idx_t * bileft  ;
	rsb_coo_idx_t * biright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;

	rsb_coo_idx_t * bjleft  ;
	rsb_coo_idx_t * bjright ;
	size_t tn=0;
	float complex * left  ;
	float complex * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;

		*biresult = *biarray;
		*bjresult = *bjarray;
		*(float complex*)result = *(float complex*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;


	bileft  = biarray;
	bjleft  = bjarray;
	biright = biarray+middle;
	bjright = bjarray+middle;
	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_float_complex_VBR
	( ileft, jleft,
		bileft, bjleft,
		left,   middle,
	        iresult  ,       jresult,
		biresult, bjresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_float_complex_VBR
	(iright, jright,
		biright, bjright,
		right, length-middle,  iresult+middle  ,jresult+middle,
		biresult+middle, bjresult+middle,
	((float complex*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);

	RSB_MEMCPY(bileft ,biresult       ,so*middle);
	RSB_MEMCPY(bjleft ,bjresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(float complex)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));

	RSB_MEMCPY(biright ,biresult+middle       ,so*(length-middle));
	RSB_MEMCPY(bjright ,bjresult+middle       ,so*(length-middle));
	RSB_MEMCPY( right, ((float complex*)result)+middle ,sizeof(float complex)*(length-middle));

	rsb_do_merge_float_complex_VBR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			bileft, biright, biresult,
			bjleft, bjright, bjresult,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_float_complex_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const float complex* left, const float complex* restrict right,  float complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our VBR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,BIEL,BJEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	biresult[result_index]=(BIEL);  \
	bjresult[result_index]=(BJEL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	biresult[result_index]=bileft[left_index];  \
	bjresult[result_index]=bjleft[left_index];  \
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	biresult[result_index]=biright[right_index];  \
	bjresult[result_index]=bjright[right_index];  \
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
	
		bileft[left_index] < biright[right_index] ||
		(	bileft[left_index] == biright[right_index]	&&
			bjleft[left_index] <= bjright[right_index]	)
		)
					{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}

rsb_err_t rsb_do_mergesort_double_complex_VBR(
	rsb_coo_idx_t *restrict iarray,
	rsb_coo_idx_t *restrict jarray,
	rsb_coo_idx_t *restrict biarray,
	rsb_coo_idx_t *restrict bjarray,
	double complex *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *restrict iresult,
	rsb_coo_idx_t *restrict jresult,
	rsb_coo_idx_t *restrict biresult,
	rsb_coo_idx_t *restrict bjresult,
	double complex *restrict result)

{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      double complex matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  double complex array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
	 * 	\param biarray  the input  block row    indices array
	 * 	\param bjarray  the input  block column indices array
	 * 	\param biresult the output block row    indices array
	 * 	\param bjresult the output block column indices array
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;

	rsb_coo_idx_t * bileft  ;
	rsb_coo_idx_t * biright ;
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;

	rsb_coo_idx_t * bjleft  ;
	rsb_coo_idx_t * bjright ;
	size_t tn=0;
	double complex * left  ;
	double complex * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;

		*biresult = *biarray;
		*bjresult = *bjarray;
		*(double complex*)result = *(double complex*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;


	bileft  = biarray;
	bjleft  = bjarray;
	biright = biarray+middle;
	bjright = bjarray+middle;
	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

/* 20121016 commented out omp usage because broke serial compilation  */
	{
	rsb_do_mergesort_double_complex_VBR
	( ileft, jleft,
		bileft, bjleft,
		left,   middle,
	        iresult  ,       jresult,
		biresult, bjresult,
		result         );

	if(tn==1)
	rsb_do_mergesort_double_complex_VBR
	(iright, jright,
		biright, bjright,
		right, length-middle,  iresult+middle  ,jresult+middle,
		biresult+middle, bjresult+middle,
	((double complex*)result)+middle  );
	}

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);

	RSB_MEMCPY(bileft ,biresult       ,so*middle);
	RSB_MEMCPY(bjleft ,bjresult       ,so*middle);
	RSB_MEMCPY(  left, result       ,sizeof(double complex)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));

	RSB_MEMCPY(biright ,biresult+middle       ,so*(length-middle));
	RSB_MEMCPY(bjright ,bjresult+middle       ,so*(length-middle));
	RSB_MEMCPY( right, ((double complex*)result)+middle ,sizeof(double complex)*(length-middle));

	rsb_do_merge_double_complex_VBR		(
			ileft,iright,iresult,
			jleft,jright,jresult,

			bileft, biright, biresult,
			bjleft, bjright, bjresult,
			left, right, result,
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}



void rsb_do_merge_double_complex_VBR(
		const rsb_coo_idx_t* restrict ileft, const rsb_coo_idx_t* restrict iright,  rsb_coo_idx_t*restrict iresult,
		const rsb_coo_idx_t* restrict jleft, const rsb_coo_idx_t* restrict jright,  rsb_coo_idx_t*restrict jresult,
const rsb_coo_idx_t * restrict bileft, const rsb_coo_idx_t * restrict biright, rsb_coo_idx_t * restrict biresult,const rsb_coo_idx_t * restrict bjleft, const rsb_coo_idx_t * restrict bjright, rsb_coo_idx_t * restrict bjresult,		const double complex* left, const double complex* restrict right,  double complex* restrict result,
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )


{
	/*!
	 * \ingroup gr_util
	 * The merge function for our VBR matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
	          ^- index
	 */


#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

#define RESULT_APPEND(IEL,JEL,BIEL,BJEL,EL)	\
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
	biresult[result_index]=(BIEL);  \
	bjresult[result_index]=(BJEL);  \
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
	biresult[result_index]=bileft[left_index];  \
	bjresult[result_index]=bjleft[left_index];  \
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
	biresult[result_index]=biright[right_index];  \
	bjresult[result_index]=bjright[right_index];  \
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
	
		bileft[left_index] < biright[right_index] ||
		(	bileft[left_index] == biright[right_index]	&&
			bjleft[left_index] <= bjright[right_index]	)
		)
					{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}


#endif /* RSB_OBSOLETE_QUARANTINE */

#ifdef __cplusplus
}
#endif  /* __cplusplus */

/* @endcond */
