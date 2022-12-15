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
/*!
 \ingroup rsb_doc_examples
 @file
 @author Michele Martone

 @brief A first C "hello RSB" example program using
        a Sparse BLAS interface and <rsb.h>.
	Uses #BLAS_duscr_begin(), #BLAS_ussp(), #BLAS_usgp(),
	#BLAS_duscr_insert_entries(), #BLAS_duscr_end(),
	#BLAS_dusget_element(),#BLAS_dusmv(),#BLAS_usds().

 \include hello-spblas.c
*/
#include <rsb.h>	/* for rsb_lib_init */
#include <blas_sparse.h>	/* Sparse BLAS on the top of librsb */
#include <stdio.h>	/* printf */

int main(const int argc, char * const argv[])
{
	/*!
	 * A Hello/Sparse BLAS program.
	 *
	 * This program shows how to use the blas_sparse.h
	 * interface correctly to:
	 *
	 * - initialize the library using #rsb_lib_init()
	 * - allocate (build) a single sparse matrix in the RSB
	 *   format using #BLAS_duscr_begin()/#BLAS_duscr_insert_entries()
	 *   /#BLAS_duscr_end()
	 * - extract one matrix element with #BLAS_dusget_element()
	 * - multiply the matrix times a vector using #BLAS_dusmv()
	 * - deallocate the matrix using #BLAS_usds() 
	 * - finalize the library using
	 *   #rsb_lib_exit(#RSB_NULL_EXIT_OPTIONS) 
	*/
#ifndef RSB_NUMERICAL_TYPE_DOUBLE   
	printf("'double' type configured out."
	" Please reconfigure the library with it and recompile.\n");
	return EXIT_SUCCESS;
#else /* RSB_NUMERICAL_TYPE_DOUBLE */
	blas_sparse_matrix A = blas_invalid_handle; /* handle for A */
	const int nnz = 4;	/* number of nonzeroes of matrix A */
	const int  nr = 3;	/* number of A's rows */
	const int  nc = 3;	/* number of A's columns */
	/* A's nonzero elements row indices (coordinates): */
#ifdef RSB_WANT_LONG_IDX_TYPE 
	const int64_t IA[] = { 0, 1, 2, 2 };
#else /* RSB_WANT_LONG_IDX_TYPE */
	const int   IA[] = { 0, 1, 2, 2 };
#endif /* RSB_WANT_LONG_IDX_TYPE */
	/* A's nonzero elements column indices (coordinates): */
#ifdef RSB_WANT_LONG_IDX_TYPE 
	const int64_t JA[] = { 0, 1, 0, 2 };
#else /* RSB_WANT_LONG_IDX_TYPE */
	const int   JA[] = { 0, 1, 0, 2 };
#endif /* RSB_WANT_LONG_IDX_TYPE */
	/* A's nonzero values (matrix coefficients): */
	double VA[] = { 11.0, 22.0, 13.0, 33.0  };
       	/* the X vector's array: */
	double X[] = { 0.0, 0.0, 0.0 };
       	/* the B vector's array: */
	const double B[] = { -1.0, -2.0, -2.0 };
       	/* the (known) result array: */
	const double AB[] = { 11.0+26.0, 44.0, 66.0+13.0 };
	/* rsb error variable: */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	int i;

	printf("Hello, RSB!\n");
	/* initialize the library */
	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) 
			!= RSB_ERR_NO_ERROR)
	{
		goto err;
	}
	printf("Correctly initialized the library.\n");

	/* initialize a matrix descriptor */
	A = BLAS_duscr_begin(nr,nc);
	if( A == blas_invalid_handle )
	{
		goto err;
	}
	
	/* specify properties (e.g.: symmetry)*/
	if( BLAS_ussp(A,blas_lower_symmetric) != 0 )
	{
		goto err;
	}

	/* get properties (e.g.: symmetry) */
	if( BLAS_usgp(A,blas_lower_symmetric) != 1 )
	{
		printf("Symmetry property non set ?!\n");
		goto err;
	}

	/* insert the nonzeroes (here, all at once) */
	if( BLAS_duscr_insert_entries(A, nnz, VA, IA, JA)
			== blas_invalid_handle)
	{
		goto err;
	}

	/* finalize (allocate) the matrix build  */
	if( BLAS_duscr_end(A) == blas_invalid_handle )
	{
		goto err;
	}
	printf("Correctly allocated a matrix.\n");

	VA[0] = 0.0;
	if( BLAS_dusget_element(A, IA[0], JA[0], &VA[0]) )
	{
		goto err;
	}

	/* a check */
	if( VA[0] != 11.0 )
	{
		goto err;
	}

	/* compute X = X + (-1) * A * B   */
	if(BLAS_dusmv(blas_no_trans,-1,A,B,1,X,1))
	{
		goto err;
	}

	for( i = 0 ; i < nc; ++i )
		if( X[i] != AB[i] )
		{
			printf("Computed SPMV result seems wrong. Terminating.\n");
			goto err;
		}
	printf("Correctly performed a SPMV.\n");

	/* deallocate matrix A */
	if( BLAS_usds(A) )
	{
		goto err;
	}
	printf("Correctly freed the matrix.\n");

	/* finalize the library */
	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!= RSB_ERR_NO_ERROR)
	{
		goto err;
	}
	printf("Correctly finalized the library.\n");
	printf("Program terminating with no error.\n");

	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	printf("Program terminating with error.\n");
	return EXIT_FAILURE;
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
}

