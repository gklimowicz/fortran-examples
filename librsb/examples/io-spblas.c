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

 @brief Example C program using the Sparse BLAS interface
        and reading from file using \ref rsb_blas_file_mtx_load(),
	\ref BLAS_usgp(), \ref BLAS_dusmv(), \ref BLAS_usds().

 \include io-spblas.c
*/
#include <rsb.h>	/* for rsb_lib_init */
#include <blas_sparse.h>
#include <stdio.h>
	
int main(const int argc, char * const argv[])
{
#ifndef RSB_NUMERICAL_TYPE_DOUBLE   
	printf("Skipping a test because of 'double' type opted out.\n");
	return EXIT_SUCCESS;
#else /* RSB_NUMERICAL_TYPE_DOUBLE */
	blas_sparse_matrix A = blas_invalid_handle;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
	const rsb_char_t * filename = argc > 1 ? argv[1] : "pd.mtx";

	printf("Hello, RSB!\n");
	if((rsb_perror(NULL,
		rsb_lib_init(RSB_NULL_INIT_OPTIONS)))!=RSB_ERR_NO_ERROR)
	{
		printf("Error while initializing the library.\n");
		goto err;
	}

	printf("Correctly initialized the library.\n");

	A = rsb_blas_file_mtx_load(filename, typecode );
	if( A == blas_invalid_handle )
	{
		printf("Error while loading matrix %s from file.\n",
				filename);
		goto err;
	}

	printf("Correctly loaded and allocated a matrix"
			" from file %s.\n",filename);

	if( BLAS_usgp(A,blas_symmetric) == 1 )
		printf("Matrix is symmetric\n");

	if( BLAS_usgp(A,blas_hermitian) == 1 )
		printf("Matrix is hermitian\n");

	printf("Now SPMV with NULL vectors will be attempted,"
			" resulting in an error (so don't worry).\n");

	if(BLAS_dusmv(blas_no_trans,-1,A,NULL,1,NULL,1))
	{
		printf("Correctly detected an error condition.\n");
		goto okerr;
	}

	printf("No error detected ?\nIf you see this line printed out,"
		" please report	as a bug, because the above NULL pointers"
		" should have been detected\n");
	return EXIT_FAILURE;
okerr:
	printf("Program correctly recovered from intentional"
			" error condition.\n");
	if(BLAS_usds(A))
	{
		printf("Error while freeing the matrix!\n");
		goto err;
	}

	printf("Correctly freed the matrix.\n");
err:
	if(rsb_perror(NULL,
		rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))!=RSB_ERR_NO_ERROR)
	{
		printf("Failed finalizing the library.\n");
		goto ferr;
	}

	printf("Correctly finalized the library.\n");
	return EXIT_SUCCESS;
ferr:
	return EXIT_FAILURE;
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
}

