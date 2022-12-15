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

 @brief C triangular solve example program.
   Uses #rsb_spsv(),#rsb_tune_spsm().
   Based on <rsb.h>.

 \include hello.c
*/
#include <rsb.h>	/* librsb header to include */
#include <stdio.h>	/* printf() */

int main(const int argc, char * const argv[])
{
	/*!
	  A Hello-RSB program.
	 
	  This program shows how to use the rsb.h interface correctly to:
	 
	  - initialize the library using #rsb_lib_init()
	  - allocate (build) a single sparse matrix in the RSB format
	    using #rsb_mtx_alloc_from_coo_const(), with implicit diagonal
	  - print information obtained via #rsb_mtx_get_info_str()
	  - multiply the triangular matrix using #rsb_spmv()
	  - autotune the matrix for rsb_spsv with #rsb_tune_spsm()
	  - solve the triangular system using #rsb_spsv()
	  - deallocate the matrix using #rsb_mtx_free() 
	  - finalize the library using #rsb_lib_exit()
	 
	  In this example, we use #RSB_DEFAULT_TYPE as matrix type.
	  This type depends on what was configured at library build time.
	 * */
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const RSB_DEFAULT_TYPE one = 1;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnzA = 7;		/* matrix nonzeroes count */
	const rsb_coo_idx_t nrA = 6;		/* matrix rows count */
	const rsb_coo_idx_t ncA = 6;		/* matrix columns count */
	/* nonzero row indices coordinates: */
	const rsb_coo_idx_t IA[] = {1,2,3,4,5,6,1};
	/* nonzero column indices coordinates: */
	const rsb_coo_idx_t JA[] = {1,2,3,4,5,6,6};
	const RSB_DEFAULT_TYPE VA[] = {1,1,1,1,1,1,1};/*values of nonzeroes*/
	RSB_DEFAULT_TYPE X[] = { 0,0,0,0,0,0 };	/* X vector's array */
	const RSB_DEFAULT_TYPE B[] = { 1,1,1,1,1,1 }; /* B */
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */
	char ib[200];
	int i;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_int_t wvat = 1; /* want verbose autotuning; see 
                  documentation of RSB_IO_WANT_VERBOSE_TUNING */

	printf("Hello, RSB!\n");
	printf("Initializing the library...\n");
	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != 
			RSB_ERR_NO_ERROR)
	{
		printf("Error initializing the library!\n");
		goto err;
	}
	printf("Correctly initialized the library.\n");

	errval = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING, &wvat );
	if( (errval) != RSB_ERR_NO_ERROR )
	{
		printf("Error setting option!\n");
		goto err;
	}

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS /* force rsb */
		| RSB_FLAG_DUPLICATES_SUM/* sum dups */
		| RSB_FLAG_UNIT_DIAG_IMPLICIT/* ask diagonal implicit */
		| RSB_FLAG_TRIANGULAR /* need triangle for spsv */
		| RSB_FLAG_FORTRAN_INDICES_INTERFACE /* treat indices as 1-based */
		, &errval);
	if((!mtxAp) || (errval != RSB_ERR_NO_ERROR))
	{
		printf("Error while allocating the matrix!\n");
		goto err;
	}
	printf("Correctly allocated a matrix with %ld nonzeroes.\n",
		(long int)nnzA);
	printf("Summary information of the matrix:\n");
	/* print out the matrix summary information  */
	rsb_mtx_get_info_str(mtxAp,"RSB_MIF_MATRIX_INFO__TO__CHAR_P",
			ib,sizeof(ib));
	printf("%s",ib);
	printf("\nMatrix printout:\n");

	rsb_file_mtx_save(mtxAp, NULL);

	if((errval = 
		rsb_spmv(RSB_TRANSPOSITION_N,&one,mtxAp,B,1,&one,X,1))
			!= RSB_ERR_NO_ERROR )
	{
		printf("Error performing a multiplication!\n");
		goto err;
	}
	printf("\nWe have a unitary vector:\n");
	rsb_file_vec_save(NULL, typecode, B, nrA);
	
	printf("\nMultiplying matrix by unitary vector we get:\n");
	rsb_file_vec_save(NULL, typecode, X, nrA);

	errval = rsb_tune_spsm(&mtxAp, NULL, NULL, 0, 0, RSB_TRANSPOSITION_N,
	       	&one, NULL, 1, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, NULL, nrA,
		 NULL, NULL, nrA);
	if( (errval) != RSB_ERR_NO_ERROR )
	{
		printf("Error performing autotuning!\n");
		goto err;
	}

	if((errval = rsb_spsv(RSB_TRANSPOSITION_N,&one,mtxAp,X,1,X,1))
			!= RSB_ERR_NO_ERROR )
	{
		printf("Error performing triangular solve!\n");
		goto err;
	}

	printf("\nBacksolving we should get a unitary vector:\n");
	rsb_file_vec_save(NULL, typecode, X, nrA);

	for(i=0;i<nrA;++i)
		if(X[i]!=one)
		{
			printf("Warning! Result vector not unitary!:\n");
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;
		}

	printf("All done.\n");
	rsb_mtx_free(mtxAp);
	printf("Correctly freed the matrix.\n");
	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!= RSB_ERR_NO_ERROR)
	{
		printf("Error finalizing the library!\n");
		goto err;
	}
	printf("Correctly finalized the library.\n");
	printf("Program terminating with no error.\n");
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	printf("Program terminating with error.\n");
	return EXIT_FAILURE;
}

