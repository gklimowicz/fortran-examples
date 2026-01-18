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

 @brief A first "hello RSB" example C program based on <rsb.h>.
   Uses #rsb_lib_set_opt(), #rsb_mtx_get_info_str().

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
	  - set library options using #rsb_lib_set_opt()
	  - revert such changes 
	  - allocate (build) a single sparse matrix in the RSB format
	    using #rsb_mtx_alloc_from_coo_const()
	  - prints information obtained via #rsb_mtx_get_info_str()
	  - multiply the matrix times a vector using #rsb_spmv()
	  - deallocate the matrix using #rsb_mtx_free() 
	  - finalize the library using #rsb_lib_exit()
	 
	  In this example, we use #RSB_DEFAULT_TYPE as matrix type.
	  This type depends on what was configured at library build time.
	 * */
	const rsb_blk_idx_t bs = RSB_DEFAULT_BLOCKING;
	const rsb_blk_idx_t brA = bs, bcA = bs;
	const RSB_DEFAULT_TYPE one = 1;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnzA = 4;		/* matrix nonzeroes count */
	const rsb_coo_idx_t nrA = 3;		/* matrix rows count */
	const rsb_coo_idx_t ncA = 3;		/* matrix columns count */
	/* nonzero row indices coordinates: */
	const rsb_coo_idx_t IA[] = {0,1,2,2};
	/* nonzero column indices coordinates: */
	const rsb_coo_idx_t JA[] = {0,1,2,2};
	const RSB_DEFAULT_TYPE VA[] = {11,22,32,1};/* values of nonzeroes */
	RSB_DEFAULT_TYPE X[] = { 0, 0, 0 };	/* X vector's array */
	const RSB_DEFAULT_TYPE B[] = { -1, -2, -5 }; /* B vector's array */
	char ib[200];
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	printf("Hello, RSB!\n");
	printf("Initializing the library...\n");
	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != 
			RSB_ERR_NO_ERROR)
	{
		printf("Error initializing the library!\n");
		goto err;
	}
	printf("Correctly initialized the library.\n");

	printf("Attempting to set the"
	       " RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE library option.\n");
	{
		rsb_int_t evi=1; 

		/* Setting a single optional library parameter. */
		errval = rsb_lib_set_opt(
			RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE, &evi);
		if(errval != RSB_ERR_NO_ERROR)
		{
			char errbuf[256];
			rsb_strerror_r(errval,&errbuf[0],sizeof(errbuf));
			printf("Failed setting the"
			" RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE"
			" library option (reason string:\n%s).\n",errbuf);
			if(errval&RSB_ERRS_UNSUPPORTED_FEATURES)
			{
			  printf("This error may be safely ignored.\n");
			}
			else
			{
			  printf("Some unexpected error occurred!\n");
			  goto err;
			}
		}
		else
		{
			printf("Setting back the "
				"RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE"
				" library option.\n");
			evi = 0;
			errval = rsb_lib_set_opt(RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE,
					&evi);
			errval = RSB_ERR_NO_ERROR;
		}
	}

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_NOFLAGS    /* default format will be chosen */
		|RSB_FLAG_DUPLICATES_SUM/* duplicates will be summed */
			,&errval);
	if((!mtxAp) || (errval != RSB_ERR_NO_ERROR))
	{
		printf("Error while allocating the matrix!\n");
		goto err;
	}
	printf("Correctly allocated a matrix.\n");
	printf("Summary information of the matrix:\n");
	/* print out the matrix summary information  */
	rsb_mtx_get_info_str(mtxAp,"RSB_MIF_MATRIX_INFO__TO__CHAR_P",
			ib,sizeof(ib));
	printf("%s",ib);
	printf("\n");

	if((errval = 
		rsb_spmv(RSB_TRANSPOSITION_N,&one,mtxAp,B,1,&one,X,1))
			!= RSB_ERR_NO_ERROR )
	{
		printf("Error performing a multiplication!\n");
		goto err;
	}
	printf("Correctly performed a SPMV.\n");
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

