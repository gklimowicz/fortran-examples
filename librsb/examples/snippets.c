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

 @brief Collection of C snippets of other examples.
        Used piecewise the documentation.
        Not intended to be read as example.
*/
#include <rsb.h>	/* librsb header to include */
#include <stdio.h>	/* printf() */
#include <ctype.h>	/* isdigit() */
#include <stdlib.h>	/* atoi() */
#include <string.h>	/* strstr() */

static int get_coo_block_snippet(struct rsb_mtx_t *mtxAp,
		rsb_nnz_idx_t nnzA )
{
	/*! [Extract one sparse matrix block] */
	rsb_coo_idx_t nzi;
	rsb_coo_idx_t *IA = NULL;
	rsb_coo_idx_t *JA = NULL;
	const rsb_coo_idx_t IREN[]={0,1,2,3};
	const rsb_coo_idx_t JREN[]={3,2,1,0};
	RSB_DEFAULT_TYPE *VA = NULL;
	const size_t so = sizeof(RSB_DEFAULT_TYPE);
	const size_t si = sizeof(rsb_coo_idx_t);
	rsb_err_t errval;
	rsb_flags_t flagsA = RSB_FLAG_NOFLAGS;
	rsb_nnz_idx_t rnz = 0;
	rsb_coo_idx_t frA=0,lrA=1; // first two rows
	rsb_coo_idx_t fcA=0,lcA=4; // 5 (all) columns

	// get the nnz count only
	errval=rsb_mtx_get_coo_block
		(mtxAp,NULL,NULL,NULL,frA,lrA,fcA,lcA,NULL,NULL,&rnz,flagsA);
	if(errval != RSB_ERR_NO_ERROR )
		goto err;

	// allocate VA, IA, JA to rnz elements
	IA = calloc(rnz, si);
	JA = calloc(rnz, si);
	VA = calloc(rnz, so);

	// get the  rnz  values then
	errval=rsb_mtx_get_coo_block
		(mtxAp,  VA,  IA,  JA,frA,lrA,fcA,lcA,NULL,NULL,NULL,flagsA);
	if(errval != RSB_ERR_NO_ERROR )
		goto err;

	for(nzi=0;nzi<rnz;++nzi)
		printf("%d/%d  %d %d -> %d\n",(int)nzi,(int)rnz,
			(int)IA[nzi],(int)JA[nzi],(int)VA[nzi]);

	// get the  rnz  values again, renumbered
	errval=rsb_mtx_get_coo_block
		(mtxAp,  VA,  IA,  JA,frA,lrA,fcA,lcA,IREN,JREN,NULL,flagsA);
	if(errval != RSB_ERR_NO_ERROR )
		goto err;

	for(nzi=0;nzi<rnz;++nzi)
		printf("%d/%d  %d %d -> %d\n",(int)nzi,(int)rnz,
			(int)IA[nzi],(int)JA[nzi],(int)VA[nzi]);

	free(VA);
	free(IA);
	free(JA);
	/*! [Extract one sparse matrix block] */
err:
	return errval;
}

static int main_backsolve(const int argc, char * const argv[])
{
	/*!
	  A Hello-RSB program.
	 
	  This program shows how to use the rsb.h interface correctly to:
	 
	  - initialize the library using #rsb_lib_init()
	  - allocate (build) a single sparse matrix in the RSB format
	    using #rsb_mtx_alloc_from_coo_const(), with implicit diagonal
	  - print information obtained via #rsb_mtx_get_info_str()
	  - multiply the triangular matrix using #rsb_spmv()
	  - solve the triangular system using #rsb_spsv()
	  - deallocate the matrix using #rsb_mtx_free() 
	  - finalize the library using #rsb_lib_exit(RSB_NULL_EXIT_OPTIONS) 
	 
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
	const rsb_coo_idx_t IA[] = {0,1,2,3,4,5,1};
	/* nonzero column indices coordinates: */
	const rsb_coo_idx_t JA[] = {0,1,2,3,4,5,5};
	const RSB_DEFAULT_TYPE VA[] = {1,1,1,1,1,1,1};/*values of nonzeroes*/
	RSB_DEFAULT_TYPE X[] = { 0,0,0,0,0,0 };	/* X vector's array */
	const RSB_DEFAULT_TYPE B[] = { 1,1,1,1,1,1 }; /* B */
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */
	char ib[200];
	int i;
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

	/*! [Allocate a matrix with triangular flags] */
	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS /* force rsb */
		| RSB_FLAG_DUPLICATES_SUM/* sum dups */
		| RSB_FLAG_UNIT_DIAG_IMPLICIT/* ask diagonal implicit */
		| RSB_FLAG_TRIANGULAR /* need triangle for spsv */
		, &errval);
	if((!mtxAp) || (errval != RSB_ERR_NO_ERROR))
	{
		printf("Error while allocating the matrix!\n");
		goto err;
	}
	printf("Correctly allocated a matrix with %ld nonzeroes.\n",
		(long int)nnzA);
	/*! [Allocate a matrix with triangular flags] */

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

	/*! [Backsolve a triangular system] */
	if((errval = rsb_spsv(RSB_TRANSPOSITION_N,&one,mtxAp,X,1,X,1))
			!= RSB_ERR_NO_ERROR )
	{
		printf("Error performing triangular solve!\n");
		goto err;
	}
	/*! [Backsolve a triangular system] */

	printf("\nBacksolving we should get a unitary vector:\n");

	/*! [Print to stdout an nrA-long numerical vector] */
	errval = rsb_file_vec_save(NULL, typecode, X, nrA);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error printing vector!\n");
		goto err;
	}
	/*! [Print to stdout an nrA-long numerical vector] */

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

static int hello_snip(const int argc, char * const argv[])
{
	/**
	  This program shows how to use the rsb.h interface correctly to:
	 
	  - initialize the library using #rsb_lib_init()
	  - set library options using #rsb_lib_set_opt()
	  - revert such changes 
	  - allocate (build) a single sparse matrix in the RSB format
	    using #rsb_mtx_alloc_from_coo_const()
	  - prints information obtained via #rsb_mtx_get_info_str()
	  - multiply the matrix times a vector using #rsb_spmv()
	  - deallocate the matrix using #rsb_mtx_free() 
	  - finalize the library using #rsb_lib_exit(RSB_NULL_EXIT_OPTIONS) 
	 
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
        /*! [Declare error codes variable] */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        /*! [Declare error codes variable] */

	printf("Hello, RSB!\n");
	printf("Initializing the library...\n");
        /*! [Initialize the library] */
	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != 
			RSB_ERR_NO_ERROR)
	{
		printf("Error initializing the library!\n");
		goto err;
	}
        /*! [Initialize the library] */
	printf("Correctly initialized the library.\n");

	printf("Attempting to set the"
	       " RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE library option.\n");
	{
		/*! [Setting a single optional library parameter] */
		rsb_int_t evi=1; 

		/* Setting a single optional library parameter. */
		errval = rsb_lib_set_opt(
			RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE, &evi);
		if(errval != RSB_ERR_NO_ERROR)
		{
		       /*! [Copy error message to string] */
			char errbuf[256];
			rsb_strerror_r(errval,&errbuf[0],sizeof(errbuf));
			printf("Failed setting the"
			" RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE"
			" library option (reason string:\n%s).\n",errbuf);
		       /*! [Copy error message to string] */
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
		/*! [Setting a single optional library parameter] */
	}

	/*! [Allocate matrix with error flags check] */
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
	/*! [Allocate matrix with error flags check] */
	printf("Correctly allocated a matrix.\n");
	printf("Summary information of the matrix:\n");
	/* print out the matrix summary information  */

	/*! [Get an info string for the matrix] */
	rsb_mtx_get_info_str(mtxAp,"RSB_MIF_MATRIX_INFO__TO__CHAR_P",
			ib,sizeof(ib));
	printf("%s",ib);
	/*! [Get an info string for the matrix] */
	printf("\n");

	/*! [Multiply a sparse matrix by a dense vector] */
	if((errval = 
		rsb_spmv(RSB_TRANSPOSITION_N,&one,mtxAp,B,1,&one,X,1))
			!= RSB_ERR_NO_ERROR )
	{
		printf("Error performing a multiplication!\n");
		goto err;
	}
	/*! [Multiply a sparse matrix by a dense vector] */

	printf("Correctly performed a SPMV.\n");

	/*! [Free a sparse matrix] */
	rsb_mtx_free(mtxAp);
	/*! [Free a sparse matrix] */

	printf("Correctly freed the matrix.\n");

        /*! [Finalize the library] */
	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!= RSB_ERR_NO_ERROR)
	{
		printf("Error finalizing the library!\n");
		goto err;
	}
        /*! [Finalize the library] */

	printf("Correctly finalized the library.\n");
	printf("Program terminating with no error.\n");
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	printf("Program terminating with error.\n");
	return EXIT_FAILURE;
}

static int vec_load_snip(const int argc, char * const argv[])
{
	struct rsb_mtx_t *mtxAp = NULL;
	const rsb_blk_idx_t brA = RSB_DEFAULT_BLOCKING,
                            bcA = RSB_DEFAULT_BLOCKING;
	rsb_nnz_idx_t nnzA = 4;
	rsb_coo_idx_t  nrA = 3;
	rsb_coo_idx_t  ncA = 3;
	const rsb_coo_idx_t    IA[] = { 0, 1, 2, 0 };
	const rsb_coo_idx_t    JA[] = { 0, 1, 2, 2 };
	const RSB_DEFAULT_TYPE VA[] = { 11, 22, 33, 13 };
	RSB_DEFAULT_TYPE XV[] = { 0,0,0,0,0,0 };
	rsb_coo_idx_t  vl = 0;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	/* library initialization */
	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS)!=RSB_ERR_NO_ERROR)
	{
		return EXIT_FAILURE;
	}

	/* allocation */

	/*! [Allocate matrix without error flags check] */
	mtxAp = rsb_mtx_alloc_from_coo_const(
			VA,IA,JA,nnzA,typecode,nrA,ncA,
			brA,bcA,RSB_FLAG_NOFLAGS,NULL);
	if(!mtxAp)
	{
		return EXIT_FAILURE;
	}
	/*! [Allocate matrix without error flags check] */

	/* printout */
	if(RSB_ERR_NO_ERROR!=(errval = rsb_file_mtx_save(mtxAp,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err;
	}
	
	/* matrix transposition */
	/*! [Clone and transpose a sparse matrix] */
	if( RSB_ERR_NO_ERROR != (errval =
		rsb_mtx_clone(&mtxAp,RSB_NUMERICAL_TYPE_SAME_TYPE,
		RSB_TRANSPOSITION_T,NULL,mtxAp,RSB_FLAG_IDENTICAL_FLAGS)))
	{
		goto err;
	}
	/*! [Clone and transpose a sparse matrix] */

	/* printout */
	if(RSB_ERR_NO_ERROR!=(errval = rsb_file_mtx_save(mtxAp,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err;
	}

	rsb_mtx_free(mtxAp);

	/* doing the same after load from file */

        /*! [Load a matrix from Matrix Market file] */
	mtxAp = rsb_file_mtx_load("pd.mtx",
		RSB_FLAG_NOFLAGS,typecode,NULL);
	if(!mtxAp)
	{
		return EXIT_FAILURE;
	}
        /*! [Load a matrix from Matrix Market file] */

	
        /*! [Print a matrix to standard output] */
	if(RSB_ERR_NO_ERROR!=(errval = rsb_file_mtx_save(mtxAp,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err;
	}
        /*! [Print a matrix to standard output] */

	/* one can see dimensions in advance, also */
	/*! [Get dimensions of sparse matrix stored in Matrix Market file] */
	if(RSB_ERR_NO_ERROR!=(errval =
		rsb_file_mtx_get_dims("pd.mtx",&nrA,&ncA,&nnzA,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err; /* may have not configured what needed */
	}
	/*! [Get dimensions of sparse matrix stored in Matrix Market file] */

	/* A matrix can be rendered to Postscript. */
	{
	/*! [Render a Sparse matrix to Postscript] */
		if(RSB_ERR_NO_ERROR!=(errval =
		rsb_mtx_rndr("pd.eps",mtxAp,512,512,RSB_MARF_EPS_B)))
			goto err;
	/*! [Render a Sparse matrix to Postscript] */
	}

	rsb_mtx_free(mtxAp);

	/*! [Load vector matrix from file] */
	/* also vectors can be loaded */
	if(RSB_ERR_NO_ERROR!=(errval = 
		rsb_file_vec_load("vf.mtx",typecode,NULL,&vl )))
		goto err;
	/* we expect vf.mtx to be 6 rows long */
	if( vl != 6 )
	{
		goto err;
	}

	if(RSB_ERR_NO_ERROR!=(errval = 
		rsb_file_vec_load("vf.mtx",typecode,XV, NULL )))
		goto err;
	/*! [Load vector matrix from file] */

	/*! [Render matrix from Matrix Market pixelmap in memory] */
	/* matrices can be rendered from file to a pixelmap as well */
	{
		unsigned char pixmap[3*2*2];

		if(RSB_ERR_NO_ERROR!=(errval =
		rsb_file_mtx_rndr(pixmap,"pd.mtx",2,2,2,RSB_MARF_RGB)))
			goto err;
	}
	/*! [Render matrix from Matrix Market pixelmap in memory] */

	if(RSB_ERR_NO_ERROR != rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
	{
		goto err;
	}
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	return EXIT_FAILURE;
}

static int tune_snip__tune_from_file(char * const filename,
                rsb_int_t wvat)
{
	struct rsb_mtx_t *mtxMp = NULL;
	/* spmv specific variables */
	const RSB_DEFAULT_TYPE alpha = 1;
	const RSB_DEFAULT_TYPE beta = 1;
       	rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
      	const rsb_coo_idx_t nrhs = 2;  /* number of right hand sides */
       	rsb_trans_t transA = RSB_TRANSPOSITION_N; /* transposition */
       	rsb_nnz_idx_t ldB = 0;
       	rsb_nnz_idx_t ldC = 0;
	/* misc variables */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t dt;
	char ib[200];
	const char*is = "RSB_MIF_MATRIX_INFO__TO__CHAR_P";
	/* misc variables */
	/* input autotuning variables */
       	rsb_int_t oitmax = 1 /*15*/;	/* auto-tune iterations */
       	rsb_time_t tmax = 0.1;	/* time per autotune operation */
	/* output autotuning variables */
	rsb_flags_t flagsA = RSB_FLAG_NOFLAGS;
	/* int ione = 1; */
	rsb_type_t typecodea [] = RSB_MATRIX_TYPE_CODES_ARRAY;
	int typecodei;

	errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS);

	if( (errval) != RSB_ERR_NO_ERROR )
		goto err;

	errval = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING, &wvat );
	
	/*
	errval = rsb_lib_set_opt(RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE, &ione);
	*/

	if( (errval) != RSB_ERR_NO_ERROR )
		goto err;

	printf("Loading matrix from file \"%s\".\n",filename);

	mtxMp = rsb_file_mtx_load(filename, flagsA, typecodea[0], &errval);

	if( (errval) != RSB_ERR_NO_ERROR )
		goto err;

	for( typecodei = 0 ; typecodei < RSB_IMPLEMENTED_TYPES; ++typecodei )
	{
		rsb_type_t typecode = typecodea[typecodei];
		struct rsb_mtx_t *mtxAp = NULL;
		struct rsb_mtx_t *mtxOp = NULL;
		rsb_real_t sf = 0.0;
       		rsb_int_t tn = 0;

		sf = 0.0;
       		tn = 0;

		printf("Considering %c clone.\n",typecode);
		
		errval = rsb_mtx_clone(&mtxAp, typecode, transA, NULL, mtxMp,
			       	flagsA);

		if( (errval) != RSB_ERR_NO_ERROR )
			goto err;

		printf("Base matrix:\n");
		rsb_mtx_get_info_str(mtxAp,is,ib,sizeof(ib));
		printf("%s\n\n",ib);

		dt = -rsb_time();
		errval = rsb_tune_spmm(NULL, &sf, &tn, oitmax, tmax, transA,
		     &alpha, mtxAp, nrhs, order, NULL, ldB, &beta, NULL, ldC);

		dt += rsb_time();
		if(tn == 0)
		printf("After %lfs, autotuning routine did not find a better"
			" threads count configuration.\n",dt);
		else
		printf("After %lfs, thread autotuning declared speedup of %lg x,"
			" when using threads count of %d.\n",dt,sf,tn);
		printf("\n");


		dt = -rsb_time();

		mtxOp = mtxAp;
		errval = rsb_tune_spmm(&mtxAp, &sf, &tn, oitmax, tmax, transA,
		       	&alpha, NULL, nrhs, order, NULL, ldB, &beta, NULL, ldC);
		if( (errval) != RSB_ERR_NO_ERROR )
			goto err;

		dt += rsb_time();
		if( mtxOp == mtxAp )
		{
			printf("After %lfs, global autotuning found old matrix optimal,"
			" with declared speedup %lg x when using %d threads\n",dt,sf,tn);
		}
		else
		{
			printf("After %lfs, global autotuning declared speedup of %lg x,"
			" when using threads count of %d and a new matrix:\n",dt,sf,tn);
			rsb_mtx_get_info_str(mtxAp,is,ib,sizeof(ib));
			printf("%s\n",ib);
		}
		printf("\n");

		/* user is expected to:
		errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn);
		and use mtxAp in SpMV.
 		*/
		rsb_mtx_free(mtxAp);
		mtxAp = NULL;
	}
	rsb_mtx_free(mtxMp);
	mtxMp = NULL;

	goto ret;
ret:
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	printf("Program terminating with error.\n");
	return EXIT_FAILURE;
}

int tune_snip__main(const int argc, char * const argv[])
{
	/*!
	 Autotuning example.
	 */
	/* matrix variables */
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */
	const int bs = RSB_DEFAULT_BLOCKING;
	rsb_coo_idx_t nrA = 5; /* number of rows */
	rsb_coo_idx_t ncA = 5; /* number of cols */
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_coo_idx_t rd = 1;/* every rd rows one is non empty */
	const rsb_coo_idx_t cd = 4;/* every cd cols one is non empty */
	rsb_nnz_idx_t nnzA = (nrA/rd)*(ncA/cd); /* nonzeroes */
	rsb_coo_idx_t*IA = NULL;
	rsb_coo_idx_t*JA = NULL;
	RSB_DEFAULT_TYPE*VA = NULL;
	/* spmv specific variables */
	const RSB_DEFAULT_TYPE alpha = 1;
	const RSB_DEFAULT_TYPE beta = 1;
	RSB_DEFAULT_TYPE*Cp = NULL;
	RSB_DEFAULT_TYPE*Bp = NULL;
       	rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
      	const rsb_coo_idx_t nrhs = 2;  /* number of right hand sides */
       	const rsb_trans_t transA = RSB_TRANSPOSITION_N; 
       	rsb_nnz_idx_t ldB = nrA;
       	rsb_nnz_idx_t ldC = ncA;
	/* misc variables */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const size_t so = sizeof(RSB_DEFAULT_TYPE);
	const size_t si = sizeof(rsb_coo_idx_t);
	rsb_time_t dt,odt;
       	rsb_int_t t;
       	const rsb_int_t tt = 100;	/* will repeat spmv tt times */
	char ib[200];
	const char*is = "RSB_MIF_MATRIX_INFO__TO__CHAR_P";
	/* misc counters */
       	rsb_coo_idx_t ci; 
	rsb_coo_idx_t ri;
	rsb_coo_idx_t ni;
	rsb_int_t nrhsi;
	/* misc variables */
	rsb_time_t etime = 0.0;
	/* input autotuning variables */
       	const rsb_int_t oitmax = 15;	/* auto-tune iterations */
       	const rsb_time_t tmax = 0.1;	/* time per autotune operation */
	/* input/output autotuning variables */
       	rsb_int_t tn = 0;	/* threads number */
	/* output autotuning variables */
	rsb_real_t sf = 0.0;	/* speedup factor obtained from auto tuning */
	const rsb_int_t wvat = 1; /* want verbose autotuning; see 
                  documentation of RSB_IO_WANT_VERBOSE_TUNING */

	if(argc > 1 && !isdigit(argv[1][0]) )
		return tune_snip__tune_from_file(argv[1],wvat);

	if(argc > 1)
	{
		nrA = ncA = atoi(argv[1]);
		if ( nrA < RSB_MIN_MATRIX_DIM || (nrA > (RSB_MAX_MATRIX_DIM) ))
			goto err;

		nnzA = (nrA/rd)*(ncA/cd);
       		ldB = nrA;
       		ldC = ncA;
	}

	printf("Creating %d x %d matrix with %d nonzeroes.\n",(int)nrA,
		(int)ncA, (int)nnzA);

	IA = calloc(nnzA, si);
	JA = calloc(nnzA, si);
	VA = calloc(nnzA, so);
	Bp = calloc(nrhs*ncA ,so);
	Cp = calloc(nrhs*nrA ,so);

	if( ! ( VA && IA && JA && Bp && Cp ) )
		goto err;

	for(nrhsi=0;nrhsi<nrhs;++nrhsi)
		for(ci=0;ci<ncA/cd;++ci)
			Bp[nrhsi*ldC+ci] = 1.0;

	for(nrhsi=0;nrhsi<nrhs;++nrhsi)
		for(ri=0;ri<nrA/rd;++ri)
			Cp[nrhsi*ldC+ri] = 1.0;

	ni = 0;

	for(ci=0;ci<ncA/cd;++ci)
		for(ri=0;ri<nrA/rd;++ri)
		{
			VA[ni] = nrA * ri + ci,
			IA[ni] = ri;
			JA[ni] = ci;
			printf("%d/%d  %d %d -> %d\n",(int)ni,(int)nnzA,
				(int)IA[ni],(int)JA[ni],(int)VA[ni]);
			ni++;
		}

	printf("Done.\n");

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))
			!= RSB_ERR_NO_ERROR) goto err;

	errval = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING, &wvat );

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,bs,bs,
		RSB_FLAG_NOFLAGS,&errval);

	/* VA, IA, JA are not necessary anymore */
	free(VA);
	free(IA);
	free(JA);
	VA = NULL;
       	IA = NULL;
       	JA = NULL;

	if((!mtxAp) || (errval != RSB_ERR_NO_ERROR))
		goto err;

	printf("Allocated matrix of %zd nonzeroes:\n",(size_t)nnzA);
	rsb_mtx_get_info_str(mtxAp,is,ib,sizeof(ib));
	printf("%s\n\n",ib);

	dt = - rsb_time();
	for(t=0;t<tt;++t)
		/* 
		   If nrhs == 1, the following is equivalent to
		   rsb_spmv(transA,&alpha,mtxAp,Bp,1,&beta,Cp,1);
		*/
		rsb_spmm(transA,&alpha,mtxAp,nrhs,order,Bp,ldB,&beta,Cp,ldC);
	dt += rsb_time();
	odt = dt;
	printf("Before auto-tuning, %d multiplications took %lfs.\n",tt,dt);

	printf("Threads autotuning (may take more than %lfs)...\n",
			oitmax*tmax);
	dt = -rsb_time();
	errval = rsb_tune_spmm(NULL, &sf, &tn, oitmax, tmax, transA,
		       	&alpha, mtxAp, nrhs, order, Bp, ldB, &beta, Cp, ldC);
	dt += rsb_time();
	if(errval != RSB_ERR_NO_ERROR)
		goto err;

	if(tn == 0)
	printf("After %lfs, autotuning routine did not find a better"
			" threads count configuration.\n",dt);
	else
	printf("After %lfs, autotuning routine declared speedup of %lg x,"
			" when using threads count of %d.\n",dt,sf,tn);

	errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn);
	if(errval != RSB_ERR_NO_ERROR)
		goto err;

	rsb_mtx_get_info_str(mtxAp,is,ib,sizeof(ib));
	printf("%s\n",ib);

	dt = -rsb_time();
	for(t=0;t<tt;++t)
		/*rsb_spmv(transA,&alpha,mtxAp,Bp,1,&beta,Cp,1);*/
		rsb_spmm(transA,&alpha,mtxAp,nrhs,order,Bp,ldB,&beta,Cp,ldC);
	dt += rsb_time();
	printf("After threads auto-tuning, %d multiplications took %lfs"
			"  --  effective speedup of %lg x\n",tt,dt,odt/dt);
	odt = dt;


	tn = 0; /* this will restore default threads count */
	errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn);
	if(errval != RSB_ERR_NO_ERROR)
		goto err;
	errval = rsb_lib_get_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn);
	if(errval != RSB_ERR_NO_ERROR)
		goto err;

	printf("Matrix autotuning (may take more than %lfs; using %d"
			" threads )...\n", oitmax*tmax, tn);

	/* A negative tn will request also threads autotuning: */
	/* tn = -tn; */

	dt = -rsb_time();
	errval = rsb_tune_spmm(&mtxAp, &sf, &tn, oitmax, tmax, transA,
		       	&alpha,  NULL, nrhs, order, Bp, ldB, &beta, Cp, ldC);
	dt += rsb_time();

	if(errval != RSB_ERR_NO_ERROR)
		goto err;

	if(tn == 0)
	printf("After %lfs, autotuning routine did not find a better"
			" threads count configuration.\n",dt);
	else
	printf("After %lfs, autotuning routine declared speedup of %lg x,"
			" when using threads count of %d.\n",dt,sf,tn);

	rsb_mtx_get_info_str(mtxAp,is,ib,sizeof(ib));
	printf("%s\n",ib);

	dt = -rsb_time();
	for(t=0;t<tt;++t)
		/*rsb_spmv(transA,&alpha,mtxAp,Bp,1,&beta,Cp,1);*/
		rsb_spmm(transA,&alpha,mtxAp,nrhs,order,Bp,ldB,&beta,Cp,ldC);
	dt += rsb_time();
	printf("After threads auto-tuning, %d multiplications took %lfs"
			"  --  further speedup of %lg x\n",tt,dt,odt/dt);

	get_coo_block_snippet(mtxAp,nnzA);

	rsb_mtx_free(mtxAp);
	free(Cp);
	free(Bp);


	errval = rsb_lib_get_opt(RSB_IO_WANT_LIBRSB_ETIME,&etime);
	if(errval == RSB_ERR_UNSUPPORTED_FEATURE)
	{
		printf("librsb timer-based profiling is not supported in "
		"this build. If you wish to have it, re-configure librsb "
	        "with its support. So you can safely ignore the error you"
		" might just have seen printed out on screen.\n");
		errval = RSB_ERR_NO_ERROR;
	}
	else
	if(etime) /* This will only work if enabled at configure time. */
		printf("Elapsed program time is %5.2lfs\n",etime);

	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!=RSB_ERR_NO_ERROR)
		goto err;
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	printf("Program terminating with error.\n");
	return EXIT_FAILURE;
}

int main_rsb_coo_cleanup1(const int argc, char * const argv[])
{
	/*! [COO cleanup 1] */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t nnzA = 4;
	const rsb_coo_idx_t  nrA = 4;
	const rsb_coo_idx_t  ncA = 4;
	rsb_coo_idx_t    IA[] = { 1, 1, 1, 2 };
	rsb_coo_idx_t    JA[] = { 1, 1, 3, 2 };
	RSB_DEFAULT_TYPE VA[] = { 1, 10, 13, 22 };
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	rsb_flags_t flagsA = RSB_FLAG_DUPLICATES_SUM | RSB_FLAG_SORTED_INPUT;
	// IA={1,1,1,2} JA={1,1,3,2} VA={1,10,13,22} nnzA=4 nrA=4 nca=4

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))
			!= RSB_ERR_NO_ERROR) goto err;

	errval = rsb_coo_cleanup(&nnzA, VA, IA, JA,
		 nnzA, nrA, ncA, typecode, flagsA );
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_coo_cleanup!\n");
		goto err;
	}
	// IA={1,1,2} JA={1,3,2} VA={11,13,22} nnzA=3 nrA=4 nca=4
	/*! [COO cleanup 1] */

	if(nnzA!=3)
	{
		printf("Unexpected nnz count out of rsb_coo_cleanup!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	if(VA[0]!=11)
	{
		printf("Unexpected VA out of rsb_coo_cleanup!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	if(IA[1]!=1 || JA[1]!=3)
	{
		printf("Unexpected IA/JA out of rsb_coo_cleanup!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	if((errval = rsb_lib_exit(RSB_NULL_INIT_OPTIONS))
			!= RSB_ERR_NO_ERROR) goto err;
err:
	return errval;
}

int main_rsb_lib_reinit(const int argc, char * const argv[])
{
	/*! [rsb_lib_reinit__rsb_lib_set_opt_str_snip] */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_initopts io;

	rsb_int_t ione={1};
	enum rsb_opt_t keys[]={RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE};
	void*values[]={&ione};
	io.action=RSB_IO_SPECIFIER_SET;
	io.keys=keys;
	io.values=values;
	io.n_pairs=1;

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))
			!= RSB_ERR_NO_ERROR) goto err;

	// won't print anything
	if((errval = rsb_lib_reinit(&io))
			!= RSB_ERR_NO_ERROR) goto err;

	// may print verbose message (depends on configure)
	if((errval = rsb_lib_reinit(NULL))
			!= RSB_ERR_NO_ERROR) goto err;

	// may print verbose message (depends on configure)
	if((errval = rsb_lib_set_opt_str(
		"RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE","0"))
			!= RSB_ERR_NO_ERROR) goto err;

	// won't print anything anymore
	if((errval = rsb_lib_exit(&io))
			!= RSB_ERR_NO_ERROR) goto err;
	/*! [rsb_lib_reinit__rsb_lib_set_opt_str_snip] */
err:
	return errval;
}

int main_rsb_coo_cleanup2(const int argc, char * const argv[])
{
	/*! [COO cleanup 2] */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t nnzA = 3;
	const rsb_coo_idx_t  nrA = 2;
	const rsb_coo_idx_t  ncA = 2;
	rsb_coo_idx_t    IA[] = { 1, 1, 1 };
	rsb_coo_idx_t    JA[] = { 2, 1, 1 };
	RSB_DEFAULT_TYPE VA[] = { 1, 2, 3 };
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_flags_t flagsA = RSB_FLAG_DUPLICATES_SUM 
			| RSB_FLAG_SORTED_INPUT
			| RSB_FLAG_FORTRAN_INDICES_INTERFACE;

        // IA={1,1,1} JA={2,1,1} VA={1,2,3} nnzA=3 nrA=2 nca=2
	errval =rsb_coo_sort(VA, IA, JA, nnzA, nrA, ncA,  typecode, flagsA);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_coo_cleanup!\n");
		goto err;
	}
        // IA={1,1,1} JA={1,1,2} VA={2,3,1} nnzA=3 nrA=2 nca=2

	errval = rsb_coo_cleanup(&nnzA, VA, IA, JA, nnzA,
		 nrA, ncA, typecode, flagsA );
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_coo_cleanup!\n");
		goto err;
	}
	// IA={1,1} JA={1,2} VA={5,1} nnzA=2 nrA=2 nca=2
	/*! [COO cleanup 2] */

	if(nnzA!=2)
	{
		printf("Unexpected nnz count out of rsb_coo_cleanup!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	if(VA[0]!=5)
	{
		printf("Unexpected VA out of rsb_coo_cleanup!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	if(IA[1]!=1 || JA[1]!=2)
	{
		printf("Unexpected IA/JA out of rsb_coo_cleanup!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
err:
	return errval;
}

static int main_rsb_mtx_alloc_from_csc_const(const int argc,
	 char * const argv[])
{
	/*! [snip__rsb_mtx_alloc_from_csc_const] */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;
	const rsb_blk_idx_t brA = RSB_DEFAULT_BLOCKING,
                            bcA = RSB_DEFAULT_BLOCKING;
	const rsb_nnz_idx_t nnzA = 4;
	const rsb_coo_idx_t  nrA = 3;
	const rsb_coo_idx_t  ncA = 3;
	const rsb_coo_idx_t    IA[] = { 0, 2, 1, 2 };
	const rsb_coo_idx_t    CP[] = { 0, 2, 3, 4 };
	const RSB_DEFAULT_TYPE VA[] = { 11, 31, 22, 33 };
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;

	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS)!=RSB_ERR_NO_ERROR)
	{
		return EXIT_FAILURE;
	}

	mtxAp = rsb_mtx_alloc_from_csc_const(
			VA,IA,CP,nnzA,typecode,nrA,ncA,
			brA,bcA,RSB_FLAG_NOFLAGS,NULL);
	if(!mtxAp)
	{
		return EXIT_FAILURE;
	}

	rsb_file_mtx_save(mtxAp, NULL);
	/*! [snip__rsb_mtx_alloc_from_csc_const] */

{
	/*! [snip__rsb_mtx_get_info] */
	rsb_real_t isopnnz;
	const enum rsb_mif_t miflags = 
		RSB_MIF_INDEX_STORAGE_IN_BYTES_PER_NNZ__TO__RSB_REAL_T;
	errval = rsb_mtx_get_info(mtxAp, miflags, &isopnnz);

	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_get_info!\n");
		goto err;
	}
	printf("RSB matrix uses %lf bytes per nnz.\n",(double)isopnnz);
	
	/*! [snip__rsb_mtx_get_info] */
}


{
	/*! [snip__rsb_mtx_upd_vals] */
        enum rsb_elopf_t elop_flags = RSB_ELOPF_NEG;
	const RSB_DEFAULT_TYPE omegap[] = {10};

	errval = rsb_mtx_upd_vals(mtxAp, elop_flags, NULL);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_upd_vals!\n");
		goto err;
	}

        elop_flags = RSB_ELOPF_MUL;
	errval = rsb_mtx_upd_vals(mtxAp, elop_flags, omegap);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_upd_vals!\n");
		goto err;
	}
	/*! [snip__rsb_mtx_upd_vals] */
}
{
	/*! [snip__rsb_mtx_get_vals] */
	const rsb_coo_idx_t    IA[] = { 2, 0, 2, 0 };
	const rsb_coo_idx_t    JA[] = { 2, 0, 0, 0 };
	RSB_DEFAULT_TYPE VA[] = { -1, -1, -1, -1 };

	errval = rsb_mtx_get_vals(mtxAp, VA, IA, JA, nnzA, RSB_FLAG_NOFLAGS);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_get_vals!\n");
		goto err;
	}
	/*! [snip__rsb_mtx_get_vals] */

	if( ! ( VA[0]==-330 && VA[1]==-110 && VA[3]==-110 && VA[2]==-310 ) )
	{
		printf("Unexpected rsb_mtx_get_vals output!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
}

{
	/*! [snip__rsb_mtx_get_rows_sparse] */
	rsb_coo_idx_t    IA[] = { 0, 0, 0, 0 };
	rsb_coo_idx_t    JA[] = { 0, 0, 0, 0 };
	RSB_DEFAULT_TYPE VA[] = { -1, -1, -1, -1 };
	rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t frA = 2, lrA = 2;
	rsb_nnz_idx_t rnz;
	RSB_DEFAULT_TYPE *alphap = NULL;

	errval = rsb_mtx_get_rows_sparse(transA, NULL, mtxAp, NULL,
		 NULL, NULL, frA, lrA, &rnz, RSB_FLAG_NOFLAGS);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_get_rows_sparse!\n");
		goto err;
	}

	printf("Rows between %d and %d have %d nnz\n",
		(int)frA,(int)lrA,(int)rnz);

	errval = rsb_mtx_get_rows_sparse(transA, alphap, mtxAp,
		 VA, IA, JA, frA, lrA, &rnz, RSB_FLAG_NOFLAGS);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_get_vals!\n");
		goto err;
	}
	/*! [snip__rsb_mtx_get_rows_sparse] */

	if( ! (rnz == 2) )
	{
		printf("Unexpected rsb_mtx_get_vals output!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
}

{
	/*! [snip__rsb_mtx_add_to_dense] */
	const rsb_nnz_idx_t ldB = 4, nrB = 3, ncB = 3;
	const rsb_bool_t rowmajorB = RSB_BOOL_TRUE;
	RSB_DEFAULT_TYPE Bp[ /*ldB*nrB*/ ] = {
		-1, -1, -1, -1,
		-1, -1, -1, -1,
		-1, -1, -1, -1
	};
	RSB_DEFAULT_TYPE *alphap = NULL;

	errval = rsb_mtx_add_to_dense(alphap, mtxAp, ldB,
		 nrB, ncB, rowmajorB, Bp);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_add_to_dense!\n");
		goto err;
	}

	/*! [snip__rsb_mtx_add_to_dense] */
	
	if( !  (
		Bp[ldB*0+0] == -111 &&
		Bp[ldB*1+1] == -221 &&
		Bp[ldB*2+0] == -311 &&
		Bp[ldB*2+2] == -331
		) 
	)
	{
		printf("Unexpected rsb_mtx_add_to_dense result!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
}

{
	/*! [snip__rsb_spmsp_to_dense] */
	const rsb_nnz_idx_t ldC = 4, nrC = 3, ncC = 3;
	const rsb_bool_t rowmajorC = RSB_BOOL_TRUE;
	RSB_DEFAULT_TYPE Cp[ /*ldC*nrC*/ ] = {
		0, 0, 0, -99,
		0, 0, 0, -99,
		0, 0, 0, -99
	};
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_trans_t transB = RSB_TRANSPOSITION_N;
	RSB_DEFAULT_TYPE *alphap = NULL;
	RSB_DEFAULT_TYPE *betap = NULL;

	errval = rsb_spmsp_to_dense(typecode, transA, alphap, mtxAp,
		 transB, betap, mtxAp , ldC, nrC, ncC, rowmajorC, Cp);
	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_spmsp_to_dense!\n");
		goto err;
	}
	/*! [snip__rsb_spmsp_to_dense] */
	/*
octave:1> A=[-110,0,0;0,-220,0;-310,0,-330]**2

		A =
		
		    12100        0        0
		        0    48400        0
		   136400        0   108900
	*/
	rsb_file_mtx_save(mtxAp, NULL);
	if( !  (
		Cp[ldC*0+0] == 12100 &&
		Cp[ldC*1+1] == 48400 &&
		Cp[ldC*2+0] == 136400 &&
		Cp[ldC*2+2] == 108900
		) 
	)
	{
		printf("Unexpected rsb_spmsp_to_dense result!\n");
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
}

{
	/*! [snip__rsb_spmsp] */
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_trans_t transB = RSB_TRANSPOSITION_N;
	RSB_DEFAULT_TYPE *alphap = NULL;
	RSB_DEFAULT_TYPE *betap = NULL;
	struct rsb_mtx_t * mtxCp = NULL;

	mtxCp = rsb_spmsp(typecode, transA, alphap, mtxAp,
		 transB, betap, mtxAp, &errval);
	if( !mtxCp )
	{
		printf("Error calling rsb_spmsp!\n");
		goto err;
	}	

	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_spmsp!\n");
		goto err;
	}
	/*! [snip__rsb_spmsp] */
	/*
octave:1> A=[-110,0,0;0,-220,0;-310,0,-330]**2

		A =
		
		    12100        0        0
		        0    48400        0
		   136400        0   108900
	*/
	rsb_file_mtx_save(mtxCp, NULL);
	rsb_mtx_free(mtxCp);
}

{
	/*! [snip__rsb_sppsp] */
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_trans_t transB = RSB_TRANSPOSITION_N;
	RSB_DEFAULT_TYPE *alphap = NULL;
	RSB_DEFAULT_TYPE *betap = NULL;
	struct rsb_mtx_t * mtxCp = NULL;

	mtxCp = rsb_sppsp(typecode, transA, alphap, mtxAp,
		 transB, betap, mtxAp, &errval);
	if( !mtxCp )
	{
		printf("Error calling rsb_sppsp!\n");
		goto err;
	}	

	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_sppsp!\n");
		goto err;
	}
	/*! [snip__rsb_sppsp] */
	rsb_file_mtx_save(mtxCp, NULL);
	rsb_mtx_free(mtxCp);
}

	rsb_mtx_free(mtxAp);
	if(RSB_ERR_NO_ERROR != rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
	{
		goto err;
	}
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	return EXIT_FAILURE;
}

static int main_rsb_strerror_r(const int argc, char * const argv[])
{
	/*! [snip__rsb_strerror_r] */
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;

	// ...
	
	if(errval != RSB_ERR_NO_ERROR)
	{
		char errbuf[256];

		rsb_strerror_r(errval,&errbuf[0],sizeof(errbuf));

		// error handling ...

		/*! [snip__rsb_strerror_r] */
		if(strstr(errbuf,"An error occurred")==NULL)
			return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

static int main_rsb_mtx_switch_to_csr(const int argc,
	 char * const argv[])
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;
	const rsb_blk_idx_t brA = RSB_DEFAULT_BLOCKING,
                            bcA = RSB_DEFAULT_BLOCKING;
	const rsb_nnz_idx_t nnzA = 4;
	const rsb_coo_idx_t  nrA = 3;
	const rsb_coo_idx_t  ncA = 3;
	rsb_coo_idx_t    IA[] = { 0, 2, 1, 2 };
	rsb_coo_idx_t    JA[] = { 0, 0, 1, 2 };
	RSB_DEFAULT_TYPE VA[] = { 11, 31, 22, 33 };
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;

	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS)!=RSB_ERR_NO_ERROR)
	{
		return EXIT_FAILURE;
	}

	mtxAp = rsb_mtx_alloc_from_coo_inplace(
			VA,IA,JA,nnzA,typecode,nrA,ncA,
			brA,bcA,RSB_FLAG_NOFLAGS,&errval);
	if(!mtxAp)
	{
		return EXIT_FAILURE;
	}

	rsb_file_mtx_save(mtxAp, NULL);
	
{
	/*! [snip__rsb_mtx_switch_to_csr] */
	rsb_coo_idx_t  *IA = NULL;
	rsb_coo_idx_t  *JA = NULL;
	RSB_DEFAULT_TYPE *VA = NULL;

	errval = rsb_mtx_switch_to_csr(mtxAp, (void**)&VA,
		&IA, &JA, RSB_FLAG_NOFLAGS);

	// NOTE: no rsb_mtx_free() necessary now..
	/*! [snip__rsb_mtx_switch_to_csr] */

	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_switch_to_csr!\n");
		goto err;
	}

	if( ! ( IA && JA && VA ) )
	{
		printf("Error calling rsb_mtx_switch_to_csr!\n");
		goto err;
	}

	if( IA[0] != 0 )
	{
		printf("Error using data from rsb_mtx_switch_to_csr!\n");
		goto err;
	}
	if( IA[nrA] != nnzA  )
	{
		printf("Error using data from rsb_mtx_switch_to_csr!\n");
		goto err;
	}
}

	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!= RSB_ERR_NO_ERROR)
	{
		goto err;
	}

	return EXIT_SUCCESS;
err:
	return EXIT_FAILURE;
}

static int main_rsb_mtx_switch_to_coo(const int argc,
	 char * const argv[])
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;
	const rsb_blk_idx_t brA = RSB_DEFAULT_BLOCKING,
                            bcA = RSB_DEFAULT_BLOCKING;
	const rsb_nnz_idx_t nnzA = 4;
	const rsb_coo_idx_t  nrA = 3;
	const rsb_coo_idx_t  ncA = 3;
	rsb_coo_idx_t    IA[] = { 0, 2, 1, 2 };
	rsb_coo_idx_t    JA[] = { 0, 0, 1, 2 };
	RSB_DEFAULT_TYPE VA[] = { 11, 31, 22, 33 };
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;

	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS)!=RSB_ERR_NO_ERROR)
	{
		return EXIT_FAILURE;
	}

	mtxAp = rsb_mtx_alloc_from_coo_inplace(
			VA,IA,JA,nnzA,typecode,nrA,ncA,
			brA,bcA,RSB_FLAG_NOFLAGS,&errval);
	if(!mtxAp)
	{
		return EXIT_FAILURE;
	}

	rsb_file_mtx_save(mtxAp, NULL);
	
{
	/*! [snip__rsb_mtx_switch_to_coo] */
	rsb_coo_idx_t  *RP = NULL;
	rsb_coo_idx_t  *JA = NULL;
	RSB_DEFAULT_TYPE *VA = NULL;

	errval = rsb_mtx_switch_to_coo(mtxAp, (void**)&VA,
		&RP, &JA, RSB_FLAG_NOFLAGS);

	// NOTE: no rsb_mtx_free() necessary now..
	/*! [snip__rsb_mtx_switch_to_coo] */

	if(errval != RSB_ERR_NO_ERROR )
	{
		printf("Error calling rsb_mtx_switch_to_coo!\n");
		goto err;
	}

	if( ! ( RP && JA && VA ) )
	{
		printf("Error calling rsb_mtx_switch_to_coo!\n");
		goto err;
	}

	if( RP[0] != 0 || RP[nnzA-1] != 2 || VA[nnzA-1] != 33  )
	{
		printf("Error using data from rsb_mtx_switch_to_coo!\n");
		goto err;
	}
}

	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!= RSB_ERR_NO_ERROR)
	{
		goto err;
	}

	return EXIT_SUCCESS;
err:
	return EXIT_FAILURE;
}

int main_rsb_psblas_trans_to_rsb_trans
 (const int argc, char * const argv[])
{
	if( rsb_psblas_trans_to_rsb_trans('N') 
		!= RSB_TRANSPOSITION_N )
		goto err;

	if( rsb_psblas_trans_to_rsb_trans('T') 
		!= RSB_TRANSPOSITION_T )
		goto err;

	if( rsb_psblas_trans_to_rsb_trans('C') 
		!= RSB_TRANSPOSITION_C )
		goto err;

	if( rsb_psblas_trans_to_rsb_trans('?') 
		!= -1 )
		goto err;

	return EXIT_SUCCESS;
err:
	return 1;
}

int main_rsb_mtx_get_prec(const int argc,
   char * const argv[])
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnzA = 7;		/* matrix nonzeroes count */
	const rsb_coo_idx_t nrA = 6;		/* matrix rows count */
	const rsb_coo_idx_t ncA = 6;		/* matrix columns count */
	/* nonzero row indices coordinates: */
	const rsb_coo_idx_t IA[] = {0,1,2,3,4,5,1};
	/* nonzero column indices coordinates: */
	const rsb_coo_idx_t JA[] = {0,1,2,3,4,5,5};
	const RSB_DEFAULT_TYPE VA[] = {11,22,33,44,55,66,16};
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/*! [snip__rsb_mtx_get_prec] */
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */
	struct rsb_mtx_t *mtxLUp [2];	/* matrix structure pointer */
	rsb_precf_t prec_flags = RSB_PRECF_ILU0;

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != 
			RSB_ERR_NO_ERROR)
	{
		printf("Error initializing the library!\n");
		goto err;
	}

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS /* force rsb */
		| RSB_FLAG_DUPLICATES_SUM /* sum dups */
		| RSB_FLAG_TRIANGULAR /* need triangle for spsv */
		, &errval);

	if((!mtxAp) || (errval != RSB_ERR_NO_ERROR))
	{
		printf("Error while allocating the matrix!\n");
		goto err;
	}

	errval = rsb_mtx_get_prec(mtxLUp,mtxAp, prec_flags, NULL);
	if( errval != RSB_ERR_NO_ERROR )
	{
		printf("Error while calling rsb_mtx_get_prec!\n");
		goto err;
	}
	// ...
	
	rsb_mtx_free(mtxLUp[0]);
	rsb_mtx_free(mtxLUp[1]);
	rsb_mtx_free(mtxAp );
	/*! [snip__rsb_mtx_get_prec] */


	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!= RSB_ERR_NO_ERROR)
	{
		goto err;
	}
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	printf("Program terminating with error.\n");
	return EXIT_FAILURE;
}

int main(const int argc, char * const argv[])
{
	/*!
	  A RSB-by-snips program.
	  */
	return EXIT_SUCCESS
		| hello_snip(argc,argv)
		| vec_load_snip(argc,argv)
		| tune_snip__main(argc,argv)
		| main_backsolve(argc,argv)
		| main_rsb_coo_cleanup1(argc,argv)
		| main_rsb_coo_cleanup2(argc,argv)
		| main_rsb_lib_reinit(argc,argv)
		| main_rsb_mtx_alloc_from_csc_const(argc,argv)
		| main_rsb_strerror_r(argc,argv)
		| main_rsb_mtx_switch_to_csr(argc,argv)
		| main_rsb_mtx_switch_to_coo(argc,argv)
		| main_rsb_psblas_trans_to_rsb_trans(argc,argv)
		| main_rsb_mtx_get_prec(argc,argv)
	;
}
