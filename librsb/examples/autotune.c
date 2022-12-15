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

 @brief C "RSB autotune" example program based on <rsb.h>.
   uses #rsb_file_mtx_load(), rsb_spmm(), rsb_tune_spmm().

 \include autotune.c
*/
#include <rsb.h>	/* librsb header to include */
#include <stdio.h>	/* printf() */
#include <ctype.h>	/* isdigit() */
#include <stdlib.h>	/* atoi() */
/* #include "rsb_internals.h" */

static int tune_from_file(char * const filename, rsb_int_t wvat)
{
	struct rsb_mtx_t *mtxMp = NULL;
	/* spmv specific variables */
	const void * alphap = NULL; // equivalent to 1
	const void * betap = NULL; // equivalent to 1
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
		     alphap, mtxAp, nrhs, order, NULL, ldB, betap, NULL, ldC);

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
			alphap, NULL, nrhs, order, NULL, ldB, betap, NULL, ldC);
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

err:
	rsb_perror(NULL,errval);
	if ( errval != RSB_ERR_NO_ERROR )
		printf("Program terminating with error.\n");
	return errval;
}

int main(const int argc, char * const argv[])
{
	/*!
	 Autotuning example.
	 */
	/* matrix variables */
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */
	const int bs = RSB_DEFAULT_BLOCKING;
	rsb_coo_idx_t nrA = 500; /* number of rows */
	rsb_coo_idx_t ncA = 500; /* number of cols */
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
	{
		errval = tune_from_file(argv[1],wvat);
		goto ret;
	}

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
			ni++;
		}

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))
			!= RSB_ERR_NO_ERROR) goto err;

	errval = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING, &wvat );
	if( (errval) != RSB_ERR_NO_ERROR )
	{
		printf("Error setting option!\n");
		goto err;
	}

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
		   rsb_spmv(transA,alphap,mtxAp,Bp,1,betap,Cp,1);
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

ret:
	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
			!=RSB_ERR_NO_ERROR)
		goto err;
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	printf("Program terminating with error.\n");
	return EXIT_FAILURE;
}
