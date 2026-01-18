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
/* @cond INNERDOC  */
/*
 * unsorted random matrix market generator
 * */
/**
 @file
 @author Michele Martone
 @brief
 A toy program generating sparse matrices.
 */

//#include <stdlib.h>	/* bsearch, calloc, malloc */
//#include <stdio.h>	/* printf */
//#include "rsb.h"
#include "rsb_common.h"
#include "rsb_internals.h"

//#define RSB_CONDITIONAL_FREE(p) {if((p))free((p));(p)=NULL;}

int g_allow_duplicates;

static rsb_err_t gen_mm(rsb_coo_idx_t r, rsb_coo_idx_t c, rsb_nnz_idx_t nnz)
{
	/*
	 * Since we could modify our implementation of bitmap data structures over time, 
	 * there will always be a testing need, so we place here this routine.
	 * */
	rsb_coo_idx_t * IA=NULL,*JA=NULL;
	float *VA;
	rsb_nnz_idx_t i,k;
	rsb_nnz_idx_t duplicates = 0;
	
	/* We generate a r x c bitmap and two integer index arrays */
	IA = rsb__calloc(sizeof(rsb_coo_idx_t)*nnz);
	JA = rsb__calloc(sizeof(rsb_coo_idx_t)*nnz);
	VA = rsb__calloc(sizeof(float)*nnz);
	if(!VA || !IA || !JA)goto err;
	/* We populate the arrays with random coordinates (please note that there could be duplicates) */
	for(k=0;k<nnz;++k) {IA[k]=rand()%r;} 
	for(k=0;k<nnz;++k) {JA[k]=rand()%c;}
	for(k=0;k<nnz;++k) {VA[k]=(float)rand();} 
	/* with a very stupid algorithm for avoiding duplicates (COULD TAKE FOREVER!) */
	duplicates = 0;
	if(!g_allow_duplicates)
	do
	{
		/* WARNING : THERE IS NO PROOF THIS COULD TERMINATE, AS IT IS DEPENDENT ON THE PSEUDORANDOM NUMER GENERATOR */
		duplicates=0;
		for(k=0;k<nnz;++k) for(i=0;i<k;++i)if(IA[i]==IA[k]) if(JA[i]==JA[k])
		{
			IA[k]=rand()%r;
			JA[k]=rand()%c;
			++duplicates;
		}
	}
	while(duplicates);

	printf("%s","%%MatrixMarket matrix coordinate real general\n");
	printf("%zd %zd %zd\n",(rsb_printf_int_t)r,(rsb_printf_int_t)c,(rsb_printf_int_t)nnz);
	for(k=0;k<nnz;++k)
	{
		printf("%6zd %6zd %20g\n",(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),VA[k]);
	}

	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(JA);
	return RSB_ERR_NO_ERROR;
err:
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(JA);
	return RSB_ERR_GENERIC_ERROR;
}

rsb_option_t options[] = {
    {"nnz",     required_argument, NULL, 'n'},  
    {"cols",     required_argument, NULL, 'c'},  
    {"rows",     required_argument, NULL, 'r'},  
    {"banded",     required_argument, NULL, 'b'},  
    {"allow-duplicates",     no_argument, NULL, 'd'},  
    {"generate-matrix",     no_argument, NULL, 'g'},  
    {0,0,0,0}
};

//size_t strlen(char *s){char *ss=s;while(*s)++s;return s-ss;}

int rsb_genmm_main(int argc,char *argv[])
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	int c;
	int              opt_index = 0;
	rsb_coo_idx_t rows=0,cols=0;
	rsb_nnz_idx_t nnz=0;
	double want_percentage=0.0;
	rsb_coo_idx_t g_want_banded = 0;
	rsb_bool_t g_want_lower = 0;
	rsb_bool_t g_diagonal = 0;
	rsb_coo_idx_t g_want_spacing=1;
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT_INTEGER;

	g_allow_duplicates = 0;
	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS) != RSB_ERR_NO_ERROR)
	{
		RSB_ERROR("initialization error!\n");
		goto err;
	}
    	for (;;)
	{
		c = rsb__getopt_long(argc, argv, "gDb:dr:c:n:ls:", options, &opt_index);
		if (c == -1)break;

		switch (c)
		{
			case 'b':
			g_want_banded = rsb__util_atoi_km10(optarg);
			g_want_banded++;/* just a trick */
			break;
			case 'D':
			g_diagonal = 1;
			break;
			case 'd':
			g_allow_duplicates = 1;
			break;
			case 'c':
			cols = rsb__util_atoi_km10(optarg);
			break;
			case 'l':
			g_want_lower = 1;
			break;
			case 'n':
			nnz = rsb__util_atoi_km10(optarg);
			if(*optarg && optarg[strlen(optarg)-1]=='%')want_percentage=nnz;
			break;
			case 'r':
			rows = rsb__util_atoi_km10(optarg);
			break;
			case 'g':
			break;
			case 's':
			g_want_spacing = rsb__util_atoi_km10(optarg);
			break;
		}
	}

	if( g_want_lower )
	{
		struct rsb_mtx_t * mtxAp;

		if(rows==0)
			rows=cols;

		if( g_want_banded != 0 )
			errval = rsb__generate_blocked_banded_mtx(rows,g_want_spacing,g_want_banded-1,0,&mtxAp,typecode);
		else
			errval = rsb__generate_blocked_banded_mtx(rows,g_want_spacing,rows-1,0,&mtxAp,typecode);
		if(!mtxAp) goto err;
		rsb__do_file_mtx_save(mtxAp,NULL);
		RSB_MTX_FREE(mtxAp);
		return 0;
	}

	if( g_want_banded != 0 )
	{
		/* we want a banded matrix */
		struct rsb_mtx_t * mtxAp;
		mtxAp = rsb__generate_banded(1, 1 , rows, cols, /*rows/4*/ g_want_banded-1 , NULL, typecode);
		/* errval = rsb__generate_blocked_banded_mtx(rows,1,g_want_banded-1,g_want_banded-1,&mtxAp,typecode); */
		if(!mtxAp) goto err;
		rsb__do_file_mtx_save(mtxAp,NULL); /* FIXME: errval = rsb_pr..*/
		RSB_MTX_FREE(mtxAp);
		return 0;
	}

	if( want_percentage )
	{
		want_percentage *= 0.01;
		want_percentage *= rows;
		want_percentage *= cols;
		nnz = want_percentage ;
		if(!nnz)++nnz;
	}

	if( g_diagonal )
	{
		/* we want a diagonal matrix */
		struct rsb_mtx_t * mtxAp = NULL;
		/* mtxAp = rsb_generate_diagonal( rows, NULL, typecode); */
		errval = rsb__generate_blocked_banded_mtx(rows,1,0,0,&mtxAp,typecode);
		if(!mtxAp) goto err;
		rsb__do_file_mtx_save(mtxAp,NULL);
		RSB_MTX_FREE(mtxAp);
		return RSB_PROGRAM_SUCCESS;
	}

	if( nnz < 1 )
	{
		fprintf(stderr,
			"usage: %s -g -r rows -c cols \n"
			"\t [ -n nonzeros [%%] ] "
			"| [ -b bandwidth ] (-b for a banded matrix with 'bandwidth' wide bandwidth)\n"
			"\t[-d ] (-d means that duplicates are allowed) !\n",
			argv[0]
			);
		return RSB_PROGRAM_ERROR;
	}
	if( nnz > rows * cols )
	{
		fprintf(stderr,"can't generate more nonzeros than rows x columns!\n");
		return RSB_PROGRAM_ERROR;
	}
	errval = gen_mm(rows,cols,nnz);
	return RSB_ERR_TO_PROGRAM_ERROR(errval | rsb_lib_exit(NULL));
err:	
	fprintf(stderr,"some error occurred during matrix generation\n");
	/* no deallocation, though */
	return RSB_PROGRAM_ERROR;
}

/*
int main(int argc,char *argv[])
{
	return rsb_genmm_main(argc,argv);
}
*/

/* @endcond */
