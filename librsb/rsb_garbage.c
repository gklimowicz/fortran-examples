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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains functions for benchmarking, integration testing, and miscellaneous.
 * 
 * FIXME : this file contains both important and obsolete code.
 * */

#include <ctype.h>	/*isdigit*/
#include <string.h>	/*memcmp, strchr*/
#include <math.h>	/*fabs*/
#include "rsb_garbage.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

#if 0
// FIXME: unused. 
static rsb_bool_t rsb_are_same_coo(
	void * VA,
	void * new_VA,
	rsb_coo_idx_t * IA, 
	rsb_coo_idx_t * new_IA, 
	rsb_coo_idx_t * JA, 
	rsb_coo_idx_t * new_JA, 
	rsb_nnz_idx_t nnz,
	size_t el_size,
	rsb_err_t * errvalp
	)
{
	/**
	 * \ingroup gr_internals

		FIXME
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!VA || !new_VA || !IA || !new_IA || !JA || !new_JA || RSB_INVALID_NNZ_INDEX(nnz))
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	if(RSB_MEMCMP(IA,new_IA,sizeof(rsb_coo_idx_t)*(nnz)))
		goto diff;

	if(RSB_MEMCMP(JA,new_JA,sizeof(rsb_coo_idx_t)*(nnz)))
		goto diff;

	if(RSB_MEMCMP(VA,new_VA,el_size*(nnz)))
		goto diff;

	return RSB_BOOL_TRUE;
diff:
	return RSB_BOOL_FALSE;
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return RSB_BOOL_MAYBE;
}
#endif

#if 0
// FIXME: unused. 
static rsb_bool_t rsb_util_is_matrix_equal_to_coo(
	struct rsb_mtx_t * mtxCp,
	void * VA,
	void * new_VA,
	rsb_coo_idx_t * IA, 
	rsb_coo_idx_t * new_IA, 
	rsb_coo_idx_t * JA, 
	rsb_coo_idx_t * new_JA, 
	rsb_nnz_idx_t nnz,
	rsb_coo_idx_t m,
	rsb_coo_idx_t k,
	size_t el_size,
	rsb_type_t typecode,
	rsb_flags_t flags,
	rsb_err_t * errvalp)
{
	/**
	 * \ingroup gr_internals
	 *
	 * FIXME
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_bool_t same = RSB_BOOL_FALSE;

	RSB_INFO("getting back matrix in coo format ... \n");
	if( (errval = rsb__do_get_coo_noalloc(mtxCp,new_VA,new_IA,new_JA,NULL,RSB_FLAG_C_INDICES_INTERFACE)) != RSB_ERR_NO_ERROR )
	{
		RSB_INFO("rsb__do_get_coo_noalloc returned with an error code\n");
		goto err;
	}

	RSB_INFO("sorting back reconstructed matrix in 1x1 blocks  ... \n");
	if((errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,m,k,typecode,flags))!=RSB_ERR_NO_ERROR)
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	if((errval = rsb__util_sort_row_major_inner(new_VA,new_IA,new_JA,nnz,m,k,typecode,flags))!=RSB_ERR_NO_ERROR)
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	same = rsb_are_same_coo( VA, new_VA, IA, new_IA, JA, new_JA, nnz, el_size, &errval);
	if(same == RSB_BOOL_MAYBE)
	{
		RSB_STDERR("error while comparing coo\n");
//		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
	if(same == RSB_BOOL_FALSE)
	{
		RSB_STDERR("mismatch\n");
		goto diff;
	}

	return 1;
diff:
	return 0;
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return -1;
}
#endif

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
/* Need a bitmap testing routine, but this is insufficient. */
int rsb__test_bitmap_driver(rsb_coo_idx_t r, rsb_coo_idx_t c)
{
	/**
	 * It will test for the library bitmap functionalities.
	 *
	 * \param r	should specify the bitmap rows count
	 * \param r	should specify the bitmap columns count
	 * 
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * FIXME : seems broken (./rsbench -Ot -b)
	 * */
	rsb_bitmap_data_t *bitmap    = NULL;
	rsb_coo_idx_t * IA=NULL,*JA = NULL;
	rsb_coo_idx_t i,k,nnz;
	rsb_nnz_idx_t counted_nnz = 0;
	rsb_nnz_idx_t duplicates  = 0;

	if(r<1 || c<1){RSB_ERROR("rsb__test_bitmap_driver(r=%d, c=%d)\n",r,c);goto err;}
	nnz=(r*c)/2;
	RSB_INFO("testing a %ld x %ld bitmap...\n",(long)r,(long)c);
	
	/* We generate a r x c bitmap and two integer index arrays */
	bitmap = rsb__allocate_bitmap(r,c);
	IA = rsb__malloc(sizeof(rsb_coo_idx_t)*nnz);
	JA = rsb__malloc(sizeof(rsb_coo_idx_t)*nnz);
	if(!bitmap || !IA || !JA)goto err;
	/* We populate the arrays with random coordinates (please note that there could be duplicates) */
	RSB_INFO("generating coefficients...\n");
	for(k=0;k<nnz;++k) {IA[k]=rsb__rand_coo_index(r);} 
	for(k=0;k<nnz;++k) {JA[k] =rsb__rand_coo_index(c);}
	/* with a very stupid algorithm for avoiding duplicates :) */
	RSB_INFO("fixing duplicates (WARNING : could take forever!)..  \n");
	duplicates = 0;
	do
	{
		/* WARNING : THERE IS NO PROOF THIS COULD TERMINATE, AS IT IS DEPENDENT ON THE PSEUDORANDOM NUMBER GENERATOR */
		duplicates=0;
		for(k=0;k<nnz;++k) for(i=0;i<k;++i)if(IA[i]==IA[k]) if(JA[i]==JA[k])
		{
			//RSB_STDOUT("k:%d i:%d ; ",k,i); RSB_STDOUT("IA[k]:%d JA[k]:%d - ",IA[k],JA[k]); RSB_STDOUT("IA[i]:%d JA[i]:%d\n",IA[i],JA[i]);
			IA[k] = rsb__rand_coo_index(r);
			JA[k]  = rsb__rand_coo_index(c);
			++duplicates;
		}
	}
	while(duplicates);
	
	/* We try to set bits in the bitmap according to the coordinate arrays */
	RSB_INFO("setting bits randomly...\n");
	for(k=0;k<nnz;++k) RSB_BITMAP_SET(bitmap,r,c,IA[k],JA[k]);

	/* We try to read bits in the bitmap according to the coordinate arrays and cross-validate */
	RSB_INFO("checking count directly...  ");

	counted_nnz=0;
	for(k=0;k<nnz;++k) if(RSB_BITMAP_GET(bitmap,r,c,IA[k],JA[k]))++counted_nnz;

#ifdef RSB_DEBUG_BITMAP
	if(counted_nnz!=nnz) {RSB_ERROR("inserted nonzeros : %d\ncounted nonzeros : %d\n",nnz,counted_nnz);goto err;}
#endif

	RSB_INFO(" ...ok\n");
	RSB_INFO("checking count indirectly...");
	counted_nnz=0;
	counted_nnz = rsb__bitmap_bit_count(bitmap,r,c);
	if(counted_nnz!=nnz-duplicates) {RSB_ERROR("inserted nonzeros : %d - %d = %d\ncounted nonzeros  : %d\n",nnz,duplicates,nnz-duplicates,counted_nnz);goto err;}
	RSB_INFO(" ...ok\n");

	/* We validate the bit count */
	/* ok, tests passed */
	rsb__free(bitmap);
	rsb__free(IA);
	rsb__free(JA);
	return 0;
err:
	RSB_CONDITIONAL_FREE(bitmap);
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	return -1;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_err_t rsb__test_gen_matrix(rsb_type_t typecode, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, int allow_duplicates)
{
	/**
	 * Generates a randomly populated matrix.
	 * FIXME : move to some macro
	 * FIXME : UNTESTED
	 */

	if( !IA || !JA || !VA)goto err;
	*IA=NULL;
	*JA=NULL;
	*VA=NULL;

	if( nnz<1 )goto err;
	if( rows < 1 || cols < 1 )goto err;

	if(RSB_SOME_ERROR( rsb__util_coo_alloc( VA, IA, JA,nnz,typecode,RSB_BOOL_TRUE)))
		goto err;

	if(rsb__test_fill_matrix_nnz(typecode, nnz, *VA ))goto err;
	if(rsb__test_fill_matrix_coords(*IA, *JA, rows, cols, nnz, allow_duplicates))goto err;

	return 0;
err:
	RSB_CONDITIONAL_FREE(*IA);
	RSB_CONDITIONAL_FREE(*JA);
	RSB_CONDITIONAL_FREE(*VA);
	return RSB_ERR_GENERIC_ERROR;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__test_fill_matrix_nnz(rsb_type_t typecode, rsb_nnz_idx_t nnz, void *VA )
{
	/**
	 * Fills a given array with random values.
	 *
	 * FIXME : move to some macro
	 * FIXME : UNTESTED
	 * */
	rsb_nnz_idx_t k;
	
	if(!VA || nnz<1 )goto err;

#ifdef RSB_NUMERICAL_TYPE_INT
	if(typecode == RSB_NUMERICAL_TYPE_INT ) for(k=0;k<nnz;++k) {((int *)(VA))[k]=(int )rand();} 
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
	if(typecode == RSB_NUMERICAL_TYPE_FLOAT ) for(k=0;k<nnz;++k) {((float *)(VA))[k]=(float )rand();} 
	else
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if(typecode == RSB_NUMERICAL_TYPE_DOUBLE) for(k=0;k<  nnz;++k) {((double*)(VA))[k]=(double)rand();} 
	else
#endif
	return -1;

	return 0;
err:
	return -1;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_coo_idx_t rsb__rand_coo_index(rsb_coo_idx_t max_plus_one)
{
	/**
	 * \ingroup gr_internals
	 */
	rsb_coo_idx_t i = (rsb_coo_idx_t )rand()%max_plus_one;
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(i));
	return i;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_blk_idx_t rsb__rand_blk_index(rsb_blk_idx_t max_plus_one)
{
	/**
	 * \ingroup gr_internals
	 */
	rsb_coo_idx_t i = (rsb_blk_idx_t )rand()%max_plus_one;
	RSB_DEBUG_ASSERT(RSB_IS_VALID_BLK_INDEX(i));
	return i;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__test_fill_matrix_coords(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t allow_duplicates)
{
	/**
	 * Since we could modify our implementation of bitmap data structures over time, 
	 * there will always be a testing need, so we place here this routine.
	 *
	 * FIXME : document this
	 * */
	rsb_nnz_idx_t i,k;
	rsb_nnz_idx_t duplicates = 0;
	
	if( !IA || !JA)goto err;
	if( rows < 1 || cols < 1 || nnz < 1 )goto err;

	/* We populate the arrays with random coordinates (please note that there could be duplicates) */
	for(k=0;k<nnz;++k) {IA[k]=rsb__rand_coo_index(rows);} 
	for(k=0;k<nnz;++k) {JA[k] =rsb__rand_coo_index(cols);}
	/* with a very stupid algorithm for avoiding duplicates (COULD TAKE FOREVER!) */
	duplicates = 0;
	if(!allow_duplicates)
	do
	{
		/* WARNING : THERE IS NO PROOF THIS COULD TERMINATE, AS IT IS DEPENDENT ON THE PSEUDORANDOM NUMER GENERATOR */
		duplicates=0;
		for(k=0;k<nnz;++k) for(i=0;i<k;++i)if(IA[i]==IA[k]) if(JA[i]==JA[k])
		{
			IA[k] = rsb__rand_coo_index(rows);
			JA[k]  = rsb__rand_coo_index(cols);
			++duplicates;
		}
	}
	while(duplicates);

	return 0;
err:
	return RSB_ERR_GENERIC_ERROR;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

int rsb__test_dump_main(const int argc,rsb_char_t *const argv[])
{
	/**
	 * \ingroup gr_internals
	 * This example main program reads in a Matrix Market file in fixed
	 * block format and dumps it out.
	 *
	 * TODO :
	 *         * this a strictly debugging function.
	 *         * it needs better documentation.
	 *         * only double precision numbers for now.
	 *         * it should go to some macro file
	 * 
	 * FIXME : this file functionality is already in test_matops.m4.
	 **/
	rsb_option_t options[] = {
	    {"block-rowsize",	required_argument, NULL, 'r'},  
	    {"block-columns",	required_argument, NULL, 'c'},  
	    {"type",		required_argument, NULL, 'T'},  
	    {"matrix-filename",	required_argument, NULL, 'f'},  
	    {"in-place-permutation",	no_argument, NULL, 'P'},  
	    {0,0,0,0}
	};

	int c;
	int opt_index = 0;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
#else
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
#endif

	const rsb_char_t * filename=NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_STORAGE_FLAGS;

    	for (;;)
	{
		c = rsb__getopt_long(argc, argv, RSB_SAMPLE_PROGRAM_OPTIONS_GET_FLAGS"T:f:", options, &opt_index);/* Flawfinder: ignore */
		if (c == -1)
			break;
		RSB_DO_FLAG_ADD(flags,rsb__sample_program_options_get_flags(c,optarg));
		switch (c)
		{
			case 'T':
			typecode = ( *optarg == ':' ? RSB_NUMERICAL_TYPE_DEFAULT : *optarg );
			break;
			case 'f':
			filename = optarg;
			break;
	    	}
	}

	if(!RSB_SOME_ERROR(errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)))
	if(filename)
	{
		struct rsb_mtx_t * mtxAp = rsb_file_mtx_load(filename, flags, typecode, &errval);
		errval = rsb__do_file_mtx_save( mtxAp, NULL );
		RSB_MASK_OUT_SOME_ERRORS(errval);
		RSB_MTX_FREE(mtxAp);
	}
	RSB_MASK_OUT_SOME_ERRORS(errval);
	RSB_DO_ERROR_CUMULATE(errval,rsb_lib_exit(RSB_NULL_EXIT_OPTIONS));
	return RSB_ERR_TO_PROGRAM_ERROR(errval);
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__test_gen_and_print_matrix(rsb_type_t typecode, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz)
{
	/**
	 * \ingroup gr_internals
	 * Allocates a whole matrix and returns it back to the caller.
	 * */
	int allow_duplicates = 0;
	if( !IA || !JA || !VA)goto err;
	*IA=NULL;
	*JA=NULL;
	*VA=NULL;
	if(nnz<1 )goto err;
	if( rows < 1 || cols < 1 )goto err;
	
	if(rsb__test_gen_matrix(typecode, IA, JA, VA, rows, cols, nnz, allow_duplicates ))goto err;
	if(rsb__test_print_coo_mm(typecode,RSB_FLAG_NOFLAGS,*IA,*JA,*VA,rows,cols,nnz,RSB_BOOL_TRUE,RSB_DEFAULT_STREAM ))goto err;

	/* Note : we do not free these arrays, but give them back to the caller */
	return 0;
err:
	RSB_CONDITIONAL_FREE(*IA);
	RSB_CONDITIONAL_FREE(*JA);
	RSB_CONDITIONAL_FREE(*VA);
	return -1;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_flags_t rsb__sample_program_options_get_flags(int c, const rsb_char_t * optarg)
{
	/**
		\ingroup gr_internals
		
		Function for demo programs letting user set sparse matrix format flags.

	 	c = rsb__getopt_long(argc, argv, ...

		\see RSB_SAMPLE_PROGRAM_OPTIONS_GET_FLAGS
	*/
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

	switch (c)
	{
			case 0x41: /* A */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_AUTO_BLOCKING);
			break;
#if 0
			case 0x45: /* E */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_DOUBLE_DETECTED_CACHE);/* FIXME : EXPERIMENTAL */;
			break;
			case 0x43: /* C */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_HALF_DETECTED_CACHE);/* FIXME : EXPERIMENTAL */;
			break;
#endif
			case 0x44: /* D */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG);/* FIXME : EXPERIMENTAL */;
			break;
			case 0x48: /* H */
			//RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES_COO);/* FIXME : EXPERIMENTAL */;
			//RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES_CSR);/* FIXME : EXPERIMENTAL */;
			RSB_ERROR("-H switch is DEPRECATED\n");
			break;
#ifdef RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES
			case 0x4C: /* L */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES);/* FIXME : EXPERIMENTAL */;
			break;
#endif
			case 0x52: /* R */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_QUAD_PARTITIONING /* FIXME : EXPERIMENTAL */ /*  | RSB_FLAG_RECURSIVE_SHRINK_BOUNDING_BOX*/);
			break;
			case 0x69: /* i */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR);
			break;
//			case 0x64: /* d */
//			RSB_DO_FLAG_ADD(flags,RSB_FLAG_SHOULD_DEBUG);/* new */
//			break;
			case 0x73: /* s */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORT_INPUT);
			break;
			case 'P': /* P */
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT);
			break;
//			case 'V': /* V */
//			RSB_STDERR("RSB_FLAG_SHOULD_BE_VERBOSE flag is not supported anymore\n");
//			//RSB_DO_FLAG_ADD(flags,RSB_FLAG_SHOULD_BE_VERBOSE);
//			break;
			case 'q': /* q */
			{
				const rsb_char_t * op=optarg;

				while(op && *op)
				{
					switch(toupper(*op))
					{
//						case 'T':
//						RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
//						break;
						case 'U':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
						break;
//						case 'L':
//						RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
//						break;
#if 0
						case 'C':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_HALF_DETECTED_CACHE);
						break;
						case 'E':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_DOUBLE_DETECTED_CACHE);
						break;
#endif
//						case 'E':
//						RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG);
//						break;
						case 'H':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES);
						break;
//						case 'H':
//						RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
//						break;
						case 'O':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES_COO);
						break;
						case 'R':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_QUAD_PARTITIONING);
						break;
						case 'T':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS);
						break;
#ifdef RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES
						case 'L':
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES);
						break;
#endif
						case 'Z': /* Z */
						RSB_WARN("Warning: Using experimental Z sorting flag.\n");
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING); // experimental only
						break;
					}
					++op;
				}
			}
			break;
			case 'F': /* F */
			/* FIXME : UNFINISHED, UNDOCUMENTED */
			{int _sf=0;
			if(strchr(optarg,0x6F))/* o */
				{RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_COO_STORAGE);_sf++;}
#ifdef RSB_FLAG_WANT_LINKED_STORAGE
			if(strchr(optarg,0x6C))/* l */ /* FIXME */
				{RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_LINKED_STORAGE);_sf++;}
#endif
			if(strchr(optarg,0x62))/* b */
				{RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);_sf++;}
			if(strchr(optarg,0x63))/* c */
				{RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER);_sf++;}
#ifdef RSB_FLAG_WANT_AUTO_MAJOR_ORDER
			if(strchr(optarg,'a'))/* c */
				{RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_AUTO_MAJOR_ORDER);_sf++;}
#endif
			if(strchr(optarg,0x76))/* v */
				{RSB_DO_FLAG_DEL(flags,RSB_FLAG_WANT_BCSS_STORAGE);_sf++;}
			if( !_sf ){ RSB_STDERR("specified an unknown matrix format (should be [b|v|l|o][c])\n");goto err;}
			}
			break;
    	}
err:
	return flags;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__oski_estimate_bcsr_fillin_from_csr(const rsb_nnz_idx_t * pntr, const rsb_coo_idx_t * indx, const rsb_coo_idx_t m, const rsb_coo_idx_t k, const rsb_nnz_idx_t nnz, rsb_fillin_t * efillinmap);

rsb_err_t rsb__oski_estimate_bcsr_fill_from_coo(/*  const*/ rsb_coo_idx_t * IA, /*const*/ rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_fillin_t * efillinmap )
{
	/*
		FIXME : we waste resources here.
	*/

	rsb_nnz_idx_t * ptr=NULL, * indx=NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void * VA=NULL;

	VA = rsb__malloc_vector(nnz,typecode);	/* ! FIXME ! */

	if(!VA)
	{
		errval = RSB_ERR_ENOMEM;
		goto err;
	}
	
	errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,m,k,typecode,RSB_FLAG_NOFLAGS);	/* ! */
	if(RSB_SOME_ERROR(errval))
		goto err;
	if(RSB_SOME_ERROR(errval = rsb__allocate_csr_arrays_from_coo_sorted(NULL,IA,JA,nnz,m,k,RSB_NUMERICAL_TYPE_INVALID_TYPE,NULL,&indx,&ptr)))
		goto err;

	rsb__oski_estimate_bcsr_fillin_from_csr(ptr, indx, m, k, nnz, efillinmap );
err:
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(ptr);
	RSB_CONDITIONAL_FREE(indx);
	
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__oski_estimate_bcsr_fillin_from_csr(const rsb_nnz_idx_t * pntr, const rsb_coo_idx_t * indx, const rsb_coo_idx_t m, const rsb_coo_idx_t k, const rsb_nnz_idx_t nnz, rsb_fillin_t * efillinmap)
{
	/*
		\ingroup gr_internals
		FIXME : UNFINISHED
		A basic implementation of OSKI's fillin blocking estimation algorithm.
	*/

	const rsb_int rua[] = RSB_ROWS_UNROLL_ARRAY;
	const rsb_int cua[] = RSB_COLUMNS_UNROLL_ARRAY;
	rsb_int ci=0,ri=0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if (RSB_WANT_EXPERIMENTAL_FILLIN_ESTIMATOR==1)
	rsb_nnz_idx_t * block_count=NULL;

	block_count = rsb__malloc(k*RSB_COLUMNS_UNROLL_ARRAY_LENGTH);
	if(!block_count)
	{
		errval = RSB_ERR_ENOMEM;
		goto err;
	}
#endif
	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
			efillinmap[ri*RSB_COLUMNS_UNROLL_ARRAY_LENGTH+ci]=RSB_REAL_ZERO;

//	RSB_STDOUT("#experimenal fillin estimator\n");

	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		rsb_coo_idx_t i=0,mi=0,/*Mi=0,*/j=0,bi=0;
//		const rsb_int fraction=m<10000?1000:m/(100*rua[ri]);
		const rsb_int fraction=1000/rua[ri];	/* highest the constant, the ligher the computation */
		rsb_nnz_idx_t num_blocks[RSB_COLUMNS_UNROLL_ARRAY_LENGTH];
		rsb_coo_idx_t last_block_index[RSB_COLUMNS_UNROLL_ARRAY_LENGTH];
		rsb_nnz_idx_t nnz_visited = 0;
//		rsb_nnz_idx_t dr=((m/rua[ri])/fraction)*rua[ri];/* FIXME */
		rsb_nnz_idx_t dr=fraction*rua[ri];/* FIXME */
//		Mi = RSB_MIN((m)/fraction,m-1);	/* FIXME */

//		dr=m;

		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
			num_blocks[ci]=0;
	
		//for(i=mi;i+rua[ri]-1<=Mi && i+rua[ri]-1<m;i+=dr)	/* FIXME */
		for(i=mi;i+rua[ri]-1<m;i+=dr)	/* FIXME */
		{
			rsb_nnz_idx_t ja[RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE];

//			RSB_STDOUT("#%zd / %zd\n",(rsb_printf_int_t)i,(rsb_printf_int_t)m);

			for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
				last_block_index[ci]=RSB_MARKER_NNZ_VALUE;

			/* for each nonzero */
			for(bi=0;bi<rua[ri];++bi)
				ja[bi]=pntr[i+bi];

#if (RSB_WANT_EXPERIMENTAL_FILLIN_ESTIMATOR==1)
			/* UNFINISHED */
			for(bi=0;bi<rua[ri];++bi)
			for(j=pntr[i+bi];j!=pntr[i+bi+1];++j)
			for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
			{
				if(!block_count(k*ci+j/cua[ci]))
				{
					block_count(k*ci+j/cua[ci])++;
					nnz_visited++;
				}
			}
			RSB_BITMAP_CLEAR(bitmap,RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE,k);
#else
			for(j=0;j<k;++j)
			for(bi=0;bi<rua[ri];++bi)
			if(pntr[i+bi]!=pntr[i+bi+1])
			{
/*				RSB_STDOUT("#%zd , %zd ? %zd , %zd \n",
					(rsb_printf_int_t)(bi),(rsb_printf_int_t)j,
					(rsb_printf_int_t)indx[ja[bi]],(rsb_printf_int_t)ja[bi]
					);*/
				if(indx[ja[bi]]==j)
				{
//					RSB_STDOUT("#%zd , %zd\n",(rsb_printf_int_t)(i+bi),(rsb_printf_int_t)j);
					nnz_visited++;

			//		RSB_STDOUT("#%zd -> ",(rsb_printf_int_t)ja[bi]);
			//		if(ja[bi]+1<pntr[i+bi+1])
					ja[bi]++;

//					RSB_STDOUT("%zd\n",(rsb_printf_int_t)ja[bi]);

					for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
					{
						if(j/cua[ci] != last_block_index[ci])
						{
							last_block_index[ci]= j/cua[ci] ;
							num_blocks[ci]++;
						}
					}
				}
			}
#endif

//			RSB_STDOUT("at the end of %zd, %zd nnz\n",(rsb_printf_int_t)i,(rsb_printf_int_t)nnz_visited);
		}
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		{
			efillinmap[ri*RSB_COLUMNS_UNROLL_ARRAY_LENGTH+ci]=(((rsb_fillin_t)num_blocks[ci])*rua[ri]*cua[ci])/((rsb_fillin_t)nnz_visited);
//			RSB_STDOUT("#nnz_visited : %d, num_blocks : %d \n",nnz_visited,num_blocks[ci]);
//			RSB_STDOUT("#%d %d %f\n",rua[ri],cua[ci],efillinmap[ri][ci]);
		}
	}

	goto err;
//	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
//		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
//			RSB_STDOUT("#%d %d %f\n",rua[ri],cua[ci],efillinmap[ri*RSB_COLUMNS_UNROLL_ARRAY_LENGTH+ci]);
err:
#if (RSB_WANT_EXPERIMENTAL_FILLIN_ESTIMATOR==1)
	RSB_CONDITIONAL_FREE(block_count);
#endif
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__do_column_expand(rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t * kp, rsb_int factor)
{
	rsb_nnz_idx_t i;
	rsb_coo_idx_t k;

	/*!
	*/

	k=*kp;

	if(factor>0)
	{
		for(i=0;i<nnz;++i)
		{
			JA[i]*=factor;
		}
		*kp=*kp*(factor);
	}
	else
	{
		factor=-factor;
#if 1
		/* mirror */
		for(i=0;i<nnz;++i)
		{
			JA[i]=(k-(JA[i]+1))*factor;
		}
		*kp=*kp*factor;
#else
		/* this has the potential of introducing duplicates, so do not use it */
		for(i=0;i<nnz;++i)
		{
			JA[i]=JA[i]/factor;
		}
		*kp=*kp/factor;
#endif
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__do_print_some_vector_stats(const void * p, rsb_type_t typecode, rsb_nnz_idx_t m, rsb_nnz_idx_t inc)
{
#if RSB_ALLOW_STDOUT
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_aligned_t errnorm[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	/* this is debug information, very cheap to include */
	rsb__util_find_min(errnorm,p,typecode,m,inc);
	RSB_STDOUT("#min:");
	RSB_DO_ERROR_CUMULATE(errval,rsb__debug_print_value(errnorm,typecode));
	RSB_STDOUT("\n");

	rsb__util_find_max(errnorm,p,typecode,m,inc);
	RSB_STDOUT("#max:");
	RSB_DO_ERROR_CUMULATE(errval,rsb__debug_print_value(errnorm,typecode));
	RSB_STDOUT("\n");

	rsb__util_vector_sum_strided(errnorm,p,typecode,m,inc);
	RSB_STDOUT("#sum:");
	RSB_DO_ERROR_CUMULATE(errval,rsb__debug_print_value(errnorm,typecode));
	RSB_STDOUT("\n");

	RSB_DO_ERROR_CUMULATE(errval,rsb__vector_norm_strided(errnorm,p,typecode,m,inc));
	RSB_STDOUT("#norm:");
	RSB_DO_ERROR_CUMULATE(errval,rsb__debug_print_value(errnorm,typecode));
	RSB_STDOUT("\n");

	RSB_DO_ERR_RETURN(errval)
#else
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif
}



#if 0
rsb_err_t _rsb_BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sU_de_uG(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup gr_kernels
	 * computes \f$y \leftarrow y + {A} \cdot x, where A \neq A^T. \f$
         * Matrix A should be blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 */

	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const rsb_coo_idx_t incx=1,incy=1;

	/*	Outer loop. Occurs on the major dimension.	*/
	for(i=0;RSB_LIKELY(i<Mdim);++i)
	{


#if 1
		register double c_0=0;				
		double *c=out+(1*(i*1));
		for(k=bpntr[i]+0,j=bindx[k];k<bpntr[i+1]-5  ;k+=6,a += rows*columns,j=bindx[k])	/* k is the index of the block */
		{
			const double *b = rhs+(1*(j*1));
			c_0+=a[(0*1)+0]*b[0];
			c_0+=a[(0*1)+1]**(rhs+(1*(bindx[k+1]*1)));
			c_0+=a[(0*1)+2]**(rhs+(1*(bindx[k+2]*1)));
			c_0+=a[(0*1)+3]**(rhs+(1*(bindx[k+3]*1)));
			c_0+=a[(0*1)+4]**(rhs+(1*(bindx[k+4]*1)));
			c_0+=a[(0*1)+5]**(rhs+(1*(bindx[k+5]*1)));
		}

		for(;k<bpntr[i+1]  ;++k,a += rows*columns,j=bindx[k])	/* k is the index of the block */
		{
			const double *b = rhs+(1*(j*1));
			c_0+=a[(0*1)+0]*b[0];
		}
		c[0]+=c_0;
#else
		for(k=bpntr[i]+0,j=bindx[k];k<bpntr[i+1]  ;k+=1,a += rows*columns,j=bindx[k])	/* k is the index of the block */
		{
			const double *b = rhs+(1*(j*1));
			double *c=out+(1*(i*1));
		register double c_0=0;				
		c_0+=a[(0*1)+0]*b[0];
			c[0]+=c_0;
		}
#endif

	}
	return 0;
}
#endif
/* @endcond */
