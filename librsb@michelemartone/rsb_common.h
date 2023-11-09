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
/* @cond INNERDOC */
/**
 * @author Michele Martone
 * @file
 * @brief Low level routines and tools for our sparse matrix formats implementations.
 */
#ifndef RSB_COMMON_H_INCLUDED
#define RSB_COMMON_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#ifdef __cplusplus
#define restrict	/* for now, the restrict keyword is not allowed in C++ */
#endif  /* __cplusplus */
/**
 *
 * VBR internals the user shouldn't use, no way.
 * This file contents are not meant to be used as an API (Application Programming Interface).
 * 
 * Manipulate this file at your own risk.
 *
 */
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* fileno ... */
#endif /* _XOPEN_SOURCE */
#include <stdlib.h>	/* bsearch, calloc, malloc */
#include <stdio.h>	/* printf */

#ifdef HAVE_CONFIG_H	/* hopefully, the only non-RSB_ prefixed symbol of ours */
#define RSB_HAVE_CONFIG_H HAVE_CONFIG_H
#endif /* HAVE_CONFIG_H */

#ifdef RSB_HAVE_CONFIG_H		/* autotools makefiles define this */
#include "rsb-config.h"		/* this should be the first include */
#endif /* RSB_HAVE_CONFIG_H */
#if RSB_WANT_OMP_RECURSIVE_KERNELS
#include <omp.h>
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#ifdef RSB_HAVE_UNISTD_H
#include <unistd.h>	/* getopt, gethostname (some system don't have it here!) */
#endif /* RSB_HAVE_UNISTD_H */
#ifdef RSB_HAVE_GETOPT_LONG
#ifdef RSB_HAVE_GETOPT_H
#include <getopt.h>	/* getopt_long is not always available (e.g.: on AIX, or any non GNU system) */
#endif /* RSB_HAVE_GETOPT_H */
           typedef struct option rsb_option_t;
           struct rsb_option_struct {
           	struct option ro;
           };
#else /* RSB_HAVE_GETOPT_LONG */
/*#undef required_argument*/
/*#undef no_argument*/
/*#undef optional_argument*/
#define required_argument	0 /* ONLY A STUB */
#define no_argument 		1 /* ONLY A STUB */
#define optional_argument	2 /* ONLY A STUB */
       extern char *optarg;
       extern int optind, opterr, optopt;

           struct rsb_option {
               const char *name;
               int         has_arg;
               int        *flag;
               int         val;
           };
           typedef struct rsb_option rsb_option_t;
           struct rsb_option_struct {
               struct rsb_option ro;
           };
/*!
 * \ingroup gr_internals
 * \brief An internal, helper structure.
 */
#endif /* RSB_HAVE_GETOPT_LONG */
#include <ctype.h>	/*isdigit*/

/*
 * Comment / uncomment the following to your likes
 */
#define RSB_WITH_MM
/*#undef RSB_WITH_MM*/

#ifndef __cplusplus
#ifdef RSB_WITH_MM
#include "rsb_mmio.h"
#endif /* RSB_WITH_MM */
#endif /* __cplusplus */


/* the NDEBUG and DEBUG symbols will affect lot of checking code */
/* these flags are NOT supported : they should be debugged :) */
/*#define NDEBUG 1*/
/*#undef  NDEBUG*/

/*
#ifdef NDEBUG
#define DEBUG 1
#endif
*/

/*!
 \internal
 For all the situations where 'int' would be used.
 */
typedef int rsb_int;

/* some ill situations could arise the need for this */
#if 0
#ifndef NULL
#define NULL ((VOID *)(0))
#endif
#endif
#define RSB_NUL '\0'

/*  #define RSB_INLINE inline	*/ /* experimental */
#define RSB_INLINE 	/* experimental */

/* if defined, less partitioning arrays allocations will occur for plain CSR and CSC matrices (EXPERIMENTAL) */
#define RSB_WANT_EXPERIMENTAL_IN_PLACE_RECURSIVE	 0	/** EXPERIMENTAL */
#define RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS 1	/** EXPERIMENTAL: should avoid extra (non BCSR-related) arrays*/
#define RSB_WANT_EXPERIMENTAL_FILLIN_ESTIMATOR 2		/* */

#define RSB_WANT_INDEX_BASED_SORT 1	/**< if defined, will enable a faster index+permutation based sorting (EXPERIMENTAL) */
#define RSB_WANT_Z_SORT 1		/**< if 1, enable Z sorting at all */
#define RSB_FORTRAN_VERBOSE_CALLS 0		/**< */
#define RSB_EXPERIMENTAL_WANT_PURE_BCSS 0	/**< EXPERIMENTAL : for BCSR, will prevent from allocating VBR arrays */
#define RSB_EXPERIMENTAL_USE_PURE_BCSS  1	/**< EXPERIMENTAL : for BCSR, will prevent from using VBR arrays   */
#define RSB_WANT_BCSC_LEAVES  0	/** */
#define RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR  1	/**< EXPERIMENTAL : for BCSR, will prevent from using VBR arrays   */

#define RSB_WANT_OMP_RECURSIVE_SPSV						1 	/**< EXPERIMENTAL  */
#define RSB_WANT_OMP_RECURSIVE_SPMV						1 	/**< EXPERIMENTAL  */
#define RSB_EXPERIMENTAL_FORCE_ROW_SUBDIVISIONS_UNTIL_CORES_NUMBER	   	0 	/**< EXPERIMENTAL  */
#define RSB_EXPERIMENTAL_ONE_SINGLE_LOCK_FOR_ALL_SUBMATRICES		   	(1&&RSB_WANT_OMP_RECURSIVE_KERNELS) 	/**< EXPERIMENTAL: incompatible with the prev.  */

#define RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT 1
#define RSB_ALLOW_EMPTY_MATRICES 1

#define RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS   	(1&&RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT ) 	/**< EXPERIMENTAL  */

#define RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS2   	(1&&RSB_WANT_OMP_RECURSIVE_KERNELS) 	/**< EXPERIMENTAL  */
/* mutually exclusive options for : RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS2   	 */
#define RSB_EXPERIMENTAL_SHOULD_TRAVERSE_WITHOUT_BLOCKING		   	1 	/**< EXPERIMENTAL  */
#define RSB_EXPERIMENTAL_ALTERNATING_SUBMATRIX_HEURISTIC			0	/**< EXPERIMENTAL  */
#define RSB_EXPERIMENTAL_WORK_BALANCING_HEURISTIC				1	/**< EXPERIMENTAL  */

#define RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_LINKED_LIST	0 	/**< EXPERIMENTAL  */
#define RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_TREE		(0) 	/**< The standard mechanism, 2-partitioned.  */
#define RSB_SHOULD_FAIL_INIT_IF_MEMHIER_DETECTION_FAILS		1 	/**<  */

#define RSB_EXPERIMENTAL_NO_SUBDIVIDE_ON_MIN_NNZ_PER_ROW_OR_COLUMN		1 	/**< Block matrix subdivision under a threshold */
#define RSB_EXPERIMENTAL_ROWS_SUBDIVIDE_TO_CORES_NUM				1 	/**< Subdivide to obtain no less matrices than cores */
#define RSB_CONST_MIN_NNZ_PER_ROW_OR_COLUMN_PER_SUBMATRIX				3	/**< Lower threshold for nnz/m or nnz/k, for subdivision, should be determined heuristically  */

#define RSB_EXPERIMENTAL_WANT_BEST_TIMES   1 	/**< EXPERIMENTAL  */
#define RSB_EXPERIMENTAL_UNLIMITED_RECURSION  		0 	/**< EXPERIMENTAL  */
#define RSB_EXPERIMENTAL_MORTON_ORDERED_RECURSION  	0 	/**< EXPERIMENTAL, UNFINISHED, just for demo purposes  */

#define RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE 		0	/**< EXPERIMENTAL */
#define RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_UNSYMMETRIC 	1	/**< EXPERIMENTAL, UNIMPLEMENTED */
#define RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY 			RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE

#define RSB_CONST_IMPOSSIBLY_BIG_TIME   		1000000000 	/**< in seconds, used when computing 'minimum' running times. any measured time interval should be less than RSB_CONST_IMPOSSIBLY_BIG_TIME.  */
#define RSB_MIN_ABOVE_INF(X,Y,MIN) RSB_MAX(RSB_MIN(X,Y),MIN)
#define RSB_CONST_TIME_FOR_MICRO_BENCHMARK   		0.1 	/**<  */
#define RSB_CONST_MIN_TIMES_FOR_MICRO_BENCHMARK   	10 	/**<  */

#define RSB_EXPERIMENTAL_EXPAND_SYMMETRIC_MATRICES_BY_DEFAULT	0	/**< EXPERIMENTAL */
#define RSB_WANT_BOUNDED_BOXES 			1
#define RSB_WANT_BOUNDED_BOXES_SPMV 			(RSB_WANT_BOUNDED_BOXES && 1)
#define RSB_WANT_BOUNDED_BOXES_SPSV 			(RSB_WANT_BOUNDED_BOXES && 1)
#define RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV	1
#define RSB_WANT_SM_TO_THREAD_MOD_MAPPING		1
#define RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPSV	1
#define RSB_WANT_RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING 0

#if defined(RSB_HAVE_GZFREAD) /* Binary input hack */
#define RSB_WANT_EXPERIMENTAL_BINARY_COO 1
#endif /* RSB_HAVE_GZFREAD */

#define RSB_WANT_BITMAP 0
	/** if RSB_WANT_BITMAP, should att to rsb_mtx_t:
	 * auxiliary data structures */
	/*struct rsb_options_t *options;*/	/* FIXME: deprecated -- will be deleted soon */

#define RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE	(/*RSB_OUT_ERR_VERBOSITY>0 &&*/ RSB_INT_ERR_VERBOSITY>0) 
#define RSB_WANT_PERFORMANCE_FILE	0
#if defined(RSB_HAVE_SIGNAL_H) /* && defined(RSB_HAVE_BITS_SIGACTION_H) */
#define RSB_WANT_ACTION 1
#else
#define RSB_WANT_ACTION 0
#endif
/* #define RSB_WANT_ACTION RSB_ALLOW_INTERNAL_GETENVS */ /* 1 */
#if RSB_WANT_ACTION
#define RSB_SHALL_QUIT ( rsb__quit_rsbench != 0 )
#define RSB_INTERNALS_RSBENCH_HEAD_DECLS extern int rsb__quit_rsbench;
void rsb__sigh(int signal);
#define RSB_SIGHR rsb__sigr();
#else /* RSB_WANT_ACTION */
#define RSB_SHALL_QUIT ( 0 != 0 )
#define RSB_INTERNALS_RSBENCH_HEAD_DECLS
#define RSB_SIGHR
#endif /* RSB_WANT_ACTION */

#include "rsb_struct.h"		/* */

/* @cond INNERDOC */
/*!
 * \internal
 * \ingroup gr_internals
 * \brief An internal, helper structure (OBSOLETE).
 *
 * A rsb_options_t structure could keep track of helper information
 * like :
 *
 * - desired nonzero pattern
 * - pointers to optimal or desired functions for certain operations
 * - ..
 * - .. data which is perfectly optional
 * */
struct rsb_options_t{
	/** An auxiliary bitmap */
	rsb_bitmap_data_t *bitmap;
	double a;
};
/* @endcond */


#if  RSB_FORTRAN_VERBOSE_CALLS
#define  RSB_FORTRAN_VERBOSE_CALL(M) RSB_ERROR(M)
#else /* RSB_FORTRAN_VERBOSE_CALLS */
#define  RSB_FORTRAN_VERBOSE_CALL(M)
#endif /* RSB_FORTRAN_VERBOSE_CALLS */


/*                                DEBUG FLAGS                                */
/*
	Enable any combination of the following flags to activate debug code.
	Do this only when debugging/developing, because it will slow down the code a lot, 
*/

#define RSB_WANT_DEBUG_PARANOID_ASSERTS 1	/**< if 1, will activate a number of assert() calls which won't change the code flow but will check for anomalous error conditions */
#if RSB_WANT_DEBUG_PARANOID_ASSERTS
#define RSB_DEBUG_ASSERT(e) assert(e)
#else /* RSB_WANT_DEBUG_PARANOID_ASSERTS */
#define RSB_DEBUG_ASSERT(e) 
#endif /* RSB_WANT_DEBUG_PARANOID_ASSERTS */

#define RSB_ASSERT(e) assert(e)		/* NOTE: in the future, could be replaced with some {if(..)goto err;} or exit()-like statement  */

/* commented out 20120915, since it was not used anyway most of the time
#ifdef DEBUG
#define RSB_DEBUG_BITMAP 1
#define RSB_DEBUG_BLOCK_STUFF 1
#define RSB_DEBUG_SORT_STUFF 1
#endif
*/

/*#define RSB_QUIET 1*/

/* FIXME : TODO : join these macros as a single debug flag */
#define RSB_MEM_DEBUG 		0	/* if 1, will make trigger printouts on allocations and deallocations */
#define RSB_QUIET_MEM_ERRORS	0	/* if 1, will not even print out fatal error conditions */

/*                            END DEBUG FLAGS                                */


/** Macros to check basic indices values validity ( FIXME : unfinished, should be much stricter )  */
#define RSB_INVALID_COO_INDEX(I)	((I)>RSB_MAX_MATRIX_DIM || (I)<0)	/* should fail only if signed */
#define RSB_INVALID_NNZ_INDEX(I)	((I)>RSB_MAX_MATRIX_NNZ)
#define RSB_INVALID_BLK_INDEX(I)	(RSB_INVALID_COO_INDEX(I))
#define RSB_INVALID_NNZ_COUNT(I)	((I)<1L || ( RSB_NNZ_ADD_OVERFLOW((I),RSB_INDEX_OF_SAFE_EXTRA) ))
#define RSB_INVALID_NNZ_COUNT_FOR_FLAGS(I,F) ((!RSB_DO_FLAG_HAS((F),RSB_FLAG_UNIT_DIAG_IMPLICIT)) && RSB_INVALID_NNZ_COUNT(I))

#define RSB_DO_FLAG_HAVE_AND(V1,V2,F) RSB_BOOL_AND(RSB_DO_FLAG_HAS(V1,F),RSB_DO_FLAG_HAS(V2,F))
#define RSB_DO_FLAG_HAVE_NAND(V1,V2,F) RSB_BOOL_NAND(RSB_DO_FLAG_HAS(V1,F),RSB_DO_FLAG_HAS(V2,F))
#define RSB_DO_FLAG_HAVE_XOR(V1,V2,F) RSB_BOOL_XOR(RSB_DO_FLAG_HAS(V1,F),RSB_DO_FLAG_HAS(V2,F))
#define RSB_DO_FLAG_HAVE_OR(V1,V2,F) RSB_BOOL_OR(RSB_DO_FLAG_HAS(V1,F),RSB_DO_FLAG_HAS(V2,F))
#define RSB_DO_FLAG_HAVE_NOR(V1,V2,F) RSB_BOOL_NOR(RSB_DO_FLAG_HAS(V1,F),RSB_DO_FLAG_HAS(V2,F))
#define RSB_DO_FLAG_SUBST(FLAGSVAR,FLAGS_OLD,FLAGS_NEW) RSB_DO_FLAG_DEL(FLAGSVAR,(FLAGS_OLD)), RSB_DO_FLAG_ADD(FLAGSVAR,(FLAGS_NEW))

#define RSB_INVALID_COO_COUNT(I)	((I)<1L || ( RSB_COO_ADD_OVERFLOW((I),RSB_INDEX_OF_SAFE_EXTRA) ))
#define RSB_IS_VALID_NNZ_COUNT(I)	(!RSB_INVALID_NNZ_COUNT(I))
#define RSB_IS_VALID_COO_INDEX(I)	(!RSB_INVALID_COO_INDEX(I))
#define RSB_IS_VALID_COO_DIM(I)		(!RSB_INVALID_COO_COUNT(I))
#define RSB_IS_VALID_BLK_INDEX(I)	(!RSB_INVALID_BLK_INDEX(I))
#define RSB_IS_VALID_NNZ_INDEX(I)	(!RSB_INVALID_NNZ_INDEX(I))
#define RSB_IS_VALID_INCX_VALUE(I)	(!RSB_INVALID_NNZ_COUNT(I))
#define RSB_ARE_VALID_MATRIX_INIT_PARS(R,C,NNZ,TYPE) (	\
	RSB_IS_VALID_NNZ_COUNT(NNZ)&&			\
	RSB_IS_VALID_COO_INDEX(R)&&			\
	RSB_IS_VALID_COO_INDEX(C)&&			\
	(!RSB_MATRIX_UNSUPPORTED_TYPE(TYPE)) )
/*#define RSB_IS_VALID_NNZ_SUM(NZ1,NZ2)	RSB_IS_VALID_NNZ_COUNT(((size_t)(NZ1))+((size_t)(NZ2)))*/
#define RSB_IS_VALID_NNZ_SUM(NZ1,NZ2)	(((size_t)(RSB_MAX_MATRIX_NNZ))>=(((size_t)(NZ1))+((size_t)(NZ2))))
#define RSB_IS_INVALID_TYPE_SIZE(TS) ((TS)<1)

#define RSB_IS_VALID_TRANS(T)  ((T)>=RSB_MIN(RSB_MIN(RSB_TRANSPOSITION_T,RSB_TRANSPOSITION_C),(RSB_TRANSPOSITION_N)) && (T)<=RSB_MAX(RSB_MAX(RSB_TRANSPOSITION_T,RSB_TRANSPOSITION_C),(RSB_TRANSPOSITION_N))) /* */
#define RSB_IS_VALID_THREAD_COUNT(C)	((C)> 0 && (C)<=RSB_CONST_MAX_SUPPORTED_CORES)	/* */
#define RSB_IS_VALID_THREAD_SPEC(C)	((C)>=0 && (C)<=RSB_CONST_MAX_SUPPORTED_CORES)	/* */

/** An initializer value for index variables. */
#define RSB_INI ((rsb_coo_idx_t)(-1))


#include <stdlib.h>		/* basic types and functions definitions */

/*
 * Bitmap stuff macros and functions.
 * Uses row major order by default.
 * By defining RSB_BITMAP_ROW_MAJOR_ORDER, row major order will be adopted.
 *
 * p.s.: please DO NOT use the following fixed macros outside the nearby macros.
 * */
/*#define RSB_BITMAP_ROW_MAJOR_ORDER 1*/
#define RSB_BITS_PER_INT  	(sizeof(rsb_bitmap_data_t)*RSB_CHAR_BIT)
#define RSB_BITS_PER_ROW(cols)  ((cols)+(RSB_BITS_PER_INT-1))
#define RSB_BYTES_PER_ROW(cols) ((cols+(RSB_CHAR_BIT-1))/RSB_CHAR_BIT)
#define RSB_INTS_PER_ROW(cols)   ((cols+((RSB_BITS_PER_INT)-1))/(RSB_BITS_PER_INT))
#define RSB_INT_IN_ROW(cols)   ((cols)/(RSB_BITS_PER_INT))

/* note that this is (logically) machine independent code */
#define RSB_SET_BIT(p,b)  (*(rsb_bitmap_data_t*)(p))=(*(rsb_bitmap_data_t*)(p)|(1<<(b)))
#define RSB_UNSET_BIT(p,b)  (*(rsb_bitmap_data_t*)(p))=(*(rsb_bitmap_data_t*)(p)&~(1<<(b)))
#define RSB_GET_BIT(p,b)  ((*(rsb_bitmap_data_t*)(p))&(1<<(b)))
#define RSB_BITMAP_POINTER(p,rw,r,c) (((rsb_bitmap_data_t*)(p))+(RSB_INTS_PER_ROW(rw)*(r)+RSB_INT_IN_ROW(c)))

#ifdef RSB_BITMAP_ROW_MAJOR_ORDER
/* Note that only a swap in the following three macros is needed to switch the storage format of our bitmap */
#define RSB_BITMAP_GET(p,rows,cols,r,c) RSB_GET_BIT((RSB_BITMAP_POINTER((p),(rows),(c),(r))),((r)%(RSB_BITS_PER_INT)))
#define RSB_BITMAP_SET(p,rows,cols,r,c) RSB_SET_BIT((RSB_BITMAP_POINTER((p),(rows),(c),(r))),((r)%(RSB_BITS_PER_INT)))
#define RSB_BITMAP_UNSET(p,rows,cols,r,c) RSB_UNSET_BIT((RSB_BITMAP_POINTER((p),(rows),(c),(r))),((r)%(RSB_BITS_PER_INT)))
#define RSB_BYTES_PER_BITMAP_(ld,d) (sizeof(rsb_bitmap_data_t)*(RSB_INTS_PER_ROW(ld)) * (d))
#define RSB_BITMAP_CLEAR(p,rows,cols)	RSB_BZERO((p),(RSB_BYTES_PER_BITMAP_(cols,rows)))
#else /* RSB_BITMAP_ROW_MAJOR_ORDER */
#define RSB_BITMAP_GET(p,rows,cols,r,c) RSB_GET_BIT((RSB_BITMAP_POINTER((p),(cols),(r),(c))),((c)%(RSB_BITS_PER_INT)))
#define RSB_BITMAP_SET(p,rows,cols,r,c) RSB_SET_BIT((RSB_BITMAP_POINTER((p),(cols),(r),(c))),((c)%(RSB_BITS_PER_INT)))
#define RSB_BITMAP_UNSET(p,rows,cols,r,c) RSB_UNSET_BIT((RSB_BITMAP_POINTER((p),(cols),(r),(c))),((c)%(RSB_BITS_PER_INT)))
#define RSB_BYTES_PER_BITMAP_(d,ld) (sizeof(rsb_bitmap_data_t)*(RSB_INTS_PER_ROW(ld)) * (d))
#define RSB_BITMAP_CLEAR(p,rows,cols)	RSB_BZERO((p),(RSB_BYTES_PER_BITMAP_(rows,cols)))
#endif /* RSB_BITMAP_ROW_MAJOR_ORDER */
#define RSB_WORDS_PER_BITMAP(ld,d) ((RSB_BYTES_PER_BITMAP_(ld,d)+(sizeof(rsb_bitmap_data_t)-1))/sizeof(rsb_bitmap_data_t))

#define RSB_BYTES_PER_BITMAP(rows,cols) RSB_BYTES_PER_BITMAP_(rows,cols)
#define RSB_BLOCK_UNSET_BIT_FOR_NNZ(IA,JA,k,M) {rsb_coo_idx_t i=(RSB_GET_BLOCK_ROW_FOR_NZ(IA,(M))); rsb_coo_idx_t j=(RSB_GET_BLOCK_COL_FOR_NZ(JA,M)); RSB_BITMAP_UNSET((M)->options->bitmap,(M)->M_b,(M)->K_b,i,j);}
#define RSB_BLOCK_SET_BIT_FOR_NNZ(IA,JA,k,M) {rsb_coo_idx_t i=(RSB_GET_BLOCK_ROW_FOR_NZ(IA+k,M)); rsb_coo_idx_t j=(RSB_GET_BLOCK_COL_FOR_NZ(JA+k,M)); RSB_BITMAP_SET((M)->options->bitmap,(M)->M_b,(M)->K_b,i,j);}
#define RSB_BLOCK_GET_BIT_FOR_NNZ(IA,JA,k,M) {rsb_coo_idx_t i=(RSB_GET_BLOCK_ROW_FOR_NZ(IA+k,M)); rsb_coo_idx_t j=(RSB_GET_BLOCK_COL_FOR_NZ(JA+k,M)); RSB_BITMAP_GET((M)->options->bitmap,(M)->M_b,(M)->K_b,i,j);}

/* 
 * Macros for one dimensional bitmaps -- easier to use.
 * */
#ifdef RSB_BITMAP_ROW_MAJOR_ORDER
#define RSB_BITVECTOR_GET(p,bits,bit)     RSB_BITMAP_GET(p,bits,1,bit,0) 
#define RSB_BITVECTOR_SET(p,bits,bit)     RSB_BITMAP_SET(p,bits,1,bit,0) 
#define RSB_BITVECTOR_UNSET(p,bits,bit)   RSB_BITMAP_UNSET(p,bits,1,bit,0)
#define RSB_BITVECTOR_CLEAR(p,bits)       RSB_BITMAP_CLEAR(p,bits,1)
#define RSB_BYTES_PER_BITVECTOR(bits)     RSB_BYTES_PER_BITMAP(bits,1)
#define RSB_WORDS_PER_BITVECTOR(bits)     RSB_WORDS_PER_BITMAP(bits,1)
#else /* RSB_BITMAP_ROW_MAJOR_ORDER */
#define RSB_BITVECTOR_GET(p,bits,bit)     RSB_BITMAP_GET(p,1,bits,0,bit) 
#define RSB_BITVECTOR_SET(p,bits,bit)     RSB_BITMAP_SET(p,1,bits,0,bit) 
#define RSB_BITVECTOR_UNSET(p,bits,bit)   RSB_BITMAP_UNSET(p,1,bits,0,bit)
#define RSB_BITVECTOR_CLEAR(p,bits)       RSB_BITMAP_CLEAR(p,1,bits)
#define RSB_BYTES_PER_BITVECTOR(bits)     RSB_BYTES_PER_BITMAP(1,bits)
#define RSB_WORDS_PER_BITVECTOR(bits)     RSB_WORDS_PER_BITMAP(1,bits)
#endif /* RSB_BITMAP_ROW_MAJOR_ORDER */

	/* given :
	 * the address of the nonzero element column index
	 * a rsb_options_t pointer
	 * will return the block column variable as it should.
	 * note that this macro will work only if o->M_b, o->K_b, o->p_r IA are properly initialized and o->p_r sorted.
	 *
	 * p.s.: note that the use of this could be avoided with a modest memory allocation...
	 * p.s.: in the following, we blindly trust that bsearch won't fail
	 * */
#define RSB_GET_BLOCK_COL_FOR_NZ_(columnidxp,cpntr,K_b) (((rsb_coo_idx_t*)bsearch((columnidxp),(cpntr),(K_b),sizeof(rsb_coo_idx_t),(rsb__nnz_coord_compar))-((cpntr))))
#define RSB_GET_BLOCK_COL_FOR_NZ(columnidxp,M)		(RSB_GET_BLOCK_COL_FOR_NZ_((columnidxp),(M)->cpntr,(M)->K_b))

#define RSB_GET_BLOCK_ROW_FOR_NZ_(rowidxp,rpntr,M_b) (((rsb_coo_idx_t*)bsearch((rowidxp)   ,(rpntr),(M_b),sizeof(rsb_coo_idx_t),(rsb__nnz_coord_compar))-((rpntr))))
#define RSB_GET_BLOCK_ROW_FOR_NZ(rowidxp   ,M)		(RSB_GET_BLOCK_ROW_FOR_NZ_((rowidxp),(M)->rpntr,(M)->M_b))

#define RSB_GET_BLOCK_MAJ_FOR_NZ_(majidxp,Mpntr,Md_b) (((rsb_coo_idx_t*)bsearch((majidxp)   ,(Mpntr),(Md_b),sizeof(rsb_coo_idx_t),(rsb__nnz_coord_compar))-((Mpntr))))
#define RSB_GET_BLOCK_MAJ_FOR_NZ(majidxp   ,M)		(RSB_GET_BLOCK_MAJ_FOR_NZ_((majidxp),(M)->Mpntr,(M)->Mdim))

#define RSB_GET_BLOCK_MIN_FOR_NZ_(minidxp,mpntr,md_b) (((rsb_coo_idx_t*)bsearch((minidxp)   ,(mpntr),(md_b),sizeof(rsb_coo_idx_t),(rsb__nnz_coord_compar))-((mpntr))))
#define RSB_GET_BLOCK_MIN_FOR_NZ(minidxp   ,M)		(RSB_GET_BLOCK_MIN_FOR_NZ_((minidxp),(M)->mpntr,(M)->mdim))

#define GET_BLOCK_FIRST_COLUMN(column,M)	((M)->cpntr[(column)])
#define GET_BLOCK_FIRST_ROW(row,M)		((M)->rpntr[ (row)  ])

#define GET_BLOCK_WIDTH(column,M)		(((M)->cpntr[(column)+1])-((M)->cpntr[(column)]))
#define GET_BLOCK_HEIGHT(row,M)			(((M)->rpntr[ (row)  +1])-((M)->rpntr[ (row)  ]))

#define GET_BLOCK_SIZE(row,column,M)	((GET_BLOCK_WIDTH((column),(M)))*(GET_BLOCK_HEIGHT((row),(M))))

/*
 * Blanks a whole matrix block.
 * */
#define RSB_BLANK_BLOCK(BP,M,BLOCKROW,BLOCKCOLUMN)				\
	RSB_BZERO( ((rsb_byte_t*)(BP)),						\
		( (M)->el_size * GET_BLOCK_SIZE((BLOCKROW),(BLOCKCOLUMN),(M)) ) );

#define RSB_INTRA_BLOCK_ROW(row,blockrow,M) ((row) - (M)->rpntr[(blockrow)])
#define RSB_INTRA_BLOCK_COLUMN(column,blockcolumn,M) ((column) - (M)->rpntr[(blockcolumn)])
#define RSB_GET_INTRA_BLOCK_OFFSET_ROW_MAJOR(row,column,blockrow,blockcolumn,M) ((( (row) - (M)->rpntr[(blockrow)]) * (GET_BLOCK_WIDTH((blockcolumn),(M))) + ( (column) - (M)->cpntr[(blockcolumn)] )) * (M)->el_size)
#define RSB_GET_INTRA_BLOCK_OFFSET_COLUMN_MAJOR(row,column,blockrow,blockcolumn,M) ((( (column) - (M)->cpntr[(blockcolumn)]) * (GET_BLOCK_HEIGHT((blockrow,(M))) + ( (row) - (M)->rpntr[(blockrow)] )) * (M)->el_size)

#define RSB_GET_INTRA_BLOCK_ROW_STRIDE(blockrow,blockcolumn,M) (GET_BLOCK_WIDTH((blockcolumn),(M)))
#define RSB_GET_INTRA_BLOCK_OFFSET(row,column,blockrow,blockcolumn,M) \
	(RSB_GET_INTRA_BLOCK_OFFSET_ROW_MAJOR(row,column,blockrow,blockcolumn,M)) 
#define RSB_GET_INTRA_BLOCK_OFFSET_TRANSPOSED(row,column,blockrow,blockcolumn,M) \
	(RSB_GET_INTRA_BLOCK_OFFSET_COLUMN_MAJOR(row,column,blockrow,blockcolumn,M)) 

/*!
 * Macros for diagonal-related comparisons.
 *
 * \code
 *  (ROW,COL)      (ROW,COL+COLS)
 *     +--------------+
 *     |              |
 *     |              |
 *     |              |
 *    ...            ...
 *     |              |
 *     |              |
 *     +--------------+
 *  (ROW+ROWS,COL)      (ROW+ROWS,COL+COLS)
 *
 *
 * 	under diagonal  	        over diagonal
 * 	intersects first at row COL     intersects first at row ROW
 * 	intersects last at row ROW+ROWS     intersects last at row COL+COLS
 *
 *	+---------------+       +-\-------------+
 *     \|               |       |  \            |
 *      \               |       |   \           |
 *      |\              |       |    \          |
 *     ...             ...     ...             ...
 *      |               |       |              \|
 *      |               |       |               \
 *      +-----\---------+	+---------------+
 * \endcode
 * */

#define RSB_POINT_QUASI_UNDER_DIAGONAL(ROW,COL) 	((ROW)>=(COL))
#define RSB_POINT_QUASI_OVER_DIAGONAL(ROW,COL)	 	((ROW)<=(COL))
#define RSB_POINT_UNDER_DIAGONAL(ROW,COL) 		((ROW)> (COL))
#define RSB_POINT_OVER_DIAGONAL(ROW,COL) 		((ROW)< (COL))

#define RSB_POINT_UNDER_SUPRA_DIAGONAL(ROW,COL,OFFSET) 	  (RSB_POINT_UNDER_DIAGONAL(((ROW)+(OFFSET)),(COL)))
#define RSB_POINT_UNDER_SUB_DIAGONAL(ROW,COL,OFFSET) 	  (RSB_POINT_UNDER_DIAGONAL(((ROW)),((COL)+(OFFSET))))
#define RSB_POINT_OVER_SUPRA_DIAGONAL(ROW,COL,OFFSET) 	  (RSB_POINT_OVER_DIAGONAL(((ROW)+(OFFSET)),(COL)))
#define RSB_POINT_OVER_SUB_DIAGONAL(ROW,COL,OFFSET) 	  (RSB_POINT_OVER_DIAGONAL((ROW),((COL)+(OFFSET))))

/* assumes COLS>=1, ROWS>=1 */
#define RSB_BLOCK_CROSSED_BY_DIAGONAL(ROW,COL,ROWS,COLS)	\
	(							\
	RSB_POINT_QUASI_UNDER_DIAGONAL((ROW)+(ROWS-1),(COL)) && 	\
	RSB_POINT_QUASI_OVER_DIAGONAL((ROW),((COL)+(COLS-1))) )
#define RSB_BLOCK_CROSSED_BY_SUPRA_DIAGONAL(ROW,COL,ROWS,COLS,OFFSET)	\
	RSB_BLOCK_CROSSED_BY_DIAGONAL((ROW)+(OFFSET),(COL),(ROWS),(COLS))
#define RSB_BLOCK_CROSSED_BY_SUB_DIAGONAL(ROW,COL,ROWS,COLS,OFFSET)	\
	RSB_BLOCK_CROSSED_BY_SUPRA_DIAGONAL(COL,ROW,COLS,ROWS,OFFSET)
#define RSB_BLOCK_CROSSED_BY_SUPRA_OR_SUB_DIAGONAL(ROW,COL,ROWS,COLS,LOFFSET,UOFFSET)	\
	(										\
	(LOFFSET)>(UOFFSET)?								\
	(RSB_BLOCK_CROSSED_BY_SUB_DIAGONAL(ROW,COL,ROWS,COLS,LOFFSET)):			\
	(RSB_BLOCK_CROSSED_BY_SUPRA_DIAGONAL(ROW,COL,ROWS,COLS,UOFFSET))	)		

/* 
 * The offset in the block to the first element which is on the diagonal 
 * The stride will be ROWS+1 or COLS+1, depending on the internal storage.
 * We here assume C storage.
 * */
#define RSB_BLOCK_DIAGONAL_OFFSET(ROW,COL,ROWS,COLS)	\
	((RSB_POINT_UNDER_DIAGONAL((ROW),(COL)))  ? ((ROW)-(COL)) : (((COL)-(ROW))*(COLS)) )
/* if the block is internally stored in Fortran, then : */
#define RSB_BLOCK_DIAGONAL_OFFSET_FORTRAN_STORED(ROW,COL,ROWS,COLS)	\
	((RSB_POINT_OVER_DIAGONAL((ROW),(COL)))  ? (((COL)-(ROW))*(ROWS)) : ((ROW)-(COL))  )

/* The following is ordering-neutral.

   +------------------->                 	
   | \   +-----+ ^^
   |   \ |     | ||
   |     \     | |v RSB_BLOCK_DIAGONAL_OFFSET_FIRST_ROW
   |     +-\---+ v  RSB_BLOCK_DIAGONAL_OFFSET_LAST_ROW
   |<--->    \      RSB_BLOCK_DIAGONAL_OFFSET
   |           \
   v  

   +------------------->                 	
   |    \        ^^
   |      \      ||
   |     +--\--+ v| RSB_BLOCK_SUPRA_DIAGONAL_OFFSET_FIRST_ROW
   |     |    \|  v RSB_BLOCK_SUPRA_DIAGONAL_OFFSET_LAST_ROW
   |     |     |\
   |     +-----+  \
   |
   v  

   +------------------->                 	
   |             ^^
   |             ||
   | \   +-----+ ||
   |   \ |     | ||
   |     \     | |v RSB_BLOCK_SUB_DIAGONAL_OFFSET_FIRST_ROW
   |     +-\---+ v  RSB_BLOCK_SUB_DIAGONAL_OFFSET_LAST_ROW
   |         \
   |           \
   v  
 */
#define RSB_BLOCK_DIAGONAL_OFFSET_FIRST_ROW(ROW,COL,ROWS,COLS)	\
	((RSB_POINT_UNDER_DIAGONAL((ROW),(COL)))  ? (ROW) :  (COL) )

#define RSB_BLOCK_DIAGONAL_OFFSET_LAST_ROW(ROW,COL,ROWS,COLS)	\
	((RSB_POINT_UNDER_DIAGONAL((((ROW)+(ROWS))-1),(((COL)+(COLS))-1)))  ? (((COL)+(COLS))-1):(((ROW)+(ROWS))-1)  )

#define RSB_BLOCK_SUPRA_DIAGONAL_OFFSET_FIRST_ROW(ROW,COL,ROWS,COLS,OFFSET)	\
	(RSB_BLOCK_DIAGONAL_OFFSET_FIRST_ROW(((ROW)+(OFFSET)),(COL),ROWS,COLS))

/*#define RSB_BLOCK_SUPRA_DIAGONAL_OFFSET_LAST_ROW(ROW,COL,ROWS,COLS,OFFSET)	*/
/*	RSB_BLOCK_DIAGONAL_OFFSET_LAST_ROW((ROW)+(OFFSET),COL,ROWS,COLS)*/

#define RSB_BLOCK_SUB_DIAGONAL_OFFSET_FIRST_ROW(ROW,COL,ROWS,COLS,OFFSET)	\
	(RSB_BLOCK_DIAGONAL_OFFSET_FIRST_ROW(ROW,(COL)+(OFFSET),ROWS,COLS))

#define RSB_BLOCK_SUB_OR_SUPRA_DIAGONAL_OFFSET_FIRST_ROW(ROW,COL,ROWS,COLS,LOFFSET,UOFFSET)	\
	((RSB_POINT_UNDER_DIAGONAL((ROW)+(UOFFSET),(COL)+(LOFFSET)))  ? (ROW) :  (COL)+(LOFFSET)-(UOFFSET) )

#define RSB_BLOCK_SUB_OR_SUPRA_DIAGONAL_OFFSET_FIRST_COL(ROW,COL,ROWS,COLS,LOFFSET,UOFFSET)	\
	RSB_BLOCK_SUB_OR_SUPRA_DIAGONAL_OFFSET_FIRST_ROW(COL,ROW,COLS,ROWS,UOFFSET,LOFFSET)

#define RSB_BLOCK_SUB_DIAGONAL_OFFSET_LAST_ROW(ROW,COL,ROWS,COLS,OFFSET)	\
	RSB_BLOCK_DIAGONAL_OFFSET_LAST_ROW((ROW),(COL)+(OFFSET),ROWS,COLS)

#define RSB_BLOCK_SUB_OR_SUPRA_DIAGONAL_OFFSET_LAST_ROW(ROW,COL,ROWS,COLS,LOFFSET,UOFFSET)	\
	((RSB_POINT_UNDER_DIAGONAL((ROW)+((ROWS)-1)+(UOFFSET),(COL)+((COLS)-1)+(LOFFSET)))  ? (COL)+((COLS)-1)+(LOFFSET)-(UOFFSET) : (ROW)+(ROWS)-1)

#define RSB_BLOCK_SUB_OR_SUPRA_DIAGONAL_OFFSET_LAST_COL(ROW,COL,ROWS,COLS,LOFFSET,UOFFSET)	\
	RSB_BLOCK_SUB_OR_SUPRA_DIAGONAL_OFFSET_LAST_ROW(COL,ROW,COLS,ROWS,UOFFSET,LOFFSET)


/* pure VBR, with no trailing structs : */

/* row major order (default) : */

#define	RSB_GET_NEXT_BLOCK_POINTER(BP,M,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	/*										\
	 * *input*									\
	 * M		should be a valid rsb_mtx_t structure pointer		\
	 * *output*									\
	 * ROWVAR	will be set to the base row    of this block			\
	 * COLVAR	will be set to the base column of this block			\
	 * BLOCKROWSVAR	will be set to the rows   count of this block			\
	 * BLOCKCOLSVAR	will be set to the column count of this block			\
	 * BP		 will be set to the current block pointer			\
	 * */										\
	++_k;										\
	if(_k>=(M)->bpntr[*_pi+1])							\
	{										\
		++*_pi;	/* new blocks row */						\
		while( (M)->bpntr[*_pi] == (M)->bpntr[*_pi+1] )		/* skipping empty rows */		\
			++*_pi;											\
	}													\
	*_pj=(M)->bindx[_k]; 						/* the current block column index  */	\
	_lastk=_k;												\
	(BLOCKROWVAR)=_i;											\
	(BLOCKCOLUMNVAR)=_j;											\
	(ROWVAR)=(M)->rpntr[_i];					/* _i is the current block row index */	\
	(COLVAR)=(M)->cpntr[_j]; 					/* the current block column index  */	\
	/*(BLOCKROWSVAR)=(M)->rpntr[_i+1]-(M)->rpntr[_i];*/ 		/* the current block rows    count */	\
	/*(BLOCKCOLSVAR)=(M)->cpntr[_j+1]-(M)->cpntr[_j];*/			/* the current block columns count */	\
	/*(BP)=(rsb_byte_t*)((M)->VA ) + (M)->el_size * (M)->indptr[_k] ; */						\
	(BLOCKROWSVAR)=GET_BLOCK_HEIGHT(_i,(M));	/* the current block rows    count */			\
	(BLOCKCOLSVAR)=GET_BLOCK_WIDTH( _j,(M)); 	/* the current block rows    count */			\
	(BP)=(rsb_byte_t*)(RSB_BLOCK_ADDRESS((M),_k));										\
	;

/* row major order : */
#define RSB_GET_FIRST_BLOCK_POINTER(BP,M,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	rsb_blk_idx_t _i=0,_j=0;										\
	rsb_blk_idx_t *_pi=NULL,*_pj=NULL;									\
	rsb_nnz_idx_t _k=0,_lastk=0;										\
	if((M)->flags&RSB_FLAG_WANT_COLUMN_MAJOR_ORDER){_pi=&_j;_pj=&_i;}else{_pi=&_i;_pj=&_j;} 		\
	while( (M)->bpntr[*_pi] == (M)->bpntr[*_pi+1] )								\
		++*_pi;												\
	_k=(M)->bpntr[*_pi]; 		/* _k is the first block index for the current row of blocks */		\
	*_pj=(M)->bindx[_k]; 						/* the current block column index  */	\
	(BLOCKROWVAR)=_i;											\
	(BLOCKCOLUMNVAR)=_j;											\
	(ROWVAR)=(M)->rpntr[_i];					/* _i is the current block row index */	\
	(COLVAR)=(M)->cpntr[_j]; 					/* the current block column index  */	\
	(BLOCKROWSVAR)=GET_BLOCK_HEIGHT(_i,(M));	/* the current block rows    count */			\
	(BLOCKCOLSVAR)=GET_BLOCK_WIDTH( _j,(M)); 	/* the current block rows    count */			\
	(BP)=(rsb_byte_t*)(RSB_BLOCK_ADDRESS((M),_k));										
#if RSB_WANT_DBC
#define RSB_GOT_LAST_BLOCK_POINTER(M)	( _lastk >= (M)->block_count )
#else
#define RSB_GOT_LAST_BLOCK_POINTER(M)	( _lastk >= (M)->nnz )
#endif


#define RSB_POINT_IN_BOX(R0,C0,RH,CW,R,C)	( ((R)>=(R0)) && (R)<((R0)+(RH)) && ((C)>=(C0)) && (C)<((C0)+(CW)) )
#define RSB_MATRIX_CONTAINS(M,R,C) 		( RSB_POINT_IN_BOX((M)->roff,(M)->coff,(M)->nr,(M)->nc,R,C) )
#define RSB_SUBMATRIX_CONTAINS_ROW(M,R) 		( RSB_POINT_IN_BOX((M)->roff,0,(M)->nr,0,R,1) )
#define RSB_SUBMATRIX_INTERSECTS_COLS(M,C0,C1) 							\
	   ( ( (M)->coff <= (C1) ) && ( (M)->coff+(M)->nc > (C0) ) )
#define RSB_SUBMATRIX_INTERSECTS_ROWS(M,R0,R1) 							\
	   ( ( (M)->roff <= (R1) ) && ( (M)->roff+(M)->nr > (R0) ) )

#define RSB_SUBMATRIX_INTERSECTS_BOX(M,R0,R1,C0,C1) 			\
	(RSB_SUBMATRIX_INTERSECTS_ROWS(M,R0,R1)&&RSB_SUBMATRIX_INTERSECTS_COLS(M,C0,C1))

#define RSB_FIND_SUBMATRIX_CONTAINING(M,R,C)	( \
	((M)->sm[0]&&RSB_MATRIX_CONTAINS((M)->sm[0],R,C)?(M)->sm[0]: \
	((M)->sm[1]&&RSB_MATRIX_CONTAINS((M)->sm[1],R,C)?(M)->sm[1]: \
	((M)->sm[2]&&RSB_MATRIX_CONTAINS((M)->sm[2],R,C)?(M)->sm[2]: \
	((M)->sm[3]&&RSB_MATRIX_CONTAINS((M)->sm[3],R,C)?(M)->sm[3]:NULL )))))

#define RSB_SUBMATRIX_INDEX(M,I,J) (M->sm[(I)*2+(J)])
/*
 * this should be fixed. we would prefer to use intrinsics here. TODO
 * */
#define RSB_FABS(x) ((x)<(0)?(-x):(x))

/*!
 * Misc macros.
 */
#define RSB_ASSIGN_IF_ZERO(VAR,VAL) if( (VAR) == 0) (VAR) = (VAL);
#define RSB_SWAP(TYPE,X,Y) {TYPE __tmp=(X);(X)=(Y);(Y)=(__tmp);}

#define RSB_SUBMATRIX_FOREACH_(matrix,submatrix,smi,smj,smk) 					\
	/*int smk;*/										\
	for(smk=0;smk<4;++smk)									\
	if( (smi=smk/2) >=0 && (smj=smk%2) >= 0 && (submatrix=(matrix)->sm[smi*2+(smj)] ) )	\
 	/* NOTE : handle with care (the 'submatrix' pointer could be NULL) */		\

#define RSB_SUBMATRIX_FOREACH_REVERSE(matrix,submatrix,smi,smj) 				\
	/*int smi,smj;*/								\
	for(smi=1;smi+1>0;--smi)/* fisrt smi, then smj, or will break spmv_uxux, ... */		\
	for(smj=1,submatrix=matrix->sm[smi*2+smj];					\
		smj+1>0;									\
		--smj,submatrix=(smi<2 && smj<2)?matrix->sm[smi*2+(smj)]:NULL)		\
 	/* NOTE : handle with care (the 'submatrix' pointer could be NULL) */		\

#define RSB_SUBMATRIX_FOREACH(MTXAP,submatrix,smi,smj) 				\
	/*int smi,smj;*/								\
	for(smi=0;smi<2;++smi)/* fisrt smi, then smj, or will break spmv_uxux, ... */		\
	for(smj=0,submatrix=MTXAP->sm[smi*2+smj];					\
		smj<2;									\
		++smj,submatrix=(smi<2 && smj<2)?MTXAP->sm[smi*2+(smj)]:NULL)		\
 	/* NOTE : handle with care (the 'submatrix' pointer could be NULL) */		\

/* The following is incorrect, as it accesses one further pointer. */
/*
#define RSB_SUBMATRIX_FOREACH_LEAF(MTXAP,submatrix,smi) 				\
	for(	(smi)=0,submatrix=(MTXAP)->all_leaf_matrices[(smi)].mtxlp;		\
		(smi)<(MTXAP)->all_leaf_matrices_n;					\
			++(smi),submatrix=(MTXAP)->all_leaf_matrices[smi].mtxlp)	\
*/

#define RSB_SUBMATRIX_FOREACH_LEAF_IDX(MTXAP,smi) 				\
	for(	smi = 0; smi < (MTXAP)->all_leaf_matrices_n ; ++smi)

/* The following is correct, even if less elegant because of the bad style of assignment. */
#define RSB_SUBMATRIX_FOREACH_LEAF(MTXAP,submatrix,smi) 				\
	for(	(smi)=0;		\
		((smi)<(MTXAP)->all_leaf_matrices_n) && ( submatrix=(MTXAP)->all_leaf_matrices[(smi)].mtxlp );	\
			++(smi))

#define RSB_SUBMATRIX_FOREACH_LEAF_PERMUTED(MTXAP,submatrix,smi,PV)				\
	for(	(smi)=0;		\
		((smi)<(MTXAP)->all_leaf_matrices_n) && ( submatrix=(MTXAP)->all_leaf_matrices[PV[(smi)]].mtxlp );	\
			++(smi))

#define RSB_SUBMATRIX_IS_ON_DIAG(matrix) 	((matrix)->roff==(matrix)->coff)
#define RSB_SUBMATRIX_IS_LOWDIAG(matrix) 	((matrix)->roff>(matrix)->coff)
#define RSB_SUBMATRIX_IS_UPPDIAG(matrix) 	((matrix)->roff<(matrix)->coff)

#define RSB_SUBMATRIX_FOREACH_DIAG_LEAF(matrix,submatrix,smi) 				\
	RSB_SUBMATRIX_FOREACH_LEAF(matrix,submatrix,smi) 				\
		if(RSB_SUBMATRIX_IS_ON_DIAG(submatrix))

#define RSB_SUBMATRIX_FOREACH_LOWDIAG_LEAF(matrix,submatrix,smi)			\
	RSB_SUBMATRIX_FOREACH_LEAF(matrix,submatrix,smi) 				\
		if(RSB_SUBMATRIX_IS_LOWDIAG(submatrix))

#define RSB_SUBMATRIX_FOREACH_UPPDIAG_LEAF(matrix,submatrix,smi) 				\
	RSB_SUBMATRIX_FOREACH_LEAF(matrix,submatrix,smi) 				\
		if(RSB_SUBMATRIX_IS_UPPDIAG(submatrix))

#define RSB_SUBMATRIX_COLS_INTERSECTION_FIRST(matrix,C)					\
	   RSB_MAX(((matrix)->coff),(C))

#define RSB_SUBMATRIX_COLS_INTERSECTION_LAST(matrix,C)					\
	   RSB_MIN(((matrix)->coff+(matrix->nc-1)),(C))	/* FIXME: requires matrix->nc > 0 */

#define RSB_SUBMATRIX_ROWS_INTERSECTION_FIRST(matrix,R)					\
	   RSB_MAX(((matrix)->roff),(R))

#define RSB_SUBMATRIX_ROWS_INTERSECTION_LAST(matrix,R)					\
	   RSB_MIN(((matrix)->roff+(matrix->nr-1)),(R))	/* FIXME: requires matrix->nr > 0 */

#define RSB_BCSS_MATRIX_FOREACH_BLOCK(matrix,blockpointer,bri,bci,blockindex,baserow,basecolumn,BR,BC)	\
	RSB_DEBUG_ASSERT((matrix)->VA);									\
	RSB_DEBUG_ASSERT((matrix)->el_size>0);								\
	RSB_DEBUG_ASSERT((matrix)->br>0 && (matrix)->bc>0);							\
	blockpointer=(matrix)->VA;									\
	for(	bri=0,											\
		baserow=(bri)*(BR);									\
		bri<(matrix)->Mdim;									\
		++bri,											\
		baserow=(bri)*(BR)									\
		)											\
	for(	bi=(matrix)->bpntr[bri],									\
		bci=(matrix)->bindx[bi],									\
		basecolumn=(bci)*(BC);								\
		bi<(matrix)->bpntr[(bri)+1];								\
		++bi,											\
		blockpointer=((rsb_byte_t*)blockpointer)+(matrix)->el_size*(BR)*(BC),		\
		bci=(matrix)->bindx[bi],									\
		baserow=(bri)*(BR),									\
		basecolumn=(bci)*(BC)								\
		)

#define RSB_BCSR_MATRIX_FOREACH_BLOCK(matrix,blockpointer,bri,bci,blockindex,baserow,basecolumn)	\
	RSB_BCSS_MATRIX_FOREACH_BLOCK(matrix,blockpointer,bri,bci,blockindex,baserow,basecolumn,matrix->br,matrix->bc)

#define RSB_BCSC_MATRIX_FOREACH_BLOCK(matrix,blockpointer,bri,bci,blockindex,baserow,basecolumn)	\
	RSB_BCSS_MATRIX_FOREACH_BLOCK(matrix,blockpointer,bci,bri,blockindex,basecolumn,baserow,matrix->bc,matrix->br)

#define RSB_CONST_ENOUGH_BYTES_FOR_ANY_TYPE RSB_MAX(RSB__MAX_SIZEOF,32)	/** max currently supported sizeof */
#define RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE (RSB_CONST_ENOUGH_BYTES_FOR_ANY_TYPE/sizeof(rsb_aligned_t))	/** should adapt this in case of need */


#define RSB_INTERNAL_FLAG_CSR_SORTING_MASK (RSB_FLAG_QUAD_PARTITIONING) /* TODO: abandon this */

#define RSB_DO_FLAGS_EXTRACT_STORAGE(F)	 ( \
		/*((F) & RSB_FLAG_WANT_LINKED_STORAGE) */ 0 | \
		((F) & RSB_FLAG_WANT_COO_STORAGE) | \
		((F) & RSB_FLAG_WANT_FIXED_BLOCKING_VBR) | \
		((F) & RSB_FLAG_WANT_BCSS_STORAGE) | \
		RSB_FLAG_NOFLAGS )

#define RSB__FLAG_HAS_UNSPECIFIED_TRIANGLE(FLAGS) (RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_TRIANGULAR) && !RSB_DO_FLAG_HAS_INTERSECTION((FLAGS),RSB_FLAG_LOWER|RSB_FLAG_UPPER))

#if 1
#define rsb_do_spmv(TRANSA,ALPHAP,MTXAP,XP,INCX,BETAP,YP,INCY)	\
       	rsb__do_spmv_general(TRANSA,ALPHAP,MTXAP,XP,INCX,BETAP,YP,INCY,(RSB_OP_FLAG_DEFAULT) RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS)
#else
rsb_err_t rsb_do_spmv(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY)
{
	rsb_err_t errval = rsb__do_spmv_general(transA,alphap,mtxAp,Xp,incX,betap,Yp,incY,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
	return errval;
}
#endif

/* We may use custom memcpy functions. */
#define RSB_MEMCPY(DST,SRC,BYTES) rsb__memcpy((DST),(SRC),(BYTES))

#define RSB_A_BZERO(ID,DOFF,NNZ,ES) \
	RSB_BZERO(((rsb_byte_t*)(ID))+(ES)*(DOFF),(ES)*(NNZ)) \

#define RSB_A_MEMCPY(ID,IS,DOFF,SOFF,NNZ,ES) \
	RSB_MEMCPY(((rsb_char_t*)(ID))+(ES)*(DOFF),((const rsb_char_t*)(IS))+(ES)*(SOFF),(ES)*(NNZ)) \

#define RSB_A_MEMMOVE(ID,IS,DOFF,SOFF,NNZ,ES) \
	RSB_MEMMOVE(((rsb_char_t*)(ID))+(ES)*(DOFF),((const rsb_char_t*)(IS))+(ES)*(SOFF),(ES)*(NNZ)) \

#define RSB_COA_MEMCPY(ID,IS,DOFF,SOFF,NNZ) \
	RSB_MEMCPY(((rsb_coo_idx_t*)(ID))+(DOFF),((const rsb_coo_idx_t*)(IS))+(SOFF),sizeof(rsb_coo_idx_t)*(NNZ)) \

#define RSB_COA_MEMCPY2H(ID,IS,DOFF,SOFF,NNZ,ADD) 					\
{											\
	rsb_nnz_idx_t RSB_DUMMY_ID=0;							\
	for(RSB_DUMMY_ID=0;RSB_DUMMY_ID<(NNZ);++RSB_DUMMY_ID)				\
		((rsb_half_idx_t*)(ID))[(DOFF)+(RSB_DUMMY_ID)]=			\
		((const rsb_coo_idx_t*)IS)[(SOFF)+(RSB_DUMMY_ID)]+(ADD);		\
}

#define RSB_COA_MEMMOVE(ID,IS,DOFF,SOFF,NNZ) \
	RSB_MEMMOVE(((rsb_coo_idx_t*)(ID))+(DOFF),((const rsb_coo_idx_t*)(IS))+(SOFF),sizeof(rsb_coo_idx_t)*(NNZ)) \

#define RSB_COO_MEMMOVE(VD,ID,JD,VS,IS,JS,DOFF,SOFF,NNZ,ES) \
	RSB_MEMMOVE(((rsb_char_t*)(VD))+(ES)*(DOFF),((const rsb_char_t*)(VS))+(ES)*(SOFF),(ES)*(NNZ)), \
	RSB_MEMMOVE(((rsb_coo_idx_t*)(ID))+(DOFF),((const rsb_coo_idx_t*)(IS))+(SOFF),sizeof(rsb_coo_idx_t)*(NNZ)), \
	RSB_MEMMOVE(((rsb_coo_idx_t*)(JD))+(DOFF),((const rsb_coo_idx_t*)(JS))+(SOFF),sizeof(rsb_coo_idx_t)*(NNZ))

#define RSB_COO_MEMCPY(VD,ID,JD,VS,IS,JS,DOFF,SOFF,NNZ,ES) \
	RSB_MEMCPY(((rsb_char_t*)(VD))+(ES)*(DOFF),((const rsb_char_t*)(VS))+(ES)*(SOFF),(ES)*(NNZ)), \
	RSB_MEMCPY(((rsb_coo_idx_t*)(ID))+(DOFF),((const rsb_coo_idx_t*)(IS))+(SOFF),sizeof(rsb_coo_idx_t)*(NNZ)), \
	RSB_MEMCPY(((rsb_coo_idx_t*)(JD))+(DOFF),((const rsb_coo_idx_t*)(JS))+(SOFF),sizeof(rsb_coo_idx_t)*(NNZ))

#define RSB_CSR_MEMCPY(VD,ID,JD,VS,IS,JS,NNZ,NR,ES) \
	RSB_MEMCPY(((rsb_char_t   *)(VD)),((const rsb_char_t   *)(VS)),(ES)*(NNZ)), \
	RSB_MEMCPY(((rsb_nnz_idx_t*)(ID)),((const rsb_nnz_idx_t*)(IS)),sizeof(rsb_nnz_idx_t)*(NR)), \
	RSB_MEMCPY(((rsb_coo_idx_t*)(JD)),((const rsb_coo_idx_t*)(JS)),sizeof(rsb_coo_idx_t)*(NNZ))

#define RSB_CSR2COO_MEMCPY(VD,ID,JD,VS,I,JS,DOFF,SOFF,NNZ,ES) \
	RSB_MEMCPY(((rsb_char_t*)(VD))+(ES)*(DOFF),((const rsb_char_t*)(VS))+(ES)*(SOFF),(ES)*(NNZ)), \
	rsb__util_coo_array_set(((rsb_coo_idx_t*)(ID))+(DOFF),(NNZ),(I)), \
	RSB_MEMCPY(((rsb_coo_idx_t*)(JD))+(DOFF),((const rsb_coo_idx_t*)(JS))+(SOFF),sizeof(rsb_coo_idx_t)*(NNZ))

#define RSB_COO_MEMCPY_parallel(VD,ID,JD,VS,IS,JS,DOFF,SOFF,NNZ,ES) \
	RSB_A_MEMCPY_parallel(VD,VS,DOFF,SOFF,NNZ,ES), \
	RSB_COA_MEMCPY_parallel(ID,IS,DOFF,SOFF,NNZ), \
	RSB_COA_MEMCPY_parallel(JD,JS,DOFF,SOFF,NNZ)

#define RSB_FCOO_ASUM(S,X,LI,UI) {rsb_coo_idx_t i; for(i=(LI);i<(UI);++i)(S)+=(X)[i];}
#define RSB_XCOO_ISET(X,  LI,UI) {rsb_coo_idx_t i; for(i=(LI);i<(UI);++i)(X)[i]=i-(LI);}
#define RSB_FCOO_ISET RSB_XCOO_ISET
#define RSB_XCOO_VSET(X,V,LI,UI) {rsb_coo_idx_t i; for(i=(LI);RSB_LIKELY((i)<(UI));++i)(X)[(i)] =(V);}
#define RSB_XCOO_VADD(X,V,LI,UI) {rsb_coo_idx_t i; for(i=(LI);RSB_LIKELY((i)<(UI));++i)(X)[(i)]+=(V);}
#define RSB_XCOO_IREN	/* TODO: to write one */

#define RSB_TYPED_OFF_PTR(TYPECODE,VA,OFF) (((rsb_byte_t*)(VA))+(((size_t)(RSB_SIZEOF(TYPECODE))*(OFF))))
#define RSB_COO_LT(I1,J1,I2,J2) ( (I1) < (I2) || ( (I1) == (I2) && ( (J1) < (J2) ) ) )
#define RSB_COO_GT(I1,J1,I2,J2) RSB_COO_LT(I2,J2,I1,J1)
#define RSB_COO_GE(I1,J1,I2,J2) ( !RSB_COO_LT(I1,J1,I2,J2) )
#define RSB_COO_EQ(I1,J1,I2,J2) ( (I1) == (I2) && ( (J1) == (J2) ) )

/*!
 * \ingroup gr_internals
 * \brief An internal, helper structure.
 */
struct rsb_memory_level_t
{
	size_t size;				/*  */
	size_t level;				/*  */
	size_t associativity;			/*  */
	size_t linesize;			/*  */
};

#define RSB_MEGABYTE (1024*1024)
#define RSB_MEGABYTE_SYM "MiB"

#define RSB_DEFAULT_STREAM stdout
#define RSB_DIR_SEPARATOR	'/'	/*  */
#define RSB_MAX_STRERRLEN  	128	/*  */
#define RSB_MAX_LINE_LENGTH  	1025	/*  */
#define RSB_MAX_COMPILE_COMMAND_LENGTH 	1025	/*  */
#define RSB_MAX_VERSION_STRING_LENGTH  	4096	/*  */
#define RSB_MAX_FILENAME_LENGTH  RSB_MAX_LINE_LENGTH	/* the maximal supported file name length (in buffers) */

#define RSB_MAX_SUPPORTED_CACHE_LEVELS 32L	/* the maximal supported height of memory hierarchy */
#define RSB_MIN_THREAD_MEMCPY_NNZ 1024		/* minimal count of nonzeros to move for a thread during parallel memcpy */
#define RSB_MIN_THREAD_BZERO_BYTES 8192		/* minimal count of nonzeros to bzero for a thread during parallel bzero */
#define RSB_MIN_THREAD_XAXPY_NNZ 256 /* 1024 */		/* minimal count of elements for a parallel vector-vector operation */
#define RSB_MIN_THREAD_SORT_NNZ 256		/* minimal count of nonzeros to sort for a thread during parallel sort */
#define RSB_MIN_NNZ_FOR_PARALLEL_ADD_TO_DENSE 100000

#define RSB_POWER_OF_2(N) (1<<(N))
#define RSB_FRAC(Q,D) (((Q)+((D)-1))/(D))
#define RSB_MIDDLE(X) RSB_FRAC(X,2)
#define RSB_IS_INTEGER_ODD(X)   ( (X)&0x01)
#define RSB_IS_INTEGER_EVEN(X)	(!RSB_IS_INTEGER_ODD(X))

#define RSB_HAVE_STREAMS RSB_HAVE_STDIO_H

#if defined(RSB_WANT_LIBRSB_STATS) && (RSB_WANT_LIBRSB_STATS> 0)
#define RSB_WANT_LIBRSB_TIMER 1
#else
#define RSB_WANT_LIBRSB_TIMER 0
#endif

#define RSB_WANT_SPMV_TRACE RSB_ALLOW_INTERNAL_GETENVS /* internal, undocumented */
#define RSB_WANT_COO2RSB_THREADS RSB_ALLOW_INTERNAL_GETENVS /* internal, undocumented */
#define RSB_WANT_CACHE_TIMER_GRANULARITY 1 /* internal, undocumented */

/*!
 * \ingroup gr_internals
 * \brief An internal, helper structure.
 */
struct rsb_session_handle_t
{
	#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	/*!
	 * A global memory counter, used for debugging purposes.
	 * */
	size_t allocated_memory;			/* total of allocated memory, in bytes */
	size_t allocations_count;		/* total number of current allocations */
	size_t allocations_cumulative;		/* cumulative number of memory allocations */
	#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	size_t min_leaf_matrix_bytes;		/*  */
	size_t avg_leaf_matrix_bytes;		/*  */
	size_t rsb_g_threads;			/* detected threads */
#if RSB_WANT_PERFORMANCE_FILE
	/*rsb_byte_t * performance_binary_dump_file;*/	/* TODO: obsolete feature */
	rsb_char_t * performance_binary_dump_file;	/* TODO: obsolete feature */
#endif /* RSB_WANT_PERFORMANCE_FILE */
	/* beginning of user settable variables declarations */
	size_t rsb_want_threads;		/* RSB_IO_WANT_EXECUTING_THREADS ; active threads (may be <> rsb_g_threads) */
	rsb_int_t asm_sort_method;		/* RSB_IO_WANT_SORT_METHOD */
	rsb_real_t subdivision_multiplier;	/* RSB_IO_WANT_SUBDIVISION_MULTIPLIER */
	rsb_int_t want_bounded_box;		/* RSB_IO_WANT_BOUNDED_BOX_COMPUTATION */
	rsb_int_t cache_blocking_method;	/* RSB_IO_WANT_CACHE_BLOCKING_METHOD */
	rsb_int_t want_outer_spmm;		/* RSB_IO_WANT_LEAF_LEVEL_MULTIVEC */
#if RSB_HAVE_STREAMS
	FILE * out_stream;			/* RSB_IO_WANT_OUTPUT_STREAM */
	FILE * error_stream;			/* RSB_IO_WANT_VERBOSE_ERRORS */
	FILE * init_stream;			/* RSB_IO_WANT_VERBOSE_INIT */
	FILE * exit_stream;			/* RSB_IO_WANT_VERBOSE_EXIT */
#endif /* RSB_HAVE_STREAMS */
	const rsb_char_t * mhis;		/* RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING ; set via rsb_lib_reinit */
#if RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE
	rsb_int_t rsb_g_verbose_interface;	/* RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE */
#endif /* RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE */
	/* end of user settable variables declarations */
	long memory_hierarchy_levels;		/*  */
	struct rsb_memory_level_t caches[RSB_MAX_SUPPORTED_CACHE_LEVELS];	/* 0,..,memory_hierarchy_levels-1*/
	rsb_bool_t rsb_g_initialized;		/*  */
#if RSB_WANT_ALLOCATOR_LIMITS
       	size_t memory_count_max;		/*  */
	size_t allocations_count_max;		/*  */
#endif /* RSB_WANT_ALLOCATOR_LIMITS */
#if RSB_WANT_LIBRSB_TIMER
	rsb_time_t etime;
#endif /* RSB_WANT_LIBRSB_TIMER */
	rsb_int_t verbose_tuning;		/*  */
#if RSB_USE_LIBRSBPP
	rsb_int_t use_rsbpp; /*  */
#endif /* RSB_USE_LIBRSBPP */
#if RSB_USE_MKL
	rsb_int_t use_mkl; /*  */
#endif /* RSB_USE_MKL */
#if RSB_WANT_SPMV_TRACE
	rsb_int_t want_spmv_trace; /* allow EPS plots as traces from SpMV / SpMM's */
#endif /* RSB_WANT_SPMV_TRACE */
#ifdef RSB_WANT_COO2RSB_THREADS
	size_t coo2rsb_threads; /* */
#endif /* RSB_WANT_COO2RSB_THREADS */
#ifdef RSB_WANT_CACHE_TIMER_GRANULARITY
	rsb_time_t timer_granularity; /* */
#endif /* RSB_WANT_CACHE_TIMER_GRANULARITY */
};

#define RSB_INTERNALS_COMMON_HEAD_DECLS extern struct rsb_session_handle_t rsb_global_session_handle;
#define RSB_DO_ERROR_CUMULATE(ERRVAL,ERRFLAG) RSB_DO_FLAG_ADD((ERRVAL),(ERRFLAG))

#define RSB_IF_NOT_NULL_CAST_TO(P,TYPE,FALLBACK) ((P)?*(TYPE*)(P):(FALLBACK))
#define RSB_SET_TO_CASTED(V,P,TYPE) {(V)=*(TYPE*)(P);}
#define RSB_IF_NOT_NULL_SET_TO_CASTED(V,P,TYPE) {if((P)!=NULL)RSB_SET_TO_CASTED(V,P,TYPE)}
#define RSB_IF_NOT_NULL_GET_TO_CASTED(V,P,TYPE) {if((P)!=NULL){*(TYPE*)(P)=(V);}}
#define RSB_IF_NOT_NULL_GET_SET_TO_CASTED(V,P,TYPE,F,ERRVAL)	{	\
	switch(F){							\
		case(RSB_IO_SPECIFIER_GET):				\
		RSB_IF_NOT_NULL_GET_TO_CASTED((V),(P),TYPE);break;		\
		case(RSB_IO_SPECIFIER_SET):				\
		RSB_IF_NOT_NULL_SET_TO_CASTED((V),(P),TYPE);break;	\
		default: RSB_DO_ERROR_CUMULATE(ERRVAL,RSB_ERR_BADARGS); }}
#define rsb__sprintf sprintf

#define RSB_BLAS_ERROR		-1	/* */
#define RSB_BLAS_NO_ERROR 	0	/* */
#define RSB_BLAS_ERROR_UNSUPPORTED   RSB_BLAS_ERROR			/* TODO: spread usage of this throughout the code */
#define RSB_BLAS_ERROR_UNIMPLEMENTED RSB_BLAS_ERROR			/* TODO: spread usage of this throughout the code */
#define RSB_BLAS_ERROR_WRONG_USGP_ARG RSB_BLAS_ERROR			/* TODO: spread usage of this throughout the code */

#define RSB_SET_IF_NOT_NULL(P,V) if((P)!=NULL)*(P)=V
#ifdef RSB_WANT_LONG_IDX_TYPE 
typedef RSB_WANT_LONG_IDX_TYPE rsb_blas_int_t;
#else /* RSB_WANT_LONG_IDX_TYPE */
typedef int                    rsb_blas_int_t;
#endif /* RSB_WANT_LONG_IDX_TYPE */
typedef double rsb_aligned_t;	/* see RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE and RSB_CONST_ENOUGH_BYTES_FOR_ANY_TYPE */


#define RSB_MASK_OUT_SOME_ERRORS(ERRVAL) {if((ERRVAL)==RSB_ERR_UNSUPPORTED_FEATURE)(ERRVAL)=RSB_ERR_NO_ERROR;}/* NOTE; this is a macro only used to prevent the test suite to complain for failing printouts when output is disabled! */

/*!
 Macros to get indices types liminal values, configuration-dependent.
*/
 #define RSB_COO_HALF_BITS_SIZE	((sizeof(rsb_half_idx_t)*RSB_CHAR_BIT)) /* new: assuming RSB_COO_HALF_BITS_SIZE	refers to rsb_half_idx_t */
 /*#define RSB_COO_HALF_BITS_SIZE	((sizeof(rsb_coo_idx_t)*RSB_CHAR_BIT)/2) */ /* old: assuming rsb_half_idx_t is half of rsb_coo_idx_t */

#define RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS {int ___foo=1;++___foo;}/* will avoid things like error: label at end of compound statement */

#if RSB_WANT_ZLIB_SUPPORT
#define RSB_FOPEN(X,Y) (FILE*)gzopen((X),(Y))
#define RSB_FCLOSE(X) gzclose((gzFile)(X))
#define RSB_GETC(X) gzgetc((gzFile)X)
#define RSB_UNGETC(X,Y) gzungetc(X,(gzFile)Y)
#define RSB_FREAD(BUF,SIZE,NITEMS,FILE) gzfread(BUF,SIZE,NITEMS,FILE)
#else /* RSB_WANT_ZLIB_SUPPORT */
#define RSB_FOPEN(X,Y) fopen((X),(Y))
#define RSB_FCLOSE(X) fclose(X)
#define RSB_GETC(X) getc(X)
#define RSB_UNGETC ungetc
#define RSB_FREAD(BUF,SIZE,NITEMS,FILE) fread(BUF,SIZE,NITEMS,FILE)
#endif /* RSB_WANT_ZLIB_SUPPORT */

#define RSB_EMPTY_FILE_FILLER \
	static void foo(void);		\
	static void baz(void);		\
	static void baz(void){foo();}	\
	static void foo(void){baz();}	/* We want: no empty translation unit; no used function.. */

#define RSB_DECLARE_COO_ARRAYS_FROM_MATRIX(IA,JA,MATRIX,TYPE) 	\
		TYPE *IA=(TYPE*)(MATRIX)->bpntr;			\
		TYPE *JA=(TYPE*)(MATRIX)->bindx;

#define RSB_DECLARE_COO_IARRAY_FROM_MATRIX(IA,MATRIX,TYPE) 	\
		TYPE *IA=(TYPE*)(MATRIX)->bpntr;

#define RSB_DECLARE_COO_JARRAY_FROM_MATRIX(JA,MATRIX,TYPE) 	\
		TYPE *JA=(TYPE*)(MATRIX)->bindx;

#define RSB_DECLARE_CSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX,PTYPE,ITYPE) 	\
		PTYPE *PA=(PTYPE*)(MATRIX)->bpntr;			\
		ITYPE *JA=(ITYPE*)(MATRIX)->bindx;

#define RSB_DECLARE_CONST_HALFCSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX) 	\
	RSB_DECLARE_CSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX,const rsb_nnz_idx_t,const rsb_half_idx_t)

#define RSB_DECLARE_HALFCSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX) 	\
	RSB_DECLARE_CSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX,rsb_nnz_idx_t,rsb_half_idx_t)

#define RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX) 	\
	RSB_DECLARE_CSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX,const rsb_nnz_idx_t,const rsb_coo_idx_t)

#define RSB_DECLARE_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX) 	\
	RSB_DECLARE_CSR_ARRAYS_FROM_MATRIX(PA,JA,MATRIX,rsb_nnz_idx_t,rsb_coo_idx_t)

#define RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,MATRIX) 	\
	RSB_DECLARE_COO_ARRAYS_FROM_MATRIX(IA,JA,MATRIX,const rsb_half_idx_t)

#define RSB_DECLARE_CONST_HALFCOO_IARRAY_FROM_MATRIX(IA,MATRIX) 	\
	RSB_DECLARE_COO_IARRAY_FROM_MATRIX(IA,MATRIX,const rsb_half_idx_t)

#define RSB_DECLARE_CONST_HALFCOO_JARRAY_FROM_MATRIX(JA,MATRIX) 	\
	RSB_DECLARE_COO_JARRAY_FROM_MATRIX(JA,MATRIX,const rsb_half_idx_t)

#define RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,MATRIX) 	\
	RSB_DECLARE_COO_ARRAYS_FROM_MATRIX(IA,JA,MATRIX,const rsb_coo_idx_t)

#define RSB_DECLARE_CONST_FULLCOO_IARRAY_FROM_MATRIX(IA,MATRIX) 	\
	RSB_DECLARE_COO_IARRAY_FROM_MATRIX(IA,MATRIX,const rsb_coo_idx_t)

#define RSB_DECLARE_CONST_FULLCOO_JARRAY_FROM_MATRIX(JA,MATRIX) 	\
	RSB_DECLARE_COO_JARRAY_FROM_MATRIX(JA,MATRIX,const rsb_coo_idx_t)

#define RSB_DECLARE_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,MATRIX) 	\
	RSB_DECLARE_COO_ARRAYS_FROM_MATRIX(IA,JA,MATRIX,rsb_coo_idx_t)

/*!
 * If the restrict keyword is supported, we use it in our declarations.
 * */
/* #ifdef restrict */
#ifdef RSB_restrict
#define RSB_RESTRICT restrict
#else /* RSB_restrict */
#define RSB_RESTRICT
#endif /* RSB_restrict */

#define RSB_VA_OFFSET_POINTER(VA,ES,OFF) 		((rsb_byte_t*)(VA)+(size_t)(ES)*(OFF))
#define RSB_VA_OFFSET_POINTER_CONST(VA,ES,OFF) 		((const rsb_byte_t*)(VA)+(size_t)(ES)*(OFF))
#define RSB_VA_MEMCMP(LVA,LOFF,RVA,ROFF,ES) 		\
	RSB_MEMCMP(RSB_VA_OFFSET_POINTER((LVA),(ES),(LOFF)),RSB_VA_OFFSET_POINTER((RVA),(ES),(ROFF)),(ES))		

/*!
 * \brief Auxiliary structure for a coo-stored matrix (usually for temporary operations).
 * */
struct rsb_coo_mtx_t{
	rsb_coo_idx_t * IA, * JA;/** row and columns indices */
	rsb_coo_idx_t nr,nc;	/** matrix (declared) nonzeros */
	rsb_nnz_idx_t nnz;	/** matrix rows, columns */
	void * VA;		/** values of data elements */
	rsb_type_t typecode;	/** as specified in the RSB_NUMERICAL_TYPE_* preprocessor symbols in rsb_types.h 	*/
};

#define RSB_INIT_COO_FROM_MTX(COOP,MTXAP)	{ \
		(COOP)->nr=(MTXAP)->nr;	\
		(COOP)->nc=(MTXAP)->nc;	\
		(COOP)->nnz=(MTXAP)->nnz;	\
		(COOP)->typecode=(MTXAP)->typecode; }

#define RSB_INIT_CXX_FROM_MTX(COOP,MTXAP)	{ \
		(COOP)->nr=(MTXAP)->nr;	\
		(COOP)->nc=(MTXAP)->nc;	\
		(COOP)->nnz=RSB_MAX((MTXAP)->nnz,1+RSB_MAX((MTXAP)->nr,(MTXAP)->nc)); \
		(COOP)->typecode=(MTXAP)->typecode; }

#define RSB_BIND_COO_TO_MTX(COOP,MTXAP)	{ \
		(COOP)->VA=(MTXAP)->VA;	\
		(COOP)->IA=(MTXAP)->bpntr;	\
		(COOP)->JA=(MTXAP)->bindx;	}

#define RSB_FLAG_ANY_SYMMETRY				(RSB_FLAG_HERMITIAN|RSB_FLAG_SYMMETRIC)
#define RSB_FLAG_ALL_STRUCTURAL_FLAGS	(RSB_FLAG_ANY_SYMMETRY|RSB_FLAG_DIAGONAL|RSB_FLAG_TRIANGULAR|RSB_FLAG_UNIT_DIAG_IMPLICIT)
#define RSB_FLAG_ALL_DUPLICATE_FLAGS	(RSB_FLAG_DUPLICATES_KEEP_LAST|RSB_FLAG_DUPLICATES_SUM)
#define RSB_DUMMY_ID		rsb_dummy_id
#define RSB_DUMMY_MTX		(NULL)
#define RSB_DEFAULT_TEST_MATRIX_FILENAME "pd.mtx"	/**< this file should always be included in the library distribution (FIXME: should enforce this) */

#define RSB_VECTORS_DIFF_DISPLAY_N 10
#define RSB_VECTORS_DIFF_DISPLAY_N_SMALL 3
#define RSB_DEFAULT_UNDEFINED_COO_VALUE 0

#define RSB_PSORT_CHUNK 10000			/* FIXME: hardcoded constants are bad (and the PGI compiler won't accept them) */
#define RSB_MINIMUM_VECOP_OMP_CHUNK 1000			/* FIXME: hardcoded constants are bad (and the PGI compiler won't accept them) */

#define RSB_BOOL_IS_POINTER_NON_NULL(P) ((P)?RSB_BOOL_TRUE:RSB_BOOL_FALSE)

#define RSB_CONDITIONAL_ERRPSET(ERRVALP,ERRVAL) {if(ERRVALP)(*(ERRVALP)=(ERRVAL));}
#define RSB_MTX_FREE(MTXAP) if(MTXAP){rsb__do_mtx_free(MTXAP);(MTXAP)=NULL;}  /* frees the matrix and nullifies the associated pointer. */

/* #define RSB_FLAGS_RSB_AGNOSTIC RSB_FLAG_FORTRAN_INDICES_INTERFACE */
/* #define RSB_FLAGS_RSB_NON_AGNOSTIC (RSB_FLAG_USE_HALFWORD_INDICES|RSB_FLAG_WANT_COO_STORAGE|RSB_FLAG_WANT_CSR_STORAGE)  --- see RSB_DO_FLAGS_EXTRACT_STORAGE(flags) */
#define RSB_FLAGS_RSB_AGNOSTIC (RSB_FLAG_FORTRAN_INDICES_INTERFACE|RSB_FLAG_UNIT_DIAG_IMPLICIT|RSB_FLAG_UPPER|RSB_FLAG_LOWER|RSB_FLAG_SORTED_INPUT|RSB_FLAG_TRIANGULAR|RSB_FLAG_SYMMETRIC|RSB_FLAG_HERMITIAN)

#define RSB_INDEX_FIT_IN_HALFWORD(I) ((I)<=RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t))
#define RSB_INDICES_FIT_IN_HALFWORD(I,J) ( RSB_INDEX_FIT_IN_HALFWORD(I) && RSB_INDEX_FIT_IN_HALFWORD(J) )

#define RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(V)  						\
	if(RSB_DO_FLAG_FILTEROUT((V),RSB_FLAGS_RSB_AGNOSTIC)==RSB_FLAG_NOFLAGS)	\
 		RSB_DO_FLAG_ADD((V),RSB_FLAG_DEFAULT_MATRIX_FLAGS);

#define RSB_DIVIDE_IN_CHUNKS(N,NTHREADS) RSB_MAX(((N)+(NTHREADS)-1)/(NTHREADS),1)
#define RSB_EXIT exit
#define RSB_DO_ERR_RETURN(ERRVAL) {return (ERRVAL);}
#define RSB_DO_MTX_RETURN(MATRIX,ERRVAL) {return (MATRIX);}
#define RSB_FLAG_UPPTRI (RSB_FLAG_UPPER|RSB_FLAG_LOWER)
#define RSB_DO_FLAG_FLIP_UPLO(V)	{\
if(RSB_DO_FLAG_HAS((V),RSB_FLAG_UPPER)) \
	RSB_DO_FLAG_DEL((V),RSB_FLAG_UPPER),RSB_DO_FLAG_ADD((V),RSB_FLAG_LOWER); \
else \
if(RSB_DO_FLAG_HAS((V),RSB_FLAG_LOWER)) \
	RSB_DO_FLAG_ADD((V),RSB_FLAG_UPPER),RSB_DO_FLAG_DEL((V),RSB_FLAG_LOWER); \
}
#define RSB__ERR_FLAG_DEL(V,F)	RSB_DO_FLAG_DEL((V),(F))
#define RSB_PERR_GOTO(LABEL,...) {RSB_ERROR(__VA_ARGS__);goto LABEL;}
#define RSB_SERR_GOTO(LABEL)     {goto LABEL;}

#define RSB_SYMMETRY_STRING(FLAGS) (RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_HERMITIAN)?"hermitian":(RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_SYMMETRIC)?"symmetric":"general"))

typedef float rsb_float_t;
#define RSB_FLOAT_ONE 1.0f 

#define RSB_RECURSION_MIN_DIM (2)
#define RSB_RECURSION_MIN_NNZ (4)

/*!
 An integer type for thread indices.
 */
typedef int rsb_thread_t;

/*! \internal  */
typedef rsb_flags_t rsb_order_t;



/** \internal \todo:OBSOLETE, REMOVE */ 
#define RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING	 	0x008000

/** \internal \todo:obsolete: FIXME  */ 
/*#define RSB_FLAG_BLOCK_ASYMMETRIC_Z_SORTED	 	0x010000*/
/** \internal \todo:temporary fix: FIXME  */ 
#define RSB_FLAG_FIX_FOR_BINARY_LOADED_MATRIX		 	0x010000

/** \internal \todo:EXPERIMENTAL */ 
#define RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR	 	0x020000

/** if set, the matrix will be partitioned with a block size chosen automatically */ 
#define RSB_FLAG_AUTO_BLOCKING				0x80000000	/* FIXME: obsolete, unsupported */


/** if set, will decide between column or row major for each (leaf) matrix (NEW: EXPERIMENTAL) */
/*#define RSB_FLAG_WANT_AUTO_MAJOR_ORDER 			0x200000	*/	/* Unsupported */


#if 0
/** if set, the matrix .. (remember to unindent it if active) */
 #define RSB_FLAG_WANT_RECURSIVELY_NON_UNIFORM_AUTO_BLOCKING 0x000200	/* experimental, but works well */
#endif /* 0 */

/** if set, the matrix will take possession of partitioning arrays p_r and p_c on input. if unset, a copy will be made	*/ 
#define RSB_FLAG_OWN_PARTITIONING_ARRAYS		0x000080	/*  */

/** if set, the blocks will be linked in some way */
/* FIXME: delete this flag */
/*#define RSB_FLAG_WANT_LINKED_STORAGE 			0x000400*/

/** if set, operating routines will check input more aggressively (may break operation)  */
#define RSB_FLAG_SHOULD_DEBUG 				0x000800

/** if set, the block partitioning will be fixed but VBR or LBR (Unsupported)	*/
#define RSB_FLAG_WANT_FIXED_BLOCKING_VBR	 	0x001000

/** if set, will mark a leaf matrix */
#define RSB_FLAG_NON_ROOT_MATRIX	 	0x100000

/** if set, the matrix will be prevent from being subdivided too much (OUTLAWED) */
/*#define RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES 		0x4000000*/

/**
 * if set, the blocks will cycle column after column.
 * if RSB_FLAG_WANT_BCSS_STORAGE is also set, the matrix storage format will be BCSC.
 * otherwise it will be VBC.
 * */
/* see rsb.h*/
/*#define RSB_FLAG_WANT_COLUMN_MAJOR_ORDER 		0x4000000*/


/** if set, the code will sort the input 			*/
#define RSB_FLAG_SORT_INPUT			0x2000000	/* FIXME: delete this flag */


/** a parameter to determine if a matrix is really 'small' or not (FIXME) */
#define RSB_EXPERIMENTAL_MIN_LEAF_ELEMENTS 		1024

/*#define RSB_FLAG_RECURSIVE_SHRINK_BOUNDING_BOX		0x40000000*/

#if 0
/* only flags left :  (remember to unindent it if active)*/
 #define RSB_FLAG_ALLOW_PARALLEL_OPERATION		0x40000000		/* NEW : UNUSED */
#endif /* 0 */

#if 0
#define RSB_FLAG_DEFAULT		 		(RSB_FLAG_DISCARD_ZEROS  /*| RSB_FLAG_WANT_BCSS_STORAGE*/ /* | RSB_FLAG_SORT_INPUT*/)
#endif /* 0 */


/**
 * \brief It is an internal structure, so beware, you should not use it.
 * \internal
 *
 * This structure will be used for keeping information about matrix partitioning.
 * It should be used primarily during matrix building, when the matrix arrays are 
 * not all allocated.
 *
 * It is an internal structure, so beware, you should not use it.
 * */
struct rsb_mtx_partitioning_info_t
{
	rsb_blk_idx_t M_b, K_b;		/**< just as in rsb_mtx_t */
	rsb_blk_idx_t br, bc;			/**< block row and column size (only if BCSR) (NEW) */
	rsb_coo_idx_t *rpntr,*cpntr;		/**< just as in rsb_mtx_t */
	
	rsb_coo_idx_t nr,nc;			/**< just as in rsb_mtx_t */
	rsb_submatrix_idx_t should_subdivide_levels;		/**< for recursive partitioning (EXPERIMENTAL) */
};


typedef signed   long rsb_long_t;

#define	RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT__VAL 0x10
#define RSB_OP_FLAG_WANT_SERIAL__VAL 0x2
#define RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT_SERIAL_VAL (RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT__VAL+RSB_OP_FLAG_WANT_SERIAL__VAL) 
/*!
 * \ingroup gr_internals
 * \brief An internal, helper enumeration.
 * \internal
 */
enum rsb_op_flags_t { 	RSB_OP_FLAG_DEFAULT=0x1, /* normal operation */
       			RSB_OP_FLAG_WANT_SERIAL=RSB_OP_FLAG_WANT_SERIAL__VAL,
		       	RSB_OP_FLAG_MAY_PARALLEL=0x4,
			RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE=0x5, /* will process only diagonal blocks */
			RSB_OP_FLAG_FAKE_LOCK=0x6, /* will perform operations with no locking (thus giving incorrect results) just to determine lock overhead */
       			RSB_OP_FLAG_WANT_PARALLEL_SORT=0x7,
       			RSB_OP_FLAG_WANT_SERIAL_SORT=0x8,
       			RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT=RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT__VAL,
       			RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT_SERIAL=RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT_SERIAL_VAL,
       			RSB_OP_FLAG_WANT_TRACE_PLOT=0x10
			};

struct rsb_optrace_t {
	/*! Submatrices operation time */
	rsb_time_t t0, t1;
	rsb_int_t th_id;
};

#define RSB_BLOCK_ROWMAJOR_ADDRESS(P,LDP,NR,NC,R,C,ES) \
	((rsb_char_t*)P)+((size_t)(ES))*((LDP)*(R)+(C))
#define RSB_BLOCK_COLMAJOR_ADDRESS(P,LDP,NR,NC,R,C,ES) \
	((rsb_char_t*)P)+((size_t)(ES))*((LDP)*(C)+(R))

#define RSB_BLOCK_X_MAJOR_REFERENCE(A,LDP,R,C,ONEIFISCOLMAJOR) \
	A[(ONEIFISCOLMAJOR)?((LDP)*(C)+(R)):((LDP)*(R)+(C))]

#define RSB_SOME_ERROR(ERRVAL) ((ERRVAL)!=RSB_ERR_NO_ERROR)

#define RSB_USE_OMP_SET_NUM_THREADS 0

#if RSB_USE_OMP_SET_NUM_THREADS
#define rsb_set_num_threads(RNT) omp_set_num_threads(RNT)
#define rsb_get_num_threads()    omp_get_num_threads()
#else
#define rsb_set_num_threads(RNT) rsb__set_num_threads(RNT)
#define rsb_get_num_threads()    rsb__set_num_threads(RSB_THREADS_GET)
#endif
#define RSB_DO_THREADS_PUSH(RNT)	{if((RNT)>0)rsb_set_num_threads(RNT); /* push */}
#define RSB_DO_THREADS_POP(RNT,ORNT)	{if((RNT)>0)rsb_set_num_threads(ORNT); /* pop */}

#if defined(RSB_WANT_RSB_NUM_THREADS) && (RSB_WANT_RSB_NUM_THREADS>1) /* Notice this branch is disabled (RSB_WANT_RSB_NUM_THREADS can't be >1), but RSB_NUM_THREADS is effective at rsb_lib_init() level */
#define RSB_NUM_THREADS_DECL	const char * rnt_str = getenv("RSB_NUM_THREADS"); rsb_int_t ornt = rsb_get_num_threads(), rnt = (rnt_str? rsb__util_atoi(rnt_str) :0);
#define RSB_NUM_THREADS_PUSH	{RSB_DO_THREADS_PUSH(rnt); /* push */}
#define RSB_NUM_THREADS_POP	{RSB_DO_THREADS_POP(rnt,ornt); /* pop */}
#else /* defined(RSB_WANT_RSB_NUM_THREADS) && (RSB_WANT_RSB_NUM_THREADS>0) */
#define RSB_NUM_THREADS_DECL
#define RSB_NUM_THREADS_PUSH
#define RSB_NUM_THREADS_POP
#endif /* defined(RSB_WANT_RSB_NUM_THREADS) && (RSB_WANT_RSB_NUM_THREADS>0) */

#define RSB__ENOUGH_NNZ_FOR_PARALLEL_SORT (RSB_MIN_THREAD_SORT_NNZ*rsb_get_num_threads())

#if defined(RSB_BLAS_WANT_EXPERIMENTAL_TUNING)
#define RSB_SPB_THREADS_PUSH	{RSB_DO_THREADS_PUSH(rnt); /* push */}
#define RSB_SPB_THREADS_POP	{RSB_DO_THREADS_PUSH(rnt,ornt); /* pop */}
#else /* defined(RSB_BLAS_WANT_EXPERIMENTAL_TUNING) */
#define RSB_SPB_THREADS_PUSH
#define RSB_SPB_THREADS_POP
#endif /* defined(RSB_BLAS_WANT_EXPERIMENTAL_TUNING) */

#define RSB_SPB_THREADS_DEFAULT 0
#define RSB_SPB_THREADS_AUTO -1
#define RSB_SPB_THR_STR_AUTO -2
#define RSB_SPB_THR_STR_AUTO_NEXTOP -3 /* TODO: need to diversify in thr.-only vs str.+thr. tuning */

/* #define RSB_PRINT_THREAD_STATS RSB_STDOUT("rsb_want_threads / rsb_g_threads / omp_get_max_threads / omp_get_num_threads / omp_get_thread_limit: %d / %d / %d / %d / %d\n",rsb_global_session_handle.rsb_want_threads,rsb_global_session_handle.rsb_g_threads,omp_get_max_threads(),omp_get_num_threads(),omp_get_thread_limit()); */
#define RSB_PRINT_THREAD_STATS RSB_STDOUT("rsb_want_threads / rsb_g_threads / omp_get_max_threads / omp_get_num_threads / omp_get_thread_limit: %d / %d / %d / %d\n",(int)rsb_global_session_handle.rsb_want_threads,(int)rsb_global_session_handle.rsb_g_threads,(int)omp_get_max_threads(),(int)omp_get_num_threads());

#define RSB_ERRMSG_NOSTREAMS "streams usage configured out."
#define RSB_ERRMSG_BADFORMAT "submatrix format unrecognized."
#define RSB_ERRMSG_NOTMTXMKT "not a Matrix Market format matrix"
#define RSB_ERRMSG_FILEOPENP "problems opening"
#define RSB_ERRMSG_PROIFAMM "problems reading or interpreting file as Matrix Market"
#define RSB_ERRMSG_FILEOPENPGZ "problems opening gzipped"
#define RSB_ERRMSG_TMXMKTBANNER "Could not process Matrix Market banner"
#define RSB_ERRMSG_BADCOO "bad input coo elements"
#define RSB_INFOMSG_SAK "is a swiss army knife for testing the library functionality and performance"

#define RSB_STDOUT_MATRIX_ESSENTIALS(M,MN,TN) RSB_STDOUT("%s\t%c\t%c\t%zd\t%zd\t%zd\t%zd",(const rsb_char_t*)rsb__basename(MN),rsb__do_get_symmetry_char(M),RSB_TRANSPOSITION_AS_CHAR(transA),(rsb_printf_int_t)(TN),(rsb_printf_int_t)(M)->nr,(rsb_printf_int_t)(M)->nc,(rsb_printf_int_t)(M)->nnz)
#define RSB_FPRINTF_MATRIX_ESSENTIALS(FD,M,MN,TN) RSB_FPRINTF(FD,"%s\t%c\t%c\t%zd\t%zd\t%zd\t%zd",(const rsb_char_t*)rsb__basename(MN),rsb__do_get_symmetry_char(M),RSB_TRANSPOSITION_AS_CHAR(transA),(rsb_printf_int_t)(TN),(rsb_printf_int_t)(M)->nr,(rsb_printf_int_t)(M)->nc,(rsb_printf_int_t)(M)->nnz)
#define RSB_FPINV(FPV) (1.0/(FPV))
#define RSB_MILLION_I 1000000
#define RSB_MILLION_F 1000000.0
#define RSB_CLEARTERM_STRING "\x1B\x4D"
#define RSB_MARKER_INT_VALUE (RSB_MAX_VALUE_FOR_TYPE(rsb_int_t)-RSB_NNZ_BLK_MAX)
/*#define RSB_MAX_SHORTIDX_MATRIX_DIM (RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t)-RSB_NNZ_BLK_MAX)*/
#define RSB_MAX_SHORTIDX_MATRIX_DIM (RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t))
#define RSB_BENCH_PROG_OPTS \
	    {"nthreads",	required_argument, NULL, 0x6E} /* n */
#define RSB_MAX_ALLOCATABLE_MEMORY_CHUNK \
((size_t)((sizeof(void*)==sizeof(unsigned int))? RSB_MAX_VALUE_FOR_TYPE(unsigned int):RSB_MAX_VALUE_FOR_TYPE(size_t)))
#define RSB_DOES_TRANSPOSE(TRANSA) ((TRANSA)!=RSB_TRANSPOSITION_N)
#define RSB_DOES_NOT_TRANSPOSE(TRANSA) (!RSB_DOES_TRANSPOSE(TRANSA))
#define RSB_DOES_CONJUGATE(TRANSA) ((TRANSA)==RSB_TRANSPOSITION_C)
#define RSB_DOES_NOT_CONJUGATE(TRANSA) (!RSB_DOES_CONJUGATE(TRANSA))
#define RSB_ELSE_IF_TRANSPOSE(NR,NC,TRANSA) (RSB_DOES_TRANSPOSE((TRANSA))?(NC):(NR))

#define RSB_ALLOW_ZERO_DIM RSB_MIN_MATRIX_DIM == 0

#if defined(RSB_WANT_OMP_RECURSIVE_KERNELS) && (RSB_WANT_OMP_RECURSIVE_KERNELS>0)
#define RSB_NT rsb_global_session_handle.rsb_g_threads
#define RSB_NTC num_threads(RSB_NT)
#define RSB__GET_MAX_THREADS() RSB_MIN(omp_get_max_threads(),RSB_NT) /* for use after rsb_lib_init */
#else
#define RSB_NT
#define RSB_NTC
#define RSB__GET_MAX_THREADS() RSB_NT /* for use after rsb_lib_init */
#endif
#define RSB_STORE_IDXSA 1

#define RSB_ASSIGN_IF_SP(DSTV,SRCP) 	\
	if ( (SRCP) != NULL )		\
		(DSTV) = *(SRCP);		/* FIXME: move this declaration elsewhere */

#define RSB_ASSIGN_IF_DP(DSTP,SRCV) 	\
	if ( (DSTP) != NULL )		\
		*(DSTP) = (SRCV);		/* FIXME: move this declaration elsewhere */

#define RSB_ASSIGN_IF(DSTP,SRCV) RSB_ASSIGN_IF_DP(DSTP,SRCV)

#define RSB_RSBENCH_MAX_MTXFILES 256
#define RSB_FAF_CHKDUP 0x01
#define RSB_FAF_CHKFNP 0x02
#define RSB_FAF_VRBSRC 0x04
#define RSB_FAF_CHKREC 0x08
#define RSB_FAF_VRBADD 0x10
#define RSB_FAF_CHKMTX 0x20
#define RSB_FAF_CHKEXS 0x40
#define RSB_FAF_CHKGSS 0x80
#define RSB_FAF_DEFAULTS RSB_FAF_CHKFNP|RSB_FAF_CHKDUP|RSB_FAF_CHKREC|RSB_FAF_VRBADD|RSB_FAF_CHKMTX|RSB_FAF_CHKEXS|RSB_FAF_CHKGSS 

#define RSB__ERR_ZERO_INF_NORM RSB_ERR_CAST(0x40000)
#define RSB__ERR_UNSUPPORTED_SYMM     RSB_ERR_UNSUPPORTED_TYPE
#define RSB__ERR_UNSUPPORTED_TRANSA   RSB_ERR_UNSUPPORTED_TYPE
#define RSB__ERR_UNSUPPORTED_IDX_TYPE RSB_ERR_UNSUPPORTED_FEATURE
#define RSB__ERR_UNSUPPORTED_DIAG     RSB_ERR_UNSUPPORTED_FEATURE
#define RSB__ERR_NO_SYM_SPSV          RSB_ERR_BADARGS
#define RSB__ERR_CANTUPDATE_DIAGI     RSB_ERR_CAST(0x80000)

#define RSB_MARF_NOFLAGS	(0x00000000)			/*!< */
#define RSB_MARF_EPS_T	(0x00000200)			/*!< #rsb_marf_t Flag value for for limiting the rsb__dump_postscript_recursion_from_mtx_t plot to the submatrices times pointed by *otv and info from *pv. */
#define RSB_MARF_EPS_O	(0x00000400)			/*!< #rsb_marf_t Combine this with RSB_MARF_EPS_L to get operands as well. */
#define RSB_MARF_EPS_NO_TEXT	(0x00000800)			/*!< #rsb_marf_t No text in EPS. */
#define RSB_MARF_LATEX_RECURSION	(0x00001000)			/*!< #rsb_marf_t Matrix recursion as LaTeX. */

rsb_err_t rsb__adddir(rsb_char_t ** filenameap, rsb_int_t * filenamenp, const rsb_char_t * matrixpath, rsb_flags_t faflags);
void rsb__setenv(const rsb_char_t * var_val);
rsb_real_t rsb__getenv_real_t(const char*envv, const rsb_real_t altv);
rsb_int_t rsb__getenv_int_t(const char*envv, const rsb_int_t altv);
const rsb_char_t * rsb__getenv_str(const char*envv, const rsb_char_t* altv);
rsb_char_t rsb__getenv_char(const char *envv, const rsb_char_t altv);

#ifdef RSB_HAVE_ASSERT_H 
#ifdef RSB_USE_ASSERT
/* ok, no extra action needed */
#else /* RSB_USE_ASSERT */
/* according to POSIX.1-2001, C89, C99, NDEBUG will cause assert to generate no code.  */
#define NDEBUG 1
#endif /* RSB_USE_ASSERT */
#include <assert.h>	/* the assert() macro */
#endif /* RSB_HAVE_ASSERT_H */

#include "rsb.h"		/* public API specification */
#include "rsb_init.h"		/* initialization functions */
#include "rsb_rec.h"		/* recursion handling functions */
#include "rsb_permute.h"	/* permutation functions */
#include "rsb_strmif.h"		
#include "rsb_srt.h"		/* sorting functions */
#include "rsb_mergesort.h"	/* sorting functions */
#include "rsb_merge.h"		/* merging functions */
#include "rsb_srtp.h"		/* parallel sorting functions */
#include "rsb_prec.h"		/* toy preconditioning */
#include "rsb_msort_up.h"	/* sorting functions, adapted from PSBLAS */
#include "rsb_unroll.h"		/* computational kernels */
#include "rsb_is.h"			/* coordinate handling functions */
#include "rsb_src.h"		/* search functions */
#include "rsb_clone.h"		/* clone functions */
#include "rsb_err.h"		/* error handling functions */
#include "rsb_internals.h"		/* */
#include "rsb_do.h"		/* */
#include "rsb_mio.h"			/* I/O functions */
#include "rsb_get.h"		/* matrix getter functions */
#include "rsb_set.h"		/* matrix setter functions */
#include "rsb_dump.h"		/* matrix dumping functions */
#include "rsb_coo.h"		/* coordinate handling functions */
#include "rsb_csr.h"		/* csr handling functions */
#include "rsb_blas_stuff.h"		/* BLAS like stuff */
#include "rsb_op.h"			/* */
#include "rsb_bio.h"		/* Binary Matrix I/O */
#include "rsb_asm.h"		/* Matrix assembly functions */
#include "rsb_coo_check.h"	/* */
#include "rsb_coo_symm.h"		/* */
#include "rsb_idx.h"		/* index manipulation */
/* #include "rsb_ftn.h"*/		/* fortran interface functions (obsolete) */
#include "rsb_libspblas_handle.h"	/*  */
#include "rsb_render.h"		/* matrix as pixmap rendering functions */
#include "rsb_eps.h"		/* matrix as (encapsulated) postscript rendering functions */
#include "rsb_gen.h"		/* matrix generating functions */
#include "rsb_sys.h"		/* system related functions */
#include "rsb_mbw.h"		/* memory benchmark related functions */
#include "rsb_limiter.h"	/*  */
#include "rsb_fpb.h"		/* floating point benchmark related functions */
#include "rsb_garbage.h"		/* misc helpers routines */
#include "rsb_pcnt.h"		/* performance counters code */
#include "rsb_perf.h"		/* performance info gathering code */
#include "rsb_pr.h"		/* performance reporting */
#include "rsb_util.h"		/* sorting and computational stuff */
#include "rsb_spmv.h"		/* sparse matrix-vector multiplication */
#include "rsb_swt.h"		/* switching format functions */
#include "rsb_lock.h"		/* */
#include "rsb_partition.h"	/* custom partitioning stuff (OBSOLETE) */
#ifndef __cplusplus
#include "rsb_krnl.h"		/* kernels rsb_krnlers */
#include "rsb_krnl_vb.h"	/* vb specific functions */
/* #include "libspblas_tests.h" */	/*  */
#include "rsb_test_accuracy.h"	/* accuracy testing functions */
#include "rsb_krnl_bcss.h"	/* bcss specific functions */
#include "rsb_krnl_bcoo_spmv_u.h"	/* bcoo specific functions */
#include "rsb_bench.h"		/* performance info gathering code (OBSOLETE) */
#include "rsb_spgemm.h"		/* sparse matrices multiplication */
#include "rsb_spgemm_csr.h"	/* sparse matrices multiplication */
#endif  /* __cplusplus */
#include "rsb_spsum_misc.h"	/* sum of Sparse Matrices */
#include "rsb_spsum.h"		/* Sum of Sparse Matrices */
#include "rsb_spsv.h"		/* */
#include "rsb_lbl.h"		/* OBSOLETE */
/* #include "rsb_experiments.h" */	/* experiments (obsolete) */
#include "rsb_coo2rec.h"		/* */
#include "rsb_rec2coo.h"		/* */
#include "rsb_rec2csr.h"		/* */
#include "rsb_csr2coo.h"		/* */
#include "rsb_cpmv.h"		/* */
#include "rsb_tune.h"		/* */

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_COMMON_H_INCLUDED */
/* @endcond */
