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
#ifndef RSB_RSB_STRUCT_H_INCLUDED
#define RSB_RSB_STRUCT_H_INCLUDED

/* @cond INNERDOC  */

/*! 
 * A type for the size (in bytes) of memory areas.  */
typedef size_t rsb_size_t;

/*!
 * A typedef for printing non negative integers.
 */
typedef rsb_size_t rsb_printf_int_t;

/*!
 * The inner datatype used for the bitmap structure.
 * */
typedef signed int rsb_bitmap_data_t;

/*! A type for specifying a sparse matrix format. */
typedef rsb_flags_t rsb_fmt_t;

/*!  A floating point numerical type for performance (MFLOPS) measurements.  */
typedef rsb_real_t rsb_perf_t;

/*!  A floating point numerical type for fillin measurements (obsolete).  */
typedef rsb_real_t rsb_fillin_t;

/*!
 An integer type for submatrix indices.
 */
typedef int rsb_submatrix_idx_t;

#define RSB_MAX_SUBM_COUNT (RSB_MAX_VALUE_FOR_TYPE(rsb_submatrix_idx_t))
#define RSB_SUBM_IDX_MARKER (RSB_MAX_SUBM_COUNT)

/*!
 \internal
 An integer type which by definition should not overflow, in most cases of interest.
 */
#ifdef RSB_WANT_LONG_IDX_TYPE 
typedef RSB_WANT_LONG_IDX_TYPE rsb_non_overflowing_t;
#else /* RSB_WANT_LONG_IDX_TYPE */
typedef int                    rsb_non_overflowing_t;
#endif /* RSB_WANT_LONG_IDX_TYPE */

#define RSB_BLK_MUL_OVERFLOW(R,C) RSB_MUL_OVERFLOW((R),(C),rsb_blk_idx_t,rsb_non_overflowing_t)
#define RSB_COO_MUL_OVERFLOW(R,C) RSB_MUL_OVERFLOW((R),(C),rsb_coo_idx_t,rsb_non_overflowing_t)
#define RSB_NNZ_MUL_OVERFLOW(R,C) RSB_MUL_OVERFLOW(R,C,rsb_nnz_idx_t,rsb_non_overflowing_t)

/*!
 A type for byte strings.
 */
typedef unsigned char rsb_byte_t;


/* @cond INNERDOC  */
/*!
 \name Macros for overflow detection in common (INTERNALS) operations. 

 They are tricky because should serve both signed and unsigned typedefs.
 The following macros should be handled with care.
 */
#define RSB_INDEX_OF_SAFE_EXTRA 2 	/*< this is the value that could be added with no overflow to indices values */
#define RSB_ADD_OVERFLOW(R,C,T) ((int)((T)(((T)(R))+((T)(C)))<((T)(R))) || (int)((T)(((T)(R))+((T)(C)))<((T)(C))))
#define RSB_MUL_OVERFLOW(R,C,T,H) ( ((C) == 0 || (R) == 0 ) ? 0 : ( (  R > C && RSB_MAX_VALUE_FOR_TYPE(T) / C < R ) || (RSB_MAX_VALUE_FOR_TYPE(T) / R < C ) ) )
#define RSB_BLK_ADD_OVERFLOW(R,C) RSB_ADD_OVERFLOW((R),(C),rsb_blk_idx_t)
#define RSB_COO_ADD_OVERFLOW(R,C) RSB_ADD_OVERFLOW((R),(C),rsb_coo_idx_t)
#define RSB_NNZ_ADD_OVERFLOW(R,C) RSB_ADD_OVERFLOW((R),(C),rsb_nnz_idx_t)

/*!
 Macros to get indices types liminal values (which we often use as markers).
*/
#define RSB_IS_UNSIGNED(T) (!RSB_IS_SIGNED(T))
#define RSB_PROBABLY_SAME_TYPES(T1,T2) ((RSB_IS_SIGNED(T1)==RSB_IS_SIGNED(T2)) && sizeof(T1)==sizeof(T2))
#define RSB_MIN_SIGNED(T) (-1 - RSB_MAX_SIGNED(T))

#define RSB_IS_VALUE_MORE_THAN_HALF_BITS_LONG(V,T)	(((V)>>RSB_COO_HALF_BITS_SIZE)>0) /*< */
#define RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(V) RSB_IS_VALUE_MORE_THAN_HALF_BITS_LONG(V,rsb_coo_idx_t) /*< */

#define RSB_MIN(X,Y) ((X)<(Y)?(X):(Y))		/*!< quick macro for minimum */
#define RSB_MAX(X,Y) ((X)<(Y)?(Y):(X))		/*!< quick macro for maximum */
#define RSB_ABS(X) ((X)<0?-(X):(X))		/*!< quick macro for abs()*/

#define RSB_FOUR 4			/*!< a constant with the number of quadrants */
/* @endcond */

/* @cond INNERDOC  */
#define RSB_REAL_ZERO 0.0		/*!< \internal internal */
#define RSB_TIME_ZERO RSB_REAL_ZERO 	/*!< \internal internal */	
#define RSB_BOOL_MAYBE	(-1) /*!< a reserved, "maybe" value for rsb_bool_t */
#define RSB_INVALID_FLAGS	(-1)		/*!< \internal internal */
#define RSB_INVALID_TRANS RSB_INVALID_FLAGS	/*!< \internal internal */
#define RSB_XOR(X,Y) 	(((X)!=0)^ ((Y)!=0))	/*!< \internal internal */
#define RSB_AND(X,Y) 	(((X)!=0)&&((Y)!=0))	/*!< \internal internal */
#define RSB_OR(X,Y) 	(((X)!=0)||((Y)!=0))	/*!< \internal internal */
#define RSB_NAND(X,Y)   (!RSB_AND(X,Y))		/*!< \internal internal */
/* @endcond */

/* @cond INNERDOC  */
#define RSB_DO_FLAG_XOR(V,F)	(V) ^= (F)	/*!< The flag variable \c V gets the logical XOR value with flag \c F. */
/* @endcond */

/* @cond INNERDOC  */
#define RSB_BOOL_XOR(X,Y) 	((X)^(Y)) /*!< A logical XOR for rsb_bool_t values. */
#define RSB_BOOL_OR(X,Y) 	((X)||(Y)) /*!< A logical OR for rsb_bool_t values. */
#define RSB_BOOL_AND(X,Y) 	((X)&&(Y)) /*!< A logical OR for rsb_bool_t values. */
#define RSB_BOOL_NOT(X) 	(!(X)) /*!< A logical NOT for rsb_bool_t values. */
#define RSB_BOOL_NAND(X,Y) 	RSB_BOOL_NOT(RSB_BOOL_AND(X,Y)) /*!< A logical NAND for rsb_bool_t values. */
#define RSB_BOOL_NOR(X,Y) 	RSB_BOOL_NOT(RSB_BOOL_OR(X,Y)) /*!< A logical NOR for rsb_bool_t values. */
/* @endcond */

/* @cond INNERDOC  */
#define RSB_WANT_DBC 0 /* want dense block count */
#define RSB_WANT_SSC 0 /* want spare (extra) sparse submatrices count */
/* @endcond */

/* @cond INNERDOC  */
/*!
 * \internal
 * \ingroup gr_internals
 * \brief An internal, helper structure (OBSOLETE).
 * \internal
 */
struct rsb_expected_info_t{
	/*! Expected fillin */
	/* FIXME : here should also be a map of expected fillin */
	rsb_fillin_t efillin;
};
/* @endcond */

/* @cond INNERDOC  */
/*!
 * \internal
 * \ingroup gr_internals
 * \brief An internal, helper structure (not for end users).
 * \internal
 */
struct rsb_translated_matrix_t
{
	struct rsb_mtx_t * mtxlp;
	rsb_submatrix_idx_t level;
	rsb_coo_idx_t	roff,coff;
	rsb_coo_idx_t	nr,nc;
};
/* @endcond */

/*!
 * \ingroup rsb_doc_matrix_assembly
 * \brief A structure for the RSB (Recursive Sparse Blocks) representation of sparse matrices.
 * \n
 * This is an opaque container for a recursive storage of COO/CSR submatrices. 
 * \n
 * The user is not supposed to manipulate this structure directly.
 * \n
 * This structure shall be only manipulated through the use of appropriate functions. 
 * \n
 * Knowledge of this structure is not required at all (in any case) use the library.
 * \see rsb_doc_matrix_assembly on how to instantiate/destroy this structure.
 * \see rsb_doc_matrix_operations for computational operations using it.
 *
 * \note: VBR and BCSR submatrices are not supported.
 */
struct rsb_mtx_t
{
	/*!
		values of matrix coefficients.
		array sized ( element_count == nnz * fillin ) * el_size (CSR,BCSR,VBR) 
	 */
	void * VA;

	/*!  bpntr[bri] points to the location of bindx of the first nonzero block entry of block row bri.
			   if the ith block row contains only zeros then bpntr[i]==bpntr[i+1] (VBR,BCSR,CSR) */
	rsb_nnz_idx_t *bpntr;

	/*!  bindx[bi] contains the block column index of the bi^th nonzero block (VBR,BCSR,CSR) */
	rsb_coo_idx_t	*bindx;	/* bindx[m->block_count] should be zero, for technical reasons (for the last 'virtual' block) */

	rsb_nnz_idx_t nnz;	/*! matrix rows, columns */
	rsb_coo_idx_t nr,nc;	/*! matrix (declared) nonzeros */
	rsb_flags_t flags; 	/*! structural flags, describing some optional features */
	rsb_blk_idx_t br, bc;	/*! block row and column size (only if BCSR) */
	rsb_type_t typecode; 	/*! as specified in the RSB_NUMERICAL_TYPE_* preprocessor symbols in types.h (See \ref matrix_type_symbols_section)	*/
	rsb_fmt_t matrix_storage; /*! as specified in the RSB_MATRIX_STORAGE_* preprocessor symbols in types.h 	*/

	/*!
		intptr[bi] points (logically: in terms of numerical elements count) to the location in VA of the (0,0) entry in the bi^th block entry (VBR).
		array sized 
	*/
	rsb_nnz_idx_t  *indptr;

	/*!  rpntr[bri] contains the row index of first row in the bri^th block row
	          ( row    partitioning indices : M_b +1 elements )  (CSR,BCSR,VBR)
	     note that rpntr[Mdim] could be more than m.
	*/
	rsb_coo_idx_t	*rpntr;

	/*!  cpntr[bcj] contains the column index of the first column in the bcj^th block column
	          ( column partitioning indices : K_b +1 elements ) (VBR) */
	rsb_coo_idx_t *cpntr;

	/*!  these are aliases for rpntr and cpntr for the major dimension (Mpntr) and minor one (mpntr) 
	 */
	rsb_coo_idx_t *mpntr,*Mpntr;

	/* int  *mpntr,*Mpntr;*/	/* aliases for rpntr and cpntr (M stays for major, m for minor) */
	
	/*! block row and column counts */
	rsb_blk_idx_t M_b, K_b;

	/*!  these are aliases for M_b and K_b for the major dimension (Mdim) and minor one (mdim) 
	 *  ifdef RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, the aliasing is swapped.
	 * */
	rsb_blk_idx_t Mdim,mdim;

#if RSB_WANT_DBC
	/*! The count of blocks (regardless their size) : <= nnz */
	rsb_nnz_idx_t block_count;
#endif

#if RSB_WANT_SSC
	/*! The spare (extra) sparse blocks count. */
	rsb_nnz_idx_t ssbc;
#endif

	/*! The overall number of elements el_size bytes each (>=nnz) */
	rsb_size_t element_count;
	
	/*! the size >= 1, in bytes, of the sparse matrix numerical elements type */
	rsb_size_t el_size;

	/*! Time needed for matrix structure analysis, during construction */
	rsb_time_t sat;

	/*! Time needed for elements insertion, during construction */
	rsb_time_t eit;

	/*! Time needed for sorting elements (if sorted), during construction */
	rsb_time_t est;

	/*! Performance estimation time */
	rsb_time_t pet;

	/*! Cooordinate cleanup time */
	rsb_time_t ect;

	/*! Coordinate partitioning time */
	rsb_time_t cpt;

	/*! Recursive sort time  */
	rsb_time_t rpt;

	/*! Total assembly time */
	rsb_time_t tat;

	/*! Submatrix pointers for recursion storage */
	struct rsb_mtx_t * sm[RSB_FOUR];

/* #if RSB_STORE_IDXSA */
	/*! Index storage amount. Temporarily here: FIXME. */
	rsb_size_t idxsa;
	/*
#else */
	/*! A structure with expectation info during construction (FIXME: this member is obsolete and will be deleted soon) */
	/* struct rsb_expected_info_t einfo; */
/* #endif */

	/*! A pointer to an array of leaf submatrices pointers (only valid on root) */
	struct rsb_translated_matrix_t * all_leaf_matrices;

	/*! The number of leaf submatrices pointers in all_leaf_matrices (only valid on root) */
	rsb_submatrix_idx_t all_leaf_matrices_n;

	/*! In a recursive representation, the offset of the submatrix (respectively, rows and columns)  */
	rsb_coo_idx_t	roff,coff;

	/*! In a recursive representation, with the RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS flag, the offset of these data arrays from the beginning of the global ones  */
	rsb_nnz_idx_t	nzoff;

	/*! In a recursive representation, broff (bcoff) is the offset of the submatrix first non empty row (column) with respect to the matrix.  */
	rsb_coo_idx_t	broff,bcoff;

	/*! In a recursive representation, bm (bk) is the last non empty row (column) in the submatrix.  */
	rsb_coo_idx_t bm,bk;
};

/*!
 * Macros for printing out summary info about a matrix.
 * Accept a valid \ref rsb_mtx_t  pointer as an argument.
 *
 * Usage example:
 * \code
 * printf(RSB_PRINTF_MTX_SUMMARY_ARGS(mtxAp));
 * \endcode
 */
#define RSB_PRINTF_FLAGS_FMT_STR  		\
					"%s"	\
					"%s"	\
					"%s"	\
					"%s"	\
					"%s"	
#define RSB_NNZ_PER_ROW_AS(MTXAP,TYPE)	((MTXAP)->nr ? ((TYPE)(MTXAP)->nnz)/(MTXAP)->nr : 0)
#define RSB_PRINTF_FLAGS_ARGS(FLAGS)  		\
	RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_UPPER)?"U":"",		\
	RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_LOWER)?"L":"",		\
	RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_TRIANGULAR)?"T":"",	\
	RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_SYMMETRIC)?"S":"",	\
	RSB_DO_FLAG_HAS((FLAGS),RSB_FLAG_HERMITIAN)?"H":""
#define RSB_PRINTF_MTX_SUMMARIZED_ARGS(PRE,MTXAP,POST)  \
			"%s(%zd x %zd)[%p]{%c} @ (%zd(%zd..%zd),%zd(%zd..%zd)) (%zd nnz, %.2lg nnz/r) flags 0x%x (coo:%d, csr:%d, hw:%d, ic:%d, fi:%d), storage: %zx, subm: %zd, symflags:'"\
					RSB_PRINTF_FLAGS_FMT_STR  \
					"'%s"	\
					, \
					(PRE), \
				(rsb_printf_int_t)(MTXAP)->nr, (rsb_printf_int_t)(MTXAP)->nc, (const void*)(MTXAP),				\
				(char)(MTXAP)->typecode,						\
				(rsb_printf_int_t)(MTXAP)->roff,						\
				(rsb_printf_int_t)(MTXAP)->broff,						\
				(rsb_printf_int_t)((MTXAP)->roff+(MTXAP)->bm),				\
				(rsb_printf_int_t)(MTXAP)->coff,						\
				(rsb_printf_int_t)(MTXAP)->bcoff,						\
				(rsb_printf_int_t)((MTXAP)->coff+(MTXAP)->bk),				\
			       	(rsb_printf_int_t)(MTXAP)->nnz,						\
			       	RSB_NNZ_PER_ROW_AS(MTXAP,double),					\
			       	(unsigned int)((MTXAP)->flags),						\
				(int)(RSB_DO_FLAG_HAS((MTXAP)->flags,RSB_FLAG_WANT_COO_STORAGE)),	\
				(int)(RSB_DO_FLAG_HAS((MTXAP)->flags,RSB_FLAG_WANT_BCSS_STORAGE)),	\
				(int)(RSB_DO_FLAG_HAS((MTXAP)->flags,RSB_FLAG_USE_HALFWORD_INDICES)),	\
				(int)(RSB_DO_FLAG_HAS((MTXAP)->flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS)),\
				(int)(RSB_DO_FLAG_HAS((MTXAP)->flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE)), 	\
				(rsb_printf_int_t)((MTXAP)->matrix_storage),					\
				(rsb_printf_int_t)((MTXAP)->all_leaf_matrices_n),				\
				RSB_PRINTF_FLAGS_ARGS((MTXAP)->flags),					\
				(POST)

#define RSB_PRINTF_MTX_SUMMARY_ARGS(MTXAP)  RSB_PRINTF_MTX_SUMMARIZED_ARGS("",MTXAP,"")

#define RSB_PRINTF_MATRIX_AT_SUMMARY_ARGS(MTXAP)  \
			"%ld x %ld, type %c, %ld nnz, %.2lg nnz/r, %ld subms, %d lsubms, %2.4lf bpnz"\
					, 								\
				(long int)(MTXAP)->nr, (long int)(MTXAP)->nc, 						\
				(MTXAP)->typecode,							\
			       	(long int)(MTXAP)->nnz,								\
			       	RSB_NNZ_PER_ROW_AS(MTXAP,double),					\
				rsb__submatrices(MTXAP),							\
				(MTXAP)->all_leaf_matrices_n,						\
				((double)rsb__get_index_storage_amount(MTXAP)) / ((MTXAP)->nnz)

#define RSB_PRINTF_MATRIX_BOUNDS_SUMMARY_ARGS(MTXAP)  \
			"(nr=%ld x nc=%ld, nnz=%ld)[%p]{type=%c} @ (nzoff=%d, roff=%d,broff=%d,bm=%d, coff=%d,bcoff=%d,bk=%d) " \
					, \
				(long int)(MTXAP)->nr, (long int)(MTXAP)->nc, (long int)(MTXAP)->nnz, (const void*)(MTXAP),				\
				(MTXAP)->typecode,						\
				(MTXAP)->nzoff,						\
				(MTXAP)->roff,						\
				(MTXAP)->broff,						\
				(MTXAP)->bm,						\
				(MTXAP)->coff,						\
				(MTXAP)->bcoff,						\
				(MTXAP)->bk

#define RSB_WANT_COO_BEGIN 1

#if RSB_WANT_COO_BEGIN
#define RSB_MTX_SET_HBDF(MTXAP) (MTXAP)->RSB_MTX_BMF=RSB_MTX_BMV
#define RSB_MTX_HBDF(MTXAP) ((MTXAP)->RSB_MTX_BMF==RSB_MTX_BMV)
#define RSB_MTX_HBDFH(MTXAP) ((MTXAP)->RSB_MTX_BDF)
#define RSB_MTX_BDF nnz
#define RSB_MTX_BMF nr
#define RSB_MTX_BMV -1
#endif
/* RSB_WANT_COO_BEGIN */

#define RSB_MTX_TRANSPOSED_ROWS(MTX,TRANSA) RSB_ELSE_IF_TRANSPOSE((MTX)->nr,(MTX)->nc,(TRANSA))
#define RSB_MTX_TRANSPOSED_COLS(MTX,TRANSA) RSB_ELSE_IF_TRANSPOSE((MTX)->nc,(MTX)->nr,(TRANSA))
#define RSB_MTX_DIAG_SIZE(MTX) RSB_MIN( (MTX)->nc,(MTX)->nr )
#define RSB_MTX_DIAG_SIZE_BLK(MTX)  RSB_MTX_DIAG_SIZE_BLK(MTX) + RSB_NNZ_BLK_MAX
#define RSB_MTX_EFF_R(MTXAP) (MTXAP->bm-((MTXAP)->broff-(MTXAP)->roff)) /* effective rows count */
#define RSB_MTX_EFF_C(MTXAP) (MTXAP->bk-((MTXAP)->bcoff-(MTXAP)->coff)) /* effective columns count */
#define RSB_ANY_MTX_DIM_ZERO(MTXAP) ((MTXAP) && (((MTXAP)->nr==0)||(MTXAP)->nc==0))

#define RSB_NNZ_OF(MTXAP) ((MTXAP)?((MTXAP)->nnz):0)
/* @endcond */

#endif
