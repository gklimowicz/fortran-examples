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
 *  \file
 *  \brief  This file declares the user interface functions and data structures for the \librsb library (see \ref rsb_doc_rsb).
 *  \author Michele Martone
 * */
/*!
 \mainpage
 
 A sparse matrix library implementing the `Recursive Sparse Blocks' (\b RSB) matrix storage.

 
 This is the documentation for the application programming interface (API)
 of the \e `\librsb' library.
 \n
 In order to use \librsb, there is no need for the user to know the RSB
 layout and algorithms: this documentation should be sufficient.
 \n
 \librsb is dual-interfaced; it supports:
 a native (`RSB') interface (with identifiers prefixed by `rsb_' or `RSB_'),
 and a (mostly complete) Sparse BLAS interface, as a wrapper around the RSB interface.
 \n
 Many computationally intensive operations are implemented with thread
 parallelism, by using OpenMP.
 \n
 Thread parallelism can be turned off at configure time, if desired, or controlled
 at execution time.
 \n
 Many of the computational kernels source code files (mostly internals) were
 automatically generated.
 \n
 This user documentation concerns the end user API only.
 \n
 
 This library was born as \b research \b software.
 \n
 Over the years it has been stabilized and its \b test \b suite expanded considerably.
 \n
 There are projects making available \librsb in Octave, Python or Julia;
 you might want to check them out as well.

 \n
 For a first approach, see the \ref rsb_doc_examples  documentation
 section, or the \ref examples_section "quick start examples" section on this page.
 
 \n
 Information about the supported matrix types and matrix operations 
 resides in the \link rsb_types.h rsb_types.h \endlink file. 

 A C/C++ user can use the native API of RSB by including the \link rsb.h rsb.h \endlink header.
 \n
 The same interface is available in Fortran via the ISO C Binding interface, specified in \link rsb.F90 rsb.F90\endlink.
 \n
 An \b experimental \b  C++ API is available in (optional) \ref rsb.hpp.
 \n
 
 The C header file for \ref rsb_doc_sparse_blas  is \link blas_sparse.h blas_sparse.h\endlink (notice there is no native C++ equivalent).

 \author Michele Martone < michelemartone AT users DOT sourceforge DOT net >
 
 Contents of the README file :
 \verbinclude README

 \anchor examples_section 

 For a quick startup, consider the following three example programs.

 The first, using the RSB interface (for native C++, please see the optional rsblib/examples/example.cpp):
 \include examples/hello.c

 The second, using the Sparse BLAS interface:
 \include examples/hello-spblas.c

 For more, see the \ref rsb_doc_examples  section.

 */

/*!
 \defgroup rsb_doc_rsb The librsb library interface (rsb.h, optional ones rsb.hpp and rsb.F90)
 \brief
 \n
 The reference documentation of the \librsb library comes in both HTML and Unix man pages formats.
 \n
 The following sections/man pages are available: \ref rsb_doc_rsb ; \ref rsb_doc_sparse_blas ; \ref rsb_doc_examples.


 In general, users of this library are interested in high performance sparse matrix computations on cache based shared memory parallel computers.
 \n
 For this, \librsb offers a native C interface (here documented) and a Fortran one (in \ref rsb.F90, equivalent to the C declaration headers from \ref rsb.h), in addition to the Sparse BLAS one (both C and Fortran, documented).
 \n
 Please refer to optional <rsb.hpp> for the C++ API.

 Configuration, build, and installation instructions are contained in the \c README file distributed in the sources archive.

 <b> Typical C program structure </b>

 \li initialize \librsb with #rsb_lib_init()
 \li (in any order)
      allocate matrices (e.g.: with #rsb_mtx_alloc_from_coo_inplace() or others);
      do any computation with them (e.g.: #rsb_spmv(), #rsb_spsv() );
      converting matrices (e.g.: with #rsb_mtx_switch_to_coo() );
      freeing matrices (#rsb_mtx_free() )
 \li finalize \librsb with #rsb_lib_exit()


 <b> Important usage notes </b>

 <b> General C program structure </b>
 Before calling any \librsb function, a program is required to initialize \librsb's internal status.
 This is done by calling #rsb_lib_init() .
 Afterwards, any \librsb function can be safely used.
 When \librsb functions are not intended to be called anymore, a program may call #rsb_lib_exit() to free any resource.
 Then, #rsb_lib_init() should be called for further usage of \librsb.

 <b> Manipulating matrices and vectors </b>
 In order to use \librsb, the user is not required to use explicitly any of \librsb's data structures: their manipulation is to be performed by \librsb functions.
 Therefore, knowledge of \librsb's matrix type (\c rsb_mtx_t) is not necessary at all: this structure is intended to be used as an opaque container.

 On the contrary, arrays for numerical vectors (or more generally, dense matrices) are expected to be managed by the user: \librsb does not furnish any specific vector type.
 Computational functions treat dense vectors/matrices are simple arrays of a specified type; see the \ref rsb_doc_examples  .

 <b> Computational functions </b>
 This library can be configured at build time to support a custom subset of numerical types.
 To keep the programming interface compact, it has been decided to not replicate the computational functions to each numerical type.
 Instead, the type is expected to be specified by the user via a type flag. 
 For instance, matrix assembly functions (e.g.: #rsb_mtx_alloc_from_coo_const() ) accept a type information and keep it stored in the matrix structure.
 Therefore, computational functions (e.g.: #rsb_spmv() ) can fetch this information from their \c rsb_mtx_t operand, and treat accordingly the other parameters (e.g.: \a alphap, \a Xp, ...).
 Mixed type operations are currently not supported.


 <b> Memory management </b>

 Matrix structures (\c rsb_mtx_t) allocated by \librsb shall be freed only via #rsb_mtx_free() .

 <b> Benchmarking </b>

 If you want to benchmark this library, there are different possibilities:
 \include ./examples/benchex.sh

 <b> Tuning and Customization </b>
 
 There are different \c ./configure  options you may look at for tuning or customizing the library.


 \defgroup rsb_doc_examples	Example programs and code 
 \brief	Examples of usage of \librsb in C (<rsb.h> and <blas_sparse.h>), C++ (optional <rsb.hpp>), Fortran (<rsb.F90> and <rsb_blas_sparse.F90>).

 	The following fully working example programs illustrate correct ways of using the library.
 
 - First example in C, using <rsb.h>: \ref examples_hello_c.
 - First example in C, using <blas_sparse.h>: \ref examples_hello_spblas_c.
 - Autotuning example in C, using <rsb.h>: \ref examples_autotune_c.
 - I/O example in C, using <blas_sparse.h>: \ref examples_io_spblas_c.
 - Example transposing a matrix in C, using <rsb.h> in C: \ref examples_transpose_c.
 - Example showing the power method in C, using <rsb.h> in C: \ref examples_power_c.
 - Example in Fortran, using <rsb_blas_sparse.F90>: \ref examples_fortran_F90.
 - Example in Fortran, using <rsb.F90>: \ref examples_fortran_rsb_fi_F90.
 - Example in C, using <rsb.h>: \ref examples_backsolve_c.
 - Misc example snippets in C, using <rsb.h>: \ref examples_snippets_c.
 - Benchmark invocation from shell script: \ref examples_bench_sh.

Once installed \librsb, the script displayed here (\ref examples/make.sh) should be sufficient to build these examples:
 \include examples/make.sh

 \anchor  examples_hello_c
          examples/hello.c:
 \include examples/hello.c

 \anchor  examples_hello_spblas_c
          examples/hello-spblas.c:
 \include examples/hello-spblas.c

 \anchor  examples_autotune_c
          examples/autotune.c:
 \include examples/autotune.c

 \anchor  examples_io_spblas_c
          examples/io-spblas.c:
 \include examples/io-spblas.c

 \anchor  examples_transpose_c
          examples/transpose.c:
 \include examples/transpose.c

 \anchor  examples_power_c
          examples/power.c:
 \include examples/power.c

 \anchor  examples_fortran_F90
          examples/fortran.F90:
 \include examples/fortran.F90

 \anchor  examples_fortran_rsb_fi_F90
          examples/fortran_rsb_fi.F90:
 \include examples/fortran_rsb_fi.F90

 \anchor  examples_backsolve_c
          examples/backsolve.c:
 \include examples/backsolve.c

 \anchor  examples_snippets_c

 \anchor  examples_bench_sh
          examples/bench.sh:
 \include examples/bench.sh

 Most of the snippets in the documentation come from examples/snippets.c.
*/

/*!
 \defgroup rsb_doc_sparse_blas The Sparse BLAS interface to librsb (blas_sparse.h, rsb_blas_sparse.F90) 
 \brief
 	A Sparse BLAS interface (see http://www.netlib.org/blas/blast-forum/) to \librsb.  Level 1 (vector-vector operations) is supported in a basic way.  Level 2 (sparse matrix-dense vector operations) is supported fully.  Level 3 (sparse matrix-dense matrix operations) is supported as a wrapper around Level 2.

We also implement a number of useful extra functions as custom extensions, giving access to other \librsb functionality.

The usage pattern of this interface matches that of the Sparse BLAS standard, exception made for the necessity of initialization/finalization of \librsb.
The Sparse BLAS interface is also available for Fortran: see \ref rsb_blas_sparse.F90.


The user should be aware of the following:
\li Because this Sparse BLAS implementation is built around \librsb, initialization with #rsb_lib_init() and finalization with #rsb_lib_exit() is necessary. Inclusion of the \c rsb.h header is necessary.
\li \librsb gives users freedom of in/out arbitrarily BLAS types support at configure/build time. Hence, while all the interface functions are always included the Sparse BLAS header file, they may return an error code. Be sure of having configured correctly the library at configure time (and see the \ref blas_sparse.h header file for types configured in the current build).
\li According to the standard, the complex type functions for C accept scalar values by reference rather than by copy; equivalent functions for other types do not do so, so this may cause confusion. Be careful.
\li Error checking is weak; so for instance, passing a function the handle of a matrix of mismatching type will not be detected as an error, although it's incorrect.
\li According to the standard, VBR and BCSR styled constructors are supported, although these are interfaces for \librsb's own matrix representation.
\li Here functions for both Fortran and C are listed. The Fortran functions are declared and documented with the C notation. We may provide a better documentation in a future release.
\li Each identifier documented here suffixed by \c _  (e.g.: #blas_susdot_()) can be used from Fortran with the name stripped by that suffix (so in this case, \c blas_susdot).
We may provide a proper fix to this inconvenience in a subsequent release.
\li Each Fortran program using \librsb's Sparse BLAS Implementation shall \c use  modules \c blas_sparse  and \c rsb.
\li Also Fortran programs have to call #rsb_lib_init() and #rsb_lib_exit() e.g.:
\verbatim
       	USE blas_sparse             ! module implementing the Sparse BLAS on the top of librsb
       	USE rsb                     ! rsb module
	...
	INTEGER :: istat            ! integer variable
	...
       	istat = rsb_lib_init(RSB_NULL_INIT_OPTIONS) ! please note that this is not part of Sparse BLAS but it is needed by librsb
	if(istat.NE.0)STOP          ! a value different than zero signals an error
	...
	! code calling Sparse BLAS routines
	...
       	istat = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS) ! please note that this is not part of Sparse BLAS but it is needed by librsb
	if(istat.NE.0)STOP          ! a value different than zero signals an error
	...
\endverbatim
	\li For Fortran, more procedures exist, although they are not documented here. According to the Sparse BLAS (http://www.netlib.org/blas/blast-forum/), for almost each subroutine whose identifier prefixed with \c blas_X (with \c X being one of S,D,C,Z), a corresponding generic modern Fortran version exists.
	Please note how not all of the certain procedures identifier prefixes include the type character. 

	E.g.:
\code
      ! the following code ('d' stays for 'double precision'):
      CALL blas_duscr_begin(nr,nc,A,istat)
      CALL blas_ussp(A,blas_lower_symmetric,istat)
      CALL blas_duscr_insert_entries(A,nnz,VA,IA,JA,istat)
      CALL blas_duscr_end(A,istat)
      CALL blas_dusmv(transT,alpha,A,X,incX,B,incB,istat) 
      CALL blas_dusds(A,istat)
      ! is equivalent to:
      CALL duscr_begin(nr,nc,A,istat) ! here, 'd' must be retained for avoiding ambiguity
      CALL ussp(A,blas_lower_symmetric,istat)
      CALL uscr_insert_entries(A,nnz,VA,IA,JA,istat)
      CALL uscr_end(A,istat)
      CALL usmv(transT,alpha,A,X,incX,B,incB,istat) 
      CALL usds(A,istat)
\endcode
*/

#ifndef RSB_RSB_H_INCLUDED
#define RSB_RSB_H_INCLUDED

#if ( defined(__cplusplus) && (__cplusplus>=201103L) )
#include <cstdlib>	/* size_t */
#else  /* __cplusplus */
#include <stdlib.h>	/* size_t */
#endif /* __cplusplus */

/*! \internal
 NOTE: user programs should never include explicitly rsb_types.h.
 */
#ifndef RSB_WANT_NO_RSB_TYPES_H
#include "rsb_types.h"
#endif /* RSB_WANT_NO_RSB_TYPES_H */

#ifdef RSB_WANT_LONG_IDX_TYPE
/* The user wants a long type for indices. */
#if ( defined(__cplusplus) && (__cplusplus>=201103L) )
#include <cstdint>
#else  /* __cplusplus */
#include <stdint.h>
#endif /* __cplusplus */
#endif /* RSB_WANT_LONG_IDX_TYPE */

#if RSB_HAVE_ANY_COMPLEX_TYPE
#ifdef __cplusplus
#include <complex>	/* std::complex */
#else
#include <complex.h>	/* ISO C99 */
#endif /* __cplusplus */
#endif /* RSB_HAVE_ANY_COMPLEX_TYPE */

#ifdef __cplusplus
extern "C" {
#endif

/*!
 \name Type definitions
 \anchor definitions_section

 These are definitions of \librsb base types.
 */
/*!
 * The block arrays index type. 
 *
 * Could be an unsigned type.
 * Should not overflow when indexing matrix blocks by block coordinates.
 * */
typedef signed int rsb_blk_idx_t;

/*!
 * The coordinate arrays index type.
 *
 * Should not overflow when indexing matrix elements by coordinates.
 * Legal values when specifying a matrix size should be within #RSB_MIN_MATRIX_DIM and #RSB_MAX_MATRIX_DIM 
 * */
#ifdef RSB_WANT_LONG_IDX_TYPE 
typedef RSB_WANT_LONG_IDX_TYPE rsb_coo_idx_t;	
#else /* RSB_WANT_LONG_IDX_TYPE */
typedef signed int rsb_coo_idx_t;	
#endif /* RSB_WANT_LONG_IDX_TYPE */

/*!
 * The nnz counter index type.
 *
 * Should not overflow when indexing matrix elements.
 * On most common archs sizeof(long)>=sizeof(int).
 * Legal values when specifying a matrix size should be within #RSB_MIN_MATRIX_NNZ and #RSB_MAX_MATRIX_NNZ 
 * */
#ifdef RSB_WANT_LONG_IDX_TYPE 
typedef RSB_WANT_LONG_IDX_TYPE rsb_nnz_idx_t;	
#else /* RSB_WANT_LONG_IDX_TYPE */
typedef signed int rsb_nnz_idx_t;	
#endif /* RSB_WANT_LONG_IDX_TYPE */

/* We would like the following typedefs to be long, but
   they should be compatible with many int functions */

/*!
 A type for specifying matrix assembly or coordinate conversions option flags.
 Should be >= 4 bytes.
 See \ref flags_section for possible values.
 */
typedef signed int rsb_flags_t;

/*!
 A type for specifying numerical type codes (See \ref matrix_type_symbols_section for a list of valid values, corresponding to \ref matrix_supported_numerical_types_section).
 */
typedef char rsb_type_t;

/*!
 A type specific for error flags.
 Should be >= 4 bytes.

 A textual description of an error value may be obtained via #rsb_strerror_r() or #rsb_perror().
 */
typedef signed int rsb_err_t; /* note that an unsigned would break the RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG macros! */

/*!
 An integer type declaration for interface functions.
 Should always be 'int'.
 */
typedef signed    int rsb_int_t;		/*!< A signed integer type */

/*! A boolean type. */
typedef rsb_flags_t rsb_bool_t;

/*!
 * The type for specifying transposition (See \ref matrix_transposition_flags_section)
 */
typedef rsb_flags_t rsb_trans_t;

/*!  A floating point numerical type for certain functions e.g. \ref rsb_tune_spmm. Not to be used for sparse matrices! For that, see \ref matrix_type_symbols_section. */
typedef double rsb_real_t;

/*!
 A type for character strings.
 */
typedef char rsb_char_t;

/*!  A floating point numerical type for time measurements with #rsb_time().  */
typedef rsb_real_t rsb_time_t;

/*!
 \name Other constants

 Other constants for some typedefs.
 */
#define RSB_BOOL_TRUE	1	/*!< A "true"  value for #rsb_bool_t. */
#define RSB_BOOL_FALSE	0 /*!< A "false" value for #rsb_bool_t. */
#define RSB_DO_FLAG_ADD(V,F)	(V) |=  (F)	/*!< The flag variable \c V gets the logical OR value with flag \c F. */
#define RSB_DO_FLAG_DEL(V,F)	(V) &= ~(F)	/*!< The flag variable \c V gets the logical NAND value with flag \c F. */
#define RSB_DO_FLAG_FILTEROUT(V,F)	((V) & ~(F))	/*!< The flag variable \c V after logical NAND value with flag \c F. */
#define RSB_DO_FLAG_FILTERONLY(V,F)	((V) & (F))	/*!< The flag variable \c V after logical AND value with flag \c F. */
#define RSB_DO_FLAG_HAS(V,F)	((((V)&(F))==(F))?RSB_BOOL_TRUE:RSB_BOOL_FALSE)	 /*!< Presence check for flag \c F. */
#define RSB_DO_FLAG_HAS_INTERSECTION(V,F)	(((V)&(F))?RSB_BOOL_TRUE:RSB_BOOL_FALSE)	/*!< Presence check for flag \c F.*/

#define RSB_DEFAULT_ROW_BLOCKING 1	/*!< Reserved for future use. */
#define RSB_DEFAULT_COL_BLOCKING 1	/*!< Reserved for future use. */
#define RSB_DEFAULT_BLOCKING 1 /*!< A safe value for column blocking (reserved for future use). */

/*  Macros to get indices types liminal values.  */
#define RSB_IS_SIGNED(T)   (((T)0) > (((T)-1)))
#define RSB_MAX_UNSIGNED(T) ((T)-1)
#define RSB_CHAR_BIT 8	/* bits per byte; if not 8, librsb compilation should fail */
#define RSB_HALF_MAX_SIGNED(T) ((T)1 << (sizeof(T)*RSB_CHAR_BIT-2))
#define RSB_MAX_SIGNED(T) (RSB_HALF_MAX_SIGNED(T) - 1 + RSB_HALF_MAX_SIGNED(T))
#define RSB_MAX_VALUE_FOR_TYPE(T) (RSB_IS_SIGNED(T)?RSB_MAX_SIGNED(T):RSB_MAX_UNSIGNED(T))

#define RSB_MIN_MATRIX_DIM 0 /*!> Minimum allowed matrix dimension. */
#define RSB_MIN_MATRIX_NNZ 0 /*!> Minimum allowed matrix nonzeroes count. */
#define RSB_NNZ_BLK_MAX 255 /* Dense block maximal allowed size (still unused, for now internal) */
#define RSB_MAX_MATRIX_DIM (RSB_MAX_VALUE_FOR_TYPE(rsb_coo_idx_t)-RSB_NNZ_BLK_MAX-255) /*!> Maximum allowed matrix dimension. */
#define RSB_MAX_MATRIX_NNZ (RSB_MAX_VALUE_FOR_TYPE(rsb_nnz_idx_t)-RSB_NNZ_BLK_MAX) /*!> Maximum allowed matrix nonzeroes count. */
#define RSB_MARKER_COO_VALUE (RSB_MAX_MATRIX_DIM+1)		/* */
#define RSB_MARKER_NNZ_VALUE (RSB_MAX_MATRIX_NNZ+1)		/* */
#define RSB_INVALID_COO_IDX_VAL ((RSB_MARKER_COO_VALUE)+1)	/*< A value which is illegal for any #rsb_coo_idx_t variable. */
#define RSB_INVALID_NNZ_IDX_VAL ((RSB_MARKER_NNZ_VALUE)+1)	/*< A value which is illegal for any #rsb_nnz_idx_t variable. */

/*! \anchor rsb_mtx_t  struct rsb_mtx_t declaration is in a separate, internal include file */

/*!
 \ingroup rsb_doc_rsb
 \name Matrix assembly flags
 \anchor flags_section

 These are flags which could be combined to specify the assembly of sparse matrices and in various matrix-related operations.
 \n
 If unsure what flags to use to a function, #RSB_FLAG_NOFLAGS shall be a good default in most cases.
 */
/*!@{*/

/*! Default storage flags. */
#define RSB_FLAG_DEFAULT_STORAGE_FLAGS		 	(RSB_FLAG_WANT_BCSS_STORAGE|RSB_FLAG_WANT_COO_STORAGE)

/*! A flag combination specifying a pure COO matrix.  */
#define RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS		 	RSB_FLAG_WANT_COO_STORAGE 

/*! A flag combination specifying a pure CSR matrix.  */
#define RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS		 	RSB_FLAG_WANT_BCSS_STORAGE

/*! A flag combination specifying a pure RSB matrix.  */
#define RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS (RSB_FLAG_QUAD_PARTITIONING|RSB_FLAG_USE_HALFWORD_INDICES|RSB_FLAG_WANT_COO_STORAGE|RSB_FLAG_WANT_BCSS_STORAGE)

/*! A flag combination specifying a matrix in a default, supported format.  */
#define RSB_FLAG_DEFAULT_MATRIX_FLAGS			RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS

/*! The null (empty) flag. */
#define RSB_FLAG_NOFLAGS		 		0x000000

/*! The identical flag (used in cloning function #rsb_mtx_clone). */
#define RSB_FLAG_IDENTICAL_FLAGS RSB_FLAG_NOFLAGS

/*! If set, the input/output coordinate indices will be assumed to be 1 based. */
#define RSB_FLAG_FORTRAN_INDICES_INTERFACE		0x000001

/*! If set, the input/output coordinate indices will be assumed to be 0 based (default). */
#define RSB_FLAG_C_INDICES_INTERFACE		0x000000

/*! If set, the matrix will internally use a half word (16 bit) type for indices. */
#define RSB_FLAG_USE_HALFWORD_INDICES	0x000002

/*! Used to specify multi-vector (dense matrix) operations. */
#define RSB_FLAG_WANT_ROW_MAJOR_ORDER 			0x000000

/*! Used to specify multi-vector (dense matrix) operations. */
#define RSB_FLAG_WANT_COLUMN_MAJOR_ORDER 		0x4000000

/*! If set, the code will assume the input nonzeroes as sorted.	*/
#define RSB_FLAG_SORTED_INPUT				0x000004

/*! If set, the matrix is considered as triangular. \see #RSB_FLAG_LOWER,#RSB_FLAG_UPPER. */ 
#define RSB_FLAG_TRIANGULAR 				0x000008

/*! If set, the matrix will be stored in as lower (triangular or symmetric). \see #RSB_FLAG_TRIANGULAR,#RSB_FLAG_SYMMETRIC,#RSB_FLAG_UPPER. */
#define RSB_FLAG_LOWER		 			0x000010

/*! If set, the matrix will be stored in as upper (triangular or symmetric). \see #RSB_FLAG_LOWER*/
#define RSB_FLAG_UPPER		 			0x000020

/*! If set, the (whole super-)matrix will not store the diagonal, which will be assumed to be unitary. */
#define RSB_FLAG_UNIT_DIAG_IMPLICIT			0x000040

/* ghost flag ( moved in a non-public header, and reserved): 0x000080	*/
/* ghost flag ( moved in a non-public header, and reserved): 0x80000000	*/

/*! If set, the matrix will use COO storage, where necessary. */
#define RSB_FLAG_WANT_COO_STORAGE		0x000100

/*! Keep the last nonzero duplicate, at matrix assembly time. */ 
#define RSB_FLAG_DUPLICATES_KEEP_LAST				0x000000

/*! The default nonzeroes duplicates handling.  */ 
#define RSB_FLAG_DUPLICATES_DEFAULT_HANDLE			0x000000

/*! Compute and keep the sum of nonzero duplicates, at matrix assembly time.  */ 
#define RSB_FLAG_DUPLICATES_SUM				0x000200

/*! If set, explicit zeros will not be inserted	\warning: this flag is active by default	*/
#define RSB_FLAG_DISCARD_ZEROS				0x000400

/* ghost flag ( moved in a non-public header, and reserved): 0x000800 */
/* ghost flag ( moved in a non-public header, and reserved): 0x001000 */

/*! If set, matrix will be organized as a quad tree of submatrices. */
#define RSB_FLAG_QUAD_PARTITIONING 			0x002000

/*! If set, the block partitioning will be fixed (BCSS: BCSR or BCSC, but no VBR).	*/
#define RSB_FLAG_WANT_BCSS_STORAGE 			0x004000

/* ghost flag ( moved in a non-public header, and reserved): 0x008000 */
/* ghost flag ( moved in a non-public header, and reserved): 0x010000 */
/* ghost flag ( moved in a non-public header, and reserved): 0x020000 */

/*! If set, matrices will be fit in the three input coo arrays, after conversion. */ 
#define RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS		0x040000

/*! \internal \todo: should remove this. */ 
#define RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT	0x080000

/* ghost flag ( moved in a non-public header, and reserved): 0x100000*/
/* ghost flag (temporarily reserved): 0x200000*/

/*! If set, the input matrix will be treated as symmetric (stored as a lower triangular one by default). \see #RSB_FLAG_LOWER,#RSB_FLAG_LOWER. */
#define RSB_FLAG_SYMMETRIC 			0x400000

/*! If set, the input matrix will be treated as symmetric hermitian (stored as a lower triangular one). \see #RSB_FLAG_LOWER,#RSB_FLAG_LOWER. */
#define RSB_FLAG_HERMITIAN 			0x800000

/*! If set, recursion on small matrices will last at least the number of active threads. */
#define RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS	0x1000000

/* ghost flag ( moved in a non-public header, and reserved): 0x2000000	*/

/*! Combined flags for a lower hermitian matrix. */
#define RSB_FLAG_LOWER_HERMITIAN			(RSB_FLAG_HERMITIAN | RSB_FLAG_LOWER)

/*! Combined flags for an upper hermitian matrix. */
#define RSB_FLAG_UPPER_HERMITIAN			(RSB_FLAG_HERMITIAN | RSB_FLAG_UPPER)

/*! Combined flags for a lower triangular matrix. */
#define RSB_FLAG_LOWER_TRIANGULAR 			(RSB_FLAG_TRIANGULAR | RSB_FLAG_LOWER)

/*! Combined flags for an upper triangular matrix. */
#define RSB_FLAG_UPPER_TRIANGULAR 			(RSB_FLAG_TRIANGULAR | RSB_FLAG_UPPER)

/*! Combined flags for a symmetric, lower-stored matrix. */

#define RSB_FLAG_LOWER_SYMMETRIC 			(RSB_FLAG_SYMMETRIC | RSB_FLAG_LOWER)

/*! Combined flags for a diagonal matrix. */
#define RSB_FLAG_DIAGONAL 				(RSB_FLAG_UPPER_TRIANGULAR | RSB_FLAG_LOWER_TRIANGULAR)

/*! Combined flags for a symmetric, upper-stored matrix. */
#define RSB_FLAG_UPPER_SYMMETRIC 			(RSB_FLAG_SYMMETRIC | RSB_FLAG_UPPER)

/*! If set, the matrix will be subdivided at a finer grain on diagonal blocks. */
#define RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG 	0x8000000

/*! If set, the input COO arrays to the assembly functions will not be freed at matrix destruction time.
  \warning Please do NOT use this flag, for the default memory allocation handling is still not specified. Instead, use the in place allocation functions: #rsb_mtx_alloc_from_csr_inplace() and #rsb_mtx_alloc_from_coo_inplace().
 */
#define RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS 		0x40000000

/* Reserved, undocumented flags. Not for use. */
#define RSB_FLAG_USE_CSR_RESERVED	0x200000

/*! \internal Combined flags for half word CSR. */
#define RSB_FLAG_USE_HALFWORD_INDICES_CSR	(RSB_FLAG_USE_HALFWORD_INDICES|RSB_FLAG_USE_CSR_RESERVED)

/*! Combined flags for half word COO. */
#define RSB_FLAG_USE_HALFWORD_INDICES_COO	(RSB_FLAG_USE_HALFWORD_INDICES|RSB_FLAG_WANT_COO_STORAGE)

/*! A combination of flags which is forbidden (so don't use it). */
#define RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES	(RSB_FLAG_USE_HALFWORD_INDICES_COO|RSB_FLAG_USE_HALFWORD_INDICES_CSR)
/*!@}*/


/*! A macro for the error code value. */
#define RSB_ERR_CAST(E) (-(E))

/* Error handling functions. */
rsb_err_t rsb_strerror_r(rsb_err_t errval, rsb_char_t * buf, size_t buflen);
rsb_err_t rsb_perror(void *stream, rsb_err_t errval);

/*! No error occurred (success). The return value that means function operation success, in most cases.   */
#define RSB_ERR_NO_ERROR		RSB_ERR_CAST(0x000)

/*! An unspecified, generic error occurred. */
#define RSB_ERR_GENERIC_ERROR		RSB_ERR_CAST(0x001)

/*! The user requested an operation which is not supported (e.g.: was opted out at build time). */
#define RSB_ERR_UNSUPPORTED_OPERATION	RSB_ERR_CAST(0x002)

/*! The user requested to use a type which is not supported (e.g.: was opted out at build time). */
#define RSB_ERR_UNSUPPORTED_TYPE	RSB_ERR_CAST(0x004)

/*! The user requested to use a matrix storage format which is not supported (e.g.: was opted out at build time). */
#define RSB_ERR_UNSUPPORTED_FORMAT	RSB_ERR_CAST(0x008)

/*! An error occurred which is not apparently caused by a user's fault (internal error). */
#define RSB_ERR_INTERNAL_ERROR		RSB_ERR_CAST(0x010)

/*! The user supplied corrupt or inconsistent data as argument. */
#define RSB_ERR_BADARGS			RSB_ERR_CAST(0x020)

/*! There is not enough dynamical memory to perform the requested operation. */
#define RSB_ERR_ENOMEM			RSB_ERR_CAST(0x040)

/*! The requested operation was not implemented yet in this code revision (but probably will be, someday). */
#define RSB_ERR_UNIMPLEMENTED_YET	RSB_ERR_CAST(0x100)

/*! The requested operation could not be executed, or index overflow will happen. */
#define RSB_ERR_LIMITS			RSB_ERR_CAST(0x200)

/*! A Fortran-specific error occurred. */
#define RSB_ERR_FORTRAN_ERROR		RSB_ERR_GENERIC_ERROR

/*! The requested feature is not available because it was opted out or not configured at build time. */
#define RSB_ERR_UNSUPPORTED_FEATURE	RSB_ERR_CAST(0x400)

/*! A file containing user set configuration was not present. */
#define RSB_ERR_NO_USER_CONFIGURATION	RSB_ERR_CAST(0x800)

/*! User-supplied data (e.g.: from file) was corrupt. */
#define RSB_ERR_CORRUPT_INPUT_DATA	RSB_ERR_CAST(0x1000)

/*! Memory hierarchy info failed to be detected. You can bypass this by setting a meaningful \c RSB_USER_SET_MEM_HIERARCHY_INFO environment variable. */
#define RSB_ERR_FAILED_MEMHIER_DETECTION	RSB_ERR_CAST(0x2000)

/*! User gave flags for an in place assembly in a copy-based function. */
#define RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS	RSB_ERR_CAST(0x4000)

/*! User requested writing to a file stream, while this feature is configured out. */
#define RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT	RSB_ERR_CAST(0x8000)

/*! User gave some input with invalid numerical data. */
#define RSB_ERR_INVALID_NUMERICAL_DATA	RSB_ERR_CAST(0x10000)

/*! Probable memory leak (user did not deallocate librsb structures before calling rsb_lib_exit()). */
#define RSB_ERR_MEMORY_LEAK	RSB_ERR_CAST(0x20000)

/*! Element not found by rsb_mtx_get_vals() or rsb_mtx_set_vals(). */
#define RSB_ERR_ELEMENT_NOT_FOUND RSB_ERR_CAST(0x40000000)

/*! Collation of "unsupported" type errors. */
#define RSB_ERRS_UNSUPPORTED_FEATURES	(RSB_ERR_UNSUPPORTED_FEATURE|RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT)

/*! Program success error code (int). */
#ifdef EXIT_SUCCESS
#define RSB_PROGRAM_SUCCESS	(EXIT_SUCCESS)
#else
#define RSB_PROGRAM_SUCCESS		(0)
#endif

/*! Program error code (int). */
#ifdef EXIT_FAILURE
#define RSB_PROGRAM_ERROR		(EXIT_FAILURE)
#else
#define RSB_PROGRAM_ERROR		(-1)
#endif

/*! Program error code (int). */
#define RSB_ERR_TO_PROGRAM_ERROR(E)	((E)==(RSB_ERR_NO_ERROR)?RSB_PROGRAM_SUCCESS:RSB_PROGRAM_ERROR)


/*! \ingroup rsb_doc_misc rsb_doc_rsb
\brief library option values for \see_lib_init_funcs. */
enum rsb_opt_t
{
/*! #RSB_IO_WANT_VERBOSE_INIT prompts for a verbose initialization of the library: messages will be written
 * to the file descriptor (\c FILE*) pointed by the value pointer when calling \ref rsb_lib_init  or \ref rsb_lib_reinit.
 */
  RSB_IO_WANT_VERBOSE_INIT =0x000001	/* (FILE*) */
,
/*! #RSB_IO_WANT_VERBOSE_EXIT prompts for a verbose finalization of the library: messages will be written
 * to the file descriptor (\c FILE*) pointed by the value pointer when calling \ref rsb_lib_exit.
 */
  RSB_IO_WANT_VERBOSE_EXIT =0x000002	/* (FILE*) */
,
/*! Specifies the default output stream. Output (debug info) info will be written  
 * to the file descriptor (\c FILE*) pointed by the value pointer.
 */
  RSB_IO_WANT_OUTPUT_STREAM =0x000003	/* (FILE*) */
,
/*! Specifies the default sorting method. Specified as a pointed integer (#rsb_int_t) number, in {[0],1}. (internal)
 */
  RSB_IO_WANT_SORT_METHOD =0x000004	/* (rsb_int_t) */
,
/*! Specifies the default cache blocking method. Specified as a pointed integer (#rsb_int_t) number, in {-1,[0],1}. (internal)
 */
  RSB_IO_WANT_CACHE_BLOCKING_METHOD =0x000005	/* (rsb_int_t) */
,
/*! Specifies a multiplier for finer (if >1.0) or coarser (if <1.0) subdivisions. Specified as a pointed (#rsb_real_t) number, in {..,[1.0],..}. (internal)
 */
  RSB_IO_WANT_SUBDIVISION_MULTIPLIER =0x000006	/* (rsb_real_t) */
,
/*! Prompts for a verbose error reporting: messages will be written
 * to the file descriptor (\c FILE*) pointed by the value pointer. Only meaningful if an
 * interface error verbosity greater than 0 was set at configure time. 
 */
  RSB_IO_WANT_VERBOSE_ERRORS =0x000007	/* (FILE*) */
,
/*! Prompts for bounded box computation, for a smoother submatrices locking; pointed #rsb_int_t in {0,[1]}. (internal).
 */
  RSB_IO_WANT_BOUNDED_BOX_COMPUTATION =0x000008	/* (rsb_int_t) */
,
/*! Specifies the number of desired executing threads; pointed #rsb_int_t in {[0],1,..}.
 */
  RSB_IO_WANT_EXECUTING_THREADS =0x000009	/* (rsb_int_t) */
,
/*! Specifies the level of interface verbosity; if setting, pointed #rsb_int_t values should be in {[0],1,..}. Support may be enabled or disabled at build time via the \c --enable-internals-error-verbosity configure option. If disabled, only getting is supported and yields -1, but setting is not supported and the #RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT error will be returned.
 */
  RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE =0x000010	/* (rsb_int_t) */
,
/*! Specifies a custom memory hierarchy info string; pointed \c const #rsb_char_t*; (may point to a NULL string pointer).
 */
  RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING =0x000011	/* (const rsb_char_t*) */
,
/*! Used for getting whether the library has been initialized (#RSB_BOOL_TRUE) or not (#RSB_BOOL_FALSE) ; pointed \c const #rsb_bool_t*; (this is NOT for general users).
 */
  RSB_IO_WANT_IS_INITIALIZED_MARKER =0x000012	/* (const rsb_bool_t*) */
,
/*! Used for getting the count of memory allocations performed by librsb employing librsb's memory allocation wrapper (if disabled, will return zero); pointed \c const \c size_t*; (this is for debugging purposes).
 */
  RSB_IO_WANT_MEM_ALLOC_CNT =0x000013	/* (const size_t*) */
,
/*! Used for getting the total amount of memory allocated by librsb employing librsb's memory allocation wrapper (if disabled, will return zero); pointed \c const \c size_t*; (this is for debugging purposes).
 */
  RSB_IO_WANT_MEM_ALLOC_TOT =0x000014	/* (const size_t*) */
,
/*! Specifies whether the default multi-vector ops shall act at a leaf level (default value of 0 is yes). Specified as a pointed integer (#rsb_int_t) number, in {-1,[0]}. (internal)
 */
  RSB_IO_WANT_LEAF_LEVEL_MULTIVEC =0x000015	/* (rsb_int_t) */
,
/*! Specifies an upper limit to the count of allocated memory areas (default value of 0 means no limit). Specified as a pointed \c size_t. \rsb_configure_memwrap
 */
  RSB_IO_WANT_MAX_MEMORY_ALLOCATIONS =0x000016	/* (size_t) */
,
/*! Specifies an upper limit to the amount of allocated memory (default value of 0 means no limit). Specified as a pointed \c size_t. \rsb_configure_memwrap
 */
  RSB_IO_WANT_MAX_MEMORY_ALLOCATED =0x000017	/* (size_t) */
,
/*! Represents time spent in librsb. Specified as a pointed #rsb_time_t. Only works if statistics collection (\c --enable-librsb-stats) was specified at configure time.
 */
  RSB_IO_WANT_LIBRSB_ETIME =0x000018	/* (rsb_time_t) */
,
/*! Auto tuning verbosity level for rsb_tune_spmm/rsb_tune_spsm. If 0, no verbosity; if 1, verbose; if 2, verbose with trace files being dumped.
 */
  RSB_IO_WANT_VERBOSE_TUNING =0x000019	/* (rsb_int_t) */
};

/*! A handy macro for invoking #rsb_lib_reinit() with a single get/set specifier.
 * An appropriate I/O flag is supplied as first parameter; a valid pointer (according to the flag) should be passed as second parameter; either #RSB_IO_SPECIFIER_SET or #RSB_IO_SPECIFIER_GET is passed as third parameter; a #rsb_err_t variable as fourth one, in order to detect any error.
 * \deprecated	This macro has been deprecated and will be removed in a future version: use #rsb_lib_set_opt or #rsb_lib_get_opt instead.
 * */
#define RSB_REINIT_SINGLE_VALUE(IOF,IOP,IOS,ERRVAL) { enum rsb_opt_t keys[]={IOF}; void*values[]={(IOP)}; struct rsb_initopts io; io.action=(IOS); io.keys=keys; io.values=values; io.n_pairs=1; ERRVAL=rsb_lib_reinit(&io); }

/*! Like #RSB_REINIT_SINGLE_VALUE, but considering \c IOP \c const.
 * \deprecated	This macro has been deprecated and will be removed in a future version: use #rsb_lib_set_opt instead.
 * */
#define RSB_REINIT_SINGLE_VALUE_C_IOP(IOF,IOP,IOS,ERRVAL) { enum rsb_opt_t keys[]={IOF}; const void*values[]={(IOP)}; struct rsb_initopts io; io.action=(IOS); io.keys=keys; (io.values)=(void**)values; io.n_pairs=1; ERRVAL=rsb_lib_reinit(&io); }

/*! A handy macro for invoking #RSB_REINIT_SINGLE_VALUE with a single set specifier.
 * An appropriate I/O flag is supplied as first parameter; a valid pointer (according to the flag) should be passed as second parameter; a #rsb_err_t variable as third one, in order to detect any error.
 * \deprecated	This macro has been deprecated and will be removed in a future version: use #rsb_lib_set_opt instead.
 * */
#define RSB_REINIT_SINGLE_VALUE_SET(IOF,IOP,ERRVAL) RSB_REINIT_SINGLE_VALUE(IOF,IOP,RSB_IO_SPECIFIER_SET,ERRVAL)

/*! A handy macro for invoking #RSB_REINIT_SINGLE_VALUE with a single get specifier.
 * An appropriate I/O flag is supplied as first parameter; a valid pointer (according to the flag) should be passed as second parameter; a #rsb_err_t variable as third one, in order to detect any error.
 * \deprecated	This macro has been deprecated and will be removed in a future version: use #rsb_lib_get_opt instead.
 * */
#define RSB_REINIT_SINGLE_VALUE_GET(IOF,IOP,ERRVAL) RSB_REINIT_SINGLE_VALUE(IOF,IOP,RSB_IO_SPECIFIER_GET,ERRVAL)


/*!
 * @brief A structure specifying library (initialization) options, to be used with the \ref rsb_lib_reinit() function.
 * \n
 *
 * The structure specifies, for \c i=0,..,n_pairs-1 , a list of (key,value)
 * pairs, stored respectively as (\c keys[i],values[i]).
 * \n
 * Each flag specifies the type and possible range of values it accepts. 
 * \n
 * The structure may he used to set or query various library parameters. 
 *
 * Example:
 * \code
 	const int max_io=10; // the number of different options we want to set
	struct rsb_initopts io={NULL,NULL,0,RSB_IO_SPECIFIER_SET},
 	*iop=&io; // pointer to the options structure
	void * io_values[max_io]; // an array of pointers to max_io different option values (we shall set)
	enum rsb_opt_t io_keys[max_io]; // an array of max_io flag values specifying the type of values we are handing over to the library
	io.keys=io_keys; // io.keys will now point to io_keys as its keys array
	io.values=io_values; // io.values will now point to io_keys as its values array
	io.n_pairs=0; // we have 0 pairs specified so far
	io.keys[io.n_pairs]=RSB_IO_WANT_BOUNDED_BOX_COMPUTATION; // the first (at index 0) option we want to specify is RSB_IO_WANT_BOUNDED_BOX_COMPUTATION
	io.values[io.n_pairs]=1; // the value we want to set the RSB_IO_WANT_BOUNDED_BOX_COMPUTATION option to
	io.n_pairs++; // io.n_pairs is set to 1: we have one option set, so even if we have (max_io-io.n_pairs) left, only the first will be read
	... // we are free to specify other option (type, value) pairs
 * \endcode
 * */
struct rsb_initopts
{
	/*! An array of value types key flags. */
	enum rsb_opt_t * keys;
	/*! An array of value pointers, as specified by each flag value. */
	void ** values;
	/*! The length of the \c keys and \c values arrays. */
	rsb_int_t n_pairs;
	/*! The action we are requesting (either one of #RSB_IO_SPECIFIER_GET or #RSB_IO_SPECIFIER_SET)*/
	rsb_int_t action;
};

#define RSB_IO_SPECIFIER_GET	1 /*!< Specifies to #RSB_REINIT_SINGLE_VALUE that a given #rsb_initopts is going to be get by the user. */
#define RSB_IO_SPECIFIER_SET	0 /*!< Specifies to #RSB_REINIT_SINGLE_VALUE that a given #rsb_initopts is going to be set by the user. */
#define RSB_NULL_INIT_OPTIONS NULL /*!<  A valid value for specifying default (null) options to #rsb_lib_init().  */
#define RSB_NULL_EXIT_OPTIONS NULL /*!<  A valid value for specifying default (null) options to #rsb_lib_exit().  */

rsb_err_t rsb_lib_init(struct rsb_initopts * iop);
rsb_err_t rsb_lib_reinit(struct rsb_initopts * iop);
rsb_err_t rsb_lib_set_opt_str(const rsb_char_t* opnp, const rsb_char_t* opvp);
rsb_err_t rsb_lib_set_opt(enum rsb_opt_t iof, const void*iop);
rsb_err_t rsb_lib_get_opt(enum rsb_opt_t iof, void*iop);
rsb_err_t rsb_lib_exit(struct rsb_initopts * iop);

struct rsb_mtx_t * rsb_mtx_alloc_from_coo_begin(rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_flags_t flagsA, rsb_err_t * errvalp);
rsb_err_t rsb_mtx_alloc_from_coo_end(struct rsb_mtx_t ** mtxApp);
struct rsb_mtx_t * rsb_mtx_alloc_from_csr_const(const void *VA, const rsb_coo_idx_t * RP, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb_mtx_alloc_from_csc_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * CP, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb_mtx_alloc_from_csr_inplace(void * VA, rsb_nnz_idx_t * RP, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb_mtx_alloc_from_coo_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb_mtx_alloc_from_coo_inplace(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t * errvalp );
rsb_err_t rsb_mtx_clone(struct rsb_mtx_t ** mtxBpp, rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_flags_t flags);
struct rsb_mtx_t * rsb_mtx_free(struct rsb_mtx_t * mtxAp);

/*! \ingroup rsb_doc_misc rsb_doc_rsb
 \brief Extraction filter flags, to be used with #rsb_mtx_get_nrm()/#rsb_mtx_get_vec(). */
enum rsb_extff_t
{
  RSB_EXTF_NORM_ONE	=0x00001001			/*!< #rsb_mtx_get_nrm() flag value for computing the one-norm. */
, RSB_EXTF_NORM_TWO	=0x00001002			/*!< #rsb_mtx_get_nrm() flag value for computing the two-norm (Frobenius norm). */
, RSB_EXTF_NORM_INF	=0x00001003			/*!< #rsb_mtx_get_nrm() flag value for computing the infinity-norm. */
, RSB_EXTF_SUMS_ROW	=0x00001004			/*!< #rsb_mtx_get_vec() flag value for computing the sum along each row. */
, RSB_EXTF_SUMS_COL	=0x00001005			/*!< #rsb_mtx_get_vec() flag value for computing the sum along each column. */
, RSB_EXTF_ASUMS_ROW	=0x00001006			/*!< #rsb_mtx_get_vec() flag value for computing the absolute values sum, along each row. */
, RSB_EXTF_ASUMS_COL	=0x00001007			/*!< #rsb_mtx_get_vec() flag value for computing the absolute values sum, along each column. */
, RSB_EXTF_DIAG		=0x00000004			/*!< #rsb_mtx_get_vec() flag value for extracting the diagonal submatrix.*/
};

typedef rsb_flags_t rsb_marf_t;					/*!< Matrix rendering flags (see \ref marf_section for possible values). */
/*!@{*/
/*!
 \ingroup rsb_doc_rsb
 \name Matrix rendering flags
 \anchor marf_section

 These are flags which could be combined to specify rendering options to #rsb_mtx_rndr and #rsb_file_mtx_rndr.
 */
#define RSB_MARF_RGB	0x00000001			/*!< #rsb_marf_t Flag value for requesting an RGB rendering of a matrix. */
#define RSB_MARF_EPS_S	0x00000010			/*!< #rsb_marf_t Flag value for requesting an Encapsulated Postscript rendering of a matrix (spy plot). */
#define RSB_MARF_EPS_B	0x00000020			/*!< #rsb_marf_t Flag value for requesting an Encapsulated Postscript rendering of a matrix (blocks plot). */
#define RSB_MARF_EPS	0x00000030			/*!< #rsb_marf_t Flag value for requesting an Encapsulated Postscript rendering of a matrix (spy plot + blocks). */
#define RSB_MARF_EPS_L	0x00000070			/*!< #rsb_marf_t Flag value for requesting an Encapsulated Postscript rendering of a matrix (spy plot + blocks + labels). */
/*!@}*/

rsb_err_t rsb_mtx_get_nrm(const struct rsb_mtx_t * mtxAp , void * Np, enum rsb_extff_t flags);
#define rsb_mtx_get_norm rsb_mtx_get_nrm /*!< \deprecated #rsb_mtx_get_norm has been deprecated: use #rsb_mtx_get_nrm . */
rsb_err_t rsb_mtx_get_vec(const struct rsb_mtx_t * mtxAp , void * Dp, enum rsb_extff_t flags);
rsb_err_t rsb_mtx_rndr(const rsb_char_t * filename, const struct rsb_mtx_t*mtxAp, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags);
rsb_err_t rsb_file_mtx_rndr(void * pmp, const rsb_char_t * filename, rsb_coo_idx_t pmlWidth, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags);
#define rsb_file_mtx_render rsb_file_mtx_rndr /*!< \deprecated #rsb_file_mtx_render has been deprecated: use #rsb_file_mtx_rndr. */

rsb_err_t rsb_spmv(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY);
rsb_err_t rsb_spmm(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC);
rsb_err_t rsb_spsv(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, const void * Xp, rsb_coo_idx_t incX, void * Yp, rsb_coo_idx_t incY);
rsb_err_t rsb_spsm(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * betap, const void * Bp, rsb_nnz_idx_t ldB, void * Cp, rsb_nnz_idx_t ldC);

rsb_err_t rsb_mtx_add_to_dense(const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_nnz_idx_t ldB, rsb_nnz_idx_t nrB, rsb_nnz_idx_t ncB, rsb_bool_t rowmajorB, void * Bp);

struct rsb_mtx_t * rsb_sppsp(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_err_t * errvalp);
struct rsb_mtx_t * rsb_spmsp(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_err_t * errvalp);
rsb_err_t rsb_spmsp_to_dense(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp , rsb_nnz_idx_t ldC, rsb_nnz_idx_t nrC, rsb_nnz_idx_t ncC, rsb_bool_t rowmajorC, void *Cp);

rsb_err_t rsb_mtx_switch_to_coo(struct rsb_mtx_t * mtxAp, void ** VAp, rsb_coo_idx_t ** IAp, rsb_coo_idx_t ** JAp, rsb_flags_t flags);
rsb_err_t rsb_mtx_switch_to_csr(struct rsb_mtx_t * mtxAp, void ** VAp, rsb_coo_idx_t ** IAp, rsb_coo_idx_t ** JAp, rsb_flags_t flags);
rsb_err_t rsb_mtx_get_coo(const struct rsb_mtx_t * mtxAp, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_flags_t flags );
rsb_err_t rsb_mtx_get_csr(rsb_type_t typecode, const struct rsb_mtx_t * mtxAp, void * VA, rsb_nnz_idx_t * RP, rsb_coo_idx_t * JA, rsb_flags_t flags );
rsb_err_t rsb_mtx_get_rows_sparse(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags );
rsb_err_t rsb_mtx_get_coo_block(const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_coo_idx_t fcA, rsb_coo_idx_t lcA, const rsb_coo_idx_t * IREN, const rsb_coo_idx_t * JREN, rsb_nnz_idx_t *rnzp, rsb_flags_t flags );

/*! \ingroup rsb_doc_misc rsb_doc_rsb
\brief Flags for getting matrix information via #rsb_mtx_get_info()/#rsb_mtx_get_info_str().
*/
enum rsb_mif_t
{
  RSB_MIF_INDEX_STORAGE_IN_BYTES__TO__SIZE_T		=0x00000001	/*!< Index storage occupation, in bytes. (size_t) */
, RSB_MIF_INDEX_STORAGE_IN_BYTES_PER_NNZ__TO__RSB_REAL_T	=0x00000002	/*!< Index storage occupation per nnz, in bytes. (#rsb_real_t) */
, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T		=0x00000004	/*!< Rows count(#rsb_coo_idx_t) */
, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T		=0x00000008	/*!< Columns count (#rsb_coo_idx_t) */
, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T		=0x00000010	/*!< Nonzeroes count (#rsb_nnz_idx_t) */
, RSB_MIF_TOTAL_SIZE__TO__SIZE_T			=0x00000020	/*!< Total size, in bytes (size_t) */
, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T			=0x00000040	/*!< Matrix flags (#rsb_flags_t) */
, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T		=0x00000080	/*!< Matrix type code (#rsb_type_t) */
, RSB_MIF_MATRIX_INFO__TO__CHAR_P			=0x00000100	/*!< Matrix info string, only for #rsb_mtx_get_info_str() (#rsb_char_t*) */
, RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T		=0x00000200	/*!< Leaf submatrices count (#rsb_blk_idx_t) */
};
rsb_err_t rsb_mtx_get_info(const struct rsb_mtx_t *mtxAp, enum rsb_mif_t miflags, void* minfop);
rsb_err_t rsb_mtx_get_info_str(const struct rsb_mtx_t *mtxAp, const rsb_char_t *mis, void* minfop, size_t buflen);

/*! \ingroup rsb_doc_misc rsb_doc_rsb
\brief Flags for specifying a particular elemental/row-wise operation with #rsb_mtx_upd_vals(). */
enum rsb_elopf_t
{
  RSB_ELOPF_MUL		=0x00000001		/*!< Elemental multiplication of the matrix by a specified scalar (usable with #rsb_mtx_upd_vals(), binary operation). */
, RSB_ELOPF_DIV		=0x00000002		/*!< Elemental division by a specified scalar (usable with #rsb_mtx_upd_vals(), binary operation). */
, RSB_ELOPF_POW		=0x00000004		/*!< Elemental power to a specified scalar (usable with #rsb_mtx_upd_vals(), binary operation). */
, RSB_ELOPF_NEG		=0x00000008		/*!< Elemental negation (usable with #rsb_mtx_upd_vals(), unary operation). */
, RSB_ELOPF_SCALE_ROWS	=0x00000010		/*!< Row    wise scaling by a specified scaling vector (usable with #rsb_mtx_upd_vals(), binary operation). */
, RSB_ELOPF_SCALE_COLS	=0x00000020		/*!< Column wise scaling by a specified scaling vector (usable with #rsb_mtx_upd_vals(), binary operation). */
, RSB_ELOPF_SCALE_ROWS_REAL	=0x00000040	/*!< Row    wise scaling by a specified scaling vector. If matrix is of a complex type, the argument is expected to be of the corresponding real type (assumed that that type has been enabled). (usable with #rsb_mtx_upd_vals(), binary operation). */
, RSB_ELOPF_SCALE_COLS_REAL	=0x00000080	/*!< Column wise scaling by a specified scaling vector. If matrix is of a complex type, the argument is expected to be of the corresponding real type (assumed that that type has been enabled). (usable with #rsb_mtx_upd_vals(), binary operation). */
};

rsb_err_t rsb_mtx_upd_vals(struct rsb_mtx_t * mtxAp, enum rsb_elopf_t elop_flags, const void * omegap);
#define rsb_mtx_upd_values rsb_mtx_upd_vals /*!< \deprecated #rsb_mtx_upd_values has been deprecated: use #rsb_mtx_upd_vals. */
typedef rsb_flags_t rsb_precf_t;	/*!< Basic preconditioner flags to be used with #rsb_mtx_get_prec(). */
#define RSB_PRECF_ILU0		0x00000001		/*!< ILU-0 preconditioner request to #rsb_mtx_get_prec(). */
rsb_err_t rsb_mtx_get_prec(void *opdp, const struct rsb_mtx_t * mtxAp, rsb_precf_t prec_flags, const void *ipdp);
#define rsb_mtx_get_preconditioner rsb_mtx_get_prec 	/*!< \deprecated #rsb_mtx_get_preconditioner has been deprecated: use #rsb_mtx_get_prec. */
rsb_err_t rsb_mtx_set_vals(struct rsb_mtx_t * mtxAp, const void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags);
#define rsb_mtx_set_values rsb_mtx_set_vals /*!< \deprecated #rsb_mtx_set_values has been deprecated: use #rsb_mtx_set_vals. */
rsb_err_t rsb_mtx_get_vals(const struct rsb_mtx_t * mtxAp, void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags);
#define rsb_mtx_get_values rsb_mtx_get_vals /*!< \deprecated #rsb_mtx_get_values has been deprecated: use #rsb_mtx_get_vals. */
rsb_err_t rsb_tune_spmm(struct rsb_mtx_t ** mtxOpp, rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC);
rsb_err_t rsb_tune_spsm(struct rsb_mtx_t ** mtxOpp, rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC);
rsb_trans_t rsb_psblas_trans_to_rsb_trans(const char psbtrans);
rsb_err_t rsb_file_mtx_save(const struct rsb_mtx_t * mtxAp, const rsb_char_t * filename);
struct rsb_mtx_t * rsb_file_mtx_load(const rsb_char_t * filename, rsb_flags_t flagsA, rsb_type_t typecode, rsb_err_t *errvalp);
rsb_err_t rsb_file_vec_load(const rsb_char_t * filename, rsb_type_t typecode, void * Yp, rsb_coo_idx_t *yvlp);
rsb_err_t rsb_file_vec_save(const rsb_char_t * filename, rsb_type_t typecode, const void * Yp, rsb_coo_idx_t yvl);
rsb_err_t rsb_file_mtx_get_dims(const rsb_char_t * filename, rsb_coo_idx_t* nrp, rsb_coo_idx_t *ncp, rsb_coo_idx_t *nzp, rsb_flags_t*flagsp);
#define rsb_file_mtx_get_dimensions rsb_file_mtx_get_dims /*!< \deprecated #rsb_file_mtx_get_dimensions has been deprecated: use #rsb_file_mtx_get_dims. */
rsb_err_t rsb_coo_sort(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA,  rsb_type_t typecode, rsb_flags_t flagsA );
rsb_err_t rsb_coo_cleanup(rsb_coo_idx_t* nnzp, void* VA, rsb_coo_idx_t* IA, rsb_coo_idx_t* JA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_type_t typecode, rsb_flags_t flagsA );
rsb_time_t rsb_time(void);

/*! \ingroup rsb_doc_misc rsb_doc_rsb
 Use #RSB_SIZEOF macro to get the size (in bytes) of a type supported by the library (e.g.: when allocating numerical vectors).
 */
#define RSB_SIZEOF(TYPE) RSB_NUMERICAL_TYPE_SIZE(TYPE)

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif	/* RSB_RSB_H_INCLUDED */
