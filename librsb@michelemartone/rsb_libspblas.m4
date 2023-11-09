dnl
dnl
include(`libspblas_macros.m4')dnl
include(`rsb_fortran_macros.m4')dnl
include(`rsb_misc.m4')dnl
dnl
/*!
        @file
        @author Michele Martone

ifdef(`ONLY_WANT_HEADERS',`dnl
	@brief  This file specifies the Sparse BLAS interface to librsb.
',`dnl
dnl /* @cond INNERDOC  */
	@brief  This file implements Sparse BLAS for librsb.
')dnl
dnl
dnl	all types        :RSB_M4_SPBLAS_MATRIX_ALL_TYPES
	Supported types  :RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES .
	Unsupported types:RSB_M4_SPBLAS_MATRIX_UNSUPPORTED_TYPES .
	Level 1 ops      :RSB_M4_SPBLAS_MATRIX_ALL_L1_MOPS .
	Level 2 ops      :RSB_M4_SPBLAS_MATRIX_ALL_L2_MOPS .
	Level 3 ops      :RSB_M4_SPBLAS_MATRIX_ALL_L3_MOPS .
*/

ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_LIBSPBLAS_H_INCLUDED
#define RSB_LIBSPBLAS_H_INCLUDED
dnl typedef int rsb_blas_int_t;
dnl
#ifndef RSB_RSB_H_INCLUDED
#error "You are using Sparse BLAS headers from librsb -- You should include <rsb.h> first!"
#endif /* RSB_RSB_H_INCLUDED */
dnl
dnl
ifelse(RSB_M4_LONG_IDX,`0',`',`dnl
/* The user wants a long type for indices. */
#if ( defined(__cplusplus) && (__cplusplus>=201103L) )
#include <cstdint>
#else  /* __cplusplus */
#include <stdint.h>
#endif /* __cplusplus */
#define rsb_blas_int_t int64_t
')dnl
dnl alternative A:
define(`rsb_blas_pname_t',`int')dnl here and in internal headers.
dnl alternative B:
dnl typedef int rsb_blas_pname_t;
dnl
dnl
',`dnl
')dnl
dnl
dnl #include "blas_sparse/blas_enum.h"
include(`blas_sparse/blas_enum.h')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl

/** the sparse matrix descriptor type */
typedef int blas_sparse_matrix;
',`dnl
#include "rsb_libspblas_handle.h"
')dnl

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

ifdef(`ONLY_WANT_HEADERS',`',`dnl
/* @cond INNERDOC  */
RSB_INTERNALS_COMMON_HEAD_DECLS
/* @endcond */
')dnl
dnl #include "blas_sparse/blas_sparse.h"
dnl #include "blas_sparse/blas_sparse_proto.h"
dnl
ifelse(RSB_M4_LONG_IDX,`0',`dnl
define(`rsb_blas_int_t',`int')dnl generate an interface using `int'; if commented, will use the `rsb_blas_int_t' typedef
',`dnl
dnl
')dnl
define(`RSB_WRAPPER_PREFIX',`rsb__wp_')`'dnl
dnl
define(`RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L1',`dnl
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
`#define' RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`ID',`0',lang) RSB_WRAPPER_PREFIX`'RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`ID',`0',lang)
`#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */'
')dnl
dnl
define(`RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L2',`dnl
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
`#define' RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`ID',`0',lang) RSB_WRAPPER_PREFIX`'RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`ID',`0',lang)
`#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */'
')dnl
dnl
define(`RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_CF',`dnl
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
`#define' RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(type,mop,tri,`ID',`0',lang) RSB_WRAPPER_PREFIX`'RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(type,mop,tri,`ID',`0',lang)
`#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */'
')dnl
dnl
define(`RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_CFNT',`dnl
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
`#define' RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(`',mop,tri,`ID',`0',lang) RSB_WRAPPER_PREFIX`'RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(`',mop,tri,`ID',`0',lang)
`#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */'
')dnl
dnl
define(`RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_EF',`dnl
#if !RSB_WITH_SPARSE_BLAS_INTERFACE
`#define' RSB_M4_SPBLAS_EXTRA_FUNCTION(type,mop,tri,`ID',`0',lang) RSB_WRAPPER_PREFIX`'RSB_M4_SPBLAS_EXTRA_FUNCTION(type,mop,tri,`ID',`0',lang)
`#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */'
')dnl
dnl
               /* Level 1 Computational Routines */
foreach(`mop',RSB_M4_SPBLAS_MATRIX_ALL_L1_MOPS,`dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_ALL_TYPES,`dnl
foreach(`tri',RSB_M4_SPBLAS_SYMMETRY_UL_CHARCODE,`dnl
foreach(`lang',RSB_M4_SPBLAS_MATRIX_ALL_LANGUAGES,`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L1`'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`function_declaration',`0',lang)dnl
',`dnl
RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L1`'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`function_definition',`0',lang)dnl
')dnl
')dnl
')dnl
')dnl

')dnl
               /* Level 2 Computational Routines */

foreach(`mop',RSB_M4_SPBLAS_MATRIX_ALL_L2_MOPS,`dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_ALL_TYPES,`dnl
foreach(`tri',RSB_M4_SPBLAS_SYMMETRY_UL_CHARCODE,`dnl
foreach(`lang',RSB_M4_SPBLAS_MATRIX_ALL_LANGUAGES,`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L2`'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`function_declaration',`0',lang)
',`dnl
RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L2`'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`function_definition',`0',lang)
')dnl
')dnl
')dnl
')dnl

')dnl
               /* Level 3 Computational Routines */

foreach(`mop',RSB_M4_SPBLAS_MATRIX_ALL_L3_MOPS,`dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_ALL_TYPES,`dnl
foreach(`tri',RSB_M4_SPBLAS_SYMMETRY_UL_CHARCODE,`dnl
foreach(`lang',RSB_M4_SPBLAS_MATRIX_ALL_LANGUAGES,`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L2`'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`function_declaration',`0',lang)
',`dnl
RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_L2`'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`function_definition',`0',lang)
')dnl
')dnl
')dnl
')dnl

')dnl
               /* Handle Management Routines */
               /*             +              */
               /* Creation Routines */
               /*             +              */
               /* Insertion Routines */
               /*             +              */
               /* Completion of Construction Routines */
               /*             +              */
               /* Matrix Property Routines */
               /*             +              */
               /* Destruction Routine */

foreach(`mop',RSB_M4_SPBLAS_MATRIX_CREATION_MOPS,`dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_ALL_TYPES,`dnl
foreach(`lang',RSB_M4_SPBLAS_MATRIX_ALL_LANGUAGES,`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_CF`'dnl
RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(type,mop,`',`function_declaration',`0',lang)`'dnl
',`dnl
RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_CF`'dnl
RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(type,mop,`',`function_definition',`0',lang)`'dnl
')dnl
')dnl
')dnl

')dnl


foreach(`mop',(`cr_end',`ds'),`dnl
foreach(`lang',RSB_M4_SPBLAS_MATRIX_ALL_LANGUAGES,`dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_CFNT`'dnl
RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(`',mop,`',`function_declaration',`0',lang)`'dnl
',`dnl
RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_CFNT`'dnl
RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(`',mop,`',`function_definition',`0',lang)`'dnl
')dnl
')dnl
')dnl

foreach(`mop',RSB_M4_SBLAS_EXTRA_INTERFACE_OPS,`dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_ALL_TYPES,`dnl
foreach(`lang',RSB_M4_SPBLAS_MATRIX_ALL_LANGUAGES,`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_EF`'dnl
RSB_M4_SPBLAS_EXTRA_FUNCTION(type,mop,`',`function_declaration',`0',lang)`'dnl
',`dnl
RSB_SPARSE_BLAS_INTERFACE_REWRAPPER_EF`'dnl
RSB_M4_SPBLAS_EXTRA_FUNCTION(type,mop,`',`function_definition',`0',lang)`'dnl
')dnl
')dnl
')dnl

')dnl



dnl
dnl
dnl

ifdef(`ONLY_WANT_HEADERS',`
`#define' BLAS_ussp RSB_WRAPPER_PREFIX`_'BLAS_ussp
`#define' BLAS_usgp RSB_WRAPPER_PREFIX`_'BLAS_usgp
rsb_blas_int_t BLAS_ussp( blas_sparse_matrix A, rsb_blas_pname_t pname );
rsb_blas_int_t BLAS_usgp( blas_sparse_matrix A, rsb_blas_pname_t pname );
',`dnl
rsb_blas_int_t BLAS_usgp( blas_sparse_matrix A, rsb_blas_pname_t pname ) /*  FIXME: temporarily here */
{
	/**
	 \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_gp_msg
	 \rsb_spblas_return_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	RSB_SPB_INTERFACE_RETURN_EXP(rsb__BLAS_usgp(A,pname))
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
`#define'  `'RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C`'blas_usgp`'RSB_M4_FORTRAN_SYMBOL_ADD_TO_C`'	RSB_WRAPPER_PREFIX`'`'RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C`'blas_usgp`'RSB_M4_FORTRAN_SYMBOL_ADD_TO_C`'
`#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */'
void `'RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C`'blas_usgp`'RSB_M4_FORTRAN_SYMBOL_ADD_TO_C`'( blas_sparse_matrix*A, rsb_blas_pname_t *pname, rsb_blas_int_t * istat ) /*  FIXME: temporarily here */
{
	/** \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_gp_msg
	 \rsb_spblas_istat_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	*istat=BLAS_usgp(*A,*pname); /*  FIXME: temporarily here */
	RSB_SPB_INTERFACE_RETURN_VOID()
}

#if !RSB_WITH_SPARSE_BLAS_INTERFACE
`#define'  `'RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C`'blas_ussp`'RSB_M4_FORTRAN_SYMBOL_ADD_TO_C`'	`'RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C`'RSB_WRAPPER_PREFIX`'blas_ussp`'RSB_M4_FORTRAN_SYMBOL_ADD_TO_C`'
`#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */'
void `'RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C`'blas_ussp`'RSB_M4_FORTRAN_SYMBOL_ADD_TO_C`'( blas_sparse_matrix*A, rsb_blas_pname_t *pname, rsb_blas_int_t * istat ) /*  FIXME: temporarily here */
{
	/**
	 \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_sp_msg
	 \rsb_spblas_istat_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	*istat=BLAS_ussp(*A,*pname); /*  FIXME: temporarily here */
	RSB_SPB_INTERFACE_RETURN_VOID()
}

rsb_blas_int_t BLAS_ussp( blas_sparse_matrix A, rsb_blas_pname_t pname ) /*  FIXME: temporarily here */
{
	/**
	 \ingroup rsb_doc_sparse_blas
	 \rsb_spblasl2_sp_msg
	 \rsb_spblas_return_msg
	 */
	RSB_SPB_INTERFACE_PREAMBLE
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_ussp(A,pname))
}

')

dnl
ifdef(`ONLY_WANT_HEADERS',`
blas_sparse_matrix rsb_blas_file_mtx_load(const rsb_char_t * filename, rsb_type_t typecode ); /* This is a librsb extension. */
',`dnl
blas_sparse_matrix rsb_blas_file_mtx_load(const rsb_char_t * filename, rsb_type_t typecode )
{
	/*!
 	\ingroup rsb_doc_sparse_blas
	Load Matrix Market matrix file of specified type to a matrix, and return Sparse BLAS handler.

	\param \rsb_filename_inp_param_msg
	\param \rsb_type_o_param_msg
	\return \rsb_spblasl2_A_msg_ftn_rt

	\n
	 */
	RSB_SPB_INTERFACE_PREAMBLE
{
	RSB_SPB_INTERFACE_RETURN_EXP( rsb__load_spblas_matrix_file_as_matrix_market(filename,typecode) )
}
}
')dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`
struct rsb_mtx_t * rsb_blas_get_mtx(blas_sparse_matrix A);
',`dnl
struct rsb_mtx_t * rsb_blas_get_mtx(blas_sparse_matrix A)
{
	/*!
 	\ingroup rsb_doc_sparse_blas
	\rsb_BLAS_get_mtx_msg

	\rsb_spblasl2_A_msg
	\return \rsbmtxpmessage_bg

	\n
	
	\rsb_BLAS_get_mtx_example
	\see_rsb_BLAS_get_mtx_msg
	\rsb_BLAS_get_mtx_msg_todo
	\rsb_BLAS_get_mtx_msg_note
	\rsb_BLAS_get_mtx_msg_warn
	 */
	RSB_SPB_INTERFACE_PREAMBLE
{
	struct rsb_mtx_t * mtxAp = NULL;
	mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	RSB_SPB_INTERFACE_RETURN_EXP( mtxAp )
}
}
')dnl

#ifdef __cplusplus
}
#endif  /* __cplusplus */

ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_LIBSPBLAS_H_INCLUDED */
')

ifdef(`ONLY_WANT_HEADERS',`',`dnl
dnl /* @endcond */
')dnl
dnl
