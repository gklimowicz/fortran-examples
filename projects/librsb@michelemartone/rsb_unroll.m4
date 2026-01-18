dnl
dnl	@author: Michele Martone
dnl
dnl	execute this script with M4 to obtain a loop unrolled functions collection
dnl	this `forloop' macro is the one in the ./examples  directory distributed with the M4 package
dnl	TODO : To use goto's and labels to code much less boundary loops!
dnl	TODO : Introduce variables for indices inside for loops: do not rely on the compiler from a point on...
dnl
dnl
dnl
dnl	The generated code will expand from here
dnl
dnl
/* @cond INNERDOC */
dnl
/**
 * @file
 * @brief
 * Unrolled kernels, for each type, operation, submatrix.
 * Right now, they are used for VBR and alike formats.
 */
dnl
/* Take care of compiling this code without loop unrolling optimizations (-fno-unroll-loops, or -ON with N<=2 on gcc) */
include(`rsb_misc.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_UNROLL_H_INCLUDED
#define RSB_UNROLL_H_INCLUDED
')
dnl 
dnl
dnl
dnl
dnl
#include "rsb_types.h"
#include "rsb_common.h"
include(`do_unroll.m4')dnl
dnl

/* NULL should be defined. */
#ifndef NULL
#define NULL ((int*)(0))
#endif /* NULL */
dnl
dnl
ifelse(dnl
RSB_M4_MEMBER(`VBR',WANT_MATRIX_STORAGE)dnl
RSB_M4_MEMBER(`VBC',WANT_MATRIX_STORAGE)dnl
,`00',`dnl
/**
 * No VBR/VBC formats compiled in.
 */
',`dnl
dnl

/* FIXME : we only want this code if VB and L stuff is in */
#if defined(RSB_MATRIX_STORAGE_LC) || defined(RSB_MATRIX_STORAGE_LR) || defined(WANT_MATRIX_VB_STORAGE)


/*!
 * This code instance has coverage for submatrices of sizes as in the 
 * cartesian (rows) x (columns) product of RSB_M4_ROWS_UNROLL x RSB_M4_COLUMNS_UNROLL.
 * 
 *  For each submatrix operation, this code instance offers a macro for
 *  dispatching a function pointer to the specialized function, regarding
 *  rows and columns unrolling factors.
 *
 * In case a submatrix of a size not in the above cartesian product is 
 * given, the dispatching macros will assign a pointer to a function
 * unrolled RSB_M4_ROWS_FALLBACK_UNROLL times on the rows and RSB_M4_COLUMNS_FALLBACK_UNROLL times on the columns.
 */

/* Function headers */

ifdef(`ONLY_WANT_HEADERS',`dnl
dnl unrollm(`mrowsu',`Mrowsu',`mcolsu',`Mcolsu',
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`unrolling',(`u',`l'),`dnl
RSB_M4_UNROLL_KERNEL(`row',`rows',rowsu,`column',`columns',colsu,`type',`h',mop,unrolling,transposition)dnl
')')')')')')')dnl
')dnl

dnl
dnl	
ifdef(`ONLY_WANT_HEADERS',,`
/* Function definitions */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`unrolling',(`u',`l'),`dnl
RSB_M4_UNROLL_KERNEL(`row',`rows',rowsu,`column',`columns',colsu,`type',,mop,unrolling,transposition)
')')')')')')') ')dnl


dnl
dnl	Function dispatch table name generating macros
dnl
pushdef(`rsb_mx_dispatch_table_name',`dnl
pushdef(`type',$1)dnl
pushdef(`unrolling',ifelse($2,`l',`l',`'))dnl FIXME : this is a temporary fix (setting to `' instead of `u')
pushdef(`optype',$3)dnl
void(*RSB_M4_PREFIX`'optype`_'RSB_M4_TYPE_CODE(type)`_'unrolling`_pointer_table'[])
RSB_M4_KERNEL_FUNCTION_ARGS(type,unrolling,optype)
popdef(`unrolling')dnl
popdef(`optype')dnl
popdef(`type')
')dnl

dnl popdef(`rsb_mx_dispatch_table_name')
dnl
dnl
')dnl	ifelse(RSB_M4_MEMBER..
dnl

dnl
dnl	Extern declarations
dnl
dnl foreach(`mop',(m,spmv_uauz),
dnl foreach(`looped',`(l,)',
dnl ifdef(`ONLY_WANT_HEADERS',foreach(`type',RSB_M4_MATRIX_TYPES, ` extern rsb_mx_dispatch_table_name(type,looped,mop); /* Dispatch table */ '))))

dnl
dnl	Function dispatch table names generating macros
dnl

dnl dnl	THE FOLLOWING WAS COMMENTED BY HAND BECAUSE IT WAS TROUBLESOME TO USE ONLY M4 TO DEACTIVATE IT.
dnl dnl THEREFORE, WHEN YOU RE-ACTIVATE THE VBR STUFF, YOU SHOULD UN-COMMENT THE FOLLOWING BLOCK! (FIXME, TODO)
dnl dnl
dnl dnl
dnl dnl	leave glombo where it is : it is a dummy parameter, or 'l' handling macros won't work :)
dnl foreach(`mop',RSB_M4_MATRIX_OPS,`
dnl foreach(`looped',(`l',`u'),`
dnl foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
dnl ifdef(`ONLY_WANT_HEADERS',,`dnl
dnl foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
dnl rsb_mx_dispatch_table_name(type,looped,mop)dnl
dnl ={ /* Dispatch table compilation */
dnl foreach(`rowsu',RSB_M4_ROWS_UNROLL, `foreach(`colsu',RSB_M4_COLUMNS_UNROLL,
dnl `	RSB_M4_KERNEL_FUNCTION_NAME(type,rowsu,colsu,looped,mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE),
dnl ')')	(void*)NULL
dnl };
dnl ') ') ') ') ')
dnl
dnl	

dnl maximum element in a list
dnl : FIXME : these macros are also defined in sort.m4. you should remove this redundance!
define(`max2',`ifelse(eval(`$1>$2'),1,`$1',`$2')')dnl
define(`maxn',`ifelse($#,1,$1,`ifelse($#,2,`max2($1,$2)',`max2($1,maxn(shift($@)))')')')dnl
define(`min2',`ifelse($2,,$1,`ifelse(eval(`$1<$2'),1,`$1',`$2')')')dnl
define(`minn',`ifelse($2,,$1,`ifelse($#,2,`min2($1,$2)',`min2($1,minn(shift($@)))')')')dnl

/**
 * Loops unroll factors.
 */
`#define RSB_MIN_ROW_UNLOOP_FACTOR'	minn(WANT_ROW_UNLOOP_FACTORS)
`#define RSB_MAX_ROW_UNLOOP_FACTOR'	maxn(WANT_ROW_UNLOOP_FACTORS)
`#define RSB_MIN_COLUMN_UNLOOP_FACTOR'	minn(WANT_COLUMN_UNLOOP_FACTORS)
`#define RSB_MAX_COLUMN_UNLOOP_FACTOR'	maxn(WANT_COLUMN_UNLOOP_FACTORS)

dnl
dnl please automate generation of such macros
dnl
dnl #if 0
dnl #define RSB_GET_MV_DOUBLE_KERNEL(rows,columns) \
dnl 	( \
dnl 	(columns)>7 \
dnl 		? \
dnl 		( \
dnl 			(rows) > 2 ? \
dnl 			rsb_mv_double_r4_c8_l: \
dnl 			rsb_mv_double_r1_c8_l \
dnl 		): \
dnl 		( \
dnl 			(columns)>1? \
dnl 			rsb_mv_double_r4_c8_l: \
dnl 			rsb_mv_double_r4_c4_l \
dnl 		) \
dnl 	)
dnl //#else
dnl #define RSB_GET_MV_DOUBLE_KERNEL(rows,columns) \
dnl ( \
dnl 			rsb_mv_double_r4_c8_l \
dnl )
dnl #endif
dnl	
dnl
dnl
dnl
dnl
ifelse(dnl
RSB_M4_MEMBER(`VBR',WANT_MATRIX_STORAGE)dnl
RSB_M4_MEMBER(`VBC',WANT_MATRIX_STORAGE)dnl
,`00',`dnl
/**
 * No VBR/VBC formats compiled in.
 */
',`dnl
dnl
dnl
dnl
#if 1
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`looped',(,),`dnl (,) or (l)
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
pushdef(`fix', `RSB_M4_KERNEL_FUNCTION_NAME(type,rowsu,colsu,looped,mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE) \
' )
DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_(mop,type,looped,`UNLOOP_R_C_PAIRS`RSB_M4_KERNEL_FUNCTION_NAME(type,RSB_M4_ROWS_FALLBACK_UNROLL,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE)'')
popdef(`fix')
')')')')
dnl
#endif
dnl	TODO : implement a macro choosing a suboptimal kernel function
dnl
dnl #undef RSB_double_spmv
dnl #define RSB_double_spmv(R,C) rsb_mv_double_r4_c4

/* FIXME : we only want this code if VB and L stuff is in */
#endif
#
dnl
')dnl	ifelse(RSB_M4_MEMBER..
dnl


ifdef(`ONLY_WANT_HEADERS',`dnl
#endif  /* RSB_UNROLL_H_INCLUDED */
')

ifdef(`ONLY_WANT_HEADERS',`',`dnl
RSB_EMPTY_FILE_FILLER
')dnl

ifdef(`ONLY_WANT_HEADERS',`
`#define RSB_FITTING_SAMPLES		'dnl
RSB_M4_FITTING_SAMPLES
')dnl
dnl
/* @endcond */
dnl

