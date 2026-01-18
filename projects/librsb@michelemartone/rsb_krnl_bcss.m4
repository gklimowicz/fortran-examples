dnl
dnl	@author: Michele Martone
dnl
dnl
/* @cond INNERDOC */
dnl
/*!
 @file
 @brief

 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block compressed sparse stripes format.
 */
dnl
include(`rsb_misc.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
include(`rsb_krnl_bcss_macros.m4')dnl
include(`rsb_krnl_vb_macros.m4')dnl FIXME : RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME
dnl
dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl				Function definitions
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_BCSS_H_INCLUDED
#define RSB_BCSS_H_INCLUDED
')dnl
dnl 
dnl
ifelse(dnl
RSB_M4_MEMBER(`BCSR',WANT_MATRIX_STORAGE)dnl
RSB_M4_MEMBER(`BCSC',WANT_MATRIX_STORAGE)dnl
,`00',`dnl
/**
 * No BCSS formats compiled in.
 */
',`dnl
dnl
dnl
dnl
dnl
#include "rsb_internals.h"
dnl
ifelse(RSB_M4_WANT_OMP_IN_KERNELS,`1',`dnl
#include <omp.h>	/* OpenMP parallelism (EXPERIMENTAL) */
')dnl
dnl
dnl
dnl Now the following macros are split across files.
dnl RSB_M4_BCSS_KERNELS((`l',`u'))
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
foreach(`unrolling',(`l',`u'),`dnl
foreach(`symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`mtype',RSB_M4_MATRIX_TYPES,`dnl
DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH_MACRO_(mop,mtype,unrolling,matrix_storage,transposition,symmetry,`UNLOOP_R_C_PAIRS`RSB_M4_BCSS_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,RSB_M4_ROWS_FALLBACK_UNROLL,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop)dnl
'')
')dnl
DOUBLE_LINEAR_KERNEL_DISPATCHER_TYPE_SEARCH_MACRO_(mop,unrolling,matrix_storage,transposition)
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
dnl
dnl
')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#endif /* RSB_BCSS_H_INCLUDED */
')dnl
dnl
/* @endcond */
dnl
