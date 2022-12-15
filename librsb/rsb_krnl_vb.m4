dnl
dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
dnl
/**
 @file
 @brief

 Performance kernels dispatching code, for each type, submatrix size, operation.
 For variable block partitioned matrices.

 */
dnl
include(`rsb_misc.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
include(`rsb_krnl_vb_macros.m4')dnl
dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl				Function definitions
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_VBR_X_H_INCLUDED
#define RSB_VBR_X_H_INCLUDED 
')dnl
dnl 
dnl
dnl
#include "rsb_internals.h"
dnl
dnl	VARIABLE BLOCK SIZE KERNELS :
dnl	:
dnl foreach(`matrix_storage',(`vbr',`block',WANT_MATRIX_VB_STORAGE),`dnl
dnl
/* Function definitions */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`matrix_storage',(WANT_MATRIX_VB_STORAGE),`dnl
foreach(`unrolling',(`l',`u'),`dnl
RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION(type,matrix_storage,,,unrolling,mop)dnl
')dnl
')dnl
')dnl
')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#endif /* RSB_VBR_X_H_INCLUDED */
')dnl
dnl
/* @endcond */
dnl
