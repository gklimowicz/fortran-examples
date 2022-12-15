dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
dnl
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block compressed sparse stripes format.
 Kernels unrolled, with no loops, for only user-specified blockings.
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
#ifndef RSB_BCSS_MISC_U_H_INCLUDED
#define RSB_BCSS_MISC_U_H_INCLUDED
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
#include "rsb.h"
dnl
ifdef(`ONLY_WANT_HEADERS',`',`
#pragma GCC visibility push(hidden)
')dnl
dnl
RSB_M4_BCSS_MISC_KERNELS((`u'))
dnl
')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#endif /* RSB_BCSS_MISC_U_H_INCLUDED */
')dnl
dnl
/* @endcond */
dnl
