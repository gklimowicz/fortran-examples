dnl
dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
dnl
/*! 
 @file
 @brief Matrix type dispatching code, for each matrix operation.
 */
dnl
include(`rsb_misc.m4')dnl
dnl
RSB_M4_HEADER_MESSAGE()dnl
dnl
include(`rsb_krnl_bcss_macros.m4')dnl
include(`rsb_krnl_linked_lists.m4')dnl
include(`rsb_krnl_macros.m4')dnl
dnl
dnl
include(`rsb_misc.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
include(`rsb_krnl_linked_lists.m4')dnl
dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl				Function definitions
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_LL_H_INCLUDED
#define RSB_LL_H_INCLUDED
')dnl
dnl 
dnl
ifelse(dnl
RSB_M4_MEMBER(`LR',WANT_MATRIX_STORAGE)dnl
RSB_M4_MEMBER(`LC',WANT_MATRIX_STORAGE)dnl
,`00',`dnl
/**
 * No linked lists formats compiled in.
 */
#include "rsb_internals.h"
RSB_EMPTY_FILE_FILLER 
',`dnl
dnl
dnl
#include "rsb_internals.h"
dnl
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`matrix_storage',(`LR',`LC'),`dnl
foreach(`unrolling',(`l',`u'),`dnl
dnl RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(want_what,type,matrix_storage,unrolling,,,mop,citype,diagonal,uplo)
RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION(type,matrix_storage,unrolling,mop)
')dnl
')dnl
')dnl
')dnl
dnl
dnl
')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#endif /* RSB_LL_H_INCLUDED */
')dnl
dnl
/* @endcond */
dnl
