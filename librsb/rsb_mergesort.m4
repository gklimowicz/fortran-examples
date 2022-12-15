dnl
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
include(`mergesort_macros.m4')dnl
dnl
/* @cond INNERDOC */
dnl
/**
 * @file
 * @brief
 * Sorting functions.
 */
RSB_M4_HEADER_MESSAGE()dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_MERGESORT_H_INCLUDED
#define RSB_MERGESORT_H_INCLUDED
')
dnl
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

ifdef(`ONLY_WANT_HEADERS',`dnl
',`dnl
dnl
dnl /* We may use custom memcpy functions. */
dnl #define RSB_MEMCPY(DST,SRC,BYTES) rsb__memcpy((DST),(SRC),(BYTES))
')dnl
dnl

dnl
RSB_M4_INCLUDE_HEADERS
dnl 

#if RSB_OBSOLETE_QUARANTINE

ifdef(`ONLY_WANT_HEADERS',`',`dnl
RSB_INTERNALS_COMMON_HEAD_DECLS
')dnl

dnl
define(`blockorientations',`(CSR,BCSR,VBR)')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
foreach(`blockoriented',blockorientations,`dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_PROTOTYPE(RSB_M4_TYPES,blockoriented);
')dnl
')dnl
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
foreach(`blockoriented',blockorientations,`dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER(RSB_M4_TYPES,blockoriented)
')dnl
')dnl
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
foreach(`blockoriented',blockorientations,`dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION(mtype,blockoriented)
RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION(mtype,blockoriented)
')dnl
')dnl
')dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
foreach(`blockoriented',blockorientations,`dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_PROTOTYPE(mtype,blockoriented);
RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_PROTOTYPE(mtype,blockoriented);
')dnl
')dnl
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE */

dnl
#ifdef __cplusplus
}
#endif  /* __cplusplus */
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_MERGESORT_H_INCLUDED */
')
dnl
/* @endcond */
dnl
