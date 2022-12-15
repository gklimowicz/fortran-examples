dnl
dnl
dnl	@author: Michele Martone
dnl
dnl
/* @cond INNERDOC */
/*! 
 @file
 @brief Matrix type dispatching code, for each matrix operation.
 */
dnl
include(`rsb_misc.m4')dnl
dnl
RSB_M4_HEADER_MESSAGE()dnl

ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_DISPATCH_H_INCLUDED
#define RSB_DISPATCH_H_INCLUDED
')dnl

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

dnl
include(`rsb_krnl_bcss_macros.m4')dnl
include(`rsb_krnl_linked_lists.m4')dnl
include(`rsb_krnl_macros.m4')dnl
dnl
dnl #include "rsb_internals.h"
#include "rsb_common.h"
dnl #include "rsb_krnl_vb.h"	/* uhm */
dnl #include "rsb_krnl_lb.h"	/* uhm */
dnl #include "rsb_krnl_bcss.h"	/* uhm */
dnl #include "rsb_krnl_bcss_u.h"	/* uhm */
dnl #include "rsb_krnl_bcss_l.h"	/* uhm */
#include "rsb_krnl_bcss_spmv_u.h"	/* uhm */
#include "rsb_krnl_bcss_spsv_u.h"	/* uhm */
#include "rsb_krnl_bcss_misc_u.h"	/* uhm */
dnl
dnl

ifdef(`ONLY_WANT_HEADERS',`dnl
dnl
',`dnl
#pragma GCC visibility push(hidden)
')dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
#if RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS
#define RSB_MTX_NOT_OK(MTXAP) ( ((MTXAP)->rpntr || (MTXAP)->cpntr) || RSB_MATRIX_UNSUPPORTED_TYPE((MTXAP)->typecode) )
#else /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
#define RSB_MTX_NOT_OK(MTXAP) RSB_MATRIX_UNSUPPORTED_TYPE((MTXAP)->typecode)
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
dnl
')dnl
dnl

ifdef(`ONLY_WANT_HEADERS',`dnl
#define	RSB_BCSR_GET_NEXT_BLOCK_POINTER(BP,mtxAp,ROWVAR,COLVAR,BLOCKROWSPAR,BLOCKCOLSPAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	/*										\
	 * *input*									\
	 * mtxAp		should be a valid rsb_mtx_t structure pointer		\
	 * BLOCKROWSPAR	should be set to the rows   count of this block			\
	 * BLOCKCOLSPAR	should be set to the column count of this block			\
	 * *output*									\
	 * ROWVAR	will be set to the base row    of this block			\
	 * COLVAR	will be set to the base column of this block			\
	 * BP		will be set to the current block pointer			\
	 * */										\
	while( (mtxAp)->bpntr[_i] == (mtxAp)->bpntr[_i+1] ) 				/* skipping empty rows */	\
	{++_i;_k=(mtxAp)->bpntr[_i];} 		/* _k is the first block index for the current row of blocks */	\
	_j=(mtxAp)->bindx[_k]; 						/* the current block column index  */	\
	_lastk=_k;	\
	(BLOCKROWVAR)=_i;	\
	(BLOCKCOLUMNVAR)=_j;	\
	(ROWVAR)=(BLOCKROWSPAR)*_i;					/* _i is the current block row index */	\
	(COLVAR)=(BLOCKCOLSPAR)*_j; 					/* the current block column index  */	\
	BP+=(mtxAp)->options->el_size*(BLOCKROWSPAR)*(BLOCKCOLSPAR);			\
	_k++; 		/* for the future macro calls */						\
	if( _k >= (mtxAp)->bpntr[_i+1] )++_i;								\
	;

#define RSB_BCSR_GET_FIRST_BLOCK_POINTER(BP,mtxAp,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	int _i=0,_j=0,_k=0,_lastk=0;									\
	(BLOCKROWSVAR)=(mtxAp)->rpntr[1]-(mtxAp)->rpntr[0];		/* _i is the current block row index */	\
	(BLOCKCOLSVAR)=(mtxAp)->cpntr[1]-(mtxAp)->cpntr[0]; 		/* the current block column index  */	\
	(BP)=(mtxAp)->VA;											\
	RSB_BCSR_GET_NEXT_BLOCK_POINTER(BP,mtxAp,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)

#if RSB_WANT_DBC
#define RSB_BCSR_GOT_LAST_BLOCK_POINTER(mtxAp)	( _lastk >= (mtxAp)->block_count )
#endif
')dnl

ifdef(`ONLY_WANT_HEADERS',`
`#define RSB_BENCHMARK_MIN_SECONDS	'dnl
RSB_M4_BENCHMARK_MIN_SECONDS
`#define RSB_BENCHMARK_MIN_RUNS		'dnl
RSB_M4_BENCHMARK_MIN_RUNS 
')dnl
dnl


foreach(`mop',RSB_M4_MATRIX_ALL_OPS,`dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION(RSB_M4_MATRIX_TYPES,mop)
')
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_TESTING_FUNCTION(RSB_M4_MATRIX_TYPES,mop)
RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION(RSB_M4_MATRIX_TYPES,mop)
')
dnl
dnl
dnl	FIXME : still not for transposed kernels
dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`mtype',RSB_M4_MATRIX_TYPES,`dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION(mop,mtype,`function_declaration')
',`dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION(mop,mtype,`function_definition')
')dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION(RSB_M4_MATRIX_TYPES,mop)
')dnl
')dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`
/* Dispatch table for type and mop specific benchmarks */
 rsb_err_t (* rsb_benchmark_dispatch_table [RSB_IMPLEMENTED_TYPES][RSB_IMPLEMENTED_MOPS]) 
   (RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION(,void,`function_args')`_pointer_table)'
` =  {' 
foreach(`mtype',RSB_M4_MATRIX_TYPES,`dnl
{
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
 RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION(mop,mtype,`function_identifier')`'dnl
ifelse(mop,RSB_M4_LAST_LIST_ELEMENT(WANT_MATRIX_OPS),` ',`,')
')dnl
}
ifelse(mtype,RSB_M4_LAST_LIST_ELEMENT(WANT_TYPES),` ',`,')
')dnl
dnl 	the following breaks xlc:
dnl        (void*)NULL	/* FIXME : is not this overflow declaration ? */
};
')dnl
dnl
dnl
foreach(`mtype',RSB_M4_MATRIX_TYPES,`dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION(mtype)
')dnl
dnl
dnl
dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
#if 0
ifdef(`ONLY_WANT_HEADERS',`dnl
RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(RSB_M4_MATRIX_TYPES,mop,`function_declaration')
',`dnl
RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(RSB_M4_MATRIX_TYPES,mop,`function_definition')
')dnl
#endif /* 0 */
ifdef(`ONLY_WANT_HEADERS',`dnl
',`dnl
#pragma GCC visibility push(default)
')dnl
RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION(mop)
')dnl
dnl
dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION`'dnl
dnl
dnl
RSB_M4_DUMP_PERFORMANCE_INFO_ARRAY_FUNCTION()
dnl
#ifdef __cplusplus
}
#endif  /* __cplusplus */

ifdef(`ONLY_WANT_HEADERS',`
#endif	/* RSB_DISPATCH_H_INCLUDED */
')



dnl
dnl NEW : FIXME
/*!
 @file
 @brief ...
 */
dnl
/* @endcond */
dnl
