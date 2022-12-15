/* @cond INNERDOC */
/*! 
 @file
 @brief 

 Matrix Operations testing code source file.
 This is NOT part of the library: only of companion programs.

 */

/*

Copyright (C) 2008-2022 Michele Martone

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
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

 */
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block coordinates format.
 Kernels unrolled, with no loops, for only user-specified blockings.
 */

/*

Copyright (C) 2008-2022 Michele Martone

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
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

 */

#ifndef RSB_TEST_MATOPS_H_INCLUDED
#define RSB_TEST_MATOPS_H_INCLUDED

/* FIXME: necessary, until we use so many #ifdefs in this program */
#include "rsb-config.h"
#include "rsb_common.h"
#include "rsb_mkl.h"
#if RSB_WANT_ARMPL
#include <armpl.h>
#endif /* RSB_WANT_MKL */

#if RSB_WITH_LIKWID
#include <likwid.h>
#define RSB_LIKWID_MARKER_INIT	{RSBENCH_STDOUT("# Initializing the LIKWID API with likwid_markerInit().\n");likwid_markerInit();}
#define RSB_LIKWID_MARKER_EXIT {RSBENCH_STDOUT("# Finalizing the LIKWID API with likwid_markerClose().\n");likwid_markerClose();}
#define RSB_LIKWID_MARKER_R_START(R) likwid_markerStartRegion(R)
#define RSB_LIKWID_MARKER_R_STOP(R) likwid_markerStopRegion(R)
#else /* RSB_WITH_LIKWID */
#define RSB_LIKWID_MARKER_INIT
#define RSB_LIKWID_MARKER_EXIT
#define RSB_LIKWID_MARKER_R_START(R)
#define RSB_LIKWID_MARKER_R_STOP(R)
#endif /* RSB_WITH_LIKWID */

#define RSB_STDOUT_FD stdout
#define RSBENCH_STDOUT( ... ) fprintf(RSB_STDOUT_FD, __VA_ARGS__ )
#define RSB_WAT_FMT_H "Ss[Xx[Tt[V[V]]]]"
#define RSB_WAT_FMT "%lfs[%dx[%dt[%c]]]"

#define RSB_MIN_ABOVE_INF(X,Y,MIN) RSB_MAX(RSB_MIN(X,Y),MIN)
#define RSB_INT_MILLION 1000000
#define RSB_REAL_MILLION 1000000.0 
enum rsb_dumpvec_enum { rsb_dumpvec_no= 0, rsb_dumpvec_res= 1, rsb_dumpvec_rhs= 2 };

#define RSB_HAVE_METIS 0 /* FIXME: unfinished */
#if RSB_HAVE_METIS
#include <metis/metis.h>
#endif /* RSB_HAVE_METIS */

#ifdef RSB_WANT_OSKI_BENCHMARKING 
#ifdef RSB_HAVE_OSKI_OSKI_H 
#include <oski/oski.h>
#else /* RSB_HAVE_OSKI_OSKI_H */
#error "you should disable oski benchmarking at configure time!"
#endif /* RSB_HAVE_OSKI_OSKI_H */
#endif /* RSB_WANT_OSKI_BENCHMARKING */
#ifdef RSB_HAVE_UNISTD_H
#include <unistd.h>
#endif /* RSB_HAVE_UNISTD_H */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#ifdef RSB_WANT_PERFORMANCE_COUNTERS
	#define RSB_PERFORMANCE_COUNTERS_DUMP_MEAN(MSG,PMSG,TIMES,PCIP) if(want_perf_counters>0){rsb_perf_counters_update(); if(PMSG)rsb_perf_counters_dump(MSG,NULL,TIMES,PCIP); rsb_perf_counters_reset();/* TEMPORARY */}
	#define RSB_PERFORMANCE_COUNTERS_DUMP(MSG,PMSG) if(want_perf_counters>1)RSB_PERFORMANCE_COUNTERS_DUMP_MEAN(MSG,PMSG,1,NULL) 
#else /* RSB_WANT_PERFORMANCE_COUNTERS */
	#define RSB_PERFORMANCE_COUNTERS_DUMP_MEAN(MSG,PMSG,TIMES,PCIP) 
	#define RSB_PERFORMANCE_COUNTERS_DUMP(MSG,PMSG)
#endif /* RSB_WANT_PERFORMANCE_COUNTERS */

#if RSB_WITH_LIKWID
#define RSB_TM_LIKWID_MARKER_R_START(R) if(want_likwid == RSB_BOOL_TRUE)RSB_LIKWID_MARKER_R_START(R)
#define RSB_TM_LIKWID_MARKER_R_STOP(R)  if(want_likwid == RSB_BOOL_TRUE)RSB_LIKWID_MARKER_R_STOP(R)
#else
#define RSB_TM_LIKWID_MARKER_R_START(R)
#define RSB_TM_LIKWID_MARKER_R_STOP(R)
#endif /* RSB_WITH_LIKWID */


#if defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS==1)
#define RSB_WANT_PERFORMANCE_COUNTERS_IN_RSBENCH 1
#else
#define RSB_WANT_PERFORMANCE_COUNTERS_IN_RSBENCH 0
#endif





int rsb_test_help_and_exit(const rsb_char_t *argv0, const rsb_option_t *o, int code);

#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	#define RSBENCH_MEM_ALLOC_INFO(LT) RSBENCH_STDOUT("( allocated_memory:%zd allocations_count:%zd allocations_cumulative:%zd)%s",rsb_global_session_handle.allocated_memory,rsb_global_session_handle.allocations_count,rsb_global_session_handle.allocations_cumulative,LT);
#else
	#define RSBENCH_MEM_ALLOC_INFO(LT)
#endif

#if RSB_WANT_ARMPL



rsb_err_t rsb__armpl_csr_spmv_bench(const void *VA, const armpl_int_t m, const armpl_int_t k, const armpl_int_t nnz, const armpl_int_t * IP, const armpl_int_t *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp,rsb_bool_t want_at);
#endif /* RSB_WANT_ARMPL */

/* one function for each of (spmv_sxsa,spsv_sxsx,mat_stats)*/
int rsb__main_block_partitioned_spmv_sxsa(const int argc, rsb_char_t * const argv[])
;
int rsb__main_block_partitioned_spsv_sxsx(const int argc, rsb_char_t * const argv[])
;
int rsb__main_block_partitioned_mat_stats(const int argc, rsb_char_t * const argv[])
;

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif	/* RSB_TEST_MATOPS_H_INCLUDED */

/* @endcond */
