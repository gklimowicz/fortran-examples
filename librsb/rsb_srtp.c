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
/* @cond INNERDOC  */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains parallel sorting functions.
 * */

#include "rsb_common.h"

#define RSB_DO_WANT_PSORT_VERBOSE 0	/* set this to >0 to print parallel sort statistics */
#define RSB_DO_WANT_PSORT_TIMING (RSB_DO_WANT_PSORT_VERBOSE+0) 	/* set this to 1 to print some statistics */
#define RSB_DO_WANT_PSORT_FASTER_BUT_RISKY 0	/* FIXME: STILL UNFINISHED  */
#define RSB_WANT_SORT_PARALLEL_BUT_SLOW 1
#define RSB_IS_ENOUGH_NNZ_FOR_PARALLEL_SORT(NNZ) ((NNZ)>=RSB__ENOUGH_NNZ_FOR_PARALLEL_SORT)

RSB_INTERNALS_COMMON_HEAD_DECLS

rsb_err_t rsb__util_sort_row_major_parallel(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t nr, rsb_coo_idx_t nc,  rsb_type_t typecode, rsb_flags_t flags)
{
	/**
	 * \ingroup gr_util
	 * TODO: should describe somewhere our technique: this is a mixed counting sort + merge sort.
	 * 
	 * FIXME: it is known that very small matrices (e.g.: 2x2 from `make tests`) are not handled, here.
	 * require RSB_IS_ENOUGH_NNZ_FOR_PARALLEL_SORT(nnz).
	*/

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	int cfi;
	const size_t el_size = RSB_SIZEOF(typecode);
	const float cfa[] = {2}; // will use threads count as a reference
//	const float cfa[] = {1}; // will use cache size as a reference
//	const float cfa[] = {-1};  // will use as much memory as possible
	const long wet = rsb_get_num_threads(); /* want executing threads */

	if( nr > RSB_MAX_MATRIX_DIM || nc > RSB_MAX_MATRIX_DIM || nnz > RSB_MAX_MATRIX_NNZ )
	{
		errval = RSB_ERR_LIMITS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(!RSB_IS_ENOUGH_NNZ_FOR_PARALLEL_SORT(nnz))
	{
		//errval = RSB_ERR_BADARGS;
		//RSB_PERR_GOTO(err,RSB_ERRM_ES);
		return rsb__util_sort_row_major_inner(VA, IA, JA, nnz, nr, nc, typecode , flags);
	}

	if(!IA || !JA || !VA)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	for(cfi=0;cfi<sizeof(cfa)/sizeof(float);++cfi)
	{
		void *W = NULL;
		rsb_char_t *IW = NULL;
		size_t ti;
		long cs,tc,tcs;
		rsb_nnz_idx_t cnnz = 0;
		size_t bnnz = 0;
		size_t wb = 0;
		const size_t fsm = rsb__sys_free_system_memory();
		const long ns=2*sizeof(rsb_coo_idx_t)+el_size;
		const long bs=nnz*ns;

#if RSB_DO_WANT_PSORT_TIMING
		rsb_time_t dt,st,mt,tt;
#endif /* RSB_DO_WANT_PSORT_TIMING */

		// compute how many bytes are necessary for an element
		// compute how many bytes are necessary for the whole processing
		if(cfa[cfi]>1)
			tcs=bs+(wet*ns);
		else if(cfa[cfi]>0)
			tcs = rsb__get_lastlevel_c_size()*cfa[cfi]*wet;/* FIXME: '*wet' is a hack just for benchmark-related issues */
		else
			tcs=fsm/2; /* could be 0 */

		if(tcs<1)
			tcs=bs;	/* could happen, for an interfacing problem */
		else
		if(fsm>0)
		{
			tcs = RSB_MIN(fsm,tcs);
		}

		// prepare a buffer
		W = rsb__malloc(tcs);

		cs=tcs/wet;
		//RSB_INFO("cache is %d bytes\n",cs);

		// compute nnz count fitting in the buffer
		cnnz=cs/ns;
		// compute total count of necessary passes
		//tc=(bs+cs-1)/cs;
		tc=(nnz+(cnnz-1))/(cnnz);
		
		wb = RSB_DO_REQUIRE_BYTES_FOR_INDEX_BASED_SORT(cnnz,nr,nc,1,1);
		IW = rsb__malloc(wb*wet);

		if(!IW || !W)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(erri,RSB_ERRM_ES);
		}	

		//RSB_INFO("there are %d nnz (%d bytes), %d times the cache (%z bytes), %d nnz per cache\n",nnz,bs,tc,cs,cnnz);
#if RSB_DO_WANT_PSORT_TIMING
		dt = rsb_time();
		st=-dt;
#endif /* RSB_DO_WANT_PSORT_TIMING */

		// this phase is potentially parallel, and is the slower one.
		// NOTE:before parallelization, one should avoid allocations during sort, or serialize in them in some way! 
		#pragma omp parallel for schedule(static,1) shared(IA,JA,VA) reduction(|:errval)   RSB_NTC  
		for(ti=0;ti<tc;++ti)
		{
			const size_t fnnz=ti*cnnz;
			const rsb_nnz_idx_t bnnz=(ti==tc-1)?(nnz-fnnz):cnnz;
//			RSB_INFO("s:%d..%d (bnnz=%d)\n",fnnz,fnnz+bnnz-1,bnnz);
			errval = rsb__util_sort_row_major_buffered(((rsb_byte_t*)VA)+el_size*fnnz,IA+fnnz,JA+fnnz,bnnz,nr,nc,typecode,flags,IW+wb*ti,wb);
			//assert(!rsb__util_is_sorted_coo_as_row_major(IA+fnnz,JA+fnnz,bnnz,typecode,NULL,flags));
			RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IA+fnnz,JA+fnnz,bnnz,typecode,NULL,flags));
		}

#if RSB_DO_WANT_PSORT_TIMING
		dt = rsb_time();
		st+=dt;
#endif /* RSB_DO_WANT_PSORT_TIMING */
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}

		// this phase is potentially parallel, too
		for(bnnz=cnnz;bnnz<nnz;)
		{
		//	size_t fnnz;
			rsb_nnz_idx_t fi, fn = ((nnz-bnnz)+(2*bnnz-1))/(2*bnnz);
			#pragma omp parallel for schedule(static,1) shared(IA,JA,VA)   RSB_NTC 
//			for(fnnz=0;fnnz<nnz-bnnz;fnnz+=2*bnnz)
			for(fi=0;fi<fn;++fi)
			{
#if RSB_WANT_OMP_RECURSIVE_KERNELS
				rsb_char_t * lW=((rsb_char_t*)W)+cs*omp_get_thread_num();
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
				rsb_char_t * lW=((rsb_char_t*)W)+cs*0;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
				const size_t fnnz=fi*2*bnnz;
				const size_t lnnz=(fnnz+2*bnnz>nnz)?(nnz-(fnnz+bnnz)):bnnz;
				void *fVA=RSB_VA_OFFSET_POINTER(VA,el_size,fnnz);
				void *fIA=IA+fnnz;
				void *fJA=JA+fnnz;
#if RSB_PS_ASSERT
				void *bVA=RSB_VA_OFFSET_POINTER(VA,el_size,(fnnz+bnnz));
				void *bIA=IA+fnnz+bnnz;
				void *bJA=JA+fnnz+bnnz;
#endif /* RSB_PS_ASSERT */
				//RSB_INFO("nr:%d..%d %d..%d (%d) (bnnz=%d) (lnnz=%d)\n",fnnz,fnnz+bnnz-1,fnnz+bnnz,fnnz+bnnz+lnnz-1,cs,bnnz,lnnz);
//				RSB_INFO("sentinel:%x %d %d\n",IA+fnnz+bnnz+lnnz,IA[fnnz+bnnz+lnnz],JA[fnnz+bnnz+lnnz]);
				RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(fIA,fJA,bnnz,typecode,NULL,flags));
				RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(bIA,bJA,lnnz,typecode,NULL,flags));
				rsb__do_util_merge_sorted_subarrays_in_place(fVA,fIA,fJA,lW,bnnz,lnnz,cs,flags,typecode);
//				RSB_INFO("sentinel:%x %d %d\n",IA+fnnz+bnnz+lnnz,IA[fnnz+bnnz+lnnz],JA[fnnz+bnnz+lnnz]);
				RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(fIA,fJA,bnnz,typecode,NULL,flags));
				RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(fIA,fJA,bnnz+lnnz,typecode,NULL,flags));
			}
			#pragma omp barrier
			bnnz *= 2;
		}
#if RSB_DO_WANT_PSORT_TIMING
		mt = - dt;
		dt = rsb_time();
		mt += dt;
		tt = mt + st;
		RSB_INFO("using %zd partitions, (sort=%.5lg+merge=%.5lg)=%.5lg, on %d threads\n",(size_t)tc,st,mt,tt,wet);
#endif /* RSB_DO_WANT_PSORT_TIMING */

//		assert(!rsb__util_is_sorted_coo_as_row_major(IA,JA,nnz,typecode,NULL,flags));
//		if(rsb__util_is_sorted_coo_as_row_major(IA,JA,nnz,typecode,NULL,flags))
//			RSB_PERR_GOTO(err,RSB_ERRM_EM);

erri:
		RSB_CONDITIONAL_FREE(W);
		RSB_CONDITIONAL_FREE(IW);
	}
err:
	RSB_DO_ERR_RETURN(errval)
}


rsb_err_t rsb__util_sort_row_major_bucket_based_parallel(void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const  rsb_type_t typecode , const rsb_flags_t flags /*, void * WA, size_t wb */)
{
	/**
		\ingroup gr_internals
		Tolerate duplicates
	       	So tolerate e.g. nr*nc<nnz, arising in sparse matrices' sum.
	*/
	rsb_nnz_idx_t psc = RSB_PSORT_CHUNK;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t mnzpr = 0;
	rsb_nnz_idx_t *PA = NULL;
	void *WA = NULL;
	rsb_coo_idx_t *iWA = NULL;
	rsb_coo_idx_t *jWA = NULL;
	rsb_coo_idx_t *nWA = NULL;
	void *vWA = NULL;
	rsb_nnz_idx_t n = 0;
	size_t el_size = RSB_SIZEOF(typecode);
	//struct rsb_mtx_partitioning_info_t pinfop;
	/* const long wet = rsb_get_num_threads();*/ /* want executing threads */
	const long wet = rsb__set_num_threads(RSB_THREADS_GET_MAX_SYS); /* want executing threads; FIXME: it seems there is a severe bug with the definition of RSB_NTC */
#if RSB_DO_WANT_PSORT_TIMING
	rsb_time_t dt,pt,st,mt,ct;
	rsb_time_t tt = - rsb_time();
#endif /* RSB_DO_WANT_PSORT_TIMING */
	const rsb_int_t ei = RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE) ? 1 : 0;

	/* Follows a workaround for possible overflows encountered when compiling with gcc versions through 11.2.0 and running
	$ ./rsbench  --limits-testing
	A stand-alone reproducer for the bug:

// # gcc -Wall -Wextra -pedantic -fopenmp -std=c11 -O0 overflow.c
#include <stdlib.h> // abort
#include <limits.h> // INT_MAX
#include <omp.h> // compile with -fopenmp
int main ()
{
  const int chunk_size = 1000;
  const int n = INT_MAX - 100; // 2147483547
  int l = 0;
  #pragma omp parallel for schedule (static,chunk_size) num_threads (1)
  for (int i = 0; i < n; ++i)
    {
      l = i;
      if(i < 0)
        abort ();
    }
  if ( l != n-1 )
    abort ();
  return 0;
}

	The workaround consists of maximizing chunk size, and seems not necessary with clang-11.0.1 or icc-19.0.5.281.
	More information about the bug on the gcc bug tracker [https://gcc.gnu.org/bugzilla/show_bug.cgi?id=102461] and
	on the GNU Octave tracker [https://savannah.gnu.org/bugs/?60042].
	 */
	if(RSB_INVALID_COO_COUNT(psc+nr))
		psc = RSB_MAX(1,( RSB_MAX_VALUE_FOR_TYPE(rsb_nnz_idx_t) - nr ) / wet);

	if(wet==0)
	{
		/* see RSB_CHECK_LIB_INIT_CHECK */
		errval = RSB_ERR_GENERIC_ERROR;
		RSB_PERR_GOTO(ret,"%s\n","No threads detected! Are you sure to have initialized the library?\n"); /* RSB_ERRM_FLI */
	}

	if(RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		RSB_PERR_GOTO(ret,"\n");
	}

	if(nnz<2)
		goto ret;

	if( nr > RSB_MAX_MATRIX_DIM || nc > RSB_MAX_MATRIX_DIM || nnz > RSB_MAX_MATRIX_NNZ )
	{
		errval = RSB_ERR_LIMITS;
		RSB_PERR_GOTO(ret,RSB_ERRM_ES);
	}

	if(RSB_MUL_OVERFLOW(sizeof(rsb_nnz_idx_t),(nr+2),size_t,rsb_non_overflowing_t)
	|| RSB_MUL_OVERFLOW(sizeof(rsb_nnz_idx_t)*2+el_size,(nnz),size_t,rsb_non_overflowing_t))
	{
		errval = RSB_ERR_LIMITS;
		RSB_PERR_GOTO(err,"sorry, allocating that much memory would cause overflows\n");
	}
	PA = rsb__calloc(sizeof(rsb_nnz_idx_t)*(nr+2));
//	WA = rsb__calloc(RSB_MAX(sizeof(rsb_coo_idx_t),el_size)*(nnz+1));
//	WA = rsb__calloc((2+3*sizeof(rsb_coo_idx_t)+el_size)*nnz);
//	WA = rsb__calloc((2*sizeof(rsb_coo_idx_t)+el_size)*nnz);
	WA = rsb__calloc_parallel((sizeof(rsb_coo_idx_t)*2+el_size)*nnz); // NEW 20101201

	if(!PA || !WA)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,"after calloc, pa=%p, wa=%p\n",PA,WA);
	}
	iWA=((rsb_coo_idx_t*) WA);
	jWA=((rsb_coo_idx_t*)iWA)+nnz;
	vWA=((rsb_coo_idx_t*)jWA)+nnz;
//	nWA=((rsb_char_t*)vWA)+(el_size*nnz);

	/* saving one head element with a trick */
	++PA;
#if RSB_DO_WANT_PSORT_TIMING
	dt = rsb_time();
	mt=-dt;
#endif /* RSB_DO_WANT_PSORT_TIMING */

#if RSB_DO_WANT_PSORT_FASTER_BUT_RISKY
# if 1
	/* FIXME: unfinished and incorrect code */
	#pragma omp parallel for schedule(static,psc) shared(PA,IA) RSB_NTC 
	for(n=0;n<nnz;++n)
		PA[IA[n]+1-ei]++;
# else
	/* actually, this code is VERY SLOW :) */
	#pragma omp parallel reduction(|:errval) shared(PA,IA) 
	{
		rsb_nnz_idx_t n;
		rsb_thr_t th_id = omp_get_thread_num();
		rsb_thr_t tn = omp_get_num_threads();

		for(n=0;RSB_LIKELY(n<nnz);++n)
			if(IA[n]%tn==th_id)
				PA[IA[n]+1-ei]++;
	}
	#pragma omp barrier
#endif
#else /* RSB_DO_WANT_PSORT_FASTER_BUT_RISKY */
	/* setting PA[i] to contain the count of elements on row i */
	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		RSB_ASSERT(IA[n]>=0);
		RSB_ASSERT(IA[n]<=nr);
		PA[IA[n]+1-ei]++;
#if RSB_DO_WANT_PSORT_VERBOSE>1
		RSB_INFO("PA[nr] = %d\n",PA[nr]);
		RSB_INFO("IA[%d] = %d   PA[%d] = %d\n",n,IA[n],IA[n]+1-ei,PA[IA[n]+1-ei]);
#endif /* RSB_DO_WANT_PSORT_VERBOSE */
	}
#endif /* RSB_DO_WANT_PSORT_FASTER_BUT_RISKY */
	/* setting PA[i] to contain the count of elements before row i */
	for(n=0;RSB_LIKELY(n<nr);++n)
	{
		PA[n+1] += PA[n];
#if RSB_DO_WANT_PSORT_VERBOSE>1
		RSB_INFO("PA[%d] = %d\n",n,PA[n]);
#endif /* RSB_DO_WANT_PSORT_VERBOSE */
	}
#if RSB_DO_WANT_PSORT_VERBOSE>1
	RSB_INFO("PA[%d] = %d\n",n,PA[n]);
#endif /* RSB_DO_WANT_PSORT_VERBOSE */

#if RSB_DO_WANT_PSORT_TIMING
	dt = rsb_time();
	mt+=dt;
	pt=-dt;
#endif /* RSB_DO_WANT_PSORT_TIMING */
	/* shuffling elements on the basis of their row 
	 * FIXME : this is the slowest part of this code
	 * its performance is largely dependent on cache lenghts and latencies.. */
	/* FIXME : parallelization of this is challenging */
	rsb__util_do_scatter_rows(vWA,iWA,jWA,VA,IA,JA,PA-ei,nnz,typecode);
	--PA; /* PA has been modified. */
#if RSB_DO_WANT_PSORT_TIMING
	dt = rsb_time();
	pt+=dt;
	st=-dt;
#endif /* RSB_DO_WANT_PSORT_TIMING */
//	RSB_COA_MEMCPY(IA,iWA,0,0,nnz);
//	RSB_COA_MEMCPY(JA,jWA,0,0,nnz);
//	RSB_A_MEMCPY(VA,vWA,0,0,nnz,el_size);
//	SB_COA_MEMCPY_parallel(IA,iWA,0,0,nnz);
//	RSB_COA_MEMCPY_parallel(JA,jWA,0,0,nnz);
//	RSB_A_MEMCPY_parallel(VA,vWA,0,0,nnz,el_size);
	/* restore the row pointers with a trick */

	RSB_ASSERT(PA[nr+1-ei]==nnz);

	/* TODO: parallelization of this ? FIXME: is this necessary ? */
	for(n=0;n<nr;++n)
		mnzpr = RSB_MAX(mnzpr,PA[n+1]-PA[n]);

	nWA = rsb__malloc(sizeof(rsb_nnz_idx_t)*(mnzpr+2)*wet);/* rsb__malloc is evil inside openmp */
	if(!nWA) { RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_ENOMEM);RSB_PERR_GOTO(err,RSB_ERRM_ES); }

	psc = RSB_MIN(psc,nr);
	/* the rows are ready to be sorted (FIXME: this is slow, and could be optimized very much) */
//	#pragma omp parallel for reduction(|:errval)
//	#pragma omp parallel for
//	#pragma omp parallel for schedule(static,10)
//	#pragma omp parallel for schedule(static,1)
//	#pragma omp parallel for schedule(static,psc) shared(iWA,jWA,vWA,nWA,PA)   num_threads(wet)
	#pragma omp parallel for schedule(static,psc) shared(iWA,jWA,vWA,nWA,PA)   RSB_NTC 
	for(n=0;n<nr;++n)
	{
		rsb_nnz_idx_t nnz1,nnz0;
#if 1
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		rsb_thread_t th_id = omp_get_thread_num();
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		rsb_thread_t th_id=0;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		rsb_nnz_idx_t tnoff=th_id*(mnzpr+2);
		nnz1=PA[n+1];
		nnz0=PA[n];
#if RSB_DO_WANT_PSORT_VERBOSE
		RSB_INFO("psort row %d/%d: nonzeros [%d .. %d/%d] on thread %d\n",(int)n,nr,(int)nnz0,(int)nnz1,(int)nnz,(int)th_id);
#endif /* RSB_DO_WANT_PSORT_VERBOSE */
		if(nnz1-nnz0<2)
			continue;/* skip empty line. TODO: could implement with macro sorting algorithms for few nnz */
		if(!RSB_SOME_ERROR(rsb__do_msort_up(nnz1-nnz0,jWA+nnz0,nWA+tnoff)))
			rsb__ip_reord(nnz1-nnz0,RSB_VA_OFFSET_POINTER(vWA,el_size,nnz0),iWA+nnz0,jWA+nnz0,nWA+tnoff,typecode);
#else
		nnz1=PA[n+1];
		nnz0=PA[n];
		rsb__do_util_sortcoo(vWA+nnz0,iWA+nnz0,jWA+nnz0,nr,nc,nnz1-nnz0,typecode,NULL,flags,NULL,0);
#endif
	}

#if RSB_DO_WANT_PSORT_TIMING
	dt = rsb_time();
	st+=dt;
	ct=-dt;
#endif /* RSB_DO_WANT_PSORT_TIMING */

	RSB_COA_MEMCPY_parallel(IA,iWA,0,0,nnz);
	RSB_COA_MEMCPY_parallel(JA,jWA,0,0,nnz);
	RSB_A_MEMCPY_parallel(VA,vWA,0,0,nnz,el_size);
#if RSB_DO_WANT_PSORT_TIMING
	dt = rsb_time();
	ct+=dt;
#endif /* RSB_DO_WANT_PSORT_TIMING */
err:
	RSB_CONDITIONAL_FREE(PA);
	RSB_CONDITIONAL_FREE(WA);
	RSB_CONDITIONAL_FREE(nWA);
#if RSB_DO_WANT_PSORT_TIMING
	dt = rsb_time();
	tt+=dt;
	RSB_INFO("tt:%lg  pt:%lg  st:%lg  mt:%lg ct:%lg\n",tt,pt,st,mt,ct);
#endif /* RSB_DO_WANT_PSORT_TIMING */
ret:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__util_sort_row_major_inner(void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const  rsb_type_t typecode , const rsb_flags_t flags /*, void * WA, size_t wb */)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( nr == 0 || nc == 0 )
		goto ret;

#if RSB_WANT_SORT_PARALLEL_BUT_SLOW
	if( rsb_global_session_handle.asm_sort_method > 0 && RSB_IS_ENOUGH_NNZ_FOR_PARALLEL_SORT(nnz) )
		/* parallel and scaling but slow */
		errval = rsb__util_sort_row_major_parallel(VA,IA,JA,nnz,nr,nc,typecode,flags);
	else
#else /* RSB_WANT_SORT_PARALLEL_BUT_SLOW */
#endif  /* RSB_WANT_SORT_PARALLEL_BUT_SLOW */
		/* not so parallel nor scaling but fast */
		errval = rsb__util_sort_row_major_bucket_based_parallel(VA,IA,JA,nnz,nr,nc,typecode,flags);

ret:
	return errval;
}
/* @endcond */
