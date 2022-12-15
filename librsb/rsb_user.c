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
 * @brief Functions dumping system information to users.
 * */

#include <unistd.h>	/* sysconf */
#include "rsb_internals.h"
#include "rsb.h"
#ifdef RSB_HAVE_LIMITS_H
#include <limits.h>	/* CHAR_BIT */
#endif /* RSB_HAVE_LIMITS_H */
#include <assert.h>	/* assert */
#ifdef RSB_HAVE_MALLOC_H
#include <malloc.h>	/* posix_memalign */
#endif /* RSB_HAVE_MALLOC_H */

#ifdef RSB_HAVE_TIMES_H
#include <sys/times.h>
#endif /* RSB_HAVE_TIMES_H */
#ifdef RSB_HAVE_SYS_SYSTEMCFG_H 
#include <sys/systemcfg.h>	/* for _H_SYSTEMCFG */
#endif /* RSB_HAVE_SYS_SYSTEMCFG_H */
#ifdef RSB_HAVE_SCHED_H
#include <sched.h>	/* sched_getaffinity; FIXME: move to sys.c */
#include "rsb-config.h"
#endif /* RSB_HAVE_SCHED_H */

RSB_INTERNALS_COMMON_HEAD_DECLS

#ifdef _H_SYSTEMCFG
#if 0
from systemcfg.h :
extern struct {
        int architecture;       /* processor architecture */
        int implementation;     /* processor implementation */
        int version;            /* processor version */
        int width;              /* width (32 || 64) */
        int ncpus;              /* 1 = UP, n = n-way MP */
        int cache_attrib;       /* L1 cache attributes (bit flags)      */
                                /* bit          0/1 meaning             */
                                /* -------------------------------------*/
                                /* 31    no cache / cache present       */
                                /* 30    separate I and D / combined    */
        int icache_size;        /* size of L1 instruction cache */
        int dcache_size;        /* size of L1 data cache */
        int icache_asc;         /* L1 instruction cache associativity */
        int dcache_asc;         /* L1 data cache associativity */
        int icache_block;       /* L1 instruction cache block size */
        int dcache_block;       /* L1 data cache block size */
        int icache_line;        /* L1 instruction cache line size */
        int dcache_line;        /* L1 data cache line size */
        int L2_cache_size;      /* size of L2 cache, 0 = No L2 cache */
        int L2_cache_asc;       /* L2 cache associativity */
        int tlb_attrib;         /* TLB attributes (bit flags)           */
                                /* bit          0/1 meaning             */
                                /* -------------------------------------*/
                                /* 31    no TLB / TLB present           */
                                /* 30    separate I and D / combined    */
        int itlb_size;          /* entries in instruction TLB */
        int dtlb_size;          /* entries in data TLB */
        int itlb_asc;           /* instruction tlb associativity */
        int dtlb_asc;           /* data tlb associativity */
        long long physmem;      /* bytes of OS available memory             */
..
}_system_configuration;
#endif

static rsb_err_t aix_sys_info()
{
	/*!
	 	\ingroup internals
	*/
	RSB_INFO("Working on an AIX system\n");
	RSB_INFO("CPU		: %ld \n",_system_configuration.ncpus);
	RSB_INFO("cache_at	:%ld \n",_system_configuration.cache_attrib);
	RSB_INFO("L1		: %ld \n",_system_configuration.dcache_size);
	RSB_INFO("L2		: %ld \n",_system_configuration.L2_cache_size);
	RSB_INFO("MEM		: %lld \n",_system_configuration.physmem);
}
#endif /* _H_SYSTEMCFG */



static rsb_err_t get_sysconf_cacheinfo( long *cpa, long *cpb, long *cpc, int cac,  int cbc,  int ccc, int cl)
{
	/*!
	 \ingroup internals
	*/
	*cpa = sysconf(cac);
	*cpb = sysconf(cbc);
	*cpc = sysconf(ccc);
	if(*cpa<1 || *cpb < 1 || *cpc < 1)
		RSB_INFO("sysconf() : no level %d cache\n",cl);
	else
	{
		RSB_INFO("sysconf() : level %d cache size %ld \n",cl,*cpc);
		RSB_INFO("sysconf() : level %d cache associativity %ld \n",cl,*cpa);
		RSB_INFO("sysconf() : level %d cache line size %ld \n",cl,*cpb);
	}
	return RSB_ERR_NO_ERROR;
}

static long rsb_max_threads(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Just a user-oriented function.
	 *
	 * \return the maximum number of available hardware threads
	 *
	 * If on AIX, we use the native solution, as sysconf() gives values with are not usable as threads.
	 * */
#ifdef _H_SYSTEMCFG
	return _system_configuration.ncpus;
#else /* _H_SYSTEMCFG */
#ifdef RSB_HAVE_SYSCONF 
	/*
	 * _SC_NPROCESSORS_ONLN : The number of processors currently online (available).
	 * _SC_NPROCESSORS_CONF : The number of processors configured.
	 */
	//return sysconf(_SC_NPROCESSORS_CONF);
	return sysconf(_SC_NPROCESSORS_ONLN);
#else /* RSB_HAVE_SYSCONF  */
	return 0;	/* this should be regarded as an error */
#endif /* RSB_HAVE_SYSCONF  */
#endif /* _H_SYSTEMCFG */
}

rsb_err_t rsb__sys_info()
{
	/*!
	 \ingroup internals
	 *
	 * A function printing out information about the system.
	 * It gives information for the user about the library configuration.
	 * It should be called after library initialization.
	 *
	 * \return an error code or RSB_ERR_NO_ERROR.
	 * TODO: move to sys.c
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if RSB_WITH_HWLOC
	{
		int i;
		for(i=1;i<4;++i)
		{
			size_t sz = rsb__get_lnc_size_hwloc(i);
			if(sz)
       				RSB_INFO("hwloc size of cache level %d: %zd\n",i,sz);
		}
	}
#endif	/* RSB_WITH_HWLOC */

       	RSB_INFO("detected max available cores/threads : %ld\n",(long int)rsb_max_threads());
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	#pragma omp parallel
	{
       	RSB_INFO("detected max OpenMP procs : %ld\n",(long int)omp_get_num_procs());
	}
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
       	RSB_INFO("detected %ld levels of cache\n",(long int)rsb__get_cache_levels_num());
	{
		int i;
		for(i=1;i <= rsb__get_cache_levels_num();++i)
		       	RSB_INFO("L%d size: %ld \n",i,(long int)rsb__get_lnc_size(i));
	}

#ifdef _H_SYSTEMCFG
	aix_sys_info();
#endif /* _H_SYSTEMCFG */

//       	RSB_INFO("LL size: %ld \n",(long int)rsb__get_lastlevel_c_size());

#ifndef RSB_HAVE_SYSCONF
        RSB_INFO("sysconf() not available\n");
#else /* RSB_HAVE_SYSCONF */
//       	RSB_INFO("detected %ld levels of cache\n",(long int)rsb__get_cache_levels_num());
#endif /* RSB_HAVE_SYSCONF */
	{
#ifdef RSB_HAVE_SYSCONF
        long int pagesize = 0;
        long int mem_pages = 0;
        size_t tot_mem = 0;
#if   defined(PAGESIZE)
        pagesize = sysconf(PAGESIZE);
#elif defined(_SC_PAGESIZE)
        pagesize = sysconf(_SC_PAGESIZE);
#elif defined(PAGE_SIZE)
        pagesize = sysconf(PAGE_SIZE);
#else /* PAGE_SIZE */
#endif /* PAGE_SIZE */
        if( pagesize)RSB_INFO("sysconf() : %ld bytes per pagesize\n",pagesize);
        if(!pagesize)RSB_INFO("sysconf() available, PAGESIZE _SC_PAGESIZE PAGE_SIZE undefined\n");

	/* 
	 _SC_AVPHYS_PAGES : The number of currently available pages of physical memory.
	 _SC_PHYS_PAGES   : The  number  of pages of physical memory.
	*/
#if   defined(_SC_PHYS_PAGES)
        mem_pages = sysconf(_SC_PHYS_PAGES);
#else /* _SC_PHYS_PAGES */
#endif /* _SC_PHYS_PAGES */
        tot_mem = (size_t)mem_pages;
	tot_mem *= (size_t)pagesize;
        if( mem_pages)RSB_INFO("sysconf() : %zu physical pages\n",(size_t)mem_pages);
        if(!mem_pages)RSB_INFO("sysconf() available, _SC_PHYS_PAGES undefined\n");
        if( mem_pages && pagesize)RSB_INFO("sysconf() : %zu bytes (%zu MB) of physical memory\n",tot_mem,(tot_mem)/(1024*1024));
#if   defined(_SC_AVPHYS_PAGES)
        RSB_INFO("sysconf() : %zu available (free) physical pages\n",(size_t)sysconf(_SC_AVPHYS_PAGES));
        RSB_INFO("sysconf() : %zu available (free) physical memory\n",(size_t)(sysconf(_SC_AVPHYS_PAGES)*pagesize));
#endif /* _SC_AVPHYS_PAGES */
#endif /* RSB_HAVE_SYSCONF */
	}
	{
#ifdef RSB_HAVE_SYSCONF 
	long int sc_nprocessors_conf;
	long int sc_nprocessors_onln;
	/*
	 * _SC_NPROCESSORS_ONLN : The number of processors currently online (available).
	 * _SC_NPROCESSORS_CONF : The number of processors configured.
	 */
	sc_nprocessors_conf = sysconf(_SC_NPROCESSORS_CONF);
	sc_nprocessors_onln = sysconf(_SC_NPROCESSORS_ONLN);
	RSB_INFO("sysconf() , processors : %ld\n",sc_nprocessors_conf);
	RSB_INFO("sysconf() , processors online : %ld\n",sc_nprocessors_onln);
#endif /* RSB_HAVE_SYSCONF  */
	}
#ifdef RSB_HAVE_SYSCONF 
	{
#ifdef _SC_LEVEL1_DCACHE_SIZE 
	long int c1a,c1b,c1c;
	c1a = sysconf(_SC_LEVEL1_DCACHE_ASSOC);
	c1b = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
	c1c = sysconf(_SC_LEVEL1_DCACHE_SIZE);
	get_sysconf_cacheinfo(&c1a,&c1b,&c1c,_SC_LEVEL1_DCACHE_ASSOC,_SC_LEVEL1_DCACHE_LINESIZE,_SC_LEVEL1_DCACHE_SIZE,1);
#else /* _SC_LEVEL1_DCACHE_SIZE */
	RSB_INFO("sysconf() implementation obsolete: no L%d cache info\n",1);
#endif /* _SC_LEVEL1_DCACHE_SIZE */
	}
	{
#ifdef _SC_LEVEL2_CACHE_SIZE 
	long int c2a,c2b,c2c;
	c2a = sysconf(_SC_LEVEL2_CACHE_ASSOC);
	c2b = sysconf(_SC_LEVEL2_CACHE_LINESIZE);
	c2c = sysconf(_SC_LEVEL2_CACHE_SIZE);
	get_sysconf_cacheinfo(&c2a,&c2b,&c2c,_SC_LEVEL2_CACHE_ASSOC,_SC_LEVEL2_CACHE_LINESIZE,_SC_LEVEL2_CACHE_SIZE,2);
#else /* _SC_LEVEL2_CACHE_SIZE */
	RSB_INFO("sysconf() implementation obsolete: no L%d cache info\n",2);
#endif /* _SC_LEVEL2_CACHE_SIZE */
	}
	{
#ifdef _SC_LEVEL3_CACHE_SIZE 
	long int c3a,c3b,c3c;
	c3a = sysconf(_SC_LEVEL3_CACHE_ASSOC);
	c3b = sysconf(_SC_LEVEL3_CACHE_LINESIZE);
	c3c = sysconf(_SC_LEVEL3_CACHE_SIZE);
	get_sysconf_cacheinfo(&c3a,&c3b,&c3c,_SC_LEVEL3_CACHE_ASSOC,_SC_LEVEL3_CACHE_LINESIZE,_SC_LEVEL3_CACHE_SIZE,3);
#else /* _SC_LEVEL3_CACHE_SIZE  */
	RSB_INFO("sysconf() implementation obsolete: no L%d cache info\n",3);
#endif /* _SC_LEVEL3_CACHE_SIZE  */
	}
	{
#ifdef _SC_LEVEL4_CACHE_SIZE 
	long int c4a,c4b,c4c;
	c4a = sysconf(_SC_LEVEL4_CACHE_ASSOC);
	c4b = sysconf(_SC_LEVEL4_CACHE_LINESIZE);
	c4c = sysconf(_SC_LEVEL4_CACHE_SIZE);
	get_sysconf_cacheinfo(&c4a,&c4b,&c4c,_SC_LEVEL4_CACHE_ASSOC,_SC_LEVEL4_CACHE_LINESIZE,_SC_LEVEL4_CACHE_SIZE,4);
#else /* _SC_LEVEL4_CACHE_SIZE */
	RSB_INFO("sysconf() implementation obsolete: no L%d cache info\n",4);
#endif /* _SC_LEVEL4_CACHE_SIZE */
	}
#endif /* RSB_HAVE_SYSCONF */
#ifdef CHAR_BIT
	/* It should happen, but it could not. */
	RSB_ASSERT(CHAR_BIT==sizeof(char)*8);

	/* It should not happen, but it could. */
	if(CHAR_BIT!=8)
	{
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
		RSB_INFO("%d bits per byte! This is catastrophic.\n",CHAR_BIT);
	}
	else
		RSB_INFO("8 bits per byte. Good.\n");
#else /* CHAR_BIT */
	RSB_INFO("We have no information on bits per byte. Beware!\n");
#endif /* CHAR_BIT */
#if 1
	{
		long cbbs = rsb__get_cache_block_byte_size();
		if(cbbs)
			RSB_STDERR("cache block size		: %ld \n",cbbs);
		else
			RSB_STDERR("cache block size unknown (detected %ld: this is a problem!)\n",cbbs);
	}
#endif
#ifdef INT_MAX
	RSB_INFO("SHRT_MAX : %hd\n",(short)SHRT_MAX);
	RSB_INFO("SHRT_MIN : %hd\n",(short)SHRT_MIN);
	RSB_INFO("USHRT_MAX : %hu\n",(unsigned short)USHRT_MAX);
	RSB_INFO("INT_MIN : %d\n",(int)INT_MIN);
	RSB_INFO("INT_MAX : %d\n",(int)INT_MAX);
	RSB_INFO("UINT_MAX : %u\n",(unsigned)UINT_MAX);
	RSB_INFO("LONG_MAX : %ld\n",(long int)LONG_MAX);
	RSB_INFO("LONG_MIN : %ld\n",(long int)LONG_MIN);
#ifdef ULONG_MAX
	RSB_INFO("ULONG_MAX : %lu\n",(long unsigned)ULONG_MAX);
#else /* ULONG_MAX */
	RSB_INFO("ULONG_MAX : undefined\n");
#endif /* ULONG_MAX */
#ifdef 	LLONG_MAX 
	RSB_INFO("LLONG_MAX : %lld\n",(long long int)LLONG_MAX);
#else /* LLONG_MAX  */
	RSB_INFO("LLONG_MAX : undefined\n");
#endif /* LLONG_MAX  */
#ifdef LLONG_MIN
	RSB_INFO("LLONG_MIN : %lld\n",(long long int)LLONG_MIN);
#else /* LLONG_MIN */
	RSB_INFO("LLONG_MIN : undefined\n");
#endif /* LLONG_MIN */
#ifdef ULLONG_MAX
	RSB_INFO("ULLONG_MAX : %llu\n",(long long unsigned)ULLONG_MAX);
#else /* ULLONG_MAX */
	RSB_INFO("ULLONG_MAX : undefined\n");
#endif /* ULLONG_MAX */
#else /* INT_MAX */
	RSB_INFO("INT_MAX : undefined\n");
#endif /* INT_MAX */
	RSB_INFO("RSB_MARKER_COO_VALUE : %llu\n",(long long unsigned)RSB_MARKER_COO_VALUE);
	RSB_INFO("RSB_MARKER_NNZ_VALUE : %llu\n",(long long unsigned)RSB_MARKER_NNZ_VALUE);
	RSB_INFO("RSB_SUBM_IDX_MARKER : %llu\n",(long long unsigned)RSB_SUBM_IDX_MARKER);
	RSB_INFO("RSB_MAX_ALLOCATABLE_MEMORY_CHUNK: %llu\n",(long long unsigned)RSB_MAX_ALLOCATABLE_MEMORY_CHUNK);

	RSB_INFO("timing min delta (if negative, don't complain with us)   : %lg s\n", rsb__timer_sanity());
	RSB_INFO("timing granularity : %lg s\n", RSB_CACHED_TIMER_GRANULARITY);
#if   defined(RSB_CFLAGS)
	RSB_INFO("CFLAGS   : %s\n",RSB_CFLAGS);
#else /* RSB_CFLAGS */
	RSB_INFO("no CFLAGS info\n");
#endif /* RSB_CFLAGS */
#if   defined(RSB_CXXFLAGS)
	RSB_INFO("CXXFLAGS : %s\n",RSB_CXXFLAGS);
#else /* RSB_CXXFLAGS */
	RSB_INFO("no CXXFLAGS info\n");
#endif /* RSB_CXXFLAGS */
#if   defined(RSB_CC)
	RSB_INFO("CC       : %s\n",RSB_CC);
#else /* RSB_CC */
	RSB_INFO("no CC info\n");
#endif /* RSB_CC */
#ifdef RSB_HAVE_SCHED_H
#ifdef RSB_HAVE_SCHED_GETAFFINITY 
#ifdef _GNU_SOURCE
	{
		size_t num_cpus = CPU_SETSIZE;
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		if(1)
		{
			int sgar = 0;

			if( (sgar = sched_getaffinity(0, num_cpus, &cpuset)) != 0 )
			{
				RSB_INFO("sched_getaffinity error : %d\n",sgar);
			}
			else
			{
				RSB_INFO("sched_getaffinity's CPU_COUNT() of set:    %d\n", CPU_COUNT_S(CPU_ALLOC_SIZE(1), &cpuset));
				RSB_INFO("sched_getaffinity runnable : %zd\n",CPU_COUNT(&cpuset));
			}
		}
	}
#endif /* _GNU_SOURCE */
#endif /* RSB_HAVE_SCHED_H */
#endif /* RSB_HAVE_SCHED_GETAFFINITY */
	{
	rsb_char_t usmhib[RSB_MAX_LINE_LENGTH];
	RSB_INFO("memhinfo : %s\n",rsb__get_mem_hierarchy_info_string(usmhib));
	}
        RSB_INFO("detected free  memory : %zd\n",(size_t)rsb__sys_free_system_memory());
        RSB_INFO("detected total memory : %zd\n",(size_t)rsb__sys_total_system_memory());

	{
		rsb_nnz_idx_t *p = NULL;
		rsb_nnz_idx_t n, i, maxtries = RSB_MAX_MATRIX_NNZ, res = 0, cookie = 0, tries, v, mintries = RSB_CONST_MIN_TIMES_FOR_MICRO_BENCHMARK;
		n = 4*(rsb__get_lastlevel_c_size()/sizeof(rsb_nnz_idx_t));
		if(n<2)goto failed;
		p = rsb__malloc(sizeof(rsb_nnz_idx_t)*n);
		if(!p)goto failed;
		v = 1;
		while(2*v <= n)v *= 2;
		--v;
		for(i=0;i<n;++i)p[i] = i;
		while((v/2)>1)
		{
			rsb_time_t mtl = RSB_CONST_IMPOSSIBLY_BIG_TIME,mtb = RSB_CONST_IMPOSSIBLY_BIG_TIME,bt,mbt = RSB_CONST_TIME_FOR_MICRO_BENCHMARK,tt=0;
			for(tt=0,tries=0;tries<mintries || (tt<mbt && tries<maxtries);++tries)
			{
				/* NOTE: we are not interested in flushing the cache, here */
				bt = - rsb_time();
				cookie += rsb__seek_nnz_idx_t(p,v,n);
				bt += rsb_time();
				mtb = RSB_MIN(mtb,bt);
				tt += bt;
				bt = - rsb_time();
				cookie += rsb__seek_nnz_idx_t_linear(p,v,n);
				bt += rsb_time();
				mtl = RSB_MIN(mtl,bt);
				tt += bt;
			}
			res = cookie;
			RSB_INFO("for array sized %ld elems, took %g s for linear search and %g s for binary search for element %ld, in %ld tries, for a total of %f s (ignore this:%ld)\n", (long int)n,mtl,mtb,(long int)v,(long int)tries,tt,(long int)res);
			v = v/2;
		}
failed:
		RSB_CONDITIONAL_FREE(p);
	}	errval = rsb__dump_system_performance_summary();
	
	goto err;
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_real_t rsb__getenv_real_t(const char*envv, const rsb_real_t altv)
{
	const char * v = rsb__getenv(envv);
	return v ? rsb__util_atof(v) : altv;
}

rsb_int_t rsb__getenv_int_t(const char*envv, const rsb_int_t altv)
{
	const char * v = rsb__getenv(envv);
	return v ? rsb__util_atoi(v) : altv;
}

const rsb_char_t * rsb__getenv_str(const char*envv, const rsb_char_t* altv)
{
	const char * v = rsb__getenv(envv);
	return v ? v : altv;
}

rsb_char_t rsb__getenv_char(const char *envv, const rsb_char_t altv)
{
	const char * v = rsb__getenv(envv);
	return v ? *v : altv;
}

/* @endcond */
