/*

Copyright (C) 2008-2021 Michele Martone

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
/**
 * @file
 * @author Michele Martone
 * @brief Memory bandwidth related (e.g.: read, write, and read-write) microbenchmarks for cache based machines.
 */

#include "rsb_common.h"
#include "rsb.h"
#include <strings.h>	/* RSB_BZERO */
#include <stdint.h> /* uintptr_t */

#define RSB__PC_STRIDE 1 /* stride for pointer chasing */

/* require sizeof(w_t) >= sizeof(void*) */
typedef uintptr_t w_t;

size_t rsb__entropy;		/**< private checksum-only variable, necessary to avoid compiler optimization of memory scan operations */

RSB_INTERNALS_COMMON_HEAD_DECLS

static const char * rsb__mbw_s2s(rsb_flags_t btype)
{
	/**
		\ingroup gr_internals
	 	\return a pointer to a const string descriptive of this particular memory measurement
	 */
	switch(btype)
	{
		case RSB_MB_INVALID:	
			return "INVALID";
			break;
		case RSB_MB_READ	:	
			return "READ";
			break;
		case RSB_MB_WRITE	:	
			return "WRITE";
			break;
		case RSB_MB_RW		:
			return "RW";
			break;
		case RSB_MB_BZERO		:
			return "BZERO";
			break;
		case RSB_MB_ZERO		:
			return "ZERO";
			break;
		case RSB_MB_MEMCPY		:
			return "MEMCPY";
			break;
		case RSB_MB_MEMCPY2		:
			return "MEMCPY2";
			break;
		case RSB_MB_MEMSET		:
			return "MEMSET";
			break;
		case RSB_MB_LINEAR_CHASE	:
			return "LINEAR_CHASE";
			break;
		case RSB_MB_MORTON_CHASE	:
			return "MORTON_CHASE";
			break;
		default:
		/* error */
		return "";
	}
}

typedef int cb_t;
typedef size_t zb_t ;

static void h2c(cb_t *bip, cb_t *bjp, zb_t bz)
{
	/**
	 * \ingroup gr_internals
	 * morton to coordinate
	 */
        int b;

        *bip=0;
        *bjp=0;
//      RSB_STDERR("-> %ld\n",bz);
	/* this would greatly benefit of bit interleaving (absent on x86) */
        for(b=0;b<8*sizeof(cb_t);++b)
        {
                /* mancano queste due righe */
                *bip|=(bz&(0x1<<(2*b+1)))>>(b+1);
                *bjp|=(bz&(0x1<<(2*b+0)))>>(b+0);
        }
}

#if 0
/* unused function */
static void c2h(cb_t bi, cb_t bj, zb_t * bzp)
{
	/** 
	 * \ingroup gr_internals
	 * coordinate to morton
	 */
        int b;
        *bzp=0;
//      RSB_STDERR("b : %d %d\t",bi,bj);
        for(b=0;b<8*sizeof(cb_t);++b)
        {
                *bzp|=(bi&(0x1<<b))<<(b+1);
                *bzp|=(bj&(0x1<<b))<<(b+0);
        }
//      RSB_STDERR("z : %9ld\n",*bzp);
}
#endif

static int morton_pointer_chasing_pattern_init(w_t *p, size_t dim)
{
	/**
	 * \ingroup gr_internals
	 *
	 * TODO: Document me.
	 * Note: requires dim==4^k for some k. If not, then only half benchmark is performed.
	 * May fix this with a tiling approach.
         */
	const int stride = RSB__PC_STRIDE;
	int words,i;
	int ni=0,oi=0,e=0 /* dim>=2^e */;
	int side=0;
	{int tmp=dim/(stride*2*sizeof(w_t));while(tmp>0){++e;tmp/=2;}e/=2;e*=2;/* e is even */}

	if(e<1)
		return -1;

	words = (1<<e);
	side = (1<<(e/2));
	oi = 0;
	ni = 0;
/*
	RSB_STDERR("morton_pointer_chasing_pattern_init\n");
	RSB_STDERR("%d\n",e);
	RSB_STDERR("should span %d bytes\n",dim);
	RSB_STDERR("%d side\n",side);
	RSB_STDERR("will span %d words \n",words);
	RSB_STDERR("will span %d bytes \n",words*sizeof(w_t));*/
	for(i=0;i<words;++i)
	{
		//int j=0;
		int nx=0,ny=0;
		/* WARNING : we could not have 8 bits per byte */
		h2c( &nx,  &ny,  i+1);
		nx=nx%side;
		ny=ny%side;
		ni=(nx+ny*side)%words;
		//RSB_STDERR("%d %d\n",nx,ny);
		//RSB_STDERR("%d\n",ni);
		//RSB_STDERR("%d %d %d %d\n",nx,ny,i+1,ni);
		*(w_t**)&p[oi*stride]=(w_t*)&p[ni*stride];
		oi=ni;
	}
	return 0;
}

static int pointer_chasing_pattern_init(w_t *p, size_t dim)
{
	/**
	 * \ingroup gr_internals
         * Initializes a memory area to perform a linear pointer chasing.
         */
	int i;
	const int stride = RSB__PC_STRIDE;
	const int words=dim/(stride * sizeof(w_t));

	for(i=1;i<=words;++i)
		*(w_t**)&p[(i-1)*stride]=(w_t*)&p[(i*stride)%words];

/*        for (i = stride; i < range; i += stride) {
                *(char **)&addr[i - stride] = (char*)&addr[i];
        }
        *(char **)&addr[i - stride] = (char*)&addr[0];*/
	return 0;
}

static int rsb__scan_cache(w_t *p, size_t dim, int should, size_t times, w_t *q)
{
	/**
	 * \ingroup gr_internals
	 * Performs a naive memory scan with side effect.
	 *
	 *  The memory scan will operate on a memory area of dim bytes.
	 *  Will touch consecutively memory locations with a stride of 
	 *  sizeof(w_t) bytes.
	 *  Since we hadn't unrolled the following loops, this should be 
	 *  aggressively unrolled in a way to remain memory bound.
	 *
	 *  Please note that if dim is less than L1/L2/L3 cache, you will
	 *  not effectively benchmark your memory subsystem, but only caches.
	 *
	 *  It is advised for p to be aligned in some way for better performance.
	 *
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * */
	w_t n = 0;
	int i,t;
	const int words = dim / sizeof(w_t);
	
	RSB_DEBUG_ASSERT(times < INT_MAX);
	if (times > INT_MAX) /* avoid t could overflow if times too large */
		goto err;

	if(should == RSB_MB_MEMCPY2 && !q)
		goto err;

	if(should == RSB_MB_FLUSH)/* we ignore times */
	{
		#pragma omp parallel for schedule(static,1) RSB_NTC
		for(i=0;i<words;++i)
			n += p[i];
	}
	else
	/*
	 * Warning : if the compiler is really, really smart, it
	 * could detect we are zeroing.
	 * That could result in an excessively high value here.
	 * */
		/* Note: xlc -O4 was much smarter than us here (maybe on times!). */
		if(should == RSB_MB_ZERO)
			for(t=0;t<times;++t)
			{
				i=0;
				for(;i+7<words;i+=8)
					p[i+0]=0,p[i+1]=0,
					p[i+2]=0,p[i+3]=0,
					p[i+4]=0,p[i+5]=0,
					p[i+6]=0,p[i+7]=0;
				for(;i<words;i++)
					p[i]=0;
			}
	else
		if(should == RSB_MB_BZERO)
			for(t=0;t<times;++t)
				RSB_BZERO(p,dim);
	/*
		WARNING : memcpy operations involve two buffers or two halves ! 
		so be careful when interpreting these results:
		memory bandwidth is double than transfer speed.
	 */
	else
		if(should == RSB_MB_MEMCPY)
			for(t=0;t<times;++t)
				memcpy(p,((char*)p)+dim/2,dim/2);
	else
		if(should == RSB_MB_MEMCPY2)
			for(t=0;t<times;++t)
				memcpy(p,q,dim);
	else
		/* Note: xlc -O4 was much smarter than us here (maybe on times!). */
		if(should == RSB_MB_WRITE)
			for(t=0;t<times;++t)
			{
				i=0;
				for(;i+7<words;i+=8)
					p[i+0]=i+0,p[i+1]=i+1,
					p[i+2]=i+2,p[i+3]=i+3,
					p[i+4]=i+4,p[i+5]=i+5,
					p[i+6]=i+6,p[i+7]=i+7;
				for(;i<words;i++)
					p[i]=i;
			}
	else
		/* Note: xlc -O4 was much smarter than us here (maybe on times!). */
		if(should == RSB_MB_READ)
			for(t=0;t<times;++t)
			{
				// double loop == loop overhead
				i=0;
				for(;i+7<words;i+=8)
					n+=p[i] +p[i+1] +p[i+2] +p[i+3] +p[i+4] +p[i+5] +p[i+6] +p[i+7];
				for(;i<words;i++)
					n+=p[i];
			}
	else
		if(should == RSB_MB_MORTON_CHASE || should == RSB_MB_LINEAR_CHASE)
			for(t=0;t<times;++t)
				for(i=0;i<words;++i)
					p=*(w_t**)p;	/* pointer chasing */
	else
		if(should == RSB_MB_RW)
			for(t=0;t<times;++t)
				for(i=0;i<words;++i)
					p[i]+=i;
	else
		if(should == RSB_MB_MEMSET)
			for(t=0;t<times;++t)
				memset(p,0x0A0B0C0D,dim);
	else
		goto err;

	// return n+*p;	/* WARNING: easily optimizable! should rather accumulate into a 'pool' of entropy! */
	rsb__entropy += n+*p;
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_GENERIC_ERROR;
}

static int rsb_mbw_area_init_and_cache_flush(size_t sz, w_t *fc, w_t *p, int btype /*, int * entropy*/, size_t times)
{
	/**
	 * \ingroup gr_internals
	 * Will init the memory area for benchmarking and then
	 * flush the cache, assuming that its size is sz.
	 */
	switch(btype)
	{
		case(RSB_MB_LINEAR_CHASE):
			return pointer_chasing_pattern_init(p, sz);
		case(RSB_MB_MORTON_CHASE):
			return morton_pointer_chasing_pattern_init(p, sz);
		default:
			return 0;
	}
	rsb__scan_cache(fc,sz,RSB_MB_FLUSH,times, NULL);	/* flush cache */
}

static rsb_time_t mbw_total_time( struct rsb_mbw_m_t *mbw_m  )
{
	/**
	 * \ingroup gr_internals
	 */
	rsb_time_t t = RSB_REAL_ZERO;
	int i;
	if(!mbw_m)
		return t;
	for(i=0;i<RSB_MB_N;++i)
	{
		t+=mbw_m->mb[i].t;
	}
	return t;
}

static rsb_err_t mbw_test( struct rsb_mbw_m_t *mbw_m  );

static rsb_err_t probe_approx_mbw( struct rsb_mbw_m_t * mbwm, rsb_time_t s )
{
	/**
	 * \ingroup gr_internals
	 * we run this quick test and return a rough estimate
	 * of the number of times the test should be performed
	 * on the given memory area to last circa s seconds.
	 * (assumes all memory tests)
	 * */
	const rsb_time_t min_time=0.1;
	size_t times=0;
	rsb_time_t t=0;/* some compilers (e.g.: pgcc) don't init variables for us :) */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( !mbwm )
		{errval = RSB_ERR_BADARGS;goto err;}
	mbwm->times=1;	/* we set times */
	if(s<min_time)
		{errval = RSB_ERR_BADARGS;goto err;}

	while(t<min_time && mbwm->times<=RSB_MAX_TIMES_T)
	{
		mbwm->times*=2;	/* we set times */
		if((errval=mbw_test(mbwm))) /* we perform benchmarking */
			goto err;
		t=mbw_total_time( mbwm  );
		if(t <= RSB_REAL_ZERO)
			{errval = RSB_ERR_INTERNAL_ERROR;goto err;}
	}
	/* times/s == mbwm.times/t */
	times=(int)(((double)mbwm->times)/t)*s;
	if(times<=0 /*overflow ?*/ /* || times < 100*/)
#ifdef INT_MAX 
	{
		times=INT_MAX;
		return 0;
	}
#else /* INT_MAX  */
		{errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#endif /* INT_MAX  */
	/* finally, we set our estimate 'times' value for s seconds benchmarking */
	mbwm->times=times;
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t mbw_test( struct rsb_mbw_m_t *mbw_m  )
{
	/**
	 * \ingroup gr_internals
	 * Run each memory benchmark.
	 * Assume existence of hardware-managed caches.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	w_t * p=NULL,*fc=NULL,*q=NULL;
	struct rsb_mbw_m_t m; 
	int i;

	if(!mbw_m)
		return RSB_ERR_BADARGS;

	m=*mbw_m;

	/* if m.times is zero, probe for an appropriate value */
	if(m.times == 0 && (errval=probe_approx_mbw(&m,1.0)))
	{
		RSB_STDERR("uhm. timing problems ?!.\n");
		rsb__do_perror(NULL,errval);
		goto errl;
	}

	/* 
	 * In order for flushing to work effectively; the flush array must be big enough.
	 * Therefore it is advised to set m.hlcs to at least the size of the bigger cache.
	 * */
	p = rsb__aligned_malloc( m.sz , m.sz );
	q = rsb__aligned_malloc( m.sz , m.sz );/* q is auxiliary */
	fc= rsb__aligned_malloc( m.hlcs , m.hlcs );

	if(!p || !fc || !q)
	{
		if(!p || !q) { RSB_STDERR("problems allocating %zd bytes.\n",(size_t)m.sz  ); }
		if(!fc) { RSB_STDERR("problems allocating %zd bytes.\n",(size_t)m.hlcs); }
		errval = RSB_ERR_GENERIC_ERROR;
		goto errl;
	}

	for(i=0;i<RSB_MB_N;++i)
	{
		m.mb[i].btype=i;	/* we set benchmark type */
		rsb_mbw_area_init_and_cache_flush(m.sz, fc, p, i/*, int * entropy*/, m.times);
		m.mb[i].t = - rsb_time();
		rsb__entropy += rsb__scan_cache(p,m.sz,i,m.times,q);	/* we perform measurement */
		m.mb[i].t += rsb_time();
	}

	// about commenting the following : DANGER
	//if(m.entropy)fprintf(stderr,"the following number is printed only for tricking the compiler optimizer, and has no other use: %d\n",entropy); /* this is essential */
	if(mbw_m)
		*mbw_m=m;
errl:
	RSB_CONDITIONAL_FREE(p);
	RSB_CONDITIONAL_FREE(fc);
	RSB_CONDITIONAL_FREE(q);
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t mbw_ratio_printf(struct rsb_mbw_m_t *h, struct rsb_mbw_m_t *l)
{
	/**
	 * \ingroup gr_internals
	 * Print the ratio in performance of two measurements.
	 */
	const double M = 1.0;
	int i;

	RSB_DEBUG_ASSERT(!(!h||!l));

	for(i=0;i<RSB_MB_N;++i)
		RSB_INFO("#%-32s ratio  %lg \n"  ,
			rsb__mbw_s2s(i),
			((((double)h->times)*h->sz)/(h->mb[i].t*M))/
			((((double)l->times)*l->sz)/(l->mb[i].t*M))
			);
	return RSB_ERR_NO_ERROR;
}

static rsb_err_t mbw_printf(struct rsb_mbw_es_t *esp, struct rsb_mbw_m_t *m, int level)
{
#define RSB_MBW_M_T(MP,I) ((((double)MP->times)*MP->sz)/(MP->mb[i].t*M))
	/**
	 * \ingroup gr_internals
		Prints out memory benchmarks results. 
	 * TODO: if esp is here, it does not print anymore :-)
	*/
	int i;
	const double M=1000000.0;

	if(!m)
		return RSB_ERR_BADARGS;

	if(!esp)
		RSB_INFO("#%-32s\tsize\tlevel\tbw(MBps)\n","size");

	for(i=0;i<RSB_MB_N;++i)
	{
		if(!esp)
			RSB_INFO("%-32s\t%zd\t%zd\t%lg\n",rsb__mbw_s2s(m->mb[i].btype),(rsb_printf_int_t)m->sz,(rsb_printf_int_t)level,RSB_MBW_M_T(m,i));

		if(esp)
			esp[i].bw=RSB_MBW_M_T(m,i),
			esp[i].sz=m->sz,
			esp[i].lvl=level,
			esp[i].mbt=m->mb[i].btype;
	}

	return RSB_ERR_NO_ERROR;
#undef RSB_MBW_M_T
}

rsb_err_t rsb__mem_hier_timings(struct rsb_mbw_cm_t * cm)
{
	/**
	 * \ingroup gr_internals
	 * Measures memory bandwidth in scanning increasingly sized arrays.
	 * These are sized and aligned like initial memory hierarchies (caches) and then more.
	 */
	int cln=0,cl;
	struct rsb_mbw_m_t mbw_m,*mbw_ms=NULL;
	long cs = 0;
	const long extra_levels = RSB_MEMSCAN_EXTRA_LEVELS;

	if( !cm )
		return RSB_ERR_BADARGS;

	cln = rsb__get_cache_levels_num();

	if(cln<1)
	{
		RSB_INFO("No information about caches, sorry\n");
		return -1;
	}

	mbw_ms = rsb__calloc((cln+extra_levels) * sizeof(*mbw_ms));

	if(!mbw_ms)
	{
		goto err;
	}

	RSB_INFO("# This test will measure times in scanning arrays sized and aligned to fit in caches.\n");
	RSB_INFO("# %d cache levels detected\n",cln);

	/* we do measure for each level in the cache hierarchy plus extra_levels */
	for(cl=1;cl<=cln+extra_levels;++cl)
	{
		/* timing for cache level cl  */

		if(cl<=cln)
			cs = rsb__get_lnc_size(cl);
		else
			cs=2*cs;

		if(cs<1)
		{
			RSB_ERROR("#uhm. overflow ?\n");
			goto err;
		}
		mbw_m.so=sizeof(w_t);
		mbw_m.sz=cs;
		mbw_m.times=0;/* mbw_test will probe for a default reasonable time */
		mbw_m.cln=cln;
		mbw_m.cl=cl;
		mbw_m.hlcs = rsb__get_lnc_size(cln);
		if(mbw_m.hlcs<1)
			goto err;

		if(RSB_SOME_ERROR(mbw_test(&mbw_m)))
			goto err;

		memcpy( &(mbw_ms[cl-1]) ,&mbw_m,sizeof(struct rsb_mbw_m_t));
	}

	cm->mb=mbw_ms;
	cm->cln=cln;
	cm->extra_level=extra_levels;
	return RSB_ERR_NO_ERROR;
err:
	RSB_CONDITIONAL_FREE(mbw_ms);
	RSB_STDERR("An error occurred during memory benchmarking.\n");
	return RSB_ERR_GENERIC_ERROR;
}

static rsb_err_t rsb__print_mem_hier_timings(struct rsb_mbw_es_t *esp, const struct rsb_mbw_cm_t * cm)
{
	/**
	 * \ingroup gr_internals
	 * Quiet if esp != NULL.
	 */
	long cl;
	const long print_ratio=1;

	if(!cm)
		return RSB_ERR_BADARGS;

	for(cl=1;cl<=cm->cln+cm->extra_level;++cl)
	{
		if(!esp)
		{
			if(cl<=cm->cln)
				RSB_INFO("#Level %ld:\n",cl);
			else
				RSB_INFO("#Level %ld (RAM) (sample size 2^%ld times the last cache size):\n",cm->cln+1,cl-cm->extra_level);
		}
		mbw_printf(esp,&cm->mb[cl-1],cl);

		if(!esp)
		if(cl>1 && print_ratio)
			if(RSB_SOME_ERROR(mbw_ratio_printf(&cm->mb[cl-1],&cm->mb[cl-2])))
			{ RSB_ERROR(RSB_ERRM_ES); /* Note: may propagate error code */ }

		if( esp ) 
			esp += RSB_MB_N;
	}
	return RSB_ERR_NO_ERROR;
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_err_t rsb_tlb_benchmark(void)
{
	/**
		Unfinished code.

		Performance of this benchmark should expose shortcomings of memory fragmentation on performance.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t sz,psz,pn,wpp,times;
	const size_t K=1024;
	w_t * p=NULL;
	w_t c=0;
	rsb_int i,j;
	rsb_time_t t;

	RSB_WARN("TLB benchmark code is unfinished!\n");
	RSB_STDERR("#TLB benchmark.\n");
	for(sz=K*K/2;sz<K*K*K;sz*=2)
	{
		double mBps = 1.0;

		psz=4096;
		pn=sz/psz;
		wpp=psz/sizeof(w_t);
		times=100;
		p = rsb__aligned_malloc( sz , sz );
		if(!p)
			goto err;
		
		t = -rsb_time();
		for(j=0;j<times;++j)
		{
			rsb_time_t ft = -rsb_time();

			errval = rsb__flush_cache(0);
			if( RSB_SOME_ERROR(errval) ) 
				goto err;
			ft += rsb_time();
			t -= ft;
			for(i=0;i<pn;++i)
			{
				c+=p[i*wpp];
			}
		}
		t += rsb_time();
		RSB__FREE(p);

		mBps*=pn*sizeof(w_t);
		mBps/=t;
		mBps*=times;
		mBps/=1024*1024;

		RSB_STDERR("#TLB timing benchmark : scanned %zd entries spaced %zd bytes across %zd bytes in %lg s (%lg MBps)\n",pn,psz,sz,t,mBps);
	}

err:
	return errval;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

static rsb_err_t rsb_indirect_scan_benchmark(long ss, long * spiffero, long times, rsb_time_t *bt)
{
	/**
	*/
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
	rsb_time_t dt,rt,lt;
	rsb_coo_idx_t *IA=NULL;		/* the array to be scanned */
	rsb_coo_idx_t acc=0;			/* accumulator */
	rsb_nnz_idx_t *IP=NULL;		/* the array setting the scan order */
	void *CA=NULL;				/* the array setting the scan order */
	long els=0,fas=0;
	long i,ab,it;

	times = RSB_MAX(1,times);
	els=ss/(sizeof(rsb_coo_idx_t)),fas=4*ss;	/* the number of elements   */
	if(els<1 || fas<1)
		{ RSB_ERROR(RSB_ERRM_ES); goto err; }
	ab=sizeof(rsb_nnz_idx_t)*els+sizeof(rsb_coo_idx_t)*els;
	IP = rsb__malloc(sizeof(rsb_nnz_idx_t)*els);
	IA = rsb__malloc(sizeof(rsb_coo_idx_t)*els);
	CA = rsb__malloc(fas);
	if(!IP){RSB_ERROR(RSB_ERRM_ES);goto erri;}
	if(!IA){RSB_ERROR(RSB_ERRM_ES);goto erri;}
	if(!CA){RSB_ERROR(RSB_ERRM_ES);goto erri;}

	// random fill
	for(i=0;i<els;++i)
		IA[i]=rand()%els;
	// first phase: random scan
	for(i=0;i<els;++i)
		IP[i]=rand()%els;
	rsb__scan_cache(CA,fas,RSB_MB_FLUSH,RSB_FLUSH_TIMES,NULL);	/* flush cache */
	dt = - rsb_time();
	for(it=0;it<times;++it)
		for(i=0;i<els;++i)
			acc+=IA[IP[i]];
	dt += rsb_time();
	rt=dt/times;

	// second phase: linear scan
	for(i=0;i<els;++i)
		IP[i]=i;
	rsb__scan_cache(CA,fas,RSB_MB_FLUSH,RSB_FLUSH_TIMES,NULL);	/* flush cache */
	dt = - rsb_time();
	for(it=0;it<times;++it)
		for(i=0;i<els;++i)
			acc+=IA[IP[i]];
	dt += rsb_time();
	lt=dt/times;
	if(spiffero)
		RSB_INFO("for %ld elements, %ld bytes, random access time: %lg, linear access time: %lg, ratio %lg (%ld times)\n",els,ab,rt,lt,rt/lt,times);
	else
		;/* tuning mode only */
	if(spiffero)
		*spiffero+=acc;
	else
	{	RSB_INFO("ignore this value: %zd\n",(size_t)acc);}
	if(bt)
		*bt=(rt+lt)*times;
	errval = RSB_ERR_NO_ERROR;
erri:
	RSB_CONDITIONAL_FREE(CA);
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(IP);
err:
	return errval;
}

	static void rsb__memory_benchmark_p1_memcpy(void)
	{
		// Phase 1 of memory benchmark: MEMCPY-like.
		const long lcs = rsb__get_lastlevel_c_size();
		long times = RSB_MEMSCAN_MIN_TIMES, tinc = 1;
		const long wet = rsb_get_num_threads();
		const long fsm = rsb__sys_free_system_memory();
		int ci;
		const size_t wss = RSB_MIN(fsm / 3, 4*wet*lcs);
		long i;
		const long els = wss/sizeof(rsb_coo_idx_t);
		rsb_coo_idx_t *IS=NULL,*ID=NULL;
		rsb_time_t bt = RSB_REAL_ZERO,dt = RSB_REAL_ZERO;

		if(wss<1)
			goto errm;
		IS = rsb__malloc(wss);
		ID = rsb__calloc(wss);
		if(!IS || !ID)
			goto errm;
		for(i=0;i<els;++i)
			IS[i]=rand()%els;

		while( times<(RSB_MEMSCAN_MAX_TIMES/2) && bt<RSB_MEMSCAN_TIME )
		{
			int it;
			times += tinc;
			dt = - rsb_time();
			for(it=0;it<tinc;++it)
				RSB_A_MEMCPY_parallel(ID,IS,0,0,wss/RSB_CHAR_BIT,RSB_CHAR_BIT);
			dt += rsb_time();
			bt+=dt;
			tinc*=2;
		}
		RSB_INFO("# entering memory benchmark, phase 1 (%ldx repeated parallel MEMCPY of %zd bytes)\n", times, wss);

		if(0)
		{	RSB_WARN("first estimate of MEMCPY on %zd bytes: %lg GB/s (%ld times in %lg s)\n",(size_t)wss,
			((((double)wss)*times)/bt)/1.e9,times,bt);
		}
		for(i=0;i<els;++i)
			IS[i]=rand()%els;
		for(ci=1;ci<=wet;++ci)
		{
			int it;

			rsb__flush_cache(0);

			rsb__set_num_threads(ci);
			dt = - rsb_time();
			for(it=0;it<times;++it)
				RSB_A_MEMCPY_parallel(ID,IS,0,0,wss/RSB_CHAR_BIT,RSB_CHAR_BIT);
			dt += rsb_time();
			bt=dt;
			RSB_WARN("%zu cores MEMCPY on %zd bytes: %lg GB/s (%ld times in %lg s)\n",(size_t)ci,wss,
				((((double)wss)*times)/bt)/1.e9,times,bt);
		}
		rsb__set_num_threads(wet);
errm:
		RSB_CONDITIONAL_FREE(IS);
		RSB_CONDITIONAL_FREE(ID);
	}

	static void rsb__memory_benchmark_p2(void)
	{
		// Phase 2 of memory benchmark.
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		const long fcs = rsb__get_first_level_c_size();
		const long lcs = rsb__get_lastlevel_c_size();
		const long rcs = lcs;
		const long fsm = rsb__sys_free_system_memory();
		long spiffero=0,times = RSB_MEMSCAN_MIN_TIMES,tinc=1,reftimes=0;
		rsb_time_t bt = RSB_REAL_ZERO,dt = RSB_REAL_ZERO;

		while( times<(RSB_MEMSCAN_MAX_TIMES/2) && bt<RSB_MEMSCAN_TIME )
		{
			times += tinc;
			errval |= rsb_indirect_scan_benchmark(rcs,NULL,tinc,&dt);
			bt+=dt;
			tinc *= 2;
		}
		RSB_INFO("# entering memory benchmark, phase 2 (%ldx repeated parallel MEMCPY of %ld bytes)\n", times, fsm);
		RSB_WARN("begin experimental indirect array scan benchmark\n");
		reftimes = times;
		RSB_WARN("autotuning done. will proceed with presumably %lg s samples\n",bt);
#define RSB_MEMSCAN_TIMES_FROM_REF(reftimes,refsize,bufsize) \
	((refsize)<(bufsize)? \
	RSB_MAX(reftimes/((bufsize)/(refsize)),RSB_MEMSCAN_MIN_TIMES): \
	RSB_MAX(((refsize)/(bufsize))*reftimes,RSB_MEMSCAN_MIN_TIMES))

		times = RSB_MEMSCAN_TIMES_FROM_REF(reftimes,rcs,fcs);
		errval|=rsb_indirect_scan_benchmark(fcs,&spiffero,times,&bt);
		times = RSB_MEMSCAN_TIMES_FROM_REF(reftimes,rcs,(lcs-fcs)/2);
		errval|=rsb_indirect_scan_benchmark(fcs+(lcs-fcs)/2,&spiffero,times,&bt);
		times = RSB_MEMSCAN_TIMES_FROM_REF(reftimes,rcs,lcs);
		errval|=rsb_indirect_scan_benchmark(lcs,&spiffero,times,&bt);
		times = RSB_MEMSCAN_TIMES_FROM_REF(reftimes,rcs,4*lcs);
		errval|=rsb_indirect_scan_benchmark(4*lcs,&spiffero,times,&bt);
		times = RSB_MEMSCAN_TIMES_FROM_REF(reftimes,rcs,16*lcs);
		errval|=rsb_indirect_scan_benchmark(RSB_MIN(fsm,16*lcs),&spiffero,times/2,&bt);
		times = RSB_MEMSCAN_TIMES_FROM_REF(reftimes,rcs,32*lcs);
		errval|=rsb_indirect_scan_benchmark(RSB_MIN(fsm,32*lcs),&spiffero,times/4,&bt);
		times = RSB_MEMSCAN_TIMES_FROM_REF(reftimes,rcs,64*lcs);
		errval|=rsb_indirect_scan_benchmark(RSB_MIN(fsm,64*lcs),&spiffero,times/4,&bt);
		RSB_INFO("#please ignore this value: %ld\n",spiffero);
		RSB_INFO("end experimental indirect array scan benchmark\n");
#undef RSB_MEMSCAN_TIMES_FROM_REF
		// Note: we ignore errval.
	}

rsb_err_t rsb__memory_benchmark(struct rsb_mbw_et_t * mbetp)
{
	/**
	 * Benchmark the memory hierarchy.
	 */
	struct rsb_mbw_cm_t cm;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_DEBUG_ASSERT(rsb_global_session_handle.rsb_g_initialized == RSB_BOOL_TRUE);

	if(mbetp)
		goto onlyforrecord;

	rsb__memory_benchmark_p1_memcpy();

	rsb__memory_benchmark_p2();
	
#if RSB_OBSOLETE_QUARANTINE_UNUSED
	rsb_tlb_benchmark(); /* Note: temporarily here (to be completed). */
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

	RSB_INFO("# entering memory benchmark, phase 3\n");
onlyforrecord:	/* start here if only taking record */

	if(RSB_SOME_ERROR(rsb__mem_hier_timings(&cm)))
		goto err;

	if( ! mbetp ) /* if struct then no print needed */
		if(RSB_SOME_ERROR(rsb__print_mem_hier_timings(NULL,&cm)))
			goto err;

	if( mbetp )
	{
		RSB_DO_ERROR_CUMULATE(errval,rsb__mbw_es_fill(mbetp,&cm));
		rsb__do_perror(NULL,errval);
	}

	RSB_CONDITIONAL_FREE(cm.mb);
	if(RSB_SOME_ERROR(errval))
		goto err;

	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_GENERIC_ERROR;
}

rsb_err_t rsb__flush_cache(size_t sz)
{
	/**
	 Flush caches by repeated memory scans.
	 */
	void * fc=NULL;
	const size_t times = RSB_MIN_CACHE_FLUSH_SCAN_TIMES;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(sz==0)
	{
		sz = rsb__get_lastlevel_c_size();
		sz = RSB_MAX(sz,2*sz);
	}
	fc = rsb__calloc(sz);
	if(fc==NULL)
		return RSB_ERR_ENOMEM;
	errval = rsb__scan_cache(fc,sz,RSB_MB_FLUSH,times,NULL);	/* flush cache */
	RSB_CONDITIONAL_FREE(fc);	
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__mbw_es_print(const struct rsb_mbw_et_t * mbetp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_int_t sni;

	if(mbetp)
	{
		// RSB_STDOUT("Record comprises %d memory benchmark samples.\n",mbetp->sn);
                // RSB_STDOUT("Record comprises %d memory benchmark samples (each record %zd bytes).\n",mbetp->sn,sizeof(mbetp->et[0]));
                // RSB_STDOUT("offset of sz : %zd bytes\n",((void*)&(mbetp->et[0].sz ))-(void*)(&mbetp->et[0]));
                // RSB_STDOUT("offset of mbt: %zd bytes\n",((void*)&(mbetp->et[0].mbt))-(void*)(&mbetp->et[0]));
                // RSB_STDOUT("offset of lvl: %zd bytes\n",((void*)&(mbetp->et[0].lvl))-(void*)(&mbetp->et[0]));
                // RSB_STDOUT("offset of bw : %zd bytes\n",((void*)&(mbetp->et[0].bw ))-(void*)(&mbetp->et[0]));
		RSB_INFO("#%-32s\tsize\tlevel\tbw(MBps)\n","");
	}

	for ( sni = 0; sni < mbetp->sn; ++sni )
	{
		printf("%-32s\t%d\t%d\t%lg\n",rsb__mbw_s2s(mbetp->et[sni].mbt),mbetp->et[sni].sz,mbetp->et[sni].lvl,mbetp->et[sni].bw);
	}
	return errval;
}

rsb_err_t rsb__mbw_es_fill(struct rsb_mbw_et_t * mbetp, const struct rsb_mbw_cm_t * cm)
{
	/*
	 * accept an unused struct rsb_mbw_et_t; allocate and fill it.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!cm)
		return RSB_ERR_BADARGS;

	if(!mbetp || cm->cln == 0)
		goto ret;

	RSB_BZERO(mbetp,sizeof(*mbetp));
	mbetp->sn=(cm->cln+cm->extra_level)*RSB_MB_N;
	mbetp->et=rsb__calloc( mbetp->sn * sizeof(struct rsb_mbw_es_t) );

	RSB_INFO("Will fill struct with %d samples...\n",mbetp->sn);

	if ( ! mbetp->et )
	{
		errval = RSB_ERR_ENOMEM;
		goto ret;
	}


	if (rsb__print_mem_hier_timings(mbetp->et,cm))
		goto ret;

	return RSB_ERR_NO_ERROR;
ret:
	return errval;
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_err_t rsb__mbw_es_clone(struct rsb_mbw_et_t * mbetcp, const struct rsb_mbw_et_t * mbetsp)
{
	// unused.
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_BZERO(mbetcp,sizeof(*mbetcp));
	if( mbetsp->et )
	{
		mbetcp->et = rsb__clone_area(mbetsp->et,sizeof(*mbetsp->et)*mbetsp->sn);
		if(mbetcp->et)
			mbetcp->sn = mbetsp->sn;
		else
			errval = RSB_ERR_ENOMEM;
	}
	return errval;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__mbw_es_free(struct rsb_mbw_et_t * mbetp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mbetp)
		goto ret;

	RSB_CONDITIONAL_FREE(mbetp->et);	
	RSB_BZERO(mbetp,sizeof(*mbetp));
ret:
	return errval;
}

/* @endcond */
