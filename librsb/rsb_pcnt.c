/*                                                                                                                            

Copyright (C) 2008-2019 Michele Martone

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
 * @brief Perfomance counters related code (old).
 * @author Michele Martone
 * */

#include "rsb_internals.h"
#include "rsb-config.h"
#if RSB_HAVE_PTHREAD_H
#include <pthread.h>
#endif /* RSB_HAVE_PTHREAD_H */

#ifdef RSB_HAVE_PAPI
static const rsb_papi_int_t rsb_papi_eventlist [] = {
#if 0
	PAPI_L1_DCA,
	PAPI_L1_DCM,
	PAPI_L1_DCH,
	PAPI_L1_DCR,
	PAPI_L1_DCW,

	PAPI_L1_TCA,
	PAPI_L1_TCR,
	PAPI_L1_TCW,
	PAPI_L1_TCH,
	PAPI_L1_TCM,

	PAPI_L2_DCA,
	PAPI_L2_DCM,
	PAPI_L2_DCH,
	PAPI_L2_DCR,
#endif /* RSB_HAVE_PAPI */
#define RSB_WANT_PAPI_P4_COUNTERS 0
#if RSB_WANT_PAPI_P4_COUNTERS
	PAPI_L2_TCH,	/* p4 */
	PAPI_L2_TCM,	/* p4 */
	
//	PAPI_PRF_DM,	/* Prefetch data instruction caused a miss */
	PAPI_L1_LDM,	/* Level 1 load misses */
#if 0
	PAPI_L1_STM,
	PAPI_L2_LDM,
	PAPI_L2_STM,

	PAPI_L2_TCM,	/* p4 */
	PAPI_L2_TCH,	/* p4 */
#endif
	PAPI_TOT_CYC,	/* p4 */
	PAPI_TOT_INS,	/* p4 */
	PAPI_TLB_IM,	/* p4 */
#if 0
	PAPI_TLB_TL,
	PAPI_FMA_INS,	/* p4 */
#endif
//	PAPI_TOT_IIS,	/* p4 */
	PAPI_FP_INS,	/* p4 */
#if 1
//	PAPI_INT_INS,	/* p4 */
//	PAPI_LD_INS,	/* p4 */
//	PAPI_SR_INS,	/* p4 */
//	PAPI_VEC_INS,	/* p4 */
#endif
	PAPI_RES_STL,	/* p4 */
#endif
//	PAPI_FP_STAL,	/* p4 */
#if 0

	PAPI_FML_INS,
	PAPI_FAD_INS,
	PAPI_FDV_INS,
	PAPI_FSQ_INS,
	PAPI_FNV_INS,
#endif
//	PAPI_FP_OPS,	/* p4 */
#define RSB_WANT_PAPI_ATOM_COUNTERS 1
#if RSB_WANT_PAPI_ATOM_COUNTERS
	/*PAPI_L1_LDM,*/	/* Level 1 load misses */
	/* PAPI_L2_TCM, */	/* Level 2 load misses */
	PAPI_L1_TCM,PAPI_L2_TCM,	/* Level 1 load misses */
#endif /* RSB_WANT_PAPI_ATOM_COUNTERS */
	0
};
#define rsb_num_hwcntrs (sizeof(rsb_papi_eventlist)/sizeof(rsb_papi_int_t)-1)
#if 0
rsb_papi_long eventvals [ rsb_num_hwcntrs ];	/* FIXME */
rsb_papi_int_t EventSet=PAPI_NULL;
float real_time=0,proc_time=0,mflips=0;
float mflops=0,ips=0,ipc=0;/* FIXME : temporary */
char descr[PAPI_MAX_STR_LEN];

rsb_papi_long flpins=0,flpops=0;
rsb_papi_long ins=0;
const PAPI_hw_info_t *hwinfo = NULL;
#else
struct rsb_papi_stuff_t{
rsb_papi_int_t eventlist [rsb_num_hwcntrs];
rsb_papi_long eventvals [ rsb_num_hwcntrs ];	/* FIXME */
rsb_papi_int_t EventSet;
float real_time,proc_time,mflips;
float mflops,ips,ipc;/* FIXME : temporary */
char descr[PAPI_MAX_STR_LEN];
rsb_papi_long flpins,flpops;
rsb_papi_long ins;
const PAPI_hw_info_t *hwinfo;
};
static struct rsb_papi_stuff_t rps;
#endif
#define rsb_papi_mode 0

static rsb_err_t rsb_perf_counters_call(void)
{
#define E(M)  {RSB_STDERR("rsb_perf_counters_call : "M);/*goto err;*/}
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if((errval=PAPI_flips(&rps.real_time, &rps.proc_time, &rps.flpins, &rps.mflips))!=PAPI_OK)
	E("problem calling PAPI_flips()\n")
	if((errval=PAPI_flops(&rps.real_time, &rps.proc_time, &rps.flpops, &rps.mflops))!=PAPI_OK)
	E("problem calling PAPI_flops()\n")
	if((errval=PAPI_ipc(&rps.real_time, &rps.proc_time, &rps.ins, &rps.ipc))!=PAPI_OK)
	E("problem calling PAPI_ipc()\n")
#undef E

	RSB_DO_ERR_RETURN(errval)
	//err:
	/* should print the error code */
	//RSB_DO_ERR_RETURN(errval)
}
#endif

#if defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS>1)
rsb_err_t rsb_perf_counters_init(void)
#ifndef RSB_HAVE_PAPI
{
	return RSB_ERR_UNSUPPORTED_FEATURE;
}
#else /* RSB_HAVE_PAPI */
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_papi_err_t perrval=PAPI_OK;
	rps.EventSet=PAPI_NULL;
	RSB_BZERO_P(&rps);
	RSB_MEMCPY(&rps.eventlist,&rsb_papi_eventlist,sizeof(rsb_papi_eventlist));

#if 0
	/* Initialize the PAPI library and get the number of counters available */
	if ((rsb_num_hwcntrs = PAPI_num_counters()) <= PAPI_OK) 
		{errval = RSB_ERR_INTERNAL_ERROR;goto err;}

	/* FIXME :  RSB_STDOUT seems flawed!  */
	RSB_STDOUT("This system has %d available counters.\n", rsb_num_hwcntrs);

	if (rsb_num_hwcntrs > 2)
		rsb_num_hwcntrs = 2;

#endif

	if(rsb_papi_mode!=0)
	{
		int i;

		if (((perrval=PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) && perrval > 0 )
			{RSB_ERROR("PAPI_library_init() failed with code %d\n",(int)perrval);errval = RSB_ERR_INTERNAL_ERROR;goto err;}

		if ((rps.hwinfo = PAPI_get_hardware_info()) == NULL)
			{RSB_ERROR("PAPI_get_hardware_info() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}

#if 0
		if ( perrval=PAPI_set_granularity(PAPI_GRN_PROCG ) != PAPI_OK )
			{RSB_ERROR(RSB_ERRM_NL);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
		if ( perrval=PAPI_set_cmp_granularity(PAPI_GRN_PROCG /*PAPI_GRN_MAX*/,0) != PAPI_OK )
			{RSB_ERROR(RSB_ERRM_NL);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#endif
		/* Start counting events */
		if ((perrval=PAPI_create_eventset(&rps.EventSet)) != PAPI_OK)
			{RSB_ERROR("PAPI_create_eventset() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}

		for (i = 0; rps.eventlist[i] != 0; i++)
		{
			PAPI_event_code_to_name(rps.eventlist[i], rps.descr);
			if((perrval=PAPI_add_event(rps.EventSet,rps.eventlist[i])) != PAPI_OK)
			{
				//RSB_STDERR("PAPI_add_event(%d=%s) failed\n",rps.eventlist[i],descr);continue;
			}
			else
			{
				//RSB_STDERR("PAPI_add_event(%d=%s) ok\n",rps.eventlist[i],descr);
			}
	//		if(PAPI_remove_event(EventSet,rps.eventlist[i]) != PAPI_OK)
	//		{RSB_STDERR("PAPI_remove_event(%d) failed\n",rps.eventlist[i]);errval = RSB_ERR_INTERNAL_ERROR;continue;}
		}

		if ((perrval=PAPI_start(rps.EventSet)) != PAPI_OK)
		{
			RSB_ERROR("PAPI_start() failed : %s\n",PAPI_strerror(perrval));
			errval = RSB_ERR_INTERNAL_ERROR;goto err;
		}

	}
	else
	{
		int num_events = PAPI_num_counters();
		if(num_events<2)
			{RSB_ERROR("PAPI_num_counters() = %d < 2 \n",num_events);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
		//else
		//	RSB_STDERR("PAPI_num_counters() = %d\n",num_events);

		perrval = PAPI_library_init(PAPI_VER_CURRENT);
		if (perrval != PAPI_VER_CURRENT) {RSB_ERROR("PAPI_library_init() failed: %x\n",perrval);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#if 0
		// perrval = PAPI_set_cmp_granularity(PAPI_GRN_THD, 0);
		perrval = PAPI_set_cmp_granularity(PAPI_GRN_PROC, 0);
		/* perrval = PAPI_set_cmp_granularity(PAPI_GRN_SYS, 0); */
		if (perrval != PAPI_OK) {RSB_ERROR("PAPI_set_cmp_granularity() failed: %x\n",perrval);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
		// perrval = PAPI_create_eventset(&EventSet);
		// if (perrval != PAPI_OK) {RSB_ERROR(" \n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#endif
		/* PAPI_set_domain( PAPI_DOM_ALL ); */
		/* Start counting events */
		if ((perrval = PAPI_start_counters(rps.eventlist, rsb_num_hwcntrs)) != PAPI_OK)
			{RSB_ERROR("PAPI_start_counters() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}

#if RSB_WANT_OMP_RECURSIVE_KERNELS
#if RSB_HAVE_PTHREAD_H
		/* FIXME: in this version of librsb, PAPI threading support is NOT complete. */
		/* omp_get_thread_num is not adequate, see doc*/
	//	if (((perrval= PAPI_thread_init(omp_get_thread_num)) != PAPI_OK) )
		//if (((perrval= PAPI_thread_init(&pthread_self)) != PAPI_OK) )
		if (((perrval = PAPI_thread_init(pthread_self)) != PAPI_OK) )
			{RSB_ERROR("PAPI_thread_init() failed with code %d\n",(int)perrval);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#endif /* RSB_HAVE_PTHREAD_H */
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */


	}

	/* TODO : finish here */
//	rsb_perf_counters_call();

	RSB_STDERR("rsb_perf_counters_init: PAPI initialization ok.\n");
	return RSB_ERR_NO_ERROR;
err:
	{
		char pes[RSB_MAX_STRERRLEN];/* Flawfinder: ignore */
		/* PAPI_perror(perrval,pes,sizeof(pes)); */
		PAPI_perror("");
		pes[RSB_MAX_STRERRLEN-1]=RSB_NUL;
		RSB_ERROR("error in rsb_perf_counters_init(\"%s\")\n",pes);
	}
	rsb__do_perror(NULL,errval);
	/* error message printout */
	RSB_DO_ERR_RETURN(errval)
}
#endif /* defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS>1) */

rsb_err_t rsb_perf_counters_reset(void)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if 0
	if (PAPI_stop_counters ( eventvals, rsb_num_hwcntrs ) != PAPI_OK)
	{RSB_STDERR("PAPI_stop_counters() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}
	if (PAPI_start_counters ( eventvals, rsb_num_hwcntrs ) != PAPI_OK)
	{RSB_STDERR("PAPI_start_counters() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#else
	if (PAPI_reset ( rps.EventSet ) != PAPI_OK)
	{RSB_STDERR("PAPI_reset() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#endif
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb_perf_counters_dump(const rsb_char_t *premsg, const rsb_char_t *postmsg, rsb_int_t tdiv, struct rsb_pci_t *pcip)
{
	/*
	  \ingroup gr_internals
	 hardware counters information update
	 */
#if RSB_ALLOW_STDOUT
	const rsb_char_t *prmstr="";
	const rsb_char_t *pomstr="";
	int i;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(premsg)prmstr=premsg;
	if(postmsg)pomstr=postmsg;

#ifdef RSB_HAVE_PAPI
	if(pcip)
		RSB_BZERO_P(pcip);
	if(pcip)
	for (i = 0; i < RSB_MIN(rsb_num_hwcntrs,RSB_PC_MAX_ITEMS); i++) /* FIXME: shall put a stringent limit */
	{
		pcip->eventvals[i]=rps.eventvals[i]/tdiv;
		pcip->eventlist[i]=rps.eventlist[i];
		PAPI_event_code_to_name(pcip->eventlist[i],pcip->eventdesc[i]);
		pcip->eventnum=i+1; /* ehm...  */
	}
#endif /* RSB_HAVE_PAPI */
	else
	/*RSB_STDOUT("summary of PAPI output: \n");*/
	for (i = 0; i < rsb_num_hwcntrs; i++)
	{	
		PAPI_event_code_to_name(rps.eventlist[i], rps.descr);
		/*RSB_STDOUT("counters values : %s%s %lld\n",prmstr,descr,eventvals[i]);*/
		RSB_STDOUT("%s%s:\t%lld%s\n",prmstr,rps.descr,rps.eventvals[i]/tdiv,pomstr);
	}
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif /* RSB_ALLOW_STDOUT */
}

rsb_err_t rsb_perf_counters_update(void)
{
	/**
	  \ingroup gr_internals
	   
	 */

	/* hardware counters information update */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	/* fixme ! this accumulates, and does not reset! */
	/* unlike said in http://icl.cs.utk.edu/projects/papi/files/html_man3/papi_read_counters.html */
	if (PAPI_read_counters ( rps.eventvals, rsb_num_hwcntrs ) != PAPI_OK)
	{RSB_STDERR("PAPI_read_counters() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}

//	if (PAPI_accum_counters( eventvals, rsb_num_hwcntrs ) != PAPI_OK)
//	{RSB_STDERR("PAPI_accum_counters() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}

err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_HAVE_PAPI */

#if defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS>1)
rsb_err_t rsb_perf_counters_finalize(void)
#ifndef RSB_HAVE_PAPI
{ return RSB_ERR_UNSUPPORTED_FEATURE; }
#else /* RSB_HAVE_PAPI */
{
	/* will finalize and maybe dump some performance info */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

//	RSB_STDERR("rsb_perf_counters_finalize..\n");

	/* WRITE ME */
	/* TODO : finish here */
//	rsb_perf_counters_call();

#if 0
	{
	RSB_STDERR("sizeof(rsb_papi_long)  %d\n",sizeof(rsb_papi_long));
	RSB_STDERR("real time          %g\n",real_time);
	RSB_STDERR("processor time     %g\n",proc_time);
	RSB_STDERR("fl.p. instructions %d\n",flpins);
	RSB_STDERR("fl.p. ops          %d\n",flpops);
	RSB_STDERR("mflips             %g\n",mflips);
	RSB_STDERR("mflops             %g\n",mflops);
	RSB_STDERR("instructions       %d\n",ins);
	RSB_STDERR("ins. per cycle     %g\n",ipc);
	RSB_STDERR("ins. per second    %g\n",ips);
	}
#endif

	RSB_BZERO(rps.eventvals,sizeof(rps.eventvals));

	//rsb_perf_counters_update();
	//rsb_perf_counters_dump();

	if(rsb_papi_mode!=0)
	{
		int code=0;
		if (PAPI_stop(rps.EventSet, rps.eventvals ) != PAPI_OK)
		{RSB_STDERR("PAPI_stop() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}

		if (PAPI_cleanup_eventset(rps.EventSet ) != PAPI_OK)
		{RSB_STDERR("PAPI_cleanup_eventset() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}

		if ((code=PAPI_destroy_eventset(&rps.EventSet)) != PAPI_OK)
		{
			RSB_STDERR("PAPI_destroy_eventset() failed : %s\n",PAPI_strerror(code));
			errval = RSB_ERR_INTERNAL_ERROR;goto err;
		}
	}
	else
	{
		if (PAPI_stop_counters ( rps.eventvals, rsb_num_hwcntrs ) != PAPI_OK)
		{RSB_STDERR("PAPI_stop_counters() failed\n");errval = RSB_ERR_INTERNAL_ERROR;goto err;}
	}


	goto ok;
err:
	RSB_STDERR("error in rsb_perf_counters_finalize()\n");
	rsb__do_perror(NULL,errval);
	/* error message printout */
ok:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_HAVE_PAPI */
#endif /* defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS>1) */

rsb_err_t rsb_hc_main()		/* preliminary */
{
	/**
	  \ingroup gr_internals
	   UNFINISHED
	   A miniprogram for preliminary play with hardware counters.
	  
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_int i,j;
	//register double f=0;
	//register float f=0;
	//register char f=0;
	register int f=0;


	errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS);	
	if(RSB_SOME_ERROR(errval))
		goto err;
	
	for(j=0;j<10;++j)
	{
		// see diff for 1.0 or 2.0 add !
		//for(i=0;i<100000;++i) f=f+1.0f;
		for(i=0;i<100000;++i) f=f+1;
		//for(i=0;i<100000;++i) f=f+2.0f;
		//rsb_perf_counters_update();
		//rsb_perf_counters_dump();
		//rsb_perf_counters_reset();
	}

	errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);	
	if(RSB_SOME_ERROR(errval))
		goto err;

	RSB_STDERR("ignore this printout :) hc: %lf\n",(double)f);

err:
	RSB_DO_ERR_RETURN(errval)
}



/* @endcond */
