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
/* @cond INNERDOC */
/**
 * @file
 * @brief 
 * @author Michele Martone
 * */
#ifndef RSB_TUNE_H_INCLUDED
#define RSB_TUNE_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb.h"
#include "rsb_internals.h"

/* FIXME: the following constants need a systematization... */
#define RSB_AUT0_TUNING_DEFAULT_TIME 10.0
#define RSB_AUT0_TUNING_SILENT 0
#define RSB_AUT0_TUNING_VERBOSE 1
#define RSB_AUT0_TUNING_QUATSCH 2
#define RSB_AT_TIME_AUTO 0.0
#define RSB_AT_NTIMES_AUTO 0
#define RSB_AT_MAX_TIME 3.0
#define RSB_TRACE_MAX_THREADS_P1 1 + RSB_CONST_MAX_SUPPORTED_THREADS
#define RSB_AT_THREADS_AUTO 0 /* see RSB_THREADS_AUTO (different meaning): this one may be used for auto+overwrite */
#define RSB_CONST_AUTO_TUNING_ROUNDS 1
#define RSB_CONST_MAX_TUNING_SUBROUNDS 5
#define RSB_CONST_AT_OP_SAMPLES_MIN  3 /* measurement steps in the autotuning (inner benchmarking loop) */
#define RSB_CONST_AT_OP_SAMPLES_MAX 10 /* measurement steps in the autotuning (inner benchmarking loop) */
#define RSB_CONST_MAX_TUNING_SAMPLES ((RSB_CONST_MAX_TUNING_ROUNDS)*(1+RSB_CONST_MAX_TUNING_SUBROUNDS))
#define RSB_AT_WANT_BEST_TIME 1 /* Autotuning shall base on 'best time' for a given (matrix, op, sampling run) */
#define RSB_DT_SAME_THREADS_TNP(TNP) ( (TNP)==NULL || *(TNP)== 0 ) /* On tnp==NULL will use the default thread count. */
#define RSB_DT_THREADS_TUNE_TNP(TNP) ( (TNP)!=NULL && *(TNP) < 0 ) /* On tnp!=NULL && *tnp<0  will probe different thread counts. */
#define RSB_DT_SPEC_THREADS_TNP(TNP) ( (TNP)!=NULL && *(TNP) > 0 ) /* On tnp!=NULL && *tnp>0  will use the given thread count. */
#define RSB_AT_DESTROYS_MTX 1 /* whether autotuning is allowed to destroy a suboptimal matrix after tuning */

#define RSB_MERGED_V11_AND_V12_AUTOTUNING 1
#if RSB_MERGED_V11_AND_V12_AUTOTUNING
#define RSB_CONST_DEF_MS_AT_AUTO_STEPS 6 /* 1.2 (merge/split clone) based autotuning */
#else
#define RSB_CONST_DEF_MS_AT_AUTO_STEPS 0 /* 1.1 (full clone) based autotuning */
#endif

#define RSB_FINISH_MIN_T(CUR_T,MIN_T,NOS) ( (MIN_T) > 0.0 ? ( (CUR_T) >  (MIN_T) ) : (NOS) )
#define RSB_FINISH_MIN_I(CUR_I,MIN_I,NOS) ( (MIN_I) > 0   ? ( (CUR_I) >= (MIN_I) ) : (NOS) )
#define RSB_REPEAT_MIN_T(CUR_T,MIN_T,NOS) ( (MIN_T) > 0.0 ? ( (CUR_T) <= (MIN_T) ) : (NOS) )
#define RSB_REPEAT_MIN_I(CUR_I,MIN_I,NOS) ( (MIN_I) > 0   ? ( (CUR_I) <  (MIN_I) ) : (NOS) )
#define RSB_REPEAT_MAX_T(CUR_T,MAX_T,NOS) ( (MAX_T) > 0.0 ? ( (CUR_T) <= (MAX_T) ) : (NOS) )
#define RSB_REPEAT_MAX_I(CUR_I,MAX_I,NOS) ( (MAX_I) > 0   ? ( (CUR_I) <  (MAX_I) ) : (NOS) )

#define RSB_WANT_REPEAT_TILL_MAXTIMES 0 /* if 0, and all specified min conditions met, terminate; if 1, require maxtimes */
#define RSB_REPEAT(CUR_T,CUR_I,MIN_T,MIN_I,MAX_T,MAX_I) ( \
	( /* conditions sufficient for continuation */ \
	  RSB_REPEAT_MIN_T((CUR_T),(MIN_T),0) || \
	  RSB_REPEAT_MIN_I((CUR_I),(MIN_I),0) \
       	) || \
	( /* conditions necessary for continuation */ \
	  (! ( RSB_FINISH_MIN_I((CUR_I),(MIN_I),1) && RSB_FINISH_MIN_T((CUR_T),(MIN_T),1) && ! RSB_WANT_REPEAT_TILL_MAXTIMES ) ) &&  \
	  RSB_REPEAT_MAX_I((CUR_I),(MAX_I),1) && \
	  RSB_REPEAT_MAX_T((CUR_T),(MAX_T),1) && \
	1 ) )

#define RSB_REPEAT_(CUR_T,CUR_I,MIN_T,MIN_I,MAX_T,MAX_I) ( \
	  RSB_REPEAT_MAX_T((CUR_T),(MAX_T),0)  \
 )

/* A substitute for obsolete RSB_REPEAT_MAX_T. 
 * To be called repeatedly after each loop.
 *
 * Input constants:
 *  IT: time before loop
 *  JF: value to use if DT is zero (e.g. due to timer caching)
 * Input variables:
 *  DT: delta time time (must be initialized to IT)
 *  TI: loop counter (must be initialized to 0)
 *  SS: square sum (must be initialized to 0)
 * Output variables:
 *  CT: current time
 *  TT: total time
 *  BT: best time
 * */
#define RSB_SAMPLE_STAT(IT,CT,DT,TT,BT,WT,SS,JF,TI) { \
		(CT) = rsb_time(); \
		(DT) = (CT) - (DT); \
		(BT) = RSB_MAX((JF),RSB_MIN((DT), (BT))); /* Initialize JF (jiffie) to the smallest dt and this will avoid having zeroes. */ \
		(WT) = RSB_MAX((DT), (WT)); \
		(TT) = (CT) - (IT); \
		(SS)+= (DT)*(DT); \
		(DT) = (CT); \
		(TI) ++; \
	}

#define RSB_SPEEDUP_TO_PCT(X) ((((double)(X))-1.0)*100.0) /* best coupled to a %6.1lf printf code; FIXME: shall use this format all thorough */
#define RSB_PCT(F,T) ((((double)(F)) / ((double)(T)))*100.0) /* FIXME: PCT -> PCNT */
#define RSB_STAT_DUMP(IT,TN,CT,DT,TT,BT,WT,SS,TI) { \
	RSB_STDOUT("%zd iterations (%d th.) took %0.4lgs; avg %0.4lgs ( +/- %6.2lf/%6.2lf %%); best %0.4lgs; worst %0.4lgs; std dev. %0.4lg (taking %s).\n", \
			(size_t)(TI),TN,TT,TT/TI, \
			RSB_PCT(((TT/TI)-(BT)),(TT/TI)), \
			RSB_PCT(((WT)-(TT/TI)),(TT/TI)), \
			BT,WT, \
			sqrt( (SS)/(TI) - (TT/TI)*(TT/TI) ),  \
			RSB_AT_WANT_BEST_TIME ? "best":"avg" \
			); \
	}

#define RSB_STAT_TAKE(IT,TN,CT,DT,TT,BT,WT,SS,TI,TSTP) if(TSTP) /* struct rsb_ts_t */ { \
		(TSTP)->avg = TT/TI;	\
		(TSTP)->min = BT;	\
		(TSTP)->max = WT;	\
		(TSTP)->sd  = sqrt( (SS)/(TI) - (TT/TI)*(TT/TI) );	\
		(TSTP)->ns = TI;	\
	}

#define RSB_STAT_DUMP_TS(TS) { RSB_STDOUT("%zd iterations took avg %0.4lgs ( +/- %6.2lf/%6.2lf %%); best %0.4lgs; worst %0.4lgs; std dev. %0.4lg.\n", (TS).ns, (TS).avg, RSB_PCT((((TS).avg)-((TS).min)),((TS).avg)), RSB_PCT((((TS).max)-((TS).avg)),((TS).avg)), (TS).min, (TS).max, (TS).sd); }

struct rsb_tattr_t /* threads auto tuning trace */
{
	rsb_time_t tpo[RSB_TRACE_MAX_THREADS_P1]; /**/
	rsb_int_t nit[RSB_TRACE_MAX_THREADS_P1]; /* number of iterations */
	rsb_real_t bpn; /* bytes per nonzeroes number */
	rsb_perf_t ofe; /* operation mflops estimate */
	rsb_perf_t btpo; /* best time per operation */
	rsb_perf_t dtpo; /* time per operation */
	rsb_time_t ttt; /* threads tuning time */
	rsb_int_t mint,maxt,optt,deft; /* min/max/optimal/default threads */
	/* what about bytes per non zero ? */
	/*const struct rsb_mtx_t*mtxAp;*/
	struct rsb_mtx_t mtxAc; /* TODO: rsb__blank_ptrs() */
	rsb_int_t vl; /* verbosity level */
};

struct rsb_attr_t /* auto tuning trace */
{
	rsb_bool_t dtr;  /* dump trace ?*/
	FILE*lfp,*pfp;	 /* log/plot file pointer */
	rsb_int_t trc,br;   /* tuning rounds count (max allowed RSB_CONST_MAX_TUNING_ROUNDS-1),best round */
	rsb_trans_t transA;
	rsb_int_t nrhs;
	rsb_int_t opid;  /* TODO: use this */
	rsb_type_t typecode;
	struct rsb_tattr_t clattr; /* competitor's library a.t.t. */
	struct rsb_tattr_t tattra[RSB_CONST_MAX_TUNING_SAMPLES]; /* threads auto tuning trace array */
       	char bname[RSB_MAX_FILENAME_LENGTH];
       	char mname[RSB_MAX_FILENAME_LENGTH];
};

void rsb__tattr_init(struct rsb_tattr_t* TTRP, const struct rsb_mtx_t*MTXAP, rsb_coo_idx_t nrA, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags, rsb_int_t nrhs);
void rsb__tattr_sets(struct rsb_tattr_t* ttrp, rsb_int_t dnt, rsb_int_t nt, rsb_time_t tpo, rsb_int_t bnt, rsb_int_t nits);
void rsb__attr_dump(struct rsb_attr_t*TTRP);

rsb_err_t rsb__tune_spxx( struct rsb_mtx_t ** mtxOpp, rsb_real_t *tsfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_int_t maxms, rsb_int_t maxss, rsb_int_t mintimes, rsb_int_t maxtimes, rsb_time_t maxt, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC, enum rsb_op_t op, rsb_int_t*epsp, rsb_time_t*otpopp, rsb_time_t*btpopp, int verbose, const char * fprfn, const char * mtxns, struct rsb_attr_t *attrp, struct rsb_ts_t*otposp, struct rsb_ts_t*btposp );
rsb_err_t rsb__do_tune_spmm(struct rsb_mtx_t ** mtxOpp, rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC);
rsb_err_t rsb__do_tune_spsm(struct rsb_mtx_t ** mtxOpp, rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC);
rsb_err_t rsb__do_bench_spxm(rsb_time_t *tpop, rsb_int_t *timesp, rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC, rsb_time_t maxdt, rsb_int_t mintimes, enum rsb_op_t op, rsb_int_t maxtimes, int verbose, rsb_int_t *tnp, struct rsb_ts_t * tstp);
rsb_err_t rsb__mtx_ms_check(struct rsb_mtx_t ** mtxOpp);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* RSB_TUNE_H_INCLUDED */

/* @endcond */
