dnl
dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
dnl
/*! 
 @file
 @brief 

 Matrix Operations testing code source file.
 This is NOT part of the library: only of companion programs.

 */
dnl
include(`rsb_misc.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
include(`rsb_krnl_macros.m4')dnl
include(`libspblas_macros.m4')dnl
dnl
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_TEST_MATOPS_H_INCLUDED
#define RSB_TEST_MATOPS_H_INCLUDED
')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
#include "rsb_test_matops.h"
')dnl
dnl

/* FIXME: necessary, until we use so many #ifdefs in this program */
#include "rsb-config.h"
#include "rsb_common.h"
#include "rsb_mkl.h"
#if RSB_WANT_ARMPL
#include <armpl.h>
#endif /* RSB_WANT_MKL */

ifdef(`ONLY_WANT_HEADERS',`dnl
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

')dnl
dnl
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
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
#define RSB_UTIL_CSR_IDX_OCCUPATION(R,C,NNZ) (sizeof(rsb_coo_idx_t)*nnz+sizeof(rsb_nnz_idx_t)*nrA)
#define RSB_UTIL_COO_IDX_OCCUPATION(R,C,NNZ) (sizeof(rsb_coo_idx_t)*2*nnz)
#define RSB_UTIL_COO_OCCUPATION(R,C,NNZ,TYPE) (RSB_UTIL_COO_IDX_OCCUPATION(R,C,NNZ)+(NNZ)*(RSB_SIZEOF(TYPE)))
#define RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH() RSB_FPRINTF_MATRIX_ESSENTIALS(stdout,mtxAp,filename,cc) 
#define RSB_DIV(Q,D) ( ( (Q)+(D)-1 ) / (D) )
')dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
RSB_INTERNALS_COMMON_HEAD_DECLS
#define RSB_NEGATED_EXAGGERATED_TUNER_TIMES -999999.0
#define RSB_MKL_APPROPRIATE_AT_TIME_SPEC(TS) ( (TS) != RSB_NEGATED_EXAGGERATED_TUNER_TIMES )
#define RSB__APPROPRIATE_AT_TIME_SPEC(TS) ( (TS) != RSB_NEGATED_EXAGGERATED_TUNER_TIMES )
RSB_INTERNALS_RSBENCH_HEAD_DECLS
#define RSBENCH_MAY_SQUIT(LABEL,ACTION) { if(RSB_SHALL_QUIT) { RSB_INFO("Terminating execution earlier due to interactive user request.\n"); ACTION; goto LABEL; } }
#define RSBENCH_MAY_TQUIT(LABEL,ACTION) { if(maxtprt > RSB_TIME_ZERO && maxtprt < rsb_time()+totprt) { RSB_INFO("Terminating execution earlier due to user set max timer of %2.3lg s.\n",maxtprt); ACTION; goto LABEL; } }
')dnl
dnl
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

ifdef(`ONLY_WANT_HEADERS',`',`dnl
#ifdef RSB_HAVE_REGEX_H 
#include <regex.h>
#endif /* RSB_HAVE_REGEX_H */
#define RSBENCH_STDERR RSB_STDERR
#define fim_strchr strchr
')dnl

#if defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS==1)
#define RSB_WANT_PERFORMANCE_COUNTERS_IN_RSBENCH 1
#else
#define RSB_WANT_PERFORMANCE_COUNTERS_IN_RSBENCH 0
#endif

ifdef(`ONLY_WANT_HEADERS',`',`dnl
static int rsb__echo_cargs(const int argc, rsb_char_t * const argv[])
{
	int argci;

	if(argc > 0)
		RSBENCH_STDOUT("# %s",argv[0]);
	for(argci=1; argci<argc; ++argci)
	{
		RSBENCH_STDOUT(" %s",argv[argci]);
	}
	RSBENCH_STDOUT("\n");
	return 0;
}
')dnl

ifdef(`ONLY_WANT_HEADERS',`',`dnl
#ifdef RSB_HAVE_REGEX_H 
static	rsb_bool_t rsb_regexp_match(const rsb_char_t*s, const rsb_char_t*r)
	{
		regex_t regex;
		const int nmatch = 1;
		regmatch_t pmatch[nmatch];
		rsb_bool_t match = RSB_BOOL_FALSE;
		const int ignorecase = 0;
		const int ignorenewlines = 0;

		if(!r || !strlen(r))
			goto ret;

		if(regcomp(&regex,r, 0 | REG_EXTENDED | (ignorecase==0?0:REG_ICASE) )!=0)
		{
			RSB_ERROR("error calling regcomp; invalid regexp: %s\n",s);
			goto ret;
		}

		if(regexec(&regex,s+0,nmatch,pmatch,0)!=REG_NOMATCH)
		{
			match = RSB_BOOL_TRUE;
		}
		regfree(&regex);
ret:
		return match;
	}
#endif /* RSB_HAVE_REGEX_H */
#
       static int rsb__cmpstringp(const void *p1, const void *p2)
       {
		return strcmp(* (char * const *) p1, * (char * const *) p2);
       }
')dnl

ifdef(`ONLY_WANT_HEADERS',`',`dnl
static void rsb__echo_timeandlabel(const char*l, const char*r, rsb_time_t *stp)
{
	rsb_time_t ct = rsb_time();

	if(stp && *stp)
		RSBENCH_STDOUT("#%s%.0lf (after %.1lfs of w.c.t.)%s",l?l:"",ct,ct-*stp,r?r:"");
	else
		RSBENCH_STDOUT("#%s%.0lf%s",l?l:"",ct,r?r:"");
	if(stp)
		*stp = ct;
}
')dnl

ifdef(`ONLY_WANT_HEADERS',`',`dnl
static void rsb__impcdstr(char * dst, const char * h, const char * pp, const char * ap, rsb_bool_t with_mkl, const char * tcrprs, const char * fnrprs, const char *t)
{
	/* There is some overlap with rsb__cat_compver and rsb__sprint_matrix_implementation_code that shall be resolved. */
	rsb_char_t buf[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];/* Flawfinder: ignore */

	rsb__cat_compver(buf);
#if RSB_WANT_MKL
	if(with_mkl)
		rsb__print_mkl_version(buf+strlen(buf), RSB_BOOL_FALSE);
#endif /* RSB_WANT_MKL */
	strcat(buf,"");
	rsb__sprintf(dst,"%s%s_%s_%.0lf_%s%s" "%s%s%s" "%s%s" "%s",pp?pp:"",h,rsb__getenv_nnr("HOSTNAME"),rsb_time(),buf,ap?ap:"",
		(tcrprs?"-":""), (tcrprs?tcrprs:""), (tcrprs?"th":""),
		(fnrprs&&*fnrprs?"-":""), (fnrprs?fnrprs:""), t );
}
')dnl

int rsb_test_help_and_exit(const rsb_char_t *argv0, const rsb_option_t *o, int code)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	    size_t i=0;

            printf("%s %s",argv0," where OPTIONS are taken from :\n");
            for(i=0;o[i].val;++i)
            {
#if RSB_HAVE_GETOPT_H
#endif /* RSB_HAVE_GETOPT_H */
		const rsb_option_t ro = o[i];
                if(ro.val<RSB_MAX_VALUE_FOR_TYPE(rsb_char_t) && isprint(ro.val)  )/* please do not swap conditions : some isprint() implementations segfault on this */
		{
                	printf("\t-%c",(rsb_char_t)(ro.val));
		}
		else
			printf("\t");
                printf("\t\t");
		if(ro.name)
	                printf("--%s",ro.name);
                switch(ro.has_arg)
		{
	                case no_argument:
	                break;
	                case required_argument:
	                printf(" <arg>");
	                break;
	                case optional_argument:
	                printf(" [=arg]");
	                break;
	                default:
        	        ;
                };
                printf("\n");
	    }
            printf("\n");
	    printf("Arguments to --want-autotune of the format \"%s\", where S is the autotuning time in seconds, X is the number of tries, T the number of starting threads, V can be either q for quiet autotuning or v for a verbose one (can be specified twice). Valid examples: 3.0s2x4tv, 3.0s2x0tq, 3.0s, 2.0s10x . See documentation of rsb_tune_spmm for a full explanation of these parameters role in auto-tuning.\n",RSB_WAT_FMT_H);
            printf("Report bugs to %s.\n",RSB_PACKAGE_BUGREPORT);
            return code;
}
')

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
#ifdef RSB_HAVE_UNISTD_H
	extern char **environ;
#endif /* RSB_HAVE_UNISTD_H */
')dnl
dnl
dnl
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	#define RSBENCH_MEM_ALLOC_INFO(LT) RSBENCH_STDOUT("( allocated_memory:%zd allocations_count:%zd allocations_cumulative:%zd)%s",rsb_global_session_handle.allocated_memory,rsb_global_session_handle.allocations_count,rsb_global_session_handle.allocations_cumulative,LT);
#else
	#define RSBENCH_MEM_ALLOC_INFO(LT)
#endif

#if RSB_WANT_ARMPL

define(`RSB_M4_RSB_TYPE_TO_ARMPL_TYPE',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(type,`double complex',`armpl_doublecomplex_t')`'dnl
ifelse(type,`float complex',`armpl_singlecomplex_t')`'dnl
ifelse(type,`long double',`long double')`'dnl
ifelse(type,`double',`double')`'dnl
ifelse(type,`float',`float')`'dnl
ifelse(type,`int',`MKL_INT')`'dnl
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
int rsb_rsb_to_armpl_trans(rsb_trans_t transA)
{
	/**
	 * \ingroup gr_internals
	 */
	switch(transA)
	{
		case(RSB_TRANSPOSITION_N):
		return ARMPL_SPARSE_OPERATION_NOTRANS;
		break;
		case(RSB_TRANSPOSITION_T):
		return ARMPL_SPARSE_OPERATION_TRANS;
		break;
		case(RSB_TRANSPOSITION_C):
		return ARMPL_SPARSE_OPERATION_CONJTRANS;
		break;
		default:
		return ARMPL_SPARSE_OPERATION_NOTRANS;
	}
}
')dnl
dnl

ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb__do_armpl_csr_spmv(armpl_spmat_t armpl_mat, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)dnl
{
	armpl_status_t info = ARMPL_STATUS_INPUT_PARAMETER_ERROR;

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
		info = `armpl_spmv_exec_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`'(rsb_rsb_to_armpl_trans(transA),*(RSB_M4_RSB_TYPE_TO_ARMPL_TYPE(mtype)*)alphap,armpl_mat,(RSB_M4_RSB_TYPE_TO_ARMPL_TYPE(mtype)*)x,*(RSB_M4_RSB_TYPE_TO_ARMPL_TYPE(mtype)*)betap,(RSB_M4_RSB_TYPE_TO_ARMPL_TYPE(mtype)*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	if ( info == ARMPL_STATUS_SUCCESS )
		return RSB_ERR_NO_ERROR;
	return RSB_ERR_INTERNAL_ERROR;
}
')dnl
dnl

rsb_err_t rsb__armpl_csr_spmv_bench(const void *VA, const armpl_int_t m, const armpl_int_t k, const armpl_int_t nnz, const armpl_int_t * IP, const armpl_int_t *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp,rsb_bool_t want_at)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
#define RSB_ARMPL_THREADS_TUNING_ODECLS					\
		const rsb_time_t tinf = RSB_CACHED_TIMER_GRANULARITY;	\
		rsb_time_t best = RSB_CONST_IMPOSSIBLY_BIG_TIME;	\
		rsb_thread_t ont = 1;		\
		rsb_thread_t nt, lnt = 1, unt = 1; \
		rsb_thread_t otn = ont;					\
		rsb_thread_t dtn = 1;

#define RSB_ARMPL_THREADS_TUNING_IDECLS									\
		rsb_time_t it = rsb_time(), ct = RSB_TIME_ZERO;	/* initial/current time */	\
		rsb_time_t dt = it, tt = RSB_TIME_ZERO; /* elapsed (delta) / total  time */	\
		rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO; /* best / worst  time */	\
		rsb_time_t ss = RSB_TIME_ZERO; /* sum of squares */				\
		rsb_time_t mint = RSB_TIME_ZERO; /* minimal time */				\
		rsb_int_t times = 0; \
		const mintimes = RSB_CONST_AT_OP_SAMPLES_MIN, maxtimes = RSB_CONST_AT_OP_SAMPLES_MAX;	\
		rsb_time_t maxt = RSB_AT_MAX_TIME/* RSB_MKL_MAX_AT_TIME*/;

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t dt, tpo;
	rsb_time_t aplat_t;
	armpl_spmat_t armpl_mat;
	const int creation_flags = 0;
	armpl_status_t info = ~ ARMPL_STATUS_SUCCESS;

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ) // Note: only SDCZ
		info = `armpl_spmat_create_csr_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`(&armpl_mat, m, k, IP, JA, VA, creation_flags)';
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	if (info!=ARMPL_STATUS_SUCCESS) RSB_ERROR("ERROR: armpl_spmat_create_csr_d returned %d\n", info);

	info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_STRUCTURE, ARMPL_SPARSE_STRUCTURE_UNSTRUCTURED);
	if (info!=ARMPL_STATUS_SUCCESS) RSB_ERROR("ERROR: armpl_spmat_hint returned %d\n", info);

	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN) )
	{
		info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_STRUCTURE, ARMPL_SPARSE_STRUCTURE_HERMITIAN);
		if (info!=ARMPL_STATUS_SUCCESS) RSB_ERROR("ERROR: armpl_spmat_hint returned %d\n", info);
	}
	else
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC) )
	{
		info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_STRUCTURE, ARMPL_SPARSE_STRUCTURE_SYMMETRIC);
		if (info!=ARMPL_STATUS_SUCCESS) RSB_ERROR("ERROR: armpl_spmat_hint returned %d\n", info);
	}

	// See also:
	// ARMPL_SPARSE_STRUCTURE_DENSE
	// ARMPL_SPARSE_STRUCTURE_TRIANGULAR
	// ARMPL_SPARSE_STRUCTURE_HPCG

	if( want_at == RSB_BOOL_TRUE )
	{
		aplat_t = -rsb_time();

		info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_SPMV_OPERATION, ARMPL_SPARSE_OPERATION_NOTRANS);
		if (info!=ARMPL_STATUS_SUCCESS) RSB_ERROR("ERROR: armpl_spmat_hint returned %d\n", info);

		info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_SPMV_INVOCATIONS, ARMPL_SPARSE_INVOCATIONS_MANY);
		if (info!=ARMPL_STATUS_SUCCESS) RSB_ERROR("ERROR: armpl_spmat_hint returned %d\n", info);

		info = armpl_spmv_optimize(armpl_mat);
		if (info!=ARMPL_STATUS_SUCCESS) RSB_ERROR("ERROR: armpl_spmv_optimize returned %d\n", info);

		aplat_t += rsb_time();
	}

	if(otnp)
	{
		RSB_ARMPL_THREADS_TUNING_ODECLS	

		for(nt=lnt;nt<=unt;++nt)
		{
			/* TODO: change number of threads here (now is current OMP thread count) */
			RSB_ARMPL_THREADS_TUNING_IDECLS
			do
			{
				errval = rsb__do_armpl_csr_spmv(armpl_mat, x, y, alphap, betap, transA, typecode, flags);
				RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,tinf,times);
			}
			while(RSB_REPEAT(ct-it,times,mint,mintimes,maxt,maxtimes));

			dt = bt;
			if( dt < best )
			{
				otn = nt;
				best = RSB_MIN_ABOVE_INF(best,dt,tinf);
				RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp);
			}
			rsb__tattr_sets(ttrp,dtn,nt,dt,otn,times);/* FIXME: if no threads tuning, shall set dtpo = btpo, as well as ttrp.optt=0 */
			if(dtn == nt) RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp+1);
		}
		ttrp->ttt += rsb_time(); /* ttrp->ttt = tt; */
		tpo = best; /* tpo = 1.0 / best; */
		*otnp = otn;
	}
	else
	{
		dt = -rsb_time();
		errval = rsb__do_armpl_csr_spmv(armpl_mat, x, y, alphap, betap, transA, typecode, flags);
		dt += rsb_time();
		/* tpo = 1.0 / dt; */
		tpo = dt;
	}
	if(tpop)
		*tpop = tpo;
	return errval;

#undef RSB_ARMPL_THREADS_TUNING_ODECLS
#undef RSB_MKL_THREADS_TUNING_IDECLS

} /* rsb__armpl_csr_spmv_bench */
')dnl
dnl

#endif /* RSB_WANT_ARMPL */

define(`RSB_M4_COMPLETE_TEST_PROGRAM_FUNCTION',`dnl
pushdef(`mop',$1)dnl
int rsb__main_block_partitioned_`'mop`'(const int argc, rsb_char_t * const argv[])
ifdef(`ONLY_WANT_HEADERS',`dnl
;',`dnl
{
#if RSB_HAVE_GETOPT_H
	/*!
	 * \ingroup gr_bench
	 * This function implements a complete program for using our variable block
	 * rows sparse matrix storage as it was a fixed block size format.
	 * It is useful for benchmark against fixed block sparse matrix codes.
	 *
	 * This function will benchmark the "mop" matrix operation.
	 * */

	/*
	 * This example main program reads in a Matrix Market file in block format and multiplies it against a unit vector.
	 **/

	const struct rsb_option_struct rsb_options[] = {
	    {{"all-flags",	no_argument, NULL, 0x51}},/* Q */
	    {{"all-formats",	no_argument, NULL, 0x616c666f}},/* alfo */
	    {{"all-blas-opts",	no_argument, NULL, 0x616c626f}},/* albo */
	    {{"all-blas-types",	no_argument, NULL, 0x616c6274}},/* albt */
	    {{"allow-any-transposition-combination",	no_argument, NULL, 0x61617463 }},/* aatc */
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	    {{"alpha",	required_argument, NULL , 0x414C}},/* AL */
')dnl
	    {{"alternate-sort",	required_argument, NULL , 0x4153}},/* AS */
	    {{"auto-blocking",	no_argument, NULL, 0x41}},/* A */
	    {{"be-verbose",		no_argument, NULL, 0x76}},	/* v */
	    {{"bench",	no_argument , NULL, 0x626e6368}},	/* bnch */
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	    {{"beta",	required_argument, NULL ,  0x4246}},/* BE */
')dnl
	    {{"block-columnsize",	required_argument, NULL, 0x63}},/* c */
	    {{"block-rowsize",   required_argument, NULL, 0x72 }},/* r */
	    {{"cache-blocking",	required_argument, NULL , 0x4342}},/* CB */
/*	    {{"cache-flush",	no_argument, NULL, 0x4343}},*/ /*   */
	    {{"chdir",	required_argument, NULL, 0x6364}},/* cd */
	    {{"column-expand",	required_argument, NULL, 0x6B}},/* k */
	    {{"compare-competitors",	no_argument, NULL, 0x6363}},/* cc */
	    {{"no-compare-competitors",	no_argument, NULL, 0x776e6d6c }},/* wnml */  /* want no mkl */
	    {{"convert",	no_argument, NULL, 0x4B}},/* K */
/*	    {{"convert",	required_argument, NULL, 0x4B}},*//* K   */
dnl	    {{"dense",	required_argument, NULL, 0xbabb0 }},   /* */
	    {{"dense",	required_argument, NULL, 0x64 }},   /* d */
	    {{"diagonal-dominance-check",	no_argument , NULL, 0x4444}},/* DD */  /* new */
dnl	    {{"dump-profile",	no_argument, NULL, 0x4F}},/* O */
	    {{"dump-n-lhs-elements",	required_argument , NULL, 0x444444}},/* DDD */  /* new */
	    {{"echo-arguments",	no_argument , NULL, 0x6563686f}},/* echo */  /* new */
ifelse(mop,`mat_stats',`dnl
	    {{"estimate-samples",		required_argument, NULL, 0x53}},	/* S */
dnl 	    {{"until-confidence",required_argument, NULL, 0x75}},	/* u */	/* dead option */
	    {{"estimate-fillin",required_argument, NULL, 0x65}},	/* e */
')dnl
	    {{"flush-cache-in-iterations",	no_argument, NULL, 0x4343}},/*  */
	    {{"impatient",	no_argument, NULL, 0x696d7061}},/* impa[tient] */
	    {{"no-flush-cache-in-iterations",	no_argument, NULL, 0x434E}},/*  */
	    {{"flush-cache-around-loop",	no_argument, NULL, 0x434343}},/*  */
	    {{"want-ancillary-execs",	no_argument, NULL, 0x767646}},/*  */
	    {{"no-want-ancillary-execs",	no_argument, NULL, 0x42767646}},/*  */
	    {{"no-flush-cache-around-loop", no_argument	, NULL, 0x43434E}},/*  */
	    {{"want-no-recursive",	no_argument, NULL, 0x776e720a}},/*  */
	    {{"want-memory-benchmark",	no_argument, NULL, 0x776d62}},/* wmb */
	    {{"want-no-memory-benchmark",	no_argument, NULL, 0x776e6d62}},/* wnmb */
	    {{"nmb",	no_argument, NULL, 0x776e6d62}},/* wnmb */
	    {{"guess-blocking",	no_argument , NULL, 0x47}},/* G */
	    {{"help",	no_argument , NULL, 0x68}},	/* h */
	    {{"ilu0",	no_argument , NULL, 0x494B55}},/* ILU */  /* new */
	    {{"inc",	required_argument, NULL, 0x696e632a }},/* inc* */
	    {{"incx",	required_argument, NULL, 0xb1bb0 }},/* */
	    {{"incy",	required_argument, NULL, 0xb1bb1 }},/* */
	    {{"in-place-assembly-experimental",	no_argument , NULL, 0x6970}},/* i */
	    {{"in-place-csr",	no_argument, NULL, 0x69}},/* i */
	    {{"in-place-permutation",	no_argument, NULL, 0x50}},   /* P */
#if RSB_WITH_LIKWID
	    {{"likwid",	no_argument, NULL, 0x6c696b77}},   /* likw */
#endif /* RSB_WITH_LIKWID */
dnl	    {{"lower",	required_argument, NULL, 0xbabb1 }},   /* */
	    {{"lower",	required_argument, NULL, 0x6c}},   /* l */
	    {{"lower-dense",	required_argument, NULL, 0x6c64}},   /* ld */
	    {{"generate-lowerband",	required_argument, NULL, 0x6c6c}},   /* ll */
	    {{"gen-lband",	required_argument, NULL, 0x6c6c}},   /* ll */
	    {{"generate-spacing",	required_argument, NULL, 0xbabb2 }},   /* */
	    {{"matrix-dump",	no_argument, NULL, 0x44044}},/* D */
	    {{"matrix-dump-graph",	required_argument , NULL, 0x44047}},/* DG */
	    {{"matrix-dump-internals",	no_argument, NULL, 0x49049}},/* I 0x0 I */
	    {{"merge-experimental",	required_argument , NULL, 0x6d656578}},/* meex */
	    {{"split-experimental",	required_argument , NULL, 0x73706578}},/* spex */
	    {{"ms-experimental",	required_argument , NULL, 0x6d736578}},/* msex */
	    {{"matrix-filename",	required_argument, NULL, 0x66}},/* f */
	    {{"matrix-sample-pcnt",	required_argument, NULL, 0x6d747873}},/* mtxs */
	    {{"matrix-storage",	required_argument, NULL, 0x46}},/* F */
	    {{"matrix-time",	no_argument, NULL, 0x4D}},/* M */  /* new */
	    {{"mem-hierarchy-info",	required_argument , NULL, 0x4D4D}},/* MM */  /* new */
	    {{"max-runtime",	required_argument , NULL, 0x6d617275}},/* maru */
	    {{"no-op",		no_argument, NULL, 0x4E}},	/* N */
	    {{"notranspose",	no_argument, NULL, 0x5051}},   /* do not transpose the operation */
	    {{"no-transpose",	no_argument, NULL, 0x5051}},   /* do not transpose the operation */
	    {{"nrhs",	required_argument, NULL, 0x6e726873}},   /* */
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),1,`dnl
	    {{"nrhs-by-rows",	no_argument, NULL, 0x726f7773}},   /* */
	    {{"by-rows",	no_argument, NULL, 0x726f7773}},   /* */
	    {{"nrhs-by-columns", no_argument, NULL, 0x636f6c73}},   /* */
	    {{"by-columns",      no_argument, NULL, 0x636f6c73}},   /* */
	    {{"nrhs-by-cols",    no_argument, NULL, 0x636f6c73}},   /* undocumented alias */
	    {{"by-cols", no_argument, NULL, 0x636f6c73}},   /* undocumented alias */
')dnl
	    {{"one-nonunit-incx-incy-nrhs-per-type",	no_argument, NULL, 0x6e697270}},   /* */
	    {RSB_BENCH_PROG_OPTS},
	    {{"oski-benchmark",	no_argument, NULL, 0x42}},/* B: only long option *//* comparative benchmarking agains OSKI */
	    {{"out-lhs",		no_argument, NULL, 0x6F6C6873}},/* o */	/* should accept an output file name, optionally */
	    {{"out-rhs",		no_argument, NULL, 0x6F6F}},/* o */	/* should accept an output file name, optionally */
	    {{"override-matrix-name",	required_argument , NULL, 0x6F6D6E}},/* omn */	
	    {{"pattern-mark",	no_argument, NULL, 0x70}},/* p */
	    {{"pre-transpose",	no_argument, NULL, 0x5454}},   /* transpose the matrix before assembly  */
	    {{"read-as-binary",		required_argument, NULL, 0x62}},/* b */
	    {{"repeat-constructor",	required_argument , NULL, 0x4A4A}},
	    {{"reuse-io-arrays",	no_argument , NULL, 0x726961}}, /* ria */
	    {{"no-reuse-io-arrays",	no_argument , NULL, 0x6e726961 }}, /* nria */
	    {{"reverse-alternate-rows",	no_argument , NULL, 0x4A4A4A}},
	    {{"generate-upperband",	required_argument, NULL, 0x7575}},   /* uu */
	    {{"gen-uband",	required_argument, NULL, 0x7575}},   /* uu */
	    {{"generate-diagonal",	required_argument, NULL, 0x6464 }},   /* dd */
	    {{"gen-diag",	required_argument, NULL, 0x6464 }},   /* dd */
	    {{"implicit-diagonal",	no_argument, NULL, 0x6964}},   /* id */
	    {{"also-implicit-diagonal",	no_argument , NULL, 0x776264}},/* wbd */
	    {{"also-symmetries",	no_argument , NULL, 0x61736875}},/* ashu */
	    {{"also-short-idx",	no_argument , NULL, 0x616c6c69}},/* alli */
	    {{"also-coo-csr",	no_argument , NULL, 0x616363}},/* acc */
	    {{"also-recursive",	no_argument , NULL, 0x61726563}},/* arec */
	    {{"zig-zag",	no_argument , NULL, 0x4A4A4A}},
	    {{"subdivision-multiplier",	required_argument, NULL , 0x534D}},/* SM */
#if RSB_WANT_BOUNDED_BOXES
	    {{"bounded-box",	required_argument, NULL , 0x4242}},/* BB */
#endif /* RSB_WANT_BOUNDED_BOXES */
ifelse(mop,`mat_stats',`dnl
	    {{"max-nnz-samples",	required_argument, NULL, 0x73}},	/* s */
',`dnl
	    {{"sort",		no_argument, NULL, 0x73}},	/* s */
')dnl
	    {{"no-leaf-multivec",	no_argument, NULL , 0x6e6c6d6d}},/* nlmm */
	    {{"with-leaf-multivec",	no_argument, NULL , 0x636c6d6d}},/* wlmm */
	    {{"setenv",	required_argument, NULL, 0x7365}},/* se */
	    {{"unsetenv", required_argument, NULL, 0x757365}},/* use */
	    {{"sort-after-load",	no_argument, NULL, 0x7373}},/* ss */
	    {{"sort-filenames-list",	no_argument, NULL, 0x73666e6c}},/* sfnl */
	    {{"no-sort-filenames-list",	no_argument, NULL, 0x6e73666e}},/* nsfn */
	    {{"skip-loading-symmetric-matrices",	 no_argument, NULL, 0x736c736d}},/* slsm */
	    {{"skip-loading-unsymmetric-matrices",no_argument, NULL, 0x736c756d}},/* slum */
	    {{"skip-loading-hermitian-matrices",no_argument, NULL, 0x736c686d}},/* slhm */
	    {{"skip-loading-not-unsymmetric-matrices",no_argument, NULL, 0x736c6e75}},/* slnu */
	    {{"skip-loading-if-more-nnz-matrices",required_argument, NULL, 0x736c6d6}},/* slmn */
	    {{"skip-loading-if-less-nnz-matrices",required_argument, NULL, 0x736c6e6e}},/* slnn */
	    {{"skip-loading-if-more-filesize-kb-matrices",required_argument, NULL, 0x736c6d73}},/* slms */
#ifdef RSB_HAVE_REGEX_H
	    {{"skip-loading-if-matching-regex",required_argument, NULL, 0x736c6d72}},/* slmr */
#endif /* RSB_HAVE_REGEX_H */
	    {{"skip-loading-if-matching-substr",required_argument, NULL, 0x736c7373}},/* slss */
	    {{"times",		required_argument, NULL, 0x74}},/* t */
	    {{"transpose-as",	required_argument, NULL, 0x5040}},   /* do transpose the operation */
	    {{"transpose",	no_argument, NULL, 0x5050}},   /* do transpose the operation */
	    {{"also-transpose",	no_argument, NULL, 0x4150}},  /* N,T: do transpose the operation after no transposition */
	    {{"all-transposes",	no_argument, NULL, 0x616c6c74}},  /* N,T,C */
	    {{"type",		required_argument, NULL, 0x54}},/* T */
	    {{"types",		required_argument, NULL, 0x54}},/* T */
	    {{"update",		no_argument, NULL, 0x55}},	/* U */
	    {{"as-unsymmetric",		no_argument, NULL, 0x5555}},	/* UU: TODO: to insert such a test in as default, in order to quantify the benefit of symmetry */
	    {{"as-symmetric",		no_argument, NULL, 0x5353}},	/* SS */
	    {{"expand-symmetry",		no_argument, NULL, 0x4553}},	/* ES */
	    {{"as-hermitian",		no_argument, NULL, 0x4853}},	/* HS */
	    {{"only-lower-triangle",		no_argument, NULL, 0x4F4C54}},	/* OLT */
   	    {{"only-upper-triangle",		no_argument, NULL, 0x4F4554}},	/* OUT */
	    {{"verbose",	no_argument , NULL, 0x56}},/* V */
	    {{"less-verbose",	no_argument , NULL, 0x6c56}},/* lV */
	    {{"want-io-only",	no_argument , NULL, 0x4949}},/* --want-io-only */
	    {{"want-nonzeroes-distplot",	no_argument, NULL, 0x776E68}},/* wnh */
	    {{"want-accuracy-test",	no_argument, NULL, 0x776174}},/* wat */
	    {{"want-getdiag-bench",	no_argument , NULL, 0x774446}},/* wde */  /* FIXME: obsolete ? */
	    {{"want-getrow-bench",	no_argument , NULL, 0x777246}},/* wre */  /* FIXME: obsolete ? */
#ifdef RSB_WANT_PERFORMANCE_COUNTERS
	    {{"want-perf-counters",	no_argument , NULL, 0x707763}},/* wpc */
#endif
	    {{"want-print-per-subm-stats",	no_argument , NULL, 0x77707373}},/* wpss */
	    {{"want-only-accuracy-test",	no_argument, NULL, 0x776F6174}},/* woat */
	    {{"want-autotune",	optional_argument, NULL, 0x7772740a}},/* wrt */
	    {{"want-no-autotune",	no_argument, NULL, 0x776e7274}},/* wnrt */
	    {{"want-no-ones-fill",	no_argument, NULL, 0x776e6f66}},/* wnof */
#if RSB_HAVE_METIS
	    {{"want-metis-reordering",	no_argument, NULL, 0x776d6272 }},/* wmbr */
#endif
	    {{"want-mkl-autotune",	optional_argument, NULL, 0x776d6174}},/* wmat */
	    {{"want-mkl-one-based-indexing",	no_argument, NULL, 0x776d6f62 }},/* wmob */
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
	    {{"mkl-inspector-super-light",	no_argument , NULL, 0x6d6b6c73}},/* mkls */
	    {{"mkl-inspector-light",	no_argument , NULL, 0x6d6b6c6c}},/* mkll */
	    {{"mkl-inspector",	no_argument , NULL, 0x6d6b6c69}},/* mkli */
	    {{"mkl-no-inspector",	no_argument , NULL, 0x6d6b6c6f}},/* mklo */
')dnl
	    {{"want-unordered-coo-test",	no_argument, NULL, 0x775563}},/* */
	    {{"with-flags",	required_argument, NULL, 0x71}},/* q */
	    {{"write-as-binary",	required_argument, NULL, 0x77 }}, /* w */
	    {{"write-as-csr",	required_argument, NULL,  0x63777273 }}, /* wcsr */
	    {{"write-performance-record",	optional_argument, NULL, 0x77707266 }}, /* write performance record file  */
	    {{"performance-record-name-append",	required_argument, NULL, 0x77707261 }}, /* ...append  */
	    {{"performance-record-name-prepend",	required_argument, NULL, 0x77707270 }}, /* ...prepend  */
	    {{"write-no-performance-record",	no_argument, NULL, 0x776e7072 }}, /* write no performance record */
	    {{"discard-read-zeros",	no_argument, NULL,  0x64697a65 }}, /* dize */
	    {{"z-sorted-coo",	no_argument, NULL , 0x7A}},/* z */
	    {{0,no_argument,0,0}}	};

	const int rsb_options_count = sizeof(rsb_options)/sizeof(rsb_option_t);
	rsb_option_t options[rsb_options_count];
	rsb_nnz_idx_t nnz = 0;/* was 0 */
	int c;
	int opt_index = 0;

	rsb_coo_idx_t *IA = NULL, *JA = NULL;
	void *VA = NULL;

	int g_estimate_matrix_construction_time = 0;
dnl	int g_dump_performance_profile = 0;
	int g_all_flags = 0;
	int g_sort_only = 0;
	int repeat_construction = 1;	/* times to call the matrix constructor (the more times, the more accurate measurements) */

	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT, typecode_old = RSB_NUMERICAL_TYPE_DEFAULT;
	rsb_int ntypecodes = 0,typecodesi;
	const rsb_int maxtypes = 2*RSB_IMPLEMENTED_TYPES;
	rsb_type_t typecodes[maxtypes+1] ;
	char * typess = NULL; /* types specification string */

	rsb_blk_idx_t br = 1;
	rsb_blk_idx_t bc = 1;
	char * bcs = NULL, *brs = NULL, *cns = NULL, *mhs = NULL;
	const char *tcrprs = NULL; // thread count rpr string
	rsb_char_t fnrprs[RSB_MAX_FILENAME_LENGTH]; //  file  name  rpr string
	rsb_blk_idx_t * brv = NULL;
	rsb_blk_idx_t * bcv = NULL;
	int brl = 0;
	int bcl = 0;
	rsb_thread_t ca_[1] = {1};
	rsb_thread_t * ca = ca_;
	rsb_thread_t cn = 1, cc = ca[0];

	int times = 100;	/* the default number of times to perform mop */
	rsb_coo_idx_t nrA = 0, ncA = 0, ndA = 0;
	int filenamen = 0, filenamei = 0;
	int want_filenamelist_sort = 1;
#define RSB_RSBENCH_STATIC_FILENAMEA 1
#if RSB_RSBENCH_STATIC_FILENAMEA
	rsb_char_t *filenamea[RSB_RSBENCH_MAX_MTXFILES];
#else
	rsb_char_t **filenamea = NULL;
#endif
	const rsb_char_t *filename = NULL;
	const rsb_char_t *filename_old = NULL;
	const rsb_char_t *usfnbuf = NULL;
	rsb_char_t *fprfn = NULL, *cprfn = NULL, *apprfn = NULL, *ppprfn = NULL; /* final/checkpoint      performance file name , append/prepend */
	rsb_char_t fprfnb[RSB_MAX_FILENAME_LENGTH], cprfnb[RSB_MAX_FILENAME_LENGTH];/* final/checkpoint      performance file name buffers */
	rsb_char_t mtrfpf[RSB_MAX_FILENAME_LENGTH];/* matrix tuning rendering file prefix */
	rsb_char_t fnbuf[RSB_MAX_FILENAME_LENGTH] = {RSB_NUL, /* ... */ }; /* matrix file name buffer */
	rsb_char_t *fnbufp[1]={&(fnbuf[0])};
	rsb_char_t *dump_graph_file=NULL;
	rsb_flags_t flags_o = RSB_FLAG_NOFLAGS|RSB_FLAG_OWN_PARTITIONING_ARRAYS;
	rsb_flags_t flags_s = RSB_FLAG_NOFLAGS; /* sorting flags */
/*	RSB_DO_FLAG_ADD(flags_o,RSB_FLAG_DISCARD_ZEROS)	;	*/ /* FIXME : EXPERIMENTAL (watch nnz count on a multi blocking run ...) */
	rsb_flags_t flagsa[128] = RSB_M4_ZEROS_ARRAY(128);
	rsb_flags_t r_flags = RSB_FLAG_NOFLAGS; /* recycling flags */
	int fn = 1, fi = 0;/* for flags */
	int tn = 1, ti = 0;/* for transposition */
	int g_debug = 0;
	int be_verbose = 0;
	int pattern_only = 0;
	int dumpout = 0;
	int dumpout_internals = 0, merge_experimental = 0, split_experimental = 0;
	int just_enter_tuning = 1;
	rsb_char_t * csr_w_filename = NULL;
	rsb_char_t * b_w_filename = NULL;
	rsb_char_t * b_r_filename = NULL;
	int dumpvec = rsb_dumpvec_no;
	struct rsb_mtx_t * mtxAp = NULL;
dnl
	rsb_err_t errval = RSB_ERR_NO_ERROR;
dnl
	rsb_blk_idx_t rua[] = RSB_ROWS_UNROLL_ARRAY;
	rsb_blk_idx_t cua[] = RSB_COLUMNS_UNROLL_ARRAY;
dnl
	int guess_blocking_test = 0;		/* guess test stuff */
	rsb_int want_column_expand = 0;
dnl
ifelse(mop,`mat_stats',`',`dnl
	rsb_perf_t bperf=0,wperf=0,cperf=0;			/* guess test stuff */
	rsb_fillin_t egfillin=0,ebfillin=0,bfillin=0,maxfillin=0;	/* guess test stuff */
	rsb_blk_idx_t bri=0,bci=0;		/* guess test stuff */
	rsb_perf_t omta = RSB_REAL_ZERO; /* op memory traffic amount */
	rsb_fillin_t fillin = RSB_REAL_ZERO;
	rsb_perf_t raw_Mflops = RSB_REAL_ZERO,true_Mflops = RSB_REAL_ZERO, true_gem_Mflops = RSB_REAL_ZERO;
	rsb_char_t buf[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];/* Flawfinder: ignore */
	rsb_fillin_t efillin = RSB_REAL_ZERO;
	rsb_perf_t eperf = RSB_REAL_ZERO;
')dnl

	rsb_bool_t should_recycle_matrix = RSB_BOOL_FALSE; /* reuse the matrix across measurements */
	rsb_bool_t should_recycle_io = RSB_BOOL_TRUE;/* reuse the input arrays */
	rsb_bool_t g_allow_any_tr_comb = RSB_BOOL_FALSE; /* allow any transposition combination */
	rsb_bool_t g_fill_va_with_ones = RSB_BOOL_TRUE;
	
ifelse(mop,`mat_stats',`dnl
	int g_estimate_fillin = 0;
	int want_percentage = 0;
	double until_confidence = 0;

	rsb_nnz_idx_t  max_nnzs = 0;
	rsb_nnz_idx_t nnzn = 10;
	rsb_nnz_idx_t * nnzs = NULL;
	size_t * element_count = NULL;
	size_t * block_count = NULL;
	//rsb_nnz_idx_t i = 0;
	struct rsb_mtx_partitioning_info_t pinfo;
dnl	struct rsb_mop_performance_info_t mpi;
')dnl
	rsb_trans_t transAo = RSB_DEFAULT_TRANSPOSITION;
	rsb_trans_t transA = RSB_DEFAULT_TRANSPOSITION;
	rsb_nnz_idx_t should_generate_dense = 0;
	rsb_nnz_idx_t should_generate_dense_nc = 0;
	rsb_nnz_idx_t should_generate_lband = -1, should_generate_uband = -1;
	rsb_nnz_idx_t want_generated_spacing = 0;
	rsb_bool_t want_only_star_scan = RSB_BOOL_FALSE;
dnl
	rsb_blk_idx_t nrhs = 1, nrhsn = 1, nrhsi = 1, nrhsl = 1;
	rsb_blk_idx_t asfii = 0, asfin = 1;
	rsb_blk_idx_t asiii = 0, asiin = 1;
	rsb_blk_idx_t accii = 0, accin = 1;
	rsb_blk_idx_t areci = 0, arecn = 1;
	rsb_int_t diagin = 1, diagii = 0;
	const char*nrhss = NULL;
	rsb_blk_idx_t *nrhsa = NULL;
dnl
	const size_t outnri = 0, rhsnri = 0; /* Could be ndA for order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER and nrhs otherwise; this way is auto. */
	rsb_nnz_idx_t n_dumpres = 0;
	rsb_nnz_idx_t n_dumprhs = 0;
	rsb_bool_t ignore_failed_fio = RSB_BOOL_TRUE; /* FIXME 20140912 experimental */
	rsb_bool_t want_convert = RSB_BOOL_FALSE;
	rsb_bool_t want_update = RSB_BOOL_FALSE;
	rsb_int_t want_impatiently_soon_pre_results = 0; /* FIXME: temporary */
	rsb_bool_t want_inner_flush = RSB_BOOL_FALSE;
	rsb_bool_t want_outer_flush = RSB_BOOL_TRUE;
	rsb_bool_t want_ancillary_execs = RSB_BOOL_FALSE;
	rsb_time_t st = RSB_TIME_ZERO;
	rsb_time_t totiot = RSB_TIME_ZERO; /* total I/O time */
	rsb_time_t totatt = RSB_TIME_ZERO; /* total ancillary tests time */ /* FIXME: is this complete ? */
	rsb_time_t totct = RSB_TIME_ZERO; /* total conversions time */ /* FIXME: is this complete ? */
	rsb_time_t tottt = RSB_TIME_ZERO; /* total tuning time */
	rsb_time_t totht = RSB_TIME_ZERO, ht = RSB_TIME_ZERO; /* total checks time */ /* FIXME: is this complete ? */
	rsb_time_t maxtprt = RSB_TIME_ZERO; /* max total program run time */
	const rsb_time_t totprt = - rsb_time(); /* total program run time */
	rsb_bool_t want_as_unsymmetric = RSB_BOOL_FALSE;
	rsb_bool_t want_as_symmetric = RSB_BOOL_FALSE;
	rsb_bool_t want_expand_symmetry = RSB_BOOL_FALSE;
	rsb_bool_t want_as_hermitian = RSB_BOOL_FALSE;
	rsb_bool_t want_only_lowtri = RSB_BOOL_FALSE;
	rsb_bool_t want_only_upptri = RSB_BOOL_FALSE;
	rsb_bool_t want_sort_after_load = RSB_BOOL_FALSE;
	rsb_bool_t want_slsm = RSB_BOOL_FALSE, want_slum = RSB_BOOL_FALSE, want_slnu = RSB_BOOL_FALSE, want_slhm = RSB_BOOL_FALSE;
	rsb_nnz_idx_t want_slmn = 0,  want_slnn = 0,  want_slms = 0;
#ifdef RSB_HAVE_REGEX_H
	const rsb_char_t * want_slmr = NULL;
#endif /* RSB_HAVE_REGEX_H */
	const rsb_char_t * want_slss = NULL;
	rsb_bool_t do_perform_ilu = RSB_BOOL_FALSE;
	rsb_bool_t do_perform_ddc = RSB_BOOL_FALSE;
	rsb_bool_t want_in_place_assembly = RSB_BOOL_FALSE;
	rsb_bool_t want_accuracy_test = 0;	/* FIXME-EXPERIMENTAL */
	rsb_bool_t want_nonzeroes_distplot = 0;	/* FIXME-EXPERIMENTAL */
	rsb_bool_t want_getdiag_bench = 0;	/* FIXME-EXPERIMENTAL */
	rsb_bool_t want_getrow_bench = 0;	/* FIXME-EXPERIMENTAL */
dnl
	rsb_coo_idx_t mib = 0; /* MKL index base (FIXME: declared here and not within RSB_WANT_MKL because CSR copy made even with no MKL) */
dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
	rsb_bool_t want_mkl_bench = RSB_BOOL_FALSE;
	rsb_bool_t want_mkl_bench_csr = RSB_BOOL_TRUE;
	rsb_bool_t want_mkl_bench_gem = RSB_BOOL_TRUE;
	rsb_bool_t want_mkl_bench_coo = RSB_BOOL_FALSE;
#if RSB_WANT_MKL_INSPECTOR
	rsb_bool_t want_mkl_inspector = RSB_BOOL_TRUE;
	rsb_int_t mkl_ie_hintlvl = 2;
#endif /* RSB_WANT_MKL_INSPECTOR */
#endif /* RSB_WANT_MKL */
')dnl
dnl
#if RSB_WANT_ARMPL
	rsb_bool_t want_armpl_bench = RSB_BOOL_FALSE;
#endif /* RSB_WANT_ARMPL */
dnl
	rsb_time_t totmt = RSB_TIME_ZERO; /* total mkl/competitors (tuning) time */
	rsb_bool_t want_perf_dump = RSB_BOOL_FALSE;
	void*rspr = NULL; /* rsb sampled performance record structure pointer */

ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	rsb_aligned_t errnorm[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t alpha[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t beta[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	void * alphap = &(alpha[0]);
	void * betap = &(beta[0]);
	rsb_blk_idx_t * alphaip = NULL; // here, int
	rsb_blk_idx_t * betaip = NULL; // here, int
	rsb_int alphan = 1, betan = 1;
	rsb_int alphai = 1, betai = 1;
	const char * betass = NULL; /* beta specification string */
	const char * alphass = NULL; /* alpha specification string */
')dnl
dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
dnl
	rsb_coo_idx_t incX = 1, incY = 1;
	rsb_blk_idx_t incXn = 1, incXi = 1;
	rsb_blk_idx_t incYn = 1, incYi = 1;
	rsb_blk_idx_t *incXa = NULL, *incYa = NULL;
	const char * dis = "1,2"; /* default inc string */
	const char * incYss = NULL; /* incY specification string */
	const char * incXss = NULL; /* incX specification string */
',`dnl
	const rsb_coo_idx_t incX = 1, incY = 1;
	const rsb_blk_idx_t incXn = 1, incYn = 1;
	const rsb_blk_idx_t incXi = 0, incYi = 0;
	const rsb_blk_idx_t *const incXa = (const rsb_blk_idx_t *const)&incX, *const incYa = (const rsb_blk_idx_t *const)&incY;
')dnl
dnl
dnl	rsb_coo_idx_t ldX = 0, ldY = 0;
	rsb_int_t want_verbose = 0;
	rsb_int_t want_verbose_tuning = 0;
	rsb_bool_t want_mbw = RSB_BOOL_FALSE;
	rsb_bool_t want_transpose = RSB_BOOL_FALSE;
	#if 1
	const int max_io = 10;
	struct rsb_initopts io={NULL,NULL,0,RSB_IO_SPECIFIER_SET},*iop=&io;
dnl	rsb_int_t preferred_sorting_method=1;
	rsb_int_t should_use_cb_method = 0;
	rsb_real_t subdivision_multiplier = 0.0;
#if RSB_WANT_BOUNDED_BOXES
	rsb_int_t want_bounded_box=1;
#endif /* RSB_WANT_BOUNDED_BOXES */
	rsb_int_t want_no_leaf_spmm=0;
	void * io_values[max_io];
	enum rsb_opt_t io_keys[max_io];
	#else /* 1 */
	struct rsb_initopts *iop = RSB_NULL_INIT_OPTIONS;
	#endif /* 1 */
	rsb_int_t should_use_alternate_sort = 0;
	rsb_bool_t reverse_odd_rows = RSB_BOOL_FALSE;
	rsb_bool_t zsort_for_coo = RSB_BOOL_FALSE;
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	rsb_bool_t want_unordered_coo_bench = RSB_BOOL_FALSE;
')dnl
dnl
#ifdef RSB_WANT_OSKI_BENCHMARKING 
	/* FIXME : unfinished */
	rsb_time_t oski_t = RSB_TIME_ZERO,oski_m_t = RSB_TIME_ZERO,oski_a_t = RSB_TIME_ZERO,oski_t_t = RSB_TIME_ZERO;
	oski_idx_t * Aptr=NULL;
	oski_idx_t * Aind=NULL;
	oski_value_t * Aval=NULL;
	oski_matrix_t A_tunable;
        oski_vecview_t x_view;
        oski_vecview_t y_view;
	void * Oval = NULL;
	rsb_coo_idx_t *OIA=NULL,*OJA=NULL;
        rsb_char_t oxform[256];
        double oalpha = 1, obeta = 0;
	rsb_bool_t want_oski_bench=0;
	#ifdef RSB_HAVE_SETENV
	rsb__setenv("OSKI_LUA_PATH",OSKI_LUA_PATH,0/* if 0, will not override. if 1, it would. */);
	#endif /* RSB_HAVE_SETENV */
#endif /* RSB_WANT_OSKI_BENCHMARKING */
	const rsb_time_t tinf = rsb__timer_granularity(); /* After rsb_lib_init one may use RSB_CACHED_TIMER_GRANULARITY; */
dnl
	rsb_aligned_t pone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
dnl
	rsb_bool_t want_likwid = RSB_BOOL_FALSE;
dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),1,`dnl
	rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
')dnl
dnl
	rsb_time_t want_autotuner = RSB_NEGATED_EXAGGERATED_TUNER_TIMES, want_mkl_autotuner = RSB_NEGATED_EXAGGERATED_TUNER_TIMES;
	rsb_bool_t want_io_only = RSB_BOOL_FALSE;
	rsb_int wat = 1;	/* want autotuning threads choice */
	rsb_int wai = 1;	/* want autotuning rounds */
	char wav = 0x56;	/* want autotuning verbose */
	int wavf = RSB_AUT0_TUNING_VERBOSE;
#ifdef RSB_WANT_PERFORMANCE_COUNTERS
	int want_perf_counters = 0;
#endif
	rsb_bool_t want_print_per_subm_stats = RSB_BOOL_FALSE;
#if RSB_HAVE_METIS
	rsb_bool_t want_wmbr = RSB_BOOL_FALSE;
#endif
	rsb_bool_t want_recursive = RSB_BOOL_TRUE;
	struct rsb_mbw_et_t mbet, *mbetp = NULL;
	rsb_real_t mtx_sample_rate = 100.0;

	{
		int roi;
		for(roi=0;roi<rsb_options_count;++roi)
			options[roi] = rsb_options[roi].ro;
	}

	RSB_DEBUG_ASSERT( fnbuf[0]==RSB_NUL );
	io.keys = io_keys;
	io.values = io_values;
	io.n_pairs = 0;

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))!=RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err,"Error while initializing the library.");
	}

	optind = 0; /* FIXME: NEW 20160708 */
    	for (;;)
	{
		c = rsb__getopt_long(argc,argv,RSB_SAMPLE_PROGRAM_OPTIONS_GET_FLAGS"b:w:BGht:f:r:c:vpn:MNS:Bk:KU" /* Flawfinder: ignore */
ifelse(mop,`mat_stats',`dnl
		"s:e"
',`dnl
		/* s is in anyway, with RSB_SAMPLE_PROGRAM_OPTIONS_GET_FLAGS */
')dnl
		"o:O:"
		, options, &opt_index);
		if (c == -1)break;

		RSB_DO_FLAG_ADD(flags_o,rsb__sample_program_options_get_flags(c,optarg));
		RSB_DO_FLAG_ADD(flags_s,flags_o & RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING);

		switch (c)
		{
			case 0x62:	/* b */
			b_r_filename = optarg;
			break;
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
			case  0xb1bb0: // incX
			incXss=optarg;
			break;
			case  0xb1bb1: // incY
			incYss=optarg;
			break;
			case  0x696e632a: // inc*
			incXss=optarg;
			incYss=optarg;
			break;
')dnl
			case  0x6970:
				RSBENCH_STDOUT("# WARNING: in place assembly is an UNFINISHED, EXPERIMENTAL feature\n");
				want_in_place_assembly = RSB_BOOL_TRUE;
			break;
			case 0x6c:
			case 0x6c64: /* lower-dense */
dnl			case 0xbabb1:
			{
				should_generate_dense = - rsb__util_atoi_km10(optarg); // FIXME ! PROBLEMS
			}
			break;
			case 0x6c696b77:
#if RSB_WITH_LIKWID
				want_likwid = RSB_BOOL_TRUE;
				dnl RSBENCH_STDOUT("Usage of the LIKWID API requested.\n");
#else /* RSB_WITH_LIKWID */
				dnl RSBENCH_STDOUT("Sorry, LIKWID has not been configured in !\n");
#endif /* RSB_WITH_LIKWID */
			break;
			case 0x6c6c:
			{
				should_generate_lband = rsb__util_atoi_km10(optarg); // FIXME ! PROBLEMS
				if(should_generate_uband==-1)should_generate_uband=0;
			}
			break;
			case 0x7575:
			{
				should_generate_uband = rsb__util_atoi_km10(optarg); // FIXME ! PROBLEMS
				if(should_generate_lband==-1)should_generate_lband=0;
			}
			break;
			case 0xbabb2:
			{
				want_generated_spacing = rsb__util_atoi_km10(optarg);
			}
			break;
			case 0x6e697270:
dnl
			want_only_star_scan = RSB_BOOL_TRUE;
			break;
dnl
			case 0x6964:
				RSB_DO_FLAG_ADD(flags_o,RSB_FLAG_UNIT_DIAG_IMPLICIT);
			break;
			case 0x6464: /* gen-diag */
				should_generate_uband = 0;
				should_generate_lband = 0;
			case 0x64: /* dense */
dnl			case 0xbabb0:
			{
dnl				/* should_generate_dense = rsb__util_atoi_km10(optarg); */  // FIXME ! PROBLEMS
dnl				int sargs = sscanf(optarg,"%dx%d",&should_generate_dense,&should_generate_dense_nc);
				const char*p= optarg, *xp="x";
				should_generate_dense = rsb__util_atoi_km10(optarg);
				while(*p && *p != *xp)
					++p;
				if(*p == *xp)
					should_generate_dense_nc = rsb__util_atoi_km10(p+1);
				if( should_generate_dense_nc == 0)
					should_generate_dense_nc = should_generate_dense;
dnl				/* RSBENCH_STDOUT("# Requested generation of a %d by %d matrix\n",should_generate_dense,should_generate_dense_nc); */
			}
			break;
			/* note that specifying multiple times -r (or -c) will overwrite old setting */
			case 0x72:/* r */
			brs=optarg;
			break;
			case 0x63: /* c */
			bcs=optarg;
			break;
			case 0x42: /* oski : B */
#ifdef RSB_WANT_OSKI_BENCHMARKING 
			want_oski_bench = RSB_BOOL_TRUE;
#else /* RSB_WANT_OSKI_BENCHMARKING */
			RSB_ERROR("Sorry, OSKI comparative benchmarking was opted out at compile time\n");
			goto err;
#endif /* RSB_WANT_OSKI_BENCHMARKING */
			break;
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
			case 0x4C: /* MKL : L */
#if RSB_WANT_MKL
			want_mkl_bench = RSB_BOOL_TRUE;
#else /* RSB_WANT_MKL */
			RSB_ERROR("Sorry, MKL comparative benchmarking was opted out at compile time\n");
			goto err;
#endif /* RSB_WANT_MKL */
			break;
			case 0x6d6b6c73: /* mkls */
#if RSB_WANT_MKL_INSPECTOR
			mkl_ie_hintlvl = 0;
#endif /* RSB_WANT_MKL_INSPECTOR */
			goto wi;
			break;
			case 0x6d6b6c6c: /* mkll */
#if RSB_WANT_MKL_INSPECTOR
			mkl_ie_hintlvl = 1;
#endif /* RSB_WANT_MKL_INSPECTOR */
			goto wi;
			break;
			case 0x6d6b6c69: /* mkli */
wi:
#if RSB_WANT_MKL_INSPECTOR
			want_mkl_inspector = RSB_BOOL_TRUE;
#else /* RSB_WANT_MKL_INSPECTOR */
			RSB_ERROR("Sorry, no MKL Inspector Interface in --- aborting!\n");
			goto err;
#endif /* RSB_WANT_MKL_INSPECTOR */
			break;
#if RSB_WANT_MKL_INSPECTOR
			case 0x6d6b6c6f: /* mklo */
				want_mkl_inspector = RSB_BOOL_FALSE;
			break;
#endif /* RSB_WANT_MKL_INSPECTOR */
')dnl
			case 0x61617463:
			g_allow_any_tr_comb = RSB_BOOL_TRUE;
			break;
			case 0x51: /* Q (do not ask me why) */
			g_all_flags = 1;
			break;
			case 0x616c6274: /* albt */
				RSBENCH_STDOUT("# EXPERIMENTAL --all-blas-types option implies types: " RSB_BLAS_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS "\n");
				typess="B"; // RSB_BLAS_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS
			break;
			case 0x616c626f: /* albo */
			RSBENCH_STDOUT("# EXPERIMENTAL --all-blas-opts option implies: --types : --inc : --alpha : --beta : :\n");
			RSBENCH_STDOUT("# EXPERIMENTAL --types :\n"); // 0x54
			typess=":";
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
			RSBENCH_STDOUT("# EXPERIMENTAL --inc :\n"); // 0x696e632a
			incXss=":";
			incYss=":";
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
			RSBENCH_STDOUT("# EXPERIMENTAL --alpha :\n"); // 0x414C
			alphass=":";
			RSBENCH_STDOUT("# EXPERIMENTAL --beta :\n"); // 0x4246
			betass=":";
')dnl
dnl
			break;
			case 0x616c666f: /* alfo */
			/* FIXME: need formal bond with these options */
			RSBENCH_STDOUT("# EXPERIMENTAL --all-formats option implies: --also-implicit-diagonal --also-symmetries --also-short-idx --also-coo-csr --also-recursive:\n");
			RSBENCH_STDOUT("# EXPERIMENTAL --also-implicit-diagonal\n"); // 0x6964: /* id */
			diagin=2;
			RSBENCH_STDOUT("# EXPERIMENTAL --also-symmetries\n"); // 0x61736875: /* ashu */
			asfin=3;
			RSBENCH_STDOUT("# EXPERIMENTAL --also-short-idx\n"); // 0x616c6c69: /* alli */
			asiin=2;
			RSBENCH_STDOUT("# EXPERIMENTAL --also-coo-csr\n"); // 0x616363: /* acc */
			accin=2;
			RSBENCH_STDOUT("# EXPERIMENTAL --also-recursive\n"); // 0x61726563: /* arec */
			arecn=2;
			break;
dnl			break;
dnl			case 0x4F: /* O */
dnl			g_dump_performance_profile=1;
			break;
			case 0x44044: /* D */
			dumpout = 1;
			break;
			case 0x5040: /*  */
				transAo = rsb__do_transposition_from_char(*optarg);	/* */
				tn = 1;
				RSBENCH_STDOUT("# Requested user-specified transposition: %c\n",(char)transAo );
			break;
			case 0x4150:
				RSBENCH_STDOUT("# Requested no transposition and transposition.\n");
				tn = 2;
			break;
			case 0x616c6c74:
				RSBENCH_STDOUT("# Requested all transposes: no transposition, transposition and conjugate transposition.\n");
				tn = 3;
			break;
			case 0x5050: /*  */
				tn = 1;
				RSBENCH_STDOUT("# Requested transposition.\n");
				transAo = rsb__do_transpose_transposition(transAo);
			break;
			case 0x5051: /*  */
				RSBENCH_STDOUT("# Requested no transposition.\n");
				tn = 1;
				transAo = RSB_TRANSPOSITION_N;
			break;
			case 0x6e726873: /*  */
			nrhss = optarg;
			if(RSB_SOME_ERROR(rsb__util_get_bx_array(nrhss,&nrhsn,&nrhsa)) || nrhsn<1)
				{RSB_ERROR(RSB_ERRM_ES);goto err;}
			break;
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),1,`dnl
			case 0x726f7773: /* --nrhs-by-rows --by-rows */
				order = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
			break;
			case 0x636f6c73: /* --nrhs-by-columns --by-columns */
				order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
			break;
')dnl
			case 0x5454: /*  */
			want_transpose = !want_transpose;
			break;
			case 0x44047: /* DG */
			dump_graph_file = optarg;
			break;
			case 0x49049: /* I */
			dumpout_internals = 1;
			break;
			case 0x6d656578: /* meex */
				merge_experimental = rsb__util_atoi(optarg);
				RSB_ASSIGN_IF_ZERO(merge_experimental,RSB_CONST_DEF_MS_AT_AUTO_STEPS);
			break;
			case 0x73706578: /* spex */
				split_experimental = rsb__util_atoi(optarg);
				RSB_ASSIGN_IF_ZERO(split_experimental,RSB_CONST_DEF_MS_AT_AUTO_STEPS);
			break;
			case 0x6d736578: /* msex */
				merge_experimental = split_experimental = rsb__util_atoi(optarg);
				RSB_ASSIGN_IF_ZERO(merge_experimental,RSB_CONST_DEF_MS_AT_AUTO_STEPS);
				RSB_ASSIGN_IF_ZERO(split_experimental,RSB_CONST_DEF_MS_AT_AUTO_STEPS);
			break;
			case 0x4444 : /* DD */
			do_perform_ddc = RSB_BOOL_TRUE;
			break;
			case 0x444444 : /* DDD */
			n_dumprhs = n_dumpres = rsb__util_atoi_km10(optarg);
			break;
			case 0x6563686f: /* echo */
			{
				rsb_int argi=0;
				if(argc>0) printf("#args: %s",argv[0]);
				for(argi=1;argi<argc;++argi)
					printf(" %s",argv[argi]);
				printf("\n");
			}
			break;
			case 0x494B55 : /* ILU */
			do_perform_ilu = RSB_BOOL_TRUE;
			break;
			case 0x696d7061: /* */
			want_impatiently_soon_pre_results = 1;
			break;
			case 0x4343: /* */
			want_inner_flush = RSB_BOOL_TRUE;
			break;
			case 0x6364: /* */
			{
				if( chdir(optarg) ) // unistd.h
					RSB_WARN("change dir to %s failed!\n", optarg);
				else
					RSBENCH_STDOUT("# chdir to %s succeeded\n", optarg);
			}
			break;
			case 0x434E: /* */
			want_inner_flush = RSB_BOOL_FALSE;
			break;
			case 0x434343: /*  */
			want_outer_flush = RSB_BOOL_TRUE;
			break;
			case 0x43434E: /*  */
			want_outer_flush = RSB_BOOL_FALSE;
			break;
			case 0x776e720a: /*  */
			want_recursive = RSB_BOOL_FALSE;
			break;
			case 0x776e6d62: /* wnmb */
			want_mbw = RSB_BOOL_FALSE;
			break;
			case 0x776d62: /* wmb */
			want_mbw = RSB_BOOL_TRUE;
			break;
			case 0x4D: /* M */
			g_estimate_matrix_construction_time=1;
			break;
ifelse(mop,`mat_stats',`dnl
			case 0x65: /* e */
			g_estimate_fillin=1;
			break;
')dnl
			case 0x7A:
			zsort_for_coo = RSB_BOOL_TRUE;
			RSBENCH_STDOUT("# WARNING: the now active Z sort feature will only apply to COO submatrices\n");
			break;
			case 0x726961:
			RSBENCH_STDOUT("# setting the reuse I/O arrays option in e.g.: type transitions\n");
			should_recycle_io = RSB_BOOL_TRUE;
			break;
			case 0x6e726961:
			RSBENCH_STDOUT("# unsetting the reuse I/O arrays option in e.g.: type transitions\n");
			should_recycle_io = RSB_BOOL_FALSE;
			break;
			case 0x4A4A4A:
			reverse_odd_rows = RSB_BOOL_TRUE;
			RSBENCH_STDOUT("# WARNING: the row reversal feature only applies to CSR submatrices, and on indices only\n");
			break;
			case 0x6F6D6E:
			usfnbuf = optarg;
			break;
			case 0x4A4A:
			repeat_construction = rsb__util_atoi_km10(optarg);
			if(repeat_construction<1)
			{
				RSB_ERROR("Constructor repetition times should be a positive number!\n");goto err;
			}
			break;
			case 0x4342: /* CB */
			should_use_cb_method = rsb__util_atoi(optarg);
			break;
			case 0x4153: /* AS */
			should_use_alternate_sort = rsb__util_atoi(optarg);
			break;
			case 0x534D: /* SM */
			subdivision_multiplier = rsb__util_atof(optarg);
			break;
#if RSB_WANT_BOUNDED_BOXES
			case 0x4242: /* BB */
			want_bounded_box = rsb__util_atoi(optarg);
			break;
#endif /* RSB_WANT_BOUNDED_BOXES */
			case 0x6e6c6d6d: /* nlmm */
			want_no_leaf_spmm = /*rsb__util_atoi(optarg)*/ -1;
			break;
			case 0x636c6d6d: /* wlmm */
#if RSB_ENABLE_INNER_NRHS_SPMV
			want_no_leaf_spmm = 0;
#else
			RSB_ERROR("Cannot activate the RSB_IO_WANT_LEAF_LEVEL_MULTIVEC option because RSB_ENABLE_INNER_NRHS_SPMV is opted out!\n");goto err;
#endif
			break;
			case 0x4D4D: /* MM */
			mhs = optarg;
			break;
			case 0x6d617275:
			maxtprt = rsb__util_atof(optarg);
			maxtprt = RSB_MAX( RSB_TIME_ZERO, maxtprt  );
			break;
			case 0x6F6C6873: /* o */
			dumpvec = rsb_dumpvec_res;
			break;
			case 0x6F6F: /* o */
			dumpvec = rsb_dumpvec_rhs;
			break;
			case 0x70: /* p */
			pattern_only = 1;
			break;
			case 0x4E: /* N */
			g_sort_only = 1;
			break;
ifelse(mop,`mat_stats',`dnl
			case 0x73: /* s */
			/* FIXME : BROKEN! */
			max_nnzs = rsb__util_atoi_km10(optarg); // or revive rsb__util_atonnz
			if(*optarg && optarg[rsb__util_strlen(optarg)-1]==0x25)want_percentage=1;/* 0x25 == % */
			break;
',`dnl	
			/* handled by rsb__sample_program_options_get_flags() */
			case 0x73: /* s */
				RSB_DEPRECATED("use of the sort flag");
				flags_o = flags_o;
			break;
')dnl
ifelse(mop,`mat_stats',`dnl
			case 0x53: /* S */
			nnzn = rsb__util_atoi_km10(optarg); // or revive rsb__util_atonnz
			if(nnzn<1){RSB_ERROR(RSB_ERRM_ES);goto err;}
			break;
')dnl
#ifdef RSB_HAVE_SETENV
			case 0x7365: /* se */
			rsb__setenv(optarg);
			break;
#elif RSB_HAVE_PUTENV
			case 0x7365: /* se */
			putenv(optarg);
			break;
#endif /* RSB_HAVE_SETENV */
#ifdef RSB_HAVE_UNSETENV
			case 0x757365: /* use */
			if(optarg)
			{
				RSBENCH_STDOUT("# Calling unsetenv() with argument %s\n",optarg);
				unsetenv(optarg);
			}
			break;
#endif /* RSB_HAVE_UNSETENV */
#ifndef RSB_HAVE_SETENV
#ifndef RSB_HAVE_PUTENV
			case 0x7365: /* se */
			RSBENCH_STDERR("neither of setenv or putenv found on your machine ! Skipping --setenv !");
			break;
#endif /* RSB_HAVE_PUTENV */
#endif /* RSB_HAVE_SETENV */
#ifndef RSB_HAVE_UNSETENV
			case 0x757365: /* use */
			RSBENCH_STDERR("unsetenv not found on your machine ! Skipping ---unsetenv !");
			break;
#endif /* RSB_HAVE_UNSETENV */
			case 0x7373: /* ss */
			want_sort_after_load = RSB_BOOL_TRUE;
			break;
			case 0x73666e6c: /* sfnl */
			want_filenamelist_sort = 1;
			break;
			case 0x6e73666e: /* nsfn */
			want_filenamelist_sort = 0;
			break;
			case 0x736c736d: /* slsm */
			want_slsm = RSB_BOOL_TRUE;
			break;
			case 0x736c756d: /* slum */
			want_slum = RSB_BOOL_TRUE;
			break;
			case 0x736c686d: /* slhm */
			want_slhm = RSB_BOOL_TRUE;
			break;
			case 0x736c6e75: /* slnu */
			want_slnu = RSB_BOOL_TRUE;
			break;
			case 0x736c6d6: /* slmn */
			want_slmn = rsb__util_atoi_km10(optarg);
			break;
			case 0x736c6e6e: /* slnn */
			want_slnn = rsb__util_atoi_km10(optarg);
			break;
			case 0x736c6d73: /* slms */
			want_slms = rsb__util_atoi_km2(optarg);
			break;
#ifdef RSB_HAVE_REGEX_H
			case 0x736c6d72: /* slmr */
			want_slmr = (optarg);
			break;
#endif /* RSB_HAVE_REGEX_H */
			case 0x736c7373: /* slss */
			want_slss = (optarg);
			break;
			case 0x74: /* t */
			times = rsb__util_atoi(optarg);
			break;
			case 0x47: /* G */
			guess_blocking_test = 1;
			break;
			case 0x54: /* T */
			typess=optarg;
			break;
			case 0x56: /* V */
			want_verbose++;
			want_verbose_tuning ++;
			break;
			case 0x6c56: /* lV */
			if(want_verbose_tuning >  0)
				want_verbose_tuning --;
			if(want_verbose_tuning == 0)
				want_verbose--;
			break;
			case 0x4949: /* II */
			want_io_only = RSB_BOOL_TRUE;
			break;
			case 0x46: /* F */
			/* F handled by rsb__sample_program_options_get_flags() */
			break;
			case 0x66: /* f */
			filename = optarg;
#if RSB_RSBENCH_STATIC_FILENAMEA
/* #define RSB_RSBENCH_ADDF(FILENAME)	if(filenamen<RSB_RSBENCH_MAX_MTXFILES)filenamea[filenamen++] = (FILENAME); else {errval = RSB_ERR_INTERNAL_ERROR;RSB_ERROR("Please increase RSB_RSBENCH_MAX_MTXFILES (%d) and recompile !!\n",RSB_RSBENCH_MAX_MTXFILES);goto err;} */
#define RSB_RSBENCH_ADDF(FILENAME, FAF)	if ( RSB_ERR_NO_ERROR != ( errval = rsb__adddir(&filenamea[0], &filenamen, (FILENAME), FAF) ) ){if(errval == RSB_ERR_LIMITS ){RSB_ERROR("Please increase RSB_RSBENCH_MAX_MTXFILES (%d) and recompile !!\n",RSB_RSBENCH_MAX_MTXFILES);}errval = RSB_ERR_INTERNAL_ERROR;goto err;}
#else
 /* FIXME: for some reason, this seems to break e.g.  ./rsbench -oa -Ob --nrhs 1,2 -f pd.mtx -f A.mtx.
    Of course this is wrong also w.r.t. rsb_calloc/rsb_lib_init, but that is not a problem.
    Using calloc / realloc does not solve the problem.  */
#define RSB_RSBENCH_ADDF(FILENAME,RSB_FAF_DEFAULTS)		if(filenamen==0) \
				filenamea = rsb__calloc(sizeof(filenamea)*(filenamen+1)); \
			else \
				filenamea = rsb__do_realloc(filenamea, sizeof(filenamea)*(filenamen+1), sizeof(filenamea)); \
			filenamea[filenamen++] = (FILENAME);
#endif
			RSB_RSBENCH_ADDF(filename,RSB_FAF_DEFAULTS) /* FIXME */
			break;
			case 0x6d747873: /* mtxs **/
				mtx_sample_rate = rsb__util_atof(optarg);
				RSBENCH_STDOUT("# Using only %lg %% of the matrix nonzeroes\n",mtx_sample_rate);
				if( mtx_sample_rate > 100 || mtx_sample_rate < 1 )
					RSB_PERR_GOTO(err,"--matrix-sample-pcnt: want arg between 1 and 100!");
			break;
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
			case 0x414C: /* AL */
				alphass=optarg;
			break;
			case 0x4246: /* BE */
				betass=optarg;
			break;
')dnl
			case 0x626e6368: /* bnch */
			RSBENCH_STDOUT("# --bench option implies -qH -R --write-performance-record --want-mkl-autotune --mkl-benchmark --types : --split-experimental %d --merge-experimental %d --also-transpose --sort-filenames-list --want-memory-benchmark\n",RSB_CONST_DEF_MS_AT_AUTO_STEPS,RSB_CONST_DEF_MS_AT_AUTO_STEPS);
			want_filenamelist_sort = 1; /* 0x73666e6c: sfnl */
			want_perf_dump = RSB_BOOL_TRUE; // 0x77707266
			want_mkl_autotuner = 1.0;
			want_mkl_autotuner = RSB_MAX(1.0,want_mkl_autotuner); // 0x776d6174
			RSB_DO_FLAG_ADD(flags_o,RSB_FLAG_QUAD_PARTITIONING); // -R
#if RSB_USE_MKL
			if(rsb__getenv_int_t("RSB_WANT_USE_MKL",1)==0)
#endif /* RSB_USE_MKL */
				RSB_DO_FLAG_ADD(flags_o,RSB_FLAG_USE_HALFWORD_INDICES); // 0x71 = -qH
			merge_experimental = split_experimental = RSB_CONST_DEF_MS_AT_AUTO_STEPS; /* msex */
			want_mbw = RSB_BOOL_TRUE; /* 0x776d62: wmb */
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
			want_mkl_bench = RSB_BOOL_TRUE; // 0x4C
#endif /* RSB_WANT_MKL */
#if RSB_WANT_ARMPL
			want_armpl_bench = RSB_BOOL_TRUE; // 0x4C
#endif /* RSB_WANT_ARMPL */
')dnl
			tn = 2; /* 0x4150 */
			typess=":"; /* 0x54 */
			optarg = NULL; goto lab_0x7772740a; // 0x7772740a
			break;
			case 0x4B: /* K */
			want_convert = RSB_BOOL_TRUE; /* FIXME: ignoring argument */
			break;
			case 0x55: /* U */
			want_update = RSB_BOOL_TRUE;
			break;
			case 0x4853: /* HS */
			want_as_hermitian = RSB_BOOL_TRUE;
			break;
			case 0x5353: /* SS */
			want_as_symmetric = RSB_BOOL_TRUE;
			break;
			case 0x4553: /* ES */
			want_expand_symmetry = RSB_BOOL_TRUE;
			break;
			case 0x5555: /* UU */
			want_as_unsymmetric = RSB_BOOL_TRUE;
			break;
			case 0x4F4C54: /* OLT */
			want_only_lowtri = RSB_BOOL_TRUE;
			break;
			case 0x4F4554: /* OUT */
			want_only_upptri = RSB_BOOL_TRUE;
			break;
			case 0x6363:
			/* this flag activates all interfaced libraries (if any) */
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
			want_mkl_bench = RSB_BOOL_TRUE;
#endif /* RSB_WANT_MKL */
')dnl
#if RSB_WANT_ARMPL
			want_armpl_bench = RSB_BOOL_TRUE;
			RSB_STDOUT("# Using ARMPL.\n");
			{
				//void armplversion(armpl_int_t *major, armpl_int_t *minor, armpl_int_t *patch, const char **tag);
				armplinfo();
			}
#endif /* RSB_WANT_ARMPL */
			break;
			case 0x776e6d6c:
			/* this flag deactivates all interfaced libraries (if any) */
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
			want_mkl_bench = RSB_BOOL_FALSE;
#endif /* RSB_WANT_MKL */
')dnl
#if RSB_WANT_ARMPL
			want_armpl_bench = RSB_BOOL_FALSE;
#endif /* RSB_WANT_ARMPL */
			break;
			case 0x6B: /* ncA */
			want_column_expand = rsb__util_atoi(optarg);
			break;
			case 0x6E: /* n */
			cns = optarg; /* cores (threads) numbers (specification) string */
			break;
ifelse(mop,`mat_stats',`dnl
			case 0x75 :	/* u */
			until_confidence = rsb__util_atof(optarg);
			break;
')dnl
			case 0x76: /* spmv_uauz */
			be_verbose = 1;
			break;
			case 0x774446:	/* wde */
			want_getdiag_bench = 1;
			break;
			case 0x776264:	/* wbd */
			RSBENCH_STDOUT("# EXPERIMENTAL --also-implicit-diagonal\n");
			diagin=2;
			break;
			case 0x61736875: /* ashu */
			RSBENCH_STDOUT("# EXPERIMENTAL --also-symmetries\n");
			asfin=3;
			break;
			case 0x616c6c69: /* alli */
			RSBENCH_STDOUT("# EXPERIMENTAL --also-short-idx\n");
			asiin=2;
			break;
			case 0x616363: /* acc */
			RSBENCH_STDOUT("# EXPERIMENTAL --also-coo-csr\n");
			accin=2;
			break;
			case 0x61726563: /* arec */
			RSBENCH_STDOUT("# EXPERIMENTAL --also-recursive\n");
			arecn=2;
			break;
			case 0x776E68:	/* wnh */
			want_nonzeroes_distplot = 1;
			break;
			case 0x777246:	/* wre */
			want_getrow_bench = 1;
			break;
#ifdef RSB_WANT_PERFORMANCE_COUNTERS
			case 0x707763:	/* wpc */
			want_perf_counters = 1; /* 1 is what user wants; 2 is for debug purposes */
			break;
#endif
			case 0x77707373:	/* wpss */
			want_print_per_subm_stats = RSB_BOOL_TRUE;
			break;
			case 0x776F6174:	/* woat */
			want_accuracy_test = 2;
			break;
			case 0x776e6f66:	/* wnof */
				g_fill_va_with_ones = RSB_BOOL_FALSE;
			break;
			case 0x776e7274:	/* wnrt */
			RSBENCH_STDOUT("# Requesting no autotuning via --want-no-autotune\n");
			want_autotuner = RSB_TIME_ZERO;
			wai=wat=0;
			want_autotuner = merge_experimental = split_experimental = RSB_NEGATED_EXAGGERATED_TUNER_TIMES;
			just_enter_tuning = 0;
			break;
			case 0x7772740a:	/* wrt */
lab_0x7772740a:		/* want_autotuner = rsb__util_atof(optarg); */
			{
				const char *optarg_ = optarg ? optarg : "";
				char wavv = 0x0;
				int sargs = sscanf(optarg_,"%lfs%dx%dt%c%c",&want_autotuner,&wai,&wat,&wav,&wavv);

				if(!*optarg_)
					sargs = 0;
				RSBENCH_STDOUT("# Passed %d arguments via autotuning string \"%s\" (an empty string requests defaults)\n",sargs,optarg_);
				if(sargs < 0)
				{
					RSBENCH_STDOUT("  Wrong autotuning string detected!\n");
					rsb_test_help_and_exit(argv[0],options, 0);
					exit(0);
				}
				switch(sargs)
				{
					case(EOF):
					case(0):
						want_autotuner = 10.0;
					case(1):
						wai = 1;
					case(2):
						wat = 0;
					case(3):
						wav = 0;
					case(4):
						wavv = 0;
					case(5):
					break;
				}
				/* RSBENCH_STDOUT("Got an autotuning string: %lfs%dx%dt%c%c\n",want_autotuner,wai,wat,wav,wavv); */
				if(toupper(wav)==0x56) /* V */
					wavf = RSB_AUT0_TUNING_VERBOSE;
				else
					wavf = RSB_AUT0_TUNING_SILENT ;
				if(toupper(wavv)==0x56) /* V */
					wavf++;
				if(toupper(wai)>RSB_CONST_MAX_TUNING_ROUNDS)
				{
					RSBENCH_STDOUT("Restricting the number of tuning round to %d (%d is too much!).\n",RSB_CONST_MAX_TUNING_ROUNDS,wai);
					wai = RSB_CONST_MAX_TUNING_ROUNDS;
				}
				RSBENCH_STDOUT("Will invoke autotuning for ~%lf s x %d rounds, specifying verbosity=%d and threads=%d. (>0 means no structure tuning; 0 means only structure tuning, <0 means tuning of both with (negated) thread count suggestion).\n",want_autotuner,wai,wavf,wat);
			}
			want_mkl_autotuner = want_autotuner;
			break;
#if RSB_HAVE_METIS
			case 0x776d6272:	/* wmbr */
			want_wmbr = RSB_BOOL_TRUE;
			break;
#endif
			case 0x776d6174:	/* wmat */
			sscanf(optarg?optarg:"","%lf",&want_mkl_autotuner);
			want_mkl_autotuner = RSB_MAX(1.0,want_mkl_autotuner); /* FIXME: actual value is unimportant as long as it is positive ! */
			break;
			case 0x776d6f62:	/* wmob */
			mib = 1;
			break;
			case 0x776174:	/* wat */
			want_accuracy_test = 1;
			break;
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
			case 0x775563:
			want_unordered_coo_bench = RSB_BOOL_TRUE;
			break;
')dnl
			case 0x767646:	/* wae */
			want_ancillary_execs = RSB_BOOL_TRUE;
			break;
			case 0x42767646:	/* nwae */
			want_ancillary_execs = RSB_BOOL_FALSE;
			break;
			case 0x77:	/* w */
			b_w_filename = optarg;
			break;
			case 0x63777273:	/* wcsr */
			csr_w_filename = optarg;
			break;
			case 0x77707266:
			fprfn = optarg ? optarg : "";
			want_perf_dump = RSB_BOOL_TRUE;
			if(fprfn && !*fprfn)
				fprfn = NULL;
			else
				RSBENCH_STDOUT("# performance record file set to: %s\n", fprfn);
				/* and tuning matrix renderings as well */
			break;
			case 0x776e7072:
			fprfn = NULL;
			want_perf_dump = RSB_BOOL_FALSE;
			break;
			case 0x77707261:
			apprfn = optarg;
			break;
			case 0x77707270:
			ppprfn = optarg;
			break;
			case 0x64697a65 :	/* dize */
			RSB_DO_FLAG_ADD(flags_o,RSB_FLAG_DISCARD_ZEROS);
			break;
			case 0x68: /* h */
			/* should use rsb_test_help_and_exit */
			RSBENCH_STDOUT(
				"%s "RSB_INFOMSG_SAK".\n"
				"You can use it to perform sparse matrix - unitary vector multiplication, "
				"specifying the blocking parameters, the times to perform multiplication.\n"
				"\n"
				"Additional debugging flags (-d, -p) are present.\n"
				"\n"
				"Usage : %s [OPTIONS]\n where OPTIONS are taken from "
				"[ -f filename ] \n"
				"[ -F matrix_storage=[b|c|bc] ] \n"
				"[ -r br ] \n"
				"[ -c bc ] \n"
				"[ -t TIMES ]\n"
				"[ -n OPENMP_THREADS ]\n"
				"[ -T ( S | D | I | C ) /* float, double, integer, character*/ ] \n"
				"[ -s /* will internally sort out nnzs */ ] \n"
				"[ -p /* will set to 1 nonzeros */ ] \n"
				"[-d /* if debugging on */]: \n"
				"[-A /* for auto-blocking */]: \n"
				"[ -h ] \n"
				"\n"
				"please note that not all of the suggested numerical types could be compiled in right now and/or work well.default is double.\n"
				"\n"
				"\n"
				"e.g.: %s -f raefsky4.mtx -t 10 -T :   # 10 times for each of the supported numerical types\n",
				argv[0],
				argv[0],
				argv[0]);
			rsb_test_help_and_exit(argv[0],options, 0);
			exit(0);
	    	}
	}

ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	if(incXss)
	{
		if(RSB_SOME_ERROR(rsb__util_get_bx_array_or_default(0x3A/*colon*/,dis,incXss,&incXn,&incXa)) || incXn<1)
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
		RSBENCH_STDOUT("# setting incX=%d\n",incXa[0]);
	}
	if(incYss)
	{
		if(RSB_SOME_ERROR(rsb__util_get_bx_array_or_default(0x3A/*colon*/,dis,incYss,&incYn,&incYa)) || incYn<1)
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
		RSBENCH_STDOUT("# setting incY=%d\n",incYa[0]);
	}
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	if(alphass)
	{
		errval = rsb__util_get_bx_array_or_default(0x3A/*colon*/,"-1,1,2",alphass,&alphan,&alphaip);
		if(RSB_SOME_ERROR(errval))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
	}
	if(betass)
	{
		errval = rsb__util_get_bx_array_or_default(0x3A/*colon*/,"0,1,2",betass,&betan,&betaip);
		if(RSB_SOME_ERROR(errval))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
	}
')dnl
dnl
	if(typess)
	{
		const char*toa = typess;
		ntypecodes=0; /* this neutralizes former -T ... option */
		/* if( *typess == 0x3A || *typess == 0x2A ) */ /* : or * aka colon or asterisk */
		if( ( ! isalpha(*typess) ) || ( strstr(typess,"all") != NULL ) )
			toa = RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS ;
		if( strstr(typess,"B") == typess )
			toa = RSB_BLAS_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS;
		for(;*toa;++toa)
		if(isalpha(*toa))
		{
			if(ntypecodes<maxtypes)
				typecodes[ntypecodes++]=typecode=toupper(*toa);
			else
			{
				RSB_ERROR("Up to %d types supported! P.s.: Use a punctuation symbol to ask for all supported types.\n",maxtypes);
				goto err;
			}
		}
		typecodes[ntypecodes] = RSB_NUL;
	}

	if( (!RSB_DO_FLAG_HAS(flags_o,RSB_FLAG_QUAD_PARTITIONING)) && want_recursive != RSB_BOOL_FALSE )
	{
		RSB_WARN("Assuming a recursive matrix structure is requested...\n");
		RSB_DO_FLAG_ADD(flags_o,RSB_FLAG_QUAD_PARTITIONING);
	}
dnl
	fnrprs[0]=RSB_NUL;
dnl	
	for (c = optind; c < argc; c++)                                                     
	{
		if( c == optind && optind == argc - 1 ) /* if only one and first */
			rsb__mtxfn_bncp(fnrprs,rsb__basename(argv[c]),0);
		RSB_RSBENCH_ADDF(argv[c],RSB_FAF_DEFAULTS)
	}
dnl
	if(filenamen>0 && want_filenamelist_sort == 1)
	{
		RSBENCH_STDOUT("# Sorting matrices list (use --no-sort-filenames-list to prevent this)\n");
		qsort(filenamea, filenamen, sizeof(char *), rsb__cmpstringp);
	}
	if( /* want_verbose > 0 && */ filenamen )
	{
		RSBENCH_STDOUT("# Using matrices:");
		for ( filenamei = 0; filenamei < filenamen ; ++filenamei )
			RSBENCH_STDOUT(" %s",rsb__basename(filenamea[filenamei]));
		RSBENCH_STDOUT("\n");
	}
dnl
	if(want_verbose > 0)
	{
		rsb_char_t cbuf[RSB_MAX_COMPILE_COMMAND_LENGTH];
		rsb__echo_timeandlabel(" beginning run at ","\n",&st);
		rsb__echo_cargs(argc, argv);
		errval = rsb__do_lib_get_info_str(0, &cbuf[0], sizeof(cbuf)-1);
		if(RSB_SOME_ERROR(errval))
			errval = RSB_ERR_NO_ERROR;
		else
			RSBENCH_STDOUT("# compiled with: %s\n",cbuf);
	}
dnl
	if(cn == 1)
		tcrprs = cns;
	if(cns)
	{
		ca = NULL;
		cn = 0;

		RSB_DEBUG_ASSERT( sizeof(int) == sizeof(rsb_blk_idx_t )  );
		RSB_DEBUG_ASSERT( sizeof(int) == sizeof(rsb_thread_t)  );

		if(cns && *cns==0x3A/*colon*/)
			errval=rsb__util_get_tn_array(cns,&cn,&ca);
		else
			errval=rsb__util_get_bx_array(cns,&cn,&ca);
		if(RSB_SOME_ERROR(errval))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
		if(cn == 1)
			tcrprs = cns;
	}
	else
	{
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		cn = 1;
		ca_[0] = RSB__GET_MAX_THREADS();
		RSBENCH_STDOUT("# User did not specify threads; assuming %d. Environment provides max %d threads; this build supports max %d.\n", cn, ca_[0], RSB_CONST_MAX_SUPPORTED_THREADS  );
		if( RSB_CONST_MAX_SUPPORTED_THREADS  < ca_[0] )
		{
			RSBENCH_STDOUT("# Warning: environment provides more threads than supported by this configuration -- expect trouble !\n");
		}
		RSBENCH_STDOUT("# User did not specify threads; assuming %d. Environment provides max %d threads; this build supports max %d.\n", cn, ca_[0], RSB_CONST_MAX_SUPPORTED_THREADS  );
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		tcrprs = rsb__getenv("OMP_NUM_THREADS");
	}
dnl
	printf("# average timer granularity: %2.3lg s\n",tinf);
	if(want_perf_dump)
	{
		if(!fprfn)
		{
			rsb_bool_t with_mkl = RSB_BOOL_FALSE;
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
			with_mkl = want_mkl_bench;
#endif /* RSB_WANT_MKL */
')dnl
			rsb__impcdstr(fprfnb,"rsbench_pr",ppprfn,apprfn,with_mkl,tcrprs,fnrprs,".rpr");
			fprfn = fprfnb;
		}

		if(fprfn)
		{
			const size_t sl = strlen(fprfn);

			strcpy(mtrfpf,fprfn);
			if( sl >=4 && strcmp(".rpr",mtrfpf+sl-4) == 0 )
				mtrfpf[sl-4]=RSB_NUL;
			strcat(mtrfpf,"-tuning-");
		}

		if(!cprfn)
			rsb__sprintf(cprfnb,"%s.tmp",fprfn),
			cprfn = cprfnb;
		printf("# Will write a final performance record to file %s and periodic checkpoints to %s\n",fprfn,cprfn);
	}
	if( maxtprt > RSB_TIME_ZERO )
		printf("# If program run time will exceed %2.3lg s, will attempt early termination.\n",maxtprt );

dnl	printf("# average OpenMP timer granularity: %lg\n",omp_get_wtick());
dnl
	RSBENCH_STDOUT("# will %s""perform ancillary tests.\n", want_ancillary_execs ?"":"NOT ");
	RSBENCH_STDOUT("# will flush cache memory: %s between each operation measurement series, and %s between each operation.\n", want_outer_flush?"":"NOT", want_inner_flush?"":"NOT");
	RSBENCH_STDOUT("# will %s any zero encountered in the matrix.\n", ( RSB_DO_FLAG_HAS(flags_o,RSB_FLAG_DISCARD_ZEROS) )?"discard":"keep");
dnl
	if( nrhsa == NULL ) nrhsa = &nrhs;
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	if( incXa == NULL ) incXa = (rsb_blk_idx_t *const)&incX;
	if( incYa == NULL ) incYa = (rsb_blk_idx_t *const)&incY;
')dnl
dnl
	if(want_likwid == RSB_BOOL_TRUE){RSB_LIKWID_MARKER_INIT;}

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if(ntypecodes==0)
		typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	if(ntypecodes==0)
	{
		typecodes[ntypecodes++] = typecode;
		typecodes[ntypecodes] = RSB_NUL;
	}

	io.n_pairs=0;
	if(should_use_alternate_sort!=0)
	{
		io.values[io.n_pairs]=&should_use_alternate_sort;
		io.keys[io.n_pairs]=RSB_IO_WANT_SORT_METHOD;
		io.n_pairs++;
	}
	if(should_use_cb_method!=0)
	{
		io.values[io.n_pairs]=&should_use_cb_method;
		io.keys[io.n_pairs]=RSB_IO_WANT_CACHE_BLOCKING_METHOD;
		io.n_pairs++;
	}
	if(mhs!=NULL)
	{
		io.values[io.n_pairs]=&mhs;
		io.keys[io.n_pairs]=RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING;
		io.n_pairs++;
	}
	if(subdivision_multiplier!=0.0)
	{
		io.values[io.n_pairs]=&subdivision_multiplier;
		io.keys[io.n_pairs]=RSB_IO_WANT_SUBDIVISION_MULTIPLIER;
		io.n_pairs++;
	}
#if RSB_WANT_BOUNDED_BOXES
	if(want_bounded_box==0)
	{
		io.values[io.n_pairs]=&want_bounded_box;
		io.keys[io.n_pairs]=RSB_IO_WANT_BOUNDED_BOX_COMPUTATION;
		io.n_pairs++;
	}
#endif /* RSB_WANT_BOUNDED_BOXES */
	if(want_no_leaf_spmm!=0)
	{
		io.values[io.n_pairs]=&want_no_leaf_spmm;
		io.keys[io.n_pairs]=RSB_IO_WANT_LEAF_LEVEL_MULTIVEC;
		io.n_pairs++;
	}

#if RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE
	{
		static rsb_int_t one=1; // static: needs to survive context leave
		io.values[io.n_pairs]=&one;
		io.keys[io.n_pairs]=RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE;
		io.n_pairs++;
	}
#endif /* RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE */

	if(want_verbose >= 2)
	{
		static FILE *sp = NULL;
		sp = RSB_DEFAULT_STREAM;
		io.values[io.n_pairs]=&sp;
		io.keys[io.n_pairs]=RSB_IO_WANT_VERBOSE_EXIT;
		io.n_pairs++;
	}

	if(want_verbose > 1)
	{
		static rsb_int_t one=1; // static: needs to survive context leave
		io.values[io.n_pairs]=&one;
		io.keys[io.n_pairs]=RSB_IO_WANT_VERBOSE_TUNING;
		io.n_pairs++;
	}

#if RSB_HAVE_STREAMS
	if(want_verbose >= 2)
	{
		static FILE *sp = NULL;
		sp = RSB_DEFAULT_STREAM;
		io.values[io.n_pairs]=&sp;
		io.keys[io.n_pairs]=RSB_IO_WANT_VERBOSE_INIT;
		io.n_pairs++;
	}
#endif /* RSB_HAVE_STREAMS */

#ifdef RSB_HAVE_UNISTD_H
{
	char **me = NULL;
	rsb_int_t rpevc = 0; /* RSB_ prefixed environment variables count */

	for(me=environ;*me;++me)
		if( strstr(*me,"RSB_") == *me )
			rpevc++;

	if( rpevc )
	{
		RSB_STDOUT("# The user specified %d RSB_ prefixed environment variables:\n",rpevc);
		for(me=environ;*me;++me)
			if( strstr(*me,"RSB_") == *me )
				RSB_STDOUT("#  export %s\n",*me);
	}
}
#endif /* RSB_HAVE_UNISTD_H */
	
#define RSB_TM_GETENV_STDOUT(VAR)						\
	if( rsb__getenv(VAR) )							\
		RSB_STDOUT("# env: export " VAR "=%s\n",rsb__getenv(VAR));	\
	else									\
		RSB_STDOUT("# env: " VAR " is not set\n");

	if(want_verbose >= -2)
	{
		RSB_TM_GETENV_STDOUT("PATH");
		RSB_TM_GETENV_STDOUT("LD_LIBRARY_PATH");
		RSB_TM_GETENV_STDOUT("HOSTNAME");
#if defined(RSB_WANT_OMP_RECURSIVE_KERNELS) && (RSB_WANT_OMP_RECURSIVE_KERNELS>0)
		RSB_TM_GETENV_STDOUT("KMP_AFFINITY");
		RSB_TM_GETENV_STDOUT("OMP_AFFINITY_FORMAT");
		RSB_TM_GETENV_STDOUT("OMP_ALLOCATOR");
		RSB_TM_GETENV_STDOUT("OMP_CANCELLATION");
		RSB_TM_GETENV_STDOUT("OMP_DEBUG");
		RSB_TM_GETENV_STDOUT("OMP_DEFAULT_DEVICE");
		RSB_TM_GETENV_STDOUT("OMP_DISPLAY_ENV");
		RSB_TM_GETENV_STDOUT("OMP_DISPLAY_AFFINITY");
		RSB_TM_GETENV_STDOUT("OMP_DYNAMIC");
		RSB_TM_GETENV_STDOUT("OMP_MAX_ACTIVE_LEVELS");
		RSB_TM_GETENV_STDOUT("OMP_MAX_TASK_PRIORITY");
		RSB_TM_GETENV_STDOUT("OMP_NESTED");
		RSB_TM_GETENV_STDOUT("OMP_NUM_THREADS");
		RSB_TM_GETENV_STDOUT("OMP_PLACES");
		RSB_TM_GETENV_STDOUT("OMP_PROC_BIND");
		RSB_TM_GETENV_STDOUT("OMP_SCHEDULE");
		RSB_TM_GETENV_STDOUT("OMP_STACKSIZE");
		RSB_TM_GETENV_STDOUT("OMP_TARGET_OFFLOAD");
		RSB_TM_GETENV_STDOUT("OMP_THREAD_LIMIT");
		RSB_TM_GETENV_STDOUT("OMP_TOOL");
		RSB_TM_GETENV_STDOUT("OMP_TOOL_LIBRARIES");
		RSB_TM_GETENV_STDOUT("OMP_WAIT_POLICY");
	//	tcrprs = rsb__set_num_threads() ;
#else
		RSB_STDOUT("# serial build: ignoring environment variables: KMP_AFFINITY OMP_PROC_BIND OMP_NUM_THREADS\n");
#endif
		RSB_TM_GETENV_STDOUT("RSB_WANT_RSBPP");
		{
#if RSB_USE_LIBRSBPP
			const int wrp = rsb__getenv_int_t("RSB_WANT_RSBPP",1);
			if(wrp==0)
				RSBENCH_STDOUT("# NOT using kernels from librsbpp (opted off via environment variable).\n");
			else
				RSBENCH_STDOUT("#     using kernels from librsbpp (default).\n");
#else
			const int wrp = rsb__getenv_int_t("RSB_WANT_RSBPP",0);
			if(wrp!=0)
			{
				RSBENCH_STDOUT("Error: RSB_WANT_RSBPP environment variable set but librsbpp configured out!\n");
				goto err;
			}
#endif /* RSB_USE_LIBRSBPP */
		}
		RSB_TM_GETENV_STDOUT("SLURM_CLUSTER_NAME");
		RSB_TM_GETENV_STDOUT("SLURM_CPUS_ON_NODE");
		RSB_TM_GETENV_STDOUT("SLURM_JOB_CPUS_PER_NODE");
		RSB_TM_GETENV_STDOUT("SLURM_JOB_ID");
		RSB_TM_GETENV_STDOUT("SLURM_JOBID");
		RSB_TM_GETENV_STDOUT("SLURM_JOB_NAME");
		RSB_TM_GETENV_STDOUT("SLURM_JOB_NUM_NODES");
		RSB_TM_GETENV_STDOUT("SLURM_JOB_PARTITION");
		RSB_TM_GETENV_STDOUT("SLURM_NPROCS");
		RSB_TM_GETENV_STDOUT("SLURM_NTASKS");
		RSB_TM_GETENV_STDOUT("SLURM_STEP_TASKS_PER_NODE");
		RSB_TM_GETENV_STDOUT("SLURM_TASKS_PER_NODE");
		if(1)
		{
			rsb_char_t hn[RSB_MAX_HOSTNAME_LEN];
			rsb__strcpy_hostname(hn);
			if(hn[0])
				RSBENCH_STDOUT("# detected hostname: %s\n",hn);
		}

		if( want_verbose >= 0 )
			RSBENCH_STDOUT("# user specified a verbosity level of %d (each --verbose occurrence counts +1)\n",want_verbose_tuning );
		else
			RSBENCH_STDOUT("# user did not specify any verbosity level (each --verbose occurrence counts +1)\n");
	}

	if((errval = rsb_lib_reinit(iop))!=RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err,"Error while reinitializing the library.");
	}

	if ( want_mbw == RSB_BOOL_TRUE ) 
	{
		rsb_time_t dt = RSB_TIME_ZERO;
		dt = - rsb_time();
		if((errval = rsb__memory_benchmark(&mbet))!=RSB_ERR_NO_ERROR)
		{
			RSB_ERROR("Error while performing bandwidth benchmark!");
			goto err;
		}
		else
			mbetp = &mbet; // FIXME: need option to turn this feature on/off
		dt += rsb_time();
		if( want_verbose >= 1 )
			RSB_STDOUT("# Memory benchmark took %.3lfs\n",dt);
	}

#ifdef RSB_HAVE_SETENV
{
	/* special environment variables set just for the sake of being saved in .rpr files */
#ifdef RSB_CC
	setenv("RSB_CC",RSB_CC,1);
#endif /* RSB_CC */
#ifdef RSB_CFLAGS
	setenv("RSB_CFLAGS",RSB_CFLAGS,1);
#endif /* RSB_CFLAGS */
#ifdef RSB_DETECTED_MEM_HIERARCHY_INFO
	setenv("RSB_DETECTED_MEM_HIERARCHY_INFO",RSB_DETECTED_MEM_HIERARCHY_INFO,1);
#endif /* RSB_DETECTED_MEM_HIERARCHY_INFO */
	{
		rsb_char_t buf[RSB_MAX_LINE_LENGTH];
		rsb_char_t usmhib[RSB_MAX_LINE_LENGTH];
		buf[0]=RSB_NUL;
		strcat(buf,rsb__get_mem_hierarchy_info_string(usmhib));
		//rsb_lib_get_opt(RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING,buf);//FIXME
		setenv("RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING",&buf[0],1);
	}
}
#endif /* RSB_HAVE_SETENV */

#if RSB_WANT_PERFORMANCE_COUNTERS_IN_RSBENCH 
	if((errval = rsb_perf_counters_init())!=RSB_ERR_NO_ERROR)
	{
		RSBENCH_STDERR("problem initializing performance counters (rsb_perf_counters_init gave %d)\n",(int)errval);
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
#endif

	if( RSB__APPROPRIATE_AT_TIME_SPEC( want_autotuner ) || RSB__APPROPRIATE_AT_TIME_SPEC( merge_experimental ) || RSB__APPROPRIATE_AT_TIME_SPEC( split_experimental ) )
	{
		RSB_STDOUT("# auto-tuning oriented output implies  times==0 iterations and sort-after-load.\n");
		times = 0;
		/* if(want_verbose>0) */
		want_impatiently_soon_pre_results = 1;
		want_sort_after_load = RSB_BOOL_TRUE;
	}
	else
	if( times < 1 )
	{
		RSB_STDOUT("# The iteration times should be specified as a positive number!\n");
		RSB_ERROR(RSB_ERRM_ES);
		goto err;
	}
	else
		RSB_STDOUT("# Will measure on times=%d iterations.\n",times);

	if( 0 == filenamen )
#if RSB_RSBENCH_STATIC_FILENAMEA
	       	filenamea[0] = &(fnbuf[0]);
#else
	       	filenamea = &fnbufp;
#endif
	filenamen = RSB_MAX(1,filenamen);

ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
	if( RSB_MKL_APPROPRIATE_AT_TIME_SPEC( want_mkl_autotuner ) )
		want_mkl_bench_csr = RSB_BOOL_FALSE;
#endif /* RSB_WANT_MKL */
')dnl

	if(want_perf_dump) 
		rsb__pr_init(&rspr, NULL, filenamen, cn, incXn, incYn, nrhsn, ntypecodes, tn, mbetp );

	for(     filenamei=0;     filenamei<filenamen+want_impatiently_soon_pre_results  ;++filenamei     )
	{
		if( /*filenamea &&*/ ( filenamea[filenamei] != filename_old) && filename_old && want_impatiently_soon_pre_results && want_perf_dump && filenamei>0 && filenamen>1) 
		{
			int filenameif = filenamei-1;
			RSBENCH_STDOUT("# ====== BEGIN Impatient results record for matrix %d/%d: %s.\n",filenamei,filenamen,rsb__basename(filename_old));
			errval = rsb__pr_dump_inner(rspr, filenamea, ca, incXa, incYa, nrhsa, typecodes, NULL,&filenameif, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, NULL, NULL);
			RSBENCH_STDOUT("# ======  END  Impatient results record for matrix %d/%d: %s.\n",filenamei,filenamen,rsb__basename(filename_old));
			if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
			if( filenameif > 0 && filenameif < filenamen-1) /* not after first and not at last */
				RSBENCH_STDOUT("# ====== BEGIN Impatient summary record for the %d/%d matrices so far.\n", filenameif+1,filenamen),
				errval = rsb__pr_dump_inner(rspr, filenamea, ca, incXa, incYa, nrhsa, typecodes, NULL, NULL,&filenameif, NULL, NULL, NULL, NULL, NULL, NULL, NULL, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, NULL, NULL),
				RSBENCH_STDOUT("# ======  END  Impatient summary record for the %d/%d matrices so far.\n", filenameif+1,filenamen);
			if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
			errval = rsb__pr_save(cprfn, rspr, filenamea, ca, incXa, incYa, nrhsa, typecodes, NULL, RSB_BOOL_TRUE);
			if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}

		if( filenamei >= filenamen )
			continue; /* temporary: only for the want_impatiently_soon_pre_results trick */

		if(filenamea[filenamei])
		{
			filename = filenamea[filenamei];
		}

		if(filenamen>1)
		{
			RSBENCH_STDOUT("# multi-file benchmarking (file %d/%d) -- now using %s\n",filenamei+1,filenamen,rsb__basename(filename));
		}

ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
		if(!alphaip)
			errval=rsb__util_get_bx_array("1",&alphan,&alphaip);
		if(!betaip)
			errval=rsb__util_get_bx_array("1",&betan,&betaip);
')dnl

	for(areci=0;areci<arecn;++areci)
	{
dnl
	for(accii=0;accii<accin;++accii)
	{
dnl
	for(asiii=0;asiii<asiin;++asiii)
	{
dnl
	for(asfii=0;asfii<asfin;++asfii)
	{
dnl
	for(diagii=0;diagii<diagin;++diagii)
	{
dnl
	for(typecodesi=0;typecodesi<ntypecodes;++typecodesi)
	{
dnl
	rsb_flags_t flags = flags_o;
	rsb_thread_t cl; /* cores number last (overrides cn for this typecode cycle) */
	typecode = typecodes[typecodesi];
dnl
	if(areci==0 && arecn > 1)
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_QUAD_PARTITIONING);
	if(areci==1 && arecn > 1)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_QUAD_PARTITIONING);
dnl
	if(accii==0 && accin > 1)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_HALFWORD_INDICES);
	if(accii==1 && accin > 1)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_USE_FULLWORD_INDICES);
dnl
	if(asiii==0 && asiin > 1)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);
	if(asiii==1 && asiin > 1)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_COO_STORAGE);
dnl
	if(asfii==0 && asfin > 1)
		want_as_unsymmetric = RSB_BOOL_TRUE, want_as_symmetric = RSB_BOOL_FALSE, want_as_hermitian = RSB_BOOL_FALSE;
	if(asfii==1)
		want_as_unsymmetric = RSB_BOOL_FALSE, want_as_symmetric = RSB_BOOL_TRUE, want_as_hermitian = RSB_BOOL_FALSE;
	if(asfii==2)
		want_as_unsymmetric = RSB_BOOL_FALSE, want_as_symmetric = RSB_BOOL_FALSE, want_as_hermitian = RSB_BOOL_TRUE;
dnl
	if(diagii>0)
		RSB_DO_FLAG_XOR(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT);
dnl
	if(ntypecodes>1)
	{
		if(want_verbose >= -2)
			RSBENCH_STDOUT("# multi-type benchmarking (%s) -- now using typecode %c (last was %c).\n",typecodes,typecode,typecode_old);
		if( RSB_MATRIX_UNSUPPORTED_TYPE ( typecode ) )
		{
			RSBENCH_STDOUT("# Skipping unsupported type \"%c\" -- please choose from \"%s\".\n",typecode,RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS );
			continue;
		}
	}

	if(want_verbose >= -2)
	{
		RSBENCH_STDOUT("# Cache block size total %ld bytes, per-thread %ld bytes\n",rsb__get_lastlevel_c_size(),rsb__get_lastlevel_c_size_per_thread()); /* as used in rsb__allocate_recursive_sparse_matrix_from_row_major_coo */
 		RSBENCH_STDOUT("# so far, program took %.3lfs of wall clock time; ancillary tests %.3lfs; I/O %.3lfs; checks %.3lfs; conversions %.3lfs; rsb/mkl tuning %.3lfs/%.3lfs ",totprt + rsb_time(),totatt,totiot,totht,totct,tottt,totmt);
		RSBENCH_MEM_ALLOC_INFO("")
		RSBENCH_STDOUT(".\n");
	}

	if(cns)
	{
		cc = ca[0];
		if( cc == 0 )
			RSBENCH_STDOUT("# Using auto threads\n");
		else
			RSBENCH_STDOUT("# Using %d threads\n",(int)cc);
		RSB_DEBUG_ASSERT( cc == 0 || RSB_IS_VALID_THREAD_COUNT ( cc ) );
	}
	cl=cn;
	if(bcs)
		if(RSB_SOME_ERROR(rsb__util_get_bx_array(bcs,&bcl,&bcv)))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
	if(brs)
		if(RSB_SOME_ERROR(rsb__util_get_bx_array(brs,&brl,&brv)))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}


#ifdef RSB_WANT_OSKI_BENCHMARKING 
	/* FIXME : note that this option is not compatible with g_sort_only .. */
        oski_Init();
#endif /* RSB_WANT_OSKI_BENCHMARKING */
	g_debug = ((flags & RSB_FLAG_SHOULD_DEBUG) != 0);

	if(g_sort_only)RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORT_INPUT);

	if( RSB_MATRIX_UNSUPPORTED_TYPE ( typecode ) )
	{
		RSBENCH_STDERR("error : please recompile with double precision floating point numbers supported! \n");
		return RSB_ERR_GENERIC_ERROR;
	}
	rsb__util_set_area_to_converted_integer(&pone[0],typecode,+1);

dnl	if(g_dump_performance_profile)
dnl	{
dnl		if((errval = rsb_do_dump_performance_record_for_op_and_type(
dnl			RSB_NUMERICAL_TYPE_INDEX_FROM_CODE(typecode), 
dnl			RSB_NUMERICAL_OP_INDEX_FROM_CODE(RSB_M4_OPTYPE_INDEX_PREPROCESSOR_SYMBOL(mop))))!=RSB_ERR_NO_ERROR)
dnl			{RSB_ERROR(RSB_ERRM_ES);goto err;}
dnl		goto done;
dnl	}

ifelse(mop,`mat_stats',`dnl
	if(until_confidence && g_estimate_fillin)
	{
		RSBENCH_STDERR("cannot perform -e functionality in one run. one at a time please..\n");
		goto err;
	}
')dnl

	if(brl<1) { /* this is a hack */ brv = rua; brl = RSB_ROWS_UNROLL_ARRAY_LENGTH;}
	if(bcl<1) { /* this is a hack */ bcv = cua; bcl = RSB_COLUMNS_UNROLL_ARRAY_LENGTH;}

	if(RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		RSBENCH_STDERR("This numerical type is not supported.\n");
		goto err;
	}

	/* CONDITIONALLY, GENERATING A MATRIX */
	if(should_generate_dense!=0)
	{
		rsb_nnz_idx_t dim = RSB_FABS(should_generate_dense);
		rsb_nnz_idx_t spacing = want_generated_spacing>1?want_generated_spacing:1;
		
		if(((should_generate_lband>-1) || (should_generate_uband>-1)) && should_generate_dense>0)
		{
			rsb__sprintf(fnbuf,"banded-%zdx%zd-%zd+%zd-%zdnz-spaced-%zd",(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)should_generate_lband,(rsb_printf_int_t)should_generate_uband,(rsb_printf_int_t)(RSB_NNZ_OF_BANDED(dim,should_generate_lband,should_generate_uband)),(rsb_printf_int_t)spacing);
		}
		else
		{
		if(want_generated_spacing>0)
		{
			if(should_generate_dense>0)
				rsb__sprintf(fnbuf,"dense-%zdx%zd-%zdnz",(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)(should_generate_dense_nc*spacing)/*dim*spacing*/,(rsb_printf_int_t)(dim*dim));
			else
				rsb__sprintf(fnbuf,"lower-%zdx%zd-%zdnz-spaced-%zd",(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)((dim*(dim-1))/2+dim),(rsb_printf_int_t)spacing);
		}
		else
		{
			if(should_generate_dense>0)
				rsb__sprintf(fnbuf,"dense-%zdx%zd-%zdnz",(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)(should_generate_dense_nc*spacing/*dim*spacing*/),(rsb_printf_int_t)(dim*should_generate_dense_nc));
			else
				rsb__sprintf(fnbuf,"lower-%zdx%zd-%zdnz",(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)(dim*spacing),(rsb_printf_int_t)((dim*(dim-1))/2+dim));
		}
		}
/*		rsb__sprintf(fnbuf,"dense-%dx%d-%dnz",dim,dim,dim*dim);*/
/*		rsb__sprintf(fnbuf,"dense-%dx%d",dim,dim);*/
		filename=&(fnbuf[0]);
		// RSB_RSBENCH_ADDF(filename,RSB_FAF_CHKDUP|RSB_FAF_CHKREC|RSB_FAF_VRBADD); /* FIXME */
	}

	if(usfnbuf)
		filename=usfnbuf;

	/* CONDITIONALLY, READING A MATRIX FROM FILE */
if(!(filename || b_r_filename))
{
	RSBENCH_STDOUT("%s (mop) : Please specify a matrix filename (with -f)\n",argv[0]);
}
else
{

	rsb_blk_idx_t M_b=0;/* was 0 */
	rsb_blk_idx_t K_b=0;
	rsb_nnz_idx_t i=0;

	rsb_coo_idx_t *p_r=NULL,*p_c=NULL;	/* FIXME : get rid of these */
	void *lhs=NULL,*rhs=NULL;
	int bcvi=0;
	int brvi=0;
	rsb_time_t frt = RSB_TIME_ZERO;

	if( filename != filename_old )
	{
		RSB_CONDITIONAL_FREE(IA);
		RSB_CONDITIONAL_FREE(JA);
		RSB_CONDITIONAL_FREE(VA);
	}
	if(!should_recycle_io) { RSB_DEBUG_ASSERT( VA == NULL ); }
	if( should_recycle_io && VA && filename == filename_old )
	{
		flags = r_flags;
dnl
		if( typecode != typecode_old )
		{
			void *VA_ = rsb__malloc_vector(nnz,typecode);
			errval = rsb__do_copy_converted_scaled(VA, VA_, NULL, typecode_old, typecode, nnz, RSB_DEFAULT_TRANSPOSITION);
			if(RSB_SOME_ERROR(errval)) { RSB_ERROR(RSB_ERRM_ES);goto err; }
			RSB_CONDITIONAL_FREE(VA);
			VA = VA_;
			RSBENCH_STDOUT("# Reusing type converted (%c->%c) arrays from last iteration instead of reloading matrix file.\n",typecode_old,typecode);
			typecode_old = typecode;
		}
		else
		{
			RSBENCH_STDOUT("# Reusing same type     (type %c) arrays from last iteration instead of reloading matrix file.\n",typecode);
		}
dnl
		goto have_va_ia_ja;
	}
dnl
	if((!should_generate_dense) && (!b_r_filename))
	{
		rsb_bool_t is_symmetric = RSB_BOOL_FALSE;
		rsb_bool_t is_hermitian = RSB_BOOL_FALSE;
		const size_t fsz = rsb__sys_filesize(filename);

		frt = - rsb_time();

ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
			{
				/* FIXME : we remove symmetry flags, for they are incompatible with triangular solve */
				RSB_DO_FLAG_DEL(flags,RSB_FLAG_SYMMETRIC);
				RSB_DO_FLAG_DEL(flags,RSB_FLAG_HERMITIAN);
			/*
				if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER))
				{
					RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER_TRIANGULAR);
				}
				else
					RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER_TRIANGULAR);
			*/
				//RSB_DO_FLAG_ADD(flags,RSB_FLAG_DISCARD_ZEROS) ;//problematic : FIXME
			}
')dnl
#ifdef RSB_HAVE_REGEX_H
		if( want_slmr && rsb_regexp_match(rsb__basename(filename),want_slmr) == RSB_BOOL_TRUE )
		{
			RSB_STDOUT("# skipping loading matrix file %s, because it matches regex /%s/.\n",filename,want_slmr);
			goto nfnm;
		}
#endif /* RSB_HAVE_REGEX_H */
		if( want_slss && ( strstr( rsb__basename(filename), want_slss ) != NULL ) )
		{
			RSB_STDOUT("# skipping loading matrix file %s, because it matches substring %s.\n",filename,want_slss);
			goto nfnm;
		}
		if(RSB_SOME_ERROR(rsb__util_mm_info_matrix_f(filename,&nrA,&ncA,&nnz,NULL,&is_symmetric,&is_hermitian,NULL,NULL,NULL,NULL)) )
		{
			RSBENCH_STDERR(RSB_ERRMSG_PROIFAMM ": %s ..\n",filename);
			if( ignore_failed_fio )
			{
				RSBENCH_STDERR("Will ignore error and continue with the following files.\n");
				errval = RSB_ERR_NO_ERROR;
				goto nfnm;
			}
			goto err;
		}
		if( want_slnu == RSB_BOOL_TRUE && ( is_hermitian || is_symmetric ) )
		{
			RSB_STDOUT("# skipping loading not unsymmetric matrix %s, as requested.\n",filename);
			goto nfnm;
		}
		if( want_slsm == RSB_BOOL_TRUE && is_symmetric )
		{
			RSB_STDOUT("# skipping loading symmetric matrix %s, as requested.\n",filename);
			goto nfnm;
		}
		if( want_slhm == RSB_BOOL_TRUE && is_hermitian )
		{
			RSB_STDOUT("# skipping loading hermitian matrix %s, as requested.\n",filename);
			goto nfnm;
		}
		if( want_slum == RSB_BOOL_TRUE && !is_symmetric )
		{
			RSB_STDOUT("# skipping loading unsymmetric matrix %s, as requested.\n",filename);
			goto nfnm;
		}
		if( want_slmn > 0 && want_slmn <  nnz )
		{
			RSB_STDOUT("# skipping loading matrix %s, having %zd > %zd allowed nonzeroes.\n",filename,(rsb_printf_int_t)nnz,(rsb_printf_int_t)want_slmn);
			goto nfnm;
		}
		if( want_slms > 0 && want_slms <= fsz / 1024 )
		{
			RSB_STDOUT("# skipping loading matrix %s, having %zd>=%zd allowed filesize (KiB).\n",filename,(size_t)fsz,(size_t)want_slms);
			goto nfnm;
		}
		if( want_slnn > 0 && want_slnn >  nnz )
		{
			RSB_STDOUT("# skipping loading matrix %s, having %zd < %zd allowed nonzeroes.\n",filename,(rsb_printf_int_t)nnz,(rsb_printf_int_t)want_slnn);
			goto nfnm;
		}
	
		RSB_STDOUT("# reading %s (%zd bytes / %zd "RSB_MEGABYTE_SYM" / %zd nnz / %zd rows / %zd columns / %zd MiB COO) as type %c...\n",rsb__basename(filename),fsz,RSB_DIV(fsz,RSB_MEGABYTE),(size_t)nnz,(size_t)nrA,(size_t)ncA,(size_t)RSB_DIV(RSB_UTIL_COO_OCCUPATION(nrA,ncA,nnz,typecode),RSB_MEGABYTE),typecode);

		if( ( nrA == ncA ) && ( nrA > 1 ) && ( want_only_lowtri || want_only_upptri ) )
			nnz += nrA;	/* the loading routine shall allocate nnz+nrA */
		else
 			nnz = 0;	/* the loading routine should determine nnz */

		totiot -= rsb_time();
		errval = rsb__util_mm_load_matrix_f(filename,&IA,&JA,&VA,&nrA,&ncA,&nnz,typecode,flags,NULL,NULL);
		totiot += rsb_time();
		if(RSB_SOME_ERROR(errval))
		{
			RSBENCH_STDERR(RSB_ERRMSG_NOTMTXMKT" : %s ..\n",filename);
			goto err;
		}
		else
		{
			rsb_bool_t is_lower = RSB_BOOL_FALSE;
			rsb_bool_t is_upper = RSB_BOOL_FALSE;
			rsb_bool_t is_vector = RSB_BOOL_FALSE;

			filename_old = filename;
			typecode_old = typecode;

			frt += rsb_time();
			RSB_STDOUT("# file input of %s took %6.2lf s (%.0lf nnz, %.0lf nnz/s ) (%.2lf MB/s ) \n",rsb__basename(filename),frt,
				(((double)nnz)),
				(((double)nnz)/frt),
				(((double)rsb__sys_filesize(filename))/(frt*RSB_INT_MILLION))
			);

			if( mtx_sample_rate < 100.0 && nnz > 0 )
			{
				const rsb_nnz_idx_t snnz = RSB_MAX( (mtx_sample_rate * ( nnz / 100.0 )), 1);
				const rsb_real_t isf = 100.0 / mtx_sample_rate;
				rsb_coo_idx_t dnzi, snzi;
				struct rsb_coo_mtx_t coo = {IA,JA,0,0,nnz,VA,typecode};

				if( snnz < 0 || snnz > nnz )
				{
					errval = RSB_ERR_BADARGS;
					RSB_PERR_GOTO(err,RSB_ERRM_BADARGS);
				}

				for( snzi=0, dnzi=0 ; snzi < nnz; (++dnzi) , snzi = dnzi*isf )
				{
					RSB_COO_MEMCPY(VA,IA,JA,VA,IA,JA,dnzi,snzi,1,RSB_SIZEOF(typecode));
				}
				RSB_COO_MEMCPY(VA,IA,JA,VA,IA,JA,snnz-1,nnz-1,1,RSB_SIZEOF(typecode));

				RSBENCH_STDOUT("# Matrix sampling: using only %zd nonzeroes out of read %zd.\n",(rsb_printf_int_t)snnz,(rsb_printf_int_t)nnz);
				if( NULL == rsb__reallocate_coo_matrix_t(&coo, snnz) )
				{
					errval = RSB_ERR_ENOMEM;
					RSB_PERR_GOTO(err,RSB_ERRM_ENOMEM);
				}
				VA = coo.VA;
				IA = coo.IA;
				JA = coo.JA;
				nnz = snnz;
			}

			if ( g_fill_va_with_ones )
			if(RSB_SOME_ERROR(errval = rsb__fill_with_ones(VA,typecode,nnz,1)))
			{
				RSB_PERR_GOTO(err,rsb__get_errstr_ptr(errval));
			}

			if ( is_symmetric || is_hermitian )
			if ( want_expand_symmetry )
			{
				// Note: might want to move this to internals, e.g. matrix ctor
				struct rsb_coo_mtx_t coo = {IA,JA,0,0,nnz,VA,typecode};

				if( RSB_NNZ_MUL_OVERFLOW(nnz,2) )
				{
					errval = RSB_ERR_INTERNAL_ERROR;
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
				}

				if( NULL == rsb__reallocate_coo_matrix_t(&coo, nnz*2) )
				{
					errval = RSB_ERR_ENOMEM;
					RSB_PERR_GOTO(err,RSB_ERRM_ENOMEM);
				}

				VA = coo.VA;
				IA = coo.IA;
				JA = coo.JA;
				RSB_COO_MEMCPY(VA,IA,JA,VA,JA,IA,nnz,0,nnz,RSB_SIZEOF(typecode)); // transposed copy
				nnz *= 2;
				RSBENCH_STDOUT("# Expanded symmetry to %zd nnz (to be cleansed of diagonal duplicates). Deleting and symmetry / hermitianness flags.\n",(size_t)nnz);
				RSB_DO_FLAG_DEL(flags,RSB_FLAG_SYMMETRIC);
				RSB_DO_FLAG_DEL(flags,RSB_FLAG_HERMITIAN);
				is_symmetric = RSB_BOOL_FALSE;
				is_hermitian = RSB_BOOL_FALSE;
			}

			if (want_io_only)
			{
				/*  */
				goto rret;
			}

			if(want_transpose)
			{
				RSB_SWAP(rsb_coo_idx_t*,IA,JA);
				RSB_SWAP(rsb_coo_idx_t,nrA,ncA);
				flags = rsb__do_flip_uplo_flags(flags);
			}

			if( nrA==ncA && nrA>1 && ( want_only_lowtri || want_only_upptri ) )
			{
				rsb_nnz_idx_t discarded = 0;
				struct rsb_coo_mtx_t coo = {IA,JA,0,0,nnz,VA,typecode};

				if( NULL == rsb__reallocate_coo_matrix_t(&coo, nnz+nrA) )
				{
					errval = RSB_ERR_ENOMEM;
					RSB_PERR_GOTO(err,RSB_ERRM_ENOMEM);
				}

				VA = coo.VA;
				IA = coo.IA;
				JA = coo.JA;

				RSBENCH_STDOUT("# excluding a triangle and forcibly adding diagonal elements (duplicates will be removed)\n");
				RSB_FCOO_ISET(IA+nnz,0,nrA); // rsb__util_coo_array_set_sequence(IA+nnz,nrA,0,1); // may want to use this instead
				RSB_FCOO_ISET(JA+nnz,0,nrA); // rsb__util_coo_array_set_sequence(JA+nnz,nrA,0,1); // may want to use this instead
				rsb__fill_with_ones(((rsb_byte_t*)VA)+RSB_SIZEOF(typecode)*nnz,typecode,nrA,1);
				nnz += nrA;

				if( want_only_lowtri )
				{
					RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER_TRIANGULAR);
					errval = rsb__weed_out_non_lowtri(VA,IA,JA,nnz,typecode,NULL,&discarded);
					RSBENCH_STDOUT("# discarded %zd non lower elements of %zd.\n",(rsb_printf_int_t)discarded,(rsb_printf_int_t)nnz);
					nnz-=discarded;
				}
				if( want_only_upptri )
				{
					RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER_TRIANGULAR);
					errval = rsb__weed_out_non_upptri(VA,IA,JA,nnz,typecode,NULL,&discarded);
					RSBENCH_STDOUT("# discarded %zd non upper elements of %zd.\n",(rsb_printf_int_t)discarded,(rsb_printf_int_t)nnz);
					nnz-=discarded;
				}

				if(RSB_SOME_ERROR(errval))
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}

			if(RSB_SOME_ERROR(rsb__util_mm_info_matrix_f(filename,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&is_lower,&is_upper,&is_vector) ))
			{
				RSBENCH_STDERR(RSB_ERRMSG_PROIFAMM ": %s ..\n",filename);
				goto err;
			}
			if( is_vector )
			{
				RSBENCH_STDERR("file %s seems to store a vector\n",filename);
				goto err;
			}
			if(RSB_BOOL_AND(want_as_unsymmetric,want_as_symmetric))
			{
				RSBENCH_STDERR("requiring both symmetric and unsymmetric flags is contradictory!\n");
				goto err;
			}
			if(want_as_unsymmetric)
			{
				is_symmetric = RSB_BOOL_FALSE;
				is_hermitian = RSB_BOOL_FALSE;
			}
			if(want_as_symmetric || want_as_hermitian)
			{
				is_symmetric = RSB_BOOL_TRUE;
				is_hermitian = RSB_BOOL_TRUE;
			}
			if(!RSB_IS_MATRIX_TYPE_COMPLEX(typecode) && is_hermitian)
			{
				RSBENCH_STDOUT("# Warning: non complex matrix with hermitian flags! Converting to symmetric!\n");
				is_hermitian = RSB_BOOL_FALSE;
				is_symmetric = RSB_BOOL_TRUE;
			}
			if( RSB_IS_MATRIX_TYPE_COMPLEX(typecode) && is_symmetric && is_hermitian )
			{
				RSBENCH_STDOUT("# Warning: complex matrix with hermitian and symmetric flags! Converting to symmetric!\n");
				is_hermitian = RSB_BOOL_TRUE;
				is_symmetric = RSB_BOOL_FALSE;
			}
			/* TODO: use rsb__flags_from_props() */
			if(is_hermitian == RSB_BOOL_TRUE && !RSB_EXPERIMENTAL_EXPAND_SYMMETRIC_MATRICES_BY_DEFAULT)
			{
dnl				RSBENCH_STDOUT("# exploiting EXPERIMENTAL symmetry\n");
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_HERMITIAN);
			}
			if(is_symmetric == RSB_BOOL_TRUE && !RSB_EXPERIMENTAL_EXPAND_SYMMETRIC_MATRICES_BY_DEFAULT)
			{
dnl				RSBENCH_STDOUT("# exploiting EXPERIMENTAL symmetry\n");
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_SYMMETRIC);
			}

			if( (!RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER)) && (!RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER)) )
			{
				/* is_upper and is_lower as declared in the matrix file */
				if(is_upper)
 					RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
				if(is_lower)
 					RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
			}
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_cleanup_nnz(VA,IA,JA,nnz,0,0,nrA,ncA,&nnz,typecode,flags)); /* NEW */
			if(RSB_SOME_ERROR(errval))
			{ RSB_ERROR(RSB_ERRM_ES); goto err; }
			if(want_sort_after_load)
			{
				rsb_time_t dt = RSB_TIME_ZERO, ct = RSB_TIME_ZERO;
				dt = - rsb_time();
dnl				//if((errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,nrA,ncA,typecode,flags))!=RSB_ERR_NO_ERROR)
				if((errval = rsb__util_sort_row_major_parallel(VA,IA,JA,nnz,nrA,ncA,typecode,flags_s))!=RSB_ERR_NO_ERROR)
				{ RSB_ERROR(RSB_ERRM_ES); goto err; }
				dt += rsb_time();
				RSBENCH_STDOUT("#pre-sorting (%zd elements) took %lg s\n",(size_t)nnz,dt);
				dt = -rsb_time();
				nnz = rsb__weed_out_duplicates (IA,JA,VA,nnz,typecode,RSB_FLAG_SORTED_INPUT);
				dt += rsb_time();
				ct = -rsb_time();
				if(RSB_SOME_ERROR(rsb__util_is_sorted_coo_as_row_major(IA,JA,nnz,typecode,NULL,RSB_FLAG_NOFLAGS)))
				{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
				ct += rsb_time();
				RSBENCH_STDOUT("#weeding duplicates (to %zd elements) took %lg s (and check, %lg s )\n",(size_t)nnz,dt,ct);
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT);
			}
#if RSB_HAVE_METIS
			if(want_wmbr)
			{
				/* FIXME: unfinished */
				rsb_coo_idx_t *perm = NULL,*iperm = NULL,*vwgt = NULL;

				perm  = rsb__calloc(sizeof(rsb_coo_idx_t)*(nrA+1));
				iperm = rsb__calloc(sizeof(rsb_coo_idx_t)*(nrA+1));
#if 1
				vwgt  = rsb__calloc(sizeof(rsb_coo_idx_t)*(nnz));
				rsb__util_coo_array_set(vwgt,nnz,0);
#else
				vwgt  = rsb__clone_area(JA,nnz*sizeof(rsb_coo_idx_t));
#endif
				if( !perm || !iperm || !vwgt )
				{
					RSB_CONDITIONAL_FREE(iperm);
					RSB_CONDITIONAL_FREE(perm);
					RSB_CONDITIONAL_FREE(vwgt);
				}
				errval = rsb__util_sort_row_major_parallel(VA,IA,JA,nnz,nrA,ncA,typecode,flags_s);
				errval = rsb__do_switch_fullword_array_to_compressed(IA,nnz,nrA);
				RSBENCH_STDOUT("Calling METIS_NodeND\n");
				/*errval = */ METIS_NodeND(&nrA,IA,JA,vwgt,NULL,perm,iperm); /* Scotch wrapper crashes on vwgt=NULL. and is void */
				RSBENCH_STDOUT("Exited  METIS_NodeND with code %d\n",errval);
				/* if(errval == METIS_OK) */
				{
					RSBENCH_STDOUT("Permuting..\n");
					errval = rsb__do_switch_compressed_array_to_fullword_coo(IA, nrA, 0, NULL);
					errval = rsb__do_permute_rows_with_coo_index( IA, perm, nnz);
					RSBENCH_STDOUT("Permuted.\n");
					/* 
					 */
					for(i=0;i<nrA;++i){ RSB_STDOUT("%d\n",perm[i]);}
				}
				RSB_CONDITIONAL_FREE(vwgt);
				RSB_CONDITIONAL_FREE(perm);
				RSB_CONDITIONAL_FREE(iperm);
			}
			
#endif /* RSB_HAVE_METIS */
		}
	}
	else
	if(should_generate_dense!=0)
	{
		rsb_nnz_idx_t dim = RSB_FABS(should_generate_dense),spacing=1;
		if(want_generated_spacing>1)
			spacing = want_generated_spacing;
		if( should_generate_dense_nc !=0 && RSB_FABS(should_generate_dense_nc) < dim )
			dim = RSB_FABS(should_generate_dense_nc);
		dim *= spacing;

		if(((should_generate_lband>-1) || (should_generate_uband>-1)) && should_generate_dense>0)
		{
			const rsb_nnz_idx_t lbw = should_generate_lband , ubw = should_generate_uband;
			nrA = ncA = dim;
			if( should_generate_dense_nc !=0 && RSB_FABS(should_generate_dense_nc) < dim )
				ncA = RSB_FABS(should_generate_dense_nc);
			RSBENCH_STDOUT("# Generating a diagonally populated matrix of %zd x %zd\n",(rsb_printf_int_t)nrA,(rsb_printf_int_t)ncA);
			errval = rsb__generate_blocked_banded_coo(dim/spacing,spacing,lbw,ubw,&IA,&JA,&VA,&nnz,typecode);
			if(RSB_SOME_ERROR(errval))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
		}
		else
		{
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
		if(should_generate_dense>0)
		{
			RSBENCH_STDOUT("Interpreting --dense as --lower-dense (full dense makes no sense for triangular solve).\n");
			should_generate_dense = -should_generate_dense;
			should_generate_dense_nc = 0;
		}
')dnl
		if(should_generate_dense>0)
		{
			RSB_DEBUG_ASSERT( should_generate_dense_nc != 0 );
			/* full dense, no diag */
			nrA = dim;
			ncA = should_generate_dense_nc * spacing;
			errval = rsb__generate_dense_full(nrA/spacing,ncA/spacing,spacing,&IA,&JA,&VA,&nnz,typecode);
			if(RSB_SOME_ERROR(errval))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
		}
		else
		{
			/* trick: lower triangular */
			nrA=ncA=dim;
			errval = rsb__generate_dense_lower_triangular_coo(dim/spacing,spacing,&IA,&JA,&VA,&nnz,typecode);
			if(RSB_SOME_ERROR(errval))
			{RSB_ERROR(RSB_ERRM_ES);goto err;}
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER); /* 20121223	*/
		}
		}

		if(want_sort_after_load)	
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT);

		if(want_as_symmetric)
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_SYMMETRIC);
		if(want_as_hermitian)
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_HERMITIAN);
	} /* should_generate_dense */
dnl
have_va_ia_ja:
dnl
	RSB_DEBUG_ASSERT( ! ( VA == NULL && !b_r_filename ) );
	RSB_DEBUG_ASSERT( ! ( IA == NULL && !b_r_filename ) );
	RSB_DEBUG_ASSERT( ! ( JA == NULL && !b_r_filename ) );
	r_flags = flags;
dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
			flags = rsb__do_detect_and_add_triangular_flags(IA,JA,nnz,flags);
dnl			RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
			if(
		(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_TRIANGULAR) && RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_TRIANGULAR)) ||
		(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_TRIANGULAR)&&!RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_TRIANGULAR))
			)
			{
				RSB_ERROR("Matrix contains both upper and lower elements ? It is not suited for mop, then!\n");
				errval = RSB_ERR_CORRUPT_INPUT_DATA;	/* uhm */
				{RSB_ERROR(RSB_ERRM_ES);goto err;}
			}
')dnl
dnl

	/* CONDITIONALLY, PROCESSING THE INPUT */
	if(!b_r_filename)
	{
		if(want_column_expand)
		{
			errval = rsb__do_column_expand(JA,nnz,&ncA,want_column_expand);
			if(RSB_SOME_ERROR(errval))
			{
				RSB_ERROR(RSB_ERRM_ES);
				goto err;
			}
		}

		if( pattern_only )
			rsb__fill_with_ones(VA,typecode,nnz,1);

		if( dumpout )
		{
			errval = rsb__test_print_coo_mm(typecode,flags,IA,JA,VA,nrA,ncA,nnz,RSB_BOOL_TRUE,RSB_DEFAULT_STREAM);
			//COO equivalent for rsb_file_mtx_save(mtxAp,NULL);
			if(RSB_SOME_ERROR(errval))
			{
				RSB_ERROR(RSB_ERRM_ES);
				goto err;
			}
			goto ret;
		}
	}
#if 1
	if(want_nonzeroes_distplot)
	{
		/* FIXME: Unfinished: printout not adequate ! */
		/* FIXME: Shall use a separate routine for this! Please regard this code as temporary */
		rsb_coo_idx_t median_m=0,median_k=0,stdd_m=0,stdd_k=0,nzp_m=nnz/nrA,nzp_k=nnz/ncA;
		rsb_coo_idx_t*idxv=NULL;
		rsb_coo_idx_t mm=0;
		rsb_nnz_idx_t cs=0;
		rsb_bool_t po = RSB_BOOL_TRUE;
		const int histres=100;
		const rsb_char_t*const pmsg="\n\nplot \"-\" using 1:2 title \"cumulative %s population (nnz)\"\n";
		RSBENCH_STDOUT("set xtics rotate\n");
		RSBENCH_STDOUT("set term postscript eps color\n");
		RSBENCH_STDOUT("set output \"%s-distplot.eps\"\n", rsb__basename(filename));
		RSBENCH_STDOUT("set multiplot layout 1,2 title \"%s (%zd x %zd, %zd nnz)\"\n", rsb__basename(filename),(rsb_printf_int_t)nrA,(rsb_printf_int_t)ncA,(rsb_printf_int_t)nnz);

		ndA = RSB_MAX(nrA,ncA);
dnl		outnri = rhsnri = ndA;

		mm=nrA<histres?1:nrA/histres;
		idxv = rsb__calloc(sizeof(rsb_coo_idx_t)*(ndA));
		if(!idxv)
			goto nohists;

		for(i=0;i<nnz;++i)
			if(IA[i] < nrA && IA[i] >= 0 )
				idxv[IA[i]]++;
		for(i=0;i<nrA;++i)
			if(median_m<nnz/2)
				{ median_m+=idxv[i]; }
			else
				{ break; }
		median_m=i; 

		RSB_STDOUT(pmsg,"rows");
		if(po) for(i=0;i<nrA;++i){ cs+=idxv[i]; if(i%mm==0)RSB_STDOUT("%ld %ld\n",(long int)i,(long int)cs);}
		RSB_STDOUT("e\n");

		mm=ncA<histres?1:ncA/histres;

		for(i=0;i<nrA;++i)
			stdd_m+=(idxv[i]-nzp_m)*(idxv[i]-nzp_m);
		stdd_m=nrA<2?0:sqrt(stdd_m/(nrA-1));


		for(i=0;i<ncA;++i)
			idxv[i]=0;

		for(i=0;i<nnz;++i)
			if(JA[i] < ncA && JA[i] >= 0 )
				idxv[JA[i]]++;
		for(i=0;i<ncA;++i)
			if(median_k<nnz/2)
				{ median_k+=idxv[i]; }
			else
				{ break; }
		median_k=i; 

		cs=0;
		RSB_STDOUT(pmsg,"columns");
		if(po) for(i=0;i<ncA;++i){ cs+=idxv[i]; if(i%mm==0)RSB_STDOUT("%ld %ld\n",(long int)i,(long int)cs);}
		RSB_STDOUT("e\n");

		for(i=0;i<ncA;++i)
			stdd_k+=(idxv[i]-nzp_k)*(idxv[i]-nzp_k);
		stdd_k=ncA<2?0:sqrt(stdd_k/(ncA-1));

		RSBENCH_STDOUT("unset multiplot\n");
		RSBENCH_STDOUT("#%%:NNZ_PER_ROW_STDDEV:");/* RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(); */
		RSBENCH_STDOUT("\t%10.0zd\n",(rsb_printf_int_t)stdd_m);
		RSBENCH_STDOUT("#%%:ROWS_MEDIAN:");/* RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(); */
		RSBENCH_STDOUT("\t%10.0g\n",((double)median_m/(double)nrA));
		RSBENCH_STDOUT("#%%:NNZ_PER_COL_STDDEV:");/* RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(); */
		RSBENCH_STDOUT("\t%10.0zd\n",(rsb_printf_int_t)stdd_k);
		RSBENCH_STDOUT("#%%:COLS_MEDIAN:");/* RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(); */
		RSBENCH_STDOUT("\t%10.0g\n",((double)median_k/(double)ncA));
nohists:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
		RSB_CONDITIONAL_FREE(idxv); RSB_CONDITIONAL_FREE(idxv);
		goto ret;
	}
	#endif /* 1 */
dnl
dnl
	/* CONDITIONALLY, PERFORMING SOME TEST ON THE INPUT */
	if(want_accuracy_test>=1)
	{
		struct rsb_coo_mtx_t coo;
		rsb__fill_coo_struct(&coo,VA,IA,JA,nrA,ncA,nnz,typecode);
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_accuracy_test(&coo,ca,cn,flags));
		if(RSB_SOME_ERROR(errval))
		{
			RSB_ERROR("accuracy based test failed!\n");
			goto err;
		}
		if(want_accuracy_test>1)
		{
			goto done;
		}
	}

		if( (flags & RSB_FLAG_QUAD_PARTITIONING) && g_all_flags==1)
		{
dnl			int di=0,hi=0,li=0,oi=0;
			int hi=0,oi=0;
			fn=0;
			for(rsb_thread_t ci=0;ci<3;++ci)
/*			for(di=0;di<2;++di)*/
			for(oi=0;oi<2;++oi)
			for(hi=0;hi<2;++hi)
/*			for(li=0;li<2;++li)*/
			{
#if 0
				flagsa[di+hi*2+li*4+ci*8]=flags;
				//RSB_DO_FLAG_ADD(flagsa[di+hi*2+li*4+ci*8],hi?RSB_FLAG_USE_HALFWORD_INDICES_COO:0);
				RSB_DO_FLAG_ADD(flagsa[di+hi*2+li*4+ci*8],hi?RSB_FLAG_USE_HALFWORD_INDICES_CSR:0);
#ifdef RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES
				RSB_DO_FLAG_ADD(flagsa[di+hi*2+li*4+ci*8],li?RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES:0);
#endif /* RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES */
				RSB_DO_FLAG_ADD(flagsa[di+hi*2+li*4+ci*8],di?RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG:0);
	
#if 0
				RSB_DO_FLAG_ADD(flagsa[di+hi*2+li*4+ci*8],ci==1?RSB_FLAG_RECURSIVE_HALF_DETECTED_CACHE:0);
				RSB_DO_FLAG_ADD(flagsa[di+hi*2+li*4+ci*8],ci==2?RSB_FLAG_RECURSIVE_DOUBLE_DETECTED_CACHE:0);
#endif /* 0 */
#else /* 0 */
				flagsa[fn]=flags;
				//RSB_DO_FLAG_ADD(flagsa[fn],li?RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES:0);
				//RSB_DO_FLAG_ADD(flagsa[fn],di?RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG:0);
				RSB_DO_FLAG_ADD(flagsa[fn],oi?RSB_FLAG_USE_HALFWORD_INDICES_COO:0);
				RSB_DO_FLAG_ADD(flagsa[fn],hi?RSB_FLAG_USE_HALFWORD_INDICES_CSR:0);
#if 0
				RSB_DO_FLAG_ADD(flagsa[fn],ci==1?RSB_FLAG_RECURSIVE_HALF_DETECTED_CACHE:0);
				RSB_DO_FLAG_ADD(flagsa[fn],ci==2?RSB_FLAG_RECURSIVE_DOUBLE_DETECTED_CACHE:0);
#endif /* 0 */
#endif /* 0 */
				++fn;
			}
		}
		else
		{
			fn=1;
			flagsa[fn-1]=flags;
		}

		if(!want_perf_dump)
		if(!( RSB_MKL_APPROPRIATE_AT_TIME_SPEC( want_autotuner ) || RSB_MKL_APPROPRIATE_AT_TIME_SPEC( merge_experimental ) || RSB_MKL_APPROPRIATE_AT_TIME_SPEC( split_experimental ) )) /* otherwise pr__set.. cannot distinguish samples */
		if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
		{
			/* adds a no-recursion flag case */
			RSB_DO_FLAG_DEL(flags,RSB_FLAG_QUAD_PARTITIONING);
/*			if(fn)*/
/*				flags=flagsa[fn-1];	*//* copy from the last */
/*			else*/
/*				flagsa[fn]=flags;	*//* impose these flags */
			for(fi=fn;fi>0;--fi)
				flagsa[fi]=flagsa[fi-1];/* shift forward */
			RSB_DO_FLAG_DEL(flagsa[0],RSB_FLAG_QUAD_PARTITIONING);
			++fn;	/* add ours */
		}

dnl
	for(     nrhsi=0;     nrhsi<nrhsn     ;++nrhsi     )
	{
	nrhs = nrhsa[nrhsi];
	if(nrhs<1)
	{
		RSBENCH_STDOUT("# WARNING: Skipping non-positive nrhs (%zd): is this a mistake ?\n",(rsb_printf_int_t)nrhs );
		continue;
	}
	if( nrhsn > 1 && nrhss )
	{
		RSBENCH_STDOUT("# multi-nrhs benchmarking (%s) -- now using nrhs %zd.\n",nrhss,(rsb_printf_int_t)nrhs);
	}
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	for(     alphai=0;     alphai<alphan;++alphai )
	{
	for(     betai=0;     betai<betan;++betai )
	{
')dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	for(     incXi=0;     incXi<incXn     ;++incXi     )
	{
	for(     incYi=0;     incYi<incYn     ;++incYi     )
	{
')dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	incX = incXa[incXi];
	incY = incYa[incYi];
	if(want_verbose >= -2)
	{
		if(incXn>1)
		{
			RSBENCH_STDOUT("# multi-incX benchmarking (%zd/%zd) -- now using incX=%zd.\n",(rsb_printf_int_t)(incXi+1),(rsb_printf_int_t)incXn,(rsb_printf_int_t)incX);
		}
		if(incYn>1)
		{
			RSBENCH_STDOUT("# multi-incY benchmarking (%zd/%zd) -- now using incY=%zd.\n",(rsb_printf_int_t)(incYi+1),(rsb_printf_int_t)incYn,(rsb_printf_int_t)incY);
		}
	}
',`dnl
dnl	
')dnl
	if(incX<1)
	{
		RSBENCH_STDOUT("# WARNING: Skipping non-positive incX (%d): is this a mistake ?\n",(int)incX );
		continue;
	}
	if(incY<1)
	{
		RSBENCH_STDOUT("# WARNING: Skipping non-positive incY (%d): is this a mistake ?\n",(int)incY );
		continue;
	}
	if( want_only_star_scan )
	if( RSB_MIN(incXi,1) + RSB_MIN(incYi,1) + RSB_MIN(nrhsi,1) > 1 ) /* two or more exceed index one */
	{
		RSBENCH_STDOUT("# Skipping a case with incX=%zd incY=%zd nrhs=%zd.\n",(rsb_printf_int_t)incX,(rsb_printf_int_t)incY,(rsb_printf_int_t)nrhs);
		goto frv;
	}

ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
dnl	if(flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)	/* this is here only for easing triangular solve benchmarking (to avoid each time ) */
dnl		transA = RSB_TRANSPOSITION_T;
dnl	else
dnl		transA = RSB_TRANSPOSITION_N;
dnl	20110412	vectors shall be transA-independent---now both transposed and untransposed operations could be executed in the same run
	if(incX!=incY)
	{
		RSB_ERROR("setting (incX=%d) != (incY=%d) in triangular solve is unsupported in this program\n",incX,incY);
		errval = RSB_ERR_BADARGS;goto err;
	}
')dnl

ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
ifelse(RSB_M4_IS_ZEROING_KERNEL_MOP(mop),1,`dnl
	if(RSB_SOME_ERROR(errval = rsb__cblas_Xscal(typecode,1,NULL,beta,1))){ RSB_ERROR(RSB_ERRM_ES);goto err;}
',`
	if(RSB_SOME_ERROR(errval = rsb__fill_with_ones(beta,typecode,1,1))){ RSB_ERROR(RSB_ERRM_ES);goto err;}
')dnl
ifelse(RSB_M4_IS_OP_ADDING_KERNEL_MOP(mop),1,`dnl
	if(RSB_SOME_ERROR(errval = rsb__fill_with_ones(alpha,typecode,1,1))){ RSB_ERROR(RSB_ERRM_ES);goto err;}
')dnl
ifelse(RSB_M4_IS_OP_NEGATING_KERNEL_MOP(mop),1,`dnl
	if(RSB_SOME_ERROR(errval = rsb__fill_with_ones(alpha,typecode,1,1))){ RSB_ERROR(RSB_ERRM_ES);goto err;}
')dnl
	/* FIXME: the following collides with the former */
	rsb__util_set_area_to_converted_integer(alphap,typecode,alphai);
	rsb__util_set_area_to_converted_integer(betap ,typecode,betai);
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	/* FIXME: the following overwrites what above */
	rsb__util_set_area_to_converted_integer(alphap,typecode,alphaip[alphai]);
	rsb__util_set_area_to_converted_integer(betap ,typecode, betaip[betai]);
	if(want_verbose >= -1)
		RSBENCH_STDOUT("# Using alpha=%d beta=%d order=%s for rsb_spmv/rsb_spsv/rsb_spmm/rsb_spsm.\n",alphaip[alphai],betaip[betai],((order==RSB_FLAG_WANT_ROW_MAJOR_ORDER)?"rows":"cols"));
')
dnl
		for(ti=0;ti<tn;++ti)
dnl		if(!((ti>=1)&&(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC)||RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN))))
		{
ifelse(mop,`mat_stats',`',`dnl
	rsb_time_t op_t = RSB_TIME_ZERO;
	rsb_time_t mct = RSB_TIME_ZERO;	/* matrix construction time */
	rsb_time_t fet = RSB_TIME_ZERO;	/* fillin estimation time */

	rsb_time_t sct = RSB_TIME_ZERO;	/* serial (if minimum number of cores is 1) matrix construction time */
	rsb_time_t pct = RSB_TIME_ZERO;	/* parallel (if maximum number of cores > 1) matrix construction time */

ifelse(mop,`infty_norm',`',`dnl
ifelse(mop,`rowssums',`',`dnl
ifelse(mop,`scale',`',`dnl
	rsb_time_t smt = RSB_TIME_ZERO;	/* serial multiplication time */
	rsb_time_t pmt = RSB_TIME_ZERO;	/* parallel multiplication time */
')dnl
')dnl
')dnl
	const rsb_int_t mintimes = RSB_CONST_AT_OP_SAMPLES_MIN, maxtimes = RSB_CONST_AT_OP_SAMPLES_MAX;

ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
	rsb_time_t sst = RSB_TIME_ZERO;	/* serial solve time */
	rsb_time_t pst = RSB_TIME_ZERO;	/* parallel solve time */
')dnl
	
	rsb_time_t sest = RSB_TIME_ZERO;	/**/
	//rsb_time_t sect = RSB_TIME_ZERO;	/**/
	rsb_time_t ssat = RSB_TIME_ZERO;	/**/
	rsb_time_t seit = RSB_TIME_ZERO;	/**/
	rsb_time_t scpt = RSB_TIME_ZERO;	/**/

	rsb_time_t mest = RSB_TIME_ZERO;	/**/
	rsb_time_t mect = RSB_TIME_ZERO;	/**/
	rsb_time_t msat = RSB_TIME_ZERO;	/**/
	rsb_time_t meit = RSB_TIME_ZERO;	/**/
	rsb_time_t mcpt = RSB_TIME_ZERO;	/**/

	rsb_time_t me_best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME, me_at_best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;     /* experimental merge */
	rsb_time_t at_best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME, at_mkl_csr_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME; /* experimental merge */
	rsb_thread_t at_mkl_csr_nt = RSB_AT_THREADS_AUTO, me_at_nt = RSB_AT_THREADS_AUTO;
#if RSB_EXPERIMENTAL_WANT_BEST_TIMES
	rsb_time_t best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_time_t base_best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;	/* for comparative benchmarking */
	rsb_time_t serial_best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;	/* for comparative benchmarking */
	rsb_time_t spmv_t = RSB_TIME_ZERO;
	rsb_time_t tot_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
	rsb_time_t spsv_d_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
	rsb_time_t spsv_spmv_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
	rsb_time_t best_spsv_spmv_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_time_t spsv_f_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
#endif
#ifdef RSB_WANT_PERFORMANCE_COUNTERS
	struct rsb_pci_t rsb_pci;
#endif /* RSB_WANT_PERFORMANCE_COUNTERS */
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
	void *M_VA=NULL; MKL_INT *M_IA=NULL,*M_JA=NULL;
	void *M_VAC=NULL; MKL_INT *M_IAC=NULL,*M_JAC=NULL;
	rsb_time_t mkl_coo2csr_time = RSB_TIME_ZERO;
	rsb_time_t mkl_coo_op_tot_time = RSB_TIME_ZERO;
	rsb_time_t mkl_coo_op_time = RSB_TIME_ZERO;
	rsb_time_t mkl_coo_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_time_t mkl_coo_op_time_best_serial = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_time_t mkl_csr_op_tot_time = RSB_TIME_ZERO;
	rsb_time_t mkl_csr_op_time = RSB_TIME_ZERO;
	rsb_time_t mkl_csr_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_time_t mkl_csr_op_time_best_serial = RSB_CONST_IMPOSSIBLY_BIG_TIME;

	rsb_time_t mkl_gem_op_tot_time = RSB_TIME_ZERO;
	rsb_time_t mkl_gem_op_time = RSB_TIME_ZERO;
	rsb_time_t mkl_gem_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_time_t mkl_gem_op_time_best_serial = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	struct rsb_ts_t btpms[2]; /* first is tuned, second is not */
	const rsb_flags_t mif = ( mib == 0 ) ? RSB_FLAG_NOFLAGS : RSB_FLAG_FORTRAN_INDICES_INTERFACE; /* MKL index flags */
#ifdef RSB_WANT_PERFORMANCE_COUNTERS
	struct rsb_pci_t mkl_coo_pci,mkl_csr_pci,mkl_gem_pci;
#endif /* RSB_WANT_PERFORMANCE_COUNTERS */
#endif /* RSB_WANT_MKL */
	struct rsb_attr_t attr;	/* this structure is rather large (100k, as of 20140223); with future parameters it shall be rather heap allocated */
	struct rsb_ts_t otpos, btpos;

	RSB_BZERO_P((&otpos));
	RSB_BZERO_P((&btpos));
')dnl
dnl
')dnl
dnl

ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
	RSB_BZERO_P((&attr));
')dnl

ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	if(want_unordered_coo_bench)
	{
		struct rsb_coo_mtx_t coo;
		int ctimes = times < 1 ? 100 : times;
		rsb_time_t unordered_coo_op_tot_time = RSB_CONST_IMPOSSIBLY_BIG_TIME, unordered_coo_op_time = RSB_CONST_IMPOSSIBLY_BIG_TIME, unordered_coo_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;

		rsb__fill_coo_struct(&coo,VA,IA,JA,nrA,ncA,nnz,typecode);
		ndA = RSB_MAX(nrA,ncA);
dnl		outnri = rhsnri = ndA;
		lhs = rsb__calloc_vector(ndA*nrhs*incY,typecode);
		rhs = rsb__calloc_vector(ndA*nrhs*incX,typecode);

		if(!lhs || !rhs)
		{
			RSB_ERROR("problems allocating vectors");
			RSB_CONDITIONAL_FREE(lhs); RSB_CONDITIONAL_FREE(rhs);
			{ errval = RSB_ERR_INTERNAL_ERROR; goto err; }
		}

		if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
		for(i=0;i<ctimes;++i)
		{
			if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
			unordered_coo_op_time = - rsb_time();
			if((errval = rsb__do_spmv_fullword_coo(&coo,flags,rhs,lhs,alphap,betap,incX,incY,transA))!=RSB_ERR_NO_ERROR) { goto erru; }
			unordered_coo_op_time += rsb_time();
			unordered_coo_op_time_best = RSB_MIN_ABOVE_INF(unordered_coo_op_time_best,unordered_coo_op_time,tinf);
			unordered_coo_op_tot_time+=unordered_coo_op_time;
		}
		if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
erru:
		RSB_CONDITIONAL_FREE(lhs); RSB_CONDITIONAL_FREE(rhs);
dnl
		if(want_verbose > 0)
		{
			struct rsb_mtx_t matrixs;

			mtxAp=&matrixs;
			rsb__init_rsb_struct_from_coo(mtxAp,&coo);
			mtxAp->flags = RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS|RSB_DO_FLAG_FILTEROUT((flags),RSB_DO_FLAGS_EXTRACT_STORAGE(flags));
			rsb__do_set_init_storage_flags(mtxAp,mtxAp->flags);
			raw_Mflops=nnz*2;
			RSBENCH_STDOUT("%%:UNORDERED_COO_PERFORMANCE:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
			RSBENCH_STDOUT("\t%10.2lf\n",((rsb_time_t)raw_Mflops)/(RSB_REAL_MILLION*unordered_coo_op_time_best));
			mtxAp=NULL;
		}
	}
dnl
')dnl

dnl
		transA = transAo;
dnl
		if(ti>0)
			transA = rsb__do_transpose_transposition(transAo);
dnl
		if(ti==2)
			transA = RSB_TRANSPOSITION_C;
dnl
		if(!  (
			( RSB_IS_MATRIX_TYPE_COMPLEX(typecode) && (ti!=0) && ( flags & RSB_FLAG_ANY_SYMMETRY ) )  ||
		       ((!RSB_IS_MATRIX_TYPE_COMPLEX(typecode))&& (ti!=0) && ( flags & RSB_FLAG_SYMMETRIC) )  ||
		       ((!RSB_IS_MATRIX_TYPE_COMPLEX(typecode))&& (ti==2) &&!( flags & RSB_FLAG_ANY_SYMMETRY) )  ||
			g_allow_any_tr_comb
		))
dnl
		if(tn>1)
		{
			RSBENCH_STDOUT("# multi-transpose benchmarking -- now using transA = %c.\n",RSB_TRANSPOSITION_AS_CHAR(transA));
		}
dnl
		if( /* transA != RSB_TRANSPOSITION_N */ ti>0 && RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC) )
		{
			RSBENCH_STDOUT("# symmetric matrix --- skipping transposed benchmarking\n");
			continue;
		}
dnl
		if(want_verbose > 0)
		{
			RSBENCH_STDOUT("# will use input matrix flags: ");
			rsb__dump_flags(flags,"",", ","\n");
		}
dnl
		for(fi=0;fi<fn;++fi)
		for(brvi=-1;brvi<brl;++brvi)
		for(bcvi=-1;bcvi<bcl;++bcvi)
#ifndef  RSB_COORDINATE_TYPE_H
		if(!(flagsa[fi] & RSB_FLAG_USE_HALFWORD_INDICES_CSR))
#endif /* RSB_COORDINATE_TYPE_H */
		for(rsb_thread_t ci=0;ci<cn;++ci)	/* here just for should_recycle_matrix */
		if(!(cn > 1 && ca[ci]>1 && !(RSB_DO_FLAG_HAS(flagsa[fi],RSB_FLAG_QUAD_PARTITIONING)))) /* no need for more than one core without recursion */
		{
ifelse(mop,`mat_stats',`',`dnl
			rsb_time_t diag_op_tot_time = RSB_TIME_ZERO;
			rsb_time_t diag_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			rsb_time_t getrow_op_tot_time = RSB_TIME_ZERO;
			rsb_time_t getrow_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
')dnl
dnl
ifelse(mop,`mat_stats',`',`dnl
			rsb_time_t diag_op_time_best_serial = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			rsb_time_t getrow_op_time_best_serial = RSB_CONST_IMPOSSIBLY_BIG_TIME;
')dnl
dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
dnl
			rsb_time_t no_lock_op_time = RSB_CONST_IMPOSSIBLY_BIG_TIME, no_lock_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME,
			serial_no_lock_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME, no_lock_op_tot_time = RSB_TIME_ZERO;
dnl
			rsb_time_t qt_op_time = RSB_CONST_IMPOSSIBLY_BIG_TIME, qt_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME,
			qt_op_tot_time = RSB_TIME_ZERO;
dnl
')dnl
			cc = ca[ci];
dnl
			should_recycle_matrix=(ci>0)?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
			/* if this is the special "vanilla CSR" run after/before recursive runs ... */
			if(rsb__set_num_threads(cc)!=cc)
			{
				RSB_ERROR("failed setting %d threads!\n",cc);
				errval = RSB_ERR_INTERNAL_ERROR;
				goto err;
			}
			else
			{
				const rsb_thread_t rtn = rsb__set_num_threads(RSB_THREADS_GET);
				RSBENCH_STDOUT("# Using %ld threads\n",(long)rtn );
				RSB_DEBUG_ASSERT( RSB_CONST_MAX_SUPPORTED_THREADS >= rtn );
			}
			flags=flagsa[fi];
			if(cn>1 && !RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
				RSB_DO_FLAG_DEL(flags,RSB_FLAG_USE_HALFWORD_INDICES);

ifelse(mop,`mat_stats',`',`dnl
			best_spsv_spmv_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			op_t = RSB_TIME_ZERO;
#if RSB_EXPERIMENTAL_WANT_BEST_TIMES
			best_t = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			spmv_t = RSB_TIME_ZERO;
			tot_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
			spsv_d_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
			spsv_spmv_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
			spsv_f_t = RSB_TIME_ZERO;	/* cumulative time (not best one)*/
#endif /* RSB_EXPERIMENTAL_WANT_BEST_TIMES */
')dnl

			if(brl>0 && bcl>0)
			{
				/* this is a trick and an unclean programming practice */
				if(brvi==-1)++brvi;
				if(bcvi==-1)++bcvi;
				br = brv[brvi];
				bc = bcv[bcvi];
			}
			else
			{	
				/* br, bc already set */
			}

#if RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS 
			/*	
			* FIXME : laziness
			*/
			dnl RSB_WARN("using RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS\n");
			if( br!=1 || bc!=1 || !rsb__util_are_flags_suitable_for_optimized_1x1_constructor(flags) )
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
ifelse(mop,`mat_stats',`dnl
			if(1) /* only in mat_stats */
',`dnl
#if RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT
			if(0)
#endif /* RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT */
')dnl
			{
				p_r = rsb__util_get_partitioning_array(br,nrA,&M_b,flags);
				p_c = rsb__util_get_partitioning_array(bc,ncA,&K_b,flags);

				if((! p_r) || (! p_c))
				{
					errval = RSB_ERR_ENOMEM;
					RSB_PERR_GOTO(erri,RSB_ERRM_ES);
				}
			}

			if(  ( br!=1 || bc!=1 || p_r || p_c ) && ( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR ))
			{
				/*  */
				RSB_WARN("WARNING : disabling in place allocation flag : it is only allowed for 1x1!\n");
				RSB_DO_FLAG_DEL(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR) ;
			}

ifelse(mop,`mat_stats',`dnl
			pinfo.M_b=M_b;
			pinfo.K_b=K_b;
			pinfo.rpntr=p_r;
			pinfo.cpntr=p_c;
')dnl


ifelse(mop,`mat_stats',`dnl
			if(max_nnzs==0)
				max_nnzs=nnz;
	if(until_confidence && g_estimate_fillin)
	{
		if( want_percentage && ( max_nnzs > 100 || max_nnzs < 1) ) 
		{RSBENCH_STDERR("given percentage = %zd ?\n",(rsb_printf_int_t)max_nnzs);goto err;}
		else
		{
			if( want_percentage ) max_nnzs =(rsb_nnz_idx_t ) (((double)nnz/100.0) *(double) max_nnzs );

			if(max_nnzs>nnz)
			{RSBENCH_STDERR("want more max_nnzs (%zd) than nonzeros (%zd) !\n",(rsb_printf_int_t)max_nnzs,(rsb_printf_int_t)nnz);goto err;}
			else
			if(max_nnzs<nnzn)
			{RSBENCH_STDERR("want max_nnzs (%zd) less than %zd ?\n",(rsb_printf_int_t)max_nnzs,(rsb_printf_int_t)nnzn);goto err;}
		}
	}

#if 0
	if(!until_confidence && !g_estimate_fillin)
	{
		{RSBENCH_STDERR("should choose an option : [ -S points] (-e)!\n");goto err;}
		goto err;
	}
#else /* 0 */
	g_estimate_fillin=1;
#endif /* 0 */
		if( until_confidence && ( until_confidence > 100 || until_confidence < 1) ) 
		{RSBENCH_STDERR("given percentage = %zd ?\n",(rsb_printf_int_t)until_confidence ); {RSB_ERROR(RSB_ERRM_ES);goto err;} ;}

			if(g_estimate_fillin)
			{
				size_t total_element_count=0;
				size_t total_block_count=0;
				rsb_fillin_t fillin;

				nnzs = rsb__calloc(nnzn * sizeof(size_t));
				element_count = rsb__calloc(nnzn * sizeof(size_t));
				block_count = rsb__calloc(nnzn * sizeof(size_t));

				if(!nnzs || !element_count || !block_count)
				{
					errval = RSB_ERR_ENOMEM;
					RSB_PERR_GOTO(erri,RSB_ERRM_ES);
				}

				for(i=1;i<=nnzn;++i) nnzs[i-1]=(max_nnzs/nnzn) * i;/* ach, integer arithmetics ! */
				nnzs[nnzn-1]=max_nnzs;
				nnzs[nnzn-1]=nnz;
	
				errval = rsb__compute_partial_fillin_for_nnz_fractions(IA, JA, nnzs, nnzn, &pinfo, element_count, block_count);
				if(RSB_SOME_ERROR(errval))
				{
					RSB_PERR_GOTO(erri,RSB_ERRM_ES);
				}

				errval = rsb__compute_partial_fillin_for_nnz_fractions(IA, JA, &nnz, 1, &pinfo, &total_element_count, &total_block_count);
				if(RSB_SOME_ERROR(errval))
				{
					RSB_PERR_GOTO(erri,RSB_ERRM_ES);
				}
				fillin = ((double)total_element_count)/((double)nnz);
	
				//RSB_STDOUT("#using %d up to %d nonzeros out of %d, we estimate the fillin as:\n",nnzs[0],nnzs[nnzn-1],nnz);
				RSBENCH_STDOUT("#matrix	rows	cols	br	bc	nnz	fillin	fraction	rel.error\n");
				for(i=0;i< nnzn;++i)
				{
					rsb_fillin_t partial_fillin=0;
/*					RSBENCH_STDOUT("#%d\n",nnzs[i]);*/
/*					RSBENCH_STDOUT("#%d / %d\n",element_count[i],total_element_count);*/
					RSBENCH_STDOUT("%s\t%zd\t%zd\t%zd\t%zd\t%zd\t%lg",filename,
					(rsb_printf_int_t)nrA,(rsb_printf_int_t)ncA,(rsb_printf_int_t)br,(rsb_printf_int_t)bc,(rsb_printf_int_t)nnz,fillin);
					//RSBENCH_STDOUT(" (%d,%d)",element_count[i],block_count[i]);
					partial_fillin = (element_count[i])/(double)(nnzs[i]);
					RSBENCH_STDOUT("\t%.3lg\t%+.3lg\n",
						((double)nnzs[i])/(double)nnz,
						(partial_fillin-fillin)/fillin
					);
				}
				//RSBENCH_STDOUT("\n");
			}


',`dnl

ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
#define RSB_WANT_SPSV_STABILITY_FIX 1
#if RSB_WANT_SPSV_STABILITY_FIX
#if 0
			/* FIXME : fix for numerical stability */
#if 0
			if(RSB_SOME_ERROR(rsb__fill_with_ones(VA,typecode,nnz,1))) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
#else /* 0 */
			/* FIXME : temporary fix */
			const double uthreshold=.0001;
			const double athreshold=10000000;
			if(RSB_SOME_ERROR(rsb__util_drop_to_zero_if_under_threshold(VA,typecode,nnz,&uthreshold))) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
			if(RSB_SOME_ERROR(rsb__util_drop_to_zero_if_above_threshold(VA,typecode,nnz,&athreshold))) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
#endif /* 0 */
#else /* 0 */
			{rsb_nnz_idx_t n;for(n=0;n<nnz;++n)if(IA[n]==JA[n])rsb__fill_with_ones(((rsb_byte_t*)VA)+RSB_SIZEOF(typecode)*n,typecode,1,1);}
#endif /* 0 */
#endif /* RSB_WANT_SPSV_STABILITY_FIX */
')dnl

			if(!mtxAp)
			{
				int mci=0;
				if(b_r_filename)
				{
					rsb_err_t errval_;
					mct = - rsb_time();
					mtxAp = rsb__load_matrix_file_as_binary(b_r_filename,&errval_);
					mct += rsb_time();
					if((RSB_SOME_ERROR(errval)) || !mtxAp )
					{
						RSB_PERR_GOTO(err,RSB_ERRM_ES);
					}
					else
					{
						nnz = mtxAp->nnz;
						nrA = mtxAp->nr;
						ncA = mtxAp->nc;
					}

					filename=b_r_filename;// for info purposes
					flags=mtxAp->flags;
				}
				else
				{
				mect=mest=msat=meit=mcpt = RSB_TIME_ZERO;	/* resetting al values */

				for(mci=0;mci<repeat_construction;++mci)
				{
					if(repeat_construction>1 && mci==0)
						RSBENCH_STDOUT("# will repeat constructor %d times\n",repeat_construction);
					mct = - rsb_time();
					if(RSB_UNLIKELY(want_in_place_assembly))
						mtxAp = rsb__do_mtx_alloc_from_coo_inplace(VA,IA,JA,nnz,typecode,nrA,ncA,br,bc,flags,&errval);
					else
						mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,nrA,ncA,br,bc,flags,&errval);
					mct += rsb_time();
					if((RSB_SOME_ERROR(errval)) || !mtxAp )
					{
						RSB_PERR_GOTO(err,RSB_ERRM_MBE);
					}

					if(mect == RSB_TIME_ZERO || mect>mtxAp->ect)
						mect=mtxAp->est;
					if(mest == RSB_TIME_ZERO || mest>mtxAp->est)
						mest=mtxAp->est;
					if(msat == RSB_TIME_ZERO || msat>mtxAp->sat)
						msat=mtxAp->sat;
					if(meit == RSB_TIME_ZERO || meit>mtxAp->eit)
						meit=mtxAp->eit;
					if(mcpt == RSB_TIME_ZERO || mcpt>mtxAp->cpt)
						mcpt=mtxAp->cpt;
					if(mci != repeat_construction-1)
					{ RSB_MTX_FREE(mtxAp);	/* we only wanted timings */ }
					else
					{
						/* we keep the mtxAp, and set best individual times */;
						mtxAp->est=mest;
						mtxAp->ect=mect;
						mtxAp->sat=msat;
						mtxAp->eit=meit;
						mtxAp->cpt=mcpt;
					}
				}
				}
				if(ci==0 && sct == RSB_TIME_ZERO)
					//sct=mct;
					sct=mtxAp->tat;
				if(ci==cn-1 && pct == RSB_TIME_ZERO)
					//pct=mct;
					pct=mtxAp->tat;
			} /* !mtxAp */
			
			if(do_perform_ddc == RSB_BOOL_TRUE)
			{
			if(rsb__is_square(mtxAp))
			{
				/* FIXME: experimental, new. should write a test with octave for this */
				void * DV = rsb__calloc_vector(mtxAp->nr,mtxAp->typecode);
				void * RS = rsb__calloc_vector(mtxAp->nr,mtxAp->typecode);
				rsb_aligned_t mtwo[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
				if(!RS||!DV) { errval = RSB_ERR_ENOMEM; goto noddc; }
				RSB_DO_ERROR_CUMULATE(errval,rsb_mtx_get_nrm(mtxAp, RS, RSB_EXTF_NORM_INF));
				rsb__util_set_area_to_converted_integer(mtwo,mtxAp->typecode,-2);
				RSB_DO_ERROR_CUMULATE(errval,rsb__dodo_getdiag(mtxAp,DV));
				RSB_DO_ERROR_CUMULATE(errval,rsb__vector_to_abs(DV,mtxAp->typecode,mtxAp->nr));
				RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xscal(mtxAp->typecode,mtxAp->nr,mtwo,DV,1));
				RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xaxpy(mtxAp->typecode,mtxAp->nr,NULL,DV,1,RS,1));
				if(rsb__util_count_negative(RS,mtxAp->typecode,mtxAp->nr)==mtxAp->nr)
					RSBENCH_STDOUT("#matrix is diagonal dominant\n");
				else
					RSBENCH_STDOUT("#matrix is not diagonal dominant\n");
				RSBENCH_STDOUT("#diagonal dominance computed in ? s\n");
noddc:
				RSB_CONDITIONAL_FREE(DV); RSB_CONDITIONAL_FREE(RS);
				if(RSB_SOME_ERROR(errval))
					goto err;
			}
			else
			{
				RSB_ERROR("input matrix is not square: cannot compute the diagonal dominance check\n");
			}
			}

			if( dump_graph_file )
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_print_matrix_stats(mtxAp,RSB_CONST_DUMP_DOT,dump_graph_file));

			if(do_perform_ilu == RSB_BOOL_TRUE)
			{
				/* FIXME: experimental */
				rsb_time_t ilut;
{
				// TODO: switch to coo, merge, and csr.
				while(rsb__submatrices(mtxAp)>1)
				{
					errval = rsb__leaves_merge(mtxAp, 0, NULL, NULL, NULL, 0, 0);
					if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES); }
				}
				errval = rsb__do_switch_recursive_matrix_to_fullword_storage(mtxAp);
				if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES); }
}
				ilut = - rsb_time();
				RSB_STDOUT("performing EXPERIMENTAL ILU-0\n");
				errval = rsb__prec_ilu0(mtxAp);//TODO: actually, only for CSR
				ilut += rsb_time();
				if(RSB_SOME_ERROR(errval))
				{
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
				}
				else
					RSB_STDOUT("performed EXPERIMENTAL ILU-0 with success in %lg s.\n",ilut);
				rsb_file_mtx_save(mtxAp,NULL);
dnl				goto ret;
			} /* do_perform_ilu */

			if(want_update && mtxAp)
			{
				rsb_time_t ct = - rsb_time();
				/* FIXME: this is update, not conversion, so it should not be here */

				errval = rsb__do_set_coo_elements(mtxAp,VA,IA,JA,nnz);
				if(RSB_SOME_ERROR(errval))
				{
					RSB_PERR_GOTO(erri,RSB_ERRM_ES);
				}
				ct += rsb_time();
				/* missing check */
				RSBENCH_STDOUT("#individual update of %zd elements in assembled RSB took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)nnz,ct,(100*ct)/mtxAp->tat);
			} /* want_update */

			if(want_convert && mtxAp)
			{
				/* benchmark conversions */
				rsb_time_t ct;
				rsb_nnz_idx_t rnz=0;
				struct rsb_coo_mtx_t coo;
				struct rsb_coo_mtx_t cor;

				coo.nnz = RSB_MAX(mtxAp->nnz,RSB_MAX(nrA,ncA));
				coo.typecode=mtxAp->typecode;
				if(rsb__allocate_coo_matrix_t(&coo)!=&coo)
				{
					RSB_ERROR(RSB_ERRM_ES);
					goto errc;
				}
				coo.nr = mtxAp->nr;
				coo.nc = mtxAp->nc;

				cor.nnz = RSB_MAX(mtxAp->nnz,RSB_MAX(nrA,ncA));
				cor.typecode=mtxAp->typecode;
				if(rsb__allocate_coo_matrix_t(&cor)!=&cor)
				{
					RSB_ERROR(RSB_ERRM_ES);
					goto errc;
				}
				cor.nr = mtxAp->nr;
				cor.nc = mtxAp->nc;

				ct = - rsb_time();
				errval = rsb__do_get_rows_sparse(RSB_TRANSPOSITION_N,NULL,mtxAp,coo.VA,coo.IA,coo.JA,0,mtxAp->nr-1,&rnz,RSB_FLAG_NOFLAGS);
				if(RSB_SOME_ERROR(errval))
				{ RSB_ERROR(RSB_ERRM_ES);goto erri;}
				ct += rsb_time();
				if(RSB_SOME_ERROR(rsb__util_is_sorted_coo_as_row_major(coo.IA,coo.JA,coo.nnz,coo.typecode,
					NULL,RSB_FLAG_NOFLAGS)))
					{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
				RSBENCH_STDOUT("#extraction of %zd elements in sorted COO took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);
				RSBENCH_STDOUT("#extraction to unsorted COO unimplemented\n");
				//RSBENCH_STDOUT("#extraction of %zd elements in unsorted COO took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);

				RSB_DO_ERROR_CUMULATE(errval,rsb_mtx_get_coo(mtxAp,cor.VA,cor.IA,cor.JA,RSB_FLAG_C_INDICES_INTERFACE));

				rsb__util_coo_array_set(coo.JA,coo.nnz,0);
				rsb_coo_sort(cor.VA,cor.IA,cor.JA,mtxAp->nnz,nrA,ncA,typecode,RSB_FLAG_NOFLAGS);
				if(RSB_SOME_ERROR(errval))
				{ RSB_ERROR(RSB_ERRM_ES);goto erri;}

				ct = - rsb_time();
				errval = rsb_mtx_get_csr(typecode,mtxAp, coo.VA, coo.IA, coo.JA,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS);
				if(RSB_SOME_ERROR(errval))
				{ RSB_ERROR(RSB_ERRM_ES);goto erri;}
				ct += rsb_time();
				for(i=0;i<mtxAp->nnz;++i)if(coo.JA[i]!=cor.JA[i]){RSB_ERROR("@%d: %d != %d!\n",i,coo.JA[i],cor.JA[i]);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
				if(RSB_SOME_ERROR(errval=rsb__csr_chk(coo.IA,coo.JA,coo.nr,coo.nc,coo.nnz,mib)))
					{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
				RSBENCH_STDOUT("#extraction of %zd elements in CSR took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);

/*				ct = - rsb_time();*/
/*				errval = rsb__do_get_coo(mtxAp,&coo.VA,&coo.IA,&coo.JA);	// FIXME : bugged ?*/
/*				if(RSB_SOME_ERROR(errval)) goto erri;*/
/*				ct += rsb_time();*/
/*				if(RSB_SOME_ERROR(rsb__util_is_sorted_coo_as_row_major(coo.IA,coo.JA,coo.nnz,coo.typecode,*/
/*					NULL,RSB_FLAG_NOFLAGS)))*/
/*				{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);} */
/*				RSBENCH_STDOUT("#extraction of %zd elements in sorted COO took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);*/

				rsb__util_coo_array_set(coo.IA,coo.nnz,0);
				rsb_coo_sort(cor.VA,cor.JA,cor.IA,mtxAp->nnz,ncA,nrA,typecode,RSB_FLAG_NOFLAGS);
				ct = - rsb_time();
				errval = rsb__do_get_csc(mtxAp,(rsb_byte_t**) &coo.VA,&coo.JA,&coo.IA);
				if(RSB_SOME_ERROR(errval))
					{goto erri;}
				ct += rsb_time();
				for(i=0;i<mtxAp->nnz;++i)if(coo.IA[i]!=cor.IA[i]){RSB_ERROR("@%d: %d != %d!\n",i,coo.IA[i],cor.IA[i]);errval = RSB_ERR_INTERNAL_ERROR;goto err;}
				if(RSB_SOME_ERROR(rsb__csc_chk(coo.JA,coo.IA,coo.nr,coo.nc,coo.nnz,mib)))
					{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
				RSBENCH_STDOUT("#extraction of %zd elements in CSC took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);

				{
					struct rsb_mtx_t * cmatrix=NULL;
					ct = - rsb_time();
					cmatrix = rsb__mtx_clone_simple(mtxAp);
					ct += rsb_time();
					if(!cmatrix){errval = RSB_ERR_ENOMEM;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					if(!rsb__mtx_chk(cmatrix))
						{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					RSB_MTX_FREE(cmatrix);
				}
				RSBENCH_STDOUT("#cloning of %zd elements took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);
				{
					struct rsb_mtx_t * cmatrix=NULL;
					cmatrix = rsb__mtx_clone_simple(mtxAp);
					if(!cmatrix){errval = RSB_ERR_ENOMEM;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					ct = - rsb_time();
					errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo(cmatrix,RSB_BOOL_FALSE);
					ct += rsb_time();
					if(!rsb__mtx_chk(cmatrix))
						{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					if(
rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(cmatrix,RSB_MATRIX_STORAGE_BCOR,RSB_FLAG_USE_HALFWORD_INDICES_CSR)
					!= rsb__terminal_recursive_matrix_count(cmatrix))
						{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}

					RSBENCH_STDOUT("#conversion of %zd elements to RCOO took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);
					RSB_MTX_FREE(cmatrix);
				}

				{
					struct rsb_mtx_t * cmatrix=NULL;
					struct rsb_coo_mtx_t icoo;
					cmatrix = rsb__mtx_clone_simple(mtxAp);
					if(!cmatrix){errval = RSB_ERR_ENOMEM;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					ct = - rsb_time();
					errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_sorted(cmatrix,&icoo);
					ct += rsb_time();

					if(RSB_SOME_ERROR(rsb__util_is_sorted_coo_as_row_major(icoo.IA,icoo.JA,icoo.nnz,icoo.typecode,NULL,RSB_FLAG_NOFLAGS)))
						{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					RSBENCH_STDOUT("#conversion of %zd elements to sorted COO took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);
					rsb__destroy_coo_matrix_t(&icoo);
				}
				
				if(!RSB_DO_TOOFEWNNZFORCSR(mtxAp->nnz,mtxAp->nr))
				{
					struct rsb_mtx_t * cmatrix=NULL;
					struct rsb_coo_mtx_t icoo;
					cmatrix = rsb__mtx_clone_simple(mtxAp);
					if(!cmatrix){errval = RSB_ERR_ENOMEM;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					ct = - rsb_time();
					errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_csr(cmatrix,&icoo);
					ct += rsb_time();
					if(RSB_SOME_ERROR(rsb__csr_chk(icoo.IA,icoo.JA,icoo.nr,icoo.nc,icoo.nnz,mib)))
						{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					RSBENCH_STDOUT("#conversion of %zd elements to CSR took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);
					rsb__destroy_coo_matrix_t(&icoo);
				}

				if(!RSB_DO_TOOFEWNNZFORCSR(mtxAp->nnz,mtxAp->nc))
				{
					struct rsb_mtx_t * cmatrix=NULL;
					struct rsb_coo_mtx_t icoo;
					cmatrix = rsb__mtx_clone_simple(mtxAp);
					if(!cmatrix){errval = RSB_ERR_ENOMEM;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					ct = - rsb_time();
					errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_csc(cmatrix,&icoo);
					ct += rsb_time();
					if(RSB_SOME_ERROR(rsb__csc_chk(icoo.JA,icoo.IA,icoo.nr,icoo.nc,icoo.nnz,mib)))
						{errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}

					RSBENCH_STDOUT("#conversion of %zd elements to CSC took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);
					rsb__destroy_coo_matrix_t(&icoo);
				}

				{
					struct rsb_mtx_t * cmatrix=NULL;
					struct rsb_coo_mtx_t icoo;
					cmatrix = rsb__mtx_clone_simple(mtxAp);
					if(!cmatrix){errval = RSB_ERR_ENOMEM;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
					ct = - rsb_time();
					errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_unsorted(cmatrix,&icoo);
					ct += rsb_time();

					RSBENCH_STDOUT("#conversion of %zd elements to unsorted COO took %2.5f s: %2.5f%% of construction time\n",(rsb_printf_int_t)rnz,ct,(100*ct)/mtxAp->tat);
					rsb__destroy_coo_matrix_t(&icoo);
				}
errc:
				rsb__destroy_coo_matrix_t(&coo);
				rsb__destroy_coo_matrix_t(&cor);
			} /* want_convert */

			if(RSB_SOME_ERROR(errval))
				RSB_PERR_GOTO(erri,RSB_ERRM_MASM);

			if(!mtxAp)
			{
				errval = RSB_ERR_INTERNAL_ERROR;
				RSB_PERR_GOTO(erri,RSB_ERRM_MASM);
			}

			if( want_verbose > 0 )
			{
				RSBENCH_STDOUT("# Constructed matrix (took %.3lfs): ", mct);
				RSBENCH_STDOUT(RSB_PRINTF_MTX_SUMMARY_ARGS(mtxAp));
				RSBENCH_STDOUT("\n");
			}

			ht = -rsb_time();
			if(!rsb__mtx_chk(mtxAp))
			{
				RSB_ERROR("matrix does not seem to be built correctly\n");
				errval = RSB_ERR_INTERNAL_ERROR;
				goto erri;
			}
			ht += rsb_time();
			totht += ht;
			if( want_verbose > 0 )
				RSBENCH_STDOUT("# matrix consistency check took %.3lfs (ok)\n",ht);

dnl			RSB_DO_ERROR_CUMULATE(errval,rsb__do_print_matrix_stats(mtxAp,RSB_CONST_DUMP_RECURSION_BRIEF,NULL));
dnl			if(RSB_SOME_ERROR(errval)){goto err;}

			if(zsort_for_coo)
				rsb__do_zsort_coo_submatrices(mtxAp);
			if(reverse_odd_rows)
				rsb__do_reverse_odd_rows(mtxAp);

			if(b_w_filename || csr_w_filename)
			{
				const char * w_filename = b_w_filename ;
				rsb_dump_flags_t dflags = RSB_CONST_DUMP_RSB;

				if(csr_w_filename)
					w_filename = csr_w_filename,
					dflags = RSB_CONST_DUMP_CSR;

				frt = -rsb_time();
				errval = rsb__do_print_matrix_stats(mtxAp,dflags,w_filename);
				frt += rsb_time();
				rsb_perror(NULL,errval);
				if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_NO_XDR); }
				RSB_STDOUT("#file output of %s took %lf s (%.0lf nnz, %.0lf nnz/s ) (%.5lf MB/s ) \n",rsb__basename(w_filename),frt,
					(((double)mtxAp->nnz)),
					(((double)mtxAp->nnz)/frt),
					(((double)rsb__sys_filesize(w_filename))/(frt*RSB_INT_MILLION))
				);
				goto ret;
			}

			if(dumpout_internals)
			{
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_print_matrix_stats(mtxAp,RSB_CONST_DUMP_RECURSION,NULL));
				if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}

			rsb__get_blocking_size(mtxAp,&br,&bc);

			if(br<0 || br <0)
			{
				RSB_ERROR("problems getting blocking size");
				goto erri;
			}

			if(!RSB_IS_VALID_NNZ_COUNT(nnz))
			{
				errval = RSB_ERR_INTERNAL_ERROR;
				RSB_PERR_GOTO(erri,RSB_ERRM_INZCL);
			}
			/* NOTE: if loading from a binary dump, we need to set nrA,ncA */
			nrA = mtxAp->nr;
			ncA = mtxAp->nc;
			ndA = RSB_MAX(nrA,ncA);
dnl			outnri = rhsnri = ndA;
dnl
dnl			ldX = (RSB_DOES_NOT_TRANSPOSE(transA) ? nrA : ncA) * incX; 	/* FIXME: still unused, e.g. in rsb__do_spmm_general */
dnl			ldY = (RSB_DOES_NOT_TRANSPOSE(transA) ? ncA : nrA) * incY;
dnl
			lhs = rsb__calloc((mtxAp->el_size*(ndA+br))*nrhs*incY);
			rhs = rsb__calloc((mtxAp->el_size*(ndA+bc))*nrhs*incX);

			if(!lhs || !rhs)
			{
				RSB_ERROR("problems allocating vectors");
				RSB_CONDITIONAL_FREE(lhs);
				RSB_CONDITIONAL_FREE(rhs);
				{ errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
			}

			if(RSB_SOME_ERROR(rsb__vectors_reinit(rhs,lhs,typecode,ndA,ndA,incX,incY))) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
dnl
dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),1,`dnl
			if( RSB__APPROPRIATE_AT_TIME_SPEC( merge_experimental ) || RSB__APPROPRIATE_AT_TIME_SPEC( split_experimental ) || just_enter_tuning ) /* FIXME: pass parameter */
			{
				struct rsb_mtx_t*mtxOp = NULL;
				int wvmbat = RSB_AUT0_TUNING_SILENT; /* wanted verbosity in merge based autotuning */
				int eps = 0; /* effective partitioning steps */
				rsb_time_t btt = RSB_TIME_ZERO; /* blocks tuning time */
				const rsb_bool_t want_auto_threads = ( merge_experimental < 0 || split_experimental < 0 ) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;
				const rsb_submatrix_idx_t maxms = RSB_ABS( merge_experimental );
				const rsb_submatrix_idx_t maxss = RSB_ABS( split_experimental );
				const int maxr = ( just_enter_tuning == 0 || (merge_experimental == 0 && split_experimental == 0 )) ? 0 : RSB_CONST_AUTO_TUNING_ROUNDS;
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
				const enum rsb_op_t op = rsb_op_spmv;
')dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
				const enum rsb_op_t op = rsb_op_spsvlt;
')dnl
dnl
				const rsb_time_t maxtime = /* RSB_AT_TIME_AUTO*/ RSB_AT_MAX_TIME;
				struct rsb_mtx_t mtxA = *mtxAp;
				const rsb_submatrix_idx_t lmn = mtxAp->all_leaf_matrices_n;

				/* please note at_mkl_csr_nt in the following... */
				if(want_auto_threads) { at_mkl_csr_nt = me_at_nt = RSB_THREADS_AUTO; }
				
				if(want_verbose >= -1)
					RSBENCH_STDOUT("RSB Sparse Blocks Autotuner invoked requesting max %d splits and max %d merges in %d rounds, threads spec.%d (specify negative values to enable threads tuning).\n",maxss,maxms,maxr,me_at_nt);

				if (want_verbose_tuning > 0)
					wvmbat = RSB_AUT0_TUNING_VERBOSE;
				if (want_verbose_tuning > 1)
					wvmbat = RSB_AUT0_TUNING_QUATSCH ;
				if (want_verbose_tuning > 2)
					wvmbat = RSB_AUT0_TUNING_QUATSCH + 1;
				btt -= rsb_time(); 

				mtxOp = mtxAp;
				errval = rsb__tune_spxx(&mtxOp,NULL,&me_at_nt,maxr,maxms,maxss,mintimes,maxtimes,maxtime,transA,alphap,NULL,nrhs,order,rhs,rhsnri,betap,lhs,outnri,op,&eps,&me_best_t,&me_at_best_t,wvmbat,mtrfpf,rsb__basename(filename),&attr,&otpos,&btpos);

				btt += rsb_time(); 
				tottt += btt;
dnl
				if(want_perf_dump) /* FIXME: shall give only values from the tuning routine */
				if(RSB_DO_FLAG_HAS(/*mtxAp->*/flags,RSB_FLAG_QUAD_PARTITIONING))
					rsb__pr_set(rspr, &mtxA, me_at_best_t<me_best_t?mtxOp:NULL, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti, transA, me_best_t, RSB_CONST_IMPOSSIBLY_BIG_TIME, me_at_best_t, RSB_CONST_IMPOSSIBLY_BIG_TIME, me_at_nt, RSB_THREADS_AUTO, btt, eps, &otpos, &btpos, NULL, NULL);
dnl
				if( mtxAp != mtxOp && mtxOp )
			 	{
#if RSB_AT_DESTROYS_MTX
					mtxAp = mtxOp;
					RSBENCH_STDOUT("First run of RSB Autotuner took %lg s  (%.3le s -> %.3le s per mop) (tuned: %d -> %d lsubm).\n",btt,otpos.min,btpos.min,(int)lmn,(int)(mtxAp->all_leaf_matrices_n));
#else  /* RSB_AT_DESTROYS_MTX */
#if 1
 					/* FIXME: this is to have mtxAp address constant. */
					errval = rsb__mtx_transplant_from_clone(&mtxAp, mtxOp);
					mtxOp = NULL;
					if(RSB_SOME_ERROR(errval)) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
#else
					RSB_MTX_FREE(mtxAp);
					mtxAp = mtxOp;
#endif
#endif /* RSB_AT_DESTROYS_MTX */
				}
				else
					RSBENCH_STDOUT("First run of RSB Autotuner took %lg s and did not change matrix.\n",btt);
dnl
			}
')dnl
dnl

dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),1,`dnl
			if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
			if(RSB__APPROPRIATE_AT_TIME_SPEC( want_autotuner ))
			{
				rsb_int_t otn = wat;
				rsb_int_t*otnp = NULL;
				rsb_real_t sf = RSB_REAL_ZERO;
				rsb_time_t att = - rsb_time();
				struct rsb_mtx_t * mtxOp = NULL;
				struct rsb_mtx_t ** mtxOpp = NULL;
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),1,`dnl
				const enum rsb_op_t op = rsb_op_spmv;
',`dnl
				const enum rsb_op_t op = rsb_op_spsvlt;
')dnl
				const rsb_submatrix_idx_t lmn = mtxAp->all_leaf_matrices_n;

				if(wat >  0)
					otnp = &otn; /* starting thread suggestion */
				if(wat == 0)
				{
					otnp = NULL; /* current thread count */
					mtxOpp = &mtxOp; /* matrix structure tuning */
				}
				if(wat <  0)
				{
					otn = -wat; /* ;-) */
					otnp = &otn; /* starting thread suggestion */
					mtxOpp = &mtxOp; /* matrix structure tuning */
				}
				RSBENCH_STDOUT("RSB Sparse Blocks Autotuner invoked requesting max %d splits and max %d merges in %d rounds, auto threads spec.\n",0,0,wai);
				errval = rsb__tune_spxx(mtxOpp, &sf, otnp, wai, 0, 0, mintimes, maxtimes, want_autotuner, transA, alphap, mtxAp, nrhs, order, rhs, rhsnri, betap, lhs, outnri, op , NULL, NULL, NULL, wavf, mtrfpf, rsb__basename(filename), &attr, &otpos, &btpos);
				att += rsb_time();
				tottt += att;
				if(mtxOpp && *mtxOpp)
				{
					RSBENCH_STDOUT("Second run of RSB Autotuner took %lg s and estimated a speedup of %lf x (%.3le s -> %.3le s per op) in new matrix (%d -> %d lsubm)\n",att,sf,otpos.min,btpos.min,(int)lmn,(int)((*mtxOpp)->all_leaf_matrices_n));
					RSBENCH_STDOUT("RSB Autotuner suggested a new matrix: freeing the old one.\n");
					RSB_MTX_FREE(mtxAp);
					mtxAp = mtxOp;
					mtxOp = NULL;
					mtxOpp = NULL;
				}
				else
					RSBENCH_STDOUT("Second run of RSB Autotuner took %lg s and estimated a speedup of %lf x (%.3le s -> %.3le s per op) in same matrix (%d -> %d lsubm)\n",att,sf,otpos.min,btpos.min,(int)lmn,(int)(mtxAp->all_leaf_matrices_n));
				if(wat && otn > 0)
				{
					/* FIXME: this breaks consistency! Shall skip further cycles!  */
					RSBENCH_STDOUT("Setting autotuning suggested thread count of %d (will skip further thread number configurations!)\n",otn);
					/* rsb__set_num_threads(otn); */
					RSB_DO_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_EXECUTING_THREADS,&otn,errval);
dnl
					if(want_ancillary_execs == RSB_BOOL_TRUE)
					if(incX == 1 && incY == 1)
					{
						totatt -= rsb_time();
						RSBENCH_STDOUT("# Post-autotuning performance recheck:\n");
						/* errval = */ rsb__do_bench_spxm(NULL,NULL,transA,alphap,mtxAp,nrhs,order,rhs,rhsnri,betap,lhs,outnri,RSB_AT_TIME_AUTO,RSB_AT_NTIMES_AUTO,op,10,RSB_AUT0_TUNING_QUATSCH,NULL,NULL); /* just for check purposes */
						totatt += rsb_time();
					}
dnl
					cc=otn;cl=ci+1;
				}
			}	/* want_autotuner */

			if(RSB_SOME_ERROR(errval)) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
')dnl
dnl
				if(n_dumpres)
				{
					RSBENCH_STDOUT("##RSB LHS %zd elements pre-peek:\n",(rsb_printf_int_t)n_dumpres);
					rsb__debug_print_vector(rhs,RSB_MIN(ndA*nrhs,n_dumpres),typecode,incX);
				}
				if(n_dumprhs)
				{
					RSBENCH_STDOUT("##RSB RHS %zd elements pre-peek:\n",(rsb_printf_int_t)n_dumprhs);
					rsb__debug_print_vector(rhs,RSB_MIN(ndA*nrhs,n_dumprhs),typecode,incX);
				}
dnl
			if ( times >= 0 ) /* benchmark of mop */
			{
dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
				/* 20140616 use this in conjunction with --dump-n-lhs-elements .. */
				for(nrhsl=0;nrhsl<nrhs;++nrhsl)
					rsb__util_set_array_to_converted_integer(((rsb_byte_t*)rhs)+mtxAp->el_size*ndA*nrhsl,mtxAp->typecode,ndA,incX,nrhsl+1),
					rsb__util_set_array_to_converted_integer(((rsb_byte_t*)lhs)+mtxAp->el_size*ndA*nrhsl,mtxAp->typecode,ndA,incY,nrhsl+1);
dnl
')dnl
dnl
				if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("POST_RSB_SPMV_",0,times,NULL);
				op_t = - rsb_time();
dnl
				RSB_TM_LIKWID_MARKER_R_START("RSB_SPMV");
				for(i=0;i<times;++i)  /* benchmark loop of mop begin */
				{
dnl				int e=0;
dnl				// FIXME : the following should be used only as a bugfix
dnl				//if(RSB_SOME_ERROR(rsb__cblas_Xscal(mtxAp->typecode,nrA+br,NULL,lhs,incY))) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
dnl //#if RSB_WANT_SPSV_STABILITY_FIX
#if 0
	{
				/* an extreme debugging measure */
				rsb_nnz_index_t ii;
				if(RSB_SOME_ERROR(rsb__cblas_Xscal(mtxAp->typecode,ndA,NULL,rhs,incX))) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
				for(ii=0;ii<nnz;++ii)rsb__util_increase_by_one(rhs,IA[ii],typecode);
				RSB_DO_ERROR_CUMULATE(errval,rsb__xcopy(lhs,rhs,0,0,mtxAp->nr,mtxAp->el_size));
	}
#else /* 0 */
				if(RSB_SOME_ERROR(rsb__fill_with_ones(rhs,mtxAp->typecode,ndA,incX))){ errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
				RSB_DO_ERROR_CUMULATE(errval,rsb__xcopy(lhs,rhs,0,0,mtxAp->nr,mtxAp->el_size));
#endif /* 0 */
')dnl
dnl
dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
#if RSB_EXPERIMENTAL_WANT_BEST_TIMES
				spsv_d_t -= rsb_time();

				if((errval = rsb__do_spsv_general(transA,alphap,mtxAp,lhs,incX,lhs,incY,RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR)
				{
					RSB_ERROR(RSB_ERRM_ES);
					goto err;
				}

				spsv_d_t += rsb_time();
				if(RSB_SOME_ERROR(rsb__fill_with_ones(rhs,mtxAp->typecode,ndA,incX))){ errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
				RSB_DO_ERROR_CUMULATE(errval,rsb__xcopy(lhs,rhs,0,0,mtxAp->nr,mtxAp->el_size));

				spsv_spmv_t -= rsb_time();
				/* y <- y + A x */
				if((errval = rsb__do_spmm_general(mtxAp,rhs,lhs,&pone[0],&pone[0],incX,incY,transA,RSB_OP_FLAG_DEFAULT,order RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR)
					goto err;
				spsv_spmv_t += rsb_time();
				best_spsv_spmv_t = RSB_MIN_ABOVE_INF(spsv_spmv_t,best_spsv_spmv_t,tinf);
				if(RSB_SOME_ERROR(rsb__fill_with_ones(rhs,mtxAp->typecode,ndA*nrhs,incX))){ errval = RSB_ERR_INTERNAL_ERROR; goto erri; } 
				RSB_DO_ERROR_CUMULATE(errval,rsb__xcopy(lhs,rhs,0,0,mtxAp->nr,mtxAp->el_size));

				spsv_f_t -= rsb_time();
				if(want_ancillary_execs == RSB_BOOL_TRUE)
				if((errval = rsb__do_spsv_general(transA,alphap,mtxAp,lhs,incX,lhs,incY,RSB_OP_FLAG_FAKE_LOCK RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR)
				{
					RSB_ERROR(RSB_ERRM_ES);
					goto err;
				}
				/* FIXME: if RSB_OUTER_NRHS_SPMV_ARGS_IDS defined to empty string, will not handle properly nrhs! */
#if 0
				if((errval = rsb__do_spsv_general(transA,alphap,mtxAp,lhs,incX,lhs,incY,RSB_OP_FLAG_DEFAULT RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR)/* `mop' is mop*/
				{
					RSB_ERROR(RSB_ERRM_ES);
					goto err;
				}
#endif
				spsv_f_t += rsb_time();
				/* if(RSB_SOME_ERROR(rsb__fill_with_ones(rhs,mtxAp->typecode,ndA,incX))){ errval = RSB_ERR_INTERNAL_ERROR; goto erri; } */
				for(nrhsl=0;nrhsl<nrhs;++nrhsl)
					rsb__util_set_array_to_converted_integer(((rsb_byte_t*)rhs)+mtxAp->el_size*ndA*nrhsl,mtxAp->typecode,ndA,incX,nrhsl+1),
					rsb__util_set_array_to_converted_integer(((rsb_byte_t*)lhs)+mtxAp->el_size*ndA*nrhsl,mtxAp->typecode,ndA,incY,nrhsl+1);
				/* RSB_DO_ERROR_CUMULATE(errval,rsb__xcopy(lhs,rhs,0,0,mtxAp->nr,mtxAp->el_size)); */
#endif /* RSB_EXPERIMENTAL_WANT_BEST_TIMES */
')dnl
dnl
#if RSB_EXPERIMENTAL_WANT_BEST_TIMES
				spmv_t = - rsb_time();
#endif /* RSB_EXPERIMENTAL_WANT_BEST_TIMES */
dnl spmv_uauz
dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
dnl
dnl	yes, we are passing the same vector (lhs)
dnl	
				if((errval = rsb_spsm(transA,alphap,mtxAp,nrhs,order,betap,lhs,incY*mtxAp->nr,lhs,incY*mtxAp->nr))!=RSB_ERR_NO_ERROR) /* benchmark -- `mop' is mop*/
				{
					RSBENCH_STDERR("[!] "RSB_ERRM_TS);
					goto erri;
				}
')dnl
dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
dnl
				RSB_PERFORMANCE_COUNTERS_DUMP("PRE_RSB_SPMV_",0);
dnl				if((e = rsb__do_spmm_general(mtxAp,rhs,lhs,alphap,betap,incX,incY,transA,RSB_OP_FLAG_DEFAULT,order RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR)
				if((errval = rsb__do_spmm_general(mtxAp,rhs,lhs,alphap,betap,incX,incY,transA,RSB_OP_FLAG_DEFAULT,order RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR) /* benchmark -- `mop' is mop */
				{
					RSBENCH_STDERR("[!] "RSB_ERRM_MV);
					goto erri;
				}
dnl				else
dnl
				RSB_PERFORMANCE_COUNTERS_DUMP("POST_RSB_SPMV_",1);
')dnl
dnl
ifelse(mop,`infty_norm',`dnl
				{
				double infinity_norm=RSB_REAL_ZERO;
				double infinity_norm_2=RSB_REAL_ZERO;

				if((errval = rsb__do_matrix_norm( mtxAp, &infinity_norm, RSB_EXTF_NORM_INF))!=RSB_ERR_NO_ERROR)
				{
					RSBENCH_STDERR("[!] some problem occurred in rsb_infinity_norm!\n");
					errval = RSB_ERR_INTERNAL_ERROR;
					goto erri;
				}
				else
				if(g_debug)
				{
					RSBENCH_STDOUT("matrix infinity norm is %lg\n",infinity_norm);
					if((errval = rsb__do_infinity_norm(mtxAp,&infinity_norm_2,1,transA))!=RSB_ERR_NO_ERROR)
					{
						RSBENCH_STDERR("[!] some problem occurred in rsb__do_infinity_norm!\n");
						errval = RSB_ERR_INTERNAL_ERROR;
						goto erri;
					}
					if(infinity_norm != infinity_norm_2)
					{
						RSB_ERROR("Mismatch while computing infinity norm : \n"
							"%lg != %lg\n",
							infinity_norm,infinity_norm_2);
						errval = RSB_ERR_INTERNAL_ERROR;
						goto erri;
					}
					RSBENCH_STDOUT("Infinity norm check passed.\n");

				}}
')dnl
ifelse(mop,`negation',`dnl
				/* FIXME: this section is obsolete and shall be removed (with the corresponding m4 macros). */
				int please_fix_RSB_M4_ARGS_TO_ACTUAL_ARGS=-1;/* here to fix negation */
				if(g_debug)
				{
					struct rsb_mtx_t * matrix2 = rsb__mtx_clone_simple(mtxAp);
					void *new_VA=NULL;
					rsb_coo_idx_t *new_IA=NULL, *new_JA=NULL;
					RSBENCH_STDOUT("cloning ...\n");

					if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc(&new_VA,&new_IA,&new_JA,nnz,typecode,RSB_BOOL_TRUE))){ RSB_ERROR("an allocation problem occurred\n"); errval = RSB_ERR_INTERNAL_ERROR;goto erri;}
		
					RSBENCH_STDOUT("getting back matrix in coo format ... \n");
					if(RSB_SOME_ERROR(rsb_mtx_get_coo(mtxAp,VA,new_IA,new_JA,RSB_FLAG_C_INDICES_INTERFACE) ))
					{
						RSB_ERROR("rsb_mtx_get_coo returned with an error code\n");
						errval = RSB_ERR_INTERNAL_ERROR;
						goto erri;
					}
					if(RSB_SOME_ERROR(rsb_mtx_get_coo(matrix2,new_VA,new_IA,new_JA,RSB_FLAG_C_INDICES_INTERFACE) ))
					{
						RSB_ERROR("rsb_mtx_get_coo returned with an error code\n");
						errval = RSB_ERR_INTERNAL_ERROR;
						goto erri;
					}
					RSBENCH_STDOUT("sorting back reconstructed matrix in blocks  ... \n");
					if( RSB_SOME_ERROR(rsb_util_sortcoo(    VA, new_IA, new_JA, nnz, typecode, M_b, K_b, p_r, p_c , flags )) ||
					    RSB_SOME_ERROR(rsb_util_sortcoo(new_VA, new_IA, new_JA, nnz, typecode, M_b, K_b, p_r, p_c , flags )))
						{errval = RSB_ERR_INTERNAL_ERROR; goto erri;}
	
					RSB_MTX_FREE(matrix2);

					if((errval = rsb_do_negation( mtxAp, please_fix_RSB_M4_ARGS_TO_ACTUAL_ARGS,transA ))!=RSB_ERR_NO_ERROR)
					{
						RSBENCH_STDERR("[!] some problem occurred in negation!\n");
						errval = RSB_ERR_INTERNAL_ERROR;
						goto erri;
					}
					else
					{
						if(RSB_MEMCMP(new_VA,VA,mtxAp->el_size*mtxAp->element_count)!=0)
						{
							RSB_ERROR("negation cross check failed!\n");
							errval = RSB_ERR_INTERNAL_ERROR;
							goto erri;
						}
						else
							RSBENCH_STDOUT("matrix negation cross check successfully passed\n");
					}

					RSB_CONDITIONAL_FREE(new_VA);
					RSB_CONDITIONAL_FREE(new_IA);
					RSB_CONDITIONAL_FREE(new_JA);
				}
				else
				{
					if((errval = rsb_do_negation( mtxAp, please_fix_RSB_M4_ARGS_TO_ACTUAL_ARGS, transA  ))!=RSB_ERR_NO_ERROR)
					{
						RSBENCH_STDERR("[!] some problem occurred in negation!\n");
						errval = RSB_ERR_INTERNAL_ERROR;
						goto erri;
					}
				
				}

')dnl
dnl				RSBENCH_STDERR("[!] Unimplemented!\n");
dnl				goto erri;
#if RSB_EXPERIMENTAL_WANT_BEST_TIMES
				spmv_t += rsb_time();
				tot_t += spmv_t;
				best_t = RSB_MIN_ABOVE_INF(spmv_t,best_t,tinf);
#endif /* RSB_EXPERIMENTAL_WANT_BEST_TIMES */
				if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
dnl
dnl
dnl
dnl
dnl
dnl	After sampling time.
dnl
dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
				if((g_debug || 1) && i==times-1)
				{
					/* this is debug information, very cheap to include */
					rsb_byte_t * out2=NULL;
					rsb_aligned_t mbetap[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
					out2 = rsb__calloc(mtxAp->el_size*(RSB_MAX(nrA,ncA)+br)*nrhs);
					if(!out2 /* || rsb__cblas_Xscal(mtxAp->typecode,nrA+br,NULL,out2,incY)*/) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
					if(RSB_SOME_ERROR(rsb__fill_with_ones(alphap,typecode,1,1))){ errval = RSB_ERR_INTERNAL_ERROR; goto err;}
					if(RSB_SOME_ERROR(rsb__fill_with_ones(mbetap,typecode,1,1))){ errval = RSB_ERR_INTERNAL_ERROR; goto err;}
					if(RSB_SOME_ERROR(rsb__cblas_Xscal(typecode,1,NULL,errnorm,1))){ errval = RSB_ERR_INTERNAL_ERROR; goto err;}
					if((errval = rsb__do_spmm_general(mtxAp,lhs,out2,alphap,mbetap,incX,incY,transA,RSB_OP_FLAG_DEFAULT,order RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR)
					{
						RSBENCH_STDERR("[!] some problem occurred in sparse matrix vector product!\n");
						goto erri;
					}
					RSB_DO_ERROR_CUMULATE(errval,rsb__sqrt_of_sum_of_fabs_diffs(rhs,out2,errnorm,typecode,nrA+br));
					RSBENCH_STDOUT("#error norm:");
					RSB_DO_ERROR_CUMULATE(errval,rsb__debug_print_value(errnorm,typecode));
					RSBENCH_STDOUT("\n");
dnl					RSBENCH_STDOUT("#computed rhs:");
dnl					RSB_DO_ERROR_CUMULATE(errval,rsb__debug_print_vector(out2,nrA,typecode,1));
					if(out2)rsb__free(out2);
				}
')
dnl
	#ifdef RSB_WANT_KERNELS_DEBUG
				/* ... */
	#endif /* RSB_WANT_KERNELS_DEBUG */
				}  /* times: benchmark loop of mop end */
				RSB_TM_LIKWID_MARKER_R_STOP("RSB_SPMV");
				RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("POST_RSB_SPMV_",1,times,&rsb_pci);
dnl
ifelse(RSB_M4_IS_SPXV_KERNEL_MOP(mop),1,`dnl
				if( want_verbose > 0 )
				if((g_debug || 1) /*&& i==times-1*/)
				{
					/* debug info */
					RSB_DO_ERROR_CUMULATE(errval,rsb__do_print_some_vector_stats(lhs,typecode,nrA,incY));
				}
')dnl
dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
dnl
			if(rsb__set_num_threads(cc)!=cc)
			{
				RSB_ERROR("failed setting %d threads!\n",cc);
				errval = RSB_ERR_INTERNAL_ERROR;
				goto err;
			}
			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
dnl
			if(want_ancillary_execs == RSB_BOOL_TRUE)
			for(i=0;i<times;++i)
			{
				if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				no_lock_op_time = - rsb_time();
				if((errval = rsb__do_spmm_general(mtxAp,rhs,lhs,alphap,betap,incX,incY,transA,RSB_OP_FLAG_FAKE_LOCK,order RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR) { goto erri; }
				no_lock_op_time += rsb_time();
				no_lock_op_time_best = RSB_MIN_ABOVE_INF(no_lock_op_time_best,no_lock_op_time,tinf);
				no_lock_op_tot_time += no_lock_op_time;
			}
			if(cc==1)serial_no_lock_op_time_best=no_lock_op_time_best;
			totatt += no_lock_op_tot_time;

			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));

			if(want_ancillary_execs == RSB_BOOL_TRUE)
			if(cc==1)
			for(i=0;i<times;++i)
			{
				if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				qt_op_time = - rsb_time();
				if((errval = rsb__do_spmm_general(mtxAp,rhs,lhs,alphap,betap,incX,incY,transA,RSB_OP_FLAG_WANT_SERIAL,order RSB_OUTER_NRHS_SPMV_ARGS_IDS))!=RSB_ERR_NO_ERROR) { goto erri; }
				qt_op_time += rsb_time();
				qt_op_time_best = RSB_MIN_ABOVE_INF(qt_op_time_best,qt_op_time,tinf);
				qt_op_tot_time += qt_op_time;
			}
			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
			totatt += qt_op_tot_time;
dnl
')dnl
dnl

dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
				if((g_debug) /*&& i==times-1*/)
				{
					rsb_byte_t * out2=NULL;
					out2=rsb__calloc(mtxAp->el_size*(RSB_MAX(nrA,ncA)+br)*nrhs);
					if(!out2 /* || rsb__cblas_Xscal(mtxAp->typecode,nrA+br,NULL,out2,incY)*/) { errval = RSB_ERR_INTERNAL_ERROR; goto erri; }

					RSB_DO_FLAG_ADD(mtxAp->flags,RSB_FLAG_SHOULD_DEBUG);
/*					rsb_spmv_uaua_testing( mtxAp, rhs, out2,transA );	*//* FIXME : INCOMPLETE */
					RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_SHOULD_DEBUG);
					/* bit-per-bit checking */
					
					rsb__util_vector_sum(errnorm,lhs,typecode,nrA);
					RSBENCH_STDOUT("#sum:");
					rsb__debug_print_vector(errnorm,1,typecode,1);
					RSBENCH_STDOUT("\n");

					if(dumpvec&rsb_dumpvec_res)/* new */
						rsb__debug_print_vectors(lhs,out2,nrA,1,1,typecode);
					
					if(dumpvec&rsb_dumpvec_res)/* new */
					{
					if(RSB_MEMCMP(lhs,out2,mtxAp->el_size*(nrA+br*0))!=0)
					{
						RSB_ERROR("sparse matrix vector product cross check failed. diff (bad,good):\n");
						rsb__debug_print_vectors_diff(lhs,out2,nrA,typecode,incY,incY,RSB_VECTORS_DIFF_DISPLAY_N);

						if(out2)
							rsb__free(out2);
						{ errval = RSB_ERR_INTERNAL_ERROR; goto erri; }
					}
					else
						RSBENCH_STDOUT("sparse matrix vector product cross check succeeded\n");
					}
					if(out2)rsb__free(out2);
				}
')dnl
dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
				if(dumpvec&rsb_dumpvec_res)
					rsb__debug_print_vector(lhs,nrA,typecode,incY);
				if(dumpvec&rsb_dumpvec_rhs)
					rsb__debug_print_vector(rhs,nrA,typecode,incX);
')dnl

				if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				if(n_dumpres)
				{
					RSBENCH_STDOUT("##RSB LHS %zd elements post-peek:\n",(rsb_printf_int_t)n_dumpres);
					rsb__debug_print_vector(lhs,RSB_MIN(ndA*nrhs,n_dumpres),typecode,incY);
				}
				if(n_dumprhs)
				{
					RSBENCH_STDOUT("##RSB RHS %zd elements post-peek:\n",(rsb_printf_int_t)n_dumprhs);
					rsb__debug_print_vector(rhs,RSB_MIN(ndA*nrhs,n_dumprhs),typecode,incY);
				}
				if(!g_sort_only)
				{
					op_t += rsb_time();
					op_t /= (double)times;
					/*
				if(RSB_WANT_VERBOSE_MESSAGES)
				{RSBENCH_STDOUT("performed %lf Mflops in %lf seconds (%lf Mflops)\n",raw_Mflops, op_t, (raw_Mflops)/(op_t));
				RSBENCH_STDOUT("raw data rate of (%lf Gbytes/sec)\n", ((double)(raw_Mflops)*(mtxAp->el_size))/(op_t*1000.0));	}*/
				/*
				if(RSB_WANT_VERBOSE_MESSAGES)
				RSBENCH_STDOUT("nonzero data rate of (%lf Gbytes/sec, or %lf Mflops)\n",
				(true_Mflops*(mtxAp->el_size))/(op_t*1000.0),
				true_Mflops/(op_t)
				);*/
				}

ifelse(mop,`mat_stats',`',`dnl
                                fillin = rsb__do_get_matrix_fillin(mtxAp);
				if(g_sort_only)
				{
				/* FIXME :
				 * please note that in this rudimentary model we take in account also the matrix creationtime.
				 */
                	                raw_Mflops= (rsb_perf_t) mtxAp->element_count;
        	                        true_Mflops=(((double)mtxAp->nnz)*log((double)mtxAp->nnz))/RSB_REAL_MILLION;
					op_t=mct;	/* our timed operation is matrix construction */
				}
				else
				{
	                                raw_Mflops = RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_IDENTIFIER(mop)(mtxAp);
	                                true_Mflops = raw_Mflops/fillin;
	                                raw_Mflops *=nrhs;
	                                true_Mflops*=nrhs;
				}

dnl spmv_uauz
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),1,`dnl
#if RSB_WANT_MKL
	if(want_mkl_bench && !(cc==1 && mkl_coo_op_time_best_serial != RSB_CONST_IMPOSSIBLY_BIG_TIME))
	{
			rsb_nnz_idx_t annz = RSB_MAX(nnz,nrA+1),rnz=0,mklnz=nnz;
			/* please note that mkl routines do not support stride */
			/* FIXME: a non monotonically-increasing order will do harm */
			mkl_coo2csr_time = RSB_TIME_ZERO;
			mkl_coo_op_tot_time = RSB_TIME_ZERO;
			mkl_coo_op_time = RSB_TIME_ZERO;
			mkl_coo_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			//mkl_coo_op_time_best_serial = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			mkl_csr_op_tot_time = RSB_TIME_ZERO;
			mkl_csr_op_time = RSB_TIME_ZERO;
			mkl_csr_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			//mkl_csr_op_time_best_serial = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			
			if(nrhs>1)
				want_mkl_bench_coo = RSB_BOOL_FALSE;/* 20130401 FIXME: this circumvents an Intel MKL bug */
#if 1
			//mkl_set_dynamic(1);
			//RSBENCH_STDOUT("MKL failed enabling dynamic thread number control\n");
			mkl_set_num_threads(cc);
			//RSBENCH_STDOUT("MKL has %d threads now\n",mkl_get_num_threads());
#else /* 1 */
			if(rsb__set_num_threads(cc)!=cc)
			{
				RSB_ERROR("failed setting %d threads!\n",cc);
				errval = RSB_ERR_INTERNAL_ERROR;
				goto err;
			}
#endif /* 1 */
dnl			/* FIXME: for some matrices, the following invocation of rsb_coo_sort() causes memory leaks */
dnl			//errval = rsb_coo_sort(VA,IA,JA,nnz,nrA,ncA,typecode,RSB_FLAG_NOFLAGS);
			if(!want_sort_after_load)
			if(!want_in_place_assembly)
			{
				errval = rsb__util_sort_row_major_parallel(VA,IA,JA,nnz,nrA,ncA,typecode,flags_s);
				mklnz = rsb__weed_out_duplicates (IA,JA,VA,nnz,typecode,RSB_FLAG_SORTED_INPUT);
				if((!RSB_IS_VALID_NNZ_COUNT(mklnz)) || (!mklnz) || (RSB_SOME_ERROR(errval)))
				{
					RSB_PERR_GOTO(err,RSB_ERRM_EM);
				}
				annz = RSB_MAX(mklnz,nrA+1);
			}
			mkl_set_num_threads(cc); // necessary, or MKL will get puzzled

		if(want_mkl_bench_coo)
		{
			totct -= rsb_time();
			errval = rsb__util_coo_alloc_copy_and_stats(&M_VA,&M_IA,&M_JA,want_in_place_assembly?NULL:VA,want_in_place_assembly?NULL:IA,want_in_place_assembly?NULL:JA,NULL,NULL,mklnz,(annz-mklnz),typecode,0,mib,RSB_FLAG_NOFLAGS,NULL);
			if(RSB_SOME_ERROR(errval)){RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_ENOMEM);goto mklerr;}
			//errval = rsb_mtx_get_coo(mtxAp,M_VA,M_IA,M_JA,flags); /* FIXME: use this */
			errval = rsb__do_get_rows_sparse(RSB_DEFAULT_TRANSPOSITION,NULL,mtxAp,M_VA,M_IA,M_JA,0,mtxAp->nr-1,&rnz,RSB_FLAG_NOFLAGS|mif);
			totct += rsb_time();
	
			if(!M_VA  || !M_IA  || !M_JA ){RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_ENOMEM);goto mklerr;}

			RSB_DO_ERROR_CUMULATE(errval,rsb__vectors_reinit(rhs,lhs,typecode,ndA,ndA,incX,incY));
			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("PRE_MKL_COO_SPXV_",0,times,NULL);
			RSB_TM_LIKWID_MARKER_R_START("MKL_COO_SPMV");
			for(i=0;i<times;++i)
			{
				if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				mkl_coo_op_time = - rsb_time();
				RSB_PERFORMANCE_COUNTERS_DUMP("PRE_MKL_COO_SPXV_",0);
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),1,`dnl
				if(nrhs>1)
					RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_coo_spmm(M_VA,nrA,ncA,nrhs,mklnz,M_IA,M_JA,rhs,rhsnri,lhs,outnri,alphap,betap,transA,typecode,flags));
				else

					RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_coo_spmv(M_VA,nrA,ncA,mklnz,M_IA,M_JA,rhs,lhs,alphap,betap,transA,typecode,flags));
')dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
				RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_coo_spsv(M_VA,nrA,ncA,mklnz,M_IA,M_JA,rhs,lhs,alphap,betap,transA,typecode,flags));
')dnl
				RSB_PERFORMANCE_COUNTERS_DUMP("POST_1KL_COO_SPXV_",1);
				mkl_coo_op_time += rsb_time();
				mkl_coo_op_time_best = RSB_MIN_ABOVE_INF(mkl_coo_op_time_best,mkl_coo_op_time,tinf);
				mkl_coo_op_tot_time+=mkl_coo_op_time;
			}
			RSB_TM_LIKWID_MARKER_R_STOP("MKL_COO_SPMV");
				RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("POST_MKL_COO_SPXV_",1,times,&mkl_coo_pci);
			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
			if(n_dumpres)
			{
				RSBENCH_STDOUT("##MKL COO LHS %zd elements post-peek:\n",(rsb_printf_int_t)n_dumpres);
				rsb__debug_print_vector(lhs,RSB_MIN(ndA*nrhs,n_dumpres),typecode,incY);
			}
			if(cc==1) 
				mkl_coo_op_time_best_serial = mkl_coo_op_time_best;

			RSB_CONDITIONAL_FREE(M_VA);
			RSB_CONDITIONAL_FREE(M_IA);
			RSB_CONDITIONAL_FREE(M_JA);
		} /* want_mkl_bench_coo */

		if(want_mkl_bench_csr || RSB_MKL_APPROPRIATE_AT_TIME_SPEC( want_mkl_autotuner ) )
		{
			totct -= rsb_time();
			errval = rsb__util_coo_alloc_copy_and_stats(&M_VAC,&M_IAC,&M_JAC,want_in_place_assembly?NULL:VA,want_in_place_assembly?NULL:IA,want_in_place_assembly?NULL:JA,NULL,NULL,mklnz,(annz-mklnz),typecode,0,mib,RSB_FLAG_NOFLAGS,NULL);
			errval = rsb_mtx_get_csr(mtxAp->typecode,mtxAp,M_VAC,M_IAC,M_JAC,flags|mif);
			totct += rsb_time();
	
			if(!M_VAC || !M_IAC || !M_JAC) {RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_ENOMEM);goto mklerr;}
				// FIXME: Missing error handling !

                        if(0)/* if want bogus contents (for debug/inspection) */
                        {
                                rsb_coo_idx_t i,npr=(mklnz+nrA-1)/nrA;
                                rsb_nnz_idx_t l;
                                M_IAC[0]=0;
                                for(i=1;i<nrA;++i)
                                        M_IAC[i]=M_IAC[i-1]+npr;
                                for(i=0;i<nrA;++i)
                                        for(l=M_IAC[i];l<M_IAC[i+1];++l)
                                                M_JAC[l]=l-M_IAC[i];
                                M_IAC[nrA]=mklnz;
                        }

			totct -= rsb_time();
			if(!want_in_place_assembly)
			{
				mkl_coo2csr_time = - rsb_time();
				RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_coo2csr(nrA,ncA,mklnz,VA,IA,JA,M_VAC,M_IAC,M_JAC,typecode,mib));
				mkl_coo2csr_time += rsb_time();
				if(RSB_SOME_ERROR(errval))
				{
					RSB_PERR_GOTO(err,"MKL conversion of COO to CSR: error! Is the input suitable ? E.g. NO duplicates ?")
				}
				if(RSB_SOME_ERROR(rsb__csr_chk(M_IAC,M_JAC,nrA,ncA,mklnz,mib)))
				{
					RSB_PERR_GOTO(err,"MKL conversion of COO to CSR: error.")
				}
			}
			else
			{
				RSB_WARN("warning : skipping MKL coo2csr conversion (user chose in-place RSB build) \n");
			}
			totct += rsb_time();
		} /* want_mkl_bench_csr || want_mkl_autotuner */

			if(n_dumpres)
			{
				RSBENCH_STDOUT("##MKL CSR LHS %zd elements pre-peek:\n",(rsb_printf_int_t)n_dumpres);
				rsb__debug_print_vector(rhs,RSB_MIN(ndA*nrhs,n_dumpres),typecode,incX);
			}			RSB_DO_ERROR_CUMULATE(errval,rsb__vectors_reinit(rhs,lhs,typecode,ndA,ndA,incX,incY));
			if(n_dumprhs)
			{
				RSBENCH_STDOUT("##MKL CSR RHS %zd elements pre-peek:\n",(rsb_printf_int_t)n_dumprhs);
				rsb__debug_print_vector(lhs,RSB_MIN(ndA*nrhs,n_dumprhs),typecode,incY);
			}			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));

			if(want_mkl_bench_csr)
			{
			RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("PRE_MKL_CSR_SPXV_",0,times,NULL);
			RSB_TM_LIKWID_MARKER_R_START("MKL_CSR_SPMV");
			for(i=0;i<times;++i)
			{
				if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				mkl_csr_op_time = - rsb_time();
				RSB_PERFORMANCE_COUNTERS_DUMP("PRE_MKL_CSR_SPXV_",0);
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),1,`dnl
				if(nrhs>1)
					RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_csr_spmm_bench(M_VAC,nrA,ncA,nrhs,mklnz,M_IAC,M_JAC,rhs,rhsnri,lhs,outnri,alphap,betap,transA,typecode,flags|mif,NULL,NULL,NULL,NULL));
				else
					RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_csr_spmv_bench(M_VAC,nrA,ncA,mklnz,M_IAC,M_JAC,rhs,lhs,alphap,betap,transA,typecode,flags,NULL,NULL,NULL /* &mkl_csr_op_time */,NULL ));
')dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
				if(nrhs>1)
					RSB_DO_ERROR_CUMULATE(errval,rsb__do_mkl_csr_spsm(M_VAC,nrA,nrhs,M_IAC,M_JAC,rhs,lhs,alphap,transA,typecode,flags,nrhs/*ldX*/,nrhs/*ldY*/));
					/* FIXME: rsb__mkl_csr_spsm_bench is there */
				else
					RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_csr_spsv_bench(M_VAC,nrA,ncA,mklnz,M_IAC,M_JAC,rhs,lhs,alphap,betap,transA,typecode,flags,NULL,NULL,NULL /* &mkl_csr_op_time */,NULL));
')dnl
				RSB_PERFORMANCE_COUNTERS_DUMP("POST_MKL_CSR_SPXV_",1);
				mkl_csr_op_time += rsb_time();
				mkl_csr_op_time_best = RSB_MIN_ABOVE_INF(mkl_csr_op_time_best,mkl_csr_op_time,tinf);
				mkl_csr_op_tot_time+=mkl_csr_op_time;
			}
			RSB_TM_LIKWID_MARKER_R_STOP("MKL_CSR_SPMV");
			RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("POST_MKL_CSR_SPXV_",1,times,&mkl_csr_pci);
			} /* want_mkl_bench_csr */
			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
			if(cc==1)mkl_csr_op_time_best_serial=mkl_csr_op_time_best;
			if(n_dumpres)
			{
				RSBENCH_STDOUT("##MKL CSR LHS %zd elements post-peek:\n",(rsb_printf_int_t)n_dumpres);
				rsb__debug_print_vector(lhs,RSB_MIN(ndA*nrhs,n_dumpres),typecode,incY);
			}
			if(n_dumprhs)
			{
				RSBENCH_STDOUT("##MKL CSR RHS %zd elements post-peek:\n",(rsb_printf_int_t)n_dumprhs);
				rsb__debug_print_vector(rhs,RSB_MIN(ndA*nrhs,n_dumprhs),typecode,incY);
			}
#if RSB_WANT_OMP_RECURSIVE_KERNELS
			if( mkl_csr_op_time_best != RSB_CONST_IMPOSSIBLY_BIG_TIME )
				RSBENCH_STDOUT("##MKL STUFF DEBUG omp_set_num_threads():%d==omp_get_num_threads():%d  bestserialcsr:%0.5lf vs bestcsr:%0.5lf\n",omp_get_num_threads(),cc,mkl_csr_op_time_best_serial,mkl_csr_op_time_best);
			if( mkl_coo_op_time_best != RSB_CONST_IMPOSSIBLY_BIG_TIME )
				RSBENCH_STDOUT("##MKL STUFF DEBUG omp_set_num_threads():%d==omp_get_num_threads():%d  bestserialcoo:%0.5lf vs bestcoo:%0.5lf\n",omp_get_num_threads(),cc,mkl_coo_op_time_best_serial,mkl_coo_op_time_best);
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */

			if( RSB_MKL_APPROPRIATE_AT_TIME_SPEC( want_mkl_autotuner ) && want_mkl_autotuner > RSB_TIME_ZERO )
			{
				rsb_time_t btime = RSB_TIME_ZERO, matt = -rsb_time();
				rsb_thread_t bthreads = at_mkl_csr_nt;
				rsb_real_t sf = RSB_REAL_ZERO;
				rsb_char_t * ops = "";

ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
				rsb__tattr_init(&(attr.clattr), NULL, nrA, mklnz, typecode, flags, nrhs);
				attr.clattr.vl = 1; /* FIXME: new */
')dnl
				RSBENCH_STDOUT("# MKL CSR %s autotuning for thread spec. %d  trans %c (0=current (=%d),<0=auto,>0=specified)\n",ops,bthreads,RSB_TRANSPOSITION_AS_CHAR(transA),cc);
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),1,`dnl
#if RSB_WANT_MKL_INSPECTOR
				if( want_mkl_inspector == RSB_BOOL_TRUE )
				{
					const sparse_layout_t layout = ( order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER ) ? SPARSE_LAYOUT_COLUMN_MAJOR : SPARSE_LAYOUT_ROW_MAJOR;
					MKL_INT ldc = 0, ldb = 0;
					RSB_DO_ERROR_CUMULATE(errval,rsb__set_ldX_for_spmm(transA, mtxAp, nrhs, order, &ldb, &ldc));
					if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES); }
					RSBENCH_STDOUT("# Using the MKL \"Inspector-executor Sparse BLAS routines\" interface %s!\n",mkl_ie_hintlvl == 2?" (aggressive memory usage option)": (mkl_ie_hintlvl == 1 ? "(no memory usage option)" : "(no optimization)") );
					if(nrhs>1)
						RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_inspector_csr_spmm_bench(M_VAC,nrA,ncA,layout,nrhs,mklnz,M_IAC,M_JAC,rhs,ldb,lhs,ldc,alphap,betap,transA,typecode,flags|mif,&bthreads,&btime,&(attr.clattr),&btpms[0],mkl_ie_hintlvl));
					else
						RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_inspector_csr_spmv_bench(M_VAC,nrA,ncA,nrhs,mklnz,M_IAC,M_JAC,rhs,rhsnri,lhs,outnri,alphap,betap,transA,typecode,flags|mif,&bthreads,&btime,&(attr.clattr),&btpms[0],mkl_ie_hintlvl));
				} /* want_mkl_inspector */
				else
#endif /* RSB_WANT_MKL_INSPECTOR */
				{
					RSBENCH_STDOUT("# Using the MKL \"Sparse BLAS routines\" interface !\n");
					if(nrhs>1)
						RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_csr_spmm_bench(M_VAC,nrA,ncA,nrhs,mklnz,M_IAC,M_JAC,rhs,rhsnri,lhs,outnri,alphap,betap,transA,typecode,flags|mif,&bthreads,&btime,&(attr.clattr),&btpms[0]));
					else
						RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_csr_spmv_bench(M_VAC,nrA,ncA,mklnz,M_IAC,M_JAC,rhs,lhs,alphap,betap,transA,typecode,flags,&bthreads,&btime,&(attr.clattr),&btpms[0]));
				}
				ops = "SPMV";
')dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
dnl				RSBENCH_STDOUT("# MKL SPSV/SPSM appears to be serial --- skipping MKL threads autotuning !\n");
#if 1
				if(nrhs>1)
					RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_csr_spsm_bench(M_VAC,nrA,nrhs,M_IAC,M_JAC,rhs,lhs,alphap,transA,typecode,flags,nrhs/*ldX*/,nrhs/*ldY*/,&bthreads,&btime,&(attr.clattr),&btpms[0]));
				else
					RSB_DO_ERROR_CUMULATE(errval,rsb__mkl_csr_spsv_bench(M_VAC,nrA,ncA,mklnz,M_IAC,M_JAC,rhs,lhs,alphap,betap,transA,typecode,flags,&bthreads,&btime,&(attr.clattr),&btpms[0]));
				ops = "SPSV";
#endif
')dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
				bthreads = bthreads ? bthreads : cc;
')dnl
				RSBENCH_STDOUT("# MKL CSR %s best threads / time / perf. were: %d / %lg / %lg\n",ops,bthreads,btime,(rsb__estimate_mflops_per_op_spmv_uaua(mtxAp)*nrhs)/btime);
				matt += rsb_time();
				RSBENCH_STDOUT("MKL CSR Autotuner took %.2lgs and estimated a speedup of %lf / %lf = %lf x (best round %d samples at %d threads)\n",matt,(attr.clattr).dtpo,(attr.clattr).btpo,(attr.clattr).dtpo/(attr.clattr).btpo,attr.clattr.nit[attr.clattr.optt],attr.clattr.optt);
				at_mkl_csr_op_time_best = btime;
				at_mkl_csr_nt = bthreads;
				mkl_csr_op_time_best = (attr.clattr).dtpo;
				totmt += matt;
				RSB_ASSERT( bthreads > 0 );
			} /* want_mkl_autotuner */

			if(want_mkl_bench_gem)
			{
				rsb_coo_idx_t gemdim=0;
			RSB_DO_ERROR_CUMULATE(errval,rsb__vectors_reinit(rhs,lhs,typecode,ndA,ndA,incX,incY));
			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
			RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("PRE_MKL_GEMV_",0,times,NULL);
			RSB_TM_LIKWID_MARKER_R_START("MKL_GEMV");
			for(i=0;i<times;++i)
			{
				if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				mkl_gem_op_time = - rsb_time();
				RSB_PERFORMANCE_COUNTERS_DUMP("PRE_MKL_GEMV_",0);
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),1,`dnl
				if(nrhs>1)
					; /* FIXME */
				/* FIXME: missing error handling */
				rsb__mkl_gemv(typecode,VA,rhs,lhs,nnz,ndA,&gemdim);
')dnl
				RSB_PERFORMANCE_COUNTERS_DUMP("POST_MKL_GEMV_",1);
				mkl_gem_op_time += rsb_time();
				mkl_gem_op_time_best = RSB_MIN_ABOVE_INF(mkl_gem_op_time_best,mkl_gem_op_time,tinf);
				mkl_gem_op_tot_time+=mkl_gem_op_time;
			}
			true_gem_Mflops=2*gemdim*gemdim;
			true_gem_Mflops/=RSB_REAL_MILLION;
			RSB_TM_LIKWID_MARKER_R_STOP("MKL_GEMV");
			RSB_PERFORMANCE_COUNTERS_DUMP_MEAN("POST_MKL_GEMV_",1,times,&mkl_gem_pci);
			if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
			if(cc==1)mkl_gem_op_time_best_serial=mkl_gem_op_time_best;
			if(n_dumpres)
			{
				RSBENCH_STDOUT("##MKL GEMX LHS %zd elements peek:\n",(rsb_printf_int_t)n_dumpres);
				rsb__debug_print_vector(lhs,RSB_MIN(ndA*nrhs,n_dumpres),typecode,incY);
			}
			} /* want_mkl_bench_gem */
mklerr:
			RSB_CONDITIONAL_FREE(M_VAC);
			RSB_CONDITIONAL_FREE(M_IAC);
			RSB_CONDITIONAL_FREE(M_JAC);
			RSB_CONDITIONAL_FREE(M_VA);
			RSB_CONDITIONAL_FREE(M_IA);
			RSB_CONDITIONAL_FREE(M_JA);
			rsb_perror(NULL,errval);
		} /* want_mkl_bench  */
#endif /* RSB_WANT_MKL */
dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),1,`dnl
#if RSB_WANT_ARMPL
		if ( want_armpl_bench )
		{
			const rsb_coo_idx_t aib = 0; /* ARMPL index base */
			void *M_VA = NULL;
			armpl_int_t *M_IA = NULL, *M_JA = NULL;
			const rsb_nnz_idx_t annz = RSB_MAX(nnz,nrA+1), armplnz = nnz;
			rsb_thread_t bthreads = cc;
			struct rsb_ts_t btpms[3]; /* first is tuned, second is not */

			rsb_time_t armpl_csr_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME, at_armpl_csr_op_time_best = RSB_CONST_IMPOSSIBLY_BIG_TIME;
			errval = rsb__util_coo_alloc_copy_and_stats(&M_VA,&M_IA,&M_JA,NULL,NULL,NULL,NULL,NULL,armplnz,(annz-armplnz),typecode,0,aib,RSB_FLAG_NOFLAGS,NULL);
			RSB_DEBUG_ASSERT( errval == RSB_ERR_NO_ERROR);
			errval = rsb_mtx_get_csr(mtxAp->typecode,mtxAp,M_VA,M_IA,M_JA,flags);
			RSB_DEBUG_ASSERT( errval == RSB_ERR_NO_ERROR);
			//rsb__armpl_csr_spmv_bench(VA, nrA, ncA, nnz, M_IA, M_JA, lhs, rhs, alphap, betap, transA, typecode, flags, NULL, NULL, NULL, NULL);
			rsb__armpl_csr_spmv_bench(VA, nrA, ncA, nnz, M_IA, M_JA, lhs, rhs, alphap, betap, transA, typecode, flags, &bthreads, &armpl_csr_op_time_best, &(attr.clattr), &btpms[1],RSB_BOOL_FALSE);
			btpms[0] = btpms[1]; // only and best
			rsb__armpl_csr_spmv_bench(VA, nrA, ncA, nnz, M_IA, M_JA, lhs, rhs, alphap, betap, transA, typecode, flags, &bthreads, &at_armpl_csr_op_time_best, &(attr.clattr), &btpms[2],RSB_BOOL_TRUE);
			btpms[1] = btpms[2]; // only and best
			at_armpl_csr_op_time_best = RSB_MIN(armpl_csr_op_time_best, at_armpl_csr_op_time_best);
			rsb__pr_set(rspr, mtxAp, NULL, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti, transA, RSB_CONST_IMPOSSIBLY_BIG_TIME, armpl_csr_op_time_best, RSB_CONST_IMPOSSIBLY_BIG_TIME, at_armpl_csr_op_time_best, RSB_THREADS_AUTO, bthreads, RSB_CONST_IMPOSSIBLY_BIG_TIME, -1, NULL, NULL, &btpms[0], &btpms[1]);
			// printf("VALS: %lf -> %lf\n",btpms[1].min,btpms[0].min);
		} /* want_armpl_bench */
#endif /* RSB_WANT_ARMPL */
')dnl
dnl
dnl
#ifdef RSB_WANT_OSKI_BENCHMARKING 
			/* FIXME : should only exist for double as type */
			if(want_oski_bench && guess_blocking_test!=2 /* g.b.t=2 is an extra run*/) 
			{

			rsb__sprintf(oxform,"return BCSR(InputMat, %zd, %zd)",(rsb_printf_int_t)br,(rsb_printf_int_t)bc);
			//rsb__sprintf(oxform,"return BCSR(InputMat, %d, %d)",1,1);
			/* FIXME : ncA and nrA are not enough : we should account for br and bc excess ! */

			Oval = rsb__clone_area(VA,nnz*mtxAp->el_size);
			OIA = rsb__clone_area(IA,nnz*sizeof(rsb_coo_idx_t));
			OJA = rsb__clone_area(JA,nnz*sizeof(rsb_coo_idx_t));

			/* we need duplicates, for we later will use VA as it is */
			if(!Oval || !OIA || !OJA)
			{
				RSB_ERROR("failed aux arrays allocation !\n");goto err;
			}

			/*
				Unfortunately, Oski does not have native BCSR constructors, but 
				rely on conversion from CSR.
				So the measured time is more than it should, but a better
				approximation than oski_CreateMatCSR only.
			*/

			oski_a_t = -rsb_time();
			if(RSB_SOME_ERROR(rsb__allocate_csr_arrays_from_coo_sorted(Oval, OIA, OJA, nnz, nrA, ncA, typecode, &Aval, &Aptr, &Aind)))
			{
				RSB_ERROR("failed csr allocation !\n");goto err;
			}
			oski_a_t += rsb_time();

			if(!Aval || !Aptr || !Aind)
			{
				RSB_ERROR("failed csr arrays allocation !\n");goto err;
			}

			oski_m_t = -rsb_time();
			A_tunable = oski_CreateMatCSR (Aptr, Aind, Aval, nrA, ncA,        /* CSR arrays */
                                // SHARE_INPUTMAT /*COPY_INPUTMAT*/,        /* "copy mode" */
				 /*SHARE_INPUTMAT*/ COPY_INPUTMAT,        /* "copy mode" */
                                 1, INDEX_ZERO_BASED);
				// we should add : INDEX_SORTED, INDEX_UNIQUE
				// 3, INDEX_ZERO_BASED, MAT_TRI_LOWER, MAT_UNIT_DIAG_IMPLICIT);
			oski_m_t += rsb_time();

		        if(A_tunable==INVALID_MAT)
                	{
				RSB_ERROR("invalid oski matrix!\n");goto err;
			}

			oski_t_t = -rsb_time();
			if( oski_ApplyMatTransforms (A_tunable, oxform) )
			{
				RSB_ERROR("invalid transform!\n");goto err;
			}
			oski_t_t += rsb_time();

			if(A_tunable==INVALID_MAT)
			{
				RSB_ERROR("invalid oski tuned matrix!\n");goto err;
			}

				/* FIXME : should error - check these steps */
			//	RSBENCH_STDOUT("# oski : ncA=%zd, nrA=%zd\n",(rsb_printf_int_t)ncA,(rsb_printf_int_t)nrA);
			        x_view = oski_CreateVecView( rhs, ncA, STRIDE_UNIT );
			        y_view = oski_CreateVecView( lhs, nrA, STRIDE_UNIT );
				if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				oski_t = - rsb_time();
				for(i=0;i<times;++i)
				{
#error FIXME: flush breaks measured time
					if(want_inner_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
					/* y <- alpha A * x + beta * y */
					if(oski_MatMult( A_tunable, OP_NORMAL, oalpha, x_view, obeta, y_view ))
					{
							RSB_ERROR("failed uuuoski_MatMult !\n");goto err;
					}
				}
				oski_t += rsb_time();
				if(want_outer_flush == RSB_BOOL_TRUE) RSB_DO_ERROR_CUMULATE(errval,rsb__flush_cache(0));
				if(n_dumpres)
					rsb__debug_print_vector(lhs,RSB_MIN(ndA*nrhs,n_dumpres),typecode,incY);
				/* FIXME */
	
				oski_DestroyMat( A_tunable );
				oski_DestroyVecView( x_view );
				oski_DestroyVecView( y_view );
				RSB_CONDITIONAL_FREE(Aptr);
				RSB_CONDITIONAL_FREE(Aind);
				RSB_CONDITIONAL_FREE(Aval);
				RSB_CONDITIONAL_FREE(Oval);
				RSB_CONDITIONAL_FREE(OJA  );
				RSB_CONDITIONAL_FREE(OIA );
				Aptr= Aind= Aval= NULL;
			} /* want_oski_bench  */
#endif /* RSB_WANT_OSKI_BENCHMARKING */
')dnl
			if(ti>0)
				want_getrow_bench=0;
			if(want_getrow_bench)
			{
				const rsb_coo_idx_t nr=1;
				void * RVA = NULL;
				rsb_coo_idx_t*RIA = NULL;
				rsb_coo_idx_t*RJA = NULL;

				if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc(&RVA,&RIA,&RJA,mtxAp->nc*nr,typecode,RSB_BOOL_FALSE))){goto errgr;}
				for(i=0;i<times;++i)
				{
					rsb_time_t getrow_op_time = RSB_TIME_ZERO;
					rsb_coo_idx_t ri=0;
					rsb_nnz_idx_t rnz=0;
					getrow_op_time = - rsb_time();
					for(ri=0;ri+nr-1<mtxAp->nr;ri+=nr)
						RSB_DO_ERROR_CUMULATE(errval,rsb_mtx_get_coo_block(mtxAp,RVA,RIA,RJA,ri,RSB_MIN(mtxAp->nc-1,ri+nr-1),0,mtxAp->nc-1,NULL,NULL,&rnz,mtxAp->flags));
					getrow_op_time += rsb_time();
					getrow_op_time_best = RSB_MIN_ABOVE_INF(getrow_op_time_best,getrow_op_time,tinf);
					getrow_op_tot_time+=getrow_op_time;
				}
				if(cc==1)getrow_op_time_best_serial=getrow_op_time_best;
errgr:
				RSB_CONDITIONAL_FREE(RVA);
				RSB_CONDITIONAL_FREE(RIA);
				RSB_CONDITIONAL_FREE(RJA);
				if(RSB_SOME_ERROR(errval))
				{goto err;}
			} /* want_getrow_bench */

			if(ti>0)
				want_getdiag_bench=0;
			if(want_getdiag_bench)
			{
				void * DV = rsb__calloc_vector(mtxAp->nr,mtxAp->typecode);
				if(!DV) { errval = RSB_ERR_ENOMEM; goto err; }
				for(i=0;i<times;++i)
				{
					rsb_time_t diag_op_time = RSB_TIME_ZERO;
					diag_op_time = - rsb_time();
					RSB_DO_ERROR_CUMULATE(errval,rsb__dodo_getdiag(mtxAp,DV));
					diag_op_time += rsb_time();
					diag_op_time_best = RSB_MIN_ABOVE_INF(diag_op_time_best,diag_op_time,tinf);
					diag_op_tot_time+=diag_op_time;
				}
				if(cc==1)diag_op_time_best_serial=diag_op_time_best;
				RSB_CONDITIONAL_FREE(DV);
				if(RSB_SOME_ERROR(errval))
				{goto err;}
			} /* want_getdiag_bench */

			if(g_sort_only)
			{
				/* single line output, ideal for benchmark data to be processed later */
				RSBENCH_STDOUT ( "%-20s	%s", rsb__basename(filename),rsb__sprint_matrix_implementation_code2(mtxAp,buf,flags));

				RSBENCH_STDOUT ( "	%.3lf	%lg",
					//raw_Mflops/op_t,	/* please note that in the sort case, it is an absolutely meaningless value */
					true_Mflops/op_t,	/* algorithmic millions of ops per second (not an accurated model)  */
					op_t/true_Mflops	/* the sorting algorithmic constant (not an accurated model) */
				);
			}
			else
			if(!g_estimate_matrix_construction_time)
			if(want_verbose >= -2)
			{
#if RSB_EXPERIMENTAL_WANT_BEST_TIMES
				if( best_t != RSB_CONST_IMPOSSIBLY_BIG_TIME )
				{
					rsb__dump_performance_record(rsb__basename(filename),mtxAp,true_Mflops/best_t,raw_Mflops/best_t,"mop",flags);
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
					if( spsv_spmv_t != RSB_TIME_ZERO )
						RSB_STDOUT ("# (extra) SpMV performance record:\n"),
						rsb__dump_performance_record(rsb__basename(filename),mtxAp,(true_Mflops/best_t)*(tot_t/spsv_spmv_t),raw_Mflops/best_t*(tot_t/spsv_spmv_t),"spmv_uaua*",flags);
')dnl
				}
#else /* RSB_EXPERIMENTAL_WANT_BEST_TIMES */
				rsb__dump_performance_record(rsb__basename(filename),mtxAp,true_Mflops/op_t,raw_Mflops/op_t,"mop",flags);
#endif /* RSB_EXPERIMENTAL_WANT_BEST_TIMES */
			}
			if(g_estimate_matrix_construction_time)
			{
				/* in this case the user asked us too for :
				   * matrix construction Mflops
				   * a ratio of the selected op time with the matrix construction time
				 */
				RSBENCH_STDOUT("\t%.3lg\t%.3lg	", ((double)nnz)/(mct*RSB_REAL_MILLION), mct/op_t);
				rsb__fprint_matrix_implementation_code(mtxAp, "mop", flags, RSB_STDOUT_FD);
				RSBENCH_STDOUT ( "\n");
			}
			omta=((double)rsb_spmv_memory_accessed_bytes(mtxAp));
			
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
			if(want_mkl_bench)
			{
			if(want_mkl_bench_coo)
			{
				RSBENCH_STDOUT ( "#MKL_COO_VS_US-SPMV:%-20s\t%s\t%10.2lf\t%10.2lf\n", rsb__basename(filename),rsb__sprint_matrix_implementation_code2(mtxAp,buf,RSB_FLAG_NOFLAGS),raw_Mflops/(mkl_coo_op_tot_time/times),raw_Mflops/op_t);
				RSBENCH_STDOUT ( "#MKL_COO2CSR2SPMV_VS_US:%-20s\t%s\t%10.2lf\t%10.2lf\n", rsb__basename(filename),rsb__sprint_matrix_implementation_code2(mtxAp,buf,RSB_FLAG_NOFLAGS),(mkl_coo2csr_time)/(mkl_csr_op_tot_time/times),-1.0);
			}
			if(want_mkl_bench_csr)
			{
				RSBENCH_STDOUT ( "#MKL_CSR_VS_US-SPMV:%-20s\t%s\t%10.2lf\t%10.2lf\n", rsb__basename(filename),rsb__sprint_matrix_implementation_code2(mtxAp,buf,RSB_FLAG_NOFLAGS),raw_Mflops/(mkl_csr_op_tot_time/times),raw_Mflops/op_t);
			}
			}
#endif /* RSB_WANT_MKL */
')dnl
#ifdef RSB_WANT_OSKI_BENCHMARKING 
			if(want_oski_bench)
			{
				RSBENCH_STDOUT ( "#OSKI_VS_US-SPMV:%-20s\t%s\t%10.2lf\t%10.2lf\n", rsb__basename(filename),rsb__sprint_matrix_implementation_code2(mtxAp,buf,RSB_FLAG_NOFLAGS),raw_Mflops/(oski_t/times),raw_Mflops/op_t);
				RSBENCH_STDOUT ( "#OSKI_VS_US-ASM~:%-20s\t%s\t%10.2lf\t%10.2lf\n", rsb__basename(filename),rsb__sprint_matrix_implementation_code2(mtxAp,buf,RSB_FLAG_NOFLAGS),oski_m_t+oski_t_t+oski_a_t,mct);
			}
#endif /* RSB_WANT_OSKI_BENCHMARKING  */
			/* WARNING : we cannot use RSB_FLAG_SORTED_INPUT in the recursive case
				     until the following routine will be able to use Z sorted values.. */
			efillin = RSB_REAL_ZERO,eperf = RSB_REAL_ZERO;

			/* FIXME : dies with ct20stif.mtx, now */
			#if 0
dnl			rsb__estimate_expected_fillin_for_blocking is under RSB_OBSOLETE_QUARANTINE_UNUSED
			RSB_WARN("warning : skipping rsb__estimate_expected_fillin_for_blocking\n");
			fet = - rsb_time();
			//rsb__estimate_expected_fillin_for_blocking(VA,IA,JA,nrA,ncA,nnz,typecode,flags/*|RSB_FLAG_SORTED_INPUT*/,br,bc,&efillin);/*TODO:thiscouldbedangerous:fixit!*/
			efillin=mtxAp->einfo.efillin;	/* NEW */
			fet += rsb_time();
			#else /* 0 */
			fet = RSB_TIME_ZERO;
			#endif /* 0 */
			rsb__estimate_expected_raw_performance_for_blocking(nrA,ncA,br,bc,nnz,typecode,flags,efillin,&eperf);

			if(cc==1)
			{
				/* we need input flags, not instantiated matrix flags (which could have not that flag )*/
				if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
					base_best_t=best_t;
				else
					serial_best_t=best_t;
			}
	
			if(want_perf_dump) 
			if(RSB_DO_FLAG_HAS(/*mtxAp->*/flags,RSB_FLAG_QUAD_PARTITIONING))
			{
#if RSB_WANT_MKL
				/* FIXME: this #if is horrible */
dnl				if( /* at_mkl_csr_op_time_best<mkl_csr_op_time_best && */ cc == at_mkl_csr_nt )
dnl					mkl_csr_op_time_best = at_mkl_csr_op_time_best = RSB_MIN( mkl_csr_op_time_best, at_mkl_csr_op_time_best); /* FIXME: this is just a parachute */
				rsb__pr_set(rspr, mtxAp/*NULL */ /* FIXME */, NULL, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti, transA, RSB_CONST_IMPOSSIBLY_BIG_TIME, mkl_csr_op_time_best, RSB_CONST_IMPOSSIBLY_BIG_TIME, at_mkl_csr_op_time_best, RSB_THREADS_AUTO, at_mkl_csr_nt, RSB_CONST_IMPOSSIBLY_BIG_TIME, -1, NULL, NULL, &btpms[1], &btpms[0]);
#endif
			}

			if( want_verbose >= -1 && times > 0 )
			{
				// TODO: update these statistics here.
#if RSB_EXPERIMENTAL_WANT_BEST_TIMES
				RSBENCH_STDOUT ( "#	%10.2lf	%10.2lf	( best, average net performance in %ld tries ); diff:%2.0lf%%\n",
					((double)true_Mflops/best_t), ((double)true_Mflops/op_t),
					(long)times,
					((((double)true_Mflops/best_t)-((double)true_Mflops/op_t))*100)/((double)true_Mflops/op_t)
				);
#endif /* RSB_EXPERIMENTAL_WANT_BEST_TIMES */
				RSBENCH_STDOUT ( "#	%10.2lf	%10.2lf	%10.2lf %10.6lf (min bw, reasonable bw, exceedingly max bw, w/r ratio) (MB/s)\n"
				     "#	%10.2lf (MB per mop) %10.2lf (rhs loads, with a variable degree of locality)\n"
				     "#	%10.2lf (MB per mop, estimated)\n"
				     "#	%10.2lf (assembly + extra to (best) mop time ratio) (%10.2lf s)\n"
				     "#	%10.2lf (assembly (p.e.+s.a.+e.i.+e.s.+...) to mop time ratio)\n"
/*				     "#	%10.2lf (performance estimation to mop time ratio)\n"*/
/*				     "#	%10.2lf (gross fillin estimation to mop time ratio)\n"*/
				     "#	%10.2lf (structure analysis to mop time ratio)\n"
				     "#	%10.2lf (elements insertion to mop time ratio)\n"
				     "#	%10.2lf (elements sorting to mop time ratio) (%10.2lf s)\n"
				     "#	%10.2lf (elements partitioning to mop time ratio)\n"
				     "#	%10.2lf (recursion sort to mop time ratio)\t%10.ld (max recursion depth)\n"
				     "#	%10.2lf	%10.2lf (nnz per row/column)\n"
					,
				((double)rsb_spmv_memory_accessed_bytes_min(mtxAp))*(1.e-6/best_t) ,
				((double)omta)*(1.e-6/best_t) ,
				((double)rsb_spmv_memory_accessed_bytes_max(mtxAp))*(1.e-6/best_t) ,
				((double)rsb_spmv_memory_accessed_bytes_wr_ratio(mtxAp)),
				((double)omta)*(1.e-6),
				(1.0>((fillin*nnz)/(br*ncA))?1.0:((fillin*nnz)/(br*ncA))),
				((double)rsb_spmv_memory_accessed_bytes_(br,bc,nrA,ncA,efillin*nnz,((efillin*nnz)/br)/bc,nrA/br,mtxAp->el_size))*(1.e-6),
				(mct)/(best_t),
				(mtxAp->tat),
				(mtxAp->tat)/(best_t),
/*				(mtxAp->pet)/(best_t),*/
/*				(fet)/(best_t),*/
				(mtxAp->sat)/(best_t),
				(mtxAp->eit)/(best_t),
				(mtxAp->est)/(best_t), (mtxAp->est),
				(mtxAp->cpt)/(best_t),
				((mtxAp->rpt)/(best_t)),((long)rsb__get_recursive_matrix_depth(mtxAp)),
				(double)nnz/nrA, (double)nnz/ncA
				);
			}
				if(RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE>1)
					RSBENCH_STDOUT ( 
					     "#	%10.2lf (estimated fillin)"
					     "#	%10.2lf (estimated fillin error)\n"
					     "#	%10.2lf (estimated raw performance)"
					     "#	%10.2lf (estimated raw performance error)\n"
					     "#	%10.2lf (estimated net performance)"
					     "#	%10.2lf (estimated net performance error)\n",
					efillin, (efillin-fillin)/fillin,
					eperf, (eperf-raw_Mflops/best_t)/(raw_Mflops/best_t),
					efillin?(eperf/efillin):-1,efillin?(((eperf/efillin)-(true_Mflops/best_t))/(true_Mflops/best_t)):-1
					);
dnl
				if(want_verbose >= -1)
				{
					RSBENCH_STDOUT( "#used index storage compared to COO:%zd vs %zd bytes (%.02lf%%) "
						,(size_t)rsb__get_index_storage_amount(mtxAp),(size_t)(sizeof(rsb_coo_idx_t)*2*nnz)
						,(100*(double)rsb__get_index_storage_amount(mtxAp))/RSB_UTIL_COO_IDX_OCCUPATION(mtxAp->nr,mtxAp->nc,mtxAp->nnz)
					);
					RSBENCH_STDOUT( "; compared to CSR:%zd vs %zd bytes (%.02lf%%)\n"
						,(size_t)rsb__get_index_storage_amount(mtxAp),
						 (size_t)(sizeof(rsb_coo_idx_t)*nnz+sizeof(rsb_nnz_idx_t)*(mtxAp->nr+1))
						,(100*(double)rsb__get_index_storage_amount(mtxAp))/RSB_UTIL_CSR_IDX_OCCUPATION(mtxAp->nr,mtxAp->nc,mtxAp->nnz)
					);
				}
dnl				if(0)//very verbose and annoying!
dnl				RSBENCH_STDOUT( "#"),rsb_do_print_nnz_per_row_for_each_submatrix(mtxAp);
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
				totatt += spsv_f_t;
				if( spsv_d_t != RSB_TIME_ZERO)
					RSBENCH_STDOUT( "#gain for spsv if we had infinite spmv-workers:%lf\n",((double)tot_t)/((double)(spsv_d_t)));
				if( spsv_spmv_t != RSB_TIME_ZERO)
					RSBENCH_STDOUT( "#spsv performance vs spmv_uaua*:%lf\n",spsv_spmv_t/tot_t);
				if( spsv_f_t != RSB_TIME_ZERO)
					RSBENCH_STDOUT( "#gain for spsv if we had no concurrent writes preventing locks at all:%lf\n",((double)tot_t)/((double)(spsv_f_t)));
							
')dnl
dnl
dnl
dnl
dnl
dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
			if(ci==0 && smt == RSB_TIME_ZERO && RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
				smt=best_spsv_spmv_t;
			if(ci==cl-1 && pmt == RSB_TIME_ZERO)
				pmt=best_spsv_spmv_t;
			if(ci==0 && sst == RSB_TIME_ZERO && RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
				sst=best_t;
			if(ci==cl-1 && pst == RSB_TIME_ZERO)
				pst=best_t;
')dnl
dnl
			rsb__attr_dump(&attr);
			RSB_BZERO_P((&attr));
dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
dnl
dnl	FIXME: and what about non SPMV ?
dnl
			if(ci==0 && smt == RSB_TIME_ZERO && RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
			{
				smt=best_t;
				sest=mest;
				//sect=mect;
				ssat=msat;
				seit=meit;
				scpt=mcpt;
			}
			if(ci==cl-1 && pmt == RSB_TIME_ZERO)
			{
				pmt=best_t;
			}
')dnl
dnl
dnl
				if(want_verbose > 0 && (RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING)||fn==1))
				{
					rsb_nnz_idx_t minnz=0,maxnz=0,avgnz=0;
					const rsb_bool_t vrpr = (times != 0) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;

					if(vrpr)
					{
						RSBENCH_STDOUT("%%:PERFORMANCE:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.2lf\n",true_Mflops/best_t);
						RSBENCH_STDOUT("\t%le\t%le\n",true_Mflops,best_t);

						RSBENCH_STDOUT("%%:OP_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.6lf\n",best_t);
					}

ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
					if( no_lock_op_time_best != RSB_CONST_IMPOSSIBLY_BIG_TIME )
					{
						RSBENCH_STDOUT("%%:FAKE_LOCK_PERFORMANCE:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.2lf\n",true_Mflops/no_lock_op_time_best);

						RSBENCH_STDOUT("%%:FAKE_LOCK_OP_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.2lf\n",no_lock_op_time_best);

						RSBENCH_STDOUT("%%:FAKE_LOCK_PERF_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.2lf\n",serial_no_lock_op_time_best/no_lock_op_time_best);
					}

					if(qt_op_time_best != RSB_CONST_IMPOSSIBLY_BIG_TIME && cc==1)
					{
						RSBENCH_STDOUT("%%:RECURSIVE_SERIAL_PERFORMANCE:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.2lf\n",true_Mflops/qt_op_time_best);

						RSBENCH_STDOUT("%%:RECURSIVE_SERIAL_OP_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.2lf\n",qt_op_time_best);
					}

')
					if(vrpr)
					{
						if( serial_best_t != RSB_CONST_IMPOSSIBLY_BIG_TIME )
							RSBENCH_STDOUT("%%:PERF_SCALING:"),RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(),
							RSBENCH_STDOUT("\t%10.2lf\n",serial_best_t/best_t);
					}

					RSBENCH_STDOUT("#%%:CONSTRUCTOR_*:SORT	SCAN	INSERT	SCAN+INSERT\n");
					RSBENCH_STDOUT("%%:CONSTRUCTOR_TIMES:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf\n",mest,msat,meit,msat+meit);

					RSBENCH_STDOUT("%%:UNSORTEDCOO2RSB_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n", mest+msat+meit);

					RSBENCH_STDOUT("%%:RSB_SUBDIVISION_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n", msat);

					RSBENCH_STDOUT("%%:RSB_SHUFFLE_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n", meit);

					RSBENCH_STDOUT("%%:ROW_MAJOR_SORT_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n", mest);

					RSBENCH_STDOUT("%%:ROW_MAJOR_SORT_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.3lf\n", sest/mest);

					RSBENCH_STDOUT("%%:SORTEDCOO2RSB_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n", msat+meit);

					RSBENCH_STDOUT("%%:ROW_MAJOR_SORT_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.3lf\n", mest/best_t);

					if(vrpr)
					{
					RSBENCH_STDOUT("%%:CLEANUP_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",mect/best_t);

					RSBENCH_STDOUT("%%:CONSTRUCTOR_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\t%10.2lf\t%10.2lf\t%10.2lf\n",mest/best_t,msat/best_t,meit/best_t,(msat+meit)/best_t);

dnl					RSBENCH_STDOUT("%%:TOTCOO2RSB_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
dnl					RSBENCH_STDOUT("\t%10.2lf\n",(mtat)/best_t);

					RSBENCH_STDOUT("%%:UNSORTEDCOO2RSB_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(msat+meit+mest)/best_t);

					RSBENCH_STDOUT("%%:SORTEDCOO2RSB_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(msat+meit)/best_t);

					RSBENCH_STDOUT("%%:RSB_SUBDIVISION_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(msat)/best_t);

					RSBENCH_STDOUT("%%:RSB_SHUFFLE_TO_MOP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(meit)/best_t);
					}

					RSBENCH_STDOUT("%%:UNSORTEDCOO2RSB_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(ssat+seit+sest)/(msat+meit+mest));

					RSBENCH_STDOUT("%%:SORTEDCOO2RSB_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(ssat+seit)/(msat+meit));

					RSBENCH_STDOUT("%%:RSB_SUBDIVISION_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(ssat)/(msat));

					RSBENCH_STDOUT("%%:RSB_SHUFFLE_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",(seit)/(meit));

					RSBENCH_STDOUT("%%:CONSTRUCTOR_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\t%10.2lf\t%10.2lf\t%10.2lf\n",sest/mest,ssat/msat,seit/meit,(ssat+seit)/(meit+msat));

					if( base_best_t != RSB_CONST_IMPOSSIBLY_BIG_TIME )
					RSBENCH_STDOUT("%%:PERF_SCALING2CSR:"),RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(),
					RSBENCH_STDOUT("\t%10.2lf\n",base_best_t/best_t);

dnl					RSBENCH_STDOUT("%%:FM_FRACTIONS:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
dnl					RSBENCH_STDOUT("\tTODO\n");

					RSBENCH_STDOUT("#%%:SM_COUNTS:	Tot	HalfwordCsr	FullwordCsr	HalfwordCoo	FullwordCoo\n");
					RSBENCH_STDOUT("%%:SM_COUNTS:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					//RSBENCH_STDOUT("\t%d\t%d\t%d\t%d\t%d\n",
					RSBENCH_STDOUT("\t%ld\t%ld\t%ld\t%ld\t%ld\n",
rsb__terminal_recursive_matrix_count(mtxAp),
rsb__terminal_recursive_matrix_count_with_storage_and_flags(mtxAp,RSB_MATRIX_STORAGE_BCSR,RSB_FLAG_USE_HALFWORD_INDICES_CSR),
rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(mtxAp,RSB_MATRIX_STORAGE_BCSR,RSB_FLAG_USE_HALFWORD_INDICES_CSR),
rsb__terminal_recursive_matrix_count_with_storage_and_flags(mtxAp,RSB_MATRIX_STORAGE_BCOR,RSB_FLAG_USE_HALFWORD_INDICES_COO),
rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(mtxAp,RSB_MATRIX_STORAGE_BCOR,RSB_FLAG_USE_HALFWORD_INDICES_COO)
						);

					RSBENCH_STDOUT("%%:SM_IDXOCCUPATIONRSBVSCOOANDCSR:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%zd\t%zd\t%zd\n",rsb__get_index_storage_amount(mtxAp),
						(size_t)RSB_UTIL_COO_IDX_OCCUPATION(mtxAp->nr,mtxAp->nc,mtxAp->nnz),
						(size_t)RSB_UTIL_CSR_IDX_OCCUPATION(mtxAp->nr,mtxAp->nc,mtxAp->nnz)
						);

					RSBENCH_STDOUT("%%:SM_IDXOCCUPATION:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%zd\n",rsb__get_index_storage_amount(mtxAp));

					RSBENCH_STDOUT("%%:SM_MEMTRAFFIC:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.0lf\n",omta);
#if 0
					/* new, elegant */
					RSBENCH_STDOUT("%%:SM_MINMAXAVGSUBMNNZ:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					{
						rsb_submatrix_idx_t i=0;
						rsb_real_t avgnz = ((rsb_real_t)mtxAp->nnz) / mtxAp->all_leaf_matrices_n;
						rsb_coo_idx_t maxnz = 0, minnz = RSB_MAX_MATRIX_NNZ ;

						for(i=0;i<mtxAp->all_leaf_matrices_n;++i)
						{
							struct rsb_mtx_t * submatrix = mtxAp->all_leaf_matrices[i].mtxlp;
							maxnz = RSB_MAX(maxnz,submatrix->nnz);
							minnz = RSB_MIN(minnz,submatrix->nnz);
						}
						RSBENCH_STDOUT(" %zd %zd %.2lf %zd\n",(rsb_printf_int_t)minnz,(rsb_printf_int_t)maxnz,avgnz,(rsb_printf_int_t)mtxAp->all_leaf_matrices_n);
					}
#else
					/* old, obsolete */
					rsb__do_compute_terminal_nnz_min_max_avg_count(mtxAp,&minnz,&maxnz,&avgnz);
					RSBENCH_STDOUT("%%:SM_MINMAXAVGNNZ:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%zd\t%zd\t%zd\n",(rsb_printf_int_t)minnz,(rsb_printf_int_t)maxnz,(rsb_printf_int_t)avgnz);
#endif

				if(want_print_per_subm_stats)
				{
					RSBENCH_STDOUT("%%:SM_NNZ_HISTOGRAM:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					if(!mtxAp->all_leaf_matrices)
						RSBENCH_STDOUT(" %zd\n",(size_t)mtxAp->nnz);
					else
					{
						rsb_submatrix_idx_t i=0;
						for(i=0;i<mtxAp->all_leaf_matrices_n;++i)
							RSBENCH_STDOUT(" %zd",(size_t)mtxAp->all_leaf_matrices[i].mtxlp->nnz);
						RSBENCH_STDOUT("\n");
					}

					RSBENCH_STDOUT("%%:SM_NNZ_PER_ROW:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					if(!mtxAp->all_leaf_matrices)
						RSBENCH_STDOUT(" %lf\n",((double)mtxAp->nnz)/mtxAp->nr);
					else
					{
						rsb_submatrix_idx_t i=0;
						for(i=0;i<mtxAp->all_leaf_matrices_n;++i)
							RSBENCH_STDOUT(" %.2lf",((double)mtxAp->all_leaf_matrices[i].mtxlp->nnz)/mtxAp->all_leaf_matrices[i].mtxlp->nr);
						RSBENCH_STDOUT("\n");
					}
				} /* want_print_per_subm_stats */

#ifdef RSB_WANT_PERFORMANCE_COUNTERS
			if(want_perf_counters)
				{
					int i;
					for(i=0;i<rsb_pci.eventnum;++i)
					{
						RSBENCH_STDOUT("%%:RSB_%s:",rsb_pci.eventdesc[i]);
						RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%zd\n",(size_t)(rsb_pci.eventvals[i]));
					}
				} /* want_perf_counters */
#endif /* RSB_WANT_PERFORMANCE_COUNTERS */
				}
dnl
			} /* times */
dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
#if RSB_WANT_MKL
				if(want_mkl_bench) /* 20110428 */
				if(want_verbose > 0 && (RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING)||fn==1))
				{
#ifdef mkl_get_version
					MKLVersion mv;
					mkl_get_version(&mv);
					RSBENCH_STDOUT("#%%:MKL %d.%d-%d, %s, %s, %s, %s\n",mv.MajorVersion,mv.MinorVersion,mv.UpdateVersion,mv.ProductStatus,mv.Build,mv.Processor,mv.Platform);
#else /* mkl_get_version */
					RSBENCH_STDOUT("#%%:MKL, version unknown\n");
#endif /* mkl_get_version */
			if(want_mkl_bench_coo)
			{
					RSBENCH_STDOUT("%%:MKL_COO_PERFORMANCE:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",true_Mflops/mkl_coo_op_time_best);

					RSBENCH_STDOUT("%%:MKL_COO_OP_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(); RSBENCH_STDOUT("\t%10.6lf\n",mkl_coo_op_time_best);

					if( mkl_coo_op_time_best_serial != RSB_CONST_IMPOSSIBLY_BIG_TIME )
					RSBENCH_STDOUT("%%:MKL_COO_PERF_SCALING:"),RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(),
					RSBENCH_STDOUT("\t%10.2lf\n",mkl_coo_op_time_best_serial/mkl_coo_op_time_best);
			}
#ifdef RSB_WANT_PERFORMANCE_COUNTERS
			if(want_perf_counters)
				{
					int i;
					for(i=0;i<mkl_csr_pci.eventnum;++i)
					{
						RSBENCH_STDOUT("%%:MKL_CSR_%s:",mkl_csr_pci.eventdesc[i]);
						RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%zd\n",mkl_csr_pci.eventvals[i]);
					}
					if(want_mkl_bench_coo)
					for(i=0;i<mkl_coo_pci.eventnum;++i)
					{
						RSBENCH_STDOUT("%%:MKL_COO_%s:",mkl_coo_pci.eventdesc[i]);
						RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%zd\n",mkl_coo_pci.eventvals[i]);
					}
				}
#endif /* RSB_WANT_PERFORMANCE_COUNTERS */
			if(want_mkl_bench_csr)
			{
					RSBENCH_STDOUT("%%:MKL_CSR_PERFORMANCE:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",true_Mflops/mkl_csr_op_time_best);

					RSBENCH_STDOUT("%%:MKL_CSR_OP_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n",mkl_csr_op_time_best);

					if( mkl_csr_op_time_best_serial != RSB_CONST_IMPOSSIBLY_BIG_TIME )
					RSBENCH_STDOUT("%%:MKL_CSR_PERF_SCALING:"),RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(),
					RSBENCH_STDOUT("\t%10.2lf\n",mkl_csr_op_time_best_serial/mkl_csr_op_time_best);
			}
			if(want_mkl_bench_gem)
			{
					RSBENCH_STDOUT("%%:MKL_GEMV_PERFORMANCE:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",true_gem_Mflops/mkl_gem_op_time_best);

					RSBENCH_STDOUT("%%:MKL_GEMV_OP_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n",mkl_gem_op_time_best);

					if( mkl_gem_op_time_best_serial != RSB_CONST_IMPOSSIBLY_BIG_TIME )
					RSBENCH_STDOUT("%%:MKL_GEMV_PERF_SCALING:"),RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH(),
					RSBENCH_STDOUT("\t%10.2lf\n",mkl_gem_op_time_best_serial/mkl_gem_op_time_best);
			}

					if( mkl_coo2csr_time != RSB_TIME_ZERO )
					{
						RSBENCH_STDOUT("%%:MKL_COO2CSR_T0_CSR_TIME:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.6lf\n",mkl_coo2csr_time);
						RSBENCH_STDOUT("%%:MKL_COO2CSR_T0_CSR_OP:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.2lf\n",mkl_coo2csr_time/mkl_csr_op_time_best);


						RSBENCH_STDOUT("%%:SORTEDCOO2RSB_VS_MKLCOO2CSR:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
						RSBENCH_STDOUT("\t%10.3lf\n", (msat+meit)/(mkl_coo2csr_time));
					}
				} /* want_mkl_bench */
#endif /* RSB_WANT_MKL */
')dnl
dnl
dnl
				if(want_getrow_bench)
				{
					const char*norsbnotice="";
					const char*rsbnotice="NORSB_";
					const char*notice=norsbnotice;

					if(want_verbose > 0 && (RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING)||fn==1))
						{}
					else
						notice = rsbnotice;

					RSBENCH_STDOUT("%%:%sGETROW_PERFORMANCE:",notice);RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",((rsb_time_t)mtxAp->nnz)/(RSB_REAL_MILLION*getrow_op_time_best));
					RSBENCH_STDOUT("%%:%sGETROW_OP_TIME:",notice);RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n",getrow_op_time_best);
					RSBENCH_STDOUT("%%:%sGETROW_TO_SPMV_OP_TIME:",notice);RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n",getrow_op_time_best/best_t);
dnl					RSBENCH_STDOUT("%%:GETROW_PERF_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
dnl					RSBENCH_STDOUT("\t%10.2lf\n",/getrow_op_time_best);

				}
dnl
dnl
				if(want_getdiag_bench)
				{
					const char*norsbnotice="";
					const char*rsbnotice="NORSB_";
					const char*notice=norsbnotice;
					if(want_verbose > 0 && (RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING)||fn==1))
						{}
					else
						notice = rsbnotice;

					RSBENCH_STDOUT("%%:%sGETDIAG_PERFORMANCE:",notice);RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.2lf\n",((rsb_time_t)mtxAp->nr)/(RSB_REAL_MILLION*diag_op_time_best));
					RSBENCH_STDOUT("%%:%sGETDIAG_OP_TIME:",notice);RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n",diag_op_time_best);
					RSBENCH_STDOUT("%%:%sGETDIAG_TO_SPMV_OP_TIME:",notice);RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
					RSBENCH_STDOUT("\t%10.6lf\n",diag_op_time_best/best_t);
dnl					RSBENCH_STDOUT("%%:GETDIAG_PERF_SCALING:");RSBENCH_STDOUT_MATRIX_ESSENTIALS_RSBENCH();
dnl					RSBENCH_STDOUT("\t%10.2lf\n",/diag_op_time_best);

				}
dnl
dnl
				if(want_verbose > -2)
					RSBENCH_STDOUT( "#\n");/* end of record */
dnl
dnl
				if(guess_blocking_test)
				{
					rsb_flags_t oflags = RSB_FLAG_NOFLAGS;
					/* TODO : should keep info of the worst, to */
					rsb_perf_t nrp=(true_Mflops/op_t),bomta = RSB_REAL_ZERO /* best op memory traffic amount */;

					if(guess_blocking_test==1)
					{
						if( nrp>RSB_REAL_ZERO && nrp>bperf)
						{
							bperf=nrp;
							bomta=omta;
							bfillin=fillin;
							ebfillin=efillin;
							bri=brvi;
							bci=bcvi;
						}
					
						if(brv[brvi]==1 && bcv[bcvi]==1)/* IF ANY! */
						{
							cperf=nrp;
						}
 
						if((nrp>RSB_REAL_ZERO && nrp<wperf) || wperf == RSB_REAL_ZERO)
						{
							wperf=nrp;
						}

						if( fillin > maxfillin )
						{
							maxfillin=fillin;
						}
					}

					if( guess_blocking_test==2) 
					{
						egfillin=efillin;
						RSBENCH_STDOUT("# GUESS DATA;  best performance was       :	%zd	%zd\n", (size_t)brv[bri], (size_t)bcv[bci] );
						RSBENCH_STDOUT("# GUESS DATA;  guessed was                :	%zd	%zd\n", (size_t)br, (size_t)bc );
						RSBENCH_STDOUT("# GUESS DATA:  performance diff from best :	%lg\n", (nrp-bperf)/bperf );
						RSBENCH_STDOUT("# GUESS DATA:  performance diff from worst:	%lg\n", (nrp-wperf)/wperf );
						if(cperf)
						RSBENCH_STDOUT("# GUESS DATA:  performance diff over CSR:	%lg\n", (nrp-cperf)/cperf );
						RSBENCH_STDOUT("# GUESS DATA:  best/guessed op matrix traffic amount:	%lg	%lg\n", bomta,omta);
						RSBENCH_STDOUT("#GUESS_TEST_:%-20s\t%20s\t%zd\t%zd\t%zd\t%zd\t%zd\t%zd\n",
							rsb__basename(filename),
							rsb__sprint_matrix_implementation_code2(mtxAp,buf,flags),
				(rsb_printf_int_t)((nrp>=bperf*.95) || (brv[bri]==br && bcv[bci]==bc)),	/* (fuzzy WIN) */
				(rsb_printf_int_t)((nrp>=bperf) || (brv[bri]==br && bcv[bci]==bc)),	/* if 1, best blocking guess (WIN) */
				(rsb_printf_int_t)(nrp>=bperf),			/* if 1, best performance guess */
				(rsb_printf_int_t)(brv[bri]==br && bcv[bci]==bc),	/* if 1, best blocking guess */
				(rsb_printf_int_t)(nrp>=cperf),	/* if 0, we lose over (our) plain CSR  */
				(rsb_printf_int_t)(nrp> wperf)	/* if 0, we performed as the worst blocking! */
							);
					flags=oflags;

					RSBENCH_STDOUT(	"#GUESS_TEST:%-20s\t%-20s"
						"\t%10.2lf"
						"\t%10.2lf"
						"\t%zd" "\t%zd"
						"\t%10.4lf" "\t%10.2lf" "\t%10.4lf" "\t%10.2lf" "\t%10.4lf" "\n"
						,
						rsb__basename(filename),
						rsb__sprint_matrix_implementation_code2(mtxAp,buf,flags),	
						/* grmflops */
						raw_Mflops/op_t,
						/* egfillin */
						egfillin,
						/* bbr */
						(rsb_printf_int_t)brv[bri],
						/* bbc */
						(rsb_printf_int_t)bcv[bci],
						/* bfillin */
						bfillin,
						/* brmflops */
						bperf*bfillin,
						/* ebfillin */
						ebfillin,
						/* csrmflops */
						cperf,
						/* maxfillin */
						maxfillin);

						flags=oflags;
					}
				

					if(brvi==brl-1 && bcvi==bcl-1 && guess_blocking_test==1)
					{
						oflags=flags;
						RSB_DO_FLAG_ADD(flags,RSB_FLAG_AUTO_BLOCKING);
						guess_blocking_test++;
						--bcvi;	/* un altro giro :) */
					}
				} /* guess_blocking_test */
dnl
dnl
')dnl
dnl
')dnl
dnl
		erri:
			if(want_in_place_assembly && mtxAp)
			{
dnl				rsb_coo_sort(VA,IA,JA,mtxAp->nnz,nrA,ncA,typecode,RSB_FLAG_NOFLAGS);
				rsb_time_t st = -rsb_time();
				errval = rsb_mtx_switch_to_coo(mtxAp,&VA,&IA,&JA,RSB_FLAG_SORTED_INPUT),mtxAp=NULL;
				st += rsb_time();
				RSBENCH_STDOUT("# rsb_mtx_switch_to_coo time: %lg.\n",st);
				if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES); }
			}
			RSB_MTX_FREE(mtxAp);
			RSB_CONDITIONAL_FREE(lhs);
			RSB_CONDITIONAL_FREE(rhs);

			RSB_CONDITIONAL_FREE(p_r);
			RSB_CONDITIONAL_FREE(p_c);
			
			if(RSB_SOME_ERROR(errval))
			{
				RSB_ERROR(RSB_ERRM_ES);goto err;
			}
			if(brl==0 || bcl==0) break;
		} /* ci : core (count) index */

			if(want_verbose > 0)
			{
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
            			RSBENCH_STDOUT("%%operation:matrix	CONSTRUCTOR[%d]	SPMV[%d]	SPMV[%d]	STSV[%d]	STSV[%d]\n",
					ca[0], ca[0], ca[cl-1], ca[0], ca[cl-1]);
            			RSBENCH_STDOUT("%%operation:%s	%lg	%lg	%lg	%lg	%lg\n",
					rsb__basename(filename),sct,smt,pmt,sst,pst);
')dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
            			RSBENCH_STDOUT("%%operation:matrix	CONSTRUCTOR[%d]	SPMV[%d]	SPMV[%d]\n",ca[0],ca[0],ca[cl-1]);
            			RSBENCH_STDOUT("%%operation:%s	%lg	%lg	%lg\n",
					rsb__basename(filename),sct,smt,pmt);
')dnl
            			RSBENCH_STDOUT("%%constructor:matrix	SORT[%d]	SCAN[%d]	SHUFFLE[%d]	INSERT[%d]\n",
					ca[0],ca[0],ca[0],ca[0]);
ifelse(mop,`mat_stats',`',`dnl
            			RSBENCH_STDOUT("%%constructor:%s	%lg	%lg	%lg	%lg\n",
					rsb__basename(filename),sest,ssat,scpt,seit);
')dnl
			}
		} /* ti (transposition index) */
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	}	/* incYi */
	}	/* incXi */
')dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	}	/* betai */
	}	/* alphai */
')dnl
	}	/* nrhsi */
	}

dnl
	if(want_verbose >= -2)
	{
 		RSBENCH_STDOUT("# so far, program took %.3lfs of wall clock time; ancillary tests %.3lfs; I/O %.3lfs; checks %.3lfs; conversions %.3lfs; rsb/mkl tuning %.3lfs/%.3lfs ",totprt + rsb_time(),totatt,totiot,totht,totct,tottt,totmt);
		RSBENCH_STDOUT(".\n"); /* FIXME: this takes too much space here ! */
	}

	if(want_verbose >= 0)
		rsb__getrusage();
dnl
done:
dnl
ifelse(mop,`mat_stats',`dnl
	RSB_CONDITIONAL_FREE(nnzs);
	RSB_CONDITIONAL_FREE(element_count );
	RSB_CONDITIONAL_FREE(block_count   );
')dnl
dnl
frv:
dnl
	if( !should_recycle_io )
	{
		RSBENCH_STDOUT("# Freeing I/O arrays.\n");
		RSB_CONDITIONAL_FREE(IA);
		RSB_CONDITIONAL_FREE(JA);
		RSB_CONDITIONAL_FREE(VA);
	}
	
	if(mtxAp && !should_recycle_matrix){RSB_MTX_FREE(mtxAp)}
	if( brv != rua ) {RSB_CONDITIONAL_FREE(brv);}
	if( bcv != cua ) {RSB_CONDITIONAL_FREE(bcv);}
dnl
		RSBENCH_MAY_SQUIT(ret,{}) /* early end of program */
		RSBENCH_MAY_TQUIT(ret,{}) /* early end of program */
dnl
	}	/* typecodesi */
	}	/* diagi */
	}	/* asfii */
	}	/* asiii */
	}	/* accii */
	}	/* areci */
nfnm:	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	}	/* filenamei */
dnl
	RSBENCH_STDOUT("# benchmarking terminated --- finalizing run.\n");
dnl
#if RSB_WANT_PERFORMANCE_COUNTERS_IN_RSBENCH 
	errval = rsb_perf_counters_finalize();
	if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES); }
#endif
dnl
ret:
dnl
	errval = RSB_ERR_NO_ERROR;
	goto rret;
dnl
err:
dnl
	rsb_perror(NULL,errval);
	errval = RSB_ERR_GENERIC_ERROR;
dnl
rret:
dnl
ifelse(mop,`mat_stats',`dnl
	RSB_CONDITIONAL_FREE(nnzs);
	RSB_CONDITIONAL_FREE(element_count );
	RSB_CONDITIONAL_FREE(block_count   );
')dnl
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(VA);
dnl	if(want_in_place_assembly && mtxAp)rsb_coo_sort(VA,IA,JA,mtxAp->nnz,nrA,ncA,typecode,RSB_FLAG_NOFLAGS);
	if(want_in_place_assembly && mtxAp)rsb_mtx_switch_to_coo(mtxAp,&VA,&IA,&JA,RSB_FLAG_SORTED_INPUT),mtxAp=NULL;
	RSB_MTX_FREE(mtxAp);
	if( brv != rua ) {RSB_CONDITIONAL_FREE(brv);}
	if( bcv != cua ) {RSB_CONDITIONAL_FREE(bcv);}
dnl
	if(want_perf_dump) 
	{
		RSBENCH_STDOUT("# ====== BEGIN Total summary record.\n");
		errval = rsb__pr_dump(rspr, filenamea, ca, incXa, incYa, nrhsa, typecodes, NULL, NULL );
		RSBENCH_STDOUT("# ======  END  Total summary record.\n");
		if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(tret,RSB_ERRM_ES);
		errval = rsb__pr_save(fprfn, rspr, filenamea, ca, incXa, incYa, nrhsa, typecodes, NULL, RSB_BOOL_TRUE);
		if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(tret,RSB_ERRM_ES);
		RSBENCH_STDOUT("# Removing the temporary record file %s.\n",cprfn);
		remove(cprfn);
	}
	if( ca  != ca_ ) {RSB_CONDITIONAL_FREE(ca);}
#if RSB_RSBENCH_STATIC_FILENAMEA
	{
		rsb_int_t filenamei;
		for ( filenamei = 0; filenamei < filenamen ; ++filenamei )
	       		if( filenamea[filenamei] != fnbufp[0] )
				free(filenamea[filenamei]); /* allocated in rsb__strdup */
	}
#else
	/* if(filenamea!=&fnbufp)RSB_CONDITIONAL_FREE(filenamea); */
	if(filenamea!=&fnbufp)free(filenamea); /* FIXME */
#endif
	if(nrhsa!=(&nrhs))RSB_CONDITIONAL_FREE(nrhsa); /* if same, nrhsa only points to incx */
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	if(incXa!=(&incX))RSB_CONDITIONAL_FREE(incXa); /* same trick */
 	if(incYa!=(&incY))RSB_CONDITIONAL_FREE(incYa); /* same trick */
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	RSB_CONDITIONAL_FREE(alphaip);
	RSB_CONDITIONAL_FREE(betaip);
')dnl
dnl
	if(want_likwid == RSB_BOOL_TRUE){RSB_LIKWID_MARKER_EXIT;} /* FIXME: and other cases ? */
dnl
	RSBENCH_MEM_ALLOC_INFO("\n");
	if(want_verbose > 0)
		rsb__echo_timeandlabel(" terminating run at ","\n",&st);
dnl
	rsb__mbw_es_free(mbetp);
dnl
	rsb__pr_free(rspr);
dnl
tret:
	if(RSB_SOME_ERROR(rsb_lib_exit(RSB_NULL_EXIT_OPTIONS)))
		return RSB_ERR_GENERIC_ERROR;
	return errval;
#else /* RSB_HAVE_GETOPT_H */
	return RSB_ERR_NO_ERROR;
#endif /* RSB_HAVE_GETOPT_H */
}
')dnl
dnl popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl
dnl define(`RSB_M4_MATRIX_META_OPS',(RSB_M4_LIST_PUSH_BACK(RSB_M4_MATRIX_META_OPS,`matrix_stats')))dnl
dnl
/* one function for each of RSB_M4_MATRIX_META_OPS_REDUCED*/
foreach(`mop',RSB_M4_MATRIX_META_OPS_REDUCED,`dnl
RSB_M4_COMPLETE_TEST_PROGRAM_FUNCTION(mop)
')dnl
dnl

dnl
#ifdef __cplusplus
}
#endif  /* __cplusplus */
ifdef(`ONLY_WANT_HEADERS',`
#endif	/* RSB_TEST_MATOPS_H_INCLUDED */
')
/* @endcond */
dnl
