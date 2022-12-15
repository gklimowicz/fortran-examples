/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains MKL interfacing functions.
 * */
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
include(`libspblas_macros.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_RSB_MKL_H_INCLUDED
#define RSB_RSB_MKL_H_INCLUDED
#include "rsb_internals.h"
#if RSB_WANT_MKL
#ifdef RSB_INCLUDE_MKL_MKL_H
#include <mkl/mkl.h>
#include <mkl/mkl_blas.h>	/* dgemm, ... */
#include <mkl/mkl_spblas.h>
#else
#include <mkl.h>
#include <mkl_blas.h>	/* dgemm, ... */
#include <mkl_spblas.h>
#endif
/* #include <mkl_types.h> */
/* #include <mkl_service.h> */ /* mkl_get_version */
dnl
#if defined(RSB_WANT_MKL) && defined(INTEL_MKL_VERSION) && ( (INTEL_MKL_VERSION) >= 110302 )
#define RSB_WANT_MKL_INSPECTOR 1
#endif /* RSB_WANT_MKL */
dnl
',`dnl
#include "rsb_mkl.h"
#if RSB_WANT_MKL
RSB_INTERNALS_COMMON_HEAD_DECLS
')
dnl

#include "rsb_tune.h" /* rsb_tattr_t */
dnl
define(`RSB_M4_RSB_TYPE_TO_MKL_TYPE',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(type,`double complex',`MKL_Complex16')`'dnl
ifelse(type,`float complex',`MKL_Complex8')`'dnl
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
define(`RSB_M4_RSB_TYPE_TO_MKL_REAL_TYPE',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(type,`double complex',`double')`'dnl
ifelse(type,`float complex',`float')`'dnl
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
define(`RSB_M4_MKL_ONE',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(mtype)),1,`dnl
(RSB_M4_RSB_TYPE_TO_MKL_TYPE(type)) {1.0,0.0} `'dnl
',`dnl
(RSB_M4_RSB_TYPE_TO_MKL_TYPE(type))(1.0)`'dnl
')dnl
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl

rsb_err_t rsb__mkl_gemv(rsb_type_t typecode, const void * Mp, const void*Bp, void*Xp, rsb_nnz_idx_t mdim, rsb_coo_idx_t vdim, rsb_coo_idx_t*udimp)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	/* FIXME: TODO: incX != 1 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const MKL_INT dim=(rsb_coo_idx_t)sqrt((double)mdim);
	const MKL_INT incX=1;
	char transA_mkl=110;
dnl ,transB=110;
	if(!Mp || !Xp || !Bp)
		goto err;
	if(dim<1 || dim>vdim)
		goto err;
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype alpha=RSB_M4_ONE(mtype), beta=RSB_M4_ONE(mtype);
		`'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`gemv'(&transA_mkl,&dim,&dim, (RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)(&alpha),(const RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)Mp,&dim,(const RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)Bp,&incX,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)&beta,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)Xp,&incX);
dnl		`'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`gemm'(&transA_mkl,&transB,&dim,&dim,&dim,&alpha,(const mtype*)Mp,&dim,(const mtype*)Bp,&dim,&beta,(mtype*)Xp,&dim);
	}
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	else 
')dnl
		errval=RSB_ERR_BADARGS;

	if(udimp)
		*udimp=dim;
err:
	return errval;
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static char rsb_rsb_to_mkl_trans(rsb_trans_t transA_mkl)
{
	/**
	 * \ingroup gr_internals
	 */
	switch(transA_mkl)
	{
		case(RSB_TRANSPOSITION_N):
		return singlequote(n);
		break;
		case(RSB_TRANSPOSITION_T):
		return singlequote(t);
		break;
		case(RSB_TRANSPOSITION_C):
		return singlequote(c);
		break;
		default:
		return singlequote(n);	// FIXME
	}
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static char rsb_rsb_to_mkl_sym(rsb_flags_t flags)
{
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC))
		return singlequote(s);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_TRIANGULAR))
		return singlequote(t);
	else
		return singlequote(g);
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static char rsb_rsb_to_mkl_upl(rsb_flags_t flags)
{
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER))
		return singlequote(l);
	else
		return singlequote(u);
}
')dnl
dnl

dnl
rsb_err_t rsb__mkl_coo_spmv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = singlequote(n); // not unit diagonal
	matdescra[3] = singlequote(c); // zero based indexing

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`coomv'(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap,matdescra,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)x,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__mkl_coo_spmm(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nrhs, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = singlequote(n); // not unit diagonal
	matdescra[3] = singlequote(c); // zero based indexing

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`coomm'(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(MKL_INT*)(&k),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap,matdescra,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)b,(MKL_INT*)(&ldb),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)c,(MKL_INT*)(&ldc));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__mkl_coo_spsv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = singlequote(n); // not unit diagonal
	matdescra[3] = singlequote(c); // zero based indexing
	/* 20101118	MKL 9.1 reference manual declares also k among the parameters */
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`coosv'(&transA_mkl,(MKL_INT*)(&m),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap,matdescra,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)x,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_mkl_csr_spmv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = singlequote(n); // not unit diagonal
	matdescra[3] = singlequote(c); // zero based indexing

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`csrmv'(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap,matdescra,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)x,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
/* The following three look weird, I know. */
#define RSB_GET_MKL_MAX_THREADS rsb__set_num_threads(RSB_THREADS_GET_MAX_SYS)
#define RSB_GET_MKL_BASE_THREADS 1 /* FIXME: no mkl_get_num_threads */
#define RSB_GET_MKL_DEFAULT_THREADS mkl_get_max_threads() /* omp_get_num_threads(); */
#define RSB_MKL_MAX_AT_TIME 1.0 /* FIXME */

#define RSB_MKL_SET_THREADS_RANGE(lnt,unt,otnp)				\
		if(RSB_DT_SAME_THREADS_TNP(otnp))			\
			lnt = unt = 0; 					\
		else							\
		{							\
			if(RSB_DT_THREADS_TUNE_TNP(otnp))		\
				; /* ok */				\
			else						\
				if(RSB_DT_SPEC_THREADS_TNP(otnp))	\
			lnt = unt = *otnp;				\
		}

#define RSB_MKL_THREADS_TUNING_ODECLS					\
		const rsb_time_t tinf = RSB_CACHED_TIMER_GRANULARITY;	\
		rsb_time_t best = RSB_CONST_IMPOSSIBLY_BIG_TIME;	\
		rsb_thread_t ont = RSB_GET_MKL_BASE_THREADS;		\
		rsb_thread_t nt, lnt = 1, unt = RSB_GET_MKL_MAX_THREADS;\
		rsb_thread_t otn = ont;					\
		rsb_thread_t dtn = RSB_GET_MKL_DEFAULT_THREADS;


#define RSB_MKL_THREADS_TUNING_IDECLS									\
			rsb_time_t it = rsb_time(), ct = RSB_TIME_ZERO;	/* initial/current time */	\
			rsb_time_t dt = it, tt = RSB_TIME_ZERO; /* elapsed (delta) / total  time */	\
			rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO; /* best / worst  time */	\
			rsb_time_t ss = RSB_TIME_ZERO; /* sum of squares */				\
			rsb_time_t mint = RSB_TIME_ZERO; /* minimal time */				\
			rsb_int_t times = 0; \
			const rsb_int_t mintimes = RSB_CONST_AT_OP_SAMPLES_MIN, maxtimes = RSB_CONST_AT_OP_SAMPLES_MAX;	\
			rsb_time_t maxt = RSB_AT_MAX_TIME/* RSB_MKL_MAX_AT_TIME*/;
')dnl
dnl

dnl
rsb_err_t rsb__mkl_csr_spmv_bench(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t dt, tpo;

	if(otnp)
	{
		RSB_MKL_THREADS_TUNING_ODECLS
		RSB_MKL_SET_THREADS_RANGE(lnt,unt,otnp)

		for(nt=lnt;nt<=unt;++nt)
		{
			RSB_MKL_THREADS_TUNING_IDECLS
			if(nt) mkl_set_num_threads(nt);

			do
			{
				errval = rsb__do_mkl_csr_spmv(VA, m, k, nnz, IP, JA, x, y, alphap, betap, transA, typecode, flags);
				RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,tinf,times);
			}
			while(RSB_REPEAT(ct-it,times,mint,mintimes,maxt,maxtimes));

			dt = bt;
			if(dt < best )
			{
				otn = nt;
				best = RSB_MIN_ABOVE_INF(best,dt,tinf);
				RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp);
			}
			rsb__tattr_sets(ttrp,dtn,nt,dt,otn,times);/* FIXME: if no threads tuning, shall set dtpo = btpo, as well as ttrp.optt=0 */
			if(dtn == nt) RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp+1);
		}
		mkl_set_num_threads(ont);
		ttrp->ttt += rsb_time(); /* ttrp->ttt = tt; */
		tpo = best; /* tpo = 1.0 / best; */
		*otnp = otn;
	}
	else
	{
		dt = -rsb_time();
		errval = rsb__do_mkl_csr_spmv(VA, m, k, nnz, IP, JA, x, y, alphap, betap, transA, typecode, flags);
		dt += rsb_time();
		/* tpo = 1.0 / dt; */
		tpo = dt;
	}
	if(tpop)
		*tpop = tpo;
	return errval;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_mkl_csr_spmm(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	MKL_INT ldb_ = n, ldc_ = n; /* for zero based indexing */

	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = singlequote(n); // not unit diagonal
	matdescra[3] = singlequote(c); // zero based indexing

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		ldb_ = k, ldc_ = m, /* for one based indexing */
		matdescra[3] = singlequote(f); // one based indexing

	#if 1
	/* n = nrhs */
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`csrmm'(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&n),(MKL_INT*)(&k),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap,matdescra,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)b,(MKL_INT*)(&ldb_),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)c,(MKL_INT*)(&ldc_));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
#endif /* 1 */
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__mkl_csr_spmm_bench(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t dt, tpo;

	if(otnp)
	{
		RSB_MKL_THREADS_TUNING_ODECLS
		RSB_MKL_SET_THREADS_RANGE(lnt,unt,otnp)

		for(nt=lnt;nt<=unt;++nt)
		{
			RSB_MKL_THREADS_TUNING_IDECLS
			if(nt) mkl_set_num_threads(nt);

			do
			{
				errval = rsb__do_mkl_csr_spmm(VA, m, k, n, nnz, IP, JA, b, ldb, c, ldc, alphap, betap, transA, typecode, flags);
				RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,tinf,times);
			}
			while(RSB_REPEAT(ct-it,times,mint,mintimes,maxt,maxtimes));

			dt = bt;
			if(dt < best )
			{
				otn = nt;
				best = RSB_MIN_ABOVE_INF(best,dt,tinf);
				RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp);
			}
			rsb__tattr_sets(ttrp,dtn,nt,dt,otn,times);/* FIXME: if no threads tuning, shall set dtpo = btpo, as well as ttrp.optt=0 */
			if(dtn == nt) RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp+1);
		}
		mkl_set_num_threads(ont);
		ttrp->ttt += rsb_time(); /* ttrp->ttt = tt; */
		tpo = best; /* tpo = 1.0 / best; */
		*otnp = otn;
	}
	else
	{
		dt = -rsb_time();
		errval = rsb__do_mkl_csr_spmm(VA, m, k, n, nnz, IP, JA, b, ldb, c, ldc, alphap, betap, transA, typecode, flags);
		dt += rsb_time();
		/* tpo = 1.0 / dt; */
		tpo = dt;
	}
	if(tpop)
		*tpop = tpo;
	return errval;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_mkl_csr_spsv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = singlequote(n); // not unit diagonal
	matdescra[3] = singlequote(c); // zero based indexing

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`csrsv'(&transA_mkl,(MKL_INT*)(&m),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap,matdescra,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)x,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_mkl_csr_spsm(const void *VA, const MKL_INT m, const MKL_INT nrhs, const MKL_INT * IP, const MKL_INT *JA, const void * b, void * c, const void *alphap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, const MKL_INT ldb, const MKL_INT ldc)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = singlequote(n); // not unit diagonal
	matdescra[3] = singlequote(c); // zero based indexing

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`csrsm'(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap,matdescra,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)b,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(int)*)&ldb,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)c,(RSB_M4_RSB_TYPE_TO_MKL_TYPE(int)*)&ldc);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__mkl_csr_spsv_bench(const void *VA, const MKL_INT m, const MKL_INT k/*, const MKL_INT n*/, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, /*const MKL_INT ldb,*/ void * c,/* const MKL_INT ldc,*/ const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t dt, tpo;

	if(otnp)
	{
		RSB_MKL_THREADS_TUNING_ODECLS
		RSB_MKL_SET_THREADS_RANGE(lnt,unt,otnp)

		for(nt=lnt;nt<=unt;++nt)
		{
			RSB_MKL_THREADS_TUNING_IDECLS
			if(nt) mkl_set_num_threads(nt);

			do
			{
				errval = rsb__do_mkl_csr_spsv(VA, m, k, nnz, IP, JA, b, c, alphap, betap, transA, typecode, flags);
				RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,tinf,times);
			}
			while(RSB_REPEAT(ct-it,times,mint,mintimes,maxt,maxtimes));

			dt = bt;
			if(dt < best )
			{
				otn = nt;
				best = RSB_MIN_ABOVE_INF(best,dt,tinf);
				RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp);
			}
			rsb__tattr_sets(ttrp,dtn,nt,dt,otn,times);/* FIXME: if no threads tuning, shall set dtpo = btpo, as well as ttrp.optt=0 */
			if(dtn == nt) RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp+1);
		}
		mkl_set_num_threads(ont);
		ttrp->ttt += rsb_time(); /* ttrp->ttt = tt; */
		tpo = best; /* tpo = 1.0 / best; */
		*otnp = otn;
	}
	else
	{
		dt = -rsb_time();
		errval = rsb__do_mkl_csr_spsv(VA, m, k, nnz, IP, JA, b, c, alphap, betap, transA, typecode, flags);
		dt += rsb_time();
		/* tpo = 1.0 / dt; */
		tpo = dt;
	}
	if(tpop)
		*tpop = tpo;
	return errval;
}
')dnl
dnl

dnl
rsb_err_t rsb__mkl_csr_spsm_bench(const void *VA, const MKL_INT m, const MKL_INT nrhs, const MKL_INT * IP, const MKL_INT *JA, const void * b, void * c, const void *alphap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, const MKL_INT ldb, const MKL_INT ldc, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t dt, tpo;

	if(otnp)
	{
		RSB_MKL_THREADS_TUNING_ODECLS
		RSB_MKL_SET_THREADS_RANGE(lnt,unt,otnp)

		for(nt=lnt;nt<=unt;++nt)
		{
			RSB_MKL_THREADS_TUNING_IDECLS
			if(nt) mkl_set_num_threads(nt);

			do
			{
				errval = rsb__do_mkl_csr_spsm(VA, m, nrhs, IP, JA, b, c, alphap, transA, typecode, flags, ldb, ldc);
				RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,tinf,times);
			}
			while(RSB_REPEAT(ct-it,times,mint,mintimes,maxt,maxtimes));

			dt = bt;
			if(dt < best )
			{
				otn = nt;
				best = RSB_MIN_ABOVE_INF(best,dt,tinf);
				RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp);
			}
			rsb__tattr_sets(ttrp,dtn,nt,dt,otn,times);/* FIXME: if no threads tuning, shall set dtpo = btpo, as well as ttrp.optt=0 */
			if(dtn == nt) RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp+1);
		}
		mkl_set_num_threads(ont);
		ttrp->ttt += rsb_time(); /* ttrp->ttt = tt; */
		tpo = best; /* tpo = 1.0 / best; */
		*otnp = otn;
	}
	else
	{
		dt = -rsb_time();
		errval = rsb__do_mkl_csr_spsm(VA, m, nrhs, IP, JA, b, c, alphap, transA, typecode, flags, ldb, ldc);
		dt += rsb_time();
		/* tpo = 1.0 / dt; */
		tpo = dt;
	}
	if(tpop)
		*tpop = tpo;
	return errval;
}
')dnl
dnl

dnl
rsb_err_t rsb__mkl_coo2csr(const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const void *IVA, const MKL_INT * IIA, const MKL_INT *IJA, const void *OVA, const MKL_INT * OIA, const MKL_INT *OJA, rsb_type_t typecode, const MKL_INT mib)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	MKL_INT info = 0;
	MKL_INT job[6];
	job[0] = 2; // coo2csr (=1, the matrix in the coordinate format is converted to the CSR;=2, the matrix in the coordinate format is converted to the CSR format, and the column indices in CSR representation are sorted in the increasing order within each row.)
	job[1] = mib; // 0 based csr
	job[2] = 0; // 0 based coo
	job[3] = 0; // ignored
	job[4] = nnz; // ignored here
	job[5] = 0; // fill all three arrays

foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	`mkl_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`csrcoo'(job,(MKL_INT*)(&m),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)OVA,(MKL_INT*)OJA,(MKL_INT*)OIA,(MKL_INT*)(&nnz),(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)(IVA),(MKL_INT*)IIA,(MKL_INT*)IJA,&info);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE	;
	if(info==0)
		return RSB_ERR_NO_ERROR;
	else
		return RSB_ERR_INTERNAL_ERROR;
}
')dnl
dnl
dnl
dnl
dnl

dnl
rsb_int_t rsb__print_mkl_version(char *pn, rsb_bool_t los)
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	size_t npc = 0;

	if(pn)
		*pn = RSB_NUL;
	else
		goto err;
{
#if RSB_WANT_MKL
#ifdef mkl_get_version
	MKLVersion mv;

	mkl_get_version(&mv);
	if(los)
		sprintf(pn,"MKL:%d.%d-%d, %s, %s, %s, %s",mv.MajorVersion,mv.MinorVersion,mv.UpdateVersion,mv.ProductStatus,mv.Build,mv.Processor,mv.Platform);
	else
		// sprintf(pn,"-mkl-%d.%d-%d",mv.MajorVersion,mv.MinorVersion,mv.UpdateVersion,mv.ProductStatus);
		sprintf(pn,"-mkl-%d.%d-%d-%s-%s",mv.MajorVersion,mv.MinorVersion,mv.UpdateVersion,mv.ProductStatus,mv.Build);
#else /* mkl_get_version */
	sprintf(pn,"MKL:version unknown.");
#endif /* mkl_get_version */
#else /* RSB_WANT_MKL */
	sprintf(pn,"MKL:not linked.");
#endif /* RSB_WANT_MKL */
	npc = strlen(pn);
}
err:
	return npc;
}
')dnl
dnl

dnl
#if RSB_WANT_MKL_INSPECTOR
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`
#define RSB_MKL_INSPECTOR_TRANSA(TRANS)			\
( (TRANS) == RSB_TRANSPOSITION_N ?  SPARSE_OPERATION_NON_TRANSPOSE : ( (TRANS) == RSB_TRANSPOSITION_T ?  SPARSE_OPERATION_TRANSPOSE : (SPARSE_OPERATION_CONJUGATE_TRANSPOSE /* FIXME: what about conjugation ? */ ) ) )

#define RSB_MKL_INSPECTOR_STCHK(SV,ERRVAL)						\
	{ if( (SV) != SPARSE_STATUS_SUCCESS ) { printf("Error invoking MKL routine !\n"); (ERRVAL) = RSB_ERR_INTERNAL_ERROR; } }

#define RSB_MKL_INSPECTOR_DESCRA(DESCRA,FLAGS,HINTLVL,DONELBL)				\
	rsb_time_t dt, tpo;							\
	sparse_operation_t mkl_transA = RSB_MKL_INSPECTOR_TRANSA(transA);	\
	struct matrix_descr descrA;						\
	sparse_matrix_t       csrA;						\
	rsb_time_t iot = RSB_TIME_ZERO;						\
	sparse_status_t stv = SPARSE_STATUS_NOT_INITIALIZED;				\
	\
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )	\
		`stv = mkl_sparse_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`_create_csr '( &csrA, SPARSE_INDEX_BASE_ZERO, m,  k,  IP, IP+1, JA, VA ); \
')dnl
	if ( stv != SPARSE_STATUS_SUCCESS ) goto DONELBL;			\
	DESCRA.type = SPARSE_MATRIX_TYPE_GENERAL;				\
	if(RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_SYMMETRIC))				\
		DESCRA.type = SPARSE_MATRIX_TYPE_SYMMETRIC,			\
		DESCRA.mode = SPARSE_FILL_MODE_LOWER,				\
		DESCRA.diag = SPARSE_DIAG_NON_UNIT;				\
	if ( HINTLVL > 0 )							\
	{									\
		stv = mkl_sparse_set_memory_hint ( csrA, HINTLVL < 2 ? ( SPARSE_MEMORY_NONE ) : SPARSE_MEMORY_AGGRESSIVE );	\
		RSB_MKL_INSPECTOR_STCHK(stv,errval)					\
		iot =- rsb_time();						\
		stv = mkl_sparse_optimize ( csrA );					\
		iot += rsb_time();						\
		RSB_MKL_INSPECTOR_STCHK(stv,errval)					\
	}									\
										\
	/* RSB_MKL_INSPECTOR_DESCRA */

')dnl
dnl

dnl
rsb_err_t rsb__mkl_inspector_csr_spmv_bench(void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, MKL_INT * IP, MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp, rsb_int_t hintlvl)
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_MKL_INSPECTOR_DESCRA(descrA,flags,hintlvl,done);

	if(otnp)
	{
		RSB_MKL_THREADS_TUNING_ODECLS
		RSB_MKL_SET_THREADS_RANGE(lnt,unt,otnp)

		for(nt=lnt;nt<=unt;++nt)
		{
			RSB_MKL_THREADS_TUNING_IDECLS
			if(nt) mkl_set_num_threads(nt);

			do
			{
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
				if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
				{
					RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype) alpha = alphap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap : RSB_M4_MKL_ONE(mtype), beta = betap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap : RSB_M4_MKL_ONE(mtype);
					stv = `mkl_sparse_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`_mv' ( mkl_transA, alpha, csrA, descrA, b, beta, c );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
					return RSB_ERR_UNSUPPORTED_TYPE;
				RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,tinf,times);
			}
			while(RSB_REPEAT(ct-it,times,mint,mintimes,maxt,maxtimes));

			dt = bt;
			if(dt < best )
			{
				otn = nt;
				best = RSB_MIN_ABOVE_INF(best,dt,tinf);
				RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp);
			}
			rsb__tattr_sets(ttrp,dtn,nt,dt,otn,times);/* FIXME: if no threads tuning, shall set dtpo = btpo, as well as ttrp.optt=0 */
			if(dtn == nt) RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp+1);
		}
		mkl_set_num_threads(ont);
done:
		ttrp->ttt += rsb_time(); /* ttrp->ttt = tt; */
		tpo = best; /* tpo = 1.0 / best; */
		*otnp = otn;
	}
	else
	{
		dt = -rsb_time();
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
		if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
		{
			RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype) alpha = alphap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap : RSB_M4_MKL_ONE(mtype), beta = betap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap : RSB_M4_MKL_ONE(mtype);
			stv = `mkl_sparse_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`_mv' ( mkl_transA, alpha, csrA, descrA, b, beta, c );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
			return RSB_ERR_UNSUPPORTED_TYPE;
		dt += rsb_time();
		/* tpo = 1.0 / dt; */
		tpo = dt;
	}
	if(tpop)
		*tpop = tpo;

	RSB_MKL_INSPECTOR_STCHK(stv,errval);
	// RSBENCH_STDOUT("# MKL Inspector SPMV time: %2.3le   optimization time: %2.3le!\n",tpo,iot);
	stv = mkl_sparse_destroy ( csrA );
	RSB_MKL_INSPECTOR_STCHK(stv,errval);
	return errval;
}
')dnl
dnl

rsb_err_t rsb__mkl_inspector_csr_spmm_bench(void *VA, const MKL_INT m, const MKL_INT k, const sparse_layout_t layout, const MKL_INT nrhs, const MKL_INT nnz, MKL_INT * IP, MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp, rsb_int_t hintlvl)
ifdef(`ONLY_WANT_HEADERS',`;',`
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_MKL_INSPECTOR_DESCRA(descrA,flags,hintlvl,done);
{
	const MKL_INT ldx = ldb;
	const MKL_INT ldy = ldc;

	if(otnp)
	{
		RSB_MKL_THREADS_TUNING_ODECLS
		RSB_MKL_SET_THREADS_RANGE(lnt,unt,otnp)

		for(nt=lnt;nt<=unt;++nt)
		{
			RSB_MKL_THREADS_TUNING_IDECLS
			if(nt) mkl_set_num_threads(nt);

			do
			{
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
				if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
				{
					RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype) alpha = alphap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap : RSB_M4_MKL_ONE(mtype), beta = betap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap : RSB_M4_MKL_ONE(mtype);
					stv = `mkl_sparse_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`_mm' ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
					return RSB_ERR_UNSUPPORTED_TYPE;
				RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,tinf,times);
			}
			while(RSB_REPEAT(ct-it,times,mint,mintimes,maxt,maxtimes));

			dt = bt;
			if(dt < best )
			{
				otn = nt;
				best = RSB_MIN_ABOVE_INF(best,dt,tinf);
				RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp);
			}
			rsb__tattr_sets(ttrp,dtn,nt,dt,otn,times);/* FIXME: if no threads tuning, shall set dtpo = btpo, as well as ttrp.optt=0 */
			if(dtn == nt) RSB_STAT_TAKE(it,otn,ct,dt,tt,bt,wt,ss,times,tstp+1);
		}
		mkl_set_num_threads(ont);
done:
		ttrp->ttt += rsb_time(); /* ttrp->ttt = tt; */
		tpo = best; /* tpo = 1.0 / best; */
		*otnp = otn;
	}
	else
	{
		dt = -rsb_time();
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
		if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
		{
			RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype) alpha = alphap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)alphap : RSB_M4_MKL_ONE(mtype), beta = betap ? *(RSB_M4_RSB_TYPE_TO_MKL_TYPE(mtype)*)betap : RSB_M4_MKL_ONE(mtype);
			stv = `mkl_sparse_'RSB_M4_SPBLAS_TYPE_CHARCODE(mtype)`_mm' ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
			return RSB_ERR_UNSUPPORTED_TYPE;
		dt += rsb_time();
		/* tpo = 1.0 / dt; */
		tpo = dt;
	}
	if(tpop)
		*tpop = tpo;

	RSB_MKL_INSPECTOR_STCHK(stv,errval);
	// RSBENCH_STDOUT("# MKL Inspector SPMV time: %2.3le   optimization time: %2.3le!\n",tpo,iot);
	stv = mkl_sparse_destroy ( csrA );
	RSB_MKL_INSPECTOR_STCHK(stv,errval);
	return errval;
}
}
')dnl
dnl

dnl
#endif /* RSB_WANT_MKL_INSPECTOR */

dnl
#endif /* RSB_WANT_MKL */
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#endif  /* RSB_RSB_MKL_H_INCLUDED */
')dnl
dnl
/* @endcond */

