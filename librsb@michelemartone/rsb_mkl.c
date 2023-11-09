/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains MKL interfacing functions.
 * */


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
#include "rsb_mkl.h"
#if RSB_WANT_MKL
RSB_INTERNALS_COMMON_HEAD_DECLS


#include "rsb_tune.h" /* rsb_tattr_t */



rsb_err_t rsb__mkl_gemv(rsb_type_t typecode, const void * Mp, const void*Bp, void*Xp, rsb_nnz_idx_t mdim, rsb_coo_idx_t vdim, rsb_coo_idx_t*udimp)
{
	/* FIXME: TODO: incX != 1 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const MKL_INT dim=(rsb_coo_idx_t)sqrt((double)mdim);
	const MKL_INT incX=1;
	char transA_mkl=110;
	if(!Mp || !Xp || !Bp)
		goto err;
	if(dim<1 || dim>vdim)
		goto err;
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
		const float alpha=((float)(1.0)), beta=((float)(1.0));
		sgemv(&transA_mkl,&dim,&dim, (float*)(&alpha),(const float*)Mp,&dim,(const float*)Bp,&incX,(float*)&beta,(float*)Xp,&incX);
	}
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	else 
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
		const double alpha=((double)(1.0)), beta=((double)(1.0));
		dgemv(&transA_mkl,&dim,&dim, (double*)(&alpha),(const double*)Mp,&dim,(const double*)Bp,&incX,(double*)&beta,(double*)Xp,&incX);
	}
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	else 
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
		const float complex alpha=((float complex)(1.0)), beta=((float complex)(1.0));
		cgemv(&transA_mkl,&dim,&dim, (MKL_Complex8*)(&alpha),(const MKL_Complex8*)Mp,&dim,(const MKL_Complex8*)Bp,&incX,(MKL_Complex8*)&beta,(MKL_Complex8*)Xp,&incX);
	}
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	else 
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
		const double complex alpha=((double complex)(1.0)), beta=((double complex)(1.0));
		zgemv(&transA_mkl,&dim,&dim, (MKL_Complex16*)(&alpha),(const MKL_Complex16*)Mp,&dim,(const MKL_Complex16*)Bp,&incX,(MKL_Complex16*)&beta,(MKL_Complex16*)Xp,&incX);
	}
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	else 
		errval=RSB_ERR_BADARGS;

	if(udimp)
		*udimp=dim;
err:
	return errval;
}

static char rsb_rsb_to_mkl_trans(rsb_trans_t transA_mkl)
{
	/**
	 * \ingroup gr_internals
	 */
	switch(transA_mkl)
	{
		case(RSB_TRANSPOSITION_N):
		return 'n';
		break;
		case(RSB_TRANSPOSITION_T):
		return 't';
		break;
		case(RSB_TRANSPOSITION_C):
		return 'c';
		break;
		default:
		return 'n';	// FIXME
	}
}

static char rsb_rsb_to_mkl_sym(rsb_flags_t flags)
{
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC))
		return 's';
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_TRIANGULAR))
		return 't';
	else
		return 'g';
}

static char rsb_rsb_to_mkl_upl(rsb_flags_t flags)
{
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER))
		return 'l';
	else
		return 'u';
}

rsb_err_t rsb__mkl_coo_spmv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = 'n'; // not unit diagonal
	matdescra[3] = 'c'; // zero based indexing

#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scoomv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(float*)alphap,matdescra,(float*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(float*)x,(float*)betap,(float*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcoomv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(double*)alphap,matdescra,(double*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(double*)x,(double*)betap,(double*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccoomv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(MKL_Complex8*)alphap,matdescra,(MKL_Complex8*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(MKL_Complex8*)x,(MKL_Complex8*)betap,(MKL_Complex8*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcoomv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(MKL_Complex16*)alphap,matdescra,(MKL_Complex16*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(MKL_Complex16*)x,(MKL_Complex16*)betap,(MKL_Complex16*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__mkl_coo_spmm(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nrhs, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = 'n'; // not unit diagonal
	matdescra[3] = 'c'; // zero based indexing

#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scoomm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(MKL_INT*)(&k),(float*)alphap,matdescra,(float*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(float*)b,(MKL_INT*)(&ldb),(float*)betap,(float*)c,(MKL_INT*)(&ldc));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcoomm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(MKL_INT*)(&k),(double*)alphap,matdescra,(double*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(double*)b,(MKL_INT*)(&ldb),(double*)betap,(double*)c,(MKL_INT*)(&ldc));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccoomm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(MKL_INT*)(&k),(MKL_Complex8*)alphap,matdescra,(MKL_Complex8*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(MKL_Complex8*)b,(MKL_INT*)(&ldb),(MKL_Complex8*)betap,(MKL_Complex8*)c,(MKL_INT*)(&ldc));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcoomm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(MKL_INT*)(&k),(MKL_Complex16*)alphap,matdescra,(MKL_Complex16*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(MKL_Complex16*)b,(MKL_INT*)(&ldb),(MKL_Complex16*)betap,(MKL_Complex16*)c,(MKL_INT*)(&ldc));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__mkl_coo_spsv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IA, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = 'n'; // not unit diagonal
	matdescra[3] = 'c'; // zero based indexing
	/* 20101118	MKL 9.1 reference manual declares also k among the parameters */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scoosv(&transA_mkl,(MKL_INT*)(&m),(float*)alphap,matdescra,(float*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(float*)x,(float*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcoosv(&transA_mkl,(MKL_INT*)(&m),(double*)alphap,matdescra,(double*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(double*)x,(double*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccoosv(&transA_mkl,(MKL_INT*)(&m),(MKL_Complex8*)alphap,matdescra,(MKL_Complex8*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(MKL_Complex8*)x,(MKL_Complex8*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcoosv(&transA_mkl,(MKL_INT*)(&m),(MKL_Complex16*)alphap,matdescra,(MKL_Complex16*)VA,(MKL_INT*)IA,(MKL_INT*)JA,(MKL_INT*)(&nnz),(MKL_Complex16*)x,(MKL_Complex16*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__do_mkl_csr_spmv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags){
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = 'n'; // not unit diagonal
	matdescra[3] = 'c'; // zero based indexing

#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scsrmv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(float*)alphap,matdescra,(float*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(float*)x,(float*)betap,(float*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcsrmv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(double*)alphap,matdescra,(double*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(double*)x,(double*)betap,(double*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccsrmv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(MKL_Complex8*)alphap,matdescra,(MKL_Complex8*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex8*)x,(MKL_Complex8*)betap,(MKL_Complex8*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcsrmv(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&k),(MKL_Complex16*)alphap,matdescra,(MKL_Complex16*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex16*)x,(MKL_Complex16*)betap,(MKL_Complex16*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

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

rsb_err_t rsb__mkl_csr_spmv_bench(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)
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

rsb_err_t rsb__do_mkl_csr_spmm(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags){
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	MKL_INT ldb_ = n, ldc_ = n; /* for zero based indexing */

	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = 'n'; // not unit diagonal
	matdescra[3] = 'c'; // zero based indexing

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		ldb_ = k, ldc_ = m, /* for one based indexing */
		matdescra[3] = 'f'; // one based indexing

	#if 1
	/* n = nrhs */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scsrmm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&n),(MKL_INT*)(&k),(float*)alphap,matdescra,(float*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(float*)b,(MKL_INT*)(&ldb_),(float*)betap,(float*)c,(MKL_INT*)(&ldc_));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcsrmm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&n),(MKL_INT*)(&k),(double*)alphap,matdescra,(double*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(double*)b,(MKL_INT*)(&ldb_),(double*)betap,(double*)c,(MKL_INT*)(&ldc_));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccsrmm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&n),(MKL_INT*)(&k),(MKL_Complex8*)alphap,matdescra,(MKL_Complex8*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex8*)b,(MKL_INT*)(&ldb_),(MKL_Complex8*)betap,(MKL_Complex8*)c,(MKL_INT*)(&ldc_));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcsrmm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&n),(MKL_INT*)(&k),(MKL_Complex16*)alphap,matdescra,(MKL_Complex16*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex16*)b,(MKL_INT*)(&ldb_),(MKL_Complex16*)betap,(MKL_Complex16*)c,(MKL_INT*)(&ldc_));
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
#endif /* 1 */
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__mkl_csr_spmm_bench(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)
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

rsb_err_t rsb__do_mkl_csr_spsv(const void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags)
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = 'n'; // not unit diagonal
	matdescra[3] = 'c'; // zero based indexing

#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scsrsv(&transA_mkl,(MKL_INT*)(&m),(float*)alphap,matdescra,(float*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(float*)x,(float*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcsrsv(&transA_mkl,(MKL_INT*)(&m),(double*)alphap,matdescra,(double*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(double*)x,(double*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccsrsv(&transA_mkl,(MKL_INT*)(&m),(MKL_Complex8*)alphap,matdescra,(MKL_Complex8*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex8*)x,(MKL_Complex8*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcsrsv(&transA_mkl,(MKL_INT*)(&m),(MKL_Complex16*)alphap,matdescra,(MKL_Complex16*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex16*)x,(MKL_Complex16*)y);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__do_mkl_csr_spsm(const void *VA, const MKL_INT m, const MKL_INT nrhs, const MKL_INT * IP, const MKL_INT *JA, const void * b, void * c, const void *alphap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, const MKL_INT ldb, const MKL_INT ldc)
{
	char transA_mkl = rsb_rsb_to_mkl_trans(transA);
	char matdescra[]={RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL,RSB_NUL};
	matdescra[0] = rsb_rsb_to_mkl_sym(flags); // general ?
	matdescra[1] = rsb_rsb_to_mkl_upl(flags); // up or lo ?
	matdescra[2] = 'n'; // not unit diagonal
	matdescra[3] = 'c'; // zero based indexing

#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scsrsm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(float*)alphap,matdescra,(float*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(float*)b,(MKL_INT*)&ldb,(float*)c,(MKL_INT*)&ldc);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcsrsm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(double*)alphap,matdescra,(double*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(double*)b,(MKL_INT*)&ldb,(double*)c,(MKL_INT*)&ldc);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccsrsm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(MKL_Complex8*)alphap,matdescra,(MKL_Complex8*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex8*)b,(MKL_INT*)&ldb,(MKL_Complex8*)c,(MKL_INT*)&ldc);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcsrsm(&transA_mkl,(MKL_INT*)(&m),(MKL_INT*)(&nrhs),(MKL_Complex16*)alphap,matdescra,(MKL_Complex16*)VA,(MKL_INT*)JA,(MKL_INT*)IP,(MKL_INT*)(IP+1),(MKL_Complex16*)b,(MKL_INT*)&ldb,(MKL_Complex16*)c,(MKL_INT*)&ldc);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__mkl_csr_spsv_bench(const void *VA, const MKL_INT m, const MKL_INT k/*, const MKL_INT n*/, const MKL_INT nnz, const MKL_INT * IP, const MKL_INT *JA, const void * b, /*const MKL_INT ldb,*/ void * c,/* const MKL_INT ldc,*/ const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)
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

rsb_err_t rsb__mkl_csr_spsm_bench(const void *VA, const MKL_INT m, const MKL_INT nrhs, const MKL_INT * IP, const MKL_INT *JA, const void * b, void * c, const void *alphap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, const MKL_INT ldb, const MKL_INT ldc, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp)
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

rsb_err_t rsb__mkl_coo2csr(const MKL_INT m, const MKL_INT k, const MKL_INT nnz, const void *IVA, const MKL_INT * IIA, const MKL_INT *IJA, const void *OVA, const MKL_INT * OIA, const MKL_INT *OJA, rsb_type_t typecode, const MKL_INT mib)
{
	MKL_INT info = 0;
	MKL_INT job[6];
	job[0] = 2; // coo2csr (=1, the matrix in the coordinate format is converted to the CSR;=2, the matrix in the coordinate format is converted to the CSR format, and the column indices in CSR representation are sorted in the increasing order within each row.)
	job[1] = mib; // 0 based csr
	job[2] = 0; // 0 based coo
	job[3] = 0; // ignored
	job[4] = nnz; // ignored here
	job[5] = 0; // fill all three arrays

#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	mkl_scsrcoo(job,(MKL_INT*)(&m),(float*)OVA,(MKL_INT*)OJA,(MKL_INT*)OIA,(MKL_INT*)(&nnz),(float*)(IVA),(MKL_INT*)IIA,(MKL_INT*)IJA,&info);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	mkl_dcsrcoo(job,(MKL_INT*)(&m),(double*)OVA,(MKL_INT*)OJA,(MKL_INT*)OIA,(MKL_INT*)(&nnz),(double*)(IVA),(MKL_INT*)IIA,(MKL_INT*)IJA,&info);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	mkl_ccsrcoo(job,(MKL_INT*)(&m),(MKL_Complex8*)OVA,(MKL_INT*)OJA,(MKL_INT*)OIA,(MKL_INT*)(&nnz),(MKL_Complex8*)(IVA),(MKL_INT*)IIA,(MKL_INT*)IJA,&info);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	mkl_zcsrcoo(job,(MKL_INT*)(&m),(MKL_Complex16*)OVA,(MKL_INT*)OJA,(MKL_INT*)OIA,(MKL_INT*)(&nnz),(MKL_Complex16*)(IVA),(MKL_INT*)IIA,(MKL_INT*)IJA,&info);
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
		return RSB_ERR_UNSUPPORTED_TYPE	;
	if(info==0)
		return RSB_ERR_NO_ERROR;
	else
		return RSB_ERR_INTERNAL_ERROR;
}

rsb_int_t rsb__print_mkl_version(char *pn, rsb_bool_t los)

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

#if RSB_WANT_MKL_INSPECTOR

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
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )	\
		stv = mkl_sparse_s_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO, m,  k,  IP, IP+1, JA, VA ); \
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )	\
		stv = mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO, m,  k,  IP, IP+1, JA, VA ); \
	if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )	\
		stv = mkl_sparse_c_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO, m,  k,  IP, IP+1, JA, VA ); \
	if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )	\
		stv = mkl_sparse_z_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO, m,  k,  IP, IP+1, JA, VA ); \
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


rsb_err_t rsb__mkl_inspector_csr_spmv_bench(void *VA, const MKL_INT m, const MKL_INT k, const MKL_INT n, const MKL_INT nnz, MKL_INT * IP, MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp, rsb_int_t hintlvl)

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
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
				if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
				{
					float alpha = alphap ? *(float*)alphap : (float)(1.0), beta = betap ? *(float*)betap : (float)(1.0);
					stv = mkl_sparse_s_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
				if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
				{
					double alpha = alphap ? *(double*)alphap : (double)(1.0), beta = betap ? *(double*)betap : (double)(1.0);
					stv = mkl_sparse_d_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
				if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
				{
					MKL_Complex8 alpha = alphap ? *(MKL_Complex8*)alphap : (MKL_Complex8) {1.0,0.0} , beta = betap ? *(MKL_Complex8*)betap : (MKL_Complex8) {1.0,0.0} ;
					stv = mkl_sparse_c_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
				if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
				{
					MKL_Complex16 alpha = alphap ? *(MKL_Complex16*)alphap : (MKL_Complex16) {1.0,0.0} , beta = betap ? *(MKL_Complex16*)betap : (MKL_Complex16) {1.0,0.0} ;
					stv = mkl_sparse_z_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
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
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
		if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
		{
			float alpha = alphap ? *(float*)alphap : (float)(1.0), beta = betap ? *(float*)betap : (float)(1.0);
			stv = mkl_sparse_s_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
		if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
		{
			double alpha = alphap ? *(double*)alphap : (double)(1.0), beta = betap ? *(double*)betap : (double)(1.0);
			stv = mkl_sparse_d_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
		if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
		{
			MKL_Complex8 alpha = alphap ? *(MKL_Complex8*)alphap : (MKL_Complex8) {1.0,0.0} , beta = betap ? *(MKL_Complex8*)betap : (MKL_Complex8) {1.0,0.0} ;
			stv = mkl_sparse_c_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
		if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
		{
			MKL_Complex16 alpha = alphap ? *(MKL_Complex16*)alphap : (MKL_Complex16) {1.0,0.0} , beta = betap ? *(MKL_Complex16*)betap : (MKL_Complex16) {1.0,0.0} ;
			stv = mkl_sparse_z_mv ( mkl_transA, alpha, csrA, descrA, b, beta, c );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
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

rsb_err_t rsb__mkl_inspector_csr_spmm_bench(void *VA, const MKL_INT m, const MKL_INT k, const sparse_layout_t layout, const MKL_INT nrhs, const MKL_INT nnz, MKL_INT * IP, MKL_INT *JA, const void * b, const MKL_INT ldb, void * c, const MKL_INT ldc, const void *alphap, const void * betap, rsb_trans_t transA, rsb_type_t typecode, rsb_flags_t flags, rsb_thread_t *otnp, rsb_time_t *tpop, struct rsb_tattr_t* ttrp, struct rsb_ts_t*tstp, rsb_int_t hintlvl)

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
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
				if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
				{
					float alpha = alphap ? *(float*)alphap : (float)(1.0), beta = betap ? *(float*)betap : (float)(1.0);
					stv = mkl_sparse_s_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
				if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
				{
					double alpha = alphap ? *(double*)alphap : (double)(1.0), beta = betap ? *(double*)betap : (double)(1.0);
					stv = mkl_sparse_d_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
				if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
				{
					MKL_Complex8 alpha = alphap ? *(MKL_Complex8*)alphap : (MKL_Complex8) {1.0,0.0} , beta = betap ? *(MKL_Complex8*)betap : (MKL_Complex8) {1.0,0.0} ;
					stv = mkl_sparse_c_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
				if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
				{
					MKL_Complex16 alpha = alphap ? *(MKL_Complex16*)alphap : (MKL_Complex16) {1.0,0.0} , beta = betap ? *(MKL_Complex16*)betap : (MKL_Complex16) {1.0,0.0} ;
					stv = mkl_sparse_z_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
				}
				else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
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
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
		if( typecode == RSB_NUMERICAL_TYPE_FLOAT  )
		{
			float alpha = alphap ? *(float*)alphap : (float)(1.0), beta = betap ? *(float*)betap : (float)(1.0);
			stv = mkl_sparse_s_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
		if( typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
		{
			double alpha = alphap ? *(double*)alphap : (double)(1.0), beta = betap ? *(double*)betap : (double)(1.0);
			stv = mkl_sparse_d_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
		if( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
		{
			MKL_Complex8 alpha = alphap ? *(MKL_Complex8*)alphap : (MKL_Complex8) {1.0,0.0} , beta = betap ? *(MKL_Complex8*)betap : (MKL_Complex8) {1.0,0.0} ;
			stv = mkl_sparse_c_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
		if( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
		{
			MKL_Complex16 alpha = alphap ? *(MKL_Complex16*)alphap : (MKL_Complex16) {1.0,0.0} , beta = betap ? *(MKL_Complex16*)betap : (MKL_Complex16) {1.0,0.0} ;
			stv = mkl_sparse_z_mm ( mkl_transA, alpha, csrA, descrA, layout, b, nrhs, ldx /*m*/, beta, c, ldy /*m*/ );
		}
		else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
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

#endif /* RSB_WANT_MKL_INSPECTOR */

#endif /* RSB_WANT_MKL */
/* @endcond */

