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
/**
 * @file
 * @brief Perfomance reporting code. This source uses rsb__getenv(), therefore should not be linked to librsb but to rsbench directly.
 * @author Michele Martone
 * */

#include <strings.h>		/* bzero */
#include "rsb_internals.h"
#include "rsb-config.h"
#include <stdint.h> /* int64_t / uint64_t */

#define rsb__strcpy(A,B) strcpy((rsb_char_t*)A,B)
#define rsb__strlen(A) strlen((rsb_char_t*)A)
#define RSB_RMEMCPY(DEST,SRC,N) RSB_MEMCPY((void*RSB_RESTRICT)(DEST),(const void*RSB_RESTRICT)(SRC),(N))
#define RSB__PR_FREE(P) {rsb__pr_free(P);(P)=NULL;}
#define RSB_XFLOPS(TIME,NRHS,CANONICAL_MATRIX_OP_FLOPS)  ( (TIME)?(((double)NRHS)*(CANONICAL_MATRIX_OP_FLOPS))/(TIME):RSB_TIME_ZERO )

/* BEGIN values for RSB_PR_SR: */
#define RSB_PRD_STYLE_TBL 0
#define RSB_PRD_STYLE_CMP 1
/* base plots */
#define RSB_PRD_STYLE_PLT_BASE 2 /* new, experimental */
#define RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB 2 /* new, experimental */
#define RSB_PRD_STYLE_PLT_SUBM_BS 3 /* new, experimental */
/* polar equivalent plots */
#define RSB_PRD_STYLE_PLT_BASE_POLAR 102 /* new, experimental; TODO: unused */
#define RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB_POLAR 103 /* new, experimental */
#define RSB_PRD_STYLE_PLT_MKL_SPEEDUP_POLAR 104 /* new, experimental */
#define RSB_PRD_STYLE_PLT_SUBM_BS_POLAR 105 /* new, experimental; TODO: unused */
/*  END  values for RSB_PR_SR. */
/* */
#define RSB_PRD_CMP_MDUMP -1
#define RSB_PRD_CMP_DFLT 0
#define RSB_PRD_CMP_DIV 1
#define RSB_PRD_CMP_DIFF 2
#define RSB_PRD_CMP_APPEND 3
#define RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH 1
#define RSB_PRD_VERSION_BAD -1

#define RSB_PRD_VERSION 1
#define RSB_PRD_WANT_MBW ( RSB_PRD_VERSION > 0 )
#define RSB_PRD_WANT_ENVIRONMENT ( RSB_PRD_VERSION > 0 ) /* FIXME: _WANT_ */
#define RSB_PRD_WANT_TIMESTAMP ( RSB_PRD_VERSION > 0 )

#define RSB_ON_IF_LEM(X,Y,ONIFLESS,ONIFEQUAL,ONIFMORE)	\
	( (X)==(Y) ? (ONIFEQUAL) : ((X)<(Y)? (ONIFLESS) : (ONIFMORE) ))
#define RSB_MAX_LABEL_LENGTH RSB_MAX_FILENAME_LENGTH

/* rsb sampled performance sample structure (internal) */
/* to keep I/O portable, don't use pointers or unportable variables within it */
/* at_-prefixed variables are those relative to autotuning samples */
struct rsb_rsps_t
{
	rsb_perf_t op_time;
	rsb_perf_t mkl_csr_op_time;
	rsb_perf_t at_op_time;
	rsb_time_t at_t;
	rsb_perf_t at_mkl_csr_op_time;
	rsb_trans_t transA;
	double cmflops; /* canonical mflops considering nrhs==1 */
	rsb_flags_t flagsA;
	rsb_submatrix_idx_t nsubm, at_nsubm;
	int64_t isa, at_isa; /* indexing space allocated */
	rsb_int_t at_cn, at_mkl_csr_cn;
	rsb_int_t uc; /* updates count */
	rsb_coo_idx_t nrA,ncA;
	rsb_nnz_idx_t nnzA;
	rsb_int_t at_eps; /* effective steps */
        struct rsb_ts_t otpos, btpos; /* dumpable with RSB_STAT_DUMP_TS */
        struct rsb_ts_t otpms, btpms;
};

/* rsb sampled performance record structure  (internal)*/
struct rsb_rspr_t
{
	rsb_int_t  filenamen,   cn,   incXn,   incYn,   nrhsn,  ntypecodes,   tn, csf /* count so far */;
        rsb_int_t filenamebl, cabl, incXabl, incYabl, nrhsabl, typecodesbl, tabl; /* ... byte length */
        /* the following shall not be saved */
        rsb_bool_t ror; /* representing only ratios */
	struct rsb_rsps_t * psa; /* performance samples array */
        struct rsb_rspra_t * rsprap; /*  */
#if RSB_PRD_WANT_ENVIRONMENT
	rsb_int_t nenvv; // number of environment variables
	rsb_int_t enoib; // environment occupation in bytes
	char * envvp;
#endif /* RSB_PRD_WANT_ENVIRONMENT */
#if RSB_PRD_WANT_MBW
	struct rsb_mbw_et_t mbet;
#endif /* RSB_PRD_WANT_MBW */
#if RSB_PRD_WANT_TIMESTAMP
	rsb_time_t tbeg,tend;
#endif /* RSB_PRD_WANT_TIMESTAMP */
};

#define RSB_PRL_TCS "pr: "
#define RSB_PRL_LCC_IE rsb__getenv("RSB_PR_WLTC") ? '%' : ( rsb__getenv_char("RSB_PR_PRL_LCC", '#') ) /*  line comment char */
#define RSB_PRL_TCS_IE rsb__getenv("RSB_PR_WLTC") ? " " : ( rsb__getenv_str("RSB_PR_PRL_TCS",RSB_PRL_TCS) ) /* table comment string */
#define RSB_PRL_ENDLSTR_IE rsb__getenv("RSB_PR_WLTC") ? "\\\\" : ( rsb__getenv_str("RSB_PR_ENDLSTR","") )
#define RSB_PRL_FSEPSTR_IE rsb__getenv("RSB_PR_WLTC") ? " & " : ( rsb__getenv_str("RSB_PR_FSEPSTR"," ") )
#define RSB_PR_NOC(RSPRP) ((RSPRP)->filenamen * (RSPRP)->cn * (RSPRP)->incXn * (RSPRP)->incYn * (RSPRP)->nrhsn * (RSPRP)->ntypecodes * (RSPRP)->tn )
#define RSB_PRC RSB_STDOUT
#define RSB_PRL RSB_PRC("%c%s",rsb_prl_lcc,RSB_PRL_TCS),RSB_PRC
#define RSB_PRT RSB_PRC("%s",rsb_prl_tcs),RSB_PRC
#define RSB_PRL_SEP RSB_STDOUT("%cpr: ======== ",rsb_prl_lcc),RSB_STDOUT
#define RSB_PRWL RSB_PRC("#pr: Warning:"),RSB_PRC

#if RSB_WANT_ARMPL
#define RSB_MKL_S "APL" /* TODO: 'mkl' substrings shall be translated similarly */
#else
#define RSB_MKL_S "MKL"
#endif

static rsb_err_t rsb__rspr_all_env(struct rsb_rspr_t * rsprp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_PRD_WANT_ENVIRONMENT
	extern char **environ;

	if( environ && environ[0] )
	{
		rsb_int_t envvi = 0;

		for ( envvi = 0; environ[envvi] ; ++envvi)
	       		rsprp->enoib += rsb__strlen(environ[envvi]);
	       	rsprp->nenvv = envvi;
		rsprp->envvp = rsb__calloc( rsprp->enoib + rsprp->nenvv );
		if( rsprp->envvp && rsprp->nenvv )
		{
			size_t co = 0, tco = 0; // char offset, total char offset
			for ( envvi = 0; envvi < rsprp->nenvv ; ++envvi)
				co = strlen(environ[envvi]),
				rsb__strcpy(rsprp->envvp+tco+envvi,environ[envvi]),
				rsprp->envvp[tco+co+envvi] = RSB_NUL,
				tco += co;

   		     	RSB_ASSERT( tco == rsprp->enoib );
		}
	}
#endif /* RSB_PRD_WANT_ENVIRONMENT */
	return errval;
}

static rsb_err_t rsb__pr_alloc(struct rsb_rspr_t ** rsprpp, const struct rsb_rspr_t * rsprcp, rsb_int_t filenamen, rsb_int_t cn, rsb_int_t incXn, rsb_int_t incYn, rsb_int_t nrhsn, rsb_int_t ntypecodes, rsb_int_t tn)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const char rsb_prl_lcc = RSB_PRL_LCC_IE;
	struct rsb_rspr_t * rsprp = NULL;
	rsb_int_t noc = 0; /* number of combinations */
	size_t ab = 0; /* allocated bytes */

	rsprp = rsb__calloc(sizeof(struct rsb_rspr_t));
	if( ! rsprp )
       	{
	       	errval = RSB_ERR_ENOMEM;
	       	goto err;
       	}

        rsprp->ror = RSB_BOOL_FALSE;

        if( rsprcp )
	{
                *rsprp = *rsprcp;
#if RSB_PRD_WANT_ENVIRONMENT
		rsprp->envvp = NULL;
#endif /* RSB_PRD_WANT_ENVIRONMENT */
		rsprp->rsprap = NULL;
#if RSB_PRD_WANT_MBW
		rsprp->mbet.et = NULL;
#endif /* RSB_PRD_WANT_MBW */
	}

	rsprp->filenamen = filenamen;
	rsprp->cn = cn;
	rsprp->incXn = incXn;
	rsprp->incYn = incYn;
	rsprp->nrhsn = nrhsn;
	rsprp->ntypecodes = ntypecodes;
	rsprp->tn = tn;

	noc = RSB_PR_NOC(rsprp);
	ab = sizeof(struct rsb_rsps_t)*noc;
	rsprp->psa = rsb__calloc(ab);
	if( ! rsprp->psa )
	{
	       	errval = RSB_ERR_ENOMEM;
	       	goto err;
       	}

	RSB_ASSIGN_IF(rsprpp,rsprp)

	RSB_PRL("allocated a performance record for %d samples (%zd bytes).\n",noc,ab);
err:
        return errval;
}

struct rsb_rspra_t /* ... record arrays */
{
        rsb_char_t** RSB_RESTRICT filenamea; rsb_int_t*ca; const rsb_int_t*incXa; const rsb_int_t*incYa; const rsb_int_t*nrhsa; const rsb_type_t*typecodes; const rsb_int_t*ta;
};

#define RSB_RPR_FILE_HDR_V0 "%RPR-0..""        ""        ""        "
#define RSB_RPR_FILE_HDR "%RPR-1..""        ""        ""        "
#define RSB_RPR_FILE_HDL 32
#define RSB_RW(ROW,PTR,SIZE,NMEMB,STREAM)                               \
        {                                                               \
                sh = (SIZE) * (NMEMB);                                  \
                if(ROW)                                                 \
                        hd = fread((PTR),(SIZE),(NMEMB),(STREAM));           \
                else                                                    \
                        hd = fwrite ((PTR),(SIZE),(NMEMB),(STREAM));         \
                hd *= (SIZE);                                           \
                if( hd != sh )						\
		{							\
	       		RSB_PERR_GOTO(err,RSB_ERRM_ES);     		\
		}							\
        }

rsb_err_t rsb__rspr_rw(void * p, size_t sop, FILE * stream, rsb_bool_t row)
{
	const rsb_err_t errval = RSB_ERR_NO_ERROR;
        size_t hd = 0, sh = 0; /* have done, should have */

	RSB_RW(row,p,sop,1,stream);
err:
	/* FIXME */
	return errval;
}

static rsb_err_t rsb__rspr_rw_env(struct rsb_rspr_t * rsprp, FILE * stream, rsb_bool_t row)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        size_t hd = 0, sh = 0; /* have done, should have */

#if RSB_PRD_WANT_ENVIRONMENT
	RSB_RW(row,&rsprp->nenvv,sizeof(rsprp->nenvv),1,stream);
	RSB_RW(row,&rsprp->enoib,sizeof(rsprp->enoib),1,stream);

	RSB_ASSERT(rsprp->enoib > 0);
	RSB_ASSERT(rsprp->nenvv > 0);

	if( rsprp->enoib <= 0 )
		goto err;

	if( row ==  RSB_PR_RD )
	{
		rsprp->envvp = rsb__calloc( rsprp->nenvv + rsprp->enoib );
		if( rsprp->envvp && rsprp->nenvv )
		{
			if( 1 != fread( rsprp->envvp, rsprp->nenvv + rsprp->enoib, 1, stream) )
		        {
         		       errval = RSB_ERR_ENOMEM;
 		               RSB_PERR_GOTO(err,RSB_ERRM_ES);
		        }
        	}
	}
	else
	if( row ==  RSB_PR_WR )
	{
		if( 1 != fwrite(((rsb_byte_t*)(rsprp->envvp)), rsprp->nenvv + rsprp->enoib, 1, stream) )
       		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
	}
	goto ret;
#endif /* RSB_PRD_WANT_ENVIRONMENT */
err:
        RSB_ERROR("%s only %zd bytes instead of %zd !\n",row?"read":"wrote",hd,sh);
ret:
	return errval;
}

static rsb_err_t rsb__rsprp_rw(struct rsb_rspr_t * rsprp, FILE * stream, rsb_bool_t row)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        size_t hd = 0, sh = 0; /* have done, should have */

	RSB_RW(row,&rsprp->filenamen,sizeof(rsprp->filenamen),1,stream);
	RSB_RW(row,&rsprp->cn,sizeof(rsprp->cn),1,stream);
	RSB_RW(row,&rsprp->incXn,sizeof(rsprp->incXn),1,stream);
	RSB_RW(row,&rsprp->incYn,sizeof(rsprp->incYn),1,stream);
	RSB_RW(row,&rsprp->nrhsn,sizeof(rsprp->nrhsn),1,stream);
	RSB_RW(row,&rsprp->ntypecodes,sizeof(rsprp->ntypecodes),1,stream);
	RSB_RW(row,&rsprp->tn,sizeof(rsprp->tn),1,stream);
	RSB_RW(row,&rsprp->csf,sizeof(rsprp->csf),1,stream);
	RSB_RW(row,&rsprp->filenamebl,sizeof(rsprp->filenamebl),1,stream);
	RSB_RW(row,&rsprp->cabl,sizeof(rsprp->cabl),1,stream);
	RSB_RW(row,&rsprp->incXabl,sizeof(rsprp->incXabl),1,stream);
	RSB_RW(row,&rsprp->incYabl,sizeof(rsprp->incYabl),1,stream);
	RSB_RW(row,&rsprp->nrhsabl,sizeof(rsprp->nrhsabl),1,stream);
	RSB_RW(row,&rsprp->typecodesbl,sizeof(rsprp->typecodesbl),1,stream);
	RSB_RW(row,&rsprp->tabl,sizeof(rsprp->tabl),1,stream);
        goto ret;
err:
        errval = RSB_ERR_INTERNAL_ERROR;
        RSB_ERROR("%s only %zd bytes instead of %zd !\n",row?"read":"wrote",hd,sh);
ret:
        return errval;
}

static rsb_err_t rsb__ts_rw(struct rsb_ts_t * tsp, FILE * stream, rsb_bool_t row)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        size_t hd = 0, sh = 0; /* have done, should have */

	RSB_RW(row,&tsp->avg,sizeof(tsp->avg),1,stream);
	RSB_RW(row,&tsp->min,sizeof(tsp->min),1,stream);
	RSB_RW(row,&tsp->max,sizeof(tsp->max),1,stream);
	RSB_RW(row,&tsp->sd ,sizeof(tsp->sd ),1,stream);
	RSB_RW(row,&tsp->ns ,sizeof(tsp->ns ),1,stream);
        goto ret;
err:
        errval = RSB_ERR_INTERNAL_ERROR;
        RSB_ERROR("%s only %zd bytes instead of %zd !\n",row?"read":"wrote",hd,sh);
ret:
        return errval;
}

static rsb_err_t rsb__psp_rw(struct rsb_rsps_t * psp, FILE * stream, rsb_bool_t row)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        size_t hd = 0, sh = 0; /* have done, should have */

	RSB_RW(row,&psp->op_time,sizeof(psp->op_time),1,stream);
	RSB_RW(row,&psp->mkl_csr_op_time,sizeof(psp->mkl_csr_op_time),1,stream);
	RSB_RW(row,&psp->at_op_time,sizeof(psp->at_op_time),1,stream);
	RSB_RW(row,&psp->at_t,sizeof(psp->at_t),1,stream);
	RSB_RW(row,&psp->at_mkl_csr_op_time,sizeof(psp->at_mkl_csr_op_time),1,stream);
	RSB_RW(row,&psp->transA,sizeof(psp->transA),1,stream);
	RSB_RW(row,&psp->cmflops,sizeof(psp->cmflops),1,stream);
	RSB_RW(row,&psp->flagsA,sizeof(psp->flagsA),1,stream);
	RSB_RW(row,&psp->nsubm,sizeof(psp->nsubm),1,stream);
        RSB_RW(row,&psp->at_nsubm,sizeof(psp->at_nsubm),1,stream);
	RSB_RW(row,&psp->isa,sizeof(psp->isa),1,stream);
        RSB_RW(row,&psp->at_isa,sizeof(psp->at_isa),1,stream);
	RSB_RW(row,&psp->at_cn,sizeof(psp->at_cn),1,stream);
        RSB_RW(row,&psp->at_mkl_csr_cn,sizeof(psp->at_mkl_csr_cn),1,stream);
	RSB_RW(row,&psp->uc,sizeof(psp->uc),1,stream);
	RSB_RW(row,&psp->nrA,sizeof(psp->nrA),1,stream);
	RSB_RW(row,&psp->ncA,sizeof(psp->ncA),1,stream);
	RSB_RW(row,&psp->nnzA,sizeof(psp->nnzA),1,stream);
	RSB_RW(row,&psp->at_eps,sizeof(psp->at_eps),1,stream);

        errval = rsb__ts_rw(&psp->otpos,stream,row);
        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(ret,RSB_ERRM_ES);
        errval = rsb__ts_rw(&psp->btpos,stream,row);
        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(ret,RSB_ERRM_ES);
        errval = rsb__ts_rw(&psp->otpms,stream,row);
        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(ret,RSB_ERRM_ES);
        errval = rsb__ts_rw(&psp->btpms,stream,row);
        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(ret,RSB_ERRM_ES);

        goto ret;
err:
        errval = RSB_ERR_INTERNAL_ERROR;
        RSB_ERROR("%s only %zd bytes instead of %zd !\n",row?"read":"wrote",hd,sh);
ret:
        return errval;
}

rsb_bool_t rsb__file_exists(const rsb_char_t * RSB_RESTRICT filename)
{
	FILE*stream = NULL;

        stream = fopen(filename,"r");
        if(stream != NULL)
        {
                fclose(stream);
                return RSB_BOOL_TRUE;
        }
        return RSB_BOOL_FALSE;
}

#if RSB_PRD_WANT_TIMESTAMP
static rsb_err_t rsb__pr_rw_time(struct rsb_rspr_t * rsprp, FILE * stream, rsb_bool_t row)
{
	const rsb_err_t errval = RSB_ERR_NO_ERROR;

	rsb__rspr_rw( &rsprp->tbeg,sizeof(rsprp->tbeg), stream, row);
	rsb__rspr_rw( &rsprp->tend,sizeof(rsprp->tend), stream, row);

	return errval;
}
#endif /* RSB_PRD_WANT_TIMESTAMP */

#if RSB_PRD_WANT_MBW
static rsb_err_t rsb__mbw_es_rw(struct rsb_mbw_et_t * mbetp, FILE * stream, rsb_bool_t row)
{
	/* FIXME: shall be elsewhere */

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_int_t sni;
	struct rsb_mbw_es_t etz;
	struct rsb_mbw_et_t mbetz;
	rsb_bool_t bogus = ( mbetp == NULL ) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;

       	RSB_BZERO_P(&mbetz);
       	RSB_BZERO_P(&etz);

	if(bogus)
        	mbetp=&mbetz;

	rsb__rspr_rw( &mbetp->sn,           sizeof(mbetp->sn   ), stream, row);

	if( row == RSB_PR_RD )
	{
		if( ! bogus )
		{
			mbetp->et = rsb__calloc( mbetp->sn*sizeof( *mbetp->et ) );
			if( mbetp->et == NULL)
				goto ret; // FIXME
		}
		else
			mbetp->et = &etz;
	}

	for ( sni = 0; sni <  mbetp->sn ; ++sni )
	{
		rsb_int_t bi = bogus ? 0 : 1;
		rsb__rspr_rw(  &(mbetp->et[sni*bi].bw ), sizeof(mbetp->et[0].bw ), stream, row);
		rsb__rspr_rw(  &(mbetp->et[sni*bi].sz ), sizeof(mbetp->et[0].sz ), stream, row);
		rsb__rspr_rw(  &(mbetp->et[sni*bi].lvl), sizeof(mbetp->et[0].lvl), stream, row);
		rsb__rspr_rw(  &(mbetp->et[sni*bi].mbt), sizeof(mbetp->et[0].mbt), stream, row);
	}
	if( bogus )
		mbetp->sn = 0;

	for ( sni = 0; sni <  mbetp->sn ; ++sni )
	{
		//printf("%d %s %d %lg\n",sni,rsb__mbw_s2s(.mbt),mbetp->et[sni].lvl,mbetp->et[sni].bw);
	}

	goto ret;
//err:
	// FIXME!
ret:
	return errval;
}
#endif /* RSB_PRD_WANT_MBW */

rsb_err_t rsb__pr_save(const rsb_char_t * RSB_RESTRICT filename, /*const*/ void * RSB_RESTRICT rsprpv,
        rsb_char_t**RSB_RESTRICT filenamea, rsb_int_t*RSB_RESTRICT ca, const rsb_int_t*RSB_RESTRICT incXa, const rsb_int_t*RSB_RESTRICT incYa, const rsb_int_t*RSB_RESTRICT nrhsa, const rsb_type_t*RSB_RESTRICT typecodes, const rsb_int_t*RSB_RESTRICT ta, rsb_bool_t can_overwrite)
{
        /*
                Saves a performace record.
                FIXME: TODO: error handling can be improved.
                TODO: join common code with rsb__pr_load .
        */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        /*const*/ struct rsb_rspr_t * rsprp = rsprpv; /* FIXME: this shall be const */
	rsb_int_t noc = RSB_PR_NOC(rsprp);
	FILE*stream = NULL;
        int filenamei;
        struct rsb_rspra_t rspra;
        struct rsb_rspra_t * rsprap = NULL;
        rsb_byte_t * bbp = NULL; /* binary blob pointer */
        size_t bbo = 0, bbl = 0, bbs = 0; /* binary blob offset/length/skip */
        rsb_int_t idx;
        const rsb_char_t * const sgntr = RSB_RPR_FILE_HDR;
	const char rsb_prl_lcc = RSB_PRL_LCC_IE;

       	if(filename == NULL)
		stream = RSB_DEFAULT_FD;
	else
        {
                if(can_overwrite == RSB_BOOL_FALSE )
                if(rsb__file_exists(filename) != RSB_BOOL_FALSE)
                {
		        RSB_WARN("File %s already exists! Refusing to overwrite.\n",filename);
		        errval = RSB_ERR_INTERNAL_ERROR;
                        RSB_PERR_GOTO(err,RSB_ERRM_ES);
                }
		stream = fopen(filename,"wb");
        }
 
        if( stream == NULL )
        {
		errval = RSB_ERR_INTERNAL_ERROR;
	        RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

        RSB_BZERO_P(&rspra);

        rsprp->filenamebl = 0;
        bbs = sizeof(rspra) + sizeof(filenamea[0])*rsprp->filenamen;
        for(     filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
                rsprp->filenamebl += rsb__strlen(filenamea[filenamei]) + 1;
        rsprp->cabl = sizeof(*rsprap->ca)*rsprp->cn;
        rsprp->incXabl = sizeof(*rsprap->incXa)*rsprp->incXn;
        rsprp->incYabl = sizeof(*rsprap->incYa)*rsprp->incYn;
        rsprp->nrhsabl = sizeof(*rsprap->nrhsa)*rsprp->nrhsn;
        rsprp->typecodesbl = sizeof(*rsprap->typecodes)*( rsprp->ntypecodes + 1 );
        rsprp->tabl = ta ? sizeof(*rsprap->ta)*rsprp->tn : 0;
        bbl = rsprp->filenamebl  + rsprp->cabl + rsprp->incXabl + rsprp->incYabl + rsprp->nrhsabl + rsprp->typecodesbl + rsprp->tabl;

        fwrite(sgntr,RSB_RPR_FILE_HDL,1,stream);;

        errval = rsb__rsprp_rw(rsprp, stream, RSB_PR_WR);
        if(RSB_SOME_ERROR(errval))
                RSB_PERR_GOTO(err,RSB_ERRM_ES);

        for(idx=0;idx<noc;++idx)
        {
                struct rsb_rsps_t*psp = &(rsprp->psa[idx]);

                errval = rsb__psp_rw(psp, stream, RSB_PR_WR);
                if(RSB_SOME_ERROR(errval))
                        RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

	RSB_ASSERT(rsprp->filenamebl);

	rsprap = rsb__calloc( bbl + bbs );
	bbp = (void*) rsprap;
        if (! bbp ) { errval = RSB_ERR_ENOMEM; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
        bbp += sizeof(rspra);
        rspra.filenamea = (void*) bbp;
        bbp += sizeof(filenamea[0])*rsprp->filenamen;
        bbp += rsprp->filenamebl;
        rspra.ca = (void*) bbp;
        bbp += rsprp->cabl;
        rspra.incXa = (void*) bbp;
        bbp += rsprp->incXabl;
        rspra.incYa = (void*) bbp;
        bbp += rsprp->incYabl;
        rspra.nrhsa = (void*) bbp;
        bbp += rsprp->nrhsabl;
        rspra.typecodes = (rsb_char_t*) bbp;
        bbp += rsprp->typecodesbl;
        if(rsprp->tabl)
                        rspra.ta = (void*) bbp;
        bbp += rsprp->tabl;

        if(!rsprap)
        {
                errval = RSB_ERR_ENOMEM;
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }
        *rsprap = rspra;

        bbp = (void*) rsprap;

        bbo = 0;
        for(    filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
        {
                rsb__strcpy(bbp+bbs+bbo,filenamea[filenamei]);
                bbo += rsb__strlen(filenamea[filenamei]) + 1;
        }
	RSB_ASSERT(bbl);
        RSB_ASSERT(bbo == rsprp->filenamebl);
        RSB_RMEMCPY(rspra.ca       ,ca       ,rsprp->cabl       );
        RSB_RMEMCPY(rspra.incXa    ,incXa    ,rsprp->incXabl    );
        RSB_RMEMCPY(rspra.incYa    ,incYa    ,rsprp->incYabl    );
        RSB_RMEMCPY(rspra.nrhsa    ,nrhsa    ,rsprp->nrhsabl    );
        RSB_RMEMCPY(rspra.typecodes,typecodes,rsprp->typecodesbl);
        if(ta)
                        RSB_RMEMCPY(rspra.ta       ,ta       ,rsprp->tabl       );

        if( 1 != fwrite(((rsb_byte_t*)(rsprap))+bbs, bbl, 1, stream) )
        {
                errval = RSB_ERR_ENOMEM;
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

#if RSB_PRD_WANT_ENVIRONMENT
	errval = rsb__rspr_rw_env(rsprp, stream, RSB_PR_WR);
        if(RSB_SOME_ERROR(errval))
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif /* RSB_PRD_WANT_ENVIRONMENT */

#if RSB_PRD_WANT_MBW
	errval = rsb__mbw_es_rw(&rsprp->mbet,stream,RSB_PR_WR);
       	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif /* RSB_PRD_WANT_MBW */

#if RSB_PRD_WANT_TIMESTAMP
	errval = rsb__pr_rw_time(rsprp, stream,RSB_PR_WR);
       	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif /* RSB_PRD_WANT_TIMESTAMP */

       	if(filename == NULL)
		;
	else
		errval = ( EOF == fclose(stream) ) ? RSB_ERR_INTERNAL_ERROR : errval;
err:
	if(!RSB_SOME_ERROR(errval))
                RSB_PRL_SEP ("Saved a performance record of %d samples to %s\n",noc,filename);

	RSB_CONDITIONAL_FREE(rsprap);
        return errval;
}

static rsb_int_t rsb__pr_idx(const void*rsprpv, rsb_int_t filenamei, rsb_int_t ci, rsb_int_t incXi, rsb_int_t incYi, rsb_int_t nrhsi, rsb_int_t ntypecodei, rsb_int_t ti)
{
	/* 
	 * compute performance record index
	 * */
	const struct rsb_rspr_t * rsprp = rsprpv;
	rsb_int_t
	off6 = 1    * rsprp->tn,
	off5 = off6 * rsprp->ntypecodes,
	off4 = off5 * rsprp->nrhsn,
	off3 = off4 * rsprp->incYn,
	off2 = off3 * rsprp->incXn,
	off1 = off2 * rsprp->cn;

	return ti + ntypecodei * off6 + nrhsi * off5 + incYi * off4 + incXi * off3 + ci * off2 + filenamei * off1;
}

static rsb_err_t rsb__pr_load(const rsb_char_t * filename, struct rsb_rspr_t ** rsprpvp)
{
        /* 
                Loads a performace record.
                Overwrites the target pointer.
                If *rsprpvp, then only the first struct will be load.

                TODO: a header, e.g.: "RPR " would do no wrong.
                TODO: a header signature check.
                TODO: error handling can be improved.
                TODO: relax the typecode checks, as values from a differently configured build may have to be be read.
         */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_rspr_t rspr;
	struct rsb_rspr_t * rsprp = NULL;
	struct rsb_rspr_t ** rsprpp = rsprpvp;
	FILE*stream = NULL;
	rsb_int_t noc = 0; /* number of combinations */
        rsb_char_t sgntr [ RSB_RPR_FILE_HDL ];
        struct rsb_rspra_t rspra;
        struct rsb_rspra_t * rsprap = NULL;
        rsb_byte_t * bbp = NULL; /* binary blob pointer */
        size_t bbo = 0, bbl = 0, bbs = 0; /* binary blob offset/length/skip */
	rsb_int_t filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti;
	rsb_int_t vtag = RSB_PRD_VERSION_BAD;
	size_t pad = 0; /* padding to 4 bytes alignment */
	size_t bbpp = 0; /* bpp, padded */

        RSB_BZERO_P(&rspr);

       	if(filename == NULL)
		stream = stdin;
	else
		stream = fopen(filename,"rb");

        if( stream == NULL )
        {
		errval = RSB_ERR_INTERNAL_ERROR;
	        RSB_PERR_GOTO(err,"Error opening performance record file \"%s\" for reading !\n",filename?filename:"stdin");
        }

        if( 1 != fread(sgntr,RSB_RPR_FILE_HDL,1,stream) )
        {
		errval = RSB_ERR_INTERNAL_ERROR;
                RSB_PERR_GOTO(err,"Unable to read header!\n");
        }

        if ( 0 == strncmp(sgntr,RSB_RPR_FILE_HDR_V0,RSB_RPR_FILE_HDL) )
		vtag = 0;
	else if( 0 == strncmp(sgntr,RSB_RPR_FILE_HDR,RSB_RPR_FILE_HDL) ) 
		vtag = 1;
	if(vtag == RSB_PRD_VERSION_BAD)
        {
                /* TODO: need support for different versions ... */
		errval = RSB_ERR_INTERNAL_ERROR;
                RSB_PERR_GOTO(err,"File decoding error!\n");
        }

        rsprp = &rspr;
        errval = rsb__rsprp_rw(rsprp, stream, RSB_PR_RD);
	if(RSB_SOME_ERROR(errval))
                RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(*rsprpp)
        {
                /* ok to return now */
	        **rsprpp = rspr;
                goto cret;
        }

        errval = rsb__pr_alloc(&rsprp, &rspr, rspr.filenamen, rspr.cn, rspr.incXn, rspr.incYn, rspr.nrhsn, rspr.ntypecodes, rspr.tn);
	if(RSB_SOME_ERROR(errval))
                RSB_PERR_GOTO(err,RSB_ERRM_ES);

        if(rsprp->csf == 0)
	{
                errval = RSB_ERR_CORRUPT_INPUT_DATA;
                RSB_PERR_GOTO(err,"Seems like the input performance record is empty!\n");
	}
        RSB_ASSERT(rsprp->csf);

	noc = RSB_PR_NOC(rsprp);

        //for(idx=0;idx<noc;++idx)
	for(     filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
	for(ci=0;ci<rsprp->cn;++ci)
	for(     incXi=0;     incXi<rsprp->incXn     ;++incXi     )
	for(     incYi=0;     incYi<rsprp->incYn     ;++incYi     )
	for(     nrhsi=0;     nrhsi<rsprp->nrhsn     ;++nrhsi     )
	for(typecodesi=0;typecodesi<rsprp->ntypecodes;++typecodesi)
	for(ti=0;ti<rsprp->tn;++ti)
	{
		const size_t idx = rsb__pr_idx(rsprp, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti);
	        struct rsb_rsps_t*psp = &(rsprp->psa[idx]);

                errval = rsb__psp_rw(psp, stream, RSB_PR_RD);
               	if(RSB_SOME_ERROR(errval))
                        RSB_PERR_GOTO(err,RSB_ERRM_ES);
                /* compatibility patch for old (pre-1.2) values of RSB_TRANSPOSITION_*  */
		if(psp->transA==0x00) psp->transA = RSB_TRANSPOSITION_N;
		if(psp->transA==0x01) psp->transA = RSB_TRANSPOSITION_T;
		if(psp->transA==0x02) psp->transA = RSB_TRANSPOSITION_C;

                if ( rsb__getenv("RSB_PR_RD_NULLIFY_FILENAMEI") ) /* proof of concept */
                {
                        int nidx = rsb__util_atoi(rsb__getenv("RSB_PR_RD_NULLIFY_FILENAMEI"));
                        if ( idx / (noc/rspr.filenamen) == nidx )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_FILENAMEI") ) /* proof of concept */
                {
                        int nidx = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_FILENAMEI"));
                        if ( idx / (noc/rspr.filenamen) != nidx )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_NULLIFY_SAMPLEIDX") ) /* proof of concept */
                {
                        int nidx = rsb__util_atoi(rsb__getenv("RSB_PR_RD_NULLIFY_SAMPLEIDX"));
                        if ( idx == nidx )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_SAMPLEIDX") ) /* proof of concept */
                {
                        int nidx = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_SAMPLEIDX"));
                        if ( idx != nidx )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_TRANSA") ) /* proof of concept */
                {
                        rsb_trans_t no_transA = (*rsb__getenv("RSB_PR_RD_RESTRICT_TRANSA"));
                        if ( psp->transA != no_transA )
			       	psp->uc = 0;
                }
		/* TODO: begin of simplifiable checks */
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NR_MIN") ) /* proof of concept */
                {
                        rsb_nnz_idx_t min_nr = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NR_MIN"));
                        if ( psp->nrA < min_nr )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NR_MAX") ) /* proof of concept */
                {
                        rsb_nnz_idx_t max_nr = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NR_MAX"));
                        if ( psp->nrA > max_nr && max_nr > 0 )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NC_MIN") ) /* proof of concept */
                {
                        rsb_nnz_idx_t min_nc = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NC_MIN"));
                        if ( psp->ncA < min_nc )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NC_MAX") ) /* proof of concept */
                {
                        rsb_nnz_idx_t max_nc = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NC_MAX"));
                        if ( psp->ncA > max_nc && max_nc > 0 )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NNZ_MIN") ) /* proof of concept */
                {
                        rsb_nnz_idx_t min_nnz = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NNZ_MIN"));
                        if ( psp->nnzA < min_nnz )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NNZ_MAX") ) /* proof of concept */
                {
                        rsb_nnz_idx_t max_nnz = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NNZ_MAX"));
                        if ( psp->nnzA > max_nnz && max_nnz > 0 )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NSUBM_MIN") ) /* proof of concept */
                {
                        rsb_nnz_idx_t min_nsubm = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NSUBM_MIN"));
                        if ( psp->nsubm > 0   && psp->nsubm    < min_nsubm )
			       	psp->uc = 0;
                        if ( psp->at_nsubm> 0 && psp->at_nsubm < min_nsubm )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NSUBM_MAX") ) /* proof of concept */
                {
                        rsb_nnz_idx_t max_nsubm = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NSUBM_MAX"));
                        if ( psp->nsubm > 0   && psp->nsubm    > max_nsubm )
			       	psp->uc = 0;
                        if ( psp->at_nsubm> 0 && psp->at_nsubm > max_nsubm )
			       	psp->uc = 0;
                }
		/* TODO: end of simplifiable checks */
                if ( rsb__getenv("RSB_PR_RD_NULLIFY_TRANSA") ) /* proof of concept */
                {
                        rsb_trans_t no_transA = (*rsb__getenv("RSB_PR_RD_NULLIFY_TRANSA"));
                        if ( psp->transA == no_transA )
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_NULLIFY_NRHSI") ) /* proof of concept */
                {
                        int nidx = rsb__util_atoi(rsb__getenv("RSB_PR_RD_NULLIFY_NRHSI"));
			if(nrhsi==nidx)
			       	psp->uc = 0;
                }
                if ( rsb__getenv("RSB_PR_RD_RESTRICT_NRHSI") ) /* proof of concept */
                {
                        int nidx = rsb__util_atoi(rsb__getenv("RSB_PR_RD_RESTRICT_NRHSI"));
			if(nrhsi!=nidx)
			       	psp->uc = 0;
                }
        }

        RSB_BZERO_P(&rspra);
	
	pad = (4-(rspr.filenamebl%4))%4; /* 0..3 */
	bbl = rspr.filenamebl + pad + rspr.cabl + rspr.incXabl + rspr.incYabl + rspr.nrhsabl + rspr.typecodesbl + rspr.tabl;
	bbs = sizeof(rspra) + sizeof(rspra.filenamea[0])*rspr.filenamen;
	
        rsprap = rsb__calloc( bbl + bbs );
        bbp = (void*) rsprap;
	if (! bbp ) { errval = RSB_ERR_ENOMEM; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
	bbp += sizeof(rspra);
	rspra.filenamea = (void*) bbp;
	bbp += sizeof(rspra.filenamea[0])*rspr.filenamen;
	bbp += rspr.filenamebl;
	bbpp = (bbp - ((rsb_byte_t*) rsprap)) - sizeof(rspra);

	bbp += pad;
	rspra.ca = (void*) bbp;
	bbp += rspr.cabl;
	rspra.incXa = (void*) bbp;
	bbp += rspr.incXabl;
	rspra.incYa = (void*) bbp;
	bbp += rspr.incYabl;
	rspra.nrhsa = (void*) bbp;
	bbp += rspr.nrhsabl;
	rspra.typecodes = (rsb_char_t*) (void*) bbp;
	bbp += rspr.typecodesbl;
	if(rspr.tabl)
	        rspra.ta = (void*) bbp;
	bbp += rspr.tabl;
	
        RSB_ASSERT(rspr.filenamebl);
        RSB_ASSERT(rspr.nrhsabl);
	*rsprap = rspra;
	
	if( 1 != fread( ((rsb_byte_t*)rsprap)+bbs, bbpp, 1, stream) )
	{
                errval = RSB_ERR_INTERNAL_ERROR;
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	if( 1 != fread( ((rsb_byte_t*)rsprap)+bbs+bbpp+pad, bbl-bbpp-pad, 1, stream) )
	{
                errval = RSB_ERR_INTERNAL_ERROR;
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	bbp = (void*) rsprap;
	bbp += bbs;
	bbo = 0;

        for(     filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
	{
                rspra.filenamea[filenamei] = (rsb_char_t*) bbp + bbo;
                bbo += rsb__strlen(bbp + bbo) + 1;
	}

	if( vtag == 0 )
	{
		/* reading a pre-RSB_PRD_VERSION==1 version file */
	}
	else
	{
#if RSB_PRD_WANT_ENVIRONMENT
		errval = rsb__rspr_rw_env(rsprp, stream, RSB_PR_RD);
	        if(RSB_SOME_ERROR(errval))
        	        RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif /* RSB_PRD_WANT_ENVIRONMENT */
#if RSB_PRD_WANT_MBW
		errval = rsb__mbw_es_rw(&(rsprp->mbet), stream, RSB_PR_RD);
       		if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif /* RSB_PRD_WANT_MBW */
#if RSB_PRD_WANT_TIMESTAMP
		errval = rsb__pr_rw_time(rsprp, stream,RSB_PR_RD);
       		if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif /* RSB_PRD_WANT_TIMESTAMP */
	}

	rsprp->rsprap = rsprap;

err:
 	if (rsprp)
 	        RSB_ASSIGN_IF(rsprpp,rsprp)

 	if(rsprpp == NULL)
                RSB__PR_FREE(rsprp); /* bogus load ... */
cret:
       	if(filename == NULL)
		;
	else
		errval = ( stream && EOF == fclose(stream) ) ? RSB_ERR_INTERNAL_ERROR : errval;

        return errval;
}

#if 0
static rsb_err_t rsb__pr_dumpfile(const rsb_char_t *filename)
{
	/* Obsoleted by rsb__pr_dumpfiles. */
        /*
                Due to subtle build configuration dependent problems, it's better to declare this as experimental.
                FIXME: error handling is insufficient.
        */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_rspr_t * rsprp = NULL;

        errval = rsb__pr_load(filename,&rsprp);
	if(RSB_SOME_ERROR(errval))
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
        RSB_PRL_SEP("\n");
	errval = rsb__pr_dump(rsprp, rsprp->rsprap->filenamea, rsprp->rsprap->ca, rsprp->rsprap->incXa, rsprp->rsprap->incYa, rsprp->rsprap->nrhsa, rsprp->rsprap->typecodes, NULL, NULL );
        RSB_PRL_SEP("\n");
err:
	RSB__PR_FREE(rsprp);
        return errval;
}
#endif

#if 0
static rsb_err_t rsb__pr_sort(struct rsb_rspr_t * rsprp)
{
	/* 
	 * TODO: this function is yet incomplete.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

        if(!rsprp->rsprap)
                goto err;

        if(rsprp->rsprap->filenamea);
        if(rsprp->rsprap->ca);
        if(rsprp->rsprap->incXa);
        if(rsprp->rsprap->incYa);
        if(rsprp->rsprap->nrhsa);
        if(rsprp->rsprap->typecodes);
        if(rsprp->rsprap->ta);

	rsb_int_t  filenamen,   cn,   incXn,   incYn,   nrhsn,  ntypecodes,   tn, csf /* count so far */;
err:
        return errval;
}
#endif

static rsb_err_t rsb__pr_cmp(/*const*/ struct rsb_rspr_t * rspr0p, const struct rsb_rspr_t * rspr1p, int wr)
{
	/* 
	 * wr = 1 ratio, 2 diff.
	 * TODO: this function is yet experimental.
	 * TODO: input error checking is missing.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_rsps_t dps; /* difference performance sample */
        rsb_int_t idx = 0, noc = 0, pr = 0, csf = 0;
	const char rsb_prl_lcc = RSB_PRL_LCC_IE;

        RSB_ASSERT(rspr0p);
        RSB_ASSERT(rspr1p);

        RSB_ASSERT(rspr0p->filenamen == rspr1p->filenamen);
        RSB_ASSERT(rspr0p->filenamebl == rspr1p->filenamebl);
        RSB_ASSERT(rspr0p->cn == rspr1p->cn);
        /* RSB_ASSERT(rspr0p->csf == rspr1p->csf); */
        RSB_ASSERT(rspr0p->csf > 0 && rspr1p->csf > 0);
        RSB_ASSERT(rspr0p->cabl == rspr1p->cabl);
        RSB_ASSERT(rspr0p->incXn == rspr1p->incXn );
        RSB_ASSERT(rspr0p->incXabl == rspr1p->incXabl );
        RSB_ASSERT(rspr0p->incYn == rspr1p->incYn );
        RSB_ASSERT(rspr0p->incYabl == rspr1p->incYabl );
        RSB_ASSERT(rspr0p->nrhsn == rspr1p->nrhsn );
        RSB_ASSERT(rspr0p->nrhsabl == rspr1p->nrhsabl );
        RSB_ASSERT(rspr0p->ntypecodes == rspr1p->ntypecodes );
        RSB_ASSERT(rspr0p->typecodesbl == rspr1p->typecodesbl );
        RSB_ASSERT(rspr0p->tn == rspr1p->tn );
        RSB_ASSERT(rspr0p->tabl == rspr1p->tabl );

	csf = noc = RSB_PR_NOC(rspr0p);

        if(noc != RSB_PR_NOC(rspr1p))
        {
		errval = RSB_ERR_INTERNAL_ERROR;
	        RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

	if( rspr0p->csf != rspr1p->csf )
	{
		csf = RSB_MIN(rspr0p->csf, rspr1p->csf);
               	RSB_PRL("Out of %d samples, one record has %d and the other %d (incomplete record ?). Limiting to the minimum of the two (EXPERIMENTAL!).\n", noc, rspr0p->csf, rspr1p->csf);
	}

        rspr0p->ror = RSB_BOOL_TRUE;

        switch(wr){
        case(RSB_PRD_CMP_DIV):
        for(idx=0;idx<csf;++idx)
        if(rspr0p->psa[idx].uc && rspr1p->psa[idx].uc)
        {
	        dps = rspr0p->psa[idx];
                dps.op_time = rspr0p->psa[idx].op_time / rspr1p->psa[idx].op_time;
                dps.mkl_csr_op_time = rspr0p->psa[idx].mkl_csr_op_time / rspr1p->psa[idx].mkl_csr_op_time;
                dps.at_mkl_csr_op_time = rspr0p->psa[idx].at_mkl_csr_op_time / rspr1p->psa[idx].at_mkl_csr_op_time;
                dps.at_t = rspr0p->psa[idx].at_t / rspr1p->psa[idx].at_t;
                dps.at_op_time = rspr0p->psa[idx].at_op_time / rspr1p->psa[idx].at_op_time ;
	        dps.at_eps = 0;
	        rspr0p->psa[idx] = dps;
                ++pr;
        }
        break;
        case(RSB_PRD_CMP_DIFF):
        for(idx=0;idx<csf;++idx)
        if(rspr0p->psa[idx].uc && rspr1p->psa[idx].uc)
        {
	        dps = rspr0p->psa[idx];
                dps.op_time = rspr0p->psa[idx].op_time - rspr1p->psa[idx].op_time;
                dps.mkl_csr_op_time = rspr0p->psa[idx].mkl_csr_op_time - rspr1p->psa[idx].mkl_csr_op_time;
                dps.at_mkl_csr_op_time = rspr0p->psa[idx].at_mkl_csr_op_time - rspr1p->psa[idx].at_mkl_csr_op_time;
                dps.at_t = rspr0p->psa[idx].at_t - rspr1p->psa[idx].at_t;
                dps.at_op_time = rspr0p->psa[idx].at_op_time - rspr1p->psa[idx].at_op_time ;
                /* and: */
                dps.at_eps = rspr0p->psa[idx].at_eps - rspr1p->psa[idx].at_eps ;
                dps.at_isa = rspr0p->psa[idx].at_isa - rspr1p->psa[idx].at_isa ;
	        rspr0p->psa[idx] = dps;
                ++pr;
        }
        break;
        default:
		errval = RSB_ERR_INTERNAL_ERROR;
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

        if(pr == 0)
	{
		RSB_WARN("No pair of samples has been found to be conformable!\n");
		// errval = RSB_ERR_BADARGS; goto err;
	}
err:
        return errval;
}

static rsb_err_t rsb__pr_merge(struct rsb_rspr_t * rspr0p, const struct rsb_rspr_t * rspr1p)
{
	/* 
         * Joins two records by appending one to another and returning a new one.
	 * TODO: this function is yet experimental.
         * TODO: shall extend the joining mechanism to arbitrarily many records. 
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_int_t noc0, noc1;
        rsb_int_t idx;
	const char rsb_prl_lcc = RSB_PRL_LCC_IE;

        RSB_PRL("Warning: joining assuming ALL parameters are conformant (except filenames)\n");

        RSB_ASSERT(rspr0p);
        RSB_ASSERT(rspr1p);

        RSB_ASSERT(rspr0p->filenamen == rspr1p->filenamen);
        RSB_ASSERT(rspr0p->filenamebl == rspr1p->filenamebl);
        RSB_ASSERT(rspr0p->csf>0);
        RSB_ASSERT(rspr1p->csf>0);
        RSB_ASSERT(rspr0p->cn == rspr1p->cn);
        RSB_ASSERT(rspr0p->cabl == rspr1p->cabl);
        RSB_ASSERT(rspr0p->incXn == rspr1p->incXn );
        RSB_ASSERT(rspr0p->incXabl == rspr1p->incXabl );
        RSB_ASSERT(rspr0p->incYn == rspr1p->incYn );
        RSB_ASSERT(rspr0p->incYabl == rspr1p->incYabl );
        RSB_ASSERT(rspr0p->nrhsn == rspr1p->nrhsn );
        RSB_ASSERT(rspr0p->nrhsabl == rspr1p->nrhsabl );
        RSB_ASSERT(rspr0p->ntypecodes == rspr1p->ntypecodes );
        RSB_ASSERT(rspr0p->typecodesbl == rspr1p->typecodesbl );
        RSB_ASSERT(rspr0p->tn == rspr1p->tn );
        RSB_ASSERT(rspr0p->tabl == rspr1p->tabl );

        if(!rspr0p->rsprap)
                goto err;
        if(!rspr1p->rsprap)
                goto err;

	noc0 = RSB_PR_NOC(rspr0p);
	noc1 = RSB_PR_NOC(rspr1p);

        if( noc0 != noc1 )
        {
                errval = RSB_ERR_INTERNAL_ERROR;
                RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

        RSB_ASSERT( noc0 == noc1 );

        for(idx=0;idx<noc0;++idx)
        {
                struct rsb_rsps_t*psp = &(rspr0p->psa[idx]);

                if( psp->uc == 0 )
                        *psp = (rspr1p->psa[idx]);
        }
err:
        return errval;
}

static rsb_err_t rsb__pr_join(struct rsb_rspr_t ** rsprpp, const struct rsb_rspr_t * rspr0p, const struct rsb_rspr_t * rspr1p)
{
	/* 
         * Joins two records by appending one to another and returning a new one.
	 * TODO: this function is yet experimental.
         * TODO: shall extend the joining mechanism to arbitrarily many records. 
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        struct rsb_rspr_t * rsprp = NULL;
        struct rsb_rspra_t rspra;
	rsb_int_t noc0, noc1;
        size_t bbo = 0, bbl = 0, bbs = 0; /* binary blob offset/length/skip */
        rsb_byte_t * bbp = NULL; /* binary blob pointer */
	rsb_int_t filenamei;
	const char rsb_prl_lcc = RSB_PRL_LCC_IE;

        RSB_PRL("Warning: joining assuming ALL parameters are conformant (except filenames)\n");

        RSB_ASSERT(rspr0p);
        RSB_ASSERT(rspr1p);

        /* RSB_ASSERT(rspr0p->filenamen == rspr1p->filenamen); */
        /* RSB_ASSERT(rspr0p->filenamebl == rspr1p->filenamebl); */
        RSB_ASSERT(rspr0p->csf>0);
        RSB_ASSERT(rspr1p->csf>0);
        RSB_ASSERT(rspr0p->cn == rspr1p->cn);
        RSB_ASSERT(rspr0p->cabl == rspr1p->cabl);
        RSB_ASSERT(rspr0p->incXn == rspr1p->incXn );
        RSB_ASSERT(rspr0p->incXabl == rspr1p->incXabl );
        RSB_ASSERT(rspr0p->incYn == rspr1p->incYn );
        RSB_ASSERT(rspr0p->incYabl == rspr1p->incYabl );
        RSB_ASSERT(rspr0p->nrhsn == rspr1p->nrhsn );
        RSB_ASSERT(rspr0p->nrhsabl == rspr1p->nrhsabl );
        RSB_ASSERT(rspr0p->ntypecodes == rspr1p->ntypecodes );
        RSB_ASSERT(rspr0p->typecodesbl == rspr1p->typecodesbl );
        RSB_ASSERT(rspr0p->tn == rspr1p->tn );
        RSB_ASSERT(rspr0p->tabl == rspr1p->tabl );

        if(!rspr0p->rsprap)
                goto err;
        if(!rspr1p->rsprap)
                goto err;

        errval = rsb__pr_alloc(&rsprp, rspr0p, rspr0p->filenamen + rspr1p->filenamen, rspr0p->cn, rspr0p->incXn, rspr0p->incYn, rspr0p->nrhsn, rspr0p->ntypecodes, rspr0p->tn);

	if(RSB_SOME_ERROR(errval))
                RSB_PERR_GOTO(err,RSB_ERRM_ES);

        rsprp->csf = rspr0p->csf + rspr1p->csf;
        rsprp->filenamen = rspr0p->filenamen + rspr1p->filenamen;
        rsprp->filenamebl = rspr0p->filenamebl + rspr1p->filenamebl;

	noc0 = RSB_PR_NOC(rspr0p);
	noc1 = RSB_PR_NOC(rspr1p);

        RSB_RMEMCPY(rsprp->psa+0*noc0,rspr0p->psa,sizeof(*rsprp->psa)*noc0);
        RSB_RMEMCPY(rsprp->psa+1*noc0,rspr1p->psa,sizeof(*rsprp->psa)*noc1);

        RSB_BZERO_P(&rspra);
	
	bbl = rsprp->filenamebl + rsprp->cabl + rsprp->incXabl + rsprp->incYabl + rsprp->nrhsabl + rsprp->typecodesbl + rsprp->tabl;
	bbs = sizeof(rspra) + sizeof(rspra.filenamea[0])*rsprp->filenamen;
	
        /* FIXME: encapsulate the follwing in a function */
        rsprp->rsprap = rsb__calloc( bbl + bbs );
        bbp = (void*) rsprp->rsprap;
	if (! bbp ) { errval = RSB_ERR_ENOMEM; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
	bbp += sizeof(rspra);
	rspra.filenamea = (void*) bbp;
	bbp += sizeof(rspra.filenamea[0])*rsprp->filenamen;
	bbp += rsprp->filenamebl;
	rspra.ca = (void*) bbp;
	bbp += rsprp->cabl;
	rspra.incXa = (void*) bbp;
	bbp += rsprp->incXabl;
	rspra.incYa = (void*) bbp;
	bbp += rsprp->incYabl;
	rspra.nrhsa = (void*) bbp;
	bbp += rsprp->nrhsabl;
	rspra.typecodes = (rsb_char_t*) (void*) bbp;
	bbp += rsprp->typecodesbl;
	if(rsprp->tabl)
	        rspra.ta = (void*) bbp;
	bbp += rsprp->tabl;
	
        RSB_ASSERT(rsprp->filenamebl);
        RSB_ASSERT(rsprp->nrhsabl);
	
	bbp = (void*) rsprp->rsprap;
	bbp += bbs;
	bbo = 0;

        for(     filenamei=0;     filenamei<rspr0p->filenamen ;++filenamei     )
	{
                rspra.filenamea[rspr0p->filenamen*0+filenamei] = (rsb_char_t*) bbp + bbo;
                rsb__strcpy(bbp+bbo,rspr0p->rsprap->filenamea[filenamei]);
                bbo += rsb__strlen(bbp + bbo) + 1;
	}

        for(     filenamei=0;     filenamei<rspr1p->filenamen ;++filenamei     )
	{
                rspra.filenamea[rspr0p->filenamen*1+filenamei] = (rsb_char_t*) bbp + bbo;
                rsb__strcpy(bbp+bbo,rspr1p->rsprap->filenamea[filenamei]);
                bbo += rsb__strlen(bbp + bbo) + 1;
	}

        RSB_ASSERT(bbo == rsprp->filenamebl);
        RSB_RMEMCPY(rspra.ca       ,rspr0p->rsprap->ca       ,rsprp->cabl       );
        RSB_RMEMCPY(rspra.incXa    ,rspr0p->rsprap->incXa    ,rsprp->incXabl    );
        RSB_RMEMCPY(rspra.incYa    ,rspr0p->rsprap->incYa    ,rsprp->incYabl    );
        RSB_RMEMCPY(rspra.nrhsa    ,rspr0p->rsprap->nrhsa    ,rsprp->nrhsabl    );
        RSB_RMEMCPY(rspra.typecodes,rspr0p->rsprap->typecodes,rsprp->typecodesbl);

        *rsprp->rsprap = rspra;
	RSB_ASSIGN_IF(rsprpp,rsprp)

        /* TODO: may also adjoin the remaining arrays, extending to cabl, incXabl, incYabl, nrhsabl, typecodesbl, tabl */
err:
        return errval;
}

static char *rsb__rindex(const char *s, int c)
{
	/* return a pointer to the last occurrence of the character c in the string s. assumes '\0' is there. */
#if RSB_HAVE_STRRCHR
	return strrchr(s,c);
#elif RSB_HAVE_RINDEX
	return rindex(s,c);
#else /* RSB_HAVE_RINDEX */
#error "Need working rindex() or strrchr()!"
#endif /* RSB_HAVE_STRRCHR */
}

rsb_err_t rsb__pr_dumpfiles(const rsb_char_t **argv, const int argc)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
        int ds = RSB_PRD_CMP_DFLT;
        struct rsb_rspr_t * rsprp = NULL;
        int argi;
	rsb_int_t noc = 0, noc0; /* number of combinations */
	char rsb_prl_lcc = RSB_PRL_LCC_IE ;
	rsb_char_t dirname[RSB_MAX_FILENAME_LENGTH];

	dirname[0] = RSB_NUL;

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))!=RSB_ERR_NO_ERROR)
		goto err;

        if( !argv || argc < 1 || strlen(argv[0]) < 1 )
        {
	        RSB_PRL("No performance record files to dump !? Please specify at least one.\n");
	        RSB_PRL("Consider further options, specifiable via environment variables:\n");
		RSB_PRL("# begin of help message\n");
		RSB_PRL("# This feature of librsb is not ufficially supported.\n");
		RSB_PRL("# threshold (expressed as ratio) between values:\n");
		RSB_PRL("RSB_CMP_THR # nearly same threshold\n");
		RSB_PRL("RSB_APE_THR # close values threshold\n");
		RSB_PRL("RSB_RLD_THR # relevant difference threshold\n");
		RSB_PRL("RSB_HUD_THR # huge difference threshold\n");
		RSB_PRL("RSB_PRD_STYLE_PLT_FMT # (if RSB_PR_SR=2) plot file format: EPS if set, PNG otherwise\n");
		RSB_PRL("RSB_PRD_STYLE_PLT_PFN # (if RSB_PR_SR=2) plot file name\n");
		RSB_PRL("RSB_PR_FSEPSTR # Field separator string\n");
		RSB_PRL("RSB_PR_ENDLSTR # End of line separator string\n");
		RSB_PRL("RSB_PR_PRL_CC  # Beginning of line comment char\n");
		RSB_PRL("RSB_PR_PRL_LCC # Line Comment Character\n");
		RSB_PRL("RSB_PR_PRL_TCS # Table Comment String\n");
		RSB_PRL("RSB_PR_WLTC # If > 0 and RSB_PR_SR=0, will emit LaTeX tables	(setting accordingly RSB_PR_PRL_LCC, RSB_PR_PRL_TCS, RSB_PR_ENDLSTR, RSB_PR_FSEPSTR); if > 1 output will be colored\n");
                RSB_PRL("RSB_PR_MULTIDUMP #  %d=dump %d=auto/append %d=ratio %d=diff %d=merge.\n",RSB_PRD_CMP_MDUMP,RSB_PRD_CMP_DFLT,RSB_PRD_CMP_DIV,RSB_PRD_CMP_DIFF,RSB_PRD_CMP_APPEND);
		RSB_PRL("RSB_PR_RD_NULLIFY_FILENAMEI # exclude a matrix' index\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_FILENAMEI # restrict to one matrix' index\n");
		RSB_PRL("RSB_PR_RD_NULLIFY_TRANSA # exclude a transposition\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_TRANSA # restrict to one transposition\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NR_MIN # restrict to min of nr\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NR_MAX # restrict to max of nr\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NC_MIN # restrict to min of nc\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NC_MAX # restrict to max of nc\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NNZ_MIN # restrict to min of nnz\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NNZ_MAX # restrict to max of nnz\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NSUBM_MIN # restrict to min of nsubm\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NSUBM_MAX # restrict to max of nsubm\n");
		RSB_PRL("RSB_PR_RD_NULLIFY_NRHSI # exclude a nrhs index\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_NRHSI # restrict to one nrhs index\n");
		RSB_PRL("RSB_PR_RD_NULLIFY_SAMPLEIDX # exclude a matrix' index\n");
		RSB_PRL("RSB_PR_RD_RESTRICT_SAMPLEIDX # restrict to one matrix' index\n");
		RSB_PRL("RSB_PR_ONLY_TOTAL_TABLE # only the total table, not the 'limited' slices\n");
		RSB_PRL("RSB_PR_SAVE_MULTIDUMP # output performance record filename\n");
		RSB_PRL("RSB_PR_SR # 0 for table output, 1 for comparison table output, 2 for plot\n");
		RSB_PRL("RSB_PR_ENV # print out environment variables\n");
		RSB_PRL("RSB_PR_MBW # print out memory bandwidth benchmark info\n");
		RSB_PRL("# end of help message\n");
		goto err;
        }
        RSB_PRL_SEP("\n");

        if(argc > 1)
        {
                RSB_PRL("You can control multiple files dump with RSB_PR_MULTIDUMP= %d=dump %d=auto/append %d=ratio %d=diff %d=merge.\n",RSB_PRD_CMP_MDUMP,RSB_PRD_CMP_DFLT,RSB_PRD_CMP_DIV,RSB_PRD_CMP_DIFF,RSB_PRD_CMP_APPEND);
        }

        if( NULL != rsb__getenv("RSB_PR_MULTIDUMP") )
                ds = rsb__util_atoi(rsb__getenv("RSB_PR_MULTIDUMP"));
        else
        {
                for(argi=0;argi<argc;++argi)
                {
                        struct rsb_rspr_t rspr;
                        rsprp = &rspr;
                        RSB_BZERO_P(&rspr);
                        errval = rsb__pr_load(argv[argi],&rsprp);
		        if(RSB_SOME_ERROR(errval))
			{
                        	rsprp = NULL;
				RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}
                        noc = RSB_PR_NOC(&rspr); /* FIXME: this is quite a weak test.. */
                        if(argi == 0) noc0 = noc;
                        if(argi  > 0 && noc != noc0) { noc = 0; break; }
                        rsprp = NULL;
                }

                if ( argc > 1 &&  noc != 0 )
                {
        	        RSB_PRL("Warning: hazarding the guess you are working with complementary performance record files, therefore attempting merging!.\n");
                        ds = RSB_PRD_CMP_APPEND;
                }
        }

        if(ds < RSB_PRD_CMP_MDUMP || ds > RSB_PRD_CMP_APPEND)
        {
		RSB_ERROR("Set RSB_PR_MULTIDUMP to a bad value !\n");
		return RSB_ERR_BADARGS;
                goto err;
        }

        if(ds >= RSB_PRD_CMP_DFLT)
        {
		const rsb_char_t *const fprfn = argv[0];
                errval = rsb__pr_load(fprfn ,&rsprp);
	        if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
                sprintf(dirname,"%s",fprfn);
		if(rsb__rindex(dirname,'.'))
			*rsb__rindex(dirname,'.') = RSB_NUL;
		else
			strcat(dirname,".dir");
        }

        if(argc > 1)
                RSB_PRL("Will display summary of %d performance records\n", argc);

        if(ds == RSB_PRD_CMP_DFLT || ds == RSB_PRD_CMP_APPEND)
        {
                for(argi=1;argi<argc;++argi)
                if(ds==RSB_PRD_CMP_DFLT)
                {
                        struct rsb_rspr_t * rspr1p = NULL, * rspr0p = NULL;
                        RSB_PRL("Will append performance records of file %d/%d: %s to that of %s.\n",argi+1,argc,argv[argi],argv[0]);
                        errval = rsb__pr_load(argv[argi],&rspr0p);
                        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                        errval = rsb__pr_join(&rspr1p, rsprp, rspr0p);
                        RSB__PR_FREE(rspr0p);
                        RSB__PR_FREE(rsprp);
                        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                        rsprp = rspr1p;
                        rspr1p = NULL;
     }
                else
                if(ds==RSB_PRD_CMP_APPEND)
                {
                        struct rsb_rspr_t * rspr0p = NULL;
                        RSB_PRL("Will merge performance records of file %d/%d: %s to that of %s.\n",argi+1,argc,argv[argi],argv[0]);
                        errval = rsb__pr_load(argv[argi],&rspr0p);
                        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                        errval = rsb__pr_merge(rsprp, rspr0p);
                        RSB__PR_FREE(rspr0p);
                        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                }
                errval = rsb__pr_dump(rsprp, rsprp->rsprap->filenamea, rsprp->rsprap->ca, rsprp->rsprap->incXa, rsprp->rsprap->incYa, rsprp->rsprap->nrhsa, rsprp->rsprap->typecodes, NULL, dirname);
                if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                
                if ( rsb__getenv("RSB_PR_SAVE_MULTIDUMP") )
                {
                        const char * of = rsb__getenv("RSB_PR_SAVE_MULTIDUMP");
                        /* errval = rsb__pr_save(of, rsprp, NULL, NULL, NULL, NULL, NULL, NULL, NULL ); */
                        errval = rsb__pr_save(of, rsprp, rsprp->rsprap->filenamea, rsprp->rsprap->ca, rsprp->rsprap->incXa, rsprp->rsprap->incYa, rsprp->rsprap->nrhsa, rsprp->rsprap->typecodes, rsprp->rsprap->ta, RSB_BOOL_FALSE);
                        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                }
        }

        if(ds == RSB_PRD_CMP_DIV || ds == RSB_PRD_CMP_DIFF)
        for(argi=1;argi<argc;++argi)
        {
                struct rsb_rspr_t * rspr0p = NULL;
                RSB_PRL_SEP("\n");
		if( ds == RSB_PRD_CMP_DIV )
                	RSB_PRL("Will compare performance records of file %d/%d: %s to that of %s (first divided by second). Warning: assuming ALL parameters are conformant\n",argi+1,argc,argv[argi],argv[0]);
		if( ds == RSB_PRD_CMP_DIFF )
                	RSB_PRL("Will compare performance records of file %d/%d: %s to that of %s (first minus second). Warning: assuming ALL parameters are conformant\n",argi+1,argc,argv[argi],argv[0]);
                errval = rsb__pr_load(argv[argi],&rspr0p);
                if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
		if( rsprp->csf != rspr0p->csf && RSB_PR_NOC(rsprp) == RSB_PR_NOC(rspr0p))
                	RSB_PRL("It seems like one of the two records is incomplete!\n");
                errval = rsb__pr_cmp(rsprp,rspr0p,ds);
                RSB__PR_FREE(rspr0p);
                if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                errval = rsb__pr_dump(rsprp, rsprp->rsprap->filenamea, rsprp->rsprap->ca, rsprp->rsprap->incXa, rsprp->rsprap->incYa, rsprp->rsprap->nrhsa, rsprp->rsprap->typecodes, NULL, NULL );
                RSB__PR_FREE(rsprp);
                if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                errval = rsb__pr_load(argv[0],&rsprp);
	        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

        if(ds == RSB_PRD_CMP_MDUMP)
        for(argi=0;argi<argc;++argi)
        {
                RSB_PRL_SEP("\n");
                RSB_PRL("Dumping performance records of file %d/%d: %s\n",argi+1,argc,argv[argi]);
                errval = rsb__pr_load(argv[argi],&rsprp);
                if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
                errval = rsb__pr_dump(rsprp, rsprp->rsprap->filenamea, rsprp->rsprap->ca, rsprp->rsprap->incXa, rsprp->rsprap->incYa, rsprp->rsprap->nrhsa, rsprp->rsprap->typecodes, NULL, argv[argi] );
                RSB__PR_FREE(rsprp);
	        if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
        }

        RSB_PRL_SEP("\n");
err:
        RSB__PR_FREE(rsprp);

	if((errval |= rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))!=RSB_ERR_NO_ERROR)
		;

        return errval;
}

static void rsb__pr_upd_time(struct rsb_rspr_t * rsprp)
{
#if RSB_PRD_WANT_TIMESTAMP
	rsprp->tend = rsb_time();
#endif /* RSB_PRD_WANT_TIMESTAMP */
}

/* performance samples recording / dumping facility for rsbench : begin */
rsb_err_t rsb__pr_init(void**rsprpv, const struct rsb_mtx_t *mtxAp, rsb_int_t filenamen, rsb_int_t cn, rsb_int_t incXn, rsb_int_t incYn, rsb_int_t nrhsn, rsb_int_t ntypecodes, rsb_int_t tn, struct rsb_mbw_et_t * mbetp)
{
	/* 
	 * initialize a performance record 
	 * */
	struct rsb_rspr_t * rsprp = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( ! rsprpv )
       	{
	       	errval = RSB_ERR_ENOMEM;
	        RSB_PERR_GOTO(err,RSB_ERRM_ES);
       	}
        errval = rsb__pr_alloc(&rsprp, NULL, filenamen, cn, incXn, incYn, nrhsn, ntypecodes, tn);
	*rsprpv = rsprp;

#if RSB_PRD_WANT_ENVIRONMENT
	errval = rsb__rspr_all_env(rsprp);
        if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif /* RSB_PRD_WANT_ENVIRONMENT */
#if RSB_PRD_WANT_TIMESTAMP
	rsprp->tbeg = rsb_time();
	rsb__pr_upd_time(rsprp);
#endif /* RSB_PRD_WANT_TIMESTAMP */

#if RSB_PRD_WANT_MBW
	if(mbetp && mbetp->sn)
	{
		// FIXME: this is a transfer
		rsprp->mbet = *mbetp;
        	RSB_BZERO_P(mbetp); // zero the source struct !
	}
#endif /* RSB_PRD_WANT_MBW */

	return RSB_ERR_NO_ERROR;
err:
	RSB__PR_FREE(rsprp);
	return errval;
}

static rsb_err_t rsb__pr_set_idx(void*rsprpv, size_t idx, const struct rsb_mtx_t *mtxAp, const struct rsb_mtx_t *at_mtxAp, rsb_trans_t transA, rsb_perf_t op_time_best, rsb_perf_t mkl_csr_op_time_best, rsb_perf_t at_op_time_best, rsb_perf_t at_mkl_csr_op_time_best, rsb_int_t at_cn, rsb_int_t at_mkl_csr_cn, rsb_time_t at_t, rsb_int_t at_eps, const struct rsb_ts_t*otposp, const struct rsb_ts_t*btposp, const struct rsb_ts_t*otpmsp, const struct rsb_ts_t*btpmsp)
{
	/* 
	 * set performance record information
	 * Note: This inner, idx-based version can be invoked by internal, index-agnostic functions.
	 * */
	struct rsb_rspr_t * rsprp = rsprpv;
	
	rsb_bool_t have_own = RSB_BOOL_FALSE;

	if( rsprp->psa[idx].at_nsubm && rsprp->psa[idx].nsubm  )
		have_own = RSB_BOOL_TRUE; /* only mkl missing */

	if(RSB_CONST_IMPOSSIBLY_BIG_TIME != mkl_csr_op_time_best)
		rsprp->psa[idx].mkl_csr_op_time = mkl_csr_op_time_best;
	if(RSB_CONST_IMPOSSIBLY_BIG_TIME != op_time_best)
		rsprp->psa[idx].op_time = op_time_best;
	if(RSB_CONST_IMPOSSIBLY_BIG_TIME != at_mkl_csr_op_time_best)
		rsprp->psa[idx].at_mkl_csr_op_time = at_mkl_csr_op_time_best;
	if(RSB_IS_VALID_THREAD_COUNT( at_mkl_csr_cn )  )
		rsprp->psa[idx].at_mkl_csr_cn = at_mkl_csr_cn;

       	RSB_ASSIGN_IF_SP(rsprp->psa[idx].btpms,btpmsp)
       	RSB_ASSIGN_IF_SP(rsprp->psa[idx].otpms,otpmsp)

	if(!have_own)
	{
		if(RSB_CONST_IMPOSSIBLY_BIG_TIME != at_op_time_best)
			rsprp->psa[idx].at_op_time = at_op_time_best;
		if(RSB_CONST_IMPOSSIBLY_BIG_TIME != at_t)
			rsprp->psa[idx].at_t = at_t;
		if(RSB_IS_VALID_THREAD_COUNT( at_cn) )
			rsprp->psa[idx].at_cn = at_cn;
		if( -1 != at_eps ) 
			rsprp->psa[idx].at_eps = at_eps;

		rsprp->psa[idx].transA = transA;

        	RSB_ASSIGN_IF_SP(rsprp->psa[idx].btpos,btposp)
        	RSB_ASSIGN_IF_SP(rsprp->psa[idx].otpos,otposp)
	}

	if(!have_own)
	if(mtxAp)
	{
		rsprp->psa[idx].cmflops = rsb__estimate_mflops_per_op_spmv_uaua(mtxAp);
		rsprp->psa[idx].flagsA = mtxAp->flags;
		rsprp->psa[idx].nsubm = mtxAp->all_leaf_matrices_n;
		rsprp->psa[idx].nrA = mtxAp->nr;
		rsprp->psa[idx].ncA = mtxAp->nc;
		rsprp->psa[idx].nnzA = mtxAp->nnz;
#if RSB_STORE_IDXSA
		rsprp->psa[idx].isa = mtxAp->idxsa;
#else
		rsprp->psa[idx].isa = rsb__get_index_storage_amount(mtxAp);
#endif
	}

	if( at_mtxAp == NULL )
		at_mtxAp = mtxAp; /* FIXME: this shall be handled better  */

	if(!have_own)
	if( at_mtxAp )
	{
#if RSB_STORE_IDXSA
		rsprp->psa[idx].at_isa = at_mtxAp->idxsa;
#else
		rsprp->psa[idx].at_isa = rsb__get_index_storage_amount(at_mtxAp);
#endif
		rsprp->psa[idx].at_nsubm = at_mtxAp->all_leaf_matrices_n;
	}

	if( 0 == rsprp->psa[idx].uc ) /* if first encounter of sample, we increment pointer */
		rsprp->csf ++;

	rsprp->psa[idx].uc ++ ;

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__pr_set(void*rsprpv, const struct rsb_mtx_t *mtxAp, const struct rsb_mtx_t *at_mtxAp, rsb_int_t filenamei, rsb_int_t ci, rsb_int_t incXi, rsb_int_t incYi, rsb_int_t nrhsi, rsb_int_t typecodesi, rsb_int_t ti, rsb_trans_t transA, rsb_perf_t op_time_best, rsb_perf_t mkl_csr_op_time_best, rsb_perf_t at_op_time_best, rsb_perf_t at_mkl_csr_op_time_best, rsb_int_t at_cn, rsb_int_t at_mkl_csr_cn, rsb_time_t at_t, rsb_int_t at_eps, const struct rsb_ts_t*otposp, const struct rsb_ts_t*btposp, const struct rsb_ts_t*otpmsp, const struct rsb_ts_t*btpmsp)
{
	/* 
	 * set performance record information
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_rspr_t * rsprp = rsprpv;
	const size_t idx = rsb__pr_idx(rsprpv, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti);
	const char rsb_prl_lcc = RSB_PRL_LCC_IE;
#if 1
	rsb__pr_upd_time(rsprp);
	RSB_PRL("updating sample at index %zd (%d^th of %d), %d^th touch for (%d,%d,%d,%d,%d,%d,%d).\n",idx+1,rsprp->csf,RSB_PR_NOC(rsprp),rsprp->psa[idx].uc,filenamei,ci,incXi,incYi,nrhsi,typecodesi,ti);
	errval = rsb__pr_set_idx(rsprpv, idx, mtxAp, at_mtxAp, transA, op_time_best, mkl_csr_op_time_best, at_op_time_best, at_mkl_csr_op_time_best, at_cn, at_mkl_csr_cn, at_t, at_eps,otposp,btposp,otpmsp,btpmsp);
	return errval;
#else
	rsb_bool_t have_own = RSB_BOOL_FALSE;

	if( rsprp->psa[idx].at_nsubm && rsprp->psa[idx].nsubm  )
		have_own = RSB_BOOL_TRUE; /* only mkl missing */

	if(RSB_CONST_IMPOSSIBLY_BIG_TIME != mkl_csr_op_time_best)
		rsprp->psa[idx].mkl_csr_op_time = mkl_csr_op_time_best;
	if(RSB_CONST_IMPOSSIBLY_BIG_TIME != op_time_best)
		rsprp->psa[idx].op_time = op_time_best;
	if(RSB_CONST_IMPOSSIBLY_BIG_TIME != at_mkl_csr_op_time_best)
		rsprp->psa[idx].at_mkl_csr_op_time = at_mkl_csr_op_time_best;
	if(RSB_IS_VALID_THREAD_COUNT( at_mkl_csr_cn )  )
		rsprp->psa[idx].at_mkl_csr_cn = at_mkl_csr_cn;

	if(!have_own)
	{
		if(RSB_CONST_IMPOSSIBLY_BIG_TIME != at_op_time_best)
			rsprp->psa[idx].at_op_time = at_op_time_best;
		if(RSB_CONST_IMPOSSIBLY_BIG_TIME != at_t)
			rsprp->psa[idx].at_t = at_t;
		if(RSB_IS_VALID_THREAD_COUNT( at_cn) )
			rsprp->psa[idx].at_cn = at_cn;
		if( -1 != at_eps ) 
			rsprp->psa[idx].at_eps = at_eps;

		rsprp->psa[idx].transA = transA;
	}

	if(!have_own)
	if(mtxAp)
	{
		rsprp->psa[idx].cmflops = rsb__estimate_mflops_per_op_spmv_uaua(mtxAp);
		rsprp->psa[idx].flagsA = mtxAp->flags;
		rsprp->psa[idx].nsubm = mtxAp->all_leaf_matrices_n;
		rsprp->psa[idx].nrA = mtxAp->nr;
		rsprp->psa[idx].ncA = mtxAp->nc;
		rsprp->psa[idx].nnzA = mtxAp->nnz;
		rsprp->psa[idx].isa = rsb__get_index_storage_amount(mtxAp);
	}

	if( at_mtxAp == NULL )
		at_mtxAp = mtxAp; /* FIXME: this shall be handled better  */

	if(!have_own)
	if( at_mtxAp )
	{
		rsprp->psa[idx].at_isa = rsb__get_index_storage_amount(at_mtxAp);
		rsprp->psa[idx].at_nsubm = at_mtxAp->all_leaf_matrices_n;
	}

	if( 0 == rsprp->psa[idx].uc ) /* if first encounter of sample, we increment pointer */
		rsprp->csf ++;

	rsprp->psa[idx].uc ++ ;
	RSB_PRL("updating sample at index %zd (%d^th of %d), %d^th touch for (%d,%d,%d,%d,%d,%d,%d).\n",idx+1,rsprp->csf,RSB_PR_NOC(rsprp),rsprp->psa[idx].uc,filenamei,ci,incXi,incYi,nrhsi,typecodesi,ti);

	return RSB_ERR_NO_ERROR;
#endif
}

#define RSB_SYMCHAR(FLAGS) ( (RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_SYMMETRIC)) ? 'S' : ( (RSB_DO_FLAG_HAS(FLAGS,RSB_FLAG_HERMITIAN)) ? 'H' : 'G') )
#define RSB_DIV_NOT_BY_ZERO(D,Q) D = ( (Q) ? (D) / (Q) : 0.0 )
#define RSB_UPD_TO_MAX(VAR,VAL) (VAR)=RSB_MAX((VAR),(VAL))
#define RSB_UPD_TO_MIN(VAR,VAL) (VAR)=RSB_MIN((VAR),(VAL))

#define RSB_UPD_AMM(WHAT,ACC,MIN,MAX) (ACC)+=(WHAT); RSB_UPD_TO_MIN((MIN),(WHAT)); RSB_UPD_TO_MAX((MAX),(WHAT));
/* #define RSB_CMP_THR 1.00  */ /* 0% */
#define RSB_CMP_THR 1.01 /* below this we ignore the difference */
#define RSB_APE_THR 1.05 /* approximately equal or small difference  */
#define RSB_RLD_THR 2.00 /* relevant difference threshold */
#define RSB_HUD_THR 10.00 /* huge difference threshold */
#define RSB_CMP_THR_EXP /* "nearly same" threshold */ rsb__getenv_real_t("RSB_CMP_THR",RSB_CMP_THR)
#define RSB_APE_THR_EXP /* "approximately close" threshold */ rsb__getenv_real_t("RSB_APE_THR",RSB_APE_THR)
#define RSB_RLD_THR_EXP /* "relevant difference" threshold */ rsb__getenv_real_t("RSB_RLD_THR",RSB_RLD_THR)
#define RSB_HUD_THR_EXP /* "huge difference" threshold */ rsb__getenv_real_t("RSB_HUD_THR",RSB_HUD_THR)

#define RSB_FSTR_THN_THR(T1,T2,CMPT) ( (T1) * (CMPT) < (T2) ) /* faster only according to a small threshold; e.g. CMPT = RSB_CMP_THR */
#define RSB_SLWR_THN_THR(T1,T2) ( (T1) * RSB_CMP_THR >=(T2) && (T1) != (T2)  ) /* slower only according to a small threshold */
#define RSB_APPROX_EQUAL(T1,T2,APPT) ( ( (T1) * (APPT) >= (T2) ) && ( (T2) * (APPT) >= (T1) ) ) /* e.g. APPT = RSB_APE_THR  */
#define RSB_BOTH_FINITE(T1,T2)  ( ((T1)!=RSB_TIME_ZERO) && ((T2)!=RSB_TIME_ZERO) )
#define RSB_BIG_DIFF(T1,T2,RLDT) ( ( (T1) * (RLDT) < (T2) ) || ( (T2) * (RLDT) < (T1) ) ) /* e.g. RLDT = RSB_RLD_THR */
#define RSB_HUGE_DIFF(T1,T2,HGDT) ( ( (T1) * (HGDT) < (T2) ) || ( (T2) * (HGDT) < (T1) ) ) /* e.g. HGDT = RSB_HUD_THR  */
#define RSB_MIN_FINITE(X,Y) ( ((X)==(Y)) ? (X) : ( RSB_MIN(RSB_MAX(X,RSB_TIME_ZERO),RSB_MAX(Y,RSB_TIME_ZERO))) )

void rsb__mtxfn_bncp(char* dst, const char*src, int lm)
{
	/* 
	 * Matrix file name base name copy.
	 * This shall NOT be a library function !
	 * */
	size_t sl = 0;
		
	if(!dst || !src)
		goto ret;

	if(lm==0)
		rsb__strcpy(dst,rsb__basename(src));
	else
	{
		/* LaTeX underscore sanitization */
		const char*sp = rsb__basename(src);
		char*dp = dst;
		while(*sp)
		{
			*dp = *sp;
			if(*dp=='_')
				dp[0] = '\\',
				dp[1] = '_',
				++dp;
			++sp;
			++dp;
		}
		*dp = *sp;
	}
	sl = rsb__strlen(dst);
	if( sl > 0 && dst[sl-1] == RSB_DIR_SEPARATOR )
		dst[sl-1] = RSB_NUL;	/* for the case we want to 'sanitize' a dir basename */
	sl = rsb__strlen(dst);
#if RSB_WANT_EXPERIMENTAL_BINARY_COO
	if( sl >11 && strcmp(dst+sl-11,".mtx.bin.gz") == 0 )
		dst[sl-11] = RSB_NUL;
	if( sl > 8 && strcmp(dst+sl-8,".mtx.bin") == 0 )
		dst[sl-8] = RSB_NUL;
#endif /* RSB_WANT_EXPERIMENTAL_BINARY_COO */
	if( sl > 7 && strcmp(dst+sl-7,".mtx.gz") == 0 )
		dst[sl-7] = RSB_NUL;
	if( sl > 4 && strcmp(dst+sl-4,".mtx") == 0 )
		dst[sl-4] = RSB_NUL;
ret:
	return;
}

static void rsb__stdout_magnitude(double qty)
{
	const char * ooms = " KMGTPEZY"; // order of magnitude string

	while( qty >= 1000.0 && ooms[1] ) // stops on NUL
		qty /= 1000.0, ooms++;

	RSB_STDOUT("%4.1lf", qty);
	if(*ooms!=' ')
		RSB_STDOUT("%c", *ooms);
}

static rsb_err_t rsb__pr_dump_sample(const void*rsprpv, rsb_char_t**RSB_RESTRICT filenamea, rsb_int_t*ca, const rsb_int_t*incXa, const rsb_int_t*incYa, const rsb_int_t*nrhsa, const rsb_type_t*typecodes, const rsb_int_t*ta, const int*filenameifp, const int*ifilenameifp, const int*cifp , const int*incXifp , const int*incYifp , const int*nrhsifp , const int*typecodefip , const int*tifp, const rsb_trans_t*tfp, rsb_flags_t flagsA, rsb_flags_t nflagsA, rsb_int_t filenamei, rsb_int_t ci, rsb_int_t incXi, rsb_int_t incYi, rsb_int_t nrhsi, rsb_int_t typecodesi, rsb_int_t ti, int rds, int wltm)
{
        /* TODO: may print a different record if( rsprp->ror == RSB_BOOL_TRUE ) */

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const struct rsb_rspr_t * rsprp = rsprpv;
	{
		const size_t idx = rsb__pr_idx(rsprpv, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti);
		const struct rsb_rsps_t*psp = &(rsprp->psa[idx]);
		const rsb_int_t nnzA = psp->nnzA, nrA = psp->nrA, ncA = psp->ncA;
		rsb_char_t fnbuf[RSB_MAX_FILENAME_LENGTH];
		rsb_perf_t brsb_op_time = psp->op_time;
		rsb_perf_t bmkl_op_time = psp->mkl_csr_op_time;

		rsb__mtxfn_bncp(fnbuf,filenamea[filenamei],wltm);
		
		if( psp->at_op_time != RSB_TIME_ZERO )
			brsb_op_time = RSB_MIN( brsb_op_time, psp->at_op_time );

		if( psp->at_mkl_csr_op_time != RSB_TIME_ZERO )
			bmkl_op_time = RSB_MIN( bmkl_op_time, psp->at_mkl_csr_op_time );
#if 0
                /* TODO: will have to start here for the detailed stats: */ 
                RSB_STAT_DUMP_TS(psp->otpos); 
                RSB_STAT_DUMP_TS(psp->btpos); 
                RSB_STAT_DUMP_TS(psp->otpms); 
                RSB_STAT_DUMP_TS(psp->btpms);
#endif
        	if(rds==RSB_PRD_STYLE_TBL)
		{
			const double rldt = RSB_RLD_THR_EXP/*, cmpt = RSB_CMP_THR_EXP, appt = RSB_APE_THR_EXP, hgdt = RSB_HUD_THR_EXP*/;
			const char * ss = RSB_PRL_FSEPSTR_IE; /* separator string */
			const char * ts = RSB_PRL_ENDLSTR_IE; /* terminator string */
			const char * rcc = "\\cellcolor{red}"; /* red colored cell */
			const char * gcc = "\\cellcolor{PaleGreen1}"; /* requires x11names */
			const char * Gcc = "\\cellcolor{green}"; /* green (greener!) colored cell */
			//const char * gcc = "\\cellcolor{green}"; /* green colored cell */
			const char * bcc = "\\cellcolor{blue}"; /* blue colored cell */
			//const char * ycc = "\\cellcolor{yellow}"; /* yellow colored cell */
			const char * ycc = "\\cellcolor{LightGoldenrod1}"; /* requires x11names */
			const char * pcc = "\\cellcolor{pink}"; /* pink colored cell */
			const char * ncc = ""; /* no color cell  (no LaTeX Markup) */
			const char * bfs = "\\bfseries "; /* bold font series */
			const char * rlns = psp->mkl_csr_op_time ? RSB_ON_IF_LEM( psp->op_time*rldt ,psp->mkl_csr_op_time, bfs , bfs, ncc) : ""; /* relevant non-autotuned speedup over mkl */
			const char * rlas = psp->at_mkl_csr_op_time ? RSB_ON_IF_LEM( psp->at_op_time*rldt ,psp->at_mkl_csr_op_time,bfs , bfs, ncc) : ""; /* relevant autotuned speedup over mkl */

			if(wltm < 2)
				rcc = gcc = Gcc = bcc = ycc = pcc = bfs = rlns = rlas = ncc; /* nullify LaTeX markup */

			RSB_STDOUT("%s", RSB_ON_IF_LEM(psp->at_op_time, psp->at_mkl_csr_op_time , gcc, ncc, rcc));/* EXPERIMENTAL */
			//RSB_STDOUT("%s", RSB_ON_IF_LEM(psp->at_op_time, psp->at_mkl_csr_op_time , RSB_ON_IF_LEM( psp->at_op_time*rldt ,psp->at_mkl_csr_op_time,Gcc , gcc, gcc), ncc, rcc));/* EXPERIMENTAL */
			RSB_STDOUT("%s",ss);
			RSB_STDOUT("%s%s%d%s%d%s",
				fnbuf,ss,
				nrA,ss,
				ncA,ss
			);
			if( rsb__util_atoi(rsb__getenv("RSB_PR_WLTC")) > 0 )
				{ rsb__stdout_magnitude((double) nnzA); RSB_STDOUT("%s",ss); }
			else
				RSB_STDOUT("%d%s", nnzA,ss);
			if(rsprp->incXn > 1 && rsprp->incYn > 1)
			RSB_STDOUT("%d%s%d%s",
				incXa[incXi],ss,
				incYa[incYi],ss
			);
			RSB_STDOUT("%d%s%c%s%c%s%c%s",
			       	nrhsa[nrhsi],ss,
			       	typecodes[typecodesi],ss,
				RSB_SYMCHAR(psp->flagsA),ss,
				RSB_TRANSPOSITION_AS_CHAR(psp->transA),ss
			);
			RSB_STDOUT(
				"%2d%s%s%2d%s%s%2d%s",
				ca[ci],ss,
				ca[ci] == psp->at_cn ? ncc : (ca[ci] / 2 >= psp->at_cn ? rcc : ycc),/* EXPERIMENTAL */
				psp->at_cn,ss,
				ca[ci] == psp->at_mkl_csr_cn ? ncc : (ca[ci] / 2 >= psp->at_mkl_csr_cn ? rcc : ycc),/* EXPERIMENTAL */
				psp->at_mkl_csr_cn,ss
				);
			RSB_STDOUT(
				"%.4lf%s%s%.4lf%s",
				((rsb_perf_t)(psp->isa))/nnzA,ss,
				RSB_ON_IF_LEM(psp->at_isa , psp->isa , gcc, ncc, pcc),/* EXPERIMENTAL */
				((rsb_perf_t)(psp->at_isa))/nnzA,ss
				);
			RSB_STDOUT(
				"%d%s%s%d%s",
				psp->nsubm,ss,
				RSB_ON_IF_LEM(psp->at_nsubm , psp->nsubm , gcc, ncc, pcc),/* EXPERIMENTAL */
				psp->at_nsubm,ss
				);
			RSB_STDOUT(
				"%.2lf%s%2.3le%s%s%2.3le%s",
                                RSB_XFLOPS(brsb_op_time,nrhsa[nrhsi],psp->cmflops),ss,
				psp->op_time,ss,
				rlns,
				psp->mkl_csr_op_time,ss
				);
			RSB_STDOUT(
				"%s%2.3le%s%s%s%2.3le%s%2.3le%s",
				RSB_ON_IF_LEM(psp->at_op_time, psp->op_time, gcc, ncc, pcc),/* EXPERIMENTAL */
				psp->at_op_time,ss,
				RSB_ON_IF_LEM(psp->at_mkl_csr_op_time, psp->mkl_csr_op_time, gcc, ncc, pcc),/* EXPERIMENTAL */
				rlas,
				psp->at_mkl_csr_op_time,ss,
				psp->at_t,ss
				);
#if RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH
{
			const size_t so = RSB_SIZEOF_BACKUP(typecodes[typecodesi]); /* size of */
			const size_t mo = (so*psp->nnzA+((rsb_perf_t)(psp->at_isa))); /* matrix occupation */
			const size_t oo = so*nrhsa[nrhsi]*(nrA+ncA); /* operands occupation */
			const size_t owt = so*nrhsa[nrhsi]*RSB_ELSE_IF_TRANSPOSE(psp->nrA,psp->ncA,psp->transA); /* operands write traffic */
			/* size_t mrt = oo + mo; */ /* minimal read traffic */
			const size_t mwt = oo + mo + owt; /* minimal read+write traffic */
			/* rsb_perf_t mrb = ((rsb_perf_t)mrt)/(psp->at_op_time*1e9);*/ /* minimal read bandwidth, GBps */
			const rsb_perf_t mwb = ((rsb_perf_t)mwt)/(psp->at_op_time*1e9); /* minimal read+write bandwidth, GBps */
			/* RSB_STDOUT( "%3.2le%s", mrb, ss); */ /* BW/RDminBWIDTH: minimal bandwidth, GB/s */
			RSB_STDOUT( "%3.2le%s", mwb, ss); /* BW/RWminBW: minimal bandwidth, GB/s */
}
{
			const size_t so = RSB_SIZEOF_BACKUP(typecodes[typecodesi]); /* size of */
			const size_t mo = (so*psp->nnzA+((rsb_perf_t)(psp->at_isa))); /* matrix occupation */
			const size_t oo = so*nrhsa[nrhsi]*(nrA+ncA); /* operands occupation */
			const size_t mt = oo + mo; /* minimal traffic */
			const rsb_perf_t om = 1e6 *(psp->cmflops * nrhsa[nrhsi]); /* operation flops */
			const rsb_perf_t bm = (1.0 / om) * mt; /* bytes per flops */
			RSB_STDOUT( "%3.2le%s", bm, ss); /* CB: code balance, bpf = bytes per flop */
}
#endif /* RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH */
			RSB_STDOUT(
				"%s%d%s%3.2le%s\n",
				psp->at_eps == 0 ? rcc : ncc,/* EXPERIMENTAL */
				psp->at_eps,ss,
				psp->cmflops*nrhsa[nrhsi],
				ts
			);
			if(wltm > 1)
				RSB_STDOUT("%%...\n");

		}

        	if(rds==RSB_PRD_STYLE_CMP)
                {
		RSB_STDOUT("%d:%s %d %d %d %d %d %d %c %c %c",
				filenamei,fnbuf, 
				nrA,ncA,nnzA,
				incXa[incXi], 
				incYa[incYi],
			       	nrhsa[nrhsi],
			       	typecodes[typecodesi],
				RSB_SYMCHAR(psp->flagsA),
				RSB_TRANSPOSITION_AS_CHAR(psp->transA)
			);
		RSB_STDOUT(" %.2lf", bmkl_op_time/ brsb_op_time);
		RSB_STDOUT(" %.2lf %.2lf %.2lf %.2lf",
                                RSB_XFLOPS(brsb_op_time,nrhsa[nrhsi],psp->cmflops),
                                (psp->        op_time/brsb_op_time),
                                RSB_XFLOPS(bmkl_op_time,nrhsa[nrhsi],psp->cmflops),
                                (psp->mkl_csr_op_time/bmkl_op_time)
			);
		RSB_STDOUT("\n");
                }
	}
	return errval;
}

static rsb_err_t rsb__pr_filter(const struct rsb_rsps_t*psp, const rsb_int_t*ta, const int*filenameifp, const int*ifilenameifp, const int*cifp , const int*incXifp , const int*incYifp , const int*nrhsifp , const int*typecodefip , const int*tifp, const rsb_trans_t*tfp, rsb_flags_t flagsA, rsb_flags_t nflagsA,
	       	rsb_int_t filenamei, rsb_int_t ci, rsb_int_t incXi, rsb_int_t incYi, rsb_int_t nrhsi, rsb_int_t typecodesi, rsb_int_t ti)
{
		if( filenameifp && ( *filenameifp != filenamei  ) )
			goto skipit;
		if(ifilenameifp && (*ifilenameifp <  filenamei  ) )
			goto skipit;
		if( cifp        && (        *cifp != ci         ) )
			goto skipit;
		if( incXifp     && (     *incXifp != incXi      ) )
			goto skipit;
		if( incYifp     && (     *incYifp != incYi      ) )
			goto skipit;
		if( nrhsifp     && (     *nrhsifp !=    nrhsi   ) )
			goto skipit;
		if( typecodefip && ( *typecodefip != typecodesi ) )
			goto skipit;
		if( ( flagsA != RSB_FLAG_NOFLAGS) && ! (psp->flagsA &  flagsA) )
			goto skipit;
		if( (nflagsA != RSB_FLAG_NOFLAGS) &&   (psp->flagsA & nflagsA) )
			goto skipit;
		if(ta)
		{
			if( tifp        && (     *tifp    !=    ti      ) )
				goto skipit;
		}
		else
		{
			if( tfp         && (      *tfp    !=    psp->transA ) )
				goto skipit;
		}

	if( psp->uc > 0 && psp->uc < 3 )
		;
	else
	{
		goto skipit;/* we skip this iteration's sample */
	}

		return RSB_ERR_NO_ERROR;
skipit:
		return RSB_ERR_GENERIC_ERROR;
}

rsb_err_t rsb__pr_dump_inner(const void*rsprpv, rsb_char_t**RSB_RESTRICT filenamea, rsb_int_t*ca, const rsb_int_t*incXa, const rsb_int_t*incYa, const rsb_int_t*nrhsa, const rsb_type_t*typecodes, const rsb_int_t*ta, const int*filenameifp, const int*ifilenameifp, const int*cifp , const int*incXifp , const int*incYifp , const int*nrhsifp , const int*typecodefip , const int*tifp, const rsb_trans_t*tfp, rsb_flags_t flagsA, rsb_flags_t nflagsA, rsb_char_t *ltag, const rsb_char_t *fprfn)
{
	/*
	 * dump a performance record, inner
         * May introduce strict checks or verbosity options.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const struct rsb_rspr_t * rsprp = rsprpv;
	rsb_int_t filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti;
	rsb_int_t nocsa = 0; /* number of considered samples  */
	rsb_int_t nocma = 0, lacma = -1; /* number of considered matrices, last considered matrix */
	rsb_int_t noats = 0, noatf = 0; /* number of auto tuning successes / failures */
	rsb_int_t ntufm = 0, ntusm = 0; /* number (of) times untuned rsb (was) faster/slower (than) mkl */
	rsb_int_t nttfm = 0, nttsm = 0; /* number (of) times  tuned  rsb (was) faster/slower (than   tuned) mkl */
	rsb_int_t nttfu = 0, nttsu = 0; /* number (of) times  tuned  rsb (was) faster/slower (than untuned) mkl */
	rsb_int_t nnzrm = 0; /* number (of) samples tuned rsb and tuned mkl were not zero  */
	rsb_int_t ntmfm = 0, ntmsm = 0; /* number (of times)  tuned mkl (was) faster/slower (than untuned) mkl */
	rsb_int_t ntasl = 0, ntasm = 0, ntase = 0; /* number (of times)  autotuned (was) subdivided less/more/equally */
	rsb_int_t ntatl = 0, ntatm = 0, ntate = 0; /* number (of times)  autotuned used  less/more/equal threads */
	rsb_int_t ntsrf = 0; /* number of times multi-nrhs rsb was faster (than rsb unstrided)  */
	rsb_int_t ntsmf = 0; /* number of times multi-nrhs mkl was faster (than mkl unstrided)  */
	rsb_int_t vscm = 0; /* valid samples containing mkl */
	double aoatsp = 0.0; /* average of auto tuning speedup (percentage) */
	double aoatsr = 0.0; /* average of auto tuning speedup (ratio) */
	double aoatac = 0.0, miatac = RSB_CONST_IMPOSSIBLY_BIG_TIME, maatac = 0.0; /* average,min,max of auto tuning amortization cost */
	double aoatam = 0.0, miatam = RSB_CONST_IMPOSSIBLY_BIG_TIME, maatam = 0.0; /* average,min,max of auto tuning amortization to (tuned) mkl (cost) */
	double aoatau = 0.0, miatau = RSB_CONST_IMPOSSIBLY_BIG_TIME, maatau = 0.0; /* average,min,max of auto tuning amortization to (untuned) mkl (cost) */
	double aoatuo = 0.0, miatuo = RSB_CONST_IMPOSSIBLY_BIG_TIME, maatuo = 0.0, toatuo = 0.0; /* average,min,max,total of auto tuning untuned ops */
	double aoatto = 0.0, miatto = RSB_CONST_IMPOSSIBLY_BIG_TIME, maatto = 0.0, toatto = 0.0; /* average,min,max,total of auto tuning  tuned  ops */
	double aouatc = 0.0, miuatc = RSB_CONST_IMPOSSIBLY_BIG_TIME, mauatc = 0.0, touatc = 0.0; /* average,min,max,total of unsuccessful auto tuning cost  */
	double aotrtm = 0.0; /* average of  tuned  ratio   with respect to mkl  */
	double aotstm = 0.0; /* average of  tuned  speedup with respect to mkl  */
	double aotstu = 0.0; /* average of  tuned  speedup with respect to (untuned) mkl  */
	double aoustm = 0.0; /* average of untuned speedup with respect to mkl  */
	double aotssm = 0.0; /* average of  tuned  speedup (of mkl to rsb, when) slower (than) mkl */
	double aoussm = 0.0; /* average of untuned speedup (of mkl to rsb, when) slower (than) mkl */
	double aotsmm = 0.0; /* average of  tuned  speedup (of mkl to ) to mkl (always) */
	double msurwm = 0.0; /* maximal speedup (of) untuned rsb w.r.t. mkl  */
	double mstrwm = 0.0; /* maximal speedup (of)   tuned rsb w.r.t. mkl  */
	double mstrwu = 0.0; /* maximal speedup (of)   tuned rsb w.r.t. untuned mkl  */
	double mstmwm = 0.0; /* maximal speedup (of)   tuned mkl w.r.t. mkl  */
	double msumwr = 0.0; /* maximal speedup (of) untuned mkl w.r.t. rsb  */
	double mstmwr = 0.0; /* maximal speedup (of)   tuned mkl w.r.t. rsb  */
	double mstrwr = 0.0; /* maximal speedup (of)   tuned rsb w.r.t. rsb  */
	double aoratt = 0.0, miratt = RSB_CONST_IMPOSSIBLY_BIG_TIME, maratt = 0.0, toratt = 0.0; /* average,min,max,total of rsb auto tuning time */
	double aosatt = 0.0, misatt = RSB_CONST_IMPOSSIBLY_BIG_TIME, masatt = 0.0, tosatt = 0.0; /* average,min,max,total of rsb auto tuning time (when   successful) */
	double aouatt = 0.0, miuatt = RSB_CONST_IMPOSSIBLY_BIG_TIME, mauatt = 0.0, touatt = 0.0; /* average,min,max,total of rsb auto tuning time (when unsuccessful) */
	/* double aomatt = 0.0, mimatt = RSB_CONST_IMPOSSIBLY_BIG_TIME, mamatt = 0.0, tomatt = 0.0;*/ /* average,min,max,total of mkl auto tuning time */
	double avrmps = 0.0, mirmps = RSB_CONST_IMPOSSIBLY_BIG_TIME, marmps = 0.0; /* average,min,max rsb (canonical) mflops per second (tuned) */
	double avRmps = 0.0, miRmps = RSB_CONST_IMPOSSIBLY_BIG_TIME, maRmps = 0.0; /* average,min,max rsb (canonical) mflops per second (untuned) */
	double avmmps = 0.0, mimmps = RSB_CONST_IMPOSSIBLY_BIG_TIME, mammps = 0.0; /* average,min,max mkl (canonical) mflops per second (tuned) */
	double avMmps = 0.0, miMmps = RSB_CONST_IMPOSSIBLY_BIG_TIME, maMmps = 0.0; /* average,min,max mkl (canonical) mflops per second (untuned) */
	double aorott = 0.0, mirott = RSB_CONST_IMPOSSIBLY_BIG_TIME, marott = 0.0, torott = 0.0; /* average,min,max,total of rsb auto operation time (tuned) */
	double aoRott = 0.0, miRott = RSB_CONST_IMPOSSIBLY_BIG_TIME, maRott = 0.0, toRott = 0.0; /* average,min,max,total of rsb auto operation time (untuned) */
	double aomott = 0.0, mimott = RSB_CONST_IMPOSSIBLY_BIG_TIME, mamott = 0.0, tomott = 0.0; /* average,min,max,total of mkl auto operation time (tuned) */
	double aoMott = 0.0, miMott = RSB_CONST_IMPOSSIBLY_BIG_TIME, maMott = 0.0, toMott = 0.0; /* average,min,max,total of mkl auto operation time (untuned) */
	double avrsmv = 0.0, mirsmv = RSB_CONST_IMPOSSIBLY_BIG_TIME, marsmv = 0.0; /* average,min,max of rsb speedup for multi-vector */
	double avmsmv = 0.0, mimsmv = RSB_CONST_IMPOSSIBLY_BIG_TIME, mamsmv = 0.0; /* average,min,max of mkl speedup for multi-vector */
	double avnzbt = 0.0, minzbt = RSB_CONST_IMPOSSIBLY_BIG_TIME, manzbt = 0.0; /* average,min,max of nnz p.s. in untuned before tuning */
	double avnzat = 0.0, minzat = RSB_CONST_IMPOSSIBLY_BIG_TIME, manzat = 0.0; /* average,min,max of nnz p.s. in successfully tuned */
	double avbybt = 0.0, mibybt = RSB_CONST_IMPOSSIBLY_BIG_TIME, mabybt = 0.0; /* average,min,max of bytes p.s. in untuned before tuning */
	double avbyat = 0.0, mibyat = RSB_CONST_IMPOSSIBLY_BIG_TIME, mabyat = 0.0; /* average,min,max of bytes p.s. in successfully tuned */
	double avbpna = 0.0, mibpna = RSB_CONST_IMPOSSIBLY_BIG_TIME, mabpna = 0.0; /* average,min,max of bytes p.nnz in untuned before tuning */
	double avbpnb = 0.0, mibpnb = RSB_CONST_IMPOSSIBLY_BIG_TIME, mabpnb = 0.0; /* average,min,max of bytes p.nnz in successfully tuned */
#if RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH
	double avlorb = 0.0, milorb = RSB_CONST_IMPOSSIBLY_BIG_TIME, malorb = 0.0; /* average,min,max (autotuned) liminal/minimal operands reading bandwidth */
	double avlowb = 0.0, milowb = RSB_CONST_IMPOSSIBLY_BIG_TIME, malowb = 0.0; /* average,min,max (autotuned) liminal/minimal operands   r/w   bandwidth */
	double avcoba = 0.0, micoba = RSB_CONST_IMPOSSIBLY_BIG_TIME, macoba = 0.0; /* average,min,max code balance */
#endif /* RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH */
	rsb_int_t noc = 0;
        rsb_int phase = 0; /* first count samples, then dump */
        double cmpt = RSB_CMP_THR_EXP, appt = RSB_APE_THR_EXP, rldt = RSB_RLD_THR_EXP, hgdt = RSB_HUD_THR_EXP;
        const rsb_int_t wdbg = 0; /* want debug */
        int rds = rsb__getenv_real_t("RSB_PR_SR",RSB_PRD_STYLE_TBL);
        int wltm = rsb__getenv_real_t("RSB_PR_WLTC",0); /* Want LaTeX tables mode */
	const char * ss = RSB_PRL_FSEPSTR_IE; /* separator string */
	const char * ts = RSB_PRL_ENDLSTR_IE; /* terminator string */
	const char rsb_prl_lcc = RSB_PRL_LCC_IE;
	const char*rsb_prl_tcs = RSB_PRL_TCS_IE;
#if RSB_PRD_WANT_MBW
	const rsb_bool_t wpptmbw = rsprp->mbet.sn > 0 ? RSB_BOOL_TRUE : RSB_BOOL_FALSE; /* want performance proportion to memory bandwidth  */
#else
	const rsb_bool_t wpptmbw = RSB_BOOL_FALSE;
#endif /* RSB_PRD_WANT_MBW */

	if(0 /* FIXME: move this notice to an outer function */)
        if(rds!=RSB_PRD_STYLE_TBL)
	{
        	RSB_PRL("Further environment variables:\n");
        	RSB_PRL("RSB_CMP_THR\n");
        	RSB_PRL("RSB_APE_THR\n");
        	RSB_PRL("RSB_RLD_THR\n");
        	RSB_PRL("RSB_HUD_THR\n");
	}

	RSB_DEBUG_ASSERT(rsprp);

	noc = RSB_PR_NOC(rsprp);
gop2:
        if(phase == 1)
	{
		if(rds==RSB_PRD_STYLE_TBL && wltm > 0)
		{
			rsb_char_t fnbuf[RSB_MAX_FILENAME_LENGTH];
			rsb__mtxfn_bncp(fnbuf,rsb__basename(ltag),1);
			if(ltag)
				RSB_PRT("\\section{Record: %s}\n",fnbuf);
		}
        	RSB_PRL("Dump from a base of %d samples (of max %d) ordered by ",rsprp->csf,noc);
        	RSB_STDOUT("(%d,%d,%d,%d,%d,%d,%d) = (%s).\n",
        		rsprp->filenamen, rsprp->cn, rsprp->incXn, rsprp->incYn,
        	       	rsprp->nrhsn, rsprp->ntypecodes, rsprp->tn,
        		"filename x cores x incX x incY x nrhs x typecode x transA");
        if(rds==RSB_PRD_STYLE_TBL)
        {
		if(wltm > 0)
			//RSB_PRT("\\begin{table}[ht]\\begin{center}\\begin{tabular}{r*{26}{r}r}\\hline\n");
			RSB_PRT("\\begin{longtabu}{r*{26}{r}r}\\hline\n");
        	if(!ifilenameifp) /* no printout of records in this mode */
        	{
			// RSB_PRL("Each sample:\n"); 
			RSB_PRT("BESTCODE%sMTX%sNR%sNC%sNNZ%s", ss,ss,ss,ss,ss);
			if(rsprp->incXn > 1 && rsprp->incYn > 1)
				RSB_PRC("INCX%sINCY%s",ss,ss);
			RSB_PRC("NRHS%sTYPE%sSYM%sTRANS%sNT%sAT-NT%sAT-" RSB_MKL_S "-NT%sBPNZ%sAT-BPNZ%sNSUBM%sAT-SUBM%sRSBBEST-MFLOPS%sOPTIME%s" RSB_MKL_S "-OPTIME%sAT-OPTIME%sAT-" RSB_MKL_S "-OPTIME%sAT-TIME%s""RWminBW-GBps%s""CB-bpf%sAT-MS%sCMFLOPS%s\n",ss,ss,ss,ss,ss,ss,ss,ss,ss, ss,ss,ss,ss,ss,ss,ss,ss,ss,ss, ss,ts); /* FIXME: RWminBW and CB depend on RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH=0 */
        	}
		if(wltm > 0)
			RSB_PRT("\\hline\n");
        }
        if(rds==RSB_PRD_STYLE_CMP)
        {
        	if(!ifilenameifp) /* no printout of records in this mode */
        	{
        		RSB_PRL("Each sample: BESTCODE MTX NR NC NNZ INCX INCY NRHS TYPE SYM TRANS " RSB_MKL_S "_OP_T/RSB_OP_T RSB_OP_T RSB_MFLOPS " RSB_MKL_S "_OP_T " RSB_MKL_S "_MFLOPS\n");
        	}
        }
	}

	for(     filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
	for(ci=0;ci<rsprp->cn;++ci)
	for(     incXi=0;     incXi<rsprp->incXn     ;++incXi     )
	for(     incYi=0;     incYi<rsprp->incYn     ;++incYi     )
	for(     nrhsi=0;     nrhsi<rsprp->nrhsn     ;++nrhsi     )
	for(typecodesi=0;typecodesi<rsprp->ntypecodes;++typecodesi)
	for(ti=0;ti<rsprp->tn;++ti)
	{
		const size_t idx = rsb__pr_idx(rsprpv, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti);
		const struct rsb_rsps_t*psp = &(rsprp->psa[idx]);
		/* rsb_int_t nnzA = psp->nnzA, nrA = psp->nrA, ncA = psp->ncA; */
		rsb_bool_t atweost = ( ( psp->nsubm != psp->at_nsubm ) && ( psp->at_nsubm != 0 ) ); /* autotuning was effective on structure */
		rsb_bool_t atweoth = ( ( ca[ci] != psp->at_cn ) && ( psp->at_cn != 0 ) ); /* autotuning was effective on threads   */
		rsb_bool_t atweoti = ( ( psp->op_time > psp->at_op_time ) && ( psp->at_op_time != RSB_TIME_ZERO ) ); /* autotuning was effective on time (this is here for testing purposes) */
		rsb_bool_t atwe = RSB_BOOL_OR/*RSB_BOOL_AND*/(atweoti,RSB_BOOL_OR(atweost,atweoth));  /* autotuning was effective; 'or' condition is for testing purposes */
#if 1
		if( RSB_ERR_NO_ERROR != rsb__pr_filter(psp, ta, filenameifp, ifilenameifp, cifp , incXifp , incYifp , nrhsifp , typecodefip , tifp, tfp, flagsA, nflagsA, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti) )
			continue;
#else
		/* obsolete code */
		if( filenameifp && ( *filenameifp != filenamei  ) )
			continue;
		if(ifilenameifp && (*ifilenameifp <  filenamei  ) )
			continue;
		if( cifp        && (        *cifp != ci         ) )
			continue;
		if( incXifp     && (     *incXifp != incXi      ) )
			continue;
		if( incYifp     && (     *incYifp != incYi      ) )
			continue;
		if( nrhsifp     && (     *nrhsifp !=    nrhsi   ) )
			continue;
		if( typecodefip && ( *typecodefip != typecodesi ) )
			continue;
		if( ( flagsA != RSB_FLAG_NOFLAGS) && ! (psp->flagsA &  flagsA) )
			continue;
		if( (nflagsA != RSB_FLAG_NOFLAGS) &&   (psp->flagsA & nflagsA) )
			continue;
		if(ta)
		{
			if( tifp        && (     *tifp    !=    ti      ) )
				continue;
		}
		else
		{
			if( tfp         && (      *tfp    !=    psp->transA ) )
				continue;
		}
#endif
		
	if( psp->uc > 2 )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_ERROR("Updates count (%d) illegal --- max is 2 (internal error)!\n",psp->uc);
	       	goto ret;
	}

	if( psp->uc > 0 && psp->uc < 3 )
	{
		if( lacma != filenamei )
			lacma = filenamei,
			nocma++;
		nocsa ++; /* we consider this iteration's sample */
	}
	else
		continue; /* we skip this iteration's sample */

        if( phase == 0)
                continue; /* we only want to evaluate nocsa */

	if( psp->mkl_csr_op_time )
		vscm++;

	if(!ifilenameifp) /* no printout of records in this mode */
	if( psp->uc > 0 && psp->uc < 3 ) /* TODO: now useless loop; merge this ... */
	{
		rsb_char_t fnbuf[RSB_MAX_FILENAME_LENGTH];
		rsb_char_t usdc,asdc,csdc = '_'; /* untuned/autotuned/comparison  score descriptor char  */
		rsb_time_t best_rsb_op_time = RSB_MIN_FINITE(psp->op_time,        psp->at_op_time        );
		rsb_time_t best_mkl_op_time = RSB_MIN_FINITE(psp->mkl_csr_op_time,psp->at_mkl_csr_op_time);

		usdc = ( psp->op_time    <    psp->mkl_csr_op_time || psp->mkl_csr_op_time == RSB_TIME_ZERO )?(psp->op_time > psp->at_op_time ? 'R' : 'r' ) : ( psp->mkl_csr_op_time > psp->at_mkl_csr_op_time ? 'M' : 'm' );
		asdc = ( psp->at_op_time < psp->at_mkl_csr_op_time || psp->at_mkl_csr_op_time == RSB_TIME_ZERO )?(psp->op_time > psp->at_op_time ? 'R' : 'r' ) : ( psp->mkl_csr_op_time > psp->at_mkl_csr_op_time ? 'M' : 'm' );
	
		if( RSB_APPROX_EQUAL(best_rsb_op_time, best_mkl_op_time, appt ) )
			csdc = '~';
		if( RSB_BOTH_FINITE( best_rsb_op_time, best_mkl_op_time) && RSB_BIG_DIFF( best_rsb_op_time, best_mkl_op_time, rldt ) )
			csdc = '.';
		if( RSB_BOTH_FINITE( best_rsb_op_time, best_mkl_op_time) && RSB_HUGE_DIFF( best_rsb_op_time, best_mkl_op_time, hgdt ) )
			csdc = '!';
/*
		if( RSB_BOTH_FINITE( best_rsb_op_time, best_mkl_op_time) && csdc == ' ')
			csdc = '*';
*/
		rsb__mtxfn_bncp(fnbuf,filenamea[filenamei],0);
        	if(rds==RSB_PRD_STYLE_TBL) /* TODO: not for plots */
		RSB_PRT("%4zd:%c%s%c%c ",idx+1, usdc, ((wltm && csdc=='_')?"\\":""), csdc, asdc);
		rsb__pr_dump_sample(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, filenameifp, ifilenameifp, cifp, incXifp , incYifp , nrhsifp , typecodefip , tifp, tfp, flagsA, nflagsA, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti, rds, wltm);
	}

		if(     /* nrhsi==0 && */ rsprp->nrhsn > 1 )
		{
			const rsb_coo_idx_t miv = rsb__util_find_min_index(nrhsa, rsprp->nrhsn);
			const rsb_coo_idx_t mav = rsb__util_find_max_index(nrhsa, rsprp->nrhsn);

			RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(mav));
			RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(miv));

			if( nrhsa[miv] < nrhsa[mav] && nrhsi != miv ) /* if min and max differ */
			{
                                rsb_coo_idx_t nrhsi0 = /*miv*/ miv;
                                rsb_coo_idx_t nrhsi1 = /*mav*/ nrhsi;
                                const size_t idx0 = rsb__pr_idx(rsprpv, filenamei, ci, incXi, incYi, nrhsi0, typecodesi, ti);
                                const size_t idx1 = rsb__pr_idx(rsprpv, filenamei, ci, incXi, incYi, nrhsi1, typecodesi, ti);
                                /* ratio of performances (canonical flops per time) */
                                RSB_UPD_AMM( (rsprp->psa[idx0].at_op_time*nrhsa[nrhsi1] )     / (         rsprp->psa[idx1].at_op_time*nrhsa[nrhsi0] ) ,avrsmv, mirsmv , marsmv);
                                RSB_UPD_AMM( (rsprp->psa[idx0].mkl_csr_op_time*nrhsa[nrhsi1]) / ( rsprp->psa[idx1].at_mkl_csr_op_time*nrhsa[nrhsi0] ) ,avmsmv, mimsmv , mamsmv);
				ntsrf++;
				ntsmf++;
	       		}
	       	}

		RSB_UPD_AMM(psp->at_t,                                               aoratt, miratt, maratt );
		/* RSB_UPD_AMM(psp->at_mkl_csr_t,                                    aomatt, mimatt, mamatt ); */
		RSB_UPD_AMM(psp->at_op_time,                                         aorott, mirott, marott );
		RSB_UPD_AMM(psp->op_time,                                            aoRott, miRott, maRott );
		RSB_UPD_AMM(psp->at_mkl_csr_op_time,                                 aomott, mimott, mamott );
		RSB_UPD_AMM(psp->mkl_csr_op_time,                                    aoMott, miMott, maMott );
		RSB_UPD_AMM((psp->cmflops * nrhsa[nrhsi]) / psp->at_op_time,         avrmps, mirmps, marmps );
		RSB_UPD_AMM((psp->cmflops * nrhsa[nrhsi]) / psp->op_time,            avRmps, miRmps, maRmps );
		RSB_UPD_AMM((psp->cmflops * nrhsa[nrhsi]) / psp->at_mkl_csr_op_time, avmmps, mimmps, mammps );
		RSB_UPD_AMM((psp->cmflops * nrhsa[nrhsi]) / psp->mkl_csr_op_time,    avMmps, miMmps, maMmps );

		if( psp->nsubm && psp->at_nsubm )
		if( psp->nsubm != psp->at_nsubm )
		if( psp->isa == psp->at_isa )
		if( psp->isa != ((double)(8.0)) * psp->nnzA && psp->isa != ((double)(4.0)) * psp->nnzA )
		{
			RSB_PRWL("both auto tuned (%zd subm) and non autotuned (%zd subm) matrices use %zd bytes (%lg bpnz) of indices --- isn't that suspect ?\n",(size_t)psp->at_nsubm,(size_t)psp->nsubm,(size_t)psp->isa,((double)psp->isa)/psp->nnzA);
		}

		if( psp->nsubm && psp->at_nsubm )
		{
			if( RSB_FSTR_THN_THR(psp->at_op_time,psp->op_time,cmpt) && atwe )
			{
				const double ratio = ( psp->op_time / psp->at_op_time );
				const rsb_type_t typecode = typecodes[typecodesi];
				size_t so = RSB_SIZEOF(typecode);

				if( so == 0 ) /* reading the record of a differently configured build */
                                {
                                        rsb_type_t gtc = toupper(typecode);

                                        so = RSB_SIZEOF_BACKUP(gtc);
                                        if(wdbg)
		                                RSB_PRWL("reading file originating from a differently configured build: guessed type code size of '%c' to be %zd.\n",gtc,so);
				        if( so == 0 ) /* reading the record of a differently configured build */
                                        {
		                                RSB_PRWL(" Warning: reading file originating from a differently configured build, unable to guess correct type size for type code '%c'.\n",gtc);
                                        }
                                }

				noats ++;
				aoatsp += RSB_SPEEDUP_TO_PCT(ratio);
				aoatsr += ratio;
				RSB_UPD_TO_MAX(mstrwr, ratio);
				RSB_UPD_AMM(psp->at_t / ( psp->op_time - psp->at_op_time),aoatac,miatac,maatac);
				/* aoatuo += ( psp->at_t / psp->op_time ); */
		                RSB_UPD_AMM(psp->at_t / psp->op_time, aoatuo, miatuo, maatuo);
				/*aoatto += ( psp->at_t / psp->at_op_time );*/
		                RSB_UPD_AMM(psp->at_t / psp->at_op_time, aoatto, miatto, maatto);
				RSB_UPD_AMM( (((double)    psp->nnzA)  / psp->nsubm   ), avnzbt, minzbt, manzbt);
				RSB_UPD_AMM( (((double)    psp->nnzA)  / psp->at_nsubm), avnzat, minzat, manzat);
				RSB_UPD_AMM( (((double)(so*psp->nnzA)) / psp->nsubm   ), avbybt, mibybt, mabybt);
				RSB_UPD_AMM( (((double)(so*psp->nnzA)) / psp->at_nsubm), avbyat, mibyat, mabyat);
				RSB_UPD_AMM( (((double)(psp->isa))     / psp->nnzA    ), avbpnb, mibpnb, mabpnb);
#if RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH
				{
			const size_t so = RSB_SIZEOF_BACKUP(typecodes[typecodesi]); /* size of */
			const size_t mo = (so*psp->nnzA+psp->at_isa); /* matrix occupation */
			const size_t oo = so*nrhsa[nrhsi]*(psp->nrA+psp->ncA); /* operands occupation */
			const size_t owt = so*nrhsa[nrhsi]*RSB_ELSE_IF_TRANSPOSE(psp->nrA,psp->ncA,psp->transA); /* operands write traffic */
			const size_t mrt = oo + mo, mwt = oo + mo + owt; /* minimal read / read+write traffic */
			const rsb_perf_t mrb = ((rsb_perf_t)mrt)/(psp->at_op_time*1e9); /* minimal read bandwidth,       GBps */
			const rsb_perf_t mwb = ((rsb_perf_t)mwt)/(psp->at_op_time*1e9); /* minimal read/write bandwidth, GBps */
			RSB_UPD_AMM( mwb, avlowb, milowb, malowb);
			RSB_UPD_AMM( mrb, avlorb, milorb, malorb);
				}
				{
			const size_t so = RSB_SIZEOF_BACKUP(typecodes[typecodesi]); /* size of */
			const size_t mo = (so*psp->nnzA+((rsb_perf_t)(psp->at_isa))); /* matrix occupation */
			const size_t oo = so*nrhsa[nrhsi]*(psp->nrA+psp->ncA); /* operands occupation */
			const size_t mt = oo + mo; /* minimal traffic */
			const rsb_perf_t om = 1e6 *(psp->cmflops * nrhsa[nrhsi]); /* operation flops */
			const rsb_perf_t bm = (1.0 / om) * mt; /* bytes per mflops */;
			RSB_UPD_AMM( bm, avcoba, micoba, macoba);
				}
#endif /* RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH */
				RSB_UPD_AMM( (((double)(psp->at_isa))  / psp->nnzA    ), avbpna, mibpna, mabpna);
		                RSB_UPD_AMM(psp->at_t,                                   aosatt, misatt, masatt );

				RSB_ASSERT(psp->at_op_time != RSB_TIME_ZERO);
				RSB_ASSERT(avbyat != RSB_TIME_ZERO);
				RSB_ASSERT(avbybt != RSB_TIME_ZERO);
				RSB_ASSERT(psp->at_cn != 0);
				RSB_ASSERT(psp->at_nsubm != 0);
				RSB_ASSERT(psp->   op_time != RSB_TIME_ZERO);
				RSB_ASSERT(psp->   nsubm != 0);

				if( psp->nsubm > psp->at_nsubm )
					ntasl++;
				else
				{
					if( psp->nsubm < psp->at_nsubm )
						ntasm++;
					else
						ntase++;
				}

				if( ca[ci] > psp->at_cn )
					ntatl++;
				else
				{
					if( ca[ci] < psp->at_cn )
						ntatm++;
					else
						ntate++;
				}
			}
			else
#if 0
			if( RSB_SLWR_THN_THR(psp->at_op_time, psp->op_time) )
#endif
			{
				/* rsb tuning unsuccess */
				noatf ++;
		                RSB_UPD_AMM(psp->at_t / psp->op_time, aouatc, miuatc, mauatc );
		                RSB_UPD_AMM(psp->at_t,                aouatt, miuatt, mauatt );
			}

			if( psp->at_op_time && psp->at_mkl_csr_op_time )
			{
				/* tuned rsb ratio over tuned mkl */
				const double ratio = ( psp->at_mkl_csr_op_time / psp->at_op_time );
				nnzrm ++;
				aotrtm += ratio;
			}

			if( RSB_FSTR_THN_THR(psp->at_op_time, psp->mkl_csr_op_time, cmpt ) )
			{
				/* tuned rsb success over untuned mkl */
				const double ratio = ( psp->mkl_csr_op_time / psp->at_op_time );
				nttfu ++;
				aotstu += ratio;
				RSB_UPD_TO_MAX(mstrwu,ratio);
				RSB_UPD_AMM(psp->at_t / (psp->mkl_csr_op_time - psp->at_op_time),aoatau, miatau, maatau);
			}
			else
			if( RSB_SLWR_THN_THR(psp->at_op_time, psp->mkl_csr_op_time ) )
				nttsu ++;

			if( RSB_FSTR_THN_THR(psp->op_time, psp->mkl_csr_op_time, cmpt) )
			{
				/* untuned rsb success over mkl */
				const double ratio = ( psp->mkl_csr_op_time / psp->op_time );
				ntufm ++;
				aoustm += ratio;
				RSB_UPD_TO_MAX(msurwm,ratio);
			}
			else
			if( RSB_SLWR_THN_THR(psp->op_time, psp->mkl_csr_op_time) )
			{
				/* untuned rsb unsuccess over mkl */
				const double ratio = ( psp->op_time / psp->mkl_csr_op_time );
			       	ntusm ++;
				aoussm += ratio;
				RSB_UPD_TO_MAX(msumwr,ratio);
			}

			if( RSB_FSTR_THN_THR(psp->at_op_time, psp->at_mkl_csr_op_time, cmpt ) )
			{
				/* tuned rsb success over mkl */
				const double ratio = ( psp->at_mkl_csr_op_time / psp->at_op_time );
				nttfm ++;
				aotstm += ratio;
				RSB_UPD_TO_MAX(mstrwm,ratio);
				RSB_UPD_AMM(psp->at_t / (psp->at_mkl_csr_op_time - psp->at_op_time),aoatam, miatam, maatam);
			}
			else
			if( RSB_SLWR_THN_THR(psp->at_op_time, psp->at_mkl_csr_op_time) )
			{
				/* tuned rsb unsuccess over mkl */
				const double ratio = ( psp->at_op_time / psp->at_mkl_csr_op_time );
			       	nttsm ++;
				aotssm += ratio;
				RSB_UPD_TO_MAX(mstmwr,ratio);
			}

		} /* ... nsubm ...  */

		if( RSB_FSTR_THN_THR(psp->at_mkl_csr_op_time, psp->mkl_csr_op_time, cmpt ) )
		{
			/* mkl tuning success */
			const double ratio = ( psp->mkl_csr_op_time / psp->at_mkl_csr_op_time );
			ntmfm ++;
			aotsmm += ratio;
			RSB_UPD_TO_MAX(mstmwm,ratio);

		}
		else
		if( RSB_SLWR_THN_THR(psp->at_mkl_csr_op_time, psp->mkl_csr_op_time ) )
			ntmsm ++;
	} /* ti, typecodesi, ... */


        if(phase == 1)
        if(rds==RSB_PRD_STYLE_TBL)
        if(nocsa > 0)
	if(wltm > 0)
	{
		rsb_char_t fnbuf[RSB_MAX_FILENAME_LENGTH];
		rsb__mtxfn_bncp(fnbuf,rsb__basename(ltag),1);
		//RSB_PRT("\\hline\\end{tabular}\\end{center}CAPTION\\end{table}\n");
		RSB_PRT("\\hline\\caption{%s}\\\\\\hline\\end{longtabu}\n",ltag?fnbuf:"...");
	}

	/* FIXME: until 20160722 here was the plot case handling */

	if( nocsa <= 0 )
	{
		RSB_PRL(" No sample (out of %d) matched the dump criteria -- skipping dump round.\n",rsprp->csf);
		goto ret; /* skip any further printout */

	}
	else 
        {
	        if( rsprp->csf != nocsa )
	        {
                        RSB_PRL(" %d samples (out of %d) matched the dump limiting criteria.\n",nocsa,rsprp->csf);
                }
                if(phase == 0)
                {
                        phase = 1;
                        nocsa = 0;
			nocma = 0, lacma = -1;
                        goto gop2;
                }
        }
        if( rsprp->ror == RSB_BOOL_TRUE )
        {
                goto ret;
        }
        
        if(rds==RSB_PRD_STYLE_TBL) if(wltm > 0) RSB_PRT("\\begin{verbatim}\n");

	if( ( noats > 0 ) || ( noatf > 0 ) )
	{

		RSB_DIV_NOT_BY_ZERO(aoatsr, noats);
		RSB_DIV_NOT_BY_ZERO(aotrtm, nnzrm);
		/* RSB_DIV_NOT_BY_ZERO(mstrwr, noats); */
		RSB_DIV_NOT_BY_ZERO(aoatsp, noats);
		RSB_DIV_NOT_BY_ZERO(aoatac, noats);
		RSB_DIV_NOT_BY_ZERO(avnzbt, noats);
		RSB_DIV_NOT_BY_ZERO(avnzat, noats);
		RSB_DIV_NOT_BY_ZERO(avbybt, noats);
		RSB_DIV_NOT_BY_ZERO(avbyat, noats);
		RSB_DIV_NOT_BY_ZERO(avbpna, noats);
		RSB_DIV_NOT_BY_ZERO(avbpnb, noats);
#if RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH
		RSB_DIV_NOT_BY_ZERO(avlorb, noats);
		RSB_DIV_NOT_BY_ZERO(avcoba, noats);
#endif /* RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH */
		toatuo = aoatuo;
		RSB_DIV_NOT_BY_ZERO(aoatuo, noats);
		toatto = aoatto;
		RSB_DIV_NOT_BY_ZERO(aoatto, noats);
		touatc = aouatc;
		RSB_DIV_NOT_BY_ZERO(aouatc, noatf);
		touatt = aouatt;
		RSB_DIV_NOT_BY_ZERO(aouatt, noatf);
		tosatt = aosatt;
		RSB_DIV_NOT_BY_ZERO(aosatt, noats);

#if (RSB_WANT_MKL || RSB_WANT_ARMPL)
		/* we may spare this line also wen all mkl measurements are void. */
		RSB_PRL("above, '~' marks that rsb and mkl are close within %lgx, '.' marks that rsb is better than mkl by >%lgx, '!' marks that rsb is better than mkl by >%lgx\n",appt,rldt,hgdt);
#endif /* RSB_WANT_MKL */
		RSB_PRL("below, we define 'successful' autotuning when speedup of %lfx is exceeded, and 'tuned' results even the ones which are same as untuned\n", cmpt);
		
	       	if(noats >  0)
		{
			RSB_PRL("rsb autotuning was successful in %5d cases (%3.2lf %%) and unsuccessful in %d cases (%3.2lf %%)\n", noats, RSB_PCT(noats,noats+noatf), noatf, RSB_PCT(noatf,noats+noatf) );
			RSB_PRL(" (in succ. cases rsb autotuning gave    avg. %5.1lf %% faster, avg. sp. ratio %5.3lfx, max sp. ratio %5.3lfx, avg. ratio %5.3lfx)\n", aoatsp, aoatsr, mstrwr, aotrtm);
			//RSB_PRL(" (in succ. cases was  avg. %5.1lf %% faster, avg. sp. ratio %5.3lf, max sp. ratio %5.3lf)\n", aoatsp, aoatsr, mstrwr );
			RSB_PRL(" (in succ. cases rsb autotuning took an avg/min/max/tot of: %5.1lf/%5.1lf/%5.1lf/%5.1lf   tuned ops)\n", aoatto, miatto, maatto, toatto);
			RSB_PRL(" (in succ. cases rsb autotuning took an avg/min/max/tot of: %5.1lf/%5.1lf/%5.1lf/%5.1lf untuned ops)\n", aoatuo, miatuo, maatuo, toatuo);
	       		RSB_PRL(" (and amortizes from untuned rsb in avg. %5.1lf, min. %5.1lf, max. %5.1lf ops)\n",aoatac,miatac,maatac);
			RSB_PRL(" (avg/min/max (avg) nnz   per subm before successful tuning were %10.0lf/%10.0lf/%10.0lf)\n", avnzbt, minzbt, manzbt );
			RSB_PRL(" (avg/min/max (avg) nnz   per subm after  successful tuning were %10.0lf/%10.0lf/%10.0lf)\n", avnzat, minzat, manzat );
			RSB_PRL(" (avg/min/max (avg) bytes per subm before successful tuning were %10.0lf/%10.0lf/%10.0lf)\n", avbybt, mibybt, mabybt );
			RSB_PRL(" (avg/min/max (avg) bytes per subm after  successful tuning were %10.0lf/%10.0lf/%10.0lf)\n", avbyat, mibyat, mabyat );
			RSB_PRL(" (avg/min/max (avg) bytes per nnz  before successful tuning were %10.3lf/%10.3lf/%10.3lf)\n", avbpnb, mibpnb, mabpnb );
#if RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH
			RSB_PRL(" (avg/min/max operands (mtx,lhs,rhs) read bandwidth lower bound  %10.3lf/%10.3lf/%10.3lf,GBps)\n", avlorb, milorb, malorb);
			RSB_PRL(" (avg/min/max operands (mtx,rhs:r;lhs:rw) bandwidth lower bound  %10.3lf/%10.3lf/%10.3lf,GBps)\n", avlowb, milowb, malowb);
			RSB_PRL(" (avg/min/max code balance (bytes read at least once per flop)   %10.3lf/%10.3lf/%10.3lf)\n", avcoba, micoba, macoba);
#endif /* RSB_PRD_WANT_CODE_BALANCE_AND_BANDWIDTH */
			RSB_PRL(" (avg/min/max (avg) bytes per nnz  after  successful tuning were %10.3lf/%10.3lf/%10.3lf)\n", avbpna, mibpna, mabpna );
			RSB_PRL(" (matrix has been subdivided  more/less/same            in resp.  %d / %d /%d cases)\n", ntasm,ntasl,ntase );
			RSB_PRL(" (matrix has used             more/less/same    threads in resp.  %d / %d /%d cases)\n", ntatm,ntatl,ntate );
		}
	}
	if(noats == 0)
		RSB_PRL("no successful rsb autotuning attempt (according to %5.3lgx threshold)\n",RSB_CMP_THR );
	if(noatf == 0)
		RSB_PRL("no unsuccessful rsb autotuning attempt (according to %5.3lgx threshold) \n",RSB_CMP_THR );
	if(noatf >  0)
		RSB_PRL("unsuccessful rsb autotuning attempts (%5d cases) took avg/min/max/tot of equivalent %5.1lf/%5.1lf/%5.1lf/%5.1lf ops\n", noatf, aouatc, miuatc, mauatc, touatc );

	RSB_DIV_NOT_BY_ZERO(aotsmm, ntmfm);
	if(vscm)
	if( ntmfm || ntmsm )
		RSB_PRL("mkl threads tuning was successful in %5d cases (avg. sp. ratio %5.3lf, max sp. ratio %5.3lf) and unsuccessful in %5d cases\n", ntmfm, aotsmm, mstmwm, ntmsm);

	RSB_DIV_NOT_BY_ZERO(aoustm, ntufm);
	RSB_DIV_NOT_BY_ZERO(aoussm, ntusm);
	RSB_DIV_NOT_BY_ZERO(aotstm, nttfm);
	RSB_DIV_NOT_BY_ZERO(aotssm, nttsm);
	RSB_DIV_NOT_BY_ZERO(aoatam, nttfm);
	RSB_DIV_NOT_BY_ZERO(aoatau, nttfu);
	RSB_DIV_NOT_BY_ZERO(aotstu, nttfu);

	if(vscm)
	{
		RSB_PRL("untuned rsb has been faster than untuned mkl %5d times",ntufm);
		if( ntufm )
			RSB_PRC(", avg. sp. %2.3lf x, max %2.3lf x",aoustm,msurwm);
		RSB_PRC("\n");

		RSB_PRL("untuned rsb has been slower than untuned mkl %5d times",ntusm);
		if( ntusm )
			RSB_PRC(", avg. sl. %2.3lf x, max %2.3lf x",aoussm,msumwr);
		RSB_PRC("\n");

		RSB_PRL("tuned   rsb has been faster than   tuned mkl %5d times",nttfm);
		if( nttfm )
			RSB_PRC(", avg. sp. %2.3lf x, max %2.3lf x",aotstm,mstrwm);
		RSB_PRC("\n");
		if( nttfm )
			RSB_PRL(" (in these cases autotuning amortizes in avg. %5.1lf, min. %5.1lf, max. %5.1lf   tuned mkl ops)\n",aoatam,miatam,maatam);

		RSB_PRL("tuned   rsb has been faster than untuned mkl %5d times",nttfu);
		if( nttfu )
			RSB_PRC(", avg. sp. %2.3lf x, max %2.3lf x",aotstu,mstrwu);
		RSB_PRC("\n");
		if( nttfu )
			RSB_PRL(" (in these cases autotuning amortizes in avg. %5.1lf, min. %5.1lf, max. %5.1lf untuned mkl ops)\n",aoatau,miatau,maatau);

		RSB_PRL("tuned   rsb has been slower than   tuned mkl %5d times",nttsm);
		if( nttsm )
			RSB_PRC(", avg. sl. %2.3lf x, max %2.3lf x",aotssm,mstmwr);
		RSB_PRC("\n");
	}
        
        torott = aorott;
        toRott = aoRott;
        tomott = aomott;
        toMott = aoMott;
	toratt = aoratt;
#if 0
	tomatt = aomatt;
	RSB_DIV_NOT_BY_ZERO(aomatt, nocsa);
#endif
	RSB_DIV_NOT_BY_ZERO(aoratt, nocsa);
	RSB_DIV_NOT_BY_ZERO(aomott, nocsa);
	RSB_DIV_NOT_BY_ZERO(aoMott, nocsa);
	RSB_DIV_NOT_BY_ZERO(aorott, nocsa);
	RSB_DIV_NOT_BY_ZERO(aoRott, nocsa);
	RSB_DIV_NOT_BY_ZERO(avrmps, nocsa);
	RSB_DIV_NOT_BY_ZERO(avRmps, nocsa);
	RSB_DIV_NOT_BY_ZERO(avmmps, nocsa);
	RSB_DIV_NOT_BY_ZERO(avMmps, nocsa);
	RSB_DIV_NOT_BY_ZERO(avrsmv, ntsrf);
	RSB_DIV_NOT_BY_ZERO(avmsmv, ntsmf);

        if(noats || noatf)
	RSB_PRL("rsb auto tuning (either succ. or uns.) time was: on avg.: %5.2lf s, min %5.2lf s, max %5.2lf s, tot %5.2lf s (%d samples)\n",aoratt,miratt,maratt,toratt,nocsa );
        if(noats)
	RSB_PRL("rsb auto tuning (   only successful  ) time was: on avg.: %5.2lf s, min %5.2lf s, max %5.2lf s, tot %5.2lf s (%d samples)\n",aosatt,misatt,masatt,tosatt,noats );
        if(noatf)
	RSB_PRL("rsb auto tuning ( only unsuccessful  ) time was: on avg.: %5.2lf s, min %5.2lf s, max %5.2lf s, tot %5.2lf s (%d samples)\n",aouatt,miuatt,mauatt,touatt,noatf );
#if 0
	if(vscm)
	RSB_PRL("mkl auto tuning (either succ. or uns.) time was: on avg.: %5.2lf s, min %5.2lf s, max %5.2lf s, tot %5.2lf s (%d samples)\n",aomatt,mimatt,mamatt,tomatt,nocsa );
#endif

	if(noats) /* TODO: noats != nocsa */
	RSB_PRL(" best tun. rsb canon. mflops were: on avg. %2.3le,  min %2.3le,  max %2.3le  (%d samples)\n",avrmps, mirmps, marmps, nocsa);

	RSB_PRL(" ref. unt. rsb canon. mflops were: on avg. %2.3le,  min %2.3le,  max %2.3le  (%d samples)\n",avRmps, miRmps, maRmps, nocsa);
	if(vscm)
	{
		RSB_PRL(" best tun. mkl canon. mflops were: on avg. %2.3le,  min %2.3le,  max %2.3le  (%d samples)\n",avmmps, mimmps, mammps, nocsa);
		RSB_PRL(" ref. unt. mkl canon. mflops were: on avg. %2.3le,  min %2.3le,  max %2.3le  (%d samples)\n",avMmps, miMmps, maMmps, nocsa);
	}

	if(noats) /* TODO: noats != nocsa */
	RSB_PRL(" best tun. rsb operation time was: on avg. %2.3les, min %2.3les, max %2.3les, tot %2.3les (%d samples)\n",aorott,mirott,marott,torott, nocsa );
	RSB_PRL(" ref. unt. rsb operation time was: on avg. %2.3les, min %2.3les, max %2.3les, tot %2.3les (%d samples)\n",aoRott,miRott,maRott,toRott, nocsa );
        /* TODO: 'tot' -> 'sum' should be more appropriate */
	if(vscm)
	{
		RSB_PRL(" best tun. mkl operation time was: on avg. %2.3les, min %2.3les, max %2.3les, tot %2.3les (%d samples)\n",aomott,mimott,mamott,tomott, nocsa );
		RSB_PRL(" ref. unt. mkl operation time was: on avg. %2.3les, min %2.3les, max %2.3les, tot %2.3les (%d samples)\n",aoMott,miMott,maMott,toMott, nocsa );
	}

#if RSB_PRD_WANT_MBW
	if(noats >  0)
	if ( wpptmbw )
	{
		const int mml = RSB_MEMSCAN_DEFAULT_LEVELS; // max measurement level
		int si;
		double mrbw = 0; // memory reference bandwidth, GBps
		double crbw = 0; // cache reference bandwidth, GBps

		for(si=0;si<rsprp->mbet.sn;++si)
			if ( rsprp->mbet.et[si].lvl == mml )
			if ( rsprp->mbet.et[si].mbt == RSB_MB_MEMSET )

				mrbw = rsprp->mbet.et[si].bw / 1000.0;

		for(si=0;si<rsprp->mbet.sn;++si)
			if ( rsprp->mbet.et[si].lvl == 1 )
			if ( rsprp->mbet.et[si].mbt == RSB_MB_MEMSET )
				crbw = rsprp->mbet.et[si].bw / 1000.0;

		if ( mrbw )
		{
			// RSB_PRL(" MEMSET bandwidth in GBps: %2.3le \n",mrbw);
			// RSB_PRL(" Extrapolated min / max read bandwidth in GBps: %2.3le %2.3le\n",milorb,malorb);
			RSB_PRL(" min / max ratio of in-memory MEMSET bandwidth to extrapolated read bandwidth ratio: %2.3le %2.3le\n",mrbw/malorb,mrbw/milorb);
			if ( mrbw < malorb || mrbw < milorb )
			{
				RSB_PRL("# Warning: extrapolated memory I/O bandwidth exceeds memory bandwidth --- is this a tiny matrix ?\n");
			}

			if ( crbw )
			{
   		     		RSB_ASSERT( crbw > mrbw );
				RSB_PRL(" in-cache to in-memory MEMSET bandwidth ratio: %2.3le\n",crbw/mrbw);
				if ( crbw < malorb || crbw < milorb )
				{
					RSB_PRL(" min / max ratio of in-cache MEMSET bandwidth to extrapolated read bandwidth ratio: %2.3le %2.3le\n",crbw/malorb,crbw/milorb);
					RSB_PRL("# Warning: extrapolated memory I/O bandwidth exceeds cache bandwidth!\n");
					// errval = RSB_ERR_INTERNAL_ERROR; // Note: maybe introduce a strict mode have an error emitted here.
					RSB_PERR_GOTO(ret,"Error: extrapolated memory I/O bandwidth exceeds cache bandwidth!\n");
				}
			}
		}
		else
			; // something must have gone wrong..
	}
#endif /* RSB_PRD_WANT_MBW */

	if(ntsrf > 0)
		RSB_PRL(" rsb nrhs-to-overall-min-rhs speed ratio was: on avg.    %2.3le x, min %2.3le x, max %2.3le x (%d samples, the non-min-nrhs ones)\n",avrsmv, mirsmv, marsmv, ntsrf);
	if(vscm) /* vscm does not properly apply here; but ntsmf alone is not enough */
	if(ntsmf)
		RSB_PRL(" mkl nrhs-to-overall-min-rhs speed ratio was: on avg.    %2.3le x, min %2.3le x, max %2.3le x (%d samples, the non-min-nrhs ones)\n",avmsmv, mimsmv, mamsmv, ntsmf);

#define RSB_PR_LOOP_EL(IVAR,CVAL,UVAL) for(IVAR=(CVAL==-1?0:CVAL);IVAR<(CVAL==-1?UVAL:(CVAL+1));++IVAR)
#define RSB_PR_LOOP(FNV,CNV,IXV,IYV,NRV,TCV,TNV) \
	RSB_PR_LOOP_EL(filenamei,FNV,rsprp->filenamen) \
	RSB_PR_LOOP_EL(ci,CNV,rsprp->cn) \
	RSB_PR_LOOP_EL(incXi,IXV,rsprp->incXn     ) \
	RSB_PR_LOOP_EL(incYi,IYV,rsprp->incYn     ) \
	RSB_PR_LOOP_EL(nrhsi,NRV,rsprp->nrhsn) \
	RSB_PR_LOOP_EL(typecodesi,TCV,rsprp->ntypecodes) \
	RSB_PR_LOOP_EL(ti,TNV,rsprp->tn) \
	{ \
	} /* TODO: still unused; complete this... */

	if(rsprp->filenamen > 0)
		;/* plot for all matrices, performance for each case */
	if(rsprp->cn > 0)
		;/* plot for all matrices, performance for increasing cores */
	if(rsprp->incXn > 0)
		;/* plot per matrix, performance for increasing incX */
	if(rsprp->incYn > 0)
		;/* plot per matrix, performance for increasing incY */
	if(rsprp->nrhsn > 0)
		;/* plot per matrix, performance for increasing nrhs */
	if(rsprp->ntypecodes > 0)
		;/* plot per matrix, performance for different types */
	if(rsprp->tn > 0)
		;/* plot per matrix, different transpositions: TODO: incorporate in the usual plots .. */

	/* TODO: for each matrix, performance + MKL performance / speedup */
	/* TODO: for each matrix, thread tuned performance to normal performance / speedup */
	/* TODO: for each matrix, structure tuned performance to normal performance / speedup */
	/* TODO: for each matrix, total tuned performance to normal performance / speedup */
	/* TODO: label each of these with timestamp;
	 *
	 * plot-DATE-mtrx-TYPE-CORES-INCX-INCY-NRHS-TRANS.eps
	 * plot-DATE-nrhs-TYPE-CORES-MTRX-INCX-INCY-TRANS.eps 
	 * plot-DATE-incx-TYPE-CORES-MTRX-INCY-NRHS-TRANS.eps 
	 * plot-DATE-incy-TYPE-CORES-MTRX-INCX-NRHS-TRANS.eps 
	 * ...
	 * */

	/* wishlist: */
	/* plot per matrix, then different indexing per nonzero */
	/* plot for all matrices, then different indexing per nonzero */
        if(rds==RSB_PRD_STYLE_TBL) if(wltm > 0) RSB_PRT("\\end{verbatim}\n");

	/* begin plot subcase */
        if(phase == 1) /* need the statistics to be ready */
        if(nocsa > 0)
        if(rds>=RSB_PRD_STYLE_PLT_BASE)
        {
		const rsb_char_t * pl = NULL;
		const rsb_char_t * ppl = "";
		rsb_char_t pfn[RSB_MAX_FILENAME_LENGTH];
		rsb_int_t nocsai = 0; /* number of considered samples, index */

        	if(nocsa < 2) /* at least two samples per plot */
                	goto ret;

                if(rsb__getenv("RSB_PRD_STYLE_PLT_PFN"))
			ppl = rsb__getenv("RSB_PRD_STYLE_PLT_PFN");

                if( rsb__util_atoi(rsb__getenv("RSB_PRD_STYLE_PLT_FMT")) )
		{
			pl = "set term postscript eps color size 2,2 noclip font \"Times-Roman,14\";";
			sprintf(pfn,"%s%s.eps",ppl,ltag?ltag:"plot");
		}
		else
		{
			pl = "set term png;";
			sprintf(pfn,"%s%s.png",ppl,ltag?ltag:"plot");
		}

        	if(rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB)
		{
       			RSB_STDOUT("%sset output '%s'; set title 'autotuning effect'; unset ytics;set yrange [0: 2];\n",pl,pfn);
			RSB_STDOUT("plot '-' using 1:2 title 'rsb' lt rgb 'red'\n");
			RSB_STDOUT("set xlabel 'speedup'\n");
			RSB_STDOUT("set ylabel ' '\n");
		}

        	if(rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB_POLAR || rds==RSB_PRD_STYLE_PLT_MKL_SPEEDUP_POLAR)
		{
			// TODO: need to modify this
			double rval, mrval, arval;
			const char * title = NULL, *xlabel = NULL, *ylabel = ltag ? ltag: "";
			char xlabelbuf[RSB_MAX_LABEL_LENGTH];

        		if(rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB_POLAR)
				//rval = RSB_MAX(ceil(aoatsr),2.0), mrval = aoatsr, // FIXME: this is average, we want max for mrval
				rval = RSB_MAX(ceil(mstrwr),2.0),
				mrval = mstrwr,
				arval = aoatsr, // FIXME: this is average, we want max for mrval
				title = "autotuning effect",
				xlabel="''";

        		if(rds==RSB_PRD_STYLE_PLT_MKL_SPEEDUP_POLAR)
				rval = RSB_MAX(ceil(mstrwm),2.0),
				arval = aotrtm,
				mrval = mstrwm, 
				title = "RSB to " RSB_MKL_S " speed ratio",
				sprintf(xlabelbuf,"\"(avg impr. is %3.2lfx, max impr. is %3.2lfx,\\n avg. ratio. is %3.2lfx)\"",aotstm,mrval,arval),
				xlabel=xlabelbuf;

       			//RSB_STDOUT("#polar plot instructions RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB_POLAR\n");
       			RSB_STDOUT("\
# PLOT BEGIN #\n\
my_rval=%lg # max ratio is %lg !\n\
my_clen=2*pi\n\
my_nsam=%d # number of samples\n\
my_nmat=%d # number of matrices\n\
my_rnge=my_rval*1.2\n\
my_title='RSB'\n\
#my_size=600 # pixels\n\
#set term png size my_size,my_size\n\
set key noinvert samplen 0.75 spacing 1 width 0 height 0 at graph 1.0,1.0\n\
set title '%s';\n\
set xlabel %s\n\
set ylabel '[%s]'\n\
#my_avg_str(x) = sprintf(\"avg: ... x\")\n\
my_avg(x) = %lg \n\
my_dir=\"%s\"\n\
%s\n\
",
	rval,mrval,nocsa,nocma,title,xlabel,ylabel,arval,
	(fprfn&&*fprfn)?fprfn:".",
	(fprfn&&*fprfn)?"system('mkdir -p '.my_dir)":"");

       			RSB_STDOUT("\
set polar\n\
#rgb_type(t) = ( t eq 'D' ) ? red : ( ( t eq 'Z' ) ? blue : (( t eq 'S' ) ? green : black )  )\n\
#my_avg(v,l) = sprintf(\"avg: %%.2f%%s\",v,l)\n\
#set grid polar min(my_clen/my_nsam,2*pi/my_maxnsam)\n\
max_nsec=36 # after this won't draw sectors\n\
#my_nsec=my_nsam # one sector per sample\n\
my_nsec=my_nmat # one sector per matrix\n\
my_pangle = ((my_clen/my_nsec)>((2*pi)/max_nsec)?(my_clen/my_nsec):2*pi)\n\
set grid polar my_pangle\n\
set grid layerdefault linetype 0 linewidth 1.0, linetype 0 linewidth 4.0\n\
set grid noxtics nomxtics noytics nomytics noztics nomztics nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics\n\
unset xtics\n\
unset ytics\n\
unset border\n\
set origin 0.0, 0.0;\n\
set rtics autofreq \n\
set rrange [ 0 : my_rnge ] noreverse nowriteback \n\
");
       			RSB_STDOUT("\
%s\nset output my_dir.'/%s' \n\
",pl,pfn);
       			RSB_STDOUT("\
my_arry = my_rval\n\
my_arrx = my_avg(-1)\n\
my_max(x,y) = ( x < y ? x : y)\n\
my_sposl(x) = ( x < 1.0 ? '(slowdown)' : '(speedup)')\n\
set arrow from 1,my_arry to my_arrx,my_arry ls 5 lw 0.4 lc rgbcolor 'black' front\n\
set arrow nohead from           my_arrx,my_arry to           my_arrx,0 lw .4 lt 0 lc 'black' front\n\
set arrow nohead from           1      ,my_arry to           1      ,0 lw .4 lt 0 lc 'black' front\n\
set label sprintf(' %%.2fx %%s',my_arrx,my_sposl(my_arrx)) at my_max(1,my_arrx),my_arry*1.05 front font 'Times-Roman,10'\n\
");
       			RSB_STDOUT("\
set yrange [-my_rval: my_rval];\n\
set xrange [-my_rval: my_rval];\n\
set multiplot\n\
my_r(x)=(x+.5)*my_clen/my_nsam\n\
my_v(x)=1*x\n\
red='#dd0000'\n\
#green='#00dd00'\n\
#black='#000000'\n\
my_i_argb(r,g,b) = 0 + 65536 * int(r) + 256 * int(g) + int(b)\n\
my_i_red=my_i_argb(255,0,0)\n\
my_i_green=my_i_argb(0,255,0)\n\
my_i_blue=my_i_argb(0,0,255)\n\
my_i_black=my_i_argb(0,0,0)\n\
my_rgb_type_s(t) = ( t eq 'D' ) ? my_i_red : (( t eq 'Z' ) ? my_i_blue: ((t eq 'S') ? my_i_green:my_i_black ))\n\
my_rgb_symm_s(s) = ( s eq 'S' ) ? my_i_red : my_i_black \n\
log2(n) = log(n)/log(2.0) # FIXME\n\
my_rgb_nrhs_s(nrhs) = ( nrhs == 1 ) ? my_i_red : my_i_black \n\
#my_rgb_nrhs_s(nrhs) = int(log2(nrhs))\n\
#my_rgb_nrhs_s(nrhs) = nrhs\n\
my_rgb_cols(nc,tc,sc) = my_rgb_type_s(stringcolumn(tc)) # type->color\n\
#my_rgb_cols(nc,tc,sc) = my_rgb_symm_s(stringcolumn(sc)) # symm->color\n\
#my_rgb_cols(nc,tc,sc) = my_rgb_nrhs_s(column(nc)) # nrhs->color\n\
#my_rgb_type_col(tc) = my_i_red # type->color\n\
plot 1 notitle with filledcurves below linetype 1 linewidth 0.000 linecolor rgb '#dddddd' \n\
plot '-' using ((my_r($2))):((my_v($1))):((my_rgb_cols(3,4,5))) title my_title lc rgbcolor variable ps 1 pt 6,\
	my_avg(-1) notitle lt 0 lc rgbcolor red\
\n\
");
		}

        	if(rds==RSB_PRD_STYLE_PLT_SUBM_BS)
		{
       			RSB_STDOUT("%sset output '%s';",pl,pfn);
		       	//RSB_STDOUT("set title 'autotuning effect'; unset ytics;set yrange [0: 3];\n");
			RSB_STDOUT("set xlabel 'bytes per submatrix'\n");
			RSB_STDOUT("set ylabel 'performance, Mflops/s'\n");
			RSB_STDOUT("set xtics rotate by -45\n");
			RSB_STDOUT("plot '-' using 1:2:3:4 with vectors title 'rsb' lt rgb 'red'\n");
		}

        	if(rds==RSB_PRD_STYLE_PLT_SUBM_BS_POLAR)
		{
			// TODO: need to modify this
       			RSB_STDOUT("#polar plot instructions RSB_PRD_STYLE_PLT_SUBM_BS_POLAR (unfinished)\n");
		}

		for(     filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
		for(ci=0;ci<rsprp->cn;++ci)
		for(     incXi=0;     incXi<rsprp->incXn     ;++incXi     )
		for(     incYi=0;     incYi<rsprp->incYn     ;++incYi     )
		for(     nrhsi=0;     nrhsi<rsprp->nrhsn     ;++nrhsi     )
		for(typecodesi=0;typecodesi<rsprp->ntypecodes;++typecodesi)
		for(ti=0;ti<rsprp->tn;++ti)
		{
			const size_t idx = rsb__pr_idx(rsprpv, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti);
			const struct rsb_rsps_t*psp = &(rsprp->psa[idx]);
			const rsb_type_t typecode = typecodes[typecodesi];
			const size_t so = RSB_SIZEOF_BACKUP(toupper(typecode));

			if( RSB_ERR_NO_ERROR != rsb__pr_filter(psp, ta, filenameifp, ifilenameifp, cifp , incXifp , incYifp , nrhsifp , typecodefip , tifp, tfp, flagsA, nflagsA, filenamei, ci, incXi, incYi, nrhsi, typecodesi, ti) )
				continue;
        		if(rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB)
        			RSB_STDOUT("%le %d\n",psp->op_time/psp->at_op_time,1);
        		if(rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB_POLAR)
        			RSB_STDOUT("%le %d\n",psp->op_time/psp->at_op_time,nocsai);
        		if(rds==RSB_PRD_STYLE_PLT_MKL_SPEEDUP_POLAR)
        			RSB_STDOUT("%le %d %d %c %c\n",psp->mkl_csr_op_time/psp->at_op_time,nocsai,nrhsa[nrhsi],typecodes[typecodesi],RSB_SYMCHAR(psp->flagsA));
        		if(rds==RSB_PRD_STYLE_PLT_SUBM_BS)
			{
				const double avmbybt = (((double)(so*psp->nnzA)) / psp->nsubm   );
				const double avmbyat = (((double)(so*psp->nnzA)) / psp->at_nsubm);
				//double avmbpnb = (((double)(psp->isa))     / psp->nnzA    );
				//double avmbpna = (((double)(psp->at_isa))     / psp->nnzA    );
				const double avmrmps = ((psp->cmflops * nrhsa[nrhsi]) / psp->at_op_time);
		     		const double avmRmps = ((psp->cmflops * nrhsa[nrhsi]) / psp->op_time);
				//double avmmmps = ((psp->cmflops * nrhsa[nrhsi]) / psp->at_mkl_csr_op_time);
			     	//double avmMmps = ((psp->cmflops * nrhsa[nrhsi]) / psp->mkl_csr_op_time);

        			RSB_STDOUT("%le %le %le %le\n",avmbybt,avmRmps,avmbyat,avmrmps-avmRmps);
        			//RSB_STDOUT("%le %d %le %d\n",avmbybt,1,avmbyat,1);
        			//RSB_STDOUT("%le %le %le %le\n",avmbybt,avmbpnb,avmbyat,(avmbpna-avmbpnb));
			}
			nocsai++;
		}
        	if(rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB || rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB_POLAR || rds==RSB_PRD_STYLE_PLT_MKL_SPEEDUP_POLAR)
        		RSB_STDOUT("e\n");
        	if(rds==RSB_PRD_STYLE_PLT_AT_SPEEDUP_RSB_POLAR || rds==RSB_PRD_STYLE_PLT_MKL_SPEEDUP_POLAR)
        		RSB_STDOUT("unset multiplot;unset label;unset arrow;\n# PLOT END\n\n");
        	if(rds==RSB_PRD_STYLE_PLT_SUBM_BS)
        		RSB_STDOUT("e\n");
		goto ret;
        }
	/* end plot subcase */
ret:
	return errval;
} /* rsb__pr_dump_inner */

rsb_err_t rsb__pr_dump(const void*rsprpv, rsb_char_t**RSB_RESTRICT filenamea, rsb_int_t*ca, const rsb_int_t*incXa, const rsb_int_t*incYa, const rsb_int_t*nrhsa, const rsb_type_t*typecodes, const rsb_int_t*ta, const rsb_char_t *fprfn)
{
	/*
	 * dump a performance record
         * TODO: use rsb__basename().
         * TODO: use a systematic combinations enumeration algorithm.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const int*filenameifp = NULL; /* filename index ? pointer */
	const int*cifp = NULL; /* [used] cores index ? pointer */
	const int*incXifp = NULL; /* incX index ? pointer */
	const int*incYifp = NULL; /* incY index ? pointer */
	const int*nrhsifp = NULL; /* nrhs index ? pointer */
	const int*typecodefip = NULL; /* typecode index ? pointer */
	const int*tifp = NULL; /* transposition index ? pointer */
	const rsb_trans_t*tfp = NULL; /* transposition ? pointer */
	const struct rsb_rspr_t * rsprp = rsprpv;
	rsb_int_t filenamei, /* ci, incXi, incYi,*/ nrhsi, typecodesi, ti = 0;
	rsb_trans_t transAa [] = { RSB_TRANSPOSITION_N, RSB_TRANSPOSITION_T, RSB_TRANSPOSITION_C };
	rsb_char_t tag[2*RSB_MAX_FILENAME_LENGTH];
	rsb_char_t bfn[RSB_MAX_FILENAME_LENGTH];
	rsb_int_t noc = 0; /* number of combinations */
	char rsb_prl_lcc = RSB_PRL_LCC_IE ;
	const char*rsb_prl_tcs = RSB_PRL_TCS_IE;
        int rds = rsb__getenv_real_t("RSB_PR_SR",RSB_PRD_STYLE_TBL);
        int wltm = rsb__getenv_real_t("RSB_PR_WLTC",0); /* Want LaTeX tables mode */
	const rsb_bool_t only_total_table = rsb__util_atoi(rsb__getenv("RSB_PR_ONLY_TOTAL_TABLE"));

        if(!rsprp)
	{
		errval = RSB_ERR_BADARGS;
	       	goto err;
	}

        if(rds==RSB_PRD_STYLE_TBL && wltm > 0)
		RSB_PRT( "\\documentclass[a1,portrait,plainsections]{sciposter} \\usepackage{longtable,tabu,url,color} \\usepackage[cm]{fullpage} \\usepackage[table,x11names]{xcolor} \\usepackage[hyperindex,bookmarks]{hyperref}%% bookmarks do not seem to work\n\\begin{document}\\title{" RSB_PACKAGE_NAME " performance, postprocessed with " RSB_PACKAGE_STRING ".}\\author{} \\begin{tiny} \\rowcolors{1}{white!80!gray}{white}\n");

	noc = RSB_PR_NOC(rsprp);

	if( only_total_table )
		goto total;

	if( filenamea )
	if( rsprp->filenamen > 1 )
	if( rsprp->filenamen < noc )
	for(     filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
	{
                rsb_int is_symm = RSB_DO_FLAG_HAS(((rsprp->psa[rsb__pr_idx(rsprpv, filenamei, 0, 0, 0, 0, 0, 0)]).flagsA ),RSB_FLAG_SYMMETRIC);
        	int etn = rsprp->tn;
		rsb__mtxfn_bncp(bfn,rsb__basename(filenamea[filenamei]),0);
		sprintf(tag,"file-%d-%s",filenamei+1,bfn);
		RSB_PRL_SEP(" Limiting to file %d/%d --- %s:\n",filenamei+1,rsprp->filenamen,filenamea[filenamei]);
		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
					  &filenamei, NULL,cifp, incXifp, incYifp, nrhsifp, typecodefip, tifp, tfp, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn);

                if(is_symm && rsprp->tn == 2 )
                        etn = 1;

        	if( etn > 1 ) /* otherwise the above dump suffices */
        	if( etn * rsprp->filenamen < noc )
        	for(     ti=0;     ti<RSB_MIN(etn,3)     ;++ti     )
        	{
        		const rsb_trans_t tf = transAa[ti];
			rsb__mtxfn_bncp(bfn,rsb__basename(filenamea[filenamei]),0);
			sprintf(tag,"file-%d-%s-transA-%c",filenamei+1,bfn,RSB_TRANSPOSITION_AS_CHAR(tf));
        		RSB_PRL_SEP(" Limiting to both file %d/%d --- %s and transA=%c:\n",filenamei+1,rsprp->filenamen,filenamea[filenamei],RSB_TRANSPOSITION_AS_CHAR(tf));
        		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
					  &filenamei, NULL,cifp, incXifp, incYifp, nrhsifp, typecodefip, &ti, &tf, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn);
        	}

	        if( typecodes )
        	if( rsprp->ntypecodes > 1 )
        	if( rsprp->ntypecodes * rsprp->filenamen < noc )
        	for(typecodesi=0;typecodesi<rsprp->ntypecodes;++typecodesi)
        	{
			rsb__mtxfn_bncp(bfn,rsb__basename(filenamea[filenamei]),0);
			sprintf(tag,"file-%d-%s-type-%c",filenamei+1,bfn,typecodes[typecodesi]);
        		RSB_PRL_SEP(" Limiting to both file %d/%d --- %s and type %c:\n",filenamei+1,rsprp->filenamen,filenamea[filenamei],typecodes[typecodesi]);
        		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
					  &filenamei, NULL,cifp, incXifp, incYifp, nrhsifp, &typecodesi, tifp, tfp, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn);
        	}
       	}

	if( filenamea )
	if( rsprp->filenamen > 1 )
	{
		rsb_flags_t sf[3] = {RSB_FLAG_NOFLAGS,RSB_FLAG_NOFLAGS,RSB_FLAG_NOFLAGS}, nf[3] = {RSB_FLAG_NOFLAGS,RSB_FLAG_NOFLAGS,RSB_FLAG_NOFLAGS};
		rsb_char_t scb[3] = {0x0,0x0,0x0};
		int sc = 0, si;
		rsb_flags_t flag;

		for(     flag = RSB_FLAG_SYMMETRIC, filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
			if(rsprp->psa[rsb__pr_idx(rsprpv, filenamei, 0, 0, 0, 0, 0, 0)].flagsA & flag)
			{
				scb[sc] = RSB_SYMCHAR(flag); nf[sc] = RSB_FLAG_NOFLAGS;  sf[sc++] = flag; break;
			}

		for(     flag = RSB_FLAG_HERMITIAN, filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
			if(rsprp->psa[rsb__pr_idx(rsprpv, filenamei, 0, 0, 0, 0, 0, 0)].flagsA & flag)
			{
				scb[sc] = RSB_SYMCHAR(flag); nf[sc] = RSB_FLAG_NOFLAGS;  sf[sc++] = flag; break;
			}

		for(     flag = (RSB_FLAG_SYMMETRIC|RSB_FLAG_HERMITIAN), filenamei=0;     filenamei<rsprp->filenamen ;++filenamei     )
			if(rsprp->psa[rsb__pr_idx(rsprpv, filenamei, 0, 0, 0, 0, 0, 0)].flagsA & flag)
				;
			else
			{
				scb[sc] = RSB_SYMCHAR(RSB_FLAG_NOFLAGS); sf[sc] = RSB_FLAG_NOFLAGS;  nf[sc++] = flag; break;
			}

                if(sc > 1) /* TODO: shall be strictier and check whether at least more than one matrix applies for each given sc */
                if(sc * rsprp->filenamen < noc)
	        for(     si = 0; si < sc ; ++si )
        	{
			sprintf(tag,"symmetry-%c",scb[si]);
        		RSB_PRL_SEP(" Limiting to symmetry %c (0x%x) \n",scb[si],sf[si]);
        		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
        					  filenameifp, NULL,cifp, incXifp, incYifp, nrhsifp, typecodefip, tifp, tfp, sf[si], nf[si], tag, fprfn);
        	}
	}

	if( typecodes )
	if( rsprp->ntypecodes > 1 )
	if( rsprp->ntypecodes < noc )
	for(typecodesi=0;typecodesi<rsprp->ntypecodes;++typecodesi)
	{
		sprintf(tag,"type-%c",typecodes[typecodesi]);
		RSB_PRL_SEP(" Limiting to type %c:\n",typecodes[typecodesi]);
		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
					  filenameifp, NULL,cifp, incXifp, incYifp, nrhsifp, &typecodesi, tifp, tfp, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn);
	}

	if( nrhsa )
	if( rsprp->nrhsn > 1 )
	if( rsprp->nrhsn < noc )
	for(     nrhsi=0;     nrhsi<rsprp->nrhsn     ;++nrhsi     )
	{
		sprintf(tag,"nrhs-%d",nrhsa[nrhsi]);
		RSB_PRL_SEP(" Limiting to nrhs=%d:\n",nrhsa[nrhsi]);
		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
					  filenameifp, NULL,cifp, incXifp, incYifp, &nrhsi, typecodefip, tifp, tfp, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn );
	}

	if( ta && rsprp->tn > 1 )
        {
	if( ta && rsprp->tn < noc )
	for(     ti=0;     ti<rsprp->tn     ;++ti     ) /** FIXME: why this case ? */
	{
		sprintf(tag,"transA-%c",RSB_TRANSPOSITION_AS_CHAR(ta[ti]));
		RSB_PRL_SEP(" Limiting to transA=%d:\n",RSB_TRANSPOSITION_AS_CHAR(ta[ti]));
		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
					  filenameifp, NULL,cifp, incXifp, incYifp, nrhsifp, typecodefip, &ti, tfp, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn );
	}
        }
	else
	if( rsprp->tn > 1 )
        {
	if( rsprp->tn < noc )
	for(     ti=0;     ti<RSB_MIN(rsprp->tn,3)     ;++ti     )
	{
		const rsb_trans_t tf = transAa[ti];
		sprintf(tag,"transA-%c",RSB_TRANSPOSITION_AS_CHAR(tf));
		RSB_PRL_SEP(" Limiting to transA=%c:\n",RSB_TRANSPOSITION_AS_CHAR(tf));
		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
					  filenameifp, NULL, cifp, incXifp, incYifp, nrhsifp, typecodefip, &ti, &tf, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn );
        	if( nrhsa )
        	if( rsprp->nrhsn > 1 )
	        if( rsprp->tn * rsprp->nrhsn < noc )
        	for(     nrhsi=0;     nrhsi<rsprp->nrhsn     ;++nrhsi     )
        	{
			sprintf(tag,"transA-%c-nrhs-%d",RSB_TRANSPOSITION_AS_CHAR(tf),nrhsa[nrhsi]);
        		RSB_PRL_SEP(" Limiting to both transA=%c and nrhs=%d:\n",RSB_TRANSPOSITION_AS_CHAR(tf),nrhsa[nrhsi]);
        		errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
        					  filenameifp, NULL,cifp, incXifp, incYifp, &nrhsi, typecodefip, &ti, &tf, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn );
        	}
	}
        }

total:
	RSB_PRL_SEP(" All results (not limiting)\n");
	sprintf(tag,"all");
	errval = rsb__pr_dump_inner(rsprpv, filenamea, ca, incXa, incYa, nrhsa, typecodes, ta, 
						  filenameifp, NULL,cifp, incXifp, incYifp, nrhsifp, typecodefip, tifp, tfp, RSB_FLAG_NOFLAGS, RSB_FLAG_NOFLAGS, tag, fprfn);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

        if(rds==RSB_PRD_STYLE_TBL && wltm > 0)
		RSB_PRT("\\end{tiny}\\end{document}\n");

#if RSB_PRD_WANT_TIMESTAMP
	RSB_PRL("Record collection took %5.2lf s.\n",rsprp->tend-rsprp->tbeg);
#endif /* RSB_PRD_WANT_TIMESTAMP */
#if RSB_PRD_WANT_MBW
	RSB_PRL("Record comprises %d memory benchmark samples (prepend RSB_PR_MBW=1 to dump this).\n",rsprp->mbet.sn);
	if ( rsprp->mbet.sn > 0 && rsb__util_atoi(rsb__getenv("RSB_PR_MBW")) == 1 )
		rsb__mbw_es_print(&rsprp->mbet);
#endif /* RSB_PRD_WANT_MBW */
#if RSB_PRD_WANT_ENVIRONMENT
	RSB_PRL("Record comprises %d environment variables in %d bytes (prepend RSB_PR_ENV=1 to dump this).\n",(int)rsprp->nenvv,(int)rsprp->enoib);
	if ( rsb__util_atoi(rsb__getenv("RSB_PR_ENV")) == 1 )
	if ( rsprp->nenvv )
	{
		rsb_int_t evo = 0;
		rsb_int_t envvi = 0;

		for( evo = 0; evo < rsprp->enoib+rsprp->nenvv ; evo += strlen(rsprp->envvp+evo)+1 )
		{
			//RSB_PRL("%s\n",rsprp->envvp+evo);
			//RSB_PRL("%5d/%5d %s\n",evo,rsprp->enoib,rsprp->envvp+evo);
			RSB_PRL("%5d/%5d %s\n",envvi,rsprp->nenvv,rsprp->envvp+evo);
			envvi++;
		}
	}
#endif /* RSB_PRD_WANT_ENVIRONMENT */
/*
	for(ci=0;ci<rsprp->cn;++ci)
	for(     incXi=0;     incXi<rsprp->incXn     ;++incXi     )
	for(     incYi=0;     incYi<rsprp->incYn     ;++incYi     )
*/
err:
	return errval;
} /* rsb__pr_dump */

rsb_err_t rsb__pr_free(void * rsprpv)
{
	/*
	 * free a performance record
	 * */
	struct rsb_rspr_t * rsprp = rsprpv;
        if(!rsprp)
                goto err;

	RSB_CONDITIONAL_FREE(rsprp->psa);
	RSB_CONDITIONAL_FREE(rsprp->rsprap);
#if RSB_PRD_WANT_MBW
	// RSB_CONDITIONAL_FREE(rsprp->mbet.et);
	rsb__mbw_es_free(&rsprp->mbet);
#endif /* RSB_PRD_WANT_MBW */
#if RSB_PRD_WANT_ENVIRONMENT
	RSB_CONDITIONAL_FREE(rsprp->envvp);
#endif /* RSB_PRD_WANT_ENVIRONMENT */
	RSB_CONDITIONAL_FREE(rsprp);
err:
	return RSB_ERR_NO_ERROR;
}

/* performance samples reporting / dumping facility for rsbench : end */
/* @endcond */
