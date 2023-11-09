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
/*!
 * @file
 * @brief  Recursive Sparse matrices assembling code.
 * @author Michele Martone
 * */
/*
 * TODO: improve this code, because there is a number of unclean practices which could break a build.
 * */
#include "rsb_common.h"

#ifndef RSB_C2R_ASSERT
//#define RSB_C2R_ASSERT(e) assert(e)		// uncomment this to use   asserts
#define RSB_C2R_ASSERT(e)			// uncomment this to avoid asserts
#else /* RSB_C2R_ASSERT */
#undef RSB_C2R_ASSERT
#define RSB_C2R_ASSERT(e) 
#endif /* RSB_C2R_ASSERT */

#define RSB_DO_ENOUGHNNZFORINDEXBASEDBUILD(M) (!RSB_DO_TOOFEWNNZFORRCSR((M)->nnz,(M)->nr))
#define RSB_C2R_IF_VERBOSE 0	/* activates output which is useful for debugging */
#define RSB_C2R_PARANOIA 0

#define RSB_MEMCPY_SMALL_GENERAL(ID,IS,DOFF,SOFF,NNZ,TYPE) \
	{ \
		TYPE*dp = ((TYPE*)(ID))+(DOFF),*ld = dp+(NNZ); \
		const register TYPE*sp = ((TYPE*)(IS))+(SOFF); \
		for(;dp<ld;++sp,++dp)*dp = *sp; \
       	}

#define RSB_C2R_WANT_MAYBE_FASTER 0

#if RSB_C2R_WANT_MAYBE_FASTER 
#define RSB_COA_MEMCPY_SMALL(ID,IS,DOFF,SOFF,NNZ) RSB_MEMCPY_SMALL_GENERAL(ID,IS,DOFF,SOFF,NNZ,rsb_coo_idx_t)
#define RSB_COA_MEMCPY_ROWSZ(ID,IS,DOFF,SOFF,NNZ) RSB_COA_MEMCPY_parallel(ID,IS,DOFF,SOFF,NNZ)
#define RSB_A_MEMCPY_SMALL(ID,IS,DOFF,SOFF,NNZ,ES) RSB_A_MEMCPY(ID,IS,DOFF,SOFF,NNZ,ES) 
#else /* RSB_C2R_WANT_MAYBE_FASTER */
#define RSB_COA_MEMCPY_SMALL(ID,IS,DOFF,SOFF,NNZ) RSB_COA_MEMCPY(ID,IS,DOFF,SOFF,NNZ) 
#define RSB_A_MEMCPY_SMALL(ID,IS,DOFF,SOFF,NNZ,ES) RSB_A_MEMCPY(ID,IS,DOFF,SOFF,NNZ,ES) 
#endif /* RSB_C2R_WANT_MAYBE_FASTER */

#define RSB_TIC(T) (T) = -rsb_time()
#define RSB_TOC(T) (T) += rsb_time()
#define RSB_TOC_TIC(T,U) {rsb_time_t t = rsb_time();(T) += t;(U) = -t;}
#define RSB_TIC_TOC(U,T) {rsb_time_t t = rsb_time();(T) += t;(U) = -t;}

#define RSB_WANT_BINSEARCH_MIN_NZPR 8	/* FIXME */
#define RSB_WANT_VERBOSE_TIMINGS 0
#define RSB_WANT_VERBOSE_SUBDIVISION 0	/* */
#define RSB_WANT_VERBOSE_SUBDIVISION2 0	/* */
#define RSB_WANT_MORE_PARALLELISM 1	/* DEFAULT ON (1) */
#define RSB_WANT_FIRST_VERSION 0	/* FIXME: EXPERIMENTAL, BY DEFAULT TURNED ON (1) */
#define RSB_WANT_LITTLE_IMPROVED 1	/* FIXME: EXPERIMENTAL, BY DEFAULT TURNED OFF (0) */
#if RSB_WANT_OMP_RECURSIVE_KERNELS
#define RSB_WANT_PARALLEL_SUBDIVISION 1	/* FIXME: EXPERIMENTAL, BY DEFAULT TURNED OFF (0) */
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#define RSB_WANT_PARALLEL_SUBDIVISION 0	/* FIXME: EXPERIMENTAL, BY DEFAULT TURNED OFF (0) */
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#define RSB_WANT_QUADRANT_QUICK_DETECT 0/* FIXME: EXPERIMENTAL, BY DEFAULT TURNED OFF (0) */
#define RSB_WANT_SUBDIVISION_FIXES_20101120 1 /* FIXME: EXPERIMENTAL, BY DEFAULT TURNED OFF (0) */
#define RSB_WANT_SUBDIVISION_FIXES_20101213 0 /* FIXME: EXPERIMENTAL (DEFAULT OFF) */
#define RSB_WANT_FIX_BUG_DISCOVERED_20121210 1	/* this bug prevents from HCSR usage */
#if RSB_WANT_VERBOSE_SUBDIVISION2
//#define RSB_MTXASM_INFO RSB_INFO	/* NEW */
#define RSB_MTXASM_INFO printf	/* NEW */
#else /* RSB_MTXASM_INFO */
#define RSB_MTXASM_INFO 	/* NEW */
#endif /* RSB_MTXASM_INFO */

#define RSB_SUBDIVISION_SKEW_MAX  (RSB_FLOAT_ONE/2.0)
#define RSB_MAX_QUADRANTS_UNBALANCE (4-1)
#define RSB_SUBDIVISION_BUG_EXTRA (4)		/* incorrect behaviour is encountered if setting this to 0, as it should (experienced on a 12 core machine). proper bugfix remains unknown to me. */

#define RSB_WANT_FASTER_EXPERIMENTAL_CONSTRUCTOR 0 /* FIXME: this is experimental code and shall be finished ! */
RSB_INTERNALS_COMMON_HEAD_DECLS

#define RSB__VERBOSE_REC2REC rsb_global_session_handle.verbose_tuning>0
#if RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE
#define RSB__VERBOSE_COO2REC RSB__VERBOSE_REC2REC || (rsb_global_session_handle.rsb_g_verbose_interface&2)
#else  /* RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE */
#define RSB__VERBOSE_COO2REC RSB__VERBOSE_REC2REC
#endif /* RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE */

#if RSB_WANT_FASTER_EXPERIMENTAL_CONSTRUCTOR
#define RSB_POW2(P) (1<<(P))
//#define RSB_MUL2(P) ((P) *= 2)
#define RSB_MUL2(P) ((P)<<=1)
//#define RSB_HALF(V) (((V)+1)/2)
#define RSB_HALF(V) (((V)+1)>>1)
#define RSB_POW4(P) (RSB_POW2(P)*RSB_POW2(P))
static inline rsb_nnz_idx_t rsb_coo_index_bit_interleave(rsb_coo_idx_t o, rsb_coo_idx_t e)
{
	/* FIXME: this is DUPLICATE code !!! */
	rsb_nnz_idx_t i = 0, O = o, E = e;
	RSB_DEBUG_ASSERT(O>=0);
	RSB_DEBUG_ASSERT(E>=0);
	if (sizeof(rsb_nnz_idx_t)==1)
	{
		E = (E | (E << 2)) & 0x33;
		E = (E | (E << 1)) & 0x55;
		O = (O | (O << 2)) & 0x33;
		O = (O | (O << 1)) & 0x55;
	}
	else
	if (sizeof(rsb_nnz_idx_t)==2)
	{
		E = (E | (E << 4)) & 0x0F0F;
		E = (E | (E << 2)) & 0x3333;
		E = (E | (E << 1)) & 0x5555;
		O = (O | (O << 4)) & 0x0F0F;
		O = (O | (O << 2)) & 0x3333;
		O = (O | (O << 1)) & 0x5555;
	}
	else
	if (sizeof(rsb_nnz_idx_t)==4)
	{
		E = (E | (E << 8)) & 0x00FF00FF;
		E = (E | (E << 4)) & 0x0F0F0F0F;
		E = (E | (E << 2)) & 0x33333333;
		E = (E | (E << 1)) & 0x55555555;
		O = (O | (O << 8)) & 0x00FF00FF;
		O = (O | (O << 4)) & 0x0F0F0F0F;
		O = (O | (O << 2)) & 0x33333333;
		O = (O | (O << 1)) & 0x55555555;
	}
	else
	if (sizeof(rsb_nnz_idx_t)==8)
	{
		E = (E | (E <<16)) & 0x0000FFFF0000FFFF;
		E = (E | (E << 8)) & 0x00FF00FF00FF00FF;
		E = (E | (E << 4)) & 0x0F0F0F0F0F0F0F0F;
		E = (E | (E << 2)) & 0x3333333333333333;
		E = (E | (E << 1)) & 0x5555555555555555;
		O = (O | (O <<16)) & 0x0000FFFF0000FFFF;
		O = (O | (O << 8)) & 0x00FF00FF00FF00FF;
		O = (O | (O << 4)) & 0x0F0F0F0F0F0F0F0F;
		O = (O | (O << 2)) & 0x3333333333333333;
		O = (O | (O << 1)) & 0x5555555555555555;
	}
	else
	{
		RSB_ERROR(RSB_ERRM_FYRYNS);
		/* FIXME : fatal! */
	}

	i = (E | (O << 1));
	RSB_DEBUG_ASSERT((i & ~-1)>=0);
	return i;
}

rsb_err_t rsb_assign_subm__(rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_coo_idx_t nlev = 4;
	/* FIXME: irm, jrm fit together in a single byte! (1byte/nnz!) */
	rsb_coo_idx_t nlevi,nlevj;
	rsb_coo_idx_t mink = RSB_POW2(nlev),maxk = mink;
	rsb_nnz_idx_t scnt[16*16], nzi;
	rsb_coo_idx_t idiv[   16];
	rsb_coo_idx_t jdiv[   16];
	rsb_nnz_idx_t nzoff = 0;
	rsb_coo_idx_t * zIA = NULL, * zJA = NULL;
	rsb_coo_idx_t * oIA = NULL, * oJA = NULL;

	RSB_BZERO_P(&scnt);
	oJA = rsb__malloc(nnz*sizeof(rsb_coo_idx_t));
	oIA = rsb__malloc(nnz*sizeof(rsb_coo_idx_t));
	zJA = rsb__malloc(nnz*sizeof(rsb_coo_idx_t));
	zIA = rsb__malloc(nnz*sizeof(rsb_coo_idx_t));
	if(!oJA || !oIA) goto err;
	if(!zJA || !zIA) goto err;
	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	{
		const rsb_coo_idx_t i = IA[nzi],j = JA[nzi];
		rsb_coo_idx_t irm = 0,jrm = 0;
		rsb_coo_idx_t sm = m,sk = k;
		rsb_coo_idx_t om = 0,ok = 0;
		rsb_int sidx;
		for(nlevi=0;nlevi<nlev;++nlevi)
		{
			rsb_coo_idx_t hm = RSB_HALF(sm),hk = RSB_HALF(sk);
			RSB_MUL2(irm);RSB_MUL2(jrm);
			if(i>=hm+om){irm += 1;sm-=hm;om += hm;}else{sm = hm;}
			if(j>=hk+ok){jrm += 1;sk-=hk;ok += hk;}else{sk = hk;}
		}
		zIA[nzi] = irm;
		zJA[nzi] = jrm;
		//sidx = 16*irm+jrm;
		sidx = rsb_coo_index_bit_interleave(irm,jrm);
		//scnt[sidx]++;
		//printf("hm:%d sm:%d\n",hm,sm); printf("hk:%d sk:%d\n",hk,sk);
		//printf("%d %d -> %d %d (%d)    at %d %d  sized %d %d\n",i,j,irm,jrm,sidx,om,ok,sm,sk);
	}
	if(0)
	for(nlevi=0;nlevi<nlev;++nlevi)
	for(nlevj=0;nlevj<nlev;++nlevj)
	{
		printf("%d %d : %d\n",nlevi,nlevj,scnt[16*nlevi+nlevj]);
	}
	for(nlevi=1;nlevi<nlev*nlev;++nlevi)
		scnt[nlev*nlev-nlevi] = scnt[nlev*nlev-nlevi-1];
	scnt[0] = 0;
	for(nlevi=1;nlevi<nlev*nlev;++nlevi)
		scnt[nlevi] += scnt[nlevi-1]; /* FIXME: shall use rsb__do_prefix_sum_coo_idx_t */
	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	{
		rsb_int sidx;
		rsb_coo_idx_t irm = 0,jrm = 0;
		irm = zIA[nzi];
	       	jrm = zJA[nzi];
		sidx = rsb_coo_index_bit_interleave(irm,jrm);
		//sidx = 0;
		oIA[ scnt[sidx]  ] = IA[nzi];
		oJA[ scnt[sidx]++] = JA[nzi];
	}
	rsb__memcpy(IA,oIA,sizeof(rsb_coo_idx_t)*nnz);
	rsb__memcpy(JA,oJA,sizeof(rsb_coo_idx_t)*nnz);
	//for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	//printf("please ignore this value: %d\n",scnt[0]);
err:
	RSB_CONDITIONAL_FREE(oJA);
	RSB_CONDITIONAL_FREE(oIA);
	RSB_CONDITIONAL_FREE(zJA);
	RSB_CONDITIONAL_FREE(zIA);
	return errval;
}

static void rsb_allocate_new__(void *RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_err_t *errvalp)
{
	rsb_time_t dt;
	long cs = rsb__get_first_level_c_size();
	//long cs = rsb__get_lastlevel_c_size();
	rsb_nnz_idx_t bnz = RSB_MIN(cs/(4*sizeof(rsb_coo_idx_t)),nnz),fnz;
	//rsb_nnz_idx_t bnz = nnz,fnz;
	if(!getenv("RSB_CB"))bnz = nnz;
	RSB_TIC(dt);
	for(fnz=0;fnz<nnz;fnz+=bnz)
	{
		rsb_assign_subm__(IA+fnz,JA+fnz,m,k,RSB_MIN(bnz,nnz-fnz));
	}
	RSB_TOC(dt);
	printf("%d cache blocks\n",(nnz+bnz-1)/bnz);
	printf("EXPERIMENTAL: processed indices at %lf Mnnz/s in %lf s\n",(RSB_FPINV(dt)*nnz)/RSB_MILLION_F,dt);
	exit(0);
}
#endif /* RSB_WANT_FASTER_EXPERIMENTAL_CONSTRUCTOR */

void rsb__do_set_in_place_submatrices_offsets(struct rsb_mtx_t *RSB_RESTRICT submatrices, rsb_submatrix_idx_t cmc, rsb_char_t *RSB_RESTRICT  VA, rsb_coo_idx_t *RSB_RESTRICT  IA, rsb_coo_idx_t *RSB_RESTRICT JA, size_t el_size)
{
	/**
		\ingroup gr_internals
		\note: if nnz==0 and diagonal implicit, this could be dangerous
	 */
	rsb_submatrix_idx_t smi;
	for(smi=0;smi<cmc;++smi)
	{
		struct rsb_mtx_t * submatrix = submatrices+smi;
		if(!RSB_DO_FLAG_HAS(submatrix->flags,RSB_FLAG_QUAD_PARTITIONING))
		{
			submatrix->bpntr = IA+submatrix->nzoff;
			submatrix->bindx = JA+submatrix->nzoff;
			submatrix->VA = ((rsb_char_t*)VA)+el_size*submatrix->nzoff;
		}
	}
}

rsb_err_t rsb__do_switch_recursive_matrix_to_fullword_storage(struct rsb_mtx_t * mtxAp)
{
	/**
		\ingroup gr_internals
		TODO: move somewhere else
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if(!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_E_MTXAP);
		return RSB_ERR_BADARGS;
	}

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
			if(submatrix)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_recursive_matrix_to_fullword_storage(submatrix));
	}
	else
	{
		if( mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR )
		{
		       	if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
			{
				rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*)(mtxAp->bpntr),mtxAp->nnz,0);
				rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*)(mtxAp->bindx),mtxAp->nnz,0);
			       	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES);
			}
		}
		else
		if( mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCSR )
		{
		       	if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
			{
				rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t*)(mtxAp->bindx),mtxAp->nnz,0);
			       	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES);
			}
		}
		else
			errval = RSB_ERR_BADARGS;
	}
	
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_switch_fresh_terminal_matrix_to_halfword_storages(struct rsb_mtx_t * mtxAp)
{
	/**
		\ingroup gr_unfinished
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_BADARGS;
	}
	else
	{
		if(RSB_C2R_IF_VERBOSE && 0)
			RSB_INFO_MATRIX_SUMMARY(mtxAp),
			RSB_INFO(" -switch.."),
			RSB_INFO("HCOO?(%d)..",rsb__do_is_candidate_for_halfword_coo(mtxAp)),
			RSB_INFO("HCSR?(%d)",rsb__do_is_candidate_for_halfword_csr(mtxAp)),
			RSB_INFO("\n");
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE))
		{
			if(rsb__do_is_candidate_for_halfword_coo(mtxAp))
			{	
				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("to halfword COO:"),RSB_INFO_MATRIX_SUMMARY(mtxAp),RSB_INFO("\n");
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_to_halfword_coo(mtxAp));
			}
			else
				RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES);
		}
		else
//		if(rsb__do_is_candidate_for_fullword_coo(mtxAp))
//		{
//			RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo(mtxAp,RSB_BOOL_FALSE));
//			RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_to_fullword_zcoo(mtxAp));
//		}
//		else
#if RSB_WANT_FIX_BUG_DISCOVERED_20121210
		if(!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE))
#else /* RSB_WANT_FIX_BUG_DISCOVERED_20121210 */
		if( RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE))
#endif /* RSB_WANT_FIX_BUG_DISCOVERED_20121210 */
		{
			if(RSB_SOME_ERROR(rsb__do_is_candidate_for_halfword_csr(mtxAp)))
			{
				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("to halfword CSR:"),RSB_INFO_MATRIX_SUMMARY(mtxAp),RSB_INFO("\n");
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_to_halfword_csr(mtxAp));
			}
			else
				RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES);
		}
		else
		//if(!rsb__is_root_matrix(mtxAp) || rsb__is_terminal_recursive_matrix(mtxAp)) /* root recursive or root nonrec. */
		//	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES);
		//else
			RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES);
		;/* for root matrices, we keep the flags, because some of the leaves MAY have it */
	}
	
	RSB_DO_ERR_RETURN(errval)
}

#if !RSB_WANT_MORE_PARALLELISM 
static rsb_err_t rsb_do_switch_fresh_recursive_matrix_to_halfword_storages(struct rsb_mtx_t * mtxAp)
{
	/**
		\ingroup gr_unfinished
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_E_MTXAP);
		return RSB_ERR_BADARGS;
	}

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
			if(submatrix)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_switch_fresh_recursive_matrix_to_halfword_storages(submatrix));
	}
	else
		errval = rsb_do_switch_fresh_terminal_matrix_to_halfword_storages(mtxAp);
	
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_MORE_PARALLELISM */

rsb_err_t rsb__check_bounds(struct rsb_mtx_t * mtxAp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(
			mtxAp->broff<0 || 
			mtxAp->bcoff<0 
			/* mtxAp->broff<mtxAp->roff+broff || mtxAp->bcoff<mtxAp->coff+bcoff ||  */
	  )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_ERROR(RSB_ERRM_BCE);
		RSB_ERROR(RSB_PRINTF_MTX_SUMMARIZED_ARGS(RSB_ERRM_BM,mtxAp,RSB_ERRM_NL));

		RSB_ERROR(RSB_PRINTF_MATRIX_BOUNDS_SUMMARY_ARGS(mtxAp)); RSB_ERROR("\n");

		RSB_ASSERT(! ( mtxAp->broff<0) );
		RSB_ASSERT(! ( mtxAp->bcoff<0) );
	/*	RSB_ASSERT(! ( mtxAp->broff<mtxAp->roff+broff) );
		RSB_ASSERT(! ( mtxAp->bcoff<mtxAp->coff+bcoff) ); */

		goto ret;
	}
	if(
			mtxAp->bm>mtxAp->nr || 
			mtxAp->bk>mtxAp->nc 
	  )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_ERROR(RSB_ERRM_BCE);
		RSB_ERROR(RSB_PRINTF_MTX_SUMMARIZED_ARGS(RSB_ERRM_BM,mtxAp,RSB_ERRM_NL));
		RSB_ERROR(RSB_PRINTF_MATRIX_BOUNDS_SUMMARY_ARGS(mtxAp)); RSB_ERROR("\n");

		RSB_ASSERT(! ( mtxAp->bm>mtxAp->nr) );
		RSB_ASSERT(! ( mtxAp->bk>mtxAp->nc) );

		goto ret;
	}
ret:
	return errval;
}

rsb_err_t rsb__compute_bounded_box(struct rsb_mtx_t * mtxAp)
{
	/*
	 * Set have to be: nr, nc.
	 * Will compute: ...
	 *
	 * TODO: make sure it does not depend on RSB_FLAG_QUAD_PARTITIONING.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t broff = RSB_INVALID_COO_IDX_VAL,bcoff = RSB_INVALID_COO_IDX_VAL,bm = RSB_INVALID_COO_IDX_VAL,bk = RSB_INVALID_COO_IDX_VAL;

	if(rsb__is_coo_matrix(mtxAp))
	{
		//rsb_nnz_idx_t nnz0 = 0,nnz1 = mtxAp->nnz;
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		{
			RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			rsb_half_idx_t li, ui;
			rsb_half_idx_t lj, uj;
			// TODO: optimize: IA is sorted
			rsb__util_find_extremal_half_index_val(IA,mtxAp->nnz,0,mtxAp->nr,&li,&ui);
			rsb__util_find_extremal_half_index_val(JA,mtxAp->nnz,0,mtxAp->nc,&lj,&uj);
			bk = 1;bk += uj; bm = 1;bm += ui; broff = li; bcoff = lj;
			RSB_DEBUG_ASSERT( mtxAp->nr>=bm && mtxAp->nc>=bk );
		}
		else
		{
			RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			rsb_coo_idx_t li, ui;
			rsb_coo_idx_t lj, uj;
			// TODO: optimize: IA is sorted
			rsb__util_find_extremal_full_index_val(IA,mtxAp->nnz,0,mtxAp->nr,&li,&ui);
			rsb__util_find_extremal_full_index_val(JA,mtxAp->nnz,0,mtxAp->nc,&lj,&uj);
			bk = 1;bk += uj; bm = 1;bm += ui; broff = li; bcoff = lj;
			RSB_DEBUG_ASSERT( mtxAp->nr>=bm && mtxAp->nc>=bk );
		}
	}
	else
	if(rsb__is_csr_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
		{
			RSB_DECLARE_CONST_HALFCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			rsb_half_idx_t lj, uj;
			rsb_coo_idx_t li, ui;
			rsb__util_find_extremal_half_index_val(JA,mtxAp->nnz,0,mtxAp->nr,&lj,&uj);
			ui = rsb__nnz_split_nnz_bsearch(PA,mtxAp->nnz,mtxAp->nr+1);
			li = rsb__nnz_split_nnz_bsearch(PA,1,mtxAp->nr+1)-1;
			bk = 1;bk += uj; bm = ui; broff = li; bcoff = lj;
			RSB_DEBUG_ASSERT( mtxAp->nr>=bm && mtxAp->nc>=bk );
		}
		else
		{
			RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			rsb_coo_idx_t lj, uj;
			rsb_coo_idx_t li, ui;
			rsb__util_find_extremal_full_index_val(JA,mtxAp->nnz,0,mtxAp->nr,&lj,&uj);
			ui = rsb__nnz_split_nnz_bsearch(PA,mtxAp->nnz,mtxAp->nr+1);
			li = rsb__nnz_split_nnz_bsearch(PA,1,mtxAp->nr+1)-1;
			bk = 1;bk += uj; bm = ui; broff = li; bcoff = lj;
			RSB_DEBUG_ASSERT( mtxAp->nr>=bm && mtxAp->nc>=bk );
		}
	}
	else
		RSB_ERROR(RSB_ERRMSG_BADFORMAT);

	mtxAp->broff = mtxAp->roff+broff;
	mtxAp->bcoff = mtxAp->coff+bcoff;
	mtxAp->bm = bm;
	mtxAp->bk = bk;

	errval = rsb__check_bounds(mtxAp);
#if 0
	RSB_INFO("bounding box of "),RSB_INFO_MATRIX_SUMMARY(mtxAp),
		RSB_INFO(": %.2f%% x  %.2f %%\n",(100.0f*(float)(bm-broff))/mtxAp->nr,(100.0f*(float)(bk-bcoff))/mtxAp->nc);
		RSB_INFO(": %d,%d %d,%d\n",mtxAp->roff+broff,mtxAp->coff+bcoff,bm,bk);
		RSB_INFO(": %d,%d %d,%d\n",mtxAp->roff+broff,mtxAp->coff+bcoff,bm,bk);
#endif
	return errval;
}

static rsb_err_t rsb_do_compute_bounded_boxes(struct rsb_mtx_t * mtxAp)
{
	/**
		\ingroup gr_unfinished
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t smi = 0;
	rsb_bool_t want_really = 0;
#if RSB_WANT_BOUNDED_BOXES
	want_really = (rsb_global_session_handle.want_bounded_box!=0);
#else /* RSB_WANT_BOUNDED_BOXES */
	mtxAp->broff = roff;
	mtxAp->bcoff = coff;
	mtxAp->bm = m;
	mtxAp->bk = k;
	goto err;
#endif /* RSB_WANT_BOUNDED_BOXES */

	if(!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_E_MTXAP);
		return RSB_ERR_BADARGS;
	}

	if(mtxAp->nnz==0)
		goto err;

	if(want_really)
	{
		if(rsb__is_terminal_recursive_matrix(mtxAp)) // fix for serial 20101206
		{
			RSB_DO_ERROR_CUMULATE(errval,rsb__compute_bounded_box(mtxAp));
			goto err;
		}
		#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC 
		for(smi=0;smi<mtxAp->all_leaf_matrices_n;++smi)
		{
			struct rsb_mtx_t * submatrix = mtxAp->all_leaf_matrices[smi].mtxlp;
			RSB_DO_ERROR_CUMULATE(errval,rsb__compute_bounded_box(submatrix));
		}
		#pragma omp barrier
	}
	else
	{
		#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC 
		for(smi=0;smi<mtxAp->all_leaf_matrices_n;++smi)
		{
			struct rsb_mtx_t * submatrix = mtxAp->all_leaf_matrices[smi].mtxlp;
			submatrix->bm = submatrix->nr;
			submatrix->bk = submatrix->nc;
			submatrix->broff = submatrix->roff;
			submatrix->bcoff = submatrix->coff;
		}
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_switch_fresh_recursive_matrix_to_halfword_storages_parallel(struct rsb_mtx_t * mtxAp)
{
	/**
		\ingroup gr_internals
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t smi = 0;

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}

	/* FIXME: 20100809 it seems that 'switching' a 0-nnz matrix overwrites something which should not be overwritten  */
	if(mtxAp->nnz==0)
		goto err;

	if(rsb__is_terminal_recursive_matrix(mtxAp)) // fix for serial 20101206
	{
		errval = rsb_do_switch_fresh_terminal_matrix_to_halfword_storages(mtxAp);
		if(RSB_SOME_ERROR(errval))
			RSB_ERROR("%s\n",rsb__get_errstr_ptr(errval));
		goto err;
	}
	#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(mtxAp) RSB_NTC 
	for(smi=0;smi<mtxAp->all_leaf_matrices_n;++smi)
	{
		struct rsb_mtx_t * submatrix = mtxAp->all_leaf_matrices[smi].mtxlp;
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_switch_fresh_terminal_matrix_to_halfword_storages(submatrix));
	}
	#pragma omp barrier
err:
	RSB_DO_ERR_RETURN(errval)
}

#if 0
static rsb_err_t rsb_do_shuffle_left_and_right_rows_inner(rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t m, rsb_coo_idx_t m0, rsb_nnz_idx_t nnz, rsb_nnz_idx_t nnz0, rsb_coo_idx_t * RSB_RESTRICT IL, rsb_coo_idx_t * IM, rsb_coo_idx_t * RSB_RESTRICT WA, size_t sz)
{
	/**
		\ingroup gr_unfinished
		FIXME: UNFINISHED, EXPERIMENTAL
	 */
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		rsb_coo_idx_t iu = m0,id = m-1;
		rsb_coo_idx_t wl = nnz,wr = 0,ns = 0,nu = 0,ws = 0,nd = nnz;

		if( sz<1 || !IA || !IM || !IL || !WA || RSB_INVALID_NNZ_INDEX(nnz) || RSB_INVALID_COO_INDEX(m) )
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
		}
		if(RSB_UNLIKELY(IL[m]!=nnz))
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
		if(iu>=id)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
		nu = IL[iu];

		while(RSB_LIKELY(iu<=id))
		{
			/* compute the left subrow length */
			ns = IM[iu]-IL[iu];
			/* shift left the left subrow */
			RSB_A_MEMMOVE(IA,IA,nu,IL[iu],ns,sz);
			/* update the counter of left subrows elements in IA */
			nu += ns;
			/* compute the right subrow length */
			ws = (IL[iu+1]-IM[iu]);
			/* buffer the right subrow */
			RSB_A_MEMCPY(WA,IA,wr,IM[iu],ws,sz);
			/* update the (complementary) counter of right subrows elements in the buffer */
			wr += ws;

			if(RSB_UNLIKELY(iu>=id))
			{
				/* skip row, as id was already done */
				++id;
				goto done;
			}
			/* compute the right subrow length */
			ns = IL[id+1]-IM[id];
			/* update the (complementary) counter of right subrows elements in IA */
			nd -= ns;
			/* shift right the right subrow */
			RSB_A_MEMMOVE(IA,IA,nd,IM[id],ns,sz);
			/* compute the left subrow length */
			ws = IM[id]-IL[id];
			/* update the counter of right subrows elements in the buffer */
			wl -= ws;
			/* buffer the left subrow */
			RSB_A_MEMCPY(WA,IA,wl,IL[id],ws,sz);

			++iu,--id;
		}
		/* IA has definitive elements, from left at  0..nu-1 and from right at (nnz-nd)..nnz-1  */
		{
			//rsb_nnz_idx_t
		}
		/* WA has definitive elements, from right at  0..wr-1 and from left at  wl..nnz-1  */
		/* it should be : nnz == nu + nnz-wl+nd-wr */
		if(nu+((nnz)-wl)!=nd-wr)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
done:
		/* compute the number of left submatrix elements in the buffer */
		ns = (nnz)-wl;
		/* copy the partial left submatrix from the buffer to the array */
		RSB_A_MEMMOVE(IA,WA,nu,wl,ns,sz);
		/* update the counter of left subrows elements in IA */
		nu += ns;
		/* compute the number of right submatrix elements in the buffer */
		ns = wr;
		/* copy the partial right submatrix from the buffer to the array */
		RSB_A_MEMMOVE(IA,WA,nu,0,ns,sz);
		/* update the counter to all subrows elements in IA (those already present, too) */
		nd -= ns;

		/* minimal coherence check */
err:
		if(RSB_UNLIKELY(nu!=nd))
		{
			RSB_ERROR("nnz=%d != nu+nd = %d; nu=%d, wl=%d, wr=%d, nd=%d\n",nnz,nu+nd,nu,wl,wr,nd);
//			RSB_ERROR("nnz=%d != nu+nnz-wl+nd = %d; nu=%d, wl=%d, wr=%d, nd=%d\n",nnz,nu+nnz+wl+nd,nu,wl,wr,nd);
			errval = RSB_ERR_INTERNAL_ERROR;
		}
		/* the buffer is empty now, and the arrays are left-right partitioned */
		RSB_DO_ERR_RETURN(errval)
}
#endif

#if 0
static rsb_err_t rsb_do_shuffle_left_and_right_rows(void * RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t m, rsb_coo_idx_t m0, rsb_nnz_idx_t nnz, rsb_nnz_idx_t nnz0, rsb_type_t typecode, rsb_coo_idx_t * RSB_RESTRICT IL, rsb_coo_idx_t * IM, rsb_coo_idx_t * RSB_RESTRICT WA)
{
	/**
		\ingroup gr_unfinished
		FIXME: UNFINISHED, EXPERIMENTAL
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const size_t sz = sizeof(rsb_coo_idx_t);
	if(!IL || !IM || !WA)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_shuffle_left_and_right_rows_inner(IA,m,m0,nnz,nnz0,IL,IM,WA,sz));
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_shuffle_left_and_right_rows_inner(JA,m,m0,nnz,nnz0,IL,IM,WA,sz));
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_shuffle_left_and_right_rows_inner(VA,m,m0,nnz,nnz0,IL,IM,WA,RSB_SIZEOF(typecode)));
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif

static rsb_err_t rsb_do_compute_vertical_split_search_only(
		const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA,
	       	rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t m, rsb_coo_idx_t k,
	       	rsb_coo_idx_t hm, rsb_coo_idx_t hk, rsb_nnz_idx_t nnz,
	       	const rsb_coo_idx_t * IB, rsb_nnz_idx_t *ulp, rsb_nnz_idx_t *urp, rsb_nnz_idx_t *llp, rsb_nnz_idx_t *lrp)
{
	/**
	\ingroup gr_unfinished
	

	 */
	rsb_nnz_idx_t ul = 0,ur = 0,ll = 0,lr = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t dnnz = 0/*,wdnnz = 0,rnnz = 0,hrnnz = 0*/;
	register rsb_coo_idx_t i;
	//rsb_nnz_idx_t nnz0 = 0;
	//rsb_coo_idx_t xroff = 0;


	if(nnz>m || 1)
	//if(nnz>m)
	{
	for(i = roff;RSB_LIKELY(i<roff+m);++i)
	{
		// offset of line i in the global line pointers array
		rsb_nnz_idx_t nnz0 = IB[i];
		// nnz1..nnz0 are the boundaries of line i
		rsb_nnz_idx_t nnz1 = IB[i+1];
		rsb_nnz_idx_t nnz2 = 0;
		// check
		RSB_C2R_ASSERT(nnz0>=IB[i]);
		RSB_C2R_ASSERT(nnz1<=IB[i+1]);
		// skip line if empty
		if(nnz1-nnz0<1)continue;
		// find first element of line i also in the submatrix
		nnz0 += rsb__nnz_split_coo_bsearch(JA+nnz0,coff,nnz1-nnz0);
		// skip line if empty in the submatrix
		if(nnz1-nnz0<1)continue;
		// find the length of the subrow i in the submatrix
		nnz1 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,coff+k,nnz1-nnz0);
		//check 
		RSB_C2R_ASSERT(JA[nnz0+0]>=coff);
		// skip line if empty in the submatrix
		if(nnz1-nnz0<1)continue;
		nnz2 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,coff+hk,nnz1-nnz0);
	       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
		RSB_C2R_ASSERT(JA[nnz0+0]>=coff);
		RSB_C2R_ASSERT(JA[nnz1-1]< coff+k);
		dnnz += nnz1-nnz0;
		if(i<roff+hm)
			ul += nnz2-nnz0,
			ur += nnz1-nnz2;
		else
			ll += nnz2-nnz0,
			lr += nnz1-nnz2;
	}
	}
	else
	{
		// FIXME: UNFINISHED
		rsb_nnz_idx_t nnz0,nnz1,n;
		//RSB_INFO("almost empty matrix !\n");
		for(n=0;n<nnz;++n)
		{
			rsb_nnz_idx_t nnz2 = 0;
			i = IA[n];
			nnz0 = IB[i];
			nnz1 = IB[i+1];
			// ...
#if 1
		// skip line if empty
		if(nnz1-nnz0<1)continue;
		// find first element of line i also in the submatrix
		nnz0 += rsb__nnz_split_coo_bsearch(JA+nnz0,coff,nnz1-nnz0);
		// skip line if empty in the submatrix
		if(nnz1-nnz0<1)continue;
		// find the length of the subrow i in the submatrix
		nnz1 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,coff+k,nnz1-nnz0);
		//check 
		RSB_C2R_ASSERT(JA[nnz0+0]>=coff);
		// skip line if empty in the submatrix
		if(nnz1-nnz0<1)continue;
		nnz2 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,coff+hk,nnz1-nnz0);
	       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
		RSB_C2R_ASSERT(JA[nnz0+0]>=coff);
		RSB_C2R_ASSERT(JA[nnz1-1]< coff+k);
		dnnz += nnz1-nnz0;
		if(i<roff+hm)
			ul += nnz2-nnz0,
			ur += nnz1-nnz2;
		else
			ll += nnz2-nnz0,
			lr += nnz1-nnz2;
#else
			if(nnz1-nnz0<1)continue;
			if(i<roff+hm)
			{
				for(;n<nnz1;++n)
					if(JA[n]>=coff+hk)
						++ur;
					else
						++ul;
			}
			else
			{
				for(;n<nnz1;++n)
					if(JA[n]>=coff+hk)
						++lr;
					else
						++ll;
			}
#endif
		}
	}
//done:
	*llp = ll;
	*lrp = lr;
	*ulp = ul;
	*urp = ur;
//err:
	RSB_DO_ERR_RETURN(errval)
}

#if RSB_WANT_OMP_RECURSIVE_KERNELS
#if !RSB_WANT_PARALLEL_SUBDIVISION 
static rsb_err_t rsb_do_compute_vertical_split(const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_coo_idx_t hm, rsb_coo_idx_t hk, rsb_nnz_idx_t nnz, rsb_coo_idx_t * IL, rsb_coo_idx_t * RSB_RESTRICT IM, rsb_coo_idx_t * IR, rsb_nnz_idx_t *ulp, rsb_nnz_idx_t *urp, rsb_nnz_idx_t *llp, rsb_nnz_idx_t *lrp)
{
	/**
	\ingroup gr_unfinished
	FIXME: UNFINISHED, EXPERIMENTAL

	Computes two arrays: IM, IL.
	IM[i], contains the index of the first element >= hk on line i,
	IL[i], contains the index of the first element on line i.
	Notes: 
       		IM[i]==IL[i+1] if no element >=hk exists
       		IM[i]==IL[i]   if no IL[i] >=hk
       		IL[i]==IL[i+1]  if line i is empty
       		IL[0]==0
       		IL[m]==nnz
		IM is valid on the 0..nr-1 range.
		IL is valid on the 0..nr range.

	TODO: blocking support
	 */
	rsb_nnz_idx_t ul = 0,ur = 0,ll = 0,lr = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	rsb_nnz_idx_t wdnnz = 0;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	rsb_nnz_idx_t dnnz = 0, rnnz = 0, hrnnz = 0;
	register rsb_coo_idx_t i;
	rsb_nnz_idx_t nnz0 = 0;

	hk += coff;

	if(IR==NULL && IM==NULL && IL==NULL)
	{
		/* FIXME: write me, for cases where we want subdivision on barely-COO matrices */
	}
	else
	if(IR==NULL && IM==NULL)
	{
		// root matrix; should compute IL
		IL[0] = 0;
//		if(nnz>100*m)// TODO
		if(0)// TODO: determine which case is faster !
		{
			/* fill the row pointers array */
//			#pragma omp parallel for reduction(+:dnnz) RSB_NTC 
//			for(i=0;i<m;++i)
			for(i=0;RSB_LIKELY(i<m);++i)
			{
				// delimit the current row
				rnnz = rsb__nnz_split_coo_bsearch(IA+dnnz,i+1,nnz-dnnz);
				/* i==m-1 || IA[dnnz+rnnz] > i */
				dnnz += rnnz;
				IL[i+1] = dnnz;
				RSB_C2R_ASSERT(rnnz>=0);
			}
		}
		else
		{
			/* for full matrices, this is faster */
			rsb_nnz_idx_t n = 0;
#if 1
			nnz0 = 0;
#if RSB_WANT_QUADRANT_QUICK_DETECT 
			/* FIXME: UNFINISHED */
			if(IB[roff]==IB[roff+hm])
			{
				RSB_INFO("upper submatrix empty\n");
				// should do something sharp
			}
			else
			if(IB[roff+hm]==IB[roff+m])
			{
				RSB_INFO("lower submatrix empty\n");
				// should do something sharp
			}
#endif /* RSB_WANT_QUADRANT_QUICK_DETECT  */
			for(i=0;RSB_LIKELY(i<m);++i)
			{
				rnnz = 0;
				for(;RSB_LIKELY(n<nnz && IA[n]==i);++n)
					++rnnz;
				IL[i+1] = nnz0+rnnz;
				nnz0 += rnnz;
			}
#else
			for(i=0;RSB_LIKELY(i<m);++i)
				IL[i] = 0;
			for(n=0;RSB_LIKELY(n<nnz);++n)
				RSB_C2R_ASSERT(IA[n]>=0 && IA[n]<m);
			for(n=0;RSB_LIKELY(n<nnz);++n)
				IL[IA[n]+1]++;
			for(i=0;RSB_LIKELY(i<m);++i)
				IL[i+1] += IL[i];
#endif
		}
		RSB_C2R_ASSERT(IL[m]==nnz);
		goto err;
	}
	else
	if(IR==NULL)
	{
		RSB_C2R_ASSERT(0);
		RSB_ASSERT(ulp);
		RSB_ASSERT(llp);
		RSB_ASSERT(urp);
		RSB_ASSERT(lrp);
		// root matrix; should compute IL
		IL[0] = 0;
		/* fill the row pointers array */
		for(i=0;RSB_LIKELY(i<m);++i)
		{
			// delimit the current row
			rnnz = rsb__nnz_split_coo_bsearch(IA+dnnz,i+1,nnz-dnnz);
			/* i==m-1 || IA[dnnz+rnnz] > i */

			IL[i+1] = dnnz+rnnz;

			if(RSB_LIKELY(dnnz+rnnz<=nnz))
			{
				// the current row is non empty
				hrnnz = rsb__nnz_split_coo_bsearch(JA+dnnz,hk,rnnz);
				if(RSB_LIKELY(hrnnz<rnnz))
					IM[i] = dnnz+hrnnz;
				else
					// all the elements are in the left submatrix
					IM[i] = dnnz+ rnnz;
			}
			else
			{
				// last row
				hrnnz = rnnz,
				IM[i] = nnz;
			}

				if(RSB_UNLIKELY(IM[i]<IL[i]))
				{
					errval = RSB_ERR_INTERNAL_ERROR;
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
				}

			// TODO: split in two cycles: 0..hm-1, hm..nr-1
			if(i<hm)
				ul += IM[i  ]-IL[i],
				ur += IL[i+1]-IM[i];
			else
				ll += IM[i  ]-IL[i],
				lr += IL[i+1]-IM[i];
			//RSB_INFO("%d @ %d~%d (%d/%d)\n",i,dnnz,dnnz+rnnz-1,hrnnz,rnnz);
			dnnz += rnnz;
		}
		IM[m] = IL[m];
	}
	else
	{
		// compute middle pointers array, using the left and right ones
		RSB_ASSERT(ulp);
		RSB_ASSERT(llp);
		RSB_ASSERT(urp);
		RSB_ASSERT(lrp);
		nnz0 = IL[0];
		/* fill the middle row pointers array */
		for(i=0;RSB_LIKELY(i<m);++i)
		{
			// delimit the current row
			rsb_nnz_idx_t il = IL[i],ir = IR[i],im;
			rnnz = ir-il;

			RSB_C2R_ASSERT(ir>=il);
			if(ir<il)
			{
				errval = RSB_ERR_INTERNAL_ERROR;
				RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}
			if(ir==il)
			{
				// empty row
				IM[i] = IR[i];
				continue;
			}
			/* i==m-1 || IA[dnnz+rnnz] > i */

			// the current row is non empty
			RSB_C2R_ASSERT(JA[il+0]>=coff  );
			RSB_C2R_ASSERT(JA[ir-1]< coff+k);

			hrnnz = rsb__nnz_split_coo_bsearch(JA+il,hk,rnnz);
			im = il+hrnnz;

			IM[i] = im;

#if RSB_C2R_PARANOIA
			if(IM[i]>IR[i])
			{
				errval = RSB_ERR_INTERNAL_ERROR;
				RSB_PERR_GOTO(err,"i=%d, %d > %d!\n",i,IM[i],IR[i]);
			}

			if(IM[i]<IL[i])
			{
				errval = RSB_ERR_INTERNAL_ERROR;
				RSB_PERR_GOTO(err,"i=%d, %d < %d!\n",i,IM[i],IL[i]);
			}

#endif /* RSB_C2R_PARANOIA */
			// TODO: split in two cycles: 0..hm-1, hm..nr-1
			if(i<hm)
				ul += im-il,
				ur += ir-im;
			else
				ll += im-il,
				lr += ir-im;
			//RSB_INFO("%d @ %d~%d (%d/%d)\n",i,dnnz,dnnz+rnnz-1,hrnnz,rnnz);
			dnnz += rnnz;
		}
		IM[m] = IL[m];
	}

		if(RSB_C2R_PARANOIA)
		{
			rsb_coo_idx_t i;
			rsb_nnz_idx_t lnz = 0,rnz = 0,tnz = 0;
			if(IR==NULL)
				IR = IL+1;

/*			if(IL[m]!=nnz0+nnz)
			{
				RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
				RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}*/

			for(i=0;RSB_LIKELY(i<m);++i)
			{
				lnz += IM[i]-IL[i];
				rnz += IR[i]-IM[i];
				tnz += IR[i]-IL[i];

				if(RSB_UNLIKELY(IM[i]<IL[i] || IL[i]>IL[i+1]))
				{
					RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
				}
			}
			if(ul+ll!=lnz || ur+lr != rnz)
			{
				RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
				RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}
			if(tnz!=nnz || (rnz+lnz)!=nnz)
			{
				RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
				RSB_PERR_GOTO(err,"tnz:%d, nnz:%d, rnz:%d, lnz:%d\n",tnz,nnz,rnz,lnz);
			}
		}
		*llp = ll;
		*lrp = lr;
		*ulp = ul;
		*urp = ur;
err:
		RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_PARALLEL_SUBDIVISION */
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */

static rsb_err_t rsb_do_compute_vertical_split_parallel(const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_coo_idx_t hm, rsb_coo_idx_t hk, rsb_nnz_idx_t nnz, rsb_coo_idx_t * IL, rsb_coo_idx_t * RSB_RESTRICT IM, rsb_coo_idx_t * IR, rsb_nnz_idx_t *ulp, rsb_nnz_idx_t *urp, rsb_nnz_idx_t *llp, rsb_nnz_idx_t *lrp)
{

	/**
	Binary search for the boundaries of each row.
	Assign threads to rows intervals.
	Perform the row pointers vector fill calling rsb_do_compute_vertical_split.
	*/
#if 0
	return rsb_do_compute_vertical_split(IA,JA,roff,coff,m,k,hm,hk,nnz,IL,IM,IR,ulp,urp,llp,lrp);
#else
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* const rsb_thread_t wet = rsb_get_num_threads(); */

	if(m<1)
		return RSB_ERR_NO_ERROR;/* TODO: limit case */
	IL[0] = 0;
	if(m==1)
		goto after;

	#pragma omp parallel RSB_NTC 
	{
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		const rsb_thread_t tn = /*wet*/ omp_get_num_threads(), tnn = RSB_MIN(tn,m);
		const rsb_thread_t th_id = omp_get_thread_num();
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		const rsb_thread_t tn = 1, tnn = RSB_MIN(tn,m);
		const rsb_thread_t th_id = 0;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		const rsb_coo_idx_t mm = ((m+tnn-1)/tnn),m0 = mm*th_id,m1 = RSB_MIN(m0+mm,m);
		rsb_coo_idx_t i;
		rsb_nnz_idx_t nnz0 = 0,nnz1 = nnz;
		rsb_nnz_idx_t n,rnnz,fnnz,lnnz;

		if(th_id>=m)
			goto nowork;
		/* binary search for the boundaries of each row  */
		nnz0 = rsb__nnz_split_coo_bsearch(IA+nnz0,m0,nnz1-nnz0);
		nnz1 = nnz0+rsb__nnz_split_coo_bsearch(IA+nnz0,m1,nnz1-nnz0);
		/* assign threads to rows intervals */
		if(nnz0>=nnz1)
		{
			//for(i=m0;RSB_LIKELY(i<m1);++i)
			//	IL[i+1] = nnz0;
			RSB_XCOO_VSET(IL,nnz0,m0+1,m1+1);
			goto nowork;
		}
		//RSB_INFO("thread %d  rows %d..%d  nnz %d..%d\n",th_id,m0,m1,nnz0,nnz1);
		/* perform the row pointers vector fill calling rsb_do_compute_vertical_split */
		//RSB_DO_ERROR_CUMULATE(errval,rsb_do_compute_vertical_split(IA+nnz0,JA+nnz0,roff+m0,coff,m1-m0,k,hm,hk,nnz1-nnz0,IL+m0,NULL,NULL,ulp,urp,llp,lrp));
		fnnz = nnz0;
		n = nnz0;
		for(i=m0;RSB_LIKELY(i<m1);++i)
		{
			if((nnz1-nnz0)/(m1-m0)<RSB_WANT_BINSEARCH_MIN_NZPR) 
			{
				rnnz = 0;
				for(;RSB_LIKELY(n<nnz1 && IA[n]==i);++n)
					++rnnz;
				IL[i+1] = nnz0+rnnz;
				nnz0 += rnnz;
			}
			else
			{
				/* TODO : should use a smarter strategy than this one */
				lnnz = fnnz+rsb__nnz_split_coo_bsearch(IA+fnnz,i+1,nnz1-fnnz);
				//RSB_INFO("%d : %d\n",i,lnnz);
				IL[i+1] = lnnz;
				fnnz = lnnz;
			}
		}
nowork:			
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
	#pragma omp barrier
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
	}
after:
	IL[m] = nnz;
	//int i; RSB_INFO(":::"); for(i=0;RSB_LIKELY(i<m+1);++i) RSB_INFO("%d ",IL[i]); RSB_INFO("\n");
	//RSB_INFO(":::"); for(i=0;RSB_LIKELY(i<nnz);++i) RSB_INFO("%d ",IA[i]); RSB_INFO("\n");
//err:
	RSB_DO_ERR_RETURN(errval)
#endif
}

#if 0
static rsb_err_t rsb_do_fill_partially_rcsr_arrays_for_later(struct rsb_mtx_t * mtxAp, 
		const rsb_coo_idx_t * IL, const rsb_coo_idx_t * IR,
		rsb_coo_idx_t * IA, rsb_coo_idx_t * JA,
		//const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
		rsb_nnz_idx_t nzoff, rsb_coo_idx_t m, rsb_coo_idx_t roff )
{
	/**
		\ingroup gr_unfinished
	 */
	mtxAp->nzoff = nzoff;
	mtxAp->bindx = IA+nzoff;
	mtxAp->bpntr = NULL;
#if RSB_WANT_FIRST_VERSION
	RSB_COA_MEMMOVE(mtxAp->bindx,IL,0,roff,m+1);
#endif /* RSB_WANT_FIRST_VERSION */
#if RSB_WANT_FIRST_VERSION
{
	rsb_nnz_idx_t i;
	for(i=mtxAp->roff;RSB_LIKELY(i<mtxAp->roff+mtxAp->nr);++i)
	{
		rsb_nnz_idx_t nnz1 = IR[i];
		rsb_nnz_idx_t nnz0 = IL[i];
//		RSB_C2R_ASSERT(IL[i-mtxAp->roff]>=IB[i]);
//		RSB_C2R_ASSERT(IL[i-mtxAp->roff]<=IB[i+1]);
//		RSB_C2R_ASSERT(IL[i-mtxAp->roff+1]>=IB[i+1]);
		RSB_C2R_ASSERT(nnz0>=IL[i]);
		RSB_C2R_ASSERT(nnz1<=IR[i]);
		if(nnz1==nnz0)continue;
//	       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
		RSB_C2R_ASSERT(JA[nnz0+0]>=mtxAp->coff);
		RSB_C2R_ASSERT(JA[nnz1-1]< mtxAp->coff+mtxAp->nc);
	}
	}
#endif /* RSB_WANT_FIRST_VERSION */



	return RSB_ERR_NO_ERROR;
}
#endif

#if 0
static rsb_err_t rsb_do_fill_rcsr_arrays_for_later(struct rsb_mtx_t * mtxAp, 
		const rsb_coo_idx_t * IL, const rsb_coo_idx_t * IR,
	       	//const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	       	rsb_coo_idx_t * IA, rsb_coo_idx_t * JA,
		rsb_nnz_idx_t nzoff, rsb_coo_idx_t m, rsb_coo_idx_t roff )
{
	/**
		\ingroup gr_unfinished
	 */
	if(!IR)
		IR = IL+1;
	mtxAp->nzoff = nzoff;
	mtxAp->bindx = IA+nzoff;
	mtxAp->bpntr = IA+nzoff+m+1;
#if RSB_C2R_WANT_MAYBE_FASTER 
	RSB_COA_MEMCPY_ROWSZ(mtxAp->bindx,IL,0,roff,m+1);
	RSB_COA_MEMCPY_ROWSZ(mtxAp->bpntr,IR,0,roff,m+1);
//	RSB_COA_MEMCPY(mtxAp->bindx,IL,0,roff,m+1);
//	RSB_COA_MEMCPY(mtxAp->bpntr,IR,0,roff,m+1);
//	RSB_COA_MEMCPY_parallel(mtxAp->bindx,IL,0,roff,m+1);
//	RSB_COA_MEMCPY_parallel(mtxAp->bpntr,IR,0,roff,m+1);
#else /* RSB_C2R_WANT_MAYBE_FASTER */
	/* are we sure  we need MEMMOVE here ? FIXME */
	RSB_COA_MEMMOVE(mtxAp->bindx,IL,0,roff,m+1);
	RSB_COA_MEMMOVE(mtxAp->bpntr,IR,0,roff,m+1);
#endif /* RSB_C2R_WANT_MAYBE_FASTER */

#if RSB_C2R_PARANOIA
	{
	rsb_nnz_idx_t i;
	for(i=0;i<mtxAp->nr;++i)
	//for(i=mtxAp->roff;i<mtxAp->roff+mtxAp->nr;++i)
	{
		rsb_nnz_idx_t nnz1 = IR[i];
		rsb_nnz_idx_t nnz0 = IL[i];
//		RSB_C2R_ASSERT(IL[i-mtxAp->roff]>=IB[i]);
//		RSB_C2R_ASSERT(IL[i-mtxAp->roff]<=IB[i+1]);
//		RSB_C2R_ASSERT(IL[i-mtxAp->roff+1]>=IB[i+1]);
		RSB_C2R_ASSERT(nnz0>=IL[i]);
		RSB_C2R_ASSERT(nnz1<=IR[i]);
		if(nnz1==nnz0)continue;
//	       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
		RSB_C2R_ASSERT(JA[nnz0+0]>=mtxAp->coff);
		RSB_C2R_ASSERT(JA[nnz1-1]< mtxAp->coff+mtxAp->nc);
	}
	}
#endif /* RSB_C2R_PARANOIA */
	return RSB_ERR_NO_ERROR;
}
#endif

static rsb_err_t rsb_do_fill_early_leaf_matrix( struct rsb_mtx_t * mtxAp, struct rsb_mtx_t * submatrix,
		const rsb_coo_idx_t * IL, const rsb_coo_idx_t * IR, 
		rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, const rsb_coo_idx_t * VA,
		//const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_coo_idx_t * VA,
		rsb_nnz_idx_t snzoff, rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k,
		rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_type_t typecode, rsb_flags_t flags )
{
	/**
		\ingroup gr_unfinished
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

//	if(!IR)
//		IR = IL+1;

#if 0
	/* 20131206 nowadays IR and IL are always NULL */
	if(!RSB_DO_TOOFEWNNZFORRCSR(nnz,m) && IR && IL)
	{
		// the matrix could be split further: we fill it with info to continue, if necessary
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_rcsr_arrays_for_later(submatrix,IL,IR,IA,JA,snzoff,m,roff));
	}
	else
	if(!RSB_DO_TOOFEWNNZFORCSR(nnz,m) && IR && IL)
	{
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_partially_rcsr_arrays_for_later(submatrix,IL,IR,IA,JA,snzoff,m,roff));
		//RSB_ERROR("nnz=%d ,m=%d ! what shall we do ?\n",nnz,m);
		RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_WANT_BCSS_STORAGE);
	}
	else
#endif
	{
		if(RSB_C2R_IF_VERBOSE)
			RSB_INFO("building a very sparse recursive matrix\n");

		/* no hope for CSR : however, full/half word COO will fit  */
		submatrix->nzoff = snzoff;
		submatrix->bindx = NULL;
		submatrix->bpntr = NULL;
		RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_WANT_BCSS_STORAGE);
		//RSB_ERROR("nnz=%d ,m=%d ! what shall we do ?\n",nnz,m);
	}
	mtxAp->sm[roff?(coff?3:2):(coff?1:0)] = submatrix;
	submatrix->roff = roff+mtxAp->roff;
	submatrix->coff = coff+mtxAp->coff;
	RSB_DO_ERROR_CUMULATE(errval,rsb__set_init_flags_and_stuff(submatrix,NULL,NULL,m,k,nnz,nnz,nnz,typecode,flags));
//err:
	RSB_DO_ERR_RETURN(errval)
}

static int rsb__compar_rcsr_matrix_leftmost_first(const void * ap, const void * bp)
{
	/**
		\ingroup gr_internals
		Compare function to be used with qsort.
		Non-recursive matrices come after recursive ones.
		If both are or none is recursive, then order first by column, then by row.
	*/
	struct rsb_mtx_t *a = *(struct rsb_mtx_t **)ap;
	struct rsb_mtx_t *b = *(struct rsb_mtx_t **)bp;
	const rsb_bool_t at = !RSB_DO_FLAG_HAS(a->flags,RSB_FLAG_QUAD_PARTITIONING);
	const rsb_bool_t bt = !RSB_DO_FLAG_HAS(b->flags,RSB_FLAG_QUAD_PARTITIONING);
	int ss = 1;	/* should swap results ? */

	if(at && !bt)
		return 1;
	if(!at && bt)
		return -1;

	if(a->coff < b->coff)
	{
		RSB_SWAP(struct rsb_mtx_t *,a,b);
		ss = -1;/* should swap results ! */
	}

	return (a->coff==b->coff)?(a->roff>b->roff?1:(a->roff<b->roff?-1:0)):1*ss;
}

#if 0
/* 20121001 unfinished code: commented */
static struct rsb_mtx_t * rsb_do_find_ffmltart(struct rsb_mtx_t ** submatricesp, rsb_submatrix_idx_t smn, struct rsb_mtx_t * submatrix, rsb_coo_idx_t off)
{
	/**
		\ingroup gr_unfinished
	*/
	rsb_submatrix_idx_t smi = 0;
	rsb_coo_idx_t coff = submatrix->coff;
	rsb_coo_idx_t roff = submatrix->roff;
	rsb_coo_idx_t m = submatrix->nr;
	rsb_coo_idx_t k = submatrix->nc;
	/* leftmost from right */
	for(smi=0;smi<smn;++smi)
		if(submatricesp[smi]->coff>=coff+k &&
			       	submatricesp[smi]->roff<=off+0 && submatricesp[smi]->roff+submatricesp[smi]->nr>off)
			return submatricesp[smi];
	/* leftmost from left, the line after */
	for(smi=0;smi<smn;++smi)
		if(submatricesp[smi]->coff<coff &&
			       	submatricesp[smi]->roff<=off+1 && submatricesp[smi]->roff+submatricesp[smi]->nr>off)
			return submatricesp[smi];
	return NULL;
}
#endif

static rsb_bool_t rsb__should_recursively_partition_matrix(
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_nnz_idx_t element_count,
	rsb_nnz_idx_t block_count,
	rsb_nnz_idx_t nnz,
	rsb_blk_idx_t Mdim,
	rsb_blk_idx_t mdim,
	rsb_coo_idx_t roff,
	rsb_coo_idx_t coff,
	rsb_flags_t flags,
	size_t el_size,
	rsb_thread_t wet
)
{
#if (RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY == RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE)
	/*
		\ingroup gr_internals
		NEW: UNFINISHED
		TODO : ERROR HANDLING, DOCS
		TODO : should partition on high nnz per row count
			although in this case if the matrix is highly 
			compact (almost banded) it won't help.
			Therefore, some sparseness statistics in the matrix constructor
			would be nice.
	*/
	//long cs = rsb__get_lastlevel_c_size();
	//long cs = rsb__get_lastlevel_c_size_per_thread();
	long cs=(rsb__get_lastlevel_c_size()/(wet>0?wet:rsb_get_num_threads()));
	rsb_bool_t sp = RSB_BOOL_FALSE;	/* should partition */
	rsb_fillin_t efillin=1.0;		/* FIXME */
	size_t smab=0;			/* spmv memory accessed bytes */
	//cs/=20000;
	//cs/=20;
	//cs/=4;
	/* FIXME */
	if(nnz<RSB_RECURSION_MIN_NNZ  || m<RSB_RECURSION_MIN_DIM  || k<RSB_RECURSION_MIN_DIM  || !RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
	{
		sp = RSB_BOOL_FALSE;
		goto done;
	}

	if(kB<1)kB=1;
	if(mB<1)mB=1;

#if 0
	if(flags & RSB_FLAG_RECURSIVE_DOUBLE_DETECTED_CACHE) 
		cs*=2;

	if(flags & RSB_FLAG_RECURSIVE_HALF_DETECTED_CACHE) 
		cs/=2;
#endif

	if( (flags & RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG) && roff == coff )
		cs/=2;
	/* this will imply a more fine grained subdivision on the diagonal 
	 * FIXME : we could use a factor different than 2 !
	 * */


	/* subdivide at least until matrix indices can be compressed */
	if((flags & RSB_FLAG_USE_HALFWORD_INDICES_CSR) && m>1 && k>1 && 
//			nnz>1
//			nnz>(cs/4)
			nnz*el_size>2*cs
			&& !rsb__do_is_candidate_size_for_halfword_csr(m,k,nnz,flags))
		return RSB_BOOL_TRUE;
	if((flags & RSB_FLAG_USE_HALFWORD_INDICES_COO) && m>1 && k>1 && 
//			nnz>1
//			nnz>(cs/4)
			nnz*el_size>2*cs
			&& !rsb__do_is_candidate_size_for_halfword_coo(m,k,flags))
		return RSB_BOOL_TRUE;

	if(cs>0)
	{
		smab = rsb_spmv_memory_accessed_bytes_(mB,kB,m,k,efillin*nnz,((efillin*nnz)/mB)/kB,m/mB,el_size);

		if( 2*smab > 3*3*cs )	/* FIXME : overflow possible */
			sp=1;
		else
		if( 
			/* FIXME! */
			(((Mdim+mdim+m+k)*sizeof(rsb_coo_idx_t))
			/(nnz*el_size)) > 8*cs
		)
			sp = RSB_BOOL_TRUE;
		else
			sp = RSB_BOOL_FALSE;
	}
	else
	{	
		/* no cache info (FIXME: there should be no section like this one) */
		if(  
			Mdim<8 || mdim<8 || m < 500 
			|| k < 500 || nnz < 200*100)
			sp = RSB_BOOL_FALSE;
		else
			sp = RSB_BOOL_TRUE;
	}
#ifdef RSB_EXPERIMENTAL_ROWS_SUBDIVIDE_TO_CORES_NUM
	/* STILL UNIMPLEMENTED */
#endif /* RSB_EXPERIMENTAL_ROWS_SUBDIVIDE_TO_CORES_NUM */
#if RSB_EXPERIMENTAL_NO_SUBDIVIDE_ON_MIN_NNZ_PER_ROW_OR_COLUMN
	if(1)
	{
		rsb_nnz_idx_t nnzpr;
		if( (flags&RSB_FLAG_WANT_COLUMN_MAJOR_ORDER) != 0 )
			nnzpr=nnz/k;
		else
			nnzpr=nnz/m;

		if( nnzpr < RSB_CONST_MIN_NNZ_PER_ROW_OR_COLUMN_PER_SUBMATRIX )
			sp = RSB_BOOL_FALSE;
	}
#endif /* RSB_EXPERIMENTAL_NO_SUBDIVIDE_ON_MIN_NNZ_PER_ROW_OR_COLUMN */
done:
	return sp;
#else /* (RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY == RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE) */
	#error "should use a RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE partitioning policy!"
	return RSB_BOOL_FALSE;
#endif /* (RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY == RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE) */
}

#if RSB_OBSOLETE_QUARANTINE
static rsb_nnz_idx_t rsb_do_copy_submatrix_coa(struct rsb_mtx_t * submatrix, void * VA, void * WA, rsb_coo_idx_t * IL, rsb_coo_idx_t * IR, size_t el_size, rsb_coo_idx_t n0, rsb_coo_idx_t i0, rsb_coo_idx_t m0)
{
	rsb_nnz_idx_t i,n;
	RSB_C2R_ASSERT(submatrix);
	RSB_C2R_ASSERT(VA);
	RSB_C2R_ASSERT(WA);
	RSB_C2R_ASSERT(IL);
	RSB_C2R_ASSERT(IR);
	RSB_C2R_ASSERT(el_size>0);
	RSB_C2R_ASSERT(n0>=0);
	RSB_C2R_ASSERT(i0<m0);
	submatrix->VA=((char*)VA)+el_size*submatrix->nzoff;
	for(n=n0,i=i0;RSB_LIKELY(i<m0);n+=IR[i]-IL[i],++i)
	{
		RSB_C2R_ASSERT(n>=0);
		RSB_A_MEMCPY_SMALL(WA,VA,submatrix->nzoff+n,IL[i],IR[i]-IL[i],el_size);
	}
	return n;
}
#endif /* RSB_OBSOLETE_QUARANTINE */


static rsb_submatrix_idx_t rsb_do_pick_largest_open_matrix(struct rsb_mtx_t ** submatricesp, rsb_submatrix_idx_t smc)
{
	rsb_submatrix_idx_t smi = 0,msmi = RSB_SUBM_IDX_MARKER;
	rsb_nnz_idx_t maxnz = 0;
	if(RSB_WANT_VERBOSE_SUBDIVISION)
		if(smc==0)
			RSB_INFO("warning: no largest open matrix among 0 matrices\n");
	for(smi=0;smi<smc;++smi)
	{
		//RSB_INFO("looking %d : %d\n",smi,submatricesp[smi]->nnz);
		/* NOTE: ">=" down here is used to cope with diagonal implicit matrices (which could have nnz==0), too */
		if(submatricesp[smi])
		if(submatricesp[smi]->nnz>=maxnz)
		{
			maxnz = submatricesp[smi]->nnz;
			msmi = smi;
		}
	}
	return msmi;
}

static rsb_thread_t rsb_get_num_coo2rec_threads(void)
{
	return
#ifdef RSB_WANT_COO2RSB_THREADS
		rsb_global_session_handle.coo2rsb_threads ? rsb_global_session_handle.coo2rsb_threads :
#endif /* RSB_WANT_COO2RSB_THREADS */
			rsb_get_num_threads();
}

static rsb_err_t rsb_do_coo2rec_subdivide_parallel(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_err_t *errvalp, struct rsb_mtx_t ** submatricesp, struct rsb_mtx_t * mtxAp, const rsb_nnz_idx_t * IB, const rsb_nnz_idx_t * IX, rsb_coo_idx_t * IT, rsb_coo_idx_t * WA, rsb_submatrix_idx_t cmc, rsb_submatrix_idx_t omc, rsb_submatrix_idx_t tmc, rsb_thread_t wet, rsb_submatrix_idx_t *cmcp)
{
	/*
	 	TODO: clean this up.
		Note that rsb__set_num_threads is outlawed here.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t el_size = RSB_SIZEOF(typecode);
	const rsb_nnz_idx_t ttlnz = nnz;		/* total nnz */
	rsb_nnz_idx_t maxnz = nnz;			/* max encountered nnz for a leaf */
	rsb_submatrix_idx_t stmc = RSB_MIN(tmc,wet);	/* submatrices total count */
	rsb_submatrix_idx_t lmc = 1;			/* leaf matrix count */
	rsb_time_t cpt = RSB_TIME_ZERO,dt = RSB_TIME_ZERO;	/* cpt overwrites mtxAp->cpt */
	const rsb_thread_t tn = rsb_get_num_coo2rec_threads();	/* threads number */
	const rsb_thread_t mtn = RSB_MAX(1,(rsb_thread_t)(rsb_global_session_handle.subdivision_multiplier*tn));
	rsb_thread_t tnn = 1;				/* threads number */
	rsb_float_t skew = ((rsb_float_t)(maxnz))/(nnz/wet);	/* if more than one, will limit scaling */
	const long cbs = rsb__get_cache_block_byte_size();
	if(RSB_WANT_VERBOSE_SUBDIVISION)
		RSB_INFO("serial substage subdivision of "),RSB_INFO_MATRIX_SUMMARY(mtxAp),RSB_INFO("\n");
again:
	#pragma omp parallel reduction(|:errval) shared(submatricesp) num_threads(tnn)
	{
		rsb_submatrix_idx_t smi = 0;
		struct rsb_mtx_t * submatrix = NULL;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		rsb_thread_t th_id = omp_get_thread_num();
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		rsb_thread_t th_id = 0;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#if RSB_WANT_VERBOSE_SUBDIVISION2
		if(th_id==0){
			RSB_MTXASM_INFO("entering %d threads in subdivision phase\n",tnn);
			RSB_MTXASM_INFO("entering %d threads in subdivision phase\n",tnn);}
#endif /* RSB_WANT_VERBOSE_SUBDIVISION2 */
iagain:
		#pragma omp critical (rsb_coo2rsbsub_crs)
		{
			smi = rsb_do_pick_largest_open_matrix(submatricesp+cmc,omc);
			if(smi != RSB_SUBM_IDX_MARKER)
			{
				smi += cmc;
				submatrix = submatricesp[smi];
				maxnz = submatrix->nnz;
				/* RSB_ASSERT(nnz>=wet); */
				skew = ((rsb_float_t)(maxnz))/((rsb_float_t)(nnz/wet));
				RSB_ASSERT(!isinf(skew));
				omc--;
				if(smi!=cmc)
				{
#if 0
			  		assert(submatricesp[smi]);
			  		assert(submatricesp[cmc]);
#endif
					RSB_SWAP(struct rsb_mtx_t *,submatricesp[smi],submatricesp[cmc]);
				}
				++cmc;
			}
			else
			{
				submatrix = NULL;
			}
	 		if(RSB_WANT_VERBOSE_SUBDIVISION)
			{
				if(submatrix)
					RSB_INFO("subdividing "),RSB_INFO_MATRIX_SUMMARY(submatrix),RSB_INFO(" (open:%d,closed:%d) for thread %d\n",omc,cmc,th_id);
				else
					RSB_INFO("no available submatrix (open:%d,closed:%d) for thread %d/%d\n",omc,cmc,th_id,tnn);
			}
		}
	if((smi)!=RSB_SUBM_IDX_MARKER)
	{
#if 0
	for(smi=0;RSB_LIKELY(omc>0);smi=cmc+omc-1)
	while(submatrix)
	if(omc>0 && ((submatrix=submatricesp[smi])!=NULL))
#endif
	{
		const rsb_coo_idx_t k = submatrix->nc;
		const rsb_coo_idx_t m = submatrix->nr;
		const rsb_coo_idx_t hk = RSB_MIDDLE(k);
		const rsb_coo_idx_t hm = RSB_MIDDLE(m);
		rsb_nnz_idx_t ul = 0,ur = 0,ll = 0,lr = 0;
		const rsb_nnz_idx_t nnz = submatrix->nnz;
		const rsb_coo_idx_t roff = submatrix->roff;
		const rsb_coo_idx_t coff = submatrix->coff;
		const rsb_flags_t flags = submatrix->flags;
		const rsb_nnz_idx_t nzoff = submatrix->nzoff;
		rsb_bool_t sqs = RSB_BOOL_FALSE;		/* should quad subdivide */
		rsb_submatrix_idx_t smc = 0;	/* submatrices count */

		if(RSB_C2R_IF_VERBOSE)
			RSB_INFO("cmc:%d omc:%d smi:%d tmc=%d stmc=%d th_id=%d\n",cmc,omc,smi,tmc,stmc,th_id);

		/* too few nonzeros for recursion (TODO: may change in the future) */
		if(RSB_DO_TOOFEWNNZFORRCSR(nnz,m))
#if RSB_WANT_SUBDIVISION_FIXES_20101120
		if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COO_STORAGE))
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */
		{
			if(RSB_C2R_IF_VERBOSE)
				RSB_INFO("matrix too sparse for RCSR: rejoining\n");
			sqs = RSB_BOOL_FALSE;
			goto nosqstest;
		}

		/* decide if the matrix is worth subdividing further (soft) */
		sqs = rsb__should_recursively_partition_matrix(0,0,m,k,0,0,nnz,m,k,roff,coff,flags,el_size,mtn);
#if RSB_WANT_SUBDIVISION_FIXES_20101120
		if(nnz<RSB_RECURSION_MIN_NNZ  || m<RSB_RECURSION_MIN_DIM  || k<RSB_RECURSION_MIN_DIM  || !RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
		{
			sqs = RSB_BOOL_FALSE;		/* a hard condition */
			goto nosqstest;
		}
		else
			if(cmc+omc<tmc)
				if(skew>RSB_SUBDIVISION_SKEW_MAX)
					sqs = RSB_BOOL_TRUE;	/* a soft condition */
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

		if(!sqs)
			if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS))
				if(wet>lmc)
					sqs = RSB_BOOL_TRUE;

		if(sqs)
		{
			rsb_bool_t awfcsr = RSB_BOOL_FALSE; /* all of the matrices will fit csr ? */
#if RSB_WANT_SUBDIVISION_FIXES_20101120
			rsb_nnz_idx_t mqnnz = RSB_MAX(RSB_MAX(ul,ur),RSB_MAX(lr,ll));
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

			/* compute the split vector */
			dt = - rsb_time();
			if((errval = rsb_do_compute_vertical_split_search_only(IA,JA,roff,coff,m,k,hm,hk,nnz,IB,&ul,&ur,&ll,&lr))!=RSB_ERR_NO_ERROR) 
				;/* goto err; */
			dt += rsb_time();
			cpt += dt;
			RSB_C2R_ASSERT(IR);
			awfcsr = ( (ul>0 && RSB_DO_TOOFEWNNZFORCSR(ul,hm))   || (ur>0 && RSB_DO_TOOFEWNNZFORCSR(ur,hm)) || (lr>0 && RSB_DO_TOOFEWNNZFORCSR(lr,m-hm)) || (ll>0 && RSB_DO_TOOFEWNNZFORCSR(ll,m-hm)))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;

			if(awfcsr) /* FIXME: misleading naming ! */ 
			{
				/* if some leaf won't fit in CSR, we don't split anymore */
				if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COO_STORAGE))
					sqs = RSB_BOOL_TRUE;
				else
					sqs = RSB_BOOL_FALSE; 
				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("no space for conversion of some leaf: rejoining ? %d\n",!sqs);
			}

#if RSB_WANT_SUBDIVISION_FIXES_20101120
			/* an alternative would be to place this test in the branch above*/
			if(	mqnnz>RSB_MAX_QUADRANTS_UNBALANCE*(nnz-mqnnz) &&
				el_size*nnz<cbs &&
				nnz < (ttlnz/wet) )
				sqs = RSB_BOOL_FALSE; 
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

			/* how many submatrices out of four ? */
			smc = (ul?1:0)+(ur?1:0)+(ll?1:0)+(lr?1:0);
			if(cmc+omc+smc>tmc)
			{	
				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("too many submatrices (%d+%d>%d: rejoining\n",cmc+omc,smc,tmc);
				sqs = RSB_BOOL_FALSE;
				goto nosqstest;
			}

#if !RSB_WANT_SUBDIVISION_FIXES_20101120
#ifdef RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES
			if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES))
				if(wet<lmc-1)
				{	
					if(RSB_C2R_IF_VERBOSE)
						RSB_INFO("RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES: rejoining\n");
					sqs = RSB_BOOL_FALSE;
				}
#endif /* RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES */
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

			if(RSB_C2R_IF_VERBOSE)
			RSB_ERROR("splitting %d/%d -> %d/%d %d/%d %d/%d %d/%d sqs? %d\n",nnz,m,ul,hm,ur,hm,ll,m-hm,lr,m-hm,sqs);
			if(ul+ur+ll+lr != nnz)
			{
				if(RSB_C2R_IF_VERBOSE)
					RSB_ERROR("%d ?= %d + %d + %d + %d = %d\n",nnz,ul,ur,ll,lr,ul+ur+ll+lr);
				RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
			}
		}
nosqstest:
		if(sqs)
		{
			/* should quad-subdivide. let's take care of indices. */
			rsb_nnz_idx_t snzoff = nzoff;
			rsb_submatrix_idx_t smci = 0;
			rsb_submatrix_idx_t smco = 0;
			struct rsb_mtx_t*isms[4] = {NULL,NULL,NULL,NULL};
			
			/*
			the index arrays are copied/linked into the quadrants
			some quadrants may seem ready for recursion, but they not result as such later on.
			they will be made leaf later on, if necessary.
			...
			*/
			RSB_C2R_ASSERT(ur>=0 && ul>=0 && lr>=0 && ll>=0);

			#pragma omp critical (rsb_coo2rsbsub_crs)
			{
				if(cmc+omc+smc+RSB_SUBDIVISION_BUG_EXTRA>tmc)
				{	
					if(RSB_C2R_IF_VERBOSE)
						RSB_INFO("too many submatrices (%d+%d>%d): rejoining\n",cmc+omc,smc,tmc);
					sqs = RSB_BOOL_FALSE;
				}
				else
				{
					lmc += smc;
					lmc -= 1;
					smco = cmc+omc;
					snzoff = nzoff;
					if(ul){ isms[0] = submatricesp[smco+smci];submatricesp[smco+smci] = NULL;smci++;}
					if(ur){ isms[1] = submatricesp[smco+smci];submatricesp[smco+smci] = NULL;smci++;}
					if(ll){ isms[2] = submatricesp[smco+smci];submatricesp[smco+smci] = NULL;smci++;}
					if(lr){ isms[3] = submatricesp[smco+smci];submatricesp[smco+smci] = NULL;smci++;}
					smci = 0;
					omc += smc;
				}
			if(sqs)	
			{
				if(ul)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,isms[0],NULL,NULL,IA,JA,VA,snzoff,ul,hm,hk,0,0,typecode,flags)), snzoff += ul,++smci;
				if(ur)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,isms[1],NULL,NULL,IA,JA,VA,snzoff,ur,hm,k-hk,0,hk,typecode,flags)), snzoff += ur,++smci;
				if(ll)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,isms[2],NULL,NULL,IA,JA,VA,snzoff,ll,m-hm,hk,hm,0,typecode,flags)), snzoff += ll,++smci;
				if(lr)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,isms[3],NULL,NULL,IA,JA,VA,snzoff,lr,m-hm,k-hk,hm,hk,typecode,flags)), snzoff += lr,++smci;
				RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
			}

			if(sqs)
			{
				smci = 0;
				if(ul){ submatricesp[smco+smci] = isms[0];smci++;}
				if(ur){ submatricesp[smco+smci] = isms[1];smci++;}
				if(ll){ submatricesp[smco+smci] = isms[2];smci++;}
				if(lr){ submatricesp[smco+smci] = isms[3];smci++;}
			}

			if(sqs)
			{
			if(snzoff-nzoff!=nnz)
			{
				/* is this a partition ? */
				RSB_ERROR("%d - %d != %d ?= %d + %d + %d + %d = %d\n",snzoff,nzoff,nnz,ul,ur,ll,lr,ul+ur+ll+lr);
				RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
			}
			if(RSB_SOME_ERROR(errval))
			{
				RSB_ERROR(RSB_ERRM_ES); /* goto err; */
			}
			RSB_DO_FLAG_ADD(submatrix->flags,RSB_FLAG_QUAD_PARTITIONING);
			RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_NON_ROOT_MATRIX);
			submatrix->bindx = NULL; submatrix->bpntr = NULL; submatrix->indptr = NULL;
			}
			}
		}
		if(!sqs)
		{
			RSB_DO_FLAG_SUBST(submatrix->flags,RSB_FLAG_QUAD_PARTITIONING,RSB_FLAG_NON_ROOT_MATRIX);
			/* selecting a format and declaring as leaf */
			if(!RSB_DO_TOOFEWNNZFORCSR(nnz,m) /*&& IR && IL*/)
			{
/*				RSB_INFO("CSR -> COO ?\n"); */
				if(RSB_DO_FLAG_HAS(submatrix->flags,RSB_FLAG_WANT_BCSS_STORAGE))
					RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_WANT_COO_STORAGE);
				if((errval = rsb__do_set_init_storage_flags(submatrix,submatrix->flags))!=RSB_ERR_NO_ERROR)
					;/* goto err; */
			}
			else
			{
/*				RSB_INFO("COO !\n"); */
				rsb_flags_t sflags = flags;
				RSB_DO_FLAG_SUBST(sflags,RSB_FLAG_WANT_BCSS_STORAGE,RSB_FLAG_WANT_COO_STORAGE);
				if((errval = rsb__do_set_init_storage_flags(submatrix,sflags))!=RSB_ERR_NO_ERROR)
					;/* goto err; */
			}	
			if(RSB_C2R_IF_VERBOSE)
				RSB_INFO("freezing %d ",smi+1),
				RSB_INFO_MATRIX_SUMMARY(submatrix),
				RSB_INFO("\n");
		}
		/* matrix is declared as 'closed'.
		   sorting, in a way smi will point to the biggest open mtxAp, which will be picked up next */
#if 0
		qsort(submatricesp+cmc,(size_t)(omc),sizeof(struct rsb_mtx_t*),& rsb_compar_rcsr_matrix_regarding_nnz);
		RSB_SWAP(struct rsb_mtx_t *,submatricesp[smi],submatricesp[cmc-1]);
#endif
	}
		/*smi = cmc+omc-1; */
	}
		if(omc>0 && cmc+omc<stmc)
			goto iagain;
#if RSB_WANT_VERBOSE_SUBDIVISION2
		if(th_id==0)
		{
			RSB_MTXASM_INFO("thread %d:terminating subdivision",th_id);
			if(omc==0)
			{RSB_MTXASM_INFO(", no more open matrices");}
			RSB_MTXASM_INFO("(closed %d= %d nodes + %d leaves, out of %d available)",cmc,cmc-lmc,lmc,tmc);
			RSB_MTXASM_INFO(",(maxnz=%d,skew=%g)",maxnz,skew);
			if(cmc+omc>=stmc)
			{RSB_MTXASM_INFO(", no room left for submatrices");}
			RSB_MTXASM_INFO(".\n");
		}
#endif
	} /* parallel */

	if(RSB_SOME_ERROR(errval))
		goto err;

	#pragma omp barrier
	if(stmc!=tmc)
	{
		stmc = tmc;
		if(RSB_WANT_VERBOSE_SUBDIVISION)
			RSB_INFO("parallel substage subdivision of "),RSB_INFO_MATRIX_SUMMARY(mtxAp),RSB_INFO("\n");
		tnn = tn;
		goto again;
	}
	else
	{
	
		if(RSB_WANT_VERBOSE_SUBDIVISION)
			RSB_INFO("parallel substage subdivision of "),RSB_INFO_MATRIX_SUMMARY(mtxAp),RSB_INFO(" not required\n");
	}
	{
 		if(RSB_WANT_VERBOSE_SUBDIVISION)
			RSB_INFO("subdivision of "),RSB_INFO_MATRIX_SUMMARY(mtxAp),RSB_INFO("complete \n");
	}
	mtxAp->cpt = cpt;

	*cmcp = cmc;
err:
	RSB_DO_ERR_RETURN(errval)
}

#if 0
static rsb_err_t rsb_do_coo2rec_subdivide(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_err_t *errvalp, struct rsb_mtx_t ** submatricesp, struct rsb_mtx_t * mtxAp, const rsb_nnz_idx_t * IB, const rsb_nnz_idx_t * IX, rsb_coo_idx_t * IT, rsb_coo_idx_t * WA, rsb_submatrix_idx_t cmc, rsb_submatrix_idx_t omc, rsb_submatrix_idx_t tmc, rsb_thread_t wet, rsb_submatrix_idx_t *cmcp)
{
	rsb_nnz_idx_t tdnnz = 0;
	rsb_submatrix_idx_t smi = 0;/* max matrix count, done matrix count, submatrix index */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t el_size = RSB_SIZEOF(typecode);
	rsb_time_t cpt = RSB_TIME_ZERO,dt = RSB_TIME_ZERO;

	for(smi=0;RSB_LIKELY(omc>0);smi=cmc+omc-1)
	{
		struct rsb_mtx_t * submatrix = submatricesp[smi];
		rsb_coo_idx_t k = submatrix->nc;
		rsb_coo_idx_t m = submatrix->nr;
		rsb_coo_idx_t hk = (k+1)/2;
		rsb_coo_idx_t hm = (m+1)/2;
		rsb_nnz_idx_t ul = 0,ur = 0,ll = 0,lr = 0;
		rsb_nnz_idx_t nnz = submatrix->nnz;
		rsb_coo_idx_t roff = submatrix->roff;
		rsb_coo_idx_t coff = submatrix->coff;
		rsb_coo_idx_t*IL = submatrix->bindx;	// IL will be hosted here
		rsb_coo_idx_t*IM = IT;			// IM will be hosted in a temporary vector
		rsb_coo_idx_t*IR = submatrix->bpntr;	// IR will be hosted here
		rsb_flags_t flags = submatrix->flags;
		rsb_nnz_idx_t nzoff = submatrix->nzoff;
		rsb_bool_t sqs = RSB_BOOL_FALSE;		// should quad subdivide
		rsb_submatrix_idx_t smc = 0;	/* submatrices count */

//		RSB_INFO("picked up %d/%d -> %d x %d, %d nnz, @ %d %d \n",smi+1,tmc,m,k,nnz,roff,coff);

		if(!IL && !RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COO_STORAGE) )
		{
			/* if there is no line pointer, we make this submatrix leaf */
			sqs = RSB_BOOL_FALSE;
			if(RSB_C2R_IF_VERBOSE)
				RSB_INFO("no split, as no line pointer found\n");
			goto nosqstest;
		}

		if(/*!IL || */!IM)
		{
			/* if this happens, this is an error */
			RSB_ERROR("IL:%p, IM:%p\n",IL,IM);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;
		}
	
		/* too few nonzeros for recursion (TODO: may change in the future) */
		if(RSB_DO_TOOFEWNNZFORRCSR(nnz,m) 
				/* 
				 * Uncommenting the following allows subdivision for very spare matrices 
				 * However, this feature is unfinished and bugful (segfault risk)
				 * */
			       /*	&& !RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COO_STORAGE) */
				)
		{
			if(RSB_C2R_IF_VERBOSE)
				RSB_INFO("matrix too sparse for RCSR: rejoining\n");

			sqs = RSB_BOOL_FALSE;
			goto nosqstest;
		}

		/* decide if the matrix is worth subdividing further */
		sqs = rsb__should_recursively_partition_matrix(0,0,m,k,0,0,nnz,m,k,roff,coff,flags,el_size,0);

		/* if we want subdivision  */

		if(!sqs)
		if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS))
		if(wet>cmc+1+smc)/* /FIXME : this may not terminate! */
			sqs = RSB_BOOL_TRUE;

		if(sqs)
		{
			rsb_bool_t awfcsr = RSB_BOOL_FALSE; /* all of the matrices will fit csr ? */

			// compute the split vector
			dt = - rsb_time();
			if((!RSB_DO_TOOFEWNNZFORCSR(nnz,m)) && IR && IL)
			{if((errval = rsb_do_compute_vertical_split(IA,JA,roff,coff,m,k,hm,hk,nnz,IL,IM,IR,&ul,&ur,&ll,&lr))!=RSB_ERR_NO_ERROR) goto err;}
			else
			{
				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("using the sparse splitter\n");
				if((errval = rsb_do_compute_vertical_split_search_only(IA,JA,roff,coff,m,k,hm,hk,nnz,IB,&ul,&ur,&ll,&lr))!=RSB_ERR_NO_ERROR) goto err;
			}
			dt += rsb_time();
			cpt += dt;
			RSB_C2R_ASSERT(IR);
			awfcsr = ( (ul>0 && RSB_DO_TOOFEWNNZFORCSR(ul,hm))   || (ur>0 && RSB_DO_TOOFEWNNZFORCSR(ur,hm)) || (lr>0 && RSB_DO_TOOFEWNNZFORCSR(lr,m-hm)) || (ll>0 && RSB_DO_TOOFEWNNZFORCSR(ll,m-hm)))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;

			// after computing the split vector, we can still resign from subdividing
			// especially if some submatrix is deemed too small and the overall submatrices count is enough
			// TODO: if(rsb__should_rejoin_small_leaf(...
			// ...
			
			/* is there room for these additional submatrices ? */
//			if( (ul>0 && RSB_DO_TOOFEWNNZFORRCSR(ul,hm))   || (ur>0 && RSB_DO_TOOFEWNNZFORRCSR(ur,hm)) || (lr>0 && RSB_DO_TOOFEWNNZFORRCSR(lr,m-hm)) || (ll>0 && RSB_DO_TOOFEWNNZFORRCSR(ll,m-hm)))
			if(awfcsr)
			{
				/* if some leaf won't fit in CSR, we don't split anymore */
				if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COO_STORAGE))
					sqs = RSB_BOOL_TRUE;
				else
					sqs = RSB_BOOL_FALSE; 

				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("no space for conversion of some leaf: rejoining ? %d\n",!sqs);
			}

			/* how many submatrices out of four ? */
			smc = (ul?1:0)+(ur?1:0)+(ll?1:0)+(lr?1:0);
			if(cmc+omc+smc>tmc)
			{	
				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("too many submatrices: rejoining\n");
				sqs = RSB_BOOL_FALSE;
			}

#ifdef RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES
			/* 
 			  if we want to avoid micro leaves, we could stop here 
 			  FIXME: we need a better criteria (for proper load balancing!)
 			*/
			if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES))
				if(wet<cmc)
				{	
					if(RSB_C2R_IF_VERBOSE)
						RSB_INFO("RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES: rejoining\n");
					sqs = RSB_BOOL_FALSE;
				}
#endif /* RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES */

			if(RSB_C2R_IF_VERBOSE)
			RSB_INFO("splitting %d/%d -> %d/%d %d/%d %d/%d %d/%d sqs? %d\n",nnz,m,ul,hm,ur,hm,ll,m-hm,lr,m-hm,sqs);
		}

nosqstest:
		omc--;
		if(smi!=cmc)
			RSB_SWAP(struct rsb_mtx_t *,submatricesp[smi],submatricesp[cmc]);
		++cmc;
		if(sqs)
		{
			/* should quad-subdivide. let's take care of indices. */
			rsb_nnz_idx_t snzoff = nzoff;
			
			// the index arrays are copied/linked into the quadrants
			// some quadrants may seem ready for recursion, but they not result as such later on.
			// they will be made leaf later on, if necessary.
			// ...
			RSB_C2R_ASSERT(ur>=0 && ul>=0 && lr>=0 && ll>=0);

#if RSB_C2R_WANT_MAYBE_FASTER 
			if(IL)
				RSB_COA_MEMCPY_ROWSZ(IX,IL,0  ,0,m+1),
				IL = IX;
			if(IR)
				RSB_COA_MEMCPY_ROWSZ(IX,IR,m+1,0,m+1),
				IR = IX+m+1;
#else /* RSB_C2R_WANT_MAYBE_FASTER */
			if(IL)
				RSB_COA_MEMMOVE(IX,IL,0  ,0,m+1),
				IL = IX;
			if(IR)
				RSB_COA_MEMMOVE(IX,IR,m+1,0,m+1),
				IR = IX+m+1;
#endif /* RSB_C2R_WANT_MAYBE_FASTER */

			if(ul)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,submatricesp[cmc+omc],IL,IM,IA,JA,VA,snzoff,ul,hm,hk,0,0,typecode,flags), ++omc, snzoff += ul);
			if(ur)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,submatricesp[cmc+omc],IM,IR,IA,JA,VA,snzoff,ur,hm,k-hk,0,hk,typecode,flags), ++omc, snzoff += ur);
			if(ll)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,submatricesp[cmc+omc],IL,IM,IA,JA,VA,snzoff,ll,m-hm,hk,hm,0,typecode,flags), ++omc, snzoff += ll);
			if(lr)
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_fill_early_leaf_matrix(submatrix,submatricesp[cmc+omc],IM,IR,IA,JA,VA,snzoff,lr,m-hm,k-hk,hm,hk,typecode,flags), ++omc, snzoff += lr);

			if(snzoff-nzoff!=nnz)
			{
				/* is this a partition ? */
				RSB_ERROR("%d - %d != %d ?= %d + %d + %d + %d = %d\n",snzoff,nzoff,nnz,ul,ur,ll,lr,ul+ur+ll+lr);
				errval = RSB_ERR_INTERNAL_ERROR; goto err;
			}
			if(RSB_SOME_ERROR(errval))
			{
				RSB_ERROR(RSB_ERRM_ES); goto err;
			}
			RSB_DO_FLAG_ADD(submatrix->flags,RSB_FLAG_QUAD_PARTITIONING);
			RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_NON_ROOT_MATRIX);
			submatrix->bindx = NULL; submatrix->bpntr = NULL; submatrix->indptr = NULL;
		}
		else
		{
			RSB_DO_FLAG_SUBST(submatrix->flags,RSB_FLAG_QUAD_PARTITIONING,RSB_FLAG_NON_ROOT_MATRIX);
			// we should decide a format, and proceed declaring it as leaf
			if(!RSB_DO_TOOFEWNNZFORRCSR(nnz,m) && IR && IL)
			{
				if(RSB_DO_FLAG_HAS(submatrix->flags,RSB_FLAG_WANT_BCSS_STORAGE))
					RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_WANT_COO_STORAGE);
#if RSB_WANT_LITTLE_IMPROVED 
				if(submatrix==mtxAp)/*  root only */
					/* FIXME: TODO: IR is NOT needed AT ALL!  */
					rsb_do_fill_rcsr_arrays_for_later(submatrix,IL,IR,IA,JA,nzoff,m,0);
#else /* RSB_WANT_LITTLE_IMPROVED */
					rsb_do_fill_rcsr_arrays_for_later(submatrix,IL,IR,IA,JA,nzoff,m,0);
#endif /* RSB_WANT_LITTLE_IMPROVED */
				if((errval = rsb__do_set_init_storage_flags(submatrix,submatrix->flags))!=RSB_ERR_NO_ERROR)
					goto err;
				submatrix->VA = VA;	// FIXME: we will place pointers to partially swapped VA, here.
			}
			else
			if(!RSB_DO_TOOFEWNNZFORCSR(nnz,m) /*&& IR && IL*/)
			{
//				RSB_INFO("CSR -> COO ?\n");
				if(RSB_DO_FLAG_HAS(submatrix->flags,RSB_FLAG_WANT_BCSS_STORAGE))
					RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_WANT_COO_STORAGE);
				if((errval = rsb__do_set_init_storage_flags(submatrix,submatrix->flags))!=RSB_ERR_NO_ERROR)
					goto err;
			}
			else
			{
//				RSB_INFO("COO !\n");
				rsb_flags_t sflags = flags;
				RSB_DO_FLAG_SUBST(sflags,RSB_FLAG_WANT_BCSS_STORAGE,RSB_FLAG_WANT_COO_STORAGE);
				if((errval = rsb__do_set_init_storage_flags(submatrix,sflags))!=RSB_ERR_NO_ERROR)
					goto err;
			}	
			if(RSB_C2R_IF_VERBOSE)
				RSB_INFO("freezing %d ",smi+1),
				RSB_INFO_MATRIX_SUMMARY(submatrix),
				RSB_INFO("\n");
		}
		// matrix is declared as 'closed'.
		// sorting, in a way smi will point to the biggest open mtxAp, which will be picked up next
		qsort(submatricesp+cmc,(size_t)(omc),sizeof(struct rsb_mtx_t*),& rsb_compar_rcsr_matrix_regarding_nnz);
		// FIXME: a priority queue would do the job, here
	}
	mtxAp->cpt = cpt;
err:
	*cmcp = cmc;
	RSB_DO_ERR_RETURN(errval)
}
#endif

static rsb_err_t rsb_do_coo2rec_shuffle(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_err_t *errvalp, struct rsb_mtx_t ** submatricesp, struct rsb_mtx_t * mtxAp, const rsb_nnz_idx_t * IB, rsb_coo_idx_t * WA, rsb_submatrix_idx_t cmc)
{
	rsb_nnz_idx_t tdnnz = 0;
	rsb_submatrix_idx_t smi = 0;/* max matrix count, done matrix count, submatrix index */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t el_size = RSB_SIZEOF(typecode);
#if RSB_WANT_VERBOSE_TIMINGS
	rsb_time_t pmt = RSB_TIME_ZERO;
#endif /* RSB_WANT_VERBOSE_TIMINGS */

	if(!VA || cmc==1)
	{
		mtxAp->VA = VA;
		goto no_va_cp;
	}

	// the following is a highly parallel phase
	#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(tdnnz,submatricesp,IB) RSB_NTC 
	for(smi=0;smi<cmc;++smi)
	{
		struct rsb_mtx_t * submatrix = submatricesp[smi];
#if RSB_OBSOLETE_QUARANTINE
		const rsb_coo_idx_t*IL = submatrix->bindx;
		const rsb_coo_idx_t*IR = submatrix->bpntr;
#endif /* RSB_OBSOLETE_QUARANTINE */
		rsb_nnz_idx_t dnnz = 0;
		rsb_coo_idx_t i;
//		if(RSB_C2R_IF_VERBOSE)
//		RSB_INFO("%d -> %d\n",smi,omp_get_thread_num());
		if(!RSB_DO_FLAG_HAS((submatrix->flags),RSB_FLAG_QUAD_PARTITIONING))
		{
#if RSB_OBSOLETE_QUARANTINE
			if(RSB_DO_ENOUGHNNZFORINDEXBASEDBUILD(submatrix) && !rsb__is_coo_matrix(submatrix) && IL && IR)
			{
				RSB_C2R_ASSERT(IL);
				RSB_C2R_ASSERT(IR);
				dnnz = rsb_do_copy_submatrix_coa(submatrix,VA,WA,IL,IR,el_size,0,0,submatrix->nr);
#if 0
				for(i=submatrix->roff;RSB_LIKELY(i<submatrix->roff+submatrix->nr);++i)
				{
					rsb_nnz_idx_t nnz1 = IR[i-submatrix->roff];
					rsb_nnz_idx_t nnz0 = IL[i-submatrix->roff];
					RSB_C2R_ASSERT(IL[i-submatrix->roff]>=IB[i]);
					RSB_C2R_ASSERT(IL[i-submatrix->roff]<=IB[i+1]);
					RSB_C2R_ASSERT(IL[i-submatrix->roff+1]>=IB[i+1]);
					RSB_C2R_ASSERT(nnz0>=IL[i-submatrix->roff]);
					RSB_C2R_ASSERT(nnz1<=IL[i-submatrix->roff+1]);
				       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					if(IB[i]==IB[i+1])continue;
					RSB_C2R_ASSERT(nnz1>=nnz0);
					if(nnz1==nnz0)continue;
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					RSB_C2R_ASSERT(JA[nnz1-1]< submatrix->coff+submatrix->nc);
				}
#endif
			}
			else
#endif /* RSB_OBSOLETE_QUARANTINE */
			if(!RSB_DO_TOOFEWNNZFORCSR(submatrix->nnz,submatrix->nr))
			{
				//rsb_coo_idx_t*IL = submatrix->bindx;
				for(i=submatrix->roff;RSB_LIKELY(i<submatrix->roff+submatrix->nr);++i)
				{
					//rsb_nnz_idx_t fel;
					// offset of line i in the global line pointers array
					rsb_nnz_idx_t nnz0 = IB[i];
					// nnz1..nnz0 are the boundaries of line i
					rsb_nnz_idx_t nnz1 = IB[i+1];
					// check
					RSB_C2R_ASSERT(nnz0>=IB[i]);
					RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					// skip line if empty
					if(nnz1-nnz0<1)continue;
					// find first element of line i also in the submatrix
					nnz0 += rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff,nnz1-nnz0);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)continue;
					// find the length of the subrow i in the submatrix
					nnz1 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff+submatrix->nc,nnz1-nnz0);
					//check 
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)continue;
					// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
					// checks
//					RSB_C2R_ASSERT(IL[i-submatrix->roff]>=IB[i]);
//					RSB_C2R_ASSERT(IL[i-submatrix->roff]<=IB[i+1]);
//					RSB_C2R_ASSERT(IL[i-submatrix->roff+1]>=IB[i+1]);
//					RSB_C2R_ASSERT(nnz0>=IL[i-submatrix->roff]);
//					RSB_C2R_ASSERT(nnz1<=IL[i-submatrix->roff+1]);
				       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					RSB_C2R_ASSERT(JA[nnz1-1]< submatrix->coff+submatrix->nc);
					// perform the copy
					RSB_A_MEMCPY_SMALL(WA,VA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0,el_size);
					//RSB_COA_MEMCPY(WA,JA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0);
					// update the actual offset in the destination array
					dnnz += nnz1-nnz0;
				}
			}
			else
			{
				//rsb_coo_idx_t*IL = submatrix->bindx;
				for(i=submatrix->roff;RSB_LIKELY(i<submatrix->roff+submatrix->nr);++i)
				{
					//rsb_nnz_idx_t fel;
					// offset of line i in the global line pointers array
					rsb_nnz_idx_t nnz0 = IB[i];
					// nnz1..nnz0 are the boundaries of line i
					rsb_nnz_idx_t nnz1 = IB[i+1];
					// check
					RSB_C2R_ASSERT(nnz0>=IB[i]);
					RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					// skip line if empty
					if(nnz1-nnz0<1)continue;
					// find first element of line i also in the submatrix
					nnz0 += rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff,nnz1-nnz0);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)continue;
					// find the length of the subrow i in the submatrix
					nnz1 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff+submatrix->nc,nnz1-nnz0);
					//check 
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)continue;
					// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
					// checks
//					RSB_C2R_ASSERT(IL[i-submatrix->roff]>=IB[i]);
//					RSB_C2R_ASSERT(IL[i-submatrix->roff]<=IB[i+1]);
//					RSB_C2R_ASSERT(IL[i-submatrix->roff+1]>=IB[i+1]);
//					RSB_C2R_ASSERT(nnz0>=IL[i-submatrix->roff]);
//					RSB_C2R_ASSERT(nnz1<=IL[i-submatrix->roff+1]);
				       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					RSB_C2R_ASSERT(JA[nnz1-1]< submatrix->coff+submatrix->nc);
					// perform the copy
					RSB_A_MEMCPY_SMALL(WA,VA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0,el_size);
					//RSB_COA_MEMCPY(WA,JA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0);
					// update the actual offset in the destination array
					dnnz += nnz1-nnz0;
				}
			}
			if(dnnz!=submatrix->nnz)
			{
				RSB_ERROR("@%d,%d: found %d, should have found %d\n",
						submatrix->roff, submatrix->coff, dnnz,submatrix->nnz);
				RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
			}
			#pragma omp critical (rsb_coo2rsb_nzinc_crs)
			{tdnnz += dnnz;}
		}
	}
//gerr:
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
#if   !defined(__xlC__)
	/* FIXME: xlc does not allow this, but we have experienced problems, without */
	#pragma omp barrier
#endif /* __xlC__ */
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_NL);
	}

	if(tdnnz!=nnz)
	{
		RSB_ERROR("found %d, should have found %d\n", tdnnz,nnz);
		errval = RSB_ERR_INTERNAL_ERROR;
	       	goto err;
	}

#if RSB_WANT_VERBOSE_TIMINGS
	pmt -= rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */
	RSB_A_MEMCPY_parallel(VA,WA,0,0,nnz,el_size);
#if RSB_WANT_VERBOSE_TIMINGS
	pmt += rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */
no_va_cp:

	tdnnz = 0;
	#pragma omp parallel for schedule(static,1) reduction(|:errval)  shared(tdnnz,submatricesp,IB) RSB_NTC 
	for(smi=0;smi<cmc;++smi)
	{
		struct rsb_mtx_t * submatrix = submatricesp[smi];
		if(!RSB_DO_FLAG_HAS(submatrix->flags,RSB_FLAG_QUAD_PARTITIONING))
		{
			rsb_coo_idx_t oll,nll;
			rsb_coo_idx_t i;
			const rsb_coo_idx_t*IR = submatrix->bpntr;
			const rsb_coo_idx_t*IL = submatrix->bindx;
			rsb_nnz_idx_t dnnz = 0;

			//if(!RSB_DO_TOOFEWNNZFORRCSR(submatrix->nnz,submatrix->nr))
			if(RSB_DO_ENOUGHNNZFORINDEXBASEDBUILD(submatrix) && !rsb__is_coo_matrix(submatrix) && IL && IR)
			{
				if(RSB_C2R_IF_VERBOSE)
				RSB_INFO("CSR:%zd/%zd:%zd..%zd\n",(rsb_printf_int_t)smi,(rsb_printf_int_t)cmc,(rsb_printf_int_t)submatrix->nzoff,(rsb_printf_int_t)(submatrix->nzoff+submatrix->nnz));
				submatrix->bpntr = IA+submatrix->nzoff;
				RSB_C2R_ASSERT(IL); RSB_C2R_ASSERT(IR);
				RSB_C2R_ASSERT(IR< IA+submatrix->nzoff+submatrix->nnz);
				RSB_C2R_ASSERT(IL< IA+submatrix->nzoff+submatrix->nnz);
				RSB_C2R_ASSERT(IR>=IA+submatrix->nzoff);
				RSB_C2R_ASSERT(IL>=IA+submatrix->nzoff);
				for(dnnz=0,i=0;RSB_LIKELY(i<submatrix->nr);dnnz += IR[i]-IL[i],++i)
					RSB_COA_MEMCPY_SMALL(WA,JA,submatrix->nzoff+dnnz,IL[i],IR[i]-IL[i]);

//				RSB_INFO("%d x %d (%d) @ %d, %d (rcsr)\n",submatrix->nr,submatrix->nc,submatrix->nnz,submatrix->roff,submatrix->coff);
				if(dnnz!=submatrix->nnz)
				{
					RSB_ERROR("@%d,%d: found %d, should have found %d\n",
							submatrix->roff, submatrix->coff, dnnz,submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
				}
				RSB_C2R_ASSERT(IL); RSB_C2R_ASSERT(IR);

				oll = IR[0]-IL[0];
				submatrix->bpntr[0] = 0;
				for(i=1;RSB_LIKELY(i<submatrix->nr);++i)
				{
					nll = IR[i]-IL[i];
					submatrix->bpntr[i] = submatrix->bpntr[i-1]+oll;
					oll = nll;
				}
				submatrix->bpntr[submatrix->nr] = submatrix->bpntr[submatrix->nr-1]+oll;
				if(submatrix->bpntr[submatrix->nr]!=submatrix->nnz)
				{
					RSB_ERROR("@%d,%d: found %d, should have found %d\n",
					submatrix->roff,submatrix->coff,submatrix->bpntr[submatrix->nr],submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
				}
			}
			else
			if(!RSB_DO_TOOFEWNNZFORCSR(submatrix->nnz,submatrix->nr) && !rsb__is_coo_matrix(submatrix))
			{
				if(RSB_C2R_IF_VERBOSE)
					RSB_INFO("CSR:%zd/%zd:%zd..%zd\n",(rsb_printf_int_t)smi,(rsb_printf_int_t)cmc,(rsb_printf_int_t)submatrix->nzoff,(rsb_printf_int_t)(submatrix->nzoff+submatrix->nnz));
				oll = 0;
				submatrix->bpntr = IA+submatrix->nzoff;
				submatrix->bpntr[0] = 0;
//				RSB_INFO("%d x %d (%d) @ %d, %d (rcsr)\n",submatrix->nr,submatrix->nc,submatrix->nnz,submatrix->roff,submatrix->coff);
				for(i=submatrix->roff;RSB_LIKELY(i<submatrix->roff+submatrix->nr);++i)
				{
					//rsb_nnz_idx_t fel;
					// offset of line i in the global line pointers array
					rsb_nnz_idx_t nnz0 = IB[i];
					// nnz1..nnz0 are the boundaries of line i
					rsb_nnz_idx_t nnz1 = IB[i+1];
					// check
					RSB_C2R_ASSERT(nnz0>=IB[i]);
					RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					// skip line if empty
					if(nnz1-nnz0<1)goto is_empty_subrow;
					// find first element of line i also in the submatrix
					nnz0 += rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff,nnz1-nnz0);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)goto is_empty_subrow;
					// find the length of the subrow i in the submatrix
					nnz1 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff+submatrix->nc,nnz1-nnz0);
					//check 
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)goto is_empty_subrow;
					// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
					// checks
				       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					RSB_C2R_ASSERT(JA[nnz1-1]< submatrix->coff+submatrix->nc);
					// convert row indices
					nll = nnz1-nnz0;
					if(i>submatrix->roff)
						submatrix->bpntr[i-submatrix->roff] = submatrix->bpntr[i-submatrix->roff-1]+oll;
					oll = nll;
					// perform the copy
					RSB_COA_MEMCPY_SMALL(WA,JA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0);
					// update the actual offset in the destination array
					dnnz += nnz1-nnz0;
					continue;
is_empty_subrow:
					// convert row indices
					nll = 0;
					if(RSB_LIKELY(i>submatrix->roff))
						submatrix->bpntr[i-submatrix->roff] = submatrix->bpntr[i-submatrix->roff-1]+oll;
					oll = nll;
				}
				submatrix->bpntr[submatrix->nr] = submatrix->bpntr[submatrix->nr-1]+oll;
				if(dnnz!=submatrix->nnz || submatrix->bpntr[submatrix->nr]!=submatrix->nnz)
				{
					RSB_ERROR("@%d,%d: found %d, and %d; should have found %d\n",
							submatrix->roff, submatrix->coff,
							dnnz, submatrix->bpntr[submatrix->nr],submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
				}
			}
			else
			{
				if(RSB_C2R_IF_VERBOSE)
				RSB_INFO("COO:%zd/%zd:%zd..%zd\n",(rsb_printf_int_t)smi,(rsb_printf_int_t)cmc,(rsb_printf_int_t)submatrix->nzoff,(rsb_printf_int_t)(submatrix->nzoff+submatrix->nnz));
				oll = 0;
				submatrix->bpntr = IA+submatrix->nzoff;
				submatrix->bpntr[0] = 0;
				for(i=submatrix->roff;RSB_LIKELY(i<submatrix->roff+submatrix->nr);++i)
				{
					//rsb_nnz_idx_t fel;
					// offset of line i in the global line pointers array
					rsb_nnz_idx_t nnz0 = IB[i];
					// nnz1..nnz0 are the boundaries of line i
					rsb_nnz_idx_t nnz1 = IB[i+1];
					// check
					RSB_C2R_ASSERT(nnz0>=IB[i]);
					RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					// skip line if empty
					if(nnz1-nnz0<1)continue;
					// find first element of line i also in the submatrix
					nnz0 += rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff,nnz1-nnz0);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)continue;
					// find the length of the subrow i in the submatrix
					nnz1 = nnz0+rsb__nnz_split_coo_bsearch(JA+nnz0,submatrix->coff+submatrix->nc,nnz1-nnz0);
					//check 
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					// skip line if empty in the submatrix
					if(nnz1-nnz0<1)continue;
					// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
					// checks
				       	RSB_C2R_ASSERT(nnz1<=IB[i+1]);
					RSB_C2R_ASSERT(JA[nnz0+0]>=submatrix->coff);
					RSB_C2R_ASSERT(JA[nnz1-1]< submatrix->coff+submatrix->nc);
					// convert row indices
					// perform the copy
					RSB_COA_MEMCPY_SMALL(WA,JA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0);
					rsb__util_coo_array_set(IA+submatrix->nzoff+dnnz,nnz1-nnz0,i-submatrix->roff);
					// update the actual offset in the destination array
					dnnz += nnz1-nnz0;
				}
				if(dnnz!=submatrix->nnz )
				{
					RSB_ERROR("@%d,%d: found %d; should have found %d\n",
							submatrix->roff, submatrix->coff, dnnz, submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
				}
			}
			rsb__util_coo_array_sub(WA+submatrix->nzoff,dnnz,submatrix->coff);
			#pragma omp critical (rsb_coo2rsb_nzinc_crs)
			{tdnnz += dnnz;}
		}
	}
	#pragma omp barrier
#if RSB_WANT_VERBOSE_TIMINGS
	pmt -= rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */
#if 1
	RSB_COA_MEMCPY_parallel(JA,WA,0,0,nnz);
#else
	RSB_COA_MEMCPY(JA,WA,0,0,nnz);
#endif

#if RSB_WANT_VERBOSE_TIMINGS
	pmt += rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_NL);
	}

	if(tdnnz!=nnz)
	{
		RSB_ERROR("found %d, should have found %d\n", tdnnz,nnz);
		errval = RSB_ERR_INTERNAL_ERROR;
	       	goto err;
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_submatrix_idx_t rsb__estimate_subm_count(const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_flags_t flags, const rsb_thread_t wet, rsb_err_t *errvalp)
{
	const long fcs = rsb__get_first_level_c_size();
	const long lcs = rsb__get_lastlevel_c_size();
	const long lcspt = rsb__get_lastlevel_c_size_per_thread();
	const long cbs = rsb__get_cache_block_byte_size();
#if RSB_WANT_SUBDIVISION_FIXES_20101213
	const long wcbs = cbs;
#else /* RSB_WANT_SUBDIVISION_FIXES_20101213 */
	const long wcbs = lcspt;
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101213 */
	rsb_submatrix_idx_t tmc = 0; /* total (max) matrix count */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( fcs > lcs || fcs<1 || cbs<1 ) /* we allow declaration of 1 level of cache only */
	{
		/* TODO : find a reasonable solution, and declare it in ./rsbench -I, which should give some diagnostic about this */
		errval = RSB_ERR_FAILED_MEMHIER_DETECTION;
		RSB_PERR_GOTO(err,"innermost cache size:%d, outermost cache size:%d, cache block size %d\n",fcs,lcs,cbs);
	}

	tmc = RSB_SUBDIVISION_BUG_EXTRA+2*(((nnz+wcbs)*(RSB_SIZEOF(typecode)+2*sizeof(rsb_coo_idx_t)))/(wcbs)); /* TODO: clean this up */
	tmc = RSB_MAX(1,(rsb_submatrix_idx_t)(rsb_global_session_handle.subdivision_multiplier*tmc));

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS))
		tmc = RSB_MAX(tmc,wet);

#if !RSB_OLD_COO_CRITERIA
	if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_QUAD_PARTITIONING))
	{
		tmc = 1; /* No quad partitioning ? then single submatrix. */
	}
#endif
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	RSB_DEBUG_ASSERT(tmc>0);
	return tmc;
}

struct rsb_mtx_t * rsb__allocate_recursive_sparse_matrix_from_row_major_coo(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_err_t *errvalp)
{
	/**
		\ingroup gr_experimental

		Assembly a complete R-CSR/R-CSC matrix in-place in the provided COO arrays.

		TODO: characterize memory usage.
		TODO: get rid of pinfop.
		TODO: interleave the matrix structs into the data arrays.
		TODO: guaranteed preallocation needed
		TODO: may continue the line of rsb_allocate_new__ and plug it here.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * submatrices = NULL;
	struct rsb_mtx_t ** submatricesp = NULL;
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_time_t hct = RSB_TIME_ZERO;
	rsb_time_t drt = RSB_TIME_ZERO;
	rsb_time_t mat = RSB_TIME_ZERO;
	rsb_time_t lst = RSB_TIME_ZERO;
	rsb_coo_idx_t * IB = NULL;
	rsb_coo_idx_t * IT = NULL;
	rsb_coo_idx_t * IX = NULL;
	rsb_coo_idx_t * WA = NULL;
	rsb_submatrix_idx_t smi = 0; /* submatrix index */
	rsb_submatrix_idx_t cmc = 0, omc = 0; /* closed matrices count, open matrices count */
	rsb_submatrix_idx_t lm = 0;           /* leaf matrices */
	rsb_time_t dt = RSB_TIME_ZERO,eit = RSB_TIME_ZERO,est = RSB_TIME_ZERO,ect = RSB_TIME_ZERO,tat = RSB_TIME_ZERO,sat = RSB_TIME_ZERO;
	const rsb_thread_t wet = rsb_get_num_coo2rec_threads(); /* want executing threads: */
	const size_t el_size = RSB_SIZEOF(typecode);
	rsb_coo_idx_t roff = 0;
	rsb_coo_idx_t coff = 0;
	rsb_nnz_idx_t dnnz = 0;
	const rsb_submatrix_idx_t tmc = rsb__estimate_subm_count(nnz,typecode,flags,wet,&errval); /* total (max) matrix count */

	if(RSB_SOME_ERROR(errval))
		goto err;

	tat = -(dt = rsb_time());
#if RSB_ALLOW_EMPTY_MATRICES
	if(!RSB_HAVE_GOOD_PARMS_FOR_EMPTY(m,k,nnz,flags))
#endif /* RSB_ALLOW_EMPTY_MATRICES */
	if(!RSB_HAVE_GOOD_PARMS_FOR_IN_PLACE_RCSR(m,k,nnz,flags))
	{
	       	errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_MDNFARTS);
	}

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC) && RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_BFSAH);
	}

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	{
	       	errval = RSB_ERR_UNIMPLEMENTED_YET;
		RSB_PERR_GOTO(err,RSB_ERRM_CMOINIY);
	}

#if 1
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,"Unexpected RSB_FLAG_FORTRAN_INDICES_INTERFACE flags here!\n");
	}
#else /* */
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE)) /* TODO: this is *slow*, speed this up */
		rsb__util_coo_array_from_fortran_indices(IA,nnz,RSB_BOOL_TRUE),
		rsb__util_coo_array_from_fortran_indices(JA,nnz,RSB_BOOL_TRUE),
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
#endif /* */

	if(RSB__VERBOSE_COO2REC)
		RSB_STDOUT("Building a matrix with %ld nnz, %ld x %ld\n",(long int)nnz,(long int)m,(long int)k);

	if(RSB_DO_FLAGS_EXTRACT_STORAGE(flags)==RSB_FLAG_NOFLAGS)
	{
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_DEFAULT_STORAGE_FLAGS);
	}

	/* TODO: may plug *here* upcoming RSB_WANT_FASTER_EXPERIMENTAL_CONSTRUCTOR stuff */

	mat = -dt;
	submatrices = rsb__calloc(sizeof(struct rsb_mtx_t )*tmc);
	submatricesp = rsb__calloc(sizeof(struct rsb_mtx_t*)*tmc);

	if(!submatrices || !submatricesp)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	mat += (dt = rsb_time());

	for(smi=0;smi<tmc;++smi)
		submatricesp[smi] = submatrices + smi;

	ect = -dt;
	if((errval = rsb__do_cleanup_nnz(VA,IA,JA,nnz,roff,coff,m,k,&nnz,typecode,flags))!=RSB_ERR_NO_ERROR)
		goto err;
	ect += (dt = rsb_time());

	est = -dt;
	if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_SORTED_INPUT))
	if((errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,m,k,typecode,flags))!=RSB_ERR_NO_ERROR)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT); /* TODO: is this needed ? */
	est += (dt = rsb_time()); /* with 'sorting' (est) we DO NOT intend also cleanup (in ect) */

	/* we need duplicates removal, and this can only take place after sorting */
	drt = -dt;
	dnnz = nnz - rsb__weed_out_duplicates(IA,JA,VA,nnz,typecode,flags);
	nnz -= dnnz;
	drt += (dt = rsb_time());
	if( RSB__VERBOSE_COO2REC )
		RSB_INFO("Duplicates check: %zd - %zd = %zd\n",(size_t)(nnz+dnnz),(size_t)dnnz,(size_t)nnz);


	/* work vectors allocation */
/*	IL = rsb__malloc(sizeof(rsb_coo_idx_t)*(m+1)); */
	mat -= dt;
	IT = rsb__malloc(sizeof(rsb_coo_idx_t)*(m+1));
	IX = rsb__malloc(sizeof(rsb_coo_idx_t)*2*(m+1));
	IB = rsb__malloc(sizeof(rsb_coo_idx_t)*(m+1));
	mat += (dt = rsb_time());
	if(/*  !IL ||*/ !IT || !IX || !IB)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	/* declaring this first matrix (smi == 0) as 'open' */
	smi = 0; omc = 1;
	/* compile in the first mtxAp, linking into to the temporary split vector */
	submatrices[smi].nzoff = 0;
	submatrices[smi].roff = roff;
	submatrices[smi].coff = coff;
	submatrices[smi].bindx = IB;
	submatrices[smi].bpntr = IB+1;
	submatrices[smi].indptr = NULL;
/*	RSB_DO_FLAG_ADD(flags,RSB_FLAG_QUAD_PARTITIONING); */
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS);
	if( (errval = rsb__set_init_flags_and_stuff(
		&submatrices[smi],NULL,NULL,m,k,nnz,nnz,nnz,typecode,flags)
				)!=RSB_ERR_NO_ERROR)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	
	if(nnz==0)
	{
		/* a special case. we copy the arrays addresses because they may be non-NULL and containing duplicate/diagonal/etc.. elements we have honoured to free, afterwards. */
		++cmc; --omc;
		mtxAp = &submatrices[0];
		RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_QUAD_PARTITIONING); /* necessary, too */
		mtxAp->bpntr = IA;
		mtxAp->bindx = JA;
		mtxAp->indptr = NULL;
		mtxAp->VA = VA;
		goto arrays_done;
	}
	else
		mtxAp = &submatrices[0];

	sat = -(dt = rsb_time());
	/* computing the first right-left pointer vectors */
#if RSB_WANT_PARALLEL_SUBDIVISION 
	if((errval = rsb_do_compute_vertical_split_parallel(IA,JA,roff,coff,m,k,0,0,nnz,IB,NULL,NULL,NULL,NULL,NULL,NULL))!=RSB_ERR_NO_ERROR)
#else /* RSB_WANT_PARALLEL_SUBDIVISION */
	if((errval = rsb_do_compute_vertical_split_parallel(IA,JA,roff,coff,m,k,0,0,nnz,IB,NULL,NULL,NULL,NULL,NULL,NULL))!=RSB_ERR_NO_ERROR)
	/*if((errval = rsb_do_compute_vertical_split(IA,JA,roff,coff,m,k,0,0,nnz,IB,NULL,NULL,NULL,NULL,NULL,NULL))!=RSB_ERR_NO_ERROR) */
#endif /* RSB_WANT_PARALLEL_SUBDIVISION */
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_C2R_IF_VERBOSE)
		RSB_INFO("beginning (%zd x %zd) @ %p with flags 0x%x (coo:%d, csr:%d), storage: 0x%x, max %zd submatrices\n",
			(rsb_printf_int_t)submatrices[smi].nr, (rsb_printf_int_t)submatrices[smi].nc, (const void*)&submatrices[smi], submatrices[smi].flags,
			RSB_DO_FLAG_HAS(submatrices[smi].flags,RSB_FLAG_WANT_COO_STORAGE),
			RSB_DO_FLAG_HAS(submatrices[smi].flags,RSB_FLAG_WANT_BCSS_STORAGE),
			submatrices[smi].matrix_storage,(rsb_printf_int_t)tmc
			);

/*	if(!RSB_WANT_MORE_PARALLELISM || (RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_QUAD_PARTITIONING))) */ /* TODO */
#if 1 /* the code is not yet ready for this */
/* #if RSB_WANT_OMP_RECURSIVE_KERNELS */
	errval = rsb_do_coo2rec_subdivide_parallel(VA,IA,JA,m,k,nnz,typecode,pinfop,flags,errvalp,submatricesp,mtxAp,IB,IX,IT,WA,cmc,omc,tmc,RSB_MAX(1,RSB_MIN(wet,nnz)),&cmc);
#else
	errval = rsb_do_coo2rec_subdivide(VA,IA,JA,m,k,nnz,typecode,pinfop,flags,errvalp,submatricesp,mtxAp,IB,IX,IT,WA,cmc,omc,tmc,wet,&cmc);
#endif
	sat += (dt = rsb_time());

	RSB_CONDITIONAL_FREE(IX);
	if(RSB_SOME_ERROR(errval))
		goto err;

	/*
	RSB_CONDITIONAL_FREE(IL);
	RSB_CONDITIONAL_FREE(IT);
       	*/
	/* WA will is needed for shuffle only (so, after IL,IM deallocation, in a way that total memory need is max(WA,IL)) */
	mat -= dt;
	WA = rsb__malloc(RSB_MAX(sizeof(rsb_coo_idx_t),el_size)*nnz);
	if(!WA)
	{
	       	errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_FAOTAFS);
	}
	mat += (dt = rsb_time());

	eit = -dt;
#if 0
	/* after symbolic partitioning is done, we are ready to shuffle all of the arrays using the temporary storage and add an intermediate node */
	RSB_INFO("assembling leaf %d -> %d x %d, %d\n",smi,submatricesp[smi]->nr,submatricesp[smi]->nc,submatricesp[smi]->nnz);
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_shuffle_left_and_right_rows(VA,IA,JA,m,0,nnz,0,typecode,IL,IM,WA));
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_shuffle_left_and_right_rows(VA,IA,JA,(m+1)/2,0,IL[(m+1)/2],0,typecode,IL,IR,WA));
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_shuffle_left_and_right_rows(VA,IA,JA,m,(m+1)/2,nnz,IL[(m+1)/2],typecode,IL,IR,WA));
	/* TODO : should use a temporary vector, here. */
#endif

	for(smi=0;smi<cmc;++smi)
		if(rsb__is_terminal_recursive_matrix(submatricesp[smi]))
			++lm;

/*	qsort(submatricesp+(cmc-lm),(size_t)(lm),sizeof(struct rsb_mtx_t*),&rsb__compar_rcsr_matrix_leftmost_first); */
	qsort(submatricesp,(size_t)(cmc),sizeof(struct rsb_mtx_t*),&rsb__compar_rcsr_matrix_leftmost_first);
	/* TODO: a priority queue would do the job, here */
	for(smi=0;smi<cmc-lm;++smi)
		if(rsb__is_terminal_recursive_matrix(submatricesp[smi]))
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ANLSMIT);
		}
	for(smi=cmc-lm;smi<cmc;++smi)
		if(!rsb__is_terminal_recursive_matrix(submatricesp[smi]))
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ALSMINT);
		}

	errval = rsb_do_coo2rec_shuffle(VA,IA,JA,m,k,nnz,typecode,pinfop,flags,errvalp,submatricesp,mtxAp,IB,WA,cmc);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_SEOWS);
	}

	rsb__do_set_in_place_submatrices_offsets(submatrices,cmc,VA,IA,JA,el_size);
	
/*	RSB_INFO("VA:%p, IA:%p, JA:%p\n",VA,IA,JA); */

	if(RSB_C2R_IF_VERBOSE)
		RSB_INFO("IA? :%p / %p\n",(const void*)IA,
			(const void*)(rsb__do_get_first_submatrix(mtxAp)->bpntr-
			(rsb__do_get_first_submatrix(mtxAp)->nr+1))
/*			rsb__do_get_first_submatrix(mtxAp)->roff-
  			((submatricesp[0])->nr+1) */
			);

	/* after shuffling, the last vectors conversion happens and we are done. */
arrays_done:
	eit += (dt = rsb_time());

	lst -= dt;
	/* first matrix is always root (even if a CSR one) */
	RSB_DO_FLAG_DEL(submatrices[0].flags,RSB_FLAG_NON_ROOT_MATRIX);
	#if RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS
	if(!(submatrices[0].flags & RSB_FLAG_NON_ROOT_MATRIX))
	{
		submatrices[0].all_leaf_matrices = NULL;
		errval = rsb__get_array_of_leaf_matrices(&submatrices[0],&(submatrices[0].all_leaf_matrices),&submatrices[0].all_leaf_matrices_n);
		if(RSB_SOME_ERROR(errval))
			goto err;
	}
	else
	{
		/* this is a non root matrix */
		submatrices[0].all_leaf_matrices = NULL;
		submatrices[0].all_leaf_matrices_n = 0;
	}
	#endif /* RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS */
	lst += (dt = rsb_time());

	hct = -dt;
	if(
		/* RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES_COO)
		       ||	RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES_CSR)
		       ||	RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COO_STORAGE)*/
		       	RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES)
		)
#if RSB_WANT_MORE_PARALLELISM 
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_switch_fresh_recursive_matrix_to_halfword_storages_parallel(mtxAp));
#else /* RSB_WANT_MORE_PARALLELISM */
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_switch_fresh_recursive_matrix_to_halfword_storages(mtxAp));
#endif /* RSB_WANT_MORE_PARALLELISM */
	else
	{
		if(RSB_C2R_IF_VERBOSE)
			RSB_INFO("no  RSB_FLAG_USE_HALFWORD_INDICES flag\n");
	}
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_compute_bounded_boxes(mtxAp));

	hct += rsb_time();

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	tat += rsb_time();

	mtxAp->sat = sat;
	mtxAp->ect = ect;
	mtxAp->eit = eit;
	mtxAp->est = est;
	mtxAp->pet = RSB_TIME_ZERO;
	mtxAp->rpt = RSB_TIME_ZERO;
	mtxAp->tat = tat;
	mtxAp->cpt = RSB_TIME_ZERO;/* cpt is contained in sat, so should not be counted here! */

	if( RSB__VERBOSE_COO2REC )
	{
		RSB_INFO(" converted COO to RSB in %.3le s (%.2lf %%)\n",tat,RSB_PCT(tat,tat));
		RSB_INFO(" analyzed arrays in %.3le s (%.2lf %%)\n",sat,RSB_PCT(sat,tat));
		RSB_INFO(" cleaned-up arrays in %.3le s (%.2lf %%)\n",ect,RSB_PCT(ect,tat));
		RSB_INFO(" deduplicated arrays in %.3le s (%.2lf %%)\n",drt,RSB_PCT(drt,tat));
		RSB_INFO(" sorted arrays in %.3le s (%.2lf %%)\n",est,RSB_PCT(est,tat));
		if( mtxAp->cpt )
			RSB_INFO(" computed partitions in %.3le s (%.2lf %%)\n",mtxAp->cpt,RSB_PCT(mtxAp->cpt,tat));
		RSB_INFO(" shuffled partitions in %.3le s (%.2lf %%)\n",eit,RSB_PCT(eit,tat));
		RSB_INFO(" memory allocations took %.3le s (%.2lf %%)\n",mat,RSB_PCT(mat,tat));
		RSB_INFO(" leafs setup took %.3le s (%.2lf %%)\n",lst,RSB_PCT(lst,tat));
		RSB_INFO(" halfword conversion took %.3le s (%.2lf %%)\n",hct,RSB_PCT(hct,tat));
	}
#if RSB_STORE_IDXSA
	mtxAp->idxsa = rsb__get_index_storage_amount(mtxAp);
#endif

#if 0
  	RSB_INFO("got %d matrices (%d leaves)\n",cmc,lm);

  	if( !rsb__mtx_chk(mtxAp) )
  	{
  		errval = RSB_ERR_INTERNAL_ERROR;
  		RSB_PERR_GOTO(derr,RSB_ERRM_NL);
  	}

	errval = rsb__do_switch_recursive_matrix_to_fullword_storage(mtxAp);

	errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo(mtxAp);
ok:
#endif
	goto noerr;
err:
	mtxAp = NULL;
	RSB_CONDITIONAL_FREE(submatrices);
noerr:
	if(RSB_SOME_ERROR(errval))
		rsb__do_perror(NULL,errval);
	RSB_CONDITIONAL_FREE(IB);
	RSB_CONDITIONAL_FREE(IT);
	RSB_CONDITIONAL_FREE(IX);
/*	RSB_CONDITIONAL_FREE(IL); */
	RSB_CONDITIONAL_FREE(WA);
	RSB_CONDITIONAL_FREE(submatricesp);
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	if( ( RSB__VERBOSE_COO2REC ) && mtxAp )
	{
		RSB_STDOUT("Built ");
		RSB_STDOUT_MATRIX_SUMMARY(mtxAp);
		RSB_STDOUT("\n");
	}
	return mtxAp;
}

rsb_err_t rsb__project_rsb_to_coo(const struct rsb_mtx_t *mtxAp, struct rsb_coo_mtx_t *coop)
{
	// copy parms and (shallow) pointers from rsb to coo
	const rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_INIT_COO_FROM_MTX(coop,mtxAp);
	RSB_BIND_COO_TO_MTX(coop,mtxAp);

	return errval;
}

/* @endcond */
