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
 * This source file contains functions for sparse matrices multiplication.
 */
/* FIXME : UNFINISHED, UNCHECKED, UNSECURED, preliminar code
 * TODO : spscatter
 *  	  should parallelize this code (not difficult; based on nproc accumulator arrays, one external omp loop)
 * */

#include "rsb_common.h"
#include "rsb_clone.h"

#define RSB_WANT_SPGEMM_MFLOPS 1
#define RSB_WANT_SPGEMM_VERBOSE (RSB_WANT_SPGEMM_MFLOPS&&0)  
#define RSB_WANT_SPGEMM_VERBOSE_PROGRAM 1
#define WANT_SPGEMM_FULL 1
// take care that having more threads will not give you benefits if memory access is slow 
#define RSB_WANT_OMP_RECURSIVE_SPGEMM_KERNELS RSB_WANT_OMP_RECURSIVE_KERNELS

#if 0
/* 20120930	unused ?! */
static rsb_err_t rsb_spgemm_sym_count_blocks(const struct rsb_mtx_t * mtxAp, const struct rsb_mtx_t * mtxBp, rsb_nnz_idx_t * cblocksp)
{
	/**
	 * \ingroup gr_unfinished
	 * Counts an upper bound estimate of the output blocks of A * B conformant sparse matrices product.
	 * It is a valid estimate for the expected work.
	 * Matrices must be pairwise CSR and CSC formats or BCSR br x bc and BCSC bc x br formats to conform.
	 * Matrices should be unsymmetric (although no error would be issues on symmetric matrices as input).
	 *
	 * ...
	 * */
	rsb_nnz_idx_t cblocks=0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t bi,bj;
	rsb_coo_idx_t cm,ck;

	if(!mtxAp || !mtxBp || !cblocksp || !rsb_are_matrices_spgemm_block_conformant(mtxAp,mtxBp))
	{errval = RSB_ERR_BADARGS;goto err;}
	
	cm=mtxAp->nr;
	ck=mtxBp->nc;

	for(bi=0;bi<cm;++bi)
	for(bj=0;bj<ck;++bj)
	{
		rsb_nnz_idx_t aro=mtxAp->bpntr[bi];
		rsb_nnz_idx_t are=mtxAp->bpntr[bi+1];
//		rsb_nnz_idx_t arb=mtxAp->bpntr[bi+1] - mtxAp->bpntr[bi];
		rsb_nnz_idx_t bro=mtxBp->bpntr[bj];
		rsb_nnz_idx_t bre=mtxBp->bpntr[bj+1];
//		rsb_nnz_idx_t bcb=mtxBp->bpntr[bj+1] - mtxBp->bpntr[bj];
		rsb_nnz_idx_t al=0,bl=0;

		if(/*  arb<bcb*/ 1 )
		{
			for(al=aro,bl=bro;al<are&&bl<bre;)
			{
				while(mtxBp->bindx[bl]<mtxAp->bindx[al] && bl<bre)
					++bl;

				/* TODO : a fast binary search codelet here */

				while(mtxBp->bindx[bl]>mtxAp->bindx[al] && al<are)
					++al;

				if(mtxAp->bindx[al]==mtxBp->bindx[bl])
					++cblocks,++al;
			}
		}
		else
		{
			/* FIXME : write me */
		}
	}
	
	*cblocksp=cblocks;
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif

RSB_INTERNALS_COMMON_HEAD_DECLS

static rsb_err_t rsb_spgemm_inner(const struct rsb_coo_mtx_t * acoo, const struct rsb_coo_mtx_t * bcoo, rsb_nnz_idx_t * cblocksp, rsb_coo_idx_t ** PA, rsb_coo_idx_t ** JA, void ** VA, size_t * opsp)
{
	/**
	 * \ingroup gr_unfinished
	 * Counts an upper bound estimate of the output blocks of A * B conformant sparse matrices product.
	 * It is a valid estimate for the expected work.
	 * Matrices must be pairwise CSR formats or CSC formats to conform.
	 * Matrices should be unsymmetric (although no error would be issues on symmetric matrices as input).
	 *
	 * CSR-CSR or CSC-CSC (and swapped A/B)
	 * ...
	 * */
	rsb_nnz_idx_t cblocks=0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t bi;
	rsb_coo_idx_t ai,aj;
//	rsb_coo_idx_t ci,cj;
	rsb_coo_idx_t cm,ck;
	/*rsb_coo_idx_t am,ak;*/
	rsb_coo_idx_t /*bm,*/bk;
	rsb_coo_idx_t al,bl;
	size_t el_size=0;
	rsb_coo_idx_t*abpntr=NULL,*bbpntr=NULL;
	rsb_coo_idx_t*abindx=NULL,*bbindx=NULL;

	void *acc=NULL;
	rsb_nnz_idx_t * p=NULL;
#if RSB_WANT_SPGEMM_MFLOPS
	rsb_nnz_idx_t ops=0;
	rsb_time_t t=0;
#endif /* RSB_WANT_SPGEMM_MFLOPS */

	if(!acoo || !bcoo || !cblocksp /* || !rsb_are_matrices_spgemm_block_conformant(acoo,bcoo)*/)
	{errval = RSB_ERR_BADARGS;goto err;}
	el_size = RSB_SIZEOF(acoo->typecode);

	abpntr=acoo->IA;
	bbpntr=bcoo->IA;
	abindx=acoo->JA;
	bbindx=bcoo->JA;

#if RSB_WANT_SPGEMM_MFLOPS
	t = - rsb_time();
#endif /* RSB_WANT_SPGEMM_MFLOPS */

	cm=acoo->nr;
	/*bm=bcoo->nr;*/
	bk=bcoo->nc;
	ck=bcoo->nc;
	/*ak=acoo->nc;*/

	acc = rsb__calloc(el_size*bk);
	if(!acc) {errval = RSB_ERR_ENOMEM;goto err;}
	p = rsb__calloc(sizeof(rsb_nnz_idx_t)*(bk+1));
	if(!p) {errval = RSB_ERR_ENOMEM;goto err;}

	if(PA && JA && VA)
		*PA = rsb__calloc(sizeof(rsb_coo_idx_t)*(cm+1));

	if(!*PA)
	{
		RSB_CONDITIONAL_FREE(*PA);
		errval = RSB_ERR_ENOMEM;
		goto err;
	}

	for(ai=0;ai<cm;++ai)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		/* rsb_nnz_idx_t arb; */
		rsb_nnz_idx_t marker;
		rsb_nnz_idx_t rcblocks=0;

		marker=cblocks+1;
		aro=abpntr[ai];
		are=abpntr[ai+1];
		/* arb=abpntr[ai+1] - abpntr[ai]; */
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro;
			rsb_nnz_idx_t bre;
			/* rsb_nnz_idx_t bcb; */
			aj=abindx[al];
			bi=aj;
			bro=bbpntr[bi];
			bre=bbpntr[bi+1];
			/* bcb=bbpntr[bi+1] - bbpntr[bi]; */
			for(bl=bro;bl<bre;++bl)
			{
				rsb_coo_idx_t bj;
				bj=bbindx[bl];
			//	RSB_STDOUT("(%d %d) x (%d %d)\n",ai,aj,bi,bj);
				if(p[bj]<marker)
					p[bj]=marker,
					rcblocks++;
				else
					;
			}
		}
		if(!RSB_IS_VALID_NNZ_SUM(cblocks,rcblocks))
		{
			RSB_ERROR("number of matrices product nnz may exceed maximum allowed.\n");
			errval = RSB_ERR_LIMITS;
			goto err;
		}
		cblocks += rcblocks;
		(*PA)[ai+1]=cblocks;
	}
	if(cblocksp)
		*cblocksp=cblocks;

	/* FIXME: and what to do if cblocks == 0 ? */
	if(PA && JA && VA)
	{
		*JA = rsb__calloc(sizeof(rsb_coo_idx_t)*cblocks);
		*VA = rsb__calloc(el_size*cblocks);
		if(!*PA || !*JA || !*VA)
		{
			RSB_CONDITIONAL_FREE(*PA);
			RSB_CONDITIONAL_FREE(*JA);
			RSB_CONDITIONAL_FREE(*VA);
			errval = RSB_ERR_ENOMEM;
			goto err;
		}
#if RSB_WANT_OMP_RECURSIVE_SPGEMM_KERNELS
	#pragma omp parallel RSB_NTC
	{
		const rsb_thread_t nt = omp_get_num_threads();
		const rsb_thread_t tn = omp_get_thread_num();
		rsb_coo_idx_t bj;
		void *acc_=NULL; rsb_nnz_idx_t *p_=NULL;
		#pragma omp critical (spgemm_alloc)
	{
#if RSB_WANT_SPGEMM_VERBOSE
		//if(!tn)RSB_STDOUT("parallel SPGEMM with %d threads\n",nt);
#endif /* RSB_WANT_SPGEMM_VERBOSE */
		if(tn)
		{
			acc_=rsb__calloc(el_size*bk);
			if(!acc_) {errval = RSB_ERR_ENOMEM;}
			p_=rsb__calloc(sizeof(rsb_nnz_idx_t)*(bk+1));
			if(!p_) {errval = RSB_ERR_ENOMEM;}
		}
		else
		{
			acc_=acc;p_=p;
			for(bj=0;bj<bk;++bj) p_[bj]=0;
		}
	}
		if(RSB_SOME_ERROR(errval))
			goto ierr;
		/* FIXME: ops in parallel section ! */
		rsb__do_util_csr_csr_sparse_mul_serial(*PA,*JA,*VA,abpntr,bbpntr,abindx,bbindx,acoo->VA,bcoo->VA,cm,ck,p_,acc_,&ops,acoo->typecode,tn,nt);
ierr:
		#pragma omp critical (spgemm_alloc)
		{ if(tn) { RSB_CONDITIONAL_FREE(acc_);RSB_CONDITIONAL_FREE(p_); } }
	}
#else /* RSB_WANT_OMP_RECURSIVE_SPGEMM_KERNELS */
	{
		rsb_coo_idx_t bj;
		for(bj=0;bj<bk;++bj)
			p[bj]=0;
		
		//rsb__do_util_csr_csr_sparse_mul_serial(*PA,*JA,*VA,acoo->bpntr,bcoo->bpntr,acoo->bindx,bcoo->bindx,acoo->VA,bcoo->VA,cm,ck,p,acc,&ops,acoo->typecode,0,1);
		rsb__do_util_csr_csr_sparse_mul_serial(*PA,*JA,*VA,acoo->IA,bcoo->IA,acoo->JA,bcoo->JA,acoo->VA,bcoo->VA,cm,ck,p,acc,&ops,acoo->typecode,0,1);
	}
#endif /* RSB_WANT_OMP_RECURSIVE_SPGEMM_KERNELS */
	}
	*PA = rsb__realloc(*PA,RSB_MAX(cblocks,RSB_MAX(1+cm,ck))*sizeof(rsb_nnz_idx_t));
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_compressed_array_to_fullword_coo(*PA,cm,0,NULL));
#if RSB_WANT_SPGEMM_MFLOPS
	ops*=2;
	t += rsb_time();
#endif /* RSB_WANT_SPGEMM_MFLOPS */
#if RSB_WANT_SPGEMM_VERBOSE
#if RSB_ALLOW_STDOUT
	RSB_STDOUT("%zd nonzeros , %zd x %zd\n",(size_t)cblocks,cm,ck);
	RSB_STDOUT("%g MFLOPS\n",((double)(ops))/(t*1000000));
#endif /* RSB_ALLOW_STDOUT */
#endif /* RSB_WANT_SPGEMM_VERBOSE */
err:
	RSB_CONDITIONAL_FREE(acc);
	RSB_CONDITIONAL_FREE(p);
	if(opsp)*opsp=ops;
	if(RSB_SOME_ERROR(errval))
	{
		RSB_CONDITIONAL_FREE(*PA);
		RSB_CONDITIONAL_FREE(*JA);
		RSB_CONDITIONAL_FREE(*VA);
	}
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_spgemm_dense_inner(rsb_coo_idx_t ldc, rsb_coo_idx_t nr, rsb_coo_idx_t nc, rsb_bool_t isccolmajor, void *cVA_, const struct rsb_coo_mtx_t * acoo, const struct rsb_coo_mtx_t * bcoo, rsb_nnz_idx_t * cblocksp, size_t * opsp)
{
	/**
	 * \ingroup gr_unfinished
	 * */
	//rsb_nnz_idx_t cblocks=0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//rsb_coo_idx_t bi;
	//rsb_coo_idx_t ai,aj;
	rsb_coo_idx_t cm,ck;
	//rsb_coo_idx_t /*am,*/ak;
	//rsb_coo_idx_t bm,bk;
	//rsb_coo_idx_t al,bl;
	//size_t el_size=0;

#if RSB_WANT_SPGEMM_MFLOPS
	rsb_nnz_idx_t ops=0;
	rsb_time_t t=0;
#endif /* RSB_WANT_SPGEMM_MFLOPS */
	rsb_coo_idx_t*abpntr=NULL,*bbpntr=NULL;
	rsb_coo_idx_t*abindx=NULL,*bbindx=NULL;

	if(!acoo || !bcoo || !cblocksp /* || !rsb_are_matrices_spgemm_block_conformant(acoo,bcoo)*/)
	{errval = RSB_ERR_BADARGS;goto err;}
	//el_size = RSB_SIZEOF(acoo->typecode);

	if( (isccolmajor && (ldc<nr)) ||  ((!isccolmajor) && (ldc<nc))
	 || (acoo->nc != bcoo->nr) || (acoo->nr > nr) || (bcoo->nc > nc) 
	  )
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}
#if RSB_WANT_SPGEMM_MFLOPS
	t = - rsb_time();
#endif /* RSB_WANT_SPGEMM_MFLOPS */

	abpntr=acoo->IA;
	bbpntr=bcoo->IA;
	abindx=acoo->JA;
	bbindx=bcoo->JA;

	cm=acoo->nr;
	//bm=bcoo->nr;
	//bk=bcoo->nc;
	ck=bcoo->nc;
	//ak=acoo->nc;

	{
		rsb_nnz_idx_t opss=0;
#if RSB_WANT_OMP_RECURSIVE_SPGEMM_KERNELS
#if RSB_WANT_SPGEMM_MFLOPS
		//size_t opss=0;
#else /* RSB_WANT_SPGEMM_MFLOPS */
		//size_t opss=NULL;
#endif /* RSB_WANT_SPGEMM_MFLOPS */
	#pragma omp parallel reduction(+:opss) RSB_NTC
	{
		const rsb_thread_t nt = omp_get_num_threads();
		const rsb_thread_t tn = omp_get_thread_num();
		//rsb_coo_idx_t bj;
		rsb__do_util_csr_csr_dense_mul_serial(ldc,nr,nc,isccolmajor,cVA_,abpntr,bbpntr,abindx,bbindx,acoo->VA,bcoo->VA,cm,ck,&opss,acoo->typecode,tn,nt);
	}
#else /* RSB_WANT_OMP_RECURSIVE_SPGEMM_KERNELS */
	{
		rsb__do_util_csr_csr_dense_mul_serial(ldc,nr,nc,isccolmajor,cVA_,abpntr,bbpntr,abindx,bbindx,acoo->VA,bcoo->VA,cm,ck,&opss,acoo->typecode,0,1);
	}
#endif /* RSB_WANT_OMP_RECURSIVE_SPGEMM_KERNELS */
		ops+=opss;
	}

#if RSB_WANT_SPGEMM_MFLOPS
	ops*=2;
	t += rsb_time();
#endif /* RSB_WANT_SPGEMM_MFLOPS */
#if RSB_WANT_SPGEMM_VERBOSE
#if RSB_ALLOW_STDOUT
	//RSB_STDOUT("%zd nonzeros , %zd x %zd\n",(size_t)cblocks,cm,ck);
	RSB_STDOUT("%g MFLOPS\n",((double)(ops))/(t*1000000));
#endif /* RSB_ALLOW_STDOUT */
#endif /* RSB_WANT_SPGEMM_VERBOSE */
	if(opsp)*opsp=ops;
err:
	RSB_DO_ERR_RETURN(errval)
}

static struct rsb_mtx_t * rsb_spgemm_tmp(rsb_type_t typecode, const struct rsb_mtx_t * mtxAp, const struct rsb_mtx_t * mtxBp, rsb_trans_t transA, rsb_trans_t transB, const  void *alphap, const void *betap, rsb_err_t *errvalp, rsb_time_t * dtp, size_t * opsp)
{
	/**
	 * \ingroup gr_unfinished
	 * FIXME: requires CSR format everywhere !
	 * TODO need to-CSR cloning, here!
	 * FIXME: missing a smart approach for handling symmetric matrices (current policy is explicit expansion) 
	 * */
	// FIXME: WON'T WORK FOR SOME VERY SPARSE MATRICES
	void * VA = NULL;
	rsb_coo_idx_t * IA = NULL, *JA = NULL;
	rsb_nnz_idx_t nnz = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS;
	struct rsb_coo_mtx_t acsr,bcsr;
	struct rsb_mtx_t *mtxCp = NULL;
	rsb_time_t dt = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_coo_idx_t cm = 0,ck = 0;
	rsb_coo_idx_t tak = 0, tbm = 0;

	RSB_BZERO_P(&acsr);
	RSB_BZERO_P(&bcsr);

	if(RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}

	if( !mtxAp || !mtxBp )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"Supplied a NULL matrix pointer.\n");
	}
	RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);

#if 0
	if( mtxAp->typecode != mtxBp->typecode )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"Matrix types do not match.\n");
	}
#else
#endif

	cm = RSB_MTX_TRANSPOSED_ROWS(mtxAp,transA);
	ck = RSB_MTX_TRANSPOSED_COLS(mtxBp,transB);
	tak = RSB_MTX_TRANSPOSED_COLS(mtxAp,transA);
	tbm = RSB_MTX_TRANSPOSED_ROWS(mtxBp,transB);

	if( (transA != RSB_TRANSPOSITION_N) || (transB != RSB_TRANSPOSITION_N) )
	{
		errval = RSB_ERR_UNIMPLEMENTED_YET;
		RSB_PERR_GOTO(err,"Transposition parameter not yet supported !\n");
	}

	if(tak!=tbm)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"Matrix sizes do not match.\n");
	}

	acsr.nr = mtxAp->nr;
	acsr.nc = mtxAp->nc;
	acsr.nnz = RSB_MAX(mtxAp->nnz,RSB_MAX(mtxAp->nr+1,mtxAp->nc+1)); /* FIXME: temporary !*/
	acsr.nnz += mtxAp->nnz+RSB_MIN(mtxAp->nr+1,mtxAp->nc+1); /* FIXME: temporary, in case of symmetry & diagonal !*/
	acsr.typecode = typecode;
	if((rsb__allocate_coo_matrix_t(&acsr)!=&acsr))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,"problem converting the A matrix\n");
	}
	//acsr.nnz=mtxAp->nnz;
	bcsr.nr = mtxBp->nr;
	bcsr.nc = mtxBp->nc;
	bcsr.nnz = RSB_MAX(mtxBp->nnz,RSB_MAX(mtxBp->nr+1,mtxBp->nc+1)); /* FIXME: temporary !*/
	bcsr.nnz += mtxBp->nnz+RSB_MIN(mtxBp->nr+1,mtxBp->nc+1); /* FIXME: temporary, in case of symmetry & diagonal !*/
	bcsr.typecode = typecode;
	if((rsb__allocate_coo_matrix_t(&bcsr)!=&bcsr))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,"problem converting the B matrix\n");
	}
	//bcsr.nnz=mtxBp->nnz;
	errval = rsb__do_get_csr(typecode,mtxAp,acsr.VA,acsr.IA,acsr.JA,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS);
	if(RSB_SOME_ERROR(errval))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,"csr extraction problems from matrix A\n");
	}
	errval = rsb__do_get_csr(typecode,mtxBp,bcsr.VA,bcsr.IA,bcsr.JA,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS);
	if(RSB_SOME_ERROR(errval))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,"csr extraction problems from matrix B\n");
	}
	acsr.nnz=acsr.IA[acsr.nr];
	bcsr.nnz=bcsr.IA[bcsr.nr];

	if(dtp)dt = - rsb_time();
	if((errval = rsb_spgemm_inner(&acsr,&bcsr,&nnz,&IA,&JA,&VA,opsp))!=RSB_ERR_NO_ERROR)
	/* FIXME: warning: allocation size may not be max(nnz,m) in mtxCp arrays, now ! */
	{
		rsb__do_perror(NULL,errval);
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	if(dtp)dt += rsb_time();
	//else rsb__test_print_coo_mm(mtxAp->typecode,mtxAp->flags,IA,JA,VA,mtxAp->nr,mtxBp->nc,nnz,RSB_BOOL_TRUE,RSB_DEFAULT_STREAM);

	mtxCp = rsb__do_mtx_alloc_from_coo_inplace(VA,IA,JA,nnz,typecode,cm,ck,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags,&errval);
	if(mtxCp)RSB_DO_FLAG_DEL(mtxCp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
	if(!mtxCp)
	{
		RSB_ERROR("Failed allocating a matrix\n.");
		if(errval == RSB_ERR_NO_ERROR)
			errval = RSB_ERR_INTERNAL_ERROR;
		RSB_CONDITIONAL_FREE(IA);
		RSB_CONDITIONAL_FREE(JA);
		RSB_CONDITIONAL_FREE(VA);
		goto err;
	}
	//else rsb__print_matrix_unsorted_coo(mtxCp);
err:
	if(dtp)*dtp=dt;
	rsb__destroy_coo_matrix_t(&acsr);
	rsb__destroy_coo_matrix_t(&bcsr);
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxCp;
}

rsb_err_t rsb__do_spgemm_to_dense(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_coo_idx_t ldc, rsb_coo_idx_t nr, rsb_coo_idx_t nc, rsb_bool_t isccolmajor, void *cVA, rsb_time_t * dtp, size_t * opsp)
{
	/**
	 * \ingroup gr_unfinished
	 * TODO: missing input validation
	 * FIXME: WON'T WORK FOR SOME VERY SPARSE MATRICES
	 * FIXME: alphap,betap,transA,transB still unused
	 * */
	rsb_nnz_idx_t nnz=0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS;
	struct rsb_coo_mtx_t acsr,bcsr;
	rsb_time_t dt=RSB_CONST_IMPOSSIBLY_BIG_TIME;

	RSB_BZERO_P(&acsr);
	RSB_BZERO_P(&bcsr);

	if( RSB_INVALID_COO_INDEX(nr) || RSB_INVALID_COO_INDEX(nc) || RSB_INVALID_COO_INDEX(ldc) || !mtxAp || !mtxBp
#if 0
	|| (mtxAp->typecode != mtxBp->typecode)
#endif
       	)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	if( (transA != RSB_TRANSPOSITION_N) || (transB != RSB_TRANSPOSITION_N) )
	{
		RSB_ERROR("Transposition parameter not yet supported !\n");
		errval = RSB_ERR_UNIMPLEMENTED_YET;
		goto err;
	}
	acsr.nr=mtxAp->nr;
	acsr.nc=mtxAp->nc;
	acsr.nnz = RSB_MAX(mtxAp->nnz,RSB_MAX(mtxAp->nr+1,mtxAp->nc+1)); /* FIXME: temporary !*/
	//acsr.typecode=mtxAp->typecode;
	acsr.typecode=typecode;
	if((rsb__allocate_coo_matrix_t(&acsr)!=&acsr))
	{
		RSB_ERROR("allocaton problem\n");
       		errval = RSB_ERR_INTERNAL_ERROR; goto err;
	}
	acsr.nnz=mtxAp->nnz;
	bcsr.nr=mtxBp->nr;
	bcsr.nc=mtxBp->nc;
	bcsr.nnz = RSB_MAX(mtxBp->nnz,RSB_MAX(mtxBp->nr+1,mtxBp->nc+1)); /* FIXME: temporary !*/
	//bcsr.typecode=mtxBp->typecode;
	bcsr.typecode=typecode;
	if((rsb__allocate_coo_matrix_t(&bcsr)!=&bcsr))
	{
		RSB_ERROR("allocaton problem\n");
       		errval = RSB_ERR_INTERNAL_ERROR; goto err;
	}
	bcsr.nnz=mtxBp->nnz;
	errval = rsb__do_get_csr(typecode,mtxAp,acsr.VA,acsr.IA,acsr.JA,flags);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_ERROR("csr extraction problems from matrix A\n");
	      	errval = RSB_ERR_INTERNAL_ERROR; goto err;
	}
	errval = rsb__do_get_csr(typecode,mtxBp,bcsr.VA,bcsr.IA,bcsr.JA,flags);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_ERROR("csr extraction problems from matrix B\n");
	      	errval = RSB_ERR_INTERNAL_ERROR; goto err;
	}

	if(dtp)dt = - rsb_time();
	if((errval = rsb_do_spgemm_dense_inner(ldc,nr,nc,isccolmajor,cVA,&acsr,&bcsr,&nnz,opsp))!=RSB_ERR_NO_ERROR)
	{
		RSB_ERROR(RSB_ERRM_ES);
		goto err;
	}
	if(dtp)dt += rsb_time();

err:
	if(dtp)*dtp=dt;
	rsb__destroy_coo_matrix_t(&acsr);
	rsb__destroy_coo_matrix_t(&bcsr);
	return errval;
}

rsb_err_t rsb__do_spgemm_test_code(const int argc, char * const argv[])
{
	/**
	 * \ingroup gr_unfinished
	 * FIXME : temporary, testing code */
#if WANT_SPGEMM_FULL
	struct rsb_mtx_t *mtxCp = NULL;
#else /* WANT_SPGEMM_FULL */
	void * CVA=NULL; rsb_coo_idx_t * CIA=NULL,*CJA=NULL;
#endif /* WANT_SPGEMM_FULL */
	//rsb_coo_idx_t m=0,k=0,nnz=0;
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	//rsb_flags_t flags = 0;
	//rsb_flags_t flags = RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS;		 	
	rsb_flags_t flags = RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const char * filename=NULL;
	const char * cfilename=NULL;
	struct rsb_mtx_t *mtxAp = NULL;
	struct rsb_mtx_t *mtxBp = NULL;
	//int br,bc;
	/* 4x4 and 1x8 blockings give various results, as a numerical side effect */
	//rsb_nnz_idx_t cblocks=0;
	rsb_trans_t transA = RSB_TRANSPOSITION_N;
	rsb_trans_t transB = RSB_TRANSPOSITION_N;
	//rsb_aligned_t errnorm[];
#if 0
	rsb_aligned_t alpha[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t beta[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_byte_t * alphap=alpha;
	rsb_byte_t * betap=beta;
#else
	rsb_byte_t * alphap=NULL;
	rsb_byte_t * betap=NULL;
#endif
	rsb_time_t rsb_spg_time = RSB_CONST_IMPOSSIBLY_BIG_TIME,csr_spg_time = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_time_t t = RSB_CONST_IMPOSSIBLY_BIG_TIME;
	rsb_thread_t nt = rsb_get_num_threads(); 
	size_t ops=0;
	//rsb_option_t options[] = { RSB_BENCH_PROG_OPTS {0,0,0,0} };
	//br=1,bc=8; br=4,bc=4; br=1,bc=1;

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORT_INPUT);
	filename = "pd.mtx";

	if(argc>=2)
		filename=argv[1];

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
		goto err;

	if((errval = rsb__do_load_matrix_file_as_matrix_market(&mtxAp,filename,flags,typecode))!=RSB_ERR_NO_ERROR)
		goto err;

	if(argc>=3 && strcmp(argv[2],argv[1]))
		errval = rsb__do_load_matrix_file_as_matrix_market(&mtxBp,argv[2],flags,typecode);
	else
		mtxBp=mtxAp;

	if(argc>=4)
		cfilename = argv[3];
	else
		cfilename = "pd.tmp.mtx";

	t = - rsb_time();
	if((mtxCp = rsb_spgemm_tmp(typecode,mtxAp,mtxBp,transA,transB,alphap,betap,&errval,&csr_spg_time,&ops))==NULL)
		goto err;

#if RSB_WANT_SPGEMM_VERBOSE
	//RSB_STDOUT("%zd nonzeros, %g s (%g Mnnz/s)\n",(size_t)cblocks,t,((double)(cblocks))/(t*1000000));
#endif /* RSB_WANT_SPGEMM_VERBOSE */

err:
	t += rsb_time();
	rsb_spg_time=t;
#if RSB_WANT_SPGEMM_VERBOSE_PROGRAM
#if RSB_ALLOW_STDOUT
	RSB_STDOUT_MATRIX_SUMMARY(mtxAp);
	RSB_STDOUT("\n * \n");
	RSB_STDOUT_MATRIX_SUMMARY(mtxBp);
	RSB_STDOUT("\n = \n");
	RSB_STDOUT_MATRIX_SUMMARY(mtxCp);
	RSB_STDOUT("\n");
	RSB_STDOUT("%%:CSR_SPGEMM_PERFORMANCE:");RSB_STDOUT_MATRIX_ESSENTIALS(mtxCp,filename,nt);RSB_STDOUT("\t%10.6lf\n",(RSB_FPINV(csr_spg_time)*ops)/RSB_MILLION_I );
	RSB_STDOUT("%%:RSB_SPGEMM_PERFORMANCE:");RSB_STDOUT_MATRIX_ESSENTIALS(mtxCp,filename,nt);RSB_STDOUT("\t%10.6lf\n",(RSB_FPINV(rsb_spg_time)*ops)/RSB_MILLION_I );
	RSB_STDOUT("%%:CSR_SPGEMM_TIME:");RSB_STDOUT_MATRIX_ESSENTIALS(mtxCp,filename,nt);RSB_STDOUT("\t%10.6lf\n",(csr_spg_time));
	RSB_STDOUT("%%:RSB_SPGEMM_TIME:");RSB_STDOUT_MATRIX_ESSENTIALS(mtxCp,filename,nt);RSB_STDOUT("\t%10.6lf\n",(rsb_spg_time));
	RSB_STDOUT("%%:SPGEMM_OPS:");RSB_STDOUT_MATRIX_ESSENTIALS(mtxCp,filename,nt);RSB_STDOUT("\t%zd\n",ops);
#endif /* RSB_ALLOW_STDOUT */
#endif /* RSB_WANT_SPGEMM_VERBOSE_PROGRAM */

	rsb__do_perror(NULL,errval);
	if(mtxBp && mtxBp!=mtxAp)
		RSB_MTX_FREE(mtxBp);
	if(mtxAp)
	       	RSB_MTX_FREE(mtxAp);

	errval = rsb_file_mtx_save(mtxCp, cfilename);
	rsb__do_perror(NULL,errval);

#if WANT_SPGEMM_FULL
	if(mtxCp)
		RSB_MTX_FREE(mtxCp)
#else /* WANT_SPGEMM_FULL */
	RSB_CONDITIONAL_FREE(CIA);
	RSB_CONDITIONAL_FREE(CJA);
	RSB_CONDITIONAL_FREE(CVA);
#endif /* WANT_SPGEMM_FULL */

	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS)) != RSB_ERR_NO_ERROR)
		;

	RSB_DO_ERR_RETURN(errval)
}

struct rsb_mtx_t * rsb__do_matrix_mul(rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_trans_t transB, const void *betap, const struct rsb_mtx_t * mtxBp, rsb_err_t * errvalp)
{
	return rsb_spgemm_tmp(typecode,mtxAp,mtxBp,transA,transB,alphap,betap,errvalp,NULL,NULL);
}

/* @endcond */
