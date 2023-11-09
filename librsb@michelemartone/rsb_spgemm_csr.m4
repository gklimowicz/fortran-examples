/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains some CSR sparse matrices multiplication code.
 * */
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_SPGEMM_COO_H_INCLUDED
#define RSB_SPGEMM_COO_H_INCLUDED
#include "rsb_internals.h"
',`dnl
#include "rsb_internals.h"
')
dnl

rsb_err_t rsb__do_util_csr_csr_sparse_mul_serial(rsb_nnz_idx_t * PA, rsb_coo_idx_t * JA, void *VA_, const rsb_nnz_idx_t *ARP, const rsb_nnz_idx_t *BRP, const rsb_coo_idx_t *AJA, const rsb_coo_idx_t *BJA, const void * aVA_, const void * bVA_, const rsb_coo_idx_t cm, const rsb_coo_idx_t ck, rsb_nnz_idx_t * p, void * acc_, rsb_nnz_idx_t * opsp , rsb_type_t typecode, const rsb_coo_idx_t afr, const rsb_coo_idx_t ars)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	rsb_nnz_idx_t cblocks=0; 
	rsb_nnz_idx_t ops=0; 
	rsb_coo_idx_t ai,aj;
	rsb_coo_idx_t al,bl,cl;
	rsb_coo_idx_t bj;
	//rsb_coo_idx_t bi;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
dnl
	{
	mtype *VA=VA_,*acc=acc_;
	const mtype * aVA=aVA_,*bVA=bVA_; ;
dnl	for(ai=0;ai<cm;++ai)
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;
		rsb_nnz_idx_t marker;

		//assert(cblocks==PA[ai]);	// this is true on the serial execution of this loop
		cblocks=PA[ai];		// this shall work even in a parallel execution of this loop (with differing acc/p arrays)
		marker=cblocks+1;
		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];
				if(p[bj]<marker)
					p[bj]=marker,
					(JA)[cblocks++]=bj,
					acc[bj] =aVA[al]*bVA[bl];
				else
					acc[bj]+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}

		for(cl=(PA)[ai];cl<(PA)[ai+1];++cl)
		{
			((mtype*)(VA))[cl]=acc[(JA)[cl]];	/* FIXME */
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl

dnl

rsb_err_t rsb__do_util_csr_csr_dense_mul_serial(rsb_coo_idx_t ldc, rsb_coo_idx_t nr, rsb_coo_idx_t nc, rsb_bool_t isccolmajor, void *cVA_, const rsb_nnz_idx_t *ARP, const rsb_nnz_idx_t *BRP, const rsb_coo_idx_t *AJA, const rsb_coo_idx_t *BJA, const void * aVA_, const void * bVA_, const rsb_coo_idx_t cm, const rsb_coo_idx_t ck, rsb_nnz_idx_t * opsp , rsb_type_t typecode, const rsb_coo_idx_t afr, const rsb_coo_idx_t ars)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	rsb_nnz_idx_t ops=0; 
	rsb_coo_idx_t ai,aj;
	rsb_coo_idx_t al,bl;
	//rsb_coo_idx_t bi;
	rsb_coo_idx_t bj;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
dnl
	{
	mtype *cVA=cVA_;
	const mtype * aVA=aVA_,*bVA=bVA_; ;
dnl	for(ai=0;ai<cm;++ai)
	for(ai=afr;ai<cm;ai+=ars)
	{
		rsb_nnz_idx_t aro;
		rsb_nnz_idx_t are;
		//rsb_nnz_idx_t arb;

		aro=ARP[ai];
		are=ARP[ai+1];
		//arb=ARP[ai+1]-ARP[ai];
		/* we start row ai of target matrix C */
		for(al=aro;al<are;++al)
		{
			rsb_nnz_idx_t bro=BRP[aj=AJA[al]];
			rsb_nnz_idx_t bre=BRP[aj+1];
/*			rsb_nnz_idx_t bcb=BRP[aj+1] - BRP[aj];*/
			for(bl=bro;bl<bre;++bl)
			{
				//bi=aj;
				bj=BJA[bl];
dnl				*(mtype*)(RSB_BLOCK_ROWMAJOR_ADDRESS(cVA,ldc,nr,nc,ai,bj,(sizeof(mtype))))+=aVA[al]*bVA[bl];

				RSB_BLOCK_X_MAJOR_REFERENCE(cVA,ldc,ai,bj,isccolmajor)+=aVA[al]*bVA[bl];
			}
/*#if RSB_WANT_SPGEMM_MFLOPS*/
			ops+=(bre-bro);
/*#endif*/
		}
	}
	if(opsp)*opsp=ops;
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl

dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#endif /* RSB_SPGEMM_COO_H_INCLUDED */
')dnl
dnl
/* @endcond */
