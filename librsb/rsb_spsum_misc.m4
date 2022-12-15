/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * */
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_SPSUM_MISC_H_INCLUDED
#define RSB_SPSUM_MISC_H_INCLUDED
#include "rsb_common.h"
',`dnl
#include "rsb_common.h"
')
dnl

rsb_err_t rsb__do_add_submatrix_to_dense(const struct rsb_mtx_t * mtxAp, const void *alphap, void * Bp, rsb_nnz_idx_t ldb, rsb_nnz_idx_t nr, rsb_nnz_idx_t nc, rsb_bool_t rowmajor)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	if(!mtxAp || !Bp || !alphap )
		goto err;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( mtxAp->typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
dnl
	{
{
	rsb_nnz_idx_t n;
	const rsb_coo_idx_t roff=mtxAp->roff, coff=mtxAp->coff;
	const mtype *VA=mtxAp->VA;

	if(rsb__is_coo_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		{
			RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(mtype*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(mtype*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
		}
		else
		{
			RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(mtype*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(mtype*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
		}
	}
	else
	if(rsb__is_csr_matrix(mtxAp))
	{
		rsb_nnz_idx_t n,i;

		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
		{
			RSB_DECLARE_CONST_HALFCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			if(rowmajor)
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(mtype*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(mtype*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
			}
		}
		else
		{
			RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			if(rowmajor)
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(mtype*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(mtype*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(mtype*)alphap)*(VA[n]);
			}		
		}
	}
	else
		RSB_ERROR(RSB_ERRM_NL);
}

	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE;
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_GENERIC_ERROR;
}
')dnl

dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
#endif /* RSB_SPSUM_MISC_H_INCLUDED */
')dnl
dnl
/* @endcond */
