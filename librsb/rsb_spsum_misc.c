/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
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
#include "rsb_common.h"


rsb_err_t rsb__do_add_submatrix_to_dense(const struct rsb_mtx_t * mtxAp, const void *alphap, void * Bp, rsb_nnz_idx_t ldb, rsb_nnz_idx_t nr, rsb_nnz_idx_t nc, rsb_bool_t rowmajor)
{
	if(!mtxAp || !Bp || !alphap )
		goto err;

#ifdef RSB_NUMERICAL_TYPE_DOUBLE 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE  )
	{
{
	rsb_nnz_idx_t n;
	const rsb_coo_idx_t roff=mtxAp->roff, coff=mtxAp->coff;
	const double *VA=mtxAp->VA;

	if(rsb__is_coo_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		{
			RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
		}
		else
		{
			RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
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
					*(double*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(double*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
			}
		}
		else
		{
			RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			if(rowmajor)
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(double*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(double*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double*)alphap)*(VA[n]);
			}		
		}
	}
	else
		RSB_ERROR(RSB_ERRM_NL);
}

	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT  )
	{
{
	rsb_nnz_idx_t n;
	const rsb_coo_idx_t roff=mtxAp->roff, coff=mtxAp->coff;
	const float *VA=mtxAp->VA;

	if(rsb__is_coo_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		{
			RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
		}
		else
		{
			RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
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
					*(float*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(float*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
			}
		}
		else
		{
			RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			if(rowmajor)
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(float*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(float*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float*)alphap)*(VA[n]);
			}		
		}
	}
	else
		RSB_ERROR(RSB_ERRM_NL);
}

	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  )
	{
{
	rsb_nnz_idx_t n;
	const rsb_coo_idx_t roff=mtxAp->roff, coff=mtxAp->coff;
	const float complex *VA=mtxAp->VA;

	if(rsb__is_coo_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		{
			RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
		}
		else
		{
			RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(float complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
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
					*(float complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(float complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
			}
		}
		else
		{
			RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			if(rowmajor)
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(float complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(float complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(float complex*)alphap)*(VA[n]);
			}		
		}
	}
	else
		RSB_ERROR(RSB_ERRM_NL);
}

	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
	if( mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  )
	{
{
	rsb_nnz_idx_t n;
	const rsb_coo_idx_t roff=mtxAp->roff, coff=mtxAp->coff;
	const double complex *VA=mtxAp->VA;

	if(rsb__is_coo_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		{
			RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
		}
		else
		{
			RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			if(rowmajor)
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
			else
			for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				*(double complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,IA[n]+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
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
					*(double complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(double complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
			}
		}
		else
		{
			RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			if(rowmajor)
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(double complex*)(RSB_BLOCK_ROWMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
			}
			else
			{
				for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
				for(n=PA[i];RSB_LIKELY(n<PA[i+1]);++n)
					*(double complex*)(RSB_BLOCK_COLMAJOR_ADDRESS(Bp,ldb,nr,nc,i+roff,JA[n]+coff,mtxAp->el_size))+=(*(double complex*)alphap)*(VA[n]);
			}		
		}
	}
	else
		RSB_ERROR(RSB_ERRM_NL);
}

	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
	return RSB_ERR_UNSUPPORTED_TYPE;
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_GENERIC_ERROR;
}

/* @endcond */
