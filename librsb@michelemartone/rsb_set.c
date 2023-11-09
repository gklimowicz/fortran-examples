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
 * @author Michele Martone
 * @brief
 * Matrix setter/getter functions.
 * */

#include "rsb_internals.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

#define RSB_ALLOW_MTX_UPD 1

static const void * rsb_do_has_coo_element_inner(const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t i, rsb_coo_idx_t j)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Return a pointer to the value at (i,j), if present, and NULL if not present.
	 * Only for CSR/CSC.
	 * */
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(!RSB_INVALID_COO_INDEX(i));
	RSB_DEBUG_ASSERT(!RSB_INVALID_COO_INDEX(j));
	RSB_DEBUG_ASSERT(rsb__is_css_matrix(mtxAp));

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		const struct rsb_mtx_t * submatrix = RSB_FIND_SUBMATRIX_CONTAINING(mtxAp,i+mtxAp->roff,j+mtxAp->coff);
		rsb_coo_idx_t moff;
		rsb_coo_idx_t koff;

		if(!submatrix)
			return NULL;
		moff = submatrix->roff-mtxAp->roff;
		koff = submatrix->coff-mtxAp->coff;
		return rsb_do_has_coo_element_inner(submatrix,i-moff,j-koff);
	}
	else
	{
		rsb_nnz_idx_t si;
		rsb_coo_idx_t Mi,mi;
		const rsb_nnz_idx_t * bpntr = mtxAp->bpntr;
		const rsb_coo_idx_t * bindx = mtxAp->bindx;

		if( mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
			Mi = j, mi = i;
		else
			Mi = i, mi = j;

		if(rsb__is_coo_matrix(mtxAp))
		{
			rsb_nnz_idx_t nnz1,nnz0,nnz = mtxAp->nnz;

			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				// delimit the current row
				const rsb_half_idx_t *IA = (const rsb_half_idx_t *)mtxAp->bpntr;
				const rsb_half_idx_t *JA = (const rsb_half_idx_t *)mtxAp->bindx;
				// we search the beginning of line Mi
				nnz0 = rsb__nnz_split_hcoo_bsearch(IA,Mi,nnz);
				// we search the end of line Mi
				nnz1 = nnz0+rsb__nnz_split_hcoo_bsearch(IA+nnz0,Mi+1,nnz-nnz0);
				if(nnz1-nnz0<1)
					return NULL;// no row Mi
				// in line Mi, we search the index of mi, if any
				nnz0 += rsb__nnz_split_hcoo_bsearch(JA+nnz0,mi,nnz1-nnz0);
				if(nnz1-nnz0<1 || JA[nnz0]!=mi)
					return NULL;// no element mi			
			}
			else
			{
				// delimit the current row
				const rsb_coo_idx_t *IA = mtxAp->bpntr;
				const rsb_coo_idx_t *JA = mtxAp->bindx;
				// we search the beginning of line Mi
				nnz0 = rsb__nnz_split_coo_bsearch(IA,Mi,nnz);
				// we search the end of line Mi
				nnz1 = nnz0+rsb__nnz_split_coo_bsearch(IA+nnz0,Mi+1,nnz-nnz0);
				if(nnz1-nnz0<1)
					return NULL;// no row Mi
				// in line Mi, we search the index of mi, if any
				nnz0 += rsb__nnz_split_coo_bsearch(JA+nnz0,mi,nnz1-nnz0);
				if(nnz1-nnz0<1 || JA[nnz0]!=mi)
					return NULL;// no element mi
			}
			return ((const rsb_char_t*)(mtxAp->VA))+mtxAp->el_size*(nnz0);
		}
		else
		{
			if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_CSR)
			       	si = rsb__seek_half_idx_t(((rsb_half_idx_t*)bindx)+bpntr[Mi],mi,bpntr[Mi+1]-bpntr[Mi]);
			else
			       	si = rsb__seek_nnz_idx_t(bindx+bpntr[Mi],mi,bpntr[Mi+1]-bpntr[Mi]);
		}

		if(si == RSB_MARKER_NNZ_VALUE)
			return NULL;
		else
			return ((const rsb_char_t*)(mtxAp->VA))+mtxAp->el_size*(bpntr[Mi]+si);
	}
}

const void * rsb__do_coo_element_inner_address(const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t i, rsb_coo_idx_t j)
{
	return rsb_do_has_coo_element_inner(mtxAp,i,j);
}

rsb_err_t rsb__do_set_coo_elements(struct rsb_mtx_t * mtxAp, const void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz)
{
	/*!
	 * \ingroup gr_internals
	 * Will continue updating even on error.
	 * Ignore diagonal-implicit matrix update attempts (no error returned).
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t n;

	if(nnz<rsb_global_session_handle.rsb_want_threads)
	{
		for(n=0;RSB_LIKELY(n<nnz);++n)
		{
			rsb_err_t lerrval = rsb__do_set_coo_element(mtxAp,((rsb_char_t*)VA)+mtxAp->el_size*n,IA[n],JA[n]);

			RSB_DO_ERROR_CUMULATE(errval,lerrval);
#if RSB_OUT_ERR_VERBOSITY<=2
			RSB__ERR_FLAG_DEL(lerrval, RSB__ERR_CANTUPDATE_DIAGI);
#endif /* RSB_OUT_ERR_VERBOSITY */
			if(RSB_SOME_ERROR(lerrval))
			{
			       	RSB_ERROR("error updating %dth element of %d: %d %d\n",n,nnz,IA[n],JA[n]);
				// RSB_PERR_GOTO(err,RSB_ERRM_ES)
			}
		}
	}
	else
	{
		#pragma omp parallel for schedule(static,1) reduction(|:errval)  RSB_NTC
		for(n=0;n<nnz;++n)
		{
			rsb_err_t lerrval = rsb__do_set_coo_element(mtxAp,((rsb_char_t*)VA)+mtxAp->el_size*n,IA[n],JA[n]);

			RSB_DO_ERROR_CUMULATE(errval,lerrval);
#if RSB_OUT_ERR_VERBOSITY<=2
			RSB__ERR_FLAG_DEL(lerrval, RSB__ERR_CANTUPDATE_DIAGI);
#endif /* RSB_OUT_ERR_VERBOSITY */
			if(RSB_SOME_ERROR(lerrval))
			{
			       	RSB_ERROR("error updating %dth element of %d: %d %d\n",n,nnz,IA[n],JA[n]);
				// RSB_PERR_GOTO(err,RSB_ERRM_ES)
			}
		}
	}
//err:
	if ( errval & RSB__ERR_CANTUPDATE_DIAGI )
	{
		RSB__ERR_FLAG_DEL(errval, RSB__ERR_CANTUPDATE_DIAGI);
#if RSB_OUT_ERR_VERBOSITY>=3
		RSB_ERROR("Ignoring diagonal-implicit matrix update attempt\n");
#endif /* RSB_OUT_ERR_VERBOSITY */
	}
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_set_coo_element(struct rsb_mtx_t * mtxAp, const void * vp, const rsb_coo_idx_t i, const rsb_coo_idx_t j)
{
	return rsb__do_upd_coo_element(mtxAp, vp, i, j, RSB_FLAG_DUPLICATES_DEFAULT_HANDLE);
}

rsb_err_t rsb__do_upd_coo_element(struct rsb_mtx_t * mtxAp, const void * vp, const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Overwrite element at (i,j), if present.
	 * If not present, return RSB_ERR_GENERIC_ERROR.
	 * If matrix is diagonal implicit, return RSB__ERR_CANTUPDATE_DIAGI.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void * OV = NULL;

	if(!mtxAp || !vp || RSB_INVALID_COO_INDEX(i) || RSB_INVALID_COO_INDEX(j))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(!rsb__is_css_matrix(mtxAp))
	{
		errval = RSB_ERR_UNIMPLEMENTED_YET;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	
	if(!RSB_MATRIX_CONTAINS(mtxAp,i,j))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if( i == j && rsb__get_diagonal_type_flag(mtxAp)==RSB_DIAGONAL_I )
	{
		errval = RSB__ERR_CANTUPDATE_DIAGI;
#if RSB_OUT_ERR_VERBOSITY>=3
		RSB_PERR_GOTO(err,RSB_ERRM_DICBU)
#else /* RSB_OUT_ERR_VERBOSITY */
		goto err;
#endif /* RSB_OUT_ERR_VERBOSITY */
	}

	OV = (void*) rsb_do_has_coo_element_inner(mtxAp,i,j);

	if(OV)
	{
#if RSB_ALLOW_MTX_UPD
		rsb_flags_t cflag = ( mtxAp->flags & RSB_FLAG_ALL_DUPLICATE_FLAGS )
		                  | (        flags & RSB_FLAG_ALL_DUPLICATE_FLAGS );

		if(RSB_DO_FLAG_HAS(cflag,RSB_FLAG_DUPLICATES_SUM))
		{ RSB_NUMERICAL_TYPE_SUM_AND_STORE_ELEMENTS(OV,vp,mtxAp->typecode);}
		else
#endif
		{ RSB_NUMERICAL_TYPE_SET_ELEMENT(OV,vp,mtxAp->typecode); }
	}
	else
		errval = RSB_ERR_ELEMENT_NOT_FOUND;

err:
	RSB_DO_ERR_RETURN(errval)
}

#if 1
static rsb_err_t rsb_do_locate_nnz_element(const struct rsb_mtx_t * mtxAp, void ** vpp, rsb_coo_idx_t*ip, rsb_coo_idx_t*jp, rsb_nnz_idx_t nzi)
{
	/*
		With rsb__do_get_nnz_element  used (experimentally) by sparsersb.
		Note: it honors Fortran indices flags, if any.
	*/
	rsb_err_t errval = RSB_ERR_BADARGS;
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;

	if(!mtxAp)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(rsb__is_terminal_recursive_matrix(mtxAp)
		       	&& ( nzi >= mtxAp->nzoff )
		       	&& ( nzi <  mtxAp->nzoff+mtxAp->nnz )
			)
	{
		rsb_nnz_idx_t lnz = nzi-mtxAp->nzoff;
		rsb_byte_t*OV = mtxAp->VA;
		OV += mtxAp->el_size * lnz;

		if( (!ip) || (!jp) )
			goto noij;

		if(rsb__is_coo_matrix(mtxAp))
		{
			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
			{
				RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(mIA,mJA,mtxAp)
				i = mIA[lnz];
				j = mJA[lnz];
			}
			else
			{
				RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(mIA,mJA,mtxAp)
				i = mIA[lnz];
				j = mJA[lnz];
			}
		}
		else
		{
			/*	Examples:
 
				For low triangular 2x2 looking for element 'lnz==0', one 
				would return row 0.

				0: [0] -> [0  ]
				1: [1] -> [0,1]
				2: [3] -> []

				For low triangular 2x2 looking for element 'lnz==1', one 
				would return row 1.

				0: [0] -> [0  ]
				1: [1] -> [0,1]
				2: [3] -> []
   
				For low triangular 2x2 looking for element 'lnz==2', one 
				would return row 1.

				0: [0] -> [0  ]
				1: [1] -> [0,1]
				2: [3] -> []

			*/
			if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_CSR)
			{
				RSB_DECLARE_CONST_HALFCSR_ARRAYS_FROM_MATRIX(mPA,mJA,mtxAp);
				j = mJA[lnz];
				i = rsb__nnz_split_nnz_bsearch(mPA,lnz,mtxAp->nr);
				if( i>0 && mPA[i] > lnz )
					--i; // when lnz not first on its row
			}
			else
			{
				RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(mPA,mJA,mtxAp);
				j = mJA[lnz];
				i = rsb__nnz_split_nnz_bsearch(mPA,lnz,mtxAp->nr);
				if( i>0 && mPA[i] > lnz )
					--i; // when lnz not first on its row
			}
		}
		if(ip)
			*ip = i+mtxAp->roff;;
		if(jp)
			*jp = j+mtxAp->coff;;
noij:
		if(vpp)*vpp = OV;
		errval = RSB_ERR_NO_ERROR;
	}
	else
	{
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix
			       	&& ( nzi >= submatrix->nzoff ) 
			       	&& ( nzi <  submatrix->nnz+submatrix->nzoff )
		  )
		{
		       	errval = rsb_do_locate_nnz_element(submatrix,vpp,ip,jp,nzi);
			if(RSB_SOME_ERROR(errval))
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_nnz_element(const struct rsb_mtx_t * mtxAp, void * vp, rsb_coo_idx_t*ip, rsb_coo_idx_t*jp, rsb_nnz_idx_t nzi)
{
	/*
		This was used by sparsersb until 1.0.6. and known as rsb_do_get_nnz_element.
		Note: it honors Fortran indices flags, if any.
	*/
	void * OV = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	errval = rsb_do_locate_nnz_element(mtxAp,&OV,ip,jp,nzi);

	if((!RSB_SOME_ERROR(errval)) && OV)
	{
		RSB_NUMERICAL_TYPE_SET_ELEMENT(vp,OV,mtxAp->typecode);
	}
	RSB_DO_ERR_RETURN(errval)
}
#endif

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb_do_set_nnz_element(const struct rsb_mtx_t * mtxAp, const void * vp, rsb_nnz_idx_t nzi)
{
	/*
		Note: it honors Fortran indices flags, if any.
		Unused for now.
	 */
	void * OV = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	errval = rsb_do_locate_nnz_element(mtxAp,&OV,NULL,NULL,nzi);
	if((!RSB_SOME_ERROR(errval)) && OV)
	{
		RSB_NUMERICAL_TYPE_SET_ELEMENT(OV,vp,mtxAp->typecode);
	}
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__do_get_coo_element(const struct rsb_mtx_t * mtxAp, void * vp, rsb_coo_idx_t i, rsb_coo_idx_t j)
{
	/*!
	 * \ingroup gr_internals
	 * Gets the element at (i,j), if present.
	 * If not present, returns RSB_ERR_GENERIC_ERROR and zeros the area.
	 * In case of a 0x0 matrix (or blank area), it won't write anything, and return success.
	 * */
	const void * OV = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if RSB_ALLOW_ZERO_DIM
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
		goto err; /* Note: skip error checks */
#endif
	if(!mtxAp || !vp || RSB_INVALID_COO_INDEX(i) || RSB_INVALID_COO_INDEX(j))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(!rsb__is_css_matrix(mtxAp))
	{
		errval = RSB_ERR_UNIMPLEMENTED_YET;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	
	if(!RSB_MATRIX_CONTAINS(mtxAp,i,j))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if( i == j && rsb__get_diagonal_type_flag(mtxAp)==RSB_DIAGONAL_I )
		return rsb__fill_with_ones(vp,mtxAp->typecode,1,1);

	OV = rsb_do_has_coo_element_inner(mtxAp,i,j);
	if(!OV)
	{
		errval = RSB_ERR_ELEMENT_NOT_FOUND;
		rsb__cblas_Xscal(mtxAp->typecode,1,NULL,vp,1);
	}
	else
		{RSB_NUMERICAL_TYPE_SET_ELEMENT(vp,OV,mtxAp->typecode);}

err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_reverse_odd_rows(struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \rsb_warn_unoptimized_msg
	 * \note Works only for csr leaves.
	 * */
	RSB_DEBUG_ASSERT(mtxAp);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		struct rsb_mtx_t * submatrix;
		rsb_submatrix_idx_t i,j;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
			rsb__do_reverse_odd_rows(submatrix);
	}
	else
	{
		if(rsb__is_coo_matrix(mtxAp))
		{
			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				//RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
				/* FIXME: not implemented yet for COO */
			}
			else
			{
				//RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
				/* FIXME: not implemented yet for COO */
			}
		}
		else
		{
			rsb_coo_idx_t i;
			rsb_nnz_idx_t ib,ie;
			if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				RSB_DECLARE_HALFCSR_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
				for(i=1;i<mtxAp->nr;i+=2)
					ib = IA[i],ie = IA[i+1],
					rsb__util_reverse_halfword_coo_array(JA+ib,ie-ib);
			}
			else
			{
				RSB_DECLARE_FULLCSR_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
				for(i=1;i<mtxAp->nr;i+=2)
					ib = IA[i],ie = IA[i+1],
					rsb__util_reverse_fullword_coo_array(JA+ib,ie-ib);
			}
		}
	}
	return RSB_ERR_NO_ERROR; 
}

rsb_err_t rsb__do_zsort_coo_submatrices(struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \rsb_warn_unoptimized_msg
	 * */
	RSB_DEBUG_ASSERT(mtxAp);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		struct rsb_mtx_t * submatrix;
		rsb_submatrix_idx_t i,j;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
    			if(submatrix)
				rsb__do_zsort_coo_submatrices(submatrix);
	}
	else
	{
		if(rsb__is_coo_matrix(mtxAp))
		{
			rsb_err_t errval = RSB_ERR_NO_ERROR;
			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				RSB_DECLARE_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_to_fullword_zcoo(mtxAp));
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_index_based_z_morton_sort(NULL,NULL,NULL,IA,JA,mtxAp->VA,mtxAp->nr,mtxAp->nc,mtxAp->nnz,mtxAp->typecode,RSB_OP_FLAG_DEFAULT));
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_switch_to_halfword_coo(mtxAp));
			}
			else
			{
				RSB_DECLARE_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_index_based_z_morton_sort(NULL,NULL,NULL,IA,JA,mtxAp->VA,mtxAp->nr,mtxAp->nc,mtxAp->nnz,mtxAp->typecode,RSB_OP_FLAG_DEFAULT));
			}
			return errval;
		}
		else
		{
			if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				/* nothing to do */
			}
			else
			{
				/* nothing to do */
			}
		}
	}
	return RSB_ERR_NO_ERROR; 
}

/* @endcond */
