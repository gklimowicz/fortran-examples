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
 * This source file contains matrix getter functions.
 * */

#include "rsb_internals.h"

#define RSB_VA_MEMCPY(VD,VS,DOFF,SOFF,NNZ,ES) \
	RSB_MEMCPY(((rsb_char_t*)(VD))+(ES)*(DOFF),((rsb_char_t*)(VS))+(ES)*(SOFF),(ES)*(NNZ))

RSB_INTERNALS_COMMON_HEAD_DECLS

rsb_err_t rsb__do_get_coo(const struct rsb_mtx_t *mtxAp, rsb_byte_t ** VA, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, rsb_flags_t flags )
{
	/*! 
	 * \ingroup gr_internals
	 *  Returns the matrix converted in a coordinate storage format.
	 *
	 * \param VA  (optional) the values array pointer, sized at least for mtxAp->nnz elements of matrix type
	 * \param IA  (optional) an integer array pointer for row    coordinates
	 * \param JA  (optional) an integer array pointer for column coordinates
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * WARNING : If any of the given arrays points to NULL, it will be allocated accordingly.
	 *
	 * The entire matrix will be returned in COO format, in the specified VA,IA,JA arrays
	 * No more than mtxAp->nnz elements will be written to in the VA, IA and JA arrays
	 * TODO: according to flags, may sort according to rows or columns!
	 * Note: the filled array could result smaller;
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t el_size;
	rsb_nnz_idx_t nnz = 0;/* FIXME */

#if RSB_ALLOW_ZERO_DIM
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
		goto err; /* FIXME: skipping further error checks */
#endif

	if( !IA || !JA || !VA || !mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	el_size = mtxAp->el_size;

	/* Note : we do not allocate all arrays at once  (but do so in freeing)! */
	if(!*VA ) *VA= rsb__malloc(el_size     * mtxAp->nnz);
	if(!*IA ) *IA= rsb__malloc(sizeof(rsb_coo_idx_t) * mtxAp->nnz);
	if(!*JA ) *JA= rsb__malloc(sizeof(rsb_coo_idx_t) * mtxAp->nnz);

	if(!*VA || !*IA || !*JA)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if( (errval = rsb__do_get_coo_noalloc(mtxAp,*VA,*IA,*JA,&nnz,flags)) == RSB_ERR_NO_ERROR)
		return RSB_ERR_NO_ERROR;
	else
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
err:
	/* Note : we free all arrays at once ! */
	RSB_CONDITIONAL_FREE(*IA);
	RSB_CONDITIONAL_FREE(*JA);
	RSB_CONDITIONAL_FREE(*VA);
	RSB_DO_ERR_RETURN(errval)
}

/* Important: does not handle COO. */
rsb_err_t rsb__do_get_row_dense(const struct rsb_mtx_t *mtxAp , void * row , rsb_blk_idx_t rowindex)
{
	/*!
	 * \ingroup gr_internals
	 * Will write entire row rowindex of matrix \c mtxAp in the row vector.
	 *
	 * \param rowindex the specified row
	 * \param row an already allocated vector of the same type as \c mtxAp.
	 * \param matrix is a valid pointer to a rsb_mtx_t structure
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * \note This function is slow.
	 * \note This function is obsolete -- it's only used by ot.c !
	 * */
#if RSB_OBSOLETE_QUARANTINE_UNUSED
	register rsb_byte_t	*bp = NULL;
	register rsb_nnz_idx_t baserow,basecolumn;
	register rsb_blk_idx_t rows,columns;
	register rsb_blk_idx_t blockrow,blockcolumn;
#endif
	size_t el_size = mtxAp->el_size;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!row || rowindex <0 || rowindex >= mtxAp->nr)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		{
			rsb_coo_idx_t moff;
			rsb_coo_idx_t koff;
			if(submatrix)
			{
				moff = submatrix->roff-mtxAp->roff;
				koff = submatrix->coff-mtxAp->coff;
			}
		if(submatrix && rowindex >= moff && rowindex < submatrix->nr+submatrix->roff )
		{
			errval = rsb__do_get_row_dense(submatrix, ((rsb_byte_t*)row)+el_size*(j*koff) , rowindex-(i*moff));
			if(RSB_SOME_ERROR(errval))
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}}
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(rowindex));	
#if RSB_EXPERIMENTAL_USE_PURE_BCSS
	if(rsb__is_bcsr_matrix(mtxAp))
	{
		/* plain BCSR case */
		rsb_blk_idx_t br, bc;
		rsb_nnz_idx_t frb,lrb,bi;

		rsb__get_blocking_size(mtxAp, &br, &bc);
	
		frb = mtxAp->bpntr[ rowindex/br   ];
		lrb = mtxAp->bpntr[(rowindex/br)+1];
		for(bi=frb;RSB_LIKELY(bi<lrb);++bi)
		{
			/* FIXME : numerical overflow possible, in br*bc and alike */
			rsb__memcpy(
				((rsb_byte_t*)row)+el_size * mtxAp->bindx[bi]*bc,
				((rsb_byte_t*)mtxAp->VA)+(bi*br*bc + bc*(rowindex - (rowindex/br)*br))*el_size,
				el_size*bc
			);
		}
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	else
#if RSB_WANT_BCSC_LEAVES  
	if(rsb__is_bcsc_matrix(mtxAp))
	{
		rsb_blk_idx_t br, bc;
		rsb_nnz_idx_t bri, bci, bi;

		rsb__get_blocking_size(mtxAp, &br, &bc);
		bri= rowindex/br ;	/* the block row of interest */


		RSB_DEBUG_ASSERT(br>0 && bc>0);
		RSB_DEBUG_ASSERT(el_size);

		for(bci=0;RSB_LIKELY(bci<mtxAp->Mdim);++bci)
		{
			rsb_byte_t * dst = ((rsb_byte_t*)row) + el_size * bc * bci;

			if((bi = rsb__seek_nnz_idx_t(mtxAp->bindx+mtxAp->bpntr[bci+0],bri,mtxAp->bpntr[bci+1]-mtxAp->bpntr[bci+0]))!=RSB_MARKER_NNZ_VALUE)
			{
				rsb_nnz_idx_t boff = rowindex-br*bri;
				bi += mtxAp->bpntr[bci+0];
				RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(bci));
				RSB_DEBUG_ASSERT(bc );
				RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(bc*bci));
				rsb__memcpy(dst,((rsb_byte_t*)mtxAp->VA)+(bi*br*bc+bc*boff)*el_size,el_size*bc);
			}
		}
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
#endif /* RSB_WANT_BCSC_LEAVES */
#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS */

	{
	RSB_GET_FIRST_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
	while(!RSB_GOT_LAST_BLOCK_POINTER(mtxAp))
	{
		rsb_blk_idx_t c;
		if(rowindex >= baserow && rowindex <baserow+rows)
		for(c=0;c<columns;++c)
		{
			rsb_byte_t*src = (rsb_byte_t*)(bp+RSB_GET_INTRA_BLOCK_OFFSET(rowindex,basecolumn+c,blockrow,blockcolumn,mtxAp));
			rsb_byte_t*dst = ((rsb_byte_t*)row) + el_size * (basecolumn+c);
			RSB_NUMERICAL_TYPE_SET_ELEMENT(dst,src,mtxAp->typecode);	/* FIXME : SLOW */
		}
		RSB_GET_NEXT_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
	}
	}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
	errval = RSB_ERR_NO_ERROR;
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_rows_nnz(const struct rsb_mtx_t *mtxAp, rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_nnz_idx_t  *rnz)
{
        /*!
	 * \ingroup gr_internals
         *
	 * FIXME : NEW, TO DOCUMENT
         * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(fr));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(lr));
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(lr< mtxAp->roff+mtxAp->nr);
	RSB_DEBUG_ASSERT(fr>=mtxAp->roff);

	// 20101023 the following breaks `make tests`, for some uninvestigated reason
	//if(fr==0 && lr==mtxAp->nr-1)
	//	return mtxAp->nnz;

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix && RSB_SUBMATRIX_INTERSECTS_ROWS(submatrix,fr,lr))
		{
			const rsb_coo_idx_t fri = RSB_SUBMATRIX_ROWS_INTERSECTION_FIRST(submatrix,fr);
			const rsb_coo_idx_t lri = RSB_SUBMATRIX_ROWS_INTERSECTION_LAST(submatrix,lr);
			errval = rsb__do_get_rows_nnz(submatrix,fri,lri,rnz);
			if(RSB_SOME_ERROR(errval))
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
	}
	else
	{
		/* 
		 * leaf matrix processing
		 * */
		const rsb_coo_idx_t fri = RSB_SUBMATRIX_ROWS_INTERSECTION_FIRST(mtxAp,fr)-mtxAp->roff;
		const rsb_coo_idx_t lri = RSB_SUBMATRIX_ROWS_INTERSECTION_LAST(mtxAp,lr) -mtxAp->roff;

		if(rsb__is_coo_matrix(mtxAp))
		{	
			rsb_nnz_idx_t nnz1, nnz0, nnz = mtxAp->nnz;
			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
			{
				RSB_DECLARE_CONST_HALFCOO_IARRAY_FROM_MATRIX(IA,mtxAp)
				// we search the beginning of line fri
				nnz0 = rsb__nnz_split_hcoo_bsearch(IA,fri,nnz);
				// we search the end of line lri
				nnz1 = nnz0+rsb__nnz_split_hcoo_bsearch(IA+nnz0,lri+1,nnz-nnz0);
				*rnz += nnz1-nnz0;
			}
			else
			{
				RSB_DECLARE_CONST_FULLCOO_IARRAY_FROM_MATRIX(IA,mtxAp)
				// we search the beginning of line fri
				nnz0 = rsb__nnz_split_coo_bsearch(IA,fri,nnz);
				// we search the end of line lri
				nnz1 = nnz0+rsb__nnz_split_coo_bsearch(IA+nnz0,lri+1,nnz-nnz0);
				*rnz += nnz1-nnz0;
			}
		}
		else
		if(rsb__is_csr_matrix(mtxAp))
		{
			*rnz += mtxAp->bpntr[lri+1]-mtxAp->bpntr[fri];
		}
		else
			errval = RSB_ERR_UNIMPLEMENTED_YET;
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_columns_sparse(const struct rsb_mtx_t *mtxAp , void * RSB_RESTRICT VA , rsb_blk_idx_t fc, rsb_blk_idx_t lc, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT CP, rsb_coo_idx_t ioff, rsb_coo_idx_t joff)
{
        /*!
	 * \ingroup gr_internals
	 *
	 * JA can be NULL.
	 * fc, lc in global matrix space.
         * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(fc));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(lc));
	RSB_DEBUG_ASSERT(mtxAp);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
#if RSB_OBSOLETE_QUARANTINE_UNUSED
		RSB_DEBUG_ASSERT(lc< mtxAp->roff+mtxAp->nr);
		RSB_DEBUG_ASSERT(fc>=mtxAp->roff);
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix && RSB_SUBMATRIX_INTERSECTS_ROWS(submatrix,fc,lc))
		{
			const rsb_coo_idx_t fci = RSB_SUBMATRIX_COLS_INTERSECTION_FIRST(submatrix,fc);
			const rsb_coo_idx_t lci = RSB_SUBMATRIX_COLS_INTERSECTION_LAST(submatrix,lc);
			errval = rsb__do_get_columns_sparse(submatrix,VA,fci,lci,IA,JA,CP,ioff,joff);
			if(RSB_SOME_ERROR(errval))
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
		errval = RSB_ERR_BADARGS; RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	else
	{
		/* 
		 * leaf matrix processing
		 *
		 * FIXME: this code could (should?) be improved, regarding its performance
		 * */
		const rsb_coo_idx_t roff = mtxAp->roff;
		const rsb_coo_idx_t coff = mtxAp->coff;
		//const rsb_coo_idx_t fci = RSB_SUBMATRIX_COLS_INTERSECTION_FIRST(mtxAp,fc);
		//const rsb_coo_idx_t lci = RSB_SUBMATRIX_COLS_INTERSECTION_LAST(mtxAp,lc);
		rsb_nnz_idx_t nnz = mtxAp->nnz/*,dnnz = 0*/;
		register rsb_coo_idx_t i;
		register rsb_nnz_idx_t n;
		//const void * MVA = mtxAp->VA;
		const rsb_coo_idx_t * bindx = mtxAp->bindx;

		if(rsb__is_coo_matrix(mtxAp))
		{
#if 0
			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_COO)
			{
				for(n=0;n<mtxAp->nnz;++n)
				{
					rsb_coo_idx_t ij = bindx[n],j = RSB_COO_HALFWORDS_VALUES_UNPACK_UJ(ij);
					rsb_nnz_idx_t idx = CP[coff+j];
					i = RSB_COO_HALFWORDS_VALUES_UNPACK_LI(ij);
					RSB_VA_MEMCPY(VA,mtxAp->VA,idx,n,1,mtxAp->el_size);
						if(IA)IA[idx] = i+roff+ioff;
					if(JA)JA[idx] = j+coff+joff;
					CP[coff+j]++;
				}
			}
			else
#endif
			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(mIA,mJA,mtxAp)
				for(n=0;RSB_LIKELY(n<nnz);++n)
				{
					rsb_half_idx_t i = mIA[n],j = mJA[n];
					rsb_nnz_idx_t idx = CP[coff+j];
					RSB_VA_MEMCPY(VA,mtxAp->VA,idx,n,1,mtxAp->el_size);
					if(IA)IA[idx] = roff+ioff+i;
					if(JA)JA[idx] = coff+joff+j;
					CP[coff+j]++;
				}
			}
			else
			{
				RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(mIA,mJA,mtxAp)
				for(n=0;RSB_LIKELY(n<nnz);++n)
				{
					rsb_coo_idx_t i = mIA[n],j = mJA[n];
					rsb_nnz_idx_t idx = CP[coff+j];
					RSB_VA_MEMCPY(VA,mtxAp->VA,idx,n,1,mtxAp->el_size);
					if(IA)IA[idx] = roff+ioff+i;
					if(JA)JA[idx] = coff+joff+j;
					CP[coff+j]++;
				}
			}
		}
		else
		{
		if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_CSR)
		{
			const rsb_half_idx_t *hbindx = (const rsb_half_idx_t *)bindx;
			for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
			for(n=mtxAp->bpntr[i];RSB_LIKELY(n<mtxAp->bpntr[i+1]);++n)
			{
				rsb_coo_idx_t j = hbindx[n];// FIXME: is this mixed types assignment correct ?
				rsb_nnz_idx_t idx = CP[coff+j];
				RSB_VA_MEMCPY(VA,mtxAp->VA,idx,n,1,mtxAp->el_size);
				if(IA)IA[idx] = i+roff+ioff;
				if(JA)JA[idx] = j+coff+joff;
				CP[coff+j]++;
			}
		}
		else
		{
			for(i=0;RSB_LIKELY(i<mtxAp->nr);++i)
			for(n=mtxAp->bpntr[i];RSB_LIKELY(n<mtxAp->bpntr[i+1]);++n)
			{
				rsb_coo_idx_t j = bindx[n];
				rsb_nnz_idx_t idx = CP[coff+j];
//				RSB_INFO("%d %d %d %d %d %d\n",roff,coff,i,j,n,idx);
				RSB_VA_MEMCPY(VA,mtxAp->VA,idx,n,1,mtxAp->el_size);
				if(IA)IA[idx] = i+roff+ioff;
				if(JA)JA[idx] = j+coff+joff;
				CP[coff+j]++;
			}
//			if(roff>0)exit(-1);
		}}
	}
err:
        RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_rows_sparse_rec(const struct rsb_mtx_t *mtxAp , void * RSB_RESTRICT VA , rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t * RSB_RESTRICT rnz, rsb_coo_idx_t ioff, rsb_coo_idx_t joff)
{
        /*!
	 * \ingroup gr_internals
	 *
	 * IA can be NULL.
	 * FIXME: shall rewrite this to be faster.
         * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(fr));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(lr));
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(lr< mtxAp->roff+mtxAp->nr);
	RSB_DEBUG_ASSERT(fr>=mtxAp->roff);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix && RSB_SUBMATRIX_INTERSECTS_ROWS(submatrix,fr,lr))
		{
			const rsb_coo_idx_t fri = RSB_SUBMATRIX_ROWS_INTERSECTION_FIRST(submatrix,fr);
			const rsb_coo_idx_t lri = RSB_SUBMATRIX_ROWS_INTERSECTION_LAST(submatrix,lr);
			errval = rsb__do_get_rows_sparse_rec(submatrix,VA,fri,lri,IA,JA,rnz,ioff,joff);
			if(RSB_SOME_ERROR(errval))
			{
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
			}
		}
	}
	else
	{
		/* 
		 * leaf matrix processing
		 *
		 * FIXME: this code could (should?) be improved, regarding its performance
		 * */
		const rsb_coo_idx_t roff = mtxAp->roff;
		const rsb_coo_idx_t coff = mtxAp->coff;
		const rsb_coo_idx_t fri = RSB_SUBMATRIX_ROWS_INTERSECTION_FIRST(mtxAp,fr)-roff;
		const rsb_coo_idx_t lri = RSB_SUBMATRIX_ROWS_INTERSECTION_LAST (mtxAp,lr)-roff;
		rsb_nnz_idx_t nnz,dnnz = 0;
		const rsb_nnz_idx_t zoff = *rnz;
		rsb_nnz_idx_t doff = 0;
		register rsb_coo_idx_t i;
		const void * MVA = mtxAp->VA;

#define RSB_CSR2COO_MEMCPY_(VD,ID,JD,VS,I,JS,DOFF,SOFF,NNZ,ES,J0) \
	{	\
		if(ID)rsb__util_coo_array_set(((rsb_coo_idx_t*)(ID))+(DOFF),(NNZ),(I)); \
		RSB_IA_MEMCPY(JD,JS,DOFF,SOFF,NNZ,J0);	\
	}

		if(rsb__is_coo_matrix(mtxAp))
		{	
			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
			{
				rsb_nnz_idx_t nnz1, nnz0, nnz = mtxAp->nnz;
				RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(mIA,mJA,mtxAp)
				// we search the beginning of line fri
				nnz0 = rsb__nnz_split_hcoo_bsearch(mIA,fri,nnz);
				// we search the end of line lri
				nnz1 = nnz0+rsb__nnz_split_hcoo_bsearch(mIA+nnz0,lri+1,nnz-nnz0);
				nnz = nnz1-nnz0;
				if(nnz>0)
				{
					RSB_DEBUG_ASSERT( nnz <= mtxAp->nc );
					/* FIXME: need specialized little functions or macros, here */
					if(IA)
						for(i=nnz0;i<nnz1;++i)
							IA[zoff+i-nnz0] = mIA[i],
							IA[zoff+i-nnz0] += mtxAp->roff+ioff;
					for(i=nnz0;i<nnz1;++i)
						JA[zoff+i-nnz0] = mJA[i],
						JA[zoff+i-nnz0] += mtxAp->coff+joff;

					RSB_VA_MEMCPY(VA,MVA,zoff,nnz0,nnz,mtxAp->el_size);
					doff = nnz0;
					dnnz = nnz;
				}
				if(nnz<0)
				{
					errval = RSB_ERR_INTERNAL_ERROR;
					RSB_PERR_GOTO(err,RSB_ERRM_ES)
				}
			}
			else
			{
				rsb_nnz_idx_t nnz1,nnz0,nnz = mtxAp->nnz;
				RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(mIA,mJA,mtxAp)
				// we search the beginning of line fri
				nnz0 = rsb__nnz_split_coo_bsearch(mIA,fri,nnz);
				// we search the end of line lri
				nnz1 = nnz0+rsb__nnz_split_coo_bsearch(mIA+nnz0,lri+1,nnz-nnz0);
				nnz = nnz1-nnz0;
				if(nnz>0)
				{
					RSB_DEBUG_ASSERT( nnz <= mtxAp->nc );
					RSB_COA_MEMCPY(JA,mJA,zoff,nnz0,nnz);
					rsb__util_coo_array_add(JA+zoff,nnz,mtxAp->coff+joff);
					if(IA)
						RSB_COA_MEMCPY(IA,mIA,zoff,nnz0,nnz),
						rsb__util_coo_array_add(IA+zoff,nnz,mtxAp->roff+ioff);
					RSB_VA_MEMCPY(VA,MVA,zoff,nnz0,nnz,mtxAp->el_size);
					doff = nnz0;
					dnnz = nnz;
				}
				if(nnz<0)
				{
					errval = RSB_ERR_INTERNAL_ERROR;
					RSB_PERR_GOTO(err,RSB_ERRM_ES)
				}
//				RSB_INFO("COO OUT (%d..%d) (%d nnz) (@ %d)\n",fri,lri,nnz,doff);
			}
		}
		else /* csr ! FIXME: why was RSB_FLAG_USE_HALFWORD_INDICES_CSR being used ?? */
		{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		// if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES_CSR))
		{
			for(i=fri;RSB_LIKELY(i<=lri);++i)
			{
				nnz = mtxAp->bpntr[i+1]-mtxAp->bpntr[i];
				if(IA)rsb__util_coo_array_set(IA+(zoff+dnnz),nnz,i+roff+ioff);
				dnnz += nnz;
			}
			RSB_MEMCPY(JA+zoff,((rsb_half_idx_t*)mtxAp->bindx)+mtxAp->bpntr[fri],
					sizeof(rsb_half_idx_t)*dnnz);
			rsb__do_switch_array_to_fullword_coo((rsb_half_idx_t *)(JA+zoff),dnnz,mtxAp->coff+joff);
			//rsb__util_coo_array_add(JA+zoff,dnnz,mtxAp->coff+joff);
			doff = mtxAp->bpntr[fri];
			RSB_VA_MEMCPY(VA,MVA,zoff,doff,dnnz,mtxAp->el_size);
		}
#if 0
		else
		if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_COO)
		{
			/* FIXME: here, need fast seek to fri and lri */
			rsb_nnz_idx_t n,fi = mtxAp->bpntr[fri],li = mtxAp->bpntr[lri+1];//FIXME : relying on bpntr is EVIL !
		//	for(n=0;n<=mtxAp->nnz;++n)
			for(n=fi;n<li;++n)
			{
				rsb_coo_idx_t ij = mtxAp->bindx[n];
				rsb_coo_idx_t j = RSB_COO_HALFWORDS_VALUES_UNPACK_UJ(ij);
				rsb_coo_idx_t i = RSB_COO_HALFWORDS_VALUES_UNPACK_LI(ij);
//				if(i<fri || i>lri) continue;/* FIXME: slow ! */
				if(IA)IA[zoff+dnnz] = i+roff+ioff;
				if(JA)JA[zoff+dnnz] = j+coff+joff;
				dnnz += 1;
			}
			doff = mtxAp->bpntr[fri];
		}
#endif
		else
		{
			for(i=fri;RSB_LIKELY(i<=lri);++i)
			{
				nnz = mtxAp->bpntr[i+1]-mtxAp->bpntr[i];
				RSB_CSR2COO_MEMCPY_(VA,IA,JA,MVA,i+roff+ioff,mtxAp->bindx,zoff+dnnz,
					mtxAp->bpntr[i],nnz,mtxAp->el_size,coff+joff);
				dnnz += nnz;
			}
			doff = mtxAp->bpntr[fri];
			RSB_VA_MEMCPY(VA,MVA,zoff,doff,dnnz,mtxAp->el_size);
		}
		}
		*rnz += dnnz;

#if 0
			if(IA && JA && (mtxAp->nnz==3))for(i=0;i<dnnz;++i)
			{
				RSB_STDOUT("at %d %d\n",1+IA[i],1+JA[i]);
			}
#endif
	}
err:
        RSB_DO_ERR_RETURN(errval)
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_get_rows_dense(const struct rsb_mtx_t *mtxAp , void * row , rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t  *rnz, rsb_flags_t flags )
{
        /*!
	 * \ingroup gr_internals
         *
	 * FIXME : rename this function
	 * FIXME : THIS IS NOT ROWS_DENSE ! IT DOES SOMETHING ELSE !
	 *
         * \note This function is slow.
         * */
        rsb_coo_idx_t i,j;
	rsb_nnz_idx_t l;
	rsb_nnz_idx_t nnz = 0,gap = 0,discarded = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_time_t gt,ct;

	if(!mtxAp || !rnz)
	{errval = RSB_ERR_BADARGS;RSB_PERR_GOTO(err,RSB_ERRM_ES)}
					

	if(fr<0 || lr>mtxAp->nr)
	{errval = RSB_ERR_BADARGS;RSB_PERR_GOTO(err,RSB_ERRM_ES)}

	if(!row || fr>lr)
	{errval = RSB_ERR_BADARGS;RSB_PERR_GOTO(err,RSB_ERRM_ES)}

	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(fr));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(lr));

	/* input fortran indices */
	if( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
	{
		lr--;
		fr--;
	}
	
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(fr));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(lr));

	/* FIXME : slow */
        for(i=fr;RSB_LIKELY(i<=lr);++i)
        for(j=0;RSB_LIKELY(j<mtxAp->nc);++j)
	{
		l = mtxAp->nc*(i-fr)+j;
		IA[l] = i;
		JA[l] = j;
	}

//	for(i=0;i<mtxAp->nc;++i)
//		printf("%d %d %lg\n",IA[i],JA[i],((double*)row)[i]);

	gt = - rsb_time();
	RSB_BZERO(row,mtxAp->el_size*mtxAp->nc*(lr-fr+1));
        for(i=fr;RSB_LIKELY(i<=lr);++i)
        {
		errval = rsb__do_get_row_dense(mtxAp,((rsb_byte_t*)row)+(mtxAp->el_size*(i-fr))*mtxAp->nc , i);
		if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	gt += rsb_time();

//	printf("!\n");
//	for(i=0;i<mtxAp->nc*(lr-fr+1);++i)
//		printf("%d %d %lg\n",IA[i],JA[i],((double*)row)[i]);
	
	/* output fortran indices */
	if( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
        for(i=fr;RSB_LIKELY(i<=lr);++i)
        for(j=0;RSB_LIKELY(j<mtxAp->nc);++j)
	{
		l = mtxAp->nc*(i-fr)+j;
		IA[l]++;
		JA[l]++;
	}

	nnz = ((lr+1)-fr)*mtxAp->nc;
	/* FIXME : (SLOW!) (do we really need this ?) 
	 * FIXME : THIS IS NOT ANYMORE ROWS_DENSE !
	 * */
	ct = - rsb_time();
	rsb__util_compact_nonzeros(row,IA,JA,nnz,mtxAp->typecode,&gap,&discarded,RSB_FLAG_NOFLAGS);
	ct += rsb_time();
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(discarded));
	RSB_BZERO(((rsb_byte_t*)row)+(mtxAp->el_size*(nnz-discarded)),mtxAp->el_size*discarded);//NEW
	//nnz -= gap;
	nnz -= discarded;

//	printf("\ncompa %lg   getrow %lg\n",ct,gt);
/*	if(mtxAp->nr==1 && mtxAp->nc==1)
	{
		printf("\nMATRIX %d:%d %lg\n",fr,lr,*((double*)mtxAp->VA));
		printf("\nGETROWDENSEFIRST: %d %d %lg\n",IA[0],JA[0],((double*)row)[0]);	
	}*/

	if(discarded)
	{
		IA[nnz] = 0;	/* FUNDAMENTAL! FIXME ! */
		JA[nnz] = 0;
	}

	*rnz = nnz;	// we notify the callee

//	for(i=0;i<nnz;++i)
//		printf("%d %d %lg\n",IA[i],JA[i],((double*)row)[i]);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nnz));
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if RSB_OBSOLETE_QUARANTINE_UNUSED
static void rsb__do_extract_nnz_from_block(void * blockpointer, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t br, rsb_coo_idx_t bc, rsb_coo_idx_t baserow, rsb_coo_idx_t basecolumn, rsb_type_t typecode, size_t el_size, rsb_nnz_idx_t *nnzp )
{
	/**
		FIXME
	*/
	rsb_coo_idx_t r,c;
	rsb_nnz_idx_t nz = 0;

	for(r=0;r<br;++r)
	for(c=0;c<bc;++c)
	{
		rsb_byte_t*src=((rsb_byte_t*)blockpointer)+(el_size*(r*bc+c));
		/*
		 * SERVE UNA NUOVA MACRO CHE INDIVIDUI PER CIASCUN TIPO SE AREA DI MEMORIA E' NNZ
		 * */
		if( RSB_IS_ELEMENT_NONZERO(src,typecode) )
		{
			rsb_byte_t * dst = ((rsb_byte_t*)(VA)) + el_size * (nz);
			RSB_NUMERICAL_TYPE_SET_ELEMENT(dst,src,typecode) /* FIXME : this is SLOW */
			(IA)[nz] = baserow   +r;
			(JA)[nz] = basecolumn+c;
			++nz;
		//	if(nz>mtxAp->nnz)goto fatal_err;/* PLACE HERE ERROR CHECKING CODE */
		}
	}
	*nnzp += nz;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__do_get_coo_noalloc(const struct rsb_mtx_t *mtxAp, rsb_byte_t * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t * nnzp, rsb_flags_t flags)
{
	/*! 
	 * \ingroup gr_internals
	 *
	 *  Returns the matrix converted in a coordinate storage format.
	 *
	 * \param VA  the values array pointer, sized at least for mtxAp->nnz elements of matrix type
	 * \param IA  an integer array pointer for row    coordinates
	 * \param JA  an integer array pointer for column coordinates
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * The entire matrix will be returned in COO format, in the specified VA,IA,JA arrays
	 * No more than mtxAp->nnz elements will be written to in the VA, IA and JA arrays
	 * 
	 * FIXME : should add an offset argument, for recursive matrices.
	 * */
	register rsb_nnz_idx_t baserow,basecolumn;
	register rsb_blk_idx_t rows,columns;
	register rsb_blk_idx_t blockrow,blockcolumn;
	register rsb_byte_t *bp;
	size_t el_size = 0;
	rsb_nnz_idx_t nz = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( !IA || !JA || !VA )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if( !mtxAp )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	el_size = mtxAp->el_size;

	{
		rsb_nnz_idx_t dnnz = 0;
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_rows_sparse(RSB_TRANSPOSITION_N,NULL,mtxAp,VA,IA,JA,0,mtxAp->nr-1,&dnnz,RSB_FLAG_NOFLAGS));
		if(nnzp)
			*nnzp = dnnz;
		if( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
			rsb__util_coo_array_to_fortran_indices_parallel(IA,dnnz),
			rsb__util_coo_array_to_fortran_indices_parallel(JA,dnnz);
		goto ret;
	}

#if 0
	// FIXME: THE FOLLOWING IS OLD AND BROKEN
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		struct rsb_mtx_t * submatrix;
		rsb_submatrix_idx_t i,j;
		rsb_nnz_idx_t nzoff = 0;

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			rsb_nnz_idx_t snnz;
			rsb_coo_idx_t moff = submatrix->roff-mtxAp->roff;
			rsb_coo_idx_t koff = submatrix->coff-mtxAp->coff;

			snnz = submatrix->nnz;

			RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_COUNT(snnz));
			RSB_DEBUG_ASSERT(RSB_IS_VALID_NNZ_INDEX(nzoff));

			errval = rsb__do_get_coo_noalloc(submatrix,VA+el_size*nzoff,IA+nzoff,JA+nzoff,nnzp,flags);
			if(RSB_SOME_ERROR(errval))
			{
				RSB_PERR_GOTO(erri,RSB_ERRM_ES)
			}
			rsb__util_coo_arrays_add(IA+nzoff,JA +nzoff, moff*i, koff*j, snnz);

			nzoff += snnz;
		}
erri:
		goto ret;
	}
#endif

	RSB_DEBUG_ASSERT(mtxAp->bpntr);
	RSB_DEBUG_ASSERT(mtxAp->indptr);
	RSB_DEBUG_ASSERT(mtxAp->bindx);

	RSB_BZERO(VA,el_size     * mtxAp->nnz);
	RSB_BZERO(IA,sizeof(rsb_coo_idx_t) * mtxAp->nnz);
	RSB_BZERO(JA,sizeof(rsb_coo_idx_t) * mtxAp->nnz);


#if !RSB_EXPERIMENTAL_USE_PURE_BCSS
	if(rsb__is_bcss_matrix(mtxAp))
	{
		rsb_nnz_idx_t bi, nz = 0, bri, bci;

		if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
		{
			RSB_BCSC_MATRIX_FOREACH_BLOCK(mtxAp,bp,bri,bci,bi,baserow,basecolumn)
				rsb__do_extract_nnz_from_block(bp, ((rsb_byte_t*)(VA))+el_size*nz, IA+nz, JA+nz, mtxAp->br, mtxAp->bc, baserow, basecolumn, mtxAp->typecode, el_size, &nz);
		}
		else
		{
			RSB_BCSR_MATRIX_FOREACH_BLOCK(mtxAp,bp,bri,bci,bi,baserow,basecolumn)
				rsb__do_extract_nnz_from_block(bp, ((rsb_byte_t*)(VA))+el_size*nz, IA+nz, JA+nz, mtxAp->br, mtxAp->bc, baserow, basecolumn, mtxAp->typecode, el_size, &nz);
		}
		return RSB_ERR_NO_ERROR;
	}
#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
	RSB_DEBUG_ASSERT(mtxAp->rpntr);
	RSB_DEBUG_ASSERT(mtxAp->cpntr);

	{
	RSB_GET_FIRST_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
	while(!RSB_GOT_LAST_BLOCK_POINTER(mtxAp))
	{
		/*
			FIXME : should better use rsb__do_extract_nnz_from_block !
		*/
		rsb_coo_idx_t r,c;
		for(r=0;r<rows;++r)
		for(c=0;c<columns;++c)
		{
			rsb_byte_t*src = (rsb_byte_t*)((double*)(bp+RSB_GET_INTRA_BLOCK_OFFSET(baserow+r,basecolumn+c,blockrow,blockcolumn,mtxAp)));
			RSB_DEBUG_ASSERT(src>=(rsb_byte_t*)bp);
			RSB_DEBUG_ASSERT(src>=(rsb_byte_t*)mtxAp->VA);
//			RSB_DEBUG_ASSERT(src<(rsb_byte_t*)mtxAp->VA+el_size*mtxAp->element_count);// intrablock struct breaks this because of extra space
			/*
			 * SERVE UNA NUOVA MACRO CHE INDIVIDUI PER CIASCUN TIPO SE AREA DI MEMORIA E' NNZ
			 * */
			if(  RSB_IS_ELEMENT_NONZERO(src,mtxAp->typecode) )
			{
				rsb_byte_t * dst = (VA) + el_size * (nz);
				RSB_NUMERICAL_TYPE_SET_ELEMENT(dst,src,mtxAp->typecode) /* FIXME : this is SLOW */
				(IA)[nz] = baserow   +r;
				(JA)[nz] = basecolumn+c;
				++nz;
				if(nz>mtxAp->nnz)goto fatal_err;/* PLACE HERE ERROR CHECKING CODE */
			}
		}
		RSB_GET_NEXT_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
	}
	}
	/* FIXME : this should happen only in verbose mode */
	if(nz<mtxAp->nnz)
		RSB_INFO("warning : counted less nonzeros (%zd) than it should (%zd) (may be zeros..)!\n",(size_t)nz,(size_t)mtxAp->nnz);

	if(nz>mtxAp->nnz)
	{
		goto fatal_err;
	}
	if(nnzp)
		*nnzp = nz;
	else
	{
		/* FIXME : WRITE ME */	
	}

	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_BADARGS;
fatal_err:
	RSB_ERROR("fatal: counted more nonzeros (%zd) than it should (%zd)!\n",(size_t)nz,(size_t)mtxAp->nnz);
	errval = RSB_ERR_BADARGS;
ret:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_rows_sums_inner(const struct rsb_mtx_t *mtxAp , rsb_byte_t * row_sums, rsb_bool_t do_testing, rsb_trans_t transA)
{
	/*
		FIXME : document
	*/
	//size_t el_size = mtxAp->el_size;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( mtxAp && rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;

		if(RSB_SOME_ERROR(errval))
		{RSB_PERR_GOTO(err,RSB_ERRM_ES)}

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			/*
			rsb_coo_idx_t off = 0;
			if(RSB_DOES_NOT_TRANSPOSE(transA))
				off = submatrix->roff-mtxAp->roff;
			else
				off = submatrix->coff-mtxAp->coff;
			*/
			errval = rsb__do_rows_sums_inner(submatrix,row_sums/*+el_size*(off)*/,do_testing,transA);
			if(RSB_SOME_ERROR(errval))
			{RSB_PERR_GOTO(err,RSB_ERRM_ES)}
		}
	}
	else
	{
#ifndef RSB_HAVE_OPTYPE_INFTY_NORM
	return RSB_ERR_UNSUPPORTED_OPERATION;
#else /* RSB_HAVE_OPTYPE_INFTY_NORM */
#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
		if(do_testing!=0 || (mtxAp->flags & RSB_FLAG_SHOULD_DEBUG))
		{
			// uhm.. FIXME : do we really need this ?
			if(RSB_SOME_ERROR(errval = rsb__rowssums_testing(mtxAp,transA,row_sums)))
			{RSB_PERR_GOTO(err,RSB_ERRM_ES)}
		}
		else
#endif /* RSB_WANT_KERNELS_DEBUG */
		{
			if(RSB_SOME_ERROR(errval = rsb__do_rowssums(mtxAp,transA,row_sums)))
			{RSB_PERR_GOTO(err,RSB_ERRM_ES)}
		}
#endif /* RSB_HAVE_OPTYPE_INFTY_NORM */
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb__do_infinity_norm_inner( const struct rsb_mtx_t *mtxAp, rsb_byte_t * row_asums, rsb_bool_t do_testing, rsb_trans_t transA)
{
	/*
	 * On each row/column i, accumulate row/column absolutes sum in row_asums[i].
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( mtxAp && rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;

		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			errval = rsb__do_infinity_norm_inner(submatrix,row_asums,do_testing,transA);
			if(RSB_SOME_ERROR(errval))
			{
				RSB_PERR_GOTO(err,RSB_ERRM_ES)
			}
		}
	}
	else
	{
#ifndef RSB_HAVE_OPTYPE_INFTY_NORM
		return RSB_ERR_UNSUPPORTED_OPERATION;
#else /* RSB_HAVE_OPTYPE_INFTY_NORM */
#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
		if(do_testing!=0 || (mtxAp->flags & RSB_FLAG_SHOULD_DEBUG))
		{
			// uhm.. FIXME : do we really need this ?
			if(RSB_SOME_ERROR(errval = rsb__infty_norm_testing(mtxAp,transA,row_asums)))
			{RSB_PERR_GOTO(err,RSB_ERRM_ES)}
		}
		else
#endif /* RSB_WANT_KERNELS_DEBUG */
		{
			if(RSB_SOME_ERROR(errval = rsb__do_infty_norm(mtxAp,transA,row_asums)))
			{RSB_PERR_GOTO(err,RSB_ERRM_ES)}
		}
#endif /* RSB_HAVE_OPTYPE_INFTY_NORM */
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_absolute_rows_sums( const struct rsb_mtx_t * mtxAp, void * d)
{
	return rsb__do_infinity_norm_inner(mtxAp,d,RSB_BOOL_FALSE,RSB_TRANSPOSITION_N);
}

rsb_err_t rsb__do_absolute_columns_sums( const struct rsb_mtx_t * mtxAp, void * d)
{
	return rsb__do_infinity_norm_inner(mtxAp,d,RSB_BOOL_FALSE,RSB_TRANSPOSITION_T);
}

rsb_err_t rsb__do_infinity_norm(const struct rsb_mtx_t *mtxAp, void * infinity_norm, const rsb_bool_t do_testing, rsb_trans_t transA)
{
	/*!
	 * \ingroup gr_internals
	 * Returns the maximum of sums of absolute values of elements in a row, over all rows.
	 *
	 * \input mtxAp is a pointer to a valid matrix.
	 * \input infinity_norm is a pointer to the location where the norm will be written
	 * also known as the maximum absolute row sum norm : ||A|| = max_i sum_j^n |a_ij|
	 *
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#ifndef RSB_HAVE_OPTYPE_INFTY_NORM
	return RSB_ERR_UNSUPPORTED_OPERATION;
#else /* RSB_HAVE_OPTYPE_INFTY_NORM */
	void * row_sums = NULL,*mrp = NULL;
	rsb_blk_idx_t maximal_row = 0;
	const rsb_coo_idx_t tm = mtxAp ? RSB_MTX_TRANSPOSED_ROWS(mtxAp,transA) : 0;

        RSB_DEBUG_ASSERT(mtxAp);

	if(!infinity_norm)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(ret,RSB_ERRM_ES);
	}

	row_sums = rsb__calloc(mtxAp->el_size*(tm+RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE));
	if(!row_sums)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(ret,RSB_ERRM_ES)
	}

	errval = rsb__do_infinity_norm_inner(mtxAp,row_sums,do_testing,transA);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(ret,RSB_ERRM_ES)
	}

	/*
	 * TODO : should use BLAS own icamax or similar to determine the maximal row sum
	 * */
	RSB_VECTOR_FIND_MAXIMAL_ELEMENT(maximal_row,row_sums,tm,mtxAp->typecode);

	if(maximal_row<0)
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(ret,RSB_ERRM_ES)
	}

	/* 
	 * that row sum is our infinity norm now
	 * */
	mrp = (((rsb_byte_t*)row_sums)+(mtxAp->el_size*maximal_row));
	if(rsb__get_diagonal_type_flag(mtxAp)==RSB_DIAGONAL_I)
		rsb__util_increase_by_one(mrp,0,mtxAp->typecode);

	RSB_NUMERICAL_TYPE_SET_ELEMENT_REAL(infinity_norm,mrp,mtxAp->typecode)
ret:
	RSB_CONDITIONAL_FREE(row_sums);
	RSB_DO_ERR_RETURN(errval)
#endif /* RSB_HAVE_OPTYPE_INFTY_NORM */
}

double rsb__do_get_matrix_fillin(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * This function returns the fillin of a given matrix.
	 * The fillin is defined as the ratio of the allocated elements 
	 * count and the original matrix nonzeros.
	 *
         * This function is expected to be quite fast (so, it's won't be a number
         * crunching routine, no matter how fancy our data structures are).
	 *
	 * \return NULL on error, othwerwise the fillin
	 * */
	if(!mtxAp)
		return 0.0;
	return ((double)mtxAp->element_count)/((double)mtxAp->nnz) ;/* numbers / nonzeros */
}

rsb_nnz_idx_t rsb__do_get_matrix_nnz(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * This function returns the number of nonzeros of a given sparse matrix.
	 *
	 * \return NULL on error, othwerwise the nonzeros
	 * */
	return (mtxAp->nnz) ;/* nonzeros */
}

rsb_long_t rsb__submatrices(const struct rsb_mtx_t * mtxAp)
{
	/*!
	 * 	\ingroup gr_internals
	 * 	Counts submatrices: either leaves or nodes.
	 */
	rsb_long_t sm = 0;
	rsb_submatrix_idx_t i,j;
	const struct rsb_mtx_t * submatrix;

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	{
		if(submatrix)
		{
			sm += rsb__submatrices(submatrix);
		}
	}
	return sm+1;
}

void rsb__get_blocking_size(const struct rsb_mtx_t * mtxAp, rsb_blk_idx_t *brp, rsb_blk_idx_t *bcp)
{
	/*!
	 *
	 * \ingroup gr_internals
	 *
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * \attention : before it used rsb_coo_idx_t !
	 * */
	RSB_ASSERT(mtxAp);
	RSB_ASSERT(!( RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_NON_ROOT_MATRIX) && !(rsb__have_fixed_blocks_matrix_flags(mtxAp->flags) /*|| rsb__is_bcss_matrix(mtxAp)*/ )));
	
#if (RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS)
	if( (mtxAp->rpntr==NULL) && ( mtxAp->cpntr==NULL) )
	{

#if RSB_EXPERIMENTAL_USE_PURE_BCSS
		RSB_ASSERT(mtxAp->br>0);
		RSB_ASSERT(mtxAp->bc>0);
		*brp= mtxAp->br;
		*bcp= mtxAp->bc;
#else /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
		*brp = 1;
		*bcp = 1;
#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
	}
	else
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
	{
		*brp = mtxAp->rpntr[1]-mtxAp->rpntr[0];	/* we assume that block_count >= 1 */
		*bcp = mtxAp->cpntr[1]-mtxAp->cpntr[0];	/* we assume that block_count >= 1 */
	}
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(*brp));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(*bcp));
	return;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__get_blocking_size_as_row_major(const struct rsb_mtx_t * mtxAp, rsb_blk_idx_t *bMp, rsb_blk_idx_t *bmp)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * */
	rsb_err_t errval = RSB_ERR_BADARGS;

	if(!mtxAp)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	
	errval = RSB_ERR_NO_ERROR;
	if( mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
		rsb__get_blocking_size(mtxAp, bmp, bMp);
	else
		rsb__get_blocking_size(mtxAp, bMp, bmp);
err:
	return errval;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
size_t rsb__do_get_max_blocksize(const struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 *	FIXME : this is for VBR and recursive matrices.
	 *      it is UNFINISHED
	 */
	rsb_blk_idx_t maxMd = 0,maxmd = 0;
	rsb_blk_idx_t i;
	size_t sz;

//	if(rsb__is_recursive_matrix(mtxAp->flags))
//		RSB_WARN("rsb__do_get_max_blocksize unfinished for recursive formats!\n");

	RSB_DEBUG_ASSERT(mtxAp);

	if(rsb__have_fixed_blocks_matrix_flags(mtxAp->flags))
		rsb__get_blocking_size_as_row_major(mtxAp, &maxMd, &maxmd);
	else
	{
		RSB_DEBUG_ASSERT(mtxAp->Mdim);
		RSB_DEBUG_ASSERT(mtxAp->mdim);
		RSB_DEBUG_ASSERT(mtxAp->Mpntr);
		RSB_DEBUG_ASSERT(mtxAp->mpntr);

		for(i=1;RSB_LIKELY(i<=mtxAp->Mdim);++i)
			if(mtxAp->Mpntr[i]-mtxAp->Mpntr[i-1]>maxMd)
				maxMd = mtxAp->Mpntr[i]-mtxAp->Mpntr[i-1];
		for(i=1;RSB_LIKELY(i<=mtxAp->mdim);++i)
			if(mtxAp->mpntr[i]-mtxAp->mpntr[i-1]>maxmd)
				maxmd = mtxAp->mpntr[i]-mtxAp->mpntr[i-1];
	}
	sz = mtxAp->el_size;
	sz *= maxmd;
	sz *= maxMd;
	return sz;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__get_physical_blocking_size(const struct rsb_mtx_t * mtxAp, rsb_blk_idx_t *brp, rsb_blk_idx_t *bcp)
{
	const rsb_err_t errval = RSB_ERR_NO_ERROR;
	if( mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
		rsb__get_blocking_size(mtxAp, bcp, brp);
	else
		rsb__get_blocking_size(mtxAp, brp, bcp);
	return errval;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_submatrix_idx_t rsb__get_recursive_matrix_depth(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_submatrix_idx_t msmd = 0,smd = 0;

	if(!mtxAp)
		return 0;

	if(!rsb__is_recursive_matrix(mtxAp->flags))
		goto norec;

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
	{
		smd = rsb__get_recursive_matrix_depth(submatrix);
		msmd = RSB_MAX(smd,msmd);
	}

	return msmd+1;
norec:
	return 0;
}

rsb_flags_t rsb__get_hermitian_flag(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 */
	return (mtxAp->flags & RSB_FLAG_HERMITIAN);
}

rsb_flags_t rsb__get_symmetry_flag(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 */
	return (mtxAp->flags & RSB_FLAG_SYMMETRIC);
}

rsb_flags_t rsb__get_diagonal_type_flag(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 */
	if(!mtxAp)
		return RSB_FLAG_NOFLAGS;

#ifdef RSB_DIAGONAL_I
	return (mtxAp->flags & RSB_FLAG_UNIT_DIAG_IMPLICIT) ? RSB_DIAGONAL_I : RSB_DIAGONAL_E ;
#else /* RSB_DIAGONAL_I */
	return RSB_DIAGONAL_E ;
#endif /* RSB_DIAGONAL_I */
}

rsb_flags_t rsb__get_symmetry_type_flag(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * Used to dispatch to a truly symmetric kernel or not.
	 * So that e.g.
         * symmetric-diagonal matrices shall be OK for rsb_spsv/rsb_spsm.
	 */
	rsb_flags_t flags;

	RSB_DEBUG_ASSERT(mtxAp);

	flags = (rsb__get_hermitian_flag(mtxAp) | rsb__get_symmetry_flag(mtxAp));

	if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_DIAGONAL))
		RSB_DO_FLAG_DEL(flags,RSB_FLAG_ANY_SYMMETRY);

	return flags;
}

RSB_INLINE rsb_coo_idx_t rsb__do_get_rows_of(const struct rsb_mtx_t *mtxAp, rsb_trans_t transA)
{
	/*!
	 * \ingroup gr_internals
	 */
	RSB_DEBUG_ASSERT(mtxAp);
	return RSB_MTX_TRANSPOSED_ROWS(mtxAp,transA);
}

RSB_INLINE rsb_coo_idx_t rsb__do_get_columns_of(const struct rsb_mtx_t *mtxAp, rsb_trans_t transA)
{
	/*!
	 * \ingroup gr_internals
	 */
	RSB_DEBUG_ASSERT(mtxAp);
	return RSB_MTX_TRANSPOSED_COLS(mtxAp,transA);
}

static rsb_err_t rsb__do_get_elements_for_all_columns_inner(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * CP)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : NEW, UNFINISHED STUB, UNTESTED
	 * TODO : to parallelize this
	 * TODO : to declare this function in some header
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(CP);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_elements_for_all_columns_inner(submatrix,CP));
	}
	else
	{
		if(!rsb__is_css_matrix(mtxAp))
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err,RSB_ERRM_ES)
		}
		if(rsb__is_coo_matrix(mtxAp))
		{
			rsb_nnz_idx_t n;
			//const rsb_coo_idx_t roff = mtxAp->roff;
			const rsb_coo_idx_t coff = mtxAp->coff;
			if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				RSB_DECLARE_CONST_HALFCOO_JARRAY_FROM_MATRIX(mJA,mtxAp)
				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					rsb_half_idx_t /*i = mIA[n],*/j = mJA[n];
					(CP)[coff+j]++;
				}
			}
			else
			{
				RSB_DECLARE_CONST_FULLCOO_JARRAY_FROM_MATRIX(mJA,mtxAp)
				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					rsb_coo_idx_t /*i = mIA[n],*/j = mJA[n];
					(CP)[coff+j]++;
				}
			}
		}
		else
#if RSB_WANT_BCSC_LEAVES  
		if(rsb__is_bcsc_matrix(mtxAp))
		{
			rsb__util_nnz_array_add_array(CP+mtxAp->coff,mtxAp->bpntr,mtxAp->nc);
		}
		else
#endif /* RSB_WANT_BCSC_LEAVES */
		if(rsb__is_bcsr_matrix(mtxAp))
		{
			rsb_nnz_idx_t n;
			rsb_coo_idx_t *bindx = mtxAp->bindx;
			//const rsb_coo_idx_t roff = mtxAp->roff;
			const rsb_coo_idx_t coff = mtxAp->coff;

#if 0
			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_COO)
			{
				for(n=0;n<mtxAp->nnz;++n)
				{
					rsb_coo_idx_t ij = bindx[n],j = RSB_COO_HALFWORDS_VALUES_UNPACK_UJ(ij);
					(CP+1)[coff+j]++;
				}
			}
#endif
			if(mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES_CSR)
			{
				const rsb_half_idx_t *hbindx = (const rsb_half_idx_t *)bindx;
				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					rsb_coo_idx_t j = hbindx[n];
					(CP)[coff+j]++;
				}
			}
			else
			{
				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					rsb_coo_idx_t j = bindx[n];
					(CP)[coff+j]++;
				}
			}
		}
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb__do_get_elements_for_each_column(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * CP)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : NEW, UNFINISHED STUB
	 * TODO : to parallelize this
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//rsb_nnz_idx_t n;
	if(!mtxAp)
       	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(CP);
	rsb__util_nnz_array_set(CP,mtxAp->nc,0);
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_elements_for_all_columns_inner(mtxAp,CP));
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb__do_get_elements_for_all_columns(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * CP)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : NEW, UNFINISHED STUB
	 * TODO : to parallelize this
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t n;
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(CP);
	rsb__util_nnz_array_set(CP,mtxAp->nc+1,0);
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_elements_for_each_column(mtxAp,CP+1));
	for(n=0;n<mtxAp->nc;++n)CP[n+1] += CP[n];
//err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb__do_get_elements_for_all_rows_inner(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * RP)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : NEW, UNFINISHED STUB
	 * TODO : to parallelize this
	 * TODO : to declare this function in some header
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(RP);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		const struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_elements_for_all_rows_inner(submatrix,RP));
	}
	else
	{
		rsb_nnz_idx_t n;
		rsb_nnz_idx_t i;
		if(rsb__is_coo_matrix(mtxAp))
		{
			const rsb_coo_idx_t roff = mtxAp->roff;
			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			{
				RSB_DECLARE_CONST_HALFCOO_IARRAY_FROM_MATRIX(mIA,mtxAp)
				for(n=0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					rsb_half_idx_t i = mIA[n];//,j = mJA[n];
					(RP)[roff+i]++;
				}
			}
			else
			{
				RSB_DECLARE_CONST_FULLCOO_IARRAY_FROM_MATRIX(mIA,mtxAp)
				for(n = 0;RSB_LIKELY(n<mtxAp->nnz);++n)
				{
					rsb_coo_idx_t i = mIA[n];//,j = mJA[n];
					(RP)[roff+i]++;
				}
			}
		}
		else
		{
//			if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
//				rsb__util_nnz_array_add_array(RP+mtxAp->roff,mtxAp->bpntr+1,mtxAp->nr-1);
//			else
//				rsb__util_nnz_array_add_array(RP+mtxAp->roff,mtxAp->bpntr+1,mtxAp->nr-1);
			for(i=0;i<mtxAp->nr;++i)
			{
				(RP)[mtxAp->roff+i] += mtxAp->bpntr[i+1]-mtxAp->bpntr[i];
			}
		}
	}
//err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb__do_get_elements_for_all_rows(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * RP)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : NEW, UNFINISHED STUB
	 * TODO : to parallelize this
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//rsb_nnz_idx_t nzi;
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(RP);

	rsb__util_nnz_array_set(RP,mtxAp->nr+1,0);
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_elements_for_all_rows_inner(mtxAp,RP+1));
	//for(nzi=0;RSB_LIKELY(nzi<mtxAp->nr);++nzi) RP[nzi+1] += RP[nzi];
	rsb_do_prefix_sum_nnz_idx_t(RP,mtxAp->nr+1);
//err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__dodo_get_csr(const struct rsb_mtx_t *mtxAp, rsb_byte_t ** VA, rsb_nnz_idx_t ** RP, rsb_coo_idx_t ** JA)
{
	/*!
	 * \ingroup gr_internals
	 */
	rsb_nnz_idx_t rnz = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(RSB_UNLIKELY(!RP || !*RP || !mtxAp || !VA || !*VA || !JA || !*JA))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

//	RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_rows_sparse_rec(mtxAp,*VA,0,mtxAp->nr-1,NULL,*JA,&rnz,0,0));
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_rows_sparse(RSB_TRANSPOSITION_N,NULL,mtxAp,*VA,NULL,*JA,0,mtxAp->nr-1,&rnz,RSB_FLAG_NOFLAGS));

	RSB_DO_ERROR_CUMULATE(errval,rsb__do_get_elements_for_all_rows(mtxAp,*RP));
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_csc(const struct rsb_mtx_t *mtxAp, rsb_byte_t ** VA, rsb_nnz_idx_t ** CPp,rsb_coo_idx_t ** IA)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : NEW, UNFINISHED STUB
	 * TODO : to parallelize this
	 * NOTE : sort in a SPSV-like order the matrices; compute the columns pointer vector; extract ;
	 * TODO : to parallelize this
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_translated_matrix_t * all_leaf_matrices = NULL;	/** NEW, EXPERIMENTAL */
	rsb_submatrix_idx_t all_leaf_matrices_n = 0;
	rsb_submatrix_idx_t * deps = NULL;	/** NEW, EXPERIMENTAL */
	rsb_nnz_idx_t n;
	//return 0;
	rsb_nnz_idx_t * CP = *CPp;

	if(!mtxAp||!VA||!CPp||!IA)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	errval = rsb__do_get_submatrices_block_for_get_csr(mtxAp,&all_leaf_matrices,&all_leaf_matrices_n/*,RSB_TRANSPOSITION_N*/);/* ! */
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
//	deps = rsb__malloc(sizeof(rsb_submatrix_idx_t)*all_leaf_matrices_n);
//	if(RSB_SOME_ERROR(errval) || !all_leaf_matrices || !deps)
//	{errval = RSB_ERR_ENOMEM;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
//
	if(0)
	for(n=0;n<all_leaf_matrices_n;++n)
	{
		all_leaf_matrices[n].roff = mtxAp->nc-(all_leaf_matrices[n].coff+all_leaf_matrices[n].mtxlp->nc*1);
		all_leaf_matrices[n].coff = mtxAp->nr-(all_leaf_matrices[n].roff+all_leaf_matrices[n].mtxlp->nr*1);
		all_leaf_matrices[n].nr = all_leaf_matrices[n].mtxlp->nc;
		all_leaf_matrices[n].nc = all_leaf_matrices[n].mtxlp->nr;
	}

	errval = rsb__do_get_elements_for_all_columns(mtxAp,CP);/*  */
	//errval = rsb__do_get_elements_for_each_column(mtxAp,CP);
	/* ... */

	for(n=0;RSB_LIKELY(n<all_leaf_matrices_n);++n)
	{
		/* extract elements ... */
		struct rsb_mtx_t *submatrix = all_leaf_matrices[n].mtxlp;
		errval |= rsb__do_get_columns_sparse(submatrix,*VA,0,mtxAp->nc-1,*IA,NULL,CP,0,0);
//		rsb__do_get_columns_sparse(submatrix,*VA,0,mtxAp->nc-1,*IA,NULL,CP,submatrix->roff,submatrix->coff);
	}
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES)

	for(n=mtxAp->nc;RSB_LIKELY(n>0);--n) CP[n] = CP[n-1];
	CP[0] = 0;

	// FIXME : should check for correctness, now
	if(RSB_UNLIKELY(CP[mtxAp->nc]!=mtxAp->nnz))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,"%d != %d\n",CP[mtxAp->nc],mtxAp->nnz)
	}

err:
	RSB_CONDITIONAL_FREE(deps);
	RSB_CONDITIONAL_FREE(all_leaf_matrices);
	return RSB_ERR_NO_ERROR;
}

#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_get_block_sparse_pattern(const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t fr, rsb_coo_idx_t lr, rsb_coo_idx_t fc, rsb_coo_idx_t lc, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, const rsb_coo_idx_t * IREN, const rsb_coo_idx_t * JREN, rsb_nnz_idx_t *rnz, rsb_flags_t flags )
{
	/*! 
	 * \ingroup gr_internals
	 * */
	return rsb__do_get_block_sparse(mtxAp,NULL,IA,JA,fr,lr,fc,lc,IREN,JREN,rnz,flags);
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

static rsb_err_t rsb__do_get_block_sparse_leaf(const struct rsb_mtx_t * mtxAp, void* OVA, rsb_coo_idx_t fr, rsb_coo_idx_t lr, rsb_coo_idx_t fc, rsb_coo_idx_t lc, rsb_coo_idx_t * OIA, rsb_coo_idx_t * OJA, const rsb_coo_idx_t * IREN, const rsb_coo_idx_t * JREN, rsb_nnz_idx_t *rnz, rsb_coo_idx_t ioff, rsb_coo_idx_t joff, rsb_flags_t flags)
{
	/*! 
	 * \ingroup gr_internals
	 *
	 * FIXME: IREN/JREN untested
	 * */
	rsb_nnz_idx_t nnz = 0;
	rsb_coo_idx_t i = 0,roff = mtxAp->roff,coff = mtxAp->coff;
	RSB_DEFAULT_TYPE *VA = mtxAp->VA;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_nnz_idx_t zoff = *rnz;// FIXME: rnz is MANDATORY 

	fr -= mtxAp->roff; lr -= mtxAp->roff; fc -= mtxAp->coff; lc -= mtxAp->coff; 

	if(!mtxAp->bindx || !mtxAp->bpntr)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	if(rsb__is_coo_matrix(mtxAp))
	{
		rsb_nnz_idx_t nnz0 = 0,nnz1 = mtxAp->nnz;
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))
		{
			RSB_DECLARE_CONST_HALFCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			nnz0 += rsb__nnz_split_hcoo_bsearch(IA+nnz0,fr,nnz1-nnz0);
			nnz1 = nnz0+rsb__nnz_split_hcoo_bsearch(IA+nnz0,lr+1,nnz1-nnz0);
			if(nnz1-nnz0)
			{
				rsb_nnz_idx_t fnz = 0, lnz = 0, rnz = 0;
				for(i=fr;RSB_LIKELY(i<=lr);++i)
				{
					fnz = lnz;
					lnz = nnz1;
					rnz = lnz-fnz;
					fnz = fnz+rsb__nnz_split_hcoo_bsearch(IA+fnz,i,rnz);
					rnz = lnz-fnz;
					lnz = fnz+rsb__nnz_split_hcoo_bsearch(IA+fnz,i+1,rnz);
					rnz = lnz-fnz;
					fnz = fnz+rsb__nnz_split_hcoo_bsearch(JA+fnz,fc,rnz);
					rnz = lnz-fnz;
					lnz = fnz+rsb__nnz_split_hcoo_bsearch(JA+fnz,lc+1,rnz);

					if(OVA)
						RSB_VA_MEMCPY(OVA,VA,zoff+nnz,fnz,rnz,mtxAp->el_size);
					if(OJA)
						RSB_IA_MEMCPY_H(OJA,JA,zoff+nnz,fnz,rnz,coff+joff);
					if(OIA)
						RSB_IA_MEMCPY_H(OIA,IA,zoff+nnz,fnz,rnz,roff+ioff);
					nnz += lnz-fnz;				}
			}
		}
		else
		{
			RSB_DECLARE_CONST_FULLCOO_ARRAYS_FROM_MATRIX(IA,JA,mtxAp)
			nnz0 += rsb__nnz_split_coo_bsearch(IA+nnz0,fr,nnz1-nnz0);
			nnz1 = nnz0+rsb__nnz_split_coo_bsearch(IA+nnz0,lr+1,nnz1-nnz0);
			if(nnz1-nnz0)
			{
				rsb_nnz_idx_t fnz = 0, lnz = 0, rnz = 0;
				for(i=fr;RSB_LIKELY(i<=lr);++i)
				{
					fnz = lnz;
					lnz = nnz1;
					rnz = lnz-fnz;
					fnz = fnz+rsb__nnz_split_coo_bsearch(IA+fnz,i,rnz);
					rnz = lnz-fnz;
					lnz = fnz+rsb__nnz_split_coo_bsearch(IA+fnz,i+1,rnz);
					rnz = lnz-fnz;
					fnz = fnz+rsb__nnz_split_coo_bsearch(JA+fnz,fc,rnz);
					rnz = lnz-fnz;
					lnz = fnz+rsb__nnz_split_coo_bsearch(JA+fnz,lc+1,rnz);

					if(OVA)
						RSB_VA_MEMCPY(OVA,VA,zoff+nnz,fnz,rnz,mtxAp->el_size);

					if(OJA)
						RSB_IA_MEMCPY(OJA,JA,zoff+nnz,fnz,rnz,coff+joff);
					if(OIA)
						RSB_IA_MEMCPY(OIA,IA,zoff+nnz,fnz,rnz,roff+ioff);
					nnz += lnz-fnz;
				}
			}
		}
	}
	else
	if(rsb__is_csr_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
		{
			RSB_DECLARE_CONST_HALFCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			for(i=fr;RSB_LIKELY(i<=lr);++i)
			{
				rsb_nnz_idx_t lnz = PA[i+1],fnz = PA[i],rnz = lnz-fnz;
				if(rnz)
				{
					lnz = fnz+rsb__nnz_split_hcoo_bsearch(JA+fnz,lc+1,rnz);
					fnz = fnz+rsb__nnz_split_hcoo_bsearch(JA+fnz,fc,rnz);
					rnz = lnz-fnz;
					if(OVA)
						RSB_VA_MEMCPY(OVA,VA,zoff+nnz,fnz,rnz,mtxAp->el_size);
					if(OJA)
						RSB_IA_MEMCPY_H(OJA,JA,zoff+nnz,fnz,rnz,coff+joff);
					if(OIA)
						rsb__util_coo_array_set(OIA+(zoff+nnz),rnz,i+roff+ioff);
					nnz += rnz;
				}
			}
		}
		else
		{
			RSB_DECLARE_CONST_FULLCSR_ARRAYS_FROM_MATRIX(PA,JA,mtxAp)
			for(i=fr;RSB_LIKELY(i<=lr);++i)
			{
				rsb_nnz_idx_t lnz = PA[i+1],fnz = PA[i],rnz = lnz-fnz;
				if(rnz)
				{
					lnz = fnz+rsb__nnz_split_coo_bsearch(JA+fnz,lc+1,rnz);
					fnz = fnz+rsb__nnz_split_coo_bsearch(JA+fnz,fc,rnz);
					rnz = lnz-fnz;
					if(OVA)
						RSB_VA_MEMCPY(OVA,VA,zoff+nnz,fnz,rnz,mtxAp->el_size);
					if(OJA)
						RSB_IA_MEMCPY(OJA,JA,zoff+nnz,fnz,rnz,coff+joff);
					if(OIA)
						rsb__util_coo_array_set(OIA+(zoff+nnz),rnz,i+roff+ioff);
					nnz += rnz;
				}
			}
		}
	}
	else
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}


	if(IREN)
		rsb__util_coo_array_renumber(OIA+zoff,IREN,nnz,flags,flags,flags);
	if(JREN)
		rsb__util_coo_array_renumber(OJA+zoff,JREN,nnz,flags,flags,flags);
err:
	if(rnz)
		*rnz += nnz;

	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_block_sparse(const struct rsb_mtx_t * mtxAp, void* VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t fr, rsb_coo_idx_t lr, rsb_coo_idx_t fc, rsb_coo_idx_t lc, const rsb_coo_idx_t * IREN, const rsb_coo_idx_t * JREN, rsb_nnz_idx_t *rnz, rsb_flags_t flags )
{
	/*! 
	 * \ingroup gr_internals
	 *
	 * FIXME: IREN/JREN untested
	 * */
	rsb_nnz_idx_t nnz = 0;
	rsb_coo_idx_t ioff = 0,joff = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const struct rsb_mtx_t * submatrix = NULL;

#if RSB_ALLOW_ZERO_DIM
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
		goto ret; /* FIXME: skipping further error checks */
#endif

	if(mtxAp == NULL)
	{
		errval = RSB_ERR_BADARGS;
		goto ret;
	}

	if( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )
		fr--,lr--,fc--,lc--,ioff++,joff++;

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t smi;
		//#pragma omp parallel reduction(|:errval,+,nnz) shared(mtxAp)  RSB_NTC
		RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
		if(RSB_SUBMATRIX_INTERSECTS_BOX(submatrix,fr,lr,fc,lc))
		{
      			const rsb_coo_idx_t fri = RSB_SUBMATRIX_ROWS_INTERSECTION_FIRST(submatrix,fr);
			const rsb_coo_idx_t lri = RSB_SUBMATRIX_ROWS_INTERSECTION_LAST(submatrix,lr);
      			const rsb_coo_idx_t fci = RSB_SUBMATRIX_COLS_INTERSECTION_FIRST(submatrix,fc);
			const rsb_coo_idx_t lci = RSB_SUBMATRIX_COLS_INTERSECTION_LAST(submatrix,lc);
			errval = rsb__do_get_block_sparse_leaf(submatrix,VA,fri,lri,fci,lci,IA,JA,IREN,JREN,&nnz,ioff,joff,flags);
		}
	}
	else
		errval = rsb__do_get_block_sparse_leaf(mtxAp,VA,fr,lr,fc,lc,IA,JA,IREN,JREN,&nnz,ioff,joff,flags);
ret:
	if(rnz)
		*rnz = nnz;

	RSB_DO_ERR_RETURN(errval)
}

rsb_nnz_idx_t rsb__do_get_block_nnz(const struct rsb_mtx_t *mtxAp, rsb_blk_idx_t fr, rsb_blk_idx_t lr, rsb_blk_idx_t fc, rsb_blk_idx_t lc, rsb_flags_t flags, rsb_err_t * errvalp)
{
	/*! 
	 * \ingroup gr_internals
	 * */
	rsb_nnz_idx_t nnz = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	errval = rsb__do_get_block_sparse(mtxAp,NULL,NULL,NULL,fr,lr,fc,lc,NULL,NULL,&nnz,flags);
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return nnz;
}

/* @endcond */
