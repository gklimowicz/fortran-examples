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
/**
 * @file
 * @author Michele Martone
 * @brief
 *
 * Matrix cloning biotech.
 * \internal
 *
 * */
#include "rsb_common.h"

//#define RSB_MTX_REASSIGN(OLD_MTXP,NEW_MTXP) {if(rsb_do_assign(NEW_MTXP,OLD_MTXP)) {RSB_PERR_GOTO(err,RSB_ERRM_ES);}
#define RSB_MTX_REASSIGN(OLD_MTXP,NEW_MTXP) { RSB_MTX_FREE(OLD_MTXP); (OLD_MTXP)=(NEW_MTXP); }

RSB_INTERNALS_COMMON_HEAD_DECLS

void * rsb__clone_area_with_extra(const void *src, size_t csize, size_t bsize, size_t esize)
{
	/*!
	 * \ingroup gr_internals
	 * (m)allocates an area of bsize+csize+esize bytes and copies there csize bytes from src, at offset bsize
	 * esize (extra size) may be zero.
	 */
	rsb_byte_t * dst = NULL;

	if(!src)
		goto ret;

	dst = rsb__malloc(csize+bsize+esize);
	if(!dst)
		goto ret;

	rsb__memcpy(dst+bsize,src,csize);
ret:
	return dst;
}

void * rsb__clone_area(const void *src, size_t size)
{
	/*!
	 * \ingroup gr_internals
	 * \param src the source data pointer
	 * \param size the amount of data to clone
	 * \return the cloned amount, or NULL in case of error
	 *
	 * allocates an area of size bytes and copies there data from src
	 * */
	void * dst = NULL;

	if(!src || size < 1)
		goto ret;

	dst = rsb__clone_area_with_extra(src,size,0,0);
ret:
	return dst;
}

rsb_err_t rsb__util_coo_alloc(void **RSB_RESTRICT VAp, rsb_coo_idx_t ** RSB_RESTRICT IAp, rsb_coo_idx_t ** RSB_RESTRICT JAp, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_bool_t do_calloc)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL,*IA_ = NULL,*JA_ = NULL;

	if( RSB_MATRIX_UNSUPPORTED_TYPE(typecode) )
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		goto err;
	}

	if(do_calloc != RSB_BOOL_TRUE)
	{
		VA_ = rsb__malloc_vector((nnz),typecode),
		IA_ = rsb__malloc(sizeof(rsb_coo_idx_t)*(nnz)),
		JA_ = rsb__malloc(sizeof(rsb_coo_idx_t)*(nnz));
	}
	else
	{
		VA_ = rsb__calloc_vector((nnz),typecode),
		IA_ = rsb__calloc(sizeof(rsb_coo_idx_t)*(nnz)),
		JA_ = rsb__calloc(sizeof(rsb_coo_idx_t)*(nnz));
	}

	if(!VA_ || !IA_ || !JA_)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ENOMEM);
	}

	*VAp = VA_;
	*IAp = IA_;
	*JAp = JA_;
	goto done;
err:
	RSB_CONDITIONAL_FREE(IA_);
	RSB_CONDITIONAL_FREE(JA_);
	RSB_CONDITIONAL_FREE(VA_);
done:
	return errval;
}

rsb_err_t rsb__util_coo_alloc_copy_and_stats(void **RSB_RESTRICT VAp, rsb_coo_idx_t ** RSB_RESTRICT IAp, rsb_coo_idx_t ** RSB_RESTRICT JAp, const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t*RSB_RESTRICT mp, rsb_coo_idx_t*RSB_RESTRICT kp, rsb_nnz_idx_t nnz, rsb_nnz_idx_t ennz, rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t iflags, rsb_flags_t*RSB_RESTRICT flagsp)
{
	/*!
	 * Copies contents of a COO arrays triple to a freshly allocated COO arrays triple.
	 * Size is assumed to be nnz+ennz.
	 * Last ennz elements are not zeroed.
	 *
	 * Flags are determined: RSB_FLAG_UPPER_TRIANGULAR, RSB_FLAG_LOWER_TRIANGULAR.
	 *
	 * TODO: May implement input sanitization or zeroes detection.
	 * TODO: Check for nnz+ennz overflow.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL;
	rsb_coo_idx_t *IA_ = NULL,*JA_ = NULL;

	errval = rsb__util_coo_alloc((void**)(&VA_),&IA_,&JA_,nnz+ennz,typecode,RSB_BOOL_FALSE);
	if(RSB_SOME_ERROR(errval))
		goto err;

	if(!VA && !IA && !JA)
	       	goto nocopy; /* it's ok: alloc only semantics */
	/* TODO: the following shall be made parallel */
	if(mp || kp || ( flagsp && RSB__FLAG_HAS_UNSPECIFIED_TRIANGLE(*flagsp)) )
		errval = rsb__util_coo_copy_and_stats(VA,IA,JA,VA_,IA_,JA_,mp,kp,nnz,typecode,offi,offo,iflags,flagsp);
	else
	{
		errval = rsb__util_coo_copy(VA,IA,JA,VA_,IA_,JA_,nnz,typecode,offi,offo);
		/* ... flags may not always be desired! */
	/*	if(flagsp)
			(*flagsp)|=rsb__util_coo_determine_uplo_flags(IA_,JA_,nnz);*/
	}
nocopy:
	*VAp = VA_;
	*IAp = IA_;
	*JAp = JA_;
	goto done;
err:
	RSB_CONDITIONAL_FREE(IA_);
	RSB_CONDITIONAL_FREE(JA_);
	RSB_CONDITIONAL_FREE(VA_);
done:
	return errval;
}

#if !RSB_WANT_SM_TO_THREAD_MOD_MAPPING
static void * rsb__clone_area_parallel(const void *src, size_t size, size_t n)
{
	/*!
	 * \ingroup gr_internals
	 * \param src the source data pointer
	 * \param size the amount of data to clone
	 * \return the cloned amount, or NULL in case of error
	 *
	 * allocates an area of size bytes and copies there data from src
	 * */
	void * dst = NULL;

	if(!src || size < 1)
		goto ret;

	dst = rsb__malloc(size*n);

	if(!dst)
		goto ret;

	RSB_A_MEMCPY_parallel(dst,src,0,0,n,size);
ret:
	return dst;
}
#endif /* !RSB_WANT_SM_TO_THREAD_MOD_MAPPING */

#if RSB_WANT_BITMAP
static void * rsb__clone_options_t(const struct rsb_options_t *o, rsb_blk_idx_t M_b, rsb_blk_idx_t K_b)
{
	/*!
	 * \ingroup gr_internals
	 * \param o a valid option structure pointer
	 * \return the input pointer
	 *
	 * clones a rsb_options_t structure, deeply
	 *
	 * p.s.: the input rsb_options_t could be NULL. in that case it won't be cloned because there will be no need of doing so.
	 * */
	struct rsb_options_t *no = NULL;

	if(!o)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	/* we allocate a new options structure */
	if(! (no = rsb__clone_area(o,sizeof(*no))))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( o->bitmap)
	{
		no->bitmap = rsb__clone_area(o->bitmap,RSB_BYTES_PER_BITMAP( M_b,K_b));
		if(!no->bitmap)
			RSB_PERR_GOTO(err,RSB_ERRM_ENOMEM);
	}
	return no;
err:
	rsb__destroy_options_t(no);
	return NULL;
}
#endif /* RSB_WANT_BITMAP */

#define RSB_ARE_MTX_COO_CONFORMANT(MTXAP,MTXBP) \
	( ( (MTXAP)->nnz == (MTXBP)->nnz ) && ( (MTXAP)->nr == (MTXBP)->nr ) && ( (MTXAP)->nc == (MTXBP)->nc ) && ( (MTXAP)->typecode == (MTXBP)->typecode ) )

rsb_err_t rsb__mtx_shift_leaf_ptrs(struct rsb_mtx_t *RSB_RESTRICT  mtxCp, const struct rsb_mtx_t *RSB_RESTRICT  mtxAp, long smc)
{
	/* 
	 * Adjust pointers displacements in the matrix tree.
	 * Note: no pointer is being referenced.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t n,si;

	for(n=0;n<smc;++n)
		if(mtxCp[n].nnz) /* If valid. FIXME: IF NOT (E.G. MERGED), SHALL BE COMPLETELY ZEROED. */
		for(si=0;si<RSB_FOUR;++si)
			if(mtxAp[n].sm[si])
			{
				RSB_PTR_SHIFT(mtxCp[n].sm[si],mtxAp,mtxCp,(struct rsb_mtx_t*));
				//mtxCp[n].sm[si] = mtxCp+(mtxAp[n].sm[si]-mtxAp);
			/*	RSB_STDOUT("%03d/%03d: %p\n",n,si,mtxCp[n].sm[si]); */
			}
	return errval;
}

#if !RSB_AT_DESTROYS_MTX
rsb_err_t rsb__mtx_transplant_from_clone(struct rsb_mtx_t ** mtxDpp, struct rsb_mtx_t * mtxSp)
{
	/* 
	 Moves the inner content of mtxSp to mtxDp.
	 Shall free mtxSp at the end and not change the outer allocation status of mtxSp.
	 Shall work even if any of the two matrices is in-place.
	 Can only work if matrices match in size (nonzeroes, rows, columns, ...).

	 Untested.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxDp = NULL;
	struct rsb_mtx_t *fsmS = NULL;
	struct rsb_mtx_t *fsmD = NULL;
	void * VS = NULL,* VD = NULL;
	rsb_coo_idx_t * IS = NULL, *JS = NULL;
	rsb_coo_idx_t * ID = NULL, *JD = NULL;
	rsb_long_t smc = 0;

	if( ( !mtxDpp) || ( !*mtxDpp) || (!mtxSp) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err, RSB_ERRM_E_MTXAP);
	}

	mtxDp = *mtxDpp;

	if( ! ( RSB_ARE_MTX_COO_CONFORMANT( mtxSp, mtxDp ) ) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}

#if 0
	RSB_STDOUT("		==== CLONING: ==== \n");
	RSB_STDOUT("will transplant: \n");
	RSB_STDOUT(RSB_PRINTF_MTX_SUMMARY_ARGS(mtxSp)),
	RSB_STDOUT("to: \n");
	RSB_STDOUT(RSB_PRINTF_MTX_SUMMARY_ARGS(mtxDp)),
	RSB_STDOUT("\n");
	RSB_STDOUT("S ip: : %x \n",(RSB_DO_FLAG_HAS(mtxSp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS)));
	RSB_STDOUT("D ip: : %x \n",(RSB_DO_FLAG_HAS(mtxDp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS)));
#endif

	fsmS = rsb__do_get_first_submatrix(mtxSp);
	fsmD = rsb__do_get_first_submatrix(mtxDp);
	VS = fsmS->VA, VD = fsmD->VA;
	IS = fsmS->bpntr, JS = fsmS->bindx;
	ID = fsmD->bpntr, JD = fsmD->bindx;

	errval = rsb__util_coo_copy(VS, IS, JS, VD, ID, JD, mtxSp->nnz, mtxSp->typecode, 0, 0);
	if(RSB_SOME_ERROR(errval))
		goto err;

	// correct mtxSp pointers recursively with the three offsets

	/* get rid of the source COO arrays */
	if(RSB_DO_FLAG_HAS(mtxSp->flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{	
		fsmS->VA = NULL;
		fsmS->bindx = NULL;
		fsmS->bpntr = NULL;
	}
	else
	{
		RSB_CONDITIONAL_FREE(fsmS->VA);
		RSB_CONDITIONAL_FREE(fsmS->bindx);
		RSB_CONDITIONAL_FREE(fsmS->bpntr);
	}
	smc = 1 + rsb__submatrices_max_ptr_diff(mtxSp);

	if(RSB_DO_FLAG_HAS(mtxDp->flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
		RSB_DO_FLAG_ADD(mtxSp->flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
	else
		; /* mtxSp will retain its original flags with these regards  */

	/* set new matrix arrays by translating the submatrix pointers */
	rsb__do_set_in_place_submatrices_offsets(mtxSp, smc, VD, ID, JD, mtxDp->el_size);

	/* mtxDp->all_leaf_matrices shall be ok... */

	/* free now unnecessary original destination matrix pointer */
	fsmD->VA = NULL;
	fsmD->bindx = NULL;
	fsmD->bpntr = NULL;
	RSB_MTX_FREE(mtxDp);

	/* overwrite the output pointer */
	mtxDp = mtxSp;
	*mtxDpp = mtxDp;
	mtxSp = NULL;

#if 0
	RSB_STDOUT("obtained: \n");
	RSB_STDOUT(RSB_PRINTF_MTX_SUMMARY_ARGS(mtxDp)),
	RSB_STDOUT("\n");
#endif
err:
	return errval;
}
#endif  /* RSB_AT_DESTROYS_MTX */

#if 0
static rsb_err_t rsb_do_assign(struct rsb_mtx_t * mtxBp, const struct rsb_mtx_t * mtxAp)
{
	rsb_submatrix_idx_t i,j;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * submatrix=NULL;

	if(!mtxBp || !mtxAp)
		goto err;

	rsb__destroy_inner(mtxBp);

	rsb__memcpy(mtxBp,mtxAp,sizeof(*mtxAp));
	rsb__init_blank_pointers(mtxBp);

	if(rsb__clone_inner(mtxAp,mtxBp)==NULL)
		goto err;

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
	{
		if((mtxBp->sm[i*2+j]=rsb__mtx_clone_simple(submatrix))==NULL)
			goto err;
	}

#if RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS
	if(mtxAp->all_leaf_matrices)
	{
		mtxBp->all_leaf_matrices=NULL;
		errval = rsb__get_array_of_leaf_matrices(mtxBp,&mtxBp->all_leaf_matrices,&mtxBp->all_leaf_matrices_n);
		if(RSB_SOME_ERROR(errval))
			goto errr;
	}
#endif /* RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS */
	goto errr;
err:
	errval = RSB_ERR_GENERIC_ERROR;
errr:
	return errval;
}
#endif

static size_t rsb__submatrices_max_ptr_diff_inner(const struct rsb_mtx_t * mtxRp, const struct rsb_mtx_t * mtxAp)
{
	/*!
	 * 	\ingroup gr_internals
	 * 	Note: this makes only sense if submatrices are allocated in one block.
	 */
	size_t md = 0;
	rsb_submatrix_idx_t i,j;
	const struct rsb_mtx_t * submatrix;

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	{
		if(submatrix)
		{
			size_t sd = rsb__submatrices_max_ptr_diff_inner(mtxRp,submatrix);

			md = RSB_MAX(md,sd);
			md = RSB_MAX(md,(submatrix-mtxRp));
		}
	}
	return md;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
size_t rsb__submatrices_max_ptr_diff_noholes(const struct rsb_mtx_t * mtxAp)
{
	size_t md = 0;
	rsb_submatrix_idx_t i,j;
	const struct rsb_mtx_t * submatrix;

	if(!mtxAp)
	{
		return 0;
	}

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	{
		if(submatrix)
		{
			const size_t sd = rsb__submatrices_max_ptr_diff(submatrix);
			md = RSB_MAX(md,sd+(submatrix-mtxAp));
		}
	}
	return md;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

size_t rsb__submatrices_max_ptr_diff(const struct rsb_mtx_t * mtxAp)
{
	/*!
	 * 	Note: makes only sense if submatrices are allocated in one block.
	 */
#if 0
	return rsb__submatrices_max_ptr_diff_noholes(mtxAp);
#else
	return rsb__submatrices_max_ptr_diff_inner(mtxAp, mtxAp);
#endif
}

void * rsb__clone_area_guided(void * RSB_RESTRICT dst, const void *RSB_RESTRICT src, size_t size, size_t nmemb, const struct rsb_mtx_t *RSB_RESTRICT mtxAp, const rsb_thread_t * RSB_RESTRICT cta, const rsb_thread_t nct, rsb_err_t * errvalp)
{
	/*
		Initializes, allocating if needed, and/or copying in parallel, using specified chunk sizes and array.
		If dst supplied, will use it, otherwise will allocate one.
		If src supplied, will use it, otherwise will only zero the arrays.
		If mtxAp == NULL, then must also be  cta == NULL && nct == 0.
		Returns either dst or the newly allocated area address. 
	*/
	rsb_err_t errval = RSB_ERR_BADARGS;

	if(size*nmemb == 0)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(mtxAp == NULL && (cta != NULL || nct != 0) )
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(dst == NULL && ( dst = rsb__malloc(size*nmemb) ) == NULL )
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ENOMEM);
	}

#if RSB_WANT_OMP_RECURSIVE_KERNELS
	if(mtxAp == NULL)
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	{
		/* master thread */
		RSB_A_MEMCPY(dst,src,0,0,nmemb,size);
		errval = RSB_ERR_NO_ERROR;
		goto done;
	}

#if RSB_WANT_OMP_RECURSIVE_KERNELS
	if(cta == NULL)
	{
		RSB_DEBUG_ASSERT((mtxAp)->all_leaf_matrices_n);
#pragma omp parallel shared(mtxAp) 
{
	rsb_submatrix_idx_t smi;
	const rsb_thread_t omt = omp_get_max_threads(), otn = omp_get_thread_num();
	/* auto: each submatrix a round robin thread */
	for(smi=0;smi<mtxAp->all_leaf_matrices_n;++smi) /* FIXME: make this an OpenMP-friendly macro */
	if( ( smi % omt ) == otn )
	{
		const struct rsb_mtx_t * submatrix = mtxAp->all_leaf_matrices[smi].mtxlp;
		const size_t off = submatrix->nzoff;
		const rsb_nnz_idx_t nnz = submatrix->nnz;

		if(src)
			RSB_A_MEMCPY(dst,src,off,off,nnz,size);
		else
			RSB_A_BZERO(dst,off,nnz,size);
	}
}
#pragma omp barrier
		errval = RSB_ERR_NO_ERROR;
		goto done;
	}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
#pragma omp parallel shared(mtxAp) RSB_NTC 
{
	/* guided: each submatrix touched by a specified thread */
	const rsb_thread_t otn = omp_get_thread_num();
	rsb_submatrix_idx_t cti; /* thread index */

	for(cti=0;cti<nct;++cti)
	if( otn == cta[cti] )
	{
		const struct rsb_mtx_t * submatrix = mtxAp->all_leaf_matrices[cti].mtxlp;
		const size_t off = submatrix->nzoff;
		const rsb_nnz_idx_t nnz = submatrix->nnz;

		if(src)
			RSB_A_MEMCPY(dst,src,off,off,nnz,size);
		else
			RSB_A_BZERO(dst,off,nnz,size);
	}
}
#pragma omp barrier
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	errval = RSB_ERR_NO_ERROR;
	goto done;
err:
	dst = NULL;
done:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return dst;
}

static struct rsb_mtx_t *rsb__mtx_clone_simple_extra(const struct rsb_mtx_t *mtxAp, rsb_submatrix_idx_t esmc)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Clones a whole matrix, retaining the same submatrices structure.
	 * TODO: need better error handling.
	 * FIXME : Unfinished: only the RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS case is handled.
	 */
	struct rsb_mtx_t *mtxCp = NULL;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
#if RSB_ALLOW_INTERNAL_GETENVS
	rsb_time_t ct = RSB_TIME_ZERO;
	rsb_time_t mact = RSB_TIME_ZERO;
#endif /* RSB_ALLOW_INTERNAL_GETENVS */

	if(!mtxAp)
	{
		RSB_PERR_GOTO(nerr,RSB_ERRM_ES);
	}

	flags = mtxAp->flags;
#if RSB_ALLOW_INTERNAL_GETENVS
	ct = -rsb_time();
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
	RSB_DO_FLAG_DEL(flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS))
	{
		if(rsb__is_root_matrix(mtxAp))
		{
			/* this is a trick: fitting the whole recursive matrix in three arrays */
			/* void *IA = NULL,*JA = NULL,*VA = NULL; */
			const rsb_nnz_idx_t nnz = mtxAp->nnz;
			/* rsb_bool_t is_bio = rsb__do_is_matrix_binary_loaded(mtxAp);*/ /* binary I/O matrix (20120930 FIXME why is this unused ?) */
			/* rsb_long_t smc = rsb__submatrices(mtxAp); */
			const rsb_long_t smc = 1 + rsb__submatrices_max_ptr_diff(mtxAp);
			const struct rsb_mtx_t *fsm = rsb__do_get_first_submatrix(mtxAp);
			rsb_err_t errval = RSB_ERR_NO_ERROR;
			void * VA = NULL;
			rsb_coo_idx_t * bpntr = NULL, * bindx = NULL;

			RSB_DEBUG_ASSERT(1+rsb__submatrices_max_ptr_diff(mtxAp)>=rsb__submatrices(mtxAp));
			
			/* RSB_STDOUT("MAX PTR DIFF: %d  SUBM COUNT:%d\n",1+rsb__submatrices_max_ptr_diff(mtxAp), rsb__submatrices(mtxAp)); */

			/* mtxCp = rsb__clone_area(mtxAp,sizeof(struct rsb_mtx_t)*(smc+esmc)); */
			mtxCp = rsb__clone_area_with_extra(mtxAp,sizeof(struct rsb_mtx_t)*(smc),0,sizeof(struct rsb_mtx_t)*(esmc));

			if(!mtxCp)
			{
				RSB_PERR_GOTO(nerr,RSB_ERRM_PAL);
			}
			
			errval = rsb__mtx_shift_leaf_ptrs(mtxCp, mtxAp, smc);

			mtxCp->all_leaf_matrices = NULL;
#if 0
			errval = rsb__get_array_of_leaf_matrices(mtxCp,&(mtxCp->all_leaf_matrices),&(mtxCp->all_leaf_matrices_n));
#else
			mtxCp->all_leaf_matrices = rsb__clone_area_with_extra(mtxAp->all_leaf_matrices,sizeof(mtxAp->all_leaf_matrices[0])*(mtxCp->all_leaf_matrices_n),0,0);
			if( mtxCp->all_leaf_matrices == NULL )
			{
				errval = RSB_ERR_ENOMEM;
			}
			else
			{
				rsb_submatrix_idx_t si;

				for(si=0;si<mtxCp->all_leaf_matrices_n;++si)
					RSB_PTR_SHIFT(mtxCp->all_leaf_matrices[si].mtxlp,mtxAp,mtxCp,(struct rsb_mtx_t*));
				for(si=0;si<mtxCp->all_leaf_matrices_n;++si)
					RSB_DO_FLAG_DEL(mtxCp->all_leaf_matrices[si].mtxlp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
				RSB_DO_FLAG_DEL(mtxCp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
			}
#endif
			if(RSB_SOME_ERROR(errval))
			{
				RSB_PERR_GOTO(nerr,RSB_ERRM_NAOL);
			}
			/* FIXME: and what when nnz==0 and VA!=NULL ? */
			mtxCp->bindx = NULL; mtxCp->bpntr = NULL; mtxCp->VA = NULL;
			if(nnz)
			{
#if RSB_ALLOW_INTERNAL_GETENVS
				mact = -rsb_time();
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
				RSB_ASSERT( fsm->el_size );
				RSB_ASSERT( fsm->bindx );
				RSB_ASSERT( fsm->bpntr );
				RSB_ASSERT( fsm->VA );

#if RSB_WANT_SM_TO_THREAD_MOD_MAPPING
			if( 1 /*rsb__util_atoi(getenv("RSB_CLONE_SERIAL")) */ )
			{
				/* rsb_thread_t nct = (rsb__util_atoi(getenv("RSB_CLONE_SERIAL"))) - 1; */
				rsb_thread_t nct = 0;
				bindx = rsb__clone_area_guided(NULL,fsm->bindx,sizeof(rsb_coo_idx_t),nnz,mtxAp,NULL,nct,&errval);
				bpntr = rsb__clone_area_guided(NULL,fsm->bpntr,sizeof(rsb_coo_idx_t),nnz,mtxAp,NULL,nct,&errval);
				VA    = rsb__clone_area_guided(NULL,fsm->VA   ,mtxAp->el_size,       nnz,mtxAp,NULL,nct,&errval);
			}
#endif
#if !RSB_WANT_SM_TO_THREAD_MOD_MAPPING
			{
				bindx = rsb__clone_area_parallel(fsm->bindx,sizeof(rsb_coo_idx_t),nnz);
				bpntr = rsb__clone_area_parallel(fsm->bpntr,sizeof(rsb_coo_idx_t),nnz);
				VA = rsb__clone_area_parallel(fsm->VA,mtxAp->el_size,nnz);
			}
#endif /* !RSB_WANT_SM_TO_THREAD_MOD_MAPPING */

#if RSB_ALLOW_INTERNAL_GETENVS
				mact += rsb_time();
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
				if(!bindx || !bpntr || !VA || !mtxAp->el_size)
				{
					RSB_ASSERT( mtxCp->el_size );
					RSB_ASSERT( bindx );
					RSB_ASSERT( bpntr );
					RSB_ASSERT( VA );
					RSB_PERR_GOTO(ierr,RSB_ERRM_NNTC);
				}
			}
#if !RSB_ALLOW_EMPTY_MATRICES
			else
			if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
			{
				// ok
			}
			else
			{
				RSB_PERR_GOTO(ierr,RSB_ERRM_NDIANN);
			}
#endif /* RSB_ALLOW_EMPTY_MATRICES */
			rsb__do_set_in_place_submatrices_offsets(mtxCp,smc,VA,bpntr,bindx,mtxCp->el_size);

			if(!smc)
				mtxCp->bindx = bindx,
				mtxCp->bpntr = bpntr,
				mtxCp->VA = VA;

			/* note: the cloned matrix won't have the is_bio property */
			goto ret;
ierr:
			RSB_ERROR(RSB_ERRM_ES);
			if(mtxCp)
			{
				RSB_CONDITIONAL_FREE(mtxCp->bpntr);
				RSB_CONDITIONAL_FREE(mtxCp->bindx);
				RSB_CONDITIONAL_FREE(mtxCp->VA);
				RSB_CONDITIONAL_FREE(mtxCp->all_leaf_matrices);
			}
			RSB_CONDITIONAL_FREE(mtxCp);
			goto ret;
		}
		else
		{
			RSB_PERR_GOTO(ret,"no cloning possible for a non root matrix!\n");
			/* no cloning for a non root */
		}
	}
	else
	{
		/* we allocate a new matrix structure */
		mtxCp = rsb__clone_area(mtxAp,sizeof(*mtxCp));
	
		if(!mtxCp)
			{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	
		rsb__init_blank_pointers(mtxCp);

		RSB_MTX_REASSIGN(mtxCp,(struct rsb_mtx_t*)mtxAp);
		goto ret;
	}
err:
	RSB_MTX_FREE(mtxCp);
nerr:
	RSB_CONDITIONAL_FREE(mtxCp);
ret:
#if RSB_ALLOW_INTERNAL_GETENVS
	ct += rsb_time();
	if( rsb__util_atoi(getenv("RSB_MTX_CLONE_STATS") ) != 0)
	if(mtxCp)
	{
		const size_t szv = rsb__get_sizeof(mtxCp);
		RSB_STDOUT("Cloned a %zd nnz, %zd bytes matrix in %0.2lgs (%0.3lg MiB/s x 2 = r+w); of which %0.2lgs for the main arrays.\n",
				(size_t)(mtxCp->nnz),szv,ct,(((rsb_time_t)szv)/ct)/RSB_MEGABYTE,mact);
	}
#endif /* RSB_ALLOW_INTERNAL_GETENV */
 	if(mtxCp)
		RSB_DO_FLAG_DEL(mtxCp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
	return mtxCp;
}

struct rsb_mtx_t *rsb__mtx_clone_simple(const struct rsb_mtx_t *mtxAp)
{
	return rsb__mtx_clone_simple_extra(mtxAp, 0);
}

rsb_err_t rsb__clone_coo(const struct rsb_mtx_t * mtxAp, rsb_trans_t transA, const void *alphap, rsb_type_t typecode, struct rsb_coo_mtx_t*dcoop, rsb_flags_t flags/*, rsb_extff_t eflags*/)
{
	/* 
	   TODO: may integrate here fortran indices handling and so on
	   TODO: missing checks for indices overflow
	   TODO: shall juggle appropriately with 'sorted' flags
	*/
	rsb_flags_t cflags = RSB_FLAG_NOFLAGS;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t dels = 0;
	rsb_coo_idx_t ioff,joff;
	struct rsb_coo_mtx_t dcoo,scoo;
	const rsb_bool_t expsymm = (RSB_DO_FLAG_HAVE_XOR(flags,mtxAp->flags,RSB_FLAG_SYMMETRIC) && !RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_DIAGONAL));
	const rsb_bool_t expherm = (RSB_DO_FLAG_HAVE_XOR(flags,mtxAp->flags,RSB_FLAG_HERMITIAN) && !RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_DIAGONAL));

	RSB_BZERO_P(&scoo);
	RSB_BZERO_P(&dcoo);
	ioff = joff = ( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )?1:0;
	scoo.nr = dcoo.nr = mtxAp->nr;
	scoo.nc = dcoo.nc = mtxAp->nc;
	dcoo.nnz = scoo.nnz = mtxAp->nnz;
	if(expsymm || expherm)
		dcoo.nnz *= 2;/* of course, this is overkill in the case of a diagonal matrix */
	if(RSB_DO_FLAG_HAVE_XOR(flags,mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
		dels = RSB_MIN(dcoo.nr,dcoo.nc);
	dcoo.nnz += dels;
	//if(dels)
	//	RSB_STDOUT("on diag: %d\n",dels);

	scoo.typecode = mtxAp->typecode;
	dcoo.typecode = typecode;
	if(dcoo.nnz>0)
	{
		if(rsb__allocate_coo_matrix_t(&dcoo)!=&dcoo)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(ierr,RSB_ERRM_PAL);
	       	}
		if(rsb__allocate_coo_matrix_t(&scoo)!=&scoo)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(ierr,RSB_ERRM_PAL);
	       	}
		errval = rsb__do_get_coo_noalloc(mtxAp,scoo.VA,dcoo.IA,dcoo.JA,NULL,/*mtxAp->*/flags);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(ierr,RSB_ERRM_NL);
		}
		errval = rsb__do_copy_converted_scaled(scoo.VA,dcoo.VA,alphap,mtxAp->typecode,typecode,mtxAp->nnz,transA);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(ierr,RSB_ERRM_NL);
	       	}
		if(expsymm || expherm)
			RSB_COO_MEMCPY(dcoo.VA,dcoo.IA,dcoo.JA,dcoo.VA,dcoo.JA,dcoo.IA,scoo.nnz,0,scoo.nnz,RSB_SIZEOF(typecode));
		if(expherm)
			rsb__util_do_conjugate(((rsb_byte_t*)(dcoo.VA))+(RSB_SIZEOF(typecode)*scoo.nnz),typecode,scoo.nnz);
		if(RSB_DO_FLAG_HAVE_XOR(flags,mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
			rsb__do_fill_with_diag(dcoo.VA,dcoo.IA,dcoo.JA,ioff,joff,dcoo.nnz-dels,typecode,dels);
		rsb__destroy_coo_matrix_t(&scoo);
		RSB_BZERO_P(&scoo);
	}

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_TRIANGULAR))
		RSB_DO_FLAG_ADD(cflags,RSB_FLAG_UPPER_TRIANGULAR);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_TRIANGULAR))
		RSB_DO_FLAG_ADD(cflags,RSB_FLAG_LOWER_TRIANGULAR);
	if(RSB_DO_FLAG_HAVE_XOR(flags,mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
		RSB_DO_FLAG_ADD(cflags,RSB_FLAG_UNIT_DIAG_IMPLICIT);
	if((cflags != RSB_FLAG_NOFLAGS) || expsymm || expherm)
	{
		rsb__util_sort_row_major_inner(dcoo.VA,dcoo.IA,dcoo.JA,dcoo.nnz,dcoo.nr,dcoo.nc,typecode,flags);
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT);
		dcoo.nnz = rsb__weed_out_duplicates(dcoo.IA,dcoo.JA,dcoo.VA,dcoo.nnz,typecode,flags);
		errval = rsb__do_cleanup_nnz(dcoo.VA,dcoo.IA,dcoo.JA,dcoo.nnz,0,0,dcoo.nr,dcoo.nc,&dcoo.nnz,dcoo.typecode,cflags); /* FIXME: are we using roff,coff well here ? */
	}
	if(RSB_SOME_ERROR(errval))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(ierr,RSB_ERRM_CP);
       	}

	if(RSB_DOES_TRANSPOSE(transA))
		rsb__transpose_coo_matrix_t(&dcoo);
	*dcoop = dcoo;
ierr:
	return errval;
}

rsb_err_t rsb__mtx_clone(struct rsb_mtx_t ** mtxBpp, rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * clones a rsb_mtx_t structure, deeply
	 * This routine may/shall be optimized in plenty of ways, in the future.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxCp = NULL;

	if( (!mtxAp) || (!mtxBpp) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"user did not supply a valid matrix pointer\n");
	}

	if( typecode != RSB_NUMERICAL_TYPE_SAME_TYPE && RSB_MATRIX_UNSUPPORTED_TYPE(typecode) )
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,"user supplied wrong flags (RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS is illegal here)\n");
	}

	if( flags == RSB_FLAG_IDENTICAL_FLAGS )
		flags = mtxAp->flags;
	if(typecode == RSB_NUMERICAL_TYPE_SAME_TYPE)
		typecode = mtxAp->typecode;
	RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);/* unnecessary here */
	RSB_DO_FLAG_DEL(flags,RSB_FLAG_NON_ROOT_MATRIX);

	/* what about RSB_FLAG_DISCARD_ZEROS ? */
	/* what about many other structural flags ? */

	if(((!alphap) || RSB_IS_ELEMENT_ONE(alphap,mtxAp->typecode)) &&
			(flags == mtxAp->flags ) &&  /* FIXME: this condition is unnecessarily strict and cause inefficiencies */
		/* RSB_DOES_NOT_TRANSPOSE(transA) && */ (typecode==mtxAp->typecode) )
	{
		if( (*mtxBpp) != mtxAp)
			mtxCp = rsb__mtx_clone_simple(mtxAp);
		else
			mtxCp = *mtxBpp;
		if( transA == RSB_TRANSPOSITION_C )
			errval = rsb__do_transpose(&mtxCp,RSB_BOOL_TRUE);
		else
		if( transA == RSB_TRANSPOSITION_T )
			errval = rsb__do_transpose(&mtxCp,RSB_BOOL_FALSE);
		if( (*mtxBpp) == mtxAp)
		{
			*mtxBpp = mtxCp;
			goto ok;
		}
	}
	else
	{
		struct rsb_coo_mtx_t dcoo;
		RSB_BZERO_P(&dcoo);
#if 0
		struct rsb_coo_mtx_t scoo;
		RSB_BZERO_P(&scoo);
		scoo.nr = dcoo.nr = mtxAp->nr;
		scoo.nc = dcoo.nc = mtxAp->nc;
		dcoo.nnz = scoo.nnz = mtxAp->nnz;
		scoo.typecode = mtxAp->typecode;
		dcoo.typecode = typecode;
		if(mtxAp->nnz>0)
		{
		if(rsb__allocate_coo_matrix_t(&dcoo)!=&dcoo)
		{
		       	errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(ierr,RSB_ERRM_PAL);
		}
		if(rsb__allocate_coo_matrix_t(&scoo)!=&scoo)
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(ierr,RSB_ERRM_PAL);
		}
		errval = rsb__do_get_coo_noalloc(mtxAp,scoo.VA,dcoo.IA,dcoo.JA,NULL,mtxAp->flags);
		rsb__do_copy_converted_scaled(scoo.VA,dcoo.VA,alphap,mtxAp->typecode,typecode,mtxAp->nnz,transA);
		rsb__destroy_coo_matrix_t(&scoo);
		RSB_BZERO_P(&scoo);
		}
		if(RSB_DOES_TRANSPOSE(transA))
			rsb__transpose_coo_matrix_t(&dcoo);
#else
		errval = rsb__clone_coo(mtxAp,transA,alphap,typecode,&dcoo,flags);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_PERR_GOTO(ierr,RSB_ERRM_NL);
	       	}
#endif
		mtxCp = rsb__do_mtx_alloc_from_coo_inplace(dcoo.VA,dcoo.IA,dcoo.JA,dcoo.nnz,dcoo.typecode,dcoo.nr,dcoo.nc,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags,NULL);
		if(mtxCp)
			RSB_DO_FLAG_DEL(mtxCp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
ierr:
		if(!mtxCp)
			rsb__destroy_coo_matrix_t(&dcoo);
	}

	if( (*mtxBpp) == NULL )
	{
		*mtxBpp = mtxCp;
	}
	else
	{
		RSB_MTX_REASSIGN(*mtxBpp,mtxCp);
	}
ok:
err:
	return errval;
}

#if 0
void * rsb__clone_inner(const struct rsb_mtx_t *mtxAp, struct rsb_mtx_t *mtxCp)
{
	/*!
	 * \ingroup gr_internals
	 * clones a rsb_mtx_t structure, deeply
	 *
	 * \param matrix valid matrix pointer (to an empty mtxAp)
	 * \param mtxCp valid matrix pointer
	 * \return a pointer to the cloned structure (mtxCp) in case of success, NULL otherwise
	 *
	 * \note matrix flags are largely ignored in this function.
	 **/

	if(!mtxAp || !mtxCp)
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}

#if RSB_WANT_BITMAP
	/* we allocate a new options structure */
	mtxCp->options = rsb__clone_options_t(mtxAp->options,mtxAp->M_b,mtxAp->K_b);

	if(! mtxCp->options && mtxAp->options )
	{RSB_PERR_GOTO(err_opt,RSB_ERRM_ES);}
#endif /* RSB_WANT_BITMAP */
	if( mtxAp->rpntr && (mtxAp->flags & RSB_FLAG_OWN_PARTITIONING_ARRAYS))
	{
		mtxCp->rpntr = rsb__clone_area(mtxAp->rpntr,sizeof(rsb_coo_idx_t)*(mtxAp->M_b+1));
		if(!mtxCp->rpntr)
		{RSB_PERR_GOTO(err_rpntr,RSB_ERRM_ES);}
	}
	else
		mtxCp->rpntr = mtxAp->rpntr;

	if( mtxAp->cpntr && (mtxAp->flags & RSB_FLAG_OWN_PARTITIONING_ARRAYS))
	{
		mtxCp->cpntr = rsb__clone_area(mtxAp->cpntr,sizeof(rsb_coo_idx_t)*(mtxAp->K_b+1));
		if(!mtxCp->cpntr)
		{RSB_PERR_GOTO(err_cpntr,RSB_ERRM_ES);}
	}
	else
		mtxCp->cpntr = mtxAp->cpntr;

#if RSB_WANT_DBC
	if( mtxAp->bindx)
	{
		mtxCp->bindx = rsb__clone_area(mtxAp->bindx,sizeof(rsb_nnz_idx_t)*(mtxAp->block_count+1));
		if(!mtxCp->bindx)
			{RSB_PERR_GOTO(err_bindx,RSB_ERRM_ES);}
	}

	if( mtxAp->indptr)
	{
		mtxCp->indptr = rsb__clone_area(mtxAp->indptr,sizeof(rsb_nnz_idx_t)*(mtxAp->block_count+1));
		if(!mtxCp->indptr)
			{RSB_PERR_GOTO(err_indptr,RSB_ERRM_ES);}
	}
#endif

	if( mtxAp->bpntr)
	{
		mtxCp->bpntr = rsb__clone_area(mtxAp->bpntr,sizeof(rsb_nnz_idx_t)*(mtxAp->Mdim+1));
		if(!mtxCp->bpntr)
			{RSB_PERR_GOTO(err_bpntr,RSB_ERRM_ES);}
	}

#if RSB_WANT_BITMAP
	if( mtxAp->VA)
	{
		mtxCp->VA = rsb__clone_area(mtxAp->VA,(RSB_TOTAL_BLOCK_BYTES(mtxAp,mtxAp->options)));
		if(!mtxCp->VA)
			{RSB_PERR_GOTO(err_va,RSB_ERRM_ES);}
	}
#endif
	goto ret;

#if RSB_WANT_BITMAP
err_va:
	if( mtxAp->VA)
		RSB_CONDITIONAL_FREE(mtxCp->VA);
#endif
err_bpntr:
	if( mtxAp->bpntr )
		RSB_CONDITIONAL_FREE(mtxCp->bpntr);
err_indptr:
	if( mtxAp->indptr )
		RSB_CONDITIONAL_FREE(mtxCp->indptr);
err_bindx:
	if( mtxAp->bindx )
		RSB_CONDITIONAL_FREE(mtxCp->bindx);
err_cpntr:
	if( mtxAp->cpntr && (mtxAp->flags & RSB_FLAG_OWN_PARTITIONING_ARRAYS))
		RSB_CONDITIONAL_FREE(mtxCp->cpntr);
err_rpntr:
	if( mtxAp->rpntr && (mtxAp->flags & RSB_FLAG_OWN_PARTITIONING_ARRAYS))
		RSB_CONDITIONAL_FREE(mtxCp->rpntr);
#if RSB_WANT_BITMAP
err_opt:
	RSB_CONDITIONAL_FREE(mtxCp->options);
#endif /* RSB_WANT_BITMAP */
err:
	mtxCp = NULL;
ret:
	return mtxCp;
}
#endif

/* @endcond */
