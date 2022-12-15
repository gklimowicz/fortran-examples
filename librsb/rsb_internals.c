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
 * Low level routines and tools for our sparse matrix formats implementations.
 * \internal
 *
 * */
#include "rsb_common.h"
#include "rsb_util.h"
#include "rsb.h"
#include "rsb_types.h"
#include "rsb_unroll.h"

#define RSB_WANT_NULL_ALLOCATED_ZERO_NNZ_COO_MATRICES_ARRAYS 1 /* 20110419 a bugfix for the nnz == 0 case vs realloc and memory counters */
#define RSB_TOKLEN 16
#define RSB_WANT_ZERO_ON_DESTROY 0	/* a useful debug option */
#define rsb__strcpy strcpy

#define RSB_ILLEGAL_FLAGS	0xFFFFFFFF
#define RSB_WANT_DOUBLE_MATRIX_FREE_DETECT 0

#if RSB_WANT_DOUBLE_MATRIX_FREE_DETECT
#define RSB_CONDITIONAL_FREE_MTXAP(MTXAP) if(MTXAP){ if( (MTXAP)->flags == RSB_ILLEGAL_FLAGS ) { RSB_ERROR("Probable attempt to free already freed matrix at %p !\n",(MTXAP)); } (MTXAP)->flags = RSB_ILLEGAL_FLAGS; RSB_CONDITIONAL_FREE(MTXAP); }
#else
#define RSB_CONDITIONAL_FREE_MTXAP(MTXAP) RSB_CONDITIONAL_FREE(MTXAP)
#endif

RSB_INTERNALS_COMMON_HEAD_DECLS

#if RSB_WANT_OMP_RECURSIVE_KERNELS
#define rsb_get_thread_num() rsb_global_session_handle.rsb_want_threads
#else  /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#define rsb_get_thread_num() (1)
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void * rsb__init_options_t(struct rsb_options_t *o)
{
	/*!
	 * \ingroup gr_internals
	 * initializes a rsb_options_t struct to default 'vanilla' values
	 * \return the input address
	 */
	if(!o)
		goto err;
	RSB_BZERO_P(o);
err:
	return o;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if RSB_WANT_BITMAP
void * rsb__destroy_options_t(struct rsb_options_t *o)
{
	/*!
	 * \ingroup gr_internals
	 * frees the given options structure and any of its allocated arrays
	 *
	 * \return the input address
	 */
	if(!o)
		goto err;
	RSB_CONDITIONAL_FREE(o->bitmap);
	RSB_CONDITIONAL_FREE(o);
err:
	return o;
}
#endif /* RSB_WANT_BITMAP */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
const void * rsb__is_valid_options_t(const struct rsb_options_t *o, rsb_coo_idx_t m, rsb_coo_idx_t k)
{
	/*!
	 * \ingroup gr_internals
	 * checks if the structure members which are to be used for 
	 * the creation of a matrix have meaningful values.
	 *
	 * \return the input address in case of success, NULL otherwise
	 * */
	if(!o)
		goto err;

	if(RSB_INVALID_COO_INDEX(m) || RSB_INVALID_COO_INDEX(k))
		goto err;

	/*
	 * someday:
	 *
	 * if(already_initialized) { 
	 * 	if(mtxAp->el_size<1)		goto err;
	 *	if(! o->bitmap)			goto err;
	 *	....
	 */
	return o;
err:
	return NULL;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void * rsb__reallocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp, rsb_nnz_idx_t nnnz)
{
	/*!
	 * \ingroup gr_internals
	 * \return On success, return the input address; on failure, NULL.
	 */
	void * check = NULL;

	if(!cmp)
		goto err;

	if( nnnz == 0 && RSB_WANT_NULL_ALLOCATED_ZERO_NNZ_COO_MATRICES_ARRAYS )
	{
		cmp->IA = NULL;
		cmp->JA = NULL;
		cmp->VA = NULL;
		goto done;
	}

	check = rsb__realloc(cmp->IA,sizeof(rsb_coo_idx_t)*(nnnz));
	if(!check)
		goto err;
	cmp->IA = check;

	check = rsb__realloc(cmp->JA,sizeof(rsb_coo_idx_t)*(nnnz));
	if(!check)
		goto err;
	cmp->JA = check;

	check = rsb__realloc_vector(cmp->VA,nnnz,cmp->typecode);
	if(!check)
		goto err;
	cmp->VA = check;

	if(!cmp->IA || !cmp->JA || !cmp->VA)
		goto cerr;
done:
	cmp->nnz = nnnz;
	return cmp;
cerr:
	/* critical error (should not happen) */
	if(cmp)
	{
		RSB_CONDITIONAL_FREE(cmp->IA);
		RSB_CONDITIONAL_FREE(cmp->JA);
		RSB_CONDITIONAL_FREE(cmp->VA);
	}
err:
	return NULL;
}

void * rsb__callocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp)
{
	return rsb__xallocate_coo_matrix_t(cmp,RSB_BOOL_TRUE,RSB_FLAG_NOFLAGS);
}

void * rsb__allocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp)
{
	return rsb__xallocate_coo_matrix_t(cmp,RSB_BOOL_FALSE,RSB_FLAG_NOFLAGS);
}

void * rsb__xallocate_coo_matrix_t(struct rsb_coo_mtx_t *cmp, rsb_bool_t want_calloc, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * Will allocate enough memory to fit either of COO or CSR.
	 * \return the input address on success, NULL on error
	 */
	size_t es, nnz, rnz;

	if(!cmp)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	rnz = nnz = cmp->nnz;
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS))
		rnz = RSB_MAX(nnz,cmp->nr+1);

        es = RSB_NUMERICAL_TYPE_SIZE(cmp->typecode);

	if(es<1)
	{
		RSB_ERROR("typecode seem wrong: 0%x\n",cmp->typecode);
		RSB_PERR_GOTO(err,RSB_ERRM_WTC);
	}
	cmp->IA = NULL;
	cmp->JA = NULL;
	cmp->VA = NULL;
	/* the above will avoid problems in case of error */

	if(want_calloc == RSB_BOOL_TRUE)
		cmp->IA = rsb__calloc(sizeof(rsb_coo_idx_t)*(rnz)),
		cmp->JA = rsb__calloc(sizeof(rsb_coo_idx_t)*(nnz)),
		cmp->VA = rsb__calloc_vector(nnz, cmp->typecode);
	else
		cmp->IA = rsb__malloc(sizeof(rsb_coo_idx_t)*(rnz)),
		cmp->JA = rsb__malloc(sizeof(rsb_coo_idx_t)*(nnz)),
		cmp->VA = rsb__malloc_vector(nnz, cmp->typecode);

	if(!cmp->IA || !cmp->JA || !cmp->VA)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	goto ret;
err:
	if(cmp)
		rsb__destroy_coo_matrix_t(cmp);
ret:
	return cmp;
}

void * rsb__destroy_coo_matrix_t(struct rsb_coo_mtx_t *cmp)
{
	/*!
	 * \ingroup gr_internals
	 * frees the given structure allocated arrays
	 *
	 * \return the input address
	 */
	if(!cmp)
		return cmp;
	RSB_CONDITIONAL_FREE(cmp->IA);
	RSB_CONDITIONAL_FREE(cmp->JA);
	RSB_CONDITIONAL_FREE(cmp->VA);
	/* RSB_BZERO_P(cmp); */
	/*
	 * Note: we do not free cmp itself.
	 */
	return cmp;
}

void * rsb__transpose_coo_matrix_t(struct rsb_coo_mtx_t *cmp)
{
	/*!
	 * \ingroup gr_internals
	 * transposes symbolically the given matrix
	 *
	 * \return the input address
	 */
	if(!cmp)
		return cmp;
	RSB_SWAP(rsb_coo_idx_t ,cmp->nr, cmp->nc );
	RSB_SWAP(rsb_coo_idx_t*,cmp->IA,cmp->JA);
	return cmp;
}

void * rsb__init_blank_pointers(struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return the input address
	 * */
	if(!mtxAp)
		return mtxAp;

	mtxAp->VA = NULL;
	mtxAp->indptr = NULL;
	mtxAp->bindx = NULL;
	mtxAp->rpntr = NULL;
	mtxAp->cpntr = NULL;
	mtxAp->bpntr = NULL;
#if RSB_WANT_BITMAP
	mtxAp->options = NULL;		/* exactly ... */
#endif /* RSB_WANT_BITMAP */
	mtxAp->mpntr = NULL;
	mtxAp->Mpntr = NULL;
	mtxAp->all_leaf_matrices = NULL;
	mtxAp->sm[0] = NULL;
	mtxAp->sm[1] = NULL;
	mtxAp->sm[2] = NULL;
	mtxAp->sm[3] = NULL;

	return mtxAp;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__fill_struct(struct rsb_mtx_t *mtxAp, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_type_t typecode, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * initializes a rsb_mtx_t struct to default 'vanilla' values
	 * \return the input address
	 * FIXME: rsb__fill_struct -> rsb__mtx_init
	 */

	if(!mtxAp)
		return RSB_ERR_GENERIC_ERROR;

	rsb__init_struct(mtxAp);/* redundant ?*/
	mtxAp->VA = VA;
	mtxAp->bindx = JA;
	mtxAp->bpntr = IA;
	mtxAp->nr = m;
	mtxAp->nc = k;
	mtxAp->flags = flags;
	mtxAp->typecode = typecode;
	
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void * rsb__fill_coo_struct(struct rsb_coo_mtx_t *mtxAp, void * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode)
{
	if(!mtxAp)
		return NULL;
	mtxAp->IA = IA;
	mtxAp->JA = JA;
	mtxAp->VA = VA;
	mtxAp->nnz = nnz;
	mtxAp->nr = m;
	mtxAp->nc = k;
	mtxAp->typecode = typecode;
	return mtxAp;
}

void * rsb__init_struct(struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * initializes a rsb_mtx_t struct to default 'vanilla' values
	 * \return the input address
	 * */
	if(!mtxAp)
		return mtxAp;

	rsb__init_blank_pointers(mtxAp);

	mtxAp->flags = RSB_FLAG_NOFLAGS ;

	mtxAp->sat = RSB_TIME_ZERO;
	mtxAp->eit = RSB_TIME_ZERO;
	mtxAp->pet = RSB_TIME_ZERO;
	mtxAp->est = RSB_TIME_ZERO;
	mtxAp->tat = RSB_TIME_ZERO;
	mtxAp->cpt = RSB_TIME_ZERO;
	mtxAp->rpt = RSB_TIME_ZERO;

	return mtxAp;
}

void * rsb__destroy_inner(struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \input mtxAp is a pointer to a valid matrix structure.
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * Deallocates the guts of a sparse matrix.
	 * (will leave the struct in an inconsistent state)
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;

	//RSB_STDOUT("destroying matrix %p\n",mtxAp);

	if(!mtxAp)
		goto ret;

	if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS))
	{
		if(rsb__is_root_matrix(mtxAp))
		{
			/* FIXME: unfinished, temporary */
			/* this is a trick: fitting the whole recursive matrix in three arrays */
			void *IA = NULL,*JA = NULL,*VA = NULL;
			const rsb_bool_t is_bio = rsb__do_is_matrix_binary_loaded(mtxAp); // binary I/O matrix
			struct rsb_mtx_t *fsm = rsb__do_get_first_submatrix(mtxAp);
			rsb_flags_t flags = mtxAp->flags;

			if(!is_bio)
				RSB_CONDITIONAL_FREE(mtxAp->all_leaf_matrices)	// ?!
			JA = fsm->bindx;
			VA = fsm->VA;
			if(!is_bio)
			{
				//IA = fsm->bpntr-(fsm->nr+1);
				IA = fsm->bpntr;
			}
			else
				IA = mtxAp;
			if(!is_bio)
				RSB_CONDITIONAL_FREE_MTXAP(mtxAp)// extra allocation
//			RSB_INFO("VA:%p, IA:%p, JA:%p\n",VA,IA,JA);
			if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
			{
				RSB_CONDITIONAL_FREE(IA);/* these arrays are allowed to be NULL, as it happens during conversions */
				RSB_CONDITIONAL_FREE(JA);
				RSB_CONDITIONAL_FREE(VA);
			}
		}
		return NULL;/* no deallocation, in this case */
	}

	if(rsb__is_recursive_matrix(mtxAp->flags))
	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	{
		if(submatrix)
		{
			rsb__do_mtx_free(submatrix);
		}
	}

	if(!(mtxAp->flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR))
	if(!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{
		RSB_CONDITIONAL_FREE(mtxAp->VA);
		RSB_CONDITIONAL_FREE(mtxAp->bindx);
		RSB_CONDITIONAL_FREE(mtxAp->bpntr);
	}

	RSB_CONDITIONAL_FREE(mtxAp->indptr);

#if RSB_WANT_BITMAP
	if(mtxAp->options)
		rsb__destroy_options_t(mtxAp->options);
#endif /* RSB_WANT_BITMAP */

	if((mtxAp->flags & RSB_FLAG_OWN_PARTITIONING_ARRAYS)!=0)
	{
		RSB_CONDITIONAL_FREE(mtxAp->rpntr);
		RSB_CONDITIONAL_FREE(mtxAp->cpntr);
	}
#if RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS
	RSB_CONDITIONAL_FREE(mtxAp->all_leaf_matrices);
#endif /* RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS */
	RSB_BZERO_P(mtxAp);/* this enforces correct usage */
ret:
	return NULL;
}

struct rsb_mtx_t * rsb__do_get_first_submatrix(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * */
	if(!mtxAp)
		return NULL;

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			return rsb__do_get_first_submatrix(submatrix);
	}
	return (struct rsb_mtx_t*)mtxAp;/* FIXME */
}

void * rsb__do_mtx_free(struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \param mtxAp a pointer to a matrix structure
	 *
	 * Will destroy a valid matrix and deallocate all of its allocated data.
	 * */
	rsb_flags_t flags;

	if(!mtxAp)
		goto ret;

	if( RSB_MTX_HBDF( mtxAp ) )
	{
		blas_sparse_matrix bmtxA = RSB_MTX_HBDFH(mtxAp);
		rsb__BLAS_Xusds(bmtxA);
		RSB_CONDITIONAL_FREE_MTXAP(mtxAp);
		goto ret;
	}

	flags = mtxAp->flags;

	rsb__destroy_inner(mtxAp);

	if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS))
	{
		if(mtxAp && RSB_WANT_ZERO_ON_DESTROY)
		{
			RSB_BZERO_P(mtxAp);
		}
		RSB_CONDITIONAL_FREE_MTXAP(mtxAp);
	}
	else
		mtxAp = NULL;

	RSB_DEBUG_ASSERT( !mtxAp );
ret:
	return mtxAp;
}


#if RSB_WANT_BITMAP
static size_t rsb__get_sizeof_options(const struct rsb_options_t *o, rsb_blk_idx_t M_b, rsb_blk_idx_t K_b)
{	
	/*!
	 * \ingroup gr_internals
	 * \return memory usage
	 * \param o a pointer to a valid rsb_options_t structure
	 *
	 * \return the amount of memory allocated for this structure, deeply
	 * */
	size_t count = 0;	

	if(!o )
		return 0;

	count += sizeof(*o);

	/* we allocate a new options structure */
	if(o->bitmap)
		count += RSB_BYTES_PER_BITMAP(M_b,K_b) ;

	return count;
}
#endif /* RSB_WANT_BITMAP */

size_t rsb__get_sizeof(const struct rsb_mtx_t *mtxAp )
{
	/*! 
	 * \ingroup gr_internals
	 * \param mtxAp a pointer to a valid rsb_mtx_t structure
	 * \return the amount of memory allocated for this structure, deeply (indices + coefficients).
	 * */
	size_t count = 0;	
	struct rsb_mtx_t * submatrix = NULL;
	rsb_submatrix_idx_t i,j;
	rsb_bool_t istrec = RSB_BOOL_FALSE;

	if(!mtxAp )
		goto err;

	istrec = rsb__is_terminal_recursive_matrix(mtxAp);

	count += sizeof(*mtxAp);

#if RSB_WANT_BITMAP
	if(mtxAp->options)
		count += rsb__get_sizeof_options(mtxAp->options,mtxAp->M_b,mtxAp->K_b);
	else
		return count;
#endif /* RSB_WANT_BITMAP */
	/* we allocate a new options structure */
	if( mtxAp->rpntr	) count += sizeof(rsb_coo_idx_t)*(mtxAp->M_b+1);
	if( mtxAp->cpntr	) count += sizeof(rsb_coo_idx_t)*(mtxAp->K_b+1);
	if(istrec)
	{
#if RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT
	if(rsb__is_coo_matrix(mtxAp))
	{
		if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			count += sizeof(rsb_half_idx_t)*(mtxAp->nnz)*2;
		else
			count += sizeof(rsb_coo_idx_t)*(mtxAp->nnz)*2;
	}
	else
	{
		if( mtxAp->flags & RSB_FLAG_USE_HALFWORD_INDICES)
			count += sizeof(rsb_half_idx_t)*(mtxAp->nnz)+sizeof(rsb_nnz_idx_t)*(mtxAp->nr+1);
		else
			count += sizeof(rsb_coo_idx_t)*(mtxAp->nnz)+sizeof(rsb_nnz_idx_t)*(mtxAp->nr+1);
	}
	if( mtxAp->VA  	) count += mtxAp->el_size*mtxAp->nnz;
	/* FIXME: missing the amount of memory allocated as extra submatrices for root, and the redundant structs array */
#else /* RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT */
	if( mtxAp->bindx	) count += sizeof(rsb_nnz_idx_t)*(mtxAp->block_count+1);
	if( mtxAp->indptr	) count += sizeof(rsb_nnz_idx_t)*(mtxAp->block_count+1);
	if( mtxAp->bpntr	) count += sizeof(rsb_nnz_idx_t)*(mtxAp->Mdim+1);
#if RSB_WANT_BITMAP
	if( mtxAp->VA  	) count += ( RSB_TOTAL_BLOCK_BYTES(mtxAp,mtxAp->options) );
#endif /* RSB_WANT_BITMAP */
#endif /* RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT */
	}
	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			count += rsb__get_sizeof(submatrix);

err:
	return count;
}

int rsb__nnz_coord_compar(const void *key, const void *am)
{	
	/*!
	 * \ingroup gr_internals
	 * A service function. NOTE:
	 *
	 * Please note that the sole use of this function is the major bottleneck during matrix creation.
	 * When thinking about optimizing matrix creation, come back here: this routine eats up to 90% 
	 * of the time required for matrix creation.
	 * */
	register const rsb_coo_idx_t *iam = am,*ik = key;
	/*!
	 * this function is used as compar in for stdlib's bsearch 
	 * on the ordered arrays p_r and p_c, and will return :
	 * -1, 0, or 1 
	 *  respectively if the key element is :
	 * less than both, in between, or greater than both *am and am[1].
	 * */
	return (*ik < iam[1] )? (*ik<*iam?-1:0):1;
}

rsb_bitmap_data_t * rsb__allocate_bitmap(rsb_blk_idx_t rows, rsb_blk_idx_t cols)
{
	/*!
	 * \ingroup gr_internals
	 * \param rows the amount of rows the bitmap should have
	 * \param rows the amount of columns the bitmap should have
	 * \return the bitmap area
	 *
	 * Allocates an area of (((cols+sizeof(rsb_bitmap_data_t)-1))/sizeof(rsb_bitmap_data_t) * rows) bytes
	 * to use as a bitmap, through the  RSB_BITMAP_SET 
	 * and  RSB_BITMAP_GET macros.
	 *
	 * This bitmap takes ( ceil(cols/sizeof(rsb_bitmap_data_t))*sizeof(rsb_bitmap_data_t)*rows ) bytes of memory.
	 * it should be roughly 1 bit for block of data
	 * or at worst
	 * ((cols/(sizeof(rsb_bitmap_data_t)*8))+1/(sizeof(rsb_bitmap_data_t)*8))/cols bits for every data block.
	 *
	 * The bitmap is return set to zero.
	 * assumes sizeof(rsb_bitmap_data_t)>1
	 * */
	if( RSB_INVALID_COO_INDEX(rows) || RSB_INVALID_COO_INDEX(cols))
		return NULL;
	return rsb__calloc(RSB_BYTES_PER_BITMAP(rows,cols));
}

rsb_bitmap_data_t * rsb__allocate_bitvector(rsb_blk_idx_t nbits)
{
	/*
	 * \ingroup gr_internals
	 * Allocates an array for \c nbits  bits, set to zero.
	 * */
#ifdef RSB_BITMAP_ROW_MAJOR_ORDER
	return rsb__allocate_bitmap(nbits,1);
#else /* RSB_BITMAP_ROW_MAJOR_ORDER */
	return rsb__allocate_bitmap(1,nbits);
#endif /* RSB_BITMAP_ROW_MAJOR_ORDER */
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_blk_idx_t rsb__bitmap_bit_count(const rsb_bitmap_data_t *bitmap, const rsb_blk_idx_t rows, const rsb_blk_idx_t cols)
{
	/*!
	 * \ingroup gr_internals
	 * \param bitmap is the bitmap data pointer
	 * \param rows is the number of rows the bitmap was allocated for
	 * \param cols is the number of columns the bitmap was allocated for
	 * \return -1 in case of error, the set bit count otherwise
	 *
	 * This function counts the bits in a bitmap
	 * Note that it is not dependent on the internal bitmap storage scheme,
	 * (column major or row major), but only as long as the bit pool is uniform.
	 *
	 * TODO : this is not core functionality, so it could be moved somewhere else.
	 * */
	register rsb_blk_idx_t i,bc = 0;
	register rsb_bitmap_data_t w;
	if(!bitmap)
		return RSB_ERR_BADARGS;
	for(i=0;i<((rsb_blk_idx_t)RSB_BYTES_PER_BITMAP(rows,cols)/sizeof(rsb_bitmap_data_t));++i)
	{
		w = bitmap[i];
		/* TODO : no stdlib functions for counting bits in integers ? */
		//b += (w&1); while(w) b += ((w /= 2)&1);
		while(w!=0) {bc += (w&1);w /= 2;}
	}
	/* warning : overflow check missing */
	return bc;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void* rsb__get_block_address( rsb_blk_idx_t blockrow, rsb_blk_idx_t blockcolumn, const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \param blockrow the block row
	 * \param blockcolumn the block column
	 * \param mtxAp a valid matrix structure pointer
	 * \return a pointer to the block itself or NULL if it is not present
 	 *
	 * A service function for getting the (blockrow,blockcolumn) block address inside the matrix.
	 *
	 * This function is SLOW, and should be used for debugging purposes only !
	 * ( it uses indirect indexing to catch elements )
	 * */
	rsb_nnz_idx_t l = 0;
	rsb_nnz_idx_t fnze = 0;
#if RSB_WANT_BITMAP
	struct rsb_options_t *o = NULL;
#endif /* RSB_WANT_BITMAP */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t offset = 0;

	if(!mtxAp)
	{errval = RSB_ERR_BADARGS;goto err;}

#if RSB_WANT_BITMAP
	o = mtxAp->options;
	if(!o)
	{errval = RSB_ERR_BADARGS;goto err;}
#endif /* RSB_WANT_BITMAP */

	if(RSB_BLK_ADD_OVERFLOW(blockrow,RSB_INDEX_OF_SAFE_EXTRA))
	{errval = RSB_ERR_LIMITS;goto err;}

	if(RSB_BLK_ADD_OVERFLOW(blockcolumn,RSB_INDEX_OF_SAFE_EXTRA))
	{errval = RSB_ERR_LIMITS;goto err;}

	/* i is the working block row */
	if( mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
	{
		if(mtxAp->bpntr[blockcolumn]==mtxAp->bpntr[blockcolumn+1])
			goto err;/* empty block column */
		fnze = mtxAp->bpntr[blockcolumn];	/* first nonzero entry in bindx */
		while(mtxAp->bindx[fnze+l]!=blockrow)++l;
	}
	else
	{
		if(mtxAp->bpntr[blockrow]==mtxAp->bpntr[blockrow+1])
			goto err;/* empty block row */
		fnze = mtxAp->bpntr[blockrow];	/* first nonzero entry in bindx */
		while(mtxAp->bindx[fnze+l]!=blockcolumn)++l;
	}

	if(RSB_NNZ_ADD_OVERFLOW(fnze,l))
	{errval = RSB_ERR_LIMITS;goto err;}

	offset = fnze+l;
	//return ((rsb_byte_t*)(mtxAp->VA)) + mtxAp->indptr[offset] * mtxAp->el_size;
	return RSB_BLOCK_ADDRESS(mtxAp,offset);
err:
	rsb__do_perror(NULL,errval);
	return NULL;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__recheck_insertion(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, const struct rsb_mtx_t *mtxAp, const struct rsb_options_t *o)
{
	/*!
	 * \ingroup gr_internals
	 * This is a very slow debug function.
	 * It should be supplied with some sparse matrix construction arrays in any order, and 
	 * a fully constructed matrix structure.
	 * 
	 * \note: Does not support block column major matrices and sorted ones.
	 *
	 * TODO  : should support more matrix formats ( e.g.: block column majer )
	 * FIXME : obsolete, very limited function
	 * 
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * */
	register rsb_blk_idx_t i,j;
	rsb_nnz_idx_t k;
	rsb_nnz_idx_t missing = 0;
	const rsb_byte_t *moff = NULL;
	const rsb_byte_t *src = NULL;

	if(!mtxAp || !o )return RSB_ERR_BADARGS;
	if(mtxAp->flags & RSB_FLAG_SORT_INPUT) return RSB_ERR_BADARGS;
	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER) return RSB_ERR_UNIMPLEMENTED_YET;
	if(! o->bitmap)  return RSB_ERR_UNSUPPORTED_OPERATION;/* when building sorted matrices, we don't create bitmaps .. should we ? */

#if RSB_WANT_BITMAP
	for(k=0;k<nnz;++k) { rsb_coo_idx_t iI = IA[k],iJ = JA[k];RSB_BLOCK_UNSET_BIT_FOR_NNZ(&iI,&iJ,k,mtxAp); }
	for(k=0;k<nnz;++k) { RSB_BLOCK_SET_BIT_FOR_NNZ(  IA,JA,k,mtxAp); }
#endif /* RSB_WANT_BITMAP */

	for(k=0;k<nnz;++k)
	{
		i = RSB_GET_BLOCK_ROW_FOR_NZ(IA+k,mtxAp);
		j = RSB_GET_BLOCK_COL_FOR_NZ(JA+k,mtxAp);
		if(!(RSB_BITMAP_GET(o->bitmap,mtxAp->M_b,mtxAp->K_b,i,j))) ++missing;
	}

	if(!missing)
		RSB_STDERR("checking structure : there are no blocks missing.\n");
	else
		RSB_STDERR("checking structure : there are %zd blocks missing\n",(size_t)missing);
	if(missing) return RSB_ERR_GENERIC_ERROR;
	for(k=0;k<nnz;++k)
	{
		i = RSB_GET_BLOCK_ROW_FOR_NZ(IA+k,mtxAp);
		j = RSB_GET_BLOCK_COL_FOR_NZ(JA+k,mtxAp);
		moff = rsb__get_block_address(i,j,mtxAp);
		if(!moff)
		{
			RSB_STDERR("critical block error on block (%d,%d).\n",i,j);
			return RSB_ERR_GENERIC_ERROR;
		}
		else
		{
			moff += RSB_GET_INTRA_BLOCK_OFFSET(IA[k],JA[k],i,j,mtxAp);
			src = VA;
			src += mtxAp->el_size*k;

			if(RSB_MEMCMP(src,moff,mtxAp->el_size))
			{
				/* may give problems when flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER */
				RSB_ERROR("critical error: %d'th nonzero (%d,%d) at (%d,%d) in block (%d,%d) is wrong!\n",
				k,
				IA[k]+1,JA[k]+1,
				RSB_INTRA_BLOCK_ROW(IA[k],i,mtxAp),RSB_INTRA_BLOCK_COLUMN(JA[k],j,mtxAp),i,j);
				/* warning : the following instruction is potentially harmful ! */
#ifdef RSB_DEBUG_BLOCK_STUFF
				RSB_STDERR("should be : 0x%x\n",*(int*)src );
				RSB_STDERR("is : 0x%x\n",*((int*)(moff)));
/*				RSB_STDERR("should be : %g\n",src );
				RSB_STDERR("is : %g\n",*((float*)(moff)));*/
#endif /* RSB_DEBUG_BLOCK_STUFF */
				return RSB_ERR_GENERIC_ERROR;
			}
		}
	}
	return RSB_ERR_NO_ERROR;	
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_is_valid_pinfo_t(const struct rsb_mtx_partitioning_info_t * pinfop)
{
	/*!
	 * \ingroup gr_internals
	 * \param pinfop should specify a partitioning info array
	 *
	 * This is a strictly debugging function, whose sole purpose is to verify
	 * the partitioning arrays contents of a rsb_mtx_partitioning_info_t structure.
	 * */
	rsb_nnz_idx_t k;
	if(pinfop->nr<1)
	{
		RSB_PERR_GOTO(err,"m == %d ?\n",pinfop->nr);
	}
	if(pinfop->nc<1)
	{
		RSB_PERR_GOTO(err,"k == %d ?\n",pinfop->nc);
	}
	if(pinfop->rpntr && pinfop->M_b<1)
	{
		RSB_PERR_GOTO(err,"M_b == %d ?\n",pinfop->M_b-1);
	}
	if(pinfop->cpntr && pinfop->K_b<1)
	{
		RSB_PERR_GOTO(err,"K_b == %d ?\n",pinfop->K_b-1);
	}

	if(pinfop->rpntr && pinfop->cpntr )
	{
		/* FIXME */
	if(pinfop->rpntr[pinfop->M_b]<=pinfop->rpntr[pinfop->M_b-1])
	{
		RSB_PERR_GOTO(err,"last (%d) rpntr element is %d <= %d\n",pinfop->M_b,pinfop->rpntr[pinfop->M_b],pinfop->rpntr[pinfop->M_b-1]);
	}
	if(pinfop->cpntr[pinfop->K_b]<=pinfop->cpntr[pinfop->K_b-1])
	{
		RSB_PERR_GOTO(err,"last (%d) cpntr element is %d <= %d\n",pinfop->K_b,pinfop->cpntr[pinfop->K_b],pinfop->cpntr[pinfop->K_b-1]);
	}

	for(k=0;k<pinfop->M_b;++k)if(pinfop->rpntr[k]<0)
	{
		RSB_PERR_GOTO(err,"bad partitioning : rpntr[%d]=%d\n",k,pinfop->rpntr[k]);
	}
	for(k=0;k<pinfop->K_b;++k)if(pinfop->cpntr[k]<0)
	{
		RSB_PERR_GOTO(err,"bad partitioning : cpntr[%d]=%d\n",k,pinfop->cpntr[k]);
	}
	for(k=0;k<pinfop->M_b;++k)if(pinfop->rpntr[k]>pinfop->nr)
	{
		RSB_PERR_GOTO(err,"bad partitioning : rpntr[%d]=%d > m==%d\n",k,pinfop->rpntr[k],pinfop->nr);
	}
	for(k=0;k<pinfop->K_b;++k)if(pinfop->cpntr[k]>pinfop->nc)
	{
		RSB_PERR_GOTO(err,"bad partitioning : cpntr[%d]=%d > k==%d\n",k,pinfop->cpntr[k],pinfop->nc);
	}
	}
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_GENERIC_ERROR;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if 0
rsb_err_t rsb__compute_partial_fillin_for_nnz_fraction(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,const  rsb_nnz_idx_t nnz, struct rsb_mtx_partitioning_info_t * pinfop, size_t * element_countp, size_t * block_countp)
{
	/*!
	 * \ingroup gr_internals
	 * see rsb__compute_partial_fillin_for_nnz_fractions
	 * */
	return rsb__compute_partial_fillin_for_nnz_fractions(IA,JA,&nnz,1,pinfop,element_countp,block_countp);
}
#endif

rsb_err_t rsb__compute_partial_fillin_for_nnz_fractions(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,const  rsb_nnz_idx_t * nnz, const rsb_nnz_idx_t nnzn, struct rsb_mtx_partitioning_info_t * pinfop, size_t * element_countp, size_t * block_countp)
{
	/*!
	 * \ingroup gr_internals
	 * \param IA is a row indices array sized nnz
	 * \param JA is a column indices array sized nnz
	 * \param nnz is the length of IA and JA
	 * \param element_countp is where the element counts will be written
	 * \param block_countp   is where the block   counts will be written
	 *
	 * Will estimate fillin for the first nnz[ni]<nnz[ni+1]<...<nnz[nnzn-1] elements, with nnz[nnzn-1] being 
	 * less than or equal to the number of elements in the IA, JA, element_countp, block_countp arrays.
	 *
	 * Note: this function performs almost no data validation.
	 * Note : this is not a service but an experimental function, and is very slow.
	 * Note : another wayof structuring thisfunction would beto make it accept a 2*nnzn sized array with lower
	 *        and upper segment indices both specified.
	 *        This would have been more flexible but would require some change in this function code.
	 * TODO : this is not core functionality, so this function could be moved elsewhere
	 * */

	rsb_bitmap_data_t * bitmap = NULL;
	size_t  element_count = 0;
	rsb_nnz_idx_t block_count = 0;
	rsb_nnz_idx_t k = 0,l = 0;/* were -1 */
	rsb_blk_idx_t i = 0,j = 0;/* were -1 */

	if( !IA || !JA || !pinfop || !element_countp || !block_countp ) goto err;
	if( nnzn < 1 ) goto err;
	if( RSB_INVALID_BLK_INDEX(pinfop->M_b) || RSB_INVALID_BLK_INDEX(pinfop->K_b) )goto err;
	if( !pinfop->rpntr  || !pinfop->cpntr  )goto err;
	bitmap = rsb__allocate_bitmap(pinfop->M_b,pinfop->K_b);
	if(!bitmap)goto err;

	#ifdef  RSB_DEBUG_INPUT
	if(rsb__do_is_valid_pinfo_t(pinfop))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	#endif /* RSB_DEBUG_INPUT */
	
	if(1)
	{
	size_t skip = 0,skipu = 0;//new
	rsb_nnz_idx_t c = 0;
	rsb_nnz_idx_t sl;

	skip = 0;
	
	skipu = (nnz[nnzn-1]/(nnzn*nnzn));
	skip = nnzn*skipu;
	for(sl=0;sl<nnzn;++sl)
		block_countp[sl] = element_countp[sl] = 0;

	#if 1
	/* An alternative way, much more stable!
	 * However, it underestimates the nnz count so it should in some manner mitigated! */
	for(l=0;l<nnzn;++l)
	{
		for(sl=0;sl<nnzn;++sl)
		{
		//RSB_INFO("#i: %d / %d  (%d)\n",sl*skip+l*skipu,sl*skip+(l+1)*skipu,nnzn);
		for(k=sl*skip+l*skipu;k<sl*skip+(l+1)*skipu;++k)
		{
	
			++c;
			/* NOTE : the following ifelse statements are for situations where m<br or k<bc  */
			if(pinfop->M_b>1)
				i = RSB_GET_BLOCK_ROW_FOR_NZ(IA+k,pinfop);
			else
			       	i = 0;
			if(pinfop->K_b>1)
				j = RSB_GET_BLOCK_COL_FOR_NZ(JA+k,pinfop);
			else
			       	j = 0;

			/* if the bit is not already set */
			if(!(RSB_BITMAP_GET(bitmap,pinfop->M_b,pinfop->K_b,i,j)))
			{
				element_count += GET_BLOCK_SIZE(i,j,pinfop);
				(block_count)++;
				RSB_BITMAP_SET(bitmap,pinfop->M_b,pinfop->K_b,i,j) ;
			}
		}
		block_countp[l] = block_count;
		element_countp[l] = element_count;
		}
	}
	l = nnzn-1;
	//RSB_INFO("#c: %d / %d (%d)..  %d -> %d\n",c,block_countp[l],nnzn,nnzn*skip,nnz[l]);
	for(k=c;k<nnz[l];++k)
	//for(k=nnzn*skip;k<nnz[l];++k)
	{
//		++c;
	
			i = RSB_GET_BLOCK_ROW_FOR_NZ(IA+k,pinfop);
			j = RSB_GET_BLOCK_COL_FOR_NZ(JA+k,pinfop);

			/* if the bit is not already set */
			if(!(RSB_BITMAP_GET(bitmap,pinfop->M_b,pinfop->K_b,i,j)))
			{
				element_countp[l] += GET_BLOCK_SIZE(i,j,pinfop);
				block_countp[l]++;
				RSB_BITMAP_SET(bitmap,pinfop->M_b,pinfop->K_b,i,j) ;
			}
	}
	//RSB_INFO("#c: %d / %d (%d)\n",c,block_countp[l],nnzn);
	#endif
	}
	else
	for(l=0;l<nnzn;++l)
	{
		rsb_nnz_idx_t li;
		if(l==0)
			li = 0;
	       	else
		       	li = nnz[l-1];/* will the first loop optimized by the compiler ? :) */

		for(k=li;k<nnz[l];++k)
		{
	
			i = RSB_GET_BLOCK_ROW_FOR_NZ(IA+k,pinfop);
			j = RSB_GET_BLOCK_COL_FOR_NZ(JA+k,pinfop);

			/* if the bit is not already set */
			if(!(RSB_BITMAP_GET(bitmap,pinfop->M_b,pinfop->K_b,i,j)))
			{
				element_count += GET_BLOCK_SIZE(i,j,pinfop);
				(block_count)++;
				RSB_BITMAP_SET(bitmap,pinfop->M_b,pinfop->K_b,i,j) ;
			}
		}
		if(block_countp)
			block_countp[l] = block_count;
		if(element_countp)
			element_countp[l] = element_count;
	}

	RSB_CONDITIONAL_FREE(bitmap);
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_ENOMEM;
}

#if RSB_WANT_BITMAP
static rsb_err_t rsb_element_block_count_and_bitmap_from_coo_partitioning(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop, size_t * element_countp, rsb_nnz_idx_t * block_countp, rsb_bitmap_data_t ** bitmapp, const rsb_flags_t flags, struct rsb_mtx_t * mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \param IA is a row indices array sized nnz
	 * \param JA is a column indices array sized nnz
	 * \param nnz is the length of IA and JA
	 * \param pinfop should specify a partitioning info array
	 * \param element_countp is where the element counts will be written
	 * \param block_countp   is where the block   counts will be written
	 *
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * WARNING : IA and JA can respectively point to columns and rows arrays instead of rows and columns,
	 * as long as the pinfop information is swapped accordingly.
	 * In this way a transposed bitmap will be allocated (note that its size could be the same. guess why..).
	 *
	 * In case of error, no bitmap will be allocated, but its pointer may be overwritten.
	 * In case of success, a bitmap structure will be allocated.
	 * */

	rsb_nnz_idx_t k = 0;
	rsb_bitmap_data_t * bitmap = NULL;
	size_t element_count = 0;
	rsb_nnz_idx_t block_count = 0;
	rsb_blk_idx_t mI = 0,MI = 0;
	const rsb_coo_idx_t * mIndx = NULL,* MIndx = NULL;

	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	{
		mIndx = IA;
		MIndx = JA;
	}
	else
	{
		mIndx = JA;
		MIndx = IA;
	}

	if( !IA || !JA || !pinfop || !element_countp || !bitmapp || !block_countp ) goto err;
	if( RSB_INVALID_NNZ_INDEX(nnz) ) goto err;
	if( RSB_INVALID_BLK_INDEX(pinfop->M_b) || RSB_INVALID_BLK_INDEX(pinfop->K_b) )goto err;
	if( !pinfop->rpntr  || !pinfop->cpntr  )goto err;

	bitmap = rsb__allocate_bitmap(mtxAp->Mdim,mtxAp->mdim);

	if(!bitmap)goto err;

	if( mtxAp->flags & RSB_FLAG_SHOULD_DEBUG )
		if(rsb__do_is_valid_pinfo_t(pinfop))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
	
	if(RSB_WANT_VERBOSE_MESSAGES)
		RSB_INFO("counting matrix blocks ..\n");

	for(k=0;RSB_LIKELY(k<nnz);++k)
	{
		/* 
		 * We count the amount of elements for each block, setting bits in
		 * our bitmap where a block should be placed, and leaving unset bits
		 * which correspond to zero blocks 
		 * */

		MI = RSB_GET_BLOCK_MAJ_FOR_NZ(MIndx+k,mtxAp);
		mI = RSB_GET_BLOCK_MIN_FOR_NZ(mIndx+k,mtxAp);

		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(MI));
		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(mI));

		if(mI>=mtxAp->mdim)
		{
			RSB_PERR_GOTO(err," j=%d >= o->K_b=%d\n ",mI,mtxAp->mdim);
		} 
		if(mI <0 )
		{
			RSB_PERR_GOTO(err," j=%d < 0 \n",mI);
		}
		if(MI>=mtxAp->Mdim)
		{
			RSB_PERR_GOTO(err," i=%d >= o->M_b=%d\n ",MI,mtxAp->Mdim);
		}
		if(MI <0 )
		{
			RSB_PERR_GOTO(err," i=%d < 0 \n",MI);
		}
		
		/* if the bit is not already set */
		if(!(RSB_BITMAP_GET(bitmap,mtxAp->Mdim,mtxAp->mdim,MI,mI)))
		{
			if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
				element_count += GET_BLOCK_SIZE(mI,MI,pinfop);
			else
				element_count += GET_BLOCK_SIZE(MI,mI,pinfop);

			(block_count)++;
			RSB_BITMAP_SET(bitmap,mtxAp->Mdim,mtxAp->mdim,MI,mI) ;
		}
	}

	if(block_count > nnz)
	{
		RSB_PERR_GOTO(err,"(mtxAp->block_count=%d >= n=%d)!\n",block_count,nnz);
	}
	
	*block_countp = block_count;
	*element_countp = element_count;
	*bitmapp = bitmap;

	return RSB_ERR_NO_ERROR;
err:
	RSB_CONDITIONAL_FREE(bitmap);
	return RSB_ERR_GENERIC_ERROR;
}
#endif /* RSB_WANT_BITMAP */

rsb_err_t rsb__do_set_init_storage_flags(struct rsb_mtx_t *mtxAp, rsb_flags_t flags)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_flags_t storage_only_flags = RSB_DO_FLAGS_EXTRACT_STORAGE(flags);

#ifdef RSB_FLAG_WANT_LINKED_STORAGE
	if(RSB_DO_FLAG_HAS(storage_only_flags,RSB_FLAG_WANT_LINKED_STORAGE))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
#ifdef RSB_MATRIX_STORAGE_LC
			mtxAp->matrix_storage = RSB_MATRIX_STORAGE_LC;
#els /* RSB_MATRIX_STORAGE_LC */e
			{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_LC */
		else
#ifdef RSB_MATRIX_STORAGE_LR
			mtxAp->matrix_storage = RSB_MATRIX_STORAGE_LR;
#else /* RSB_MATRIX_STORAGE_LR */
			{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif  /* RSB_MATRIX_STORAGE_LR */
	}
	else
#endif /* RSB_FLAG_WANT_LINKED_STORAGE */
	{
		if(RSB_DO_FLAG_HAS(storage_only_flags,RSB_FLAG_WANT_COO_STORAGE))
		{
			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
#ifdef RSB_MATRIX_STORAGE_BCOC
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCOC;
#else /* RSB_MATRIX_STORAGE_BCOC */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_BCOC */
			else
#ifdef RSB_MATRIX_STORAGE_BCOR
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCOR;
#else /* RSB_MATRIX_STORAGE_BCOR */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_BCOR */
		}
		else //FIXME: this switch could cohexist with CSS, and be processed later or be ignored (in the old constructor)
		if(RSB_DO_FLAG_HAS(storage_only_flags,RSB_FLAG_WANT_COO_STORAGE))
		{
			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
#ifdef RSB_MATRIX_STORAGE_BCOC
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCOC;
#else /* RSB_MATRIX_STORAGE_BCOC */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_BCOC */
			else
#ifdef RSB_MATRIX_STORAGE_BCOR
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCOR;
#else /* RSB_MATRIX_STORAGE_BCOR */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_BCOR */
		}
		else
		if(RSB_DO_FLAG_HAS(storage_only_flags,RSB_FLAG_WANT_BCSS_STORAGE))
		{
			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
#ifdef RSB_MATRIX_STORAGE_BCSC
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCSC;
#else /* RSB_MATRIX_STORAGE_BCSC */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_BCSC */
			else
#ifdef RSB_MATRIX_STORAGE_BCSR
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_BCSR;
#else /* RSB_MATRIX_STORAGE_BCSR */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_BCSR */
		}
		else
		if(RSB_DO_FLAG_HAS(storage_only_flags,RSB_FLAG_WANT_FIXED_BLOCKING_VBR))
		{
			if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
#ifdef RSB_MATRIX_STORAGE_VBC
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_VBC;
#else /* RSB_MATRIX_STORAGE_VBC */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_VBC */
			else
#ifdef RSB_MATRIX_STORAGE_VBR
				mtxAp->matrix_storage = RSB_MATRIX_STORAGE_VBR;
#else /* RSB_MATRIX_STORAGE_VBR */
				{errval = RSB_ERR_UNSUPPORTED_FORMAT;goto err;}
#endif /* RSB_MATRIX_STORAGE_VBR */
		}
		else
		{
			/* undetermined format or a merge of formats (happens on recursive matrices during construction) */
			mtxAp->matrix_storage = storage_only_flags;
		}
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__set_init_flags_and_stuff( struct rsb_mtx_t *mtxAp, struct rsb_options_t * o, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_nnz_idx_t block_count, rsb_nnz_idx_t element_count, rsb_type_t typecode, rsb_flags_t flags )
{
	/**
	 * \ingroup gr_internals
	 *
	 * This inner service function sets some flags and variables during matrix construction.
	 *
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if !RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS 
	if( /* !o */ /* FIXME: o disabled lately || */
	 !pinfop 
	)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS  */

	mtxAp->flags = flags;
	mtxAp->typecode =typecode;
	mtxAp->el_size = RSB_NUMERICAL_TYPE_SIZE(typecode);
#if RSB_WANT_BITMAP
	mtxAp->options = NULL;	// we ignore o
//	mtxAp->options = o;
#endif /* RSB_WANT_BITMAP */
	mtxAp->nnz = nnz;
	mtxAp->element_count = element_count;
#if RSB_WANT_DBC
	mtxAp->block_count = block_count;
#endif

	if(pinfop)
	{
		mtxAp->M_b = pinfop->M_b;
		mtxAp->K_b = pinfop->K_b;
		mtxAp->nr = pinfop->nr;
		mtxAp->nc = pinfop->nc;
		mtxAp->rpntr = pinfop->rpntr;
		mtxAp->cpntr = pinfop->cpntr;
//#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
		mtxAp->br = pinfop->br;
		mtxAp->bc = pinfop->bc;
//#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR */
	}
	else
	{
		mtxAp->M_b = m;
		mtxAp->K_b = k;
		mtxAp->nr = m;
		mtxAp->nc = k;
		mtxAp->rpntr = NULL;
		mtxAp->cpntr = NULL;
//#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
		mtxAp->br = 1;
		mtxAp->bc = 1;
//#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR */
		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(mtxAp->br+1));
		RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(mtxAp->bc+1));
	}
	if(RSB_IS_INVALID_TYPE_SIZE(mtxAp->el_size = rsb__do_sizeof(mtxAp->typecode)))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if((errval = rsb__do_set_init_storage_flags(mtxAp,flags))!=RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(mtxAp->br==1 && mtxAp->bc==1)
	{
		mtxAp->M_b = mtxAp->nr;
		mtxAp->K_b = mtxAp->nc;
	}

	/* setting aliases */
	if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	{
		mtxAp->Mdim = mtxAp->K_b;
		mtxAp->mdim = mtxAp->M_b;
		mtxAp->mpntr = mtxAp->rpntr;
		mtxAp->Mpntr = mtxAp->cpntr;
	}
	else
	{
		mtxAp->Mdim = mtxAp->M_b;
		mtxAp->mdim = mtxAp->K_b;
		mtxAp->Mpntr = mtxAp->rpntr;
		mtxAp->mpntr = mtxAp->cpntr;
	}

//	RSB_DEBUG_ASSERT(mtxAp->Mdim);
//	RSB_DEBUG_ASSERT(mtxAp->mdim);
//	RSB_DEBUG_ASSERT(mtxAp->rpntr);
//	RSB_DEBUG_ASSERT(mtxAp->cpntr);
	RSB_DEBUG_ASSERT(mtxAp->el_size);

	return RSB_ERR_NO_ERROR;
err:
	RSB_DO_ERR_RETURN(errval)
}

#if 0
rsb_err_t rsb_dump_matrix ( const struct rsb_mtx_t *mtxAp )
{
	/*!
	 * \ingroup gr_internals
	 *
	 *  \param mtxAp is a valid matrix structure pointer
	 *  \param diagonal is an array sized as min(mtxAp->nr,mtxAp->nc) which on exit will contain the diagonal elements.
	 *  \return -1 in case of error, 0 otherwise
	 *
	 * FIXME : UNTESTED AND UNDOCUMENTED AND UNFINISHED
	 * FIXME : USE rsb_print_matrix and delete this ?
	 * */
	register rsb_nnz_idx_t baserow,basecolumn;
	register rsb_blk_idx_t rows,columns;
	register rsb_blk_idx_t blockrow,blockcolumn;
	register rsb_byte_t *bp;

	RSB_INFO("%% [!] TESTING CODE !\n");
	RSB_INFO("%%rows:%d columns:%d blocks:%d\n",mtxAp->nr,mtxAp->nc,mtxAp->block_count);
	RSB_INFO("%d %d %d\n", mtxAp->nr,mtxAp->nc,mtxAp->nnz);
	RSB_GET_FIRST_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
	while(!RSB_GOT_LAST_BLOCK_POINTER(mtxAp))
	{
		rsb_coo_idx_t r,c;
		/*
		 * FIXME
		 * */
//		RSB_INFO("%x \n", bp) ;
//		RSB_INFO("_k : %d %d   ", _k,_lastk) ;
//		RSB_INFO("%d %d ", baserow,basecolumn) ;
		RSB_INFO("%d %d\n", rows,columns) ;
#if 1
		for(r=0;r<rows;++r)
		for(c=0;c<columns;++c)
		{
				RSB_INFO("%d %d %lg\n", baserow+r,basecolumn+c,((double*)bp)[columns*r+c]) ;
		}
#endif
		RSB_GET_NEXT_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
	}

	return RSB_ERR_NO_ERROR;
}
#endif

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_insert_sorted( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop)
{
	/*!
	 * \ingroup gr_internals
	 * Inserts in the matrix structures the specified coo elements, sorted accordingly to the specified BCSR blocking.
	 */

	rsb_coo_idx_t blockrows = 0;
	rsb_coo_idx_t blockcolumns = 0;
	rsb_coo_idx_t baserow = 0;
	rsb_coo_idx_t basecolumn = 0;
	rsb_nnz_idx_t *indptr = mtxAp->indptr;
	const rsb_coo_idx_t *Mpntr = NULL;
	const rsb_coo_idx_t *mpntr = NULL;
	const rsb_coo_idx_t *MIndx = NULL;
	const rsb_coo_idx_t *mIndx = NULL;
	rsb_blk_idx_t mI = 0,MI = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;	/* will index a nnz sized array */
	rsb_nnz_idx_t K = 0;
	rsb_byte_t*dst = NULL;
	rsb_byte_t*src = NULL;

	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	{
		mpntr = pinfop->rpntr;
		Mpntr = pinfop->cpntr;
		mIndx = IA;
		MIndx = JA;
	}
	else
	{
		Mpntr = pinfop->rpntr;
		mpntr = pinfop->cpntr;
		MIndx = IA;
		mIndx = JA;
	}

	k = mI = MI = 0;K = 0;
	blockrows = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
	while( MIndx[k] >= Mpntr[MI+1] )++MI;	/* skipping 'preceding' block rows .. */
	while( mIndx[k] >= mpntr[mI+1] )++mI;	/* skipping 'preceding' block columns .. */
	baserow = Mpntr[MI];
	basecolumn = mpntr[mI];
	mtxAp->bindx [ K ] = mI;			/* a 'new' block */
	indptr[ K+1 ]=indptr[ K  ] + blockrows * blockcolumns;

#ifdef RSB_FLAG_WANT_LINKED_STORAGE
	if( rsb__have_linked_storage(mtxAp->flags) )
	{
		if(RSB_WANT_VERBOSE_MESSAGES)
			RSB_INFO("initializing linked lists stuff.\n");
		if(RSB_UNLIKELY(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),mI,MI,blockcolumns,blockrows,basecolumn,baserow)
		else
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),MI,mI,blockrows,blockcolumns,baserow,basecolumn)
	}
#endif /* RSB_FLAG_WANT_LINKED_STORAGE */

/*
	dst = mtxAp->VA;
	dst += RSB_BLOCK_OFFSET(mtxAp,K);
	{rsb_blk_idx_t ibo = 0;
	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
		ibo = RSB_GET_INTRA_BLOCK_OFFSET(mIndx[k],MIndx[k],mI,MI,mtxAp) ;
	else
		ibo = RSB_GET_INTRA_BLOCK_OFFSET(MIndx[k],mIndx[k],MI,mI,mtxAp) ;
	dst += ibo;}
	src = ((rsb_byte_t*)VA) + mtxAp->el_size * k;
	RSB_NUMERICAL_TYPE_SET_ELEMENT(dst,src,mtxAp->typecode);*/

	while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%d : (%d %d) is not ok\n",k, MIndx[k]+1,mIndx[k]+1);
			RSB_STDERR("(minor dim. index %d < base row %d)\n",MIndx[k] , baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif /* DEBUG */

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
			while( mIndx[k] >= mpntr[mI+1] )++mI;
			blockcolumns = mpntr[mI+1] - mpntr[mI];
			basecolumn = mpntr[mI];

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];
			}
			else
			{
				/* same block row  */
			}
			++K;
			mtxAp->bindx [ K ] = mI;			/* a 'new' block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
#ifdef RSB_FLAG_WANT_LINKED_STORAGE
			if( rsb__have_linked_storage(mtxAp->flags) )
			{
				if(RSB_UNLIKELY(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),mI,MI,blockcolumns,blockrows,basecolumn,baserow)
				else
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),MI,mI,blockrows,blockcolumns,baserow,basecolumn)
			}
#endif /* RSB_FLAG_WANT_LINKED_STORAGE */
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
			while( MIndx[k] >= Mpntr[MI+1] )++MI;
			blockrows = Mpntr[MI+1] - Mpntr[MI];
			baserow = Mpntr[MI];

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;
				while( mIndx[k] >= mpntr[mI+1] )++mI;
				blockcolumns = mpntr[mI+1] - mpntr[mI];
				basecolumn = mpntr[mI];
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			mtxAp->bindx [ K ] = mI;			/* a 'new' block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
#ifdef RSB_FLAG_WANT_LINKED_STORAGE
			if( rsb__have_linked_storage(mtxAp->flags) )
			{
				if(RSB_UNLIKELY(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),mI,MI,blockcolumns,blockrows,basecolumn,baserow)
				else
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),MI,mI,blockrows,blockcolumns,baserow,basecolumn)
			}
#endif /* RSB_FLAG_WANT_LINKED_STORAGE */
		}
		else
		{
			/* same block row for sure */
		}
		dst = mtxAp->VA;
		dst += RSB_BLOCK_OFFSET(mtxAp,K);
		{
		rsb_nnz_idx_t ibo = 0;
		if(RSB_UNLIKELY(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
			ibo = RSB_GET_INTRA_BLOCK_OFFSET(mIndx[k],MIndx[k],mI,MI,mtxAp) ;
		else
			ibo = RSB_GET_INTRA_BLOCK_OFFSET(MIndx[k],mIndx[k],MI,mI,mtxAp) ;
		dst += ibo;
		}
		//RSB_ERROR("%d %d %d\n",((rsb_byte_t*)dst)-((rsb_byte_t*)mtxAp->VA),MIndx[k],mIndx[k]);
		src = ((rsb_byte_t*)VA) + mtxAp->el_size * k;
		RSB_NUMERICAL_TYPE_SET_ELEMENT(dst,src,mtxAp->typecode);
		++k;
	}

	if(nnz)++K;
	mtxAp->bindx[K] = 0;	// the first element off the 'working' bindx should be set to a safe value

	if(mtxAp->flags & RSB_FLAG_SHOULD_DEBUG)
	if( K != mtxAp->block_count )
	{
		RSB_ERROR("K is %zd ! should be %zd (block count)!\n",(size_t)K,(size_t)mtxAp->block_count);
		RSB_STDERR("nnz : %zd\n",(size_t)nnz);
		RSB_STDERR("k : %zd\n",(size_t)k);
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_account_sorted( struct rsb_mtx_t * mtxAp, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_nnz_idx_t * elements_per_block_row, rsb_nnz_idx_t * blocks_per_block_row)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * The routine for inserting BCSR sorted coo elements in a fresh matrix.
	 * It is not optimized, and it should be used as a debug resort when tuning optimized ones.
	 *
	 * FIXME : this code is deprecated in favour of rsb__do_account_sorted_optimized
	 * FIXME : does not support lots of flags!
	 */
	rsb_coo_idx_t blockrows = 0;
	rsb_coo_idx_t blockcolumns = 0;
	rsb_coo_idx_t baserow = 0;
	rsb_coo_idx_t basecolumn = 0;
	const rsb_coo_idx_t *Mpntr = NULL;
	const rsb_coo_idx_t *mpntr = NULL;
	const rsb_coo_idx_t *MIndx = NULL;
	const rsb_coo_idx_t *mIndx = NULL;
	rsb_blk_idx_t mI = 0,MI = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;	/* will index a nnz sized array */
	rsb_nnz_idx_t K = 0;
	k = mI = MI = K=0;

	if( ! IA || ! JA 
#if !RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS 
	|| !pinfop
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
	)
	{
		errval = RSB_ERR_BADARGS;goto err;
	}

	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	{
		mpntr = pinfop->rpntr;
		Mpntr = pinfop->cpntr;
		mIndx = IA;
		MIndx = JA;
	}
	else
	{
		Mpntr = pinfop->rpntr;
		mpntr = pinfop->cpntr;
		MIndx = IA;
		mIndx = JA;
	}

	while( MIndx[k] >= Mpntr[MI+1] )++MI;	/* skipping 'preceding' block rows .. */
	while( mIndx[k] >= mpntr[mI+1] )++mI;	/* skipping 'preceding' block columns .. */
	blockrows = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
	baserow = Mpntr[MI];
	basecolumn = mpntr[mI];
	elements_per_block_row[MI*0] += blockrows * blockcolumns;
	blocks_per_block_row[MI] += 1;

	while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%d : (%d %d) is not ok\n",k, MIndx[k]+1,mIndx[k]+1);
			RSB_STDERR("(minor dim. index %d < base row %d)\n",MIndx[k] , baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;
		}
#endif /* DEBUG */

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
			while( mIndx[k] >= mpntr[mI+1] )++mI;
			blockcolumns = mpntr[mI+1] - mpntr[mI];
			basecolumn = mpntr[mI];

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
				while( MIndx[k] >= Mpntr[MI+1] )++MI;
				blockrows = Mpntr[MI+1] - Mpntr[MI];
				baserow = Mpntr[MI];
			}
			else
			{
				/* same block row  */
			}
			elements_per_block_row[MI*0] += blockrows * blockcolumns;
			blocks_per_block_row[MI] += 1;
			++K;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
			while( MIndx[k] >= Mpntr[MI+1] )++MI;
			blockrows = Mpntr[MI+1] - Mpntr[MI];
			baserow = Mpntr[MI];

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;
				while( mIndx[k] >= mpntr[mI+1] )++mI;
				blockcolumns = mpntr[mI+1] - mpntr[mI];
				basecolumn = mpntr[mI];
			}
			else
			{
				/* new row block, same column  */
			}
			/* get rid of this var : elements_per_block_row */
			elements_per_block_row[MI*0] += blockrows * blockcolumns;
			blocks_per_block_row[MI] += 1;
			++K;
		}
		else
		{
			/* same block row for sure */
		}
		++k;
	}
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
struct rsb_mtx_t * rsb__allocate_css_from_coo_sorted( void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t m, rsb_coo_idx_t k, struct rsb_options_t * o, rsb_type_t typecode, rsb_flags_t flags, rsb_err_t *errvalp)
{
	/*!
	 * \ingroup gr_internals
	 * The routine for matrix building from sorted coo elements. CSR only.
	 *
	 * FIXME : EXPERIMENTAL, UNFINISHED
	 * FIXME : FIX THIS FUNCTION TO ALLOCATE CSR/CSC EVEN IF NOT IN PLACE
	 * */
//	rsb_nnz_idx_t n = 0;
	rsb_nnz_idx_t * elements_per_block_row = NULL;
	rsb_time_t t = RSB_TIME_ZERO;
	struct rsb_mtx_t *mtxAp = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* const rsb_coo_idx_t *MIndx = NULL, *mIndx = NULL; */
	rsb_blk_idx_t MI = 0;

	if(!errvalp)
		return NULL;

	if(!( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR ))
	{
		return NULL;
	}
	if(flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
		return NULL;/* FIXME : only csr now */

	pinfop = NULL;/* FIXME */

	if(!o) {errval = RSB_ERR_BADARGS;goto err;}
	rsb__init_struct(mtxAp = rsb__calloc(sizeof(*mtxAp)));
	if(!mtxAp){errval = RSB_ERR_ENOMEM;goto err;}
	if((errval = rsb__set_init_flags_and_stuff(mtxAp,o,pinfop,m,k,0,0,0,typecode,flags))!=RSB_ERR_NO_ERROR)goto err;
	/*
	if(flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	       	mIndx = IA, MIndx = JA;
       	else
	       	MIndx = IA, mIndx = JA;
	*/
	elements_per_block_row = rsb__calloc(sizeof(rsb_nnz_idx_t)*(1+mtxAp->Mdim));

	if(!elements_per_block_row)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	t = - rsb_time();

	RSB_DEBUG_ASSERT(rsb__util_is_sorted_coo_as_row_major(IA,JA,nnz,typecode,pinfop,flags)==RSB_ERR_NO_ERROR);
	RSB_DEBUG_ASSERT(rsb__util_are_valid_coo_arrays(IA,JA,nnz)==RSB_ERR_NO_ERROR);

	errval = rsb__do_account_sorted_optimized(mtxAp,IA,JA,m,k,nnz,NULL,elements_per_block_row,NULL);
	mtxAp->block_count = 0;

	mtxAp->block_count = nnz;

	t += rsb_time();
	mtxAp->sat = t;

	mtxAp->indptr = rsb__malloc(sizeof(rsb_nnz_idx_t)*(mtxAp->block_count+1));
	mtxAp->bindx = JA;	/* ok, done :) FIXME : should type - convert ... */

	if(!mtxAp->bindx ){ errval = RSB_ERR_ENOMEM; goto err;}
	if(!mtxAp->indptr){ errval = RSB_ERR_ENOMEM; goto err;}

	mtxAp->indptr[0] = 0;/* */

	mtxAp->bpntr = IA;
	mtxAp->bpntr [0] = 0;
//	mtxAp->bpntr = NULL;
	for(MI=0;MI<mtxAp->Mdim;++MI) mtxAp->bpntr[MI+1]= mtxAp->bpntr[MI]+ elements_per_block_row[MI];

#if RSB_WANT_BITMAP
	mtxAp->options = o ;
#endif /* RSB_WANT_BITMAP */
	mtxAp->nnz = nnz;
	mtxAp->element_count = nnz;
	mtxAp->VA = VA;
	t = - rsb_time();
	errval = rsb__do_insert_sorted_optimized(mtxAp,VA,IA,JA,nnz,NULL);
	if(RSB_SOME_ERROR(errval))
		goto err;
	t += rsb_time();
	mtxAp->eit = t;
	RSB_CONDITIONAL_FREE(elements_per_block_row);
	return mtxAp;
err:
	RSB_STDERR("rsb__allocate_from_coo_sorted:\n");
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	if(mtxAp)
		rsb__do_mtx_free(mtxAp);	/* destroys all of the internals of matrix */
	RSB_CONDITIONAL_FREE(elements_per_block_row);
	return NULL;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
#if RSB_WANT_DBC
struct rsb_mtx_t * rsb__allocate_from_coo_sorted( const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t m, rsb_coo_idx_t k, struct rsb_options_t * o, rsb_type_t typecode, rsb_flags_t flags, rsb_err_t *errvalp)
{
	/*!
	 * \ingroup gr_internals
	 * The routine for matrix building from sorted coo elements.
	 * This function requires the coefficients to be sorted accordingly to the inter block ordering policy.
	 *
	 * \param pinfop is the pointer to a rsb_mtx_partitioning_info_t structure with partitioning information.
	 *
	 * \note : should behave well with the flags: RSB_FLAG_WANT_FIXED_BLOCKING_VBR RSB_FLAG_SORTED_INPUT
	 * \note : this function should be optimized and tested thoroughly.
	 * */

	struct rsb_mtx_t *mtxAp = NULL;
	rsb_blk_idx_t MI = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* rsb_coo_idx_t blockcolumns = 0; */

	const rsb_coo_idx_t *MIndx = NULL;

	rsb_nnz_idx_t * elements_per_block_row = NULL;
	rsb_nnz_idx_t * blocks_per_block_row = NULL;	/* per major dimension .. */
	size_t element_count = 0;

	rsb_time_t t = RSB_TIME_ZERO;

#if !RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS 
	if(!pinfop)
	{
		errval = RSB_ERR_BADARGS;goto err;
	}
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
	if( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR )
	{errval = RSB_ERR_BADARGS;goto err;}

//	if(!o)
//	{
//		errval = RSB_ERR_BADARGS;goto err;
//	}

	blocks_per_block_row = NULL;	/* per major dimension .. */
	element_count = 0;

	rsb__init_struct(mtxAp = rsb__calloc(sizeof(struct rsb_mtx_t)));
	if(!mtxAp){errval = RSB_ERR_ENOMEM;goto err;}

	if((errval = rsb__set_init_flags_and_stuff(mtxAp,o,pinfop,m,k,0,0,0,typecode,flags))!=RSB_ERR_NO_ERROR)
		goto err;
	if(flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
		MIndx = JA;
	else
		MIndx = IA;

	/* FIXME : elements_per_block_row can be replaced with a single variable */
	elements_per_block_row = rsb__calloc(sizeof(rsb_nnz_idx_t)*(1+mtxAp->Mdim));
	blocks_per_block_row = rsb__calloc(sizeof(rsb_nnz_idx_t)*(1+mtxAp->Mdim));

	if(!blocks_per_block_row  )
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	if(!elements_per_block_row)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	t = - rsb_time();
	blocks_per_block_row++;/* we increment the pointer for 1 element (we will use this array as bpntr later)*/
	errval = rsb__do_account_sorted_optimized(mtxAp,IA,JA,m,k,nnz,pinfop,elements_per_block_row,blocks_per_block_row);
	if(RSB_SOME_ERROR(errval))
		goto err;

	if(nnz==0)blocks_per_block_row[0] = 0;/* handling the degenerate nnz == 0 case (e.g.: unit diag) */
	mtxAp->block_count = 0;
	element_count = 0;

	for(MI=0;MI<mtxAp->Mdim;++MI)mtxAp->block_count += blocks_per_block_row  [MI];
	for(MI=0;MI<mtxAp->Mdim;++MI)element_count += elements_per_block_row[MI];

	t += rsb_time();
	mtxAp->sat = t;
	if(RSB_WANT_VERBOSE_MESSAGES)
	RSB_STDERR("matrix creation phase 1 (accounting) : %lf seconds \n", t);

	mtxAp->indptr = rsb__malloc(sizeof(rsb_nnz_idx_t)*(mtxAp->block_count+1));
	mtxAp->bindx = rsb__malloc(sizeof(rsb_nnz_idx_t)*(mtxAp->block_count+1));

	if(!mtxAp->bindx ){ errval = RSB_ERR_ENOMEM; goto err;}
	if(!mtxAp->indptr){ errval = RSB_ERR_ENOMEM; goto err;}

	mtxAp->indptr[0] = 0;/* */

	mtxAp->bpntr = (--blocks_per_block_row);	/* :) */
	for(MI=0;MI<mtxAp->Mdim;++MI)
		mtxAp->bpntr[MI+1] += mtxAp->bpntr[MI];	/* in this way bpntr[i] has the count of blocks before row i */
	mtxAp->bpntr [0] = 0;
	blocks_per_block_row = NULL;		/* it will be freed with the matrix */

	/* second pass : we have allocated the needed arrays and are ready to fill in data structures */

#if RSB_WANT_BITMAP
	mtxAp->options = o ;
#endif /* RSB_WANT_BITMAP */
	mtxAp->nnz = nnz;
	mtxAp->element_count = element_count;
	//mtxAp->block_count = block_count;

	if(mtxAp->block_count > nnz)
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_STDERR("more blocks (%zd) than nonzeros (%zd) ?could be a bug!\n",(size_t)mtxAp->block_count,(size_t)nnz);
		goto err;
	}

	mtxAp->bpntr[0] = 0;
	mtxAp->VA = rsb__malloc( RSB_TOTAL_BLOCK_BYTES(mtxAp,o));

	if(RSB_WANT_VERBOSE_MESSAGES)
		RSB_INFO("allocating %zd bytes.\n",(size_t)RSB_TOTAL_BLOCK_BYTES(mtxAp,o) );
		
	if(!mtxAp->VA)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_STDERR("had problems allocating %zd bytes.\n",(size_t)RSB_TOTAL_BLOCK_BYTES(mtxAp,o));
		goto err;
	}
	//	k = 0;/* nnz index */
	
	t = - rsb_time();

	/* the following code could run parallel with some work */
	errval = rsb__do_insert_sorted_optimized( mtxAp, VA, IA, JA, nnz, pinfop);

	if(RSB_SOME_ERROR(errval))
		goto err;
	t += rsb_time();
	mtxAp->eit = t;
	if(RSB_WANT_VERBOSE_MESSAGES)
		RSB_STDERR("matrix creation phase 2 (insertion) : %lf seconds \n", mtxAp->eit);
#if 0
	if((flags & RSB_FLAG_SHOULD_DEBUG) && 0)
	{
		register rsb_coo_idx_t	baserow,basecolumn,rows,columns;
		register rsb_blk_idx_t	blockrow,blockcolumn;
		register rsb_byte_t*bp;

		/* FIXME : will fail if pure bcsr */
		RSB_GET_FIRST_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
		if(0 /* super paranoia */)
		while(!RSB_GOT_LAST_BLOCK_POINTER(mtxAp))
		{
			RSB_INFO("%zd / %zd  ; block (%zd %zd)/(%zd %zd) base : (%zd %zd) size : (%zd %zd)\n",
			(rsb_printf_int_t)_lastk,(rsb_printf_int_t)mtxAp->block_count,
			(rsb_printf_int_t)blockrow,(rsb_printf_int_t)blockcolumns,
			(rsb_printf_int_t)pinfop->M_b,(rsb_printf_int_t)pinfop->K_b,
			(rsb_printf_int_t)baserow,(rsb_printf_int_t)basecolumn,(rsb_printf_int_t)rows,(rsb_printf_int_t)columns);
			RSB_GET_NEXT_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
		}
	}
#endif

	RSB_CONDITIONAL_FREE(elements_per_block_row);
	return mtxAp;
err:
	RSB_STDERR("rsb__allocate_from_coo_sorted:\n");
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	if(mtxAp)
		rsb__do_mtx_free(mtxAp);	/* destroys all of the internals of matrix */
	if(blocks_per_block_row != MIndx )
		RSB_CONDITIONAL_FREE(blocks_per_block_row );
	RSB_CONDITIONAL_FREE(elements_per_block_row);
	return NULL;
}
#endif /* RSB_WANT_DBC */
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if RSB__USE_MTX_PARTITIONING_INFO_T
rsb_err_t rsb__do_get_blocking_from_pinfo(const struct rsb_mtx_partitioning_info_t * pinfop, rsb_flags_t flags, rsb_blk_idx_t *mbp, rsb_blk_idx_t *kbp)
{
	if( ( flags & RSB_FLAG_WANT_BCSS_STORAGE ) || ( flags & RSB_FLAG_WANT_FIXED_BLOCKING_VBR ) )
	{
		if( pinfop && pinfop->cpntr && pinfop->rpntr )
		{
			/* FIXME : experimental */
			*kbp = pinfop->cpntr[1]-pinfop->cpntr[0];
			*mbp = pinfop->rpntr[1]-pinfop->rpntr[0];
		}
		else
		if(pinfop)
		{
#if RSB_EXPERIMENTAL_USE_PURE_BCSS
			*mbp = pinfop->br;
			*kbp = pinfop->bc;
#else /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
			*kbp = *mbp = -1;
#endif /* RSB_EXPERIMENTAL_USE_PURE_BCSS */
		}
		else
		{
			*kbp = *mbp = 1;
		}
		RSB_DEBUG_ASSERT(*kbp>=1);
		RSB_DEBUG_ASSERT(*mbp>=1);
		if( *kbp<1 || *mbp <1 )
		{
			return RSB_ERR_BADARGS;
		}
	}
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB__USE_MTX_PARTITIONING_INFO_T */

size_t rsb__util_strlen(const rsb_char_t *s)
{
	/*!
	 * \ingroup gr_internals
	 */
	return strlen(s);/* Flawfinder: ignore */
}

#if 0
static int rsb_util_sprintf(rsb_char_t *str, const rsb_char_t *format, ...)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME :  BUGGY
	 * */
        va_list ap;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	va_start(ap,format);
	errval = rsb__sprintf(str,format,ap);/* Flawfinder: ignore */
	va_end(ap);
	RSB_DO_ERR_RETURN(errval)
}
#endif

rsb_char_t *rsb__util_strcat(rsb_char_t *dest, const rsb_char_t *src)
{
	/*!
	 * A wrapper.
	 * \ingroup gr_internals
	 */
	return strcat(dest,src); /* Flawfinder: ignore */
}

const rsb_char_t * rsb__sprint_matrix_implementation_code2(const struct rsb_mtx_t *mtxAp, rsb_char_t * buf, rsb_flags_t inflags)
{
	/*!
	 * \ingroup gr_internals
	 *  FIXME : missing error handling 
	 * buf be at least RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH chars long.
	 */
	const rsb_char_t sep[] = "\t";
	rsb_blk_idx_t br,bc;
	if(!mtxAp) return NULL;
	buf[0] = '\0';

	rsb__get_blocking_size(mtxAp,&br,&bc);

	/* NOTE : assumes BCSR or takes into account only the first blocks */
	rsb__sprintf(buf+rsb__util_strlen(buf),"%ld%s%ld%s",(long)mtxAp->nr,sep,(long)mtxAp->nc,sep);
	rsb__sprintf(buf+rsb__util_strlen(buf),"%ld%s%ld%s",(long)br,sep,(long)bc,sep);
	rsb__sprintf(buf+rsb__util_strlen(buf),"%zd%s%lg",(size_t)rsb__do_get_matrix_nnz(mtxAp),sep,rsb__do_get_matrix_fillin(mtxAp));

	return buf;
}

rsb_char_t rsb__do_get_symmetry_char(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 */
	if(rsb__is_symmetric(mtxAp))
		return 'S';
	else
	if(rsb__is_hermitian(mtxAp))
		return 'H';
	else
		return 'G';
}

static const rsb_char_t * rsb_do_get_symmetry_string(const struct rsb_mtx_t *mtxAp, rsb_char_t * auxbuf)
{
	/*!
	 * \ingroup gr_internals
	 */
	const rsb_char_t * const s = "Symmetric";
	const rsb_char_t * const g = "General";
	const rsb_char_t * const h = "Hermitian";

	if(rsb__is_symmetric(mtxAp))
		rsb__strcpy(auxbuf,s);
	else
	if(rsb__is_hermitian(mtxAp))
		rsb__strcpy(auxbuf,h);
	else
		rsb__strcpy(auxbuf,g);
	return auxbuf;
}

static const rsb_char_t * rsb__sprint_matrix_implementation_code(const struct rsb_mtx_t *mtxAp, const rsb_char_t * op, rsb_flags_t inflags, rsb_char_t * buf);

rsb_err_t rsb__fprint_matrix_implementation_code(const struct rsb_mtx_t *mtxAp, const rsb_char_t * op, rsb_flags_t inflags, FILE*fd)
{
	const rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_char_t buf[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];/* Flawfinder: ignore */
	
	fprintf( fd, "%s", rsb__sprint_matrix_implementation_code(mtxAp,op,inflags,buf));
	return errval;
}

void rsb__cat_compver(rsb_char_t * buf)
{
	/* FIXME : fix	rsb_util_sprintf and use it ! */
#if defined(__INTEL_COMPILER)
	/* icc 10.10 is ok */
	rsb__sprintf(buf,"intel-%d",__INTEL_COMPILER);
#elif   defined(__xlC__)
	/* ok on sp5 */
	rsb__sprintf(buf,"xlc-%d",__xlC__);
#elif   defined(__PGI)
	/* pgcc-7.0.4 is ok */
	rsb__sprintf(buf,"pgcc-%d.%d.%d",__PGIC__,__PGIC_MINOR__,__PGIC_PATCHLEVEL__);
#elif   defined(__clang__)
	rsb__sprintf(buf,"clang-%d.%d.%d",__clang_major__,__clang_minor__,__clang_patchlevel__);
	/* see also __clang_version__ */
#elif   defined(__FUJITSU)
	rsb__sprintf(buf,"fujitsu-%d.%d.%d"/*"(%s)"*/,__FCC_major__,__FCC_minor__,__FCC_patchlevel__/*,__FCC_version__*/);
#elif   defined(__GNUC__)
	rsb__sprintf(buf,"gcc-%d.%d",__GNUC__,__GNUC_MINOR__);
#elif defined(__SUNPRO_CC)
	rsb__sprintf(buf,"sun-%d",__SUNPRO_CC);
#else /* __SUNPRO_CC */
	rsb__util_strcat(buf,"CC?");
#endif /* __SUNPRO_CC */
}

static const rsb_char_t * rsb__sprint_matrix_implementation_code(const struct rsb_mtx_t *mtxAp, const rsb_char_t * op, rsb_flags_t inflags, rsb_char_t * buf)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Gives back a matrix implementation info string.
	 * NOTE: for more consistency, we should translate any \t (TAB) char in ' '.
	 * NOTE: it will give some more info, too..
	 *
	 * \return a static string pointer on correct operation, NULL otherwise
	 * */
	
	rsb_char_t sep[] = "/";
	const rsb_char_t * csp;
	rsb_long_t sm = 0;
	rsb_long_t tsm = 0;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_blk_idx_t br = 0,bc = 0;
	rsb_char_t auxbuf[RSB_TOKLEN];

	if(!mtxAp)
		return NULL;

	rsb__get_blocking_size(mtxAp,&br,&bc);

	flags = mtxAp->flags|inflags;

	sm = rsb__submatrices(mtxAp);
	tsm = rsb__terminal_recursive_matrix_count(mtxAp);

	buf[0] = '\0';

	/* FIXME : DANGER */

	if(1)
	{
		long hcoo = rsb__terminal_recursive_matrix_count_with_storage_and_flags(mtxAp,RSB_MATRIX_STORAGE_BCOR,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
		long hcsr = rsb__terminal_recursive_matrix_count_with_storage_and_flags(mtxAp,RSB_MATRIX_STORAGE_BCSR,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
		long fcoo = rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(mtxAp,RSB_MATRIX_STORAGE_BCOR,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
		long fcsr = rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(mtxAp,RSB_MATRIX_STORAGE_BCSR,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
		long kinds = 0;
		if(hcoo)++kinds;
		if(hcsr)++kinds;
		if(fcoo)++kinds;
		if(fcsr)++kinds;
#if 0
		if(fcoo==0 && hcoo==0)
			rsb__util_strcat(buf,"CSR");
		else
		if(fcsr==0 && hcsr==0)
			rsb__util_strcat(buf,"COO");
		else
#endif
			rsb__util_strcat(buf,"RSB");
	}
	else
	{
	if(rsb__is_recursive_matrix(flags))
		rsb__util_strcat(buf,"R");

#ifdef RSB_MATRIX_STORAGE_BCOR
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCOR)
	{
		if(br==1&&bc==1)
			rsb__util_strcat(buf,"COR");
		else
			rsb__util_strcat(buf,RSB_MATRIX_STORAGE_BCOR_STRING);
	}
	else
#endif /* RSB_MATRIX_STORAGE_BCOR */
#ifdef RSB_MATRIX_STORAGE_BCOC
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCOC)
	{
		if(br==1&&bc==1)
			rsb__util_strcat(buf,"COC");
		else
			rsb__util_strcat(buf,RSB_MATRIX_STORAGE_BCOC_STRING);
	}
	else
#endif /* RSB_MATRIX_STORAGE_BCOC */
#ifdef RSB_MATRIX_STORAGE_BCSR
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSR)
	{
		if(br==1&&bc==1)
			rsb__util_strcat(buf,"CSR");
		else
			rsb__util_strcat(buf,RSB_MATRIX_STORAGE_BCSR_STRING);
	}
	else
#endif /* RSB_MATRIX_STORAGE_BCSR */
#ifdef RSB_MATRIX_STORAGE_BCSC
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_BCSC)
	{
		if(br==1&&bc==1)
			rsb__util_strcat(buf,"CSC");
		else
			rsb__util_strcat(buf,RSB_MATRIX_STORAGE_BCSC_STRING);
	}
	else
#endif /* RSB_MATRIX_STORAGE_BCSC */
#ifdef RSB_MATRIX_STORAGE_VBR 
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_VBR )
		rsb__util_strcat(buf,RSB_MATRIX_STORAGE_VBR_STRING);
	else
#endif /* RSB_MATRIX_STORAGE_VBR */
#ifdef RSB_MATRIX_STORAGE_VBC
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_VBC )
		rsb__util_strcat(buf,RSB_MATRIX_STORAGE_VBC_STRING);
	else
#endif /* RSB_MATRIX_STORAGE_VBC */
#ifdef RSB_MATRIX_STORAGE_LR
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_LR )
		rsb__util_strcat(buf,"LBLR");
	else
#endif /* RSB_MATRIX_STORAGE_LR */
#ifdef RSB_MATRIX_STORAGE_LC
	if(mtxAp->matrix_storage & RSB_MATRIX_STORAGE_LC )
		rsb__util_strcat(buf,"LBLC");
	else
#endif /* RSB_MATRIX_STORAGE_LC */
		return NULL;
	}
	{
//	if(sm>=1 && rsb__is_recursive_matrix(mtxAp))/* NEW */ /* FIXME : rsb__is_recursive_matrix() seems plagued by indeterminism! */
		if(sm>=1 /*&& rsb__is_recursive_matrix(mtxAp)*/)/* NEW */
			rsb__sprintf(buf+rsb__util_strlen(buf),"(@:%ld/%ld;%3.1lf%%diagnz;%3.1lf%%diagblk)",sm,tsm,
					(((double)rsb__get_diagonal_elements_count(mtxAp)*100)/(mtxAp->nnz)),
					(((double)rsb__get_diagonal_submatrices_count(mtxAp)*100)/(tsm))
					);
	}
#if 1
	/* uhm. this refers to inter block ordering. */
	if( mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
		rsb__util_strcat(buf,"-C");
	else
		rsb__util_strcat(buf,"-R");
#endif
	rsb__util_strcat(buf,sep);
	rsb__util_strcat(buf,"RowMajor");
	rsb__util_strcat(buf,sep);
	rsb__util_strcat(buf,rsb_do_get_symmetry_string(mtxAp,auxbuf));
	rsb__util_strcat(buf,sep);
	rsb__util_strcat(buf,op?op:"");
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_AUTO_BLOCKING))
		rsb__util_strcat(buf,"-AutoBlocking");
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR))
		rsb__util_strcat(buf,"-InPlace");
#if 0
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES_COO))
		rsb__sprintf(buf+rsb__util_strlen(buf),"-SwitchToHalfwordCoo:(%ld~%ld)"
		,rsb__terminal_recursive_matrix_count_with_flags(mtxAp,RSB_FLAG_USE_HALFWORD_INDICES_COO)
		,tsm);
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES_CSR))
		rsb__sprintf(buf+rsb__util_strlen(buf),"-SwitchToHalfwordCsr:(%ld~%ld)"
		,rsb__terminal_recursive_matrix_count_with_flags_but(mtxAp,RSB_FLAG_USE_HALFWORD_INDICES_CSR,RSB_FLAG_USE_HALFWORD_INDICES_COO)
		,tsm);
#else
	{
		long hcoo = rsb__terminal_recursive_matrix_count_with_storage_and_flags(mtxAp,RSB_MATRIX_STORAGE_BCOR,RSB_FLAG_USE_HALFWORD_INDICES_COO);
		long hcsr = rsb__terminal_recursive_matrix_count_with_storage_and_flags(mtxAp,RSB_MATRIX_STORAGE_BCSR,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
		long fcoo = rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(mtxAp,RSB_MATRIX_STORAGE_BCOR,RSB_FLAG_USE_HALFWORD_INDICES_COO);
		long fcsr = rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(mtxAp,RSB_MATRIX_STORAGE_BCSR,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
		rsb__sprintf(buf+rsb__util_strlen(buf),"-HalfwordCsr:(%ld~%ld)",hcsr,tsm);
		rsb__sprintf(buf+rsb__util_strlen(buf),"-FullwordCsr:(%ld~%ld)",fcsr,tsm);
		rsb__sprintf(buf+rsb__util_strlen(buf),"-HalfwordCoo:(%ld~%ld)",hcoo,tsm);
		rsb__sprintf(buf+rsb__util_strlen(buf),"-FullwordCoo:(%ld~%ld)",fcoo,tsm);
	}
#endif
#if 0
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_RECURSIVE_HALF_DETECTED_CACHE))
		rsb__util_strcat(buf,"-BlockForHalfCache");
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_RECURSIVE_DOUBLE_DETECTED_CACHE))
		rsb__util_strcat(buf,"-BlockForDoubleCache");
#endif
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG))
		rsb__util_strcat(buf,"-ExtraDiagonalSubdivisions");
#ifdef RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES))
		rsb__util_strcat(buf,"-NoMicroLeafs");
#endif /* RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES */

	rsb__util_strcat(buf,sep);
	RSB_NUMERICAL_TYPE_STRING(csp,mtxAp->typecode);
	rsb__util_strcat(buf,csp);

	if(1)
	{
		rsb_thread_t ncores = 0;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
                #pragma omp parallel RSB_NTC
                if(omp_get_thread_num()==0)
                {
                        ncores = omp_get_num_threads();
                }
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		ncores = ncores?ncores:1;
		rsb__util_strcat(buf,sep);
		rsb__sprintf(buf+rsb__util_strlen(buf),"cores:%d",ncores);
	}
	/* see http://gasnet.cs.berkeley.edu/dist/other/portable_platform.h */
	rsb__util_strcat(buf,sep);
	/* NOTE : on some systems, __GNUC__ is defined even under icc! (therefore we switched precedence) */

#if 0
	/* FIXME : fix	rsb_util_sprintf and use it ! */
#if defined(__INTEL_COMPILER)
	/* icc 10.10 is ok */
	rsb__sprintf(buf+rsb__util_strlen(buf),"intel-%d",__INTEL_COMPILER);
#elif   defined(__xlC__)
	/* ok on sp5 */
	rsb__sprintf(buf+rsb__util_strlen(buf),"xlc-%d",__xlC__);
#elif   defined(__PGI)
	/* pgcc-7.0.4 is ok */
	rsb__sprintf(buf+rsb__util_strlen(buf),"pgcc-%d.%d.%d",__PGIC__,__PGIC_MINOR__,__PGIC_PATCHLEVEL__);
#elif   defined(__GNUC__)
	rsb__sprintf(buf+rsb__util_strlen(buf),"gcc-%d.%d",__GNUC__,__GNUC_MINOR__);
#elif defined(__SUNPRO_CC)
	rsb__sprintf(buf+rsb__util_strlen(buf),"sun-%d",__SUNPRO_CC);
#else /* __SUNPRO_CC */
	rsb__util_strcat(buf,"CC?");
#endif /* __SUNPRO_CC */
#else
	rsb__cat_compver(buf+rsb__util_strlen(buf));
#endif
	rsb__util_strcat(buf,sep);
	/* still missing CXX case */
#if   defined(CFLAGS)
	//rsb_util_sprintf(buf+rsb__util_strlen(buf),"%s",CFLAGS);
	rsb__sprintf(buf+rsb__util_strlen(buf),"%s",CFLAGS);
#else /* CFLAGS */
	rsb__util_strcat(buf,"");
#endif /* CFLAGS */
	/* NEW */
	rsb__util_strcat(buf,sep);
	rsb__sprintf(buf+rsb__util_strlen(buf),"sizeof(nnz_idx_t):%zd,",sizeof(rsb_nnz_idx_t));
	rsb__sprintf(buf+rsb__util_strlen(buf),"sizeof(coo_idx_t):%zd,",sizeof(rsb_coo_idx_t));
	rsb__sprintf(buf+rsb__util_strlen(buf),"sizeof(blk_idx_t):%zd",sizeof(rsb_blk_idx_t));

	/* NEW */
	rsb__util_strcat(buf,sep);
	rsb__sprintf(buf+rsb__util_strlen(buf),"idx_storage:%zd-idx_storage_in_csr:%zd-idx_storage_in_coo:%zd"
		,(size_t)rsb__get_index_storage_amount(mtxAp)
		,((size_t)mtxAp->nnz)*sizeof(rsb_coo_idx_t)+((size_t)mtxAp->Mdim+1)*sizeof(rsb_nnz_idx_t)
		,((size_t)mtxAp->nnz)*sizeof(rsb_coo_idx_t)*2
		);

	rsb__util_strcat(buf,sep);
#ifdef RSB_PACKAGE_VERSION 
	rsb__sprintf(buf+rsb__util_strlen(buf),"version:%s",RSB_PACKAGE_VERSION);
#endif /* RSB_PACKAGE_VERSION */
	rsb__util_strcat(buf,sep);
	rsb__util_strcat(buf,"memhinfo:[");
	{rsb_char_t usmhib[RSB_MAX_LINE_LENGTH];
	rsb__util_strcat(buf,rsb__get_mem_hierarchy_info_string(usmhib));}
	rsb__util_strcat(buf,"]");

	if ( rsb__getenv_int_t("RSB_USE_HOSTNAME", 1) )
	{
		rsb__util_strcat(buf,sep);
		rsb__strcpy_hostname(buf);
	}
	return buf;
}

rsb_err_t rsb__util_get_tn_array(const rsb_char_t* optarg, int* bxlp, rsb_blk_idx_t **bxvp)
{
	// fill bxvp to 1,2,..2^K,maxt
	int bxl = 0;
	rsb_blk_idx_t * bxv = NULL;
	int mint = 1,maxt = 1;

	if(*optarg==':')
	{
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		maxt = rsb_global_session_handle.rsb_want_threads;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		while(mint<=maxt)
			mint *= 2, ++bxl;
		if( mint > maxt && mint/2 != maxt )
			bxl++;
		mint = 1;
	}
	else
		goto err;

	bxv = rsb__malloc(sizeof(rsb_blk_idx_t)*(size_t)bxl);
	if(!bxv)
		goto err;
	bxl = 0;

	if(*optarg==':')
	{
		while( mint <= maxt )
			bxv[bxl++] = mint, mint *= 2;
		if( bxv[bxl-1] != maxt )
			bxv[bxl++] = maxt;
	}
	
	if(*bxvp)
		rsb__free(*bxvp);
	*bxvp = bxv;
	*bxlp = bxl;
	
	return RSB_ERR_NO_ERROR;
err:
	RSB_CONDITIONAL_FREE(bxv);
	rsb__do_perror(NULL,RSB_ERR_GENERIC_ERROR);
	return RSB_ERR_GENERIC_ERROR;
}

rsb_err_t rsb__util_get_bx_array_or_default(const rsb_char_t defsym, const rsb_char_t* defarg, const rsb_char_t* optarg, int* bxlp, rsb_blk_idx_t **bxvp)
{
	if(!optarg || *optarg==defsym)
		return rsb__util_get_bx_array(defarg,bxlp,bxvp);
	else
		return rsb__util_get_bx_array(optarg,bxlp,bxvp);
}

rsb_err_t rsb__util_get_bx_array(const rsb_char_t* optarg, int* bxlp, rsb_blk_idx_t **bxvp)
{
	/*!
	   	\ingroup gr_internals
	  
		Extract block row and block column sizes from user data codified in optarg.
		It can accept zeroes as well.
		\param bxlp will be set to the number of the desired block sizes.
		\param bxvp will be set to an array (of dimension *bxlp) allocated with rsb__malloc() with block sizes.
	 */
	int bxl = 0;
	rsb_blk_idx_t * bxv = NULL;
	const rsb_char_t*c = optarg;

	if(!bxlp || !bxvp)
		return RSB_ERR_BADARGS;

	do
	{
		int nul = 0;
		while(*c!=nul && !isdigit(*c))
			++c;
		if(isdigit(*c))
			bxl++;
		while( isdigit(*c))
			++c;
	} while(*c);

	bxv = rsb__malloc(sizeof(rsb_blk_idx_t)*(size_t)bxl);
	if(!bxv)
		goto err;
	bxl = 0;
	c = optarg;

	do
	{
		int nul = 0,ci;

		ci = rsb__util_atoi_km10(c);
		if(ci<0 && c[0]!='-')
			goto err;
		if(*c == '-')
			++c;
		if(isdigit(*c))
			bxv[bxl++] = (rsb_blk_idx_t)ci;
		while( isdigit(*c))
			++c;
		while(*c!=nul && !(isdigit(*c) || *c == '-'))
			++c; /* e.g. 'K' in "1K" or '-' in "1,-1K" */
	} while(*c);
	
	if(*bxvp)
		rsb__free(*bxvp);
	*bxvp = bxv;
	*bxlp = bxl;
	
	return RSB_ERR_NO_ERROR;
err:
	RSB_CONDITIONAL_FREE(bxv);
	rsb__do_perror(NULL,RSB_ERR_GENERIC_ERROR);
	return RSB_ERR_GENERIC_ERROR;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_long_t rsb__util_atol(const rsb_char_t *nptr)
{
	/*!
	  	\ingroup gr_internals
	 */
	return atol(nptr);/* Flawfinder: ignore */
}

rsb_nnz_idx_t rsb__util_atonnz(const rsb_char_t * optarg)
{
	/*!
		\ingroup gr_internals

	   	Will parse a single rsb_nnz_idx_t number.
	 	\param optarg

		\note : may overflow
		\warning : may overflow
	 */
	rsb_long_t i = rsb__util_atol(optarg);/*Flawfinder: ignore */
	if(i<0)
		i = 0;
	return (rsb_nnz_idx_t)i;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_real_t rsb__util_atof(const rsb_char_t *nptr)
{
	/*!
	  	\ingroup gr_internals
	 */
	return atof(nptr);/* Flawfinder: ignore */
}

int rsb__util_atoi(const rsb_char_t *nptr)
{
	/*!
	  	\ingroup gr_internals
	 */
	int n = 0;

	if(nptr)
		n = atoi(nptr);/* Flawfinder: ignore */

	return n;
}

static int rsb__util_atoi_kmX(const rsb_char_t *nptr, int base)
{
	/*!
	  	\ingroup gr_internals
	 */
	int v = rsb__util_atoi(nptr);

	if(!nptr || !base)
	{
		RSB_ERROR(RSB_ERRM_BADARGS);
		goto ret;
	}

	if(*nptr=='-')
		++nptr;
	while(isdigit(*nptr))
		++nptr;
	if(*nptr && tolower(*nptr)=='g')
		v *= base * base * base;
	if(*nptr && tolower(*nptr)=='m')
		v *= base * base;
	if(*nptr && tolower(*nptr)=='k')
		v *= base;
ret:
	return v;
}

int rsb__util_atoi_km2(const rsb_char_t *nptr)
{
	return rsb__util_atoi_kmX(nptr, 1024);
}

int rsb__util_atoi_km10(const rsb_char_t *nptr)
{
	// TODO: check rsb__util_atonnz
	return rsb__util_atoi_kmX(nptr, 1000);
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__copy_css_arrays(const void *iVA, const rsb_coo_idx_t * iINDX, const rsb_coo_idx_t * iXA, const rsb_nnz_idx_t nnz, rsb_coo_idx_t X, rsb_type_t typecode, void *oVA, rsb_coo_idx_t * oINDX, rsb_nnz_idx_t * oXA)
{
	if(!iVA || !iINDX || !iXA || RSB_INVALID_COO_INDEX(X) || RSB_INVALID_NNZ_INDEX(nnz) || !oVA || !oINDX || !oXA)
		return RSB_ERR_BADARGS;
	RSB_CSR_MEMCPY(oVA,oINDX,oXA,iVA,iINDX,iXA,nnz,X,RSB_SIZEOF(typecode));
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_WANT_OSKI_BENCHMARKING 
rsb_err_t rsb__allocate_csc_arrays_from_coo_sorted(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_type_t typecode, void **VAp, rsb_coo_idx_t ** indxp, rsb_nnz_idx_t ** indptrp)
{
	return rsb__allocate_csr_arrays_from_coo_sorted(VA, JA, IA, nnz, k, m, typecode, VAp, indxp, indptrp);
}

rsb_err_t rsb__allocate_csr_arrays_from_coo_sorted(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_type_t typecode, void **VAp, rsb_coo_idx_t ** indxp, rsb_nnz_idx_t ** indptrp)
{
	/*!
	 	\ingroup gr_internals
		
		FIXME : UNFINISHED, UNTESTED, NEW, and SLOW
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t * indx = NULL,MI = 0;
	rsb_coo_idx_t * indpntr = NULL;
	void *cVA = NULL;

	if(!indxp || !indptrp || !JA || !IA) /* VA is optional */
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	indx = rsb__clone_area(JA,sizeof(rsb_coo_idx_t)*(nnz));	/* FIXME ! BAD ! */
	indpntr = rsb__calloc(sizeof(rsb_nnz_idx_t)*(m+1));

	if(!indx || !indpntr)
	{
		errval = RSB_ERR_ENOMEM;
		goto err;
	}

	if(VA)
	{
		cVA = rsb__clone_area(VA,RSB_SIZEOF(typecode)*(nnz));
		if(!cVA)
		{
			errval = RSB_ERR_ENOMEM;
			goto err;
		}
	}

	errval = rsb__do_account_sorted_optimized_css(IA,JA,m,k,nnz,indpntr+1,NULL);

	if(RSB_SOME_ERROR(errval))
		goto err;

	for(MI=0;MI<m;++MI)
		indpntr[MI+1] += indpntr[MI];

	if(VAp)*VAp = cVA;
	*indxp = indx;
	*indptrp = indpntr;

	goto ok;
err:
	RSB_CONDITIONAL_FREE(cVA);
	RSB_CONDITIONAL_FREE(indx);
	RSB_CONDITIONAL_FREE(indpntr);
ok:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_OSKI_BENCHMARKING */

rsb_err_t rsb__print_configuration_string(const char *pn, rsb_char_t * cs, rsb_bool_t wci)
{
	/*!
	 	\ingroup gr_internals
	*/
	/*
 		note : an interleaved 
			RSB_INFO(
			#ifdef FOO
			"string"
			# endif
			)
			style of string composition breaks xlc, so we don't use it, in the following.
		FIXME: pn should be reasonably short, as this routine does NOT check for buffer overflow. for this same reason this function shall remain INTERNAL always.	
 	*/
	rsb_char_t buf[RSB_MAX_VERSION_STRING_LENGTH];
	const rsb_char_t * const sep = " ";
	const rsb_char_t * const nl = "\n";

	rsb__strcpy(buf,"");
	if(!cs)
	{
		// TODO: missing error handling
	       	goto done;
	}
#if 0
	rsb__sprintf(buf,"%s version: %d.%d.%d\n",pn?pn:"librsb",RSB_LIBRSB_VER_MAJOR,RSB_LIBRSB_VER_MINOR,RSB_LIBRSB_VER_RELEASE);
#else
	rsb__sprintf(buf,"%s version: %s\n",pn?pn:"librsb",RSB_LIBRSB_VER_STRING);
#endif
	if(wci == RSB_BOOL_FALSE)
	{
		rsb__sprintf(cs,"%s",buf);
            	rsb__sprintf(cs+strlen(cs),"%s.\n\n",RSB_COPYRIGHT_STRING);
            	rsb__sprintf(cs+strlen(cs),"Written by %s.\n",RSB_PACKAGE_BUGREPORT);
		goto done;
	}
	rsb__util_strcat(buf,"format switches:");
#ifdef RSB_MATRIX_STORAGE_BCSR_STRING 
	rsb__util_strcat(buf,"br");
	rsb__util_strcat(buf,sep);
#endif /* RSB_MATRIX_STORAGE_BCSR_STRING */
#ifdef RSB_MATRIX_STORAGE_BCSC_STRING 
	rsb__util_strcat(buf,"bc");
	rsb__util_strcat(buf,sep);
#endif /* RSB_MATRIX_STORAGE_BCSC_STRING */
#ifdef RSB_MATRIX_STORAGE_VBR_STRING 
	rsb__util_strcat(buf,"vr");
	rsb__util_strcat(buf,sep);
#endif /* RSB_MATRIX_STORAGE_VBR_STRING */
#ifdef RSB_MATRIX_STORAGE_VBC_STRING 
	rsb__util_strcat(buf,"vc");
	rsb__util_strcat(buf,sep);
#endif /* RSB_MATRIX_STORAGE_VBC_STRING */
#ifdef RSB_MATRIX_STORAGE_LC_STRING 
	rsb__util_strcat(buf,"lc");
	rsb__util_strcat(buf,sep);
#endif /* RSB_MATRIX_STORAGE_LC_STRING */
#ifdef RSB_MATRIX_STORAGE_LR_STRING 
	rsb__util_strcat(buf,"lr");
	rsb__util_strcat(buf,sep);
#endif /* RSB_MATRIX_STORAGE_LR_STRING */
	rsb__util_strcat(buf,nl);
	rsb__util_strcat(buf,"ops:");
	rsb__util_strcat(buf,RSB_M4_MATRIX_META_OPS_STRING);
	rsb__util_strcat(buf,nl);

	rsb__util_strcat(buf,"types:");
	rsb__util_strcat(buf,RSB_M4_MATRIX_TYPES_STRING);
	rsb__util_strcat(buf,nl);
	rsb__util_strcat(buf,"type char codes:");
	rsb__util_strcat(buf,RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS );
	rsb__util_strcat(buf,nl);
	rsb__util_strcat(buf,"types count:");
	rsb__sprintf(buf+strlen(buf),"%d",RSB_IMPLEMENTED_TYPES);
	rsb__util_strcat(buf,nl);
	rsb__util_strcat(buf,"transposition codes:");
	rsb__util_strcat(buf,RSB_TRANSPOSITIONS_PREPROCESSOR_SYMBOLS );
	rsb__util_strcat(buf,nl);

	rsb__util_strcat(buf,"restrict keyword is: ");
#ifdef RSB_restrict
	rsb__util_strcat(buf,"on" );
#else /* RSB_restrict */
	rsb__util_strcat(buf,"off" );
#endif /* RSB_restrict */
	rsb__util_strcat(buf,nl);

	rsb__util_strcat(buf,"row unrolls:");
	rsb__util_strcat(buf,RSB_M4_WANT_COLUMN_UNLOOP_FACTORS_STRING);
	rsb__util_strcat(buf,nl);
	rsb__util_strcat(buf,"column unrolls:");
	rsb__util_strcat(buf,RSB_M4_WANT_ROW_UNLOOP_FACTORS_STRING	);
	rsb__util_strcat(buf,nl);
	rsb__util_strcat(buf,"reference benchmark sample minimum time (seconds):%lg\n");
	rsb__util_strcat(buf,"reference benchmark sample minimum runs:%zd\n");
	rsb__util_strcat(buf,"maximal configured block size:%zd\n");
#ifdef RSB_WANT_OSKI_BENCHMARKING 
	rsb__util_strcat(buf,"oski comparative benchmarking enabled\n");
#endif /* RSB_WANT_OSKI_BENCHMARKING */
	rsb__util_strcat(buf,"sizeof(rsb_nnz_idx_t):%zd\n");
	rsb__util_strcat(buf,"sizeof(rsb_coo_idx_t):%zd\n");
	rsb__util_strcat(buf,"sizeof(rsb_blk_idx_t):%zd\n");
	//rsb__util_strcat(buf,"sizeof(rsb_thread_t):%zd\n");
	rsb__util_strcat(buf,"sizeof(size_t):%zd\n");
	rsb__util_strcat(buf,"sizeof(struct rsb_mtx_t):%zd\n");
	rsb__util_strcat(buf,"sizeof(struct rsb_blas_sparse_matrix_t):%zd\n");
	rsb__util_strcat(buf,"sizeof(struct rsb_coo_mtx_t):%zd\n");
	rsb__util_strcat(buf,"RSB_MAX_MATRIX_DIM:%zd\n");
	rsb__util_strcat(buf,"RSB_MAX_MATRIX_NNZ:%zd\n");
	rsb__util_strcat(buf,"RSB_CONST_MAX_SUPPORTED_CORES:%zd\n");
	rsb__util_strcat(buf,"RSB_BLAS_MATRICES_MAX:%zd\n");
	rsb__util_strcat(buf,"RSB_CONST_MIN_NNZ_PER_ROW_FOR_COO_SWITCH:%zd\n");

	rsb__util_strcat(buf,"RSB_USER_SET_MEM_HIERARCHY_INFO:%s\n");
	rsb__util_strcat(buf,"RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t):%zd\n");
	rsb__util_strcat(buf,"RSB_IOLEVEL:%d\n");
	//RSB_INFO(
	rsb__sprintf(cs,
		buf,
		RSB_BENCHMARK_MIN_SECONDS,
		(rsb_printf_int_t)RSB_BENCHMARK_MIN_RUNS,
		(rsb_printf_int_t)RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE,
		sizeof(rsb_nnz_idx_t),
		sizeof(rsb_coo_idx_t),
		sizeof(rsb_blk_idx_t),
		//sizeof(rsb_thread_t),
		sizeof(size_t),
		sizeof(struct rsb_mtx_t),
		sizeof(struct rsb_blas_sparse_matrix_t),
		sizeof(struct rsb_coo_mtx_t),
		(rsb_printf_int_t)RSB_MAX_MATRIX_DIM,
		(rsb_printf_int_t)RSB_MAX_MATRIX_NNZ,
		(rsb_printf_int_t)RSB_CONST_MAX_SUPPORTED_CORES,
		(rsb_printf_int_t)RSB_BLAS_MATRICES_MAX,
		(rsb_printf_int_t)RSB_CONST_MIN_NNZ_PER_ROW_FOR_COO_SWITCH
		,(rsb_printf_int_t)rsb__init_get_mem_hierarchy_info_string(RSB_BOOL_FALSE)?rsb__init_get_mem_hierarchy_info_string(RSB_BOOL_FALSE):NULL
		,RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t)
		,RSB_IOLEVEL 
	);
done:
	return RSB_ERR_NO_ERROR;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_blk_idx_t rsb__recursive_middle_block_index(rsb_blk_idx_t i)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the index of the split point index.
	 */
#if RSB_EXPERIMENTAL_MORTON_ORDERED_RECURSION  
	rsb_blk_idx_t s = 0;
	while( (1<<s+1) < i)
		++s;
	return (1<<s);
#else /* RSB_EXPERIMENTAL_MORTON_ORDERED_RECURSION */
#if 1
	return (i+1)/2;
#else
	int p = 0;
	while(i>>(p+1) && i > (1<<(p+1)))
		++p;
//	RSB_INFO("%d %d %d\n",i,p,(1<<p));
	return (1<<p);
#endif
#endif /* RSB_EXPERIMENTAL_MORTON_ORDERED_RECURSION */
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__recursive_middle_index(const struct rsb_mtx_partitioning_info_t * pinfop, rsb_coo_idx_t * M_bp, rsb_coo_idx_t * K_bp )
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the index of ..
	 */
	rsb_blk_idx_t mBi;
	rsb_blk_idx_t	kBi;
	if(!pinfop || !M_bp || !K_bp)
		return RSB_ERR_BADARGS;
	mBi = rsb__recursive_middle_block_index(pinfop->M_b);
	kBi = rsb__recursive_middle_block_index(pinfop->K_b);
	if(pinfop->rpntr)
		*M_bp = pinfop->rpntr[rsb__recursive_middle_block_index(mBi)];
	else
		*M_bp = rsb__recursive_middle_block_index(pinfop->nr);
	if(pinfop->cpntr)
		*K_bp = pinfop->cpntr[rsb__recursive_middle_block_index(kBi)];
	else
		*K_bp = rsb__recursive_middle_block_index(pinfop->nc);
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__recursive_split_point_parms_get(
		const struct rsb_mtx_partitioning_info_t * pinfop,
		rsb_coo_idx_t * moff, rsb_coo_idx_t * koff)
{
	/*!
		\ingroup gr_internals

		FIXME: this is function is OBSOLETE
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!pinfop || !pinfop->rpntr || !pinfop->cpntr)
	{
		errval = RSB_ERR_BADARGS;
	}
	RSB_DEBUG_ASSERT(moff);
	RSB_DEBUG_ASSERT(koff);
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(pinfop->M_b));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(pinfop->K_b));
	
	*moff = pinfop->rpntr[ rsb__recursive_middle_block_index(pinfop->M_b) ];
	*koff = pinfop->cpntr[ rsb__recursive_middle_block_index(pinfop->K_b) ];

	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(*moff));
	RSB_DEBUG_ASSERT(RSB_IS_VALID_COO_INDEX(*koff));

	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_long_t rsb__terminal_recursive_matrix_count(const struct rsb_mtx_t *mtxAp)
{
	/*!
	 * \ingroup gr_internals
	 * \return the count of leaf (terminal) matrices
	 *
	 * TODO : change this function type!
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_long_t smc = 0;

	if(!mtxAp)
	{smc = 0;goto done;}

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{smc = 1;goto done;}

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
		smc += rsb__terminal_recursive_matrix_count(submatrix);
done:
	return smc;
}

rsb_err_t rsb__do_compute_terminal_nnz_min_max_avg_count(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * minnz, rsb_nnz_idx_t * maxnz, rsb_nnz_idx_t * avgnz)
{
//	struct rsb_mtx_t * submatrix = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if(minnz)
		*minnz = RSB_MAX_MATRIX_NNZ;
	if(maxnz)
		*maxnz = 0;
	if(avgnz)
		*avgnz = 0;
	errval = rsb__do_compute_terminal_nnz_min_max_count(mtxAp,minnz,maxnz);
	if(avgnz)
		*avgnz = mtxAp->nnz/rsb__terminal_recursive_matrix_count(mtxAp);
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_compute_terminal_nnz_min_max_count(const struct rsb_mtx_t *mtxAp, rsb_nnz_idx_t * minnz, rsb_nnz_idx_t * maxnz)
{
	/*!
	 * \ingroup gr_internals
	 * \return ...
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_ES);
		errval = RSB_ERR_BADARGS;
	}

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		if(minnz)
//			RSB_STDOUT_MATRIX_SUMMARY(mtxAp), RSB_INFO(" <- MIN (from %d)\n",*minnz),
			*minnz = RSB_MIN(*minnz,mtxAp->nnz);
		if(maxnz)
//			RSB_STDOUT_MATRIX_SUMMARY(mtxAp), RSB_INFO(" <- MAX (from %d)\n",*maxnz),
			*maxnz = RSB_MAX(*maxnz,mtxAp->nnz);
	}
	else
	{
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_compute_terminal_nnz_min_max_count(submatrix,minnz,maxnz));
	}
//done:
	RSB_DO_ERR_RETURN(errval)
}

rsb_long_t rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(const struct rsb_mtx_t *mtxAp, rsb_fmt_t matrix_storage, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * \return the count of leaf (terminal) matrices
	 *
	 * TODO : change this function type!
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_long_t smc = 0;

	if(!mtxAp)
	{smc=0;goto done;}

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		if((!RSB_DO_FLAG_HAS(mtxAp->flags,flags)) && (mtxAp->matrix_storage==matrix_storage))
			smc = 1;
		goto done;
	}

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
		smc += rsb__terminal_recursive_matrix_count_with_storage_and_no_flags(submatrix,matrix_storage,flags);
done:
	return smc;
}

rsb_long_t rsb__terminal_recursive_matrix_count_with_storage_and_flags(const struct rsb_mtx_t *mtxAp, rsb_fmt_t matrix_storage, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * \return the count of leaf (terminal) matrices
	 *
	 * TODO : change this function type!
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_long_t smc = 0;

	if(!mtxAp)
	{smc = 0;goto done;}

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,flags) && (mtxAp->matrix_storage==matrix_storage))
			smc = 1;
		goto done;
	}

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
		smc += rsb__terminal_recursive_matrix_count_with_storage_and_flags(submatrix,matrix_storage,flags);
done:
	return smc;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_long_t rsb__terminal_recursive_matrix_count_with_flags_but(const struct rsb_mtx_t *mtxAp, rsb_flags_t flags, rsb_flags_t nflags)
{
	/*!
	 * \ingroup gr_internals
	 * \return the count of leaf (terminal) matrices
	 *
	 * TODO : change this function type!
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_long_t smc = 0;

	if(!mtxAp)
	{smc = 0;goto done;}

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,flags) && !RSB_DO_FLAG_HAS(mtxAp->flags,nflags))
			smc = 1;
		goto done;
	}

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
		smc += rsb__terminal_recursive_matrix_count_with_flags_but(submatrix,flags,nflags);
done:
	return smc;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_long_t rsb__terminal_recursive_matrix_count_with_flags(const struct rsb_mtx_t *mtxAp, rsb_flags_t flags)
{
	/*!
	 * \ingroup gr_internals
	 * \return the count of leaf (terminal) matrices
	 *
	 * TODO : change this function type!
	 */
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_long_t smc = 0;

	if(!mtxAp)
	{smc = 0;goto done;}

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		if(RSB_DO_FLAG_HAS(mtxAp->flags,flags))
			smc = 1;
		goto done;
	}

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
		smc += rsb__terminal_recursive_matrix_count_with_flags(submatrix,flags);
done:
	return smc;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_trans_t rsb__do_transposition_from_char(rsb_char_t tc)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : document this
	 */
	if(tolower(tc)=='t')
		return RSB_TRANSPOSITION_T;
	else
	if(tolower(tc)=='n')
		return RSB_TRANSPOSITION_N;
	else
	if(tolower(tc)=='c')
		return RSB_TRANSPOSITION_C;
	else
		return RSB_INVALID_TRANS;
}

rsb_trans_t rsb__do_transpose_transposition(rsb_trans_t transA)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : document this
	 */
	if(RSB_DOES_NOT_TRANSPOSE(transA))
		return RSB_TRANSPOSITION_T;
	if(transA == RSB_TRANSPOSITION_T)
		return RSB_TRANSPOSITION_N;
	if(transA == RSB_TRANSPOSITION_C)
		return RSB_TRANSPOSITION_N;
	return transA;
}

#if 0
rsb_err_t rsb_spmm_inner(const struct rsb_mtx_t * mtxAp, const void * mrhs, void *mout, rsb_int_t bstride, rsb_int_t cstride, rsb_int_t nrhs, rsb_trans_t transA)
{
#ifdef RSB_HAVE_OPTYPE_SPMM_AZ
	/*!
	 * \ingroup gr_internals
	 * fixme */

	size_t el_size = mtxAp->el_size;

	if( rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix = NULL;
		rsb_coo_idx_t mB = (mtxAp->rpntr[rsb__recursive_middle_block_index(mtxAp->M_b)]);
		rsb_coo_idx_t kB = (mtxAp->cpntr[rsb__recursive_middle_block_index(mtxAp->K_b)]);

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			rsb_coo_idx_t moff,koff;

			moff = i*mB;
			koff = j*kB;
	
			rsb_spmm_inner(submatrix,((rsb_byte_t*)mrhs)+koff*el_size,((rsb_byte_t*)mout)+moff*el_size,bstride,cstride,nrhs,transA);
		}
		return RSB_ERR_NO_ERROR;
	}
	else
		return rsb__do_spmm_az(mtxAp,mrhs,mout,bstride,cstride,nrhs,transA);/*FIXME*/
#else /* RSB_HAVE_OPTYPE_SPMM_AZ */
	return RSB_ERR_UNSUPPORTED_OPERATION;
#endif /* RSB_HAVE_OPTYPE_SPMM_AZ */
}
#endif

rsb_flags_t rsb__do_flip_uplo_flags(rsb_flags_t flags)
{
	/**
	 * \ingroup gr_internals
	 */
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER))
		RSB_DO_FLAG_SUBST(flags,RSB_FLAG_UPPER,RSB_FLAG_LOWER);
	else
		if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER))
			RSB_DO_FLAG_SUBST(flags,RSB_FLAG_LOWER,RSB_FLAG_UPPER);
	return flags;
}

rsb_flags_t rsb__do_detect_and_add_triangular_flags(rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	/**
	 * \ingroup gr_internals
	 */
	rsb_nnz_idx_t i;
	if(!IA || !JA || RSB_INVALID_NNZ_INDEX(nnz))
		return flags;
	/* FIXME: this code could be optimized a great deal, by introducing a state machine like scan. */
	for(i=0;i<nnz;++i)
	{
		if(IA[i]==JA[i])
			continue;
		if(IA[i]>JA[i])
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
		else
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
	}
	if(RSB_NAND(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER),RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER)))
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
		
	/* it could be the case that both RSB_FLAG_LOWER and RSB_FLAG_UPPER flags get caught */
	return flags;
}

rsb_err_t rsb__do_load_matrix_file_as_matrix_market(struct rsb_mtx_t ** mtxApp, const rsb_char_t * filename, rsb_flags_t flags, rsb_type_t typecode)
{
	/*!
	 * FIXME: and typecode check ?
	 * FIXME: UNFINISHED
	 */
	/** \ingroup gr_unfinished */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	void *VA = NULL;
	rsb_coo_idx_t *IA = NULL, *JA = NULL;

	if(!mtxApp)
	{
		RSB_ERROR(RSB_ERRM_E_MTXAP);
		errval = RSB_ERR_BADARGS;
	}
	else
	if(!filename)
	{
		RSB_ERROR(RSB_ERRM_NPSFF);
		errval = RSB_ERR_BADARGS;
	}
	else
	{
		rsb_coo_idx_t m = 0,k = 0;
		rsb_nnz_idx_t nnz = 0;
#define RSB_20120309_FIX 1
#if RSB_20120309_FIX
		rsb_bool_t is_symmetric = RSB_BOOL_FALSE;
		rsb_bool_t is_hermitian = RSB_BOOL_FALSE;
		rsb_bool_t is_lower = RSB_BOOL_FALSE;
		rsb_bool_t is_upper = RSB_BOOL_FALSE;
		rsb_bool_t is_vector = RSB_BOOL_FALSE;
		/* FIXME: shall update test_matops.c accordingly. */
		/* FIXME: need limits checks! */
		if(RSB_SOME_ERROR(rsb__util_mm_info_matrix_f(filename,&m,&k,&nnz,NULL,&is_symmetric,&is_hermitian,NULL,&is_lower,&is_upper,&is_vector) ) || is_vector )
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err,RSB_ERRMSG_NOTMTXMKT" : %s ..\n",filename);
		}
		else
		{
			if( is_symmetric == RSB_BOOL_TRUE ) RSB_DO_FLAG_ADD(flags,RSB_FLAG_SYMMETRIC);
			if( is_hermitian == RSB_BOOL_TRUE ) RSB_DO_FLAG_ADD(flags,RSB_FLAG_HERMITIAN);
			if(is_upper) RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
			if(is_lower) RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
		}

		if( m==k && m>1 /* && want_only_lowtri*/ )
			nnz += m;	/* the loading routine shall allocate nnz+m */
		else
 			nnz = 0;	/* the loading routine should determine nnz */
#endif /* RSB_20120309_FIX */
		if((errval = rsb__util_mm_load_matrix_f(filename,&IA,&JA,&VA,&m,&k,&nnz,typecode,flags,NULL,NULL))!=RSB_ERR_NO_ERROR)
		{
			RSB_ERROR(RSB_ERRM_ES);
			rsb__do_perror(NULL,errval);
			goto err;
		}
		if((mtxAp = rsb__do_mtx_alloc_from_coo_inplace(VA,IA,JA,nnz,typecode,m,k,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags,&errval))==NULL)
		{
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
			// FIXME: incomplete error handling
		}
		if(mtxAp)
		{
			rsb_submatrix_idx_t smi;
			struct rsb_mtx_t * submatrix;

			RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
				RSB_DO_FLAG_DEL(submatrix->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
			RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
		}
		*mtxApp = mtxAp;
		goto ok;
	}
err:
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(VA);
ok:
	RSB_DO_ERR_RETURN(errval)
}

rsb_bool_t rsb__are_coo_matrices_equal(const struct rsb_coo_mtx_t *cm1, const struct rsb_coo_mtx_t *cm2)
{
	/* this is a debug routine and may print internal stuff */
	const rsb_bool_t no = RSB_BOOL_FALSE;
	const rsb_bool_t yes = RSB_BOOL_TRUE;
	rsb_nnz_idx_t nnz;
#if (RSB_WANT_VERBOSE_MESSAGES || 1)
#define	RSB_GOTO_DIFFERING	{ RSB_PERR_GOTO(differing,RSB_ERRM_ES);}
#else /* RSB_WANT_VERBOSE_MESSAGES */
#define	RSB_GOTO_DIFFERING	{goto differing;}
#endif /* RSB_WANT_VERBOSE_MESSAGES */
	if( cm1 ==  cm2)
		goto equal;
	if(!cm1 || !cm2)
		RSB_GOTO_DIFFERING
	nnz = cm1->nnz;
	if( RSB_INVALID_NNZ_INDEX(nnz) )
		RSB_GOTO_DIFFERING
	if(cm1->nr!= cm2->nr)
		RSB_GOTO_DIFFERING
	if(cm1->nc!= cm2->nc)
		RSB_GOTO_DIFFERING
	if(cm1->nnz!= cm2->nnz)
		RSB_GOTO_DIFFERING
	if(cm1->typecode!= cm2->typecode)
		RSB_GOTO_DIFFERING
	//if((!cm1->IA)&&(!cm2->IA)) return no;
	//if((!cm1->IA)||(!cm2->IA)) return no;
	//else
		if((cm1->IA)&&(cm2->IA))
		if(RSB_MEMCMP(cm1->IA,cm2->IA,sizeof(rsb_coo_idx_t)*nnz))
			RSB_GOTO_DIFFERING
	//if((!cm1->JA)&&(!cm2->JA)) return no;
	//if((!cm1->JA)||(!cm2->JA)) return no;
	//else
		if((cm1->JA)&&(cm2->JA))
		if(RSB_MEMCMP(cm1->JA,cm2->JA,sizeof(rsb_coo_idx_t)*nnz))
			RSB_GOTO_DIFFERING
	//if((!cm1->VA)&&(!cm2->VA)) return no;
	//if((!cm1->VA)||(!cm2->VA)) return no;
	//else
		if((cm1->VA)&&(cm2->VA))
		{
#if 1
			if(rsb__do_are_same(cm1->VA,cm2->VA, nnz,cm1->typecode, 1, 1))
				RSB_GOTO_DIFFERING
#else
			if(RSB_MEMCMP(cm1->VA,cm2->VA,RSB_SIZEOF(cm1->typecode)*nnz)) /* This is too strict: for it, -0.0 != 0.0 */
				RSB_GOTO_DIFFERING
#endif
		}


equal:
	return yes;
differing:
#if RSB_ALLOW_STDOUT
#if (RSB_WANT_VERBOSE_MESSAGES || 1)
	if(cm1)RSB_STDOUT_COO_MATRIX_SUMMARY(cm1);
	if(cm2)RSB_STDOUT_COO_MATRIX_SUMMARY(cm2);
#endif /* RSB_WANT_VERBOSE_MESSAGES */
#endif /* RSB_ALLOW_STDOUT */
	return no;
#undef	RSB_GOTO_DIFFERING
}

static rsb_bool_t rsb__is_coo_matrix_empty(const struct rsb_coo_mtx_t *cm, rsb_flags_t flags)
{
	if(!cm)
		return RSB_BOOL_FALSE;
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
		return RSB_BOOL_FALSE;
	if(cm->nnz==0)
		return RSB_BOOL_TRUE;
	return (rsb__check_for_nonzeros(cm->VA,cm->nnz,cm->typecode)==0)?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
}

rsb_bool_t rsb__are_coo_matrices_both_empty(const struct rsb_coo_mtx_t *cm1, rsb_flags_t flags1, const struct rsb_coo_mtx_t *cm2, rsb_flags_t flags2)
{
	rsb_bool_t im1e = RSB_BOOL_FALSE,im2e = RSB_BOOL_FALSE,abme = RSB_BOOL_FALSE;
	im1e = rsb__is_coo_matrix_empty(cm1,flags1);
	im2e = rsb__is_coo_matrix_empty(cm2,flags2);
	abme = RSB_BOOL_OR(im1e,im2e);
	return abme;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__are_coo_matrices_equal_or_both_empty(const struct rsb_coo_mtx_t *cm1, rsb_flags_t flags1, const struct rsb_coo_mtx_t *cm2, rsb_flags_t flags2)
{
	rsb_bool_t acme = rsb__are_coo_matrices_equal(cm1,cm2);
	if(acme)
		;
	else
		acme = rsb__are_coo_matrices_both_empty(cm1,flags1,cm2,flags2);
	return acme;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__get_row_dense(const struct rsb_mtx_t * mtxAp, void* row, rsb_coo_idx_t i )
{
	/*!
	 * \ingroup gr_mops
	 * Will write entire row i of matrix matrix in the row vector.
	 *
	 * \param i the specified row
	 * \param row an already allocated vector of the same type as the matrix
	 * \param mtxAp is a valid pointer to a rsb_mtx_t structure
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * FIXME: this function is unfinished.
	 *
	 * \note This function is obsolete -- it's only used by ot.c !
	 * */
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
	if(!mtxAp)
		return RSB_ERR_BADARGS;

	RSB_BZERO( row , mtxAp->el_size * mtxAp->nc);

	return rsb__do_get_row_dense(mtxAp, row, i );
}

#if 0
rsb_err_t rsb_spmv_unua(const struct rsb_mtx_t * mtxAp, const void * x, void * y, rsb_trans_t transA)
{
	/*!
	 * \ingroup gr_mops
	 * computes \f$y \leftarrow y - op(A) \cdot x \f$
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * 
	 * */
	rsb_aligned_t mone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb__util_set_area_to_converted_integer(&mone[0],mtxAp->typecode,-1);
	return rsb__do_spmv_general(transA,&mone[0],mtxAp,x,1,NULL,y,1,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
}

rsb_err_t rsb_spmv_uaua(const struct rsb_mtx_t * mtxAp, const void * rhs, void * out, rsb_trans_t transA)
{
	/*!
	 * \ingroup gr_mops
	 * computes \f$y \leftarrow y + op(A) \cdot x \f$
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * 
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	errval = rsb__do_spmv_general(transA,NULL,mtxAp,rhs,1,NULL,out,1,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb_spmv_uauz(const struct rsb_mtx_t * mtxAp, const void * rhs, void * out, rsb_trans_t transA)
{
	/*!
	 * \ingroup gr_mops
	 * computes \f$y \leftarrow op(A) \cdot x \f$
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * */

	/* FIXME : TEMPORARY */
	rsb_aligned_t zero[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];

	if(!mtxAp)
		return RSB_ERR_BADARGS;
	rsb__util_set_area_to_converted_integer(&zero[0],mtxAp->typecode,0);
	return rsb__do_spmv_general(transA,NULL,mtxAp,rhs,1,&zero[0],out,1,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
}
#endif

#if 0
static rsb_err_t rsb_spmv_sxsx(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_trans_t transA, rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{
	/*!
	 * \ingroup gr_mops
	 * computes \f$y \leftarrow \beta \cdot y + \alpha\cdot A\cdot x\f$
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 */
	if(!alphap || !betap)
		return RSB_ERR_BADARGS;
	return rsb__do_spmv_general(transA,alphap,mtxAp,x,incx,betap,y,incy,RSB_OP_FLAG_DEFAULT RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
}
#endif



struct rsb_mtx_t * rsb__load_matrix_file_as_binary(const rsb_char_t * filename, rsb_err_t *errvalp)
{
	/*!
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;

	if(!errvalp || !filename)
	{
		errval = RSB_ERR_BADARGS;
	}
	else
		errval = rsb__do_load_matrix_file_as_binary(&mtxAp,filename);
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

rsb_err_t rsb__do_spsm(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * betap, const void * b, rsb_nnz_idx_t ldb, void * c, rsb_nnz_idx_t ldc)
{
	/*!
	 * \ingroup gr_internals
	 * When x is a multivector with nrhs elements, b elements having stride ldb and c elements having stride ldc
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 * */
	 /* FIXME : and error detection ? **/
	 /* FIXME : UNTESTED, UNFINISHED **/
	rsb_coo_idx_t l = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if RSB_ALLOW_ZERO_DIM
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
	{
		goto err; /* FIXME: skipping further checks */
	}
#endif

	if(!mtxAp || !b || !ldb || !c || !ldc || !nrhs /*  transA*/ || !alphap || !betap ||
		( order != RSB_FLAG_WANT_COLUMN_MAJOR_ORDER && order != RSB_FLAG_WANT_ROW_MAJOR_ORDER ) )
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}


	if( 0
#if RSB_ENABLE_INNER_NRHS_SPSV
		|| ( rsb_global_session_handle.want_outer_spmm==0 )
#endif
	  ) /* 0 == yes TODO: need want_outer_spsm here */
	{
		size_t outnri = ldc, rhsnri = ldb;

		if(order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spsv_general(transT,alphap,mtxAp,b,   1,c,   1,RSB_OP_FLAG_DEFAULT RSB_INNER_NRHS_SPSV_ARGS_IDS));
		else
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spsv_general(transT,alphap,mtxAp,b,nrhs,c,nrhs,RSB_OP_FLAG_DEFAULT RSB_INNER_NRHS_SPSV_ARGS_IDS));
	}
	else
	{
		if(order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			/*  column major */
			for(l=0;l<nrhs;++l)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spsv(transT,alphap,mtxAp, ((const rsb_byte_t*)b)+(mtxAp->el_size*ldb)*l,1, ((rsb_byte_t*)c)+(mtxAp->el_size*ldc)*l,1));
		else
			/*  row major */
			for(l=0;l<nrhs;++l)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spsv(transT,alphap,mtxAp, ((const rsb_byte_t*)b)+(mtxAp->el_size*l),ldb, ((rsb_byte_t*)c)+(mtxAp->el_size*l),ldc));
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__util_sort_row_major_buffered(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k,  rsb_type_t typecode , rsb_flags_t flags, void * WA, size_t wb )
{
	/*!

	   Will sort as CSR (or CSC) the given coefficients.
	   Will ignore any fancy sorting flags.
	    
	   \param \rsb_wr_va_ia_ja_desc_msg
	   \param \rsb_flags_inp_param_msg
	   \param \rsb_nnz_inp_param_msg
	   \param \rsb_nrows_inp_param_msg
	   \param \rsb_ncols_inp_param_msg
	   \param \rsb_type_param_msg
	   \return \rsberrcodemsg
	*/
	// FIXME : should handle error conditions
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//struct rsb_mtx_partitioning_info_t pinfop;

	if(RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;goto err;
	}
#if 0
	pinfop.rpntr = rsb__util_get_partitioning_array( 1, m , &pinfop.M_b, flags), 
	pinfop.cpntr = rsb__util_get_partitioning_array( 1, k , &pinfop.K_b, flags),
	rsb__pinfo_init( &pinfop, pinfop.M_b, pinfop.K_b, pinfop.rpntr, pinfop.cpntr, m,k);/* FIXME : is this ok ? */
	errval = rsb__do_util_sortcoo(VA,IA,JA,m,k,nnz,typecode,&pinfop,flags,NULL,0);
	RSB_CONDITIONAL_FREE(pinfop.rpntr);
	RSB_CONDITIONAL_FREE(pinfop.cpntr);
#else
	//Remove internal support for RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT):
	//RSB_DO_FLAG_SUBST(flags,RSB_INTERNAL_FLAG_CSR_SORTING_MASK,RSB_FLAG_WANT_BCSS_STORAGE|RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT);
	RSB_DO_FLAG_SUBST(flags,RSB_INTERNAL_FLAG_CSR_SORTING_MASK,RSB_FLAG_WANT_BCSS_STORAGE);
	errval = rsb__do_util_sortcoo(VA,IA,JA,m,k,nnz,typecode,NULL,flags,WA,wb);
#endif
err:
	RSB_DO_ERR_RETURN(errval)
}

const rsb_char_t *rsb__basename(const rsb_char_t *path)
{
	rsb_int_t sl;

	if(!path)
		return path;
	sl = rsb__util_strlen(path);
	while(sl>0 && path[sl-1]==RSB_DIR_SEPARATOR)
		--sl;		/* not to get stuck on trailing slashes */
	while(sl>0 && path[sl-1]!=RSB_DIR_SEPARATOR)
		--sl;
	return path+sl;
}

rsb_err_t rsb__do_set_elements(struct rsb_mtx_t * mtxAp, const void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	rsb_coo_idx_t k;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t ifo = ( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )?1:0;

	if(!IA || !VA || !JA || !mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

#if RSB_WANT_COO_BEGIN 
	if( RSB_MTX_HBDF( mtxAp) )
	{
		errval = rsb__BLAS_Xuscr_insert_entries(RSB_MTX_HBDFH(mtxAp),nnz,VA,IA,JA);
		goto err;
	}
#endif /* RSB_WANT_COO_BEGIN */
	for(k=0;k<nnz;++k)
		errval |= rsb__do_upd_coo_element(mtxAp,((const rsb_char_t*)VA)+mtxAp->el_size*k,IA[k]-ifo,JA[k]-ifo,flags);
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__set_ldX_for_spmm(rsb_trans_t transA, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, rsb_nnz_idx_t * ldBp, rsb_nnz_idx_t * ldCp)
{
	if ( order == RSB_FLAG_WANT_ROW_MAJOR_ORDER )
	{
		if( !*ldBp )
			*ldBp = nrhs;
		if( !*ldCp )
			*ldCp = nrhs;
	}
	else
	if ( order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
	{
		if( !*ldBp )
		{
			if ( transA == RSB_TRANSPOSITION_N )
				*ldBp = mtxAp->nc;
			else
				*ldBp = mtxAp->nr;
		}
		if( !*ldCp )
		{
			if ( transA == RSB_TRANSPOSITION_N )
				*ldCp = mtxAp->nr;
			else
				*ldCp = mtxAp->nc;
		}
	}
	else
		goto err;

	if( !(*ldBp) || (!*ldCp) )
		goto err;

	if ( order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
	{
		if( *ldCp < rsb__do_get_rows_of(mtxAp,transA) || *ldBp < rsb__do_get_columns_of(mtxAp,transA) )
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	else
	if ( order == RSB_FLAG_WANT_ROW_MAJOR_ORDER )
	{
		if( *ldCp < nrhs || *ldBp < nrhs )
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_BADARGS;
}

rsb_err_t rsb__do_spmm(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * b, rsb_nnz_idx_t ldb, const void * betap, void * c, rsb_nnz_idx_t ldc, enum rsb_op_flags_t op_flags)
{
	/*!
	   \ingroup rsb_doc_matrix_operations

	   FIXME: this function name is WRONG

	   when x is a multivector with nrhs elements, b elements having stride ldb and c elements having stride ldc

	   \return \rsberrcodemsg
	 */
	 /* FIXME : and error detection ?
	  * e.g.: 
	    if(order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER && ldb<mtxAp->nr && transA=...)
	  * **/
	 /* TODO: incx,incy  **/
	rsb_coo_idx_t l = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_coo_idx_t incx = 1,incy = 1;
	const rsb_byte_t * bp = b;
       	rsb_byte_t * cp = c;
	rsb_aligned_t pone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];

#if RSB_ALLOW_ZERO_DIM 
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
		goto err; /* FIXME: skipping checks on ldB, ldC, op_flags*/
#endif

	if(!mtxAp || !b || !c || !nrhs /*  transA*/ ||
		( order != RSB_FLAG_WANT_COLUMN_MAJOR_ORDER && order != RSB_FLAG_WANT_ROW_MAJOR_ORDER ) )
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	if(!alphap || !betap)
	{
		rsb__fill_with_ones(&pone[0],mtxAp->typecode,1,1);
		if(!alphap)
			alphap = &pone[0];
		if(!betap)
			betap = &pone[0];
	}

	errval = rsb__set_ldX_for_spmm(transA, mtxAp, nrhs, order, &ldb, &ldc);
	if(RSB_SOME_ERROR(errval))
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

#ifdef RSB_HAVE_OPTYPE_SPMM_AZ
	/*  
		return ...
	 SPMM_AZ is not yet complete
	 */
#endif /* RSB_HAVE_OPTYPE_SPMM_AZ */
	if( rsb_global_session_handle.want_outer_spmm==0 ) /* 0 == yes */
	{
		/* inner loop: fast */
		if(order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transA,alphap,mtxAp, bp,incx, betap, cp,incy,op_flags,nrhs,ldc, ldb  ));
		else
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transA,alphap,mtxAp, bp,ldb , betap, cp, ldc ,op_flags,nrhs,incx,incy));
	}
	else
	{
		/* outer loop: slow */
		if(order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			for(l=0;l<nrhs;++l)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transA,alphap,mtxAp,bp+(mtxAp->el_size*ldb)*l,incx, betap, cp+(mtxAp->el_size*ldc)*l,incy,op_flags RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS));
		else
			for(l=0;l<nrhs;++l)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transA,alphap,mtxAp,bp+(mtxAp->el_size*l  )  ,ldb , betap, cp+(mtxAp->el_size*l  )  ,ldc ,op_flags RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS));
	}
err:
	RSB_DO_ERR_RETURN(errval)

}

rsb_err_t rsb__do_spmm_general(const struct rsb_mtx_t * mtxAp, const void * b, void * c, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA, enum rsb_op_flags_t op_flags, rsb_flags_t order,const rsb_int_t nrhs, const size_t outnri, const size_t rhsnri)
{
	/* */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(nrhs != 1)
		errval = rsb__do_spmm(transA, alphap, mtxAp, nrhs, order, b, rhsnri, betap, c, outnri, op_flags);
	else
		errval = rsb__do_spmv_general(transA, alphap, mtxAp, b, incx, betap, c, incy, op_flags RSB_OUTER_NRHS_SPMV_ARGS_IDS);
	return errval;
}

rsb_err_t rsb__do_transpose(struct rsb_mtx_t ** mtxApp, rsb_bool_t want_conj)
{ 
	// TODO: and what to to if data arrays are externally allocated ?
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t*tmatrix = NULL;
	struct rsb_mtx_t*mtxAp = NULL;
	struct rsb_coo_mtx_t coo;
	struct rsb_mtx_t *fm = NULL;

	if(!mtxApp || !*mtxApp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err, RSB_ERRM_E_MTXAP);
	}
	mtxAp = *mtxApp;
	RSB_DO_FLAG_FLIP_UPLO(mtxAp->flags);
	RSB_INIT_CXX_FROM_MTX(&coo,mtxAp);
	if(rsb__allocate_coo_matrix_t(&coo)!=&coo)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	errval = rsb__do_get_rows_sparse(RSB_TRANSPOSITION_N,NULL,mtxAp,coo.VA,coo.IA,coo.JA,0,mtxAp->nr-1,&coo.nnz,RSB_FLAG_NOFLAGS);
	if(RSB_SOME_ERROR(errval))
	{
		goto err;
	}
	fm = rsb__do_get_first_submatrix(mtxAp);
	if(!fm)
		goto err; // FIXME
	if(want_conj)
		errval = rsb__util_do_conjugate(coo.VA,coo.typecode,coo.nnz);
	if(RSB_SOME_ERROR(errval))goto err;
	RSB_SWAP(rsb_coo_idx_t*,coo.IA,coo.JA);
	RSB_SWAP(rsb_coo_idx_t,coo.nr,coo.nc);
	RSB_COO_MEMCPY_parallel(fm->VA,fm->bpntr,fm->bindx,coo.VA,coo.IA,coo.JA,0,0,coo.nnz,fm->el_size);
 	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_SORTED_INPUT);
	tmatrix = rsb__mtx_alloc_inner(fm->VA,fm->bpntr,fm->bindx,coo.nnz,0,0,coo.typecode,coo.nr,coo.nc,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,mtxAp->flags,&errval);
 	RSB_DO_FLAG_ADD(mtxAp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
	rsb__destroy_inner(mtxAp);
	rsb__destroy_coo_matrix_t(&coo);
	*mtxApp = tmatrix;
	
	//RSB_ERROR("FIXME: should implement tranpose functionality !");
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_elements(const struct rsb_mtx_t * mtxAp, void * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags)
{
	rsb_coo_idx_t k;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t ifo = ( flags & RSB_FLAG_FORTRAN_INDICES_INTERFACE )?1:0;

	if(!IA || !VA || !JA || !mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	for(k=0;k<nnz;++k)
		errval |= rsb__do_get_coo_element(mtxAp,((rsb_char_t*)VA)+mtxAp->el_size*k,IA[k]-ifo,JA[k]-ifo);
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_set_initopt_as_string(const rsb_char_t *opn, const rsb_char_t *arg)
{
	/* FIXME: document me */
	return rsb__stropts_set(opn,arg);
}

rsb_err_t rsb__do_lib_get_info_str(int what, rsb_char_t* sbuf, size_t buflen)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t rl = buflen;
	
	if( sbuf == NULL )
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	sbuf[0] = RSB_NUL;
#ifdef RSB_CC    
	snprintf(sbuf+(buflen-rl),rl,"CC=%s ",    RSB_CC    );
#else /* RSB_CC */
	errval |= RSB_ERR_GENERIC_ERROR;
#endif /* RSB_CC */
	rl -= strlen(sbuf);
#ifdef RSB_CFLAGS
	snprintf(sbuf+(buflen-rl),rl,"CFLAGS=%s",RSB_CFLAGS);
#else /* RSB_CFLAGS */
	errval |= RSB_ERR_GENERIC_ERROR;
#endif /* RSB_CFLAGS */
err:
	return errval;
}

struct rsb_mtx_t * rsb__do_mtx_alloc_from_coo_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp)

{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL;
	rsb_coo_idx_t *IA_ = NULL, *JA_ = NULL;
	struct rsb_mtx_t * mtxAp = NULL;

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{
		errval = /*RSB_ERR_BADARGS|*/RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS;
		RSB_PERR_GOTO(err,RSB_ERRM_CNHEAF);
	}

	RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(flags);

	if(nnzA>0)
	{
		rsb_coo_idx_t offi = 0;

		if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
			offi = 1, RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
		errval = rsb__util_coo_alloc_copy_and_stats(&VA_,&IA_,&JA_,VA,IA,JA,nrA?NULL:&nrA,ncA?NULL:&ncA,nnzA,0,typecode,offi,0,RSB_FLAG_NOFLAGS,&flags);

		if(!VA_ || !IA_ || !JA_)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_E_VIJ);
		}
	}
	else
	{
#if !RSB_ALLOW_EMPTY_MATRICES
		/* FIXUP CASE FOR 0-NNZ MATRICES AND IMPLICIT DIAGONAL */
		if(RSB_INVALID_NNZ_COUNT_FOR_FLAGS(nnzA,flags))
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err,RSB_ERRM_CBAEM);
		}
#endif /* RSB_ALLOW_EMPTY_MATRICES */
	}
	RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(flags);

	mtxAp = rsb__mtx_alloc_inner(VA_,IA_,JA_,nnzA,0,0,typecode,nrA,ncA,brA,bcA,flags,&errval);
	if(mtxAp && errval == RSB_ERR_NO_ERROR)
		goto ok;
	/* FIXME: and if !matrix but errval ? */
err:
	RSB_CONDITIONAL_FREE(IA_);
	RSB_CONDITIONAL_FREE(JA_);
	RSB_CONDITIONAL_FREE(VA_);
ok:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

struct rsb_mtx_t * rsb__do_mtx_alloc_from_coo_inplace(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_coo_idx_t roff = 0,coff = 0;

	RSB_ASSERT(VA || nnzA == 0 );
	RSB_ASSERT(IA || nnzA == 0 );
	RSB_ASSERT(JA || nnzA == 0 );

	if(RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}
	RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(flags);
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		roff = -1,coff = -1, RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
	mtxAp = rsb__mtx_alloc_inner(VA,IA,JA,nnzA,roff,coff,typecode,nrA,ncA,brA,bcA,flags,&errval);
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

struct rsb_mtx_t * rsb__do_mtx_alloc_from_csr_const(const void *VA, const rsb_coo_idx_t * RP, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL;
	rsb_coo_idx_t *IA_ = NULL,*JA_ = NULL;
	struct rsb_mtx_t * mtxAp = NULL;
	size_t cnnz = RSB_MAX(nnzA,nrA+1), cis = sizeof(rsb_coo_idx_t),nvs = RSB_SIZEOF(typecode);
	rsb_coo_idx_t roff = 0,coff = 0;

	RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(flags);

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{
		errval = /*RSB_ERR_BADARGS|*/RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS;
		RSB_PERR_GOTO(err,RSB_ERRM_BFEANS);
	}

	IA_ = rsb__clone_area_with_extra(RP,cis*(nrA+1),0,cis*(cnnz-(nrA+1)));
	JA_ = rsb__clone_area_with_extra(JA,cis*(nnzA ),0,cis*(cnnz-nnzA));
	VA_ = rsb__clone_area_with_extra(VA,nvs*(nnzA ),0,nvs*(cnnz-nnzA));

	if(!VA_ || !IA_ || !JA_)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_E_VIJ);
	}
	errval = rsb__util_uncompress_row_pointers_array((rsb_coo_idx_t*)RP,nrA,flags,RSB_FLAG_C_INDICES_INTERFACE,IA_);

	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}
	if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		coff = -1, RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);

	RSB_DEBUG_ASSERT(roff>=-1 && coff>=-1); /* for Fortran */
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT);
	mtxAp = rsb__mtx_alloc_inner(VA_,IA_,JA_,nnzA,roff,coff,typecode,nrA,ncA,brA,bcA,flags,&errval);
err:
	if( RSB_SOME_ERROR(errval) && ( VA_ || IA_ || JA_) )
	{
		RSB_CONDITIONAL_FREE(IA_);
		RSB_CONDITIONAL_FREE(JA_);
		RSB_CONDITIONAL_FREE(VA_);
	}
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

struct rsb_mtx_t * rsb__do_mtx_alloc_from_csc_const(const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * CP, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t * errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL;
	rsb_coo_idx_t *IA_ = NULL,*JA_ = NULL;
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_nnz_idx_t maxdim = RSB_MAX(nnzA,RSB_MAX(nrA+1,ncA+1));

	RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(flags);

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{
		errval = /*RSB_ERR_BADARGS|*/RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}

	if(nnzA>0)
	{
		rsb_time_t dt;
		rsb_coo_idx_t offi = 0;

		if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc(&VA_,&IA_,&JA_,maxdim,typecode,RSB_BOOL_FALSE)))
		{
			RSB_PERR_GOTO(err,RSB_ERRM_EM);
		}
		if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE))
			offi = 1, RSB_DO_FLAG_DEL(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
		//errval=
		dt = - rsb_time();
		rsb__util_csc2csr(VA,IA,CP,VA_,IA_,JA_,nrA,ncA,nnzA,typecode,offi,0,&flags);/* FIXME: assembly shall give the user chance to pass offo and offi */
		dt += rsb_time();
		//printf("csc 2 csr took %lg s\n",dt);
	}
	else
	{
	}
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT);
	mtxAp = rsb_mtx_alloc_from_csr_inplace (VA_,IA_,JA_,nnzA,typecode,nrA,ncA,brA,bcA,flags,errvalp);
	if(mtxAp)
		RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS);
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_spata(const void *alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY)
{
	/* FIXME: untested; details to be finished ! */
	rsb_err_t errval = RSB_ERR_UNIMPLEMENTED_YET;
	rsb_byte_t *Zp = NULL;

	rsb_time_t mvndt = RSB_TIME_ZERO, mvtdt = RSB_TIME_ZERO, atadt = RSB_TIME_ZERO;

if(getenv("WANT_FAKE_ATA")) /* FIXME: temporary */
{
	Zp = rsb__calloc(RSB_SIZEOF(mtxAp->typecode)*mtxAp->nr);
	if(!Zp){errval = RSB_ERR_ENOMEM; RSB_PERR_GOTO(err,RSB_ERRM_EM); }
	if(!getenv("WANT_ATA_SKIP_SPMV")) /* FIXME: temporary */
	{
		mvndt = -rsb_time();
		errval = rsb_spmv(RSB_TRANSPOSITION_N,alphap,mtxAp,Xp,incX,/*betap*/NULL,Zp,incY);
		mvndt += rsb_time();
	}
	if(RSB_SOME_ERROR(errval)) goto err;
	if(!getenv("WANT_ATA_SKIP_SPMVT")) /* FIXME: temporary */
	{
		mvtdt = -rsb_time();
		errval = rsb_spmv(RSB_TRANSPOSITION_T,alphap,mtxAp,Zp,incX,/*betap*/NULL,Yp,incY);
		mvtdt += rsb_time();
	}
	if(RSB_SOME_ERROR(errval)) goto err;
	errval = RSB_ERR_NO_ERROR;
	if(RSB_SOME_ERROR(errval)) goto err;
}
	else
{
	struct rsb_translated_matrix_t * all_leaf_matrices=NULL;
	struct rsb_spmv_lock_struct_t zlock,ylock;
	rsb_trans_t transA = RSB_TRANSPOSITION_N;
	rsb_submatrix_idx_t all_leaf_matrices_n=0;

	atadt = -rsb_time();

	Zp = rsb__calloc(RSB_SIZEOF(mtxAp->typecode)*mtxAp->nr);

	if(!Zp){errval = RSB_ERR_ENOMEM; RSB_PERR_GOTO(err,RSB_ERRM_EM); }
	errval = rsb__do_get_submatrices_block_for_get_csr(mtxAp,&all_leaf_matrices,&all_leaf_matrices_n);
	if(RSB_SOME_ERROR(errval)) goto err;
	errval = rsb__sort_array_of_leaf_matrices(NULL,all_leaf_matrices, all_leaf_matrices_n, rsb_op_ata );
	if(RSB_SOME_ERROR(errval)) goto err;
	errval = rsb_do_spmv_lock_init(&zlock,rsb_global_session_handle.rsb_want_threads,all_leaf_matrices_n,mtxAp,RSB_OP_FLAG_DEFAULT,transA,Zp,incy);
	if(RSB_SOME_ERROR(errval)) goto err;
	errval = rsb_do_spmv_lock_init(&ylock,rsb_global_session_handle.rsb_want_threads,all_leaf_matrices_n,mtxAp,RSB_OP_FLAG_DEFAULT,rsb__do_transpose_transposition(transA),y,incy);
	if(RSB_SOME_ERROR(errval)) goto err;

	#pragma omp parallel reduction(|:errval) shared(zlock,all_leaf_matrices,mtxAp)  RSB_NTC 
{
	rsb_thr_t th_id = rsb_get_thread_num();
	rsb_submatrix_idx_t n=0;
	rsb_submatrix_idx_t zdm=0;
	rsb_submatrix_idx_t ydm=0;

	if(th_id >= rsb_global_session_handle.rsb_want_threads)
		goto skip;

	if(th_id>=all_leaf_matrices_n)
		goto skip;
again:
	for(n=0;RSB_LIKELY(n<all_leaf_matrices_n);++n)
	{
		const struct rsb_mtx_t *submatrix=all_leaf_matrices[n].mtxlp;
		int nort = 0;
		rsb_coo_idx_t oincy=incY;
		int nrhs = 1;
		const size_t outnri=mtxAp->nr,rhsnri=mtxAp->nr;
		char *ov=Yp; /* FIXME: unused */

		nort = 0;
		#pragma omp critical (rsb_sp_ata_crs)
		{
			if(!RSB_BITVECTOR_GET(zlock.bmap,zlock.subms,n))
			{
		       		nort=(rsb_do_spmv_lock_get(&zlock,th_id,submatrix->broff,submatrix->bm,submatrix->bcoff,submatrix->bk,n,transA,&ov,&oincy)==RSB_BOOL_TRUE)?1:0;
			}
			else
			if(!RSB_BITVECTOR_GET(ylock.bmap,ylock.subms,n))
			{
				rsb_submatrix_idx_t in;

				for(in=0;in<n;++in)
					/* if(!RSB_BITVECTOR_GET(zlock.bmap,zlock.subms,in)) */
					if( all_leaf_matrices[n  ].mtxlp->roff < all_leaf_matrices[in].mtxlp->roff+all_leaf_matrices[in].mtxlp->bm && (!RSB_BITVECTOR_GET(zlock.bmap,zlock.subms,in))   )
						goto noreq;
				for(in=n+1;in<all_leaf_matrices_n;++in)
					if( all_leaf_matrices[n  ].mtxlp->roff < all_leaf_matrices[in].mtxlp->roff+all_leaf_matrices[in].mtxlp->bm && (!RSB_BITVECTOR_GET(zlock.bmap,zlock.subms,in))   )
						goto noreq;
		       		nort=(rsb_do_spmv_lock_get(&ylock,th_id,submatrix->broff,submatrix->bm,submatrix->bcoff,submatrix->bk,n,rsb__do_transpose_transposition(transA),&ov,&oincy)==RSB_BOOL_TRUE)?2:0;
noreq:	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
			}
	       	}

		if(nort>0)
		{
			const char * offx=NULL; char *offy=NULL;
			const size_t scoff=submatrix->coff-mtxAp->coff;
			const size_t sroff=submatrix->roff-mtxAp->roff;

			if(nort==1)
if(!getenv("WANT_ATA_SKIP_SPMV")) /* FIXME: temporary */
			{
				offy=((char*)Zp)+(mtxAp->el_size*sroff)*oincy,offx=((const char*)Xp)+(mtxAp->el_size*scoff)*incX;
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_non_recursive(submatrix,offx,offy,alphap,/*betap*/NULL,incX,oincy,(nort == 1 ? transA : rsb__do_transpose_transposition(transA) ) RSB_INNER_NRHS_SPMV_ARGS_IDS));
			}

			if(nort==2)
if(!getenv("WANT_ATA_SKIP_SPMVT")) /* FIXME: temporary */
			{
				offy=((char*)Yp)+(mtxAp->el_size*sroff)*oincy,offx=((const char*)Zp)+(mtxAp->el_size*scoff)*incX;
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_non_recursive(submatrix,offx,offy,alphap,/*betap*/NULL,incX,oincy,(nort == 1 ? transA : rsb__do_transpose_transposition(transA) ) RSB_INNER_NRHS_SPMV_ARGS_IDS));
			}

                     	#pragma omp critical (rsb_sp_ata_crs)
			{
				if(nort == 1)
				{ rsb_do_spmv_lock_release(&zlock,th_id,ov);RSB_DO_SPMV_LOCK_DM_INC(zlock); }
				if(nort == 2)
				{ rsb_do_spmv_lock_release(&ylock,th_id,ov);RSB_DO_SPMV_LOCK_DM_INC(ylock); }
			}
		}
		nort = 0;
	}
		#pragma omp critical (rsb_spmv_crs)
		{ zdm = RSB_DO_SPMV_LOCK_DM(zlock); ydm = RSB_DO_SPMV_LOCK_DM(ylock); }
		if(ydm<all_leaf_matrices_n || zdm<all_leaf_matrices_n)
			goto again;
skip:	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
}
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_spmv_lock_free(&zlock));
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_spmv_lock_free(&ylock));
	RSB_CONDITIONAL_FREE(all_leaf_matrices);
err:	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	errval = RSB_ERR_NO_ERROR;
	atadt += rsb_time();
}

if(getenv("WANT_ATAV")) /* FIXME: temporary */
{
	printf("ATA times: TOT:%6.3les ATA:%6.3les MVN:%6.3les MVT:%6.3les\n",atadt+mvndt+mvtdt,atadt,mvndt,mvtdt);
}
	RSB_CONDITIONAL_FREE(Zp);
	return errval;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

void rsb__debug_print_flags(rsb_flags_t flags)
{
	/* easily callable during debugging */
	RSB_STDOUT(RSB__FLAGS_PRINTF_ARGS(flags));
}

rsb_err_t rsb__print_matrix_stats(const struct rsb_mtx_t * mtxAp)
{
	/* easily callable during debugging */
	return rsb__do_print_matrix_stats(mtxAp, 
		RSB_CONST_DUMP_MIF_MTX_INFO|
		RSB_CONST_DUMP_RECURSION_BRIEF|
		RSB_CONST_DUMP_RECURSION|
		RSB_CONST_DUMP_FLAGS
		,NULL);
}
/* @endcond */
