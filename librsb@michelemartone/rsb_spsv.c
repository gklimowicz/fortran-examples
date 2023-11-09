/*                                                                                                                            

Copyright (C) 2008-2020 Michele Martone

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
 * This source file contains functions for sparse recursive multicore triangular solve.
 */
/*
 * FIXME: the submatrices sorting routines are buggy.
 * */
#include "rsb_internals.h"		/* */
#include "rsb_lock.h"		/* */
#include "rsb_spsv.h"		/* */

#define RSB_WANT_VERBOSE_SPSV	0
/* #define RSB_CBLAS_X_SCAL_SPSV rsb__cblas_Xscal */
#define RSB_CBLAS_X_SCAL_SPSV rsb__cblas_Xscal_parallel

RSB_INTERNALS_COMMON_HEAD_DECLS

static rsb_err_t rsb_do_spsv_terminal(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transl RSB_INNER_NRHS_SPSV_ARGS)
{
	/**
	  	\ingroup gr_internals
		Entry function for SPSV.
		alphap can be NULL.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp || !y || !x || transl == RSB_INVALID_FLAGS || !rsb__is_square(mtxAp) || !RSB_IS_VALID_INCX_VALUE(incx) || !RSB_IS_VALID_INCX_VALUE(incy))
	{
		errval = RSB_ERR_BADARGS;
		goto ret;
	}

	/*
		FIXME : should handle alphap in a more specialized fashion.
	*/

#if 0
	if(betap && !RSB_IS_ELEMENT_ONE(betap,mtxAp->typecode))
	{
		if(incy>1)
			rsb__cblas_Xscal(mtxAp->typecode,rsb__do_get_columns_of(mtxAp,transl),betap,y,incy);
		else
		{
			/* if should zero the output vector */
			if(RSB_IS_ELEMENT_ZERO(betap,mtxAp->typecode))
				rsb__cblas_Xscal(mtxAp->typecode,rsb__do_get_columns_of(mtxAp,transl),NULL,y,incy);
			else
			/* if should scale the output vector */
			if(!RSB_IS_ELEMENT_ONE(betap,mtxAp->typecode))
				rsb_vector_scale(y,betap,mtxAp->typecode,rsb__do_get_columns_of(mtxAp,transl));
		}
	}
#endif


#if RSB_ENABLE_INNER_NRHS_SPSV
	{
	const size_t lenx=(mtxAp->el_size*rhsnri);
	const size_t leny=(mtxAp->el_size*outnri);
	rsb_int_t nrhsi=0;
	for(nrhsi=0;nrhsi<nrhs;++nrhsi)
#endif
	{
#if RSB_ENABLE_INNER_NRHS_SPSV
		void      *out=((      rsb_byte_t*)y)+(leny*nrhsi);
		const void*rhs=((const rsb_byte_t*)x)+(lenx*nrhsi);
#else
		void      *out=((      rsb_byte_t*)y);
		const void*rhs=((const rsb_byte_t*)x);
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
	if(!alphap || RSB_IS_ELEMENT_ONE(alphap,mtxAp->typecode))
	{
		if(incy==1 && incx==1)
			errval = rsb__do_spsv_uxua(mtxAp,rhs,out,transl);
		else
			errval = rsb__do_spsv_sxsx(mtxAp,rhs,y,alphap,incx,incy,transl);
	}
	else
		errval = rsb__do_spsv_sxsx(mtxAp,rhs,out,alphap,incx,incy,transl);
	}
#if RSB_ENABLE_INNER_NRHS_SPSV
	}
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
ret:
	return errval;
}

static rsb_err_t rsb_do_spsv_recursive_serial(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transl, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPSV_ARGS)
{
	/**
	  	\ingroup gr_internals
	 *
	 *	FIXME : document
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_aligned_t mone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t pone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_int_t nrhsi = 0;
	rsb__util_set_area_to_converted_integer(&mone[0],mtxAp->typecode,-1);
	rsb__util_set_area_to_converted_integer(&pone[0],mtxAp->typecode,+1);

	if(mtxAp->roff == mtxAp->coff)
	if(rsb__is_root_matrix(mtxAp))
#if RSB_ENABLE_INNER_NRHS_SPSV
	for (nrhsi=0;nrhsi<nrhs;++nrhsi)
#endif /* RSB_ENABLE_INNER_NRHS_SPSV */
		RSB_CBLAS_X_SCAL_SPSV(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transl),alphap,RSB_TYPED_OFF_PTR(mtxAp->typecode,y,nrhsi*(outnri)*incy),incy);

	if( rsb__is_recursive_matrix(mtxAp->flags))
#if 1
	{
		void*offy=NULL; const void *offx=NULL;
//		rsb_coo_idx_t scoff=submatrix->coff; rsb_coo_idx_t sroff=submatrix->roff;
		rsb_bool_t isupp = rsb__is_upper_triangle(mtxAp->flags);
		rsb_bool_t istrans = (RSB_DOES_NOT_TRANSPOSE(transl))?0:1;
		rsb_coo_idx_t half;
//		offy=((rsb_byte_t*)y)+(mtxAp->el_size*sroff)*incy,offx=((const rsb_byte_t*)x)+(mtxAp->el_size*scoff)*incx;
		if(mtxAp->coff!=mtxAp->roff)
		{	
			RSB_ERROR("!\n");
			errval = RSB_ERR_BADARGS;goto err;
		}
		if(!(RSB_SUBMATRIX_INDEX(mtxAp,0,0)) || !(RSB_SUBMATRIX_INDEX(mtxAp,1,1)))
		{
			RSB_ERROR("@ %d %d and with null diagonal elements, %p %p %p %p\n",mtxAp->roff,mtxAp->coff,
					mtxAp->sm[0],
					mtxAp->sm[1],
					mtxAp->sm[2],
					mtxAp->sm[3]
					);
			errval = RSB_ERR_BADARGS;goto err;
		}
		half = RSB_SUBMATRIX_INDEX(mtxAp,1,1)->roff-mtxAp->roff;
		offy=((      rsb_byte_t*)y)+(mtxAp->el_size*half)*incy;
		offx=((const rsb_byte_t*)x)+(mtxAp->el_size*half)*incx;
	
		switch(isupp)
		{
		case RSB_BOOL_TRUE:
		switch(istrans)
		{
		case RSB_BOOL_TRUE:
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,0,0),x,y,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			if(RSB_SUBMATRIX_INDEX(mtxAp,0,1) && op_flags != RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transl,&mone[0],RSB_SUBMATRIX_INDEX(mtxAp,0,1),offx,incx,NULL,y,incy,RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT_SERIAL RSB_INNER_NRHS_SPSV_ARGS_IDS));
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,1,1),offx,offy,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
		break;
		default:
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,1,1),offx,offy,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			if(RSB_SUBMATRIX_INDEX(mtxAp,0,1) && op_flags != RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transl,&mone[0],RSB_SUBMATRIX_INDEX(mtxAp,0,1),offx,incx,NULL,y,incy,RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT_SERIAL RSB_INNER_NRHS_SPSV_ARGS_IDS));
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,0,0),x,y,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
		break;
		}
		break;
		case RSB_BOOL_FALSE:
		switch(istrans)
		{
		case RSB_BOOL_TRUE:
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,1,1),offx,offy,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			if(RSB_SUBMATRIX_INDEX(mtxAp,1,0) && op_flags != RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transl,&mone[0],RSB_SUBMATRIX_INDEX(mtxAp,1,0),x,incx,NULL,offy,incy,RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT_SERIAL RSB_INNER_NRHS_SPSV_ARGS_IDS));
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,0,0),x,y,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
		break;
		default:
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,0,0),x,y,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			if(RSB_SUBMATRIX_INDEX(mtxAp,1,0) && op_flags != RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE)
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_general(transl,&mone[0],RSB_SUBMATRIX_INDEX(mtxAp,1,0),x,incx,NULL,offy,incy,RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT_SERIAL RSB_INNER_NRHS_SPSV_ARGS_IDS));
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(RSB_SUBMATRIX_INDEX(mtxAp,1,1),offx,offy,&pone[0],incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
		break;
		}
		break;
		default:
			errval = RSB_ERR_INTERNAL_ERROR;goto err;
		break;
		}
	}
#else
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix=NULL;
		void*offy=NULL; const void *offx=NULL;

		if( RSB_DOES_NOT_TRANSPOSE( transl ) )
		{
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
//			RSB_STDOUT("%d %d \n",i,j);
			rsb_coo_idx_t scoff=submatrix->coff; rsb_coo_idx_t sroff=submatrix->roff;

			RSB_DEBUG_ASSERT(scoff>=0);
			RSB_DEBUG_ASSERT(sroff>=0);

			//RSB_ERROR("-> 0x%p %d %d (%d) (%d)\n",submatrix,submatrix->roff,submatrix->coff,submatrix->nnz, rsb__is_recursive_matrix(mtxAp->flags));

				offy=((rsb_byte_t*)y)+(mtxAp->el_size*sroff)*incy,offx=((const rsb_byte_t*)x)+(mtxAp->el_size*scoff)*incx;
			if(submatrix->roff==submatrix->coff && (i==(RSB_DOES_NOT_TRANSPOSE(transl))?0:1) )
			{
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(submatrix,x,y,alphap,incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			}
			else
			if(submatrix->roff==submatrix->coff && (i==(RSB_DOES_NOT_TRANSPOSE(transl))?1:0) )
			{
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(submatrix,x,y,alphap,incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			}
			else
			//if(i==((RSB_DOES_NOT_TRANSPOSE(transl))?1:0))
			if(i==1 && op_flags != RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE)
			{
			//	RSB_STDOUT("offx %g offy %g\n",*(double*)offx,*(double*)offy);
//				RSB_STDOUT("spmv %d %d\n",submatrix->roff,submatrix->coff);
				/* FIXME : DOES NOT TAKE INTO ACCOUNT INCX,INCY */
				/* transposition is not relevant, as long as we work with square matrices everywhere */
			//	if(y!=x)
			//		RSB_DO_ERROR_CUMULATE(errval,rsb__xcopy(offy,x,0,0,sroff,mtxAp->el_size));
				// FIXME : the following lines should be equivalent, but they aren't . why ?
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_recursive_serial(submatrix,offx,offy,&mone[0],NULL,1,1,transl RSB_INNER_NRHS_SPMV_ARGS_IDS));
			//	RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_recursive_parallel(submatrix,offx,offy,&mone[0],NULL,1,1,transl));
				//RSB_DO_ERROR_CUMULATE(errval,rsb_spmv_unua(submatrix,offx,offy,transl));
			//	RSB_STDOUT("offx %g offy %g\n",*(double*)offx,*(double*)offy);
			}
			if(RSB_SOME_ERROR(errval))
				goto err;
		}}
		else
		{
		RSB_SUBMATRIX_FOREACH_REVERSE(mtxAp,submatrix,i,j)
		{
		//	RSB_STDOUT("%d %d \n",i,j);
		if(submatrix)
		{
//			RSB_STDOUT("%d %d \n",i,j);
			rsb_coo_idx_t scoff=submatrix->coff; rsb_coo_idx_t sroff=submatrix->roff;

			RSB_DEBUG_ASSERT(scoff>=0);
			RSB_DEBUG_ASSERT(sroff>=0);

			offy=((rsb_byte_t*)y)+(mtxAp->el_size*sroff)*incy,offx=((const rsb_byte_t*)x)+(mtxAp->el_size*scoff)*incx;
			if(submatrix->roff==submatrix->coff && (i==(RSB_DOES_NOT_TRANSPOSE(transl)))?0:1) )
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(submatrix,x,y,alphap,incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			else
			if(submatrix->roff==submatrix->coff && (i==RSB_DOES_NOT_TRANSPOSE(transl)))?1:0) )
				RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_recursive_serial(submatrix,x,y,alphap,incx,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS));
			else
			if(i==1 && op_flags != RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE)
			{
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_recursive_serial(submatrix,offx,offy,&mone[0],NULL,1,1,transl RSB_INNER_NRHS_SPMV_ARGS_IDS));
			}
			if(errval != RSB_ERR_NO_ERROR)goto err;
		}}
		}
	}
#endif
	else
	{
		void*offy=NULL; const void *offx=NULL;
		rsb_coo_idx_t scoff=0;
		rsb_coo_idx_t sroff=0;
		RSB_DEBUG_ASSERT(scoff>=0);
		RSB_DEBUG_ASSERT(sroff>=0);
		offy=((rsb_byte_t*)y)+(mtxAp->el_size*sroff)*incy,offx=((const rsb_byte_t*)x)+(mtxAp->el_size*scoff)*incx;
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_terminal(mtxAp,offx,offy,alphap,incx,incy,transl RSB_OUTER_NRHS_SPSV_ARGS_IDS));
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_submatrices_block_for_get_csr(const struct rsb_mtx_t * mtxAp, struct rsb_translated_matrix_t ** all_leaf_matricesp, rsb_submatrix_idx_t * all_leaf_matrices_np)
{
	/**	
	 * \ingroup gr_internals
	 * FIXME: rename : csr -> csc
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t all_leaf_matrices_n=0;
	struct rsb_translated_matrix_t * all_leaf_matrices=NULL;

	all_leaf_matrices_n=mtxAp->all_leaf_matrices_n;
	all_leaf_matrices = rsb__clone_area(mtxAp->all_leaf_matrices,sizeof(struct rsb_translated_matrix_t)*all_leaf_matrices_n);
	errval = rsb__sort_array_of_leaf_matrices(NULL,all_leaf_matrices,all_leaf_matrices_n,rsb_op_get_csr);

	*all_leaf_matrices_np=all_leaf_matrices_n;
	*all_leaf_matricesp=all_leaf_matrices;
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_submatrices_for_ussv( const struct rsb_mtx_t * mtxAp, struct rsb_translated_matrix_t ** all_leaf_matricesp, rsb_submatrix_idx_t * all_leaf_matrices_np, rsb_trans_t transT)
{
	/**
	  	\ingroup gr_internals
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t all_leaf_matrices_n=0;
	struct rsb_translated_matrix_t * all_leaf_matrices=NULL;

	all_leaf_matrices_n=mtxAp->all_leaf_matrices_n;
	all_leaf_matrices = rsb__clone_area(mtxAp->all_leaf_matrices,sizeof(struct rsb_translated_matrix_t)*all_leaf_matrices_n);
	rsb__submatrices_exclude_nontriangular(all_leaf_matrices,&all_leaf_matrices_n,mtxAp);
	errval = rsb__sort_array_of_leaf_matrices_for_ussv(mtxAp,all_leaf_matrices,all_leaf_matrices_n,transT);

	*all_leaf_matrices_np=all_leaf_matrices_n;
	*all_leaf_matricesp=all_leaf_matrices;
	RSB_DO_ERR_RETURN(errval)
}

void rsb__submatrices_exclude_nontriangular(struct rsb_translated_matrix_t * all_leaf_matrices, rsb_submatrix_idx_t * all_leaf_matrices_np, const struct rsb_mtx_t * mtxAp)
{
	/**
	  	\ingroup gr_internals
	*/
	rsb_submatrix_idx_t n, all_leaf_matrices_n=0;
	RSB_DEBUG_ASSERT(mtxAp);
	RSB_DEBUG_ASSERT(all_leaf_matrices);
	RSB_DEBUG_ASSERT(all_leaf_matrices_np);
	if(rsb__is_upper_triangle(mtxAp->flags))
	{
		for(n=0;n<mtxAp->all_leaf_matrices_n;++n)
			if (mtxAp->all_leaf_matrices[n].roff<=mtxAp->all_leaf_matrices[n].coff)
				all_leaf_matrices[all_leaf_matrices_n++]=mtxAp->all_leaf_matrices[n];
	}
	else
	{
		for(n=0;n<mtxAp->all_leaf_matrices_n;++n)
			if (mtxAp->all_leaf_matrices[n].roff>=mtxAp->all_leaf_matrices[n].coff)
				all_leaf_matrices[all_leaf_matrices_n++]=mtxAp->all_leaf_matrices[n];
	}
//	if(all_leaf_matrices_n<mtxAp->all_leaf_matrices_n)
//	;
	//;RSB_STDOUT("FIX : discarded %d upper diagonal matrices out of %d \n",mtxAp->all_leaf_matrices_n-all_leaf_matrices_n,mtxAp->all_leaf_matrices_n);
	*all_leaf_matrices_np=all_leaf_matrices_n;
}

#if RSB_WANT_OMP_RECURSIVE_KERNELS
static rsb_err_t rsb__do_spsv_uxua_recursive_parallel(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transl, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPSV_ARGS)
{
	/**
	  	\ingroup gr_internals
	 * triangular solve for recursive
	 *
	 * each submatrix starting at (i,j) depends on :
	   
	   i if i=j
	   
	   1 2 3 4 5 6
	  +-+---------+
	  +-+-+       | 1
	  +-+-+-+     | 2
	  +-+-+-+-+   | 3
	  +-+-+-+-+-+ | 4
	  +-+-+-+-+-+-+ 5
	  +-+-+-+-+-+-+ 6

	  The active rows in the spmv are locked with an interval information.
	  Since the rows interval active in trsv is not under spmv, there is no need for another lock.

	  We need a ...

	  alphap NULL means alpha = 1

	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
#if	RSB_WANT_VERBOSE_SPSV
	rsb_time_t sv_time = RSB_CONST_IMPOSSIBLY_BIG_TIME,dp_time = RSB_CONST_IMPOSSIBLY_BIG_TIME;
#endif /* RSB_WANT_VERBOSE_SPSV */
	struct rsb_translated_matrix_t * all_leaf_matrices=NULL;	/** NEW, EXPERIMENTAL */
	struct rsb_rows_lock_struct_t lock;
	rsb_submatrix_idx_t * deps=NULL;	/** NEW, EXPERIMENTAL */
	rsb_submatrix_idx_t all_leaf_matrices_n=0;
	rsb_bool_t backdeps,isupptri;
	rsb_aligned_t alpha_inv[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t mone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb_aligned_t pone[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
	rsb__util_set_area_to_converted_integer(&mone[0],mtxAp->typecode,-1);
	rsb__util_set_area_to_converted_integer(&pone[0],mtxAp->typecode,+1);

	RSB_BZERO_P(&lock);
	if( !mtxAp)
	{errval = RSB_ERR_GENERIC_ERROR;goto err;}
	if( !rsb__is_recursive_matrix(mtxAp->flags) || !mtxAp)
	{errval = RSB_ERR_GENERIC_ERROR;goto err;}
	if(x!=y)/* FIXME */
	{errval = RSB_ERR_GENERIC_ERROR;goto err;}
	
	isupptri = rsb__is_upper_triangle(mtxAp->flags);
	backdeps = RSB_BOOL_XOR(RSB_DOES_TRANSPOSE(transl),isupptri);
	rsb__util_set_area_to_negated_fraction(alpha_inv,alphap,mtxAp->typecode);

#if	RSB_WANT_VERBOSE_SPSV
	dp_time = rsb_time();
#endif /* RSB_WANT_VERBOSE_SPSV */

	errval = rsb__do_get_submatrices_for_ussv(mtxAp,&all_leaf_matrices,&all_leaf_matrices_n,transl);
	deps = rsb__malloc(sizeof(rsb_submatrix_idx_t)*all_leaf_matrices_n);
	if(RSB_SOME_ERROR(errval) || !all_leaf_matrices || !deps)
	{errval = RSB_ERR_ENOMEM;goto err;}

#if 0
	{
		/* printout */
		rsb_submatrix_idx_t n;
			for(n=0;n<all_leaf_matrices_n;++n)
				RSB_STDOUT("got %d/%d:%d~%d,%d~%d\n",
					n,all_leaf_matrices_n,
					all_leaf_matrices[n].roff,all_leaf_matrices[n].mtxlp->nr+all_leaf_matrices[n].roff,
					all_leaf_matrices[n].coff,all_leaf_matrices[n].mtxlp->nc+all_leaf_matrices[n].coff);
	}
#endif
#if 0
	{
	rsb_submatrix_idx_t n;
	rsb_coo_idx_t s=0;
	for(n=0;n<all_leaf_matrices_n;++n)
		if(all_leaf_matrices[n].roff==all_leaf_matrices[n].coff)
		{
			s+=all_leaf_matrices[n].mtxlp->nr;
//			RSB_STDOUT("%d/%d [%d~%d,%d~%d] (on diag)\n",n,all_leaf_matrices_n,
//				all_leaf_matrices[n].roff, all_leaf_matrices[n].roff+all_leaf_matrices[n].mtxlp->nr-1,
//				all_leaf_matrices[n].coff, all_leaf_matrices[n].coff+all_leaf_matrices[n].mtxlp->nc-1);
		}
		else
			;
//			RSB_STDOUT("%d/%d [%d~%d,%d~%d] (not on diag)\n",n,all_leaf_matrices_n,
//				all_leaf_matrices[n].roff, all_leaf_matrices[n].roff+all_leaf_matrices[n].mtxlp->nr-1,
//				all_leaf_matrices[n].coff, all_leaf_matrices[n].coff+all_leaf_matrices[n].mtxlp->nc-1);
	if(mtxAp->nr != s)
	{	RSB_STDOUT("FATAL : sum of diagonal matrices rows %d != %d \n",s,mtxAp->nr);
		goto err;
}
	}
#endif
#if 0
	if(RSB_DOES_TRANSPOSE(transl))
	{
		rsb_submatrix_idx_t n; for(n=0;n<all_leaf_matrices_n;++n) {	
			RSB_SWAP(rsb_coo_idx_t,all_leaf_matrices[n].roff,all_leaf_matrices[n].coff);
		}
	}
#endif

	if(0)
	if(RSB_DOES_TRANSPOSE(transl))
	{
		rsb_submatrix_idx_t n;
		for(n=0;n<all_leaf_matrices_n;++n)
		{
			all_leaf_matrices[n].roff=mtxAp->nc-(all_leaf_matrices[n].coff+all_leaf_matrices[n].mtxlp->nc*1);
			all_leaf_matrices[n].coff=mtxAp->nr-(all_leaf_matrices[n].roff+all_leaf_matrices[n].mtxlp->nr*1);
			all_leaf_matrices[n].nr=all_leaf_matrices[n].mtxlp->nc;
			all_leaf_matrices[n].nc=all_leaf_matrices[n].mtxlp->nr;
		}
	}

	if(errval != RSB_ERR_NO_ERROR)
		goto err;
#if 0
	{
		/* printout */
		rsb_submatrix_idx_t n;
			for(n=0;n<all_leaf_matrices_n;++n)
				RSB_STDOUT("got %d/%d:%d~%d,%d~%d\n",
					n,all_leaf_matrices_n,
//					all_leaf_matrices[n].mtxlp->roff,all_leaf_matrices[n].mtxlp->nr+all_leaf_matrices[n].mtxlp->roff,
//					all_leaf_matrices[n].mtxlp->coff,all_leaf_matrices[n].mtxlp->nc+all_leaf_matrices[n].mtxlp->coff);
					all_leaf_matrices[n].roff,all_leaf_matrices[n].nr+all_leaf_matrices[n].roff,
					all_leaf_matrices[n].coff,all_leaf_matrices[n].nc+all_leaf_matrices[n].coff);
	}
#endif
#if 1
	{
		rsb_submatrix_idx_t n; 
		for(n=1;n<all_leaf_matrices_n;++n)
		{
			rsb_submatrix_idx_t np=n-1;
			if(all_leaf_matrices[n].roff==all_leaf_matrices[n].coff)
			{
				while(np>0 && all_leaf_matrices[np].roff!=all_leaf_matrices[np].coff)
					--np;
				//for(;np<n;++np)
				//	RSB_STDOUT("%d vs %d ? %d\n",n,np,rsb__compar_rcsr_matrix_for_spsvl(all_leaf_matrices+np,all_leaf_matrices+n));
			}
			else
			{
				while(np>0
				 && !(
					(  all_leaf_matrices[np].roff==all_leaf_matrices[np].coff ) /*&&
					( (all_leaf_matrices[np].coff+ all_leaf_matrices[np].mtxlp->nc)>=
					  (all_leaf_matrices[n ].coff+ all_leaf_matrices[n ].mtxlp->nc) )*/
					))
					--np;
			}
			deps[n]=np;
//#define RSB__TRSV_OUT__ 1
			if(RSB__TRSV_OUT__)RSB_STDOUT("dep: %d ->  %d\n",n,np);
		}
	}
#endif
#if 0
	if(RSB_DOES_TRANSPOSE(transl))
	{
		rsb_submatrix_idx_t n; for(n=0;n<all_leaf_matrices_n;++n) {	
			RSB_SWAP(rsb_coo_idx_t,all_leaf_matrices[n].roff,all_leaf_matrices[n].coff);
		}
#if 1
		for(n=0;n<all_leaf_matrices_n/2;++n) {	RSB_SWAP(struct rsb_translated_matrix_t,all_leaf_matrices[n],all_leaf_matrices[(all_leaf_matrices_n-1)-n]);}
	}
#endif
#endif
#if 0
	{
	rsb_submatrix_idx_t n;
	for(n=0;n<all_leaf_matrices_n;++n)
		if(all_leaf_matrices[n].roff<all_leaf_matrices[n].coff)
		{
			RSB_ERROR("all_leaf_matrices[n].roff<all_leaf_matrices[n].coff (%d<%d) in a lower triangular matrix!\n",
			all_leaf_matrices[n].roff,all_leaf_matrices[n].coff);
			{errval = RSB_ERR_GENERIC_ERROR;goto err;}
		}
	}
#endif
#if 0
	{
	rsb_submatrix_idx_t n,ad=0,d=0,pad=0,pda=0;
//	int cppad=0,ncppad=0;
	for(n=0;n<mtxAp->all_leaf_matrices_n;++n)
	{
		if(all_leaf_matrices[n].roff==all_leaf_matrices[n].coff)
		{
			rsb_submatrix_idx_t nn;
			rsb_submatrix_idx_t dr=all_leaf_matrices[n].roff+all_leaf_matrices[n].mtxlp->nr;
			++d;
			for(nn=n;nn<mtxAp->all_leaf_matrices_n;++nn)
			if(
				all_leaf_matrices[nn].roff >= dr &&
				all_leaf_matrices[nn].roff > all_leaf_matrices[nn].coff &&
				all_leaf_matrices[nn].coff + all_leaf_matrices[nn].mtxlp->nc <= dr &&
				1
			)
			{
				pad++;	/* parallelizable anti diagonal */
			}
			for(nn=n+1;nn<mtxAp->all_leaf_matrices_n && 
				all_leaf_matrices[nn].roff > all_leaf_matrices[nn].coff ;++nn)
				pda++;
		}
	}
		RSB_STDOUT(
			"    diagonal blocks : %d\n"
			"antidiagonal blocks : %d\n"
			"prediagonal blocks : %d\n"
	//		"antidiagonal blocks on the critical path : %d\n"
	//		"antidiagonal blocks not on the critical path : %d\n"
			,d,pad,pda//,cppad,ncppad
			);
	}
#endif

	errval = rsb__do_lock_init(&lock,rsb_global_session_handle.rsb_want_threads,all_leaf_matrices_n,mtxAp,op_flags);
	if(errval != RSB_ERR_NO_ERROR)
		goto err;


{
#if RSB_ENABLE_INNER_NRHS_SPSV
	const size_t leny=(mtxAp->el_size*outnri);
	rsb_int_t nrhsi=0;
	for(nrhsi=0;nrhsi<nrhs;++nrhsi)
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
	{
#if RSB_ENABLE_INNER_NRHS_SPSV
		void      *out=((      rsb_byte_t*)y)+(leny*nrhsi);
#else /* RSB_ENABLE_INNER_NRHS_SPMV */
		void      *out=((      rsb_byte_t*)y);
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
		RSB_CBLAS_X_SCAL_SPSV(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transl),alphap,out,incy);

#if 1
	/* corner triangle solve */
	if(all_leaf_matrices_n)
	{
		rsb_submatrix_idx_t n=0;
		void * offy=((rsb_byte_t*)out)+(mtxAp->el_size*all_leaf_matrices[n].roff)*incy;
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_terminal(all_leaf_matrices[n].mtxlp,offy,offy,&pone[0],incy,incy,transl RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS));
		rsb__do_lock_get(&lock,0,all_leaf_matrices[n].mtxlp->nr+all_leaf_matrices[n].roff,all_leaf_matrices[n].mtxlp->nr,all_leaf_matrices[n].coff,all_leaf_matrices[n].mtxlp->nc,n,transl);
		rsb__do_lock_release(&lock,0);
		lock.dm=1; /* first matrix processed */
		if(!backdeps)
			lock.dr=all_leaf_matrices[n].mtxlp->nr+all_leaf_matrices[n].roff;
		else
			lock.dr=all_leaf_matrices[n].roff;
	}
#else
	if(all_leaf_matrices_n)
	{
		if(backdeps)
			lock.dr=all_leaf_matrices[0].mtxlp->nr+all_leaf_matrices[0].roff;
		else
			lock.dr=all_leaf_matrices[0].mtxlp->nr;
	}
#endif
	}
}

#if 	RSB_WANT_VERBOSE_SPSV
	sv_time = rsb_time();
	dp_time=-dp_time+sv_time;
#endif /* RSB_WANT_VERBOSE_SPSV */
	#pragma omp parallel reduction(|:errval) shared(lock,all_leaf_matrices,mtxAp)  RSB_NTC 
	{
	const rsb_thr_t th_id = omp_get_thread_num();
	rsb_submatrix_idx_t n=0;
	rsb_submatrix_idx_t dm=0;
	#pragma omp barrier

	if(th_id >= rsb_global_session_handle.rsb_want_threads)
		goto skip;

	#pragma omp critical (rsb_spsv_crs)
	{ dm=lock.dm; }

again:
	for(n=0;n<all_leaf_matrices_n;++n)
	//if(!RSB_BITMAP_GET(lock.bmap,1,lock.subms,0,n))
	{
		struct rsb_mtx_t *submatrix=all_leaf_matrices[n].mtxlp;
		const rsb_byte_t* trhs=((rsb_byte_t*)y)+mtxAp->el_size*all_leaf_matrices[n].coff*incx;
		rsb_byte_t* tout=((rsb_byte_t*)y)+mtxAp->el_size*all_leaf_matrices[n].roff*incy;
		enum rsb_op_t op = rsb_op_nop;
		rsb_coo_idx_t roff=all_leaf_matrices[n].roff;
		rsb_coo_idx_t coff=all_leaf_matrices[n].coff;

//#define RSB__TRSV_OUT__ 1
		#pragma omp critical (rsb_spsv_crs)
		{
		if( th_id==0 )
		{
			if(	
			 	(
				 	(( backdeps) && coff+all_leaf_matrices[n].mtxlp->nc==lock.dr ) ||
				 	((!backdeps) && roff==lock.dr ) 
				)
				&& roff==coff  && 
					(!RSB_BITMAP_GET(lock.bmap,1,lock.subms,0,n))
				)
				{
					rsb_submatrix_idx_t np=deps[n];
					while(RSB_BITMAP_GET(lock.bmap,1,lock.subms,0,np))
					{
						++np;
						if(RSB__TRSV_OUT__) RSB_STDOUT("%d -> %d\n",n,np);
					}
					if(np==n && (rsb__do_lock_get(&lock,th_id,roff,all_leaf_matrices[n].mtxlp->nr,coff,all_leaf_matrices[n].mtxlp->nc,n,transl)==RSB_BOOL_TRUE))
						op = rsb_op_spsvl;
				}
		}

		if(op == rsb_op_nop)
		{
//			if(RSB__TRSV_OUT)RSB_STDOUT("%d@%d %d %d %d %d\n",n,th_id,op,all_leaf_matrices[n].roff,all_leaf_matrices[n].coff,omp_get_num_threads());
			if(
			 	(((!backdeps) && 
				(
				 (roff >=lock.dr && !isupptri &&
				roff != coff &&
				coff + all_leaf_matrices[n].mtxlp->nc <= lock.dr) ||
				 (coff >=lock.dr && isupptri &&
				roff != coff &&
				roff + all_leaf_matrices[n].mtxlp->nr <= lock.dr)
				)
			       	) 
				 ||
			  	(( backdeps) && 
				(( isupptri &&(coff >= lock.dr)) ||
				((!isupptri)&&(roff >= lock.dr))) &&
			       	roff != coff //&&
				//coff + all_leaf_matrices[n].mtxlp->nc <=lock.dr 
				))
			       	&& 
#if RSB_WANT_BOUNDED_BOXES_SPSV
				(rsb__do_lock_get(&lock,th_id,all_leaf_matrices[n].mtxlp->broff,all_leaf_matrices[n].mtxlp->bm,all_leaf_matrices[n].mtxlp->bcoff,all_leaf_matrices[n].mtxlp->bk,n,transl)==RSB_BOOL_TRUE)&&
#else /* RSB_WANT_BOUNDED_BOXES_SPSV */
				(rsb__do_lock_get(&lock,th_id,all_leaf_matrices[n].roff,all_leaf_matrices[n].mtxlp->nr,all_leaf_matrices[n].coff,all_leaf_matrices[n].mtxlp->nc,n,transl)==RSB_BOOL_TRUE)&&
#endif /* RSB_WANT_BOUNDED_BOXES_SPSV */
				1
				)
				op = rsb_op_spmv ;
		}
		}
//			if(RSB__TRSV_OUT)RSB_STDOUT("%d@%d %d %d %d %d\n",n,th_id,op,all_leaf_matrices[n].roff,all_leaf_matrices[n].coff,omp_get_num_threads());
		if(RSB__TRSV_OUT__)dm=lock.dm;
		if(RSB__TRSV_OUT__)
		if(
				(!RSB_BITMAP_GET(lock.bmap,1,lock.subms,0,n) && op == rsb_op_nop)||
				( RSB_BITMAP_GET(lock.bmap,1,lock.subms,0,n) && op != rsb_op_nop)
				)
		RSB_STDOUT("%zd/%zd [%zd~%zd,%zd~%zd] on th.%zd -> op %zd (dr:%zd) (done:%zd)\n",(rsb_printf_int_t)n,(rsb_printf_int_t)all_leaf_matrices_n,
				(rsb_printf_int_t)all_leaf_matrices[n].roff, (rsb_printf_int_t)(all_leaf_matrices[n].roff+submatrix->nr-1),
				(rsb_printf_int_t)all_leaf_matrices[n].coff, (rsb_printf_int_t)(all_leaf_matrices[n].coff+submatrix->nc-1),
				(rsb_printf_int_t)th_id,(rsb_printf_int_t)op,(rsb_printf_int_t)lock.dr,(rsb_printf_int_t)dm);

		switch(op){
		case rsb_op_spsvl:
		{
			/* diagonal blocks */
			if(RSB__TRSV_OUT__)RSB_STDOUT("spsv on %zd on %zd \n",(rsb_printf_int_t)n,(rsb_printf_int_t)th_id);
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_terminal(submatrix,trhs,tout,&pone[0],incx,incy,transl RSB_OUTER_NRHS_SPSV_ARGS_IDS));
			//RSB_DO_ERROR_CUMULATE(errval,rsb_do_spsv_terminal(submatrix,trhs,tout,alphap,incx,incy,transl RSB_OUTER_NRHS_SPSV_ARGS_IDS));
                       	#pragma omp critical (rsb_spsv_crs)
			{
				if(!backdeps)
			       		lock.dr=submatrix->nr+roff;
				else
			       		lock.dr=submatrix->roff;
				rsb__do_lock_release(&lock,th_id); ++lock.dm; dm=lock.dm;
			}
		}
		break;
		case rsb_op_spmv:
		{
			/* antidiagonal blocks */
			if(RSB__TRSV_OUT__)RSB_STDOUT("spmv on %zd on %zd \n",(rsb_printf_int_t)n,(rsb_printf_int_t)th_id);
			//RSB_DO_ERROR_CUMULATE(errval,rsb_spmv_unua(submatrix,trhs,tout,transl));
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_non_recursive(submatrix,trhs,tout,&mone[0],NULL,incx,incy,transl RSB_OUTER_NRHS_SPSV_ARGS_IDS));
                       	#pragma omp critical (rsb_spsv_crs)
			{rsb__do_lock_release(&lock,th_id);++lock.dm;dm=lock.dm;}
		}
		break;
		case rsb_op_nop:
		{
		}
		break;
		default: //rsb_op_spsvlt, rsb_op_spsvu, rsb_op_spsvut, rsb_op_get_csr
			RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
		}
		//if(errval != RSB_ERR_NO_ERROR)
		//	break;
		//if(op != rsb_op_nop)
	
	}
	
	#pragma omp critical (rsb_spsv_crs)
	{ dm=lock.dm; }
	if(RSB__TRSV_OUT__)RSB_STDOUT("on thread %d : done %d/%d \n",th_id,lock.dm,all_leaf_matrices_n);
	if(dm<all_leaf_matrices_n
#if RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPSV
			&& ((all_leaf_matrices_n-dm)>th_id)
#endif /* RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPSV */
	)goto again;
skip:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
		/* FIXME : could place a barrier here. */
		#pragma omp barrier
		/* now we can leave the parallel region safely */
	}
err:
	RSB_CONDITIONAL_FREE(all_leaf_matrices);
	RSB_CONDITIONAL_FREE(deps);
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_lock_free(&lock));
#if	RSB_WANT_VERBOSE_SPSV
	sv_time=-sv_time+rsb_time();
	RSB_INFO("SPSV: solve time:%lg   deps time:%lg\n",sv_time,dp_time);
#endif /* RSB_WANT_VERBOSE_SPSV */
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	errval = RSB_ERR_UNSUPPORTED_OPERATION;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */

rsb_err_t rsb__do_spsv_general(rsb_trans_t transl, const void * alphap, const struct rsb_mtx_t * mtxAp, const void * x, rsb_coo_idx_t incx, void * y, rsb_coo_idx_t incy, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPSV_ARGS)
{
	/**
	  	\ingroup gr_internals
	 	computes \f$y \leftarrow \alpha op(A)^{-1} \cdot y \f$

		Entry function for SPSV.
		alphap can be NULL.

		FIXME : incx and incy check should be stricter
		FIXME : x!=y are disabled
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp || !y || !x || transl == RSB_INVALID_FLAGS || !rsb__is_square(mtxAp) || incx<1 || incy<1)
	{
		errval = RSB_ERR_BADARGS;goto err;
	}
	if(!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_TRIANGULAR))
	{
		errval = RSB_ERR_BADARGS;goto err;
	}
	if(RSB__FLAG_HAS_UNSPECIFIED_TRIANGLE(mtxAp->flags))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_BOH_TRI);
	}
#if 1
	if(x!=y)
	{
			// FIXME: should parallelize this
		rsb_int_t nrhsi = 0;
#if RSB_ENABLE_INNER_NRHS_SPSV
		for (nrhsi=0;nrhsi<nrhs;++nrhsi)
#endif /* RSB_ENABLE_INNER_NRHS_SPSV */
		{
			//RSB_DO_ERROR_CUMULATE(errval,rsb__xcopy_strided_typed(y,x,0,0,mtxAp->nr,mtxAp->typecode,incy,incx));
			RSB_DO_ERROR_CUMULATE(errval,rsb__cblas_Xcopy(mtxAp->typecode,mtxAp->nr,RSB_TYPED_OFF_PTR(mtxAp->typecode,x,nrhsi*(rhsnri)*incx),incx,RSB_TYPED_OFF_PTR(mtxAp->typecode,y,nrhsi*(outnri)*incy),incy));
		}
	}
#endif

	if(op_flags == RSB_OP_FLAG_INFINITE_PARALLELISM_EMULATE)
	{
		errval = rsb_do_spsv_recursive_serial(mtxAp,y,y,alphap,incy,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS);
		goto done;
	}
	if( !rsb__is_recursive_matrix(mtxAp->flags))
		errval = rsb_do_spsv_terminal(mtxAp,y,y,alphap,incy,incy,transl RSB_OUTER_NRHS_SPSV_ARGS_IDS);
	else
	{
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		if(op_flags == RSB_OP_FLAG_WANT_SERIAL)
			errval = rsb_do_spsv_recursive_serial(mtxAp,y,y,alphap,incy,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS);
		else
			errval = rsb__do_spsv_uxua_recursive_parallel(mtxAp,y,y,alphap,incy,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS);
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		errval = rsb_do_spsv_recursive_serial(mtxAp,y,y,alphap,incy,incy,transl,op_flags RSB_INNER_NRHS_SPSV_ARGS_IDS);
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	}
	goto done;
done:
#if 0
                {
			/* FIXME : won't work with when incx or incy is not 1 */
			rsb_aligned_t checksum[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
			RSB_CBLAS_X_SCAL_SPSV(mtxAp->typecode,1,NULL,checksum,1);
			rsb_nnz_idx_t n;
			rsb__util_vector_sum(checksum,y,mtxAp->typecode,mtxAp->nr);
			RSB_STDOUT("#spsv checksum:\n");
			rsb__debug_print_value(checksum,typecode);
			RSB_STDOUT("\n");
                }
#endif
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_spsv(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, const void * Xp, rsb_coo_idx_t incX, void * Yp, rsb_coo_idx_t incY)
{
	enum rsb_op_flags_t op_flags = RSB_OP_FLAG_DEFAULT;
#if RSB_ALLOW_INTERNAL_GETENVS
	op_flags = rsb__getenv_int_t("RSB_SPSV_OP_FLAG", op_flags);
#endif /* RSB_ALLOW_INTERNAL_GETENVS*/
	return rsb__do_spsv_general(transT,alphap,mtxTp,Xp,incX,Yp,incY,op_flags RSB_DEFAULT_OUTER_NRHS_SPSV_ARGS	);
}

/* @endcond */
