/* @cond INNERDOC */
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block coordinates format.
 Kernels unrolled, with no loops, for only user-specified blockings.
 */

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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains functions for sparse recursive multicore matrix vector multiplication.
 */
/*
 * FIXME: many beta-related ops are NOT parallel, and this is BAD.
 *
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



#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "rsb_internals.h"
#include "rsb_lock.h"
#if RSB_USE_LIBRSBPP
#include "librsbpp.h"
#endif /* RSB_USE_LIBRSBPP */
#if RSB_USE_MKL
#include "rsb_mkl.h"
#endif /* RSB_USE_MKL */
RSB_INTERNALS_COMMON_HEAD_DECLS
#define RSB_WANT_OUTER_SPMM_BETA_SCALE 1 /* perform beta scaling already in rsb__do_spmv_general */
#define RSB_TRY_SPMM_BETA_SCALE_CONTIG 1 /* if right hand side multivecs contiguous, loop once */
/* FIXME: to move these macros to one header in order to avoid any identifier clash */
/* #define RSB_CBLAS_X_SCAL_SPMV rsb__cblas_Xscal */
#define RSB__FOREACH_NRHS(NRHSI, NRHS) for (NRHSI=0;NRHSI<NRHS;++NRHSI)
#define RSB_CBLAS_X_SCAL_SPMV(TYPECODE,N,ALPHAP,A,STRIDE) rsb__cblas_Xscal_parallel((TYPECODE),(N),(ALPHAP),(A),(STRIDE))
#if RSB_ENABLE_INNER_NRHS_SPMV
#define RSB_CBLAS_X_SCAL_SPMM(TYPECODE,N,ALPHAP,A,STRIDE) 					\
{	/* FIXME: this interacts with RSB_INNER_NRHS_SPMV_ARGS */					\
	rsb_int_t nrhsi = 0;										\
	const rsb_bool_t row_major = ( nrhs>1 && incy>=nrhs ) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;		\
	const rsb_bool_t is_compact = ( row_major ? (STRIDE)==(nrhs) : (N)==(outnri) ) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE; \
	if ( is_compact && RSB_TRY_SPMM_BETA_SCALE_CONTIG )					\
		RSB_CBLAS_X_SCAL_SPMV(TYPECODE, (nrhs)*(N), ALPHAP, A, 1); 				\
	else 												\
	RSB__FOREACH_NRHS(nrhsi,nrhs)									\
	{												\
		RSB_CBLAS_X_SCAL_SPMV(TYPECODE, N, ALPHAP, RSB_TYPED_OFF_PTR(TYPECODE,A,nrhsi*(outnri)), STRIDE); 	\
	}												\
}	/* RSB_CBLAS_X_SCAL_SPMM */
#else  /* RSB_ENABLE_INNER_NRHS_SPMV */
#define RSB_CBLAS_X_SCAL_SPMM(TYPECODE,N,ALPHAP,A,STRIDE) RSB_CBLAS_X_SCAL_SPMV(TYPECODE,N,ALPHAP,A,STRIDE) 
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */

rsb_err_t rsb__do_spmv_non_recursive(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA RSB_INNER_NRHS_SPMV_ARGS)
{
	/**
	  	\ingroup gr_internals
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_bool_t nostride = ( incx == 1 && incy == 1 )?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
	const rsb_bool_t should_scale_y = ( betap && !RSB_IS_ELEMENT_ONE( betap,mtxAp->typecode))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
	const rsb_bool_t use_alpha_one = (!alphap || RSB_IS_ELEMENT_ONE(alphap,mtxAp->typecode))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
	const rsb_bool_t use_y_zeroing_kernel = (should_scale_y && RSB_IS_ELEMENT_ZERO(betap,mtxAp->typecode) && nostride && use_alpha_one)?RSB_BOOL_TRUE:RSB_BOOL_FALSE;

	/*
		FIXME : should handle beta in a more specialized fashion.
			should also handle more specialized alphap cases.
	*/
#if RSB_ENABLE_INNER_NRHS_SPMV
	/*const size_t outtot=0,rhstot=0;
	const size_t outnri=0,rhsnri=0;
	const rsb_int_t nrhs=1;*/
 	/* FIXME: the above should be specified from argument */
	const size_t lenx=(mtxAp->el_size*rhsnri);
	const size_t leny=(mtxAp->el_size*outnri);
	rsb_int_t nrhsi=0;

#if RSB_USE_LIBRSBPP
	if(rsb_global_session_handle.use_rsbpp)
	if(! should_scale_y )
	if( (incx == 1 && incy == 1 && nrhs == 1) || (nrhs > 1) ) /* nrhs>1 && ldx=incx>=nrhs && ldy=incy>=nrhs is row major case */
	if(!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
	if( mtxAp->nnz > 0 )
	{
		const size_t scoff=mtxAp->coff;
		const size_t sroff=mtxAp->roff;
		void * oy=((char*)y)-(mtxAp->el_size*sroff)*incy;
		const void * ox=((const char*)x)-(mtxAp->el_size*scoff)*incx;
		const rsb_coo_idx_t ldX = (incx==1) ? rhsnri : incx;
		const rsb_coo_idx_t ldY = (incy==1) ? outnri : incy;
		const rsb_bool_t order = ( nrhs>1 && incx>=nrhs && incy>=nrhs ) ? RSB_FLAG_WANT_ROW_MAJOR_ORDER : RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;

		RSB_ASSERT(scoff>=0);
		RSB_ASSERT(sroff>=0);
		RSB_ASSERT(mtxAp->VA);
		RSB_ASSERT(mtxAp->bpntr);
		RSB_ASSERT(mtxAp->bindx);

		if( ! ( order == RSB_FLAG_WANT_ROW_MAJOR_ORDER && transA != RSB_TRANSPOSITION_N ) )
		if( rsb__is_coo_matrix(mtxAp) )
		{
			errval = rsbpp_coo_spmm(mtxAp->typecode, mtxAp->flags, mtxAp->nnz, mtxAp->nr ,mtxAp->nc, mtxAp->VA, mtxAp->bpntr, mtxAp->bindx,nrhs,ldX,ox,ldY,oy, alphap, incx, incy, transA, mtxAp->roff, mtxAp->coff, order );
			return errval;
		}

		if( ! ( order == RSB_FLAG_WANT_ROW_MAJOR_ORDER && transA != RSB_TRANSPOSITION_N ) )
		if( rsb__is_csr_matrix(mtxAp) )
		{
			errval = rsbpp_csr_spmm(mtxAp->typecode, mtxAp->flags, mtxAp->nnz, mtxAp->nr, mtxAp->nc, mtxAp->VA, mtxAp->bpntr, mtxAp->bindx,nrhs,ldX,ox,ldY,oy, alphap, incx, incy, transA, mtxAp->roff, mtxAp->coff, order );
			return errval;
		}
	}
#endif /* RSB_USE_LIBRSBPP */

#if RSB_USE_MKL
	if(rsb_global_session_handle.use_mkl)
	if( nrhs == 1 ) /* need to achieve equivalence in ldb/ldc in rsb__do_mkl_csr_spmm */
	if( sizeof(MKL_INT) == sizeof(rsb_coo_idx_t) && sizeof(MKL_INT) == sizeof(rsb_nnz_idx_t) )
	if(! should_scale_y )
	if( incx == 1 && incy == 1 )
	if(!RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
	if(!RSB_DO_FLAG_HAS(mtxAp->flags,(RSB_FLAG_USE_HALFWORD_INDICES)))
	if(!RSB_DO_FLAG_HAS_INTERSECTION(mtxAp->flags,(RSB_FLAG_ANY_SYMMETRY)))
	if( transA == RSB_TRANSPOSITION_N )
	if( mtxAp->nnz > 0 )
	{
		const rsb_coo_idx_t ldX = (incx==1) ? rhsnri : incx;
		const rsb_coo_idx_t ldY = (incy==1) ? outnri : incy;
		rsb_aligned_t alpha[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
		rsb_aligned_t beta[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];

		rsb__util_set_area_to_converted_integer(alpha,mtxAp->typecode,1);
		rsb__util_set_area_to_converted_integer(beta ,mtxAp->typecode,1);

		if( rsb__is_coo_matrix(mtxAp) )
		{
			if( nrhs == 1 )
				errval = rsb__mkl_coo_spmv(mtxAp->VA, mtxAp->nr, mtxAp->nc, mtxAp->nnz, mtxAp->bpntr, mtxAp->bindx, x, y, alphap?alphap:alpha, betap?betap:beta, transA, mtxAp->typecode, mtxAp->flags);
			else
				errval = rsb__mkl_coo_spmm(mtxAp->VA, mtxAp->nr, mtxAp->nc, nrhs, mtxAp->nnz, mtxAp->bpntr, mtxAp->bindx, x, ldX, y, ldY, alphap?alphap:alpha, betap?betap:beta, transA, mtxAp->typecode, mtxAp->flags);
			return errval;
		}

		if( rsb__is_csr_matrix(mtxAp) )
		{
			if( nrhs == 1 )
				errval = rsb__do_mkl_csr_spmv(mtxAp->VA, mtxAp->nr, mtxAp->nc, mtxAp->nnz, mtxAp->bpntr, mtxAp->bindx, x, y, alphap?alphap:alpha, betap?betap:beta, transA, mtxAp->typecode, mtxAp->flags);
			else
				errval = rsb__do_mkl_csr_spmm(mtxAp->VA, mtxAp->nr, mtxAp->nc, nrhs, mtxAp->nnz, mtxAp->bpntr, mtxAp->bindx, x, ldX, y, ldY, alphap?alphap:alpha, betap?betap:beta, transA, mtxAp->typecode, mtxAp->flags);
			return errval;
		}
	}
#endif /* RSB_USE_MKL */

	for(nrhsi=0;nrhsi<nrhs;++nrhsi)
	{
		void      *out=((      rsb_byte_t*)y)+(leny*nrhsi);
		const void*rhs=((const rsb_byte_t*)x)+(lenx*nrhsi);
#else /* RSB_ENABLE_INNER_NRHS_SPMV */
		void      *out=((      rsb_byte_t*)y);
		const void*rhs=((const rsb_byte_t*)x);
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */

	if(should_scale_y && !use_y_zeroing_kernel)
		RSB_CBLAS_X_SCAL_SPMV(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transA),betap,out,incy);
	/* no beta specified counts as beta=1, and so no scaling is needed */

	if(use_alpha_one)
	{
		/* no alpha specified counts as alpha=1 */
		if(nostride)
		{
			if(use_y_zeroing_kernel)
				/* y <- a * x  */
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_uauz(mtxAp,rhs,out,transA));
			else
				/* y <- y + a * x  */
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_uaua(mtxAp,rhs,out,transA));
		}
		else
			/* y <- a * x  , with stride */
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_sasa(mtxAp,rhs,out,incx,incy,transA));
	}
	else
	{
		if(nostride)
		{
			/* y <- - a * x  */
			if(RSB_IS_ELEMENT_MINUS_ONE(alphap,mtxAp->typecode))
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_unua(mtxAp,rhs,out,transA));
			/* y <- alpha * a * x  */
			else
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_uxua(mtxAp,rhs,out,alphap,transA));
		}
		else
			/* y <- alpha * a * x  , with stride */
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_sxsa(mtxAp,rhs,out,alphap,incx,incy,transA));
	}
#if RSB_ENABLE_INNER_NRHS_SPMV
	}
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
	RSB_DO_ERR_RETURN(errval)
}
#if RSB_WANT_SPMV_TRACE
#define RSB_SPMV_MAX_VS_ON_STACK 1024
#define RSB_SPMV_VS_ON_STACK(N) ((N)>(RSB_SPMV_MAX_VS_ON_STACK)?(0):(N))
#define RSB_SPMV_VS_DECL					\
	rsb_submatrix_idx_t *mivr=NULL, mivi=0;			\
	struct rsb_optrace_t *mvta=NULL; /* matrix visit times array */	\
	rsb_submatrix_idx_t smivr[RSB_SPMV_VS_ON_STACK(all_leaf_matrices_n)];	\
	struct rsb_optrace_t smvta[RSB_SPMV_VS_ON_STACK(all_leaf_matrices_n)];
#define RSB_SPMV_VS_ALLOC(MTXAP,ERRVAL,ERRL)			\
	if ( rsb_global_session_handle.want_spmv_trace != 0 )	\
		op_flags |= RSB_OP_FLAG_WANT_TRACE_PLOT;	\
	if(op_flags & RSB_OP_FLAG_WANT_TRACE_PLOT)		\
	{							\
		if (RSB_SPMV_VS_ON_STACK(all_leaf_matrices_n))	\
		{							\
			mvta=&smvta[0];	\
			mivr=&smivr[0];	\
		}			\
		else			\
		{			\
			mivr=rsb__calloc(sizeof(rsb_submatrix_idx_t)*(all_leaf_matrices_n));	\
			if(!mivr){ERRVAL=RSB_ERR_ENOMEM;goto ERRL;}	\
			mvta=rsb__calloc(sizeof(struct rsb_optrace_t)*(all_leaf_matrices_n));	\
			if(!mvta){ERRVAL=RSB_ERR_ENOMEM;goto ERRL;}	\
		}			\
	}
#define RSB_SPMV_VS_MARK_PRE(MI)				\
	if(op_flags & RSB_OP_FLAG_WANT_TRACE_PLOT){ mivr[mivi++]=(MI); mvta[(MI)].t0=rsb_time(); }
#define RSB_SPMV_VS_MARK_POST(MI)	\
	if(op_flags & RSB_OP_FLAG_WANT_TRACE_PLOT){ mvta[(MI)].t1=+rsb_time(); mvta[(MI)].th_id=th_id; }
#define RSB_SPMV_VS_DUMP(MTXAP)				\
	if(op_flags & RSB_OP_FLAG_WANT_TRACE_PLOT)	\
	{						\
		 /* for(mivi=0;mivi<((MTXAP)->all_leaf_matrices_n);++mivi) printf("%d ",(int)mivr[mivi]);printf("\n"); */ 	\
		const rsb_time_t t0 = mvta[mivr[0]].t0; \
		for(mivi=0;mivi<((MTXAP)->all_leaf_matrices_n);++mivi) \
			mvta[mivi].t0-=t0, mvta[mivi].t1-=t0; 	\
		rsb__dump_postscript_recursion_from_mtx_t(NULL,"spmv-times.eps",(MTXAP),1,1,RSB_DEFAULT_MATRIX_RENDERING_ROWS,RSB_DEFAULT_MATRIX_RENDERING_COLS,RSB_MARF_EPS_T,0,1,0,mivr,mvta,NULL);\
		rsb__dump_postscript_recursion_from_mtx_t(NULL,"spmv-trace.eps",(MTXAP),1,1,RSB_DEFAULT_MATRIX_RENDERING_ROWS,RSB_DEFAULT_MATRIX_RENDERING_COLS,RSB_FLAG_NOFLAGS|RSB_MARF_EPS_L,0,1,0,mivr,mvta,NULL);\
		rsb__dump_multiple_recursion_postscript_from_mtx_t("spmv-frame-",(MTXAP),1,1,RSB_DEFAULT_MATRIX_RENDERING_ROWS,RSB_DEFAULT_MATRIX_RENDERING_COLS,RSB_MARF_EPS_L|RSB_MARF_EPS_O,0,1,0,mivr,mvta);\
	}
#define RSB_SPMV_VS_DEALLOC		\
	if (all_leaf_matrices_n>RSB_SPMV_MAX_VS_ON_STACK)	\
	{							\
		RSB_CONDITIONAL_FREE(mivr);			\
		RSB_CONDITIONAL_FREE(mvta);			\
	}
#else /* 0 */
#define RSB_SPMV_VS_DECL
#define RSB_SPMV_VS_ALLOC(MTXAP,ERRVAL,ERRL)
#define RSB_SPMV_VS_MARK_PRE(MI)
#define RSB_SPMV_VS_MARK_POST(MI)
#define RSB_SPMV_VS_DUMP(MTXAP)
#define RSB_SPMV_VS_DEALLOC
#endif /* 0 */

rsb_err_t rsb__do_spmv_recursive_parallel(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPMV_ARGS)
{
	/**
	  	\ingroup gr_internals
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	const struct rsb_translated_matrix_t * all_leaf_matrices=NULL;
	struct rsb_spmv_lock_struct_t lock;
	const rsb_submatrix_idx_t all_leaf_matrices_n=mtxAp->all_leaf_matrices_n;
	const size_t el_size = mtxAp->el_size;
	RSB_SPMV_VS_DECL


	if(!rsb__is_recursive_matrix(mtxAp->flags))
		return rsb__do_spmv_non_recursive(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS);

	all_leaf_matrices  =mtxAp->all_leaf_matrices;

	if(!all_leaf_matrices || all_leaf_matrices_n<1)
	{errval = RSB_ERR_ENOMEM;goto err;}

	errval = rsb_do_spmv_lock_init(&lock,rsb_global_session_handle.rsb_want_threads,all_leaf_matrices_n,mtxAp,op_flags,transA,y,incy);
	if(RSB_SOME_ERROR(errval))
		goto err;

#if 0
	if(betap && !RSB_IS_ELEMENT_ONE(betap,mtxAp->typecode))
	{
	#pragma omp parallel shared(y,mtxAp,rsb_global_session_handle)  RSB_NTC 
	{
		rsb_nnz_idx_t tdim = rsb__do_get_rows_of(mtxAp,transA),dim,chunk;
		rsb_char_t * yy=y;
		rsb_thr_t th_id = omp_get_thread_num();
		if(th_id >= rsb_global_session_handle.rsb_want_threads)
			goto scaled;
		chunk=tdim/rsb_global_session_handle.rsb_want_threads;
		yy+=th_id*chunk*incy;
		if(th_id == rsb_global_session_handle.rsb_want_threads-1)
			dim=tdim-th_id*chunk;
		else		
			dim=chunk;

		if(RSB_IS_ELEMENT_ZERO(betap,mtxAp->typecode))
			RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,dim,NULL,yy,incy);
		else
			RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,dim,betap,yy,incy);
		scaled:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	}
	}
	#pragma omp barrier
#else /* 0 */
	if(betap && !RSB_IS_ELEMENT_ONE(betap,mtxAp->typecode))
		RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transA),betap,y,incy);
#endif /* 0 */
	RSB_SPMV_VS_ALLOC(mtxAp,errval,err)
	#pragma omp parallel reduction(|:errval) shared(lock,all_leaf_matrices,mtxAp)  RSB_NTC 
{
	const size_t want_threads = rsb_global_session_handle.rsb_want_threads;
	const size_t max_threads = RSB_MIN(all_leaf_matrices_n,want_threads);
	const rsb_thr_t th_id = omp_get_thread_num();
	rsb_submatrix_idx_t n = 0, dm = 0;

	if(RSB_UNLIKELY(th_id >= max_threads))
		goto skip;
again:
	for(n=0;RSB_LIKELY(n<all_leaf_matrices_n);++n)
#if RSB_WANT_SM_TO_THREAD_MOD_MAPPING && !RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV 
	if( (n % max_threads) == th_id )
#endif /* RSB_WANT_SM_TO_THREAD_MOD_MAPPING */
	{
		const struct rsb_mtx_t * const submatrix = all_leaf_matrices[n].mtxlp;
		rsb_bool_t gomv = RSB_BOOL_FALSE;
#if RSB_WANT_SPMV_WITH_REDUCE
		rsb_coo_idx_t rh,r0;
		char *ov = y;
		rsb_coo_idx_t oincy = incy;
#else
		char * const ov = y;
		const rsb_coo_idx_t oincy = incy;
#endif /* RSB_WANT_SPMV_WITH_REDUCE */
		#pragma omp critical (rsb_spmv_crs)
		{
#if RSB_WANT_BOUNDED_BOXES_SPMV
			const rsb_coo_idx_t roff = submatrix->broff;
			const rsb_coo_idx_t nr = RSB_MTX_EFF_R(submatrix);
			const rsb_coo_idx_t coff = submatrix->bcoff;
			const rsb_coo_idx_t nc = RSB_MTX_EFF_C(submatrix);
#else /* RSB_WANT_BOUNDED_BOXES_SPMV */
			const rsb_coo_idx_t roff = submatrix->roff;
			const rsb_coo_idx_t nr = submatrix->nr;
			const rsb_coo_idx_t coff = submatrix->coff;
			const rsb_coo_idx_t nc = submatrix->nc;
#endif /* RSB_WANT_BOUNDED_BOXES_SPMV */
			gomv=(rsb_do_spmv_lock_get(&lock,th_id,roff,nr,coff,nc,n,transA,&ov,&oincy)==RSB_BOOL_TRUE);
			if(gomv==RSB_BOOL_TRUE){RSB_SPMV_VS_MARK_PRE(n);}
		}
		if(gomv == RSB_BOOL_TRUE)
		{
			const size_t scoff = submatrix->coff-mtxAp->coff;
			const size_t sroff = submatrix->roff-mtxAp->roff;
			const char * const offx = ((const char*)x)+(el_size*scoff)*incx;
			char * const offy = ((char*) ov)+(el_size*sroff)*oincy;

			RSB_ASSERT(scoff>=0);
			RSB_ASSERT(sroff>=0);
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_non_recursive(submatrix,offx,offy,alphap,NULL,incx,oincy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS));
                       	#pragma omp critical (rsb_spmv_crs)
			{rsb_do_spmv_lock_release(&lock,th_id,ov);RSB_DO_SPMV_LOCK_DM_INC(lock);}
			RSB_SPMV_VS_MARK_POST(n);
		}
#if RSB_WANT_SPMV_WITH_REDUCE
		if(gomv == RSB_BOOL_ALMOST_TRUE)
		{
                       	#pragma omp critical (rsb_spmv_crs)
			{rsb__do_pick_candidate_interval_for_reduce(&lock,th_id,&ov,&r0,&rh);}

			if(ov && ov!=y)
			{
				rsb__vectors_left_sum_reduce_and_zero(y,ov,mtxAp->typecode,rh,oincy,r0);
                       		#pragma omp critical (rsb_spmv_crs)
                       		{ rsb__do_release_candidate_interval_for_reduce(&lock,th_id,ov,r0,rh);}
			}
		}
#endif /* RSB_WANT_SPMV_WITH_REDUCE*/
	}
		#pragma omp critical (rsb_spmv_crs)
		{ dm = RSB_DO_SPMV_LOCK_DM(lock); }
		if(dm<all_leaf_matrices_n
#if RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV
			&& ((all_leaf_matrices_n-dm)>th_id)
#endif	/* RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV */
		)goto again;
skip:
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	/* done */
}
	RSB_SPMV_VS_DUMP(mtxAp)
err:
	RSB_SPMV_VS_DEALLOC
#if   !defined(__xlC__)
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	/* FIXME: xlc does not allow this, but we have experienced problems, without */
	#pragma omp barrier
#endif /* __xlC__ */
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_spmv_lock_free(&lock));
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	errval = RSB_ERR_UNIMPLEMENTED_YET;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_spmv_recursive_serial(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA RSB_INNER_NRHS_SPMV_ARGS)
{
	/**
	  	\ingroup gr_internals
		This function does not offer result vector accumulation in case of diagonal implicit matrices.
	*/
	struct rsb_mtx_t * submatrix=NULL;
	rsb_submatrix_idx_t i,j;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_bool_t should_scale_y = ( betap && !RSB_IS_ELEMENT_ONE( betap,mtxAp->typecode) ) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;

	if( rsb__is_recursive_matrix(mtxAp->flags) )
	{
		if( should_scale_y )
		{
			RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transA),betap,y,incy);
		}

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			const size_t scoff=submatrix->coff-mtxAp->coff;
			const size_t sroff=submatrix->roff-mtxAp->roff;
			void       *offy=((char      *)y)+(mtxAp->el_size*sroff)*incy;
			const void *offx=((const char*)x)+(mtxAp->el_size*scoff)*incx;
			RSB_ASSERT(scoff>=0);
			RSB_ASSERT(sroff>=0);
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_recursive_serial(submatrix,offx,offy,alphap,NULL,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS));
		}
	}
	else
	{
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_non_recursive(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS));
	}
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_spmv_general(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, const void * x, rsb_coo_idx_t incx, const void * betap, void * y, rsb_coo_idx_t incy, enum rsb_op_flags_t op_flags RSB_OUTER_NRHS_SPMV_ARGS)
{
	/**
	  	\ingroup gr_internals
	*/
	rsb_err_t errval = RSB_ERR_BADARGS;
	const rsb_type_t typecode = mtxAp->typecode;

#if RSB_ALLOW_ZERO_DIM
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
	{
		errval = RSB_ERR_NO_ERROR;
		goto err; /* FIXME: skipping further checks */
	}
#endif
	if(x==y)
		goto err;

	if(incx<1 || incy<1)
		goto err;

/*	we tolerate NULL alhap and betap */
#if 0
	if(!alphap || !betap)
		goto err;
#endif /* 0 */

	if(!mtxAp || !x || !y || transA == RSB_INVALID_FLAGS)
		goto err;

#if RSB_WANT_OUTER_SPMM_BETA_SCALE
	if(betap && !RSB_IS_ELEMENT_ONE(betap,typecode))
	{
		RSB_CBLAS_X_SCAL_SPMM(typecode,rsb__do_get_rows_of(mtxAp,transA),betap,y,incy);
		betap = NULL;
	}
#endif /* RSB_WANT_OUTER_SPMM_BETA_SCALE */

#if RSB_WANT_OMP_RECURSIVE_KERNELS
	if(RSB_LIKELY(op_flags != RSB_OP_FLAG_WANT_SERIAL))
	{
		RSB_NUM_THREADS_DECL
		RSB_NUM_THREADS_PUSH
		errval = rsb__do_spmv_recursive_parallel(mtxAp,x,y,alphap,betap,incx,incy,transA,op_flags RSB_OUTER_NRHS_SPMV_ARGS_IDS	);
		RSB_NUM_THREADS_POP
	}
	else
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		errval = rsb__do_spmv_recursive_serial(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS);

	/* Note: the RSB_OP_FLAG_FAKE_LOCK case is handled by rsb__do_spmv_recursive_parallel */
	goto done;
done:
	if(!RSB_UNLIKELY(op_flags&RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT))// NEW: fix for odd spsv/diagonal implicit/no-parallel cases
	if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
	{
		const rsb_coo_idx_t ndy = RSB_MIN(mtxAp->nr,mtxAp->nc);
		const int row_major = ( nrhs > 1 && (incx >= nrhs || incy >= nrhs ) );
		const rsb_nnz_idx_t ldY = rsb__do_get_rows_of(mtxAp,transA);
		const rsb_nnz_idx_t ldX = rsb__do_get_columns_of(mtxAp,transA);
		rsb_int_t nrhsi, di;

		if( row_major )
		for (di=0;di<ndy;++di)
			rsb__cblas_Xaxpy(typecode,nrhs,alphap
				,RSB_TYPED_OFF_PTR(typecode,x,incx*di)
				,1
				,RSB_TYPED_OFF_PTR(typecode,y,incy*di)
				,1);
		else
		for (nrhsi=0;nrhsi<nrhs;++nrhsi)
			rsb__BLAS_Xaxpy_parallel(ndy,alphap
				,RSB_TYPED_OFF_PTR(typecode,y,incy*nrhsi*ldY)
				,incy
				,RSB_TYPED_OFF_PTR(typecode,x,incx*nrhsi*ldX)
				,incx,typecode);
	}
err:
	RSB_DO_ERR_RETURN(errval)
}


#ifdef __cplusplus
}
#endif  /* __cplusplus */

/* @endcond */
