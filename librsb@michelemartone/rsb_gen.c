/*                                                                                                                            

Copyright (C) 2008-2019 Michele Martone

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
 * Matrix generating functions.
 * */

#include "rsb_internals.h"

rsb_err_t rsb__generate_blocked_banded_coo(rsb_nnz_idx_t dim, rsb_nnz_idx_t spacing, rsb_nnz_idx_t lbw, rsb_nnz_idx_t ubw, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_nnz_idx_t *nnzp, rsb_type_t typecode)
{
	/* 
	 * duuu
	 * lduuu
	 * llduuu
	 *  llduu
	 *   lldu
	 *    lld
	 *
	 * assuming lbw<dim an ubw<dim,
	 *
	 * there are dim 'd' type nonzeros
	 * there are lbw 'l' type nonzeros
	 * there are ubw*dim - (ubw*(ubw-1))/2 'u' type nonzeros
	 * there are lbw*dim - (lbw*(lbw-1))/2 'u' type nonzeros
	 *
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t nnz = RSB_NNZ_OF_BANDED(dim,lbw,ubw);
	rsb_coo_idx_t ri=0,ci=0,nzi=0;
	//rsb_time_t dt;

	/* overflow is possible here */
	if(RSB_INVALID_NNZ_COUNT(spacing)){errval = RSB_ERR_BADARGS;goto err;}
	if(RSB_INVALID_NNZ_COUNT(dim*spacing)){errval = RSB_ERR_BADARGS;goto err;}
	if(RSB_INVALID_NNZ_COUNT(dim)){errval = RSB_ERR_BADARGS;goto err;}
	if(lbw>0)if(RSB_INVALID_NNZ_COUNT(lbw)){errval = RSB_ERR_BADARGS;goto err;}
	if(ubw>0)if(RSB_INVALID_NNZ_COUNT(ubw)){errval = RSB_ERR_BADARGS;goto err;}
	if((ubw>=dim)||(lbw>=dim)){errval = RSB_ERR_BADARGS;goto err;}
	if(!VA || !JA || !IA || !nnzp) {errval = RSB_ERR_BADARGS;goto err;}
	if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc( VA, IA, JA,nnz,typecode,RSB_BOOL_FALSE))){goto err;}
	//dt = - rsb_time();
	for(ri=0;ri<lbw;++ri)
	for(ci=0;ci<RSB_MIN(ri+1+ubw,dim);++ci)
	{
		(*IA)[nzi]=ri;
		(*JA)[nzi]=ci;
		++nzi;
	}
	
	for(ri=lbw;ri<(dim-ubw);++ri)
	for(ci = ri-lbw;ci<1+ri+ubw;++ci)
	{
		(*IA)[nzi]=ri;
		(*JA)[nzi]=ci;
		++nzi;
	}
		
	for(ri = RSB_MAX(lbw,dim-ubw);ri<dim;++ri)
	for(ci = ri-lbw;ci<dim;++ci)
	{
		(*IA)[nzi]=ri;
		(*JA)[nzi]=ci;
		++nzi;
	}
	
	//dt += rsb_time();
	//printf("TIME: %lg, %d\n",dt,omp_get_num_threads());
	*nnzp=nnz;
	if((errval = rsb__fill_with_ones(*VA,typecode,nnz,1))!=RSB_ERR_NO_ERROR)goto err;
	if(spacing>1)
		rsb__util_coo_arrays_mul(*IA,*JA,spacing,spacing,nnz);

	goto ok;
err:
	RSB_CONDITIONAL_FREE(*VA);
	RSB_CONDITIONAL_FREE(*IA);
	RSB_CONDITIONAL_FREE(*JA);
ok:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__generate_dense_full(rsb_nnz_idx_t dim_r, rsb_nnz_idx_t dim_c, rsb_nnz_idx_t spacing, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void ** VA, rsb_nnz_idx_t *nnzp, rsb_type_t typecode)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t nnz = dim_r*dim_c,/*n=0,*/lda=dim_c;
	rsb_coo_idx_t ri=0,ci=0;
	//rsb_time_t dt;

	/* FIXME : overflow possible  */
	if(!VA || !JA || !IA || !nnzp) {errval = RSB_ERR_BADARGS;goto err;}
	if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc( VA, IA, JA,nnz,typecode,RSB_BOOL_FALSE))){goto err;}
	//dt = - rsb_time();
	for(ri=0;ri<dim_r;++ri){
	for(ci=0;ci<dim_c;++ci)
	{
		(*IA)[lda*ri+ci]=ri;
		(*JA)[lda*ri+ci]=ci;
	}}
	/* n=dim*dim; */
	//dt += rsb_time();
	//printf("TIME: %lg, %d\n",dt,omp_get_num_threads());
	*nnzp=nnz;
	if((errval = rsb__fill_with_ones(*VA,typecode,nnz,1))!=RSB_ERR_NO_ERROR)goto err;
	if(spacing>1)
		rsb__util_coo_arrays_mul(*IA,*JA,spacing,spacing,nnz);

	goto ok;
err:
	RSB_CONDITIONAL_FREE(*VA);
	RSB_CONDITIONAL_FREE(*IA);
	RSB_CONDITIONAL_FREE(*JA);
ok:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__generate_dense_lower_triangular_coo(rsb_nnz_idx_t dim, rsb_nnz_idx_t spacing, rsb_coo_idx_t ** IAp, rsb_coo_idx_t ** JAp, void ** VAp, rsb_nnz_idx_t *nnzp, rsb_type_t typecode)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t nnz = (dim%2)?(((dim+1)/2)*dim):((dim/2)*(dim+1)),n=0;
	rsb_coo_idx_t ri = 0, ci = 0;

	/* overflow is possible here! */
	if(!VAp || !JAp || !IAp || !nnzp)
       	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	if(nnz == 0)
		goto skalloc; /* tolerated corner case */

	if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc( VAp, IAp, JAp,nnz,typecode,RSB_BOOL_FALSE)))
	{
		goto err;
	}
skalloc:	
	for(ri=0;ri < dim;++ri)
	for(ci=0;ci <= ri;++ci)
	{
		(*IAp)[n]=ri;
		(*JAp)[n]=ci;
		++n;
	}
	*nnzp=nnz;
	if((errval = rsb__fill_with_ones(*VAp,typecode,nnz,1))!=RSB_ERR_NO_ERROR)
		goto err;
	if(spacing>1)
		rsb__util_coo_arrays_mul(*IAp,*JAp,spacing,spacing,nnz);

	goto ok;
err:
	RSB_CONDITIONAL_FREE(*VAp);
	RSB_CONDITIONAL_FREE(*IAp);
	RSB_CONDITIONAL_FREE(*JAp);
ok:
	RSB_DO_ERR_RETURN(errval)
}

#if RSB_WANT_DO_LOCK_TEST
struct rsb_mtx_t * rsb__generate_dense_lower_triangular(const rsb_coo_idx_t dim, double * timep, rsb_type_t typecode)
{
	void * VA=NULL;
	rsb_coo_idx_t * IA=NULL;
	rsb_coo_idx_t * JA=NULL;
	rsb_nnz_idx_t nnz = RSB_MARKER_NNZ_VALUE;
	struct rsb_mtx_t * mtxAp=NULL;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_MATRIX_FLAGS;
	rsb_time_t time;
	const rsb_blk_idx_t br=1,bc=1;
	rsb_blk_idx_t m=dim,k=dim;

	errval = rsb__generate_dense_lower_triangular_coo(dim,1,&IA,&JA,&VA,&nnz,typecode);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	time = - rsb_time();
	mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,m,k,br,bc,flags,&errval);
	time += rsb_time();
	if(!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_ES);
		rsb__do_perror(NULL,errval);
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	return mtxAp;
err:
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	return NULL;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

struct rsb_mtx_t * rsb__generate_banded(const rsb_blk_idx_t br, const rsb_blk_idx_t bc, const rsb_coo_idx_t rows, const rsb_coo_idx_t cols, rsb_coo_idx_t bw, double * timep, rsb_type_t typecode)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Generate and return a sparse banded matrix.
         *
         * Note that it is not performance optimized:
         * this will result in longer benchmarking times.
         *
	 * \return a matrix pointer, NULL in case of failure.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void * VA=NULL;
	rsb_coo_idx_t * IA=NULL;
	rsb_coo_idx_t * JA=NULL;
	struct rsb_mtx_t * mtxAp=NULL;
	rsb_coo_idx_t lbw=0;
	rsb_coo_idx_t rbw=0;
	rsb_coo_idx_t ri=0,ci=0;

	rsb_blk_idx_t rob;
	rsb_blk_idx_t cob;
	rsb_nnz_idx_t blockcount=0;
	size_t blockcount_=0;
	rsb_nnz_idx_t nnz=0;
	size_t  r=0,c=0;
	rsb_flags_t flags = RSB_FLAG_DEFAULT_STORAGE_FLAGS;
	rsb_time_t time;

	if( (bw*bc) > cols )
	{
		RSB_ERROR("too much bandwidth..\n");
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	if(rows != cols)
	{
		RSB_ERROR("matrix is not square..\n");
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	if(rows != cols)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	lbw=bw>0?(bw-1)/2:0;
	rbw=bw-lbw;
	ri=0;ci=0;/* NOTE : danger : it was signed and -1 before */

	/*!
	 * Will create a blocked banded matrix with bandwidth expressed
	 * in number of blocks, as in the experiments from [ijhcp_07] 
	 * (Buttari, Eijkhout, Langou, Filippone 2007) article from 
	 * International Journal of High Performance Computing Applications 2007; 21; 467
	 * 
	 * \code
	 * +*+-------+
	 * +*+*+     |  Example for bandwidth of 1 block
	 * | +*+     |
	 * |    ...  |
	 * |       +*+
	 * +-------+*+
	 *
	 * +*+*+-----+
	 * +*+*+*+   |  Example for bandwidth of 2 blocks
	 * | +*+*+   |
	 * |    ...  |
	 * |       +*+
	 * +-------+*+
	 *
	 * +*+*+-----+
	 * +*+*+*+   |  Example for bandwidth of 3 blocks
	 * +*+*+*+   |
	 * |    ...  |
	 * |     +*+*+
	 * +-----+*+*+
	 * \endcode
	 **/
		
	rob = rows/br;	/* rows of blocks */
	cob=cols/bc;	/* cols of blocks */

	for(ri=0;ri<rob;++ri)
	for(ci = ri-lbw;ci <= ri+rbw;++ci)
		if(ci>=0 && ci<cob )
			++blockcount_;	/* yes, we waste resources, but we are in hurry. FIXME */
	blockcount=blockcount_;
		
	if(blockcount<=0 || (((size_t)blockcount_)!=blockcount))
	{errval = RSB_ERR_INTERNAL_ERROR;goto err;}/* overflow */

	nnz = blockcount * ( br * bc );

	if(nnz<=0 || nnz<blockcount){errval = RSB_ERR_INTERNAL_ERROR;goto err;}/* overflow */

	if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc(&VA,&IA,&JA,nnz,typecode,RSB_BOOL_FALSE))){goto err;}

	for(ri=0;ri<rob;++ri)
	for(ci = ri-lbw;ci <= ri+rbw;++ci)
	if(ci>=0 && ci<cob )
	for( r=0; r< br;++ r)
	for( c=0; c< bc;++ c)
	{
		*IA = ri*br+r;
		*JA=ci*bc+c;
		++IA;++JA;
	}
	if((errval = rsb__fill_with_ones(VA,typecode,nnz,1))!=RSB_ERR_NO_ERROR)
		goto err;

	IA-=nnz;
	JA-=nnz;
		
//	p_r = rsb__util_get_partitioning_array(br,rows,&M_b,flags);
//	p_c = rsb__util_get_partitioning_array(bc,cols,&K_b,flags);

//	if(! p_r || !p_c) {errval = RSB_ERR_ENOMEM;goto err;}

	time = - rsb_time();
	mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,rows,cols,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,flags,&errval);
	time += rsb_time();
	if(timep)*timep=time;

	if(!mtxAp || (RSB_SOME_ERROR(errval))) {errval = RSB_ERR_ENOMEM;goto err;}

	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	return mtxAp;
err:
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
//	RSB_CONDITIONAL_FREE(p_r);
//	RSB_CONDITIONAL_FREE(p_c);
	RSB_MTX_FREE(mtxAp);
	return NULL;
}

void rsb__do_fill_with_diag(void *VA, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, rsb_coo_idx_t ioff, rsb_coo_idx_t joff, rsb_nnz_idx_t nzoff, rsb_type_t typecode, rsb_nnz_idx_t nnz)
{
	rsb_nnz_idx_t nzi;
	void *dVA=((char*)VA)+((size_t)RSB_SIZEOF(typecode))*nzoff;

	if(VA)
		rsb__fill_with_ones(dVA,typecode,nnz,1);

	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
	{
		IA[nzoff+nzi]=nzi+ioff;
		JA[nzoff+nzi]=nzi+joff;
	}
}

#if 0
struct rsb_mtx_t * rsb_generate_diagonal(const rsb_coo_idx_t rows, double * timep, rsb_type_t typecode)
{
	/**
	 * \ingroup gr_internals
	 *
	 * FIXME : untested, undocumented
 	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void * VA = NULL;
	rsb_nnz_idx_t nnz   = rows, cols = rows;
	rsb_coo_idx_t * IA   = NULL;
	rsb_coo_idx_t * JA   = NULL;
	rsb_coo_idx_t *p_r=NULL,*p_c=NULL;
	struct rsb_mtx_t * mtxAp=NULL;
	//rsb_flags_t flags = RSB_FLAG_OWN_PARTITIONING_ARRAYS | RSB_FLAG_WANT_FIXED_BLOCKING_VBR;
	rsb_time_t time;
	rsb_blk_idx_t M_b=0,K_b=0;
	rsb_flags_t flags = RSB_FLAG_OWN_PARTITIONING_ARRAYS;
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);
//	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER);	/* experimental */

	if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc(&VA,&IA,&JA,nnz,typecode,RSB_BOOL_FALSE))){goto err;}
	p_r = rsb__util_get_partitioning_array( 1, rows , &M_b, flags);
	p_c = rsb__util_get_partitioning_array( 1, cols , &K_b, flags);
	if(! p_r || !p_c) {errval = RSB_ERR_ENOMEM;goto err;}

	rsb__do_fill_with_diag(NULL,IA,JA,0,0,0,typecode,nnz);
	if((errval = rsb__fill_with_ones(VA,typecode,nnz,1))!=RSB_ERR_NO_ERROR)goto err;

	time = - rsb_time();
	mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,rows,cols,M_b,K_b,flags,&errval);
	time += rsb_time();

	if(timep)*timep=time;

	if(!mtxAp || (RSB_SOME_ERROR(errval))) {errval = RSB_ERR_ENOMEM;goto err;}

	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE( IA );
	RSB_CONDITIONAL_FREE( JA );
	return mtxAp;
err:
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE( IA );
	RSB_CONDITIONAL_FREE( JA );
	RSB_CONDITIONAL_FREE(p_r);
	RSB_CONDITIONAL_FREE(p_c);
	RSB_MTX_FREE(mtxAp);
	return NULL;
}
#endif

struct rsb_mtx_t * rsb__generate_blocked_banded(const rsb_blk_idx_t br, const rsb_blk_idx_t bc, const rsb_coo_idx_t rows, const rsb_coo_idx_t cols, const rsb_coo_idx_t bw, double * timep, rsb_type_t typecode,rsb_bool_t want_lowtri)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * A function which generates and returns a sparse banded matrix.
         *
         * Note that it is not performance optimized, and this will
         * result in longer benchmarking times.
         *
	 * \return a matrix pointer, NULL in case of failure.
	 */

	rsb_err_t errval = RSB_ERR_NO_ERROR;

	void * VA = NULL;
	rsb_coo_idx_t * IA   = NULL;
	rsb_coo_idx_t * JA   = NULL;
	rsb_coo_idx_t *p_r=NULL,*p_c=NULL;
	struct rsb_mtx_t * mtxAp=NULL;
	rsb_coo_idx_t lbw=(bw-1)/2;
	rsb_coo_idx_t rbw= bw-lbw;
	rsb_coo_idx_t ri=0,ci=0;/* NOTE : danger : it was signed and -1 before */
	size_t rob;
	size_t cob;
	rsb_nnz_idx_t nnz = 0;
	size_t blockcount = 0;
	rsb_blk_idx_t M_b=0,K_b=0;
	size_t  r =0, c =0;
	rsb_flags_t flags = RSB_FLAG_OWN_PARTITIONING_ARRAYS | RSB_FLAG_WANT_BCSS_STORAGE;
	rsb_time_t time;

	if( bw*bc >= cols ){RSB_ERROR("too much bandwidth..\n");errval = RSB_ERR_BADARGS;goto err;}
	if(rows != cols){RSB_ERROR("matrix is not square..\n");errval = RSB_ERR_BADARGS;goto err;}
	if(rows != cols){errval = RSB_ERR_BADARGS;goto err;}
	/*!
	 * Will create a blocked banded matrix with bandwidth expressed
	 * in number of blocks, as in the experiments from [ijhcp_07] 
	 * (Buttari, Eijkhout, Langou, Filippone 2007) article from 
	 * International Journal of High Performance Computing Applications 2007; 21; 467
	 * 
	 * \code
	 * +*+-------+
	 * +*+*+     |  Example for bandwidth of 1 block
	 * | +*+     |
	 * |    ...  |
	 * |       +*+
	 * +-------+*+
	 *
	 * +*+*+-----+
	 * +*+*+*+   |  Example for bandwidth of 2 blocks
	 * | +*+*+   |
	 * |    ...  |
	 * |       +*+
	 * +-------+*+
	 *
	 * +*+*+-----+
	 * +*+*+*+   |  Example for bandwidth of 3 blocks
	 * +*+*+*+   |
	 * |    ...  |
	 * |     +*+*+
	 * +-----+*+*+
	 * \endcode
	 **/

	rob = rows/br;	/* rows of blocks */
	cob=cols/bc;	/* cols of blocks */

	for(ri=0;ri<rob;++ri)
	for(ci = ri-lbw;ci <= ri+rbw;++ci)
	if(ci>=0 && ci<cob )
		++blockcount;	/* yes, we waste resources, but we are in hurry. FIXME */

	if(blockcount<=0){errval = RSB_ERR_INTERNAL_ERROR;goto err;}/* overflow */

	nnz=blockcount*(br*bc);

	if(nnz<=0){errval = RSB_ERR_INTERNAL_ERROR;goto err;}/* overflow */

	if(RSB_SOME_ERROR(errval = rsb__util_coo_alloc(&VA,&IA,&JA,nnz,typecode,RSB_BOOL_FALSE))){goto err;}
	if(want_lowtri)
	{
		nnz=0;
		/* FIXME: a dirty hack : will result in half the nonzeros :/ */
		for(ri=0;ri<rob;++ri)
		for(ci = ri-lbw;ci <= ri+rbw;++ci)
		if(ci>=0 && ci<cob )
		for( r=0; r< br;++ r)
		for( c=0; c< bc;++ c)
		{
			*IA   = ri*br+r;
			*JA   = ci*bc+c;
			if(*IA>=*JA)
				++IA,++JA,++nnz;
		}
	}
	else
	{
		for(ri=0;ri<rob;++ri)
		for(ci = ri-lbw;ci <= ri+rbw;++ci)
		if(ci>=0 && ci<cob )
		for( r=0; r< br;++ r)
		for( c=0; c< bc;++ c)
		{
			*IA   = ri*br+r;
			*JA   = ci*bc+c;
			++IA;++JA;
		}
	}
	
	IA -= nnz;
	JA -= nnz;
	if(rsb__fill_with_ones(VA,typecode,nnz,1)) { errval = RSB_ERR_INTERNAL_ERROR; goto err; }
	
	p_r = rsb__util_get_partitioning_array(br,rows,&M_b,flags);
	p_c = rsb__util_get_partitioning_array(bc,cols,&K_b,flags);

	if(! p_r || !p_c) {errval = RSB_ERR_ENOMEM;goto err;}

	time = - rsb_time();
	mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,rows,cols,M_b,K_b,flags,&errval);
	time += rsb_time();
	if(timep)*timep=time;

	if(!mtxAp) {goto err;}
	
	goto ret;
err:
	RSB_MTX_FREE(mtxAp);
ret:
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(p_r);
	RSB_CONDITIONAL_FREE(p_c);
	return mtxAp;
}

rsb_err_t rsb__generate_blocked_banded_mtx(rsb_nnz_idx_t dim, rsb_nnz_idx_t spacing, rsb_nnz_idx_t lbw, rsb_nnz_idx_t ubw, struct rsb_mtx_t ** mtxApp, rsb_type_t typecode)
{
	rsb_err_t errval = RSB_ERR_BADARGS;
	void * VA = NULL;
	rsb_nnz_idx_t nnz = dim;
	rsb_coo_idx_t rows = dim, cols = dim;
	rsb_coo_idx_t * IA = NULL;
	rsb_coo_idx_t * JA = NULL;

	if(!mtxApp)
		goto ret; 

	errval = rsb__generate_blocked_banded_coo(dim, spacing, lbw, ubw, &IA, &JA, &VA, &nnz, typecode);
	if(RSB_SOME_ERROR(errval))
	{
		goto ret;
	}
	*mtxApp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,rows,cols,RSB_DEFAULT_ROW_BLOCKING,RSB_DEFAULT_COL_BLOCKING,RSB_FLAG_DEFAULT_MATRIX_FLAGS,&errval);
ret:
	return errval;
}

/* @endcond */
