/* @cond INNERDOC */
/*! 
 @file
 @brief Matrix type dispatching code, for each matrix operation.
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

#ifndef RSB_DISPATCH_H_INCLUDED
#define RSB_DISPATCH_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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
#include "rsb_common.h"
#include "rsb_krnl_bcss_spmv_u.h"	/* uhm */
#include "rsb_krnl_bcss_spsv_u.h"	/* uhm */
#include "rsb_krnl_bcss_misc_u.h"	/* uhm */



#define	RSB_BCSR_GET_NEXT_BLOCK_POINTER(BP,mtxAp,ROWVAR,COLVAR,BLOCKROWSPAR,BLOCKCOLSPAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	/*										\
	 * *input*									\
	 * mtxAp		should be a valid rsb_mtx_t structure pointer		\
	 * BLOCKROWSPAR	should be set to the rows   count of this block			\
	 * BLOCKCOLSPAR	should be set to the column count of this block			\
	 * *output*									\
	 * ROWVAR	will be set to the base row    of this block			\
	 * COLVAR	will be set to the base column of this block			\
	 * BP		will be set to the current block pointer			\
	 * */										\
	while( (mtxAp)->bpntr[_i] == (mtxAp)->bpntr[_i+1] ) 				/* skipping empty rows */	\
	{++_i;_k=(mtxAp)->bpntr[_i];} 		/* _k is the first block index for the current row of blocks */	\
	_j=(mtxAp)->bindx[_k]; 						/* the current block column index  */	\
	_lastk=_k;	\
	(BLOCKROWVAR)=_i;	\
	(BLOCKCOLUMNVAR)=_j;	\
	(ROWVAR)=(BLOCKROWSPAR)*_i;					/* _i is the current block row index */	\
	(COLVAR)=(BLOCKCOLSPAR)*_j; 					/* the current block column index  */	\
	BP+=(mtxAp)->options->el_size*(BLOCKROWSPAR)*(BLOCKCOLSPAR);			\
	_k++; 		/* for the future macro calls */						\
	if( _k >= (mtxAp)->bpntr[_i+1] )++_i;								\
	;

#define RSB_BCSR_GET_FIRST_BLOCK_POINTER(BP,mtxAp,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	int _i=0,_j=0,_k=0,_lastk=0;									\
	(BLOCKROWSVAR)=(mtxAp)->rpntr[1]-(mtxAp)->rpntr[0];		/* _i is the current block row index */	\
	(BLOCKCOLSVAR)=(mtxAp)->cpntr[1]-(mtxAp)->cpntr[0]; 		/* the current block column index  */	\
	(BP)=(mtxAp)->VA;											\
	RSB_BCSR_GET_NEXT_BLOCK_POINTER(BP,mtxAp,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)

#if RSB_WANT_DBC
#define RSB_BCSR_GOT_LAST_BLOCK_POINTER(mtxAp)	( _lastk >= (mtxAp)->block_count )
#endif


#define RSB_BENCHMARK_MIN_SECONDS	/*0.5*/1.0
#define RSB_BENCHMARK_MIN_RUNS		/*5*/10 


rsb_err_t rsb__do_spmv_uaua(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_spmv_uauz(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_spmv_uxua(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);

rsb_err_t rsb__do_spmv_unua(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_spmv_sasa(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_spsv_uxua(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_spmv_sxsa(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_spsv_sxsx(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_infty_norm(const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_rowssums(const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_scale(struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);


#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spmv_uaua_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spmv_uaua(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spmv_uauz_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spmv_uauz(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spmv_uxua_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spmv_uxua(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spmv_unua_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spmv_unua(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spmv_sasa_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spmv_sasa(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spsv_uxua_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spsv_uxua(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spmv_sxsa_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spmv_sxsa(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__spsv_sxsx_testing(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_spsv_sxsx(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__infty_norm_testing(const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_infty_norm(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__rowssums_testing(const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_rowssums(double * elapsed_time, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t rsb__scale_testing(struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);
#endif /* RSB_WANT_KERNELS_DEBUG */

rsb_err_t rsb_do_time_scale(double * elapsed_time, struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);


rsb_err_t rsb__do_fullrangebenchmark_double_spmv_uaua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spmv_uaua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spmv_uaua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spmv_uaua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spmv_uaua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spmv_uaua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spmv_uaua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spmv_uaua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_spmv_uauz(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spmv_uauz(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spmv_uauz(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spmv_uauz(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spmv_uauz(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spmv_uauz(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spmv_uauz(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spmv_uauz(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_spmv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spmv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spmv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spmv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spmv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spmv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spmv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spmv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_spmv_unua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spmv_unua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spmv_unua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spmv_unua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spmv_unua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spmv_unua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spmv_unua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spmv_unua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_spmv_sasa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spmv_sasa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spmv_sasa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spmv_sasa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spmv_sasa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spmv_sasa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spmv_sasa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spmv_sasa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_spsv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spsv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spsv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spsv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spsv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spsv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spsv_uxua(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spsv_uxua(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_spmv_sxsa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spmv_sxsa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spmv_sxsa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spmv_sxsa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spmv_sxsa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spmv_sxsa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spmv_sxsa(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spmv_sxsa(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_spsv_sxsx(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_spsv_sxsx(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_spsv_sxsx(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_spsv_sxsx(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_spsv_sxsx(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_spsv_sxsx(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_spsv_sxsx(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_spsv_sxsx(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);

rsb_err_t rsb__do_fullrangebenchmark_double_infty_norm(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_infty_norm(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_float_infty_norm(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_infty_norm(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_infty_norm(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_infty_norm(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_infty_norm(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_infty_norm(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_double_rowssums(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_rowssums(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_float_rowssums(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_rowssums(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_rowssums(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_rowssums(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_rowssums(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_rowssums(double * total_elapsed_time, double * m_flops, const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);

rsb_err_t rsb__do_fullrangebenchmark_double_scale(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_scale(double * total_elapsed_time, double * m_flops, struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);

rsb_err_t rsb__do_fullrangebenchmark_float_scale(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_scale(double * total_elapsed_time, double * m_flops, struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);

rsb_err_t rsb__do_fullrangebenchmark_float_complex_scale(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_float_complex_scale(double * total_elapsed_time, double * m_flops, struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);

rsb_err_t rsb__do_fullrangebenchmark_double_complex_scale(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags);

rsb_err_t rsb__do_benchmark_double_complex_scale(double * total_elapsed_time, double * m_flops, struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);





#if 0
rsb_err_t rsb__do_spmv_uaua_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spmv_uaua(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_spmv_uauz_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spmv_uauz(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_spmv_uxua_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spmv_uxua(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_spmv_unua_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spmv_unua(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_spmv_sasa_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spmv_sasa(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_spsv_uxua_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spsv_uxua(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_spmv_sxsa_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spmv_sxsa(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_spsv_sxsx_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const void * restrict rhs, void * restrict out,const void * alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy,const rsb_trans_t transA);
#endif /* 0 */
double rsb__estimate_mflops_per_op_spsv_sxsx(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_infty_norm_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);
#endif /* 0 */
double rsb__estimate_mflops_per_op_infty_norm(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_rowssums_with_macros_vbr(const struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,void * row_sums);
#endif /* 0 */
double rsb__estimate_mflops_per_op_rowssums(const struct rsb_mtx_t * mtxAp);

#if 0
rsb_err_t rsb__do_scale_with_macros_vbr(struct rsb_mtx_t * mtxAp,const rsb_trans_t transA,const void * scale_factors);
#endif /* 0 */
double rsb__estimate_mflops_per_op_scale(const struct rsb_mtx_t * mtxAp);

rsb_err_t rsb__do_completebenchmark(const int argc, char *const argv[]);
rsb_err_t rsb__dump_performance_array(const char * an, const double*array);

#ifdef __cplusplus
}
#endif  /* __cplusplus */


#endif	/* RSB_DISPATCH_H_INCLUDED */




/*!
 @file
 @brief ...
 */
/* @endcond */
