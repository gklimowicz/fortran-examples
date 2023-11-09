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
/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains locks for sparse recursive multicore operations.
 * */

#ifndef RSB_LOCK_H_INCLUDED
#define RSB_LOCK_H_INCLUDED

#include "rsb_internals.h"

#define RSB__TRSV_OUT  0
#define RSB__TRSV_OUT_ 0
#define RSB__TRSV_OUT__ 0
#define RSB_WANT_DO_LOCK_TEST 0 /* compile and use lock test function (broken, especially with assertions) */


#define RSB_CONST_MIN_SUPPORTED_CORES 	1
#define RSB_CONST_MAX_SUPPORTED_CORES 	RSB_CONST_MAX_SUPPORTED_THREADS /* The maximum number of cores (TODO: support any number of cores) */
#define RSB_CONST_MAX_SUPPORTED_TEMPORARY_VECTORS RSB_CONST_MAX_SUPPORTED_CORES

#define RSB__MAX_BITMAP_SUBMS_ON_STACK (4096) /* large matrices (say, >>1K subms) may justify a per-spmv/spsv malloc */

typedef int rsb_thr_t;

/*!
 * \ingroup gr_internals
 * \brief An internal, helper structure.
 */
struct rsb_rows_lock_struct_t
{
	/* FIXME : EXPERIMENTAL,NEW  */
	/* FIXME : THE LOCK SHULD BE SIZED PROPORTIONALLY TO THE MATRIX, INSTEAD !  */
	rsb_coo_idx_t coresrowf[RSB_CONST_MAX_SUPPORTED_CORES];	/*  first locked row, for each thread */
	rsb_coo_idx_t coresrowl[RSB_CONST_MAX_SUPPORTED_CORES];	/*  last  locked row, for each thread */
	rsb_coo_idx_t corescolf[RSB_CONST_MAX_SUPPORTED_CORES];	/*  first locked col, for each thread */
	rsb_coo_idx_t corescoll[RSB_CONST_MAX_SUPPORTED_CORES];	/*  last  locked col, for each thread */
	rsb_bitmap_data_t * bmap;	/* done matrices bitmap */
#if RSB__MAX_BITMAP_SUBMS_ON_STACK > 0
	rsb_bitmap_data_t bos[RSB_BYTES_PER_BITVECTOR(RSB__MAX_BITMAP_SUBMS_ON_STACK)];	/* bmap on stack; note that this makes the struct not shallow copyable  */
#endif
	rsb_submatrix_idx_t subms;	/* all matrices count */
	rsb_submatrix_idx_t dm;	/* done matrices count */
	rsb_submatrix_idx_t dr;	/* last done row */
	rsb_int_t nt;				/* number of threads */
	rsb_bool_t want_symlock;	/* symmetrical lock -- will lock both row and column region of output vector */
	rsb_bool_t want_fake_lock;	/* fake lock -- will allow concurrent writes (debug only) */
};

/*!
 * \ingroup gr_internals
 * \brief An internal, helper structure.
 */
struct rsb_bti_lock_struct
{
	rsb_coo_idx_t mvleaves;	/* maximal vertical leaves (>=itl) (2**(nlevels)) */
	rsb_coo_idx_t nlevels;	/* number of subdivisions  */
	rsb_coo_idx_t bsz;		/* (=2*mvleaves-1)*/
	rsb_coo_idx_t itl;		/* lock interval total length (e.g.: matrix dimension) */
	rsb_bitmap_data_t * bmap;	/* done intervals bitmap */
	rsb_bitmap_data_t * tmap;	/* tainted intervals bitmap */
};

/*!
 * \ingroup gr_internals
 * \brief An internal, helper structure.
 */
struct rsb_mv_lock_t
{
	/** 
	 * NEW: EXPERIMENTAL
	 * */
	struct rsb_rows_lock_struct_t olock;				/* output vector lock */
	struct rsb_bti_lock_struct locks[RSB_CONST_MAX_SUPPORTED_TEMPORARY_VECTORS];	/* it has no sense to have more locks than cores */
	size_t el_size;							/* numerical element size */
	rsb_type_t typecode;						/* type code */
	rsb_coo_idx_t nv;						/* number of vectors  */
	rsb_char_t * mv[RSB_CONST_MAX_SUPPORTED_TEMPORARY_VECTORS];		/* multiple vectors */
	rsb_char_t * ov;						/* master (output) vector */
	rsb_coo_idx_t itl;						/* interval total length */
	rsb_submatrix_idx_t last_subm[RSB_CONST_MAX_SUPPORTED_CORES];	/* last (tried unsuccessfully) matrix, per thread */
	rsb_coo_idx_t   in[RSB_CONST_MAX_SUPPORTED_CORES];		/* interval index, non transposed */
	rsb_coo_idx_t   it[RSB_CONST_MAX_SUPPORTED_CORES];		/* interval index, transposed */
	rsb_coo_idx_t   incov;					/* FIXME: NEW */
	rsb_trans_t	transA;						/* FIXME: NEW */
/*	rsb_bitmap_data_t ir[RSB_WORDS_PER_BITVECTOR(RSB_CONST_MAX_SUPPORTED_CORES)];	*/	/* is reducing ? */
};

#define RSB_WANT_SPMV_WITH_REDUCE 0

#if !RSB_WANT_SPMV_WITH_REDUCE
#define RSB_BOOL_ALMOST_TRUE 2 /* :) */
#define rsb_spmv_lock_struct_t rsb_rows_lock_struct_t
#define rsb_do_spmv_lock_init(LOCK,NT,SUMBS,MATRIX,OPFLAGS,TRANSA,OV,IO) rsb__do_lock_init(LOCK,NT,SUMBS,MATRIX,OPFLAGS)
#define rsb_do_spmv_lock_free(LOCK) rsb__do_lock_free(LOCK)
#define rsb_do_spmv_lock_release(LOCK,THID,OV) rsb__do_lock_release(LOCK,THID)
#define rsb_do_spmv_lock_get(LOCK,THID,ROFF,M,COFF,K,SUBM,TRANSA,OV,OI) rsb__do_lock_get(LOCK,THID,ROFF,M,COFF,K,SUBM,TRANSA)
#define RSB_DO_SPMV_LOCK_DM(LOCK) ((LOCK).dm)
#define RSB_DO_SPMV_LOCK_DM_INC(LOCK) ((LOCK).dm)++
#else
#define RSB_BOOL_ALMOST_TRUE 2 /* :) */
#define rsb_spmv_lock_struct_t rsb_mv_lock_t
#define rsb_do_spmv_lock_init(LOCK,NT,SUMBS,MATRIX,OPFLAGS,TRANSA,OV,IO) rsb__do_mv_lock_init(LOCK,NT,SUMBS,MATRIX,OPFLAGS,TRANSA,OV,IO)
#define rsb_do_spmv_lock_free(LOCK) rsb__do_mv_lock_free(LOCK)
#define rsb_do_spmv_lock_release(LOCK,THID,OV) rsb__do_mv_lock_release(LOCK,THID,OV)
#define rsb_do_spmv_lock_get(LOCK,THID,ROFF,M,COFF,K,SUBM,TRANSA,OV,OI) rsb__do_mv_lock_get(LOCK,THID,ROFF,M,COFF,K,SUBM,TRANSA,OV,OI)
#define RSB_DO_SPMV_LOCK_DM(LOCK) ((LOCK).olock.dm)
#define RSB_DO_SPMV_LOCK_DM_INC(LOCK) ((LOCK).olock.dm)++
#endif /* RSB_WANT_SPMV_WITH_REDUCE */

rsb_err_t rsb__do_mv_lock_init(struct rsb_mv_lock_t *lock, rsb_int_t num_threads, rsb_submatrix_idx_t subms, const struct rsb_mtx_t * mtxAp, enum rsb_op_flags_t op_flags, rsb_trans_t transA, rsb_char_t * ov, rsb_coo_idx_t incov);
rsb_err_t rsb__do_mv_lock_free(struct rsb_mv_lock_t *lock);
rsb_err_t rsb__do_mv_lock_release(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_char_t *ov);
rsb_bool_t rsb__do_mv_lock_get(struct rsb_mv_lock_t *lock ,rsb_thr_t th_id, rsb_coo_idx_t roff, rsb_coo_idx_t m, rsb_coo_idx_t coff, rsb_coo_idx_t k, rsb_submatrix_idx_t subm, rsb_trans_t transA, rsb_char_t **ov, rsb_coo_idx_t *incov);
rsb_err_t rsb__do_pick_candidate_interval_for_reduce(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_char_t ** ov, rsb_coo_idx_t * roff, rsb_coo_idx_t * m);
rsb_err_t rsb__do_release_candidate_interval_for_reduce(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_char_t *ov, rsb_coo_idx_t roff, rsb_coo_idx_t m);

rsb_bool_t rsb__do_lock_release(struct rsb_rows_lock_struct_t *lock, rsb_thr_t th_id);
rsb_bool_t rsb__do_lock_get(struct rsb_rows_lock_struct_t *lock,rsb_thr_t th_id, rsb_coo_idx_t roff, rsb_coo_idx_t m, rsb_coo_idx_t coff, rsb_coo_idx_t k, rsb_submatrix_idx_t subm, rsb_trans_t transA);
rsb_err_t rsb__do_lock_init(struct rsb_rows_lock_struct_t *lock, rsb_int_t num_threads, rsb_submatrix_idx_t subms, const struct rsb_mtx_t * mtxAp, enum rsb_op_flags_t op_flags);
rsb_err_t rsb__do_lock_free(struct rsb_rows_lock_struct_t *lock);
#if 0
rsb_err_t rsb__do_lock_test(void);
#endif

#endif /* RSB_LOCK_H_INCLUDED */

/* @endcond */
