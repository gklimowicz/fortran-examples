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
 * This source file contains locks for sparse recursive multicore operations.
 * */
#include "rsb_lock.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

/*
 TODO: one shall reduce the external interface, e.g. to a single rsb__lock function.
*/

rsb_bool_t rsb__do_lock_release(struct rsb_rows_lock_struct_t *lock, rsb_thr_t th_id)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	if(RSB__TRSV_OUT_)RSB_INFO("thread %zd releases  %zd %zd\n",(rsb_printf_int_t)th_id,(rsb_printf_int_t)lock->coresrowf[th_id],(rsb_printf_int_t)lock->coresrowl[th_id]);
	lock->corescoll[th_id]=RSB_MARKER_COO_VALUE;
	lock->corescolf[th_id]=RSB_MARKER_COO_VALUE;
	lock->coresrowl[th_id]=RSB_MARKER_COO_VALUE;
	lock->coresrowf[th_id]=RSB_MARKER_COO_VALUE;
	return RSB_BOOL_TRUE;
}

#if RSB_WANT_DO_LOCK_TEST
static RSB_INLINE rsb_bool_t rsb_do_lock_check_if_matrix_done(const struct rsb_rows_lock_struct_t *lock, rsb_submatrix_idx_t subm)
{
	/**
	 * 	\ingroup gr_internals
	 *  */
	if(RSB_BITVECTOR_GET(lock->bmap,lock->subms,subm))
		return RSB_BOOL_TRUE;
	else
		return RSB_BOOL_FALSE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

static RSB_INLINE rsb_bool_t rsb_do_lock_check_interval(const struct rsb_rows_lock_struct_t *lock, rsb_thr_t th_id, rsb_coo_idx_t roff, rsb_coo_idx_t m, rsb_coo_idx_t coff, rsb_coo_idx_t k, rsb_trans_t transA)
{
	/**
	 * 	\ingroup gr_internals
	 *  */
	rsb_thr_t tn;
	rsb_bool_t want_both=(lock->want_symlock == RSB_BOOL_TRUE);

	if(want_both)
	{
		for(tn=0;tn<lock->nt; ++tn)
		if( tn!=th_id && (
                           ((lock->coresrowf[tn] >= roff) && (lock->coresrowf[tn] < roff+m))
			|| ((lock->coresrowf[tn] <= roff) && (lock->coresrowl[tn]+1 > roff))
			|| ((lock->corescolf[tn] >= coff) && (lock->corescolf[tn] < coff+k))
			|| ((lock->corescolf[tn] <= coff) && (lock->corescoll[tn]+1 > coff))

                        || ((lock->coresrowf[tn] >= coff) && (lock->coresrowf[tn] < coff+k))
			|| ((lock->coresrowf[tn] <= coff) && (lock->coresrowl[tn]+1 > coff))
			|| ((lock->corescolf[tn] >= roff) && (lock->corescolf[tn] < roff+m))
			|| ((lock->corescolf[tn] <= roff) && (lock->corescoll[tn]+1 > roff))
			))
		{
			if(RSB__TRSV_OUT_)RSB_INFO("%zd %zd blocks %zd %zd\n",(rsb_printf_int_t)lock->coresrowf[tn],(rsb_printf_int_t)lock->coresrowl[tn],(rsb_printf_int_t)roff,(rsb_printf_int_t)m);
			goto l_false;
		}
	}
	else
	{
		if((RSB_DOES_NOT_TRANSPOSE(transA)) || want_both)
		for(tn=0;tn<lock->nt; ++tn)
		if( tn!=th_id
			&& (((lock->coresrowf[tn] >= roff) && (lock->coresrowf[tn] < roff+m))
			|| ((lock->coresrowf[tn] <= roff) && (lock->coresrowl[tn]+1 > roff))))
		{
			if(RSB__TRSV_OUT_)RSB_INFO("%zd %zd blocks %zd %zd\n",(rsb_printf_int_t)lock->coresrowf[tn],(rsb_printf_int_t)lock->coresrowl[tn],(rsb_printf_int_t)roff,(rsb_printf_int_t)m);
				goto l_false;
		}
	
		if(RSB_DOES_TRANSPOSE(transA) || want_both)
		for(tn=0;tn<lock->nt; ++tn)
		if( tn!=th_id
			&& (((lock->corescolf[tn] >= coff) && (lock->corescolf[tn] < coff+k))
			|| ((lock->corescolf[tn] <= coff) && (lock->corescoll[tn]+1 > coff))))
		{
			if(RSB__TRSV_OUT_)RSB_INFO("%zd %zd blocks %zd %zd\n",(rsb_printf_int_t)lock->coresrowf[tn],(rsb_printf_int_t)lock->coresrowl[tn],(rsb_printf_int_t)coff,(rsb_printf_int_t)k);
			goto l_false;
		}
	}
	return RSB_BOOL_TRUE;
l_false:
	return RSB_BOOL_FALSE;
}

	/* sets only the interval info for a given thread */
#define RSB_DO_LOCK_INTERVALS(LOCK,TH_ID,R0,R,C0,C) \
	(LOCK)->coresrowf[(TH_ID)]=(R0), (LOCK)->coresrowl[(TH_ID)]=(R0)+((R)-1), \
	(LOCK)->corescolf[(TH_ID)]=(C0), (LOCK)->corescoll[(TH_ID)]=(C0)+((C)-1)

#define RSB_DO_LOCK_INTERVAL(LOCK,TH_ID,R0,R) \
	(LOCK)->coresrowf[(TH_ID)]=(R0), (LOCK)->coresrowl[(TH_ID)]=(R0)+((R)-1), \
	(LOCK)->corescolf[(TH_ID)]=(R0), (LOCK)->corescoll[(TH_ID)]=(R0)+((R)-1)	/* FIXME: is there a reason for redundance ? */

/* FIXME: actually, this is the interval +1  */
#define RSB_GET_LOCK_INTERVAL_W(LOCK,TH_ID,R0,R1) \
	(R0)=(LOCK)->coresrowf[(TH_ID)], (R1)=(LOCK)->coresrowl[(TH_ID)]+1
#define RSB_GET_LOCK_INTERVAL_L(LOCK,TH_ID,R0,R1) \
	(R0)=(LOCK)->coresrolf[(TH_ID)], (R1)=(LOCK)->coresroll[(TH_ID)]+1

#if 0
#define RSB_GET_LOCK_INTERVALS(LOCK,TH_ID,R0,R,C0,C) \
	(R0)=(LOCK)->coresrowf[(TH_ID)], \
	(C0)=(LOCK)->corescolf[(TH_ID)], \
	(R)=(LOCK)->coresrowl[(TH_ID)]-(R0)+1, \
	(C)=(LOCK)->corescoll[(TH_ID)]-(C0)+1
#endif

rsb_bool_t rsb__do_lock_get(struct rsb_rows_lock_struct_t *lock, rsb_thr_t th_id, rsb_coo_idx_t roff, rsb_coo_idx_t m, rsb_coo_idx_t coff, rsb_coo_idx_t k, rsb_submatrix_idx_t subm, rsb_trans_t transA)
{
	/**
	 * 	\ingroup gr_internals
	 *  */
#if 0
	if(th_id)
	if(RSB__TRSV_OUT_)RSB_INFO("blocked by %p %d @ %d .. %d\n",lock->bmap,lock->subms,th_id,subm);
#endif
	
	if(RSB_BITVECTOR_GET(lock->bmap,lock->subms,subm))
		goto l_false;

	if(lock->want_fake_lock == RSB_BOOL_TRUE)
		goto l_true;	/* debug only : no locked rows check */

	if(!rsb_do_lock_check_interval(lock,th_id,roff,m,coff,k,transA))
		goto l_false;

	RSB_DO_LOCK_INTERVALS(lock,th_id,roff,m,coff,k);

	if(RSB__TRSV_OUT_)RSB_INFO("thread %zd locks  %zd %zd with matrix %zd\n",(rsb_printf_int_t)th_id,(rsb_printf_int_t)lock->coresrowf[th_id],(rsb_printf_int_t)lock->coresrowl[th_id],(rsb_printf_int_t)subm);
l_true:
	/* 
	 * WARNING : this does not mean that the matrix is 'done'.
	 * It only means that the matrix is now assigned to some core, and it will be processed soon.
	 * The guarantee that the matrix is done will be given us only by the lock-if this matrix
	 * is marked AND its row (or column) interval is free, then the matrix is done (in SPSV/SPMV).
	 * */
	RSB_BITVECTOR_SET(lock->bmap,lock->subms,subm);
	return RSB_BOOL_TRUE;
l_false:
	return RSB_BOOL_FALSE;
}

rsb_err_t rsb__do_lock_init(struct rsb_rows_lock_struct_t *lock, rsb_int_t num_threads, rsb_submatrix_idx_t subms, const struct rsb_mtx_t * mtxAp, enum rsb_op_flags_t op_flags)
{
	/** 
	 * 	\ingroup gr_internals
	 *
	 * 	Initializes lock data.
	 * 	Up to RSB__MAX_BITMAP_SUBMS_ON_STACK it won't allocate the bitmap but use own.
	 * 	So this excludes shallow copy.
	 * 	If RSB__MAX_BITMAP_SUBMS_ON_STACK is zero, shallow copy OK and will always allocate the bitmap.
	 * */
	rsb_int tn;

	RSB_ASSERT(mtxAp);
	RSB_ASSERT(lock);
	RSB_BZERO_P(lock);

	lock->nt=num_threads;
	for(tn=0;tn<RSB_CONST_MAX_SUPPORTED_CORES; ++tn)
		lock->corescolf[tn]=RSB_MARKER_COO_VALUE, lock->corescoll[tn]=RSB_MARKER_COO_VALUE,
		lock->coresrowf[tn]=RSB_MARKER_COO_VALUE, lock->coresrowl[tn]=RSB_MARKER_COO_VALUE;
	lock->dm=0;
	lock->subms=subms;
	lock->want_symlock = rsb__is_not_unsymmetric(mtxAp);
	lock->want_fake_lock=(op_flags == RSB_OP_FLAG_FAKE_LOCK);

#if RSB__MAX_BITMAP_SUBMS_ON_STACK > 0
	if ( subms <= RSB__MAX_BITMAP_SUBMS_ON_STACK )
		lock->bmap = (rsb_bitmap_data_t*) & (lock->bos);
	else
#endif
		lock->bmap = rsb__allocate_bitvector(subms);

	return (lock->bmap!=NULL)?RSB_ERR_NO_ERROR:RSB_ERR_ENOMEM;
}

rsb_err_t rsb__do_lock_free(struct rsb_rows_lock_struct_t *lock)
{
	/** 
	 * 	\ingroup gr_internals
	 * */
	RSB_ASSERT(lock);
	RSB_ASSERT(lock->subms > 0);

#if RSB__MAX_BITMAP_SUBMS_ON_STACK > 0
	if( lock->subms > RSB__MAX_BITMAP_SUBMS_ON_STACK )
#endif
		RSB_CONDITIONAL_FREE(lock->bmap);

	return RSB_ERR_NO_ERROR;
}

/*  BEGIN EXPERIMENTAL CODE */

#if RSB_WANT_DO_LOCK_TEST
size_t static rsb_do_log2(size_t n)
{
	/*!
	 * \ingroup gr_internals
	 * FIXME : document this
	 */
	size_t res = 0;

	while(n /= 2)
		++res;
	return res;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#define RSB_MULTINT_BY_TWO(X)   ((X)<<1)	/* FIXME: this is not portable */
#if RSB_WANT_DO_LOCK_TEST
#define RSB_UPPER_BOUNDING_LOG2(X) (rsb_do_log2(rsb__nearest_power_of_two(X)))
#endif
#define RSB_LOUD_BTILS_TESTING 0 /*  */
#define RSB_LOUD_MVL_TESTING 0   /* multivector lock   */
#define RSB_LOUD_MVR_TESTING 0   /* multivector reduce */
#define RSB_INHIBIT_MULTIVECTOR 1   /* multivector reduce */
#define RSB_INHIBIT_REDUCE 0   /* multivector reduce */

#if RSB_WANT_DO_LOCK_TEST
static rsb_err_t rsb_do_btils_init(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t itl, rsb_coo_idx_t nlevels)
{
	/**
	 * 	\ingroup gr_internals
	 * 	Initializes a lock structure.
	 * 	The input structure shall be freshly instantiated or freed.
	 * 	In case of error, it is safe but not required to call rsb_do_btils_free() to free it.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!lock || nlevels<0)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_BZERO_P(lock);
	lock->bmap=NULL;
	lock->nlevels=nlevels;
	lock->itl=itl;
	lock->mvleaves = RSB_POWER_OF_2(nlevels);
	lock->bsz=(2*lock->mvleaves-1);
	/* FIXME: need a check on nlevels */
	lock->bmap = rsb__allocate_bitvector(lock->bsz);
	lock->tmap = rsb__allocate_bitvector(lock->bsz);
	if(!lock->bmap || !lock->tmap)
	{
		RSB_CONDITIONAL_FREE(lock->bmap);
		RSB_CONDITIONAL_FREE(lock->tmap);
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
err:
	return errval;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_err_t rsb_do_btils_free(struct rsb_bti_lock_struct * lock)
{
	/** 
	 * 	\ingroup gr_internals
	 * 	Frees a lock structure.
	 * 	The input structure shall be initialized with success.
	 * */
	if(!lock)
		return RSB_ERR_BADARGS;
	RSB_CONDITIONAL_FREE(lock->bmap);
	RSB_CONDITIONAL_FREE(lock->tmap);
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_coo_idx_t rsb_do_rindex_to_lindex(rsb_coo_idx_t r0, rsb_coo_idx_t r1, rsb_coo_idx_t n, rsb_coo_idx_t nlevels)
{
	/** 
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t l0=0,l1=0,doffset=1,offset=0;
	rsb_coo_idx_t n0=n,n1=n;
	rsb_int i,delta=0;
	if(nlevels<1)
	{
		return 0;
	}
	if(r1==n1)
		l1=2,r1=0;
	for(i=0;i<nlevels;++i)
	{
		rsb_coo_idx_t m0=RSB_MIDDLE(n0);
		rsb_coo_idx_t m1=RSB_MIDDLE(n1);

		if(r0>=m0)
			r0-=m0,++l0,n0-=m0;
		else
			n0=m0;
		if(r1>=m1)
			r1-=m1,++l1,n1-=m1;
		else
			n1=m1;	

		if(i<nlevels-1)
			l0*=2,l1*=2;
	}
#if 0
  	RSB_INFO("%d!\n",l1-l0);
#endif
	delta=l1-l0;
	l0=l0/(l1-l0);
	offset = RSB_POWER_OF_2(nlevels)-1;
	doffset = RSB_POWER_OF_2(nlevels-1);
	for( ;delta>1;delta/=2,doffset/=2)
		offset-=doffset;
		
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("@ bit %d + %d\n",l0,offset);
	return offset+l0;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TES
static rsb_bool_t rsb_do_btils_lock_update_tmap(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t i)
{
	rsb_coo_idx_t iu,il,ii;
	/* we taint the vector: after a lock, it will mark the interval as tainted 
	 * (as opposed to the cases where an untainted vector is unlocked (e.g.: after a reduce))  */
	RSB_BITVECTOR_SET(lock->tmap,lock->bsz,i);
	/* did already any ancestor taint all way up ? */
	for(iu=(i-1)/2;iu>0;iu=(iu-1)/2)
	{
		/* TODO: could speed up a little while by inverting the visit order */
		if(RSB_LOUD_BTILS_TESTING)
			RSB_INFO("updating tmap\n");
		if(RSB_BITVECTOR_GET(lock->tmap,lock->bsz,iu))
			goto l_done;
	}
	/* no ancestor tainted all way up */
	RSB_ASSERT(iu==0);
	if(RSB_BITVECTOR_GET(lock->tmap,lock->bsz,iu))
		goto l_done;
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("reducing taint map:\n"),rsb__do_dump_bitmap(lock->tmap,1,lock->bsz),RSB_INFO("\n");

	/* we look for neighbor leaves needing collapse, at any upper level  */
	while(i>0)
	{
		il=2*((i-1)/2  )+1;
		iu=2*((i-1)/2+1)+1;
		for(ii=il;ii<iu;++ii)
			if(!RSB_BITVECTOR_GET(lock->tmap,lock->bsz,ii))
				goto skip;/* The sibling interval is not tainted: we may stop merging here.
				    Pay attention: some descendant of ours may have still its bit set despite 
				    at this level we are done with merging: thus that bit would be obsolete,
				    and it could be possible for it to remain.
				    This does not cause harm, so we don't force bit-clear to lower nodes, here.
			           */
		/* merge the current subtree */
		for(ii=il;ii<iu;++ii)
			RSB_BITVECTOR_UNSET(lock->tmap,lock->bsz,ii);
		/* collapse to the upper node */
		i=(i-1)/2;
		RSB_BITVECTOR_SET(lock->tmap,lock->bsz,i);
		continue;
skip:
		i=(i-1)/2;
	}
	/* the taint map is done */
	
l_done:
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("taint map:\n"),rsb__do_dump_bitmap(lock->tmap,1,lock->bsz),RSB_INFO("\n");
	return RSB_BOOL_TRUE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static RSB_INLINE rsb_bool_t rsb_do_btils_lock_probe_inner(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t i)
{
	/**
	 * */
	rsb_coo_idx_t iu,il;
	rsb_coo_idx_t ili=2,ilii;
	RSB_ASSERT(lock);
	RSB_ASSERT(i>=0);

	if(RSB_BITVECTOR_GET(lock->bmap,lock->bsz,i))
		goto l_false;

	for(iu=(i-1)/2;iu>0;iu=(iu-1)/2)
	{
#if 0
		if(1) RSB_INFO("checking bit .. %d:%d\n",iu,RSB_BOOL_TRUE == RSB_BITVECTOR_GET(lock->bmap,lock->bsz,iu));
#endif
		if(RSB_BITVECTOR_GET(lock->bmap,lock->bsz,iu))
			goto l_false;
	}
	if(RSB_BITVECTOR_GET(lock->bmap,lock->bsz,iu))/* iu==0 */
		goto l_false;
	for(il=2*i+1;il<lock->bsz;il=2*il+1,ili*=2)
	{
		for(ilii=0;ilii<ili;++ilii)
			if(RSB_BITVECTOR_GET(lock->bmap,lock->bsz,il+ilii))
				goto l_false;
	}
	
	return RSB_BOOL_TRUE;
l_false:
	return RSB_BOOL_FALSE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_bool_t rsb_do_btils_lock_probe(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t m0, rsb_coo_idx_t m1, rsb_coo_idx_t *ip)
{
	/**
	 * */
	rsb_coo_idx_t i;
	RSB_ASSERT(lock);
	RSB_ASSERT(ip);

	i = rsb_do_rindex_to_lindex(m0,m1,lock->itl,lock->nlevels);
	if(!rsb_do_btils_lock_probe_inner(lock,i))
		goto l_false;
	*ip=i;
	return RSB_BOOL_TRUE;
l_false:
	return RSB_BOOL_FALSE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_bool_t rsb_do_btils_lock_get_sym(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t m0, rsb_coo_idx_t m1, rsb_coo_idx_t k0, rsb_coo_idx_t k1, rsb_trans_t transA, rsb_coo_idx_t *ip, rsb_coo_idx_t *jp)
{
	/** 
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t i,j;
	if(!rsb_do_btils_lock_probe(lock,m0,m1,&i))
		goto l_false;
	j=i;
	if((m0!=k0) && (m1!=k1))
		if(!rsb_do_btils_lock_probe(lock,k0,k1,&j))
			goto l_false;

	RSB_BITVECTOR_SET(lock->bmap,lock->bsz,i);
	if(i!=j)
		RSB_BITVECTOR_SET(lock->bmap,lock->bsz,j);
	if(RSB_LOUD_BTILS_TESTING)
		rsb__do_dump_bitmap(lock->bmap,1,lock->bsz),RSB_INFO(" (%d)\n",lock->bsz);

	/* we're going to lock up to i */
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("(nlev=%d)(%d .. %d) -> %d ok\n",lock->nlevels,m0,m1,i);

	/* TODO: update the taint vector accordingly */
	rsb_do_btils_lock_update_tmap(lock,i);
	if(i!=j)
		rsb_do_btils_lock_update_tmap(lock,j);

	if(RSB_DOES_TRANSPOSE(transA))
		RSB_SWAP(rsb_coo_idx_t,i,j);

	*ip=i;
	*jp=j;

	return RSB_BOOL_TRUE;
l_false:
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("(nlev=%d)(%d .. %d) -> (%d %d) busy \n",lock->nlevels,m0,m1,i,j);
	return RSB_BOOL_FALSE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_bool_t rsb_do_btils_lock_get(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t m0, rsb_coo_idx_t m1, rsb_trans_t transA, rsb_coo_idx_t *ip, rsb_coo_idx_t *jp)
{
	/** 
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t i = RSB_MARKER_COO_VALUE,j = RSB_MARKER_COO_VALUE;
	if(!rsb_do_btils_lock_probe(lock,m0,m1,&i))
		goto l_false;

	RSB_BITVECTOR_SET(lock->bmap,lock->bsz,i);
	if(RSB_LOUD_BTILS_TESTING)
		rsb__do_dump_bitmap(lock->bmap,1,lock->bsz),RSB_INFO(" (%d)\n",lock->bsz);

	/* we're going to lock up to i */
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("(nlev=%d)(%d .. %d) -> %d ok\n",lock->nlevels,m0,m1,i);

	/* TODO: update the taint vector accordingly */
	rsb_do_btils_lock_update_tmap(lock,i);
	
	if(RSB_DOES_TRANSPOSE(transA))
		RSB_SWAP(rsb_coo_idx_t,i,j);
	
	*ip=i,*jp=j;

	return RSB_BOOL_TRUE;
l_false:
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("(nlev=%d)(%d .. %d) -> %d busy \n",lock->nlevels,m0,m1,i);
	return RSB_BOOL_FALSE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_err_t rsb_do_get_interval_info_from_btils_lock(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t i, rsb_coo_idx_t *m0p, rsb_coo_idx_t * m1p)
{
	/** 
	 * 	\ingroup gr_internals
	 * 	FIXME: unfinished
	 * */
	rsb_coo_idx_t m0=0,m1=lock->itl,h=lock->itl,iu,l=0,ii=0,pot=1, nl=0;

	for(iu=i;iu>0;iu=(iu-1)/2)
	{
		ii = RSB_MULTINT_BY_TWO(ii);
		if(RSB_IS_INTEGER_EVEN(iu))
			++ii;
		++nl;
	}

	for(l=0;l<nl;++l)
	{
		if(ii&pot)
			m0+=RSB_MIDDLE(h),
			h=h-RSB_MIDDLE(h);
		else
			m1-=h-RSB_MIDDLE(h),
			h = RSB_MIDDLE(h);
#if 0
		RSB_INFO("BIBO: ii=%d m0=%d, h=%d, pot=%d\n",ii,m0,h,pot);
#endif
		pot = RSB_MULTINT_BY_TWO(pot);
	}
	*m0p=m0;
	*m1p=m1;
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static void RSB_INLINE rsb__do_btils_lock_release_inner(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t i)
{
	/** 
	 * 	\ingroup gr_internals
	 * 	FIXME: does this call free one or two intervals ?
	 * */
	if(RSB_LOUD_BTILS_TESTING)
		rsb__do_dump_bitmap(lock->bmap,1,lock->bsz),RSB_INFO(" (%d)\n",lock->bsz);
	RSB_BITVECTOR_UNSET(lock->bmap,lock->bsz,i);
	if(RSB_LOUD_BTILS_TESTING)
	{
		rsb_coo_idx_t m0,m1;
		rsb_do_get_interval_info_from_btils_lock(lock,i,&m0,&m1);
		RSB_INFO("freeing (%d .. %d)\n",m0,m1),
		rsb__do_dump_bitmap(lock->bmap,1,lock->bsz),RSB_INFO(" (%d)\n",lock->bsz);
	}
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_err_t rsb_do_btils_lock_release(struct rsb_bti_lock_struct * lock, rsb_coo_idx_t m0, rsb_coo_idx_t m1)
{
	/** 
	 * 	\ingroup gr_internals
	 * 	FIXME: does this call free one or two intervals ?
	 * 	FIXME: deprecated
	 * */
	rsb_coo_idx_t i;

	if(!lock)
		return RSB_ERR_BADARGS;

	i = rsb_do_rindex_to_lindex(m0,m1,lock->itl,lock->nlevels);
	RSB_ASSERT(i>=0);
	rsb__do_btils_lock_release_inner(lock,i);
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#define RSB_MV_OFFSET(LOCK,INDEX,OFFSET) \
	((((rsb_char_t *)((LOCK)->mv[INDEX]))) +(LOCK)->el_size*(OFFSET))

#if RSB_WANT_DO_LOCK_TEST
static rsb_err_t rsb__do_mv_lock_release_single(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_char_t *ov)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t nvi = RSB_MARKER_COO_VALUE;

	/* in the case the locked vector was the master one */
	if(ov==lock->ov)
	{
		if(RSB_LOUD_MVL_TESTING)
			RSB_INFO("releasing master vector from thread %d\n",th_id);
		rsb__do_lock_release(&(lock->olock),th_id);
		goto ok;
	}

	/* in the case the locked vector was not the master one */
	for(nvi=0;nvi<lock->nv;++nvi)
		if( (ov >= RSB_MV_OFFSET(lock,nvi,0)) && (ov<RSB_MV_OFFSET(lock,nvi+1,0)))
		{
			struct rsb_bti_lock_struct * vlock=&(lock->locks[nvi]);
			/* we localized the vector. now we shall see if that interval was locked */
			if(RSB_LOUD_MVL_TESTING)
			{
				RSB_INFO("releasing vector %d from thread %d\n",nvi,th_id);
/*  			if(RSB_BITVECTOR_GET(vlock->bmap,vlock->bsz,rsb_do_rindex_to_lindex(roff,roff+m,vlock->itl,vlock->nlevels)))
				RSB_INFO("freeing interval %d .. %d on vector %d (thread %d)\n",roff,roff+m,nvi,th_id);
				else
				{RSB_INFO("guessed pointer !?\n");goto failure;}*/
			}
			if(lock->it[th_id]!=RSB_MARKER_COO_VALUE)
			{
				if(RSB_LOUD_MVL_TESTING) RSB_INFO("releasing inner\n");
				rsb__do_btils_lock_release_inner(vlock,lock->it[th_id]);
			}
			if(lock->in[th_id]!=RSB_MARKER_COO_VALUE)
			{
				if(RSB_LOUD_MVL_TESTING) RSB_INFO("releasing inner\n");
				rsb__do_btils_lock_release_inner(vlock,lock->in[th_id]);
			}
			lock->it[th_id]=RSB_MARKER_COO_VALUE;
			lock->in[th_id]=RSB_MARKER_COO_VALUE;
			goto ok;
		}
#if 0
failure:
	if(RSB_LOUD_MVL_TESTING)
		RSB_INFO("did not find a vector to release for thread %d\n",th_id);
	return RSB_ERR_GENERIC_ERROR;
#endif
ok:
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
rsb_err_t rsb__do_mv_lock_release(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_char_t *ov)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
#if RSB_LOUD_MVL_TESTING
#if 0
	rsb_coo_idx_t roff = RSB_MARKER_COO_VALUE,m = RSB_MARKER_COO_VALUE,coff = RSB_MARKER_COO_VALUE,k = RSB_MARKER_COO_VALUE;
#endif
#endif /* RSB_LOUD_MVL_TESTING */
	rsb_bool_t is_reduce_only = RSB_BOOL_TRUE;/* FIXME: ?? */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_ASSERT(lock);
	RSB_ASSERT(ov);
	RSB_ASSERT(th_id>=0);
#if 0
	RSB_GET_LOCK_INTERVALS(&(lock->olock),th_id,roff,m,coff,k);
#endif
	errval = rsb__do_mv_lock_release_single(lock,th_id,ov);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(failure,RSB_ERRM_ES);
	}
	if(is_reduce_only)
		goto reduce_ok;
	else
	       	goto ok;
reduce_ok:
#if 0
	RSB_ASSERT(roff>=0); RSB_ASSERT(m>=0); RSB_ASSERT(th_id>=0);
	if(RSB_LOUD_MVL_TESTING)
		RSB_INFO("freeing interval %d .. %d on master vector (thread %d) \n",roff,roff+m,th_id);
	/* it may still be the master vector, or a random pointer :) */
	rsb__do_lock_release(&(lock->olock),th_id);
#endif
	goto ok;
failure:
#if 0
  	RSB_ASSERT(roff>=0); RSB_ASSERT(m>=0); RSB_ASSERT(th_id>=0);
  	if(RSB_LOUD_MVL_TESTING)
  		RSB_INFO("not freeing interval %d .. %d\n",roff,roff+m);
#endif
ok:
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_bool_t rsb_do_is_bitmap_blank(rsb_bitmap_data_t *bmap, rsb_coo_idx_t r, rsb_coo_idx_t c)
{
	/* *
	 * 	\ingroup gr_internals
	 * 	FIXME: new, untested
	 * */
	size_t bs = RSB_WORDS_PER_BITMAP(r,c);
	rsb_coo_idx_t i;
	for(i=0;i<bs;++i)
	{
		if(bmap[i])
			return RSB_BOOL_FALSE;
	}
	return RSB_BOOL_TRUE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_bool_t rsb_do_is_bitvector_blank(rsb_bitmap_data_t *bmap, rsb_coo_idx_t c)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	return rsb_do_is_bitmap_blank(bmap,1,c);
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_bool_t rsb__do_mv_lock_is_used(struct rsb_mv_lock_t *lock)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t nvi;
	for(nvi=0;nvi<lock->nv;++nvi)
	{
		struct rsb_bti_lock_struct * vlock=&(lock->locks[nvi]);
		if(!rsb_do_is_bitvector_blank(vlock->bmap,vlock->bsz))
			return RSB_BOOL_TRUE;
	}
	return RSB_BOOL_FALSE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
static rsb_bool_t rsb__do_mv_lock_is_tainted(struct rsb_mv_lock_t *lock)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t nvi;
	for(nvi=0;nvi<lock->nv;++nvi)
	{
		struct rsb_bti_lock_struct * vlock=&(lock->locks[nvi]);
		if(!rsb_do_is_bitvector_blank(vlock->tmap,vlock->bsz))
			return RSB_BOOL_TRUE;
	}
	return RSB_BOOL_FALSE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
rsb_err_t rsb__do_mv_lock_free(struct rsb_mv_lock_t *lock)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t nvi;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if !RSB_INHIBIT_MULTIVECTOR
	rsb_bool_t tainted = rsb__do_mv_lock_is_tainted(lock);
#endif /* RSB_INHIBIT_MULTIVECTOR */
	if(RSB_LOUD_MVL_TESTING)
	{
		if(lock->nv)
			RSB_INFO("taint maps:\n");
		for(nvi=0;nvi<lock->nv;++nvi)
			rsb__do_dump_bitmap(lock->locks[nvi].tmap,1,lock->locks[nvi].bsz),RSB_INFO("\n");
	}

	/* FIXME: TODO: reduce all of the vectors, here. */
#if !RSB_INHIBIT_MULTIVECTOR
#if RSB_INHIBIT_REDUCE
	/* no reduce. this will produce wrong results, of course */
#else /* RSB_INHIBIT_REDUCE */
	if(rsb__do_mv_lock_is_used(lock))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,"no vector should not be in use before reducing!");
	}
#if 0
	/* this approach is likely to be faster for high nnz/row cases */
	if(RSB_LOUD_MVR_TESTING)
		RSB_INFO("summing up (%d) vectors to the master (strided %d)\n",lock->nv,lock->incov);
	if(RSB_LOUD_MVR_TESTING)
		RSB_INFO("on master vector:\n"),
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_print_some_vector_stats(lock->ov,lock->typecode,lock->itl));
	for(nvi=0;nvi<lock->nv;++nvi)
	{
		if(RSB_LOUD_MVR_TESTING)
			RSB_INFO("on vector %d:\n",nvi),
			RSB_DO_ERROR_CUMULATE(errval,rsb__do_print_some_vector_stats(lock->mv[nvi],lock->typecode,lock->itl));
		rsb__vectors_left_sum_reduce_and_zero(lock->ov,lock->mv[nvi],lock->typecode,lock->itl,lock->incov,0);
	}
	if(RSB_LOUD_MVR_TESTING)
		RSB_INFO("\n");
#else
#if 1
	#pragma omp parallel shared(tainted) RSB_NTC
	{
		rsb_thr_t th_id = omp_get_thread_num();
		rsb_coo_idx_t oincy=lock->incov,rh,r0;
		rsb_char_t *ov=NULL;

		if(th_id>=lock->nv)
			goto skip;
		if(th_id >= rsb_global_session_handle.rsb_want_threads)
			goto skip;

		while(tainted)
		{

			if(RSB_LOUD_MVR_TESTING)
			{
				if(lock->nv)
					RSB_INFO("taint maps:\n");
				for(nvi=0;nvi<lock->nv;++nvi)
					rsb__do_dump_bitmap(lock->locks[nvi].tmap,1,lock->locks[nvi].bsz),RSB_INFO(" (%d)\n",rsb__do_mv_lock_is_tainted(lock));
				if(lock->nv)
					RSB_INFO("use maps:\n");
				for(nvi=0;nvi<lock->nv;++nvi)
					rsb__do_dump_bitmap(lock->locks[nvi].bmap,1,lock->locks[nvi].bsz),RSB_INFO(" (%d)\n",rsb__do_mv_lock_is_tainted(lock));
			}

			ov=lock->ov;
			#pragma omp critical (rsb_lock_crs)
			{ rsb__do_pick_candidate_interval_for_reduce(lock,th_id,&ov,&r0,&rh); }
	
			if(ov && ov!=lock->ov)
			{
				if(RSB_LOUD_MVR_TESTING)
					RSB_INFO("%d .. %d (incov = %d)\n",r0,rh,oincy);
				rsb__vectors_left_sum_reduce_and_zero(lock->ov,ov,lock->typecode,rh,oincy,r0);/*wrong ?*/
#if 0
				rsb__vectors_left_sum_reduce_and_zero(lock->ov,ov,lock->typecode,lock->itl,oincy,0);/*~works*/
				rsb__vectors_left_sum_reduce_and_zero(lock->ov,ov,lock->typecode,lock->itl,lock->incov,0);/*~works*/
#endif
	                     	#pragma omp critical (rsb_lock_crs)
	                   	{ rsb__do_release_candidate_interval_for_reduce(lock,th_id,ov,r0,rh); }
			}
			#pragma omp critical (rsb_lock_crs)
			{ tainted = rsb__do_mv_lock_is_tainted(lock); }
	}
skip:
		#pragma omp barrier
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	}
#else
	if(RSB_LOUD_MVR_TESTING)
	{
		if(lock->nv)
			RSB_INFO("taint maps:\n");
		for(nvi=0;nvi<lock->nv;++nvi)
			rsb__do_dump_bitmap(lock->locks[nvi].tmap,1,lock->locks[nvi].bsz),RSB_INFO("\n");
	}
	/* serial approach, for debugging purposes (very slow) ; it should be used to debug the rest */
	for(nvi=0;nvi<lock->nv;++nvi)
		rsb__vectors_left_sum_reduce_and_zero(lock->ov,lock->mv[nvi],lock->typecode,lock->itl,lock->incov,0);
#endif
#endif
#endif /* RSB_INHIBIT_REDUCE */
#endif /* RSB_INHIBIT_MULTIVECTOR */
	goto nosync;
nosync:
	if(!lock)
	{	
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	for(nvi=0;nvi<lock->nv;++nvi)
		RSB_CONDITIONAL_FREE(lock->mv[nvi]);

	for(;nvi>=0;--nvi)
		rsb_do_btils_free(&(lock->locks[nvi]));

	RSB_DO_ERROR_CUMULATE(errval,rsb__do_lock_free(&(lock->olock)));
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
rsb_err_t rsb__do_mv_lock_init(struct rsb_mv_lock_t *lock, rsb_int_t num_threads, rsb_submatrix_idx_t subms, const struct rsb_mtx_t * mtxAp, enum rsb_op_flags_t op_flags, rsb_trans_t transA, rsb_char_t * ov, rsb_coo_idx_t incov)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t nvi,nlevels;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_int_t th_id=0;
	rsb_coo_idx_t tn;

	if(!lock || !mtxAp)
		return RSB_ERR_BADARGS;

	RSB_BZERO_P(lock);
	errval = rsb__do_lock_init(&(lock->olock),num_threads,subms,mtxAp,op_flags);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err0,RSB_ERRM_ES);
	}
	/* FIXME: we need a policy for this */
#if RSB_INHIBIT_MULTIVECTOR
	lock->nv=0;	/* FIXME: for debugging purposes */
#else
#if 0
  	lock->nv = RSB_MIN(num_threads-1,((mtxAp->nnz+1)/(4*mtxAp->nr+1)));
  	lock->nv = RSB_MIN(num_threads-1,1);/*FIXME: temporary */
  	lock->nv = RSB_MIN(num_threads-1,rsb__submatrices(mtxAp));/*FIXME: temporary */
  	lock->nv=1;	/* FIXME: for debugging purposes */
#endif
	lock->nv = RSB_MIN(num_threads-1,mtxAp->all_leaf_matrices_n-1);/* FIXME: temporary */
#endif /* RSB_INHIBIT_MULTIVECTOR */
	if(RSB_LOUD_MVR_TESTING)
		RSB_INFO("Will use %d temporary vectors for %d threads\n",lock->nv,num_threads);

	RSB_ASSERT(lock->nv<RSB_CONST_MAX_SUPPORTED_CORES);
	lock->el_size=mtxAp->el_size;
	lock->typecode=mtxAp->typecode;
	lock->itl = rsb__do_get_rows_of(mtxAp,transA);
	lock->ov=ov;
	lock->incov=incov;
	lock->transA=transA;
	nlevels = rsb__get_recursive_matrix_depth(mtxAp);
	for(tn=0;tn<RSB_CONST_MAX_SUPPORTED_CORES; ++tn)
		lock->it[tn]=
		lock->in[tn]=RSB_MARKER_COO_VALUE;
	for(nvi=0;nvi<lock->nv;++nvi)
		if((errval = rsb_do_btils_init(&(lock->locks[nvi]),lock->itl,nlevels))!=RSB_ERR_NO_ERROR)
		{
			RSB_PERR_GOTO(err1,RSB_ERRM_ES);
		}
	/* time to allocate the temporary vectors */
	for(nvi=0;nvi<lock->nv;++nvi)
		if((lock->mv[nvi]=rsb__calloc(lock->el_size*lock->itl))==NULL)
		{
			RSB_PERR_GOTO(err2,RSB_ERRM_ES);
		}
	for(th_id=0;th_id<num_threads;++th_id)
		lock->last_subm[th_id]=RSB_SUBM_IDX_MARKER;
	/* the multivector lock is allocated. nice! */
	
	return RSB_ERR_NO_ERROR;
err2:
	for(nvi=0;nvi<lock->nv;++nvi)
		RSB_CONDITIONAL_FREE(lock->mv[nvi]);
err1:
	for(nvi=0;nvi<lock->nv;++nvi)
		rsb_do_btils_free(&(lock->locks[nvi]));
	RSB_DO_ERROR_CUMULATE(errval,rsb__do_lock_free(&(lock->olock)));
err0:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
rsb_bool_t rsb__do_mv_lock_get(struct rsb_mv_lock_t *lock ,rsb_thr_t th_id, rsb_coo_idx_t roff, rsb_coo_idx_t m, rsb_coo_idx_t coff, rsb_coo_idx_t k, rsb_submatrix_idx_t subm, rsb_trans_t transA, rsb_char_t **ov, rsb_coo_idx_t *incov)
{
	/* *
	 * 	\ingroup gr_internals
	 * */
	rsb_coo_idx_t nvi = RSB_MARKER_COO_VALUE;
	rsb_bool_t was_looping=(lock->last_subm[th_id]==subm);
	rsb_coo_idx_t i,j;
	if(!ov)
		return RSB_BOOL_FALSE;
	if(rsb_do_lock_check_if_matrix_done(&(lock->olock),subm))
	{
		if(RSB_LOUD_MVL_TESTING)
			RSB_INFO("matrix %d is already locked (for thread %d)\n",subm,th_id);
		/* if the thread was looping on this mtxAp, there's no reason to do so anymore (mtxAp locked or done) */
		if(was_looping)
			lock->last_subm[th_id]=RSB_SUBM_IDX_MARKER;

		return RSB_BOOL_FALSE;/* nothing to do: matrix done */
	}
	/* first, we try to get a lock on the master vector */
	if(rsb__do_lock_get(&(lock->olock),th_id,roff,m,coff,k,subm,transA))
	{
		if(RSB_LOUD_MVL_TESTING)
			RSB_INFO("locking matrix %d [%d...%d) to thread %d on master vector\n",subm,roff,roff+m,th_id);
		goto found;
	}
	/* if the master vector was not available, we check if this thread was in a loop on this matrix */
	if(!was_looping)
	{
		/* it was not looping on this submatrix */
		if(lock->last_subm[th_id]==RSB_SUBM_IDX_MARKER)
		{
			/* it was not looping at all; 
			 * now, if the thread will be back here with the value unchanged, the loop will be detected */
			lock->last_subm[th_id]=subm;
			if(RSB_LOUD_MVL_TESTING)
				RSB_INFO("not locking matrix %d to thread %d : waiting for a loop\n",subm,th_id);
			return RSB_BOOL_ALMOST_TRUE;
		}
		else
			;/*  the thread is looping on another submatrix: let it loop there, then */
		return RSB_BOOL_FALSE;
	}

	/* the thread was looping, and then it has the right to use a temporary vector (if any) */
	if(RSB_DOES_TRANSPOSE(transA))
	{ RSB_SWAP(rsb_coo_idx_t,k,m);RSB_SWAP(rsb_coo_idx_t,coff,roff); } /* FIXME: a dirty trick */
	if((lock->olock.want_symlock == RSB_BOOL_TRUE))
	{
		for(nvi=0;nvi<lock->nv;++nvi)
			if(rsb_do_btils_lock_get_sym(&(lock->locks[nvi]),roff,m+roff,coff,k+coff,transA,&i,&j))
			{
				lock->it[th_id]=j,lock->in[th_id]=i;
				goto found;
			}
	}
	else
	{
		for(nvi=0;nvi<lock->nv;++nvi)
			if(rsb_do_btils_lock_get(&(lock->locks[nvi]),roff,m+roff,transA,&i,&j))/* FIXME */
			{
				lock->in[th_id]=i,lock->it[th_id]=RSB_MARKER_COO_VALUE;
				goto found;
			}
	}
	/* TODO:  are we sure that pushing the thread for looping is always the best thing ? */
	/* TODO: implement here the task of picking up some vector and "reducing" it (not here, but in a "returned" signalling)! */
	return RSB_BOOL_FALSE;
found:
	/* found a temporary vector to perform computation */
	if(RSB_LOUD_MVL_TESTING)
		if(nvi != RSB_MARKER_COO_VALUE)
			RSB_INFO("locking interval %d .. %d on vector %d to thread %d\n",roff,roff+m,nvi,th_id);

	lock->last_subm[th_id]=RSB_SUBM_IDX_MARKER;
	if(nvi == RSB_MARKER_COO_VALUE)
		*ov=lock->ov,			/* we'll work on the master vector */
		*incov=lock->incov;			/* unchanged */
	else
		*ov = RSB_MV_OFFSET(lock,nvi,0),	/* we'll work on an auxiliary vector */
		*incov=1;
	RSB_BITVECTOR_SET(lock->olock.bmap,lock->olock.subms,subm);
	return RSB_BOOL_TRUE;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_SPMV_WITH_REDUCE
rsb_err_t rsb__do_release_candidate_interval_for_reduce(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_char_t *ov, rsb_coo_idx_t roff, rsb_coo_idx_t m)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* retrieve working interval */
	/* ... */
	rsb_coo_idx_t m0=0,m1=0;
	RSB_GET_LOCK_INTERVAL_W(&(lock->olock),th_id,m0,m1);
	rsb__do_lock_release(&(lock->olock),th_id);
	rsb__do_mv_lock_release(lock,th_id,ov);
	if(RSB_LOUD_MVR_TESTING)
		RSB_INFO("releasing reduce interval %d .. %d from thread %d\n",m0,m1,th_id);
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_SPMV_WITH_REDUCE */

#if RSB_WANT_SPMV_WITH_REDUCE
rsb_err_t rsb__do_pick_candidate_interval_for_reduce(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_char_t ** ov, rsb_coo_idx_t * roff, rsb_coo_idx_t * m)
{
	/*
		pick an interval which is free on both the master vector AND some vector v_i, and candidate for reducing
		we begin with the last temporary vector first
	*/
	rsb_coo_idx_t nvi;
	rsb_coo_idx_t i;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!lock)
		goto err;
	/* we start looking for a vector to reduce with the last one */
	for(nvi=lock->nv-1;nvi>=0;--nvi)
	{
		struct rsb_bti_lock_struct * vlock=&(lock->locks[nvi]);
		if(RSB_LOUD_MVR_TESTING)
			RSB_INFO("looking for tainted subvectors in temporary vector %d\n",nvi),
			RSB_INFO("taint map: "),rsb__do_dump_bitmap(vlock->tmap,1,vlock->bsz),RSB_INFO("\n"),
			RSB_INFO("use   map: "),rsb__do_dump_bitmap(vlock->bmap,1,vlock->bsz),RSB_INFO("\n");
		for(i=0;i<vlock->bsz;++i)
		{
			if(RSB_BITVECTOR_GET(vlock->tmap,vlock->bsz,i) && rsb_do_btils_lock_probe_inner(vlock,i))
			{
				/* let's see if the master vector has available subvector i (TODO) */
				/* first step is to obtain the bounds of the lock */
				rsb_coo_idx_t m0,m1;
				rsb_bool_t goir = RSB_BOOL_FALSE;

				errval = rsb_do_get_interval_info_from_btils_lock(vlock,i,&m0,&m1);
				if(RSB_LOUD_MVR_TESTING)
					RSB_INFO("temporary vector %d is tainted at interval %d\n",nvi,i);

				if(RSB_LOUD_MVR_TESTING)
					RSB_INFO("let's see if the master vector has available subvector %d at [%d .. %d]... ",nvi,m0,m1);
				goir = rsb_do_lock_check_interval(&(lock->olock),th_id,m0,m1-m0,m0,m1-m0,lock->transA);
				if(RSB_LOUD_MVR_TESTING)
					if(!goir)
						RSB_INFO("no\n");
				if(goir)
				{
					if(RSB_LOUD_MVR_TESTING)
						RSB_INFO("yes. will lock now [%d .. %d].\n",m0,m1);
					/* The interval is free on both subvectors.
					 *  We lock the interval on both vectors, then.
					 *  */
					/* mark the interval as not tainted anymore, but locked, on nvi  */
					RSB_BITVECTOR_SET(vlock->bmap,vlock->bsz,i);
					RSB_BITVECTOR_UNSET(vlock->tmap,vlock->bsz,i);
					/* mark the interval as locked, on the master vector  */
					RSB_DO_LOCK_INTERVAL(&(lock->olock),th_id,m0,m1-m0);/* FIXME */

					lock->it[th_id]=i;/* FIXME */
/*					RSB_BITVECTOR_SET(lock->ir,RSB_CONST_MAX_SUPPORTED_CORES,th_id); */
					
					/* let's give the vector info */
					*ov = RSB_MV_OFFSET(lock,nvi,0);
					*roff=m0;
					*m=m1-m0;
					goto done;
				}
			}
#if 0
			else
			if(RSB_LOUD_MVR_TESTING)
					RSB_INFO("interval %d in vector %d is not available\n",i,nvi);
#endif
		}
	}
	if(RSB_LOUD_MVR_TESTING)
		RSB_INFO("there are no available taint vectors.\n");
	/* no tainted subvectors found, or no free ones found. */
	/* in this case, we could look for a common subvector to both some v_i and v_j, i<j, and reduce it into v_i */
	goto done;
err:
	errval = RSB_ERR_BADARGS;
done:
	return errval;
}
#endif /* RSB_WANT_SPMV_WITH_REDUCE */

#if RSB_WANT_DO_LOCK_TEST
rsb_err_t rsb__do_perform_partial_reduce(struct rsb_mv_lock_t *lock, rsb_thr_t th_id, rsb_trans_t transA, rsb_coo_idx_t incv)
{
	/* FIXME: this is only an example routine */
	rsb_char_t * ov=NULL;
	rsb_coo_idx_t roff, m;
	
	rsb__do_pick_candidate_interval_for_reduce(lock,th_id,&ov,&roff,&m);

	if(!ov)
		goto done;
	RSB_ASSERT(lock->ov);
	if(RSB_LOUD_BTILS_TESTING)
		RSB_INFO("on thread %d about to reduce %d .. %d\n",th_id,roff,m+roff);
	rsb__vectors_left_sum_reduce_and_zero(lock->ov,ov,lock->typecode,m,incv,roff);
	/* perform reduce here */
	rsb__do_release_candidate_interval_for_reduce(lock,th_id,ov,roff,m);
	/*
	(with no symmetry or transposition, here)
	lock it for both vectors
	reduce the corresponding subvector (via sum), and zero it on the v_i vector
	update v_i's taint vector accordingly
	release the lock
	*/
done:
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_WANT_DO_LOCK_TEST */

#if RSB_WANT_DO_LOCK_TEST
rsb_err_t rsb__do_lock_test()
{
	/** 
	 * 	\ingroup gr_internals
	 * 	FIXME: NEW, UNFINISHED
	 **/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_bti_lock_struct lock;
       	rsb_coo_idx_t itl=8,nlevels = RSB_UPPER_BOUNDING_LOG2(itl);
	rsb_int in,it,i;
       	rsb_trans_t transA = RSB_TRANSPOSITION_N;

	RSB_INFO("LOCK TEST: BEGIN\n");

	RSB_ASSERT(nlevels==3);
	if((errval = rsb_do_btils_init(&lock,itl,nlevels))!=RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,2,4,transA,&in,&it));
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,4,8,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,0,4,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,0,8,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,2,3,transA,&in,&it));
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,0,1,transA,&in,&it));
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,1,2,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,0,2,transA,&in,&it));
	rsb_do_btils_lock_release(&lock,1,2);
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,0,2,transA,&in,&it));
	rsb_do_btils_lock_release(&lock,0,1);
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,0,2,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,2,3,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,3,4,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,4,5,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,6,7,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,7,8,transA,&in,&it));
	rsb_do_btils_lock_release(&lock,4,8);
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,4,5,transA,&in,&it));
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,6,7,transA,&in,&it));
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,7,8,transA,&in,&it));
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,7,8,transA,&in,&it));
	rsb_do_btils_free(&lock);

	itl=10,nlevels = RSB_UPPER_BOUNDING_LOG2(itl);/* 4 */
	RSB_ASSERT(nlevels==4);
	if((errval = rsb_do_btils_init(&lock,itl,nlevels))!=RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,0,itl,transA,&in,&it));
	for(i=0;i<itl;++i)
		RSB_ASSERT(!rsb_do_btils_lock_get(&lock,i,i+1,transA,&in,&it));
	for(i=0;i<RSB_MIDDLE(itl);++i)
		RSB_ASSERT(!rsb_do_btils_lock_get(&lock,2*i,2*i+1,transA,&in,&it));
	rsb_do_btils_lock_release(&lock,0,itl);
	for(i=0;i<RSB_MIDDLE(itl);++i)
		RSB_ASSERT( rsb_do_btils_lock_get(&lock,2*i,2*i+1,transA,&in,&it)),
		RSB_ASSERT(!rsb_do_btils_lock_get(&lock,2*i,2*i+1,transA,&in,&it));
	for(i=0;i<RSB_MIDDLE(itl);++i)
		rsb_do_btils_lock_release(&lock,2*i,2*i+1);
	RSB_ASSERT( rsb_do_btils_lock_get(&lock,0,RSB_MIDDLE(itl),transA,&in,&it)),
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,0,RSB_MIDDLE(itl),transA,&in,&it)),
	RSB_ASSERT(!rsb_do_btils_lock_get(&lock,0,RSB_MIDDLE(itl),transA,&in,&it)),
	rsb_do_btils_free(&lock);

	/*
	 * TODO:
	 * need symmetry and transposition support.
	 * need routines for reducing the temporary vectors after 'failed' double loops
	 * */
	RSB_INFO("BINARY LOCK TEST OK\n");
	//RSB_INFO("binary tree lock test ok\n");
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
       	struct rsb_mv_lock_t lock;
	rsb_int_t num_threads=4,th_id=0;
	struct rsb_mtx_t *mtxAp=NULL;
	struct rsb_mtx_t *submatrix=NULL;
	enum rsb_op_flags_t op_flags = RSB_OP_FLAG_DEFAULT;
	rsb_submatrix_idx_t si=0;
	rsb_char_t * y=NULL,*oy=NULL,*oY=NULL,*oX=NULL;
       	rsb_submatrix_idx_t subms=0;
	rsb_coo_idx_t incv=1;
	mtxAp = rsb__generate_dense_lower_triangular(2000,NULL,RSB_NUMERICAL_TYPE_DEFAULT);

	if(!mtxAp)
	{ RSB_PERR_GOTO(erri,RSB_ERRM_ES); }
       	subms=mtxAp->all_leaf_matrices_n;
	y = rsb__calloc(mtxAp->el_size*mtxAp->nr*incv);
	if(!y)
	{ RSB_PERR_GOTO(erri,RSB_ERRM_ES); }
	if((errval = rsb__do_mv_lock_init(&lock,num_threads,subms,mtxAp,op_flags,transA,y,incv))!=RSB_ERR_NO_ERROR)
	{ RSB_PERR_GOTO(erri,RSB_ERRM_ES); }
	
	submatrix=mtxAp->all_leaf_matrices[si].mtxlp;
	RSB_ASSERT(rsb__do_mv_lock_get(&lock,th_id,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si,transA,&oy,&incv));
RSB_ASSERT(!rsb__do_mv_lock_get(&lock,th_id+1,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si+1,transA,&oY,&incv));
RSB_ASSERT( rsb__do_mv_lock_get(&lock,th_id+1,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si+1,transA,&oY,&incv));
RSB_ASSERT(!rsb__do_mv_lock_get(&lock,th_id,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si,transA,&oy,&incv));
RSB_ASSERT(!rsb__do_mv_lock_get(&lock,th_id+2,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si+2,transA,&oX,&incv));
RSB_ASSERT( rsb__do_mv_lock_get(&lock,th_id+2,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si+2,transA,&oX,&incv));
RSB_ASSERT(!rsb__do_mv_lock_get(&lock,th_id+3,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si+3,transA,&oy,&incv));
RSB_ASSERT(!rsb__do_mv_lock_get(&lock,th_id+3,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,si+3,transA,&oy,&incv));
	RSB_ASSERT(!rsb__do_mv_lock_release(&lock,th_id,oy));
	RSB_ASSERT(!rsb__do_mv_lock_release(&lock,th_id+2,oY));
	RSB_ASSERT(!rsb__do_mv_lock_release(&lock,th_id+1,oX));
	RSB_ASSERT(!rsb__do_mv_lock_release(&lock,th_id,oy));	/* harmless duplicate */
	RSB_ASSERT(!rsb__do_mv_lock_release(&lock,th_id+1,oX));	/* harmless duplicate */

	rsb__do_perform_partial_reduce(&lock,th_id,transA,incv);
	rsb__do_perform_partial_reduce(&lock,th_id+1,transA,incv);

	/* 
		The following idea was inspired by Frigo's 'reducers & hyperobjects' paper.
		To support it, we could extend the rsb_bool_t to handle a trivalent logic:
		RSB_BOOL_FALSE=0, RSB_BOOL_TRUE=1, RSB_BOOL_ALMOST=2
		When encountering RSB_BOOL_ALMOST, the thread would perform the reducing strategy.
		After a detected loop, the lock (which could effectively turned out into a scheduler)
	       	would "propose" the thread to perform a "partial reduce".

	 * */
#if RSB_INT_ERR_VERBOSITY==1
	/* TODO: this way of handling things forces 'incx' to be 1, then */
	RSB_INFO("FIXME: missing handling of after-reduce release.\n");
#endif /* RSB_INT_ERR_VERBOSITY */
	goto oki;
erri:
	//RSB_INFO("binary tree based multi-lock test problems..\n");
	RSB_INFO("LOCK TEST: END (PROBLEMS)\n");
oki:
	RSB_MTX_FREE(mtxAp);
	RSB_CONDITIONAL_FREE(y);
	rsb__do_mv_lock_free(&lock);
}
	goto ok;
ok:
	//RSB_INFO("binary tree based multi-lock test ok\n");
	RSB_INFO("BINARY LOCK TEST: END (SUCCESS)\n");
	return RSB_ERR_NO_ERROR;
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_DO_LOCK_TEST */

/*  END EXPERIMENTAL CODE */
/* @endcond */
