dnl
dnl
dnl	@author: Michele Martone
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
/* @cond INNERDOC */
/**
 * @file
 * @brief Auxiliary functions.
 */
RSB_M4_HEADER_MESSAGE()dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_MERGE_H_INCLUDED
#define RSB_MERGE_H_INCLUDED
')
dnl
dnl
#include "rsb_common.h"
dnl
#ifndef RSB_PS_ASSERT
/*#define RSB_PS_ASSERT(e) assert(e)	*/ 	/* uncomment this to use   asserts */
#define RSB_PS_ASSERT(e)			/* uncomment this to avoid asserts */
#else /* RSB_PS_ASSERT */
#undef RSB_PS_ASSERT
#define RSB_PS_ASSERT(e) 
#endif /* RSB_PS_ASSERT */
dnl
dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
dnl

ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb__do_util_merge_sorted_subarrays_in_place_`'RSB_M4_TYPE_CODE(mtype)`'(
		mtype *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_char_t * W,
		rsb_nnz_idx_t annz, rsb_nnz_idx_t bnnz,
	       	size_t wsize, rsb_flags_t flags
		)
{
dnl
pushdef(`typecode',RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))`'dnl
dnl
	/**
	 * \ingroup gr_util
	 *       A         B
	 * +----------+----------+
	 *  <- annz -> <- bnnz ->
	 *
	 *        W
	 * +------------+
	 *  <-  wnnz  ->
	 *
	 * Merges an array containing two ordered subarrays A and B, 
	 * sized respectively with annz and bnnz elements, using a
	 * swap area sized wsize bytes. 
	 *
	 * NOTE: this is NOT an optimized code, just a naive one to have this functionality working.
	 */
	rsb_int_t wpasses;
	rsb_nnz_idx_t wnnz,nnz=annz+bnnz;
	mtype *VW=NULL;
       	rsb_coo_idx_t * IW=NULL;
       	rsb_coo_idx_t * JW=NULL;
	mtype *VB=NULL;
       	rsb_coo_idx_t * IB=NULL;
       	rsb_coo_idx_t * JB=NULL;
	size_t el_size=sizeof(mtype);
	int step;
	rsb_nnz_idx_t aoff=0,boff=0,woff=0;

	wnnz=wsize/(el_size+2*sizeof(rsb_coo_idx_t));
	VW=(mtype*)W;
	W+=el_size*wnnz;
	IW=(rsb_coo_idx_t*)W;
	W+=sizeof(rsb_coo_idx_t)*wnnz;
	JW=(rsb_coo_idx_t*)W;

	VB=VA+annz;
	IB=IA+annz;
	JB=JA+annz;

	wpasses=(annz+bnnz+wnnz-1)/wnnz;

#define RSB_COO_MOVE(VD,ID,JD,VS,IS,JS,doff,soff) \
		VD[doff]=VS[soff], \
		ID[doff]=IS[soff], \
		JD[doff]=JS[soff],++soff,++doff

#define RSB_CMP_COO_LESS_THAN(IA,JA,IB,JB,aoff,boff) \
		IA[aoff]<IB[boff] || ( IA[aoff]==IB[boff] && JA[aoff] < JB[boff] )

/*	RSB_STDOUT(" * \n");*/
/*	RSB_STDOUT("wsize=%d steps:%d wnnz=%d nnz=%d\n",wsize,wpasses,wnnz,nnz);*/

/*	RSB_STDOUT("SSSSSSsentinel:%x %d %d\n",IA+annz+bnnz,IA[annz+bnnz],JA[annz+bnnz]);*/
	for(step=0;step<wpasses;++step)
	{
		rsb_nnz_idx_t cnnz;
		if(step==wpasses-1)
			wnnz=nnz-step*wnnz;

		cnnz=wnnz;
		cnnz = RSB_MIN(cnnz,annz);
		cnnz = RSB_MIN(cnnz,bnnz);
/*		RSB_STDOUT("step:%d wnnz=%d annz=%d bnnz=%d cnnz=%d\n",step,wnnz,annz,bnnz,cnnz);*/
		/* merge wnnz elements from A and B in W */
		woff=0;
		aoff=boff=0;
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IA,JA,annz,typecode,NULL,flags));
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IB,JB,bnnz,typecode,NULL,flags));
	/*	RSB_STDOUT("SSSSsentinel:%x %d %d\n",IB+bnnz,IB[bnnz],JB[bnnz]); */
		while(woff<wnnz && aoff<annz && boff<bnnz)
		{
			if(RSB_CMP_COO_LESS_THAN(IA,JA,IB,JB,aoff,boff))
				RSB_COO_MOVE(VW,IW,JW,VA,IA,JA,woff,aoff);
			else
				RSB_COO_MOVE(VW,IW,JW,VB,IB,JB,woff,boff);
		}
/*		RSB_STDOUT("aoff=%d boff=%d woff=%d\n",aoff,boff,woff);*/
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IW,JW,woff,typecode,NULL,flags));
		if(woff<wnnz)
		{
			if(aoff==annz)
				RSB_COO_MEMMOVE(VW,IW,JW,VB,IB,JB,woff,boff,wnnz-woff,el_size),boff+=(wnnz-woff);
			else
			if(boff==bnnz)
				RSB_COO_MEMMOVE(VW,IW,JW,VA,IA,JA,woff,aoff,wnnz-woff,el_size),aoff+=(wnnz-woff);
		}
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IW,JW,wnnz,typecode,NULL,flags));
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IB,JB,bnnz,typecode,NULL,flags));
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IA,JA,annz,typecode,NULL,flags));
/*		RSB_STDOUT("aoff:%d boff=%d wnnz=%d annz=%d\n",aoff,boff,wnnz,annz);*/
		/* memmove A boff places forward */
		bnnz-=boff;
		annz-=aoff;
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IB,JB,bnnz,typecode,NULL,flags));
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IA,JA,annz,typecode,NULL,flags));
/*		RSB_STDOUT("SSSSsentinel:%x %d %d\n",IB+boff+bnnz,IB[boff+bnnz],JB[boff+bnnz]);*/
		RSB_COO_MEMMOVE(VA,IA,JA,VA,IA,JA,wnnz,aoff,annz,el_size);
/*		RSB_STDOUT("PSSSsentinel:%x %d %d\n",IB+boff+bnnz,IB[boff+bnnz],JB[boff+bnnz]);*/
		VB+=boff;
		IB+=boff;
		JB+=boff;
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IB,JB,bnnz,typecode,NULL,flags));
		RSB_COO_MEMMOVE(VA,IA,JA,VW,IW,JW,0,0,wnnz,el_size);
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IA,JA,wnnz,typecode,NULL,flags));
		RSB_PS_ASSERT(!rsb__util_is_sorted_coo_as_row_major(IA,JA,annz,typecode,NULL,flags));
		VA+=wnnz;
		IA+=wnnz;
		JA+=wnnz;
	}
	return RSB_ERR_NO_ERROR;
	#undef RSB_COO_MOVE
	#undef RSB_CMP_COO_LESS_THAN
dnl
popdef(`typecode')`'dnl
dnl
}
')dnl
')dnl
dnl
rsb_err_t rsb__do_util_merge_sorted_subarrays_in_place(
		void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_char_t * W,
		rsb_nnz_idx_t annz, rsb_nnz_idx_t bnnz,
	       	size_t wsize, rsb_flags_t flags, rsb_type_t typecode
		)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	switch(typecode)
	{
foreach(`mtype',RSB_M4_TYPES,`dnl
		case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype):
			return rsb__do_util_merge_sorted_subarrays_in_place_`'RSB_M4_TYPE_CODE(mtype)`'(VA,IA,JA,W,annz,bnnz,wsize,flags);
		break;
')dnl
		default :
			return RSB_ERR_UNSUPPORTED_TYPE;
	}
}
')dnl
dnl
dnl

dnl
rsb_err_t rsb__do_util_merge_sorted_subarrays_in_place_test(void)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;

foreach(`mtype',RSB_M4_TYPES,`dnl
{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const rsb_nnz_idx_t nnzA = 3, nnzB = 3;
	const rsb_nnz_idx_t nnzW = 3;
	rsb_coo_idx_t IA[] = {0,1,2,3,4,5};
	rsb_coo_idx_t JA[] = {0,1,2,3,4,5};
	mtype VA[] = {1,2,3,1,2,3};
	const rsb_coo_idx_t IR[] = {0,1,2,3,4,5};
	const rsb_coo_idx_t JR[] = {0,1,2,3,4,5};
	const mtype VR[] = {1,2,3,1,2,3};
	rsb_byte_t WA[2*nnzW*sizeof(IA[0])+nnzW*sizeof(VA[0])];
	const size_t wsize = sizeof(WA);
	rsb_flags_t flags = RSB_FLAG_WANT_ROW_MAJOR_ORDER;

	errval = rsb__do_util_merge_sorted_subarrays_in_place(VA, IA, JA, (void*)WA, nnzA, nnzB, wsize, flags, typecode);
 
	if(RSB_MEMCMP(IA,IR,(nnzA+nnzB)*sizeof(IA[0])))
		errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(JA,JR,(nnzA+nnzB)*sizeof(JA[0])))
		errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(VA,VR,(nnzA+nnzB)*sizeof(VA[0])))
		if(RSB_SOME_ERROR(rsb__do_are_same(VA, VR, nnzA+nnzB,typecode, 1, 1)))
			errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const rsb_nnz_idx_t nnzA = 3, nnzB = 3;
	const rsb_nnz_idx_t nnzW = 3;
	rsb_coo_idx_t IA[] = {0,1,2,3,4,5};
	rsb_coo_idx_t JA[] = {3,4,5,5,0,1};
	mtype VA[] = {1,2,3,1,2,3};
	const rsb_coo_idx_t IR[] = {0,1,2,3,4,5};
	const rsb_coo_idx_t JR[] = {3,4,5,5,0,1};
	const mtype VR[] = {1,2,3,1,2,3};
	rsb_byte_t WA[2*nnzW*sizeof(IA[0])+nnzW*sizeof(VA[0])];
	const size_t wsize = sizeof(WA);
	rsb_flags_t flags = RSB_FLAG_WANT_ROW_MAJOR_ORDER;

	errval = rsb__do_util_merge_sorted_subarrays_in_place(VA, IA, JA, (void*)WA, nnzA, nnzB, wsize, flags, typecode);
 
	if(RSB_MEMCMP(IA,IR,(nnzA+nnzB)*sizeof(IA[0])))
		errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(JA,JR,(nnzA+nnzB)*sizeof(JA[0])))
		errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(VA,VR,(nnzA+nnzB)*sizeof(VA[0])))
		if(RSB_SOME_ERROR(rsb__do_are_same(VA, VR, nnzA+nnzB,typecode, 1, 1)))
			errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}
')dnl
err:
	return errval;
}
')dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_MERGE_H_INCLUDED */
')
dnl


/* @endcond */
dnl
