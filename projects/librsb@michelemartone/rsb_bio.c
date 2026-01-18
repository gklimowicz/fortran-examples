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
/* @cond INNERDOC  */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains matrix binary I/O functions.
 * */
/*
 * RSB_HAVE_RPC_XDR_H -> RSB_WANT_XDR_SUPPORT
 * */

#include "rsb_common.h"
#if RSB_WANT_XDR_SUPPORT
//#ifndef RSB_HAVE_RPC_XDR_H
#if RSB_HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include <rpc/xdr.h>

#pragma GCC visibility push(hidden)

/*                                                                                                                   
                                                |  5  | >0 && <RSB_VSL|             50>RSB_VSL                           |  */
#define RSB_BINARY_SPARSE_MATRIX_FILE_SIGNATURE "%RSB-"
#ifndef RSB_PACKAGE_VERSION
#define RSB_PACKAGE_VERSION "?"
#endif
#define RSB_BINARY_SPARSE_MATRIX_FILE_HEADER RSB_BINARY_SPARSE_MATRIX_FILE_SIGNATURE""RSB_PACKAGE_VERSION"                                                  "
#define RSB_BINARY_SPARSE_MATRIX_FILE_SIGNATURE_LEN   5
#define RSB_BINARY_SPARSE_MATRIX_FILE_HEADER_LEN     32	/* the first RSB_BINARY_SPARSE_MATRIX_FILE_HEADER_LEN bytes of RSB_BINARY_SPARSE_MATRIX_FILE_HEADER are written at the beginning of the file */

#define RSB_XDR_ERROR {errval = RSB_ERR_INTERNAL_ERROR;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
#define RSB_XDR_LOAD_ERROR RSB_XDR_ERROR 
#define RSB_XDR_SAVE_ERROR RSB_XDR_ERROR 
#define RSB_XDR_SAVE_TRY(EXP) if(EXP!=1)RSB_XDR_ERROR 

#ifdef RSB_WANT_LONG_IDX_TYPE 
#define RSB_XDR_COO_IDX xdr_int64_t
#define RSB_XDR_NNZ_IDX xdr_int64_t
#else /* RSB_WANT_LONG_IDX_TYPE */
#define RSB_XDR_COO_IDX xdr_int32_t
#define RSB_XDR_NNZ_IDX xdr_int32_t
#endif /* RSB_WANT_LONG_IDX_TYPE */


rsb_err_t rsb__do_bindump_init(void)
{
	/*!
	 * \ingroup gr_bio
	 *
	 * \return RSB_ERR_NO_ERROR if the binary dumping of matrices is supported, or RSB_ERR_UNSUPPORTED_FEATURE.
	 */
#if RSB_WANT_XDR_SUPPORT
	rsb_err_t errval = RSB_ERR_UNSUPPORTED_FEATURE;
#else
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#endif /* RSB_WANT_XDR_SUPPORT */
	RSB_DO_ERR_RETURN(errval)
}



#if RSB_WANT_XDR_SUPPORT
static rsb_err_t rsb_do_rw_matrix_dimensions_xdr(struct rsb_mtx_t * mtxAp, XDR *xdrsp)
{
	/*!
	 * \ingroup gr_bio
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	uint64_t el_size_ = mtxAp->el_size;
	uint64_t element_count_ = mtxAp->element_count;

	RSB_XDR_SAVE_TRY(RSB_XDR_NNZ_IDX(xdrsp,&(mtxAp->nnz)));
	RSB_XDR_SAVE_TRY(xdr_uint64_t(xdrsp,&(el_size_)));
	RSB_XDR_SAVE_TRY(xdr_uint64_t(xdrsp,&(element_count_)));
#if RSB_WANT_DBC
	RSB_XDR_SAVE_TRY(RSB_XDR_NNZ_IDX(xdrsp,&(mtxAp->block_count)));
#endif
	RSB_XDR_SAVE_TRY(xdr_int32_t(xdrsp,&(mtxAp->all_leaf_matrices_n)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->nr)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->nc)));
	RSB_XDR_SAVE_TRY(xdr_char(xdrsp,&(mtxAp->typecode)));
	RSB_XDR_SAVE_TRY(xdr_int32_t(xdrsp,&(mtxAp->flags)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->roff)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->coff)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->bm)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->bk)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->broff)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->bcoff)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->roff)));
	RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,&(mtxAp->coff)));
	RSB_XDR_SAVE_TRY(RSB_XDR_NNZ_IDX(xdrsp,&(mtxAp->nzoff)));
	mtxAp->element_count = element_count_;
	mtxAp->el_size = el_size_;
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_rw_matrix_times_xdr(struct rsb_mtx_t * mtxAp, XDR *xdrsp)
{
	/*!
	 * \ingroup gr_bio
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_XDR_SAVE_TRY(xdr_double(xdrsp,&(mtxAp->sat)));
	RSB_XDR_SAVE_TRY(xdr_double(xdrsp,&(mtxAp->eit)));
	RSB_XDR_SAVE_TRY(xdr_double(xdrsp,&(mtxAp->est)));
	RSB_XDR_SAVE_TRY(xdr_double(xdrsp,&(mtxAp->pet)));
	RSB_XDR_SAVE_TRY(xdr_double(xdrsp,&(mtxAp->cpt)));
	RSB_XDR_SAVE_TRY(xdr_double(xdrsp,&(mtxAp->rpt)));
	RSB_XDR_SAVE_TRY(xdr_double(xdrsp,&(mtxAp->tat)));
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif

#if RSB_WANT_XDR_SUPPORT
static rsb_err_t rsb_do_rw_matrix_struct_xdr(struct rsb_mtx_t * mtxAp, struct rsb_mtx_t ** smp, XDR *xdrsp, const rsb_char_t rw)
{
	/*!
	 * \ingroup gr_bio
	 * FIXME: this code will break much of library's configurability
	 * to fix this, should use an intermediate struct before saving.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	errval = rsb_do_rw_matrix_dimensions_xdr(mtxAp,xdrsp);
	if(RSB_SOME_ERROR(errval))
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	errval = rsb_do_rw_matrix_times_xdr(mtxAp,xdrsp);
	if(RSB_SOME_ERROR(errval))
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(rw=='w')
	{
		uint32_t submatrices = 0;
		submatrices = (mtxAp->sm[0]!=NULL)*1+ (mtxAp->sm[1]!=NULL)*2+ (mtxAp->sm[2]!=NULL)*4+ (mtxAp->sm[3]!=NULL)*8;
		RSB_XDR_SAVE_TRY(xdr_uint32_t(xdrsp,&(submatrices)));
	}
	else
	if(rw=='r')
	{
		uint32_t submatrices = 0;
		int i;
		struct rsb_mtx_t * sm = *smp;
		RSB_XDR_SAVE_TRY(xdr_uint32_t(xdrsp,&(submatrices)));
		for(i=0;i<RSB_FOUR;++i)
		if(submatrices&(1<<i))
		{
			mtxAp->sm[i] = sm++;
		}
		else
			mtxAp->sm[i] = NULL;
		*smp = sm;
		RSB_DO_FLAG_ADD(mtxAp->flags,RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS);
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_rw_matrix_xdr_ia(struct rsb_mtx_t * mtxAp, struct rsb_mtx_t ** smp, rsb_nnz_idx_t *rnnzp, XDR *xdrsp, const rsb_char_t rw)
{
	/*!
	 * \ingroup gr_bio
	 * FIXME: this code will break much of library's configurability
	 * to fix this, should use an intermediate struct before saving.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if(rw!='r' &&  rw!= 'w') {	errval = RSB_ERR_INTERNAL_ERROR; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
	if(rw=='r')
	{
	}
	else
	{
		/* FIXME: write me */
	}
	errval = rsb_do_rw_matrix_struct_xdr(mtxAp,smp,xdrsp,rw);
#if RSB_WANT_DBC
	errval = rsb__set_init_flags_and_stuff(mtxAp,NULL,NULL,mtxAp->nr,mtxAp->nc,mtxAp->nnz,mtxAp->block_count,mtxAp->element_count,mtxAp->typecode,mtxAp->flags);
#else
	errval = rsb__set_init_flags_and_stuff(mtxAp,NULL,NULL,mtxAp->nr,mtxAp->nc,mtxAp->nnz,                 0,mtxAp->element_count,mtxAp->typecode,mtxAp->flags);
#endif
	if(RSB_SOME_ERROR(errval))
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i = 0,j = 0;
		struct rsb_mtx_t * submatrix = NULL;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_rw_matrix_xdr_ia(submatrix,smp,rnnzp,xdrsp,rw));
	}
	else
	{
		rsb_nnz_idx_t n = 0;
		rsb_nnz_idx_t * bpntr = NULL;
		if(smp && rw=='r')
			bpntr = (rsb_nnz_idx_t*)(*smp),
			mtxAp->bpntr = bpntr;

		if(mtxAp->bpntr)
		for(n=0;n<mtxAp->Mdim+1;++n)
		RSB_XDR_SAVE_TRY(RSB_XDR_NNZ_IDX(xdrsp,mtxAp->bpntr+n));

		*rnnzp += mtxAp->nnz;
		if(smp && rw=='r')
			*smp = (struct rsb_mtx_t*)(bpntr+mtxAp->Mdim+1);
	/*  
		Now it remains:
		void * VA;
		rsb_nnz_idx_t  *indptr;
		rsb_coo_idx_t	*bindx;
		rsb_coo_idx_t	*rpntr;
		rsb_coo_idx_t *cpntr;
		rsb_coo_idx_t *mpntr,*Mpntr;
		rsb_nnz_idx_t *bpntr;
		struct rsb_options_t *options;
	
		struct rsb_mtx_t * sm[RSB_FOUR];
		struct rsb_expected_info_t einfo;
		struct rsb_translated_matrix_t * all_leaf_matrices;
	 */
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_rw_matrix_xdr_ja(struct rsb_mtx_t * mtxAp, rsb_coo_idx_t * JA, rsb_nnz_idx_t *rnnzp, XDR *xdrsp, const rsb_char_t rw)
{
	/*!
	 * \ingroup gr_bio
	 * FIXME: this code will break much of library's configurability
	 * to fix this, should use an intermediate struct before saving.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if(rw!='r' &&  rw!= 'w') {	errval = RSB_ERR_INTERNAL_ERROR; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
	if(RSB_SOME_ERROR(errval))
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_rw_matrix_xdr_ja(submatrix,JA,rnnzp,xdrsp,rw));
	}
	else
	{
		rsb_nnz_idx_t n;
		if(rw=='r')
			mtxAp->bindx = JA+*rnnzp;
		/* should dump the rnnz VA and JA elements */
		for(n=0;n<mtxAp->nnz;++n)
			RSB_XDR_SAVE_TRY(RSB_XDR_COO_IDX(xdrsp,mtxAp->bindx+n));
		*rnnzp += mtxAp->nnz;
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_rw_matrix_xdr_va(struct rsb_mtx_t * mtxAp, rsb_char_t * VA, rsb_nnz_idx_t *rnnzp, XDR *xdrsp, const rsb_char_t rw)
{
	/*!
	 * \ingroup gr_bio
	 * FIXME: this code will break much of library's configurability
	 * to fix this, should use an intermediate struct before saving.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if(rw!='r' &&  rw!= 'w') {	errval = RSB_ERR_INTERNAL_ERROR; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
	if(rw=='r')
		mtxAp->VA = VA+mtxAp->el_size**rnnzp;
	if(RSB_SOME_ERROR(errval))
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_rw_matrix_xdr_va(submatrix,VA,rnnzp,xdrsp,rw));
	}
	else
	{
		rsb_nnz_idx_t n;
		if(0)
			;
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
		else
		if(mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE)
		{
			for(n=0;n<mtxAp->nnz;++n)
				RSB_XDR_SAVE_TRY(xdr_double(xdrsp,((double*)(mtxAp->VA))+n));
		}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT
		else
		if(mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT)
		{
			for(n=0;n<mtxAp->nnz;++n)
				RSB_XDR_SAVE_TRY(xdr_float(xdrsp,((float*)(mtxAp->VA))+n));
		}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
		else
		if(mtxAp->typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX)
		{
			for(n=0;n<2*mtxAp->nnz;++n)
				RSB_XDR_SAVE_TRY(xdr_double(xdrsp,((double*)(mtxAp->VA))+n));
		}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
		else
		if(mtxAp->typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX)
		{
			for(n=0;n<2*mtxAp->nnz;++n)
				RSB_XDR_SAVE_TRY(xdr_float(xdrsp,((float*)(mtxAp->VA))+n));
		}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
		else
		{
			/* TODO: if you have a new type, complete here */
			errval = RSB_ERR_UNSUPPORTED_TYPE;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
		*rnnzp += mtxAp->nnz;
	}
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_XDR_SUPPORT */

static rsb_err_t rsb_do_compute_total_bytes_for_binary_dump_recursive(const struct rsb_mtx_t * mtxAp, uint64_t * ia_size, uint64_t * ja_size, uint64_t * va_size)
{
	/*!
	 * \ingroup gr_bio
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i,j;
		struct rsb_mtx_t * submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_compute_total_bytes_for_binary_dump_recursive(submatrix,ia_size,ja_size,va_size));
	}

	*ia_size += sizeof(struct rsb_mtx_t)+sizeof(rsb_nnz_idx_t)*(mtxAp->Mdim+1);
	*ja_size += 0;
	*va_size += 0;

	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_compute_total_bytes_for_binary_dump(const struct rsb_mtx_t * mtxAp, uint64_t * ia_size, uint64_t * ja_size, uint64_t * va_size)
{
	/*!
	 * \ingroup gr_bio
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	*ia_size = 0;
	*ja_size = 0;
	*va_size = 0;
	*ia_size += sizeof(struct rsb_translated_matrix_t)*(rsb__terminal_recursive_matrix_count(mtxAp));
	*ja_size += sizeof(rsb_coo_idx_t)*(mtxAp->nnz);// FIXME: nnz or nnz+1 (locally in each mtxAp, of course) ?
	*va_size += mtxAp->el_size*mtxAp->nnz;
	errval = rsb_do_compute_total_bytes_for_binary_dump_recursive(mtxAp,ia_size,ja_size,va_size);
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_load_matrix_file_as_binary(struct rsb_mtx_t ** mtxApp, const rsb_char_t * filename)
{
	/*!
	 * \ingroup gr_bio
	 */
#if RSB_WANT_XDR_SUPPORT
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	uint64_t ia_size = 0,ja_size = 0,va_size = 0;
	FILE *fd = NULL;
	void *IA = NULL, *JA = NULL, *VA = NULL;
	XDR xdrs;
	rsb_time_t lt;
	struct rsb_mtx_t * mtxAp = NULL, *smp = NULL;
	rsb_nnz_idx_t rnnz = 0;
	const char * signature[RSB_BINARY_SPARSE_MATRIX_FILE_HEADER_LEN];
	u_int slen = RSB_BINARY_SPARSE_MATRIX_FILE_HEADER_LEN;

	RSB_IO_NOTICE("binary loading file %s..\n",filename);
	lt = - rsb_time();

	fd = fopen(filename,"r");
	xdrstdio_create(&xdrs,fd,XDR_DECODE); 
	RSB_XDR_SAVE_TRY(fread(signature,slen,1,fd));
//	if(RSB_MEMCMP(signature,RSB_BINARY_SPARSE_MATRIX_FILE_HEADER,RSB_BINARY_SPARSE_MATRIX_FILE_HEADER_LEN))
	if(RSB_MEMCMP(signature,RSB_BINARY_SPARSE_MATRIX_FILE_SIGNATURE,RSB_BINARY_SPARSE_MATRIX_FILE_SIGNATURE_LEN))
	{
		RSB_IO_ERROR("wrong file signature (not beginning with %s): skipping..\n",RSB_BINARY_SPARSE_MATRIX_FILE_SIGNATURE);
		goto ierr;
	
	}
	RSB_XDR_SAVE_TRY(xdr_uint64_t(&xdrs,&(ia_size)));
	RSB_XDR_SAVE_TRY(xdr_uint64_t(&xdrs,&(ja_size)));
	RSB_XDR_SAVE_TRY(xdr_uint64_t(&xdrs,&(va_size)));
	/* FIXME: should validate input, here */
	IA = rsb__calloc(ia_size);
	JA = rsb__calloc(ja_size);
	VA = rsb__calloc(va_size);
	if(!IA || !JA || !VA)
		goto ierr;	/* FIXME: err should close streams */
	mtxAp = IA; // FIXME
	smp = mtxAp+1;

	errval = rsb_do_rw_matrix_xdr_ia(mtxAp,&smp,&rnnz,&xdrs,'r');
	if(RSB_SOME_ERROR(errval))
		goto ierr;
	if(rnnz!=mtxAp->nnz)
	{
		RSB_IO_ERROR("error : read %ld instead of %ld nnz!\n",(long int)rnnz,(long int)mtxAp->nnz);
		errval = RSB_ERR_GENERIC_ERROR; goto ierr;
	}


	rnnz = 0;
	errval = rsb_do_rw_matrix_xdr_ja(mtxAp,JA,&rnnz,&xdrs,'r');
	if(rnnz!=mtxAp->nnz)
	{
		RSB_IO_ERROR("error : read %ld instead of %ld nnz!\n",(long int)rnnz,(long int)mtxAp->nnz);
		errval = RSB_ERR_GENERIC_ERROR; goto ierr;
	}
	rnnz = 0;
	errval = rsb_do_rw_matrix_xdr_va(mtxAp,VA,&rnnz,&xdrs,'r');
	if(rnnz!=mtxAp->nnz)
	{
		RSB_IO_ERROR("error : read %ld instead of %ld nnz!\n",(long int)rnnz,(long int)mtxAp->nnz);
		errval = RSB_ERR_GENERIC_ERROR; goto ierr;
	}

	mtxAp->all_leaf_matrices = (struct rsb_translated_matrix_t*)smp;
	errval = rsb__get_array_of_leaf_matrices(mtxAp,&mtxAp->all_leaf_matrices,NULL);
	if(RSB_SOME_ERROR(errval))
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
//	rsb__do_get_first_submatrix(mtxAp)->bindx = JA;
//	rsb__do_get_first_submatrix(mtxAp)->VA = VA;		// FIXME: temporarily here

	// place a check here
	if(!rsb__mtx_chk(mtxAp))
	{
		errval = RSB_ERR_CORRUPT_INPUT_DATA;
		RSB_PERR_GOTO(err,RSB_ERRM_EM);
	}
ierr:
	xdr_destroy(&xdrs);
	if(fclose(fd)!=0)
	{
		// NOTE: we ignore this error
	}
	/* FIXME: the matrix should be validated, now, before returning */
	if(!mtxAp)
	{RSB_PERR_GOTO(err,RSB_ERRM_ES);}
	*mtxApp = mtxAp;
	RSB_DO_FLAG_ADD(mtxAp->flags,RSB_FLAG_FIX_FOR_BINARY_LOADED_MATRIX);

	lt += rsb_time();
	RSB_IO_NOTICE("#ia_size %lld..\n",(long long int)ia_size);
	RSB_IO_NOTICE("#ja_size %lld..\n",(long long int)ja_size);
	RSB_IO_NOTICE("#va_size %lld..\n",(long long int)va_size);
	RSB_IO_NOTICE("#binary loading file %s succeeded and took %lf s (%.0f nnz/s).\n",filename,lt,(1.0/(lt/mtxAp->nnz)));
	/* FIXME : this is debug info */
//	RSB_DO_ERROR_CUMULATE(errval,rsb__do_print_matrix_stats(mtxAp,RSB_CONST_DUMP_TIMES|RSB_CONST_DUMP_DIMENSIONS|RSB_CONST_DUMP_RECURSION,NULL));
	goto ret;
err:
	/* FIXME : missing error handling */
	RSB_MTX_FREE(mtxAp);
#else /* RSB_WANT_XDR_SUPPORT */
	rsb_err_t errval = RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_WANT_XDR_SUPPORT */
ret:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_save_matrix_file_as_binary(const struct rsb_mtx_t * mtxAp, FILE * fd)
{
	/*!
	 * \ingroup gr_bio
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_WANT_XDR_SUPPORT
	uint64_t ia_size = 0,ja_size = 0,va_size = 0;
	XDR xdrs;
	rsb_nnz_idx_t rnnz = 0;
	const char * signature = RSB_BINARY_SPARSE_MATRIX_FILE_HEADER;
	u_int slen = RSB_BINARY_SPARSE_MATRIX_FILE_HEADER_LEN;

	if(!mtxAp || !fd)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	// file signature dump
	if(!rsb__mtx_chk(mtxAp))
	{
		errval = RSB_ERR_CORRUPT_INPUT_DATA;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	errval = rsb_do_compute_total_bytes_for_binary_dump(mtxAp,&ia_size,&ja_size,&va_size);
	xdrstdio_create(&xdrs,fd,XDR_ENCODE); 
	RSB_XDR_SAVE_TRY(fwrite(signature,slen,1,fd));
	RSB_XDR_SAVE_TRY(xdr_uint64_t(&xdrs,&(ia_size)));
	RSB_XDR_SAVE_TRY(xdr_uint64_t(&xdrs,&(ja_size)));
	RSB_XDR_SAVE_TRY(xdr_uint64_t(&xdrs,&(va_size)));
	rnnz = 0;
	errval = rsb_do_rw_matrix_xdr_ia((struct rsb_mtx_t*)mtxAp,NULL,&rnnz,&xdrs,'w');
	if(rnnz!=mtxAp->nnz)
	{
		RSB_IO_ERROR("error : wrote %ld instead of %ld nnz!\n",(long int)rnnz,(long int)mtxAp->nnz);
		errval = RSB_ERR_GENERIC_ERROR; goto ierr;
	}
	rnnz = 0;
	errval = rsb_do_rw_matrix_xdr_ja((struct rsb_mtx_t*)mtxAp,NULL,&rnnz,&xdrs,'w');
	if(rnnz!=mtxAp->nnz)
	{
		RSB_IO_ERROR("error : wrote %ld instead of %ld nnz!\n",(long int)rnnz,(long int)mtxAp->nnz);
		errval = RSB_ERR_GENERIC_ERROR; goto ierr;
	}
	rnnz = 0;
	errval = rsb_do_rw_matrix_xdr_va((struct rsb_mtx_t*)mtxAp,NULL,&rnnz,&xdrs,'w');
	if(rnnz!=mtxAp->nnz)
	{
		RSB_IO_ERROR("error : read %ld instead of %ld nnz!\n",(long int)rnnz,(long int)mtxAp->nnz);
		errval = RSB_ERR_GENERIC_ERROR; goto ierr;
	}
ierr:
	xdr_destroy(&xdrs);
#else /* RSB_WANT_XDR_SUPPORT */
	rsb_err_t errval = RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_WANT_XDR_SUPPORT */
err:
	RSB_DO_ERR_RETURN(errval)
}

#else /* RSB_WANT_XDR_SUPPORT */
rsb_err_t rsb__do_bindump_init(void){return RSB_ERR_UNSUPPORTED_FEATURE;}
rsb_err_t rsb__do_load_matrix_file_as_binary(struct rsb_mtx_t ** mtxApp, const rsb_char_t * filename){return RSB_ERR_UNSUPPORTED_FEATURE;}
rsb_err_t rsb__do_save_matrix_file_as_binary(const struct rsb_mtx_t * mtxAp, FILE * fd){return RSB_ERR_UNSUPPORTED_FEATURE;}
#endif /* RSB_WANT_XDR_SUPPORT */
/* @endcond */
