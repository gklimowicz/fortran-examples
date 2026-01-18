/* @cond INNERDOC */
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block compressed sparse stripes format.
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
#include "rsb_internals.h"
#include "rsb.h"

#pragma GCC visibility push(hidden)
rsb_err_t rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*1)];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*1));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*1));
			const double aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * rhs, double * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dE_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[fk];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dE_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=VA[lk-1];
		if(aa == ((double)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double bb_0=rhs[(1*i*(incx))];
		double ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double *b=out + (1*(j*(incx)));
			double *c=&ax_0;
{	{

		register double c_0 = ((double)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double *c_0=out+(1*(i*(incy)));
			const double aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dI_uU(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dI_uL(const double * restrict VA, const double * restrict rhs, double * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double *a=VA;
	const double alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double ax_0;
		const double aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*1)];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*1));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*1));
			const float aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * rhs, float * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dE_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[fk];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dE_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=VA[lk-1];
		if(aa == ((float)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float bb_0=rhs[(1*i*(incx))];
		float ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float *b=out + (1*(j*(incx)));
			float *c=&ax_0;
{	{

		register float c_0 = ((float)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float *c_0=out+(1*(i*(incy)));
			const float aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dI_uU(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dI_uL(const float * restrict VA, const float * restrict rhs, float * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float *a=VA;
	const float alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float ax_0;
		const float aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*1)];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*1));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*1));
			const float complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * rhs, float complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_float_complex_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conjf(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dE_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[fk];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dE_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=VA[lk-1];
		if(aa == ((float complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const float complex bb_0=rhs[(1*i*(incx))];
		float complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const float complex *b=out + (1*(j*(incx)));
			float complex *c=&ax_0;
{	{

		register float complex c_0 = ((float complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			float complex *c_0=out+(1*(i*(incy)));
			const float complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dI_uU(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dI_uL(const float complex * restrict VA, const float complex * restrict rhs, float complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const float complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type float complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const float complex *a=VA;
	const float complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_float_complex_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		float complex ax_0;
		const float complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conjf(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*1)];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*1));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*1));
			const double complex aa=1;
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^T}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=*a*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * rhs, double complex * out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags)
{

	/**
	 * \ingroup rsb_doc_kernels
	 * Computes \f$y \leftarrow {A^H}^{-1} \cdot x, where A \neq A^T. \f$
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_uxua_double_complex_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*1)]/=aa;
		ax_0=out[1*(i*1)];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*1)]-=conj(*a)*ax_0;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dE_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+1;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a -= rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dE_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-1  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		if(lk-fk>0)
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
			a += rows*columns;
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dE_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dE_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[fk];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		a += rows*columns;
		for(k=fk+1,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dE_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal explicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dE_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=VA[lk-1];
		if(aa == ((double complex)(0)))return RSB_ERR_INVALID_NUMERICAL_DATA;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-1,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dI_uU\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tN_r1_c1_uu_sU_dI_uL\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		const double complex bb_0=rhs[(1*i*(incx))];
		double complex ax_0;
		ax_0=0;
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

			const double complex *b=out + (1*(j*(incx)));
			double complex *c=&ax_0;
{	{

		register double complex c_0 = ((double complex)(0));
				

		c_0 += a[(0*1)+0]*b[0];
			c[0]+= c_0 ;
	}	
}
		}
		{
			/* the last element (which for a lower triangular solve is on the diagonal)*/
			/* Lx=y ; x_0=y_0/L_1_1  */
			double complex *c_0=out+(1*(i*(incy)));
			const double complex aa=1;
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
		}
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tT_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=*a*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_coo_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_coo_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_C__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}


rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dI_uU(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0,j=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dI_uU\n");
	for(i=br;RSB_LIKELY(i<bc);++i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];
		
		for(k=fk+0,j=bindx[k];k<lk-0  ;++k,a += rows*columns,j=bindx[k])
		{

		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dI_uL(const double complex * restrict VA, const double complex * restrict rhs, double complex * restrict out, const rsb_coo_idx_t  Mdim,const rsb_coo_idx_t  mdim,const rsb_half_idx_t * restrict bindx,const rsb_nnz_idx_t * restrict bpntr,const rsb_nnz_idx_t *restrict indptr,const rsb_coo_idx_t * restrict rpntr,const rsb_coo_idx_t * restrict cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const double complex * restrict alphap,rsb_coo_idx_t incx, rsb_coo_idx_t incy)
{

	/**
	 * \ingroup rsb_doc_kernels
         * A blocked 1 x 1, stored in BCSR format, diagonal implicit, of type double complex, with rsb_half_idx_t column indices.
	 * \return \rsb_errval_inp_param_msg
	 */


	register rsb_coo_idx_t i=0;
	register rsb_nnz_idx_t k=0;
	const register rsb_coo_idx_t columns=1,rows=1;
	const double complex *a=VA;
	const double complex alpha=*alphap;	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in rsb__BCSR_spsv_sxsx_double_complex_H__tC_r1_c1_uu_sU_dI_uL\n");

	for(i=Mdim-1; RSB_LIKELY((i+1)>0);--i)
	{
		const rsb_nnz_idx_t fk=bpntr[i],lk=bpntr[i+1];
		double complex ax_0;
		const double complex aa=1;

		out[1*(i*(incy))]/=aa;
		ax_0=out[1*(i*(incx))];

		for(k=lk-1-0,a=VA+k;k+1>=fk+1+0;--k,a -= rows*columns)
		{
			const rsb_coo_idx_t j=bindx[k];
		out[1*(j*(incy))]-=conj(*a)*ax_0;
		}
		out[1*(i*(incy))]*=alpha;
	}

	return RSB_ERR_NO_ERROR;
}



/* @endcond */
