/*
Copyright (C) 2021-2022 Michele Martone

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
#ifndef RSBPP_COO_HPP
#define RSBPP_COO_HPP

rsb_err_t rsbpp_coo_spmm(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const void* VA, const void* IA, const void* JA, rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const void * rhs, const rsb_coo_idx_t ldY, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, const rsb_flags_t order)
{
	const bool by_rows = (order == RSB_FLAG_WANT_ROW_MAJOR_ORDER);

	typecode = toupper(typecode);

	if( ! is_valid_trans(transA) )
		return RSB_ERR_BADARGS;

	if(nrhs < 1)
		return RSB_ERR_BADARGS;

	if( ( incy  > 1 || incx  > 1 ) && nrhs == 1 )
		return RSB_ERR_UNSUPPORTED_OPERATION;

	if( ( incy != 1 || incx != 1 ) && nrhs == 1 )
		return RSB_ERR_BADARGS;

	if ( RSB_DO_FLAG_HAS(flags , RSB_FLAG_USE_HALFWORD_INDICES) )
	{
#if RSB_NUMERICAL_TYPE_INT 
		if ( typecode == RSB_NUMERICAL_TYPE_INT )
			return rsbpp_coo_spmx<int,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const int*>(VA), reinterpret_cast<const rsb_half_idx_t*>(IA), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const int*>(rhs), ldY, reinterpret_cast<int*>(out), reinterpret_cast<const int*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_INT */
#if RSB_NUMERICAL_TYPE_DOUBLE 
		if ( typecode == RSB_NUMERICAL_TYPE_DOUBLE )
			return rsbpp_coo_spmx<double,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const double*>(VA), reinterpret_cast<const rsb_half_idx_t*>(IA), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const double*>(rhs), ldY, reinterpret_cast<double*>(out), reinterpret_cast<const double*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
#if RSB_NUMERICAL_TYPE_FLOAT 
		if ( typecode == RSB_NUMERICAL_TYPE_FLOAT )
			return rsbpp_coo_spmx<float,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const float*>(VA), reinterpret_cast<const rsb_half_idx_t*>(IA), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const float*>(rhs), ldY, reinterpret_cast<float*>(out), reinterpret_cast<const float*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_FLOAT */
#if RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX 
		if ( typecode == RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX )
			return rsbpp_coo_spmx<std::complex<long double>,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const std::complex<long double>*>(VA), reinterpret_cast<const rsb_half_idx_t*>(IA), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const std::complex<long double>*>(rhs), ldY, reinterpret_cast<std::complex<long double>*>(out), reinterpret_cast<const std::complex<long double>*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX */
#if RSB_NUMERICAL_TYPE_LONG_DOUBLE 
		if ( typecode == RSB_NUMERICAL_TYPE_LONG_DOUBLE )
			return rsbpp_coo_spmx<long double,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const long double*>(VA), reinterpret_cast<const rsb_half_idx_t*>(IA), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const long double*>(rhs), ldY, reinterpret_cast<long double*>(out), reinterpret_cast<const long double*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_LONG_DOUBLE */
#if RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
		if ( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
			return rsbpp_coo_spmx<std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const std::complex<double>*>(VA), reinterpret_cast<const rsb_half_idx_t*>(IA), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const std::complex<double>*>(rhs), ldY, reinterpret_cast<std::complex<double>*>(out), reinterpret_cast<const std::complex<double>*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
#if RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
		if ( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
			return rsbpp_coo_spmx<std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const std::complex<float>*>(VA), reinterpret_cast<const rsb_half_idx_t*>(IA), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const std::complex<float>*>(rhs), ldY, reinterpret_cast<std::complex<float>*>(out), reinterpret_cast<const std::complex<float>*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
	}
	else
	{
#if RSB_NUMERICAL_TYPE_INT 
		if ( typecode == RSB_NUMERICAL_TYPE_INT )
			return rsbpp_coo_spmx<int,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const int*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IA), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const int*>(rhs), ldY, reinterpret_cast<int*>(out), reinterpret_cast<const int*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_INT */
#if RSB_NUMERICAL_TYPE_DOUBLE 
		if ( typecode == RSB_NUMERICAL_TYPE_DOUBLE )
			return rsbpp_coo_spmx<double,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const double*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IA), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const double*>(rhs), ldY, reinterpret_cast<double*>(out), reinterpret_cast<const double*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
#if RSB_NUMERICAL_TYPE_FLOAT 
		if ( typecode == RSB_NUMERICAL_TYPE_FLOAT )
			return rsbpp_coo_spmx<float,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const float*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IA), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const float*>(rhs), ldY, reinterpret_cast<float*>(out), reinterpret_cast<const float*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_FLOAT */
#if RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX 
		if ( typecode == RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX )
			return rsbpp_coo_spmx<std::complex<long double>,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const std::complex<long double>*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IA), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const std::complex<long double>*>(rhs), ldY, reinterpret_cast<std::complex<long double>*>(out), reinterpret_cast<const std::complex<long double>*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX */
#if RSB_NUMERICAL_TYPE_LONG_DOUBLE 
		if ( typecode == RSB_NUMERICAL_TYPE_LONG_DOUBLE )
			return rsbpp_coo_spmx<long double,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const long double*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IA), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const long double*>(rhs), ldY, reinterpret_cast<long double*>(out), reinterpret_cast<const long double*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_LONG_DOUBLE */
#if RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX 
		if ( typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX )
			return rsbpp_coo_spmx<std::complex<double>,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const std::complex<double>*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IA), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const std::complex<double>*>(rhs), ldY, reinterpret_cast<std::complex<double>*>(out), reinterpret_cast<const std::complex<double>*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
#if RSB_NUMERICAL_TYPE_FLOAT_COMPLEX 
		if ( typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX )
			return rsbpp_coo_spmx<std::complex<float>,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const std::complex<float>*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IA), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const std::complex<float>*>(rhs), ldY, reinterpret_cast<std::complex<float>*>(out), reinterpret_cast<const std::complex<float>*>(alphap), incx, incy, transA, roff, coff, by_rows);
#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
	}
	return RSB_ERR_UNSUPPORTED_TYPE;
}

rsb_err_t rsbpp_coo_spmv(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const void* VA, const void* IA, const void* JA, const void * rhs, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff)
{
	const rsb_coo_idx_t nrhs=1;
	const rsb_coo_idx_t ldX=0, ldY=0;
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;

	return rsbpp_coo_spmm(typecode, flags, nnz, nr, nc, VA, IA, JA, nrhs, ldX, rhs, ldY, out, alphap, incx, incy, transA, roff, coff, order);
}

#endif /* RSBPP_COO_HPP */
