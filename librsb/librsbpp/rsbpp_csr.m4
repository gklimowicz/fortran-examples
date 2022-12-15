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
dnl
include(`rsbpp.m4')dnl
dnl
dnl
#ifndef RSBPP_CSR_HPP
#define RSBPP_CSR_HPP

rsb_err_t rsbpp_csr_spmm(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const void* VA, const void* IP, const void* JA, rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const void * rhs, const rsb_coo_idx_t ldY, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, const rsb_flags_t order)
{
	const bool by_rows = (order == RSB_FLAG_WANT_ROW_MAJOR_ORDER);

	typecode = toupper(typecode);

	if(nrhs < 1)
		return RSB_ERR_BADARGS;

	if(nr < 1)
		return RSB_ERR_BADARGS;

	if( ! is_valid_trans(transA) )
		return RSB_ERR_BADARGS;

	if( ( incy  > 1 || incx  > 1 ) && nrhs == 1 )
		return RSB_ERR_UNSUPPORTED_OPERATION;

	if( ( incy != 1 || incx != 1 ) && nrhs == 1 )
		return RSB_ERR_BADARGS;

	if ( RSB_DO_FLAG_HAS(flags , RSB_FLAG_USE_HALFWORD_INDICES) )
	{
foreach(`NT',(`RSBPP_M4_IMPLEMENTED_TYPES'),`dnl
pushdef(`CXX_NT',RSB_M4_MATRIX_TYPE_AS_CXX(NT))`'dnl
`#if 'RSB_M4_CHOPTRAILINGSPACES(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(NT))
		if ( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(NT))
			return rsbpp_csr_spmx<CXX_NT,rsb_coo_idx_t,rsb_half_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const CXX_NT*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IP), reinterpret_cast<const rsb_half_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const CXX_NT*>(rhs), ldY, reinterpret_cast<CXX_NT*>(out), reinterpret_cast<const CXX_NT*>(alphap), incx, incy, transA, roff, coff, by_rows);
`#endif '/* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(NT)*/
popdef(`CXX_NT')`'dnl
')dnl
	}
	else
	{
foreach(`NT',(`RSBPP_M4_IMPLEMENTED_TYPES'),`dnl
pushdef(`CXX_NT',RSB_M4_MATRIX_TYPE_AS_CXX(NT))`'dnl
`#if 'RSB_M4_CHOPTRAILINGSPACES(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(NT))
		if ( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(NT))
			return rsbpp_csr_spmx<CXX_NT,rsb_coo_idx_t>(flags,nnz,nr,nc, reinterpret_cast<const CXX_NT*>(VA), reinterpret_cast<const rsb_coo_idx_t*>(IP), reinterpret_cast<const rsb_coo_idx_t*>(JA), nrhs, ldX, reinterpret_cast<const CXX_NT*>(rhs), ldY, reinterpret_cast<CXX_NT*>(out), reinterpret_cast<const CXX_NT*>(alphap), incx, incy, transA, roff, coff, by_rows);
`#endif '/* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(NT)*/
popdef(`CXX_NT')`'dnl
')dnl
	}
	return RSB_ERR_UNSUPPORTED_TYPE;
}

rsb_err_t rsbpp_csr_spmv(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const void* VA, const void* IP, const void* JA, const void * rhs, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff)
{
	const rsb_coo_idx_t nrhs=1;
	const rsb_coo_idx_t ldX=0, ldY=0;
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;

	return rsbpp_csr_spmm(typecode, flags, nnz, nr, nc, VA, IP, JA, nrhs, ldX, rhs, ldY, out, alphap, incx, incy, transA, roff, coff, order);
}

#endif /* RSBPP_CSR_HPP */
