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
#ifndef RSBCK_HPP_INCLUDED
#define RSBCK_HPP_INCLUDED

#include "rsbpp.hpp"

// <int,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<int,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const int* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const int* rhs, const rsb_coo_idx_t ldY, int* out, const int* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <double,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<double,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<float,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<long double>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<long double>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<long double>* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<long double>* rhs, const rsb_coo_idx_t ldY, std::complex<long double>* out, const std::complex<long double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <long double,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<long double,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const long double* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const long double* rhs, const rsb_coo_idx_t ldY, long double* out, const long double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <int,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<int,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const int* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const int* rhs, const rsb_coo_idx_t ldY, int* out, const int* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <double,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<double,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<float,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<long double>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<long double>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<long double>* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<long double>* rhs, const rsb_coo_idx_t ldY, std::complex<long double>* out, const std::complex<long double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <long double,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<long double,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const long double* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const long double* rhs, const rsb_coo_idx_t ldY, long double* out, const long double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);

// <int,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<int,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const int* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const int* rhs, const rsb_coo_idx_t ldY, int* out, const int* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <double,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<double,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<float,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<long double>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<long double>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<long double>* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<long double>* rhs, const rsb_coo_idx_t ldY, std::complex<long double>* out, const std::complex<long double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <long double,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<long double,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const long double* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const long double* rhs, const rsb_coo_idx_t ldY, long double* out, const long double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <int,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<int,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const int* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const int* rhs, const rsb_coo_idx_t ldY, int* out, const int* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <double,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<double,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<float,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<long double>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<long double>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<long double>* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<long double>* rhs, const rsb_coo_idx_t ldY, std::complex<long double>* out, const std::complex<long double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <long double,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<long double,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const long double* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const long double* rhs, const rsb_coo_idx_t ldY, long double* out, const long double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);

#endif /* RSBCK_HPP_INCLUDED */
