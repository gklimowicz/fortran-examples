/*
Copyright (C) 2020-2022 Michele Martone

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
#include "rsbpp.hpp"
#include "rsbck.hpp"
#include "librsbpp.h"
#include <cassert>

#if defined(RSB_NUMERICAL_TYPE_DOUBLE) && defined(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX) && defined(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX) && defined(RSB_NUMERICAL_TYPE_FLOAT) && RSBPP_WANT_ALL

static
void test_Coo_BuildsEmptyFromVectors(void)
{
	std::vector<int> IA{},JA{};
	std::vector<double> VA{};
#if USE_CXX17
	Coo m{VA,IA,JA};
#else
	Coo<int,double> m{VA,IA,JA};
#endif

	assert(1 == m.nr());
	assert(1 == m.nc());
	assert(0 == m.nnz());
}

static
void test_Coo_BuildsNonEmptyFromVectors(void)
{
	std::vector<int> IA{0},JA{0};
	std::vector<double> VA{+11};
#if USE_CXX17
	Coo m{VA,IA,JA};
#else
	Coo<int,double> m{VA,IA,JA};
#endif

	assert(1 == m.nr());
	assert(1 == m.nc());
	assert(1 == m.nnz());
}

static
void test_Csr_BuildsNonEmptyFromVectors(void)
{
	std::vector<int> IA{0},JA{0};
	std::vector<double> VA{+11};
	Csr<int,double> m{VA,IA,JA};

	assert(1 == m.nr());
	assert(1 == m.nc());
	assert(1 == m.nnz());
}

static
void test_Coo_IsEmptyInstantiable(void)
{
	Coo<int,double> m{};

	assert(0 == m.nr());
	assert(0 == m.nc());
	assert(0 == m.nnz());
}

static
void test_Coo_SpMV_Symmetric(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{+11};
	const std::vector<NT> X{21};
	std::vector<NT> Y{+11};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=1;
	const IT ncA=1;
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t roff{0};
	const rsb_coo_idx_t coff{0};

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( +11 + 2 * +11 * 21 == Y[0]);

	auto ret2 = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert(     0 == ret2 );
	assert(+473 + 2 * +11 * 21 == Y[0]);

	auto ret3 = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );
	assert(  0 == ret3 );

	auto ret4 = rsbpp_coo_spmv('X',flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );
	assert( RSB_ERR_UNSUPPORTED_TYPE == ret4 );
}

static
void test_Coo_SpMV_Complex_Unsymmetric_Conjugated_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,-22,-23,+24};
	std::vector<NT> Y{-11,+13};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA{4};
	const IT ncA{2};
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_C;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11           ) == Y[0]);
	assert(NT(+13) + NT(0,+2 *-+44 * +24) == Y[1]);
}

static
void test_Csr_SpMV_Complex_Unsymmetric_Conjugated_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IP{0,1},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,-22,-23,+24};
	std::vector<NT> Y{-11,+13};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA{1};
	const IT ncA{1};
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_C;

	auto ret = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11           ) == Y[0]);
	assert(NT(+13) + NT(0,+2 *-+44 * +24) == Y[1]);
}

static
void test_Coo_SpMV_Complex_Symmetric_Transposed_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA{4};
	const IT ncA{2};
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_T;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 * +44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 * +44 * +22) == Y[3]);
}

static
void test_Coo_SpMV_Complex_Symmetric_Conjugated_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA{1};
	const IT ncA{1};
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_C;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 * -44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 * -44 * +22) == Y[3]);
}

static
void test_Csr_SpMV_Complex_Symmetric_Transposed_At(void)
{
	using IT = rsb_coo_idx_t;
	using CT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<CT> JA{0};
	const std::vector<IT> IP{0,1};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA=1;
	const IT ncA=1;
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_T;

	auto ret = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 * +44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 * +44 * +22) == Y[3]);
}

static
void test_Csr_SpMV_Complex_Symmetric_Conjugated_At(void)
{
	using IT = rsb_coo_idx_t;
	using CT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<CT> JA{0};
	const std::vector<IT> IP{0,1};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA=1;
	const IT ncA=1;
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_C;

	auto ret = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 * -44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 * -44 * +22) == Y[3]);
}

static
void test_Coo_SpMV_Symmetric_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{+11};
	const std::vector<NT> X{-21,+22};
	std::vector<NT> Y{-11,3};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	const IT roff{1};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(-11            == Y[0]);
	assert( +3 + 2 * +11 * +22 == Y[1]);

	auto ret2 = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret2 );
	assert(-11             == Y[0]);
	assert(+487 + 2 * +11 * +22 == Y[1]);
}

static
void test_Coo_SpMV_Complex_Hermitian_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA=4;
	const IT ncA=4;
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_HERMITIAN};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 *-+44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 * +44 * +22) == Y[3]);
}

static
void test_Csr_SpMV_Complex_Hermitian_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IP{0,1},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA=1;
	const IT ncA=1;
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_HERMITIAN};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	auto ret = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 *-+44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 * +44 * +22) == Y[3]);
}

static
void test_Coo_SpMV_Complex_Hermitian_Transposed_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA=4;
	const IT ncA=4;
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_HERMITIAN};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_T;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 * +44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 *-+44 * +22) == Y[3]);
}

static
void test_Csr_SpMV_Complex_Hermitian_Transposed_At(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IP{0,1},JA{0};
	const std::vector<NT> VA{{0,+44}};
	const std::vector<NT> X{-21,+22,-23,+24};
	std::vector<NT> Y{-11,+12,-13,+14};
	const NT alpha {+2,0};
	const IT nnzA=VA.size();
	const IT nrA=4;
	const IT ncA=4;
	const IT roff{3};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_HERMITIAN};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_T;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(NT(-11)                 == Y[0]);
	assert(NT(+12) + NT(0,+2 * +44 * +24) == Y[1]);
	assert(NT(-13)                 == Y[2]);
	assert(NT(+14) + NT(0,+2 *-+44 * +22) == Y[3]);
}

static
void test_Coo_SpMV_Symmetric_At_Short(void)
{
	using IT = rsb_coo_idx_t;
	using CT = rsb_half_idx_t;
	using NT = double;
	const std::vector<CT> IA{0},JA{0};
	const std::vector<NT> VA{+44};
	const std::vector<NT> X{-11,+22};
	std::vector<NT> Y{-11,+12};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	const IT roff{1};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC|RSB_FLAG_USE_HALFWORD_INDICES};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(-11            == Y[0]);
	assert(+12 + 2 * +44 * +22 == Y[1]);

	auto ret2 = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret2 );
	assert(-11            == Y[0]);
	assert(+1948 + 2 * +44 * +22 == Y[1]);
}

static
void test_Csr_SpMV_Symmetric_At_Short(void)
{
	using IT = rsb_coo_idx_t;
	using CT = rsb_half_idx_t;
	using NT = double;
	const std::vector<IT> IP{0,1,1};
	const std::vector<CT> JA{0};
	const std::vector<NT> VA{+44};
	const std::vector<NT> X{-11,+22};
	std::vector<NT> Y{-11,+12};
	const NT alpha {+2};
	const IT nnzA=2;
	const IT nrA=1;
	const IT ncA=1;
	const IT roff{1};
	const IT coff{1};
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC|RSB_FLAG_USE_HALFWORD_INDICES};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	auto ret = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret );
	assert(-11            == Y[0]);
	assert(+12 + 2 * +44 * +22 == Y[1]);

	auto ret2 = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff);

	assert( RSB_ERR_NO_ERROR == ret2 );
	assert(-11            == Y[0]);
	assert(+1948 + 2 * +44 * +22 == Y[1]);
}

static
void test_Coo_SpMV_Complex(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IA{0},JA{0};
	const std::vector<NT> VA{+44};
	const std::vector<NT> X{+21};
	std::vector<NT> Y{+31};
	const NT alpha {+2,+1};
	const IT nnzA=VA.size();
	const IT nrA=1;
	const IT ncA=1;
	rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t roff{0}, coff{0};

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( NT{+31}  + (NT{+2,+1}) * NT{+44} * NT{+21} == Y[0]);
}

static
void test_Coo_SpMV_Complex_Hermitian(void)
{
	using IT = rsb_coo_idx_t;
	using NT = std::complex<double>;
	const std::vector<IT> IA{0},JA{1};
	const std::vector<NT> VA{{+4,+1}};
	const std::vector<NT> X{+1,+2};
	std::vector<NT> Y{-1,-2};
	const NT alpha {+3,+0};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_HERMITIAN};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t roff{0}, coff{0};

	auto ret1 = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert( RSB_ERR_NO_ERROR == ret1 );
	assert( NT{-1}  + (NT{3} * NT{+4,+1}) * NT{+2} == Y[0]);
	assert( NT{-2}  + (NT{3} * NT{+4,-1}) * NT{+1} == Y[1]);
}

static
void test_Coo_SpMV_No_Inc(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IA{0,0},JA{0,1};
	const std::vector<NT> VA{11,12};
	const std::vector<NT> X{+1,-1};
	std::vector<NT> Y{+3,-3};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t roff{0}, coff{0};

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,2*incx,incy,transA,roff,coff );
	assert( RSB_ERR_UNSUPPORTED_OPERATION == ret );
}

static
void test_Coo_SpMV_Unsymmetric(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IA{0,0},JA{0,1};
	const std::vector<NT> VA{11,12};
	const std::vector<NT> X{+1,-1};
	std::vector<NT> Y{+3,-3};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t roff{0}, coff{0};

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( 2 == Y.size() );

	assert(+3 + 2 * 11 * +1 + 2 * 12 * -1 == Y[0]);
	assert(-3 + 2 *  0 * +1 + 2 *  0 * -1 == Y[1]);
}

static
void test_Csr_SpMV_Unsymmetric(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IP{0,2,2},JA{0,1};
	const std::vector<NT> VA{+11,+12};
	const std::vector<NT> X{+21,+22};
	std::vector<NT> Y{+11,-12};
	const NT alpha {+2};
	const IT nnzA=2;
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t roff{0}, coff{0};

	auto ret = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( 2 == Y.size() );

	assert(+11 + 2 * +11 * +21 + 2 * +12 * +22 == Y[0]);
	assert(-12 + 2 *   0  * +21 + 2 *  0 * +22 == Y[1]);
}

static
void test_Csr_SpMM_Unsymmetric(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IP{0,2,2},JA{0,1};
	const std::vector<NT> VA{+11,+12};
	const std::vector<NT> X{+21,+22,+23,+24};
	std::vector<NT> Y{+11,-12,+13,-14};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_coo_idx_t nrhs = 2;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t ldX{2};
	const rsb_coo_idx_t ldY{2};
	const rsb_coo_idx_t roff{0};
	const rsb_coo_idx_t coff{0};
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;

	auto ret = rsbpp_csr_spmm(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),nrhs,ldX,X.data(),ldY,Y.data(),&alpha,incx,incy,transA,roff,coff,order );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( 4 == Y.size() );

	assert(+11 + 2 * +11 * +21 + 2 * +12 * +22 == Y[0]);
	assert(-12 + 2 *   0 * +21 + 2 *   0 * -22 == Y[1]);
	assert(+13 + 2 * +11 * +23 + 2 * +12 * +24 == Y[2]);
	assert(-14 + 2 *   0 * +23 + 2 *   0 * -24 == Y[3]);
}

static
void test_Csr_SpMM_br_Unsymmetric(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IP{0,2,2},JA{0,1};
	const std::vector<NT> VA{+11,+12};
	//const std::vector<NT> X{+21,+22,+23,+24};
	//std::vector<NT> Y{+11,-12,+13,-14};
	const std::vector<NT> X{+21,+23,+22,+24};
	std::vector<NT> Y{+11,+13,-12,-14};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t nrhs = 2;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t ldX{2};
	const rsb_coo_idx_t ldY{2};

	using CT = rsb_coo_idx_t;
	CT roff {0};
	CT coff {0};

	const rsbpp_CsrP<NT,IT,CT> csrp {nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),roff,coff};
	const IT bidx{0};
	const IT eidx{nrA};

	using RT = decltype(X);
	using OT = decltype(Y);
	constexpr int uf = 1;

	auto ret = spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,nrhs>(csrp, flags, transA, ldX, X, ldY, Y, 0, alpha, bidx, eidx);

	assert( RSB_ERR_NO_ERROR == ret );
	assert( 4 == Y.size() );

	assert(+11 + 2 * +11 * +21 + 2 * +12 * +22 == Y[0]);
	assert(-12 + 2 *   0 * +21 + 2 *   0 * -22 == Y[2]);
	assert(+13 + 2 * +11 * +23 + 2 * +12 * +24 == Y[1]);
	assert(-14 + 2 *   0 * +23 + 2 *   0 * -24 == Y[3]);
}

static
void test_Csr_SpMM_Unsymmetric_Short_At(void)
{
	using CT = short int;
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IP{0,2,2};
	const std::vector<CT> JA{0,1};
	const std::vector<NT> VA{+21,+22};
	const std::vector<NT> X{+21,+22,+23,+24};
	const NT alpha {+2};
	const IT nnzA=2;
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS|RSB_FLAG_USE_HALFWORD_INDICES};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_coo_idx_t nrhs = 2;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t coff{0};
	const rsb_coo_idx_t ldX{2};
	const rsb_coo_idx_t roff{1};
	const rsb_coo_idx_t ldY{2+roff};
	std::vector<NT> Y{0,-11,+12,0,-13,+14};
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;

	auto ret = rsbpp_csr_spmm(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),nrhs,ldX,X.data(),ldY,Y.data(),&alpha,incx,incy,transA,roff,coff,order );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( 6 == Y.size() );

	assert(  0                                         == Y[0+0]);
	assert(-11 + alpha * +21 * +21 + alpha * +22 * +22 == Y[1+0]);
	assert(+12                                         == Y[1+1]);
	assert(  0                                         == Y[1+2]);
	assert(-13 + alpha * +21 * +23 + alpha * +22 * +24 == Y[1+3]);
	assert( 14                                         == Y[1+4]);
}

static
void test_Coo_SpMM_Unsymmetric(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IA{0,0},JA{0,1};
	const std::vector<NT> VA{11,12};
	const std::vector<NT> X{+1,-1,+2,-2};
	std::vector<NT> Y{+3,-3};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const IT nrA{2};
	const IT ncA{2};
	const rsb_coo_idx_t ldX{ncA};
	const rsb_coo_idx_t ldY{nrA};
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
	const rsb_coo_idx_t roff{0};
	const rsb_coo_idx_t coff{0};

	auto ret1 = rsbpp_coo_spmm(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),   1,ldX,X.data(),ldY,Y.data(),&alpha,incx,incy,transA,roff,coff,order );

	assert(  0 == ret1 );
	assert( 2 == Y.size() );

	assert(+3 + 2 * 11 * +1 + 2 * 12 * -1 == Y[0]);
	assert(-3 + 2 *  0 * +1 + 2 *  0 * -1 == Y[1]);


	const IT nrhs{2};
	Y = {+3,-3,+3,-3};
	auto ret2 = rsbpp_coo_spmm(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),nrhs,ldX,X.data(),ldY,Y.data(),&alpha,incx,incy,transA,roff,coff,order );

	assert(  0 == ret2 );
	assert( 4 == Y.size() );

	assert(+3 + 2 * 11 * +1 + 2 * 12 * -1 == Y[0]);
	assert(-3 + 2 *  0 * +1 + 2 *  0 * -1 == Y[1]);

	assert(+3 + 2 * 11 * +2 + 2 * 12 * -2 == Y[2]);
	assert(-3 + 2 *  0 * +2 + 2 *  0 * -2 == Y[3]);
}

static
void test_Coo_SpMV_Unsymmetric_Transposed(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IA{0,1},JA{0,0};
	const std::vector<NT> VA{+11,+21};
	const std::vector<NT> X{+21,+22};
	std::vector<NT> Y{+31,-32};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_T;
	const rsb_coo_idx_t roff{0};
	const rsb_coo_idx_t coff{0};

	auto ret = rsbpp_coo_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IA.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( 2 == Y.size() );

	assert(+31 + 2 * +11 * +21 + 2 * +21 * +22 == Y[0]);
	assert(-32 + 2 *   0 * +21 + 2 *   0 * +22 == Y[1]);
}

void test_Csr_SpMV_Unsymmetric_Transposed(void)
{
	using IT = rsb_coo_idx_t;
	using NT = double;
	const std::vector<IT> IP{0,1,2},JA{0,0};
	const std::vector<NT> VA{+11,+21};
	const std::vector<NT> X{+21,+22};
	std::vector<NT> Y{+31,-32};
	const NT alpha {+2};
	const IT nnzA=VA.size();
	const IT nrA=2;
	const IT ncA=2;
	rsb_flags_t flags {RSB_FLAG_NOFLAGS};
	const rsb_coo_idx_t incx = 1, incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_T;
	const rsb_coo_idx_t roff{0}, coff{0};

	auto ret = rsbpp_csr_spmv(RSB_NUMERICAL_TYPE_DOUBLE,flags,nnzA,nrA,ncA,VA.data(),IP.data(),JA.data(),X.data(),Y.data(),&alpha,incx,incy,transA,roff,coff );

	assert( RSB_ERR_NO_ERROR == ret );
	assert( 2 == Y.size() );

	assert(+31 + 2 * +11 * +21 + 2 * +21 * +22 == Y[0]);
	assert(-32 + 2 *   0 * +21 + 2 *   0 * +22 == Y[1]);
}

int main()
{
	test_Coo_BuildsEmptyFromVectors();
	test_Coo_BuildsNonEmptyFromVectors();
	test_Coo_IsEmptyInstantiable();
	test_Coo_SpMV_Symmetric();
	test_Coo_SpMV_Unsymmetric();
	test_Coo_SpMM_Unsymmetric();
	test_Coo_SpMV_Unsymmetric_Transposed();
	test_Coo_SpMV_Complex();
	test_Coo_SpMV_No_Inc();
	test_Coo_SpMV_Symmetric_At();
	test_Coo_SpMV_Complex_Symmetric_Transposed_At();
	test_Coo_SpMV_Complex_Symmetric_Conjugated_At();
	test_Coo_SpMV_Complex_Hermitian_At();
	test_Coo_SpMV_Complex_Hermitian_Transposed_At();
	test_Coo_SpMV_Complex_Unsymmetric_Conjugated_At();
	test_Coo_SpMV_Symmetric_At_Short();
	test_Coo_SpMV_Complex_Hermitian();

	test_Csr_BuildsNonEmptyFromVectors();
	test_Csr_SpMV_Unsymmetric();
	test_Csr_SpMM_Unsymmetric();
	test_Csr_SpMM_br_Unsymmetric();
	test_Csr_SpMM_Unsymmetric_Short_At();
	test_Csr_SpMV_Symmetric_At_Short();
	test_Csr_SpMV_Complex_Symmetric_Transposed_At();
	test_Csr_SpMV_Complex_Symmetric_Conjugated_At();
	test_Csr_SpMV_Unsymmetric_Transposed();
	test_Csr_SpMV_Complex_Unsymmetric_Conjugated_At();
	test_Csr_SpMV_Complex_Hermitian_At();
	test_Csr_SpMV_Complex_Hermitian_Transposed_At();
}
#else
int main()
{
}
#endif
