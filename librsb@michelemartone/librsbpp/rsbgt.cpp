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
#include "config.h"

#if (!RSB_HAVE_GTEST || !RSBPP_HAS_RSB_H)
int main()
{
}
#else

#include "rsb.h"
#include <librsbpp.h>
#include <vector>
#include <complex>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#if defined(__llvm__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments" /* for gtest under clang (yes, with '#pragma GCC') */
#endif

typedef unsigned short int rsbpp_half_idx_t;

namespace rsb_internal
{
const rsb_flags_t RSBPP_INVALID_FLAGS = -1;
template<typename NT> constexpr rsb_type_t rsb_type_t_for(void) { return RSB_NUMERICAL_TYPE_INVALID_TYPE; }
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
template<> constexpr rsb_type_t rsb_type_t_for<double> (void) { return RSB_NUMERICAL_TYPE_DOUBLE; }
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
template<> constexpr rsb_type_t rsb_type_t_for<long double> (void) { return RSB_NUMERICAL_TYPE_LONG_DOUBLE; }
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
template<> constexpr rsb_type_t rsb_type_t_for<std::complex<long double>> (void) { return RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX; }
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
template<> constexpr rsb_type_t rsb_type_t_for<std::complex<double>> (void) { return RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX; }
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
template<> constexpr rsb_type_t rsb_type_t_for<std::complex<float>> (void) { return RSB_NUMERICAL_TYPE_FLOAT_COMPLEX; }
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
template<> constexpr rsb_type_t rsb_type_t_for<float> (void) { return RSB_NUMERICAL_TYPE_FLOAT; }
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
template<> constexpr rsb_type_t rsb_type_t_for<int> (void) { return RSB_NUMERICAL_TYPE_INT; }
#endif
template<typename IT> constexpr rsb_flags_t rsb_flags_t_for (void) { return RSBPP_INVALID_FLAGS; }
template<> constexpr rsb_flags_t rsb_flags_t_for<rsb_coo_idx_t> (void) { return RSB_FLAG_NOFLAGS; }
template<> constexpr rsb_flags_t rsb_flags_t_for<rsbpp_half_idx_t> (void) { return RSB_FLAG_USE_HALFWORD_INDICES; }
}

template<typename T>
typename std::enable_if< std::is_scalar<T>::value,T>::type rsbppt_i(void) {return {0};}

template<typename T>
typename std::enable_if<!std::is_scalar<T>::value,T>::type rsbppt_i(void) {return {0,1};}

template<typename T>
const std::vector<T> plus_imag(std::vector<T> v) { for (auto &e:v) e += rsbppt_i<T>(); return v; }

template <typename NT>
class TestSymmetricCoo: public ::testing::Test {
 public:
  static constexpr rsb_flags_t flags { RSB_FLAG_NOFLAGS };
  static constexpr rsb_type_t typecode { rsb_internal::rsb_type_t_for<NT>() };
  const rsb_nnz_idx_t nnzA { 1 };
  static constexpr rsb_coo_idx_t nrA { 1 };
  static constexpr rsb_coo_idx_t ncA { 1 };
  const std::vector<NT> VA { std::vector<NT>(nnzA,4) };
  const std::vector<rsb_coo_idx_t> IA { nrA-1 };
  const std::vector<rsb_coo_idx_t> JA { nrA-1 };
  const std::vector<NT> X { std::vector<NT>(nrA,2) };
  std::vector<NT> Y { std::vector<NT>(nrA,3) };
  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  const rsb_coo_idx_t roff { 0 };
  const rsb_coo_idx_t coff { 0 };
};

using BLASTypes = ::testing::Types<
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
	long double,
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
	std::complex<long double>,
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
	int,
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	std::complex<double>,
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	std::complex<float>,
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	double,
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
	float
#else
	RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE
#endif
>;
TYPED_TEST_CASE(TestSymmetricCoo, BLASTypes);

TYPED_TEST(TestSymmetricCoo, SpMV) {
  const rsb_err_t errval = rsbpp_coo_spmv(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->X.data(), this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( TypeParam{ 3 + 2 * 4 * 2 }, ::testing::Eq( (this->Y)[0] ) );
}

template <typename NT, typename IT=rsb_coo_idx_t>
class TestUpperTriangle {
 public:
  const rsb_nnz_idx_t nnzA { 3 };
  static constexpr rsb_coo_idx_t nrA { 2 };
  static constexpr rsb_coo_idx_t ncA { 2 };
  const std::vector<NT> VA { 11, 12, 22 };
  const std::vector<rsb_coo_idx_t> IP { 0, 2, 3 };
  const std::vector<IT> IA { 0, 0, 1 };
  const std::vector<IT> JA { 0, 1, 1 };
  static constexpr rsb_type_t typecode { rsb_internal::rsb_type_t_for<NT>() };
};

template <typename NT>
class TestUnsymmetric:  public TestUpperTriangle<NT>, public ::testing::Test {
 public:
  static constexpr rsb_flags_t flags { RSB_FLAG_NOFLAGS };
  const std::vector<NT> X { std::vector<NT>(this->nrA,2) };
  std::vector<NT> Y { std::vector<NT>(this->nrA,3) };
  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  const rsb_coo_idx_t roff { 0 };
  const rsb_coo_idx_t coff { 0 };
  const std::vector<NT> R { 95, 91 };
};

TYPED_TEST_CASE(TestUnsymmetric, BLASTypes);

TYPED_TEST(TestUnsymmetric, CsrSpMV) {
  const rsb_err_t errval = rsbpp_csr_spmv(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->X.data(), this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->R)), ::testing::Eq(true) );
}

TYPED_TEST(TestUnsymmetric, CooSpMV) {
  const rsb_err_t errval = rsbpp_coo_spmv(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->X.data(), this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->R)), ::testing::Eq(true) );
}

template <typename NT>
class TestSymmetricSpMV:  public TestUpperTriangle<NT>, public ::testing::Test {
 public:
  static constexpr rsb_flags_t flags { RSB_FLAG_SYMMETRIC };
  const std::vector<NT> X { std::vector<NT>(this->nrA,2) };
  std::vector<NT> Y { std::vector<NT>(this->nrA,3) };
  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  const rsb_coo_idx_t roff { 0 };
  const rsb_coo_idx_t coff { 0 };
  const std::vector<NT> R { 95, 139 };
};

TYPED_TEST_CASE(TestSymmetricSpMV, BLASTypes);

TYPED_TEST(TestSymmetricSpMV, Coo) {
  const rsb_err_t errval = rsbpp_coo_spmv(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->X.data(), this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->R)), ::testing::Eq(true) );
}

TYPED_TEST(TestSymmetricSpMV, Csr) {
  const rsb_err_t errval = rsbpp_csr_spmv(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->X.data(), this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->R)), ::testing::Eq(true) );
}


template <typename IT>
std::vector<IT> rsbpp_bc2br(const std::vector<IT> & a, const size_t nrhs, const size_t lda, const size_t len) {
	std::vector<IT> r {a};

	for (auto ni=0ULL; ni < nrhs; ++ni)
		for (auto li=0ULL; li < len; ++li)
			r[nrhs*li+ni] = a[lda*ni+li];
	return r;
}

class Test_rsbpp_bc2br: public ::testing::Test {
 public:
  const std::vector<int> RbC { 3, 3, 95, 139, 3, 3, 95, 139 };
  const std::vector<int> RbR { 3, 3, 3, 3, 95, 95, 139, 139 };
};

TEST_F(Test_rsbpp_bc2br,basic) {
  const std::vector<int> r { rsbpp_bc2br(this->RbC,2,4,4) };
  EXPECT_THAT( std::equal(std::begin(r), std::end(r), std::begin(this->RbR)), ::testing::Eq(true) );
}


template <typename NT>
class TestSymmetricMultivec:  public TestUpperTriangle<NT>, public ::testing::Test {
 public:
  static constexpr rsb_flags_t flags { RSB_FLAG_SYMMETRIC };
  const rsb_coo_idx_t nrhs {2};
  const rsb_coo_idx_t ld { 4 };
  const std::vector<NT> X { std::vector<NT>(nrhs*this->ld,2) };
  std::vector<NT> Y { std::vector<NT>(nrhs*this->ld,3) };
  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  const rsb_coo_idx_t ldX { ld };
  const rsb_coo_idx_t ldY { ld };
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
};

template <typename NT>
class TestSymmetricMultivecDiagBlock:  public TestSymmetricMultivec<NT> {
 public:
  const rsb_coo_idx_t roff { 2 };
  const rsb_coo_idx_t coff { 2 };
  const std::vector<NT> RbC { 3, 3, 95, 139, 3, 3, 95, 139 };
  const std::vector<NT> RbR { rsbpp_bc2br(this->RbC,this->nrhs,this->ldY,this->ldY) };
};

TYPED_TEST_CASE(TestSymmetricMultivecDiagBlock, BLASTypes);

TYPED_TEST(TestSymmetricMultivecDiagBlock, CsrSpMMByCols) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
  const rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbC)), ::testing::Eq(true) );
}

TYPED_TEST(TestSymmetricMultivecDiagBlock, CooSpMMByCols) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
  const rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbC)), ::testing::Eq(true) );
}

TYPED_TEST(TestSymmetricMultivecDiagBlock, CsrSpMMByRows) {
  const rsb_coo_idx_t ldX { this->ncA + this->coff };
  const rsb_coo_idx_t ldY { this->nrA + this->roff };
  const rsb_coo_idx_t incx { this->nrhs };
  const rsb_coo_idx_t incy { this->nrhs };
  const rsb_flags_t by_rows = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
  rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, ldX, this->X.data(), ldY, this->Y.data(), &this->alpha, incx, incy, this->transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}

TYPED_TEST(TestSymmetricMultivecDiagBlock, CooSpMMByRows_transC) {
  const rsb_trans_t transA { RSB_TRANSPOSITION_C };
  const rsb_coo_idx_t ldX { this->ncA + this->coff };
  const rsb_coo_idx_t ldY { this->nrA + this->roff };
  const rsb_coo_idx_t incx { this->nrhs };
  const rsb_coo_idx_t incy { this->nrhs };
  const rsb_flags_t by_rows = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
  rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, ldX, this->X.data(), ldY, this->Y.data(), &this->alpha, incx, incy,       transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}

TYPED_TEST(TestSymmetricMultivecDiagBlock, CsrSpMMByRows_transC) {
  const rsb_trans_t transA { RSB_TRANSPOSITION_C };
  const rsb_coo_idx_t ldX { this->ncA + this->coff };
  const rsb_coo_idx_t ldY { this->nrA + this->roff };
  const rsb_coo_idx_t incx { this->nrhs };
  const rsb_coo_idx_t incy { this->nrhs };
  const rsb_flags_t by_rows = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
  rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, ldX, this->X.data(), ldY, this->Y.data(), &this->alpha, incx, incy,       transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}

TYPED_TEST(TestSymmetricMultivecDiagBlock, CooSpMMByRows) {
  const rsb_coo_idx_t ldX { this->ncA + this->coff };
  const rsb_coo_idx_t ldY { this->nrA + this->roff };
  const rsb_coo_idx_t incx { this->nrhs };
  const rsb_coo_idx_t incy { this->nrhs };
  const rsb_flags_t by_rows = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
  rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, ldX, this->X.data(), ldY, this->Y.data(), &this->alpha, incx, incy, this->transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}


template <typename NT>
class TestSymmetricMultivecNonDiagBlock:  public TestSymmetricMultivec<NT> {
 public:
  const rsb_coo_idx_t roff { 0 };
  const rsb_coo_idx_t coff { 2 };
  const std::vector<NT> RbC { 95, 91, 47, 139, 95, 91, 47, 139 };
  const std::vector<NT> RbR { rsbpp_bc2br(this->RbC,this->nrhs,this->ldY,this->ldY) };
};

TYPED_TEST_CASE(TestSymmetricMultivecNonDiagBlock, BLASTypes);

TYPED_TEST(TestSymmetricMultivecNonDiagBlock, CooSpMMByRows) {
  const rsb_coo_idx_t ldX { this->ncA + this->coff };
  const rsb_coo_idx_t ldY { this->nrA + this->roff };
  const rsb_coo_idx_t incx { this->nrhs };
  const rsb_coo_idx_t incy { this->nrhs };
  const rsb_flags_t by_rows = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
  rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, ldX, this->X.data(), ldY, this->Y.data(), &this->alpha, incx, incy, this->transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}


template <typename NT>
class TestUnsymmetricMultivec:  public TestUpperTriangle<NT>, public ::testing::Test {
 public:
  static constexpr rsb_flags_t flags { RSB_FLAG_NOFLAGS };
  const rsb_coo_idx_t nrhs {2};
  const rsb_coo_idx_t roff { 2 };
  const rsb_coo_idx_t coff { 2 };
  const std::vector<NT> X { std::vector<NT>(nrhs*(this->ncA+coff),2) };
  std::vector<NT> Y { std::vector<NT>(nrhs*(this->nrA+roff),3) };
  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  const rsb_coo_idx_t ldX { this->ncA + coff };
  const rsb_coo_idx_t ldY { this->nrA + roff };
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  const std::vector<NT> R { 3, 3, 95, 91, 3, 3, 95, 91 };
};

TYPED_TEST_CASE(TestUnsymmetricMultivec, BLASTypes);

TYPED_TEST(TestUnsymmetricMultivec, SpMM_Multi_NRHS) {
  // Note: this is to be expanded systematically
  for ( const rsb_flags_t by_rows : { RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, RSB_FLAG_WANT_ROW_MAJOR_ORDER} )
  for ( const rsb_coo_idx_t nrhs : { 1,2,3,4,8,16,32,50,100 } )
  {
        const std::vector<TypeParam> X (nrhs*(this->ncA+this->coff),2);
        std::vector<TypeParam> Y { std::vector<TypeParam>(nrhs*(this->nrA+this->roff),3) };
        rsb_err_t errval { RSB_ERR_NO_ERROR };

        errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(),       nrhs, this->ldX,       X.data(), this->ldY,       Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, by_rows);
        EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

        errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(),       nrhs, this->ldX,       X.data(), this->ldY,       Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, by_rows);
        EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  }
}

TYPED_TEST(TestUnsymmetricMultivec, CooSpMMByCols) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
  const rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->R)), ::testing::Eq(true) );
}

TYPED_TEST(TestUnsymmetricMultivec, CsrSpMMByCols) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
  const rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->R)), ::testing::Eq(true) );
}

TYPED_TEST(TestUnsymmetricMultivec, CooTransSpMMByCols) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
  const rsb_trans_t transA { RSB_TRANSPOSITION_T };
  const std::vector<TypeParam> R { 3, 3, 47, 139, 3, 3, 47, 139 };
  const rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(R)), ::testing::Eq(true) );
}

TYPED_TEST(TestUnsymmetricMultivec, CsrTransSpMMByCols) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
  const rsb_trans_t transA { RSB_TRANSPOSITION_T };
  const std::vector<TypeParam> R { 3, 3, 47, 139, 3, 3, 47, 139 };
  const rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, transA, this->roff, this->coff, by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(R)), ::testing::Eq(true) );
}

TYPED_TEST(TestUnsymmetricMultivec, CooTransSpMMByRows) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
  const rsb_coo_idx_t ldX { this->nrA + this->roff };
  const rsb_coo_idx_t ldY { this->ncA + this->coff };
  const rsb_coo_idx_t incx { this->nrhs };
  const rsb_coo_idx_t incy { this->nrhs };
  const std::vector<TypeParam> R { 3, 3, 3, 3, 47, 47, 139, 139 };
  const std::vector<rsb_trans_t> transAa {RSB_TRANSPOSITION_T, RSB_TRANSPOSITION_C };
  for (const auto transA: transAa) {
    auto Y {this->Y};
    const rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, ldX, this->X.data(), ldY, Y.data(), &this->alpha, incx, incy, transA, this->roff, this->coff, by_rows);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
    EXPECT_THAT( std::equal(std::begin(Y), std::end(Y), std::begin(R)), ::testing::Eq(true) );
  }
}

TYPED_TEST(TestUnsymmetricMultivec, CsrTransSpMMByRows) {
  const rsb_flags_t by_rows = RSB_FLAG_WANT_ROW_MAJOR_ORDER;
  const rsb_coo_idx_t ldX { this->nrA + this->roff };
  const rsb_coo_idx_t ldY { this->ncA + this->coff };
  const rsb_coo_idx_t incx { this->nrhs };
  const rsb_coo_idx_t incy { this->nrhs };
  const std::vector<rsb_trans_t> transAa {RSB_TRANSPOSITION_T, RSB_TRANSPOSITION_C };
  const std::vector<TypeParam> R { 3, 3, 3, 3, 47, 47, 139, 139 };
  for (const auto transA: transAa) {
    auto Y {this->Y};
    const rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, ldX, this->X.data(), ldY, Y.data(), &this->alpha, incx, incy, transA, this->roff, this->coff, by_rows);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
    EXPECT_THAT( std::equal(std::begin(Y), std::end(Y), std::begin(R)), ::testing::Eq(true) );
  }
}

template <typename NT>
class TestUpperTriangleComplex {
 public:
  const rsb_nnz_idx_t nnzA { 3 };
  static constexpr rsb_coo_idx_t nrA { 2 };
  static constexpr rsb_coo_idx_t ncA { 2 };
  const std::vector<NT> VA { 11, { 12, -1 }, 22 };
  const std::vector<rsb_coo_idx_t> IP { 0, 2, 3 };
  const std::vector<rsb_coo_idx_t> IA { 0, 0, 1 };
  const std::vector<rsb_coo_idx_t> JA { 0, 1, 1 };
  const std::vector<rsb_coo_idx_t>TIP { 0, 1, 3 };
  const std::vector<rsb_coo_idx_t>TIA { 0, 1, 1 };
  const std::vector<rsb_coo_idx_t>TJA { 0, 0, 1 };
  static constexpr rsb_type_t typecode { rsb_internal::rsb_type_t_for<NT>() };
};

template <typename NT>
class TestHermitianMultivecByCols:  public TestUpperTriangleComplex <NT>, public ::testing::Test {
 public:
  static constexpr rsb_flags_t flags { RSB_FLAG_HERMITIAN };
  static constexpr rsb_coo_idx_t nrhs {2};
  const rsb_coo_idx_t roff { 2 };
  const rsb_coo_idx_t coff { 2 };
  const std::vector<NT> X { std::vector<NT>(this->nrhs*(this->ncA+coff),2) };
  std::vector<NT> Y { std::vector<NT>(this->nrhs*(this->nrA+roff),3) };
  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  const rsb_coo_idx_t ldX { this->ncA + coff };
  const rsb_coo_idx_t ldY { this->nrA + roff };
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  const std::vector<NT> RbC { 3, 3, {95,-4}, {139,+4}, 3, 3, {95,-4}, {139,+4} };
  const std::vector<NT> RbR { rsbpp_bc2br(this->RbC,this->nrhs,this->ldY,this->ldY) };
  static constexpr rsb_flags_t by_rows { RSB_FLAG_WANT_COLUMN_MAJOR_ORDER };
};

#if defined(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX) && defined(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX)
using BLASTypesComplex = ::testing::Types<std::complex<double>,std::complex<float>>;
#else
using BLASTypesComplex = ::testing::Types<>;
#endif

TYPED_TEST_CASE(TestHermitianMultivecByCols, BLASTypesComplex );

TYPED_TEST(TestHermitianMultivecByCols, CooSpMMByCols) {
  const rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbC)), ::testing::Eq(true) );
}

TYPED_TEST(TestHermitianMultivecByCols, CsrSpMMByCols) {
  const rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbC)), ::testing::Eq(true) );
}

template <typename NT>
class TestHermitianMultivecByRows:  public TestUpperTriangleComplex <NT>, public ::testing::Test {
 public:
  static constexpr rsb_flags_t flags { RSB_FLAG_HERMITIAN };
  static constexpr rsb_coo_idx_t nrhs {2};
  const rsb_coo_idx_t roff { 2 };
  const rsb_coo_idx_t coff { 2 };
  const std::vector<NT> X { std::vector<NT>(this->nrhs*(this->ncA+coff),2) };
  std::vector<NT> Y { std::vector<NT>(this->nrhs*(this->nrA+roff),3) };
  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { nrhs };
  static constexpr rsb_coo_idx_t incy { nrhs };
  const rsb_coo_idx_t ldX { this->ncA + this->coff };
  const rsb_coo_idx_t ldY { this->nrA + this->roff };
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  static constexpr rsb_flags_t by_rows { RSB_FLAG_WANT_ROW_MAJOR_ORDER };
  const std::vector<NT> RbC { 3, 3, {95,-4}, {139,+4}, 3, 3, {95,-4}, {139,+4} };
  const std::vector<NT> RbR { rsbpp_bc2br(this->RbC,this->nrhs,this->ldY,this->ldY) };
};

TYPED_TEST_CASE(TestHermitianMultivecByRows, BLASTypesComplex );

TYPED_TEST(TestHermitianMultivecByRows, CooSpMMByRows) {
  const rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}

TYPED_TEST(TestHermitianMultivecByRows, CsrSpMMByRows) {
  const rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}

TYPED_TEST(TestHermitianMultivecByRows, CsrSpMMByRowsTransT) {
  const rsb_trans_t transA { RSB_TRANSPOSITION_T };
  const rsb_err_t errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->TIP.data(), this->TJA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, transA, this->roff, this->coff, this->by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}

TYPED_TEST(TestHermitianMultivecByRows, CooSpMMByRowsTransT) {
  const rsb_trans_t transA { RSB_TRANSPOSITION_T };
  const rsb_err_t errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->TIA.data(), this->TJA.data(), this->nrhs, this->ldX, this->X.data(), this->ldY, this->Y.data(), &this->alpha, this->incx, this->incy, transA, this->roff, this->coff, this->by_rows);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::equal(std::begin(this->Y), std::end(this->Y), std::begin(this->RbR)), ::testing::Eq(true) );
}


template <int roff,int coff,int lio,int transA,int extra_flags=RSB_FLAG_NOFLAGS>
struct SpMM_EnumerableParms {
	static constexpr rsb_coo_idx_t Roff{roff};
	static constexpr rsb_coo_idx_t Coff{coff};
	static constexpr rsb_flags_t Ops_layout { lio ? RSB_FLAG_WANT_COLUMN_MAJOR_ORDER : RSB_FLAG_WANT_ROW_MAJOR_ORDER };
	static constexpr rsb_trans_t TransA { transA == 'T' ? RSB_TRANSPOSITION_T : ( transA == 'C' ? RSB_TRANSPOSITION_C : RSB_TRANSPOSITION_N ) };
	static constexpr rsb_flags_t ExtraFlags { extra_flags };
};

using ParmsByRowsTransN = SpMM_EnumerableParms <0,0,0,'N'>;
using ParmsByColsTransN = SpMM_EnumerableParms <0,0,1,'N'>;

using ParmsByRowsTransT = SpMM_EnumerableParms <0,0,0,'T'>;
using ParmsByColsTransT = SpMM_EnumerableParms <0,0,1,'T'>;

using ParmsByRowsTransC = SpMM_EnumerableParms <0,0,0,'C'>;
using ParmsByColsTransC = SpMM_EnumerableParms <0,0,1,'C'>;

#ifdef RSB_NUMERICAL_TYPE_FLOAT
using Index_and_BLAS_Types_Default = std::tuple<float,rsbpp_half_idx_t,ParmsByColsTransN>; // default
#else
using Index_and_BLAS_Types_Default = std::tuple<RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE,rsbpp_half_idx_t,ParmsByColsTransN>; // default
#endif

using Index_and_BLAS_Types = ::testing::Types<
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
	std::tuple<long double,rsb_coo_idx_t,ParmsByRowsTransN>,
	std::tuple<long double,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<long double,rsb_coo_idx_t,ParmsByColsTransN>,
	std::tuple<long double,rsbpp_half_idx_t,ParmsByColsTransN>,
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
	std::tuple<std::complex<long double>,rsb_coo_idx_t,ParmsByRowsTransN>,
	std::tuple<std::complex<long double>,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<std::complex<long double>,rsb_coo_idx_t,ParmsByColsTransN>,
	std::tuple<std::complex<long double>,rsbpp_half_idx_t,ParmsByColsTransN>,
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
	std::tuple<int,rsb_coo_idx_t,ParmsByRowsTransN>,
	std::tuple<int,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<int,rsb_coo_idx_t,ParmsByColsTransN>,
	std::tuple<int,rsbpp_half_idx_t,ParmsByColsTransN>,
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	std::tuple<std::complex<double>,rsb_coo_idx_t,ParmsByRowsTransN>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<std::complex<double>,rsb_coo_idx_t,ParmsByColsTransN>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,ParmsByColsTransN>,
	// variation:
	std::tuple<std::complex<double>,rsb_coo_idx_t,ParmsByRowsTransC>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,ParmsByRowsTransC>,
	std::tuple<std::complex<double>,rsb_coo_idx_t,ParmsByColsTransC>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,ParmsByColsTransC>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<2,0,0,'C'>>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<2,0,1,'C'>>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<2,0,1,'C',RSB_FLAG_SYMMETRIC>>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<0,0,0,'N',RSB_FLAG_SYMMETRIC>>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<0,0,1,'N',RSB_FLAG_SYMMETRIC>>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<2,0,1,'N',RSB_FLAG_SYMMETRIC>>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<0,0,1,'T',RSB_FLAG_SYMMETRIC>>,
	std::tuple<std::complex<double>,rsbpp_half_idx_t,SpMM_EnumerableParms<2,0,1,'T',RSB_FLAG_HERMITIAN>>,
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	std::tuple<double,rsb_coo_idx_t,ParmsByRowsTransN>,
	std::tuple<double,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<double,rsb_coo_idx_t,ParmsByColsTransN>,
	std::tuple<double,rsbpp_half_idx_t,ParmsByColsTransN>,
	// variation:
	std::tuple<double,rsb_coo_idx_t,ParmsByRowsTransT>,
	std::tuple<double,rsbpp_half_idx_t,ParmsByRowsTransT>,
	std::tuple<double,rsb_coo_idx_t,ParmsByColsTransT>,
	std::tuple<double,rsbpp_half_idx_t,ParmsByColsTransT>,
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	std::tuple<std::complex<float>,rsb_coo_idx_t,ParmsByRowsTransN>,
	std::tuple<std::complex<float>,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<std::complex<float>,rsb_coo_idx_t,ParmsByColsTransN>,
	std::tuple<std::complex<float>,rsbpp_half_idx_t,ParmsByColsTransN>,
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
	std::tuple<float,rsb_coo_idx_t,ParmsByRowsTransN>,
	std::tuple<float,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<float,rsb_coo_idx_t,ParmsByColsTransN>,
	std::tuple<float,rsbpp_half_idx_t,ParmsByColsTransN>,
#endif
	Index_and_BLAS_Types_Default
>;

template <typename NT, typename IT=rsb_coo_idx_t>
class SquareTestMatrix {
 public:
  const rsb_nnz_idx_t nnzA { 7 };
  static constexpr rsb_coo_idx_t nrA { 7 };
  static constexpr rsb_coo_idx_t ncA { 7 };
  const std::vector<NT> VA { plus_imag<NT>({11, 12, 13, 14, 15, 16, 17}) };
  const std::vector<rsb_coo_idx_t> IP { 0, 7, 7, 7, 7, 7, 7, 7 };
  const std::vector<IT> IA { 0, 0, 0, 0, 0, 0, 0 };
  const std::vector<IT> JA { 0, 1, 2, 3, 4, 5, 6 };
  static constexpr rsb_type_t typecode { rsb_internal::rsb_type_t_for<NT>() };
};

template <typename TV>
class TestSpMM_OfAnyTypeSquare: public SquareTestMatrix <std::tuple_element_t<0,TV>,std::tuple_element_t<1,TV>>, public ::testing::Test {
 public:
  using NT = typename std::tuple_element_t<0,TV>;
  using IT = typename std::tuple_element_t<1,TV>;
  using Parms = typename std::tuple_element_t<2,TV>;
  using MtxP = SquareTestMatrix <NT,IT>;
  static constexpr rsb_flags_t flags { rsb_internal::rsb_flags_t_for<IT>() | Parms::ExtraFlags };

  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  static constexpr rsb_coo_idx_t roff { Parms::Roff };
  static constexpr rsb_coo_idx_t coff { Parms::Coff };
  static constexpr auto by_roc { Parms::Ops_layout };
  static constexpr rsb_trans_t transA { Parms::TransA };
  static constexpr auto nrO { transA == RSB_TRANSPOSITION_N ? MtxP::nrA : MtxP::ncA };
  static constexpr auto ncO { transA == RSB_TRANSPOSITION_N ? MtxP::ncA : MtxP::nrA };
  static constexpr auto ld { std::max(roff,coff) + std::max(MtxP::nrA,MtxP::ncA) };
};

TYPED_TEST_CASE(TestSpMM_OfAnyTypeSquare, Index_and_BLAS_Types);

TYPED_TEST(TestSpMM_OfAnyTypeSquare, SpMM) {
  using NT = typename TestFixture::NT;

  for ( const NT alpha : { 1,2 } )
  for ( const rsb_coo_idx_t nrhs : { 1,2,3,4,8,16,32,50,100 } )
  {
    const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const std::vector<NT> X (nrhs*(this->ld),2);
    std::vector<NT> Y1 { std::vector<NT>(nrhs*(this->ld),3) };
    std::vector<NT> Y2 { Y1 };
    rsb_err_t errval { RSB_ERR_NO_ERROR };

    errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y1.data(), &alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y2.data(), &alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    EXPECT_THAT( std::equal(std::begin(Y1), std::end(Y1), std::begin(Y2)), ::testing::Eq(true) );
  }
}


template <typename TV>
class TestSpMM_OfAnyTypeUpper: public TestUpperTriangle<std::tuple_element_t<0,TV>,std::tuple_element_t<1,TV>>, public ::testing::Test {
 public:
  using NT = typename std::tuple_element_t<0,TV>;
  using IT = typename std::tuple_element_t<1,TV>;
  using Parms = typename std::tuple_element_t<2,TV>;
  using MtxP = TestUpperTriangle<NT,IT>;
  static constexpr rsb_flags_t flags { rsb_internal::rsb_flags_t_for<IT>() | Parms::ExtraFlags };

  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  static constexpr rsb_coo_idx_t roff { Parms::Roff };
  static constexpr rsb_coo_idx_t coff { Parms::Coff };
  static constexpr auto by_roc { Parms::Ops_layout };
  static constexpr rsb_trans_t transA { Parms::TransA };
  static constexpr auto nrO { transA == RSB_TRANSPOSITION_N ? MtxP::nrA : MtxP::ncA };
  static constexpr auto ncO { transA == RSB_TRANSPOSITION_N ? MtxP::ncA : MtxP::nrA };
  static constexpr auto ld { std::max(roff,coff) + std::max(MtxP::nrA,MtxP::ncA) };
};

TYPED_TEST_CASE(TestSpMM_OfAnyTypeUpper, Index_and_BLAS_Types);

TYPED_TEST(TestSpMM_OfAnyTypeUpper, SpMM) {
  using NT = typename TestFixture::NT;

  for ( const rsb_coo_idx_t nrhs : { 1,2,3,4,8,16,32,50,100 } )
  {
    const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const std::vector<NT> X (nrhs*(this->ld),2);
    std::vector<NT> Y1 { std::vector<NT>(nrhs*(this->ld),3) };
    std::vector<NT> Y2 { Y1 };
    rsb_err_t errval { RSB_ERR_NO_ERROR };

    errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y1.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y2.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    EXPECT_THAT( std::equal(std::begin(Y1), std::end(Y1), std::begin(Y2)), ::testing::Eq(true) );
  }
}


template <typename NT, typename IT=rsb_coo_idx_t>
class TestUpperTriangle2 {
 public:
  const rsb_nnz_idx_t nnzA { 4 };
  static constexpr rsb_coo_idx_t nrA { 4 };
  static constexpr rsb_coo_idx_t ncA { 4 };
  const std::vector<NT> VA { plus_imag<NT>({ 12, 13, 14, 44 }) };
  const std::vector<rsb_coo_idx_t> IP { 0, 3, 3, 3, 4 };
  const std::vector<IT> IA { 0, 0, 0, 3 };
  const std::vector<IT> JA { 1, 2, 3, 3 };
  static constexpr rsb_type_t typecode { rsb_internal::rsb_type_t_for<NT>() };
};

template <typename TV>
class TestSpMM_OfAnyTypeUpper2: public TestUpperTriangle2<std::tuple_element_t<0,TV>,std::tuple_element_t<1,TV>>, public ::testing::Test {
 public:
  using NT = typename std::tuple_element_t<0,TV>;
  using IT = typename std::tuple_element_t<1,TV>;
  using Parms = typename std::tuple_element_t<2,TV>;
  using MtxP = TestUpperTriangle2<NT,IT>;
  static constexpr rsb_flags_t flags { rsb_internal::rsb_flags_t_for<IT>() | Parms::ExtraFlags };

  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  static constexpr rsb_coo_idx_t roff { Parms::Roff };
  static constexpr rsb_coo_idx_t coff { Parms::Coff };
  static constexpr auto by_roc { Parms::Ops_layout };
  static constexpr rsb_trans_t transA { Parms::TransA };
  static constexpr auto nrO { transA == RSB_TRANSPOSITION_N ? MtxP::nrA : MtxP::ncA };
  static constexpr auto ncO { transA == RSB_TRANSPOSITION_N ? MtxP::ncA : MtxP::nrA };
  static constexpr auto ld { std::max(roff,coff) + std::max(MtxP::nrA,MtxP::ncA) };
};

TYPED_TEST_CASE(TestSpMM_OfAnyTypeUpper2, Index_and_BLAS_Types);

TYPED_TEST(TestSpMM_OfAnyTypeUpper2, SpMM) {
  using NT = typename TestFixture::NT;

  for ( const rsb_coo_idx_t nrhs : { 1,2,3,4,8,16,32,50,100 } )
  {
    const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const std::vector<NT> X (nrhs*(this->ld),2);
    std::vector<NT> Y1 { std::vector<NT>(nrhs*(this->ld),3) };
    std::vector<NT> Y2 { Y1 };
    rsb_err_t errval { RSB_ERR_NO_ERROR };

    errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y1.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y2.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    EXPECT_THAT( std::equal(std::begin(Y1), std::end(Y1), std::begin(Y2)), ::testing::Eq(true) );
  }
}


template <typename NT, typename IT=rsb_coo_idx_t>
class TestLowerTriangle2 {
 public:
  const rsb_nnz_idx_t nnzA { 5 };
  static constexpr rsb_coo_idx_t nrA { 5 };
  static constexpr rsb_coo_idx_t ncA { 5 };
  const std::vector<NT> VA { plus_imag<NT>({ 11, 42, 43, 44, 53 }) };
  const std::vector<rsb_coo_idx_t> IP { 0, 1, 1, 1, 4, 5 };
  const std::vector<IT> IA { 0, 3, 3, 3, 4 };
  const std::vector<IT> JA { 0, 1, 2, 3, 2 };
  static constexpr rsb_type_t typecode { rsb_internal::rsb_type_t_for<NT>() };
};

template <typename TV>
class TestSpMM_OfAnyTypeLower2: public TestLowerTriangle2<std::tuple_element_t<0,TV>,std::tuple_element_t<1,TV>>, public ::testing::Test {
 public:
  using NT = typename std::tuple_element_t<0,TV>;
  using IT = typename std::tuple_element_t<1,TV>;
  using Parms = typename std::tuple_element_t<2,TV>;
  using MtxP = TestLowerTriangle2<NT,IT>;
  static constexpr rsb_flags_t flags { rsb_internal::rsb_flags_t_for<IT>() | Parms::ExtraFlags | RSB_FLAG_LOWER};

  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  static constexpr rsb_coo_idx_t roff { Parms::Roff };
  static constexpr rsb_coo_idx_t coff { Parms::Coff };
  static constexpr auto by_roc { Parms::Ops_layout };
  static constexpr rsb_trans_t transA { Parms::TransA };
  static constexpr auto nrO { transA == RSB_TRANSPOSITION_N ? MtxP::nrA : MtxP::ncA };
  static constexpr auto ncO { transA == RSB_TRANSPOSITION_N ? MtxP::ncA : MtxP::nrA };
  static constexpr auto ld { std::max(roff,coff) + std::max(MtxP::nrA,MtxP::ncA) };
};

TYPED_TEST_CASE(TestSpMM_OfAnyTypeLower2, Index_and_BLAS_Types);

TYPED_TEST(TestSpMM_OfAnyTypeLower2, SpMM) {
  using NT = typename TestFixture::NT;

  for ( const rsb_coo_idx_t nrhs : { 1,2,3,4,8,16,32,50,100 } )
  {
    const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const std::vector<NT> X (nrhs*(this->ld),2);
    std::vector<NT> Y1 { std::vector<NT>(nrhs*(this->ld),3) };
    std::vector<NT> Y2 { Y1 };
    rsb_err_t errval { RSB_ERR_NO_ERROR };

    errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y1.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y2.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    EXPECT_THAT( std::equal(std::begin(Y1), std::end(Y1), std::begin(Y2)), ::testing::Eq(true) );
  }
}


template <typename TV>
struct TestMatrixOneByOne {
  using NT = typename std::tuple_element_t<0,TV>;
  using IT = typename std::tuple_element_t<1,TV>;

  static constexpr rsb_flags_t flags { rsb_internal::rsb_flags_t_for<IT>() };
  static constexpr rsb_type_t typecode { rsb_internal::rsb_type_t_for<NT>() };
  const rsb_nnz_idx_t nnzA { 1 };
  static constexpr rsb_coo_idx_t nrA { 1 };
  static constexpr rsb_coo_idx_t ncA { 1 };
  const std::vector<NT> VA { plus_imag<NT>(std::vector<NT>(nnzA,4)) };
  const std::vector<IT> IA { nrA-1 };
  const std::vector<rsb_coo_idx_t> IP { 0,1 };
  const std::vector<IT> JA { nrA-1 };
};

template <typename TV>
class TestSpMM_OfAnyType: public TestMatrixOneByOne<TV>, public ::testing::Test {
 public:
  using NT = typename std::tuple_element_t<0,TV>;
  using IT = typename std::tuple_element_t<1,TV>;
  using Parms = typename std::tuple_element_t<2,TV>;
  using MtxP = TestMatrixOneByOne<TV>;

  const NT alpha {2};
  static constexpr rsb_coo_idx_t incx { 1 };
  static constexpr rsb_coo_idx_t incy { 1 };
  static constexpr rsb_coo_idx_t roff { Parms::Roff };
  static constexpr rsb_coo_idx_t coff { Parms::Coff };
  static constexpr auto by_roc { Parms::Ops_layout };
  static constexpr rsb_trans_t transA { Parms::TransA };
  static constexpr auto nrO { transA == RSB_TRANSPOSITION_N ? MtxP::nrA : MtxP::ncA };
  static constexpr auto ncO { transA == RSB_TRANSPOSITION_N ? MtxP::ncA : MtxP::nrA };
  static constexpr auto ld { std::max(roff,coff) + std::max(MtxP::nrA,MtxP::ncA) };
};

TYPED_TEST_CASE(TestSpMM_OfAnyType, Index_and_BLAS_Types);

#include <rsbpp.hpp> /* TODO: move following tests (using rsbpp_csr_spmx,rsbpp_coo_spmx,Coo) to a different program! */

TYPED_TEST(TestSpMM_OfAnyType, SpMM) {
  using NT = typename TestFixture::NT;

  for ( const rsb_coo_idx_t nrhs : { 1,2,3,4,8,16,32,50,100 } )
  {
    const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ld };
    const std::vector<NT> X (nrhs*(this->ld),2);
    std::vector<NT> Y1 { std::vector<NT>(nrhs*(this->ld),3) };
    std::vector<NT> Y2 { Y1 };
    rsb_err_t errval { RSB_ERR_NO_ERROR };

    errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y1.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y2.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );

    EXPECT_THAT( std::equal(std::begin(Y1), std::end(Y1), std::begin(Y2)), ::testing::Eq(true) );

#if !RSBPP_WANT_MINIMAL_LIBRSBPP
   if( this->roff == 0 && this->coff == 0 )
     if( this->by_roc != RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
     {
       std::vector<NT> Y3 { std::vector<NT>(nrhs*(this->ld),3) };
       using IT = typename TestFixture::IT;
       const NT beta {1};
       const Coo<IT,NT> coo(this->VA,this->IA,this->JA,this->flags);
       const rsb_err_t errval = coo.spmm(this->transA, &this->alpha, nrhs, X.data(), &beta, Y3.data());
       EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_NO_ERROR ) );
       EXPECT_THAT( std::equal(std::begin(Y1), std::end(Y1), std::begin(Y3)), ::testing::Eq(true) );
     }
#endif
  }
}

template <typename TV>
struct ShortTestSpMM_OfAnyType:
  public TestSpMM_OfAnyType<TV> {
  using NT = typename TestSpMM_OfAnyType<TV>::NT;
  using IT = typename TestSpMM_OfAnyType<TV>::IT;
  using Parms = typename TestSpMM_OfAnyType<TV>::Parms;
};

using Index_and_BLAS_Types_for_Errors = ::testing::Types<
#ifdef RSB_NUMERICAL_TYPE_FLOAT
	std::tuple<float,rsbpp_half_idx_t,ParmsByRowsTransN>,
	std::tuple<float,rsb_coo_idx_t,ParmsByColsTransN>,
#endif
	Index_and_BLAS_Types_Default
>;

TYPED_TEST_CASE(ShortTestSpMM_OfAnyType, Index_and_BLAS_Types_for_Errors);

TYPED_TEST(ShortTestSpMM_OfAnyType, SpMM_UnsupportedIncXYError) {
  using NT = typename TestFixture::NT;
  static constexpr rsb_coo_idx_t incx { 2 }; // note: an override
  static constexpr rsb_coo_idx_t incy { 2 }; // note: an override
  static constexpr rsb_coo_idx_t nrhs { 1 };

  const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ncO };
  const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->nrO };
  const std::vector<NT> X (nrhs*(this->ncO),2);
  std::vector<NT> Y { std::vector<NT>(nrhs*(this->nrO),3) };
  rsb_err_t errval { RSB_ERR_NO_ERROR };

  errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, incx, incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_UNSUPPORTED_OPERATION ) );

  errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, incx, incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_UNSUPPORTED_OPERATION ) );
}

TYPED_TEST(ShortTestSpMM_OfAnyType, SpMM_NegativeIncXYError) {
  using NT = typename TestFixture::NT;
  static constexpr rsb_coo_idx_t incx { -1 }; // note: an override
  static constexpr rsb_coo_idx_t incy { -1 }; // note: an override
  static constexpr rsb_coo_idx_t nrhs { 1 };

  const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ncO };
  const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->nrO };
  const std::vector<NT> X (nrhs*(this->ncO),2);
  std::vector<NT> Y { std::vector<NT>(nrhs*(this->nrO),3) };
  rsb_err_t errval { RSB_ERR_NO_ERROR };

  errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, incx, incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );

  errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, incx, incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );
}

TYPED_TEST(ShortTestSpMM_OfAnyType, SpMM_BadNRHS_Error) {
  using NT = typename TestFixture::NT;
  static constexpr rsb_coo_idx_t nrhs { 0 }; // note: an override

  const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ncO };
  const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->nrO };
  const std::vector<NT> X (nrhs*(this->ncO),2);
  std::vector<NT> Y { std::vector<NT>(nrhs*(this->nrO),3) };
  rsb_err_t errval { RSB_ERR_NO_ERROR };

  errval = rsbpp_coo_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );

  errval = rsbpp_csr_spmm(this->typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );
}

TYPED_TEST(ShortTestSpMM_OfAnyType, SpMM_Bad_Type_Error) {
  using NT = typename TestFixture::NT;
  static constexpr rsb_coo_idx_t nrhs { 1 };
  static constexpr rsb_type_t typecode { RSB_NUMERICAL_TYPE_INVALID_TYPE }; // note: an override

  const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ncO };
  const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->nrO };
  const std::vector<NT> X (nrhs*(this->ncO),2);
  std::vector<NT> Y { std::vector<NT>(nrhs*(this->nrO),3) };
  rsb_err_t errval { RSB_ERR_NO_ERROR };

  errval = rsbpp_coo_spmm(typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_UNSUPPORTED_TYPE ) );

  errval = rsbpp_csr_spmm(typecode, this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, this->transA, this->roff, this->coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_UNSUPPORTED_TYPE ) );
}

TYPED_TEST(ShortTestSpMM_OfAnyType, SpMM_Symmetric_non_Square_Error)
 {
  // TODO: reduce further
  using NT = typename TestFixture::NT;
  static constexpr rsb_coo_idx_t nrhs { 1 };

  const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ncO };
  const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->nrO };
  const std::vector<NT> X (nrhs*(this->ncO),2);
  std::vector<NT> Y { std::vector<NT>(nrhs*(this->nrO),3) };
  rsb_err_t errval { RSB_ERR_NO_ERROR };
  const rsb_coo_idx_t roff { 0 };
  const rsb_coo_idx_t coff { 0 };

  errval = rsbpp_coo_spmm(this->typecode, this->flags | RSB_FLAG_SYMMETRIC, this->nnzA, this->nrA+1, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, this->transA, roff, coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );

  errval = rsbpp_csr_spmm(this->typecode, this->flags | RSB_FLAG_SYMMETRIC, this->nnzA, this->nrA+1, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, this->transA, roff, coff, this->by_roc);
  EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );
}


TYPED_TEST(ShortTestSpMM_OfAnyType, SpMM_BadTransaError) {
  using IT = typename TestFixture::IT;
  using NT = typename TestFixture::NT;
  static constexpr rsb_coo_idx_t nrhs { 1 };

  const rsb_coo_idx_t ldX { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->ncO };
  const rsb_coo_idx_t ldY { this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER ? nrhs : this->nrO };
  const std::vector<NT> X (nrhs*(this->ncO),2);
  std::vector<NT> Y { std::vector<NT>(nrhs*(this->nrO),3) };
  rsb_err_t errval { RSB_ERR_NO_ERROR };
  const rsb_trans_t transA { RSB_TRANSPOSITION_INVALID };

  errval = rsbpp_csr_spmx<NT,rsb_coo_idx_t,IT>(this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IP.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, transA, this->roff, this->coff, (this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER) );
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );

  errval = rsbpp_coo_spmx<NT,rsb_coo_idx_t,IT>(this->flags, this->nnzA, this->nrA, this->ncA, this->VA.data(), this->IA.data(), this->JA.data(), nrhs, ldX, X.data(), ldY, Y.data(), &this->alpha, this->incx, this->incy, transA, this->roff, this->coff, (this->by_roc == RSB_FLAG_WANT_ROW_MAJOR_ORDER) );
    EXPECT_THAT( errval, ::testing::Eq( RSB_ERR_BADARGS ) );
}
#endif
