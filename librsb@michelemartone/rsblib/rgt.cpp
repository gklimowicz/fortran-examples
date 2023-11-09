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

#if !RSB_HAVE_GTEST

int main()
{
}

#else

#if RSBP_WANT_CPP20
#include <array>
#endif
#include <vector>
#include <list>
#include <numeric>
#include <algorithm> // fill
#include <limits>
#include <regex>
#include <complex> // no guarantee rsb.hpp provides this
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#define RSBP_WANT_FS (__cplusplus >= 201703)
#if RSBP_WANT_FS
#include <filesystem>
#endif
#define RSBP_TESTING_ONLY 1
#define RSBP_NO_NAMESPACE 1
#include <rsb.hpp>
#if _OPENMP
#include <omp.h>
#endif

#if defined(__llvm__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments" /* for gtest under clang (yes, with '#pragma GCC') */
#endif

using namespace ::testing;

using DefaultType = RSB_DEFAULT_TYPE /* provided by rsb.h */;

struct NoInitTest : public testing::Test {
  void mtxInit(void) const {
    const RsbMatrix<DefaultType> mtx (nullptr,nullptr,nullptr,0);
  }
};

TEST_F(NoInitTest,OkWithLibRSB_Init)
{
  RsbLib rsblib;
  EXPECT_NO_THROW ( this->mtxInit(); );
}

TEST_F(NoInitTest,CatchMissingLibRSB_Init)
{
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW ( this->mtxInit(); );
#else
  EXPECT_ANY_THROW ( this->mtxInit(); );
#endif
}

TEST_F(NoInitTest,CatchAfterLibRSB_Exit)
{
  {
     RsbLib rsblib;
  }
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW ( this->mtxInit(); );
#else
  EXPECT_ANY_THROW ( this->mtxInit(); );
#endif
}


class InitTest : public testing::Test {
 protected:
  RsbLib rsblib;
};

TEST_F(InitTest,meminfo)
{
  EXPECT_THAT( rsblib.meminfo() , Eq( 0 ) );
}

TEST_F(InitTest,reinit_tolerate_no_option)
{
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW ( rsblib.reinit(nullptr ) );
#else
  EXPECT_THAT(      rsblib.reinit<rsb_err_t>(nullptr ) , Eq( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(InitTest,set_opt_str)
{
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW ( rsblib.set_opt_str("RSB_IO_WANT_VERBOSE_TUNING","0") );
#else
  EXPECT_THAT(      rsblib.set_opt_str<rsb_err_t>("RSB_IO_WANT_VERBOSE_TUNING","0") , Eq( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(InitTest,set_opt)
{
  const rsb_int_t v {0};
  rsb_int_t r {-1};
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW ( rsblib.set_opt(RSB_IO_WANT_VERBOSE_TUNING,&v) );
  EXPECT_NO_THROW ( rsblib.get_opt(RSB_IO_WANT_VERBOSE_TUNING,&r) );
#else
  EXPECT_THAT(      rsblib.set_opt<rsb_err_t>(RSB_IO_WANT_VERBOSE_TUNING,&v), Eq(RSB_ERR_NO_ERROR) );
  EXPECT_THAT(      rsblib.get_opt<rsb_err_t>(RSB_IO_WANT_VERBOSE_TUNING,&r), Eq(RSB_ERR_NO_ERROR) );
#endif
  EXPECT_THAT(      r ,  Eq(v) );
}

TEST_F(InitTest,get_opt)
{
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW ( rsblib.get_opt(RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING) );
#else
  EXPECT_NO_THROW ( rsblib.get_opt(RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING) );
#endif
}

TEST_F(InitTest,str_get_opt_err)
{
  rsb_string_t ov;
  const enum rsb_opt_t u_iof {RSB_IO_WANT_VERBOSE_TUNING}; // unsupported flags

#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW ( ov = rsblib.get_opt(u_iof) );
#else
  EXPECT_NO_THROW ( ov = rsblib.get_opt(u_iof) );
#endif
}

TEST_F(InitTest,set_num_threads)
{
#if _OPENMP
  const auto mot = omp_get_max_threads();
#endif
  EXPECT_NO_THROW ( rsblib.set_num_threads( 2) );
  EXPECT_THAT(      rsblib.get_num_threads(), Eq(2) );
#if _OPENMP
  EXPECT_THAT(             omp_get_max_threads(), Eq(mot) );
#endif
  EXPECT_NO_THROW ( rsblib.set_num_threads( 1) );
  EXPECT_THAT(      rsblib.get_num_threads(), Eq(1) );
#if _OPENMP
  EXPECT_THAT(             omp_get_max_threads(), Eq(mot) );
#endif
  EXPECT_NO_THROW ( rsblib.set_num_threads( 0) );
  EXPECT_THAT(      rsblib.get_num_threads(), Ne(0) );
#if _OPENMP
  EXPECT_THAT(             omp_get_max_threads(), Eq(mot) );
#endif
}

TEST_F(InitTest,tolerate_double_initialize)
{
  RsbLib rsblib_double;
}

TEST_F(InitTest,second_verbose_initialize)
{
  RsbLib rsblib_again(true);
}

template <class IT>
const std::vector<IT> ia_to_pa(const std::vector<IT> IA, IT nrA) {
	// indices array to pointers array
	std::vector<IT> PA(nrA+1,0);
	for (auto i: IA)
		PA[i+1]++;
	for (IT i=0; i < nrA; i++)
		PA[i+1]+=PA[i];
	return PA;
}

class LowerTest : public testing::Test {
 protected:
  const RsbLib rsblib;
  static constexpr rsb_nnz_idx_t nnzA { 7 };
  const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 }, incX { 1 }, incY { 1 };
  const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
  const std::vector<rsb_coo_idx_t> PA { ia_to_pa(IA,nrA) };
  const std::vector<double> VA {1,1,1,1,1,1,2};
  const std::vector<double> X {1,1,1,1,1,1};
  std::vector<double> Y{0,0,0,0,0,0};
  const double alpha {2}, beta {1};
  rsb_int_t tn {0};
  RsbMatrix<double> mtx{IA.data(),JA.data(),VA.data(),nnzA,RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS | RSB_FLAG_DUPLICATES_SUM};
  rsb_err_t errval { RSB_ERR_NO_ERROR };
};

TEST_F(LowerTest,getters)
{
  EXPECT_THAT( mtx.nnz(), Eq( nnzA ) );
  EXPECT_THAT( mtx.rows(), Eq( nrA ) );
  EXPECT_THAT( mtx.cols(), Eq( ncA ) );
  EXPECT_THAT( mtx.blocks(), Gt( 0 ) );
}

TEST_F(LowerTest,normOne)
{
  EXPECT_THAT( mtx.normOne() , Eq(3.0) );
}

TEST_F(LowerTest,has_triangular_flags)
{
  EXPECT_THAT( mtx.rsbflags() & RSB_FLAG_TRIANGULAR, Eq(RSB_FLAG_TRIANGULAR) );
}

TEST_F(LowerTest,normInf)
{
  EXPECT_THAT( mtx.normInf() , Eq(3.0) );
}

TEST_F(LowerTest,_info)
{
  const rsb_string_t info = mtx._info();
  const rsb_string_t sbs {"(6 x 6)"};
  EXPECT_THAT( info.size(), Gt(20) );
  EXPECT_THAT( info.size(), Lt(200) );
  EXPECT_THAT( std::search(info.begin(), info.end(), std::begin(sbs), std::end(sbs)) !=info.end(), Eq(true) );
  std::cmatch m;
#if RSB_LIBRSB_VER < 10300
  EXPECT_THAT( std::regex_search(info.c_str(), m, std::regex {"(6 x 6).*(7 nnz, 1.2 nnz/r).*flags 0x204639e.*(coo:1, csr:1, hw:1, ic:1)"}), Eq(true) );
#else
  EXPECT_THAT( std::regex_search(info.c_str(), m, std::regex {"(6 x 6).*(7 nnz, 1.2 nnz/r).*flags 0x204639e.*(coo:1, csr:1, hw:1, ic:1, fi:0)"}), Eq(true) );
#endif
  EXPECT_THAT( mtx.get_info_str("RSB_MIF_MATRIX_INFO__TO__CHAR_P"), Eq(info) );
}

TEST_F(LowerTest,get_info_Rows)
{
  rsb_coo_idx_t nr{};
  errval = mtx.get_info<rsb_err_t>(RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T,&nr);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( nr , Eq( mtx.rows() ) );
}

TEST_F(LowerTest,get_info_Cols)
{
  rsb_coo_idx_t nc{};
  errval = mtx.get_info<rsb_err_t>(RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T,&nc);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( nc , Eq( mtx.cols() ) );
}

TEST_F(LowerTest,get_info_Nnz)
{
  rsb_coo_idx_t nnz{};
  errval = mtx.get_info<rsb_err_t>(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T,&nnz);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( nnz, Eq( mtx.nnz() ) );
}

TEST_F(LowerTest,rsbflags)
{
  rsb_flags_t flags {};
  flags = mtx.rsbflags();
  EXPECT_THAT( RSB_DO_FLAG_HAS(flags, RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS), Eq( 1 ) );
  EXPECT_THAT( RSB_DO_FLAG_HAS(flags, RSB_FLAG_DUPLICATES_SUM), Eq( 1 ) );
}

TEST_F(LowerTest,get_info_coo_t_Rows)
{
  const rsb_coo_idx_t n{mtx.get_info_coo_t(RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T)};
  EXPECT_THAT( n, Eq( mtx.rows() ) );
}

TEST_F(LowerTest,get_info_coo_t_Cols)
{
  const rsb_coo_idx_t n{mtx.get_info_coo_t(RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T)};
  EXPECT_THAT( n, Eq( mtx.cols() ) );
}

TEST_F(LowerTest,get_info_size_t)
{
  const size_t n{mtx.get_info_size_t(RSB_MIF_INDEX_STORAGE_IN_BYTES__TO__SIZE_T)};
  EXPECT_THAT( n, Eq( mtx._get_index_storage_bytes() ) );
  EXPECT_THAT( n, Le( mtx.nnz() * 2 * sizeof(rsb_coo_idx_t) ) );
  EXPECT_THAT( n + 2 * sizeof(double) * mtx.nnz(), Eq( mtx._get_storage_bytes() ) );
}

TEST_F(LowerTest,get_info_blk_t)
{
  const rsb_blk_idx_t n{mtx.get_info_blk_t(RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T)};
  EXPECT_THAT( n, Eq( mtx.blocks() ) );
}

TEST_F(LowerTest,get_info_nnz_t)
{
  const rsb_nnz_idx_t n{mtx.get_info_nnz_t(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T)};
  EXPECT_THAT( n, Eq( mtx.nnz() ) );
}

TEST_F(LowerTest,get_flags_t)
{
  const rsb_flags_t f{mtx.get_flags_t(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T)};
  EXPECT_THAT( f, Eq( mtx.rsbflags() ) );
}

TEST_F(LowerTest,get_type_t)
{
  EXPECT_THAT( mtx.rsbtype(), Eq( mtx.get_type_t(RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T) ) );
}

TEST_F(LowerTest,tune_spmm_with_vectors)
{
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  errval = mtx.tune_spmm<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X.data(),ncA,&beta,Y.data(),nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
}

TEST_F(LowerTest,tune_spmm_with_vectors_and_refs)
{
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  rsb_real_t sf;
  errval = mtx.tune_spmm<rsb_err_t>(sf,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X.data(),beta,Y.data());
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
}

TEST_F(LowerTest,tune_spmm_without_vectors)
{
  errval = mtx.tune_spmm<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(LowerTest,SpMV)
{
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, X.data(), incX, &beta, Y.data(), incY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
}

TEST_F(LowerTest,SpSV)
{
  // same as SpMV
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, X.data(), incX, &beta, Y.data(), incY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );

  auto salpha {1/alpha};
  errval = mtx.spsv<rsb_err_t>(RSB_TRANSPOSITION_N, &salpha, Y.data(), incX, Y.data(), incY);
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 6 ) );
}

#if RSBP_WANT_CPP20
TEST_F(LowerTest,SpMV_Span_Ok)
{
  auto R{X};
  EXPECT_THAT( std::accumulate(R.begin(),R.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, alpha, R, beta, Y);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(R.begin(),R.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
}

TEST_F(LowerTest,SpMV_Span_Error_Catch)
{
  const std::vector<double> R {1};
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (          mtx.spmv(RSB_TRANSPOSITION_N, alpha, R, beta, Y); );
#else
  EXPECT_NO_THROW  ( errval = mtx.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, alpha, R, beta, Y); );
  EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
#endif

  const std::vector<double> L {1,1,1,1,1,1,1};
  errval = mtx.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, alpha, L, beta, Y);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );

  const rsb_trans_t transA = RSBP_INVALID_TRANS_CHAR;
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (          mtx.spmv(transA, alpha, L, beta, Y); );
#else
  EXPECT_NO_THROW  ( errval = mtx.spmv<rsb_err_t>(transA, alpha, L, beta, Y); );
  EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
#endif
}
#endif /* RSBP_WANT_CPP20 */

TEST_F(LowerTest,SpMM)
{
  const rsb_nnz_idx_t ldX{nrA}, ldY{ncA};
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmm<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, nrhs, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, X.data(), ldX, &beta, Y.data(), ldY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
}

TEST_F(LowerTest,SpMM_Alpha_and_Beta_By_Value)
{
  const rsb_nnz_idx_t ldX{nrA}, ldY{ncA};
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmm<rsb_err_t>(RSB_TRANSPOSITION_N, alpha, nrhs, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, X.data(), ldX, beta, Y.data(), ldY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
}

#if RSBP_WANT_IMPLICIT_ORDER
TEST_F(LowerTest,SpMM_Alpha_and_Beta_By_Value_By_Columns)
{
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmm<rsb_err_t>(RSB_TRANSPOSITION_N, alpha, nrhs, X.data(), beta, Y.data());
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
}

TEST_F(LowerTest,SpMM_Alpha_and_Beta_By_Reference_By_Columns)
{
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmm<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, nrhs, X.data(), &beta, Y.data());
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
}

TEST_F(LowerTest,SpMM_NoAlphaNoBetaNoTrans)
{
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx.spmm<rsb_err_t>(Y.data(), X.data(), nrhs, false );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 8 ) );
}
#endif /* RSBP_WANT_IMPLICIT_ORDER */

TEST_F(LowerTest,SpMV_Error_Catch)
{
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (          mtx.spmv(RSB_TRANSPOSITION_N, &alpha, nullptr, incX, &beta, nullptr, incY) );
#else
  errval = mtx.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, nullptr, incX, &beta, nullptr, incY);
  EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(LowerTest,SpMM_Error_Catch)
{
  const rsb_nnz_idx_t ldX{nrA}, ldY{ncA};
  const rsb_trans_t transA = RSBP_INVALID_TRANS_CHAR;
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (          mtx.spmm(transA, &alpha, nrhs, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, X.data(), ldX, &beta, Y.data(), ldY) );
#else
  errval = mtx.spmm<rsb_err_t>(transA, &alpha, nrhs, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, X.data(), ldX, &beta, Y.data(), ldY);
  EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(LowerTest,get_vec)
{
  std::vector<double> D {0,0,0,0,0,0};
  const std::vector<double> O {1,1,1,1,1,1};
  errval = mtx.get_vec<rsb_err_t>(D.data(), RSB_EXTF_DIAG);
  EXPECT_THAT( std::equal(std::begin(D), std::end(D), std::begin(O)), Eq(true) );
}

#if RSB_LIBRSB_VER >= 10300
TEST_F(LowerTest,get_vec_error)
{
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (
           mtx.get_vec(nullptr , RSB_EXTF_DIAG);
  );
#else
  errval = mtx.get_vec<rsb_err_t>(nullptr , RSB_EXTF_DIAG);
#endif
}
#endif

TEST_F(LowerTest,vals__get_vals)
{
  const std::vector<rsb_coo_idx_t> IA {0,5}, JA{0,5};
  const std::vector<double> VA {100,-99};
  const std::vector<double> VC {101,-98};
  std::vector<double> VB {0,0};
  const rsb_nnz_idx_t nnzB {2};
  errval = mtx.set_vals<rsb_err_t>(VA.data(), IA.data(), JA.data(), nnzB, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  errval = mtx.get_vals<rsb_err_t>(VB.data(), IA.data(), JA.data(), nnzB, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( RSB_DO_FLAG_HAS(mtx.rsbflags(), RSB_FLAG_DUPLICATES_SUM), Eq( 1 ) );
  EXPECT_THAT( std::equal(std::begin(VC), std::end(VC), std::begin(VB)), Eq(true) );
}

TEST_F(LowerTest,upd_vals)
{
  errval = mtx.upd_vals<rsb_err_t>(RSB_ELOPF_NEG,nullptr);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  std::vector<double> D {0,0,0,0,0,0};
  const std::vector<double> O {-1,-1,-1,-1,-1,-1};
  errval = mtx.get_vec<rsb_err_t>(D.data(), RSB_EXTF_DIAG);
  EXPECT_THAT( std::equal(std::begin(D), std::end(D), std::begin(O)), Eq(true) );
}

TEST_F(LowerTest,upd_vals_pow)
{
  const double y {2.0};
  errval = mtx.upd_vals<rsb_err_t>(RSB_ELOPF_POW,y);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( mtx.get_val(1,0), Eq( 4.0 ) );
}

TEST_F(LowerTest,get_val)
{
  EXPECT_THAT( mtx.get_val(1,0), Eq( 2.0 ) );
  EXPECT_THAT( mtx.get_val(5,1), Eq( 0.0 ) ); // Note: non-existent nonzero
}

TEST_F(LowerTest,get_rows_sparse)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  rsb_nnz_idx_t rnz {0};
  errval = mtx.get_rows_sparse<rsb_err_t>(transA, nullptr, nullptr, nullptr, nullptr, 0, nrA-1, &rnz, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( rnz, Eq( mtx.nnz() ) );
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (
  mtx.get_rows_sparse(transA, nullptr, nullptr, nullptr, nullptr, 0, nrA-1, nullptr, RSB_FLAG_NOFLAGS );
  );
#else
  errval = mtx.get_rows_sparse<rsb_err_t>(transA, nullptr, nullptr, nullptr, nullptr, 0, nrA-1, nullptr, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
#endif
}

#if RSBP_WANT_CPP20
TEST_F(LowerTest,get_coo_span_ok)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  std::vector<rsb_coo_idx_t> IA (nnzA);
  std::vector<rsb_coo_idx_t> JA (nnzA);
  std::array<double,nnzA> VA;
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW (
  mtx.get_coo(transA, VA, IA, JA, RSB_FLAG_NOFLAGS );
  );
#else
  errval = mtx.get_coo<rsb_err_t>(transA, VA, IA, JA, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
#endif
}
#endif

TEST_F(LowerTest,get_coo_ok)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  std::vector<rsb_coo_idx_t> IA (nnzA);
  std::vector<rsb_coo_idx_t> JA (nnzA);
  std::vector<double> VA (nnzA);
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW (
  mtx.get_coo(transA, VA.data(), IA.data(), JA.data(), RSB_FLAG_NOFLAGS );
  );
#else
  errval = mtx.get_coo<rsb_err_t>(transA, VA.data(), IA.data(), JA.data(), RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(LowerTest,get_coo_err_trans)
{
  std::vector<rsb_coo_idx_t> IA (nnzA);
  std::vector<rsb_coo_idx_t> JA (nnzA);
  std::vector<double> VA (nnzA);
  const rsb_trans_t transA { RSB_TRANSPOSITION_T };
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (
  mtx.get_coo(transA, VA.data(), IA.data(), JA.data(), RSB_FLAG_NOFLAGS);
  );
#else
  errval = mtx.get_coo<rsb_err_t>(transA, VA.data(), IA.data(), JA.data(), RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(LowerTest,get_coo_err_null)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (
  mtx.get_coo(transA, nullptr, nullptr, nullptr, RSB_FLAG_NOFLAGS);
  );
#else
  errval = mtx.get_coo<rsb_err_t>(transA, nullptr, nullptr, nullptr, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(LowerTest,get_csr_ok)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  std::vector<rsb_coo_idx_t> JA (nnzA);
  std::vector<rsb_coo_idx_t> PA (nrA+1);
  std::vector<double> VA (nnzA);
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW (
  mtx.get_csr(transA, VA.data(), PA.data(), JA.data(), RSB_FLAG_NOFLAGS );
  );
#else
  errval = mtx.get_csr<rsb_err_t>(transA, VA.data(), PA.data(), JA.data(), RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
#endif
}

#if RSBP_WANT_CPP20
TEST_F(LowerTest,get_csr_span_ok)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  std::vector<rsb_coo_idx_t> JA (nnzA);
  std::vector<rsb_coo_idx_t> PA (nrA+1);
  std::array<double,nnzA> VA;
#ifndef RSBP_NOTHROW
  EXPECT_NO_THROW (
  mtx.get_csr(transA, VA, PA, JA, RSB_FLAG_NOFLAGS );
  );
#else
  errval = mtx.get_csr<rsb_err_t>(transA, VA, PA, JA, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
#endif
}
#endif

TEST_F(LowerTest,get_csr_err_trans)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_T };
  std::vector<rsb_coo_idx_t> JA (nnzA);
  std::vector<rsb_coo_idx_t> PA (ncA+1);
  std::vector<double> VA (nnzA);
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (
  mtx.get_csr(transA, VA.data(), PA.data(), JA.data(), RSB_FLAG_NOFLAGS);
  );
#else
  errval = mtx.get_csr<rsb_err_t>(transA, VA.data(), PA.data(), JA.data(), RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Eq( RSB_ERR_UNIMPLEMENTED_YET ) );
#endif
}

TEST_F(LowerTest,get_csr_err_null)
{
  const rsb_trans_t transA { RSB_TRANSPOSITION_N };
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (
  mtx.get_csr(transA, nullptr, nullptr, nullptr, RSB_FLAG_NOFLAGS);
  );
#else
  errval = mtx.get_csr<rsb_err_t>(transA, nullptr, nullptr, nullptr, RSB_FLAG_NOFLAGS );
  EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
#endif
}

#if RSBP_WANT_FS
TEST_F(LowerTest,rndr)
{
  const char * pf = "rgt.eps";
  EXPECT_THAT( std::filesystem::is_other(pf) , Eq( false ) );
  errval = mtx.rndr<rsb_err_t>(pf);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::filesystem::exists(pf) , Eq( true ) );
  EXPECT_THAT( std::filesystem::file_size(pf) , Gt( 0 ) );
  std::filesystem::remove_all(pf);
  EXPECT_THAT( std::filesystem::exists(pf) , Eq( false ) );
}
#endif

TEST_F(LowerTest,compare_same)
{
  const auto mtx2 { mtx };
  EXPECT_THAT( ( mtx  == mtx2 ) , Eq( true ) );
  EXPECT_THAT( ( mtx2 == mtx  ) , Eq( true ) );
  EXPECT_THAT( ( mtx  != mtx2 ) , Eq( false ) );
  EXPECT_THAT( ( mtx2 != mtx  ) , Eq( false ) );
  EXPECT_THAT( ( mtx2._equal_to(mtx)  ) , Eq( true) );
  EXPECT_THAT( ( mtx._equal_to(mtx2)  ) , Eq( true) );
}

TEST_F(LowerTest,compare_different)
{
  const double y {2.0};
  auto mtx2 { mtx };
  errval = mtx2.upd_vals<rsb_err_t>(RSB_ELOPF_POW,y);
  EXPECT_THAT( ( mtx  != mtx2 ) , Eq( true ) );
  EXPECT_THAT( ( mtx2 != mtx  ) , Eq( true ) );
  EXPECT_THAT( ( mtx  == mtx2 ) , Eq( false ) );
  EXPECT_THAT( ( mtx2 == mtx  ) , Eq( false ) );
  EXPECT_THAT( ( mtx2._equal_to(mtx)  ) , Eq( false ) );
  EXPECT_THAT( ( mtx._equal_to(mtx2)  ) , Eq( false ) );
}


class ThreadsTuningTest: public LowerTest {
 protected:
  const RsbMatrix<double> mtx{LowerTest::mtx};
};

TEST_F(ThreadsTuningTest, find_threads)
{
  rsb_int_t tn {-1};
  errval = mtx.tune_spmm_threads<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( tn    , Gt( 0 ) );
}

TEST_F(ThreadsTuningTest, auto_threads)
{
  rsb_int_t tn {0};
  errval = mtx.tune_spmm_threads<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( tn    , Eq( 0 ) );
}

TEST_F(ThreadsTuningTest, fixed_threads)
{
  rsb_int_t tn {1};
  errval = mtx.tune_spmm_threads<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( tn    , Eq( 1 ) );
}

TEST_F(ThreadsTuningTest, spsm_on_auto_triangular_flags)
{
  rsb_int_t tn {1};
  
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW (
  mtx.tune_spsm_threads(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA)
  );
#else
  errval = mtx.tune_spsm_threads<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA), Eq( RSB_ERR_NO_ERROR );
#endif
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}

class LowerTriangularSpSV_Test : public testing::Test {
 protected:
  const RsbLib rsblib;
  const rsb_nnz_idx_t nnzA { 7 };
  const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 }, incX { 1 }, incY { 1 };
  const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
  const std::vector<double> VA {1,1,1,1,1,1,2};
  const std::vector<double> X {1,1,1,1,1,1};
  std::vector<double> Y{0,0,0,0,0,0};
  const double alpha {2}, beta {1};
  RsbMatrix<double> mtx{IA.data(),JA.data(),VA.data(),nnzA,RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS | RSB_FLAG_TRIANGULAR};
  rsb_err_t errval { RSB_ERR_NO_ERROR };
};

TEST_F(LowerTriangularSpSV_Test,has_triangular_flags)
{
  EXPECT_THAT( mtx.rsbflags() & RSB_FLAG_TRIANGULAR, Eq(RSB_FLAG_TRIANGULAR) );
}

TEST_F(LowerTriangularSpSV_Test, find_threads)
{
  rsb_int_t tn {-1};
  errval = mtx.tune_spsm_threads<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( tn    , Gt( 0 ) );
}

TEST_F(LowerTriangularSpSV_Test, auto_threads)
{
  rsb_int_t tn {0};
  errval = mtx.tune_spsm_threads<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( tn    , Eq( 0 ) );
}

TEST_F(LowerTriangularSpSV_Test, fixed_threads_spsm)
{
  rsb_int_t tn {1};
  errval = mtx.tune_spsm_threads<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,nullptr,ncA,&beta,nullptr,nrA);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( tn    , Eq( 1 ) );
}


class Types_Test: public testing::Test {
 protected:
  const RsbLib rsblib;
};

TEST_F(Types_Test,isComplex)
{
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
  const RsbMatrix<double> dmtx(1, 1);
  EXPECT_THAT( dmtx._is_complex() , Eq( false ) );
#endif

#ifdef RSB_NUMERICAL_TYPE_FLOAT
  const RsbMatrix<float> smtx(1, 1);
  EXPECT_THAT( smtx._is_complex() , Eq( false ) );
#endif

#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
  const RsbMatrix<std::complex<double>> zmtx(1, 1);
  EXPECT_THAT( zmtx._is_complex() , Eq( true) );
#endif

#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
  const RsbMatrix<std::complex<float>> cmtx(1, 1);
  EXPECT_THAT( cmtx._is_complex() , Eq( true) );
#endif
}

class MatrixConstructors_Test: public LowerTest {
 protected:
 rsb_mtx_t * get_mtxAp(void) { return mtx.mtxAp_; }
};

TEST_F(MatrixConstructors_Test,CopyAndSpMV)
{
  const RsbMatrix<double> mtx2 {mtx};
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx2.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, X.data(), incX, &beta, Y.data(), incY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
  EXPECT_THAT( get_mtxAp() , Ne(nullptr) );
}

TEST_F(MatrixConstructors_Test,Move1AndSpMV)
{
  const RsbMatrix<double> mtx2 {std::move(mtx)};
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx2.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, X.data(), incX, &beta, Y.data(), incY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
  EXPECT_THAT( get_mtxAp() , Eq(nullptr) );
}

TEST_F(MatrixConstructors_Test,Copy2AndSpMV)
{
  const RsbMatrix<double> mtx2 {std::move(mtx),true};
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq( 0 ) );
  errval = mtx2.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, X.data(), incX, &beta, Y.data(), incY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( std::accumulate(X.begin(),X.end(),double{}), Eq( 6 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Ne( 0 ) );
  EXPECT_THAT( std::accumulate(Y.begin(),Y.end(),double{}), Eq(16 ) );
  EXPECT_THAT( get_mtxAp() , Ne(nullptr) );
}

TEST_F(MatrixConstructors_Test,CopyNoTransAndTestEqualOperator)
{
  const RsbMatrix<double> mtx2 {mtx,false};
  EXPECT_THAT( mtx == mtx2, Eq(true) );
  EXPECT_THAT( mtx != mtx2, Eq(false) );
  EXPECT_THAT( get_mtxAp() , Ne(nullptr) );
}

TEST_F(MatrixConstructors_Test,CopyTransAndTestEqualOperator)
{
  const RsbMatrix<double> mtx2 {mtx,true};
  EXPECT_THAT( mtx == mtx2, Eq(false) );
  EXPECT_THAT( mtx != mtx2, Eq(true) );
  EXPECT_THAT( get_mtxAp() , Ne(nullptr) );
}

TEST_F(MatrixConstructors_Test,Move)
{
  const RsbMatrix<double> mtx2 {std::move(mtx)};
  EXPECT_THAT( mtx2.mtxAp_ , Ne(nullptr) );
  EXPECT_THAT(  get_mtxAp() , Eq(nullptr) );
  EXPECT_THAT(  get_mtxAp() != mtx2.mtxAp_, Eq(true) );
}

#ifdef RSB_NUMERICAL_TYPE_FLOAT
TEST_F(MatrixConstructors_Test,InternalsMoveBadType)
{
  const rsb_nnz_idx_t nnzB { 1 };
  const std::vector<rsb_coo_idx_t> IB {0}, JB {0};
  const std::vector<float> VB {1};
  RsbMatrix<float> mtx2{ IB.data(),JB.data(),VB.data(),nnzB,RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS };
  EXPECT_THAT( mtx2.mtxAp_ , Ne(nullptr) );
  auto mtxBp = mtx2.mtxAp_;
  mtx2.mtxAp_ = nullptr;
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW ( mtx.use_mtx_ptr<void>(mtxBp); );
#else
  EXPECT_THAT( mtx.use_mtx_ptr<rsb_err_t>(mtxBp), Ne( RSB_ERR_NO_ERROR) );
#endif
}
#endif

#ifdef RSBP_TESTING_ONLY
TEST_F(MatrixConstructors_Test,InternalsCtorTypeOk)
{
  const RsbMatrix<double> mtx2{ mtx.mtxAp_ };
  mtx.mtxAp_ = nullptr;
  EXPECT_THAT( mtx2.mtxAp_, Ne(nullptr) );
}

TEST_F(MatrixConstructors_Test,InternalsCtorTypeBad)
{
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW ( RsbMatrix<float> mtx2{ mtx.mtxAp_ } );
#else
  EXPECT_ANY_THROW ( RsbMatrix<float> mtx2{ mtx.mtxAp_ } );
#endif
}
#endif

TEST_F(MatrixConstructors_Test,MoveAssignment)
{
  const RsbMatrix<double> mtx2 = {std::move(mtx)};
  EXPECT_THAT( mtx2.mtxAp_ , Ne(nullptr) );
  EXPECT_THAT(  get_mtxAp() , Eq(nullptr) );
  EXPECT_THAT(  get_mtxAp() != mtx2.mtxAp_, Eq(true) );
  EXPECT_ANY_THROW ( if(mtx != mtx2){} ); // mtx is now invalid
  EXPECT_ANY_THROW ( if(mtx == mtx2){} ); // mtx is now invalid
}

#ifndef RSBP_NOTHROW
TEST_F(MatrixConstructors_Test,NoMoveBug)
{
  EXPECT_THAT(  mtx  == mtx  , Eq(true) );
  const RsbMatrix<double> mtx2 = {std::move(mtx)};
  if (true)
  {EXPECT_THAT(  mtx2 == mtx2, Eq(true) );}
  else
  {EXPECT_ANY_THROW ( if(mtx2 == mtx2){} );}
}
#endif

TEST_F(MatrixConstructors_Test,clone)
{
  rsb_mtx_t * mtxBp {};
  mtx.clone(&mtxBp);
  EXPECT_THAT( mtxBp, Ne(nullptr) );
  rsb_mtx_free (mtxBp);
}

TEST_F(MatrixConstructors_Test,CtorFromNull)
{
  const rsb_nnz_idx_t nnzB { 1 }; // nullptr is error only if nnzB>0
  const std::vector<rsb_coo_idx_t> IB {}, JB {};
  const std::vector<double> VB {};
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW ( RsbMatrix<double> mtx2( IB.data(),JB.data(),VB.data(),nnzB,RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS ); );
#else
  EXPECT_ANY_THROW ( RsbMatrix<double> mtx2( IB.data(),JB.data(),VB.data(),nnzB,RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS ); );
#endif
}

#if RSBP_WANT_CPP20
TEST_F(MatrixConstructors_Test,CtorFromSpan)
{
  const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
  const std::vector<double> VA {1,1,1,1,1,1,2};
  RsbMatrix<double> mtx(IA,JA,VA);

  EXPECT_THAT( mtx.nnz(), Eq( 7 ) );
  EXPECT_THAT( mtx.rows(), Eq( 6 ) );
  EXPECT_THAT( mtx.cols(), Eq( 6 ) );
  EXPECT_THAT( mtx.blocks(), Gt( 0 ) );
}
#endif

TEST_F(MatrixConstructors_Test,CooAllowsNoExternalFlags)
{
  // related error: RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS
  EXPECT_ANY_THROW ( RsbMatrix<double> mtx(IA.data(),JA.data(),VA.data(),nnzA,RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS); );
}

TEST_F(MatrixConstructors_Test,FromCsr)
{
  RsbMatrix<double> mtx(nrA,PA.data(),JA.data(),VA.data());
  EXPECT_THAT( mtx.nnz(), Eq( 7 ) );
}

TEST_F(MatrixConstructors_Test,BeginAddClose_set_val)
{
  RsbMatrix<double> mtx2 {nrA, ncA};
  for (auto l = 0; l < nnzA ; ++l)
  	mtx2.set_val(VA[l],IA[l],JA[l]);
  mtx2.close();
  assert ( mtx2.nnz() == nnzA );
  errval = mtx2.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, X.data(), incX, &beta, Y.data(), incY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(MatrixConstructors_Test,BeginAddCloseClose)
{
  RsbMatrix<double> mtx2 {nrA, ncA};
  for (auto l = 0; l < nnzA ; ++l)
  	mtx2.set_val(VA[l],IA[l],JA[l]);
  mtx2.close();
#ifndef RSBP_NOTHROW
  EXPECT_ANY_THROW ( mtx2.close(); );
  EXPECT_NO_THROW ( mtx2.nnz() );
  EXPECT_THAT( mtx2.nnz(), Eq( nnzA ) );
#else
  EXPECT_NO_THROW ( mtx2.close<rsb_err_t>(); ); // bad: shall be no throw
  EXPECT_NO_THROW ( mtx2.nnz() );
  EXPECT_THAT( mtx2.nnz(), Eq( nnzA ) );
#endif
}

TEST_F(MatrixConstructors_Test,BeginAddClose_set_vals)
{
  RsbMatrix<double> mtx2 {nrA, ncA};
  mtx2._add(IA[nnzA-1],JA[nnzA-1],VA[nnzA-1]); // TODO: deprecated set_val
  mtx2.set_val(VA[nnzA-1],IA[nnzA-1],JA[nnzA-1]);
  mtx2.set_vals(VA.data(),IA.data(),JA.data(),nnzA-1,RSB_FLAG_NOFLAGS);
  mtx2.close();
  assert ( mtx2.nnz() == nnzA );
  errval = mtx2.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, &alpha, X.data(), incX, &beta, Y.data(), incY);
  EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}

template <typename NT>
class TestSymmetric: public testing::Test {
 public:
  const RsbLib rsblib;
  const NT iv {1};
  static constexpr rsb_nnz_idx_t nnzA { 1 };
  static constexpr rsb_coo_idx_t nrA { 6 };
  const std::vector<rsb_coo_idx_t> IA { nrA-1 };
  const std::vector<rsb_coo_idx_t> JA { nrA-1 };
  const NT alpha {2}, beta {1};
  const std::vector<NT> VA { std::vector<NT>(nrA,iv) };
  const std::vector<NT> X { std::vector<NT>(nrA,iv) };
  std::vector<NT> Yn { std::vector<NT>(nrA,iv) };
  std::vector<NT> Yt { std::vector<NT>(nrA,iv) };
  const RsbMatrix<NT> mtx { RsbMatrix<NT>(IA.data(),JA.data(),VA.data(),nnzA,RSB_FLAG_SYMMETRIC) };
};

#ifdef RSB_MATRIX_TYPES_LIST_CXX
using AllBLASTypes = ::testing::Types<RSB_MATRIX_TYPES_LIST_CXX>;
#else
using AllBLASTypes = ::testing::Types<
#if RSB_NUMERICAL_TYPE_INT
	int,
#endif
#if RSB_NUMERICAL_TYPE_LONG_DOUBLE
	long double,
#endif
	double,float,
#if RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
	std::complex<long double>,
#endif
	std::complex<double>,std::complex<float>
	>;
#endif

TYPED_TEST_CASE(TestSymmetric, AllBLASTypes);

TYPED_TEST(TestSymmetric, SpMV) {
  this->mtx.spmv(RSB_TRANSPOSITION_N, this->alpha, this->X.data(), this->beta, this->Yn.data());
  EXPECT_THAT( std::equal(std::begin(this->Yn), std::end(this->Yn), std::begin(this->Yt)), Ne(true) );
  this->mtx.spmv(RSB_TRANSPOSITION_T, this->alpha, this->X.data(), this->beta, this->Yt.data());
  EXPECT_THAT( std::equal(std::begin(this->Yn), std::end(this->Yn), std::begin(this->Yt)), Eq(true) );
  EXPECT_THAT( this->Yn[this->nrA-1], Eq(this->beta * this->iv + this->VA[1] * this->alpha * this->iv ) );
}

template <typename NT>
class TestHermitian: public testing::Test {
 public:
  const RsbLib rsblib;
  const NT iv {1};
  static constexpr rsb_nnz_idx_t nnzA { 2 };
  static constexpr rsb_coo_idx_t nrA { 6 };
  static constexpr rsb_coo_idx_t ncA { 6 };
  // NOTE: What if matrix has complex on diagonal ? Behaviour shall be defined there.
  const std::vector<rsb_coo_idx_t> IA { nrA-1, nrA-1 };
  const std::vector<rsb_coo_idx_t> JA { ncA-2, ncA-1 };
  const NT alpha {2}, beta {1};
  const std::vector<NT> VA { {+1,+1},+1 };
  const std::vector<NT> X { std::vector<NT>(ncA,iv) };
  std::vector<NT> Yn { std::vector<NT>(nrA,iv) };
  std::vector<NT> Yt { std::vector<NT>(nrA,iv) };
  const RsbMatrix<NT> mtx { RsbMatrix<NT>(IA.data(),JA.data(),VA.data(),nnzA,RSB_FLAG_HERMITIAN) };
};

#if defined(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX) && defined(RSB_NUMERICAL_TYPE_FLOAT_COMPLEX)
using BLASCmplxTypes = ::testing::Types<std::complex<double>,std::complex<float>>;
#else
using BLASCmplxTypes = ::testing::Types<>;
#endif
TYPED_TEST_CASE(TestHermitian, BLASCmplxTypes );

TYPED_TEST(TestHermitian, SpMV) {
  this->mtx.spmv(RSB_TRANSPOSITION_N, this->alpha, this->X.data(), this->beta, this->Yn.data());
  EXPECT_THAT( std::equal(std::begin(this->Yn), std::end(this->Yn), std::begin(this->Yt)), Ne(true) );
  this->mtx.spmv(RSB_TRANSPOSITION_C, this->alpha, this->X.data(), this->beta, this->Yt.data());
  EXPECT_THAT( std::equal(std::begin(this->Yn), std::end(this->Yn), std::begin(this->Yt)), Eq(true) );
  EXPECT_THAT( this->Yn[this->nrA-2], Eq(this->alpha * ( this->X[this->nrA-2] * conj(this->VA[this->nnzA-2])) + this->beta * this->iv ) );
}

#if RSBP_WANT_CPP20
class SpanInterfaces_Test: public testing::Test {
 protected:
  const RsbLib rsblib;
  static const rsb_nnz_idx_t nnzA { 7 };
  static const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 };
  const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1};
  const rsb_coo_idx_t JA [nnzA] {0,1,2,3,4,5,0};
  const std::vector<double> VA {1,1,1,1,1,1,2};
  const std::vector<double> X{1,1,1,1,1,1};
  std::array<double,nrA> Y;
  const double alpha {2}, beta {1};
  rsb_int_t tn {0};
};

TEST_F(SpanInterfaces_Test,SpMM)
{
  RsbMatrix<double> mtx(IA,JA,VA,nnzA);
#ifndef RSBP_NOTHROW
  mtx.file_save(nullptr);
  mtx.tune_spmm(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X,ncA,beta,Y,nrA);
  mtx.spmv(RSB_TRANSPOSITION_N, alpha, X, beta, Y);
#else
  EXPECT_THAT( mtx.file_save<rsb_err_t>(nullptr), Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( mtx.tune_spmm<rsb_err_t>(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X,ncA,beta,Y,nrA), Eq( RSB_ERR_NO_ERROR ) );
  EXPECT_THAT( mtx.spmv<rsb_err_t>(RSB_TRANSPOSITION_N, alpha, X, beta, Y), Eq( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(SpanInterfaces_Test,SpanSpMM)
{
  RsbMatrix<double> mtx(IA,JA,VA,nnzA);
const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
#ifndef RSBP_NOTHROW
  mtx.spmm(RSB_TRANSPOSITION_N,alpha,nrhs,order,X,beta,Y);
#else
  EXPECT_THAT( mtx.spmm<rsb_err_t>(RSB_TRANSPOSITION_N,alpha,nrhs,order,X,beta,Y), Eq( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(SpanInterfaces_Test,TuneSpMM)
{
  RsbMatrix<double> mtx(IA,JA,VA,nnzA);
  rsb_real_t sf{};
#ifndef RSBP_NOTHROW
  mtx.tune_spmm(sf,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X,ncA,beta,Y,nrA);
#else
  EXPECT_THAT( mtx.tune_spmm<rsb_err_t>(sf,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X,ncA,beta,Y,nrA), Eq( RSB_ERR_NO_ERROR ) );
#endif
}

TEST_F(SpanInterfaces_Test,TuneSpMM_short)
{
  RsbMatrix<double> mtx(IA,JA,VA,nnzA);
  rsb_real_t sf{};
#ifndef RSBP_NOTHROW
  mtx.tune_spmm(sf,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X,beta,Y);
#else
  EXPECT_THAT( mtx.tune_spmm<rsb_err_t>(sf,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X,beta,Y), Eq( RSB_ERR_NO_ERROR ) );
#endif
}
#endif /* RSBP_WANT_CPP20 */


template <typename NT>
class TestEmpty: public testing::Test {
 public:
  const RsbLib rsblib;
  const NT iv {1};
  const rsb_nnz_idx_t nnzA { 0 };
  const rsb_coo_idx_t nrA { 0 };
  const rsb_coo_idx_t ncA { 0 };
  const RsbMatrix<NT> mtx { RsbMatrix<NT>(nullptr,nullptr,nullptr,nnzA) };
};

TYPED_TEST_CASE(TestEmpty, BLASCmplxTypes );

TYPED_TEST(TestEmpty, base) {
  EXPECT_THAT( this->mtx.nnz(), Eq(this->nnzA) );
  EXPECT_THAT( this->mtx.rows(), Eq(this->nrA) );
  EXPECT_THAT( this->mtx.cols(), Eq(this->ncA) );
}

#endif
