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
#include "rsb-config.h"

#if !RSB_WANT_GTEST

int main()
{
}

#else

#include <vector>
#include <string> // std::to_string
#include <stdexcept>
#include <numeric> // std::iota
#include <type_traits>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "rsb.h"
#include "rsb_common.h"
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
#include <complex>
#endif

using namespace ::testing;

std::string mtxPath(const char * mtxFile) {
	/* srcdir is exported by the Makefile */
	std::string filename {};
	if (getenv("srcdir"))
		filename += getenv("srcdir"),
		filename += "/";
	filename += mtxFile;
	return filename;
}
#define MTXPATH(MTXFILE) mtxPath(MTXFILE).c_str()

rsb_err_t on_error(rsb_err_t errval) {
	const auto buflen {256};
	char buf[buflen];

	rsb_strerror_r(errval, buf, buflen);
	if ( (errval) != RSB_ERR_NO_ERROR )
		throw  std::runtime_error(buf);
	return errval;
}

TEST(General,on_an_error)
{
	EXPECT_ANY_THROW( on_error(RSB_ERR_BADARGS) );
}

TEST(General,on_no_error)
{
	EXPECT_NO_THROW( on_error(RSB_ERR_NO_ERROR) );
}

struct RsbLib {
	RsbLib (void) {
		on_error( rsb_lib_init( RSB_NULL_INIT_OPTIONS ) );
	}
	~RsbLib (void) {
		on_error( rsb_lib_exit( RSB_NULL_EXIT_OPTIONS ) );
	}
};


class Test_Init: public ::testing::Test {
	protected:
};


TEST(Test_Init,init_basic)
{
	EXPECT_NO_THROW( {RsbLib rsblib;} );
}


#if RSB_HAVE_SETENV
TEST(Test_Init,init_check_rsb_num_threads)
{
	auto dthf = [] (void) {
		RsbLib rsblib;
		rsb_err_t errval { RSB_ERR_NO_ERROR };
		rsb_int_t wvt {};
		RSB_REINIT_SINGLE_VALUE_GET(RSB_IO_WANT_EXECUTING_THREADS,&wvt,errval);
		return wvt;
	};
	const auto dthc = dthf();
	const auto wtnt = RSB_MAX( (dthc+1) % ( RSB_CONST_MAX_SUPPORTED_THREADS + 1 ), 1 );
	if ( dthc != wtnt ) {
#if defined(RSB_WANT_RSB_NUM_THREADS) && (RSB_WANT_RSB_NUM_THREADS>0) && (wtnt <= RSB_CONST_MAX_SUPPORTED_THREADS) && RSB_HAVE_GETENV
		rsb_err_t errval { RSB_ERR_NO_ERROR };
		rsb_int_t wvt {};

		const auto * rntbv = getenv("RSB_NUM_THREADS");
		setenv("RSB_NUM_THREADS",std::to_string(wtnt).c_str(),1);

		const RsbLib rsblib;
		RSB_REINIT_SINGLE_VALUE_GET(RSB_IO_WANT_EXECUTING_THREADS,&wvt,errval);
		EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
		// RSB_NUM_THREADS has influence on default RSB_IO_WANT_EXECUTING_THREADS
		EXPECT_THAT(   wvt , Eq( wtnt ) );
		EXPECT_THAT(   wvt , Ne( dthc ) );
		if(rntbv)
			setenv("RSB_NUM_THREADS",rntbv,1); // restore
		else
			unsetenv("RSB_NUM_THREADS");
#endif
	}
}
#endif /* RSB_HAVE_SETENV */

TEST(Test_Init,test__rsb_strerror_on_no_error)
{
	const auto buflen {256};
	char buf[buflen];
	const char c = '!';
	buf[0] = c;
	buf[1] = RSB_NUL;
	const rsb_err_t errval = rsb_strerror_r(RSB_ERR_NO_ERROR, buf, buflen);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( buf[0], Eq(c) );
	EXPECT_THAT( strlen(buf), Eq(1) );
}

TEST(Test_Init,test__rsb_strerror_on_errors)
{
	const auto buflen {256};
	char buf[buflen];
	const char c = '!';
	buf[0] = c;
	const rsb_err_t errvals = RSB_ERR_GENERIC_ERROR|RSB_ERR_BADARGS;
	const rsb_err_t errval = rsb_strerror_r(errvals, buf, buflen);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( errvals, RSB_ERR_GENERIC_ERROR );
	EXPECT_THAT( buf[0], Ne(c) );
	EXPECT_THAT( strlen(buf), Gt(10) );
	EXPECT_THAT( strlen(buf), Lt(256) );
}

TEST(Test_Init,test__rsb_strerror_on_error)
{
	const auto buflen {256};
	char buf[buflen];
	const char c = '!';
	buf[0] = c;
	const rsb_err_t errval = rsb_strerror_r(RSB_ERR_GENERIC_ERROR, buf, buflen);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( buf[0], Ne(c) );
	EXPECT_THAT( strlen(buf), Gt(10) );
	EXPECT_THAT( strlen(buf), Lt(256) );
}


class Test_After_Init: public ::testing::Test {
	protected:
	RsbLib rsblib;
};

TEST(Test_After_Init,init_too_many_threads)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };

	rsb_int_t wvt = {RSB_CONST_MAX_SUPPORTED_THREADS + 1}; // setting wrong threads can lead to inconsistent checks checks
	RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_EXECUTING_THREADS,&wvt,errval);
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );

	RSB_REINIT_SINGLE_VALUE_GET(RSB_IO_WANT_EXECUTING_THREADS,&wvt,errval);
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
	EXPECT_THAT(   wvt , Eq( RSB_CONST_MAX_SUPPORTED_THREADS + 1 ) );
}

TEST(Test_After_Init,set_mem_hierarchy_info_0)
{
	RSB_INTERNALS_COMMON_HEAD_DECLS
	EXPECT_THAT(   rsb__set_mem_hierarchy_info(""), Eq( RSB_ERR_NO_ERROR ) );
}

TEST(Test_After_Init,set_mem_hierarchy_info_1)
{
	RSB_INTERNALS_COMMON_HEAD_DECLS
	EXPECT_THAT(   rsb__set_mem_hierarchy_info("L3:12/64/3072K,L2:8/64/256K,L1:8/64/32K"), Eq( RSB_ERR_NO_ERROR ) );
	EXPECT_THAT(   rsb_global_session_handle.memory_hierarchy_levels, Eq(3) );
}

TEST(Test_After_Init,set_mem_hierarchy_info_2)
{
	RSB_INTERNALS_COMMON_HEAD_DECLS
	// TODO: need error code here
	EXPECT_THAT(   rsb__set_mem_hierarchy_info("L2:8/64/256K;L1:8/64/32K"), Eq( RSB_ERR_NO_ERROR ) );
	EXPECT_THAT(   rsb_global_session_handle.memory_hierarchy_levels, Eq(3) );
}

TEST(Test_After_Init,set_mem_hierarchy_info_3)
{
	RSB_INTERNALS_COMMON_HEAD_DECLS
	// TODO: need error code here
	EXPECT_THAT(   rsb__set_mem_hierarchy_info("GARBAGE"), Eq( RSB_ERR_NO_ERROR ) );
	EXPECT_THAT(   rsb_global_session_handle.memory_hierarchy_levels, Eq(0) );
}


struct LowerCooWithDupsAndZeros {
	const rsb_nnz_idx_t nnzA { 9 };
	const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const std::vector<rsb_coo_idx_t> IA {1,1,0,1,2,3,4,5,0}, JA {0,0,0,1,2,3,4,5,1};
	const std::vector<RSB_DEFAULT_TYPE> VA {2,4,1,1,1,1,1,1,0};
	const rsb_blk_idx_t brA = 0, bcA = 0;
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
};


class testLowerCooWithDups: public ::testing::Test, public ::LowerCooWithDupsAndZeros {
	protected:
	RsbLib rsblib { };
	const rsb_flags_t flagsA { RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS | RSB_FLAG_DUPLICATES_SUM };
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval) };

	~testLowerCooWithDups (void) {
		rsb_mtx_free(this->mtxAp_);
	}
};

TEST_F(testLowerCooWithDups,get_vals)
{
	RSB_DEFAULT_TYPE v {-1};
	const rsb_coo_idx_t i {1}, j {0};

	errval = rsb_mtx_get_vals(mtxAp_, &v, &i, &j, 1, RSB_FLAG_NOFLAGS);
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
	EXPECT_THAT( v, Eq( VA[0]+VA[1] ) );
	EXPECT_THAT( mtxAp_->nnz, Eq( nnzA-1 ) ); // with zero, without duplicate
}


class testLowerCooWithZeros: public ::testing::Test, public ::LowerCooWithDupsAndZeros {
	protected:
	RsbLib rsblib { };
	const rsb_flags_t flagsA { RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS | RSB_FLAG_DISCARD_ZEROS};
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval) };

	~testLowerCooWithZeros (void) {
		rsb_mtx_free(this->mtxAp_);
	}
};

TEST_F(testLowerCooWithZeros,get_vals)
{
	RSB_DEFAULT_TYPE v {-1};
	const rsb_coo_idx_t i {1}, j {0};

	errval = rsb_mtx_get_vals(mtxAp_, &v, &i, &j, 1, RSB_FLAG_NOFLAGS);
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
	EXPECT_THAT( v, Eq( VA[1] ) ); // last
	EXPECT_THAT( mtxAp_->nnz, Eq( nnzA-2 ) ); // without zero, without duplicate
}


class testLowerAsRsb {
	protected:
	const rsb_nnz_idx_t nnzA { 7 };
	const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
	const std::vector<RSB_DEFAULT_TYPE> VA {1,1,1,1,1,1,2};
	const std::vector<RSB_DEFAULT_TYPE> X {1,1,1,1,1,1};
	std::vector<RSB_DEFAULT_TYPE> Y{0,0,0,0,0,0};
	const RSB_DEFAULT_TYPE alpha {2}, beta {1};
	rsb_int_t tn {0};
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp_ {nullptr};
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_flags_t flagsA { RSB_FLAG_NOFLAGS };
	RsbLib rsblib { };

	testLowerAsRsb (void) {
		this->mtxAp_ = rsb_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA, flagsA, &errval);
	}

	~testLowerAsRsb (void) {
		rsb_mtx_free(this->mtxAp_);
	}
};

class testLower : public testLowerAsRsb, public ::testing::Test {
	protected:
};

TEST_F(testLower,any)
{
	EXPECT_THAT( this->mtxAp_, Ne( nullptr ) );
}

TEST_F(testLower,begin_wrong_dim)
{
	const rsb_mtx_t * mtxAp {rsb_mtx_alloc_from_coo_begin(nnzA, typecode, RSB_MAX_MATRIX_DIM+1, RSB_MAX_MATRIX_DIM+1, flagsA, &errval)};
	EXPECT_THAT( mtxAp, Eq( nullptr ) );
}

TEST_F(testLower,begin_wrong_type)
{
	const rsb_mtx_t * mtxAp {rsb_mtx_alloc_from_coo_begin(nnzA, RSB_NUMERICAL_TYPE_INVALID_TYPE, nrA, ncA, flagsA, &errval)};
	EXPECT_THAT( mtxAp, Eq( nullptr ) );
}

TEST_F(testLower,begin_wrong_flags)
{
	const rsb_mtx_t * mtxAp {rsb_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA, flagsA | RSB_FLAG_SYMMETRIC | RSB_FLAG_HERMITIAN, &errval)};
	EXPECT_THAT( mtxAp, Eq( nullptr ) );
}

TEST_F(testLower,end)
{
	EXPECT_THAT( rsb__do_mtx_alloc_from_coo_end(&this->mtxAp_), Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(testLower,end_twice)
{
	EXPECT_THAT( rsb__do_mtx_alloc_from_coo_end(&this->mtxAp_), Eq( RSB_ERR_NO_ERROR ) );
	EXPECT_THAT( rsb__do_mtx_alloc_from_coo_end(&this->mtxAp_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testLower,end_wrong_input)
{
	struct rsb_mtx_t mtxA;
	struct rsb_mtx_t * mtxAp {&mtxA};
	RSB_MTX_SET_HBDF(mtxAp);
	EXPECT_THAT( rsb__do_mtx_alloc_from_coo_end(&mtxAp), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testLower,end_null_p)
{
	struct rsb_mtx_t * mtxAp {nullptr};
	EXPECT_THAT( rsb__do_mtx_alloc_from_coo_end(&mtxAp), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testLower,end_null)
{
	EXPECT_THAT( rsb__do_mtx_alloc_from_coo_end(nullptr), Eq( RSB_ERR_BADARGS ) );
}


template <class IT>
std::vector<IT> dense_ia(const IT nr, const IT nc) {
	std::vector<IT> v(nr*nc);
	for (IT i = 0; i < nr * nc ; ++i )
		v[i] = i / nc;
	return v;
}

template <class IT>
std::vector<IT> dense_ja(const IT nr, const IT nc) {
	std::vector<IT> v(nr*nc);
	for (IT i = 0; i < nr * nc ; ++i )
		v[i] = i % nc;
	return v;
}

template <class NT, class IT>
std::vector<NT> dense_va(const IT nr, const IT nc) {
	std::vector<NT> v(nr*nc,1.0);
	return v;
}

static int rsbg_clearenv(void) {
#if RSB_ALLOW_INTERNAL_GETENVS
#if RSB_HAVE_SETENV
	unsetenv("RSB_SPLIT_SF");
	setenv("RSB_WANT_COO2RSB_THREADS","4",1);
#endif /* RSB_HAVE_SETENV */
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
	return 0;
}

class testSquare : public ::testing::Test {
	protected:
	const rsb_coo_idx_t nrA { 16 }, ncA { 16 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_nnz_idx_t nnzA { nrA * ncA };
	const std::vector<rsb_coo_idx_t> IA { dense_ia(nrA,ncA) };
	const std::vector<rsb_coo_idx_t> JA { dense_ja(nrA,ncA) };
	const std::vector<RSB_DEFAULT_TYPE> VA { dense_va<RSB_DEFAULT_TYPE>(nrA,ncA) };
	const rsb_blk_idx_t brA = 0, bcA = 0;
	struct rsb_mtx_t * mtxAp_ {nullptr};
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	rsb_flags_t flagsA { RSB_FLAG_NOFLAGS };
	int env { rsbg_clearenv() }; // with side-effect on RsbLib, to be declared after.
	RsbLib rsblib { };

	testSquare (void) {
		rsb_err_t errval { RSB_ERR_NO_ERROR };
		this->mtxAp_ = rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval);
		on_error( errval );
		rsbg_clearenv();
	}

	~testSquare (void) override {
		rsb_mtx_free(this->mtxAp_);
		rsbg_clearenv();
	}
};

TEST_F(testSquare,no_in_place_switch)
{
	void * coo_VA{};
	rsb_coo_idx_t * coo_IA{}, * coo_JA{};
	const rsb_err_t errval = rsb_mtx_switch_to_coo(this->mtxAp_,&coo_VA,&coo_IA,&coo_JA,RSB_FLAG_NOFLAGS );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(testSquare,mtx)
{
	EXPECT_THAT( this->mtxAp_, Ne( nullptr ) );
}

TEST_F(testSquare,nnz)
{
	EXPECT_THAT( this->mtxAp_->nnz, Eq( this->nnzA ) );
}

TEST_F(testSquare,symmetry)
{
	EXPECT_THAT( ( this->mtxAp_->flags & RSB_FLAG_SYMMETRIC ), Eq( RSB_FLAG_NOFLAGS ) );
}

TEST_F(testSquare,merge)
{
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	auto errval = rsb__leaves_merge_multiple(this->mtxAp_, NULL, NULL, NULL, 0, 1);
	EXPECT_NO_THROW( on_error( errval ) );
	const auto smn1 = mtxAp_->all_leaf_matrices_n;
	EXPECT_THAT( smn0, Gt( smn1) );
}

TEST_F(testSquare,realloc_with_spare_err)
{
	struct rsb_mtx_t * old_mtxAp {this->mtxAp_};
	const rsb_err_t errval = rsb__mtx_realloc_with_spare_leaves(&old_mtxAp, 0);
	EXPECT_THAT( errval, Eq( RSB_ERR_BADARGS ) );
	EXPECT_THAT( old_mtxAp , Eq( nullptr ) );
}

TEST_F(testSquare,realloc_with_spare_in_place)
{
	const rsb_err_t errval = rsb__mtx_realloc_with_spare_leaves(&this->mtxAp_, RSB_TMP_OVERALLOC_MTX*rsb__submatrices(this->mtxAp_));
	EXPECT_NO_THROW( on_error( errval ) );
	const auto smn1 = rsb__submatrices(this->mtxAp_);
}

TEST_F(testSquare,realloc_with_spare_does_realloc)
{
	struct rsb_mtx_t * old_mtxAp { this->mtxAp_ };
	const rsb_err_t errval = rsb__mtx_realloc_with_spare_leaves(&old_mtxAp , RSB_TMP_OVERALLOC_MTX*rsb__submatrices(this->mtxAp_));
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( old_mtxAp , Ne( nullptr ) ); // may realloc to same address, but not NULL
	this->mtxAp_ = old_mtxAp;
}

TEST_F(testSquare,noLeafMultivec)
{
	const rsb_coo_idx_t nrhs {2};
	const std::vector<RSB_DEFAULT_TYPE> X(ncA*nrhs*2,1);
	std::vector<RSB_DEFAULT_TYPE> Y(nrA*nrhs*2,1);
	std::vector<RSB_DEFAULT_TYPE> R(nrA*nrhs,1);
	const rsb_trans_t transA { RSB_DEFAULT_TRANSPOSITION };
	const rsb_nnz_idx_t ldB { 0 };
	const rsb_nnz_idx_t ldC { 0 };
	const rsb_nnz_idx_t inc { 1 };
	RSB_DEFAULT_TYPE alpha {1};
	RSB_DEFAULT_TYPE beta {0};
	rsb_int_t wllmv = {1};
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_LEAF_LEVEL_MULTIVEC ,&wllmv,errval);
	EXPECT_NO_THROW( on_error( errval ) );
	rsb_blk_idx_t lc = -1;
	EXPECT_NO_THROW(rsb_mtx_get_info(mtxAp_, RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T, &lc));
	errval = rsb_spmm(transA, &alpha, mtxAp_, nrhs, RSB_FLAG_WANT_ROW_MAJOR_ORDER, X.data(), ldB, &beta, Y.data(), ldC);
	EXPECT_THAT( errval, Ne( RSB_ERR_UNSUPPORTED_TYPE ) );
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
	errval = rsb_spmv(transA, &alpha, mtxAp_, X.data(), inc, &beta, R.data(), inc);
	EXPECT_THAT( errval, Ne( RSB_ERR_UNSUPPORTED_TYPE ) );
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
	if(errval == RSB_ERR_NO_ERROR)
	for(int i = 0; i < nrA; ++i)
	{
		EXPECT_THAT( Y[i*nrhs+1], Eq( R[i]) );
	}
	wllmv = {0};
	RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_LEAF_LEVEL_MULTIVEC ,&wllmv,errval);
	EXPECT_NO_THROW( on_error( errval ) );
}

class testSquareSplit : public testSquare {
	protected:
	testSquareSplit (void) {
		rsb_err_t errval { RSB_ERR_NO_ERROR };
		errval = rsb__mtx_realloc_with_spare_leaves(&this->mtxAp_, RSB_TMP_OVERALLOC_MTX*rsb__submatrices(this->mtxAp_));
		on_error( errval );
	}
};

TEST_F(testSquareSplit,split)
{
	const int verbosity = 0; // 0..3
	const rsb_submatrix_idx_t manp {0};
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	const auto errval = rsb__mtx_split(this->mtxAp_, manp, NULL, NULL, NULL, verbosity, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	const auto smn1 = mtxAp_->all_leaf_matrices_n;
	EXPECT_THAT( smn1, Gt( smn0) );
	const double mfts = 0.5; /* matrices fraction to split */
	EXPECT_THAT( smn1, Ge( floor( smn0 * mfts ) + floor ( smn0 * mfts ) * 4 ) );
}

#if RSB_HAVE_SETENV
#if RSB_ALLOW_INTERNAL_GETENVS
TEST_F(testSquareSplit,split_two)
{
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	setenv("RSB_SPLIT_SF","2",1);
	auto errval = rsb__mtx_split(this->mtxAp_, 0, NULL, NULL, NULL, 0, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxAp_->all_leaf_matrices_n, Eq( smn0 + ( 4 - 1 ) * 2 ) );
}

TEST_F(testSquareSplit,split_two_verbose)
{
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	setenv("RSB_SPLIT_SF","2",1);
	auto errval = rsb__mtx_split(this->mtxAp_, 0, NULL, NULL, NULL, 3, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxAp_->all_leaf_matrices_n, Eq( smn0 + ( 4 - 1 ) * 2 ) );
}

TEST_F(testSquareSplit,split_half)
{
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	setenv("RSB_SPLIT_SF","0.5",1);
	auto errval = rsb__mtx_split(this->mtxAp_, 0, NULL, NULL, NULL, 0, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxAp_->all_leaf_matrices_n, Eq( smn0 + ( 4 - 1 ) * smn0 / 2 ) );
}

TEST_F(testSquareSplit,split_fraction)
{
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	setenv("RSB_SPLIT_SF","3/12",1);
	auto errval = rsb__mtx_split(this->mtxAp_, 0, NULL, NULL, NULL, 0, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxAp_->all_leaf_matrices_n, Eq( smn0 + ( 4 - 1 ) * ( smn0 / 4 ) ) );
}

TEST_F(testSquareSplit,split_percentage)
{
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	setenv("RSB_SPLIT_SF","12.5%",1);
	auto errval = rsb__mtx_split(this->mtxAp_, 0, NULL, NULL, NULL, 0, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxAp_->all_leaf_matrices_n, Eq( smn0 + ( 4 - 1 ) * ( smn0 / 8 ) ) );
}

TEST_F(testSquareSplit,split_default_on_wrong_env)
{
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	setenv("RSB_SPLIT_SF","-frmrmkmr",1);
	auto errval = rsb__mtx_split(this->mtxAp_, 0, NULL, NULL, NULL, 0, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxAp_->all_leaf_matrices_n, Eq( smn0 + ( 4 - 1 ) * ( smn0 / 2 ) ) );
}
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
#endif /* RSB_HAVE_SETENV */

TEST_F(testSquareSplit,split_manp)
{
	const int verbosity = 0; // 0..3
	const rsb_submatrix_idx_t manp {1}; // subdivide only 1
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	const auto errval = rsb__mtx_split(this->mtxAp_, manp, NULL, NULL, NULL, verbosity, 0);
	EXPECT_NO_THROW( on_error( errval ) );
	const auto smn1 = mtxAp_->all_leaf_matrices_n;
	EXPECT_THAT( smn1, Eq( smn0 + 4 - 1  ) );
}

TEST_F(testSquareSplit,CheckMergeSplit)
{
	struct rsb_mtx_t ** mtxOpp { &this->mtxAp_ };
	const rsb_err_t errval = rsb__mtx_ms_check(mtxOpp);
	EXPECT_NO_THROW( on_error( errval ) );
}


class testSquareTune : public testSquare {
	protected:
	const rsb_int_t maxr { 1 }; // one round: tuning possible
	const rsb_time_t maxt { 0 };
	const rsb_coo_idx_t nrhs { 1 };
	const rsb_trans_t transA { RSB_DEFAULT_TRANSPOSITION };
	const rsb_flags_t order { RSB_FLAG_WANT_COLUMN_MAJOR_ORDER };
	const rsb_nnz_idx_t ldB { 0 };
	const rsb_nnz_idx_t ldC { 0 };

	testSquareTune(void) {
		rsb_err_t errval { RSB_ERR_NO_ERROR };
		errval = rsb__mtx_realloc_with_spare_leaves(&this->mtxAp_, RSB_TMP_OVERALLOC_MTX*rsb__submatrices(this->mtxAp_));
		on_error( errval );
	}
};

TEST_F(testSquareTune,noRounds)
{
	const rsb_int_t maxr = 0; // no rounds: only performance-measuring benchmark
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	rsb_real_t sf {-1};

	const rsb_err_t errval = rsb_tune_spmm(NULL, &sf, NULL, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_NO_THROW( on_error( errval ) );
	const auto smn1 = mtxAp_->all_leaf_matrices_n;
	EXPECT_THAT( smn1, Eq( smn0 ) );
	EXPECT_THAT( sf, Eq( 1 ) ); // no tuning => no speedup
}

TEST_F(testSquareTune,ZeroNrhsIsOk)
{
	const rsb_int_t maxr = 0; // no rounds: only performance-measuring benchmark
	const rsb_coo_idx_t nrhs = 0; // defaults to 1
	const auto smn0 = mtxAp_->all_leaf_matrices_n;
	rsb_real_t sf {-1};

	const rsb_err_t errval = rsb_tune_spmm(NULL, &sf, NULL, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_NO_THROW( on_error( errval ) );
	const auto smn1 = mtxAp_->all_leaf_matrices_n;
	EXPECT_THAT( smn1, Eq( smn0 ) );
	EXPECT_THAT( sf, Eq( 1 ) );
}

TEST_F(testSquareTune,ValidDefault)
{
	const rsb_err_t errval = rsb_tune_spmm(NULL, NULL, NULL, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);
	EXPECT_NO_THROW( on_error( errval ) );
}

TEST_F(testSquareTune,ErrorInvalidTrans)
{
	const rsb_trans_t transA = RSB_INVALID_TRANS;
	const rsb_err_t errval = rsb_tune_spmm(NULL, NULL, NULL, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);
	EXPECT_ANY_THROW( on_error( errval ) );
}

TEST_F(testSquareTune,ErrorNoMatrix)
{
	const rsb_mtx_t * mtxAp { nullptr };
	struct rsb_mtx_t ** mtxOpp { nullptr };
	const rsb_err_t errval = rsb_tune_spmm(mtxOpp, NULL, NULL, maxr, maxt, transA, NULL, mtxAp, nrhs, order, NULL, ldB, NULL, NULL, ldC);
	EXPECT_ANY_THROW( on_error( errval ) );
}

TEST_F(testSquareTune,ErrorBothMatrixArgsForbidden)
{
	struct rsb_mtx_t ** mtxOpp { &mtxAp_ };
	const rsb_err_t errval = rsb_tune_spmm(mtxOpp, NULL, NULL, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);
	EXPECT_ANY_THROW( on_error( errval ) );
}

TEST_F(testSquareTune,ErrorDuplicateMatrix)
{
	struct rsb_mtx_t ** mtxOpp { &mtxAp_ };
	const rsb_err_t errval = rsb_tune_spmm(mtxOpp, NULL, NULL, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);
	EXPECT_ANY_THROW( on_error( errval ) );
}

TEST_F(testSquareTune,oneRound_RowMajor)
{
	const rsb_flags_t order { RSB_FLAG_WANT_ROW_MAJOR_ORDER };
	struct rsb_mtx_t * mtxOp { mtxAp_ };
	rsb_real_t sf {-1};

	const rsb_err_t errval = rsb_tune_spmm(&mtxOp, &sf, NULL, maxr, maxt, transA, NULL, NULL, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxOp, Ne( nullptr ) );
	if(mtxOp != mtxAp_)
	{
		EXPECT_THAT( sf, Gt( 1.0 ) );
		mtxAp_ = mtxOp; // mtxAp_ might have been freed
	}
	else
	{
		EXPECT_THAT( sf, Eq( 1.0 ) );
	}
}

TEST_F(testSquareTune,oneRound)
{
	struct rsb_mtx_t * mtxOp { mtxAp_ };
	rsb_real_t sf {-1};

	const rsb_err_t errval = rsb_tune_spmm(&mtxOp, &sf, NULL, maxr, maxt, transA, NULL, NULL, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxOp, Ne( nullptr ) );
	if(mtxOp != mtxAp_)
	{
		EXPECT_THAT( sf, Gt( 1.0 ) );
		mtxAp_ = mtxOp; // mtxAp_ might have been freed
	}
	else
	{
		EXPECT_THAT( sf, Eq( 1.0 ) );
	}
}

TEST_F(testSquareTune,with_excessive_threads)
{
	struct rsb_mtx_t * mtxOp { mtxAp_ };
	rsb_int_t tn = RSB_CONST_MAX_SUPPORTED_THREADS + 1;
	const rsb_int_t otn { tn };
	rsb_real_t sf {-1};
	rsb_int_t wvt = {0}; // set to >0 for verbose
	rsb_err_t errval;

	RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_VERBOSE_TUNING,&wvt,errval);
	EXPECT_NO_THROW( on_error( errval ) );
	errval = rsb_tune_spmm(&mtxOp, &sf, &tn, maxr, maxt, transA, NULL, NULL, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_THAT( tn, Eq( otn ) );
	EXPECT_THAT( errval, Eq( RSB_ERR_BADARGS	) );
	EXPECT_THAT( mtxOp, Eq( mtxAp_ ) );
	EXPECT_THAT( sf, Eq( -1.0 ) );
}

TEST_F(testSquareTune,with_specific_threads)
{
	struct rsb_mtx_t * mtxOp { mtxAp_ };
	rsb_int_t tn = RSB_MIN(2,rsb_get_num_threads()); // want specific threads count (as long as not illegal)
	const rsb_int_t otn { tn };
	rsb_real_t sf {-1};
	rsb_int_t wvt = {0}; // set to >0 for verbose
	rsb_err_t errval;

	RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_VERBOSE_TUNING,&wvt,errval);
	EXPECT_NO_THROW( on_error( errval ) );
	errval = rsb_tune_spmm(&mtxOp, &sf, &tn, maxr, maxt, transA, NULL, NULL, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_THAT( tn, Eq( otn ) );
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( mtxOp, Ne( nullptr ) );
	if(mtxOp != mtxAp_)
	{
		EXPECT_THAT( sf, Gt( 1.0 ) );
		mtxAp_ = mtxOp; // mtxAp_ might have been freed
	}
	else
	{
		EXPECT_THAT( sf, Eq( 1.0 ) );
	}
}

TEST_F(testSquareTune,find_threads)
{
	rsb_int_t tn = -1; // find threads
	rsb_real_t sf {-1};
	rsb_int_t wvt = {0}; // set to >0 for verbose
	rsb_err_t errval;

	RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_VERBOSE_TUNING,&wvt,errval);

	errval = rsb_tune_spmm(NULL, &sf, &tn, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( tn, Gt( 0 ) );
}

TEST_F(testSquareTune,sample_only)
{
	rsb_int_t tn = 0; // default threads
	rsb_real_t sf {-1};
	rsb_int_t wvt = {0}; // set to >0 for verbose
	rsb_err_t errval;

	RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_VERBOSE_TUNING,&wvt,errval);

	errval = rsb_tune_spmm(NULL, &sf, &tn, maxr, maxt, transA, NULL, mtxAp_, nrhs, order, NULL, ldB, NULL, NULL, ldC);

	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( tn, Eq( 0 ) );
}

TEST_F(testSquareTune,switch_recursive_in_place_matrix_to_in_place_coo_sorted)
{
	struct rsb_coo_mtx_t coo;
	const rsb_err_t errval = rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_sorted(mtxAp_, &coo);
	EXPECT_NO_THROW( on_error( errval ) );
	mtxAp_ = NULL;
	rsb__destroy_coo_matrix_t(&coo);
}


class testCsr2Coo: public ::testing::Test {
	protected:
	const rsb_coo_idx_t nr { 4 };
	const rsb_nnz_idx_t nnz { 5 };
};

TEST_F(testCsr2Coo,switch_compressed_array_to_fullword_coo)
{
	std::vector<rsb_coo_idx_t> PA { 0,1,2,3,nnz  };
	const std::vector<rsb_coo_idx_t> IA { 0,1,2,3,3 };
	const rsb_err_t errval = rsb__do_switch_compressed_array_to_fullword_coo(PA.data(), nr, 0, PA.data());
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( PA[4], Eq( IA[4]) );
}

TEST_F(testCsr2Coo,switch_fullword_array_to_compressed)
{
	std::vector<rsb_coo_idx_t> IA { 0,1,2,3,3 };
	const std::vector<rsb_coo_idx_t> PA { 0,1,2,3,nnz  };
	const rsb_err_t errval = rsb__do_switch_fullword_array_to_compressed(IA.data(), nnz, nr);
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( PA[4], Eq( IA[4]) );
}


class testCsr2Rsb: public ::testing::Test {
	protected:
	const rsb_coo_idx_t nr { 4 };
	const rsb_coo_idx_t br { 1 };
	const rsb_nnz_idx_t nnz { 5 };
	const std::vector<rsb_coo_idx_t> JA { 0,1,2,3,3 };
	const std::vector<rsb_coo_idx_t> PA { 0,1,2,3,nnz  };
	const std::vector<RSB_DEFAULT_TYPE> VA { 1,1,2,3,4  };
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_flags_t flags { RSB_FLAG_NOFLAGS };
	RsbLib rsblib { };
};

TEST_F(testCsr2Rsb,ok)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp = rsb_mtx_alloc_from_csr_const(VA.data(), PA.data(), JA.data(), nnz, typecode, nr, nr, br, br, flags, &errval);
	EXPECT_THAT( mtxAp, Ne( nullptr ) );
	EXPECT_NO_THROW( on_error( errval ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testCsr2Rsb,EnoMem)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp {};

	mtxAp = rsb_mtx_alloc_from_csr_const(NULL, PA.data(), JA.data(), nnz, typecode, nr, nr, br, br, flags, &errval);
	EXPECT_THAT( mtxAp, Eq( nullptr ) );
	EXPECT_THAT( errval, Eq(RSB_ERR_ENOMEM) );
}

TEST_F(testCsr2Rsb,badFlags)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	const rsb_flags_t flags { RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS };
	struct rsb_mtx_t * mtxAp = rsb_mtx_alloc_from_csr_const(VA.data(), PA.data(), JA.data(), nnz, typecode, nr, nr, br, br, flags, &errval);
	EXPECT_THAT( mtxAp, Eq( nullptr ) );
	EXPECT_THAT( errval, Eq( RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS ) );
}

TEST_F(testCsr2Rsb,badSymFlags)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	const rsb_flags_t flags { RSB_FLAG_SYMMETRIC | RSB_FLAG_HERMITIAN };
	struct rsb_mtx_t * mtxAp = rsb_mtx_alloc_from_csr_const(VA.data(), PA.data(), JA.data(), nnz, typecode, nr, nr, br, br, flags, &errval);
	EXPECT_THAT( mtxAp, Eq( nullptr ) );
	EXPECT_THAT( errval, Eq( RSB_ERR_BADARGS ) );
}


class testWeedDuplicates: public ::testing::Test {
	protected:
};

TEST_F(testWeedDuplicates,fromUnsorted)
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 3+1;
	rsb_coo_idx_t IA[] = {3,2,2,1};
	rsb_coo_idx_t JA[] = {2,1,1,0};
	RSB_DEFAULT_TYPE VA[] = {9,3,2,0};
	const rsb_nnz_idx_t nnnz = rsb__weed_out_duplicates(IA, JA, VA, nnz, typecode, RSB_FLAG_NOFLAGS);

	EXPECT_THAT( nnnz , Eq( 3 ) );
	EXPECT_THAT( IA[2], Eq( 1 ) );
	EXPECT_THAT( JA[2], Eq( 0 ) );
	EXPECT_THAT( VA[2], Eq( 0 ) );
}

TEST_F(testWeedDuplicates,fromSorted)
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 3+1;
	rsb_coo_idx_t IA[] = {1,2,2,3};
	rsb_coo_idx_t JA[] = {0,1,1,2};
	RSB_DEFAULT_TYPE VA[] = {0,3,2,9};
	const rsb_nnz_idx_t nnnz = rsb__weed_out_duplicates(IA, JA, VA, nnz, typecode, RSB_FLAG_SORTED_INPUT);

	EXPECT_THAT( nnnz , Eq( 3 ) );
	EXPECT_THAT( IA[2], Eq( 3 ) );
	EXPECT_THAT( JA[2], Eq( 2 ) );
	EXPECT_THAT( VA[2], Eq( 9 ) );
}


template <typename NT=RSB_DEFAULT_TYPE, int TC=RSB_NUMERICAL_TYPE_DEFAULT>
struct testRsbAsCxxBaseT {
	const rsb_coo_idx_t nrA { 16 }, ncA { 16 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_nnz_idx_t nnzA { nrA * ncA };
	const std::vector<rsb_coo_idx_t> IA { dense_ia(nrA,ncA) };
	const std::vector<rsb_coo_idx_t> JA { dense_ja(nrA,ncA) };
	const std::vector<NT> VA { dense_va<NT>(nrA,ncA) };
	const rsb_blk_idx_t brA = 0, bcA = 0;
	const rsb_type_t typecode { TC };
	RsbLib rsblib { };
};

class testRsbAsCxxBase {
	protected:
	const rsb_coo_idx_t nrA { 16 }, ncA { 16 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_nnz_idx_t nnzA { nrA * ncA };
	const std::vector<rsb_coo_idx_t> IA { dense_ia(nrA,ncA) };
	const std::vector<rsb_coo_idx_t> JA { dense_ja(nrA,ncA) };
	const std::vector<RSB_DEFAULT_TYPE> VA { dense_va<RSB_DEFAULT_TYPE>(nrA,ncA) };
	const rsb_blk_idx_t brA = 0, bcA = 0;
	rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	RsbLib rsblib { };
};


class testRsbAsCxx : public testRsbAsCxxBase, public ::testing::Test {
	protected:
};

TEST_F(testRsbAsCxx,mtxCreateRsb)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_NOFLAGS, &errval) };
	EXPECT_THAT( mtxAp, Ne( nullptr ) );
	EXPECT_THAT( rsb__submatrices(mtxAp), Gt( 1 ) );
	rsb_mtx_free(mtxAp);
}

TEST_F(testRsbAsCxx,mtxCreateCsr)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS, &errval) };
	EXPECT_THAT( mtxAp, Ne( nullptr ) );
	EXPECT_THAT( rsb__submatrices(mtxAp), 1 );
	rsb_mtx_free(mtxAp);
}

TEST_F(testRsbAsCxx,mtxCreateCoo)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS, &errval) };
	EXPECT_THAT( mtxAp, Ne( nullptr ) );
	EXPECT_THAT( rsb__submatrices(mtxAp), 1 );
	rsb_mtx_free(mtxAp);
}

TEST_F(testRsbAsCxx,mtxCreateHalfwordCsr)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_USE_HALFWORD_INDICES_CSR, &errval) };
	EXPECT_THAT( mtxAp, Ne( nullptr ) );
	EXPECT_THAT( rsb__submatrices(mtxAp), 1 );
	rsb_mtx_free(mtxAp);
}

TEST_F(testRsbAsCxx,mtxCreateHalfwordCoo)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_USE_HALFWORD_INDICES_COO, &errval) };
	EXPECT_THAT( mtxAp, Ne( nullptr ) );
	EXPECT_THAT( rsb__submatrices(mtxAp), 1 );
	rsb_mtx_free(mtxAp);
}


class testRsbUpLoFlagsNonSquareWide : public ::testing::Test {
	protected:
	const rsb_coo_idx_t nrA { 4 }, ncA { 10 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_nnz_idx_t nnzA { nrA * ncA };
	const std::vector<rsb_coo_idx_t> IA { dense_ia(nrA,ncA) };
	const std::vector<rsb_coo_idx_t> JA { dense_ja(nrA,ncA) };
	const std::vector<RSB_DEFAULT_TYPE> VA { dense_va<RSB_DEFAULT_TYPE>(nrA,ncA) };
	const rsb_blk_idx_t brA = 0, bcA = 0;
	rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	RsbLib rsblib { };
};

TEST_F(testRsbUpLoFlagsNonSquareWide, NoLower)
{
	// giving lower flags implies lower triangle presence, not exclusivity
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_LOWER, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( this->nnzA) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( 0 ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlagsNonSquareWide, NonSquareRectangularLowTri)
{
	// both dimensions are given, triangular flags are given, non-lower nonzeroes are removed
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_LOWER | RSB_FLAG_TRIANGULAR, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( 10 ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlagsNonSquareWide, NonSquareSymmetricAutoSquare)
{
	// both dimensions are detected, symmetric triangular flags are given, non-lower nonzeroes are removed
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, 0*nrA, 0*ncA, brA, bcA, RSB_FLAG_LOWER | RSB_FLAG_TRIANGULAR | RSB_FLAG_SYMMETRIC, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( 10 ) );
	EXPECT_THAT( mtxAp->nr, Ne( nrA ) );
	EXPECT_THAT( mtxAp->nc, Eq( ncA) );
	EXPECT_THAT( mtxAp->nr, Eq( mtxAp->nc) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_SYMMETRIC), Eq( RSB_FLAG_SYMMETRIC) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlagsNonSquareWide, NonSquareSymmetricIfNoDimentionsDetection)
{
	// one dimension is detected, symmetric triangular flags are given, non-lower nonzeroes are removed
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, 1*nrA, 0*ncA, brA, bcA, RSB_FLAG_LOWER | RSB_FLAG_TRIANGULAR | RSB_FLAG_SYMMETRIC, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( 10 ) );
	EXPECT_THAT( mtxAp->nr, Eq( this->nrA ) );
	EXPECT_THAT( mtxAp->nc, Eq( this->ncA) );
	EXPECT_THAT( mtxAp->nr, Ne( mtxAp->nc) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_SYMMETRIC), Eq( RSB_FLAG_SYMMETRIC) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlagsNonSquareWide, NoNonSquareSymmetricIfUnspecifiedTriangle)
{
	// dimensions are given, diag flags are detected, non-diag nonzeroes are removed
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, 1*nrA, 1*ncA, brA, bcA, RSB_FLAG_TRIANGULAR | RSB_FLAG_SYMMETRIC, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( 4 ) );
	EXPECT_THAT( mtxAp->nr, Eq( this->nrA ) );
	EXPECT_THAT( mtxAp->nc, Eq( this->ncA) );
	EXPECT_THAT( mtxAp->nr, Ne( mtxAp->nc) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_DIAGONAL), Eq( RSB_FLAG_DIAGONAL) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_SYMMETRIC), Eq( RSB_FLAG_SYMMETRIC) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}


class testRsbUpLoFlags : public ::testing::Test {
	protected:
	const rsb_coo_idx_t nrA { 10 }, ncA { 10 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_nnz_idx_t nnzA { nrA * ncA };
	const rsb_nnz_idx_t halfSquareNnzA { (nrA*(this->ncA+1)/2) };
	const std::vector<rsb_coo_idx_t> IA { dense_ia(nrA,ncA) };
	const std::vector<rsb_coo_idx_t> JA { dense_ja(nrA,ncA) };
	const std::vector<RSB_DEFAULT_TYPE> VA { dense_va<RSB_DEFAULT_TYPE>(nrA,ncA) };
	const rsb_blk_idx_t brA = 0, bcA = 0;
	rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	RsbLib rsblib { };
};

TEST_F(testRsbUpLoFlags,Lower)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_LOWER, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( this->nnzA) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( 0 ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlags,Upper)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_UPPER, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( this->nnzA) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( 0 ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_UPPER), Eq( RSB_FLAG_UPPER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlags,LowerTriangular)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_LOWER_TRIANGULAR, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( this->halfSquareNnzA ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlags,UpperTriangular)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_UPPER_TRIANGULAR, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( this->halfSquareNnzA ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_UPPER), Eq( RSB_FLAG_UPPER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlags,TriangularSymmetricIsDiagonal)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_TRIANGULAR | RSB_FLAG_SYMMETRIC, &errval) };
	EXPECT_THAT( mtxAp->nnz, Eq( this->nrA ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_SYMMETRIC), Eq( RSB_FLAG_SYMMETRIC) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_UPPER), Eq( RSB_FLAG_UPPER ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlags,TriangularUnsymmetricIsDiagonal)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_TRIANGULAR, &errval) };
	EXPECT_THAT( mtxAp->nnz, Lt( nnzA ) );
	EXPECT_THAT( mtxAp->nnz, Eq( this->nrA ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_UPPER), Eq( RSB_FLAG_UPPER ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_DIAGONAL), Eq( RSB_FLAG_DIAGONAL) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlags,UpperLowerTriangularSymmetricIsDiagonal)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_UPPER_TRIANGULAR | RSB_FLAG_LOWER_TRIANGULAR | RSB_FLAG_SYMMETRIC, &errval) };
	EXPECT_THAT( mtxAp->nnz, Lt( nnzA ) );
	EXPECT_THAT( mtxAp->nnz, Eq( this->nrA ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_SYMMETRIC), Eq( RSB_FLAG_SYMMETRIC) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_UPPER), Eq( RSB_FLAG_UPPER ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_DIAGONAL), Eq( RSB_FLAG_DIAGONAL) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}

TEST_F(testRsbUpLoFlags,UpperLowerTriangularIsDiagonal)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_UPPER_TRIANGULAR | RSB_FLAG_LOWER_TRIANGULAR, &errval) };
	EXPECT_THAT( mtxAp->nnz, Lt( nnzA ) );
	EXPECT_THAT( mtxAp->nnz, Eq( this->nrA ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_TRIANGULAR), Eq( RSB_FLAG_TRIANGULAR) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_UPPER), Eq( RSB_FLAG_UPPER ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_LOWER), Eq( RSB_FLAG_LOWER ) );
	EXPECT_THAT( (mtxAp->flags & RSB_FLAG_DIAGONAL), Eq( RSB_FLAG_DIAGONAL) );
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
}


class testCsr: public testRsbAsCxx {
	protected:
	rsb_coo_idx_t ib_ {0};
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS, NULL) };
	struct rsb_coo_mtx_t coo_;
	testCsr(void) { on_error(rsb__project_rsb_to_coo(mtxAp_, &coo_)); }
	~testCsr(void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(testCsr,csr_chk_Ok)
{
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,coo_.nr,coo_.nc,coo_.nnz,ib_), Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(testCsr,csr_chk_BadArgsNnz)
{
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,coo_.nr,coo_.nc,-1,ib_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCsr,csr_chk_BadArgsNr)
{
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,-1,coo_.nc,coo_.nnz,ib_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCsr,csr_chk_BadArgsJA)
{
	EXPECT_THAT( rsb__csr_chk(coo_.IA,NULL,coo_.nr,coo_.nc,coo_.nnz,ib_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCsr,csr_chk_BadArgsPA_Nnz)
{
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,coo_.nr,coo_.nc,coo_.nnz+1,ib_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCsr,csr_chk_BadArgsPA_0)
{
	coo_.IA[0] = coo_.IA[1] + 1;
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,coo_.nr,coo_.nc,coo_.nnz,ib_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCsr,csr_chk_BadArgsPA_1)
{
	coo_.IA[coo_.nr/2]++;
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,coo_.nr,coo_.nc,coo_.nnz,ib_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCsr,csr_chk_BadArgsJA_0)
{
	const rsb_nnz_idx_t nzi = coo_.IA[coo_.nr/2] + coo_.nc/2; // middle row, middle column
	coo_.JA[nzi] = coo_.JA[nzi+1];
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,coo_.nr,coo_.nc,coo_.nnz,ib_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCsr,csr_chk_BadArgsJA_1)
{
	const rsb_nnz_idx_t nzi = coo_.IA[coo_.nr/2] - 1; // pre-middle row, last column
	coo_.JA[nzi] = coo_.nc+ib_;
	EXPECT_THAT( rsb__csr_chk(coo_.IA,coo_.JA,coo_.nr,coo_.nc,coo_.nnz,ib_), Eq( RSB_ERR_BADARGS ) );
}


class testAddCooToDense: public testRsbAsCxx {
	protected:
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS, NULL) };
	struct rsb_coo_mtx_t coo_;
	testAddCooToDense(void) { on_error(rsb__project_rsb_to_coo(mtxAp_, &coo_)); }
	~testAddCooToDense(void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(testAddCooToDense,RowsMajor)
{
	const RSB_DEFAULT_TYPE alpha {1};
	const rsb_nnz_idx_t ldB {ncA};
	const rsb_nnz_idx_t nrB {nrA};
	const rsb_nnz_idx_t ncB {ncA};
	const rsb_bool_t rowmajorB {true};
	std::vector<RSB_DEFAULT_TYPE> B { dense_va<RSB_DEFAULT_TYPE>(nrB,ldB) };
	const rsb_err_t errval = rsb_mtx_add_to_dense(&alpha, mtxAp_, ldB, nrB, ncB, rowmajorB, B.data());
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(testAddCooToDense,ColsMajor)
{
	const RSB_DEFAULT_TYPE alpha {1};
	const rsb_nnz_idx_t ldB {ncA};
	const rsb_nnz_idx_t nrB {nrA};
	const rsb_nnz_idx_t ncB {ncA};
	const rsb_bool_t rowmajorB {false};
	std::vector<RSB_DEFAULT_TYPE> B { dense_va<RSB_DEFAULT_TYPE>(nrB,ldB) };
	const rsb_err_t errval = rsb_mtx_add_to_dense(&alpha, mtxAp_, ldB, nrB, ncB, rowmajorB, B.data());
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}


TEST_F(testAddCooToDense,BadArgsError)
{
	const RSB_DEFAULT_TYPE alpha {1};
	const rsb_nnz_idx_t ldB {ncA};
	const rsb_nnz_idx_t nrB {nrA};
	const rsb_nnz_idx_t ncB {ncA};
	const rsb_bool_t rowmajorB {false};
	std::vector<RSB_DEFAULT_TYPE> B { dense_va<RSB_DEFAULT_TYPE>(nrB,ldB) };
	const rsb_err_t errval = rsb_mtx_add_to_dense(&alpha, nullptr, ldB, nrB, ncB, rowmajorB, B.data());
	EXPECT_THAT( errval, Eq( RSB_ERR_BADARGS ) );
}


class testCoo: public testRsbAsCxx {
	protected:
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS, NULL) };
	struct rsb_coo_mtx_t coo_;
	testCoo(void) { on_error(rsb__project_rsb_to_coo(mtxAp_, &coo_)); }
	~testCoo(void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(testCoo,coo_chk_if_Css)
{
	EXPECT_THAT( mtxAp_, Ne( nullptr ) );
	EXPECT_THAT( rsb__is_css_matrix(mtxAp_ ), Eq( RSB_BOOL_TRUE) );
}

TEST_F(testCoo,coo_chk_err_ok)
{
	const rsb_err_t errval = rsb__mtx_chk(mtxAp_);
	EXPECT_THAT( errval, Eq( RSB_BOOL_TRUE ) );
}

TEST_F(testCoo,coo_chk_err_baddims)
{
	const auto nrA = mtxAp_->nr;
	mtxAp_->nr *= -1;
	const rsb_bool_t errval = rsb__mtx_chk(mtxAp_);
	mtxAp_->nr *= -1;
	EXPECT_THAT( errval, Eq( RSB_BOOL_FALSE ) );
}

TEST_F(testCoo,coo_chk_Valid)
{
	EXPECT_THAT( rsb__util_is_valid_coo_struct(&coo_), Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(testCoo,coo_chk_Invalid)
{
	coo_.nnz = -1;
	EXPECT_THAT( rsb__util_is_valid_coo_struct(&coo_), Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(testCoo,coo_chk_InvalidType)
{
	coo_.typecode = RSB_NUMERICAL_TYPE_INVALID_TYPE;
	EXPECT_THAT( rsb__util_is_valid_coo_struct(&coo_), Eq( RSB_ERR_UNSUPPORTED_TYPE ) );
}

TEST_F(testCoo,coo_chk_InvalidNnz)
{
	coo_.nnz = RSB_MAX_MATRIX_NNZ+1;
	EXPECT_THAT( rsb__util_is_valid_coo_struct(&coo_), Eq( RSB_ERR_BADARGS ) );
}

TEST_F(testCoo,coo_chk_Ok)
{
	EXPECT_THAT( rsb__util_is_valid_coo_array (coo_.IA,coo_.nnz), Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(testCoo,coo_chk_NotOk)
{
	coo_.IA[0]=-1;
	EXPECT_THAT( rsb__util_is_valid_coo_array (coo_.IA,coo_.nnz), Eq( RSB_ERR_GENERIC_ERROR ) );
}

TEST_F(testCoo,SortedRowMajor)
{
	EXPECT_THAT( rsb__util_is_sorted_coo_as_row_major(coo_.IA, coo_.JA, coo_.nnz, coo_.typecode, NULL, RSB_FLAG_WANT_ROW_MAJOR_ORDER), Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(testCoo,SortedRowMajorEmpty)
{
	EXPECT_THAT( rsb__util_is_sorted_coo_as_row_major(coo_.IA, coo_.JA, 0, coo_.typecode, NULL, RSB_FLAG_WANT_ROW_MAJOR_ORDER), Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(testCoo,SortedRowMajorNotColMajor)
{
	EXPECT_THAT( rsb__util_is_sorted_coo_as_row_major(coo_.IA, coo_.JA, coo_.nnz, coo_.typecode, NULL, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER), Eq(RSB_ERR_GENERIC_ERROR) );
}

TEST_F(testCoo,SortedRowMajorIsTransposedColMajor)
{
	EXPECT_THAT( rsb__util_is_sorted_coo_as_row_major(coo_.JA, coo_.IA, coo_.nnz, coo_.typecode, NULL, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER), Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(testCoo,SortedRowMajorBadArgs)
{
	EXPECT_THAT( rsb__util_is_sorted_coo_as_row_major(coo_.IA, coo_.JA, -1, coo_.typecode, NULL, RSB_FLAG_WANT_ROW_MAJOR_ORDER), Eq(RSB_ERR_BADARGS) );
}

TEST_F(testCoo,SortedRowMajorBadType)
{
	EXPECT_THAT( rsb__util_is_sorted_coo_as_row_major(coo_.IA, coo_.JA, coo_.nnz, RSB_NUMERICAL_TYPE_INVALID_TYPE, NULL, RSB_FLAG_WANT_ROW_MAJOR_ORDER), Eq(RSB_ERR_UNSUPPORTED_TYPE) );
}

TEST_F(testCoo,coo_sort)
{
	EXPECT_THAT( rsb_coo_sort(coo_.VA, coo_.IA, coo_.JA, coo_.nr, coo_.nc, coo_.nnz, coo_.typecode, RSB_FLAG_NOFLAGS ), Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(testCoo,coo_sort_zero)
{
	EXPECT_THAT( rsb_coo_sort(coo_.VA, coo_.IA, coo_.JA, coo_.nnz, coo_.nr*0, coo_.nc*0, coo_.typecode, RSB_FLAG_NOFLAGS ), Eq(RSB_ERR_NO_ERROR) );
}


class test__do_util_sortcoo: public ::testing::Test {
	protected:
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	std::vector<RSB_DEFAULT_TYPE> VA {1,1,1,1,1,1,2};
	std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1};
	std::vector<rsb_coo_idx_t> JA {0,1,2,3,4,5,0};
	const rsb_nnz_idx_t nnz { 7 };
	const rsb_coo_idx_t nr {6}, nc {6};
	const struct rsb_mtx_partitioning_info_t * pinfop {};
	const rsb_flags_t flags {RSB_FLAG_WANT_BCSS_STORAGE};
	void * WA {};
	const size_t wb {};
	RsbLib rsblib { };
};

TEST_F(test__do_util_sortcoo,ok_by_rows)
{
	EXPECT_THAT( rsb__do_util_sortcoo(VA.data(), IA.data(), JA.data(), nr, nc, nnz, typecode, pinfop, flags, WA, wb), Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( IA[1], Eq(1) );
	EXPECT_THAT( JA[1], Eq(0) );
	EXPECT_THAT( VA[1],    2  );
	EXPECT_THAT( IA[nnz-1], Eq(5) );
	EXPECT_THAT( JA[nnz-1], Eq(5) );
	EXPECT_THAT( VA[nnz-1],    1  );
}

TEST_F(test__do_util_sortcoo,ok_by_rows_dbg)
{
	const rsb_flags_t flags { this->flags | RSB_FLAG_SHOULD_DEBUG };

	EXPECT_THAT( rsb__do_util_sortcoo(VA.data(), IA.data(), JA.data(), nr, nc, nnz, typecode, pinfop, flags, WA, wb), Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( IA[1], Eq(1) );
	EXPECT_THAT( JA[1], Eq(0) );
	EXPECT_THAT( VA[1],    2  );
	EXPECT_THAT( IA[nnz-1], Eq(5) );
	EXPECT_THAT( JA[nnz-1], Eq(5) );
	EXPECT_THAT( VA[nnz-1],    1  );
}

TEST_F(test__do_util_sortcoo,todo_by_cols)
{
	const rsb_flags_t flags { RSB_FLAG_WANT_BCSS_STORAGE | RSB_FLAG_WANT_COLUMN_MAJOR_ORDER };

	EXPECT_THAT( rsb__do_util_sortcoo(VA.data(), IA.data(), JA.data(), nr, nc, nnz, typecode, pinfop, flags, WA, wb), Eq(RSB_ERR_UNSUPPORTED_FEATURE) );
}

TEST_F(test__do_util_sortcoo, ZeroNNZ)
{
	EXPECT_THAT( rsb__do_util_sortcoo(VA.data(), IA.data(), JA.data(), nr, nc,   0, typecode, pinfop, flags, WA, wb), Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( IA[nnz-1], Eq(1) );
	EXPECT_THAT( JA[nnz-1], Eq(0) );
	EXPECT_THAT( VA[nnz-1],    2  );
}

TEST_F(test__do_util_sortcoo, ipps)
{
	const rsb_flags_t flags {RSB_FLAG_WANT_BCSS_STORAGE | RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT};

	EXPECT_THAT( rsb__do_util_sortcoo(VA.data(), IA.data(), JA.data(), nr, nc, nnz, typecode, pinfop, flags, WA, wb), Eq(RSB_ERR_UNIMPLEMENTED_YET) );
}

TEST_F(test__do_util_sortcoo, BadVA)
{
	EXPECT_THAT( rsb__do_util_sortcoo(     NULL, IA.data(), JA.data(), nr, nc, nnz, typecode, pinfop, flags, WA, wb), Eq(RSB_ERR_BADARGS) );
}


class small_low_tri_for_sort_tests {
	protected:
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	std::vector<RSB_DEFAULT_TYPE> VA {1,1,1,1,1,1,2};
	std::vector<RSB_DEFAULT_TYPE> rVA {VA};
	std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1};
	std::vector<rsb_coo_idx_t> rIA {IA};
	std::vector<rsb_coo_idx_t> JA {0,1,2,3,4,5,0};
	std::vector<rsb_coo_idx_t> rJA {JA};
	const rsb_nnz_idx_t nnz { 7 };
	const rsb_coo_idx_t nr {6}, nc {6};
	const struct rsb_mtx_partitioning_info_t * pinfop {};
	const rsb_flags_t flags {RSB_FLAG_WANT_BCSS_STORAGE};
	void * WA {};
	const size_t wb {};
	const rsb_coo_idx_t br {1}, bc {1};
	const enum rsb_op_flags_t op_flags {RSB_OP_FLAG_DEFAULT};
	RsbLib rsblib { };
};

class test__do_index_based_bcsr_sort: public ::testing::Test, protected small_low_tri_for_sort_tests  {
	/* */
};

TEST_F(test__do_index_based_bcsr_sort,ok_by_rows)
{
	const rsb_err_t errval = rsb__do_index_based_bcsr_sort( IA.data(), JA.data(), VA.data(), rIA.data(), rJA.data(), rVA.data(), nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rIA[1], Eq(1) );
	EXPECT_THAT( rJA[1], Eq(0) );
	EXPECT_THAT( rVA[1],    2  );
	EXPECT_THAT( rIA[nnz-1], Eq(5) );
	EXPECT_THAT( rJA[nnz-1], Eq(5) );
	EXPECT_THAT( rVA[nnz-1],    1  );
}

TEST_F(test__do_index_based_bcsr_sort,one_nnz)
{
	const rsb_nnz_idx_t nnz { 1 };
	const rsb_err_t errval = rsb__do_index_based_bcsr_sort( IA.data(), JA.data(), VA.data(), rIA.data(), rJA.data(), rVA.data(), nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test__do_index_based_bcsr_sort,z_obsoleted)
{
	const rsb_flags_t flags {RSB_FLAG_OBSOLETE_BLOCK_ASYMMETRIC_Z_SORTING};
	const rsb_err_t errval = rsb__do_index_based_bcsr_sort( IA.data(), JA.data(), VA.data(), rIA.data(), rJA.data(), rVA.data(), nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) ); // unsupported
}

TEST_F(test__do_index_based_bcsr_sort,badargs)
{
	const rsb_coo_idx_t nr {0};
	const rsb_err_t errval = rsb__do_index_based_bcsr_sort( IA.data(), JA.data(), VA.data(), rIA.data(), rJA.data(), rVA.data(), nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) ); // unsupported
}

TEST_F(test__do_index_based_bcsr_sort,in_place)
{
	const rsb_flags_t flags {RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT};
	const rsb_err_t errval = rsb__do_index_based_bcsr_sort( IA.data(), JA.data(), VA.data(), rIA.data(), rJA.data(), rVA.data(), nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rIA[1], Eq(1) );
	EXPECT_THAT( rJA[1], Eq(0) );
	EXPECT_THAT( rVA[1],    2  );
	EXPECT_THAT( rIA[nnz-1], Eq(5) );
	EXPECT_THAT( rJA[nnz-1], Eq(5) );
	EXPECT_THAT( rVA[nnz-1],    1  );
}

TEST_F(test__do_index_based_bcsr_sort,two_passes)
{
	const rsb_coo_idx_t nr {2*RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t)};
	const rsb_coo_idx_t nc {2*RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t)};
	const rsb_err_t errval = rsb__do_index_based_bcsr_sort( IA.data(), JA.data(), VA.data(), rIA.data(), rJA.data(), rVA.data(), nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rJA[1], Eq(0) );
	EXPECT_THAT( rVA[1],    2  );
	EXPECT_THAT( rIA[nnz-1], Eq(5) );
	EXPECT_THAT( rJA[nnz-1], Eq(5) );
	EXPECT_THAT( rVA[nnz-1],    1  );
}

TEST_F(test__do_index_based_bcsr_sort,in_place_two_passes)
{
	const rsb_flags_t flags {RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT};
	const rsb_coo_idx_t nr {2*RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t)};
	const rsb_coo_idx_t nc {2*RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t)};
	const rsb_err_t errval = rsb__do_index_based_bcsr_sort( IA.data(), JA.data(), VA.data(), rIA.data(), rJA.data(), rVA.data(), nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
#if 0 /* buggy */
	EXPECT_THAT( rIA[1], Eq(1) );
	EXPECT_THAT( rJA[1], Eq(1) );
	EXPECT_THAT( rVA[1],    1  );
#else /* correct */
	EXPECT_THAT( rJA[1], Eq(0) );
	EXPECT_THAT( rVA[1],    2  );
#endif
	EXPECT_THAT( rIA[nnz-1], Eq(5) );
	EXPECT_THAT( rJA[nnz-1], Eq(5) );
	EXPECT_THAT( rVA[nnz-1],    1  );
}


class test__rsb__util_sort_row_major_bucket_based_parallel: public ::testing::Test, protected small_low_tri_for_sort_tests  {
	/* */
};

TEST_F(test__rsb__util_sort_row_major_bucket_based_parallel,ok)
{
	const rsb_err_t errval = rsb__util_sort_row_major_bucket_based_parallel(rVA.data(), rIA.data(), rJA.data(), nnz, nr, nc, typecode, flags);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test__rsb__util_sort_row_major_bucket_based_parallel,wrongType)
{
	const rsb_err_t errval = rsb__util_sort_row_major_bucket_based_parallel(rVA.data(), rIA.data(), rJA.data(), nnz, nr, nc, RSB_NUMERICAL_TYPE_INVALID_TYPE, flags);
	EXPECT_THAT( errval, Eq(RSB_ERR_UNSUPPORTED_TYPE) );
}

TEST_F(test__rsb__util_sort_row_major_bucket_based_parallel,overflow)
{
	const rsb_err_t errval = rsb__util_sort_row_major_bucket_based_parallel(rVA.data(), rIA.data(), rJA.data(), RSB_MAX_MATRIX_NNZ+1, nr, nc, typecode, flags);
	EXPECT_THAT( errval, Eq(RSB_ERR_LIMITS) );
}


class test__rsb__util_sort_row_major_parallel: public ::testing::Test, protected small_low_tri_for_sort_tests  {
	/* */
};

TEST_F(test__rsb__util_sort_row_major_parallel,ok)
{
	const rsb_err_t errval = rsb__util_sort_row_major_parallel(rVA.data(), rIA.data(), rJA.data(), nnz, nr, nc, typecode, flags);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test__rsb__util_sort_row_major_parallel,wrongType)
{
	const rsb_err_t errval = rsb__util_sort_row_major_parallel(rVA.data(), rIA.data(), rJA.data(), nnz, nr, nc, RSB_NUMERICAL_TYPE_INVALID_TYPE, flags);
	EXPECT_THAT( errval, Eq(RSB_ERR_UNSUPPORTED_TYPE) );
}

TEST_F(test__rsb__util_sort_row_major_parallel,overflow)
{
	const rsb_err_t errval = rsb__util_sort_row_major_parallel(rVA.data(), rIA.data(), rJA.data(), RSB_MAX_MATRIX_NNZ+1, nr, nc, typecode, flags);
	EXPECT_THAT( errval, Eq(RSB_ERR_LIMITS) );
}

TEST_F(test__rsb__util_sort_row_major_parallel,nullargs)
{
	const rsb_nnz_idx_t nnz {RSB__ENOUGH_NNZ_FOR_PARALLEL_SORT};
	const rsb_err_t errval = rsb__util_sort_row_major_parallel(NULL, NULL, NULL, nnz, nr, nc, typecode, flags);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}


class test__file_mtx: public ::testing::Test {
	protected:
	RsbLib rsblib { };
};

TEST_F(test__file_mtx, noTriFlags)
{
	rsb_coo_idx_t nr {}, nc {};
	rsb_coo_idx_t nz {};
	rsb_flags_t flags;
	const rsb_err_t errval = rsb_file_mtx_get_dims(MTXPATH("us.mtx"), &nr, &nc, &nz, &flags);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( nr, 6 );
	EXPECT_THAT( nc, 6 );
	EXPECT_THAT( nz, 21 );
	// no triangular flags returned yet
	EXPECT_THAT( (flags & RSB_FLAG_TRIANGULAR), Eq(RSB_FLAG_NOFLAGS) );
}

TEST_F(test__file_mtx, noTriFlagsUppTri)
{
	rsb_coo_idx_t nr {}, nc {};
	rsb_coo_idx_t nz {};
	rsb_flags_t flags;
	const rsb_err_t errval = rsb_file_mtx_get_dims(MTXPATH("ut.mtx"), &nr, &nc, &nz, &flags);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( nr, 6 );
	EXPECT_THAT( nc, 6 );
	EXPECT_THAT( nz, 21 );
	// no triangular flags returned yet
	EXPECT_THAT( (flags & RSB_FLAG_TRIANGULAR), Eq(RSB_FLAG_NOFLAGS) );
	EXPECT_THAT( (flags & RSB_FLAG_SYMMETRIC), Ne(RSB_FLAG_SYMMETRIC) );
	EXPECT_THAT( (flags & RSB_FLAG_LOWER), Ne(RSB_FLAG_LOWER) );
	EXPECT_THAT( (flags & RSB_FLAG_UPPER), Ne(RSB_FLAG_UPPER) );
}

TEST_F(test__file_mtx, noTriFlagsSym)
{
	rsb_coo_idx_t nr {}, nc {};
	rsb_coo_idx_t nz {};
	rsb_flags_t flags;
	const rsb_err_t errval = rsb_file_mtx_get_dims(MTXPATH("us.mtx"), &nr, &nc, &nz, &flags);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( nr, 6 );
	EXPECT_THAT( nc, 6 );
	EXPECT_THAT( nz, 21 );
	// no triangular flags returned yet
	EXPECT_THAT( (flags & RSB_FLAG_TRIANGULAR), Eq(RSB_FLAG_NOFLAGS) );
	EXPECT_THAT( (flags & RSB_FLAG_SYMMETRIC), Eq(RSB_FLAG_SYMMETRIC) );
	EXPECT_THAT( (flags & RSB_FLAG_LOWER), Ne(RSB_FLAG_LOWER) );
	EXPECT_THAT( (flags & RSB_FLAG_UPPER), Ne(RSB_FLAG_UPPER) );
}

TEST_F(test__file_mtx, noTriFlagsCh)
{
	rsb_coo_idx_t nr {}, nc {};
	rsb_coo_idx_t nz {};
	rsb_flags_t flags;
	const rsb_err_t errval = rsb_file_mtx_get_dims(MTXPATH("ch.mtx"), &nr, &nc, &nz, &flags);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( nr, 6 );
	EXPECT_THAT( nc, 6 );
	EXPECT_THAT( nz, 21 );
	// no triangular flags returned yet
	EXPECT_THAT( (flags & RSB_FLAG_TRIANGULAR), Eq(RSB_FLAG_NOFLAGS) );
	EXPECT_THAT( (flags & RSB_FLAG_HERMITIAN), Eq(RSB_FLAG_HERMITIAN) );
	EXPECT_THAT( (flags & RSB_FLAG_LOWER), Ne(RSB_FLAG_LOWER) );
	EXPECT_THAT( (flags & RSB_FLAG_UPPER), Ne(RSB_FLAG_UPPER) );
}


class test__clone_area: public ::testing::Test {
	protected:
	RsbLib rsblib { };
};

TEST_F(test__clone_area, ok)
{
	const void * p = rsb__clone_area(nullptr, 0);
	EXPECT_THAT( p, Eq(nullptr) );
}


class test__util_coo_alloc: public ::testing::Test {
	protected:
	RsbLib rsblib { };
};

TEST_F(test__util_coo_alloc, bad_type)
{
	const rsb_err_t errval = rsb__util_coo_alloc(nullptr, nullptr, nullptr, 0, RSB_NUMERICAL_TYPE_INVALID_TYPE, RSB_BOOL_TRUE);
	EXPECT_THAT( errval, Eq(RSB_ERR_UNSUPPORTED_TYPE) );
}


class test__clone_guided: public ::testing::Test {
	protected:
	RsbLib rsblib { };
	const struct rsb_mtx_t * mtxAp {nullptr};
	struct rsb_mtx_t src;
	struct rsb_mtx_t dst;
	size_t nmemb {1};
	size_t size {sizeof(src)};
	const rsb_thread_t * cta {nullptr};
	const rsb_thread_t nct {0};
};

TEST_F(test__clone_guided, ok)
{
	rsb_err_t errval {};

	void * p = rsb__clone_area_guided(nullptr, &src, size, nmemb, mtxAp, cta, nct, &errval);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( p, Ne(nullptr) );
	RSB_CONDITIONAL_FREE(p);
}

TEST_F(test__clone_guided, err_zero)
{
	rsb_err_t errval {};
	size_t size {0};
	size_t nmemb {0};

	void * p = rsb__clone_area_guided(nullptr, &src, size, nmemb, mtxAp, cta, nct, & errval);

	EXPECT_THAT( p, Eq(nullptr) );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
	RSB_CONDITIONAL_FREE(p);
}

TEST_F(test__clone_guided, err_nct)
{
	rsb_err_t errval {};
	const rsb_thread_t nct {1};

	void * p = rsb__clone_area_guided(nullptr, &src, size, nmemb, mtxAp, cta, nct, & errval);

	EXPECT_THAT( p, Eq(nullptr) );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
	RSB_CONDITIONAL_FREE(p);
}


class testNnzIdx : public ::testing::Test {
	protected:
};

TEST_F(test__clone_guided, rsb__util_find_coo_max_index_val)
{
	const std::vector<rsb_nnz_idx_t> IA {0,1,2,3,4,5,1};

	EXPECT_THAT( rsb__util_find_coo_max_index_val(IA.data(), IA.size()), Eq(5) );
}


class test_col_exp: public ::testing::Test {
	protected:
};

TEST_F(test_col_exp, zero)
{
	std::vector<rsb_coo_idx_t> JA {0,1,2,3,4,5,1};
	rsb_coo_idx_t K {1} ;
	const rsb_err_t errval = rsb__do_column_expand(JA.data(), JA.size(), &K, 0 );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( K , Eq(0) );
	EXPECT_THAT( JA[JA.size()-1] , Eq(0) );
}

TEST_F(test_col_exp, pone)
{
	std::vector<rsb_coo_idx_t> JA {0,1,2,3,4,5,1};
	rsb_coo_idx_t K {1} ;
	const rsb_err_t errval = rsb__do_column_expand(JA.data(), JA.size(), &K, +1);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( K , Eq(+1) );
	EXPECT_THAT( JA[JA.size()-1] , Eq(+1) );
}

TEST_F(test_col_exp, mone)
{
	std::vector<rsb_coo_idx_t> JA {0,1,2,3,4,5,1};
	rsb_coo_idx_t K {2} ;
	const rsb_err_t errval = rsb__do_column_expand(JA.data(), JA.size(), &K, -2);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( K , Eq(+4) );
	EXPECT_THAT( JA[JA.size()-1] , Eq(0) );
}


class testRsbAsCxxBaseNonConst {
	protected:
	const rsb_coo_idx_t nrA { 16 }, ncA { 16 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_nnz_idx_t nnzA { nrA * ncA };
	std::vector<rsb_coo_idx_t> IA { dense_ia(nrA,ncA) };
	std::vector<rsb_coo_idx_t> JA { dense_ja(nrA,ncA) };
	std::vector<RSB_DEFAULT_TYPE> VA { dense_va<RSB_DEFAULT_TYPE>(nrA,ncA) };
	const rsb_blk_idx_t brA = 0, bcA = 0;
	rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	RsbLib rsblib { };
};

class testGetRows : public testRsbAsCxxBaseNonConst, public ::testing::Test {
protected:
	const rsb_nnz_idx_t rnzA {ncA};
	std::vector<RSB_DEFAULT_TYPE> coo_VA{std::vector<RSB_DEFAULT_TYPE>(rnzA)};
	std::vector<rsb_coo_idx_t> coo_IA{std::vector<rsb_coo_idx_t>(rnzA)};
	std::vector<rsb_coo_idx_t> coo_JA{std::vector<rsb_coo_idx_t>(rnzA)};
	const rsb_flags_t coo_flagsA { RSB_FLAG_NOFLAGS };
	const RSB_DEFAULT_TYPE alpha {2};
	const rsb_coo_idx_t frA {0}, lrA {0};
	rsb_nnz_idx_t rnz {0};
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_NOFLAGS, NULL) };
	const rsb_trans_t transA { RSB_DEFAULT_TRANSPOSITION };
	~testGetRows(void) override { rsb_mtx_free(mtxAp_); }
};

TEST_F(testGetRows, not_null)
{
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(this->transA, &alpha, nullptr, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA, lrA, &rnz, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(testGetRows, get_rnz)
{
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(this->transA, &alpha, this->mtxAp_, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA, lrA, &rnz, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rnz, Eq(this->ncA) );
}

TEST_F(testGetRows, get_coo_wrong)
{
	rsb_nnz_idx_t * rnzp {nullptr};
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(this->transA, &alpha, this->mtxAp_, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA, lrA, rnzp, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(testGetRows, get_coo_wrong_idx)
{
	const rsb_coo_idx_t frA {1};
	const rsb_coo_idx_t lrA {0};
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(this->transA, &alpha, this->mtxAp_, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA, lrA, &rnz, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(testGetRows, get_coo_wrong_frA)
{
	const rsb_coo_idx_t frA {-1};
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(this->transA, &alpha, this->mtxAp_, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA, lrA, &rnz, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(testGetRows, get_coo)
{
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(this->transA, &alpha, this->mtxAp_, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA, lrA, &rnz, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rnz, Eq(this->ncA) );
	EXPECT_THAT( coo_IA[0], Eq(0) );
	EXPECT_THAT( coo_JA[0], Eq(0) );
}

TEST_F(testGetRows, get_coo_trans)
{
	const rsb_trans_t transA { RSB_TRANSPOSITION_C };
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(transA, &alpha, this->mtxAp_, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA, lrA, &rnz, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rnz, Eq(this->ncA) );
	EXPECT_THAT( coo_IA[ncA-1], Eq(ncA-1) );
	EXPECT_THAT( coo_JA[ncA-1], Eq( 0   ) );
}

TEST_F(testGetRows, get_coo_fidx)
{
	const rsb_flags_t coo_flagsA { RSB_FLAG_FORTRAN_INDICES_INTERFACE };
	const rsb_err_t errval = rsb_mtx_get_rows_sparse(this->transA, &alpha, this->mtxAp_, this->coo_VA.data(), this->coo_IA.data(), this->coo_JA.data(), frA+1, lrA+1, &rnz, coo_flagsA);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rnz, Eq(this->ncA) );
	EXPECT_THAT( coo_IA[0], Eq(1) );
	EXPECT_THAT( coo_JA[0], Eq(1) );
}


class testGetCooInPlace : public testRsbAsCxxBaseNonConst, public ::testing::Test {
protected:
	void * coo_VA{};
	rsb_coo_idx_t * coo_IA{}, * coo_JA{};
	const rsb_flags_t coo_flagsA_ { RSB_FLAG_NOFLAGS };
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_inplace(this->VA.data(), this->IA.data(), this->JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_NOFLAGS, NULL) };

	~testGetCooInPlace(void) override { rsb_mtx_free(mtxAp_); }
};

TEST_F(testGetCooInPlace, not_any_matrix)
{
	const rsb_err_t errval = rsb_mtx_switch_to_coo(this->mtxAp_,nullptr,nullptr,nullptr,coo_flagsA_);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(testGetCooInPlace, ok)
{
	const rsb_err_t errval = rsb_mtx_switch_to_coo(this->mtxAp_,&coo_VA,&coo_IA,&coo_JA,coo_flagsA_);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	this->mtxAp_ = NULL; // note: not elegant, but this is the official way
}


#if RSB_WITH_SPARSE_BLAS_INTERFACE
class testSparseBLAS_general : public ::testing::Test {
	RsbLib rsblib { };
};

TEST_F(testSparseBLAS_general, destroy_once)
{
	const rsb_err_t errval = rsb__BLAS_is_type_supported('i');
	EXPECT_THAT( errval, Eq(RSB_ERR_UNSUPPORTED_TYPE) );
}

TEST_F(testSparseBLAS_general, rsb__BLAS_handle_free_just_so)
{
	const auto retval = rsb__BLAS_handle_free(0);

	EXPECT_THAT( retval, Eq(RSB_BLAS_INVALID_VAL) );
}

TEST_F(testSparseBLAS_general, begin_bad)
{
	const auto handle = rsb__BLAS_Xuscr_begin(-1,-1,RSB_NUMERICAL_TYPE_FIRST_BLAS);

	EXPECT_THAT( handle, Eq(RSB_BLAS_INVALID_VAL) );
}

TEST_F(testSparseBLAS_general, begin_bad_rbp)
{
	const rsb_coo_idx_t cbp [] = {1};
	const rsb_coo_idx_t * rbp {nullptr};
	const blas_sparse_matrix A = rsb__BLAS_new_matrix_begin(1, 1, 1, RSB_NUMERICAL_TYPE_FIRST_BLAS, 1, 1, rbp, cbp);

	EXPECT_THAT( A, Eq(RSB_BLAS_INVALID_VAL) );
}

TEST_F(testSparseBLAS_general, begin_ok_xbp)
{
	const rsb_coo_idx_t rbp [] = {1};
	const rsb_coo_idx_t cbp [] = {1};
	const blas_sparse_matrix A = rsb__BLAS_new_matrix_begin(1, 1, 1, RSB_NUMERICAL_TYPE_FIRST_BLAS, 0, 0, rbp, cbp);

	EXPECT_THAT( A, Ne(RSB_BLAS_INVALID_VAL) );
	BLAS_usds(A);
}

TEST_F(testSparseBLAS_general, begin_ok)
{
	const rsb_coo_idx_t * rbp {nullptr}, * cbp {nullptr};
	const blas_sparse_matrix A = rsb__BLAS_new_matrix_begin(1, 1, 1, RSB_NUMERICAL_TYPE_FIRST_BLAS, 0, 0, rbp, cbp);

	EXPECT_THAT( A, Ne(RSB_BLAS_INVALID_VAL) );
	BLAS_usds(A);
}

TEST_F(testSparseBLAS_general, leak_handled_by_librsb)
{
	const rsb_coo_idx_t * rbp {nullptr}, * cbp {nullptr};
	const blas_sparse_matrix A = rsb__BLAS_new_matrix_begin(1, 1, 1, RSB_NUMERICAL_TYPE_FIRST_BLAS, 0, 0, rbp, cbp);

	EXPECT_THAT( A, Ne(RSB_BLAS_INVALID_VAL) );
	// BLAS_usds handled internally at rsb_lib_exit
	rsb_lib_exit( RSB_NULL_EXIT_OPTIONS ); // and remember rsb_lib_exit can be executed many times
}


class testSparseBLAS : public ::testing::Test {
	protected:
	const rsb_blas_int_t nrA { 16 }, ncA { 16 }, nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_blas_int_t nnzA { nrA * ncA };
	const std::vector<rsb_blas_int_t> IA { dense_ia(nrA,ncA) };
	const std::vector<rsb_blas_int_t> JA { dense_ja(nrA,ncA) };
	const std::vector<RSB_DEFAULT_TYPE> VA { dense_va<RSB_DEFAULT_TYPE>(nrA,ncA) };
	const rsb_blas_int_t brA = 0, bcA = 0;
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	RsbLib rsblib { };

	blas_sparse_matrix A { blas_invalid_handle };

	testSparseBLAS (void) {
		// note: we're omitting error checking here (bad)
		A = rsb__BLAS_Xuscr_begin(nrA,ncA,typecode);
		rsb__BLAS_Xuscr_insert_entries(A,nnzA,VA.data(),IA.data(),JA.data());
		BLAS_duscr_end(A);
	}

	~testSparseBLAS (void) override { BLAS_usds(A); }
};

TEST_F(testSparseBLAS, destroy_once)
{
	auto errval = BLAS_usds(A);
	EXPECT_THAT( errval, Eq(RSB_BLAS_NO_ERROR) );
	this->A = blas_invalid_handle;
	errval = BLAS_usds(A);
	EXPECT_THAT( errval, Eq(RSB_BLAS_ERROR) );
}

TEST_F(testSparseBLAS, handles_free)
{
	// pretty intruding: shall restructure
	BLAS_usds(A);
	rsb_err_t errval = rsb__BLAS_handles_free();
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__BLAS_handles_free(); // multiple calls OK
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	this->A = blas_invalid_handle;
}

TEST_F(testSparseBLAS, free_piecewise)
{
	// pretty intruding: shall restructure
	rsb_mtx_t * mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
	EXPECT_THAT( rsb__BLAS_handle_free(this->A), Eq(this->A) );
	this->A = blas_invalid_handle;
	const rsb_err_t errval = rsb__BLAS_handles_free();
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(testSparseBLAS, rsb__BLAS_matrix_retrieve_freed_index)
{
	BLAS_usds(A);
	const blas_sparse_matrix handle = rsb__BLAS_handle_free(A);
	EXPECT_THAT( handle, Eq(RSB_BLAS_INVALID_VAL) );
	this->A = blas_invalid_handle;
}
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */


class testGetCsr : public testRsbAsCxxBase, public ::testing::Test {
protected:
	std::vector<rsb_nnz_idx_t> csr_RP_ {std::vector<rsb_nnz_idx_t>(nrA+1)};
	std::vector<rsb_nnz_idx_t> csr_JA_ {std::vector<rsb_nnz_idx_t>(nnzA)};
	std::vector<double> csr_VA_ {std::vector<double>(nnzA)};
	const rsb_flags_t csr_flagsA { RSB_FLAG_NOFLAGS };
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, RSB_FLAG_NOFLAGS, NULL) };

	~testGetCsr(void) override { rsb_mtx_free(mtxAp_); }
};

TEST_F(testGetCsr, ok)
{
	EXPECT_THAT( rsb__do_get_csr(this->typecode,this->mtxAp_,(rsb_byte_t*)this->csr_VA_.data(),this->csr_RP_.data(),this->csr_JA_.data(),csr_flagsA), Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(testGetCsr, bad_mtx)
{
	EXPECT_THAT( rsb__do_get_csr(this->typecode,nullptr,(rsb_byte_t*)this->csr_VA_.data(),this->csr_RP_.data(),this->csr_JA_.data(),csr_flagsA), Eq(RSB_ERR_BADARGS) );
}

TEST_F(testGetCsr, wrong_type_handling)
{
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_INVALID_TYPE };
	EXPECT_THAT( rsb__do_get_csr(typecode,this->mtxAp_,(rsb_byte_t*)this->csr_VA_.data(),this->csr_RP_.data(),this->csr_JA_.data(),csr_flagsA), Eq(RSB_ERR_UNSUPPORTED_TYPE) );
}

#if defined(RSB_NUMERICAL_TYPE_FLOAT) && (RSB_NUMERICAL_TYPE_FLOAT != RSB_NUMERICAL_TYPE_DEFAULT)
TEST_F(testGetCsr, ok_type_convert)
{
	// one assumes RSB_NUMERICAL_TYPE_FLOAT != typecode for this test to make sense
	const rsb_type_t csr_typecode { RSB_NUMERICAL_TYPE_FLOAT };
	std::vector<float> csr_VA{std::vector<float>(nnzA,-1)};

	EXPECT_THAT( rsb__do_get_csr(csr_typecode,this->mtxAp_,(rsb_byte_t*)csr_VA.data(),this->csr_RP_.data(),this->csr_JA_.data(),csr_flagsA), Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( csr_VA[nnzA-1], Eq(VA[nnzA-1]) );
}
#endif


class test_get_vec : public testRsbAsCxx {
	protected:
	std::vector<RSB_DEFAULT_TYPE> vec { std::vector<RSB_DEFAULT_TYPE>(std::max(this->nrA,this->ncA)) };
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_NOFLAGS, NULL) };
	test_get_vec(void) { }
	~test_get_vec(void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(test_get_vec, ok)
{
	const rsb_err_t errval = rsb_mtx_get_vec(mtxAp_, vec.data(), RSB_EXTF_SUMS_ROW);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test_get_vec, bad_mtx)
{
	const rsb_err_t errval = rsb_mtx_get_vec(nullptr, vec.data(), RSB_EXTF_SUMS_ROW);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test_get_vec, bad_op)
{
	const rsb_err_t errval = rsb_mtx_get_vec(mtxAp_, vec.data(),  RSB_EXTF_NORM_ONE);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test_get_vec, bad_vec)
{
	const rsb_err_t errval = rsb_mtx_get_vec(mtxAp_, nullptr, RSB_EXTF_SUMS_ROW);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}


class test_get_nrm : public testRsbAsCxx {
	protected:
	std::vector<RSB_DEFAULT_TYPE> vec { std::vector<RSB_DEFAULT_TYPE>(std::max(this->nrA,this->ncA)) };
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_NOFLAGS, NULL) };
	test_get_nrm(void) { }
	~test_get_nrm(void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(test_get_nrm, ok)
{
	const rsb_err_t errval = rsb_mtx_get_nrm(mtxAp_, vec.data(), RSB_EXTF_NORM_ONE);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test_get_nrm, bad_mtx)
{
	const rsb_err_t errval = rsb_mtx_get_nrm(nullptr, vec.data(), RSB_EXTF_NORM_ONE);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test_get_nrm, bad_op)
{
	const rsb_err_t errval = rsb_mtx_get_nrm(mtxAp_, vec.data(), RSB_EXTF_SUMS_ROW);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test_get_nrm, bad_vec)
{
	const rsb_err_t errval = rsb_mtx_get_nrm(mtxAp_, nullptr, RSB_EXTF_NORM_ONE);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}


class test__get_coo_element : public testRsbAsCxx {
	protected:
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_NOFLAGS, NULL) };
	test__get_coo_element (void) { }
	~test__get_coo_element (void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(test__get_coo_element, ok)
{
	RSB_DEFAULT_TYPE val { 99 };
	rsb_err_t errval = rsb__do_get_coo_element(mtxAp_ , &val, 0, 0 );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	errval = rsb__do_get_coo_element(mtxAp_, &val, this->nrA-1, this->ncA-1 );
	for (auto i=0;i<this->nrA;++i)
		for (auto j=0;j<this->ncA;++j)
			errval |= rsb__do_get_coo_element(mtxAp_, &val, i, j);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test__get_coo_element, zero)
{
	RSB_DEFAULT_TYPE val { 99 };
	struct rsb_mtx_t mtxA;
	RSB_BZERO(&mtxA,sizeof(mtxA));
	const auto errval = rsb__do_get_coo_element(&mtxA, &val, 0, 0 );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( val , Eq(99) );
}

TEST_F(test__get_coo_element, error_non_css_matrix)
{
	RSB_DEFAULT_TYPE val {};
	struct rsb_mtx_t mtxA;
	rsb_err_t errval { RSB_ERR_NO_ERROR };

	RSB_BZERO(&mtxA,sizeof(mtxA));
	mtxA.nr = mtxA.nc = 2;
	mtxA.br = mtxA.bc = 2;
	errval = rsb__do_get_coo_element(&mtxA, &val, 0, 0 );
	EXPECT_THAT( errval, Eq(RSB_ERR_UNIMPLEMENTED_YET) );
}

TEST_F(test__get_coo_element, error_NULL)
{
	RSB_DEFAULT_TYPE val;
	rsb_err_t errval = rsb__do_get_coo_element(NULL , &val, 1, 1 );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test__get_coo_element, error_bad_index)
{
	RSB_DEFAULT_TYPE val;
	const auto errval = rsb__do_get_coo_element(NULL , &val, this->nrA, this->ncA );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}


struct test_upd_vals : public testRsbAsCxx {
	std::vector<RSB_DEFAULT_TYPE> vec { std::vector<RSB_DEFAULT_TYPE>(std::max(this->nrA,this->ncA)) };
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_NOFLAGS, NULL) };
	const RSB_DEFAULT_TYPE omega { 3 };
	test_upd_vals(void) { }
	~test_upd_vals(void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(test_upd_vals, ok)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	RSB_DEFAULT_TYPE v0 { }, v1 { };

	errval = rsb__do_get_coo_element(mtxAp_, &v0, this->nrA/2, this->ncA/2 );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb_mtx_upd_vals(mtxAp_ , RSB_ELOPF_MUL, &omega );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_get_coo_element(mtxAp_, &v1, this->nrA/2, this->ncA/2 );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	EXPECT_THAT( v1, Ne(v0) );
	EXPECT_THAT( v1, Eq(v0*omega) );
}

TEST_F(test_upd_vals, bad_mtx)
{
	const rsb_err_t errval = rsb_mtx_upd_vals(nullptr, RSB_ELOPF_MUL, &omega );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test_upd_vals, bad_ptr)
{
	const rsb_err_t errval = rsb_mtx_upd_vals(mtxAp_ , RSB_ELOPF_MUL, nullptr );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test_upd_vals, bad_op)
{
	const rsb_err_t errval = rsb_mtx_upd_vals(mtxAp_, static_cast<enum rsb_elopf_t>(-RSB_ELOPF_MUL), &omega );
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}


#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
struct test_upd_vals_binop_complex:
		public testRsbAsCxxBaseT<std::complex<double >,RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX>,
		public ::testing::Test {
	using NT = std::complex<double>;
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(this->VA.data(), this->IA.data(), this->JA.data(), this->nnzA, this->typecode, this->nrA, this->ncA, this->brA, this->bcA, RSB_FLAG_NOFLAGS, NULL) };
	test_upd_vals_binop_complex(void) { }
	~test_upd_vals_binop_complex(void) { rsb_mtx_free(mtxAp_); }
};

TEST_F(test_upd_vals_binop_complex, ok_real)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	std::vector<double> op (nrA); // notice: not NT, but its real counterpart
	std::iota(op.begin(),op.end(),1); // 1,2,..,nrA
	const rsb_coo_idx_t i {this->nrA/2};
	const rsb_coo_idx_t j {2};
	NT v0 { }, v1 { };
	const NT omega_i { op[i] }, omega_j { op[j] };

	errval = rsb__do_get_coo_element(mtxAp_, &v0, i, j );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_elemental_binop(mtxAp_, RSB_ELOPF_SCALE_ROWS_REAL, op.data());
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_get_coo_element(mtxAp_, &v1, i, j );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	EXPECT_THAT( v1, Ne(v0) );
	EXPECT_THAT( v1, Eq(v0*omega_i) );

	errval = rsb__do_elemental_binop(mtxAp_, RSB_ELOPF_SCALE_COLS_REAL, op.data());
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_get_coo_element(mtxAp_, &v1, i, j );
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	EXPECT_THAT( v1, Eq(v0*omega_i*omega_j) );
}

TEST_F(test_upd_vals_binop_complex, error)
{
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	NT op {};
	errval = rsb__do_elemental_binop(mtxAp_, RSB_ELOPF_NEG, &op);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}
#endif


static int foo;
class test__rndr_misc: public testRsbAsCxxBase, public ::testing::Test {
	protected:
	void * pmp {&pmp};
	const char * filename {""};
	const rsb_coo_idx_t pmlWidth {2};
	const rsb_coo_idx_t pmWidth {2};
	const rsb_coo_idx_t pmHeight {2};
	const rsb_marf_t rflags {RSB_MARF_RGB};
};

TEST_F(test__rndr_misc, err_nullptrs)
{
	const char * filename {nullptr};
	const rsb_err_t errval = rsb__do_file_mtx_rndr(pmp, filename, pmlWidth, pmWidth, pmHeight, rflags);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test__rndr_misc, err_file)
{
	const rsb_err_t errval = rsb__do_file_mtx_rndr(pmp, filename, pmlWidth, pmWidth, pmHeight, rflags);
	EXPECT_THAT( errval, Eq(RSB_ERR_GENERIC_ERROR) );
}

TEST_F(test__rndr_misc, err_sz)
{
	const rsb_coo_idx_t pmlWidth {1};
	const rsb_err_t errval = rsb__do_file_mtx_rndr(pmp, filename, pmlWidth, pmWidth, pmHeight, rflags);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test__rndr_misc, err_flags)
{
	const rsb_marf_t rflags {RSB_MARF_EPS};
	const rsb_err_t errval = rsb__do_file_mtx_rndr(pmp, filename, pmlWidth, pmWidth, pmHeight, rflags);
	EXPECT_THAT( errval, Eq(RSB_ERR_UNIMPLEMENTED_YET) );
}

#if __gnu_linux__
TEST_F(test__rndr_misc, write_to_empty_filename)
{
	const rsb_marf_t rflags {RSB_MARF_EPS_L};
	struct rsb_mtx_t * mtxAp_ { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), 1+0*nnzA, typecode, 1+0*nrA, 1+0*ncA, brA, bcA, RSB_FLAG_NOFLAGS, NULL) };
	const rsb_err_t errval = rsb_mtx_rndr("", mtxAp_, pmWidth, pmHeight, rflags);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	rsb_mtx_free(mtxAp_);
}
#endif

TEST_F(test__rndr_misc, err_flags_others)
{
	const rsb_marf_t rflags {RSB_MARF_EPS_L};
	const rsb_err_t errval = rsb__do_file_mtx_rndr(pmp, filename, pmlWidth, pmWidth, pmHeight, rflags);
	EXPECT_THAT( errval, Eq(RSB_ERR_UNIMPLEMENTED_YET) );
}


class test__vec_misc: public ::testing::Test {
	protected:
	RsbLib rsblib { };
};

TEST_F(test__vec_misc, err_but_save_file)
{
	const RSB_DEFAULT_TYPE n {1};
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	// well, I expect the file to open/close, but the -1 length to skip any writing
	// TODO: might want to test this differently
	const rsb_err_t errval = rsb_file_vec_save("/devv/stdout", typecode, &n, -1);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(test__vec_misc, err_load)
{
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_err_t errval = rsb__do_load_vector_file_as_matrix_market(nullptr, typecode, nullptr, nullptr);
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
TEST_F(test__vec_misc, are_similar)
{
	using NT = typename std::complex<float>;
	const std::vector<NT> a { 1.00000, 1.00000, 1.00000 };
	const std::vector<NT> b { 1.00009, 2.00000, 1.00009 };
	const std::vector<NT> c { 1.00011, 2.00000, 1.00011 };
	rsb_err_t errval {};

	errval = rsb__do_are_similar(a.data(), a.data(), 1, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 1);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar(a.data(), b.data(), 1, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 1);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar(a.data(), c.data(), 1, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 1);
	EXPECT_THAT( errval, Ne(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar(a.data(), b.data(), 2, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 2);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar(a.data(), a.data(), 3, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 2, 2);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar(a.data(), b.data(), 2, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 2, 2);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar(a.data(), c.data(), 2, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 2, 2);
	EXPECT_THAT( errval, Ne(RSB_ERR_NO_ERROR) );
}
#endif

#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
TEST_F(test__vec_misc, are_similar_parametric)
{
	using NT = typename std::complex<float>;
	const std::vector<NT> a { 1.00000, 1.00000, 1.00000 };
	const std::vector<NT> b { 1.00001, 1.00001, 1.00000 };
	const std::vector<NT> c { 1.00010, 1.00010, 1.00000 };
	rsb_err_t errval {};

	errval = rsb__do_are_similar_parametric(a.data(), a.data(), 1, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 1, 0);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar_parametric(a.data(), b.data(), 1, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 1, 0);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar_parametric(a.data(), c.data(), 1, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 1, 0);
	EXPECT_THAT( errval, Ne(RSB_ERR_NO_ERROR) );

	errval = rsb__do_are_similar_parametric(a.data(), c.data(), 1, RSB_NUMERICAL_TYPE_FLOAT_COMPLEX, 1, 1, 1);
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}
#endif


class test__err: public ::testing::Test {
	protected:
	RsbLib rsblib { };
};

TEST_F(test__err, perror)
{
	const rsb_err_t errval = rsb__do_perror(stdout, RSB_ERR_CAST(0x80000000));
	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test__err, anycodeok)
{
	rsb_err_t berrval = 0x1;

	for ( size_t s = 1; s < (sizeof(berrval)*RSB_CHAR_BIT) ; ++s, berrval = berrval * 2 )
	{
		const rsb_err_t errval = RSB_ERR_CAST(berrval);
		const char * es = rsb__get_errstr_ptr(errval);

		EXPECT_THAT( es, Ne(nullptr) );
	}
}

TEST_F(test__err, rsb__do_strerror_r_ok)
{
	const rsb_err_t ierrval = RSB_ERR_NO_ERROR ;
	const rsb_err_t errval = rsb__do_strerror_r(ierrval, nullptr, 1);

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
}

TEST_F(test__err, rsb__do_strerror_r_err)
{
	const rsb_err_t ierrval = RSB_ERR_BADARGS;
	const rsb_err_t errval = rsb__do_strerror_r(ierrval, nullptr, 1);

	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}


struct rsb_mtx_t * get_square(rsb_coo_idx_t nrA = 4, rsb_coo_idx_t ncA = 4 )
{
	const rsb_coo_idx_t nrhs { 1 }, incX { 1 }, incY { 1 };
	const rsb_nnz_idx_t nnzA { nrA * ncA };
	const std::vector<rsb_coo_idx_t> IA { dense_ia(nrA,ncA) };
	const std::vector<rsb_coo_idx_t> JA { dense_ja(nrA,ncA) };
	const std::vector<RSB_DEFAULT_TYPE> VA { dense_va<RSB_DEFAULT_TYPE>(nrA,ncA) };
	const rsb_blk_idx_t brA = 0, bcA = 0;
	rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	rsb_flags_t flagsA { RSB_FLAG_NOFLAGS };
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval) };
	return mtxAp;
}

class testSparseSum: public ::testing::Test {
	protected:
	RsbLib rsblib { };
	struct rsb_mtx_t * mtxAp_ { get_square() };
	struct rsb_mtx_t * mtxBp_ { get_square() };

	testSparseSum(void) {
	}

	~testSparseSum(void) override {
		rsb_mtx_free(this->mtxAp_);
		rsb_mtx_free(this->mtxBp_);
	}
};

TEST_F(testSparseSum, rsb_sppsp_bad_mtx_arg)
{
	rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_trans_t transA { RSB_DEFAULT_TRANSPOSITION };
	const void *alphap {};
	const rsb_trans_t transB { RSB_DEFAULT_TRANSPOSITION };
	const void *betap {};
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxCp { rsb_sppsp(typecode, transA, alphap, nullptr, transB, betap, mtxBp_, &errval) };

	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
}

TEST_F(testSparseSum, rsb_sppsp)
{
	rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_trans_t transA { RSB_DEFAULT_TRANSPOSITION };
	const void *alphap {};
	const rsb_trans_t transB { RSB_DEFAULT_TRANSPOSITION };
	const void *betap {};
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxCp { rsb_sppsp(typecode, transA, alphap, mtxAp_, transB, betap, mtxBp_, &errval) };

	EXPECT_THAT( errval, Eq(RSB_ERR_NO_ERROR) );
	EXPECT_THAT( rsb_mtx_free(mtxCp), Eq(nullptr) );
}

class testSparseSumMismatch: public ::testing::Test {
	protected:
	RsbLib rsblib { };
	testSparseSumMismatch(void) {
	}

	~testSparseSumMismatch(void) override {
	}
};

TEST_F(testSparseSumMismatch, rsb_sppsp_bad_dims)
{
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_trans_t transA { RSB_DEFAULT_TRANSPOSITION };
	const void *alphap {};
	const rsb_trans_t transB { RSB_DEFAULT_TRANSPOSITION };
	const void *betap {};
	rsb_err_t errval { RSB_ERR_NO_ERROR };
	struct rsb_mtx_t * mtxAp { get_square(3,4) };
	struct rsb_mtx_t * mtxBp { get_square(4,3) };

	struct rsb_mtx_t * mtxCp{ rsb_sppsp(typecode, transA, alphap, mtxAp, transB, betap, mtxBp, &errval) };
	EXPECT_THAT( errval, Eq(RSB_ERR_BADARGS) );
	EXPECT_THAT( mtxCp, Eq(nullptr) );

	EXPECT_THAT( rsb_mtx_free(mtxAp), Eq(nullptr) );
	EXPECT_THAT( rsb_mtx_free(mtxBp), Eq(nullptr) );
}


class test_filesize: public ::testing::Test {
};

TEST_F(test_filesize, ok) {
	auto sz = rsb__sys_filesize(MTXPATH("pd.mtx.bin"));
	EXPECT_THAT( sz, Gt(0) );
}

TEST_F(test_filesize, bad) {
	auto sz = rsb__sys_filesize(nullptr);
	EXPECT_THAT( sz, Eq(0) );
}


#if RSB_WANT_LONG_IDX_TYPE /* pd.mtx.bin & co have not been produced with long indices */
#if RSB_WANT_EXPERIMENTAL_BINARY_COO
#if RSB_WANT_ZLIB_SUPPORT
#include <zlib.h>
class testBinaryIO: public ::testing::Test {
	protected:
	RsbLib rsblib { };
	FILE *fd;
	FILE *ngzfd { NULL };
	int s;
	int cc;
	rsb_coo_idx_t nr{}, nc{};
	rsb_nnz_idx_t nnz{};

	void init_binary_IO(const char * fn = MTXPATH("pd.mtx.bin"))
	{
		typedef char MM_typecode[4];
		fd = (FILE*)gzopen(fn,"r");
		if ( !fd )
			throw  std::runtime_error("!"); // avoid segfault on missing test file
		MM_typecode matcode;
		rsb_coo_idx_t innz = 0;	/* FIXME */
		s = rsb__mm_read_banner(NULL ,fd,(char(*)[])(&matcode));
		rsb__mm_read_mtx_crd_size(NULL,fd,&nr,&nc,&nnz);
		cc = rsb_getc(NULL,fd);
	}
	testBinaryIO(void)
	{
		init_binary_IO();
	}
	void deinit(void) {
		RSB_FCLOSE(fd);
	}
	~testBinaryIO(void) override {
		deinit();
	}
};

TEST_F(testBinaryIO, ok_mini)
{
	EXPECT_THAT( s , Eq(0) );
	EXPECT_THAT( cc , Eq('B') );
}

TEST_F(testBinaryIO, ok)
{
	std::vector<rsb_coo_idx_t> IA(1+nnz);
	std::vector<rsb_coo_idx_t> JA(1+nnz);
	std::vector<RSB_DEFAULT_TYPE> VA(1+nnz); // need >0 nnz to avoid passing NULL pointers
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_err_t errval = rsb__read_coo_bin_fd(NULL, fd, IA.data(), JA.data(), VA.data(), nr, nc, nnz, typecode);
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}

TEST_F(testBinaryIO, bad_type)
{
	std::vector<rsb_coo_idx_t> IA(1+nnz);
	std::vector<rsb_coo_idx_t> JA(1+nnz);
	std::vector<RSB_DEFAULT_TYPE> VA(1+nnz);
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_INVALID_TYPE };
	const rsb_err_t errval = rsb__read_coo_bin_fd(NULL, fd, IA.data(), JA.data(), VA.data(), nr, nc, nnz, typecode);
	EXPECT_THAT( errval, Eq( RSB_ERR_UNSUPPORTED_TYPE ) );
}

TEST_F(testBinaryIO, bad_file)
{
	deinit();
	init_binary_IO(MTXPATH("pd.mtx.bad.bin"));

	std::vector<rsb_coo_idx_t> IA(nnz);
	std::vector<rsb_coo_idx_t> JA(nnz);
	std::vector<RSB_DEFAULT_TYPE> VA(nnz);
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_INVALID_TYPE };
	const rsb_err_t errval = rsb__read_coo_bin_fd(NULL, fd, IA.data(), JA.data(), VA.data(), nr, nc, nnz, typecode);
	EXPECT_THAT( errval, Eq( RSB_ERR_INTERNAL_ERROR ) );
}

TEST_F(testBinaryIO, ok_one)
{
	deinit();
	init_binary_IO(MTXPATH("pd.mtx.one.bin"));
	// note: this test is not portable.

	std::vector<rsb_coo_idx_t> IA(nnz);
	std::vector<rsb_coo_idx_t> JA(nnz);
	std::vector<RSB_DEFAULT_TYPE> VA(nnz);
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_DEFAULT };
	const rsb_err_t errval = rsb__read_coo_bin_fd(NULL, fd, IA.data(), JA.data(), VA.data(), nr, nc, nnz, typecode);
	EXPECT_THAT( errval, Eq( RSB_ERR_NO_ERROR ) );
}
#endif
#endif /* RSB_WANT_EXPERIMENTAL_BINARY_COO */
#endif /* RSB_WANT_LONG_IDX_TYPE */


class test__timers: public ::testing::Test {
	protected:
	RsbLib rsblib { };
};

TEST_F(test__timers, repeat_small_chunks_hit_maxtimes)
{
	// see rsb__do_bench_spxm
	rsb_time_t ct = RSB_TIME_ZERO, tt = RSB_TIME_ZERO;
	rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO;
	rsb_time_t ss = RSB_TIME_ZERO;
	rsb_int_t times = 0;
	const rsb_time_t jf = 1e-15;
	const rsb_time_t mindt = RSB_TIME_ZERO;
	const rsb_time_t maxdt = RSB_TIME_ZERO;
	const rsb_int_t mintimes = 2;
	const rsb_int_t maxtimes = 3;
	const rsb_time_t it = rsb_time();
	rsb_time_t dt = it;

	do
	{
		/* do nothing, take little time, hit maxtimes */
		RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,jf,times);
	}
	while(RSB_REPEAT(ct-it,times,mindt,mintimes,maxdt,maxtimes));
#if RSB_WANT_REPEAT_TILL_MAXTIMES
	EXPECT_THAT( times, Ge( mintimes ) );
	EXPECT_THAT( times, Eq( maxtimes ) );
#else /* RSB_WANT_REPEAT_TILL_MAXTIMES */
	EXPECT_THAT( times, Eq( mintimes ) );
#endif /* RSB_WANT_REPEAT_TILL_MAXTIMES */
}

TEST_F(test__timers, repeat_small_chunks_hit_maxtimes_minmax)
{
	// see rsb__do_bench_spxm
	rsb_time_t ct = RSB_TIME_ZERO, tt = RSB_TIME_ZERO;
	rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO;
	rsb_time_t ss = RSB_TIME_ZERO;
	rsb_int_t times = 0;
	const rsb_time_t jf = 1e-15;
	const rsb_time_t mindt = 0.0001;
	const rsb_time_t maxdt = 0.001;
	const rsb_int_t mintimes = 2;
	const rsb_int_t maxtimes = 3;
	const rsb_time_t it = rsb_time();
	rsb_time_t dt = it;

	do
	{
		// each iteration takes more than mindt, less than maxdt
		// will hit limit of maxtimes
		while ( rsb_time()-dt <= mindt )
			;
		RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,jf,times);
	}
	while(RSB_REPEAT(ct-it,times,mindt,mintimes,maxdt,maxtimes));

	EXPECT_THAT( tt   , Ge( mindt) );
	EXPECT_THAT( tt   , Lt( maxdt) );
#if RSB_WANT_REPEAT_TILL_MAXTIMES
	EXPECT_THAT( times, Ge( mintimes ) );
	EXPECT_THAT( times, Eq( maxtimes ) );
#else /* RSB_WANT_REPEAT_TILL_MAXTIMES */
	EXPECT_THAT( times, Eq( mintimes ) );
#endif /* RSB_WANT_REPEAT_TILL_MAXTIMES */
}

TEST_F(test__timers, repeat_small_chunks_hit_maxtimes_no_mindt)
{
	// see rsb__do_bench_spxm
	rsb_time_t ct = RSB_TIME_ZERO, tt = RSB_TIME_ZERO;
	rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO;
	rsb_time_t ss = RSB_TIME_ZERO;
	rsb_int_t times = 0;
	const rsb_time_t jf = 1e-15;
	const rsb_time_t mindt = RSB_TIME_ZERO;
	const rsb_time_t maxdt = 0.001;
	const rsb_int_t mintimes = 2;
	const rsb_int_t maxtimes = 3;
	const rsb_time_t it = rsb_time();
	rsb_time_t dt = it;

	do
	{
		// each iteration is very small
		while ( rsb_time()-dt <= jf )
			;
		RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,jf,times);
	}
	while(RSB_REPEAT(ct-it,times,mindt,mintimes,maxdt,maxtimes));

	EXPECT_THAT( tt   , Gt( mindt) );
	EXPECT_THAT( tt   , Lt( maxdt) );
#if RSB_WANT_REPEAT_TILL_MAXTIMES
	EXPECT_THAT( times, Gt( mintimes ) );
	EXPECT_THAT( times, Eq( maxtimes ) );
#else /* RSB_WANT_REPEAT_TILL_MAXTIMES */
	EXPECT_THAT( times, Eq( mintimes ) );
#endif /* RSB_WANT_REPEAT_TILL_MAXTIMES */
}

TEST_F(test__timers, repeat_small_chunks_mintime_hit_maxdt)
{
	// see rsb__do_bench_spxm
	rsb_time_t ct = RSB_TIME_ZERO, tt = RSB_TIME_ZERO;
	rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO;
	rsb_time_t ss = RSB_TIME_ZERO;
	rsb_int_t times = 0;
	const rsb_time_t jf = 1e-15;
	const rsb_time_t mindt = 0.0001;
	const rsb_time_t maxdt = 0.0002;
	const rsb_int_t mintimes = 2;
	const rsb_int_t maxtimes = 3;
	const rsb_time_t it = rsb_time();
	rsb_time_t dt = it;

	do
	{
		// each cycle take more than mindt, with small maxdt
		// hit limit of maxdt
		// won't hit limit of maxtimes
		while ( rsb_time()-dt <= mindt )
			;
		RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,jf,times);
	}
	while(RSB_REPEAT(ct-it,times,mindt,mintimes,maxdt,maxtimes));

	EXPECT_THAT( tt   , Ge( mindt    ) );
	EXPECT_THAT( tt   , Ge( maxdt) );
	EXPECT_THAT( times, Ge( mintimes ) );
	EXPECT_THAT( times, Le( maxtimes ) );
}

TEST_F(test__timers, repeat_big_chunk_hit_maxdt)
{
	// see rsb__do_bench_spxm
	rsb_time_t ct = RSB_TIME_ZERO, tt = RSB_TIME_ZERO;
	rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO;
	rsb_time_t ss = RSB_TIME_ZERO;
	rsb_int_t times = 0;
	const rsb_time_t jf = 1e-15;
	const rsb_time_t mindt = RSB_TIME_ZERO;
	const rsb_time_t maxdt = 0.001;
	const rsb_int_t mintimes = 2;
	const rsb_int_t maxtimes = 3;
	const rsb_time_t it = rsb_time();
	rsb_time_t dt = it;

	do
	{
		while ( rsb_time()-dt < 2 * maxdt )
			;
		// will hit maxdt in first iteration
		RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,jf,times);
	}
	while(RSB_REPEAT(ct-it,times,mindt,mintimes,maxdt,maxtimes));
	EXPECT_THAT( tt   , Ge( maxdt    ) );
	EXPECT_THAT( times, Ge( mintimes ) );
	EXPECT_THAT( times, Le( maxtimes ) );
}

TEST_F(test__timers, repeat_big_chunks_mintime_hit_maxdt)
{
	// see rsb__do_bench_spxm
	rsb_time_t ct = RSB_TIME_ZERO, tt = RSB_TIME_ZERO;
	rsb_time_t bt = RSB_CONST_IMPOSSIBLY_BIG_TIME, wt = RSB_TIME_ZERO;
	rsb_time_t ss = RSB_TIME_ZERO;
	rsb_int_t times = 0;
	const rsb_time_t jf = 1e-15;
	const rsb_time_t mindt = 0.0001;
	const rsb_time_t maxdt = 0.001;
	const rsb_int_t mintimes = 2;
	const rsb_int_t maxtimes = 3;
	const rsb_time_t it = rsb_time();
	rsb_time_t dt = it;

	do
	{
		while ( rsb_time()-dt <= maxdt )
			;
		// guaranteed to surpass maxdt in first iteration, but will need another
		RSB_SAMPLE_STAT(it,ct,dt,tt,bt,wt,ss,jf,times);
	}
	while(RSB_REPEAT(ct-it,times,mindt,mintimes,maxdt,maxtimes));

	EXPECT_THAT( tt   , Ge( mindt) );
	EXPECT_THAT( tt   , Ge( maxdt) );
	EXPECT_THAT( times, Ge( mintimes ) );
	EXPECT_THAT( times, Lt( maxtimes ) );
}

#if defined(LIBRSBPP_H_INCLUDED)
class RsbPP: public ::testing::Test {
	protected:
	RsbLib rsblib { };
	/* deliberately bad combination, to be checked in inner kernel */
	const rsb_flags_t flags {RSB_FLAG_SYMMETRIC};
	const rsb_nnz_idx_t nr {1};
	const rsb_nnz_idx_t nc {2}; // diagonal-aligned symmetric block is non supported
	const rsb_coo_idx_t roff {0};
	const rsb_coo_idx_t coff {0};
	/* these other flags are ok */
	const rsb_nnz_idx_t nnz {0};
	const void* VA {};
	const void* IA {};
	const void* IP {};
	const void* JA {};
	const rsb_coo_idx_t nrhs {1};
	const rsb_coo_idx_t ldX {1};
	const void * rhs {};
	const rsb_coo_idx_t ldY {1};
	void * out {};
	const void * alphap {};
	const rsb_coo_idx_t incx {1};
	const rsb_coo_idx_t incy {1};
	const rsb_trans_t transA {RSB_TRANSPOSITION_N };
	const rsb_flags_t order {RSB_FLAG_WANT_ROW_MAJOR_ORDER};
};

TEST_F(RsbPP, all_types_ok)
{
	for ( rsb_type_t typecode : RSB_MATRIX_TYPE_CODES_ARRAY )
	{
		const rsb_err_t errval = rsbpp_coo_spmm(typecode, flags, nnz, nr, nc, VA, IA, JA, nrhs, ldX, rhs, ldY, out, alphap, incx, incy, transA, roff, coff, order); EXPECT_THAT( errval, Ne( RSB_ERR_UNSUPPORTED_TYPE ) );
		EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
	}

	for ( rsb_type_t typecode : RSB_MATRIX_TYPE_CODES_ARRAY )
	{
		const rsb_err_t errval = rsbpp_csr_spmv(typecode, flags, nnz, nr, nc, VA, IP, JA, rhs, out, alphap, incx, incy, transA, roff, coff);
		EXPECT_THAT( errval, Ne( RSB_ERR_UNSUPPORTED_TYPE ) );
		EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
	}
}

TEST_F(RsbPP, bad_type_not_ok)
{
	const rsb_type_t typecode { RSB_NUMERICAL_TYPE_INVALID_TYPE };
	{
		const rsb_err_t errval = rsbpp_coo_spmm(typecode, flags, nnz, nr, nc, VA, IA, JA, nrhs, ldX, rhs, ldY, out, alphap, incx, incy, transA, roff, coff, order);
		EXPECT_THAT( errval, Eq( RSB_ERR_UNSUPPORTED_TYPE ) );
		EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
	}

	{
		const rsb_err_t errval = rsbpp_csr_spmv(typecode, flags, nnz, nr, nc, VA, IP, JA, rhs, out, alphap, incx, incy, transA, roff, coff);
		EXPECT_THAT( errval, Eq( RSB_ERR_UNSUPPORTED_TYPE ) );
		EXPECT_THAT( errval, Ne( RSB_ERR_NO_ERROR ) );
	}
}

#endif

#endif /* RSB_WANT_GTEST */
