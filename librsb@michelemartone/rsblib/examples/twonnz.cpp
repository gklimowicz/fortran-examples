/*

Copyright (C) 2020-2021 Michele Martone

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
/*!
 \ingroup rsb_doc_examples
 @file
 @author Michele Martone

 @brief C++ example based on <rsb.hpp> measuring RsbMatrix.spmm() performance of a matrix with only two elements; this is is effectively measuring performance of result vector scaling.

 \include twonnz.cpp
 */
#include <complex>
#include <array>
#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include <rsb.hpp>

using namespace ::rsb;

template <typename nt_t>
void bench(const rsb_coo_idx_t n, const rsb_coo_idx_t nrhs, const nt_t alpha, const nt_t beta) {
  	RsbLib rsblib;
	const rsb_int_t rnt { rsblib.get_num_threads() };
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	rsb_time_t dt;
	dt = -rsb_time();
	const std::vector<rsb_coo_idx_t> IA{0,n-1};
	const std::vector<rsb_coo_idx_t> JA{0,n-1};
	const std::vector<nt_t> VA{1,1};
  	const RsbMatrix<nt_t> mtx(IA.data(),JA.data(),VA.data(),VA.size());
	dt += rsb_time();
	const rsb_real_t mtxocc = mtx._get_storage_bytes();
	const rsb_real_t opocc = sizeof(nt_t)*nrhs*(mtx.cols()+mtx.rows());
	const rsb_nnz_idx_t nnzA {mtx.nnz()};
	const nt_t pone = 1;

	/*! [snip__RsbMatrix_dims] */
	std::cout << "# Matrix sized " << mtx.rows() << "x" << mtx.cols() << ", " << nnzA << " nnz  built in " << dt << " s and occupies " << mtxocc << " bytes " << std::endl;
	/*! [snip__RsbMatrix_dims] */

	for ( rsb_flags_t order : { RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,
#if defined(RSB_LIBRSB_VER_DATE) && (RSB_LIBRSB_VER_DATE) /* since 1.2.0.10 */
		RSB_FLAG_WANT_ROW_MAJOR_ORDER
#endif
	} )
	{
		std::vector<rsb_time_t> dta;
		const long long nnz_ops = ( mtx._is_complex() ? 6 : 1 );
		const auto flops_c = ( beta != pone ? nnz_ops*n*nrhs : 0LL ); // note we ignore the matrix
		const char oc = ( order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER ? 'C' : 'R' );

		for ( const auto nt : {1,rnt} )
		{
			const rsb_nnz_idx_t ldB {};
			const rsb_nnz_idx_t ldC {};
			const std::vector<nt_t> B(nrhs*mtx.cols(),1.);
			std::vector<nt_t> C(nrhs*mtx.rows(),0.);

			rsblib.set_num_threads(nt);
			mtx.spmm(transA,alpha,nrhs,order,B.data(),ldB,beta,C.data(),ldC); // caches warmup

			dt = -rsb_time();
			mtx.spmm(transA,alpha,nrhs,order,B.data(),ldB,beta,C.data(),ldC);
			dt += rsb_time();
			dta.push_back(dt);
		}
		/*! [snip__RsbMatrix_rsbtype] */
		std::cout << "# type=" << mtx.rsbtype() << " nt=1," << rnt << " n=" << n << " nrhs=" << nrhs << " order=" << oc << " alpha=" << alpha << " beta=" << beta << " dt=" << dta[0] << ".." << dta[1] << " spmm-scalability=" << dta[0]/dta[1] << " nnz/s=" << nnzA/dta[0] << ".." << nnzA/dta[1] << " flops=" << flops_c/dta[0] << ".." << flops_c/dta[1] << " occ.=" << opocc << " " << std::endl;
		/*! [snip__RsbMatrix_rsbtype] */
	}
}

auto main(const int argc, char * argv[]) -> int
{
	const rsb_coo_idx_t n { argc > 1 ? std::stoi(argv[1]) : 1000000};
	const rsb_coo_idx_t nrhs { argc > 2 ? std::stoi(argv[2]) : 1};
	const int alpha { argc > 3 ? std::stoi(argv[3]) : 1};
	const int beta { argc > 4 ? std::stoi(argv[4]) : 1};
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	bench<std::complex<double>>(n,nrhs,alpha,beta);
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	bench<double>(n,nrhs,alpha,beta);
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	bench<std::complex<float>>(n,nrhs,alpha,beta);
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
	bench<float>(n,nrhs,alpha,beta);
#endif
}
