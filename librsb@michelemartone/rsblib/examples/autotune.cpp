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
/*!
 \ingroup rsb_doc_examples
 @file
 @author Michele Martone

 @brief C++ example based on <rsb.hpp> for performance-benchmarking a matrix from file using RsbMatrix.tune_spmm(), for various right-hand side counts.

 \include autotune.cpp
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
void bench(const std::string filename, rsb_flags_t order) {
	RsbLib rsblib;
	rsb_int_t rnt { rsblib.get_num_threads() };
	const rsb_int_t one {1};

	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	std::cout << "Librsb initialized with " << rnt << " threads." << std::endl;

	RsbMatrix<nt_t> mtx(filename.c_str());
	const rsb_flags_t flagsA {mtx.rsbflags()};
	const bool is_sym = flagsA & ( RSB_FLAG_SYMMETRIC|RSB_FLAG_HERMITIAN ) ? 1 : 0;
	const long long nnz_ops = ( is_sym ? 2 : 1 ) * ( mtx._is_complex() ? 8 : 2 );
	const rsb_real_t mtxocc = mtx._get_storage_bytes();

	std::cout << "Read matrix " << std::quoted(filename) << " : " << mtx._info() << std::endl << std::endl;
	std::cout << "Matrix occupies " << mtxocc << " bytes " << std::endl;

	for ( rsb_coo_idx_t nrhs : {1,2,4,50,100} )
	{
		const rsb_real_t opocc = sizeof(nt_t)*nrhs*(mtx.cols()+mtx.rows());
		std::cout << std::endl << "Operands occupy " << opocc << " bytes now ( with " << nrhs << " nrhs )" << std::endl;
		const rsb_nnz_idx_t nnzA {mtx.nnz()};
		const char oc = ( order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER ? 'C' : 'R' );
		const rsb_nnz_idx_t ldB {};
		const rsb_nnz_idx_t ldC {};
		const std::vector<nt_t> B(nrhs*mtx.cols(),1.);
		std::vector<nt_t> C(nrhs*mtx.rows(),0.);

		const nt_t alpha { 1.}, beta { 1.};
		rsb_time_t dt, tt;
		const rsb_blk_idx_t nsmA{mtx.blocks()};

		mtx.spmm(transA,&alpha,nrhs,order,B.data(),ldB,&beta,C.data(),ldC); // caches warmup

		dt = -rsb_time();
		mtx.spmm(transA,&alpha,nrhs,order,B.data(),ldB,&beta,C.data(),ldC);

		dt += rsb_time();

		const auto flops_u = (nnz_ops*nnzA*nrhs)/dt;
		std::cout << "rsb_spmm-" << nrhs << "-" << oc << " took " << dt << " s, for " << nnzA/dt << " nnz/s, " << flops_u << " flops/s\n";

		if(true)
		{
			rsb_real_t sf = 1.0;
			const rsb_int_t maxr = 10;
			const rsb_time_t tmax = 20;
			rsb_int_t tn {0};

			// std::cout << "Turning on autotuning before each nrhs count..." << std::endl;

			if(false)
				rsblib.set_opt(RSB_IO_WANT_VERBOSE_TUNING, &one);

			/*! [snip__tune_spmm] */
			tt = -rsb_time();
			mtx.tune_spmm(&sf,&tn,maxr,tmax,transA,&alpha,nrhs,order,B.data(),ldB,&beta,C.data(),ldC);
			tt += rsb_time();

			auto nnsmA {mtx.blocks()};
			std::cout << "Tuning took " << tt << " s ( " << tt / dt << " ops ) and changed " << nsmA << " to " << nnsmA << " blocks" << std::endl;

			mtx.spmm(transA,&alpha,nrhs,order,B.data(),ldB,&beta,C.data(),ldC); // caches warmup
			/*! [snip__tune_spmm] */
			dt = -rsb_time();
			mtx.spmm(transA,&alpha,nrhs,order,B.data(),ldB,&beta,C.data(),ldC);
			dt += rsb_time();

			const auto flops_o = (nnz_ops*nnzA*nrhs)/dt;
			std::cout << "rsb_spmm-" << nrhs << "-" << oc << " took " << dt << " s, for " << nnzA/dt << " nnz/s, " << flops_o << " flops/s\n";
			if ( sf > 1.0 && flops_o > flops_u  )
				std::cout << "Tuning brought a " << std::max(flops_o / flops_u,sf) << " x speedup" << std::endl;
			else
				std::cout << "Tuning brought no speedup" << std::endl;
		}
	}

	std::cout << "Done." << std::endl;
}

auto main(const int argc, char * argv[]) -> int
{
	const std::string filename{ argc > 1 ? argv[1] : "../A.mtx"};

	for ( rsb_flags_t order : { RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,
#if defined(RSB_LIBRSB_VER_DATE) && (RSB_LIBRSB_VER_DATE) /* since 1.2.0.10 */
		RSB_FLAG_WANT_ROW_MAJOR_ORDER
#endif
	} )
	{
#ifdef RSB_NUMERICAL_TYPE_FLOAT
		bench<float>(filename,order);
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
		bench<double>(filename,order);
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
		bench<std::complex<float>>(filename,order);
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
		bench<std::complex<double>>(filename,order);
#endif
	}
}
