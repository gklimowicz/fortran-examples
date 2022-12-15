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

 @brief C++ example based on <rsb.hpp> for performance-benchmarking a matrix from file using RsbMatrix.spmm(), for various right-hand side counts.

 \include bench.cpp
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
void bench(const std::string filename) {
	const RsbLib rsblib(true);
	const auto mhis { rsblib.get_opt(RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING) };
	rsb_int_t rnt { rsblib.get_num_threads() };

	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	std::cout << "# Librsb initialized with " << rnt << " threads." << std::endl;
	if(mhis.size())
		std::cout << "# detected custom memory hierarchy info: " << mhis << std::endl;

	rsb_time_t dt;
	dt = -rsb_time();
  	RsbMatrix<nt_t> mtx(filename.c_str());
	dt += rsb_time();
       	const rsb_flags_t flagsA {mtx.rsbflags()};
	const bool is_sym = flagsA & ( RSB_FLAG_SYMMETRIC|RSB_FLAG_HERMITIAN ) ? 1 : 0;
	const long long nnz_ops = ( is_sym ? 2 : 1 ) * ( mtx._is_complex() ? 8 : 2 );
	const rsb_real_t mtxocc = mtx._get_storage_bytes();
	const std::string mtxname {filename.begin()+(filename.find_last_of("/\\")?filename.find_last_of("/\\")+1:0),filename.begin()+filename.rfind(".mtx")};

	/*! [snip__RsbMatrix_upd_vals] */
	mtx.upd_vals(RSB_ELOPF_POW,nt_t{0.0}); // set matrix values to ones
	/*! [snip__RsbMatrix_upd_vals] */

	std::cout << "# Read matrix " << std::quoted(mtxname) << " : " << mtx._info() << std::endl << "#" << std::endl;
	std::cout << "# Matrix file " << std::quoted(filename) << " read in " << dt << " s " << std::endl;
	std::cout << "# Matrix occupies " << mtxocc << " bytes " << std::endl;

	for ( rsb_coo_idx_t nrhs : {1,2,4,50,100} ) {
		const rsb_real_t opocc = sizeof(nt_t)*nrhs*(mtx.cols()+mtx.rows());
		std::cout << "#" << std::endl << "# Operands occupy " << opocc << " bytes now ( with " << nrhs << " nrhs )" << std::endl;
	for ( rsb_flags_t order : { RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,
#if defined(RSB_LIBRSB_VER_DATE) && (RSB_LIBRSB_VER_DATE) /* since 1.2.0.10 */
		RSB_FLAG_WANT_ROW_MAJOR_ORDER
#endif
	} )
	{
		const rsb_nnz_idx_t nnzA {mtx.nnz()};
		const char oc = ( order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER ? 'C' : 'R' );
		const rsb_nnz_idx_t ldB {};
		const rsb_nnz_idx_t ldC {};
		const std::vector<nt_t> B(nrhs*mtx.cols(),1.);
		std::vector<nt_t> C(nrhs*mtx.rows(),0.);
		const nt_t alpha { 1.}, beta { 1.};
		const std::string tag { std::string() + mtx.rsbtype() + ':' + std::to_string(nrhs) + '-' + oc + ':' + mtxname };
		const rsb_int_t minits {2}, maxits {5};
		rsb_int_t its {};
		dt = std::numeric_limits<decltype(dt)>::max();
		const rsb_time_t mindt {.3}, t0 {rsb_time()};
		do {
			const rsb_time_t t1 {rsb_time()};
			/*! [snip__spmm] */
			mtx.spmm(transA,&alpha,nrhs,order,B.data(),ldB,&beta,C.data(),ldC);
			/*! [snip__spmm] */
			dt = std::min( dt, rsb_time()-t1 );
			++its;
		} while ( its < minits || ( rsb_time()-t0 < mindt && its < maxits) ); // repeat at least minits times, at least mindt s and maximally maxits times
		const auto flops_u = (nnz_ops*nnzA*nrhs) / dt;
		std::cout << "# rsb_spmm-" << nrhs << "-" << oc << " took " << dt << " s, for " << nnzA/dt << " nnz/s, " << flops_u << " flops/s on " << its << " samples\n";
		std::cout << tag << "\t" << flops_u << std::endl;
	}
	}
	std::cout << "# Done." << std::endl;
}

auto main(const int argc, char * argv[]) -> int
{
	const std::string filename{ argc > 1 ? argv[1] : "../A.mtx"};

	std::cout << "MTX"<< "\t" << "FLOPS" << std::endl;
#ifdef RSB_NUMERICAL_TYPE_FLOAT
	bench<float>(filename);
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	bench<double>(filename);
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	bench<std::complex<float>>(filename);
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	bench<std::complex<double>>(filename);
#endif
}
