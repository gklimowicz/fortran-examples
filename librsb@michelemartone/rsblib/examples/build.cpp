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

 @brief C++ example based on <rsb.hpp> using timings for common matrix operations on RsbMatrix: RsbMatrix.get_coo(), rsb_coo_sort(), rsb_time().

 \include build.cpp
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
	RsbLib rsblib;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_int_t nt { rsblib.get_num_threads() };
	const RsbMatrix<nt_t> mtx(filename.c_str());
	std::vector<nt_t> VA(mtx.nnz());
	std::vector<rsb_coo_idx_t> IA(mtx.nnz());
	std::vector<rsb_coo_idx_t> JA(mtx.nnz());
	rsb_time_t dt;


	dt = -rsb_time();
	const RsbMatrix<nt_t> mtx2{ mtx };
	dt += rsb_time();

	std::cout << "copy-" << nt << " took " << dt << std::endl;


	mtx.get_coo(transA,VA.data(),IA.data(),JA.data(),RSB_FLAG_NOFLAGS);
	dt += rsb_time();
	std::cout << "get_coo-" << nt << " took " << dt << std::endl;

	dt = -rsb_time();
	rsb_coo_sort(VA.data(), JA.data(), IA.data(), VA.size(), mtx.cols(), mtx.rows(),  mtx.rsbtype(), RSB_FLAG_NOFLAGS);
	dt += rsb_time();
	std::cout << "rsb_coo_sort-T-" << nt << " took " << dt << std::endl;

	dt = -rsb_time();
	rsb_coo_sort(VA.data(), IA.data(), JA.data(), VA.size(), mtx.rows(), mtx.cols(),  mtx.rsbtype(), RSB_FLAG_NOFLAGS);
	dt += rsb_time();
	std::cout << "rsb_coo_sort-N-" << nt << " took " << dt << std::endl;
}

auto main(const int argc, char * argv[]) -> int
{
	const std::string filename{ argc > 1 ? argv[1] : "../A.mtx"};

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

