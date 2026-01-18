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

 @brief C++ example based on <rsb.hpp> of invoking the autotuner on an RsbMatrix matrix with RsbMatrix.tune_spmm() and rendering it with RsbMatrix.rndr().

 \include render.cpp
 */
#include <iostream>
#include <rsb.hpp>
#include <vector>
#include <string>
#include <sstream>
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace ::rsb;

template <typename NT=RSB_DEFAULT_TYPE>
void render_test(std::string filename="../A.mtx")
{
	RsbLib rsblib;
  	RsbMatrix<NT> mtx(filename.c_str());

#if defined(_OPENMP)
  	auto nt = omp_get_max_threads();
#else
	auto nt = 1;
#endif
	
	for ( ; nt > 0 ; nt /= 2 )
	{
		const rsb_trans_t transA { RSB_TRANSPOSITION_N };
  		rsblib.set_num_threads(nt);
		rsb_real_t sf{};
		rsb_int_t maxr = 0;
		rsb_time_t maxt = 0;
		const rsb_coo_idx_t nrhs {1};
		const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
		std::ostringstream oss;
		oss << filename << "-" << nt << "th.eps";
		const auto psfilename = oss.str();

		mtx.tune_spmm(&sf,nullptr,maxr,maxt,transA, nullptr, nrhs, order, nullptr, {}, nullptr, nullptr, {});
		/*! [snip__RsbMatrix_rndr] */
		mtx.rndr(psfilename.c_str());
		/*! [snip__RsbMatrix_rndr] */
	}
}

auto main(const int argc, char * argv[]) -> int
{
	if(argc > 1)
		render_test(argv[1]);
	else
		render_test();
}
