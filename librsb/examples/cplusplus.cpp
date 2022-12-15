/*

Copyright (C) 2008-2021 Michele Martone

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

 @brief C++ program using the C <rsb.h> interface.
  Uses #rsb_tune_spmm(), #rsb_spmv().

 \see See examples/span.cpp for a C++ example based on the native C++ interface of <rsb.hpp>.

 \include cplusplus.cpp
*/
#include <rsb.h>	/* librsb header to include */
#include <vector>
#include <stdexcept>

int main(void)
{
	/* You can use the C interface of rsb.h in C++ as well. */
	const rsb_nnz_idx_t nnzA { 7 };
	const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 };
	const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
	const std::vector<RSB_DEFAULT_TYPE> VA {1,1,1,1,1,1,2}, X(ncA,1);
	std::vector<RSB_DEFAULT_TYPE> Y(nrA,0);
	const RSB_DEFAULT_TYPE alpha {2}, beta {1};
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT; // see rsb_types.h
	rsb_int_t tn {0};
	rsb_real_t sf {0}; // speedup factor
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
		throw  std::runtime_error("failure running rsb_lib_init!");

	auto mtxAp = rsb_mtx_alloc_from_coo_const(
		VA.data(),IA.data(),JA.data(),nnzA,typecode,nrA,ncA,1,1,
		RSB_FLAG_NOFLAGS    /* default format will be chosen */
		|RSB_FLAG_DUPLICATES_SUM/* duplicates will be summed */
			,&errval);
	if ( ! mtxAp || errval != RSB_ERR_NO_ERROR)
		throw  std::runtime_error("failure running rsb_mtx_alloc_from_coo_const!");

	errval = rsb_file_mtx_save(mtxAp,nullptr); // rsb_file_mtx_save
	if ( ! mtxAp )
		throw  std::runtime_error("failure running rsb_file_mtx_save!");

	errval = rsb_tune_spmm(&mtxAp,&sf,&tn,0,0.0,RSB_TRANSPOSITION_N,&alpha,nullptr,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X.data(),ncA,&beta,Y.data(),nrA);
	if ( ! mtxAp )
		throw  std::runtime_error("failure running rsb_tune_spmm!");

	errval = rsb_spmv(RSB_TRANSPOSITION_N, &alpha, mtxAp, X.data(), 1, &beta, Y.data(), 1); // rsb_spmv
	if ( ! mtxAp )
		throw  std::runtime_error("failure running rsb_spmv!");

	rsb_mtx_free(mtxAp);

	if((errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS)) != RSB_ERR_NO_ERROR)
		throw  std::runtime_error("failure running rsb_lib_exit!");
}
