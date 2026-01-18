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

 @brief C++ example based on <rsb.hpp> illustrating use of RsbMatrix.file_save(), and std::span-based versions of RsbMatrix.tune_spmm(), RsbMatrix.spmv().

   Using a \librsb  program via \rsblib does not differ conceptually much \librsb.

   Errors caught by \librsb shall not go unnoticed and trigger an exception instead.

   Memory management of matrices and the library state itself follow the usual C++ RAII rules:
   the \c mtx object is freed first via RsbMatrix's destructor;
   then \librsb is finalized via RsbLib()'s destructor .

 \include span.cpp
 */
#include <rsb.hpp>
#include <vector>
#include <array>
using namespace ::rsb;
#if defined(RSBP_WANT_CPP20) && defined(RSB_NUMERICAL_TYPE_DOUBLE)
auto main() -> int {
  /*! [snip__span_RsbMatrix] */
  RsbLib rsblib;
  const rsb_nnz_idx_t nnzA { 7 };
  const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 };
  const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1};
  const rsb_coo_idx_t JA [] = {0,1,2,3,4,5,0};
  const std::vector<double> VA {1,1,1,1,1,1,2}, X(ncA,1);
  std::array<double,nrA> Y;
  const double alpha {2}, beta {1};
  rsb_int_t tn {0};

  RsbMatrix<double> mtx(IA,JA,VA,nnzA);
  /*! [snip__span_RsbMatrix] */

  mtx.file_save(nullptr);
  /*! [snip__RsbMatrix_spmv] */
  mtx.tune_spmm(nullptr,&tn,0,0.0,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X,ncA,beta,Y,nrA);
  mtx.spmv(RSB_TRANSPOSITION_N, alpha, X, beta, Y);
  /*! [snip__RsbMatrix_spmv] */
}
#else /* RSBP_WANT_CPP20 */
int main () {}
#endif /* RSBP_WANT_CPP20 */
