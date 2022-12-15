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
 \file
 @author Michele Martone

 @brief C++ example based on <rsb.hpp> using RsbMatrix.spmm().

 Using a \librsb  program via \rsblib does not differ conceptually much \librsb.

 Errors caught by \librsb shall not go unnoticed and trigger an exception instead.

 Memory management of matrices and the library state itself follow the usual C++ RAII rules:
 the \c mtx object is freed first via RsbMatrix's destructor;
 then \librsb is finalized via RsbLib()'s destructor .

 \include example.cpp
 */
#include <rsb.hpp>
using namespace rsb;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE

#if RSBP_WANT_CPP20
#include <vector>
#include <array>

/* If your compiler is C++20 compatible, the std::span-based interface is available like: */
auto main() -> int {
  RsbLib rsblib;
  const rsb_nnz_idx_t nnzA { 7 };
  const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 2 };
  const rsb_coo_idx_t IA [nnzA] {0,1,2,3,4,5,1};
  std::array<rsb_coo_idx_t,7> JA {0,1,2,3,4,5,0};
  std::vector<double> VA {1,1,1,1,1,1,2}, X(nrhs*ncA,1);
  std::vector<double> Y(nrhs*nrA,0);
  const double alpha {2}, beta {1};
  rsb_int_t tn {0};
  rsb_real_t sf {0}; // speedup factor (tune_spmm output)
  const rsb_flags_t order {RSB_FLAG_WANT_COLUMN_MAJOR_ORDER};

  // IA,JA,VA are respectively from a C array, std::vector, std::array;
  // here using the C++20's std::span interface:
  RsbMatrix<double> mtx(IA,JA,VA);

  mtx.file_save(); // rsb_file_mtx_save
  mtx.tune_spmm(sf,RSB_TRANSPOSITION_N,alpha,nrhs,order,X,beta,Y); // rsb_tune_spmm
  mtx.spmm(RSB_TRANSPOSITION_N,alpha,nrhs,order,X,beta,Y); // rsb_spmv
}
#else
#include <vector>
/* The pointer-based interface is available as well: */
auto main() -> int {
  RsbLib rsblib;
  const rsb_nnz_idx_t nnzA { 7 };
  const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 };
  const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
  const std::vector<double> VA {1,1,1,1,1,1,2}, X(ncA,1);
  std::vector<double> Y(nrA,0);
  const double alpha {2}, beta {1};
  rsb_int_t tn {0};
  rsb_real_t sf {0}; // speedup factor (tune_spmm output)

  RsbMatrix<double> mtx(IA.data(),JA.data(),VA.data(),nnzA);

  mtx.file_save(nullptr); // rsb_file_mtx_save
  mtx.tune_spmm(&sf,&tn,0,0.0,RSB_TRANSPOSITION_N,alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,X.data(),ncA,beta,Y.data(),nrA); // rsb_tune_spmm
  mtx.spmv(RSB_TRANSPOSITION_N, alpha, X.data(), beta, Y.data()); // rsb_spmv
}
#endif
#else
auto main() { }
#endif
