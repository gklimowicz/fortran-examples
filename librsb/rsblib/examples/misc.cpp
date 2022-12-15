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
/*!
 \ingroup rsb_doc_examples
 \file
 @author Michele Martone

 @brief C++ example based on <rsb.hpp> showing various RsbMatrix operations.

 \include misc.cpp
 */
#include <rsb.hpp>
#include <vector>
#include <cassert>
using namespace rsb;
auto main() -> int {
  RsbLib rsblib;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
  const rsb_nnz_idx_t nnzA { 7 };
  const rsb_coo_idx_t nrA { 6 }, ncA { 6 };
  const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
  const std::vector<double> VA {1,1,1,1,1,1,2}, X(ncA,1);
  std::vector<double> Y(nrA,0);

  RsbMatrix<double> mtx1(IA.data(),JA.data(),VA.data(),nnzA);
  RsbMatrix<double> mtx2(IA.data(),JA.data(),VA.data(),nnzA);

  // matrices can be compared:
  /*! [snip__cmp_RsbMatrix] */
  assert(   mtx1 == mtx2  );
  assert( !(mtx1 != mtx2) );
  /*! [snip__cmp_RsbMatrix] */

  // matrices internals can be std::move'd:
  /*! [snip__move_RsbMatrix] */
  assert( mtx1.nnz() == nnzA );
  RsbMatrix<double> mtx3 { std::move(mtx1) };
  assert( mtx3.nnz() == nnzA );
  /*! [snip__move_RsbMatrix] */


  {
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif
    /*! [snip__RsbMatrix_asm_coo] */
    const rsb_nnz_idx_t nnzA { 7 };
    const rsb_coo_idx_t nrA { 6 }, ncA { 6 }, nrhs { 1 };
    const std::vector<rsb_coo_idx_t> IA {0,1,2,3,4,5,1}, JA {0,1,2,3,4,5,0};
    const std::vector<double> VA {1,1,1,1,1,1,2}, X(ncA,1);
    RsbMatrix<double> mtx(IA.data(),JA.data(),VA.data(),nnzA);
    /*! [snip__RsbMatrix_asm_coo] */

    /*! [snip__RsbMatrix_file_save_stdout] */
    mtx.file_save(); // print to stdout
    /*! [snip__RsbMatrix_file_save_stdout] */
  }
#endif
 }
