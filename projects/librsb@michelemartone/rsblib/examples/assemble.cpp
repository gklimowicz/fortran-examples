/*

Copyright (C) 2021-2021 Michele Martone

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

 @brief C++ example based on <rsb.hpp> assembling RsbMatrix by pieces.

 \include assemble.cpp
 */
#include <rsb.hpp>
#include <vector>
#include <cassert>
#include <array>

using namespace rsb;

auto main() -> int {
  RsbLib rsblib;
  /*! [snip__RsbMatrix_Assemble] */
  const rsb_coo_idx_t nrA { 4 }, ncA { 4 };
  RsbMatrix<double> mtx(nrA,ncA); // begin matrix assembly

  // insert elements of a tridiagonal matrix, one by one
  for (auto i = 0; i < nrA; ++i )
    for (auto j = i-1; j <= i+1; ++j )
      if ( i >= 0 && i < nrA )
        if ( j >= 0 && j < ncA )
          mtx.set_val((i+1)*100+(j+1),i,j); // add entry

  mtx.close(); // finish matrix assembly
  assert(mtx.nnz() == 3 * nrA - 2);
  /*! [snip__RsbMatrix_Assemble] */

  mtx.file_save(nullptr); // print out
}
