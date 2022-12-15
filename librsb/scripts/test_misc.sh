#!/bin/bash
#
# Copyright (C) 2008-2020 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

set -e
set -x


# bug: symmetric upper triangle ignored by rsb_blas_file_mtx_load()
cat > crash.mtx << EOF
%%MatrixMarket matrix coordinate real symmetric
2 2 1
1 2 1
EOF
cat > tp.c << EOF
#include "rsb.h"
#include "blas_sparse.h"
  blas_sparse_matrix rsb__load_spblas_matrix_file_as_matrix_market(const rsb_char_t * filename, rsb_type_t typecode );
int main()
{
  rsb_lib_init(NULL);
  rsb_file_mtx_save(rsb_blas_get_mtx(rsb__load_spblas_matrix_file_as_matrix_market("crash.mtx",RSB_NUMERICAL_TYPE_DEFAULT)),NULL);
}
EOF
`librsb-config --cc` `librsb-config --I_opts` -o tp tp.c `librsb-config --static --ldflags --extra_libs`
ldd ./tp
./tp | grep ^2.2.1$

# alternative
#test -x examples/fortran
#cat > pd.mtx << EOF
#%%MatrixMarket matrix coordinate real symmetric
#2 2 1
#1 2 1
#EOF
#examples/fortran | grep Read.matrix.*1
#true



