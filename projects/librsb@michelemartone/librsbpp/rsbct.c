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
#if RSBPP_HAS_RSB_H
#include "rsb.h"
#include "librsbpp.h"
#include <assert.h>
#if defined(RSB_NUMERICAL_TYPE_DOUBLE)
static
void test_Coo_SpMV_Symmetric(void)
{
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_type_t typecode = 'D';
	const rsb_nnz_idx_t nnzA = 1;
	const rsb_coo_idx_t nrA = 1;
	const rsb_coo_idx_t ncA = 1;
	const double VA [] = {4};
	const rsb_coo_idx_t IA [] = {0};
	const rsb_coo_idx_t JA [] = {0};
	const double rhs [] = {2};
	double out [] = {3};
	const double alphap [] = {2};
	rsb_coo_idx_t incx = 1;
	rsb_coo_idx_t incy = 1;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	const rsb_coo_idx_t roff = 0;
	const rsb_coo_idx_t coff = 0;

	int ret = rsbpp_coo_spmv(typecode, flags, nnzA, nrA, ncA, VA, IA, JA, rhs, out, alphap, incx, incy, transA, roff, coff);
	assert(ret == 0);
	assert(3 + 2 * 4 * 2 == out[0]);
}
int main(void)
{
	test_Coo_SpMV_Symmetric();
	return 0;
}
#else
int main(void)
{
}
#endif
#else /* RSBPP_HAS_RSB_H */
int main(void)
{
}
#endif /* RSBPP_HAS_RSB_H */
