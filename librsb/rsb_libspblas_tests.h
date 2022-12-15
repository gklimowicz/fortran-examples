/*                                                                                                                            

Copyright (C) 2008-2020 Michele Martone

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
/* @cond INNERDOC */
/**
 * @file
 * @author Michele Martone
 * @brief  Sparse BLAS interface testing code
 * */
#ifndef LIBSPBLAS_TESTS_H_INCLUDED
#define LIBSPBLAS_TESTS_H_INCLUDED
#include "rsb_common.h"
struct rsb_tester_options_t{
	rsb_time_t mtt; /* maximal test time */
	rsb_bool_t rrm; /* require recursive matrices (error otherwise) */
	rsb_bool_t tur; /* test until recursive */
	rsb_bool_t wqt; /* want quiet testing */
	rsb_bool_t wqc; /* want quiet conditionally (on no tty) */
	rsb_bool_t wcs; /* want clear screen */
	rsb_bool_t wvr; /* want verbose reporting */
	rsb_bool_t wnl; /* want no limits testing */
};
rsb_err_t rsb_blas_tester_options_init(struct rsb_tester_options_t * top);
rsb_err_t rsb_blas_mini_tester(const struct rsb_tester_options_t * top);
rsb_err_t rsb_blas_bigger_matrices_tester(const struct rsb_tester_options_t * top);
rsb_err_t rsb_blas_limit_cases_tester(void);
rsb_err_t rsb_blas_runtime_limits_tester(void);
#endif /* LIBSPBLAS_TESTS_H_INCLUDED */
/* @endcond */
