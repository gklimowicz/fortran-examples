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
/* @cond INNERDOC */

/**
 * @file
 * @author Michele Martone
 * @brief Sorting routines.
 */

#ifndef RSB_MSORT_UP_H_INCLUDED
#define RSB_MSORT_UP_H_INCLUDED

#include "rsb_internals.h"

rsb_err_t rsb__do_msort_up(rsb_nnz_idx_t n, const rsb_nnz_idx_t * RSB_RESTRICT k, rsb_nnz_idx_t * RSB_RESTRICT l);
rsb_err_t rsb__do_msort_up2coo(rsb_nnz_idx_t n, const rsb_coo_idx_t * k, rsb_nnz_idx_t * l);

#endif /* RSB_MSORT_UP_H_INCLUDED */
/* @endcond */
