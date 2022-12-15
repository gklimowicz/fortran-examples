/*                                                                                                                            

Copyright (C) 2008-2015 Michele Martone

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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains functions for CSR handling.
 * */

#ifndef RSB_CSR_H_INCLUDED
#define RSB_CSR_H_INCLUDED

#include "rsb_internals.h"

rsb_err_t rsb__csr_chk(const rsb_nnz_idx_t * RSB_RESTRICT IP, const rsb_coo_idx_t * RSB_RESTRICT JA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t ib);
rsb_err_t rsb__csc_chk(const rsb_nnz_idx_t * RSB_RESTRICT IP, const rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_nnz_idx_t nnzA, rsb_coo_idx_t ib);


#endif /* RSB_CSR_H_INCLUDED */
/* @endcond */
