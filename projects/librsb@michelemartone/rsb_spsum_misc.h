/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * */

/*

Copyright (C) 2008-2022 Michele Martone

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
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

 */

#ifndef RSB_SPSUM_MISC_H_INCLUDED
#define RSB_SPSUM_MISC_H_INCLUDED
#include "rsb_common.h"


rsb_err_t rsb__do_add_submatrix_to_dense(const struct rsb_mtx_t * mtxAp, const void *alphap, void * Bp, rsb_nnz_idx_t ldb, rsb_nnz_idx_t nr, rsb_nnz_idx_t nc, rsb_bool_t rowmajor)
;
#endif /* RSB_SPSUM_MISC_H_INCLUDED */
/* @endcond */
