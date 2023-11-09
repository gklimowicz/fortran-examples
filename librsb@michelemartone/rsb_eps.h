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
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains postscript rendering functions.
 * */
#ifndef RSB_EPS_H_INCLUDED
#define RSB_EPS_H_INCLUDED
#include "rsb_internals.h"
/* rsb_err_t rsb__dump_postscript_from_matrix(const char * filename, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_bool_t all_nnz); */
rsb_err_t rsb__dump_postscript_recursion_from_matrix(const char * filename, rsb_marf_t rflags, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_flags_t flags, rsb_bool_t want_blocks, rsb_bool_t z_dump , rsb_bool_t want_nonzeros, rsb_bool_t want_recursion, rsb_type_t typecode);
rsb_err_t rsb__dump_postscript_recursion_from_mtx_t(FILE*fd, const char * filename, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_marf_t rflags, rsb_bool_t want_blocks, rsb_bool_t z_dump, rsb_bool_t want_nonzeros, const rsb_submatrix_idx_t *pv, const struct rsb_optrace_t * otv, const rsb_bitmap_data_t * smb);
rsb_err_t rsb__dump_multiple_recursion_postscript_from_mtx_t(const char * filename, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_marf_t rflags, rsb_bool_t want_blocks, rsb_bool_t z_dump, rsb_bool_t want_nonzeros, const rsb_submatrix_idx_t *pv, const struct rsb_optrace_t * otv );
#define RSB_POSTSCRIPT_DUMP_COMMENT(FD,C) {RSB_FPRINTF(FD,"%%%% ");RSB_FPRINTF(FD,C);RSB_FPRINTF(FD,"\n");}
rsb_err_t rsb__do_mtx_render(const char * filename, const struct rsb_mtx_t*mtxAp, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags);
#endif /* RSB_EPS_H_INCLUDED */
/* @endcond */
