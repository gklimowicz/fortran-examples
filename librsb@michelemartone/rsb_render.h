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
 * This source file contains pixmap rendering functions.
 * */

#ifndef RSB_RENDER_H_INCLUDED
#define RSB_RENDER_H_INCLUDED

#include "rsb_internals.h"	/* rsb_coo_mtx_t */

rsb_err_t rsb__do_get_pixmap_RGB_from_matrix(const char * filename, void * pixmap, int width, int height);
rsb_err_t rsb__mtx_as_pixmap_resize(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_nnz_idx_t *rnnz, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_coo_idx_t p_rows, rsb_coo_idx_t p_cols, rsb_type_t typecode, rsb_flags_t render_flags);
rsb_err_t rsb__do_print_postscript_header(FILE*fd, rsb_marf_t rflags, int width, int height, float csw, float csh);

#define RSB_DEFAULT_MATRIX_RENDERING_ROWS 512
#define RSB_OPS_RENDERING_EXTRA_FRAC 4
#define RSB_DEFAULT_MATRIX_RENDERING_COLS 512

#endif /* RSB_RENDER_H_INCLUDED */

/* @endcond */
