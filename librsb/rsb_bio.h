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
 * This source file contains matrix binary I/O functions.
 * */

#ifndef RSB_BIO_H_INCLUDED
#define RSB_BIO_H_INCLUDED
#include "rsb_common.h"
rsb_err_t rsb__do_bindump_init(void);
rsb_err_t rsb__do_load_matrix_file_as_binary(struct rsb_mtx_t ** mtxApp, const rsb_char_t * filename);
rsb_err_t rsb__do_save_matrix_file_as_binary(const struct rsb_mtx_t * mtxAp, FILE * fd);
#endif /* RSB_BIO_H_INCLUDED */
/* @endcond */
