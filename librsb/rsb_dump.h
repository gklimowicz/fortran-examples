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
 * @brief Matrix info dumping code
 * @author Michele Martone
 * */

#ifndef RSB_DUMP_H_INCLUDED
#define RSB_DUMP_H_INCLUDED

/*#include "rsb_common.h"*/
#include "rsb.h"

#define RSB_DEFAULT_DUMPFILENAME "/dev/stdout"
#define RSB_DEFAULT_FD stdout
typedef rsb_flags_t rsb_dump_flags_t;
#define RSB_CONST_DUMP_RECURSION	0x00000001
#define RSB_CONST_DUMP_INTERNALS	0x00000002	/* NOTE: free */
#define RSB_CONST_DUMP_TIMES		0x00000004
#define RSB_CONST_DUMP_BLOCKS		0x00000008	/* NOTE: free */
#define RSB_CONST_DUMP_MATRIX_MARKET	0x00000010
#define RSB_CONST_DUMP_DIMENSIONS	0x00000020
#define RSB_CONST_DUMP_COO		0x00000040	/* NOTE: unused */
#define RSB_CONST_DUMP_RECURSION_BRIEF	0x00000080
#define RSB_CONST_DUMP_OCTAVE_STYLE	0x00000100	/* FIXME: unfinished */
#define RSB_CONST_DUMP_CSR		0x00000200
#define RSB_CONST_DUMP_RSB		0x00000400
#define RSB_CONST_DUMP_DOT		0x00000800
#define RSB_CONST_DUMP_MIF_MTX_INFO     0x00001000	/* see RSB_MIF_MATRIX_INFO__TO__CHAR_P */
#define RSB_CONST_DUMP_FLAGS		0x00002000
#define RSB_CONST_DUMP_DEFAULT		0x00000000
rsb_err_t rsb__do_print_matrix_stats(const struct rsb_mtx_t *mtxAp, rsb_dump_flags_t flags, const rsb_char_t*filename);
rsb_err_t rsb__do_dump_bitmap(const rsb_bitmap_data_t * bmap, size_t w, size_t h);
#endif /* RSB_DUMP_H_INCLUDED */
/* @endcond */
