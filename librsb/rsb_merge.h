/* @cond INNERDOC */
/**
 * @file
 * @brief Auxiliary functions.
 */

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


#ifndef RSB_MERGE_H_INCLUDED
#define RSB_MERGE_H_INCLUDED

#include "rsb_common.h"
#ifndef RSB_PS_ASSERT
/*#define RSB_PS_ASSERT(e) assert(e)	*/ 	/* uncomment this to use   asserts */
#define RSB_PS_ASSERT(e)			/* uncomment this to avoid asserts */
#else /* RSB_PS_ASSERT */
#undef RSB_PS_ASSERT
#define RSB_PS_ASSERT(e) 
#endif /* RSB_PS_ASSERT */




rsb_err_t rsb__do_util_merge_sorted_subarrays_in_place(
		void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_char_t * W,
		rsb_nnz_idx_t annz, rsb_nnz_idx_t bnnz,
	       	size_t wsize, rsb_flags_t flags, rsb_type_t typecode
		)
;
rsb_err_t rsb__do_util_merge_sorted_subarrays_in_place_test(void)
;
#endif /* RSB_MERGE_H_INCLUDED */



/* @endcond */
