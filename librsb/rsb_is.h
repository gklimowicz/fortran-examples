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
 * This source file contains mtxAp info getter functions.
 * */

#ifndef RSB_IS_H_INCLUDED
#define RSB_IS_H_INCLUDED

#include "rsb_internals.h"

#define RSB__IS_SUBM_MISPLACED(SUBM,MTXAP) ( ((MTXAP)->roff>(SUBM)->roff) || ((MTXAP)->coff>(SUBM)->coff) || ((MTXAP)->nr<(SUBM)->nr)     || ((MTXAP)->nc<(SUBM)->nc) )

rsb_bool_t rsb__is_coo_matrix(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_square(const struct rsb_mtx_t *mtxAp);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__is_fixed_block_matrix(const struct rsb_mtx_t *mtxAp);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_bool_t rsb__is_css_matrix(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_bcsr_matrix(const struct rsb_mtx_t *mtxAp);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__is_bcss_matrix(const struct rsb_mtx_t *mtxAp);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_bool_t rsb__is_bcsc_matrix(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_recursive_matrix(rsb_flags_t flags);
rsb_bool_t rsb__is_terminal_recursive_matrix(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_csr_matrix(const struct rsb_mtx_t *mtxAp);
#ifdef RSB_FLAG_WANT_LINKED_STORAGE
rsb_bool_t rsb__have_linked_storage(const rsb_flags_t flags);
#endif
rsb_bool_t rsb__have_fixed_blocks_matrix_flags(rsb_flags_t flags);
rsb_bool_t rsb__util_are_flags_suitable_for_optimized_1x1_constructor(rsb_flags_t flags);
rsb_bool_t rsb__is_symmetric(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_not_unsymmetric(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_root_matrix(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_hermitian(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__is_lower_triangle(rsb_flags_t flags);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_bool_t rsb__is_triangle(rsb_flags_t flags);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/
rsb_bool_t rsb__is_upper_triangle(rsb_flags_t flags);
rsb_bool_t rsb__mtx_chk(const struct rsb_mtx_t *mtxAp);
rsb_bool_t rsb__do_is_matrix_binary_loaded(const struct rsb_mtx_t * mtxAp);

#endif /* RSB_IS_H_INCLUDED */
/* @endcond */
