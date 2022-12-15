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
/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief Functions for sparse format switch.
 * */

#ifndef RSB_SWT_H_INCLUDED
#define RSB_SWT_H_INCLUDED
#define RSB_MATRIX_STORAGE_AUTO 0x0	/* TODO: move to rsb_types.h */
#include "rsb_internals.h"		/* */
/*#define RSB_CONST_MIN_NNZ_PER_ROW_FOR_COO_SWITCH 4*/
#define RSB_CONST_MIN_NNZ_PER_ROW_FOR_COO_SWITCH 2

typedef unsigned short int rsb_half_idx_t;

rsb_bool_t rsb__do_is_candidate_size_for_halfword_coo(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_flags_t flags);
rsb_bool_t rsb__do_is_candidate_for_halfword_coo(const struct rsb_mtx_t * mtxAp);
rsb_err_t rsb__do_switch_to_halfword_coo(struct rsb_mtx_t * mtxAp);
rsb_err_t rsb__do_switch_to_fullword_zcoo(struct rsb_mtx_t * mtxAp);

rsb_err_t rsb__do_is_candidate_size_for_halfword(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_flags_t flags);
rsb_err_t rsb__do_is_candidate_size_for_halfword_csr(rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_flags_t flags);
rsb_err_t rsb__do_is_candidate_for_halfword_csr(const struct rsb_mtx_t * mtxAp);
#define RSB_FLAG_USE_FULLWORD_INDICES	0x00000000
rsb_err_t rsb__do_switch_leaf(struct rsb_mtx_t * mtxAp, rsb_fmt_t matrix_storage, rsb_flags_t flags, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t *TA);
rsb_err_t rsb__do_switch_to_halfword_csr(struct rsb_mtx_t * mtxAp);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_switch_to_fullword_csr(struct rsb_mtx_t * mtxAp);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/

/* #define RSB_COO_HALFWORDS_VALUES_PACK(LI,LJ)	((LJ)|((LI)<<RSB_COO_HALF_BITS_SIZE))*/ /* logical row and column index pack   (note: unused) */
/* #define RSB_COO_HALFWORDS_VALUES_UNPACK_LI(LIJ)	((((LIJ))>>RSB_COO_HALF_BITS_SIZE)&~((-1)<<RSB_COO_HALF_BITS_SIZE))*/	/* logical row index unpack (note: unused) */
#define RSB_COO_HALFWORDS_VALUES_UNPACK_UJ(LIJ)	(((LIJ) &~((~0U)<<RSB_COO_HALF_BITS_SIZE)))	/* logical row index unpack (note: shall get rid of this) */
void rsb__do_switch_array_to_halfword_coo(rsb_coo_idx_t  *p, rsb_nnz_idx_t n, const rsb_half_idx_t off);
void rsb__do_switch_array_to_fullword_coo(rsb_half_idx_t *p, rsb_nnz_idx_t n, const rsb_coo_idx_t off);


#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__do_is_candidate_for_fullword_coo(const struct rsb_mtx_t * mtxAp);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED*/
#endif /* RSB_SWT_H_INCLUDED */
/* @endcond */
