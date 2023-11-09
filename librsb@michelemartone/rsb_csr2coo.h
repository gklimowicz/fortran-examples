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
 * @brief CSR to COO conversion code
 * @author Michele Martone
 * */
#ifndef RSB_CSR2COO_H_INCLUDED
#define RSB_CSR2COO_H_INCLUDED
#include "rsb_common.h"

void rsb__do_prefix_sum_coo_idx_t(rsb_nnz_idx_t *IA, rsb_nnz_idx_t nnz);
#define rsb_do_prefix_sum_nnz_idx_t rsb__do_prefix_sum_coo_idx_t
rsb_err_t rsb__do_switch_compressed_array_to_fullword_coo(rsb_nnz_idx_t *RSB_RESTRICT IP, rsb_nnz_idx_t m, rsb_coo_idx_t off, rsb_coo_idx_t *RSB_RESTRICT TA);
rsb_nnz_idx_t rsb__do_copy_lowtri_from_csr_to_coo(const struct rsb_coo_mtx_t *csrp, struct rsb_coo_mtx_t *coop);
rsb_nnz_idx_t rsb__do_copy_upptri_from_csr_to_coo(const struct rsb_coo_mtx_t *csrp, struct rsb_coo_mtx_t *coop);
rsb_nnz_idx_t rsb__do_count_upptri_in_csr(const struct rsb_coo_mtx_t *csrp);
rsb_nnz_idx_t rsb__do_count_tri_in_csr(const struct rsb_coo_mtx_t *csrp, rsb_nnz_idx_t *lnzp, rsb_nnz_idx_t *unzp);
rsb_nnz_idx_t rsb__do_copy_tri_from_csr_to_coo(const struct rsb_coo_mtx_t *csrp, struct rsb_coo_mtx_t *lcoop, struct rsb_coo_mtx_t *ucoop);
rsb_err_t rsb__do_switch_fullword_array_to_compressed(rsb_nnz_idx_t *IA, rsb_nnz_idx_t nnz, rsb_nnz_idx_t m);
#endif /* RSB_CSR2COO_H_INCLUDED */
/* @endcond */
