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
 * @author Michele Martone
 * @brief Code for matrix format conversion. 
 * */
#ifndef RSB_REC2COO_H_INCLUDED
#define RSB_REC2COO_H_INCLUDED
#include "rsb_common.h"
rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_rcoo(struct rsb_mtx_t * mtxAp, rsb_bool_t do_shift);
rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_sorted(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop);
rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_coo_unsorted(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop);
#endif /* RSB_REC2COO_H_INCLUDED */
/* @endcond */
