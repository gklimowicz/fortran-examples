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
 * @brief Code for matrix format conversion. 
 * @author Michele Martone
 * */
#ifndef RSB_REC2CSR_H_INCLUDED
#define RSB_REC2CSR_H_INCLUDED
#include "rsb_common.h"
rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_csr(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop);
rsb_err_t rsb__do_switch_recursive_in_place_matrix_to_in_place_csc(struct rsb_mtx_t * mtxAp, struct rsb_coo_mtx_t * coop);
#endif /* RSB_REC2CSR_H_INCLUDED */
/* @endcond */
