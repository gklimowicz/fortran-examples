/*                                                                                                                            

Copyright (C) 2008-2019 Michele Martone

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
 * @brief Copy/Move primitives 
 * @author Michele Martone
 * */
#ifndef RSB_CPMV_H_INCLUDED
#define RSB_CPMV_H_INCLUDED
#include "rsb_common.h"
/*void RSB_A_MEMCPY_parallel(rsb_char_t * ID, const rsb_char_t * IS, size_t DOFF, size_t SOFF, size_t NNZ, size_t ES);*/
/*void RSB_COA_MEMCPY_parallel(rsb_char_t * ID, const rsb_char_t * IS, size_t DOFF, size_t SOFF, size_t NNZ);*/
#define RSB_A_MEMCPY_parallel rsb__a_memcpy_parallel
#define RSB_COA_MEMCPY_parallel rsb__coa_memcpy_parallel
#define RSB_BZERO_parallel rsb__bzero_parallel
void RSB_A_MEMCPY_parallel(void * RSB_RESTRICT ID, const void * RSB_RESTRICT IS, size_t DOFF, size_t SOFF, size_t NNZ, size_t ES);
void RSB_COA_MEMCPY_parallel(void* ID, const void* IS, size_t DOFF, size_t SOFF, size_t NNZ);
void RSB_BZERO_parallel(void * p, size_t n);
/* @endcond */
#endif /* RSB_CPMV_H_INCLUDED */
