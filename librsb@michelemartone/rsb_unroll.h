/* @cond INNERDOC */
/**
 * @file
 * @brief
 * Unrolled kernels, for each type, operation, submatrix.
 * Right now, they are used for VBR and alike formats.
 */
/* Take care of compiling this code without loop unrolling optimizations (-fno-unroll-loops, or -ON with N<=2 on gcc) */

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
#ifndef RSB_UNROLL_H_INCLUDED
#define RSB_UNROLL_H_INCLUDED

#include "rsb_types.h"
#include "rsb_common.h"

/* NULL should be defined. */
#ifndef NULL
#define NULL ((int*)(0))
#endif /* NULL */
/**
 * No VBR/VBC formats compiled in.
 */





/**
 * Loops unroll factors.
 */
#define RSB_MIN_ROW_UNLOOP_FACTOR	1
#define RSB_MAX_ROW_UNLOOP_FACTOR	1
#define RSB_MIN_COLUMN_UNLOOP_FACTOR	1
#define RSB_MAX_COLUMN_UNLOOP_FACTOR	1

/**
 * No VBR/VBC formats compiled in.
 */


#endif  /* RSB_UNROLL_H_INCLUDED */




#define RSB_FITTING_SAMPLES		/*12 8*/4
/* @endcond */

