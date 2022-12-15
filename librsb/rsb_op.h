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
 * This source file contains generic operation representation structures.
 * */

#ifndef RSB_OP_H_INCLUDED
#define RSB_OP_H_INCLUDED

#include "rsb_common.h"

	/* FIXME: NEW */
/*!
 * \ingroup gr_internals
 * \brief An internal, helper enumeration.
 */
enum rsb_opname_t{
            rsb_opn_spmv = 1,
            rsb_opn_spsv = 2,
            rsb_opn_scal = 3,	/* NOTE: this is a non const matrix op */
            rsb_opn_inrm = 4,
            rsb_opn_nop = 0
};

/*!
 * \ingroup gr_internals
 * \brief An internal, helper structure.
 */
struct rsb_c_mop_t
{
	/**
	 * const matrix op, of the form:
	 * y <- alpha A^(trans)  * x + beta * y
	 * y <- alpha A^(-trans) * y
	 * ...
	 */
	enum rsb_opname_t op;	/* the operation at hand				*/
	const struct rsb_mtx_t * matrix;	/* the operand matrix 			*/
	const void * alphap;	/* result vector post-scaling				*/
	const void * betap;	/* output vector pre-scaling (only if y!=x)	 	*/
	void * y;		/* (input) output vector 				*/
	const void * x;		/* input vector 					*/
	size_t nrhs;		/* number for right hand sides -- applies to both x,y 	*/
	rsb_trans_t trans;	/* transposition parameter 				*/
};

#endif /* RSB_OP_H_INCLUDED */
/* @endcond */

