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
/**
 * @file
 * @author Michele Martone
 * @brief Floating point microbenchmarks.
 */

#ifndef RSB_FPB_H_INCLUDED
#define RSB_FPB_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "rsb.h"		/* public API specification */

#define RSB_FPBENCH_TIME  2.0	/* min time for performing a floating point performance test on a type array  */
#define RSB_FPBENCH_MULTITYPE_TIME  ((RSB_FPBENCH_TIME)*(RSB_IMPLEMENTED_TYPES))	/* min time for performing a floating point performance test  */
rsb_err_t rsb__fp_benchmark(void);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_FPB_H_INCLUDED */
/* @endcond */
