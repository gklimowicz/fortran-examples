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
 * @author Michele Martone
 * @brief Timing/limiting mechanisms.
 */

#ifndef RSB_LIMITER_H_INCLUDED
#define RSB_LIMITER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "rsb.h"		/* public API specification */
#define RSB_TIMES_MAX 1000000000
#define RSB_TIMES_MIN 0
#define RSB_TIMES_ZERO 0
#define RSB_TIME_MAX RSB_CONST_IMPOSSIBLY_BIG_TIME
#define RSB_TIME_MIN RSB_TIME_ZERO 
struct rsb_limiter
{
	/*rsb_bool_t is_time_based;*/
	rsb_time_t t0,t1;
	rsb_time_t max_time;
	rsb__times_t max_times;
	rsb__times_t times;
};
rsb_err_t rsb__limiter_init(struct rsb_limiter* lsp, const rsb_time_t max_time, const rsb__times_t max_times);
rsb_err_t rsb__limiter_init_from_str(struct rsb_limiter* lsp, const char *tls);
rsb_err_t rsb__limiter_step(struct rsb_limiter* lsp);
rsb_bool_t rsb__limiter_done(const struct rsb_limiter* lsp);
rsb_bool_t rsb__limiter_continue(const struct rsb_limiter* lsp);
rsb_err_t rsb__limiter_info(const struct rsb_limiter* lsp);
#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_LIMITER_H_INCLUDED */
/* @endcond */
