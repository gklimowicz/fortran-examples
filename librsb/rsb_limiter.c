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
/* @cond INNERDOC  */
/**
 * @file
 * @author Michele Martone
 * @brief Timing/limiting mechanisms.
 */

#include "rsb_common.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

rsb_err_t rsb__limiter_init(struct rsb_limiter* lsp, const rsb_time_t max_time, const rsb__times_t max_times) 
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!lsp) { errval = RSB_ERR_BADARGS;goto err; }
	if( max_time  <  RSB_TIME_MIN ) { errval = RSB_ERR_BADARGS;goto err; }
	if( max_time  >  RSB_TIME_MAX ) { errval = RSB_ERR_BADARGS;goto err; }
	/*if( max_times < RSB_TIMES_MIN ) { errval = RSB_ERR_BADARGS;goto err; }*/
	if( max_times > RSB_TIMES_MAX ) { errval = RSB_ERR_BADARGS;goto err; }
	RSB_BZERO_P(lsp);
	lsp->max_time  = max_time;
	if( lsp->max_time  > RSB_TIME_ZERO )
		lsp->t0  = rsb_time();
	else
		lsp->t0  = RSB_TIME_ZERO;
	lsp->t1 = lsp->t0;
	lsp->max_times = max_times;
	lsp->times = RSB_TIMES_ZERO;
err:
	return errval;
}

rsb_err_t rsb__limiter_init_from_str(struct rsb_limiter* lsp, const char *tls)
{
	/* e.g.: tls = "4000" ; tls = "10s" */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_limiter lst;

	if(!tls || !*tls)
       		goto err;
	if(!lsp)
       		goto err;
	if(strstr(tls,"s")!=NULL)
	{
		lst.max_time  = rsb__util_atof(tls);
		lst.max_times = RSB_TIMES_ZERO;
	}
	else
	{
		lst.max_times = rsb__util_atoi(tls);
		lst.max_time  = RSB_TIME_ZERO;
	}
	errval = rsb__limiter_init(lsp,lst.max_time,lst.max_times); 
	goto ret;
err:
	errval = RSB_ERR_BADARGS;
ret:
	return errval;
}
 
rsb_err_t rsb__limiter_step(struct rsb_limiter* lsp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!lsp) { errval = RSB_ERR_BADARGS;goto err; }
	lsp->times++;
	if(lsp->max_time>RSB_TIME_ZERO)
		lsp->t1=rsb_time();
	/* FIXME: WRITE ME */
err:
	return errval;
}
 
rsb_bool_t rsb__limiter_done(const struct rsb_limiter* lsp)
{
	rsb_bool_t done = RSB_BOOL_TRUE;

	if( !lsp ) { goto err; }
	if( lsp->max_times > RSB_TIMES_ZERO && lsp->times >= lsp->max_times ) goto err;
	if( lsp->max_time  > RSB_TIME_ZERO  && (lsp->t1-lsp->t0) >= lsp->max_time ) goto err;
	if( lsp->max_times == RSB_TIMES_ZERO && lsp->max_time == RSB_TIME_ZERO ) goto err;
	done = RSB_BOOL_FALSE;
	/* FIXME: WRITE ME */
err:
	return done;
}

rsb_bool_t rsb__limiter_continue(const struct rsb_limiter* lsp)
{
	return RSB_BOOL_NOT(rsb__limiter_done(lsp));
}

rsb_err_t rsb__limiter_info(const struct rsb_limiter* lsp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_char_t*const tis="Timer info: ";

	if(!lsp) { errval = RSB_ERR_BADARGS;goto err; }
	if( lsp->max_time > RSB_TIME_ZERO )
		RSB_INFO("%s%lf / %lf seconds, %ld iterations.\n",tis,lsp->t1-lsp->t0,lsp->max_time,(long int)(lsp->times));
	else
	if( lsp->max_times > RSB_TIMES_ZERO )
		RSB_INFO("%s%d / %d iterations.\n",tis,(int)(lsp->times),(int)(lsp->max_times));
err:
	return errval;
}

/* @endcond */
