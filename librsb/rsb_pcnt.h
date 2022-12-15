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
 * @brief Perfomance tuning or measuring code.
 * @author Michele Martone
 * */

#ifndef RSB_PCNT_H_INCLUDED
#define RSB_PCNT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include "rsb_internals.h"
#include "rsb_perf.h"

#ifdef RSB_HAVE_PAPI
#include <papi.h>		/*  http://icl.cs.utk.edu/papi/ */
typedef  long_long rsb_papi_long; /* long_long is a typedef originating in the papi headers */
typedef int rsb_papi_int_t;
typedef int rsb_papi_err_t;
#define RSB_PC_MAX_ITEMS 3
#endif /* RSB_HAVE_PAPI */
struct rsb_pci_t
{
	int eventnum;
#ifdef RSB_HAVE_PAPI
	rsb_papi_int_t eventlist[RSB_PC_MAX_ITEMS];
	rsb_papi_long eventvals[RSB_PC_MAX_ITEMS];
	char          eventdesc[RSB_PC_MAX_ITEMS][PAPI_MAX_STR_LEN];
#endif /* RSB_HAVE_PAPI */
};

#if defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS>1)
rsb_err_t rsb_perf_counters_init(void);
rsb_err_t rsb_perf_counters_finalize(void);
#endif /* defined(RSB_WANT_PERFORMANCE_COUNTERS) && (RSB_WANT_PERFORMANCE_COUNTERS>1) */
rsb_err_t rsb_perf_counters_dump(const rsb_char_t *premsg, const rsb_char_t *postmsg, rsb_int_t tdiv, struct rsb_pci_t *pcip);

rsb_err_t rsb_hc_main(void);	/* preliminary */

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_PCNT_H_INCLUDED */
/* @endcond */
