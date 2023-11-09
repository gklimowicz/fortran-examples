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
 * @brief Initialization code.
 * @author Michele Martone
 * */
#ifndef RSB_INIT_H_INCLUDED
#define RSB_INIT_H_INCLUDED

#define RSB_DO_REINIT_SINGLE_VALUE(IOF,IOP,IOS,ERRVAL) { enum rsb_opt_t keys[]={IOF}; void*values[]={(IOP)}; struct rsb_initopts io; io.action=(IOS); io.keys=keys; io.values=values; io.n_pairs=1; ERRVAL=rsb__do_reinit(&io); }
#define RSB_DO_REINIT_SINGLE_VALUE_C_IOP(IOF,IOP,IOS,ERRVAL) { enum rsb_opt_t keys[]={IOF}; const void*values[]={(IOP)}; struct rsb_initopts io; io.action=(IOS); io.keys=keys; (io.values)=(void**)values; io.n_pairs=1; ERRVAL=rsb__do_reinit(&io); }
#define RSB_DO_REINIT_SINGLE_VALUE_SET(IOF,IOP,ERRVAL) RSB_DO_REINIT_SINGLE_VALUE(IOF,IOP,RSB_IO_SPECIFIER_SET,ERRVAL)
#define RSB_DO_REINIT_SINGLE_VALUE_GET(IOF,IOP,ERRVAL) RSB_DO_REINIT_SINGLE_VALUE(IOF,IOP,RSB_IO_SPECIFIER_GET,ERRVAL)

rsb_err_t rsb__init_mem_hierarchy_info(void);
rsb_err_t rsb__set_mem_hierarchy_info(const rsb_char_t * mhi);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_mem_hierarchy_info(void);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_err_t rsb__init_check_for_constants_correctness(void);
rsb_err_t rsb__init_check_for_system_constants_correctness(void);
const rsb_char_t * rsb__init_get_mem_hierarchy_info_string(rsb_bool_t verbose);
const rsb_char_t * rsb__get_mem_hierarchy_info_string(rsb_char_t *usmhib);
rsb_err_t rsb__do_init(struct rsb_initopts * io);
rsb_err_t rsb__do_reinit(struct rsb_initopts * io);
rsb_err_t rsb__do_exit(void);
#endif /* RSB_INIT_H_INCLUDED */
/* @endcond */
