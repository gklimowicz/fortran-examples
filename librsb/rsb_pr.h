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
 * @brief Perfomance reporting code.
 * @author Michele Martone
 * */

#ifndef RSB_PR_H_INCLUDED
#define RSB_PR_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include "rsb_internals.h"

struct rsb_ts_t /* time statistics struct; dumpable with RSB_STAT_DUMP_TS */
{
	rsb_time_t avg,min,max,sd;
	rsb_int_t ns;
};

rsb_err_t rsb__pr_init(void**rsprpv, const struct rsb_mtx_t *mtxAp, rsb_int_t filenamen, rsb_int_t cn, rsb_int_t incXn, rsb_int_t incYn, rsb_int_t nrhsn, rsb_int_t ntypecodes, rsb_int_t tn, struct rsb_mbw_et_t * mbetp);
rsb_err_t rsb__pr_set(void*rsprpv, const struct rsb_mtx_t *mtxAp, const struct rsb_mtx_t *at_mtxAp, rsb_int_t filenamei, rsb_int_t ci, rsb_int_t incXi, rsb_int_t incYi, rsb_int_t nrhsi, rsb_int_t typecodesi, rsb_int_t ti, rsb_trans_t transA, rsb_perf_t op_time_best, rsb_perf_t mkl_csr_op_time_best, rsb_perf_t at_op_time_best, rsb_perf_t at_mkl_csr_op_time_best, rsb_int_t at_cn, rsb_int_t at_mkl_csr_cn, rsb_time_t at_t, rsb_int_t at_eps, const struct rsb_ts_t*otposp, const struct rsb_ts_t*btposp, const struct rsb_ts_t*otpmsp, const struct rsb_ts_t*btpmsp);
rsb_err_t rsb__pr_dump(const void*rsprpv, rsb_char_t**RSB_RESTRICT filenamea, rsb_int_t*ca, const rsb_int_t*incXa, const rsb_int_t*incYa, const rsb_int_t*nrhsa, const rsb_type_t*typecodes, const rsb_int_t *ta, const rsb_char_t *fprfn);
rsb_err_t rsb__pr_dump_inner(const void*rsprpv, rsb_char_t**RSB_RESTRICT filenamea, rsb_int_t*ca, const rsb_int_t*incXa, const rsb_int_t*incYa, const rsb_int_t*nrhsa, const rsb_type_t*typecodes, const rsb_int_t*ta, const int*filenameifp, const int*ifilenameifp, const int*cifp , const int*incXifp , const int*incYifp , const int*nrhsifp , const int*typecodefip , const int*tifp, const rsb_trans_t*tfp, rsb_flags_t flagsA, rsb_flags_t nflagsA, rsb_char_t *ltag, const rsb_char_t *fprfn);
rsb_err_t rsb__pr_free(void*rsprpv);
rsb_err_t rsb__pr_dumpfiles(const rsb_char_t **argv, const int argc);
rsb_err_t rsb__pr_save(const rsb_char_t * RSB_RESTRICT filename, /*const*/ void * RSB_RESTRICT rsprpv, rsb_char_t**RSB_RESTRICT filenamea, rsb_int_t*RSB_RESTRICT ca, const rsb_int_t*RSB_RESTRICT incXa, const rsb_int_t*RSB_RESTRICT incYa, const rsb_int_t*RSB_RESTRICT nrhsa, const rsb_type_t*RSB_RESTRICT typecodes, const rsb_int_t*RSB_RESTRICT ta, rsb_bool_t can_overwrite);
rsb_bool_t rsb__file_exists(const rsb_char_t * RSB_RESTRICT filename);

rsb_err_t rsb__rspr_rw(void * p, size_t sop, FILE * stream, rsb_bool_t row);
void rsb__mtxfn_bncp(char* dst, const char*src, int lm);
#define RSB_PR_WR RSB_BOOL_FALSE
#define RSB_PR_RD RSB_BOOL_TRUE

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_PR_H_INCLUDED */

/* @endcond */
