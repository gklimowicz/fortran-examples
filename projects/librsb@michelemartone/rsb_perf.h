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
 * @brief Perfomance tuning or measuring code.
 * @author Michele Martone
 * */

#ifndef RSB_PERF_H_INCLUDED
#define RSB_PERF_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include "rsb_internals.h"

/*! \brief Performance info for a single matrix operation, all possible unrollings. */
struct rsb_mop_performance_info_t
{
	/* TODO : should we add 'flags' and 'runs' field ? */

	/** some matrix info */
	size_t rows,cols,nnz,element_count;

	/** training matrix info */
	rsb_flags_t flags,storage;
	rsb_type_t typecode;

	/** millions of floating point operations */
	double m_flops[RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH];

	/** millions of effective floating point operations (==m_flops/fillin) */
	double e_mflops[RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH];

	/** fillin */
	double fillin[RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH];

	/** time in seconds */
	double seconds[RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH];
};
/*! \brief Reference performance info for a single matrix operation, all possible unrollings. */
struct rsb_mop_reference_performance_info_t
{
	/** performance info per fitting sample */
	struct rsb_mop_performance_info_t pipfs[RSB_FITTING_SAMPLES];

	/** blocks per row density              */
	double                        blocks_per_row[RSB_FITTING_SAMPLES];

	/** alpha, beta, gamma parameterization as in the accels experimental setup*/
	double alpha[RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH];
	double beta [RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH];
	double gamma[RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH];
};
/*! \brief Performance info for multiple matrix operations. */
struct rsb_mops_performance_info_t
{
	/** performance info per matrix operation */
	struct rsb_mop_performance_info_t	pipmo[RSB_IMPLEMENTED_META_MOPS];
};
/*! \brief Reference performance info for multiple matrix operations. */
struct rsb_mops_reference_performance_info_t
{
	/** performance info per matrix operation */
	struct rsb_mop_reference_performance_info_t pipmo[RSB_IMPLEMENTED_META_MOPS];
};
/*! \brief Global performance info for all matrix operations and types. */
struct rsb_global_performance_info_t
{
	/** global performance info */
	struct rsb_mops_performance_info_t	gpi[RSB_IMPLEMENTED_TYPES];
};
/*! \brief Global reference performance info for all matrix operations and types. */
struct rsb_global_reference_performance_info_t
{
	rsb_bool_t initialized; /** if not zero, measurements should be considered valid **/
	/** global performance info */
	struct rsb_mops_reference_performance_info_t gpi[RSB_IMPLEMENTED_TYPES];
};

rsb_err_t rsb__print_mop_reference_performance_info_header(void);	/* temporary */
rsb_err_t rsb__print_mop_reference_performance_info(const struct rsb_mop_reference_performance_info_t *pi, char *s);	/* temporary */
rsb_err_t rsb__dump_global_performance_info(const struct rsb_global_performance_info_t *gpip);	/* temporary */
#if RSB_WANT_PERFORMANCE_FILE
//rsb_err_t rsb__dump_global_reference_performance_info(const struct rsb_global_reference_performance_info_t *gpip);	/* temporary */
rsb_err_t rsb__save_global_reference_performance_info(const struct rsb_global_reference_performance_info_t *gpip);
#endif /* RSB_WANT_PERFORMANCE_FILE */
rsb_err_t rsb__dump_performance_info(const struct rsb_mop_performance_info_t * pi, const char * pid);
rsb_err_t rsb__dump_performance_info_line(const struct rsb_mop_performance_info_t * pi);/* new */
rsb_err_t rsb__dump_current_global_reference_performance_info(void);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__pinfo_init(struct rsb_mtx_partitioning_info_t * pinfop,
	rsb_blk_idx_t M_b, rsb_blk_idx_t K_b,
	rsb_coo_idx_t *rpntr,rsb_coo_idx_t *cpntr,
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_blk_idx_t br, rsb_blk_idx_t bc);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_err_t rsb__dump_system_performance_summary(void);

#define RSB_PERFORMANCE_BINARY_DUMP_FILE "rsb_performance_profile.bin"

rsb_err_t rsb__perf_init(void);
rsb_err_t rsb__perf_exit(void);
size_t rsb_spmv_memory_accessed_bytes(const struct rsb_mtx_t * mtxAp);
size_t rsb_spmv_memory_accessed_bytes_min(const struct rsb_mtx_t * mtxAp);
size_t rsb_spmv_memory_accessed_bytes_max(const struct rsb_mtx_t * mtxAp);
double rsb_spmv_memory_accessed_bytes_wr_ratio(const struct rsb_mtx_t * mtxAp);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_nnz_idx_t rsb__fillin_estimation_nnz_count(
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, 
	const  rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags, rsb_int nprobes
);/* NEW */
rsb_err_t rsb__estimate_expected_fillin_for_blocking(
	const void * VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, 
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	const  rsb_nnz_idx_t nnz, rsb_type_t typecode, 
	rsb_flags_t flags,
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	double *efillinp);/* NEW */
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__estimate_expected_raw_performance_for_blocking(
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	const  rsb_nnz_idx_t nnz, rsb_type_t typecode, 
	rsb_flags_t flags,
	double efillin,
	double*eperf);

size_t rsb_spmv_memory_accessed_bytes_(
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_nnz_idx_t element_count,
	rsb_nnz_idx_t block_count,
	rsb_blk_idx_t Mdim,
	size_t el_size
);
rsb_err_t rsb__dump_performance_record(const char * s, const struct rsb_mtx_t * mtxAp, rsb_real_t rsb_NMflops_ps, rsb_real_t rsb_RMflops_ps, const char *op, rsb_flags_t inflags);
FILE *rsb__util_fopen(const char *path, const char *mode);

/* 
	NEW, EXPERIMENTAL
	TODO : should be user-specified
 */
#define RSB_FIRST_FITTING_SAMPLE_BW_MIN 10
#define RSB_FIRST_FITTING_SAMPLE_BW_MAX 100

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_PERF_H_INCLUDED */
/* @endcond */
