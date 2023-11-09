/*                                                                                                                            

Copyright (C) 2008-2020 Michele Martone

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
 * @brief Memory bandwidth related (e.g.: read, write, and read-write) microbenchmarks.
 */

#ifndef RSB_MBW_H_INCLUDED
#define RSB_MBW_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "rsb.h"		/* public API specification */

enum{
RSB_MB_INVALID		=0x99,
RSB_MB_READ		=0x00,
RSB_MB_WRITE		=0x01,
RSB_MB_RW		=0x02,
RSB_MB_BZERO		=0x03,
RSB_MB_ZERO		=0x04,
RSB_MB_MEMSET		=0x05,
RSB_MB_MEMCPY		=0x06,
RSB_MB_MEMCPY2		=0x07,
RSB_MB_LINEAR_CHASE	=0x08,
RSB_MB_MORTON_CHASE	=0x09,

RSB_MB_N		=0x0A,
RSB_MB_FLUSH		=0xAAAA
};

/**<
 * The available memory tests.
 */

#define RSB_MIN_CACHE_FLUSH_SCAN_TIMES 2

/*!
 * \ingroup gr_internals
 * \brief a memory bandwidth single measurement
 */
struct rsb_mbw_sm_t
{
	rsb_time_t t;			/**< time, in seconds */
	rsb_flags_t btype;		/**< measurement type */
};

#define RSB_MAX_TIMES_T INT_MAX

/*!
 * \ingroup gr_internals
 */
typedef size_t rsb__times_t;

/*!
 * \ingroup gr_internals
 * \brief a memory bandwidth benchmark record for a level
 */
struct rsb_mbw_m_t
{
	size_t so;			/**< sizeof probed word type */
	rsb__times_t times,sz;		/**< number of iterations of scanning a sz bytes wide area */
	struct rsb_mbw_sm_t mb[RSB_MB_N];	/**< measurements */
	long cln;			/**< number of cache levels */
	long cl;			/**< cache level for this measurement */
	long hlcs;			/**< higher level cache size */
};

/*!
 * \ingroup gr_internals
 * \brief  a complete memory bandwidth benchmark record
 */
struct rsb_mbw_cm_t
{
	struct rsb_mbw_m_t * mb;		/**< a memory bandwidth benchmark record for each level+extra_level */
	long cln;			/**< number of cache levels */
	long extra_level;		/**< number of additional measurements */
};

/*!
 * \ingroup gr_internals
 * \brief  export sample type
 */
struct rsb_mbw_es_t
{
	double bw; // bandwidth
	rsb_int_t sz,lvl; // size, level
	rsb_flags_t mbt; // test type
};

/*!
 * \ingroup gr_internals
 * \brief  export table type
 */
struct rsb_mbw_et_t
{
	struct rsb_mbw_es_t * et; // export sample type
	rsb_int_t sn; // samples number
};

#define RSB_FLUSH_TIMES 10	/* minimum number of scans for a "cache flush" intended array scan */
#define RSB_MEMSCAN_MIN_TIMES  2/* minimum number of scans for an array during a memory-bandwidth benchmarking */
#define RSB_MEMSCAN_MAX_TIMES RSB_MAX_SIGNED(int) /* minimum number of scans for an array during a memory-bandwidth benchmarking */
#define RSB_MEMSCAN_TIME  1.0	/* max allowable time scanning of an array during a memory-bandwidth benchmark */
#define RSB_MEMSCAN_EXTRA_LEVELS  2	/*  */
#define RSB_MEMSCAN_DEFAULT_LEVELS rsb__get_cache_levels_num() + RSB_MEMSCAN_EXTRA_LEVELS

rsb_err_t rsb__mem_hier_timings(struct rsb_mbw_cm_t * cm);
/* rsb_err_t rsb__print_mem_hier_timings(struct rsb_mbw_cm_t * cm); */
rsb_err_t rsb__memory_benchmark(struct rsb_mbw_et_t * mbetp);
rsb_err_t rsb__flush_cache(size_t fs);

rsb_err_t rsb__mbw_es_fill(struct rsb_mbw_et_t * mbet, const struct rsb_mbw_cm_t * cm);
rsb_err_t rsb__mbw_es_free(struct rsb_mbw_et_t * mbet);
rsb_err_t rsb__mbw_es_print(const struct rsb_mbw_et_t * mbetp);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__mbw_es_clone(struct rsb_mbw_et_t * mbetcp, const struct rsb_mbw_et_t * mbetsp);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_MBW_H_INCLUDED */
/* @endcond */
