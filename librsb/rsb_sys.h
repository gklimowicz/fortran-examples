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
/* @cond INNERDOC */
/**
 * @file
 * @brief System, or standard library related functions.
 * @author Michele Martone
 * */
#ifndef RSB_SYS_H_INCLUDED
#define RSB_SYS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdlib.h>
#include <time.h>	/* clock, difftime */
#include <sys/time.h>	/* timeval,gettimeofday */
#include <stdio.h>	/* printf */
#include <strings.h>	/* formerly bzero (now using memset) */
#include <string.h>	/* memset, memcpy */
#include <stdarg.h>	/* vprintf, va_start, va_end */
#if RSB_HAVE_LIMITS_H 
#include <limits.h>	/* CHAR_BIT */
#endif /* RSB_HAVE_LIMITS_H  */
#include "rsb_perf.h"

void * rsb__malloc(size_t size);
void * rsb__calloc(size_t n);
void * rsb__calloc_parallel(size_t n);
rsb_time_t rsb__do_time(void);
rsb_time_t rsb__timer_granularity(void );
void * rsb__free(void *rsb_data);
void * rsb__realloc(void *rsb_data, size_t size);
void * rsb__do_realloc(void *rsb_data, size_t size, size_t alignment);
void * rsb__aligned_malloc(size_t size, size_t alignment);
rsb_err_t rsb__sys_info(void);
long rsb__get_l1c_size(void);
#if RSB_WANT_EXPERIMENTS_CODE
long rsb__get_l2c_size(void);
#endif /* RSB_WANT_EXPERIMENTS_CODE */
#if RSB_OBSOLETE_QUARANTINE_UNUSED
long rsb__get_l3c_size(void);
long rsb__get_l4c_size(void);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
long rsb__get_lastlevel_c_size(void);
long rsb__get_first_level_c_size(void);
rsb_int_t rsb__get_cache_levels_num(void);
long rsb__get_lnc_size(int n);
#if RSB_OBSOLETE_QUARANTINE_UNUSED
long rsb__know_cache_sizes(void);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_err_t rsb__sys_init(void);
/* long rsb_want_executing_threads(void); */

#ifndef RSB_CONDITIONAL_FREE_FAKE
#define RSB_CONDITIONAL_FREE_FAKE(p) {if((p))/*rsb__free((p))*/;(p)=NULL;}
#endif /* RSB_CONDITIONAL_FREE_FAKE */
/* A useful macro */

#if defined(RSB_MEM_DBG) && (RSB_MEM_DBG>1) 
#define RSB__FREE(p) {RSB_ERRLN("About to free memory.\n");rsb__free(p);} 
#else /* RSB_MEM_DBG */
#define RSB__FREE(p) {                                    ;rsb__free(p);} 
#endif /* RSB_MEM_DBG */

#if RSB_WANT_CACHE_TIMER_GRANULARITY
#define RSB_CACHED_TIMER_GRANULARITY rsb_global_session_handle.timer_granularity
#else /* RSB_WANT_CACHE_TIMER_GRANULARITY */
#define RSB_CACHED_TIMER_GRANULARITY rsb__timer_granularity()
#endif /* RSB_WANT_CACHE_TIMER_GRANULARITY */

#ifndef RSB_CONDITIONAL_FREE
#define RSB_CONDITIONAL_FREE(p) {if((p)){RSB__FREE((p));(p)=NULL;}} /* frees and nullifies the associated pointer. */
/*#define RSB_CONDITIONAL_FREE(p) RSB_CONDITIONAL_FREE_FAKE(p) */
#endif /* RSB_CONDITIONAL_FREE */

extern int rsb__error(const char * format, ...);
void rsb__g_rsb_memory_counter_init(void);
size_t rsb__get_g_rsb_allocations_count(void);

/* ... and what about __PRETTY_FUNCTION__ ? */
#ifndef __func__
#ifdef __STDC_VERSION__
#if __STDC_VERSION__ < 199901L
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else /* __GNUC__  */
#  define __func__ "<unknown>"
# endif /* __GNUC__  */
#endif /* __STDC_VERSION__ */
#else /* __STDC_VERSION__ */
# define __func__ "<unknown>"
#endif /* __STDC_VERSION__ */
#endif /* __func__ */

/* '1' is better than '{;}' in situations like : if(..)RSB_INFO(...);else RSB_INFO(...); */
#define RSB_NULL_EXPRESSION_FOR_ZEN_HAPPINESS 1
#define RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS RSB_NULL_EXPRESSION_FOR_ZEN_HAPPINESS

#define RSB_FFL_PRINTF printf("In %s located in %20s:%d :\n",__func__,__FILE__,__LINE__)

#define RSB_ERRLN( ... ) \
	RSB_FFL_PRINTF,	\
	rsb__error( __VA_ARGS__ )

#if RSB_INT_ERR_VERBOSITY==1

#define RSB_DEPRECATED( ... ) \
	RSB_ERRLN( __VA_ARGS__ ), \
	printf(" is DEPRECATED !!\n") \

#define RSB_ERROR( ... ) RSB_ERRLN( __VA_ARGS__ )

#define RSB_OCTAVE_ERROR( ... ) \
	RSB_ERROR("ERROR:"), RSB_ERROR( __VA_ARGS__ ),octave_failed_tests++;

#else /* RSB_INT_ERR_VERBOSITY */
#define RSB_DEPRECATED( ... )  RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS 
#define RSB_ERROR( ... ) RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS 
#define RSB_OCTAVE_ERROR( ... ) RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS 
#endif /* RSB_INT_ERR_VERBOSITY */

#define RSB_IOLEVEL RSB_WANT_IO_LEVEL 
/*#define RSB_IOLEVEL 7*/ /* 0 is minimum, 7 is maximum */
#define RSB_ALLOW_FPRINTF (RSB_IOLEVEL&4)
#define RSB_ALLOW_STDOUT  (RSB_IOLEVEL&1)
#define RSB_ALLOW_STDERR  (RSB_IOLEVEL&2)

#if RSB_ALLOW_STDERR
/* WARNING : calling this without arguments causes segfaults! */
#define RSB_STDERR( ... ) fprintf(stderr, __VA_ARGS__ )
#else /* RSB_ALLOW_STDERR */
#define RSB_STDERR( ... ) RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS 
#endif /* RSB_ALLOW_STDERR */
#define RSB_IO_ERROR RSB_STDERR
#define RSB_IO_NOTICE RSB_STDERR

#if RSB_ALLOW_STDOUT
/* explicit standard output printout */
#define RSB_STDOUT( ... ) fprintf(stdout, __VA_ARGS__ )

/* */
/*#define RSB_DEBUGINFO( ... ) printf("%s @ %10d (%s):\n",__FILE__,__LINE__,__func__),RSB_STDOUT(__VA_ARGS__)*/

/** RSB_WARN is used in not-yet-implemented-feature-but-skip-error-triggering situations.  */
#define RSB_WARN( ... ) \
	RSB_STDOUT("%s\n#","#*****************************************************************************"),\
	RSB_STDOUT( __VA_ARGS__ ),\
	RSB_STDOUT("%s\n","#*****************************************************************************")
#else
#define RSB_STDOUT( ... ) RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS 
#define RSB_WARN( ... ) RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS 
#endif /* RSB_ALLOW_STDOUT */
#define RSB_WARN_CRITICAL( ... ) \
	RSB_STDERR("%s%s%s\n","#********************",RSB_HEADER_VERSION_STRING,"**********************************"),\
	RSB_STDERR( __VA_ARGS__ ),\
	RSB_STDERR("\n%s\n","#*****************************************************************************")


#if (RSB_WANT_IO_LEVEL==0)
#define RSB_QUIET 1
#endif /* RSB_WANT_IO_LEVEL */

/* RSB_INFO is the stream of informative messages which are user requested and expected (that is, not errors). */
#ifdef RSB_QUIET
#define RSB_INFO( ... ) RSB_NULL_COMMA_STATEMENT_FOR_ZEN_HAPPINESS 	
#else /* RSB_QUIET */
#define RSB_INFO( ... ) ((rsb_global_session_handle.out_stream)?fprintf(rsb_global_session_handle.out_stream, __VA_ARGS__ ):RSB_NULL_EXPRESSION_FOR_ZEN_HAPPINESS)
#endif /* RSB_QUIET */
/* RSB_FPRINTF is just a tool */
#define RSB_FPRINTF( ... ) fprintf( __VA_ARGS__ )

#if   defined(__GNUC__)
	/* GCC */
        #define RSB_UNLIKELY(expr) __builtin_expect(!!(expr),0)
        #define RSB_LIKELY(expr)   __builtin_expect(!!(expr),1)
        #define RSB_ALIGNED __attribute__((aligned (sizeof(double)*sizeof(unsigned char))))
/*        #define RSB_ALIGNED __attribute__((aligned (64)))	*/
#else /* __GNUC__ */
        #define RSB_UNLIKELY(expr)  (expr)
        #define RSB_LIKELY(expr)   (expr)
        #define RSB_ALIGNED
#endif /* __GNUC__ */
#define RSB_PERFORMANCE_BINARY_DUMP_FILE_SIGNATURE_MAX_CHARS 128
#define RSB_PERFORMANCE_BINARY_DUMP_FILE_SIGNATURE \
"this is a non portable performance dump file, dude........\x40\x40\x40\x40"

#define RSB_TIMER_GRANULARITY_TEST_TIMES (10000) /* call timer RSB_TIMER_GRANULARITY_TEST_TIMES times */
#define RSB_TIMER_SANITY_TEST_TIMES (1024)
#define RSB_MIN_ALLOWED_CACHE_BLOCK_SIZE (1024)	/* in bytes */
#define RSB_MAX_ALLOWED_CACHE_BLOCK_SIZE ((1024)*(1024)*(1024))	/* in bytes */

#define RSB_BZERO(b,len) (memset((b), '\0', (len)), (void) 0) /* recommendation from IEEE Std 1003.1 since bzero has been made legacy */

#define RSB_BZERO_P(P) RSB_BZERO(P,sizeof(*(P)))

#define RSB_MEMMOVE memmove	/**< we have the chance of using a custom memmove function in this way */

#define RSB_MEMCMP memcmp	/**< we have the chance of using a custom memcmp function in this way */

int rsb__getopt_long( int argc, char * const argv[], const char *optstring, const rsb_option_t *longopts, int *longindex);
#define rsb__numerical_memcpy(TYPECODE,DST,DOFF,SRC,SOFF,N) {size_t es = RSB_NUMERICAL_TYPE_SIZE(TYPECODE); 	\
	rsb__memcpy( 					\
			((rsb_byte_t*)(DST)+es*(DOFF)) ,	\
			((const rsb_byte_t*)(SRC)+es*(SOFF)) ,	\
			es*(N));	\
} 		/* see rsb__xcopy(DST,SRC,DOFF,SOFF,N,size_t el_size) */
void *rsb__memcpy(void *RSB_RESTRICT dest, const void *RSB_RESTRICT src, size_t n);

rsb_time_t rsb__timer_sanity(void);
size_t rsb__sys_free_system_memory(void);
size_t rsb__sys_total_system_memory(void);
long rsb__get_lastlevel_c_size_per_thread(void);
long rsb__get_cache_block_byte_size(void);
/* long rsb_set_executing_threads(long tn); */
rsb_err_t rsb__print_memory_allocation_info(void);
rsb_err_t rsb__printf_memory_allocation_info(FILE *os);
rsb_err_t rsb__lock_as_memory_resident(rsb_bool_t dolock);
int rsb__fileno(FILE *stream);
#if defined(RSB_HAVE_EXECINFO_H) && RSB_OUT_ERR_VERBOSITY>=2
void rsb__print_trace (void);
#else  /* RSB_HAVE_EXECINFO_H */
#define rsb__print_trace RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
#endif /* RSB_HAVE_EXECINFO_H */
rsb_err_t rsb__getrusage(void);
const rsb_char_t * rsb__getenv(const rsb_char_t * name);
const rsb_char_t * rsb__getenv_nnr(const rsb_char_t * name);
#if RSB_WITH_HWLOC
long rsb__get_lnc_size_hwloc(int n);
#endif	/* RSB_WITH_HWLOC */
#define RSB_MAX_HOSTNAME_LEN 64
void rsb__strcpy_hostname(rsb_char_t * buf);

#define RSB_THREADS_GET_MAX_LIB	-5
#define RSB_THREADS_GET_MAX_SYS	-4
#define RSB_THREADS_GET_MAX	-3
#define RSB_THREADS_GET		-2
#define RSB_THREADS_AUTO	-1
rsb_int_t rsb__set_num_threads(rsb_int_t tn);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_SYS_H_INCLUDED */

/* @endcond */
