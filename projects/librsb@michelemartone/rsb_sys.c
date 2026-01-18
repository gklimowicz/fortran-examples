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
/* @cond INNERDOC  */
/*!
 * @file
 * @author Michele Martone
 * @brief System, or standard library related functions.
 * */

#include "rsb_common.h"
#include <unistd.h>	/* sysconf */
#include "rsb_internals.h"
#include "rsb.h"
#ifdef RSB_HAVE_LIMITS_H
#include <limits.h>	/* CHAR_BIT */
#endif /* RSB_HAVE_LIMITS_H */
#include <assert.h>	/* assert */
#ifdef RSB_HAVE_STDLIB_H
#include <stdlib.h>	/* posix_memalign */
#endif /* RSB_HAVE_STDLIB_H */
#ifdef RSB_HAVE_MALLOC_H
#include <malloc.h>	/* posix_memalign */
#endif /* RSB_HAVE_MALLOC_H */
#ifdef RSB_HAVE_SYS_SYSTEMCFG_H 
#include <sys/systemcfg.h>	/* for _H_SYSTEMCFG */
#endif /* RSB_HAVE_SYS_SYSTEMCFG_H  */
#ifdef RSB_HAVE_SYS_MMAN_H
#if RSB_HAVE_SYS_MMAN_H
#include <sys/mman.h>	/* for mlockall */
#endif /* RSB_HAVE_SYS_MMAN_H */
#endif /* RSB_HAVE_SYS_MMAN_H */
#if defined(RSB_HAVE_DMALLOC_H) && defined(RSB_WANT_DMALLOC) && (RSB_WANT_DMALLOC!=0)
#include <dmalloc.h>	/* a debug library */
#endif /* defined(RSB_HAVE_DMALLOC_H) && defined(RSB_WANT_DMALLOC) && (RSB_WANT_DMALLOC!=0) */
#if defined(RSB_HAVE_DUMA_H) && defined(RSB_WANT_DUMA) && (RSB_WANT_DUMA!=0)
#include <duma.h>	/* a debug library */
#endif /* defined(RSB_HAVE_DUMA_H) && defined(RSB_WANT_DUMA) && (RSB_WANT_DUMA!=0) */
/* #include <stdio.h> */	/* fileno */
#if RSB_WITH_HWLOC
#include <hwloc.h>
#endif	/* RSB_WITH_HWLOC */
#ifdef RSB_HAVE_EXECINFO_H
#include <execinfo.h>
#endif /* RSB_HAVE_EXECINFO_H */
#include <stdint.h> /* int64_t / uint64_t / uintptr_t */

/* set the following to 0 to get some real fun */
#define RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION 0
#if RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION
/* FIXME: need a true random number generator interface */
#define RSB_SHOULD_RANDOMLY_FAIL (rsb__rand_coo_index(1349)==42)
#else /* RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION */
#define RSB_SHOULD_RANDOMLY_FAIL 0
#endif /* RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION */

#define RSB_SUSPECT_ALLOCATION_SIZE 1024*1024*16 /* i.e. of notable size. */
#define RSB_MEM_DEBUG 0 /* verbosity */
#define RSB_MEM_DEBUG_REALLOC 0 /* verbosity */
#define RSB_CHEAP_DEBUG 0 /* cheap extra checks on suspect memory wrapper related values */
#define RSB_DEBUG_MARKER_AFTER_FREE RSB_MEM_DBG /* double free protection */
#define RSB_DEBUG_SHRED_AFTER_FREE RSB_MEM_DBG /* protection against of re-use of freed areas */
#define RSB_SHRED_BYTE /* 0xFF */ 0xF0 /* byte value to use when shredding memory */
#define RSB_SHRED_WORD ( RSB_SHRED_BYTE | ( RSB_SHRED_BYTE<<8 ) | ( RSB_SHRED_BYTE<<16 ) | ( RSB_SHRED_BYTE<<24 ) )
#define RSB_FREE_MARKER 0xDEADBEEF /* ( 3735928559) */ /* a marker for detecting double free */
#define RSB_OVW_MARKER  0xBEEFBABE /* (-1091585346) */ /* a marker for detecting accidental overwrites */
#if RSB_MEM_DBG 
#define RSB_MW_ODMO (-3) /* memory wrapper overwrite detection marker offset (set 0 to deactivate, set to a non-*MW* overlapping value to activate) */
#else /* RSB_MEM_DBG */
#define RSB_MW_ODMO ( 0) /* memory wrapper overwrite detection marker offset (set 0 to deactivate, set to a non-*MW* overlapping value to activate) */
#endif /* RSB_MEM_DBG */
#define RSB_MW_SHMO (-1) /* memory wrapper shift marker offset */
#define RSB_MW_SZMO (-2) /* memory wrapper size  marker offset */
#define RSB_MW_ESLC (-RSB_MIN(RSB_MIN(RSB_MW_ODMO,RSB_MW_ODMO),RSB_MIN(RSB_MW_SZMO,RSB_MW_SHMO))) /* memory wrapper extra 'sizeof' locations count */

#if RSB_DEBUG_SHRED_AFTER_FREE
/* prevent previously shredded and freed, and now re-allocated *data to cause false alarm if unused before next free */
#define RSB_UNSHRED_FRESH(P,SIZE) if((P) && (SIZE)>0) memset((P), ~RSB_SHRED_BYTE, (SIZE));
#else
#define RSB_UNSHRED_FRESH(P,SIZE)
#endif /* RSB_DEBUG_SHRED_AFTER_FREE */

#define RSB_MD_ASSERT RSB_ASSERT  /* memory debug assert macro */

#if RSB_DEBUG_SHRED_AFTER_FREE
#include <string.h>	/* memset */
#endif

#if RSB_USE_GETRUSAGE
#ifdef RSB_HAVE_SYS_RESOURCE_H
#include <sys/resource.h>	/* getrusage */
#endif
#ifdef RSB_HAVE_SYS_TIME_H
#include <sys/time.h>	/* getrusage */
#endif
#endif

#ifdef RSB_HAVE_SYS_UTSNAME_H 
#include <sys/utsname.h>	/* uname */
#endif /* RSB_HAVE_SYS_UTSNAME_H  */

RSB_INTERNALS_COMMON_HEAD_DECLS

#if RSB_WANT_ALLOCATOR_LIMITS
#define RSB_ALLOC_MEMAAA_LIMIT rsb_global_session_handle.memory_count_max
#define RSB_ALLOC_MEMAAC_LIMIT rsb_global_session_handle.allocations_count_max
#define RSB_ALLOC_LIMITS_TRESPASSED(AAA,AAC) RSB_UNLIKELY( ( ( RSB_ALLOC_MEMAAA_LIMIT > 0 && rsb_global_session_handle.allocated_memory+(AAA) >= RSB_ALLOC_MEMAAA_LIMIT ) || ( RSB_ALLOC_MEMAAC_LIMIT > 0 && rsb_global_session_handle.allocations_count+(AAC) > RSB_ALLOC_MEMAAC_LIMIT ) ) ? 1 : 0 )
#else /* RSB_WANT_ALLOCATOR_LIMITS */
#define RSB_ALLOC_LIMITS_TRESPASSED(AAA,AAC) 0
#endif /* RSB_WANT_ALLOCATOR_LIMITS */
#define RSB_ZERO_BYTE_ALLOC_CHECK	0

#define RSB_TIME_SET_THREADS 0

void rsb__g_rsb_memory_counter_init(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * Memory counter reset.
	 */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	rsb_global_session_handle.allocated_memory=0;
	rsb_global_session_handle.allocations_count=0;
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}

size_t rsb__get_g_rsb_memory_count(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * A mere accessor function.
	 */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	return rsb_global_session_handle.allocated_memory;
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	return 0;
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}

size_t rsb__get_g_rsb_allocations_count(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * A mere accessor function.
	 */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	return rsb_global_session_handle.allocations_count;
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	return 0;
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}

#if 1

#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
static void * rsb_aligned_free(void *p)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \param p a generic pointer allocated with rsb__aligned_malloc
	 * \return the input pointer, in case of correct operation, NULL in case of error
	 *
	 * frees a memory area previously allocated with rsb__aligned_malloc,
	 * and returns the pointer in case of successfull free,
	 * or NULL in case of suspect error.
	 * */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	size_t size;
	size_t shift;
	if( p == NULL )
		return p;
	#if RSB_DEBUG_SHRED_AFTER_FREE
	if( (((uintptr_t)p) & RSB_SHRED_WORD ) == RSB_SHRED_WORD ) /* shred-area read pointer detection */
	{
		RSB_STDERR("Warning: it is likely that pointer %p is invalid and was read from a previously freed area. Expect a crash now.\n",p);
		RSB_MD_ASSERT(0);
	}
	#endif /* RSB_DEBUG_SHRED_AFTER_FREE */
	size  = ((size_t*)p)[RSB_MW_SZMO];
	#if RSB_MEM_DEBUG
	RSB_STDERR("freeing   %zu bytes at %p (in hex:0x%0zx bytes)\n",size,p,size);
	#endif
	shift = ((size_t*)p)[RSB_MW_SHMO];
	#if RSB_DEBUG_MARKER_AFTER_FREE
	if(size == RSB_FREE_MARKER || shift == RSB_FREE_MARKER)
	{
		RSB_STDERR("Warning: it is almost certain that memory at %p has been already deallocated! Expect a crash now.\n",p);
		RSB_MD_ASSERT(0);
	}
	((size_t*)p)[RSB_MW_SZMO] = RSB_FREE_MARKER;
	((size_t*)p)[RSB_MW_SHMO] = RSB_FREE_MARKER;
	#endif /* RSB_DEBUG_MARKER_AFTER_FREE */
	if( RSB_MW_ODMO && ((size_t*)p)[RSB_MW_ODMO] != RSB_OVW_MARKER ) 
	{
		RSB_STDERR("Warning: memory at %p has been overwritten (marker value is %0llx instead of 0x%0llx) ! Expect crashes.\n",p,(long long unsigned)((size_t*)p)[RSB_MW_ODMO],(long long unsigned)RSB_OVW_MARKER );
		RSB_MD_ASSERT(0);
	}
	#if RSB_DEBUG_SHRED_AFTER_FREE
	if( ( size >= sizeof(rsb_int_t) ) && ( ((uintptr_t)p) == RSB_SHRED_WORD ) ) /* shredded area re-free detection */
	{
		RSB_STDERR("Warning: it is possible that %zd-byte area at %p was recently freed (points to {0x%x, ...). May expect a crash now.\n",size,p,*(rsb_int_t*)p);
		/* RSB_MD_ASSERT(0); */
	}
	memset(p, RSB_SHRED_BYTE, size);
	#endif /* RSB_DEBUG_SHRED_AFTER_FREE */
	p = (( char *)p)-(shift);
	p = ((size_t*)p)-RSB_MW_ESLC;	/* we make room for the markers */
#pragma omp atomic
	rsb_global_session_handle.allocated_memory -= size;
#pragma omp atomic
	rsb_global_session_handle.allocations_count--;
	free(p);
	return p;
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	free(p);
	return p;
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */

void * rsb__realloc(void *rsb_data, size_t size)
{
#ifdef RSB_WANT_DOUBLE_ALIGNED
	return rsb__do_realloc(rsb_data,size,sizeof(double)*2);
#else /* RSB_WANT_DOUBLE_ALIGNED */
	return rsb__do_realloc(rsb_data,size,1);
#endif /* RSB_WANT_DOUBLE_ALIGNED */
}

void * rsb__do_realloc(void *rsb_data, size_t size, size_t alignment)
{
	void * p = rsb_data;
	size_t extra = 0;
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	size_t osize;
	size_t shift;
	extra = sizeof(size_t)*RSB_MW_ESLC+alignment;

	if(p==NULL)
		return p;
	osize  = ((size_t*)p)[RSB_MW_SZMO];
	#if ( RSB_MEM_DEBUG || RSB_MEM_DEBUG_REALLOC )
	RSB_STDERR("reallocating from %zu to %zu bytes (%+zd) at 0x%p (in hex: to 0x%0zx bytes)\n",osize,size,(size-osize),p,size);
	#endif
	shift = ((size_t*)p)[RSB_MW_SHMO];
	p = (( char *)p)-(shift);
	p = ((size_t*)p)-RSB_MW_ESLC;	/* we make room for the markers */
	/* no free was performed, since extra>0 */
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	if(size==0) /* a free shall be performed */
		extra=0;

#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	if(size==0 && p) /* a free shall be performed */
		rsb_global_session_handle.allocations_count--;

	if(size>osize)
		rsb_global_session_handle.allocated_memory+=size-osize;
	else
	if(p) /* a free shall be performed */
		rsb_global_session_handle.allocated_memory-=osize-size;
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	p = realloc(p,size+extra);/* if freeing, either p or NULL will be returned */

	if(!p)
		return p;/* failure */

#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	/*!
	 * \ingroup gr_internals
	 * TODO : no way to make survive the alignment ?
	 * */
#if 1
	/* restoring back allocation info (and potentially losing alignment, because we have to keep the same shift!) */
	p = ((size_t*)p)+RSB_MW_ESLC;	/* we make room for markers */
	p = (( char *)p)+(shift);
	/* to restore alignment, should perform a memmove */
	((size_t*)p)[RSB_MW_SHMO] = shift;
	((size_t*)p)[RSB_MW_SZMO] = size;
	if( RSB_MW_ODMO ) ((size_t*)p)[RSB_MW_ODMO] = RSB_OVW_MARKER; 
#else
	/* bugful way */

	{
		size_t off;
		off=((size_t)(((size_t*)p)+2))%(alignment); /* to the return address from ((size_t*)p)+2 */
		shift = (alignment-off);
		p = ((size_t*)p)+2;	/* we make room for two markers */
		p = (( char *)p)+(shift);
		((size_t*)p)[RSB_MW_SHMO] = shift;
		((size_t*)p)[RSB_MW_SZMO] = size;
		if( RSB_MW_ODMO ) ((size_t*)p)[RSB_MW_ODMO] = RSB_OVW_MARKER; 
	}
#endif
	return p;
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	return p; /* success */
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	return NULL;
}

void * rsb__aligned_malloc(size_t size, size_t alignment)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * allocates size bytes and an integer, in a way to keep track of the allocated chunk size
	 * \param size is the amount of needed bytes to allocate
	 *
	 * the returned area will be alignment bytes aligned.
	 * this area should be deallocated with rsb_aligned_free.
	 *
	 * note that although the address will be aligned as asked, 
	 * there will be 100% alignment and contiguity guarantee only
	 * on machines with no virtual memory (on linux there seems not to be a user space contiguos
	 * memory allocator like kernel's vmalloc).
	 *
	 * in case of lack of all of RSB_DISABLE_ALLOCATOR_WRAPPER, RSB_HAVE_POSIX_MEMALIGN and RSB_HAVE_MEMALIGN, malloc() will be used.
	 * */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	void * p;
	const size_t extra = sizeof(size_t)*RSB_MW_ESLC+alignment;
	size_t off;
	size_t shift;

	if(RSB_ALLOC_LIMITS_TRESPASSED(size+extra,1))
	{
		rsb__print_memory_allocation_info(); 
		return NULL;
	}

#if RSB_ZERO_BYTE_ALLOC_CHECK
	if(size == 0 && extra == 0)
	{
		RSB_ERROR(RSB_ERRM_ZSM);
	}
#endif
	p = malloc( size + extra );
	if(!p)
		return p;/* failure */
	RSB_DEBUG_ASSERT(rsb_global_session_handle.allocated_memory	<=(rsb_global_session_handle.allocated_memory+size  ));
	/* the following could trigger during very very long/big runs, so .. */
	RSB_DEBUG_ASSERT((rsb_global_session_handle.allocations_count	< (rsb_global_session_handle.allocations_count+1))
			/* ... this line should fix that cases */
			|| (rsb_global_session_handle.allocations_count+1)==0);

#pragma omp atomic
	rsb_global_session_handle.allocated_memory+=size;
#pragma omp atomic
	rsb_global_session_handle.allocations_count++;
#pragma omp atomic
	rsb_global_session_handle.allocations_cumulative++;
	off=((size_t)(((size_t*)p)+RSB_MW_ESLC))%(alignment); /* to the return address from ((size_t*)p)+RSB_MW_ESLC */
	shift = (alignment-off);
	p = ((size_t*)p)+RSB_MW_ESLC;	/* we make room for the markers */
	p = (( char *)p)+(shift);
	((size_t*)p)[RSB_MW_SHMO] = shift;
	((size_t*)p)[RSB_MW_SZMO] = size;
	if( RSB_MW_ODMO )
		((size_t*)p)[RSB_MW_ODMO] = RSB_OVW_MARKER; 
	RSB_UNSHRED_FRESH(p,size);
	return p;
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	void * p = NULL;
	#if RSB_HAVE_POSIX_MEMALIGN 
        size_t ca = sizeof(void*); /* corrected alignment */
        while(ca<alignment)
                ca*=2;
        alignment = ca; /* "The address  of  the  allocated  memory  will be a multiple of alignment, which must be a power of two and a multiple of sizeof(void *)." */                                                                                      
	if(posix_memalign(&p, alignment, size))
		;/* failure dealt elseqhere */
	else
		; /* success */
	#elif RSB_HAVE_MEMALIGN 
	p = memalign( alignment, size);
	#else
	p = malloc(size); /* no platform support for aligned alloc */
	#endif
	RSB_UNSHRED_FRESH(p,size);
	return p;
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}

void * rsb__free(void *p)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * see rsb_aligned_free
	 * */
	#if RSB_MEM_DEBUG
	RSB_STDERR("freeing %p\n",p);
	#endif
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	return rsb_aligned_free(p);
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	free(p);
	return p;
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}

void * rsb__malloc(size_t size)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * see rsb__aligned_malloc
	 * */
	void * data=NULL;

#if RSB_ZERO_BYTE_ALLOC_CHECK
	if(size == 0)
	{
		RSB_ERROR(RSB_ERRM_ZSM);
	}
#endif
	if(size >= RSB_MAX_ALLOCATABLE_MEMORY_CHUNK)
	{	
		#if RSB_MEM_DEBUG
		RSB_STDERR("cannot allocate %zu bytes since it is more than the maximum allowed %zu\n",size,RSB_MAX_ALLOCATABLE_MEMORY_CHUNK);
		#endif
		return data;
	}
#ifdef RSB_WANT_DOUBLE_ALIGNED
	data = rsb__aligned_malloc(size,sizeof(double)*2);
#else /* RSB_WANT_DOUBLE_ALIGNED */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	data = rsb__aligned_malloc(size,1);
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	data = malloc(size);
	RSB_UNSHRED_FRESH(data,size);
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
#endif /* RSB_WANT_DOUBLE_ALIGNED */
	#if RSB_MEM_DEBUG
	RSB_STDERR("allocated %zu bytes to %p (in hex:0x%0zx bytes)\n",size,data,size);
	#endif /* RSB_MEM_DEBUG */

#if RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION
	if(data)if(RSB_SHOULD_RANDOMLY_FAIL){rsb__free(data);data=NULL;}
#endif /* RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION */
	return data;
}

#else
/* BEGIN OF DEAD CODE */
void * rsb__free(void *rsb_data)
{
	/*!
	 * TODO: DELETE THIS DEAD CODE
	 *
	 * \param rsb_data a generic pointer allocated with rsb__malloc
	 * \return the input pointer, in case of correct operation, NULL in case of error
	 *
	 * deletes a memory area previously allocated with rsb__malloc,
	 * and returns the pointer in case of successfull free,
	 * or NULL in case of suspect error.
	 * */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	size_t size;
	if(!rsb_data)return rsb_data;

#ifdef RSB_WANT_DOUBLE_ALIGNED
	/*!
	 * This is a trick : we take the offset written in the byte just 
	 * behind the allocated area and use it to restore the original
	 * [m|c]alloc-ated area.
	 * */
	size_t off = ((unsigned rsb_byte_t*)rsb_data)[RSB_MW_SHMO];
	/* we decode back the whole allocated area size */
	size=((size_t*)((unsigned rsb_byte_t*)(rsb_data)-(0x10-off)))[RSB_MW_SHMO];
	/*
	RSB_STDERR("freeing address after offset %d ... \n",off);
	RSB_STDERR("freeing  %d bytes ... \n",size);
	*/
#else /* RSB_WANT_DOUBLE_ALIGNED */
	/* we decode back the whole allocated area size */
	size=*(((size_t*)(rsb_data))RSB_MW_SHMO);
#endif /* RSB_WANT_DOUBLE_ALIGNED */

	#if RSB_CHEAP_DEBUG
	if( size < 0 || size > RSB_SUSPECT_ALLOCATION_SIZE )/* this is a BAD sign and a warning should be issued */
	{
		RSB_ERROR("WARNING : pointer x%08x[%d] contains (has been overwritten with ?) a suspect value : %08x==%d", rsb_data,RSB_MW_SHMO,size,size );
	gaa
		return NULL;
	}
	#endif /* RSB_CHEAP_DEBUG */
	/* We update the global memory bookkeeping counter */
#pragma omp atomic
	rsb_global_session_handle.allocated_memory-=size;
#pragma omp atomic
	rsb_global_session_handle.allocations_count--;
#ifdef RSB_WANT_DOUBLE_ALIGNED
	/* We use the offset to restore the original pointer to free. */
	free((  size_t*)((unsigned rsb_byte_t*)(rsb_data)-(0x10-off))RSB_MW_SHMO);
#else /* RSB_WANT_DOUBLE_ALIGNED */
	free((((size_t*)(rsb_data))RSB_MW_SHMO));
#endif /* RSB_WANT_DOUBLE_ALIGNED */
	return rsb_data;
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	/* RSB_DISABLE_ALLOCATOR_WRAPPER undefined */
#ifdef RSB_WANT_DOUBLE_ALIGNED
	return free(rsb_data);
#else /* RSB_WANT_DOUBLE_ALIGNED */
	return free(rsb_data);
#endif /* RSB_WANT_DOUBLE_ALIGNED */
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}

void * rsb__malloc(size_t size)
{
	/*!
	 * TODO: DELETE THIS DEAD CODE
	 * 
	 * (c)allocates size bytes and an integer, in a way to keep track of the allocated chunk size
	 * \param size is the amount of needed bytes to allocate
	 *
	 * if RSB_WANT_DOUBLE_ALIGNED is defined, the returned area will be double (64 bits) aligned.
	 * this area should be deallocated with rsb__free.
	 * */
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
	void * p;
#ifdef RSB_HAVE_POSIX_MEMALIGN
	 /* we could integrate/replace with posix_memalign, memalign or continue using our custom code */
#endif /* RSB_HAVE_POSIX_MEMALIGN */
#ifdef RSB_WANT_DOUBLE_ALIGNED
	/*
	 * This is an explicit trick to give the user a double aligned memory area.
	 * To achieve this, we allocate one extra double element, and a byte (four, really, for congruency reasons)
	 * and write in one  there information on the shift amount.
	 *
	 * Of course; we count that sizeof(size_t)>sizeof(char).
	 * */
	size_t extra = sizeof(double)+sizeof(size_t)*2;
#else /* RSB_WANT_DOUBLE_ALIGNED */
	/*
	 * We allocate one size_t element for storing allocation information (for explicit memory leaking checking).
	 * */
	size_t extra = sizeof(size_t);
#endif /* RSB_WANT_DOUBLE_ALIGNED */
	if(RSB_ALLOC_LIMITS_TRESPASSED(size,1))
	{rsb__print_memory_allocation_info(); return NULL;}
	p = calloc( size + extra, 1 );
	if(!p)return p;
	*(size_t*)p=size;/* note : not size + extra */
#pragma omp atomic
	rsb_global_session_handle.allocated_memory+=size;
#pragma omp atomic
	rsb_global_session_handle.allocations_count++;
#pragma omp atomic
	rsb_global_session_handle.allocations_cumulative++;
#ifdef RSB_WANT_DOUBLE_ALIGNED
	/*
	 * WARNING : We determine the current 64 bits p alignment, so we are interested in the last 5 bits,
	 * really.
	 * DANGER  : is this portable ?
	 * */
	size_t off=(((size_t)p)+sizeof(size_t))&0xF;	/* can be 0 ... F */
	/*
	RSB_STDERR("allocated totally %d bytes \n",size+extra);
	RSB_STDERR("allocated %d ... \n",size);
	RSB_STDERR("allocation offset %d ... \n",off);
	*/
	((unsigned rsb_byte_t*)(p))[sizeof(size_t)+((0x10-off)RSB_MW_SHMO)]=(unsigned char)off;/* will be used to compute back the base pointer allocated */
	return ((unsigned rsb_byte_t*)p)+(0x10-off)+sizeof(size_t);
#else /* RSB_WANT_DOUBLE_ALIGNED */
	return ((size_t*)p)+1;
#endif /* RSB_WANT_DOUBLE_ALIGNED */
#else /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	/* RSB_DISABLE_ALLOCATOR_WRAPPER undefined */
#ifdef RSB_WANT_DOUBLE_ALIGNED
	return calloc(size,1);
#else /* RSB_WANT_DOUBLE_ALIGNED */
	return calloc(size,1);
#endif /* RSB_WANT_DOUBLE_ALIGNED */
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
}
/* END OF DEAD CODE */
#endif /* 1 */

rsb_time_t rsb__do_time(void)
{
	/*!
	   \ingroup gr_internals

	   Returns a current relative time in seconds.
	   The user should rely on this function only for time difference computations.
	 */
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	/* return omp_get_wtime(); */ /* for future use */
#endif
	/* SVr4, 4.3BSD.  POSIX.1-2001 */
	/* ( could also use psb_wtime or mpi_wtime )*/
	/* FIXME : gettimeofday() gives pessimistic estimates ! */
	/* FIXME : gettimeofday() could be in time.h or sys/times.h */
	/* FIXME : timer sanity is of paramount inportance ! Should check for its sanity at startup! */
#if defined(RSB_HAVE_GETTIMEOFDAY)
	register double t = RSB_REAL_ZERO;
	struct timeval tv1;
	gettimeofday(&tv1, NULL);
	t  =  (double)(tv1.tv_sec) + ((double)(tv1.tv_usec))*1.e-6;
	return t;
#elif defined(RSB_HAVE_TIMES) && defined(RSB_HAVE_SYSCONF) && defined(_SC_CLK_TCK)
	/* POSIX.1 */
	struct tms buffer;
	times(&buffer);
	return ( (rsb_time_t) ((clock_t)buffer.tms_utime) ) / ( (rsb_time_t)sysconf(_SC_CLK_TCK) );
#else /* defined(RSB_HAVE_TIMES) && defined(RSB_HAVE_SYSCONF) && defined(_SC_CLK_TCK) */
#error("You should better find timing routine, dude.\n")
	return -1;/* this is bad */
#endif /* defined(RSB_HAVE_TIMES) && defined(RSB_HAVE_SYSCONF) && defined(_SC_CLK_TCK) */
}

rsb_time_t rsb__timer_sanity(void)
{
	/*!
		\ingroup gr_internals
		
		Could the timer lead to negative time intervals ? (we found it can happen, sadly)
		\return the minimum interval length after a bunch of timer calls

		TODO : could we do something about this ?
	*/
	rsb_time_t md,d;
	int i;

	md = - rsb_time();
	md += rsb_time();
	for(i=0;i<RSB_TIMER_SANITY_TEST_TIMES;++i)
	{
		d = - rsb_time();
		d += rsb_time();
		md=d<md?d:md;
	}
	return md;
}

rsb_time_t rsb__timer_granularity(void)
{
	/*!
	   \ingroup gr_internals

	 * Estimate the granularity of the timing function (that is, its overhead) by measuring it.
         * This value is important to set minimal benchmarking times.
	 * No guarantee of return if the timing function is broken.
	 * */
	register double t = RSB_TIME_ZERO, t0 = RSB_TIME_ZERO;
	const int times = RSB_TIMER_GRANULARITY_TEST_TIMES;
	register int i = times;

#if 0
	i = times;
	/* results in the following two code snippets differ. it should be due to numerical roundoff. */
	while(i--)
	{
		t -= rsb_time();
		/* no op, only call overhead */
		t += rsb_time();
	}
	return t/(times*2);
#else
	/* this is more accurate (in particular: slower) but could be optimized out without an accumulator cookie (FIXME) */
	t0 = rsb_time();
	--i;
	t = -t0;

	while(i--)
	{
		/* no op, only call overhead */
		rsb_time();
		rsb_time();
	}

	t += rsb_time();
	
	t = t/(times*2);

	if(t <= RSB_TIME_ZERO)
		goto so_fast;
	else
		goto ret;
so_fast: /* FIXME: No guarantee of return with a broken timing function. */
	while( ( t = rsb_time() ) <= t0)
		;
	t -= t0;
ret:
	return t;
#endif
}

rsb_err_t rsb__printf_memory_allocation_info(FILE *os)
{
#ifndef RSB_DISABLE_ALLOCATOR_WRAPPER
/*!
 \ingroup gr_internals

 * A global memory counter, used for debugging purposes.
 * */
	RSB_FPRINTF(os,"rsb_global_session_handle.allocated_memory       \t:%zu\n",(rsb_printf_int_t)rsb_global_session_handle.allocated_memory);
	RSB_FPRINTF(os,"rsb_global_session_handle.allocations_count  \t:%zu\n",(rsb_printf_int_t)rsb_global_session_handle.allocations_count);
	RSB_FPRINTF(os,"rsb_global_session_handle.allocations_cumulative \t:%zu\n",(rsb_printf_int_t)rsb_global_session_handle.allocations_cumulative);
	return RSB_ERR_NO_ERROR; 
#endif /* RSB_DISABLE_ALLOCATOR_WRAPPER */
	return RSB_ERR_UNSUPPORTED_FEATURE;
}

rsb_err_t rsb__print_memory_allocation_info(void)
{
	rsb_err_t errval = RSB_ERR_UNSUPPORTED_FEATURE;
#if RSB_ALLOW_STDOUT
	errval = rsb__printf_memory_allocation_info(stdout);
#endif /* RSB_ALLOW_STDOUT */
	return errval;
}

void * rsb__calloc(size_t n)
{
	/*!
	 * \ingroup gr_internals
	 * Allocates an amount of n bytes set to zero.
	 *
	 * \param n is the amount of bytes to allocates
	 * \return the newly allocated area.
	 *
	 * This memory area should be freed with rsb__free.
	 * */
	void * p = rsb__malloc(n);
	if(p)
		RSB_BZERO(p,n);
#if(!RSB_QUIET_MEM_ERRORS)
	/* TODO : message diagnostics should be optable out, or debug levels-based */
        else
	{
                RSB_ERROR("cannot allocate %zu bytes!\n",n);
		rsb__print_memory_allocation_info();
	}
#endif /* (!RSB_QUIET_MEM_ERRORS) */
	/* should be ((int*)p)[RSB_MW_SHMO]==n */
	return p;
}

void * rsb__calloc_parallel(size_t n)
{
	/*!
	 * \ingroup gr_internals
	 */
	void *p = rsb__calloc(n);

	if(p)
		RSB_BZERO_parallel(p,n);
	return p;
}

int rsb__error(const char * format, ...)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \param as the first printf argument.
	 * \return an error code or 0
	 *
	 * For now, a wrapper around printf.
	 * It will print given arguments on stdout.
	 * */
	va_list ap;
	int rc=0;

	va_start(ap,format);
#ifndef RSB_QUIET
	rc = vprintf(format,ap);
#endif /* RSB_QUIET */
	va_end(ap);
	return rc;
}

long rsb__get_lnc_size(int n)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \returns the nth level data cache size if > 0, -1 on error, 0 if no such level cache
	 * Cache levels start from 1.
	 * */
	long cs=0;

	if(rsb_global_session_handle.memory_hierarchy_levels>0)
		return rsb_global_session_handle.caches[n].size;

	switch(n)
	{
		case 1:
		{
#ifdef RSB_HAVE_SYSCONF 
#ifdef _SC_LEVEL1_DCACHE_SIZE
			cs=sysconf(_SC_LEVEL1_DCACHE_SIZE);
#endif /* _SC_LEVEL1_DCACHE_SIZE */
#endif /* RSB_HAVE_SYSCONF  */
#ifdef _H_SYSTEMCFG
			cs=_system_configuration.dcache_size;
#endif /* _H_SYSTEMCFG */
#if RSB_WITH_HWLOC
			if(cs == 0) cs = rsb__get_lnc_size_hwloc(n);
#endif	/* RSB_WITH_HWLOC */
		}
		break;
		case 2:
		{
#ifdef RSB_HAVE_SYSCONF 
#ifdef _SC_LEVEL2_CACHE_SIZE
			cs=sysconf(_SC_LEVEL2_CACHE_SIZE);
#endif /* _SC_LEVEL2_CACHE_SIZE */
#endif /* RSB_HAVE_SYSCONF */
#ifdef _H_SYSTEMCFG
			cs=_system_configuration.L2_cache_size;
#endif /* _H_SYSTEMCFG */
#if RSB_WITH_HWLOC
			if(cs == 0) cs = rsb__get_lnc_size_hwloc(n);
#endif	/* RSB_WITH_HWLOC */
		}
		break;
		case 3:
		{
#ifdef RSB_HAVE_SYSCONF 
#ifdef _SC_LEVEL3_CACHE_SIZE
			cs=sysconf(_SC_LEVEL3_CACHE_SIZE);
#endif /* _SC_LEVEL3_CACHE_SIZE */
#endif /* RSB_HAVE_SYSCONF */
#ifdef _H_SYSTEMCFG
	//		cs=_system_configuration.L3_cache_size; // Does not exist :(
#endif /* _H_SYSTEMCFG */
#if RSB_WITH_HWLOC
			if(cs == 0) cs = rsb__get_lnc_size_hwloc(n);
#endif	/* RSB_WITH_HWLOC */
		}
		break;
		default :
		/* For now, we don't handle more cache levels */
		cs=-1;
	}
	cs=cs<0?0:cs;
	return cs;
}

long rsb__get_l1c_size(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the first level data cache size.
	 * \note see rsb__get_lnc_size
	 * */
	return rsb__get_lnc_size(1);
}

#if RSB_WANT_EXPERIMENTS_CODE
long rsb__get_l2c_size(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the second level data cache size.
	 * \note see rsb__get_lnc_size
	 * */
	return rsb__get_lnc_size(2);
}
#endif /* RSB_WANT_EXPERIMENTS_CODE */

#if RSB_OBSOLETE_QUARANTINE_UNUSED
long rsb__get_l3c_size(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the third level data cache size.
	 * \note see rsb__get_lnc_size
	 * */
	return rsb__get_lnc_size(3);
}

long rsb__get_l4c_size(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the fourth level data cache size.
	 * \note see rsb__get_lnc_size
	 * */
	return rsb__get_lnc_size(4);
}

long rsb__know_cache_sizes(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return do we know the caches sizes ?
	 * */
	return rsb__get_cache_levels_num() > 0;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

long rsb__get_first_level_c_size(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the first level data cache size
	 * or zero if not known or not available.
	 * */
	rsb_int_t cln = rsb__get_cache_levels_num();

	if(rsb_global_session_handle.memory_hierarchy_levels>0)
		return rsb_global_session_handle.caches[1].size;

	if(cln>0)
		return rsb__get_lnc_size(1);
	else
		return 0;
}

long rsb__get_lastlevel_c_size(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the last level data cache size
	 * or zero if not known or not available.
	 * */
	rsb_int_t cln = rsb__get_cache_levels_num();

	if(rsb_global_session_handle.memory_hierarchy_levels>0)
		return rsb_global_session_handle.caches[rsb_global_session_handle.memory_hierarchy_levels].size;

	if(cln>0)
		return rsb__get_lnc_size(cln);
	else
		return 0;
}

long rsb__get_cache_block_byte_size(void)
{
	/*!
	 * \ingroup gr_internals
	 * */
	long flc = RSB_MIN(RSB_MAX(1,rsb_global_session_handle.memory_hierarchy_levels),2);
	long cbs = RSB_MIN(rsb__get_lnc_size(flc),rsb__get_lastlevel_c_size_per_thread());

	switch(rsb_global_session_handle.cache_blocking_method)
	{
		case -1: return cbs/2; break;
		case  1: return cbs*2; break;
		default:
		case  0: return cbs;
	}
}

static long rsb_want_executing_threads(void)
{
	/*!
	  	\ingroup gr_internals
		Will always return a value 1 <= N <= RSB_CONST_MAX_SUPPORTED_CORES
		FIXME: make so that 
		rsb_want_executing_threads() == rsb_set_executing_threads()
		or rather write
	       	rsb__set_num_threads(RSB_THREADS_GET)
	*/
	long wt = 1;

#if RSB_WANT_OMP_RECURSIVE_KERNELS
	wt = rsb_global_session_handle.rsb_want_threads;
	if(wt<RSB_CONST_MIN_SUPPORTED_CORES)
		return RSB_CONST_MIN_SUPPORTED_CORES;
	wt = RSB_MIN(wt,RSB_CONST_MAX_SUPPORTED_CORES);
#endif
	return wt;
}

long rsb__get_lastlevel_c_size_per_thread(void)
{
	/*!
	 * \ingroup gr_internals
	 *
	 * \return the last level data cache size divided by the number of active threads
	 * or zero if not known or not available.
	 * */
	return (rsb__get_lastlevel_c_size()/rsb_want_executing_threads());
//	return rsb__get_lastlevel_c_size()/sqrt(rsb_want_executing_threads());
}

rsb_int_t rsb__get_cache_levels_num(void)
{
	/*!
	 \ingroup internals
	 \return the count of cache levels if >0, -1 if unknown, 0 if no caches
	*/
	long cs,l=1;

	for(l=1;l<RSB_MAX_SUPPORTED_CACHE_LEVELS;++l)
	{
		cs = rsb__get_lnc_size(l);
		if(!cs){--l;break;}
	}
	return l;
}

int rsb__getopt_long(
	int argc, char * const argv[], const char *optstring,
	const rsb_option_t *longopts, int *longindex)
{
	/*!
	   \ingroup internals
	  
	   A compatibility wrapper.
	 */
#ifdef RSB_HAVE_GETOPT_LONG
	return getopt_long(argc,argv,optstring,longopts,longindex);
#else /* RSB_HAVE_GETOPT_LONG */
#ifdef RSB_HAVE_GETOPT
	return getopt(argc,argv,optstring);	/* a remedy */
#else /* RSB_HAVE_GETOPT */
#error  Neither of getopt() nor getopt_long() detected!
#endif /* RSB_HAVE_GETOPT */
#endif /* RSB_HAVE_GETOPT_LONG */
}

rsb_err_t rsb__sys_init(void)
{
	/*!
	 \ingroup internals

	 should check some system related init stuff,
	 to prevent nasty errors during execution.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(sizeof(rsb_err_t)<4)
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);

	if(sizeof(rsb_flags_t)<4)
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);

	/*! \todo : we could check for 'base' overflow cases to define 
	 	    some limit cases..
	 */
	if(sizeof(rsb_coo_idx_t)>sizeof(rsb_nnz_idx_t))
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);

	if(sizeof(rsb_blk_idx_t)>sizeof(rsb_coo_idx_t))
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);

	if(rsb__get_l1c_size()<=0)
		rsb_global_session_handle.min_leaf_matrix_bytes = RSB_EXPERIMENTAL_MIN_LEAF_ELEMENTS*sizeof(double);
	else
		rsb_global_session_handle.min_leaf_matrix_bytes = rsb__get_l1c_size();

#ifdef CHAR_BIT	/* limits.h */
	if( RSB_CHAR_BIT != CHAR_BIT )
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_INTERNAL_ERROR);
#endif /* CHAR_BIT */

	if(rsb__get_lastlevel_c_size()<0)
		rsb_global_session_handle.avg_leaf_matrix_bytes=4*RSB_EXPERIMENTAL_MIN_LEAF_ELEMENTS*sizeof(double);
	else
		rsb_global_session_handle.avg_leaf_matrix_bytes = rsb__get_lastlevel_c_size()*2;
#if RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION
	{
		rsb_time_t tseed = rsb_time();
		unsigned int uiseed=(*(unsigned int*)(&tseed));
		RSB_WARN("#Starting library with enabled malloc fault injection.\n# Initializing with random seed of value: %u\n",uiseed);
		/* In this way, the user may introduce faults by recompiling the code. TODO: rsb_lib_init() based seed passing. */
		srand(uiseed);
	}
#endif /* RSB_WANT_RANDOM_MALLOC_FAULT_INJECTION */

	RSB_DO_ERR_RETURN(errval)
}

void *rsb__memcpy(void *RSB_RESTRICT dest, const void *RSB_RESTRICT src, size_t n)
{
	/*!
	  	\ingroup gr_internals

		Say you want to use a custom memcpy function
		or perform some statistics measurement.

		Well, this is the place to hack.
		TODO: use the 'restrict' keyword.
	*/
#if 0
	{
		register unsigned char*dp=NULL;
		register const unsigned char*sp=NULL;
		for(dp=dest,sp=src;RSB_LIKELY(dp<dest+n);++dp,++sp)
			*dp=*sp;
	}
#else
	return memcpy(dest,src,n);
#endif
}

size_t rsb__sys_free_system_memory(void)
{
	/*!
	  	\ingroup gr_internals

		System free memory, or 0.
	*/
        size_t free_mem=0;
        long int pagesize =0;
        long int mem_pages=0;

#ifdef RSB_HAVE_SYSCONF
#if   defined(PAGESIZE)
        pagesize=sysconf(PAGESIZE);
#elif defined(_SC_PAGESIZE)
        pagesize=sysconf(_SC_PAGESIZE);
#elif defined(PAGE_SIZE)
        pagesize=sysconf(PAGE_SIZE);
#else /* PAGESIZE */
#endif /* PAGESIZE */
#if   defined(_SC_AVPHYS_PAGES)
        mem_pages=sysconf(_SC_AVPHYS_PAGES);
#endif /* _SC_AVPHYS_PAGES */
#endif /* RSB_HAVE_SYSCONF */
	if(pagesize<1 || mem_pages<1)
		free_mem=0;
	else
		free_mem=((size_t)pagesize)*((size_t)mem_pages);
	return free_mem;
}

size_t rsb__sys_total_system_memory(void)
{
	/*!
	  	\ingroup gr_internals
	*/
        size_t tot_mem=0;
        long int pagesize =0;
        long int mem_pages=0;

#ifdef RSB_HAVE_SYSCONF
#if   defined(PAGESIZE)
        pagesize=sysconf(PAGESIZE);
#elif defined(_SC_PAGESIZE)
        pagesize=sysconf(_SC_PAGESIZE);
#elif defined(PAGE_SIZE)
        pagesize=sysconf(PAGE_SIZE);
#else /* PAGE_SIZE */
#endif /* PAGE_SIZE */
#if   defined(_SC_PHYS_PAGES)
        mem_pages = sysconf(_SC_PHYS_PAGES);
#endif /* _SC_PHYS_PAGES */
#endif /* RSB_HAVE_SYSCONF */
	if(pagesize<1 || mem_pages<1)
		tot_mem=0;
	else
		tot_mem=((size_t)pagesize)*((size_t)mem_pages);
	return tot_mem;
}

static long rsb_set_executing_threads(long tn)
{
	/*!
	 	FIXME: new
	  	\ingroup gr_internals
	*/
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	/* multi threaded case */
	if(tn > RSB_CONST_MAX_SUPPORTED_CORES)
	{
		RSB_ERROR("cannot set %ld threads: a maximum of %ld is supported\n",tn,RSB_CONST_MAX_SUPPORTED_CORES);
		return RSB_CONST_MIN_SUPPORTED_CORES;
	}
	tn = RSB_MIN(tn,RSB_CONST_MAX_SUPPORTED_CORES);
	if(tn < RSB_CONST_MIN_SUPPORTED_CORES)
	{
		/* a value < 0 means the user wants the threads count to be set automatically */
	//	return 1;
		tn = rsb_global_session_handle.rsb_g_threads;
	}

#if (RSB_TIME_SET_THREADS==1)
	rsb_time_t dt = -rsb_time();
#endif
#if (RSB_USE_OMP_SET_NUM_THREADS==1)
	omp_set_num_threads(tn);
#endif
#if (RSB_TIME_SET_THREADS==1)
	dt += rsb_time();
	RSB_STDOUT("setting threads (%d) took %lf s\n",tn,dt);
#endif
	/* FIXME : 20101111 on my GNU box, the following does not return tn, but seems to have effect. Weird */
	//tn = omp_get_num_threads();
	rsb_global_session_handle.rsb_want_threads = tn;/* FIXME : a hack */
	return tn;
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	/* single threaded case */
	return RSB_CONST_MIN_SUPPORTED_CORES;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
}

rsb_err_t rsb__lock_as_memory_resident(rsb_bool_t dolock)
{
	/* FIXME: need mechanisms to get/set/restore this setting */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_HAVE_MLOCKALL
	int retval = 0;
	/* "mlockall() locks all pages mapped into the address space of the calling process. "*/
	if(dolock)
	{
		retval = mlockall(MCL_FUTURE|MCL_CURRENT); 
	}
	else
	{
		retval = munlockall(); 
	}
	if(retval)
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		/* FIXME: shall introduce RSB_ERR_SYSCALL_ERROR */
		RSB_ERROR(RSB_ERRM_FCOVMU);
	}
#else /* RSB_HAVE_MLOCKALL */
	RSB_ERROR(RSB_ERRM_COVMUINS);
	errval = RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_HAVE_MLOCKALL */
	return errval; /* what about retval ? */
}

int rsb__fileno(FILE *stream)
{
#ifndef RSB_HAVE_FILENO
	return -1;
#else /* RSB_HAVE_FILENO */
	return fileno(stream);
#endif /* RSB_HAVE_FILENO */
}

rsb_int_t rsb__set_num_threads(rsb_int_t tn)
{
	/*!
	   	\ingroup rsb_doc_library rsb_doc_rsb 

	 	Gets and/or sets number of librsb running threads.
		If \a tn is RSB_THREADS_AUTO, sets threads count to a default value.
		If \a tn is RSB_THREADS_GET, only returns the count of active threads.
		If \a tn is RSB_THREADS_GET_MAX, returns max count of supported threads.
	 	\todo: shall promote this as a rsb.h function and make all similar ones static, and put them in thread.c !
		\return: number of running threads.
	 */

	long rtn=0;

	switch(tn)
	{
		case(RSB_THREADS_GET):
		rtn = rsb_want_executing_threads();
		break;
		case(RSB_THREADS_AUTO):
		tn=0;
		break;
		case(RSB_THREADS_GET_MAX_SYS):
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		/* rtn = omp_get_max_threads(); */
		/* rtn = omp_get_thread_limit(); */
		rtn = rsb_global_session_handle.rsb_g_threads;
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		rtn = 1;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		break;
		case(RSB_THREADS_GET_MAX_LIB):
		rtn = RSB_CONST_MAX_SUPPORTED_CORES;
		break;
		case(RSB_THREADS_GET_MAX):
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		rtn = omp_get_max_threads();
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		rtn = RSB_MIN(rtn,RSB_CONST_MAX_SUPPORTED_CORES);
		break;
		default:
		rtn = rsb_set_executing_threads(tn);
	}
	RSB_DEBUG_ASSERT(rtn>0);
	return (rsb_int_t) rtn;
}

#if defined(RSB_HAVE_EXECINFO_H) && RSB_OUT_ERR_VERBOSITY>=2
void rsb__print_trace(void)
{
	/* according to a glibc docs example */
	void *array[10];
	size_t size;
	char **strings;
	size_t i;
	size = backtrace (array, 10);
	strings = backtrace_symbols (array, size);
	printf ("Obtained %zd stack frames.\n", size);
	for (i = 0; i < size; i++)
		printf ("%s\n", strings[i]);
	free (strings);
}
#endif /* RSB_HAVE_EXECINFO_H */

#if RSB_USE_GETRUSAGE
static rsb_time_t rsb_rstv(struct timeval*tvp)
{
        register double t = 0.0;
        t  =  (double)(tvp->tv_sec) + ((double)(tvp->tv_usec))*1.e-6;
        return t;
}

#define RSB_K 1024
rsb_err_t rsb__getrusage(void)
{
	/*
	 * Shall work independently from rsb_lib_init/rsb_lib_exit.
	 * */
	struct rusage usage;
	int gru = getrusage(RUSAGE_SELF,&usage);

	RSB_STDOUT("getrusage() stats:\n");
	/*("ru_ixrss : %ld (integral shared memory size)\n",usage.ru_ixrss);*/
#ifdef __APPLE__
	/* ru_maxrss not necessarily available; unit might also differ... */
#else
	RSB_STDOUT("ru_maxrss: %ld (maximum resident set size -- MB)\n",usage.ru_maxrss / RSB_K);
#endif
	RSB_STDOUT("ru_stime : %0.4lgs (system CPU time used)\n",rsb_rstv(&usage.ru_stime));
	RSB_STDOUT("ru_utime : %0.4lgs (user CPU time used)\n",rsb_rstv(&usage.ru_utime));
#if 0
	RSB_STDOUT("ru_utime : %0.4lg (user page faults (hard page faults))\n",rsb_rstv(&usage.ru_majflt));
	RSB_STDOUT("ru_utime : %0.4lg (page reclaims (soft page faults))\n",rsb_rstv(&usage.ru_minflt));
#endif

	return gru == 0 ? RSB_ERR_NO_ERROR : RSB_ERR_GENERIC_ERROR;
}
#else /* RSB_USE_GETRUSAGE */
rsb_err_t rsb__getrusage(void)
{
	return RSB_ERR_NO_ERROR;
}
#endif /* RSB_USE_GETRUSAGE */

const rsb_char_t * rsb__getenv(const rsb_char_t * name)
{
	const rsb_char_t * evv = NULL;
	
#ifdef RSB_HAVE_GETENV
	evv = getenv(name);
#endif /* RSB_HAVE_GETENV */

	return evv;
}

const rsb_char_t * rsb__getenv_nnr(const rsb_char_t * name)
{
	RSB_DEBUG_ASSERT( name != NULL );
	return rsb__getenv(name) ? rsb__getenv(name) : name + strlen(name);
}

#if RSB_WITH_HWLOC
#if HWLOC_API_VERSION >= 0x00020000
static hwloc_uint64_t rsb__hwloc_max_cs_at_level(hwloc_topology_t topology, hwloc_obj_t obj, unsigned int depth, unsigned int cd)
{
	/* for rsb__get_lnc_size_hwloc */
    	unsigned int i, cl;
	hwloc_uint64_t cs=0;

	switch (obj->type)
	{
		case(HWLOC_OBJ_L1CACHE): cl=1; break;
		case(HWLOC_OBJ_L2CACHE): cl=2; break;
		case(HWLOC_OBJ_L3CACHE): cl=3; break;
		case(HWLOC_OBJ_L4CACHE): cl=4; break;
		case(HWLOC_OBJ_L5CACHE): cl=5; break;
		default: cl=0;
	}
	if(cl && cl==cd)
		cs = obj->attr->cache.size;

	for (i = 0; i < obj->arity; i++)
	{
		hwloc_uint64_t ccs = rsb__hwloc_max_cs_at_level(topology, obj->children[i], depth + 1, cd);
    		cs = ccs > cs ? ccs : cs;
    	}
	return cs;
}
#endif /* HWLOC_API_VERSION >= 0x00020000 */
#endif /* RSB_WITH_HWLOC */

#if RSB_WITH_HWLOC
long rsb__get_lnc_size_hwloc(int n)
{
	/* Gets cache size using hwloc.h. EXPERIMENTAL */
	long size = 0;
#if RSB_WITH_HWLOC
	hwloc_uint64_t cs;
	hwloc_topology_t topology;
	hwloc_obj_t obj;
	hwloc_topology_init(&topology);
	hwloc_topology_load(topology);
#if HWLOC_API_VERSION >= 0x00020000
	cs = rsb__hwloc_max_cs_at_level(topology, hwloc_get_root_obj(topology), 0, n);
	size = (long) cs;
#else /* HWLOC_API_VERSION */
	int levels = 0;
	for (obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_PU, 0); obj; obj = obj->parent)
#if HWLOC_API_VERSION >= 0x20000
		if (obj->type >= HWLOC_OBJ_L1CACHE && obj->type <= HWLOC_OBJ_L5CACHE)
#else
		if (obj->type == HWLOC_OBJ_CACHE)
#endif
			if(++levels == n)
        			size = obj->attr->cache.size;
#endif	/* HWLOC_API_VERSION */
    	hwloc_topology_destroy(topology);
#endif	/* RSB_WITH_HWLOC */
	return size;
}
#endif	/* RSB_WITH_HWLOC */

void rsb__strcpy_hostname(rsb_char_t * buf)
{
	char * setaname = NULL;
	const size_t len = RSB_MAX_HOSTNAME_LEN;
	char name[len+1];
	name[0] = name[len] = RSB_NUL;
#if 0 /* gethostname is POSIX.1-2001 but not C99 */
#if RSB_HAVE_GETHOSTNAME
	if ( setaname == NULL && 0 == gethostname(name, len) && name[0] )
		setaname = &name[0];
#endif /* RSB_HAVE_GETHOSTNAME */
#endif /* 0 */
#ifdef RSB_HAVE_SYS_UTSNAME_H 
	if ( setaname == NULL )
	{
		struct utsname un;
		if(uname(&un)==0)
		{
			strncpy(name,un.nodename,len);
			setaname = &name[0];
		}
#if 0
           struct utsname {
               char sysname[];
               char nodename[];
               char release[];
               char version[];
               char machine[];
           #ifdef _GNU_SOURCE
               char domainname[];
           #endif /* _GNU_SOURCE */
           };
#endif
	}
#endif /* RSB_HAVE_SYS_UTSNAME_H */
	strcpy(buf,name);
}

/* @endcond */
