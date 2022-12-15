/*

Copyright (C) 2008-2021 Michele Martone

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
 * @brief
 * This source file contains matrix I/O functions.
 * */

// FIXME: this code is messy and unclean

#include "rsb_internals.h"
#ifdef RSB_HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif /* RSB_HAVE_SYS_STAT_H */
#if RSB_WANT_ZLIB_SUPPORT
#include <zlib.h>
#endif /* RSB_WANT_ZLIB_SUPPORT */
#include <stdio.h>
#define RSB_MMIOH_CL 4 
#define RSB_20120321_IOBUFFERING 0
#if RSB_WANT_OMPIO_SUPPORT
#include "rsb_ompio.h"
#endif /* RSB_WANT_OMPIO_SUPPORT */
#define RSB_PATCH_GILLES_20130906 0
#if RSB_PATCH_GILLES_20130906
//#define _GNU_SOURCE
#include <asm/fcntl.h>
#include <unistd.h>		/* O_DIRECT */
#include <sys/types.h>          /* See NOTES */
#include <sys/stat.h>
#endif

#define rsb_util_sort_column_major(VA,IA,JA,nnz,nr,nc,typecode,flags) rsb__util_sort_row_major_inner(VA,JA,IA,nnz,nc,nr,typecode,flags)
#if 1
#define RSB_IO_VERBOSE_MSG(NZI,NNZ)
#else
#define RSB_IO_VERBOSE_GRNLRT RSB_MILLION_I
// #define RSB_IO_VERBOSE_GRNLRT 100000
#define RSB_IO_VERBOSE_MSG(NZI,NNZ) if((NZI)%RSB_IO_VERBOSE_GRNLRT==0) RSB_STDERR("%s%dM/%dM\n",RSB_CLEARTERM_STRING,(NZI)/RSB_IO_VERBOSE_GRNLRT,(NNZ)/RSB_IO_VERBOSE_GRNLRT )
#endif

#define RSB_FILE_ALLOW_LOAD_EMPTY_PATTERN 1 /* 20140324 */

RSB_INTERNALS_COMMON_HEAD_DECLS

#if RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE
#define RSB__VERBOSE_IO_ERRORS rsb_global_session_handle.rsb_g_verbose_interface
#define RSB__VERBOSE_IO_WARNINGS rsb_global_session_handle.rsb_g_verbose_interface ? RSB_BOOL_TRUE : RSB_BOOL_FALSE
#else /* RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE */
#define RSB__VERBOSE_IO_ERRORS 0
#define RSB__VERBOSE_IO_WARNINGS RSB_BOOL_FALSE
#endif /* RSB_WANT_DEBUG_VERBOSE_INTERFACE_NOTICE */

#ifdef RSB_WITH_MM
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__util_mm_load_coo_matrix(const char *filename, struct rsb_coo_mtx_t * cmp)
{
	/**
	 * \ingroup gr_internals
	 *
	 * Loads in a matrix in unsorted COO format. (new)
	 *
	 * FIXME : UNTESTED
	 *
	 * \note used by experiment.c files
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t cm;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

	if(!cmp || !filename)
	{
		errval = RSB_ERR_BADARGS;
	       	RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	RSB_BZERO_P(&cm);
	cm.typecode = cmp->typecode;	// should be specified
	errval = rsb__util_mm_load_matrix_f(filename, &cm.IA, &cm.JA,&cm.VA , &cm.nr, &cm.nc, &cm.nnz , cm.typecode, flags, NULL, NULL);
	if(RSB_SOME_ERROR(errval))
		goto err;
	rsb__memcpy(cmp,&cm,sizeof(cm));
err:
	rsb__do_perror(NULL,errval);
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__util_mm_info_matrix_f(const char *fn,  rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t *typecode, rsb_bool_t * is_symmetric, rsb_bool_t * is_hermitian, rsb_bool_t * is_pattern, rsb_bool_t * is_lower, rsb_bool_t * is_upper , rsb_bool_t * is_vector)
{
	/*! 
	 *  \ingroup internals
	 *  FIXME : does not return cleanly in case of errors
	 * */

  	FILE * fd = NULL;
	rsb_coo_idx_t innz = 0;	/* FIXME */
	rsb_coo_idx_t _m = 0,_k = 0;
	char matcode[RSB_MMIOH_CL]; // !?
	rsb_bool_t is_vector_ = RSB_BOOL_FALSE;
	int rbc = 0;

	if(nnz)*nnz = RSB_MARKER_NNZ_VALUE ;
	if(m)*m = RSB_MARKER_COO_VALUE;
	if(k)*k = RSB_MARKER_COO_VALUE;
	/* TODO: needs to define some RSB_NUMERICAL_COMPLEX_TYPE_DEFAULT macro and use it here, in case of a complex matrix */
	if(typecode && RSB_MATRIX_UNSUPPORTED_TYPE(*typecode))*typecode = RSB_NUMERICAL_TYPE_DEFAULT;

	if(!fn)
	{
		return RSB_ERR_BADARGS;
	}

	if(fd==NULL)
	if ((fd = RSB_FOPEN(fn, "r")) == NULL)
	{
		if( RSB__VERBOSE_IO_ERRORS )
			RSB_STDERR("Failed opening file: %s\n",fn);
		return RSB_ERR_GENERIC_ERROR;
	}

#if RSB_WANT_ZLIB_SUPPORT
	if ((rbc=rsb__mm_read_banner(NULL,fd,&(matcode))) != 0)
#else
	if ((rbc=rsb__mm_read_banner(fd,NULL,&(matcode))) != 0)
#endif
	{
		if( rbc == MM_LIKELY_GZIPPED_FILE )
        		RSB_STDERR("Trying to load a gzipped file as matrix without gzip decoder in ?!\n");
		if( RSB__VERBOSE_IO_ERRORS )
        		RSB_STDERR("Could not process Matrix Market banner.\n");
		RSB_FCLOSE(fd);
		return RSB_ERR_GENERIC_ERROR;
	}

	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */

	if ( !rsb_mm_is_matrix(matcode) )
	{
    		char *mstr = rsb__mm_typecode_to_str(matcode);
		RSB_STDERR("Sorry, this application does not support "
				"Matrix Market type: [%s]\n", mstr?mstr:"?");
		RSB_FCLOSE(fd);
	       	if(mstr)
			free(mstr);
 		return RSB_ERR_UNSUPPORTED_TYPE;
	}

	is_vector_ = (rsb_mm_is_sparse(matcode))?RSB_BOOL_FALSE:RSB_BOOL_TRUE;

	/* find out size of sparse matrix .... */

#if RSB_WANT_ZLIB_SUPPORT
	if( ((is_vector_) && (rsb__mm_read_mtx_array_size(NULL,fd,&_m,&_k) !=0)) || ((!is_vector_) && (rsb__mm_read_mtx_crd_size(NULL,fd,&_m,&_k,&innz)) !=0) )
#else
	if( ((is_vector_) && (rsb__mm_read_mtx_array_size(fd,NULL,&_m,&_k) !=0)) || ((!is_vector_) && (rsb__mm_read_mtx_crd_size(fd,NULL,&_m,&_k,&innz)) !=0) )
#endif
	{
		RSB_FCLOSE(fd);
        	return RSB_ERR_GENERIC_ERROR;
	}
	if(m)*m = (rsb_coo_idx_t)_m;
	if(k)*k = (rsb_coo_idx_t)_k;

	if(is_vector)
	{
		*is_vector = is_vector_;
	}

	if(is_pattern)
	{
		if(rsb_mm_is_pattern(matcode))
			*is_pattern = RSB_BOOL_TRUE;
		else
			*is_pattern = RSB_BOOL_FALSE;
	}

	if(is_symmetric)
	{
		if(rsb_mm_is_symmetric(matcode))
			*is_symmetric = RSB_BOOL_TRUE;
		else
			*is_symmetric = RSB_BOOL_FALSE;
	}

	if(is_hermitian)
	{
		if(rsb_mm_is_hermitian(matcode))
			*is_hermitian = RSB_BOOL_TRUE;
		else
			*is_hermitian = RSB_BOOL_FALSE;
	}

	if(is_lower)
	{
		// TODO
	}

	if(is_upper)
	{
		// TODO
	}

	if(m && k)
	if (((int)*m != _m)||((int)*k != _k))
	{
		/* overflow */
		RSB_IO_ERROR("overflow error while reading matrix dimensions.\n");
        	return RSB_ERR_INTERNAL_ERROR;
	}

	if(is_vector_)
		innz = _m*_k; /* 20120904 */
	if(nnz)
	{
		*nnz = innz;
#if RSB_FILE_ALLOW_LOAD_EMPTY_PATTERN
		if(*nnz<0)
#else
		if(*nnz<1)
#endif
			return RSB_ERR_GENERIC_ERROR;
	}

	RSB_FCLOSE(fd);
	return RSB_ERR_NO_ERROR;
}

static int rsb_zfscanf(FILE * fd,const char * fs,rsb_coo_idx_t *IV, rsb_coo_idx_t *JV, void * VAR, void * VAI, void * ngzfd)
{
	/**
	 *  \ingroup internals
	 * */
	if(fd)
	{
		if((!IV) && (!JV))
		{
			if(VAI)
				return fscanf(fd,fs,VAR,VAI);
			else
				return fscanf(fd,fs,VAR);
		}
		else
		{
			if(VAI)
				return fscanf(fd,fs,IV,JV,VAR,VAI);
			if(VAR)
				return fscanf(fd,fs,IV,JV,VAR);
			else
				return fscanf(fd,fs,IV,JV);
		}
	}
	else
	{
#if RSB_WANT_ZLIB_SUPPORT
		return rsb__fscanf(ngzfd,fs,IV,JV,VAR,VAI);
#else /* RSB_WANT_ZLIB_SUPPORT */
		return 0; /* error case */
#endif /* RSB_WANT_ZLIB_SUPPORT */
	}
}

int rsb__fscanf(FILE * fd,const char * fs,rsb_coo_idx_t *IV, rsb_coo_idx_t *JV, void * VAR, void * VAI)
{
	/**
	 *  \ingroup internals
	 *  FIXME: error handling is missing
	 * */
#if RSB_WANT_ZLIB_SUPPORT
	char line[MM_MAX_LINE_LENGTH];
	gzgets((gzFile)fd,line,MM_MAX_LINE_LENGTH);

	if((!IV) && (!JV))
	{
		if(VAI)
			return sscanf(line,fs,VAR,VAI);
		if(VAR)
			return sscanf(line,fs,VAR);
		else
			return 0;
	}
	else
	{
		if(VAI)
			return sscanf(line,fs,IV,JV,VAR,VAI);
		if(VAR)
			return sscanf(line,fs,IV,JV,VAR);
		else
			return sscanf(line,fs,IV,JV);
	}
#else /* RSB_WANT_ZLIB_SUPPORT */
	if((!IV) && (!JV))
	{
		if(VAI)
			return fscanf(fd,fs,VAR,VAI);
		if(VAR)
			return fscanf(fd,fs,VAR);
		else
			return 0;
	}
	else
	{
		if(VAI)
			return fscanf(fd,fs,IV,JV,VAR,VAI);
		if(VAR)
			return fscanf(fd,fs,IV,JV,VAR);
		else
			return fscanf(fd,fs,IV,JV);
	}
#endif /* RSB_WANT_ZLIB_SUPPORT */
}

char * rsb__fgets(char* RSB_RESTRICT buf, int len, FILE * RSB_RESTRICT fd)
{
	/**
	 *  \ingroup internals
	 * */
#if RSB_WANT_ZLIB_SUPPORT
	return gzgets((gzFile)fd,buf,len);
#else /* RSB_WANT_ZLIB_SUPPORT */
	return fgets(buf,len,fd);
#endif /* RSB_WANT_ZLIB_SUPPORT */
}

rsb_err_t rsb__util_mm_load_vector_f(const char *fn, void **VA, rsb_nnz_idx_t *nnz, rsb_type_t typecode)
{
	/* FIXME: in perpective, need stride, C/Fortran order, etc ... */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	errval = rsb__util_mm_load_matrix_f(fn,NULL,NULL,VA,NULL,NULL,nnz,typecode,RSB_FLAG_NOFLAGS,NULL,NULL);
	return errval;
}

#if RSB_WANT_EXPERIMENTAL_BINARY_COO
int rsb_getc(FILE *fd, void *ngzfd)
{
	/**
	 * \ingroup gr_internals
	 * */
	return fd ? getc(fd) : RSB_GETC(ngzfd);
}

int rsb_ungetc(int cc, FILE *fd, void *ngzfd)
{
	/**
	 * \ingroup gr_internals
	 * */
	return fd ? ungetc(cc,fd) : RSB_UNGETC(cc,ngzfd);
}

size_t rsb_fread(void* buf, size_t size, size_t nitems, void * fd, void *ngzfd)
{
	/**
	 * \ingroup gr_internals
	 * */
	return fd ?
		fread(buf, size, nitems, fd)
		:
		RSB_FREAD(buf, size, nitems, ngzfd);
}

rsb_err_t rsb__read_coo_bin_fd(FILE *fd, void *ngzfd, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, void *VA, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const rsb_nnz_idx_t nnz, const rsb_type_t typecode)
{
	/**
	 * \ingroup gr_internals
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	int c = rsb_getc(fd, ngzfd);
	int l = c;
	const size_t nmemb = nnz;
	void * TA = NULL;
	rsb_nnz_idx_t i;

	while( ( c = rsb_getc(fd, ngzfd) ) )
		l = c;
	if( RSB_MATRIX_UNSUPPORTED_TYPE(l) )
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		RSB_PERR_GOTO(err,RSB_ERRM_UNSUPPORTED_TYPE);
	}
	if( l != typecode )
	{
		TA = rsb__malloc_vector(nmemb,l);
		if ( !TA )
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err,RSB_ERRM_ENOMEM);
		}
	}
	else
		TA = VA;
	errval = RSB_ERR_INTERNAL_ERROR;
	if( rsb_fread( IA,sizeof(rsb_coo_idx_t),nmemb,fd,ngzfd) != nmemb )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if( rsb_fread( JA,sizeof(rsb_coo_idx_t),nmemb,fd,ngzfd) != nmemb )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if( rsb_fread( TA,RSB_SIZEOF(l),nmemb,fd,ngzfd) != nmemb )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = RSB_ERR_NO_ERROR;
	for ( i=0; i < nnz; i++ )
	{
		(IA)[i]--, (JA)[i]--; /* adjust from 1-based to 0-based (may revive rsb__util_nnz_array_from_fortran_indices for this ?) */
		if ( IA[i] < 0 || IA[i] > nr || JA[i] < 0 || JA[i] > nc )
		{
			// TODO: move this check to be called from rsb__util_mm_load_matrix_f()
			errval = RSB_ERR_CORRUPT_INPUT_DATA;
			RSB_PERR_GOTO(err,RSB_ERRM_CORRUPT_INPUT_DATA);
		}
	}

	if( l != typecode )
		errval = rsb__do_copy_converted_scaled(TA, VA, NULL, l, typecode, nnz, RSB_TRANSPOSITION_N);
err:
	if( l != typecode )
		RSB_CONDITIONAL_FREE(TA);
	return errval;
}
#endif /* RSB_WANT_EXPERIMENTAL_BINARY_COO */

#ifdef RSB_WANT_LONG_IDX_TYPE 
#define RSB_IJPC "%lld %lld"
#else /* RSB_WANT_LONG_IDX_TYPE */
#define RSB_IJPC   "%d %d"
#endif /* RSB_WANT_LONG_IDX_TYPE */
rsb_err_t rsb__util_mm_load_matrix_f(const char *fn, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void **VA, rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t typecode, rsb_flags_t flags, rsb_bool_t *is_lowerp, rsb_bool_t *is_upperp)
{
	/**
	 *  \ingroup internals
	 *  This function reads in a single Matrix Market format matrix drom a specified filename.
	 *  The matrix will be loaded in coordinate storage format, and IA,JA, VA arrays will 
	 *  be allocated for this purpose, here.
	 *
	 * \param fn	should contain a valid matrix file name
	 * \param IA	should point to a pointer which will be allocated here to contain the elements row values
	 * \param JA	should point to a pointer which will be allocated here to contain the elements column values
	 * \param VA	should point to a pointer which will be allocated here to contain the elements column values
	 * \param m	should point to the matrix rows count variable, which will be set in this function
	 * \param k	should point to the matrix columns count variable, which will be set in this function
	 * \param nnz	should point to the matrix nonzero count variable, which will be set in this function (shall be initialized to zero in advance!)
	 * \param type	should specify a valid numerical type to convert the read data in. see rsb.h for this.
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 *
	 * Notes: 
	 *
	 *	There is no guarantee that the loaded matrix will be free from duplicate values.
	 *	The specified numerical type should be enabled at library generaton time.
	 *
	 * FIXME : lots of ths function's code should be generated from macros.
	 * FIXME : error handling is awful
	 * FIXME : otype stands for 'original' or 'output' ?
	 * TODO  : some option to specify a pattern-only load
	 * FIXME : in a future version, should not allocate if pointers not NULL
	 * TODO: may print (at least for rsbench) out how many discarded duplicates, how many zeroes, etc etc
	 * TODO: may support different styles for tolerating / detecting e.g. reading a vector file instead a matrix one.
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t re = 0;/* read elements */
	FILE *fd = NULL;
//#if RSB_WANT_ZLIB_SUPPORT
	FILE *ngzfd = NULL;
//#endif
	const rsb_bool_t wv = RSB__VERBOSE_IO_WARNINGS;
	rsb_nnz_idx_t i = 0;
	rsb_nnz_idx_t annz = 0;/* allocated nnz */
	rsb_bool_t is_symmetric = RSB_BOOL_FALSE,is_lower = RSB_BOOL_FALSE,is_upper = RSB_BOOL_FALSE,is_pattern = RSB_BOOL_FALSE, is_hermitian = RSB_BOOL_FALSE, is_vector = RSB_BOOL_FALSE;/* FIXME : no expansion support for hermitian */
	rsb_bool_t aja = RSB_BOOL_FALSE, ava = RSB_BOOL_FALSE, aia = RSB_BOOL_FALSE;
	rsb_flags_t otype = typecode;/* original type */
	rsb_time_t frt = 0;
	char matcode[RSB_MMIOH_CL]; // !?
#if RSB_WANT_ZLIB_SUPPORT
	rsb_bool_t is_gz = RSB_BOOL_FALSE;
#endif /* RSB_WANT_ZLIB_SUPPORT */

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	double  **dval = NULL;
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT
	float **fval = NULL;
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */
	#ifdef RSB_NUMERICAL_TYPE_CHAR
	char **cval = NULL;
	#endif /* RSB_NUMERICAL_TYPE_CHAR */
	#ifdef RSB_NUMERICAL_TYPE_INT
	int  **ival = NULL;
	#endif /* RSB_NUMERICAL_TYPE_INT */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	float complex  **zval = NULL;
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	double complex  **Zval = NULL;
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
	rsb_coo_idx_t innz = 0;
	rsb_coo_idx_t _m = 0,_k = 0;
#if RSB_20120321_IOBUFFERING
	char*iobuf = NULL;
	size_t iobbs = 16*1024*1024;
#endif /* RSB_20120321_IOBUFFERING */

	frt = - rsb_time();

	if ( RSB_MATRIX_UNSUPPORTED_TYPE(typecode) )return RSB_ERR_UNSUPPORTED_TYPE	;

	//if(!VA || !JA || !IA) return RSB_ERR_BADARGS;
	if( (!(VA && JA && IA)) && (!(VA && (!JA) && (!IA))) ) return RSB_ERR_BADARGS;
	if(IA)aia = RSB_BOOL_IS_POINTER_NON_NULL(*IA);
	if(JA)aja = RSB_BOOL_IS_POINTER_NON_NULL(*JA);
	if(VA)ava = RSB_BOOL_IS_POINTER_NON_NULL(*VA);
	//if(!nnz || !k || !m) return RSB_ERR_BADARGS;
	if((!(nnz && k && m)) && (!(nnz && (!k) && (!m)))) return RSB_ERR_BADARGS;
	if(*nnz)
		annz = *nnz;// if user set

	errval = rsb__util_mm_info_matrix_f(fn,m,k,nnz,&typecode,&is_symmetric,&is_hermitian,&is_pattern,&is_lower,&is_upper,&is_vector);
	if(RSB_SOME_ERROR(errval))
		goto prerr;
	if(annz==0)
		annz = *nnz;
	if(annz<*nnz)
	{
		RSB_ERROR("user-set array size (%d) does not fit actual input (%d)\n",(int)*nnz,(int)annz);
		errval = RSB_ERR_BADARGS;
		goto prerr;
	}
	
	if(is_pattern)
	#ifdef RSB_NUMERICAL_TYPE_PATTERN
		typecode = RSB_NUMERICAL_TYPE_PATTERN;
	#else /* RSB_NUMERICAL_TYPE_PATTERN */
		return RSB_ERR_UNSUPPORTED_FEATURE;
	#endif /* RSB_NUMERICAL_TYPE_PATTERN */

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if (otype == RSB_NUMERICAL_TYPE_DOUBLE)dval = (double**)(VA);
	else
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT
	if (otype == RSB_NUMERICAL_TYPE_FLOAT)fval = (float **)(VA);
	else
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */
	#ifdef RSB_NUMERICAL_TYPE_INT
	if (otype == RSB_NUMERICAL_TYPE_INT)ival = (int **)(VA);
	else
	#endif /* RSB_NUMERICAL_TYPE_INT */
	#ifdef RSB_NUMERICAL_TYPE_CHAR
	if (otype == RSB_NUMERICAL_TYPE_CHAR)cval = (char **)(VA);
	else
	#endif /* RSB_NUMERICAL_TYPE_CHAR */
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	if (otype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX)Zval = (double complex **)(VA);
	else
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	if (otype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX)zval = (float complex **)(VA);
	else
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_PATTERN
	if (otype == RSB_NUMERICAL_TYPE_PATTERN){/* nothing to do */}
	else
	#endif /* RSB_NUMERICAL_TYPE_PATTERN */
	return RSB_ERR_UNSUPPORTED_TYPE	;

	if(IA)
	*IA = NULL;
	if(JA)
	*JA = NULL;
  	
	{
		// THIS IS A CODE DUPLICATION ..
#if RSB_WANT_ZLIB_SUPPORT
	is_gz = RSB_BOOL_FALSE;
	if ((fd = (FILE*)gzopen(fn,"r")) != NULL)
	{
		if(gzdirect((gzFile)fd))
			;
		else
			is_gz = RSB_BOOL_TRUE;
		RSB_FCLOSE(fd);
		fd=NULL;
	}
#endif /* RSB_WANT_ZLIB_SUPPORT */

#if RSB_WANT_ZLIB_SUPPORT
	if(is_gz)
	{
		if((fd = (FILE*)gzopen(fn,"r")) == NULL)
		{
			/* TODO: the following code is not robust, shall fix it. */
#ifndef RSB_HAVE_DUP
#error Functions 'dup/fileno' is not present ? Reconfigure without Z library then!
			const int fnum = 0;/* */
			ngzfd = NULL;
#else /* RSB_HAVE_DUP */
			const int fnum = dup( rsb__fileno(fd) );
			ngzfd = (FILE*)gzdopen(fnum,"r");
#endif /* RSB_HAVE_DUP */
			if(!ngzfd)
			{
				RSB_FCLOSE(fd);
				RSB_ERROR(RSB_ERRMSG_FILEOPENPGZ"\n");
				return RSB_ERR_GENERIC_ERROR;
			}
			else
				;
		}
		else
			ngzfd = fd, fd = NULL;
	}
	else
#endif /* RSB_WANT_ZLIB_SUPPORT */
#if !RSB_PATCH_GILLES_20130906
	if ((fd = fopen(fn,"r")) == NULL)
#else /* RSB_PATCH_GILLES_20130906 */
{
		int _fd;
		if ((_fd = open(fn,O_RDONLY|O_DIRECT)) < 0)
		{
			RSB_STDERR(RSB_ERRMSG_FILEOPENP" %s\n",fn);
			return RSB_ERR_GENERIC_ERROR;
		}
		else if ((fd = fdopen(_fd,"r")) == NULL) 
#endif /* RSB_PATCH_GILLES_20130906 */
	{
		RSB_STDERR(RSB_ERRMSG_FILEOPENP" %s\n",fn);
		return RSB_ERR_GENERIC_ERROR;
	}
	else
		;//ngzfd = fd;
#if RSB_20120321_IOBUFFERING
	//iobbs = BUFSIZ;
	if(iobbs>0)
	if(((iobuf = rsb__malloc(iobbs))==NULL) 
			//|| (0!=setvbuf(ngzfd,NULL,_IOLBF,0))// line buffering: slow
			//|| (0!=setvbuf(ngzfd,NULL,_IONBF,0))// no buffering: super-slow
			|| (0!=setvbuf(ngzfd,iobuf,_IOFBF,iobbs)) // seem to not work
				)
	{
		RSB_STDERR("problems setting up a buffer for file ""%s\n",fn);
	}
	//setbuf(ngzfd,iobbs);
	//setbuffer(ngzfd,iobuf,iobbs);
#endif /* RSB_20120321_IOBUFFERING */
#if RSB_PATCH_GILLES_20130906
}
#endif /* RSB_PATCH_GILLES_20130906 */

	if (rsb__mm_read_banner(fd,ngzfd,&(matcode)) != 0)
	{
        	RSB_STDERR(RSB_ERRMSG_TMXMKTBANNER".\n");
		RSB_FCLOSE(fd);
		return RSB_ERR_GENERIC_ERROR;
	}

	if( ((is_vector) && (rsb__mm_read_mtx_array_size(fd,ngzfd,&_m,&_k) !=0)) || ((!is_vector) && (rsb__mm_read_mtx_crd_size(fd,ngzfd,&_m,&_k,&innz)) !=0) )
	{
		RSB_FCLOSE(fd);
        	return RSB_ERR_GENERIC_ERROR;
	}
	if(m)
	*m = (rsb_coo_idx_t)_m;
	if(k)
	*k = (rsb_coo_idx_t)_k;
	}

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if (otype == RSB_NUMERICAL_TYPE_DOUBLE)*dval = *VA?*VA: rsb__calloc(sizeof(double)*annz);
	else
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT
	if (otype == RSB_NUMERICAL_TYPE_FLOAT) *fval =*VA?*VA: rsb__calloc(sizeof(float) *annz);
	else
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */
	#ifdef RSB_NUMERICAL_TYPE_INT
	if (otype == RSB_NUMERICAL_TYPE_INT) *ival =*VA?*VA: rsb__calloc(sizeof(int) *annz);
	else
	#endif /* RSB_NUMERICAL_TYPE_INT */
	#ifdef RSB_NUMERICAL_TYPE_CHAR
	if (otype == RSB_NUMERICAL_TYPE_CHAR) *cval =*VA?*VA: rsb__calloc(sizeof(char) *annz);
	else
	#endif /* RSB_NUMERICAL_TYPE_CHAR */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	if (otype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX) *zval =*VA?*VA: rsb__calloc(sizeof(float complex) *annz);
	else
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	if (otype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX) *Zval =*VA?*VA: rsb__calloc(sizeof(double complex) *annz);
	else
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_PATTERN
	if (otype == RSB_NUMERICAL_TYPE_PATTERN){/* no zval allocation (TODO : document this) */}
	else
	#endif /* RSB_NUMERICAL_TYPE_PATTERN */
	{errval = RSB_ERR_UNSUPPORTED_TYPE;RSB_PERR_GOTO(err,RSB_ERRM_ES);}

	if(IA && JA)
	{
		if( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR )
		{
			RSB_WARN("RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR will break with non matching index types\n");
			/* 	Allocating slightly oversized arrays.
				FIXME : potential type/size mismatches here ! 
			*/
			if(!*IA)
				*IA  = rsb__calloc(sizeof(rsb_nnz_idx_t)*((annz>*m?annz:*m)+1));
			if(!*JA)
				*JA   = rsb__calloc(sizeof(rsb_nnz_idx_t)*(annz+1));
		}
		else
		{
			if(!*IA)
			*IA   = rsb__calloc(sizeof(rsb_coo_idx_t)*annz);
			if(!*JA)
			*JA   = rsb__calloc(sizeof(rsb_coo_idx_t)*annz);
		}
    		if( !*IA || !*JA)
			goto err;
	}

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if (otype == RSB_NUMERICAL_TYPE_DOUBLE) if(!dval) {RSB_PERR_GOTO(err,RSB_ERRM_ES)}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT
	if (otype == RSB_NUMERICAL_TYPE_FLOAT ) if(!fval) {RSB_PERR_GOTO(err,RSB_ERRM_ES)}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */
	#ifdef RSB_NUMERICAL_TYPE_INT
	if (otype == RSB_NUMERICAL_TYPE_INT   ) if(!ival) {RSB_PERR_GOTO(err,RSB_ERRM_ES)}
	#endif /* RSB_NUMERICAL_TYPE_INT */
	#ifdef RSB_NUMERICAL_TYPE_CHAR
	if (otype == RSB_NUMERICAL_TYPE_CHAR  ) if(!cval) {RSB_PERR_GOTO(err,RSB_ERRM_ES)}
	#endif /* RSB_NUMERICAL_TYPE_CHAR */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	if (otype == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  ) if(!zval) {RSB_PERR_GOTO(err,RSB_ERRM_ES)}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	if (otype == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  ) if(!Zval) {RSB_PERR_GOTO(err,RSB_ERRM_ES)}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_PATTERN
	/* :) */
	#endif /* RSB_NUMERICAL_TYPE_PATTERN */

    	if( IA && JA)
		goto full_scan;

#if RSB_WANT_OMPIO_SUPPORT
	{errval = RSB_ERR_UNIMPLEMENTED_YET;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
#endif /* RSB_WANT_OMPIO_SUPPORT */

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if (typecode == RSB_NUMERICAL_TYPE_DOUBLE)
	{
		double iv;
		if(rsb_mm_is_complex(matcode))
		for (i=0; i<*nnz; i++)
		{
			re += (rsb_zfscanf(fd,"%lg %lg\n",NULL,NULL,*dval+i,&(iv),ngzfd)==2);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
		else
		for (i=0; i<*nnz; i++)
		{
			re += (rsb_zfscanf(fd,"%lg\n",NULL,NULL,*dval+i,NULL,ngzfd)==1);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
	}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

	#ifdef RSB_NUMERICAL_TYPE_FLOAT
	if (typecode == RSB_NUMERICAL_TYPE_FLOAT)
	{
		float iv;
		if(rsb_mm_is_complex(matcode))
		for (i=0; i<*nnz; i++)
		{
			re += (rsb_zfscanf(fd,"%g %g\n",NULL,NULL,*fval+i,&(iv),ngzfd)==2);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
		else
		for (i=0; i<*nnz; i++)
		{
			re += (rsb_zfscanf(fd, "%g\n",NULL,NULL,*fval+i,NULL,ngzfd)==1);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
	}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */

	#ifdef RSB_NUMERICAL_TYPE_INT
	if (typecode == RSB_NUMERICAL_TYPE_INT)
	for (i=0; i<*nnz; i++)
	{
		double fv,iv;
		if(rsb_mm_is_complex(matcode))
		for (i=0; i<*nnz; i++)
		{
			re += (rsb_zfscanf(fd,"%lg %lg\n",NULL,NULL,&(fv),&(iv),ngzfd)==2);
			(*ival)[i] = (int)fv;
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
		else
		for (i=0; i<*nnz; i++)
		{
			re += (rsb_zfscanf(fd,"%lg\n",NULL,NULL,&fv,NULL,ngzfd)==1);
			(*ival)[i] = (int)fv;
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
	}
	#endif /* RSB_NUMERICAL_TYPE_INT */

	#ifdef RSB_NUMERICAL_TYPE_CHAR
	if (typecode == RSB_NUMERICAL_TYPE_CHAR)
	for (i=0; i<*nnz; i++)
	{
		double fv;
		re += (rsb_zfscanf(fd,"%g\n",NULL,NULL,&fv,NULL,ngzfd)==1);
		(*cval)[i] = (char)fv;
		RSB_IO_VERBOSE_MSG(i,*nnz);
	}
	#endif /* RSB_NUMERICAL_TYPE_CHAR */

	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	if (typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX)
	{
		if(rsb_mm_is_complex(matcode))
		for (i=0; i<*nnz; i++)
		{
			float rv,iv;
			re += (rsb_zfscanf(fd,"%g %g\n",NULL,NULL,&(rv),&(iv),ngzfd)==2);
			(*zval)[i] = (rv + I * iv);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
		else
		for (i=0; i<*nnz; i++)
		{
			float rv;
			re += (rsb_zfscanf(fd,"%g\n",NULL,NULL,&(rv),NULL,ngzfd)==1);
			(*zval)[i] = (rv + I * 0);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
	}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	if (typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX)
	{
		if(rsb_mm_is_complex(matcode))
		for (i=0; i<*nnz; i++)
		{
			double rv,iv;
			re += (rsb_zfscanf(fd,"%lg %lg\n",NULL,NULL,&(rv),&(iv),ngzfd)==2);
			(*Zval)[i] = (rv + I * iv);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
		else
		for (i=0; i<*nnz; i++)
		{
			double rv;
			re += (rsb_zfscanf(fd,"%lg\n",NULL,NULL,&(rv),NULL,ngzfd)==1);
			(*Zval)[i] = (rv + I * 0);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
	}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
	goto scan_done;
full_scan:

#if RSB_WANT_EXPERIMENTAL_BINARY_COO
	if( IA && JA && VA && *IA && *JA && *VA && (fd || ngzfd) )
	{
		const int cc = rsb_getc(fd, ngzfd);

		if(cc == 'B') // Binary..
		{
			if(RSB_SOME_ERROR(errval = rsb__read_coo_bin_fd(fd, ngzfd, *IA, *JA, *VA, *m, *k , *nnz, otype)))
			{
				RSB_PERR_GOTO(err,RSB_ERRM_ES);
			}
			else
			{
				re = *nnz;
				goto scan_done;
			}
		}
		else
			rsb_ungetc(cc,fd,ngzfd);
	}
#endif /* RSB_WANT_EXPERIMENTAL_BINARY_COO */

	/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
	/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
	/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if (typecode == RSB_NUMERICAL_TYPE_DOUBLE)
	{
		if(rsb_mm_is_complex(matcode))
#if RSB_WANT_OMPIO_SUPPORT
		{errval = RSB_ERR_UNIMPLEMENTED_YET;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
#else /* RSB_WANT_OMPIO_SUPPORT */
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			double iv;
			re += (rsb_zfscanf(fd,RSB_IJPC " %lg %lg\n",&iI,&iJ,*dval+i,&(iv),ngzfd)==4);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
	        	(*JA)[i]--;
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
#endif /* RSB_WANT_OMPIO_SUPPORT */
		else
#if RSB_WANT_OMPIO_SUPPORT
			rsb__ompio_DOUBLE(nnz,fd,ngzfd,dval,IA,JA,&re);
#else /* RSB_WANT_OMPIO_SUPPORT */
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			re += (rsb_zfscanf(fd,RSB_IJPC " %lg\n",&iI,&iJ,*dval+i,NULL,ngzfd)==3);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
        		(*JA)[i]--;
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
#endif /* RSB_WANT_OMPIO_SUPPORT */
	}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

	#ifdef RSB_NUMERICAL_TYPE_FLOAT
	if (typecode == RSB_NUMERICAL_TYPE_FLOAT)
	{
		if(rsb_mm_is_complex(matcode))
#if RSB_WANT_OMPIO_SUPPORT
		{errval = RSB_ERR_UNIMPLEMENTED_YET;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
#else /* RSB_WANT_OMPIO_SUPPORT */
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			float iv;
			re += (rsb_zfscanf(fd,RSB_IJPC " %g %g\n",&iI,&iJ,*fval+i,&(iv),ngzfd)==4);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
	        	(*JA)[i]--;
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
#endif /* RSB_WANT_OMPIO_SUPPORT */
		else
#if RSB_WANT_OMPIO_SUPPORT
			rsb_ompio_FLOAT(nnz,fd,ngzfd,dval,IA,JA,&re);
#else /* RSB_WANT_OMPIO_SUPPORT */
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			re += (rsb_zfscanf(fd,RSB_IJPC " %g\n",&iI,&iJ,*fval+i,NULL,ngzfd)==3);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
	        	(*JA)[i]--;
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
#endif /* RSB_WANT_OMPIO_SUPPORT */
	}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */

	#ifdef RSB_NUMERICAL_TYPE_INT
	if (typecode == RSB_NUMERICAL_TYPE_INT)
	{
#if RSB_WANT_OMPIO_SUPPORT
			rsb_ompio_INT(nnz,fd,ngzfd,dval,IA,JA,&re);
#else /* RSB_WANT_OMPIO_SUPPORT */
	if(rsb_mm_is_complex(matcode))
#if RSB_WANT_OMPIO_SUPPORT
	{errval = RSB_ERR_UNIMPLEMENTED_YET;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
#else /* RSB_WANT_OMPIO_SUPPORT */
	for (i=0; i<*nnz; i++)
	{
		rsb_coo_idx_t iI,iJ;
		double fv,iv;
		re += (rsb_zfscanf(fd,RSB_IJPC " %lg %lg\n",&iI,&iJ,&fv,&iv,ngzfd)==4);
		(*IA)[i] = iI;
		(*JA)[i] = iJ;
        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
        	(*JA)[i]--;
		(*ival)[i] = (int)fv;
		RSB_IO_VERBOSE_MSG(i,*nnz);
	}
#endif /* RSB_WANT_OMPIO_SUPPORT */
	else
	for (i=0; i<*nnz; i++)
	{
		rsb_coo_idx_t iI,iJ;
		double fv;
		re += (rsb_zfscanf(fd,RSB_IJPC " %lg\n",&iI,&iJ,&fv,NULL,ngzfd)==3);
		(*IA)[i] = iI;
		(*JA)[i] = iJ;
        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
        	(*JA)[i]--;
		(*ival)[i] = (int)fv;
		RSB_IO_VERBOSE_MSG(i,*nnz);
	}
#endif /* RSB_WANT_OMPIO_SUPPORT */
	}
	#endif /* RSB_NUMERICAL_TYPE_INT */

	#ifdef RSB_NUMERICAL_TYPE_CHAR
	if (typecode == RSB_NUMERICAL_TYPE_CHAR)
#if RSB_WANT_OMPIO_SUPPORT
			rsb_ompio_CHAR(nnz,fd,ngzfd,dval,IA,JA,&re);
#else /* RSB_WANT_OMPIO_SUPPORT */
	for (i=0; i<*nnz; i++)
	{
		rsb_coo_idx_t iI,iJ;
		double fv;
		re += (rsb_zfscanf(fd,RSB_IJPC " %g\n",&iI,&iJ,&fv,NULL,ngzfd)==3);
		(*IA)[i] = iI;
		(*JA)[i] = iJ;
        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
        	(*JA)[i]--;
		(*cval)[i] = (char)fv;
		RSB_IO_VERBOSE_MSG(i,*nnz);
	}
#endif /* RSB_WANT_OMPIO_SUPPORT */
	#endif /* RSB_NUMERICAL_TYPE_INT */

	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	if (typecode == RSB_NUMERICAL_TYPE_FLOAT_COMPLEX)
	{
		if(rsb_mm_is_complex(matcode))
#if RSB_WANT_OMPIO_SUPPORT
			rsb_ompio_FLOAT_COMPLEX(nnz,fd,ngzfd,dval,IA,JA,&re);
#else /* RSB_WANT_OMPIO_SUPPORT */
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			float rv,iv;
			re += (rsb_zfscanf(fd,RSB_IJPC " %g %g\n",&iI,&iJ,&(rv),&(iv),ngzfd)==4);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
	        	(*JA)[i]--;
			(*zval)[i] = (rv + I * iv);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
#endif /* RSB_WANT_OMPIO_SUPPORT */
		else
#if RSB_WANT_OMPIO_SUPPORT
		{errval = RSB_ERR_UNIMPLEMENTED_YET;RSB_PERR_GOTO(err,RSB_ERRM_ES);}
#else /* RSB_WANT_OMPIO_SUPPORT */
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			float rv;
			re += (rsb_zfscanf(fd,RSB_IJPC " %g\n",&iI,&iJ,&(rv),NULL,ngzfd)==3);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
	        	(*JA)[i]--;
			(*zval)[i] = (rv + I * 0);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
#endif /* RSB_WANT_OMPIO_SUPPORT */
	}
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */

	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	if (typecode == RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX)
	{
		if(rsb_mm_is_complex(matcode))
#if RSB_WANT_OMPIO_SUPPORT
			rsb__ompio_DOUBLE_COMPLEX(nnz,fd,ngzfd,dval,IA,JA,&re);
#else /* RSB_WANT_OMPIO_SUPPORT */
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			double rv,iv;
			re += (rsb_zfscanf(fd,RSB_IJPC " %lg %lg\n",&iI,&iJ,&(rv),&(iv),ngzfd)==4);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
	        	(*JA)[i]--;
			(*Zval)[i] = (rv + I * iv);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
#endif /* RSB_WANT_OMPIO_SUPPORT */
		else
		for (i=0; i<*nnz; i++)
		{
			rsb_coo_idx_t iI,iJ;
			double rv;
			re += (rsb_zfscanf(fd,RSB_IJPC " %lg\n",&iI,&iJ,&(rv),NULL,ngzfd)==3);
			(*IA)[i] = iI;
			(*JA)[i] = iJ;
	        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
	        	(*JA)[i]--;
			(*Zval)[i] = (rv + I * 0);
			RSB_IO_VERBOSE_MSG(i,*nnz);
		}
	}
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */

	#ifdef RSB_NUMERICAL_TYPE_PATTERN
	if (typecode == RSB_NUMERICAL_TYPE_PATTERN)
#if RSB_WANT_OMPIO_SUPPORT
			rsb__ompio_PATTERN(nnz,fd,ngzfd,IA,JA,&re);
#else /* RSB_WANT_OMPIO_SUPPORT */
	for (i=0; i<*nnz; i++)
	{
		rsb_coo_idx_t iI,iJ;
		re += (rsb_zfscanf(fd,RSB_IJPC "\n",&iI,&iJ,NULL,NULL,ngzfd)==2);
		(*IA)[i] = iI;
		(*JA)[i] = iJ;
        	(*IA)[i]--;  /* adjust from 1-based to 0-based */
        	(*JA)[i]--;
		RSB_IO_VERBOSE_MSG(i,*nnz);
	}
#endif /* RSB_WANT_OMPIO_SUPPORT */
	#endif /* RSB_NUMERICAL_TYPE_PATTERN */
	if( is_lowerp || is_upperp )
	{
		rsb_flags_t flags = RSB_FLAG_NOFLAGS;

		flags = rsb__do_detect_and_add_triangular_flags(*IA,*JA,*nnz,flags);
		if( is_lowerp )
			*is_lowerp = RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER);
		if( is_upperp )
			*is_upperp = RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER);
	}

	#ifdef RSB_NUMERICAL_TYPE_PATTERN	
	if(typecode == RSB_NUMERICAL_TYPE_PATTERN && otype!=typecode && VA)
		rsb__fill_with_ones(*VA,otype,*nnz,1);
	#endif /* RSB_NUMERICAL_TYPE_PATTERN */
scan_done:
	if ( fd !=stdin || ngzfd )
	{
		if(ngzfd)
			RSB_FCLOSE(ngzfd);
		else
			fclose(fd);
	}
	if(re!=*nnz)
	{
		/* FIXME : this can happen when reading as double a complex matrix, now. */
		RSB_STDERR("read only %zu out of %zd matrix elements (incomplete or not a matrix file?)!\n",(rsb_printf_int_t)re,(rsb_printf_int_t)(*nnz));
		goto err;
	}

	if( _m != _k && is_symmetric )
	{
		RSB_STDERR("matrix declared as symmetric but not square!\n");
		goto err;
	}

	if(!(IA && JA))
		goto afterpmtxchecks;
	if(is_symmetric && RSB_EXPERIMENTAL_EXPAND_SYMMETRIC_MATRICES_BY_DEFAULT)
	{
		/* FIXME: this breaks reallocation tricks !! */
		if( rsb__reallocate_with_symmetry(IA,JA,VA,nnz,(otype)) )
		{
			RSB_STDERR("problems handling matrix symmetry!\n");
			goto err;
		}
	}

	if(is_symmetric && !RSB_EXPERIMENTAL_EXPAND_SYMMETRIC_MATRICES_BY_DEFAULT)
	{
		rsb_bool_t has_diagonal_elements = RSB_BOOL_FALSE;

		if(rsb__util_coo_check_if_triangle_non_empty(*IA,*JA,*nnz,RSB_FLAG_UPPER))
		{
			if(wv)
				RSB_STDERR("#converting upper to lower triangle..\n");
			rsb__util_coo_upper_to_lower_symmetric(*IA,*JA,*nnz);
			if(is_lower)
				is_lower = RSB_BOOL_TRUE;
			if( is_lowerp )
				*is_lowerp = RSB_BOOL_TRUE;
			if( is_upperp )
				*is_upperp = RSB_BOOL_FALSE;
		}

		if(rsb__util_coo_check_if_triangle_non_empty(*IA,*JA,*nnz,RSB_FLAG_UPPER))
		{
			RSB_STDERR("input declared as symmetric, but it is unsymmetric!\n");
			goto err;
		}

		errval = rsb__util_coo_check_if_has_diagonal_elements(*IA,*JA,*nnz,*m,&has_diagonal_elements,wv);
		if(RSB_SOME_ERROR(errval))
		{
			RSB_STDERR("error while checking diagonal elements!\n");
			goto err;
		}

		if(!has_diagonal_elements)
		if(wv)
		{
			RSB_STDERR("Input has missing elements on the diagonal.\n");
		}
	}
afterpmtxchecks:
	frt += rsb_time();
	/* this should be a disk read only routine, not an output reporting one */
	if(0)
		RSB_IO_NOTICE("file I/O took %lf s (%lf nnz, %lf nnz/s ) \n",frt,((double)*nnz),(((double)*nnz)/frt));

	if(IA && JA)
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER))
		RSB_SWAP(rsb_coo_idx_t*,*IA,*JA);

	errval = RSB_ERR_NO_ERROR;
prerr:
	goto ret;
err:
#if RSB_20120321_IOBUFFERING
	RSB_CONDITIONAL_FREE(iobuf);
#endif /* RSB_20120321_IOBUFFERING */
	/* FIXME: this will free also already allocated arrays */
	if(!aia) if(IA)
	RSB_CONDITIONAL_FREE(*IA);
	if(!aja) if(JA)
	RSB_CONDITIONAL_FREE(*JA);
	if(!ava){
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	if(dval) RSB_CONDITIONAL_FREE(*dval);
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT
	if(fval) RSB_CONDITIONAL_FREE(*fval);
	#endif /* RSB_NUMERICAL_TYPE_FLOAT */
	#ifdef RSB_NUMERICAL_TYPE_INT
	if(ival) RSB_CONDITIONAL_FREE(*ival);
	#endif /* RSB_NUMERICAL_TYPE_INT */
	#ifdef RSB_NUMERICAL_TYPE_CHAR
	if(cval) RSB_CONDITIONAL_FREE(*cval);
	#endif /* RSB_NUMERICAL_TYPE_CHAR */
	#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
	if(zval) RSB_CONDITIONAL_FREE(*zval);
	#endif /* RSB_NUMERICAL_TYPE_FLOAT_COMPLEX */
	#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
	if(Zval) RSB_CONDITIONAL_FREE(*Zval);
	#endif /* RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX */
	}
	errval = RSB_ERR_GENERIC_ERROR;
ret:
	RSB_DO_ERR_RETURN(errval);
}
#undef RSB_IJPC
#endif /* RSB_WITH_MM */

rsb_err_t rsb__do_util_get_matrix_dimensions(const char * filename, size_t * cols, size_t * rows, size_t * nnzp, rsb_flags_t*flagsp)
{
	/**
	 * \ingroup gr_internals
	 *
	 * FIXME : needs error handling
	 */
	rsb_coo_idx_t m,k;
	rsb_nnz_idx_t nnz;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	rsb_bool_t is_symmetric = RSB_BOOL_FALSE;
	rsb_bool_t is_hermitian = RSB_BOOL_FALSE;
	rsb_bool_t is_pattern = RSB_BOOL_FALSE;
	rsb_bool_t is_lower = RSB_BOOL_FALSE;
	rsb_bool_t is_upper = RSB_BOOL_FALSE;
	rsb_bool_t is_vector = RSB_BOOL_FALSE;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	errval = rsb__util_mm_info_matrix_f(filename,&m,&k,&nnz,&typecode,&is_symmetric,&is_hermitian,&is_pattern,&is_lower,&is_upper,&is_vector);
	if(cols)
		*cols = k;
	if(rows)
		*rows = m;
	if(nnzp)
		*nnzp = nnz;
	if(is_symmetric)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_SYMMETRIC);
	if(is_hermitian)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_HERMITIAN);
	/* if(is_pattern) ... */
	if(is_lower)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
	if(is_upper)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
	if(flagsp)
		*(flagsp) = flags;
	RSB_DO_ERR_RETURN(errval)
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__util_mm_load_matrix_f_as_csr(const char *filename, rsb_nnz_idx_t ** INDX, rsb_coo_idx_t ** JA, void **VA, rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t typecode, rsb_flags_t flags/*, rsb_bool_t *is_lowerp, rsb_bool_t *is_upperp*/)
{
	/**
	 * \ingroup gr_internals
	 * FIXME : should optimize
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t coo;

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);
	*INDX = NULL;
	*JA = NULL;
	*VA = NULL;
	RSB_BZERO_P(&coo);
	coo.typecode = typecode;
       	errval = rsb__util_mm_load_coo_matrix(filename,&coo);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb__util_sort_row_major_inner(coo.VA,coo.IA,coo.JA,coo.nnz,coo.nr,coo.nc,coo.typecode,flags);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb__allocate_csr_arrays_from_coo_sorted(coo.VA,coo.IA,coo.JA,coo.nnz,coo.nr,coo.nc,coo.typecode,VA,JA,INDX);
	if(RSB_SOME_ERROR(errval))
		goto err;
	*m = coo.nr; *k = coo.nc; *nnz = coo.nnz;
err:
	rsb__destroy_coo_matrix_t(&coo);
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__util_mm_load_matrix_f_as_csc(const char *filename, rsb_nnz_idx_t ** INDX, rsb_coo_idx_t ** IA, void **VA, rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t typecode, rsb_flags_t flags/*, rsb_bool_t *is_lowerp, rsb_bool_t *is_upperp*/)
{
	/** 
	 * \ingroup gr_internals
	 * FIXME : should optimize
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t coo;

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);
	*INDX = NULL;*IA = NULL;*VA = NULL;
	RSB_BZERO_P(&coo);
	coo.typecode = typecode;
       	errval = rsb__util_mm_load_coo_matrix(filename,&coo);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb_util_sort_column_major(coo.VA,coo.IA,coo.JA,coo.nnz,coo.nr,coo.nc,coo.typecode,flags);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb__allocate_csc_arrays_from_coo_sorted(coo.VA,coo.IA,coo.JA,coo.nnz,coo.nr,coo.nc,coo.typecode,VA,IA,INDX);
	if(RSB_SOME_ERROR(errval))
		goto err;
	*m = coo.nr; *k = coo.nc; *nnz = coo.nnz;
err:
	rsb__destroy_coo_matrix_t(&coo);
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__util_mm_fill_arrays_for_csc(const char *filename, rsb_nnz_idx_t * INDX, rsb_coo_idx_t * IA, void *VA, rsb_type_t typecode, rsb_flags_t flags)
{
	/** 
	 * \ingroup gr_internals
	 * FIXME : should optimize
	 * */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_coo_mtx_t coo;
	rsb_nnz_idx_t *iINDX = NULL;
	rsb_coo_idx_t *iIA = NULL;
	void *iVA = NULL;

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);
	if(!INDX || !IA || !VA)
		return RSB_ERR_BADARGS;

	RSB_BZERO_P(&coo);
	coo.typecode = typecode;
       	errval = rsb__util_mm_load_coo_matrix(filename,&coo);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb_util_sort_column_major(coo.VA,coo.IA,coo.JA,coo.nnz,coo.nr,coo.nc,coo.typecode,flags);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb__allocate_csc_arrays_from_coo_sorted(coo.VA,coo.IA,coo.JA,coo.nnz,coo.nr,coo.nc,coo.typecode,&iVA,&iIA,&iINDX);
	if(RSB_SOME_ERROR(errval))
		goto err;
	errval = rsb__copy_css_arrays(iVA,iINDX,iIA,coo.nnz,coo.nc,typecode,VA,INDX,IA);
	if(RSB_SOME_ERROR(errval))
		goto err;
err:
	rsb__destroy_coo_matrix_t(&coo);
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

size_t rsb__sys_filesize(const char *filename)
{
	/**
	 * file size, in bytes
	 * FIXME
	 * TODO : to sys.c
	 * TODO: rsb__do_util_get_matrix_dimensions shall invoke this.
	 * */
#ifdef RSB_HAVE_SYS_STAT_H
	struct stat ss;
#endif /* RSB_HAVE_SYS_STAT_H */
	if(!filename)
		goto err;
#ifdef RSB_HAVE_SYS_STAT_H
	stat(filename,&ss);
	return ss.st_size;
#else /* RSB_HAVE_SYS_STAT_H */
#endif /* RSB_HAVE_SYS_STAT_H */
err:
	return 0;
}

rsb_err_t rsb__do_file_mtx_get_dims(const char * filename, rsb_coo_idx_t* nrp, rsb_coo_idx_t *ncp, rsb_coo_idx_t *nzp, rsb_flags_t*flagsp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t nrA = 0,ncA = 0,nzA = 0;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

	errval = rsb__do_util_get_matrix_dimensions(filename,&ncA,&nrA,&nzA,&flags);

	if(RSB_SOME_ERROR(errval))
		goto ret;

	/* RSB_STDOUT("%d / %d  %d / %d  %d / %d\n",nrA,RSB_MAX_MATRIX_DIM,ncA,RSB_MAX_MATRIX_DIM,nnzA,RSB_MAX_MATRIX_NNZ); */
	if( nrp && RSB_INVALID_COO_INDEX(nrA) )
		errval |= RSB_ERR_LIMITS;
	if( ncp && RSB_INVALID_COO_INDEX(ncA) )
		errval |= RSB_ERR_LIMITS;
	if( nzp && RSB_INVALID_NNZ_INDEX(nzA) )
		errval |= RSB_ERR_LIMITS;
	if(nrp)
		*nrp = (rsb_coo_idx_t)nrA;
	if(ncp)
		*ncp = (rsb_coo_idx_t)ncA;
	if(nzp)
		*nzp = (rsb_nnz_idx_t)nzA;
	if(flagsp)
		*flagsp = (rsb_flags_t)flags;
	/* FIXME: need overflow check here */
ret:
	return errval;
}

/* @endcond */
