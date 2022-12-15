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
 @file
 @author Michele Martone
 @brief
 This is the main program used to benchmark and test our library.
 This should be the swiss army knife program for our library.
 */
/*
  This is not an example program: to build it, one needs all of the internal library headers.
 */

#include "rsb_common.h"
#include "rsb.h"
#include "rsb_test_matops.h"
#include "rsb_failure_tests.h"
#include "rsb_internals.h"
#if RSB_WITH_SPARSE_BLAS_INTERFACE 
#include "rsb_libspblas_handle.h"
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE  */
#include "rsb_libspblas_tests.h"
//#include "rsb-config.h"
#if RSB_WANT_ACTION
#include <signal.h>
#if defined(RSB_WANT_ACTION_SIGNAL)
#else /* defined(RSB_WANT_ACTION_SIGNAL) */
#include <bits/sigaction.h>
#endif /* defined(RSB_WANT_ACTION_SIGNAL) */
#endif /* RSB_WANT_ACTION */
/* #include <string.h> */	/* strdup is here, but not in C99 */

#define RSB_WANT_PERMISSIVE_RSBENCH 1
#define RSB_SHALL_UPDATE_COMPLETEBENCHS 0 /* Old, pretty complete complete coverage reference benchmark, also useful as test. */
#define RSB_WANT_REDUCED_RSB_M4_MATRIX_META_OPS 1 


#if RSB_WITH_LIKWID
#define RSB_RSBENCH_EXEC(FEXP) {rsb_err_t errval;RSB_LIKWID_MARKER_INIT;errval=RSB_ERR_TO_PROGRAM_ERROR(FEXP);RSB_LIKWID_MARKER_EXIT;return errval;}
#else /* RSB_WITH_LIKWID */
#define RSB_RSBENCH_EXEC(FEXP) {return RSB_ERR_TO_PROGRAM_ERROR(FEXP);}
#endif /* RSB_WITH_LIKWID */

RSB_INTERNALS_COMMON_HEAD_DECLS /* for RSB_INFO */
#if RSB_WANT_ACTION
	int rsb__quit_rsbench;
#if defined(RSB_WANT_ACTION_SIGNAL)
#else /* defined(RSB_WANT_ACTION_SIGNAL) */
	struct sigaction rsb_osa;
#endif /* defined(RSB_WANT_ACTION_SIGNAL) */

void rsb__sigh(int signal)
{
	/* TODO: extend this mechanism optionally to the library itself. */

	if( rsb__quit_rsbench == 0 )
	{
		RSBENCH_STDOUT("\n");
		RSBENCH_STDOUT("====================================================\n");
		RSBENCH_STDOUT("Caught signal %d: will terminate as soon as possible.\n",signal);
		RSBENCH_STDOUT("  ( next time won't catch the signal anymore ).\n");
		RSBENCH_STDOUT("====================================================\n");
		RSBENCH_STDOUT("\n");
		rsb__quit_rsbench++;
	}
	else
	if( rsb__quit_rsbench == 1 )
	{
#if defined(RSB_WANT_ACTION_SIGNAL)
#else /* defined(RSB_WANT_ACTION_SIGNAL) */
		sigaction(SIGINT,&rsb_osa,NULL);
#endif /* defined(RSB_WANT_ACTION_SIGNAL) */
	}
}

void rsb__sigr(void)
{
	rsb__quit_rsbench = 0;
	{
#if RSB_WANT_ACTION_SIGNAL
		/* signal() is part of C99 */
		signal(SIGINT,&rsb__sigh); /* not to be called from a threaded environment ... */
#else /* RSB_WANT_ACTION_SIGNAL */
		/* sigaction() is part of POSIX, not part of C99 */
		struct sigaction act;
		RSB_BZERO_P(&act);
		RSB_BZERO_P(&rsb_osa);
		act.sa_handler  = rsb__sigh;
		sigemptyset(&act.sa_mask);
    		sigaction(SIGINT, &act,&rsb_osa);
/*
		sigaction(SIGUSR1, &act, &rsb_osa);
		sigaction(SIGUSR2, &act, &rsb_osa);

		sigaction(SIGQUIT,&act,&rsb_osa);
		sigaction(SIGTERM,&act,&rsb_osa);

		sigaction(SIGABRT,&act,&rsb_osa);
		sigaction(SIGTSTP,&act,&rsb_osa);

		sigaction(SIGBUS, &act,&rsb_osa);
		sigaction(SIGILL, &act,&rsb_osa);
	    	sigaction(SIGSEGV,&act,&rsb_osa);
*/
#endif /* RSB_WANT_ACTION_SIGNAL */
	}
}
#endif /* RSB_WANT_ACTION */

const char * rsb__strstr_end(const char * haystack, const char * needle)
{
	size_t nl = strlen(needle);
	size_t hl = strlen(haystack);

	if( hl < nl)
		return NULL;

	return strstr(haystack+(hl-nl),needle);
}

static rsb_err_t rsb__matrix_market_filename(const rsb_char_t * filename, rsb_flags_t faflags)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( ! filename )
	{
		errval = RSB_ERR_GENERIC_ERROR;
		goto ret;
	}

	if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_CHKFNP ) )
	{
#if RSB_WANT_EXPERIMENTAL_BINARY_COO
		if( rsb__strstr_end(filename,".mtx.bin") )
			goto ret;
		if( rsb__strstr_end(filename,".mtx.bin.gz") )
			goto ret;
		if( rsb__strstr_end(filename,".mtx.gz.bin") )
			goto ret;
		if( rsb__strstr_end(filename,".mtx.gz.bin.gz") )
			goto ret;
#endif /* RSB_WANT_EXPERIMENTAL_BINARY_COO */
		if( rsb__strstr_end(filename,".mtx") )
			goto ret;
		if( rsb__strstr_end(filename,".mtx.gz") )
			goto ret;
		errval = RSB_ERR_GENERIC_ERROR;
	}
ret:
	return errval;
}

static rsb_bool_t rsb__fn_same_mtx_basename( const rsb_char_t * matrixpath_1, const rsb_char_t * matrixpath_2)
{
	const rsb_char_t * p_1=matrixpath_1, * p_2=matrixpath_2;

	if( ! ( *matrixpath_1 && *matrixpath_2 ))
		goto ret_false;

	while( *p_1 && *p_2 && *p_1 == *p_2 )
		++p_1,++p_2;

	if ( *p_1 == *p_2 )
		goto ret_true;

	if ( *p_1 && 0 == strcmp(p_1,".gz") )
		goto ret_true;

	if ( *p_2 && 0 == strcmp(p_2,".gz") )
		goto ret_true;

ret_false:
	return RSB_BOOL_FALSE;
ret_true:
	return RSB_BOOL_TRUE;
}

static rsb_bool_t rsb__check_mtx_vec_file(const rsb_char_t * filename)
{
	// FIXME: need is_vector .. 
	rsb_bool_t is_vector = RSB_BOOL_FALSE;

	if(RSB_SOME_ERROR(rsb__util_mm_info_matrix_f(filename,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&is_vector) ))
		goto err;
	if( is_vector )
		return RSB_BOOL_TRUE;
err:
	return RSB_BOOL_FALSE;
}

static rsb_bool_t rsb__check_presence(rsb_char_t ** filenameap, rsb_int_t filenamen, const rsb_char_t * matrixpath, rsb_flags_t faflags)
{
	/* note that this is not meant to be efficient on really long lists  */
	rsb_bool_t there = RSB_BOOL_TRUE;
	rsb_int_t filenamei = 0;

	for ( filenamei = 0; filenamei < filenamen ; ++filenamei )
	{
		if( rsb__fn_same_mtx_basename(rsb__basename(filenameap[filenamei]),rsb__basename(matrixpath)) )
		{
			// FIXME
			if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_VRBADD ) )
				RSB_STDOUT("Skipping duplicate matrix file: %s\n",matrixpath);
			goto ret;
		}
	}
	there = RSB_BOOL_FALSE;
ret:
	return there;
}

static rsb_char_t * rsb__strdup(const rsb_char_t * s)
{
	// NB: strdup is not C99
	// using malloc: this is meant to be called before rsb_lib_init().
	size_t sl = strlen(s);
	rsb_char_t * d = NULL;

	if(sl && ( d = malloc(sl+1) ) != NULL )
		return strcpy(d,s);
	return d;
}

static rsb_err_t rsb__add_file(rsb_char_t ** filenameap, rsb_int_t * filenamenp, const rsb_char_t * matrixpath, rsb_flags_t faflags)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_CHKDUP ) )
		if( rsb__check_presence(filenameap, *filenamenp, matrixpath, faflags) == RSB_BOOL_TRUE )
		{
			goto ret;
		}

	if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_CHKMTX ) )
		if( rsb__check_mtx_vec_file(matrixpath) == RSB_BOOL_TRUE )
		{
			RSB_STDOUT("Skipping non matrix file: %s\n",matrixpath);
			goto ret;
		}

	if( *filenamenp >= RSB_RSBENCH_MAX_MTXFILES )
	{
		errval = RSB_ERR_LIMITS;
		RSB_STDERR("Reached a limit of %d files! \n",RSB_RSBENCH_MAX_MTXFILES);
		goto ret;
	}
	
	if ( RSB_ERR_NO_ERROR == ( errval = rsb__matrix_market_filename(matrixpath,faflags) ) )
	{
		char aname[RSB_MAX_FILENAME_LENGTH];
		snprintf(aname,sizeof(aname),"%s",matrixpath);

		if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_CHKGSS ) )
		if( rsb__file_exists(aname) != RSB_BOOL_TRUE )
		{
#if RSB_WANT_ZLIB_SUPPORT
			if( rsb__strstr_end(aname,".mtx") )
				snprintf(aname,sizeof(aname),"%s.gz",matrixpath);
			else
#endif /* RSB_WANT_ZLIB_SUPPORT */
			if( rsb__strstr_end(aname,".mtx.gz") )
				aname[strlen(aname)-3]=RSB_NUL;
			if( rsb__file_exists(aname) )
				RSB_STDOUT("Given filename %s doesn't exist; assuming you meant %s, which exists!\n",matrixpath,aname);
			else
				snprintf(aname,sizeof(aname),"%s",matrixpath);
		}

		if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_CHKEXS ) )
		if( ! rsb__file_exists(aname) )
		{
			RSB_STDOUT("Skipping non-matrix filepath: %s\n",aname);
			goto ret;
		} /* might extend mechanism more ... */

		if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_VRBADD ) )
			RSB_STDOUT("Adding matrix file: %s\n",aname);
		// RSB_STDOUT("Adding %d-thfile %s\n",*filenamenp,aname);
		filenameap[(*filenamenp)++] = rsb__strdup( aname); //FIXME: this is a leak :-)
		// filenameap[(*filenamenp)++] = aname;
		errval = RSB_ERR_NO_ERROR;
	}
	else
	{
		RSB_STDOUT("Not adding %s: neither a directory, nor named as Matrix Market file.\n",matrixpath);
		errval = RSB_ERR_NO_ERROR;
	}
ret:
	return errval;
}

#if defined(RSB_HAVE_DIRENT_H) && defined(RSB_HAVE_SYS_TYPES_H)

#define _SVID_SOURCE
#include <dirent.h>	// FIXME: move out of here
#include <sys/types.h>	// FIXME: move out of here
#define RSB_USE_SCANDIR 0	// FIXME: unfortunately scandir is POSIX.1-2008 and not C99 (and separate compilation is not yet in librsb)

static rsb_err_t rsb__adddir_(rsb_char_t ** filenameap, rsb_int_t * filenamenp, const rsb_char_t * matrixpath, rsb_flags_t faflags)
{
	/* FIXME: unfinished, experimental */
	/* FIXME: declarations in rsb_common.h */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_USE_SCANDIR
	struct dirent **namelist=NULL;
	int n;
#endif /* RSB_USE_SCANDIR */
	DIR *dir = NULL;
	struct dirent *de = NULL;

	if(!filenameap || !filenamenp || !matrixpath)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	faflags|=RSB_FAF_VRBSRC;
	if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_VRBSRC ) )
		RSB_STDOUT("Will try %s\n",matrixpath);

	if ( ! RSB_DO_FLAG_HAS(faflags, RSB_FAF_CHKREC ) )
		goto addasfile;

	if ( NULL == ( dir = opendir(matrixpath) ) )
		goto addasfile;

	if ( RSB_DO_FLAG_HAS(faflags, RSB_FAF_VRBSRC ) )
		RSB_STDOUT("Will recurse in %s\n",matrixpath);

#if  RSB_USE_SCANDIR
	closedir(dir);

	if( (n=scandir(matrixpath, &namelist, NULL, alphasort)) > 0)
	{

		while(n--)
		{
			rsb_char_t dirname[RSB_MAX_FILENAME_LENGTH];

			de = namelist[n];
			sprintf(dirname,"%s%s%s",matrixpath,"/",de->d_name);

			if( strcmp(de->d_name,".") != 0 && strcmp(de->d_name,".." ) != 0 )
				errval = rsb__adddir_(filenameap, filenamenp, &dirname[0], faflags|RSB_FAF_CHKFNP );
		}
		free(namelist);
	}
#else /* RSB_USE_SCANDIR */
	while( ( de = readdir(dir) ) != NULL && errval != RSB_ERR_LIMITS )
	{
		rsb_char_t dirname[RSB_MAX_FILENAME_LENGTH];

		sprintf(dirname,"%s%s%s",matrixpath,"/",de->d_name);

		if( de->d_name[0] != '.' ) /* avoid hidden paths */
		if( strcmp(de->d_name,".") != 0 && strcmp(de->d_name,".." ) != 0 )
			errval = rsb__adddir_(filenameap, filenamenp, &dirname[0], faflags|RSB_FAF_CHKFNP );
	}

	closedir(dir);
#endif /* RSB_USE_SCANDIR */
	goto err;
addasfile:
	errval = rsb__add_file(filenameap, filenamenp, matrixpath, faflags);
err:
	return errval;
}

rsb_err_t rsb__adddir(rsb_char_t ** filenameap, rsb_int_t * filenamenp, const rsb_char_t * matrixpath, rsb_flags_t faflags)
{
	return rsb__adddir_(filenameap, filenamenp, matrixpath, faflags);
}
#else /* RSB_HAVE_SYS_TYPES_H && RSB_HAVE_DIRENT_H */
rsb_err_t rsb__adddir(rsb_char_t ** filenameap, rsb_int_t * filenamenp, const rsb_char_t * matrixpath, rsb_flags_t faflags)
{
	return rsb__add_file(filenameap, filenamenp, matrixpath, faflags);
}
#endif /* RSB_HAVE_SYS_TYPES_H && RSB_HAVE_DIRENT_H */

rsb_err_t rsb__print_configuration_string_rsbench(const char *pn, rsb_char_t * cs, rsb_bool_t wci)
{
	/* TODO: output buffer length check */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!cs)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}
	errval = rsb__print_configuration_string(pn, cs, wci);
	if(wci == RSB_BOOL_FALSE)
		goto err;

#if RSB_USE_LIBRSBPP
	{
		extern void rsb_rsbpp_dummy(void); // rsb_dummy.cpp
		rsb_rsbpp_dummy();
	}
	sprintf(cs+strlen(cs),"LIBRSBPP support: on.\n");
#else /* RSB_USE_LIBRSBPP */
	sprintf(cs+strlen(cs),"LIBRSBPP support: off.\n");
#endif /* RSB_USE_LIBRSBPP */

#if RSB_WANT_MKL
	rsb__print_mkl_version(cs+strlen(cs), RSB_BOOL_TRUE);
	sprintf(cs+strlen(cs),"\n");
#else /* RSB_WANT_MKL */
	sprintf(cs+strlen(cs),"MKL support: off.\n");
#endif /* RSB_WANT_MKL */

#if RSB_WANT_OMP_RECURSIVE_KERNELS
	sprintf(cs+strlen(cs),"OpenMP support: on.\n");
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	sprintf(cs+strlen(cs),"OpenMP support: off.\n");
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */

#if RSB_WANT_ARMPL
	sprintf(cs+strlen(cs),"ARMPL support: on.\n");
#else /* RSB_WANT_ARMPL */
	sprintf(cs+strlen(cs),"ARMPL support: off.\n");
#endif /* RSB_WANT_ARMPL */

#if RSB_WANT_XDR_SUPPORT
	sprintf(cs+strlen(cs),"XDR support: on.\n");
#else /* RSB_WANT_XDR_SUPPORT */
	sprintf(cs+strlen(cs),"XDR support: off.\n");
#endif /* RSB_WANT_XDR_SUPPORT */

#if RSB_WANT_ZLIB_SUPPORT
	sprintf(cs+strlen(cs),"ZLIB support: on.\n");
#else /* RSB_WANT_ZLIB_SUPPORT */
	sprintf(cs+strlen(cs),"ZLIB support: off.\n");
#endif /* RSB_WANT_ZLIB_SUPPORT */
#if RSB_WANT_EXPERIMENTAL_BINARY_COO
	sprintf(cs+strlen(cs),"Binary I/O Matrix Market hack: on.\n");
#else /* RSB_WANT_EXPERIMENTAL_BINARY_COO */
	sprintf(cs+strlen(cs),"Binary I/O Matrix Market hack: off.\n");
#endif /* RSB_WANT_EXPERIMENTAL_BINARY_COO */
#if RSB_USE_ASSERT
	sprintf(cs+strlen(cs),"Assertions: on.\n");
#else /* RSB_USE_ASSERT */
	sprintf(cs+strlen(cs),"Assertions: off.\n");
#endif /* RSB_USE_ASSERT */
#if RSB_ALLOW_INTERNAL_GETENVS
	sprintf(cs+strlen(cs),"Internal environment variables: on.\n");
#else /* RSB_ALLOW_INTERNAL_GETENVS */
	sprintf(cs+strlen(cs),"Internal environment variables: off.\n");
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
err:
	return errval;
}

#if RSB_HAVE_SETENV
void rsb__setenv(const rsb_char_t * var_val)
{
	rsb_char_t * esp = NULL;
	rsb_char_t * ava = rsb__strdup(var_val);

	if( ava && *ava && ( esp = strchr(ava,'=') ) != NULL )
	{
		const rsb_char_t *name=ava, *value=esp+1;
		*esp=RSB_NUL;
		setenv(name, value, /*overwrite*/1);
		RSBENCH_STDOUT("# Calling setenv() with arguments %s and %s\n",name,value);
	}
	if( ava )
		free(ava);
}
#endif /* RSB_HAVE_SETENV */

static int rsb__main_help(const int argc, char * const argv[], int default_program_operation, const char * program_codes, struct option *options)
{
			const char * pbn = rsb__basename(argv[0]);

			//RSB_STDOUT(
			printf(
				/*"[OBSOLETE DOCUMENTATION] \n"*/
				"Usage: %s [--bench] [OPTIONS] \n"
				"  or:  %s [ -o OPCODE] [ -O {subprogram-code}] [ {subprogram-specific-arguments} ] \n"
				"%s "RSB_INFOMSG_SAK"."
				"\n"
				"\n"
				//"\tOne may choose {option} among:\n"
				//"\t-I for getting system information and some micro benchmarking\n"
				"\t\n"
				 "Choose {subprogram-code} among:\n\n"
				"\tr for the reference benchmark (will produce a machine specific file)\n\n"
				"\tc for the complete benchmark\n\n"
				"\te for the matrix experimentation code\n\n"
				"\td for a single matrix dumpout\n\n"
				"\tb for the (current, going to be obsoleted) benchmark\n\n"
				"\tt for some matrix construction tests\n\n"
				"\to obsolete, will soon be removed\n"
				"\n"
				 "{subprogram-specific-arguments} will be available from the subprograms.\n\n"
				"\te.g.: %s      -O b -h   will show the current benchmark subprogram's options\n\n"
				"\te.g.: %s -o a -O b -h   will show the spmv     benchmark subprogram's options\n\n"
				"\te.g.: %s -o n -O b -h   will show the negation benchmark subprogram's options\n\n"
//				"\te.g.: %s -o A -O b    will run all of the benchmark programs.\n"
				"\nThe default {subprogram-code} is '%c'\n"
				"\n\tWith OPCODE among '%s'\n"/* TODO: fix this description, as it is too laconic. */
				"\n"
				,pbn
				,pbn
				,pbn
				,pbn
				,pbn
				,pbn
				,default_program_operation
				,program_codes
				);
			if(options)
				rsb_test_help_and_exit(pbn,options,0);
	return 0;
}

static rsb_err_t rsb__internals_check(void)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	
	int bxl;
	rsb_blk_idx_t *bxv = NULL;

	RSB_INFO("INTERNALS TEST: BEGIN (IGNORE THE ERROR PRINTOUT HERE BELOW, IT'S PART OF THE TEST)\n");

	if(rsb__util_atoi("-2K")!=-2)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(rsb__util_atoi_km2("1K")!=1024)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(rsb__util_atoi_km10("-2K")!=-2000)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("0",&bxl,&bxv);
	if(!bxv || bxl!=1 || bxv[0]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("0,90",&bxl,&bxv);
	if(!bxv || bxl!=2 || bxv[0]!=0 || bxv[1]!=90)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("1",&bxl,&bxv);
	if(!bxv || bxl!=1 || bxv[0]!=1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("1K",&bxl,&bxv);
	if(!bxv || bxl!=1 || bxv[0]!=1000)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("1,2",&bxl,&bxv);
	if(!bxv || bxl!=2 || bxv[0]!=1 || bxv[1]!=2)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("1k,2",&bxl,&bxv);
	if(!bxv || bxl!=2 || bxv[0]!=1000 || bxv[1]!=2)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("9,2M,3",&bxl,&bxv);
	if(!bxv || bxl!=3 || bxv[0]!=9 || bxv[1]!=2000000 || bxv[2]!=3)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("9,2M",&bxl,&bxv);
	if(!bxv || bxl!=2 || bxv[0]!=9 || bxv[1]!=2000000)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("-9,2M",&bxl,&bxv);
	if(!bxv || bxl!=2 || bxv[0]!=-9 || bxv[1]!=2000000)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_bx_array("9,-2M",&bxl,&bxv);
	if(!bxv || bxl!=2 || bxv[0]!=9 || bxv[1]!=-2000000)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL; bxl=0;
	errval = rsb__util_get_tn_array(":",&bxl,&bxv);
	if(!bxv || bxl<1 || bxv[0]<1 || bxv[bxl-1]<bxv[0] || bxv[bxl-1]>RSB_CONST_MAX_SUPPORTED_THREADS)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	RSB_CONDITIONAL_FREE(bxv);

	bxv = NULL;
	bxl = 0;
	errval = rsb__util_get_tn_array("1,2",&bxl,&bxv);
	if( bxl || bxv || errval == RSB_ERR_NO_ERROR )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	errval = RSB_ERR_NO_ERROR;

	// unit test for rsb__do_transposition_from_char
	errval = RSB_ERR_INTERNAL_ERROR;
	if(rsb__do_transposition_from_char('N')!=RSB_TRANSPOSITION_N)
		goto err;
	if(rsb__do_transposition_from_char('n')!=RSB_TRANSPOSITION_N)
		goto err;
	if(rsb__do_transposition_from_char('T')!=RSB_TRANSPOSITION_T)
		goto err;
	if(rsb__do_transposition_from_char('t')!=RSB_TRANSPOSITION_T)
		goto err;
	if(rsb__do_transposition_from_char('C')!=RSB_TRANSPOSITION_C)
		goto err;
	if(rsb__do_transposition_from_char('c')!=RSB_TRANSPOSITION_C)
		goto err;
	if(rsb__do_transposition_from_char('H')!=RSB_INVALID_TRANS)
		goto err;
	if(rsb__do_transposition_from_char('h')!=RSB_INVALID_TRANS)
		goto err;
	errval = RSB_ERR_NO_ERROR;


	// unit test for rsb__do_flip_uplo_flags
	errval = RSB_ERR_INTERNAL_ERROR;
	if(rsb__do_flip_uplo_flags(RSB_FLAG_UPPER)!=RSB_FLAG_LOWER)
		goto err;
	if(rsb__do_flip_uplo_flags(RSB_FLAG_LOWER)!=RSB_FLAG_UPPER)
		goto err;
	if(rsb__do_flip_uplo_flags(RSB_FLAG_UPPER_TRIANGULAR)!=RSB_FLAG_LOWER_TRIANGULAR)
		goto err;
	if(rsb__do_flip_uplo_flags(RSB_FLAG_LOWER_TRIANGULAR)!=RSB_FLAG_UPPER_TRIANGULAR)
		goto err;
	errval = RSB_ERR_NO_ERROR;

	RSB_INFO("INTERNALS TEST: END\n");
	goto ret;
err:
	if(!RSB_SOME_ERROR(errval))
		errval = RSB_ERR_INTERNAL_ERROR;
	rsb_perror(NULL,errval);
ret:
	RSB_CONDITIONAL_FREE(bxv);
	return errval;
}

static rsb_err_t rsb__chk_srt2p(void)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;

	RSBENCH_STDOUT("SORT CHECK: BEGIN\n");

#define RSB_EXPECT_2P_SORT_BUG 0
{
	// trigger double pass
	const struct rsb_mtx_partitioning_info_t * pinfop = NULL;
	const rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	const enum rsb_op_flags_t op_flags = RSB_OP_FLAG_DEFAULT;
	const int bs = RSB_DEFAULT_BLOCKING;
	const size_t wb = 0;
	const int br = bs, bc = bs;
	void * WA = NULL;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 2;
	const rsb_coo_idx_t nr = 1000000;
	const rsb_coo_idx_t nc = 1000000;
	rsb_coo_idx_t IA[] = { 463325, 417887};
	rsb_coo_idx_t JA[] = { 463231, 417880};
	rsb_coo_idx_t rIA[] = { 463325, 417887};
	rsb_coo_idx_t rJA[] = { 463231, 417880};
	RSB_DEFAULT_TYPE VA[] = {1,1};
	RSB_DEFAULT_TYPE rVA[] = {1,1};

	errval = rsb__do_index_based_bcsr_sort( IA, JA, VA, rIA, rJA, rVA, nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);
#if RSB_EXPECT_2P_SORT_BUG
	if( rIA[0] == 417887)
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rJA[0] == 417880)
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
#else
	if( rIA[0] != 417887)
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rJA[0] != 417880)
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
#endif
}

#define RSB_EXPECT_IP2P_SORT_BUG 0
{
	// trigger double pass
	const struct rsb_mtx_partitioning_info_t * pinfop = NULL;
	const rsb_flags_t flags = RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT;
	const enum rsb_op_flags_t op_flags = RSB_OP_FLAG_DEFAULT;
	const int bs = RSB_DEFAULT_BLOCKING;
	const size_t wb = 0;
	const int br = bs, bc = bs;
	void * WA = NULL;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 7;
	const rsb_coo_idx_t nr = 2*RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t);
	const rsb_coo_idx_t nc = 2*RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t);
	rsb_coo_idx_t rIA[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t rJA[] = {0,1,2,3,4,5,0};
	RSB_DEFAULT_TYPE rVA[] = {1,1,1,1,1,1,2};
	rsb_coo_idx_t IA[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t JA[] = {0,1,2,3,4,5,0};
	RSB_DEFAULT_TYPE VA[] = {1,1,1,1,1,1,2};

	errval = rsb__do_index_based_bcsr_sort( IA, JA, VA, rIA, rJA, rVA, nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

#if RSB_EXPECT_IP2P_SORT_BUG
	if( rIA[1] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rJA[1] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rVA[1] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
#else
	if( rIA[1] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rJA[1] != 0 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rVA[1] != 2 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
#endif
}

{
	// do not trigger double pass
	const struct rsb_mtx_partitioning_info_t * pinfop = NULL;
	const rsb_flags_t flags = RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT;
	const enum rsb_op_flags_t op_flags = RSB_OP_FLAG_DEFAULT;
	const int bs = RSB_DEFAULT_BLOCKING;
	const size_t wb = 0;
	const int br = bs, bc = bs;
	void * WA = NULL;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 7;
	const rsb_coo_idx_t nr = 6;
	const rsb_coo_idx_t nc = 6;
	rsb_coo_idx_t IA[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t JA[] = {0,1,2,3,4,5,0};
	RSB_DEFAULT_TYPE VA[] = {1,1,1,1,1,1,2};
	rsb_coo_idx_t rIA[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t rJA[] = {0,1,2,3,4,5,0};
	RSB_DEFAULT_TYPE rVA[] = {1,1,1,1,1,1,2};

	errval = rsb__do_index_based_bcsr_sort( IA, JA, VA, rIA, rJA, rVA, nr, nc, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	if( rIA[1] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rJA[1] != 0 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	if( rVA[1] != 2 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
}

	RSBENCH_STDOUT("SORT CHECK: END\n");
	//errval = RSB_ERR_NO_ERROR;
goto err;
merr:
	errval = RSB_ERR_INTERNAL_ERROR;
err:
	if(RSB_SOME_ERROR(errval))
		RSBENCH_STDOUT("SORT CHECK: FAIL\n");
	return errval;
}

static rsb_err_t rsb__chk_srt(const struct rsb_tester_options_t * top)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	int smi;
	const int sma[3] = {0,1,3}; // only 0 is default; 1,3 are internal and experimental
	rsb_coo_idx_t dim;
	const rsb_int_t asm_sort_method = rsb_global_session_handle.asm_sort_method;

// when 0, is defaults
// when 1, leads to rsb__do_index_based_bcsr_sort (two passes)
// when 3, leads to rsb__do_index_based_bcsr_sort irrespective
for(dim=30000;dim<100001;dim+=70000)
for(smi=0;smi<sizeof(sma)/sizeof(sma[0]);++smi)
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = RSB__ENOUGH_NNZ_FOR_PARALLEL_SORT;
	const rsb_coo_idx_t nrA = dim;
	const rsb_coo_idx_t ncA = dim;
	const rsb_coo_idx_t maxi = RSB_MIN(nnz,dim-1); // max i
	const rsb_coo_idx_t ioff = 0, joff = 0, nzoff = 0;
	rsb_nnz_idx_t n;
	rsb_coo_idx_t *IA = NULL;
	rsb_coo_idx_t *JA = NULL;
	RSB_DEFAULT_TYPE *VA = NULL;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

	if ( rsb__util_coo_alloc((void**)&VA, &IA, &JA, nnz, typecode, RSB_BOOL_TRUE) )
		goto err;

	rsb__do_fill_with_diag(VA, IA, JA, ioff, joff, nzoff, typecode, nnz);
	for ( n = 0 ; n < nnz -1  ; ++n )
		IA[n] = RSB_MIN(IA[n], maxi),
		JA[n] = RSB_MIN(JA[n], maxi);
	IA[nnz-1] = maxi;
	JA[nnz-1] = maxi;
	rsb__util_reverse_fullword_coo_array(IA,nnz);
	rsb__util_reverse_fullword_coo_array(JA,nnz);

	rsb_global_session_handle.asm_sort_method = sma[smi];

	errval = rsb_coo_sort(VA, IA, JA, nnz, nrA, ncA, typecode, flags );

	if( IA[nnz-1] != maxi )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(ierr,"Error: found a %ld where should be a %ld\n",(long)IA[nnz-1],(long)(maxi));
	}


	if( JA[nnz-1] != maxi )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(ierr,"Error: found a %ld where should be a %ld\n",(long)JA[nnz-1],(long)(maxi));
	}

ierr:
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(IA);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}
	goto ret;
err:
	errval = RSB_ERR_INTERNAL_ERROR;
ret:
	rsb_global_session_handle.asm_sort_method = asm_sort_method; // restore
	return errval;
}

static rsb_err_t rsb__chk_gen(const struct rsb_tester_options_t * top)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 3;
	rsb_coo_idx_t IA[] = {3,2,1};
	rsb_coo_idx_t JA[] = {2,1,0};
	RSB_DEFAULT_TYPE VA[] = {9,3,0};
	const RSB_DEFAULT_TYPE VAr[] = {1,1,1};
	const rsb_coo_idx_t IAr[] = {1,2,3};
	const rsb_coo_idx_t JAr[] = {2,3,4};
	rsb_nnz_idx_t l [] = {0,0,0,1,1};
	rsb_coo_idx_t ioff = 1, joff = 2;
	rsb_nnz_idx_t nzoff = 0;

	rsb__do_fill_with_diag(VA, IA, JA, ioff, joff, nzoff, typecode, nnz);
	//if(RSB_SOME_ERROR(errval))
	//	RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(VA,VAr,(nnz)*sizeof(VA[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_MEMCMP(IA,IAr,(nnz)*sizeof(IA[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_MEMCMP(JA,JAr,(nnz)*sizeof(JA[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}
err:
	return errval;
}

static rsb_err_t rsb__chk_sort(const struct rsb_tester_options_t * top)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
{
	const rsb_nnz_idx_t n = 3;
	const rsb_nnz_idx_t k [] = {3,2,1};
	const rsb_nnz_idx_t r [] = {3,0,1,2,0};
	rsb_nnz_idx_t l [] = {0,0,0,1,1};

	errval = rsb__do_msort_up(n, k, l);

	if(RSB_MEMCMP(l,r,(n+2)*sizeof(l[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_nnz_idx_t n = 3;
	const rsb_coo_idx_t k [] = {3,2,1,1,3,2};
	const rsb_nnz_idx_t r [] = {2,3,1,0,0};
	rsb_nnz_idx_t l [] = {0,0,0,1,1};

	errval = rsb__do_msort_up2coo(n, k, l);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(l,r,(n+2)*sizeof(l[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_nnz_idx_t n = 3;
	const rsb_coo_idx_t k [] = {1,1,3,2,2,3};
	const rsb_nnz_idx_t r [] = {1,3,0,2,0};
	rsb_nnz_idx_t l [] = {0,0,0,1,1};

	errval = rsb__do_msort_up2coo(n, k, l);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(l,r,(n+2)*sizeof(l[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_nnz_idx_t n = 3;
	const rsb_coo_idx_t k [] = {2,3,3,2,1,1};
	const rsb_nnz_idx_t r [] = {3,2,0,1,0};
	rsb_nnz_idx_t l [] = {0,0,0,1,1};

	errval = rsb__do_msort_up2coo(n, k, l);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(l,r,(n+2)*sizeof(l[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_nnz_idx_t n = 4;
	const rsb_coo_idx_t k [] = {0,0,3,3,1,1,2,2};
	const rsb_nnz_idx_t r [] = {1,3,0,4,2,0};
	rsb_nnz_idx_t l [] = {0,0,0,1,1,0};

	errval = rsb__do_msort_up2coo(n, k, l);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(l,r,(n+2)*sizeof(l[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_nnz_idx_t n = 4;
	const rsb_coo_idx_t k [] = {2,2,1,1,0,0,3,3};
	const rsb_nnz_idx_t r [] = {3,4,1,2,0,0};
	rsb_nnz_idx_t l [] = {0,0,0,1,1,0};

	errval = rsb__do_msort_up2coo(n, k, l);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(l,r,(n+2)*sizeof(l[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_nnz_idx_t n = 3;
	const rsb_coo_idx_t k [] = {1,1,2,3,3,2};
	rsb_nnz_idx_t l [] = {0,0,0,1,1};

	errval = rsb__do_msort_up2coo(n, k, l);

	if(errval != RSB_ERR_BADARGS)
		RSB_PERR_GOTO(err,RSB_ERRM_ES); // already sorted input expected to give error
}

{
	const rsb_nnz_idx_t n = 4;
	const rsb_coo_idx_t k [] = {11,11,5,5,4,4,0,0};
	const rsb_nnz_idx_t r [] = {4,0,1,2,3,0};
	rsb_nnz_idx_t l [] = {0,0,0,1,1,0};

	errval = rsb__do_msort_up2coo(n, k, l);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(RSB_MEMCMP(l,r,(n+2)*sizeof(l[0])))
		errval = RSB_ERR_INTERNAL_ERROR;

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

	errval = rsb__do_util_merge_sorted_subarrays_in_place_test();
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);


{
	// unit test for rsb__do_index_based_bcsr_sort
	const rsb_coo_idx_t ma[] = { 5, RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t)/2, RSB_MAX_VALUE_FOR_TYPE(rsb_half_idx_t) };
	rsb_coo_idx_t mi;

for(mi=0;mi<sizeof(ma)/sizeof(ma[0]);++mi)
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 3;
	RSB_DEFAULT_TYPE rVA[] = {24,13,32};
	rsb_coo_idx_t rIA[] = {2,1,3};
	rsb_coo_idx_t rJA[] = {4,3,2};
	rsb_coo_idx_t IA[] = {1,2,3};
	rsb_coo_idx_t JA[] = {3,4,2};
	RSB_DEFAULT_TYPE VA[] = {13,24,32};
	const rsb_coo_idx_t m = ma[mi], k = m;
	const rsb_coo_idx_t br = 1, bc = 1;
	const rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	const enum rsb_op_flags_t op_flags = RSB_OP_FLAG_DEFAULT;
	void * WA = NULL;
	const size_t wb = 0;
	rsb_nnz_idx_t nzi = 0;

	errval = rsb__do_index_based_bcsr_sort( IA, JA, VA, rIA, rJA, rVA,
		m, k, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	for(nzi=0;nzi<nnz;++nzi)
		if( ! (
			IA[nzi] == rIA[nzi] && 
			JA[nzi] == rJA[nzi] && 
			VA[nzi] == rVA[nzi] 
			) )
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
}

for(mi=0;mi<sizeof(ma)/sizeof(ma[0]);++mi)
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 6;
	RSB_DEFAULT_TYPE rVA[] = {24,21,11,13,33,32};
	rsb_coo_idx_t rIA[] = {2,2,1,1,3,3};
	rsb_coo_idx_t rJA[] = {4,1,1,3,3,2};
	rsb_coo_idx_t IA[] = {1,1,2,2,3,3};
	rsb_coo_idx_t JA[] = {1,3,1,4,2,3};
	RSB_DEFAULT_TYPE VA[] = {11,13,21,24,32,33};
	const rsb_coo_idx_t m = ma[mi], k = m;
	const rsb_coo_idx_t br = 1, bc = 1;
	const rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	const enum rsb_op_flags_t op_flags = RSB_OP_FLAG_DEFAULT;
	void * WA = NULL;
	const size_t wb = 0;
	rsb_nnz_idx_t nzi = 0;

	errval = rsb__do_index_based_bcsr_sort( IA, JA, VA, rIA, rJA, rVA,
		m, k, br, bc, nnz, typecode, flags, op_flags, WA, wb);

	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	for(nzi=0;nzi<nnz;++nzi)
		if( ! (
			IA[nzi] == rIA[nzi] && 
			JA[nzi] == rJA[nzi] && 
			VA[nzi] == rVA[nzi] 
			) )
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		}
}
}

err:
	return errval;
}

static rsb_err_t rsb__chk_asm(const struct rsb_tester_options_t * top)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */

	RSB_INFO("MATRIX ASSEMBLY FLAGS TEST: BEGIN\n");
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnzA = 7;
	const rsb_coo_idx_t nrA = 6;
	const rsb_coo_idx_t ncA = 6;
	const rsb_coo_idx_t IT[] = {0,1,2,3,4,5,1};
	const rsb_coo_idx_t JU[] = {0,1,2,3,4,5,5};
	const rsb_coo_idx_t JL[] = {0,1,2,3,4,5,0};
	const RSB_DEFAULT_TYPE VT[] = {1,1,1,1,1,1,2};
	char ib[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];

	if(top->wvr==RSB_BOOL_TRUE)
		RSB_INFO("Check rsb_mtx_alloc_from_coo_const for triangular\n");
	mtxAp = rsb_mtx_alloc_from_coo_const(
		VT,IT,JL,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
		| RSB_FLAG_DUPLICATES_SUM
		| RSB_FLAG_TRIANGULAR
		, &errval);
	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(top->wvr==RSB_BOOL_TRUE)
		rsb__print_matrix_stats(mtxAp);
	if( mtxAp->flags != (
		RSB_FLAG_USE_HALFWORD_INDICES |
		RSB_FLAG_SORTED_INPUT |
		RSB_FLAG_TRIANGULAR |
		RSB_FLAG_LOWER |
		RSB_FLAG_WANT_COO_STORAGE |
		RSB_FLAG_DUPLICATES_SUM |
		RSB_FLAG_QUAD_PARTITIONING |
		RSB_FLAG_WANT_BCSS_STORAGE |
		RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS |
		RSB_FLAG_OWN_PARTITIONING_ARRAYS |
		RSB_FLAG_SORT_INPUT
	))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_MTX_FREE(mtxAp);

	if(top->wvr==RSB_BOOL_TRUE)
		RSB_INFO("Check rsb_mtx_alloc_from_coo_const implicit diagonal\n");
	mtxAp = rsb_mtx_alloc_from_coo_const(
		VT,IT,JL,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
		| RSB_FLAG_DUPLICATES_SUM
		| RSB_FLAG_UNIT_DIAG_IMPLICIT
		| RSB_FLAG_TRIANGULAR
		, &errval);
	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(top->wvr==RSB_BOOL_TRUE)
		rsb__print_matrix_stats(mtxAp);
	if( mtxAp->flags != (
		RSB_FLAG_USE_HALFWORD_INDICES |
		RSB_FLAG_SORTED_INPUT |
		RSB_FLAG_TRIANGULAR |
		RSB_FLAG_LOWER |
		RSB_FLAG_UNIT_DIAG_IMPLICIT |
		RSB_FLAG_WANT_COO_STORAGE |
		RSB_FLAG_DUPLICATES_SUM |
		RSB_FLAG_WANT_BCSS_STORAGE |
		RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS |
		RSB_FLAG_OWN_PARTITIONING_ARRAYS |
		RSB_FLAG_SORT_INPUT
	))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_MTX_FREE(mtxAp);

	if(top->wvr==RSB_BOOL_TRUE)
		RSB_INFO("Check rsb_mtx_alloc_from_coo_const implicit diagonal\n");
	mtxAp = rsb_mtx_alloc_from_coo_const(
		VT,IT,JU,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
		| RSB_FLAG_DUPLICATES_SUM
		| RSB_FLAG_UNIT_DIAG_IMPLICIT
		| RSB_FLAG_TRIANGULAR
		, &errval);
	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(top->wvr==RSB_BOOL_TRUE)
		rsb__print_matrix_stats(mtxAp);
	if( mtxAp->flags != (
		RSB_FLAG_USE_HALFWORD_INDICES |
		RSB_FLAG_SORTED_INPUT |
		RSB_FLAG_TRIANGULAR |
		RSB_FLAG_UPPER |
		RSB_FLAG_UNIT_DIAG_IMPLICIT |
		RSB_FLAG_WANT_COO_STORAGE |
		RSB_FLAG_DUPLICATES_SUM |
		RSB_FLAG_WANT_BCSS_STORAGE |
		RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS |
		RSB_FLAG_OWN_PARTITIONING_ARRAYS |
		RSB_FLAG_SORT_INPUT
	))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_MTX_FREE(mtxAp);

	if(top->wvr==RSB_BOOL_TRUE)
		RSB_INFO("Check rsb_mtx_alloc_from_coo_begin for triangular\n");
	mtxAp = rsb_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
		| RSB_FLAG_DUPLICATES_SUM
		| RSB_FLAG_TRIANGULAR
		, &errval);
	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = rsb_mtx_set_vals(mtxAp, VT, IT, JL, nnzA, RSB_FLAG_NOFLAGS);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = rsb_mtx_alloc_from_coo_end(&mtxAp);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(top->wvr==RSB_BOOL_TRUE)
		rsb__print_matrix_stats(mtxAp);
	if(mtxAp->flags != (
		RSB_FLAG_USE_HALFWORD_INDICES |
		RSB_FLAG_SORTED_INPUT |
		RSB_FLAG_TRIANGULAR |
		RSB_FLAG_LOWER |
		RSB_FLAG_WANT_COO_STORAGE |
		RSB_FLAG_QUAD_PARTITIONING |
		RSB_FLAG_WANT_BCSS_STORAGE |
		RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS |
		RSB_FLAG_OWN_PARTITIONING_ARRAYS |
		RSB_FLAG_SORT_INPUT
	))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}	RSB_MTX_FREE(mtxAp);

	RSB_MTX_FREE(mtxAp);

	if(top->wvr==RSB_BOOL_TRUE)
		RSB_INFO("Check rsb_mtx_alloc_from_coo_begin implicit diagonal\n");
	mtxAp = rsb_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS /* force rsb */
		| RSB_FLAG_DUPLICATES_SUM/* sum dups */
		| RSB_FLAG_UNIT_DIAG_IMPLICIT/* ask diagonal implicit */
		| RSB_FLAG_TRIANGULAR /* */
		, &errval);
	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = rsb_mtx_set_vals(mtxAp, VT, IT, JL, nnzA, RSB_FLAG_NOFLAGS);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = rsb_mtx_alloc_from_coo_end(&mtxAp);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(top->wvr==RSB_BOOL_TRUE)
		rsb__print_matrix_stats(mtxAp);
	if(mtxAp->flags != (
		RSB_FLAG_USE_HALFWORD_INDICES |
		RSB_FLAG_SORTED_INPUT |
		RSB_FLAG_TRIANGULAR |
		RSB_FLAG_LOWER |
		RSB_FLAG_UNIT_DIAG_IMPLICIT |
		RSB_FLAG_WANT_COO_STORAGE |
		RSB_FLAG_WANT_BCSS_STORAGE |
		RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS |
		RSB_FLAG_OWN_PARTITIONING_ARRAYS |
		RSB_FLAG_SORT_INPUT
	))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}	RSB_MTX_FREE(mtxAp);
	RSB_MTX_FREE(mtxAp);
}

{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnzA = 7;
	const rsb_coo_idx_t nrA = 6;
	const rsb_coo_idx_t ncA = 6;
	rsb_coo_idx_t IT[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t JU[] = {0,1,2,3,4,5,5};
	rsb_coo_idx_t JL[] = {0,1,2,3,4,5,0};
	RSB_DEFAULT_TYPE VT[] = {1,1,1,1,1,1,2};
	char ib[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];

	if(top->wvr==RSB_BOOL_TRUE)
		RSB_INFO("Check rsb_mtx_alloc_from_coo_inplace for triangular\n");
	mtxAp = rsb_mtx_alloc_from_coo_inplace(
		VT,IT,JL,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
		| RSB_FLAG_DUPLICATES_SUM
		| RSB_FLAG_TRIANGULAR
		, &errval);
	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(top->wvr==RSB_BOOL_TRUE)
		rsb__print_matrix_stats(mtxAp);
	if(mtxAp->flags != (
		RSB_FLAG_USE_HALFWORD_INDICES |
		RSB_FLAG_SORTED_INPUT |
		RSB_FLAG_TRIANGULAR |
		RSB_FLAG_LOWER |
		RSB_FLAG_WANT_COO_STORAGE |
		RSB_FLAG_DUPLICATES_SUM |
		RSB_FLAG_QUAD_PARTITIONING |
		RSB_FLAG_WANT_BCSS_STORAGE |
		RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS |
		RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS |
		RSB_FLAG_OWN_PARTITIONING_ARRAYS |
		RSB_FLAG_SORT_INPUT
	))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_MTX_FREE(mtxAp);
}
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnzA = 7;
	const rsb_coo_idx_t nrA = 6;
	const rsb_coo_idx_t ncA = 6;
	rsb_coo_idx_t IT[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t JU[] = {0,1,2,3,4,5,5};
	rsb_coo_idx_t JL[] = {0,1,2,3,4,5,0};
	RSB_DEFAULT_TYPE VT[] = {1,1,1,1,1,1,2};
	char ib[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];

	if(top->wvr==RSB_BOOL_TRUE)
		RSB_INFO("Check rsb_mtx_alloc_from_coo_inplace for implicit diagonal\n");
	mtxAp = rsb_mtx_alloc_from_coo_inplace(
		VT,IT,JL,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
		| RSB_FLAG_DUPLICATES_SUM
		| RSB_FLAG_UNIT_DIAG_IMPLICIT
		| RSB_FLAG_TRIANGULAR
		, &errval);
	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(top->wvr==RSB_BOOL_TRUE)
		rsb__print_matrix_stats(mtxAp);
	if(mtxAp->flags != (
		RSB_FLAG_USE_HALFWORD_INDICES |
		RSB_FLAG_SORTED_INPUT |
		RSB_FLAG_TRIANGULAR |
		RSB_FLAG_LOWER |
		RSB_FLAG_UNIT_DIAG_IMPLICIT |
		RSB_FLAG_WANT_COO_STORAGE |
		RSB_FLAG_DUPLICATES_SUM |
		RSB_FLAG_WANT_BCSS_STORAGE |
		RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS |
		RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS |
		RSB_FLAG_OWN_PARTITIONING_ARRAYS |
		RSB_FLAG_SORT_INPUT
	))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
	RSB_MTX_FREE(mtxAp);
}
	RSB_INFO("MATRIX ASSEMBLY FLAGS TEST: END\n");
	errval = RSB_ERR_NO_ERROR;
err:
	if(RSB_SOME_ERROR(errval))
		RSB_INFO("MATRIX ASSEMBLY FLAGS TEST: FAIL\n");
	RSB_MTX_FREE(mtxAp);
	return errval;
}

static rsb_err_t rsb__chk_trtr(const struct rsb_tester_options_t * top)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */

	RSB_INFO("REGRESSION TEST: BEGIN\n");

#define RSB_EXPECT_TRTR_BUG_COO 0

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
	const rsb_nnz_idx_t nnzA = 3;
	const rsb_coo_idx_t incX = 1, incY = 1;
	const rsb_coo_idx_t nrA = 3;
	const rsb_coo_idx_t ncA = 1;
	const rsb_coo_idx_t IA[] = {0,1,2};
	const rsb_coo_idx_t JA[] = {0,0,0};
	const double VA[] = {1,1,1};
	const double X[] = {1,1,1};
	double Y[] = {1,1,1};
	const double alpha = 1.0;
	const double beta = 0.0;
	const rsb_trans_t transA = RSB_TRANSPOSITION_C;

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES
		, &errval);

	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( Y[0] != 1 || Y[1] != 1 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = rsb_spmv(transA, &alpha, mtxAp, X, incX, &beta, Y, incY);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#if RSB_EXPECT_TRTR_BUG_COO
	if( Y[0] != 3 || Y[1] != 0 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	// bug: Y[0:2] have been zeroed
#else
	if( Y[0] != 3 || Y[2] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	// ok : only Y[0] should have been zeroed
#endif

	RSB_MTX_FREE(mtxAp);
}
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
	const rsb_nnz_idx_t nnzA = 3;
	const rsb_coo_idx_t incX = 1, incY = 1;
	const rsb_coo_idx_t nrA = 1;
	const rsb_coo_idx_t ncA = 3;
	const rsb_coo_idx_t IA[] = {0,0,0};
	const rsb_coo_idx_t JA[] = {0,1,2};
	const double VA[] = {1,1,1};
	const double X[] = {1,1,1};
	double Y[] = {1,1,1};
	const double alpha = 1.0;
	const double beta = 0.0;
	const rsb_trans_t transA = RSB_TRANSPOSITION_C;

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES
		, &errval);

	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( Y[0] != 1 || Y[1] != 1 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = rsb_spmv(transA, &alpha, mtxAp, X, incX, &beta, Y, incY);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#if RSB_EXPECT_TRTR_BUG_COO
	if( Y[0] != 1 || Y[1] != 2 ) // bug: Y[1] should have been zeroed before add
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
#else
	if( Y[0] != 1 || Y[1] != 1 ) // ok : Y[1] should have been zeroed before add
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
#endif

	RSB_MTX_FREE(mtxAp);
}
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
	const rsb_nnz_idx_t nnzA = 4;
	const rsb_coo_idx_t incX = 1, incY = 1;
	const rsb_coo_idx_t nrA = 3;
	const rsb_coo_idx_t ncA = 2;
	const rsb_coo_idx_t IA[] = {0,1,2,2};
	const rsb_coo_idx_t JA[] = {0,0,0,1};
	const double VA[] = {1,1,1,1};
	const double X[] = {1,1,1};
	double Y[] = {1,1,1};
	const double alpha = 1.0;
	const double beta = 0.0;
	const rsb_trans_t transA = RSB_TRANSPOSITION_C;

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES
		, &errval);

	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( Y[0] != 1 || Y[1] != 1 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	errval = rsb_spmv(transA, &alpha, mtxAp, X, incX, &beta, Y, incY);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#if RSB_EXPECT_TRTR_BUG_COO
	if( Y[0] != 3 || Y[2] != 0 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	// bug: Y[0:2] have been zeroed, instead of Y[0:1]
#else
	if( Y[0] != 3 || Y[2] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	// ok : only Y[0:1] should have been zeroed
#endif

	RSB_MTX_FREE(mtxAp);
}
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
#define RSB_EXPECT_RECT_IMPDIA_NRHS_BUG 0
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
	const rsb_nnz_idx_t nnzA = 2;
	const rsb_coo_idx_t incX = 1, incY = 1;
	const rsb_coo_idx_t nrA = 2;
	const rsb_coo_idx_t ncA = 1;
	const rsb_coo_idx_t IA[] = {0,1};
	const rsb_coo_idx_t JA[] = {0,0};
	const rsb_coo_idx_t nrhs = 2;
	const double VA[] = {1,1};
	const double X[] = {1,1};
	double Y[] = {0,0,0,0};
	const double alpha = 1.0;
	const double beta = 1.0;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_UNIT_DIAG_IMPLICIT
		, &errval);

	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	errval = rsb_spmm(transA, &alpha, mtxAp, nrhs, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, X, 1, &beta, Y, 2);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#if RSB_EXPECT_RECT_IMPDIA_NRHS_BUG
	if( Y[1] != 2 || Y[3] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	// bug: Y[1] has been added an extra 1
#else
	if( Y[1] != 1 || Y[3] != 1 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	// ok : Y[1] has been added no extra 1
#endif

	RSB_MTX_FREE(mtxAp);
}
#undef RSB_EXPECT_RECT_IMPDIA_NRHS_BUG
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
#define RSB_EXPECT_TUNE_BROKEN_FOR_ROWS_MAJOR 0
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
	const rsb_nnz_idx_t nnzA = 2;
	const rsb_coo_idx_t incX = 1, incY = 1;
	const rsb_coo_idx_t nrA = 3; // overzealous rsb_tune_spmm bug triggers on e.g. ldB<nrA (like when RSB_FLAG_WANT_ROW_MAJOR_ORDER)
	const rsb_coo_idx_t ncA = 1;
	const rsb_coo_idx_t IA[] = {0,1};
	const rsb_coo_idx_t JA[] = {0,0};
	const rsb_coo_idx_t nrhs = 2;
	const rsb_coo_idx_t ldB = nrhs;
	const rsb_coo_idx_t ldC = nrhs;
	const double VA[] = {1,1};
	const double B[] = {1,1};
	double C[] = {0,0,0,0};
	const double alpha = 1.0;
	const double beta = 1.0;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_UNIT_DIAG_IMPLICIT
		, &errval);

	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	errval = rsb_tune_spmm(NULL, NULL, NULL, 0.0, 0.0, transA, &alpha, mtxAp, nrhs, RSB_FLAG_WANT_ROW_MAJOR_ORDER, B, ldB, &beta, C, ldC);

#if RSB_EXPECT_TUNE_BROKEN_FOR_ROWS_MAJOR
	if(!RSB_SOME_ERROR(errval))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
#else
	if( RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
#endif

	RSB_MTX_FREE(mtxAp);
}
#undef RSB_EXPECT_TUNE_BROKEN_FOR_ROWS_MAJOR
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

#ifdef RSB_NUMERICAL_TYPE_DOUBLE
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
	const rsb_nnz_idx_t nnzA = 1;
	const rsb_coo_idx_t nrA = 1;
	const rsb_coo_idx_t ncA = 1;
	const rsb_coo_idx_t IA[] = {0};
	const rsb_coo_idx_t JA[] = {0};
	const rsb_coo_idx_t nrhs = 2;
	const rsb_coo_idx_t ldB = nrhs;
	const rsb_coo_idx_t ldC = nrhs + 1;
	const double VA[] = {1};
	const double B[] = {1,2,3,4};
	double C[] = {0,0,0,0};
	const double alpha = 1.0;
	const double beta = 1.0;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	int i;

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_UNIT_DIAG_IMPLICIT,
		&errval);

	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	errval = rsb_spmm(transA, &alpha, mtxAp, nrhs, RSB_FLAG_WANT_ROW_MAJOR_ORDER, B, ldB, &beta, C, ldC);

	if( C[0] != 1 || C[1] != 2 || C[2] != 0 || C[3] != 0 )
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);

	RSB_MTX_FREE(mtxAp);
}
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */

	RSB_INFO("REGRESSION TEST: END\n");
	errval = RSB_ERR_NO_ERROR;
goto err;
merr:
	errval = RSB_ERR_INTERNAL_ERROR;
err:
	if(RSB_SOME_ERROR(errval))
		RSB_INFO("REGRESSION TEST: FAIL\n");
	RSB_MTX_FREE(mtxAp);
	return errval;
}

static rsb_err_t rsb__chk_coo(void)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 3+1;
	rsb_coo_idx_t IA[] = {3,2,2,1};
	rsb_coo_idx_t JA[] = {2,1,1,0};
	RSB_DEFAULT_TYPE VA[] = {9,3,2,0};

	if(3 != rsb__weed_out_duplicates(IA, JA, VA, nnz, typecode, RSB_FLAG_NOFLAGS))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[2]!=1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[2]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[2]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 4;
	const rsb_nnz_idx_t fd = 0;
	rsb_nnz_idx_t moved = 0, moves = 0;
	rsb_coo_idx_t IA[] = {3,2,1,0                   };
	rsb_coo_idx_t JA[] = {RSB_MARKER_COO_VALUE,2,1,2};
	RSB_DEFAULT_TYPE VA[] = {9,3,2,0};

	errval = rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode),fd,&moved,&moves);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[0]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[0]!=2)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[0]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(moved != 1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(moves != 1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 4;
	const rsb_nnz_idx_t fd = 2;
	rsb_nnz_idx_t moved = 0, moves = 0;
	rsb_coo_idx_t IA[] = {3,2,                   1,1};
	rsb_coo_idx_t JA[] = {2,1,RSB_MARKER_COO_VALUE,0};
	RSB_DEFAULT_TYPE VA[] = {9,3,2,0};

	errval = rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode),fd,&moved,&moves);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[2]!=1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[2]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[2]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 6;
	const rsb_nnz_idx_t fd = 1;
	rsb_nnz_idx_t moved = 0, moves = 0;
	rsb_coo_idx_t IA[] = {3,1,2,2,                   1,1};
	rsb_coo_idx_t JA[] = {2,4,2,2,RSB_MARKER_COO_VALUE,0};
	RSB_DEFAULT_TYPE VA[] = {9,3,2,2,2,0};

	errval = rsb__util_compact_marked_coo_array(IA,JA,VA,nnz,RSB_NUMERICAL_TYPE_SIZE(typecode),fd,&moved,&moves);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[3]!=1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[3]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[3]!=0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}


{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 7;
	rsb_coo_idx_t IA[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t JA[] = {0,1,2,3,4,5,9};
	RSB_DEFAULT_TYPE VA[] = {0,0,0,0,0,0,8};
	rsb_nnz_idx_t gap = -1;
	rsb_nnz_idx_t discarded = 0;

	errval = rsb__weed_out_diagonal(VA, IA, JA, nnz, typecode, &gap, &discarded);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[0]!=1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[0]!=9)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[0]!=8)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(discarded!=6)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(gap!=-1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 7;
	rsb_coo_idx_t IA[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t JA[] = {0,1,2,3,4,5,9};
	RSB_DEFAULT_TYPE VA[] = {1,1,1,1,1,1,9};
	rsb_nnz_idx_t gap = -1;
	rsb_nnz_idx_t discarded = 0;

	errval = rsb__weed_out_diagonal(VA, IA, JA, nnz, typecode, &gap, &discarded);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[0]!=1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[0]!=9)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[0]!=9)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(discarded!=6)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(gap!=-1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 7;
	rsb_coo_idx_t IA[] = {0,1,2,3,4,5,1};
	rsb_coo_idx_t JA[] = {0,1,2,3,4,5,9};
	RSB_DEFAULT_TYPE VA[] = {1,0,0,0,0,0,9};
	rsb_nnz_idx_t gap = -1;
	rsb_nnz_idx_t discarded = 0;

	errval = rsb__weed_out_diagonal(VA, IA, JA, nnz, typecode, &gap, &discarded);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[0]!=1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[0]!=9)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[0]!=9)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(discarded!=6)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(gap!=-1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnz = 6;
	rsb_coo_idx_t IA[] = {1,2,3,2,4,6};
	rsb_coo_idx_t JA[] = {1,1,2,3,5,6};
	RSB_DEFAULT_TYPE VA[] = {1,2,3,4,5,6};
	rsb_nnz_idx_t gap = -1;
	rsb_nnz_idx_t discarded = 0;

	errval = rsb__weed_out_non_lowtri(VA, IA, JA, nnz, typecode, &gap, &discarded);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(IA[3]!=6)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(JA[3]!=6)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(VA[3]!=6)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(discarded!=2)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	if(gap!=-1)
		RSB_PERR_GOTO(err,RSB_ERRM_ES); // FIXME: need handling
}
	errval = RSB_ERR_NO_ERROR;
	goto ret;
err:
	errval = RSB_ERR_INTERNAL_ERROR;
ret:
	return errval;
}

static rsb_err_t rsb__chk_csr2coo(void)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;

{
	const rsb_coo_idx_t x = -1;
	const rsb_coo_idx_t nr = 3;
	const rsb_coo_idx_t nz = 5;
	const rsb_coo_idx_t off = 0;
	rsb_coo_idx_t * TA = NULL;
	const rsb_coo_idx_t RP[] = {0,0,2,2,2};
	rsb_coo_idx_t IA[] = {0,2,2,nz,x};
	rsb_coo_idx_t i;

	errval = rsb__do_switch_compressed_array_to_fullword_coo(IA, nr, off, TA);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(ret,RSB_ERRM_ES);

	for(i=0;i<nz;++i)
		if(RP[i]!=IA[i])
			RSB_PERR_GOTO(err,RSB_ERRM_ES);

}

{
	const rsb_coo_idx_t x = -1;
	const rsb_coo_idx_t nr = 3;
	const rsb_coo_idx_t nz = 5;
	const rsb_coo_idx_t off = 0;
	rsb_coo_idx_t TA[nr+1];
	const rsb_coo_idx_t RP[] = {0,0,2,2,2};
	rsb_coo_idx_t IA[] = {0,2,2,nz,x};
	rsb_coo_idx_t i;

	errval = rsb__do_switch_compressed_array_to_fullword_coo(IA, nr, off, TA);
	if(RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(ret,RSB_ERRM_ES);

	for(i=0;i<nz;++i)
		if(RP[i]!=IA[i])
			RSB_PERR_GOTO(err,RSB_ERRM_ES);

}

{
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_coo_idx_t nr = 2000;
	const rsb_coo_idx_t mj = nr;
	const rsb_coo_idx_t nc = 2 * mj;
	rsb_coo_idx_t *IA = NULL;
	rsb_coo_idx_t *JA = NULL;
	rsb_coo_idx_t *IB = NULL;
	rsb_coo_idx_t *JB = NULL;
	RSB_DEFAULT_TYPE *VA = NULL;
	RSB_DEFAULT_TYPE *VB = NULL;
	const rsb_nnz_idx_t nnzA = nr * mj + (nr*(1+mj))/2;
	const rsb_nnz_idx_t nzoffA = 0, nzoffB = 0;
	const rsb_coo_idx_t iadd = 0, jadd = - mj;
 	rsb_nnz_idx_t nzl, nzr;
 	rsb_nnz_idx_t i,j,nzi=0;

	errval = RSB_ERR_INTERNAL_ERROR;

	if ( rsb__util_coo_alloc((void**)&VA, &IA, &JA, nnzA, typecode, RSB_BOOL_TRUE) )
		RSB_PERR_GOTO(ierr,RSB_ERRM_ES);

	if ( rsb__util_coo_alloc((void**)&VB, &IB, &JB, nnzA, typecode, RSB_BOOL_TRUE) )
		RSB_PERR_GOTO(ierr,RSB_ERRM_ES);

	for(i=0;i<nr;++i)
		for(j=0;j<mj+i+1;++j)
		{
			IA[nzi] = i;
			JA[nzi] = j;
			++nzi;
		}

	if ( nzi != nnzA )
		RSB_PERR_GOTO(ierr,"%ld != %ld\n",(long)nzi,(long)nnzA);

	rsb__coo_to_lr(VB, IB, JB, VA, IA, JA, mj, nnzA, nzoffB, nzoffA, &nzl, &nzr, iadd, jadd, typecode);

	nzi = 0;
	for(i=0;i<nr;++i)
		for(j=0;j<mj;++j)
		{
			if ( IA[nzi] != i ) 
				RSB_PERR_GOTO(ierr,"@%ld: %ld != %ld\n",(long)nzi,(long)i,(long)IA[nzi]);
			if ( JA[nzi] != j ) 
				RSB_PERR_GOTO(ierr,"@%ld: %ld != %ld\n",(long)nzi,(long)j,(long)JA[nzi]);
			++nzi;
		}

	for(i=0;i<nr;++i)
		for(j=0;j<i+1;++j)
		{
			if ( IA[nzi] != i ) 
				RSB_PERR_GOTO(ierr,"@%ld: %ld != %ld\n",(long)nzi,(long)i,(long)IA[nzi]);
			if ( JA[nzi] != j ) 
				RSB_PERR_GOTO(ierr,"@%ld: %ld != %ld\n",(long)nzi,(long)j,(long)JA[nzi]);
			++nzi;
		}

	if ( nzi != nnzA )
		RSB_PERR_GOTO(ierr,"%ld != %ld\n",(long)nzi,(long)nnzA);

	errval = RSB_ERR_NO_ERROR;
ierr:
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(IB);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(JB);
	RSB_CONDITIONAL_FREE(VA);
	RSB_CONDITIONAL_FREE(VB);
	if( errval == RSB_ERR_NO_ERROR )
		goto ret;
	else
		goto err;
}

err:
	errval = RSB_ERR_INTERNAL_ERROR;
ret:
	return errval;
}

static rsb_err_t rsb__chk_sppsp(void)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;	/* matrix structure pointer */
	struct rsb_mtx_t *mtxCp = NULL;	/* matrix structure pointer */

	RSB_INFO("MATRIX SUMS TEST: BEGIN\n");
{
	const int bs = RSB_DEFAULT_BLOCKING;
	const int brA = bs, bcA = bs;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	const rsb_nnz_idx_t nnzA = 6;
	const rsb_coo_idx_t nrA = 4;
	const rsb_coo_idx_t ncA = 4;
	const rsb_coo_idx_t IT[] = {1,1,1,2,2,4};
	const rsb_coo_idx_t JL[] = {1,3,4,1,2,4};
	const RSB_DEFAULT_TYPE VT[] = {11,13,14,-21,22,44};

	const void *alphap = NULL;
	const void *betap = NULL;
	struct rsb_mtx_t *mtxBp = NULL;
	const rsb_trans_t transA = RSB_TRANSPOSITION_T;
	const rsb_trans_t transB = RSB_TRANSPOSITION_T;

	mtxAp = rsb_mtx_alloc_from_coo_const(
		VT,IT,JL,nnzA,typecode,nrA,ncA,brA,bcA,
		RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
		| RSB_FLAG_DUPLICATES_SUM
		| RSB_FLAG_FORTRAN_INDICES_INTERFACE
		, &errval);

	if((!mtxAp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	mtxBp = mtxAp;
	mtxCp = rsb_sppsp(typecode, transA, alphap, mtxAp, transB, betap, mtxBp, &errval);

	if((!mtxCp) || RSB_SOME_ERROR(errval))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(mtxAp->nnz != mtxCp->nnz)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}
	RSB_INFO("MATRIX SUMS TEST: END\n");
	goto ret;
err:
	RSB_INFO("MATRIX SUMS TEST: FAIL\n");
	errval = RSB_ERR_INTERNAL_ERROR;
ret:
	RSB_MTX_FREE(mtxAp);
	RSB_MTX_FREE(mtxCp);
	return errval;
}

static rsb_err_t rsb__chk_swt(void)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;

	if( !RSB_INDICES_FIT_IN_HALFWORD(0, 0))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( !RSB_INDICES_FIT_IN_HALFWORD(40000, 40000))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(  RSB_INDICES_FIT_IN_HALFWORD(70000, 70000))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(  RSB_INDICES_FIT_IN_HALFWORD(70000, 0))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(  RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(1) )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(  RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(32000) )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(  RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(40000) )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( !RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(70000) )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if(!!RSB_IS_COO_VALUE_MORE_THAN_HALF_BITS_LONG(-1))
		RSB_PERR_GOTO(err,RSB_ERRM_ES); // only meant for positives

	if( rsb__do_is_candidate_size_for_halfword_coo(1,1,RSB_FLAG_USE_HALFWORD_INDICES_COO) != RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_coo(1,1,RSB_FLAG_USE_HALFWORD_INDICES) == RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_coo(70000,70000,RSB_FLAG_USE_HALFWORD_INDICES_COO) == RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_coo(70000,70000,RSB_FLAG_USE_HALFWORD_INDICES) == RSB_BOOL_TRUE )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(1,1,1,RSB_FLAG_USE_HALFWORD_INDICES_COO) != RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES); // note: this may seem counter-intuitive

	if( rsb__do_is_candidate_size_for_halfword_csr(1,1,1,RSB_FLAG_USE_HALFWORD_INDICES_CSR) != RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(70000,70000,1,RSB_FLAG_USE_HALFWORD_INDICES_COO) == RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(70000,7000 ,1,RSB_FLAG_USE_HALFWORD_INDICES_COO) != RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(70000,70000,1,RSB_FLAG_USE_HALFWORD_INDICES_CSR) == RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(70000,7000 ,1,RSB_FLAG_USE_HALFWORD_INDICES_CSR) != RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(70000,70000,1,RSB_FLAG_USE_HALFWORD_INDICES) == RSB_BOOL_TRUE )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(70000,70000,1,RSB_FLAG_USE_HALFWORD_INDICES_CSR) == RSB_BOOL_TRUE )
	         RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__do_is_candidate_size_for_halfword_csr(70000,70000,1,RSB_FLAG_USE_HALFWORD_INDICES) == RSB_BOOL_TRUE )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__is_recursive_matrix(RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS) != RSB_BOOL_TRUE )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__is_recursive_matrix(RSB_FLAG_DEFAULT_STORAGE_FLAGS) == RSB_BOOL_TRUE )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__is_recursive_matrix(RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS) == RSB_BOOL_TRUE )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__is_recursive_matrix(RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS) == RSB_BOOL_TRUE )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	errval = RSB_ERR_NO_ERROR;
err:
	return errval;
}

static rsb_err_t rsb__chk_idx(void)
{
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;

{
	const rsb_int_t A[] = {99,2,1,-20};
	const rsb_int_t n = 4;

	if( rsb__util_find_min_index_val(A,n)!= -20 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__util_find_max_index_val(A,n)!=  99 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__util_find_min_index(A,n)!= 3 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__util_find_max_index(A,n)!= 0 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	const rsb_coo_idx_t A[] = {RSB_MARKER_COO_VALUE,1,1,-20,30,20};
	const rsb_nnz_idx_t n = 4;

	if( rsb__util_find_coo_val_idx(A,n, RSB_MARKER_COO_VALUE) != 0 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__util_find_coo_val_idx(A,n, 20) != n )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__util_find_coo_val_idx(A,n,-20) != n-1 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( rsb__util_find_coo_val_idx(A,n,  1) != 1 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

{
	// unit test rsb__util_hcoo_array_add
	rsb_half_idx_t p[] = {99,2,1,100};
	const rsb_nnz_idx_t n = 4;

	rsb__util_hcoo_array_add(p, n-1, 99);

	if( p[2] != 100 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	if( p[3] != 100 )
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
}

	errval = RSB_ERR_NO_ERROR;
err:
	return errval;
}

int rsb_genmm_main(int argc,char *argv[]);
int rsb_mtx_ls_main(const int argc, char * argv[], rsb_bool_t want_latex);

int main(const int argc, char * argv[])
{
	rsb_option_t options[] = {
	    {"help",			no_argument, NULL, 'h' },
	    {"bench",			no_argument, NULL, 0x626e6368 }, /* bnch */
	    {"matrix-operation",	required_argument, NULL, 'o' },
	    {"subprogram-operation",	required_argument, NULL, 'O' },
	    {"information",		no_argument, NULL, 'I' },
	    {"configuration",		no_argument, NULL, 'C' },
	    {"hardware-counters",	no_argument, NULL, 'H' },
	    {"memory-benchmark",	no_argument, NULL, 'M' },
	    {"experiments",		no_argument, NULL, 'e' },
	    {"version",			no_argument, NULL, 'v' },
	    {"blas-testing",		no_argument, NULL, 'B' },
	    {"quick-blas-testing",		required_argument, NULL, 'Q' },
	    {"error-testing",		required_argument, NULL, 'E' },
	    {"fp-bench",		no_argument, NULL, 'F' },
	    {"transpose-test",		no_argument, NULL, 't' },
	    {"limits-testing",		no_argument, NULL, 0x6c696d74 },
	    {"guess-blocking",		no_argument, NULL, 'G' },	/* will pass guess parameters here some day (FIXME: obsolete) */
	    {"generate-matrix",		no_argument, NULL, 'g' }, /* should be synced to rsb_genmm_main */
	    {"plot-matrix",		no_argument, NULL,  0x50505050},/* should be synced to rsb__dump_postscript */
	    {"matrix-ls",		no_argument, NULL,  0x006D6C73},/* synced to rsb_mtx_ls_main */
	    {"matrix-ls-latex",		no_argument, NULL,  0x6c736c61},/* synced to rsb_mtx_ls_main */
	    {"matrix-print",		required_argument, NULL, 'P' },
	    {"read-performance-record",		required_argument, NULL,  0x7270720a},/*  */
	    {"help-read-performance-record",		no_argument, NULL,  0x72707268},/*  */
	    {"setenv",	required_argument, NULL, 0x7365},/* se */  
   {0,no_argument,0,0}
	};

	/*
	 * NOTE: this implies that unless an argument reset mechanism is implemented here,
	 * the o and O options will be forwarded to the host program!
	 * */
	//const char default_operation='v';
	const char default_operation='a';
	char operation=default_operation;
	//const char default_program_operation='r';
	const char default_program_operation='b';
	char program_operation=default_program_operation;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
#if RSB_WANT_REDUCED_RSB_M4_MATRIX_META_OPS
	const char * program_codes = "a"
#else /* RSB_WANT_REDUCED_RSB_M4_MATRIX_META_OPS */
	const char * program_codes = "avms"
#endif /* RSB_WANT_REDUCED_RSB_M4_MATRIX_META_OPS */
#if !RSB_SHALL_UPDATE_COMPLETEBENCHS
		"c"
#endif /* RSB_SHALL_UPDATE_COMPLETEBENCHS */
#ifdef RSB_OPTYPE_INDEX_SPSV_UXUA
		"t"
#endif
		"inS";
	int program_not_chosen=1;
	struct rsb_tester_options_t to;
	rsb_char_t cs[RSB_MAX_VERSION_STRING_LENGTH];
	int c;
	int opt_index = 0, aoc = 0; /* acceptable options count */
	// extern int optind; //unistd.h

	rsb_blas_tester_options_init(&to);

	for (;program_not_chosen;)
	{
		c = rsb__getopt_long(argc,argv,"CP:"
				"gBGvMHIho:O:"
/* #if RSB_WANT_EXPERIMENTS_CODE
				"e"
#endif */ /* RSB_WANT_EXPERIMENTS_CODE */
				"FQ:E:",options,&opt_index);
		if (c == -1)break;
		switch (c)
		{
			case 0x626e6368: /* bnch */
				RSBENCH_STDOUT("# --bench flag chosen. This is short for e.g. -oa -Ob and ...\n");
				program_operation='b';
				program_not_chosen=0;
				goto qos;
			break;
			case 'F':
				/*
				 * Floating point mini-benchmark
				 * */
				return (rsb_lib_init(RSB_NULL_INIT_OPTIONS) == RSB_ERR_NO_ERROR && rsb__fp_benchmark() == RSB_ERR_NO_ERROR) ?RSB_PROGRAM_ERROR:RSB_PROGRAM_SUCCESS;
			break;
			case 'G':
				/*
				 * Sparse GEMM preliminary code.
				 * TODO: remove this temporary case, as it may break other functionality with the G flag.
				 * */
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__do_spgemm_test_code(argc-1,argv+1));
			case 'g':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb_genmm_main(argc,argv));
			break;
			case 0x006D6C73:
				return RSB_ERR_TO_PROGRAM_ERROR(rsb_mtx_ls_main(argc,argv,RSB_BOOL_FALSE));
			break;
			case 0x6c736c61:
				return RSB_ERR_TO_PROGRAM_ERROR(rsb_mtx_ls_main(argc,argv,RSB_BOOL_TRUE));
			break;
			case 0x7365:
#if RSB_HAVE_SETENV
				rsb__setenv(optarg);
				aoc+=2;
#else /* RSB_HAVE_SETENV */
				RSB_STDERR("No setenv() built in -- terminating.\n");
				return RSB_ERR_TO_PROGRAM_ERROR(RSB_ERR_GENERIC_ERROR);
#endif /* RSB_HAVE_SETENV */
			break;
			case 0x7270720a:
			case 0x72707268:
			goto qos;
			break;
			case 0x50505050:
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__dump_postscript(argc,argv));
			break;
			case 0x6c696d74:
			{
				if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
					goto berr;
				RSB_DO_ERROR_CUMULATE(errval,rsb_blas_limit_cases_tester());
				if((!RSB_WANT_PERMISSIVE_RSBENCH) && RSB_SOME_ERROR(errval))goto ferr;
				return RSB_ERR_TO_PROGRAM_ERROR(errval);
			}
			break;
			case 'E':
			{
				if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
					goto berr;
				RSB_DO_ERROR_CUMULATE(errval,rsb_blas_failure_tester(optarg));
				if((!RSB_WANT_PERMISSIVE_RSBENCH) && RSB_SOME_ERROR(errval))goto ferr;
				return RSB_ERR_TO_PROGRAM_ERROR(errval);
			}
			case 'Q':
				to.mtt = rsb__util_atof(optarg);
				if(strstr(optarg,"R")!=NULL)to.rrm=RSB_BOOL_TRUE;
				if(strstr(optarg,"U")!=NULL)to.tur=RSB_BOOL_TRUE;
				if(strstr(optarg,"Q")!=NULL)to.wqt=RSB_BOOL_TRUE;
				if(strstr(optarg,"v")!=NULL)to.wvr=RSB_BOOL_TRUE;
				if(strstr(optarg,"q")!=NULL)to.wqc=RSB_BOOL_TRUE;
				if(strstr(optarg,"C")!=NULL)to.wcs=RSB_BOOL_TRUE;
				if(strstr(optarg,"L")!=NULL)to.wnl=RSB_BOOL_TRUE;
			case 'B':
			RSB_SIGHR
			/* Sparse BLAS test.  */
//#if RSB_WITH_SPARSE_BLAS_INTERFACE 
#if 1
			{
				if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
					goto berr;
#if RSB_ALLOW_INTERNAL_GETENVS
				if(getenv("RSB_RSBENCH_BBMB") && rsb__util_atoi(getenv("RSB_RSBENCH_BBMB")) )
					goto bbmb;
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
				RSB_DO_ERROR_CUMULATE(errval,rsb__internals_check());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_sppsp());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_csr2coo());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_swt());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_idx());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_coo());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_asm(&to));
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_trtr(&to));
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_sort(&to));
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_gen(&to));
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_srt(&to));
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_srt2p());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__chk_permute());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				RSB_DO_ERROR_CUMULATE(errval,rsb__util_testall());
				if(/*(!RSB_WANT_PERMISSIVE_RSBENCH) && */RSB_SOME_ERROR(errval))goto ferr;
				if(to.wnl!=RSB_BOOL_TRUE)
					RSB_DO_ERROR_CUMULATE(errval,rsb_blas_runtime_limits_tester());
				if((!RSB_WANT_PERMISSIVE_RSBENCH) && RSB_SOME_ERROR(errval))goto ferr;
#if RSB_WANT_DO_LOCK_TEST
				/* may move this to a separate place one day */
				RSB_DO_ERROR_CUMULATE(errval,rsb__do_lock_test());
#endif /* RSB_WANT_DO_LOCK_TEST */

				RSB_DO_ERROR_CUMULATE(errval,rsb_blas_mini_tester(&to));
				if((!RSB_WANT_PERMISSIVE_RSBENCH) && RSB_SOME_ERROR(errval))goto ferr;
				//RSB_LIKWID_MARKER_INIT;
				//RSB_LIKWID_MARKER_R_START("RSB-QUICKTEST");
#if RSB_ALLOW_INTERNAL_GETENVS
bbmb:
#endif /* RSB_ALLOW_INTERNAL_GETENVS */
				RSB_DO_ERROR_CUMULATE(errval,rsb_blas_bigger_matrices_tester(&to));/* TODO: options should be passed here */
				//RSB_LIKWID_MARKER_R_STOP("RSB-QUICKTEST");
				//RSB_LIKWID_MARKER_EXIT;
				if((!RSB_WANT_PERMISSIVE_RSBENCH) && RSB_SOME_ERROR(errval))
					goto ferr;

				RSB_DO_ERROR_CUMULATE(errval,rsb_lib_exit(RSB_NULL_EXIT_OPTIONS));
				if(RSB_SOME_ERROR(errval))
				{
					RSB_PERR_GOTO(berr,RSB_ERRM_ES);
				}
				goto ferr; /* ferr will rsb_lib_exit again, and it's expected to work. */
			}
#else
				RSB_STDERR("no Sparse BLAS interface built.\n");
				//return -1;
				return 0;
#endif
			break;
			case 'C':
				if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
					RSB_STDERR(RSB_ERRM_SILTC);
				errval |= rsb__print_configuration_string_rsbench(argv[0],cs,RSB_BOOL_TRUE);
				printf("%s",cs);
				goto verr;
			break;
			case 'v':
				if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
					RSB_STDERR(RSB_ERRM_SILTC);
				errval |= rsb__print_configuration_string_rsbench(argv[0],cs,RSB_BOOL_FALSE);
				printf("%s",cs);
				goto verr;
			break;
			case 'M':
				if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__memory_benchmark(NULL));
			break;
			case 'H':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb_hc_main());		/* preliminary */
			break;
			case 'I':
				rsb_lib_init(RSB_NULL_INIT_OPTIONS);
			return
				RSB_ERR_TO_PROGRAM_ERROR(rsb_perror(NULL,rsb__sys_info()));
			break;
			/*
#if RSB_WANT_EXPERIMENTS_CODE
			case 'e':
				return rsb_exp_bcsr_guess_experiments(argc,argv);
			break;
#endif */ /* RSB_WANT_EXPERIMENTS_CODE */
			case 'P':
		{
			errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS);
			if(RSB_SOME_ERROR(errval))
				goto err;
			else
			{
				struct rsb_mtx_t * mtxAp = rsb_file_mtx_load(optarg,RSB_FLAG_DEFAULT_MATRIX_FLAGS,RSB_NUMERICAL_TYPE_DEFAULT,NULL);
				errval = rsb_file_mtx_save(mtxAp,NULL);
				RSB_MTX_FREE(mtxAp);
			}
			RSB_MASK_OUT_SOME_ERRORS(errval)
			goto verr;
		}
			break;
			case 'o':
				operation=*optarg;
			break;
			case 'O':
				program_operation=*optarg;
				program_not_chosen=0;
			break;
			case 'h':
				/* getsubopt may come in help here */
				return rsb__main_help(argc, argv,default_program_operation,program_codes,options);
			break;
			/*
			case 't':
				return rsb__main_transpose(argc,argv);
			break;	    	
			*/
			default:
			{
			}
		}
	}

qos:	/* quit option selection */

	if(c == 0x72707268)
	{
		errval = rsb__pr_dumpfiles(NULL,0);
		return RSB_ERR_TO_PROGRAM_ERROR(errval);
	}
	if(c == 0x7270720a)
	{
/*
		if(argc == 3)
			return rsb__pr_dumpfile(optarg);
		if(argc == 3)
			return rsb__pr_dumpfiles(&optarg,1);
 */

		if(argc-aoc >= 3)
		{
			const int RSB__PR_DUMP_MAXARGS = 1024; /* TODO: temporary */
			const rsb_char_t*fna[RSB__PR_DUMP_MAXARGS];
			int i;
			for(i=2+aoc;i<RSB_MIN(RSB__PR_DUMP_MAXARGS,argc);++i)
				fna[i-2-aoc] = argv[i];
			errval = rsb__pr_dumpfiles(fna,i-2-aoc);
			return RSB_ERR_TO_PROGRAM_ERROR(errval);
		}
	}

	if(program_not_chosen)
		return rsb__main_help(argc, argv,default_program_operation,program_codes,options);

	switch (program_operation)	/* O */
	{
		case 'r':
	{
			/*
			 * A benchmark to compute (machine,compiled) reference performance values.
			 * */
			errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS);
			if(errval != RSB_ERR_NO_ERROR)
				goto verr;
			errval = rsb__do_referencebenchmark(); /* primarily kept for test purposes */
			RSB_MASK_OUT_SOME_ERRORS(errval)
			goto verr;
	}
		break;
		case 'R':
			/**/
			/*
			 * Dump current (hardcoded) performance info without computing anything.
			 * */
			errval = rsb__dump_current_global_reference_performance_info(); /* FIXME: probably obsolete */
			RSB_MASK_OUT_SOME_ERRORS(errval)
			goto verr;
		break;
#if !RSB_SHALL_UPDATE_COMPLETEBENCHS
		case 'c':
			errval = rsb__do_completebenchmark(argc,argv);
			RSB_MASK_OUT_SOME_ERRORS(errval)
			goto verr;
		break;
#endif /* RSB_SHALL_UPDATE_COMPLETEBENCHS */
		case 'd':
			/*
			 * A single matrix dump (almost useless).
			 * */
			return rsb__test_dump_main(argc,argv); /* FIXME: probably obsolete */
		break;
		case 'e':
			/*
			 * The matrix experimentation code.
			 * */
			RSB_STDERR("this option was obsoleted by -oS -Ob\n");
			/*return rsb_test_main_block_partitioned_matrix_stats(argc,argv); */ /* FIXME: probably obsolete */
			return -1;
		break;
		case 'b':
		{
			/*
			 * The current reference benchmark.
			 * */
			RSB_SIGHR
			switch(operation)	/* o */
			{
				case 'a':
#ifdef RSB_HAVE_OPTYPE_SPMV_UAUA
				//return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_spmv_uaua(argc,argv));
#endif /* RSB_HAVE_OPTYPE_SPMV_UAUA */
#ifdef RSB_HAVE_OPTYPE_SPMV_SXSA
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_spmv_sxsa(argc,argv));
#endif /* RSB_HAVE_OPTYPE_SPMV_UAUA */
				break;
#if !RSB_WANT_REDUCED_RSB_M4_MATRIX_META_OPS
#ifdef RSB_HAVE_OPTYPE_SPMV_UAUZ
				case 'v':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_spmv_uauz(argc,argv));
				break;
#endif /* RSB_HAVE_OPTYPE_SPMV_UAUZ */
#ifdef RSB_HAVE_OPTYPE_SPMM_AZ
				case 'm':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_spmm_az(argc,argv));
				break;
#endif /* RSB_HAVE_OPTYPE_SPMM_AZ */
#ifdef RSB_HAVE_OPTYPE_SCALE
				case 's':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_scale(argc,argv));
				break;
#endif /* RSB_HAVE_OPTYPE_SCALE */
#ifdef RSB_HAVE_OPTYPE_SPMV_UXUX
				case 'c':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_spmv_uxux(argc,argv));
				break;
#endif /* RSB_HAVE_OPTYPE_SPMV_UXUX */
#ifdef RSB_HAVE_OPTYPE_INFTY_NORM
				case 'i':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_infty_norm(argc,argv));
				break;
#endif /* RSB_HAVE_OPTYPE_INFTY_NORM */
#ifdef RSB_HAVE_OPTYPE_NEGATION
				case 'n':
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_negation(argc,argv));
				break;
#endif /* RSB_HAVE_OPTYPE_NEGATION */
#endif /* RSB_WANT_REDUCED_RSB_M4_MATRIX_META_OPS */
				case 't':
#ifdef RSB_OPTYPE_INDEX_SPSV_UXUA
				//return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_spsv_uxua(argc,argv));
#endif /* RSB_OPTYPE_INDEX_SPSV_UXUA */
#ifdef RSB_OPTYPE_INDEX_SPSV_SXSX
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_spsv_sxsx(argc,argv));
#endif /* RSB_OPTYPE_INDEX_SPSV_UXUA */
				break;
#if 1	/* this is a special case */
				case 'S':
				//return rsb__main_block_partitioned_sort_only(argc,argv);//old
				return RSB_ERR_TO_PROGRAM_ERROR(rsb__main_block_partitioned_mat_stats(argc,argv));//new
				break;
#endif
				default:
				RSB_STDERR(
					"You did not choose a correct operation code.\n"
					"Choose one among %s.\n",program_codes
					);
				errval = RSB_ERR_UNSUPPORTED_OPERATION;
				RSB_DO_ERR_RETURN(errval)
			}
		}
		break;
#if 0
		case 't': /* to reintegrate, add 't' to program_codes */
			/*
			 * A whole matrix repartitioning test.
			 * */
			return RSB_ERR_TO_PROGRAM_ERROR(rsb_test_main_block_partitioned_construction_test(argc,argv));
		break;
#endif
		default:
			RSB_STDERR("You did not choose an action. See help:\n");
			return rsb__main_help(argc, argv,default_program_operation,program_codes,NULL);
		return RSB_PROGRAM_SUCCESS;
    	}
	goto err;
ferr:
	/* rsb__getrusage(); */
	RSB_DO_ERROR_CUMULATE(errval,rsb_lib_exit(RSB_NULL_EXIT_OPTIONS));
berr:
	if(RSB_SOME_ERROR(errval))
		rsb__do_perror(NULL,errval);
verr:
	return RSB_ERR_TO_PROGRAM_ERROR(errval);
err:
	if(RSB_SOME_ERROR(errval))
		rsb__do_perror(NULL,errval);
	return RSB_PROGRAM_ERROR;
}

/* @endcond */
