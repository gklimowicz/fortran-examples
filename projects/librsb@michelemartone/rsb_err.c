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
/**
 * @file
 * @author Michele Martone
 * @brief
 * */
#include "rsb_common.h"
#include "rsb_util.h"
#include "rsb.h"

#define rsb__snprintf snprintf

RSB_INTERNALS_COMMON_HEAD_DECLS

const char *rsb__get_errstr_ptr(rsb_err_t errval)
{
	/* return pointer to const error message string area */
	const char *s = NULL;

	switch(errval)
	{
		case RSB_ERR_GENERIC_ERROR:
			s = RSB_ERRM_GENERIC_ERROR ;
		break;
		case RSB_ERR_UNSUPPORTED_OPERATION:
			s = RSB_ERRM_UNSUPPORTED_OPERATION;
		break;
		case RSB_ERR_UNSUPPORTED_TYPE:
			s = RSB_ERRM_UNSUPPORTED_TYPE;
		break;
		case RSB_ERR_UNSUPPORTED_FORMAT:
			s = RSB_ERRM_UNSUPPORTED_FORMAT;
		break;
		case RSB_ERR_INTERNAL_ERROR:
			s = RSB_ERRM_INTERNAL_ERROR;
		break;
		case RSB_ERR_BADARGS:
			s = RSB_ERRM_BADARGS;
		break;
		case RSB_ERR_ENOMEM:
			s = RSB_ERRM_ENOMEM;
		break;
		case RSB_ERR_UNIMPLEMENTED_YET:
			s = RSB_ERRM_UNIMPLEMENTED_YET;
		break;
		case RSB_ERR_LIMITS:
			s = RSB_ERRM_LIMITS;
		break;
		case RSB_ERR_NO_USER_CONFIGURATION:
			s = RSB_ERRM_NO_USER_CONFIGURATION;
		break;
		case RSB_ERR_CORRUPT_INPUT_DATA:
			s = RSB_ERRM_CORRUPT_INPUT_DATA;
		break;
		case RSB_ERR_FAILED_MEMHIER_DETECTION:
			s = RSB_ERRM_FAILED_MEMHIER_DETECTION;
		break;
		case RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS:
			s = RSB_ERRM_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS;
		break;
		case RSB_ERR_UNSUPPORTED_FEATURE:
			s = RSB_ERRM_UNSUPPORTED_FEATURE;
		break;
		case RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT:
			s = RSB_ERRM_NO_STREAM_OUTPUT_CONFIGURED_OUT;
		break;
		case RSB_ERR_INVALID_NUMERICAL_DATA:
			s = RSB_ERRM_INVALID_NUMERICAL_DATA;
		break;
		case RSB_ERR_MEMORY_LEAK:
			s = RSB_ERRM_MEMORY_LEAK;
		break;
		case RSB_ERR_ELEMENT_NOT_FOUND:
			s = RSB_ERRM_ELEMENT_NOT_FOUND;
		break;
		/*
		case RSB_ERR_FORTRAN_ERROR:
		s = "A Fortran-specific error occurred.";
		break;
		*/
		default:
		s = RSB_ERRM_ES;
		// might also consider "Unknown error code (%x)"
	}
	RSB_DEBUG_ASSERT(s);
	return s;
}

rsb_err_t rsb__do_strerror_r(rsb_err_t errval, rsb_char_t * buf, size_t buflen)
{
	/* TODO: what if buflen is not enough ? shall report this somehow. */
	rsb_char_t*sbuf = buf;
	const rsb_char_t *s = "No error occurred (success). The return value that means function operation success, in most cases.\n";

	if( errval == RSB_ERR_NO_ERROR )
		goto err;

	if( buf == NULL)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	s = rsb__get_errstr_ptr(errval);

	errval = RSB_ERR_NO_ERROR;
	rsb__snprintf(sbuf,buflen,"%s",s);
err:
	return errval;
}

rsb_err_t rsb__do_perror(FILE *stream, rsb_err_t errval)
{
	/*!
	 * \ingroup gr_internals
	 * Stateless function.
	 */
	rsb_char_t sbuf[RSB_MAX_STRERRLEN];

	if( errval == RSB_ERR_NO_ERROR )
		goto err;

	rsb__do_strerror_r(errval,sbuf,sizeof(sbuf)/sizeof(sbuf[0]));
	 
	if(stream)
		fprintf(stream,"ERROR 0x%x : %s\n",(unsigned int)errval,sbuf);
	else
		RSB_STDERR("ERROR 0x%x : %s\n",(unsigned int)errval,sbuf);
err:
	RSB_DO_ERR_RETURN(RSB_ERR_NO_ERROR)
}

/* @endcond */
