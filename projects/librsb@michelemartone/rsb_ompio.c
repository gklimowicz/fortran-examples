/* @cond INNERDOC */
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block coordinates format.
 Kernels unrolled, with no loops, for only user-specified blockings.
 */

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
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

 */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file is an adaptation of Gilles Gouaillardet OpenMP+fgets_unlocked suggested implementation
 * This code is still experimental and untested.
 */
#include "rsb_internals.h"
#include "rsb_lock.h"

#ifdef RSB_HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif /* RSB_HAVE_SYS_STAT_H */
#if RSB_WANT_ZLIB_SUPPORT
#include <zlib.h>
#endif /* RSB_WANT_ZLIB_SUPPORT */
#include <stdio.h>
#include "rsb_internals.h"
#include "rsb_ompio.h"
RSB_INTERNALS_COMMON_HEAD_DECLS
#if RSB_WANT_OMPIO_SUPPORT
void rsb_ompio_DOUBLE (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,double**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re){
	rsb_nnz_idx_t i=0;
# pragma omp parallel RSB_NTC
{
	int index = 0, j=0, toread=0;
	size_t res = 0;
	char line[RSB_MAX_FGETS_LINES][MM_MAX_LINE_LENGTH];
	while (index < *nnz) {
#pragma omp critical
{
		index=i;
		if (i < *nnz) {
			toread = ((*nnz-index)>RSB_MAX_FGETS_LINES)?RSB_MAX_FGETS_LINES:(*nnz-index);
			for (j=0;j<toread;j++)
#if RSB_WANT_ZLIB_SUPPORT
				if (!ngzfd)
					gzgets(fd,line[j],MM_MAX_LINE_LENGTH);
                                else
#endif /* RSB_WANT_ZLIB_SUPPORT */
	   				fgets_unlocked(line[j],MM_MAX_LINE_LENGTH,ngzfd);
			i += toread;
		}
}
		if (index < *nnz) {
			for(j=0;j<toread;j++,index++) {
			int iI,iJ;
			char * p1, *p2;
			p1=line[j];

			iI=strtol(p1,&p2,10);
			iJ=strtol(p2,&p1,10);
*(*dval+index)=strtod(p1,NULL);			res += 1;
			(*IA)[index]=(rsb_coo_idx_t)iI;
			(*JA)[index]=(rsb_coo_idx_t)iJ;
	       		(*IA)[index]--;  // adjust from 1-based to 0-based
        		(*JA)[index]--;
			}
		}
	}
# pragma omp critical
{
	*_re+=res;
}
}
}

void rsb_ompio_FLOAT (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,float**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re){
	rsb_nnz_idx_t i=0;
# pragma omp parallel RSB_NTC
{
	int index = 0, j=0, toread=0;
	size_t res = 0;
	char line[RSB_MAX_FGETS_LINES][MM_MAX_LINE_LENGTH];
	while (index < *nnz) {
#pragma omp critical
{
		index=i;
		if (i < *nnz) {
			toread = ((*nnz-index)>RSB_MAX_FGETS_LINES)?RSB_MAX_FGETS_LINES:(*nnz-index);
			for (j=0;j<toread;j++)
#if RSB_WANT_ZLIB_SUPPORT
				if (!ngzfd)
					gzgets(fd,line[j],MM_MAX_LINE_LENGTH);
                                else
#endif /* RSB_WANT_ZLIB_SUPPORT */
	   				fgets_unlocked(line[j],MM_MAX_LINE_LENGTH,ngzfd);
			i += toread;
		}
}
		if (index < *nnz) {
			for(j=0;j<toread;j++,index++) {
			int iI,iJ;
			char * p1, *p2;
			p1=line[j];

			iI=strtol(p1,&p2,10);
			iJ=strtol(p2,&p1,10);
*(*dval+index)=strtof(p1,NULL);			res += 1;
			(*IA)[index]=(rsb_coo_idx_t)iI;
			(*JA)[index]=(rsb_coo_idx_t)iJ;
	       		(*IA)[index]--;  // adjust from 1-based to 0-based
        		(*JA)[index]--;
			}
		}
	}
# pragma omp critical
{
	*_re+=res;
}
}
}

void rsb_ompio_FLOAT_COMPLEX (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,float complex**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re){
	rsb_nnz_idx_t i=0;
# pragma omp parallel RSB_NTC
{
	int index = 0, j=0, toread=0;
	size_t res = 0;
	char line[RSB_MAX_FGETS_LINES][MM_MAX_LINE_LENGTH];
	while (index < *nnz) {
#pragma omp critical
{
		index=i;
		if (i < *nnz) {
			toread = ((*nnz-index)>RSB_MAX_FGETS_LINES)?RSB_MAX_FGETS_LINES:(*nnz-index);
			for (j=0;j<toread;j++)
#if RSB_WANT_ZLIB_SUPPORT
				if (!ngzfd)
					gzgets(fd,line[j],MM_MAX_LINE_LENGTH);
                                else
#endif /* RSB_WANT_ZLIB_SUPPORT */
	   				fgets_unlocked(line[j],MM_MAX_LINE_LENGTH,ngzfd);
			i += toread;
		}
}
		if (index < *nnz) {
			for(j=0;j<toread;j++,index++) {
			int iI,iJ;
			char * p1, *p2;
			p1=line[j];

			iI=strtol(p1,&p2,10);
			iJ=strtol(p2,&p1,10);
*(((double*)(*dval))+2*index+0)=strtof(p1,&p2);
*(((double*)(*dval))+2*index+1)=strtof(p2,&p1);			res += 1;
			(*IA)[index]=(rsb_coo_idx_t)iI;
			(*JA)[index]=(rsb_coo_idx_t)iJ;
	       		(*IA)[index]--;  // adjust from 1-based to 0-based
        		(*JA)[index]--;
			}
		}
	}
# pragma omp critical
{
	*_re+=res;
}
}
}

void rsb_ompio_DOUBLE_COMPLEX (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd,double complex**dval, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re){
	rsb_nnz_idx_t i=0;
# pragma omp parallel RSB_NTC
{
	int index = 0, j=0, toread=0;
	size_t res = 0;
	char line[RSB_MAX_FGETS_LINES][MM_MAX_LINE_LENGTH];
	while (index < *nnz) {
#pragma omp critical
{
		index=i;
		if (i < *nnz) {
			toread = ((*nnz-index)>RSB_MAX_FGETS_LINES)?RSB_MAX_FGETS_LINES:(*nnz-index);
			for (j=0;j<toread;j++)
#if RSB_WANT_ZLIB_SUPPORT
				if (!ngzfd)
					gzgets(fd,line[j],MM_MAX_LINE_LENGTH);
                                else
#endif /* RSB_WANT_ZLIB_SUPPORT */
	   				fgets_unlocked(line[j],MM_MAX_LINE_LENGTH,ngzfd);
			i += toread;
		}
}
		if (index < *nnz) {
			for(j=0;j<toread;j++,index++) {
			int iI,iJ;
			char * p1, *p2;
			p1=line[j];

			iI=strtol(p1,&p2,10);
			iJ=strtol(p2,&p1,10);
*(((double*)(*dval))+2*index+0)=strtod(p1,&p2);
*(((double*)(*dval))+2*index+1)=strtod(p2,&p1);			res += 1;
			(*IA)[index]=(rsb_coo_idx_t)iI;
			(*JA)[index]=(rsb_coo_idx_t)iJ;
	       		(*IA)[index]--;  // adjust from 1-based to 0-based
        		(*JA)[index]--;
			}
		}
	}
# pragma omp critical
{
	*_re+=res;
}
}
}

void rsb_ompio_PATTERN (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re){
	rsb_nnz_idx_t i=0;
# pragma omp parallel RSB_NTC
{
	int index = 0, j=0, toread=0;
	size_t res = 0;
	char line[RSB_MAX_FGETS_LINES][MM_MAX_LINE_LENGTH];
	while (index < *nnz) {
#pragma omp critical
{
		index=i;
		if (i < *nnz) {
			toread = ((*nnz-index)>RSB_MAX_FGETS_LINES)?RSB_MAX_FGETS_LINES:(*nnz-index);
			for (j=0;j<toread;j++)
#if RSB_WANT_ZLIB_SUPPORT
				if (!ngzfd)
					gzgets(fd,line[j],MM_MAX_LINE_LENGTH);
                                else
#endif /* RSB_WANT_ZLIB_SUPPORT */
	   				fgets_unlocked(line[j],MM_MAX_LINE_LENGTH,ngzfd);
			i += toread;
		}
}
		if (index < *nnz) {
			for(j=0;j<toread;j++,index++) {
			int iI,iJ;
			char * p1, *p2;
			p1=line[j];

			iI=strtol(p1,&p2,10);
			iJ=strtol(p2,&p1,10);
			res += 1;
			(*IA)[index]=(rsb_coo_idx_t)iI;
			(*JA)[index]=(rsb_coo_idx_t)iJ;
	       		(*IA)[index]--;  // adjust from 1-based to 0-based
        		(*JA)[index]--;
			}
		}
	}
# pragma omp critical
{
	*_re+=res;
}
}
}


#endif	/* RSB_WANT_OMPIO_SUPPORT */

/* @endcond */
