dnl
dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
include(`rsb_krnl_macros.m4')dnl
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file is an adaptation of Gilles Gouaillardet OpenMP+fgets_unlocked suggested implementation
 * This code is still experimental and untested.
 */
dnl RSB_M4_HEADER_MESSAGE()dnl
dnl
dnl
#include "rsb_internals.h"
#include "rsb_lock.h"
dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_OMPIO_H_INCLUDED
#define RSB_OMPIO_H_INCLUDED
')
dnl
dnl 
ifdef(`ONLY_WANT_HEADERS',`dnl
#define RSB_MAX_FGETS_LINES 1536
',`dnl
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
')dnl
dnl
`#if RSB_WANT_OMPIO_SUPPORT'
dnl
foreach(`mtype',(WANT_TYPES,pattern),`dnl
void rsb_ompio_`'touppercase(RSB_M4_CHOPSPACES(mtype))`' (rsb_nnz_idx_t *nnz, FILE * fd, FILE * ngzfd`'ifelse(mtype,`pattern',`',`,'mtype`**dval'), rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, size_t *_re)`'dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
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
ifelse(mtype,`pattern',`',`dnl
`'dnl
ifelse(mtype,`double',`*(*dval+index)=strtod(p1,NULL);',`dnl
ifelse(mtype,`float',`*(*dval+index)=strtof(p1,NULL);',`dnl
ifelse(mtype,`double complex',`dnl
*(((double*)(*dval))+2*index+0)=strtod(p1,&p2);
*(((double*)(*dval))+2*index+1)=strtod(p2,&p1);dnl
',`dnl
ifelse(mtype,`float complex',`dnl
*(((double*)(*dval))+2*index+0)=strtof(p1,&p2);
*(((double*)(*dval))+2*index+1)=strtof(p2,&p1);dnl
',`dnl
ifelse(mtype,`int',`*(*dval+index)=strtol(p1,NULL);',`dnl
ifelse(mtype,`char',`*(*dval+index)=strtol(p1,NULL);',`dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
')dnl
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

')dnl
')dnl
dnl

`#endif'	/* RSB_WANT_OMPIO_SUPPORT */
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
')dnl
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif	/* RSB_OMPIO_H_INCLUDED */
')
dnl
/* @endcond */
dnl
