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
/* @cond INNERDOC  */
/**
 * @file
 * @author Michele Martone
 * @brief Matrix Market IA/O library for ANSI C.
 */
/*
 *   See http://math.nist.gov/MatrixMarket for details.
 *
 *   FIXME : remove dangling printf's
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "rsb_mmio.h"
#include "rsb_mio.h"

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static int rsb__mm_is_valid(MM_typecode matcode)
{
    if (!rsb_mm_is_matrix(matcode)) return 0;
    if (rsb_mm_is_dense(matcode) && rsb_mm_is_pattern(matcode)) return 0;
    if (rsb_mm_is_real(matcode) && rsb_mm_is_hermitian(matcode)) return 0;
    if (rsb_mm_is_pattern(matcode) && (rsb_mm_is_hermitian(matcode) || 
                rsb_mm_is_skew(matcode))) return 0;
    return 1;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

//int rsb__mm_read_banner(FILE *f, MM_typecode *matcode)
int rsb__mm_read_banner(FILE *f, FILE * ngzfd, MM_typecode *matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtxs[MM_MAX_TOKEN_LENGTH]; 
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;


    rsb_mm_clear_typecode(matcode);  

    if(ngzfd)
    {
    if (((char*)rsb__fgets(line,MM_MAX_LINE_LENGTH,ngzfd)) == NULL) /* stupid cast for PGI hairiness */
        return MM_PREMATURE_EOF;
    }
    else
    {
    		if (((char*)fgets(line,MM_MAX_LINE_LENGTH,f)) == NULL) /* stupid cast for PGI hairiness */
        return MM_PREMATURE_EOF;
    }

    if (sscanf(line, "%s %s %s %s %s", banner, mtxs, crd, data_type, 
        storage_scheme) != 5)
    {
	if ( line[0]==0x1f && line[1]==(char)0x8b )
        	return MM_LIKELY_GZIPPED_FILE;
        return MM_PREMATURE_EOF;
    }

    for (p=mtxs; *p!='\0'; *p=tolower(*p),p++);  /* convert to lower case */
    for (p=crd; *p!='\0'; *p=tolower(*p),p++);  
    for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtxs, MM_MTX_STR) != 0)
        return  MM_UNSUPPORTED_TYPE;
    rsb_mm_set_matrix(matcode);


    /* second field describes whether this is a sparse matrix (in coordinate
            storgae) or a dense array */


    if (strcmp(crd, MM_SPARSE_STR) == 0)
        rsb_mm_set_sparse(matcode);
    else
    if (strcmp(crd, MM_DENSE_STR) == 0)
            rsb_mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
    

    /* third field */

    if (strcmp(data_type, MM_REAL_STR) == 0)
        rsb_mm_set_real(matcode);
    else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        rsb_mm_set_complex(matcode);
    else
    if (strcmp(data_type, MM_PATTERN_STR) == 0)
        rsb_mm_set_pattern(matcode);
    else
    if (strcmp(data_type, MM_INT_STR) == 0)
        rsb_mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
    

    /* fourth field */

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        rsb_mm_set_general(matcode);
    else
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        rsb_mm_set_symmetric(matcode);
    else
    if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        rsb_mm_set_hermitian(matcode);
    else
    if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        rsb_mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
        

    return 0;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)
{
    if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
        return MM_COULD_NOT_WRITE_FILE;
    else 
        return 0;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

//int rsb__mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz )
int rsb__mm_read_mtx_crd_size(FILE *f, FILE * ngzfd, rsb_coo_idx_t *M, rsb_coo_idx_t *N, rsb_coo_idx_t *nz)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;
    long int lM,lN,lnz;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;

    /* now continue scanning until you reach the end-of-comments */
    do 
    {
	    if(ngzfd)
	    {
    		if (((char*)rsb__fgets(line,MM_MAX_LINE_LENGTH,ngzfd)) == NULL) /* stupid cast for PGI hairiness */
	            return MM_PREMATURE_EOF;
	    }
	    else
	    {
    		if (((char*)fgets(line,MM_MAX_LINE_LENGTH,f)) == NULL) /* stupid cast for PGI hairiness */
	            return MM_PREMATURE_EOF;
	    }
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%ld %ld %ld", &lM, &lN, &lnz) == 3)
    {
        *M=lM;
        *N=lN;
        *nz=lnz;
        return 0;
    }
        
    else
    do
    { 
#ifdef RSB_WANT_LONG_IDX_TYPE 
	    if(ngzfd)
        num_items_read = rsb__fscanf(ngzfd,"%ld %ld %ld",&lM,&lN,&lnz,NULL);
	    else
        num_items_read = fscanf(f,"%ld %ld %ld",&lM,&lN,&lnz);
	*M=lM,*N=lN,*nz=lnz;
#else /* RSB_WANT_LONG_IDX_TYPE */
	    if(ngzfd)
        num_items_read = rsb__fscanf(ngzfd,"%d %d %d",M,N,nz,NULL);
	    else
        num_items_read = fscanf(f,"%d %d %d",M,N,nz); 
#endif /* RSB_WANT_LONG_IDX_TYPE */
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 3);

    return 0;
}


//int rsb__mm_read_mtx_array_size(FILE *f, int *M, int *N)
int rsb__mm_read_mtx_array_size(FILE *f, FILE *ngzfd, rsb_coo_idx_t*M, rsb_coo_idx_t*N)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;
    /* set return null parameter values, in case we exit with errors */
#ifdef RSB_WANT_LONG_IDX_TYPE
    long int lM,lN;
#endif /* RSB_WANT_LONG_IDX_TYPE */

    *M = *N = 0;
	
    /* now continue scanning until you reach the end-of-comments */
    if(ngzfd)
    {
    do 
    {
        if ((char*)rsb__fgets(line,MM_MAX_LINE_LENGTH,ngzfd) == NULL) 
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
#ifdef RSB_WANT_LONG_IDX_TYPE 
    if (sscanf(line, "%ld %ld", &lM, &lN) == 2)
    {
	*M=lM;
	*N=lN;
        return 0;
    }
#else /* RSB_WANT_LONG_IDX_TYPE */
    if (sscanf(line, "%d %d", M, N) == 2)
        return 0;
#endif /* RSB_WANT_LONG_IDX_TYPE */
        
    else /* we have a blank line */
    do
    { 
#ifdef RSB_WANT_LONG_IDX_TYPE 
        num_items_read = rsb__fscanf(ngzfd, "%ld %ld", M, N, NULL, NULL);
#else /* RSB_WANT_LONG_IDX_TYPE */
        num_items_read = rsb__fscanf(ngzfd, "%d %d", M, N, NULL, NULL);
#endif /* RSB_WANT_LONG_IDX_TYPE */
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 2);
    }
    else
    if(f)
    {
    do 
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
#ifdef RSB_WANT_LONG_IDX_TYPE 
    if (sscanf(line, "%ld %ld", &lM, &lN) == 2)
    {
        *M=lM;
        *N=lN;
        return 0;
    }
#else /* RSB_WANT_LONG_IDX_TYPE */
    if (sscanf(line, "%d %d", M, N) == 2)
        return 0;
#endif /* RSB_WANT_LONG_IDX_TYPE */
        
    else /* we have a blank line */
    do
    { 
        /* blank lines in header are tolerated */
#ifdef RSB_WANT_LONG_IDX_TYPE 
        num_items_read = fscanf(f, "%ld %ld", &lM, &lN);
        *M=lM;
        *N=lN;
#else /* RSB_WANT_LONG_IDX_TYPE */
        num_items_read = fscanf(f, "%d %d", M, N);
#endif /* RSB_WANT_LONG_IDX_TYPE */
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 2);
    }

    return 0;
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__mm_write_mtx_array_size(FILE *f, int M, int N)
{
    if (fprintf(f, "%d %d\n", M, N) != 2)
        return MM_COULD_NOT_WRITE_FILE;
    else 
        return 0;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */



/*-------------------------------------------------------------------------*/

/******************************************************************/
/* use when IA[], JA[], and VA[]JA, and VA[] are already allocated */
/******************************************************************/

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, rsb_coo_idx_t IA[], rsb_coo_idx_t JA[],
        double VA[], MM_typecode matcode)
{
    int i;

    if (rsb_mm_is_complex(matcode))
    {
#ifdef RSB_WANT_LONG_IDX_TYPE 
        for (i=0; i<nz; i++)
            if (fscanf(f, "%ld %ld %lg %lg", &IA[i], &JA[i], &VA[2*i], &VA[2*i+1])
                != 4) return MM_PREMATURE_EOF;
#else /* RSB_WANT_LONG_IDX_TYPE */
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d %lg %lg", &IA[i], &JA[i], &VA[2*i], &VA[2*i+1])
                != 4) return MM_PREMATURE_EOF;
#endif /* RSB_WANT_LONG_IDX_TYPE */
    }
    else if (rsb_mm_is_real(matcode))
    {
#ifdef RSB_WANT_LONG_IDX_TYPE 
        for (i=0; i<nz; i++)
        {
            if (fscanf(f, "%ld %ld %lg\n", &IA[i], &JA[i], &VA[i])
                != 3) return MM_PREMATURE_EOF;

        }
#else /* RSB_WANT_LONG_IDX_TYPE */
        for (i=0; i<nz; i++)
        {
            if (fscanf(f, "%d %d %lg\n", &IA[i], &JA[i], &VA[i])
                != 3) return MM_PREMATURE_EOF;

        }
#endif /* RSB_WANT_LONG_IDX_TYPE */
    }

    else if (rsb_mm_is_pattern(matcode))
    {
#ifdef RSB_WANT_LONG_IDX_TYPE 
        for (i=0; i<nz; i++)
            if (fscanf(f, "%ld %ld", &IA[i], &JA[i])
                != 2) return MM_PREMATURE_EOF;
#else /* RSB_WANT_LONG_IDX_TYPE */
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d", &IA[i], &JA[i])
                != 2) return MM_PREMATURE_EOF;
#endif /* RSB_WANT_LONG_IDX_TYPE */
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static int rsb__mm_read_mtx_crd_entry(FILE *f, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA,
        double *real, double *imag, MM_typecode matcode)
{
    int iI,iJ;
    if (rsb_mm_is_complex(matcode))
    {
            if (fscanf(f, "%d %d %lg %lg", &iI, &iJ, real, imag)
                != 4) return MM_PREMATURE_EOF;
		IA[0]=(rsb_coo_idx_t)iI;
		JA[0]=(rsb_coo_idx_t)iJ;
    }
    else if (rsb_mm_is_real(matcode))
    {
            if (fscanf(f, "%d %d %lg\n", &iI, &iJ, real)
                != 3) return MM_PREMATURE_EOF;
		IA[0]=(rsb_coo_idx_t)iI;
		JA[0]=(rsb_coo_idx_t)iJ;
    }

    else if (rsb_mm_is_pattern(matcode))
    {
            if (fscanf(f, "%d %d", &iI, &iJ ) != 2) return MM_PREMATURE_EOF;
		IA[0]=(rsb_coo_idx_t)iI;
		JA[0]=(rsb_coo_idx_t)iJ;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */


/************************************************************************
    rsb__mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                        type code, e.g. 'MCRS'

                        if matrix is complex, values[] is of size 2*nz,
                            (nz pairs of real/imaginary values)
************************************************************************/

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static int rsb__mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, rsb_coo_idx_t **IA, rsb_coo_idx_t **JA, 
        double **VA, MM_typecode *matcode)
{
    int ret_code;
    FILE *f;

    if (strcmp(fname, "stdin") == 0) f=stdin;
    else
    if ((f = fopen(fname, "r")) == NULL)
        return MM_COULD_NOT_READ_FILE;


    if ((ret_code = rsb__mm_read_banner(NULL,f, matcode)) != 0)
        return ret_code;

    if (!(rsb__mm_is_valid(*matcode) && rsb_mm_is_sparse(*matcode) && 
            rsb_mm_is_matrix(*matcode)))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = rsb__mm_read_mtx_crd_size(f,NULL,M,N,nz)) != 0)
        return ret_code;


    *IA = (rsb_coo_idx_t *)  malloc(*nz * sizeof(rsb_coo_idx_t));
    *JA = (rsb_coo_idx_t *)  malloc(*nz * sizeof(rsb_coo_idx_t));
    *VA = NULL;

    if (rsb_mm_is_complex(*matcode))
    {
        *VA = (double *) malloc(*nz * 2 * sizeof(double));
        ret_code = rsb__mm_read_mtx_crd_data(f, *M, *N, *nz, *IA, *JA, *VA, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }
    else if (rsb_mm_is_real(*matcode))
    {
        *VA = (double *) malloc(*nz * sizeof(double));
        ret_code = rsb__mm_read_mtx_crd_data(f, *M, *N, *nz, *IA, *JA, *VA, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    else if (rsb_mm_is_pattern(*matcode))
    {
        ret_code = rsb__mm_read_mtx_crd_data(f, *M, *N, *nz, *IA, *JA, *VA, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    if (f != stdin) fclose(f);
    return 0;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__mm_write_banner(FILE *f, MM_typecode matcode)
{
    char *str = rsb__mm_typecode_to_str(matcode);
    int ret_code;

    ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
    free(str);
    if (ret_code !=2 )
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__mm_write_mtx_crd(char fname[], int M, int N, int nz, int IA[], int JA[],
        double VA[], MM_typecode matcode)
{
    FILE *f;
    int i;

    if (strcmp(fname, "stdout") == 0) 
        f = stdout;
    else
    if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;
    
    /* print banner followed by typecode */
    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", rsb__mm_typecode_to_str(matcode));

    /* print matrix sizes and nonzeros */
    fprintf(f, "%d %d %d\n", M, N, nz);

    /* print values */
    if (rsb_mm_is_pattern(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d\n", IA[i], JA[i]);
    else
    if (rsb_mm_is_real(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g\n", IA[i], JA[i], VA[i]);
    else
    if (rsb_mm_is_complex(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g %20.16g\n", IA[i], JA[i], VA[2*i], 
                        VA[2*i+1]);
    else
    {
        if (f != stdout) fclose(f);
        return MM_UNSUPPORTED_TYPE;
    }

    if (f !=stdout) fclose(f);

    return 0;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
  

/**
*  Create a new copy of a string s.  strdup() is a POSIX routine, 
*  not part of ANSI C, so it is included here.
*  Used by rsb__mm_typecode_to_str().
*/
static char *rsb__mm_strdup(const char *s)
{
	size_t len = strlen(s);
	char *s2 = (char *) malloc((len+1)*sizeof(char));
	return strcpy(s2, s);
}

char  *rsb__mm_typecode_to_str(MM_typecode matcode)
{
    char buffer[MM_MAX_LINE_LENGTH];
    char *types[4];

    /* check for MTX type */
    if (rsb_mm_is_matrix(matcode)) 
        types[0] = MM_MTX_STR;

    /* check for CRD or ARR matrix */
    if (rsb_mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else
    if (rsb_mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (rsb_mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else
    if (rsb_mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else
    if (rsb_mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else
    if (rsb_mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;


    /* check for symmetry type */
    if (rsb_mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else
    if (rsb_mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else 
    if (rsb_mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else 
    if (rsb_mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
    return rsb__mm_strdup(buffer);

}
/* @endcond */
