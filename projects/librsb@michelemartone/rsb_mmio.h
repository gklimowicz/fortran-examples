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
 * @brief
 *   Matrix Market I/O library for ANSI C.
 *   See http://math.nist.gov/MatrixMarket for details.
 * @author Michele Martone
 */

#ifndef MM_IO_H_INCLUDED
#define MM_IO_H_INCLUDED

#include "rsb.h"	/* just for some type compatibility */

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[];

char *rsb__mm_typecode_to_str(MM_typecode matcode);

int rsb__mm_read_banner(FILE *f, FILE * ngzfd, MM_typecode *matcode);
int rsb__mm_read_mtx_crd_size(FILE *f, FILE * ngzfd, rsb_coo_idx_t *M, rsb_coo_idx_t *N, rsb_coo_idx_t *nz);
int rsb__mm_read_mtx_array_size(FILE *f, FILE *ngzfd, rsb_coo_idx_t*M, rsb_coo_idx_t*N);

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__mm_write_banner(FILE *f, MM_typecode matcode);
int rsb__mm_write_mtx_crd_size(FILE *f, int M, int N, int nz);
int rsb__mm_write_mtx_array_size(FILE *f, int M, int N);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */


/********************* MM_typecode query fucntions ***************************/

#define rsb_mm_is_matrix(typecode)	((typecode)[0]=='M')

#define rsb_mm_is_sparse(typecode)	((typecode)[1]=='C')
#define rsb_mm_is_coordinate(typecode)((typecode)[1]=='C')
#define rsb_mm_is_dense(typecode)	((typecode)[1]=='A')
#define rsb_mm_is_array(typecode)	((typecode)[1]=='A')

#define rsb_mm_is_complex(typecode)	((typecode)[2]=='C')
#define rsb_mm_is_real(typecode)		((typecode)[2]=='R')
#define rsb_mm_is_pattern(typecode)	((typecode)[2]=='P')
#define rsb_mm_is_integer(typecode) ((typecode)[2]=='I')

#define rsb_mm_is_symmetric(typecode)((typecode)[3]=='S')
#define rsb_mm_is_general(typecode)	((typecode)[3]=='G')
#define rsb_mm_is_skew(typecode)	((typecode)[3]=='K')
#define rsb_mm_is_hermitian(typecode)((typecode)[3]=='H')

int rsb__mm_is_valid(MM_typecode matcode);		/* too complex for a macro */


/********************* MM_typecode modify fucntions ***************************/

#define rsb_mm_set_matrix(typecode)	((*typecode)[0]='M')
#define rsb_mm_set_coordinate(typecode)	((*typecode)[1]='C')
#define rsb_mm_set_array(typecode)	((*typecode)[1]='A')
#define rsb_mm_set_dense(typecode)	rsb_mm_set_array(typecode)
#define rsb_mm_set_sparse(typecode)	rsb_mm_set_coordinate(typecode)

#define rsb_mm_set_complex(typecode)((*typecode)[2]='C')
#define rsb_mm_set_real(typecode)	((*typecode)[2]='R')
#define rsb_mm_set_pattern(typecode)((*typecode)[2]='P')
#define rsb_mm_set_integer(typecode)((*typecode)[2]='I')


#define rsb_mm_set_symmetric(typecode)((*typecode)[3]='S')
#define rsb_mm_set_general(typecode)((*typecode)[3]='G')
#define rsb_mm_set_skew(typecode)	((*typecode)[3]='K')
#define rsb_mm_set_hermitian(typecode)((*typecode)[3]='H')

#define rsb_mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
									(*typecode)[2]=' ',(*typecode)[3]='G')

#define rsb_mm_initialize_typecode(typecode) rsb_mm_clear_typecode(typecode)


/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17
#define MM_LIKELY_GZIPPED_FILE		18


/******************** Matrix Market internal definitions ********************

   MM_matrix_typecode: 4-character sequence

				    ojbect 		sparse/   	data        storage 
						  		dense     	type        scheme

   string position:	 [0]        [1]			[2]         [3]

   Matrix typecode:  M(atrix)  C(oord)		R(eal)   	G(eneral)
						        A(array)	C(omplex)   H(ermitian)
											P(attern)   S(ymmetric)
								    		I(nteger)	K(kew)

 ***********************************************************************/

#define MM_MTX_STR		"matrix"
#define MM_ARRAY_STR	"array"
#define MM_DENSE_STR	"array"
#define MM_COORDINATE_STR "coordinate" 
#define MM_SPARSE_STR	"coordinate"
#define MM_COMPLEX_STR	"complex"
#define MM_REAL_STR		"real"
#define MM_INT_STR		"integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR		"symmetric"
#define MM_HERM_STR		"hermitian"
#define MM_SKEW_STR		"skew-symmetric"
#define MM_PATTERN_STR  "pattern"


/*  high level routines */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
int rsb__mm_write_mtx_crd(char fname[], int M, int N, int nz, int IA[], int JA[],
		 double VA[], MM_typecode matcode);
int rsb__mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, rsb_coo_idx_t IA[], rsb_coo_idx_t JA[],
		double VA[], MM_typecode matcode);
int rsb__mm_read_mtx_crd_entry(FILE *f, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, double *real, double *img, MM_typecode matcode);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
#endif /* MM_IO_H_INCLUDED */
/* @endcond */
