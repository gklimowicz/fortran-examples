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
/*!
 @file
 @author Michele Martone

 @brief A toy <rsb.h>-based C program implementing the power
        method for computing matrix eigenvalues. Uses #rsb_spmv().
 \ingroup rsb_doc_examples

 \include power.c
*/

#include <stdio.h>	// printf
#include <math.h>	// sqrt
#include <stdlib.h>	// calloc
#include <rsb.h>

int main(const int argc, char * const argv[])
{
	const int WANT_VERBOSE = 0;
	struct rsb_mtx_t *mtxAp = NULL;
	const int bs = RSB_DEFAULT_BLOCKING;
	const int br = bs, bc = bs; /* bs x bs blocked */
	const rsb_nnz_idx_t nnzA = 4;
	const rsb_coo_idx_t  nrA = 3;
	const rsb_coo_idx_t  ncA = 3;
	rsb_int_t it = 0;
	const rsb_int_t maxit = 100;
	const rsb_coo_idx_t    IA[] = { 0, 1, 2, 0 };
	const rsb_coo_idx_t    JA[] = { 0, 1, 2, 2 };
	const RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE VA[] = { 11, 22, 33, 13 };
	const RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE ZERO = 0;
	int i;
	rsb_err_t errval = 0;

	RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE norm = 0.0, /* nu */
		oldnorm = 1.0, /* oldnorm */
		*b1 = NULL, *b2 = NULL,
		*bnow = NULL, *bnext = NULL;/* b1 and b2 aliases */
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_FIRST_BLAS;
	size_t ds = 0;
       	/* tolerance */
	const RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE tol = 1e-14;

	/* library initialization */
	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS)!=RSB_ERR_NO_ERROR)
		return EXIT_FAILURE;

	/* allocation */
	mtxAp = rsb_mtx_alloc_from_coo_const(VA,IA,JA,nnzA,
			typecode,nrA,ncA,br,bc,RSB_FLAG_NOFLAGS,NULL);
	if(!mtxAp)
		return EXIT_FAILURE;

	ds = (nrA)*sizeof(RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE);
	b1 = calloc(1,ds);
	b2 = calloc(1,ds);

	if(! (b1 && b2))
	{
		errval = RSB_ERR_ENOMEM;
		goto err;
	}

	for( i = 0; i < nrA; ++i )
		b1[i] = 1;

	bnow = b1, bnext = b2;/* b,b' */

	while( fabs(norm-oldnorm) > tol && it<maxit )
	{
		++ it;
		oldnorm = norm;
		/* b'<-Ab */
		if(( rsb_spmv(RSB_TRANSPOSITION_N,NULL,mtxAp,bnow,
			1,&ZERO,bnext,1)) != RSB_ERR_NO_ERROR )
			goto err;
		/* nu<-||Ab||^2 */
		norm = 0;
		for(i=0;i<nrA;++i) 
			norm += bnext[i]*bnext[i];
		/* nu<-||Ab|| */
		norm = sqrt(norm);
		norm = 1.0/norm;
		/* b'<- Ab / ||Ab|| */
		for(i=0;i<nrA;++i)
		       	bnext[i] *= norm;
		norm = 1.0/norm;
		printf("it:%d norm:%lg norm diff:%lg\n",it,norm,norm-oldnorm);

		{void *tmp=bnow;bnow=bnext;bnext=tmp;/* pointers swap */}
		if(WANT_VERBOSE)
		{
			printf("norm:%lg\n",norm);
			if(isinf(norm))
			/* isinf is a C99 feature (need correct
			 * compilation flags) */
				goto err;

			for(i=0;i<2;++i)
				printf("x[%d]=%lg\n",i,((double*)bnext)[i]);
		}
	}
	/* the biggest eigenvalue should be in bnow */

	rsb_mtx_free(mtxAp);
	free(b1);
	free(b2);
	if(rsb_lib_exit(RSB_NULL_EXIT_OPTIONS)!=RSB_ERR_NO_ERROR)
		goto err;
	if( it == maxit )
	{
	       	printf("ERROR: hit iterations limit without convergence!");
	       	errval=RSB_ERR_GENERIC_ERROR;
       	}
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	return EXIT_FAILURE;
}

