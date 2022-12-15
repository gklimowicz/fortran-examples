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

 @brief A toy <rsb.h>-based C program showing instantiation,
 transposition and other operations on a single matrix.
 Uses \ref rsb_mtx_clone(), \ref rsb_file_mtx_save(),
  \ref rsb_file_mtx_get_dims(), \ref rsb_file_mtx_load().

 \ingroup rsb_doc_examples

 \include transpose.c
*/
#include <rsb.h>
#include <stdio.h>	/* printf */

int main(const int argc, char * const argv[])
{
	struct rsb_mtx_t *mtxAp = NULL;
	const rsb_blk_idx_t brA = RSB_DEFAULT_BLOCKING,
                            bcA = RSB_DEFAULT_BLOCKING;
	rsb_nnz_idx_t nnzA = 4;
	rsb_coo_idx_t  nrA = 3;
	rsb_coo_idx_t  ncA = 3;
	const rsb_coo_idx_t    IA[] = { 0, 1, 2, 0 };
	const rsb_coo_idx_t    JA[] = { 0, 1, 2, 2 };
	const RSB_DEFAULT_TYPE VA[] = { 11, 22, 33, 13 };
	RSB_DEFAULT_TYPE XV[] = { 0,0,0,0,0,0 };
	rsb_coo_idx_t  vl = 0;
	const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	/* library initialization */
	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS)!=RSB_ERR_NO_ERROR)
	{
		return EXIT_FAILURE;
	}

	/* allocation */
	mtxAp = rsb_mtx_alloc_from_coo_const(
			VA,IA,JA,nnzA,typecode,nrA,ncA,
			brA,bcA,RSB_FLAG_NOFLAGS,NULL);
	if(!mtxAp)
	{
		return EXIT_FAILURE;
	}

	/* printout */
	if(RSB_ERR_NO_ERROR!=(errval = rsb_file_mtx_save(mtxAp,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err;
	}
	
	/* matrix transposition */
	if( RSB_ERR_NO_ERROR != (errval =
		rsb_mtx_clone(&mtxAp,RSB_NUMERICAL_TYPE_SAME_TYPE,
		RSB_TRANSPOSITION_T,NULL,mtxAp,RSB_FLAG_IDENTICAL_FLAGS)))
	{
		goto err;
	}

	/* printout */
	if(RSB_ERR_NO_ERROR!=(errval = rsb_file_mtx_save(mtxAp,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err;
	}

	rsb_mtx_free(mtxAp);

	/* doing the same after load from file */
	mtxAp = rsb_file_mtx_load("pd.mtx",
		RSB_FLAG_NOFLAGS,typecode,NULL);
	if(!mtxAp)
	{
		return EXIT_FAILURE;
	}

	/* printout */
	if(RSB_ERR_NO_ERROR!=(errval = rsb_file_mtx_save(mtxAp,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err;
	}

	/* one can see dimensions in advance, also */
	if(RSB_ERR_NO_ERROR!=(errval =
		rsb_file_mtx_get_dims("pd.mtx",&nrA,&ncA,&nnzA,NULL)))
	{
		if(errval != RSB_ERR_UNSUPPORTED_FEATURE)
			goto err;
	}

	/* A matrix can be rendered to Postscript. */
	{
		if(RSB_ERR_NO_ERROR!=(errval =
		rsb_mtx_rndr("pd.eps",mtxAp,512,512,RSB_MARF_EPS_B)))
			goto err;
	}

	rsb_mtx_free(mtxAp);

	/* also vectors can be loaded */
	if(RSB_ERR_NO_ERROR!=(errval = 
		rsb_file_vec_load("vf.mtx",typecode,NULL,&vl )))
		goto err;
	/* we expect vf.mtx to be 6 rows long */
	if( vl != 6 )
	{
		goto err;
	}

	if(RSB_ERR_NO_ERROR!=(errval = 
		rsb_file_vec_load("vf.mtx",typecode,XV, NULL )))
		goto err;

	/* matrices can be rendered from file to a pixelmap as well */
	{
		unsigned char pixmap[3*2*2];

		if(RSB_ERR_NO_ERROR!=(errval =
		rsb_file_mtx_rndr(pixmap,"pd.mtx",2,2,2,RSB_MARF_RGB)))
			goto err;
	}

	if(RSB_ERR_NO_ERROR != rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))
	{
		goto err;
	}
	return EXIT_SUCCESS;
err:
	rsb_perror(NULL,errval);
	return EXIT_FAILURE;
}

