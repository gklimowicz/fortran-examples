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
 * @brief Matrix info dumping code
 * @author Michele Martone
 *
 * TODO: move other similar functions here.
 * */

#include "rsb_common.h"

#define RSB_CONST_DUMP_DEFAULT_INNER	(RSB_CONST_DUMP_RECURSION | /*RSB_CONST_DUMP_INTERNALS |*/ RSB_CONST_DUMP_TIMES | RSB_CONST_DUMP_DIMENSIONS)

RSB_INTERNALS_COMMON_HEAD_DECLS


static rsb_err_t rsb_do_dump_internals_brief(const struct rsb_mtx_t *mtxAp)
{
	/**
		\ingroup gr_internals
	 */
#if RSB_ALLOW_STDOUT
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//rsb_submatrix_idx_t i,j;
	//struct rsb_mtx_t * submatrix=NULL;
	rsb_submatrix_idx_t smi=0;

	if(!mtxAp)
	{
		return RSB_ERR_BADARGS;
	}

	for(smi=0;smi<mtxAp->all_leaf_matrices_n;++smi)
	{
		RSB_STDOUT_MATRIX_SUMMARY((mtxAp->all_leaf_matrices[smi]).mtxlp);
		RSB_INFO("\n");
	}
	RSB_DO_ERR_RETURN(errval)
#else
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif
}

static rsb_err_t rsb_do_print_matrix_t(const struct rsb_mtx_t *mtxAp, FILE * stream, rsb_dump_flags_t flags)
{
	/**
	 * \ingroup gr_internals
	 * This is a slow debug function to print out a Matrix Market matrix out of the argument mtxAp.
	 * */
#if RSB_ALLOW_STDOUT
	struct rsb_coo_mtx_t coo;
	rsb_flags_t aflags = RSB_FLAG_NOFLAGS;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_bool_t want_calloc = RSB_BOOL_FALSE;
	if(!mtxAp || (!stream && stream != RSB_DEFAULT_STREAM))
	{
		return RSB_ERR_BADARGS;
	}

	RSB_INIT_COO_FROM_MTX(&coo,mtxAp);

	if(flags&RSB_CONST_DUMP_CSR)
		aflags=RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS, want_calloc = RSB_BOOL_TRUE;
	if(rsb__xallocate_coo_matrix_t(&coo,want_calloc,aflags)!=&coo)
       	{
	       	RSB_ERROR(RSB_ERRM_ES); 
		goto err; 
	}
	errval = rsb__do_get_coo(mtxAp,(rsb_byte_t**)(&coo.VA),&coo.IA,&coo.JA,RSB_FLAG_NOFLAGS);
	if(RSB_SOME_ERROR(errval)){RSB_ERROR(RSB_ERRM_ES);goto merr;}
	if(flags&RSB_CONST_DUMP_MATRIX_MARKET || flags&RSB_CONST_DUMP_OCTAVE_STYLE)
		errval = rsb__test_print_coo_mm(mtxAp->typecode,mtxAp->flags,coo.IA,coo.JA,coo.VA,coo.nr,coo.nc,coo.nnz,RSB_BOOL_TRUE,stream);
	if(flags&RSB_CONST_DUMP_CSR)
	{
		errval = rsb__do_switch_fullword_array_to_compressed(coo.IA,coo.nnz,coo.nr);
		if(RSB_SOME_ERROR(errval)){RSB_ERROR(RSB_ERRM_ES);goto merr;}
		errval = rsb__test_print_csr(mtxAp->typecode,mtxAp->flags,coo.IA,coo.JA,coo.VA,coo.nr,coo.nc,coo.nnz,RSB_BOOL_TRUE,stream);
		if(RSB_SOME_ERROR(errval)){RSB_ERROR(RSB_ERRM_ES);goto merr;}
	}
merr:
	rsb__destroy_coo_matrix_t(&coo);
err:
	return errval;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif
}

static rsb_err_t rsb_do_dump_graphviz_dot_graph_do_file_inner(const struct rsb_mtx_t *mtxAp, FILE * fd)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const struct rsb_mtx_t * submatrix = NULL; 
	rsb_submatrix_idx_t i,j;

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}
	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
	{

//#define RSB_FPRINTF_MATRIX_NODE_SUMMARY(FD,M) RSB_FPRINTF_MATRIX_SUMMARY(FD,M)
//#define RSB_FPRINTF_MATRIX_NODE_SUMMARY(FD,M) RSB_FPRINTF(FD,"%dx%d@%d,%d:%d", (M)->nr, (M)->nc, (M)->roff, (M)->coff, (M)->nnz)
#define RSB_FPRINTF_MATRIX_NODE_SUMMARY(FD,M) RSB_FPRINTF(FD,"%zdx%zd\\n@%zd,%zd\\n:%zd(%s)", (rsb_printf_int_t)(M)->nr, (rsb_printf_int_t)(M)->nc, (rsb_printf_int_t)(M)->roff, (rsb_printf_int_t)(M)->coff, (rsb_printf_int_t)(M)->nnz,(rsb__is_recursive_matrix((M)->flags))?("*"):( \
((RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_WANT_COO_STORAGE))? \
(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES)?"HCSR":"CSR"): \
(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES)?"HCOO":"COO")) \
			))

		RSB_FPRINTF(fd,"\"");
		RSB_FPRINTF_MATRIX_NODE_SUMMARY(fd,mtxAp);
		RSB_FPRINTF(fd,"\" -> \"");
		RSB_FPRINTF_MATRIX_NODE_SUMMARY(fd,submatrix);
		RSB_FPRINTF(fd,"\"\n");
		errval = rsb_do_dump_graphviz_dot_graph_do_file_inner(submatrix,fd);
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_dump_graphviz_dot_graph_do_file(const struct rsb_mtx_t *mtxAp, FILE * fd)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}

	RSB_FPRINTF(fd,"/* example usage: dot -Tps filename.dot > filename.ps */\n");
	RSB_FPRINTF(fd,"digraph matrix {\n" 
				"quadtree=TRUE;\n"
			       	"ratio=1.4;\n");
	errval = rsb_do_dump_graphviz_dot_graph_do_file_inner(mtxAp,fd);
	RSB_FPRINTF(fd,"}\n");
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_dump_internals(const struct rsb_mtx_t *mtxAp)
{
	/**
		\ingroup gr_internals
	 */
#if RSB_ALLOW_STDOUT
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;

	if(!mtxAp)
	{
		return RSB_ERR_BADARGS;
	}

	if(rsb__is_root_matrix(mtxAp))
		RSB_STDOUT(	"#R %zd x %zd, %zd nnz (%zd bytes), "
				"%zd index space for bytes, %zd bytes for %zd structs (%zd of which are on the diagonal) "
				"(%3.2lg%% of nnz are on the diagonal) "
				"\n"
			,(size_t)mtxAp->nr
			,(size_t)mtxAp->nc
			,(size_t)mtxAp->nnz
			,((size_t)mtxAp->nnz)*RSB_SIZEOF(mtxAp->typecode)
			,(size_t)rsb__get_index_storage_amount(mtxAp)
			,((size_t)rsb__terminal_recursive_matrix_count(mtxAp))*sizeof(struct rsb_mtx_t)
			,((size_t)rsb__terminal_recursive_matrix_count(mtxAp))
			,(size_t)rsb__get_diagonal_submatrices_count(mtxAp)
			,((((double)rsb__get_diagonal_elements_count(mtxAp))*100)/(mtxAp->nnz))
			);

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		RSB_STDOUT("#T at %zd %zd, %zd x %zd, %zd nnz (%3.2lg%%)\n"
				,(size_t)mtxAp->roff
				,(size_t)mtxAp->coff
				,(size_t)mtxAp->nr
				,(size_t)mtxAp->nc
				,(size_t)mtxAp->nnz
				,((((double)mtxAp->nnz)*100)/(mtxAp->nr))/(mtxAp->nc)
				);
	}
	else
	{
		RSB_STDOUT("#N at %zd %zd, %zd x %zd, %zd nnz (%3.2lg%%)\n"
				,(size_t)mtxAp->roff
				,(size_t)mtxAp->coff
				,(size_t)mtxAp->nr
				,(size_t)mtxAp->nc
				,(size_t)mtxAp->nnz
				,((((double)mtxAp->nnz)*100)/(mtxAp->nr))/(mtxAp->nc)
				);

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_dump_internals(submatrix));
	}
	RSB_DO_ERR_RETURN(errval)
#else
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif
}

rsb_err_t rsb__do_print_matrix_stats(const struct rsb_mtx_t *mtxAp, rsb_dump_flags_t flags, const rsb_char_t*filename)
{
	/**
		\ingroup gr_internals
		Print matrix stats.
		Use stdout is default on NULL filename.
		NOTE: filename is being used only for RSB_CONST_DUMP_MATRIX_MARKET.
		NOTE: not using errno on possible fopen/fclose's failure.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	FILE *stream = NULL;

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}

	if(filename == NULL)
		stream = RSB_DEFAULT_STREAM;
	else
	{
		stream = rsb__util_fopen(filename,"w");
		if(!stream)
		{
			errval = RSB_ERR_GENERIC_ERROR;
			RSB_PERR_GOTO(err,"problems opening %s!\n",filename);
		}
	}

	if(flags == RSB_CONST_DUMP_DEFAULT)
		flags = RSB_CONST_DUMP_DEFAULT_INNER;

	if(flags&RSB_CONST_DUMP_MIF_MTX_INFO)
	{
#if RSB_ALLOW_STDOUT
		RSB_STDOUT(RSB_PRINTF_MTX_SUMMARY_ARGS(mtxAp)),
		RSB_STDOUT("\n");
#else
		errval = RSB_ERR_UNSUPPORTED_FEATURE;
		goto err;
#endif
	}

	if(flags&RSB_CONST_DUMP_FLAGS)
		rsb__dump_flags(mtxAp->flags,""," |\n","\n");

	if(flags&RSB_CONST_DUMP_RECURSION_BRIEF)
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_dump_internals_brief(mtxAp));

	if(flags&RSB_CONST_DUMP_RECURSION)
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_dump_internals(mtxAp));

	if(flags&RSB_CONST_DUMP_DIMENSIONS)
	{
#if RSB_ALLOW_STDOUT
		RSB_STDOUT(
				"m : %zd\n"
				"k : %zd\n"
				"submatrices : %zd\n"
				,(rsb_printf_int_t)(mtxAp->nr)
				,(rsb_printf_int_t)(mtxAp->nc)
				,(rsb_printf_int_t)(mtxAp->all_leaf_matrices_n)
				);
#else
		errval = RSB_ERR_UNSUPPORTED_FEATURE;
		goto err;
#endif
	}

	if(flags&RSB_CONST_DUMP_TIMES)
	{
#if RSB_ALLOW_STDOUT
		RSB_STDOUT(
				"assembly : %10.2lf s\n"
				"perf.est.:%10.2lf s\n"
				"str.anal.:%10.2lf s\n"
				"el.ins.  :%10.2lf s\n"
				"el.sort. :%10.2lf s\n"
				"el.part. :%10.2lf s\n"
				,
				(mtxAp->tat),
				(mtxAp->pet),
				(mtxAp->sat),
				(mtxAp->eit),
				(mtxAp->est),
				(mtxAp->cpt));
#else
		errval = RSB_ERR_UNSUPPORTED_FEATURE;
		goto err;
#endif
	}
	
#if 0
	if(flags&RSB_CONST_DUMP_INTERNALS)
	{
		errval = RSB_ERR_UNSUPPORTED_FEATURE;
	}
#endif

	if(flags&RSB_CONST_DUMP_BLOCKS)
		;

	if(flags&RSB_CONST_DUMP_MATRIX_MARKET || flags&RSB_CONST_DUMP_OCTAVE_STYLE
			|| flags&RSB_CONST_DUMP_CSR )
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_print_matrix_t(mtxAp,stream,flags));
#if 0
	else if(flags&RSB_CONST_DUMP_COO)
		/* May implement Octave/Matlab dump style, here. */
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_print_matrix_t_inner(mtxAp,0,0));
#endif

	if(flags&RSB_CONST_DUMP_RSB)
	{
		rsb_time_t st = - rsb_time();
		RSB_DO_ERROR_CUMULATE(errval,rsb__do_save_matrix_file_as_binary(mtxAp,stream));
		if(RSB_SOME_ERROR(errval)) RSB_PERR_GOTO(err,RSB_ERRM_ES);
		st += rsb_time();
		RSB_IO_NOTICE("#binary saving file %s succeeded and took %lf s (%.0f nnz/s).\n",filename,st,(1.0/(st/mtxAp->nnz)));
	}

	if(flags&RSB_CONST_DUMP_DOT)
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_dump_graphviz_dot_graph_do_file(mtxAp,stream));
err:
	if( stream && stream != RSB_DEFAULT_STREAM )
		fclose(stream);
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_dump_bitmap(const rsb_bitmap_data_t * bmap, size_t w, size_t h)
{
#if RSB_ALLOW_STDOUT
	size_t wi,hi;
	if(!bmap)
		return RSB_ERR_BADARGS;
	for(wi=0;wi<w;++wi)
		for(hi=0;hi<h;++hi)
			RSB_STDOUT("%c",RSB_BITMAP_GET(bmap,w,h,wi,hi)?'1':'0');
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif
}

/* @endcond */
