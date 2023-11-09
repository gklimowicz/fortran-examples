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
 * This source file contains pixmap rendering functions.
 * */
/*
 * this code is EXPERIMENTAL and UNFINISHED
 *
 * TODO : with very little effort, we could introduce
 *        - a mipmap based rsb_mtx_mipmap_t struct
 * 	  - column major, too
 * 	  - remove extra headers dumpout
 * */

#include "rsb_internals.h"	/* rsb_coo_mtx_t */

rsb_err_t rsb__do_print_postscript_header(FILE*fd, rsb_marf_t rflags, int width, int height, float csw, float csh)
{
	/**
	   \ingroup gr_internals
	   Prints an Encapsulated Postscript header for an area originating at (0 0),
	   bound by (width,height), and defining 'csquare', a square macro..

	   FIXME: this duplicates functionality from rsb_eps.c.
	*/
	const int ex = (RSB_DO_FLAG_HAS(rflags,RSB_MARF_EPS_O)) ? width / RSB_OPS_RENDERING_EXTRA_FRAC : 0; // extra x
#if RSB_ALLOW_STDOUT
	RSB_FPRINTF(fd,
"%%!PS-Adobe-3.0 EPSF-3.0\n"
"%%%%Creator: "RSB_PACKAGE_STRING"\n"
"%%%%Title: matrix plot\n"
"%%%%CreationDate: \n"
"%%%%DocumentData: Clean7Bit\n"
"%%%%Origin: 0 0\n"
"%%%%BoundingBox: %d 0 %d %d\n"
"%%%%LanguageLevel: 2 \n"
"%%%%Pages: 1\n"
"%%%%Page: 1 1\n"
"\n"
"/csquare {\n"
"        newpath\n"
"        0 0 moveto\n"
"        0 -%g rlineto\n"
"        %g 0 rlineto\n"
"        0 %g rlineto\n"
"        closepath\n"
"       setrgbcolor\n"
"        fill\n"
"} def\n"
"\n"
"0 0 moveto\n"
"\n",
(int)(-ex),
(int)(width+ex), (int)height,
csw,csh,csw
);
		RSB_FPRINTF(fd,"save /$LIBRSB_DICT 3 dict def $LIBRSB_DICT begin /M {moveto} bind def /Z {gsave currentpoint lineto %g setlinewidth 1 setlinecap stroke grestore} bind def /D {M Z} bind def /K {0.5 0.5 setrgbcolor} bind def\n",1.0);
		RSB_FPRINTF(fd,"/R {rlineto} bind def\n");
		RSB_FPRINTF(fd,"/N {newpath} bind def\n");
		RSB_FPRINTF(fd,"/L {lineto} bind def\n");
		RSB_FPRINTF(fd,"/C {closepath} bind def\n");
		RSB_FPRINTF(fd,"/SLW {setlinewidth} bind def\n");
		RSB_FPRINTF(fd,"/SRGB {setrgbcolor} bind def\n");
		RSB_FPRINTF(fd,"/SCF {scalefont} bind def\n");
		RSB_FPRINTF(fd,"/SF {setfont} bind def\n");
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

static rsb_err_t rsb_get_pixmap_RGB_from_coo(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, void * pixmap, rsb_coo_idx_t p_rows, rsb_coo_idx_t p_cols, int br, int bc, rsb_flags_t render_flags)
{
	/**
	 * \ingroup gr_internals
	 *
	 * Will fill the specified pixmap (assumed RGB, sized sizeof(char)*p_cols*rows ) with 
	 * non pixels for nonzeros and black pixels for zeros.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	char * dst;
//	rsb_coo_idx_t i,j;
	rsb_nnz_idx_t k;
	size_t sz;
	size_t bpp=3;
	char fgc=0x00,bgc=0x0;

	if(!pixmap || !IA || !JA)
		return RSB_ERR_BADARGS;
	/* should check I and JA and rows and cols and nnz */
//	if(cols>p_cols)
//		return RSB_ERR_BADARGS;
	
	/* DANGER : overflow is possible : FIXME */
	if( RSB_COO_ADD_OVERFLOW(p_rows, p_cols) )
	{
		errval = RSB_ERR_LIMITS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	dst=pixmap;
	sz =p_cols;
	sz*=p_rows;
	sz*=bpp;

	/* 
	 * NOTE : inserting unsorted coefficients is SLOW!
	 * wouldn't be faster to copy and sort ? 
	 * FIXME
	 */
	//RSB_STDOUT("br %d bc %d\n",br,bc);

	if(render_flags & 0x01)
	{
		bgc=~bgc;
	}
	fgc=~bgc;

	if(br==bc && br==1)
	{
		memset(dst,bgc,sz);
		for(k=0;k<nnz;++k)
		{
			dst[bpp*(IA[k]*p_cols+JA[k])+0]=fgc;
			dst[bpp*(IA[k]*p_cols+JA[k])+1]=fgc;
			dst[bpp*(IA[k]*p_cols+JA[k])+2]=fgc;
		}
	}
	else
	{
		memset(dst,bgc,sz);
		for(k=0;k<nnz;++k)
		{
			if(IA[k]/br < 0 || IA[k]/br>p_rows-1 || JA[k]/bc < 0 || JA[k]/bc>p_cols-1)
				RSB_ERROR("I %ld JA %ld %ld %ld\n",(long)IA[k]/br,(long)JA[k]/bc,(long)p_rows,(long)p_cols);
			dst[bpp*((IA[k]/br)*p_cols+JA[k]/bc)+0]=fgc;
			dst[bpp*((IA[k]/br)*p_cols+JA[k]/bc)+1]=fgc;
			dst[bpp*((IA[k]/br)*p_cols+JA[k]/bc)+2]=fgc;
		}
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__do_get_pixmap_RGB_from_matrix(const char * filename, void * pixmap, int width, int height)
{
	/**
	 * \ingroup gr_internals

	 * FIXME : needs error handling

	 * This function is experimentally used to render the sparse matrix;
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t *IA=NULL, *JA=NULL;
	void *VA=NULL;
	rsb_coo_idx_t m=0,k=0;
	rsb_nnz_idx_t nnz=0;
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
#else /* RSB_NUMERICAL_TYPE_DOUBLE */
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_flags_t render_flags=0x01; /*  */
	int br=1,bc=1;
	rsb_time_t t=0;

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE);
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT);

	if(!filename || !pixmap)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	t = - rsb_time();
	if((errval = rsb__util_mm_load_matrix_f(filename,&IA,&JA,&VA,&m,&k,&nnz,typecode,flags,NULL,NULL)))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	t += rsb_time();
//	RSB_STDOUT("load : %lg\n",t);
	
/*	sorting doesn't seem to speed up sparse matrices rendering 	*/
	t = - rsb_time();
//	if((errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,m,k,typecode,flags)))
//	{
//		RSB_PERR_GOTO(err,RSB_ERRM_ES)
//	}
	t += rsb_time();
//	RSB_STDOUT("sort : %lg\n",t);

	if(width <k)
		bc=(k+width-1)/width ;

	if(height<m)
		br=(m+height-1)/height;

	if(br<1 || bc<1)
	{
		// errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	/* rsb__mtx_as_pixmap_resize is optional */
	if( (errval = rsb__mtx_as_pixmap_resize(VA, IA, JA, nnz, &nnz, m, k, height, width, typecode, flags)))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	else
	{
		m=height;
		k=width;
		br=bc=1;
	}

	t = - rsb_time();
	if( (errval = rsb_get_pixmap_RGB_from_coo(IA, JA, nnz, m, k, pixmap, height, width, br, bc, render_flags)) )
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}
	t += rsb_time();
//	RSB_STDOUT("render : %lg\n",t);

err:
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(VA);

	RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb__mtx_as_pixmap_resize(void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_nnz_idx_t *rnnz, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_coo_idx_t p_rows, rsb_coo_idx_t p_cols, rsb_type_t typecode, rsb_flags_t render_flags)
{
	/*
	 * \ingroup gr_internals
	 * User shall be allowed to provide RSB_FLAG_SORTED_INPUT even on unsorted input.
	 * However, in that case only contiguous duplicates will be catched.
	 * FIXME : untested
	 * FIXME : missing comments and error handling.
	 * FIXME : missing input sanitizing
	 */

	double rf,cf;/* row and column factors */
	rsb_nnz_idx_t n;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	rf=((double)p_rows)/(double)m;
	cf=((double)p_cols)/(double)k;

	for(n=0;n<nnz;++n)
	{
		IA[n]= (rsb_coo_idx_t)(((double)IA[n]) * rf);
		JA[n]= (rsb_coo_idx_t)(((double)JA[n]) * cf);
	}
	m=p_rows;
	k=p_cols;
	
	render_flags &= RSB_FLAG_SORTED_INPUT;
	if(!RSB_DO_FLAG_HAS(render_flags,RSB_FLAG_SORTED_INPUT))
	if((errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,m,k,typecode,render_flags/*RSB_FLAG_NOFLAGS*/))!=RSB_ERR_NO_ERROR)
	{
		errval = RSB_ERR_GENERIC_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES)
	}

	*rnnz = rsb__weed_out_duplicates(IA,JA,VA,nnz,typecode,RSB_FLAG_DUPLICATES_SUM/*RSB_FLAG_DUPLICATES_DEFAULT_HANDLE*/|RSB_FLAG_SORTED_INPUT);

	/* missing error handling here */	
err:
	return errval;
}

/* @endcond */
