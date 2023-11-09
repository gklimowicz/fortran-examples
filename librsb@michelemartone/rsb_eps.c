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
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains postscript rendering functions.
 * */

#include "rsb_internals.h"

RSB_INTERNALS_COMMON_HEAD_DECLS

typedef float rsb_rf_t; /* real float */
#define RSB_RR_LEVELS 16

#define RSB_PRINTF_MATRIX_TIME_ARGS(MTXAP)  \
					"ect: %5.2le""  "	\
					"est: %5.2le""  "	\
					"sat: %5.2le""  "	\
					"eit: %5.2le""  "	\
					"cpt: %5.2le"	\
					, \
					(MTXAP)->ect, (MTXAP)->est, (MTXAP)->sat, (MTXAP)->eit, (MTXAP)->cpt

#define RSB_EPS_TRSL(DX,DY,SX,SY,XS,YS,MTXAP) DX = XS*(/*(MTXAP)->nc-*/(SX)); DY = YS*((MTXAP)->nr-(SY)); /* translate for EPS */
#define RSB_MTX_LAR(MTXAP) ((MTXAP)->roff+(MTXAP)->bm) /* last absolute row */
#define RSB_MTX_LAC(MTXAP) ((MTXAP)->coff+(MTXAP)->bk) /* last absolute column */
#define RSB_MTX_LER(MTXAP) ((MTXAP)->broff-(MTXAP)->roff) /* local empty rows */
#define RSB_MTX_LEC(MTXAP) ((MTXAP)->bcoff-(MTXAP)->coff) /* local empty columns */
#define RSB_EPS_NEWPATH "N" /* newpath */
#define RSB_EPS_MOVETO "M" /* moveto */
#define RSB_EPS_LINETO "L" /* lineto */
#define RSB_EPS_RLINETO "R" /* rlineto */
#define RSB_EPS_SCALEFONT "SCF" /* scalefont */
#define RSB_EPS_SETFONT "SF" /* setfont */
#define RSB_EPS_SETRGB "SRGB" /* setrgbcolor */
#define RSB_EPS_LHS_COLOR " 1 0 1 " /* */
#define RSB_EPS_RHS_COLOR " .5 1 .5 " /* */
#define RSB_EPS_SETLINEWIDTH "SLW" /* setlinewidth */
#define RSB_EPS_CLOSEPATH "C" /* closepath */

static rsb_err_t render_ps_box(FILE*fd, int r0, int c0, int dr, int dc, rsb_coo_idx_t orows, rsb_rf_t xs, rsb_rf_t ys, rsb_rf_t r, rsb_rf_t g, rsb_rf_t b)
{
	/**
	   \ingroup gr_internals
	   Prints out a box in the postscript language.

	   \param r0 the box base row
	   \param c0 the box base column
	   \param dr the box height
	   \param dc the box width
	   \param orows the box offset row
	   \param xs the scale on x (rows)
	   \param ys the scale on y (column)
	   \param r red value
	   \param g green value
	   \param b blue value
	 */
#if RSB_ALLOW_STDOUT
		RSB_FPRINTF(fd,
			"newpath\n"
			"%g %g "RSB_EPS_MOVETO"\n"
			"%g %d "RSB_EPS_RLINETO"\n"
			"%d %g "RSB_EPS_RLINETO"\n"
			"%g %d "RSB_EPS_RLINETO"\n"
			"%d %d "RSB_EPS_RLINETO"\n"
			RSB_EPS_CLOSEPATH"\n"
			"%g %g %g "RSB_EPS_SETRGB"\n"
			"1 "RSB_EPS_SETLINEWIDTH"\n"
			"stroke\n\n"
			,
			c0*xs,(orows-r0)*ys,
			 (dc)*xs,  0,
			0,  -(dr)*ys,
			-(dc)*xs, 0,
/*			c0*xs,(orows-r0)*ys,
			 (submatrix->nc)*xs,  0,
			0,  -(submatrix->nr)*ys,
			-(submatrix->nc)*xs, 0,*/
			0, 0,
			r, g, b
			);
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

static void rsb__render_eps_centered_text(FILE*fd, const rsb_rf_t xs, const rsb_rf_t ys, const struct rsb_mtx_t * submatrix, const struct rsb_mtx_t * mtxAp, const char * fstr, rsb_int_t li, rsb_int_t ln)
{
	double fs;
	rsb_rf_t ox,oy;

	RSB_EPS_TRSL(ox,oy,submatrix->coff,submatrix->roff+submatrix->nr/2,xs,ys,mtxAp);
	fs = ( xs * submatrix->nc ) /  ( strlen(fstr) ) * 1.3 ;
	oy -= (li-ln/2)*fs;
	RSB_FPRINTF(fd,"/Courier-Bold findfont %g "RSB_EPS_SCALEFONT" "RSB_EPS_SETFONT" %g %g "RSB_EPS_MOVETO" (%s) 0 0 0 "RSB_EPS_SETRGB" show\n",fs,ox,oy,fstr);
}

static rsb_err_t rsb_dump_postscript_z_curve(FILE*fd, const rsb_marf_t rflags, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_rf_t xs, rsb_rf_t ys, rsb_coo_idx_t orows, int level, rsb_submatrix_idx_t * p, const rsb_submatrix_idx_t * pv, const struct rsb_optrace_t *otv )
{
	/**
	   \ingroup gr_internals
	 */
#if RSB_ALLOW_STDOUT
	rsb_submatrix_idx_t i,j;
	const int levels = RSB_RR_LEVELS;
	const int want_eb=0;/* effective boundaries (FIXME: option currently unfinished) */
	const rsb_bool_t ezl = RSB_DO_FLAG_HAS(rflags,RSB_MARF_EPS_L);

	if(!mtxAp)
	{
		goto err;
	}

	if(level>=levels-1)
		level=levels;

#if 1
	if(pv)
	{
		struct rsb_mtx_t * submatrix=NULL;
		rsb_submatrix_idx_t smi=0;
		RSB_SUBMATRIX_FOREACH_LEAF_PERMUTED(mtxAp,submatrix,smi,pv)
		{
			rsb_rf_t fcoff;
			rsb_rf_t froff;
			rsb_rf_t fnc;
			rsb_rf_t fnr;

			const rsb_rf_t shade= .8 - (.8*smi)/(mtxAp->all_leaf_matrices_n);
			char fstr[RSB_MAX_STRERRLEN];

			if( ezl )
				fcoff=(rsb_rf_t)submatrix->bcoff,
				froff=(rsb_rf_t)submatrix->broff,
				fnc=(rsb_rf_t)RSB_MTX_EFF_C(submatrix),
				fnr=(rsb_rf_t)RSB_MTX_EFF_R(submatrix);
			else
				fcoff=(rsb_rf_t)submatrix->coff,
				froff=(rsb_rf_t)submatrix->roff,
				fnc=(rsb_rf_t)submatrix->nc,
				fnr=(rsb_rf_t)submatrix->nr;

			rsb_dump_postscript_z_curve(fd,rflags,submatrix,submatrix->roff,submatrix->coff,xs,ys,orows,level+1,p,NULL,NULL);

			if(otv)
			{
				const struct rsb_optrace_t ot = otv[pv[smi]];
				const rsb_submatrix_idx_t smc = mtxAp->all_leaf_matrices_n;
				rsb_submatrix_idx_t i, co = 1; /* concurrent operations */

				for(i=smi+1; i < smc && otv[pv[i]].t0 < ot.t1 ;++i)
						co++;
				for(i=smi  ; i > 1 && otv[pv[i-1]].t1 > ot.t0 ;--i)
						co++;
				sprintf(fstr," %d/%d[%d]:%0.1es",1+smi,smc,co,ot.t1-ot.t0);
			}
			else
				sprintf(fstr," %d/%d",1+smi,mtxAp->all_leaf_matrices_n);

			rsb__render_eps_centered_text(fd, xs, ys, submatrix, mtxAp, fstr, 2, 2);
#if 0
			if(want_eb)
			{
				fcoff=(rsb_rf_t)submatrix->bcoff;
				froff=(rsb_rf_t)submatrix->broff;
				fnc=(rsb_rf_t)submatrix->bm;
				fnr=(rsb_rf_t)submatrix->bk;
			}
#endif
			if(smi<mtxAp->all_leaf_matrices_n-1)
				RSB_FPRINTF(fd,"%g %g %g "RSB_EPS_SETRGB" 1 "RSB_EPS_SETLINEWIDTH" stroke "RSB_EPS_NEWPATH" %g %g "RSB_EPS_MOVETO"\n",
					shade,shade,1.0, (fcoff*xs+fnc*(xs/2)), ((-froff+orows))*ys-(fnr)*(ys/2));
		}
		goto ret;
	}
#endif

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		rsb_coo_idx_t scoff;
		rsb_coo_idx_t sroff;
		rsb_coo_idx_t snc;
		rsb_coo_idx_t snr;
		const int want_sc = 1;/* want submatrix comments */

		if( ezl )
			scoff=mtxAp->bcoff,
			sroff=mtxAp->broff,
			snc=RSB_MTX_EFF_C(mtxAp),
			snr=RSB_MTX_EFF_R(mtxAp);
		else
			scoff=coff,
			sroff=roff,
			snc=mtxAp->nc,
			snr=mtxAp->nr;

#if 0
		if(want_eb)
		{
			sroff=roff+RSB_MTX_LER(mtxAp);
			scoff=coff+RSB_MTX_LEC(mtxAp);
			snr=mtxAp->bm;
			snc=mtxAp->bk;
		}
#endif

		if(want_sc)
		RSB_FPRINTF(fd,"%% matrix at %ld %ld, level %ld, xs %g, ys %g, orows %ld\n",(long int)sroff,(long int)scoff,(long int)level,xs,ys,(long int)orows);
		if(*p>0)
			RSB_FPRINTF(fd, "%g %g "RSB_EPS_LINETO"\n" , scoff*xs+snc*(xs/2), ((rsb_rf_t)(orows-sroff))*ys-((rsb_rf_t)snr)*(ys/2));
		else
			RSB_FPRINTF(fd, "%g %g "RSB_EPS_MOVETO"\n" , scoff*xs+snc*(xs/2), ((rsb_rf_t)(orows-sroff))*ys-((rsb_rf_t)snr)*(ys/2));
		++*p;
	}
	else
	{
		struct rsb_mtx_t * submatrix=NULL;
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
//			rsb_coo_idx_t scoff=submatrix->coff;
//			rsb_coo_idx_t sroff=submatrix->roff;
			rsb_coo_idx_t snc=submatrix->nc;
			rsb_coo_idx_t snr=submatrix->nr;

			if(0)
			if(want_eb)
			{
		//		scoff=submatrix->bcoff;
		//		sroff=submatrix->broff;
		//		snr=submatrix->bm;
		//		snc=submatrix->bk;
			}

			rsb_dump_postscript_z_curve(fd,rflags,submatrix, roff+(i?(mtxAp->nr-snr):0), coff+(j?mtxAp->nc-snc:0),xs,ys,orows, level+1,p,NULL,NULL);
		}
	}
ret:
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_GENERIC_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

static rsb_err_t rsb__dump_postscript_ussv_order_curve(const struct rsb_mtx_t * mtxAp, rsb_rf_t xs, rsb_rf_t ys, rsb_submatrix_idx_t * p)
{
	/**
	   \ingroup gr_internals
	   NEW, EXPERIMENTAL
	 */
#if RSB_ALLOW_STDOUT
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_translated_matrix_t * all_leaf_matrices=NULL;
	rsb_submatrix_idx_t all_leaf_matrices_n=0,n;
	FILE*fd = RSB_DEFAULT_FD;

	if(!mtxAp)
	{
		goto err;
	}

	errval = rsb__do_get_submatrices_for_ussv(mtxAp,&all_leaf_matrices,&all_leaf_matrices_n,RSB_TRANSPOSITION_N);
	if(!all_leaf_matrices || RSB_SOME_ERROR(errval))
	{
		errval = RSB_ERR_ENOMEM;
	       	RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(RSB_SOME_ERROR(errval))
	{ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	RSB_FPRINTF(fd,"%%%%there are %d leaves, %d for ussv\n",mtxAp->all_leaf_matrices_n,all_leaf_matrices_n);
#if 1
	for(n=0;n<all_leaf_matrices_n;++n)
	{
		rsb_coo_idx_t rows=all_leaf_matrices[n].nr;
		rsb_coo_idx_t cols=all_leaf_matrices[n].nc;
		rsb_coo_idx_t roff=all_leaf_matrices[n].roff;
		rsb_coo_idx_t coff=all_leaf_matrices[n].coff;
		rsb_coo_idx_t my=mtxAp->nr-((roff+rows/2));
		rsb_coo_idx_t mx=(coff+cols/2);
		rsb_rf_t mys=my;
		rsb_rf_t mxs=mx;
		mys*=ys/mtxAp->nc;
		mxs*=xs/mtxAp->nr;
//		my/=mtxAp->cols;
		RSB_FPRINTF(fd,"%% matrix %ld at %ld %ld, %ld x %ld \n",(long int)n,(long int)roff,(long int)coff,(long int)rows,(long int)cols);
		if(*p>0)
			RSB_FPRINTF(fd, "%g %g "RSB_EPS_LINETO"\n" , mxs, mys);
		else
			RSB_FPRINTF(fd, "%g %g "RSB_EPS_MOVETO"\n" , mxs, mys);
		++*p;
	}
#endif
err:
	RSB_CONDITIONAL_FREE(all_leaf_matrices);
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif /* RSB_ALLOW_STDOUT */
}

int rsb__dump_postscript(const int argc, char * const argv[])
{
	/**
	   \ingroup gr_internals
	*/
	rsb_marf_t rflags = RSB_MARF_NOFLAGS;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_option_t options[] = {
	    {"matrix-filename",	required_argument, NULL, 0x66},/* f */  
	    {"dump-recursion",	no_argument, NULL, 'd'},/* NEW */
	    {"width",required_argument	, NULL, 0x5757},/* W */
	    {"height",required_argument	, NULL, 0x4848}, /* H */
	    {"no-submatrix-format-labels",no_argument	, NULL, 0x6e73666c}, /* --no-submatrix-format-labels */
	    {"latex",no_argument	, NULL, 0x7465780a}, /* --latex */
	    {"auto-size",no_argument	, NULL, 'a'},
	    {"block-dump",no_argument	, NULL, 'B'},
	    {"nonzeros-dump",no_argument, NULL, 'N'},
	    {"block-rowsize",	required_argument, NULL, 0x72 },/* r */
	    {"block-columnsize",	required_argument, NULL, 0x63},/* c */  
	    {"z-dump",	no_argument, NULL, 'z'},
	    {"ussv-dump",	no_argument, NULL, 'S'},
	    RSB_BENCH_PROG_OPTS,
	    {0,0,0,0}
	};

#ifdef RSB_NUMERICAL_TYPE_FLOAT
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_FLOAT;
#else /* RSB_NUMERICAL_TYPE_FLOAT */
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE ;
#else /* RSB_NUMERICAL_TYPE_DOUBLE */
	rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
#endif /* RSB_NUMERICAL_TYPE_DOUBLE */
#endif /* RSB_NUMERICAL_TYPE_FLOAT */
	rsb_blk_idx_t br=1;
	rsb_blk_idx_t bc=1;

	const char * filename=NULL;
	int c,w = RSB_DEFAULT_MATRIX_RENDERING_COLS,h = RSB_DEFAULT_MATRIX_RENDERING_ROWS;
	int opt_index = 0;
	int dump_recursion=0;
	int g_auto_size=0;
	rsb_bool_t want_blocks=0,want_nonzeros=0,z_dump=0;

	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

	if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))!=RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

    	for (;;)
	{
		c = rsb__getopt_long(argc,argv,RSB_SAMPLE_PROGRAM_OPTIONS_GET_FLAGS"ar:c:df:BNzSn:"/*"W:H:"*/,options,&opt_index);
		if (c == -1)
			break;
		RSB_DO_FLAG_ADD(flags,rsb__sample_program_options_get_flags(c,optarg));	/* FIXME : NEW */
		switch (c)
		{
			case 'r':
			br = rsb__util_atoi_km10(optarg);
			if(br<1) { errval = RSB_ERR_BADARGS; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
			break;
			case 'c':
			bc = rsb__util_atoi_km10(optarg);
			if(br<1) { errval = RSB_ERR_BADARGS; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
			break;
			case 'f':
			filename = optarg;
			break;
			case 'N':
			want_nonzeros=1;
			break;
			case 'B':
			want_blocks=1;
			break;
			case 'a':
			g_auto_size=1;
			break;
			case 'd':
			dump_recursion=1;
			break;
			case 'S':
			z_dump=2;
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER_TRIANGULAR);
			break;
			case 'z':
			z_dump=1;
			break;
			case 0x6e73666c:
			rflags = RSB_MARF_EPS_NO_TEXT;
			break;
			case 0x7465780a:
			rflags = RSB_MARF_LATEX_RECURSION;
			break;
 			case 0x4848:
			h = rsb__util_atoi_km10(optarg);
			if(h<1) { errval = RSB_ERR_BADARGS; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
			break;
			case 0x5757:
			w = rsb__util_atoi_km10(optarg);
			if(w<1) { errval = RSB_ERR_BADARGS; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
			break;
			case 'T':
			typecode = toupper(*optarg);
			break;
			case 'n':
#if 1
			rsb__set_num_threads(rsb__util_atoi(optarg));
#else
			{
				rsb_thread_t ca_[1]={1};
				rsb_thread_t * ca=ca_;
				rsb_thread_t cn=1,ci=0;
				ca=NULL;cn=0;
				if(RSB_SOME_ERROR(errval = rsb__util_get_bx_array(optarg,&cn,&ca)))
				{
				       	RSB_PERR_GOTO(err,RSB_ERRM_ES); 
				}
			}
#endif
			default:
			;
	    	}
	}
	
	if(!filename)
	{
		const char*usagestring=" -aRzd -f pd.mtx";
		//errval = RSB_ERR_BADARGS;
		RSB_INFO("Did not specify a matrix file.\n");
		RSB_INFO("Usage example: %s %s\n",argv[0],usagestring);
		//RSB_ERROR(RSB_ERRM_ES);
		goto err;
	}

	RSB_DO_FLAG_DEL(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR); /* FIXME : breaks  -r1 -c1 -Fbr -aRzD */
	if(g_auto_size)
	{
		/* rescale smartly to reflect the desired area and keep proportions (but lose both dimensions) */
		size_t cols,rows;
		rsb_rf_t area=1,p;
		rsb_flags_t flags = RSB_FLAG_NOFLAGS;

		if(RSB_SOME_ERROR(rsb__do_util_get_matrix_dimensions(filename, &cols, &rows, NULL, &flags)))
			goto err;

		area*=w;
		area*=h;
		p=((rsb_rf_t)cols)/((rsb_rf_t)rows);
		h=sqrt(area/p);
		w=h*p;
	}

	if(!dump_recursion)
		want_nonzeros = 1;

	errval = rsb__dump_postscript_recursion_from_matrix(filename,rflags,br,bc,w,h,flags,want_blocks,z_dump,want_nonzeros,dump_recursion,typecode);
err:
	RSB_MASK_OUT_SOME_ERRORS(errval)
	rsb__do_perror(NULL,errval);
	return RSB_ERR_TO_PROGRAM_ERROR(errval);
}

static rsb_err_t rsb__dump_block_rectangles(FILE*fd, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_rf_t xs, rsb_rf_t ys, rsb_coo_idx_t orows, int level)
{
	/**
	 */
#if RSB_ALLOW_STDOUT
	rsb_submatrix_idx_t i,j;
	struct rsb_mtx_t * submatrix = NULL;
	const int levels = RSB_RR_LEVELS;
	const int want_eb = 0;/* want effective boundaries (working) */
	const int want_sc = 1;/* want submatrix comments */
	rsb_rf_t shade;

	if(!mtxAp)
	{
		goto err;
	}

	if(level>=levels-1)
		level=levels;
	shade = 1.0*(level)/levels;

	if(want_sc)RSB_FPRINTF(fd,"%% matrix at %ld %ld, level %ld\n",(long int)roff,(long int)coff,(long int)level);
	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		rsb_coo_idx_t eroff=mtxAp->broff,ecoff=mtxAp->bcoff;
		rsb_coo_idx_t er=RSB_MTX_EFF_R(mtxAp),ec=RSB_MTX_EFF_C(mtxAp);
		if(want_eb==0)
			eroff=mtxAp->roff,ecoff=mtxAp->coff,er=mtxAp->nr,ec=mtxAp->nc;
		if(want_sc)RSB_FPRINTF(fd,"%% terminal matrix at %ld %ld\n",(long int)roff,(long int)coff);
		render_ps_box(fd,eroff, ecoff, er, ec, orows, xs, ys, shade, shade, shade);
	}

	RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
	if(submatrix)
	{
		rsb__dump_block_rectangles(fd,submatrix, roff+(i?(mtxAp->nr-submatrix->nr):0), coff+(j?mtxAp->nc-submatrix->nc:0),xs,ys,orows, level+1);
	}

	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_GENERIC_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

static rsb_err_t rsb__dump_operation_trace(FILE*fd, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, const rsb_submatrix_idx_t *pv, const struct rsb_optrace_t * otv )
{
	/* Serve RSB_MARF_EPS_T */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t i, n = mtxAp->all_leaf_matrices_n;
	rsb_thread_t mth_id = 0;
	rsb_time_t t1 = RSB_TIME_ZERO;

	const char * str = "\
\n\
/T0 { T0A  0   get } bind def %% enter time\n\
/DT { T1 T0 sub } bind def %% delta time\n\
/N2S { 10 12 string cvrs } bind def %% time (number) to string\n\
/SDT { 0 T1A {add} forall  0 T0A {add} forall  sub } bind def %% exit times' sum minus enter times' sum = total execution times\n\
\n\
/X0 { BX 0.2  mul } bind def %% left border coordinate\n\
/X1 { BX 0.95 mul } bind def %% right border coordinate\n\
/AW { X1 X0 sub } bind def %% available width\n\
/SF { AW DT div } bind def %% scale factor\n\
/LC { T0A length } bind def %% lines (per-submatrix) count\n\
/SC { LC TC add } bind def %% total lines (segments) count\n\
/FS { 10 } bind def %% font size\n\
/DY { BY FS sub SC 1 add div } bind def %% DY=BY/(SC+1)\n\
/Y0 { DY 2 mul } bind def\n\
/LCM1 { LC 1 sub } bind def %% LC-1\n\
\n\
0 0 moveto %% bottom lower left\n\
X0 Y0 rmoveto\n\
\n\
/DRWLN %% draw line\n\
{\n\
gsave %% save current graphics state\n\
X0 -1.2 div 0 rmoveto\n\
DY FS gt { /Courier-Bold findfont FS scalefont setfont LLBL 0 0 0 setrgbcolor show } {} ifelse\n\
grestore %% restore current graphics statere\n\
T0A IDX get                 SF mul 0 rmoveto\n\
T1A IDX get T0A IDX get sub SF mul 0 rlineto\n\
} bind def\n\
\n\
/COFF { DY TC 0.5 add mul } bind def\n\
\n\
0 1 LCM1 { \n\
  0 index /IDX exch def %% for /IDX i def\n\
  /XOFF { X0 } bind def %% X offset\n\
  /TOFF { DY TC mul } bind def %% thread offset\n\
  /YOFF { DY 1.5 IDX add mul FS add TOFF add } bind def %% Y offset\n\
  /LLBL { MDA IDX get } bind def %% left label\n\
  XOFF YOFF moveto %% move to start offset\n\
  DRWLN %% this references IDX\n\
\n\
  /TIDX { TIA IDX get } bind def %% thread indxex\n\
  /TOFF { DY TIDX mul } bind def %% thread offset\n\
  /YOFF { DY -.5 mul TOFF add } bind def %% Y offset (fixed + variable)\n\
  /LLBL { TIDX N2S } bind def %% left label\n\
  XOFF YOFF moveto %% move to start offset\n\
  DRWLN %% this references IDX\n\
} for %% for i=0...LCM1\n\
\n\
X0 COFF moveto %% caption coordinates\n\
\n\
/Courier-Bold findfont FS scalefont setfont (Elapsed: ) 0 0 0 setrgbcolor show\n\
/Courier-Bold findfont FS scalefont setfont  DT N2S 0 0 0 setrgbcolor show\n\
/Courier-Bold findfont FS scalefont setfont (s  Total:) 0 0 0 setrgbcolor show\n\
/Courier-Bold findfont FS scalefont setfont SDT N2S 0 0 0 setrgbcolor show\n\
/Courier-Bold findfont FS scalefont setfont (s  Speedup:) 0 0 0 setrgbcolor show\n\
/Courier-Bold findfont FS scalefont setfont SDT DT div N2S 0 0 0 setrgbcolor show\n\
/Courier-Bold findfont FS scalefont setfont (x) 0 0 0 setrgbcolor show\n\
\n\
stroke";
	if (! ( otv && pv ) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	RSB_FPRINTF(fd,"/T0A [");
	for(i=0;i<n;++i)
		RSB_FPRINTF(fd,"%g ",otv[pv[i]].t0);
	RSB_FPRINTF(fd,"] def %% enter times array\n");

	RSB_FPRINTF(fd,"/T1A [");
	for(i=0;i<n;++i)
		RSB_FPRINTF(fd,"%g ",otv[pv[i]].t1);
	RSB_FPRINTF(fd,"] def %% exit times array\n");

	for(i=0;i<n;++i)
		t1 = RSB_MAX(t1,otv[i].t1);
	RSB_FPRINTF(fd,"/T1 { %g } bind def %% exit time\n",t1);

	for(i=0;i<n;++i)
		mth_id = RSB_MAX(mth_id,otv[i].th_id);
	RSB_FPRINTF(fd,"/TC { %d } bind def %% threads count\n",(int)(mth_id+1));

	RSB_FPRINTF(fd,"/TIA [");
	for(i=0;i<n;++i)
		RSB_FPRINTF(fd,"%d ",(int)(otv[pv[i]].th_id+1));
	RSB_FPRINTF(fd,"] def %% thread id array\n");

	RSB_FPRINTF(fd,"/MDA [");
	for(i=0;i<n;++i)
		RSB_FPRINTF(fd,"(M:% 3d Th:% 2d) ",(int)(pv[i]+1),(int)otv[pv[i]].th_id);
	RSB_FPRINTF(fd,"] def %% exit times array\n");

	RSB_FPRINTF(fd,"\n\
%% Fixed part\n\
/BX { %d } bind def %% BoundingBox\n\
/BY { %d } bind def %% BoundingBox\n",(int)br,(int)bc);
	RSB_FPRINTF(fd,"%s\n",str);
err:
	return errval;
}

struct rsb_op_event_t {
	rsb_bool_t ev;
	rsb_submatrix_idx_t smi;
	rsb_time_t t;
};


static const rsb_bool_t Enter = RSB_BOOL_TRUE;

static int rsb__compar_op_event_t(const void * ap, const void * bp)
{
	struct rsb_op_event_t a = *(struct rsb_op_event_t*) ap;
	struct rsb_op_event_t b = *(struct rsb_op_event_t*) bp;

	return a.t > b.t ? 1 :
		(a.t < b.t ? -1 : ( a.ev == Enter ) ); /* avoid timer clash by ordering exit first */
}

rsb_err_t rsb__dump_multiple_recursion_postscript_from_mtx_t(const char * basename, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_marf_t rflags, rsb_bool_t want_blocks, rsb_bool_t z_dump, rsb_bool_t want_nonzeros, const rsb_submatrix_idx_t *pv, const struct rsb_optrace_t * otv )
{
	// TODO: missing unit test
	rsb_submatrix_idx_t i;
	const rsb_submatrix_idx_t n = mtxAp->all_leaf_matrices_n;
	rsb_bitmap_data_t bos[RSB_BYTES_PER_BITVECTOR(n)];
	struct rsb_op_event_t evts[2*n];
	const rsb_bool_t Exit = RSB_BOOL_FALSE;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	char filename[RSB_MAX_FILENAME_LENGTH];

	RSB_BZERO(bos,sizeof(bos));

	if (n<1)
		return errval;

	if ( strlen(basename) > (sizeof(RSB_MAX_FILENAME_LENGTH)-10) ) // for sprintf afterwards
		return errval;

	z_dump = RSB_BOOL_FALSE;

	for(i=0;i<n;++i)
	{
		evts[2*i+0].ev=Enter;
		evts[2*i+0].smi=pv[i];
		evts[2*i+0].t=otv[pv[i]].t0;
		evts[2*i+1].ev=Exit;
		evts[2*i+1].smi=pv[i];
		evts[2*i+1].t=otv[pv[i]].t1;
	}

	qsort( evts, 2*n, sizeof(struct rsb_op_event_t), &rsb__compar_op_event_t);

	for(i=0;i<2*n;++i)
	{
		sprintf(filename,"%s%04d.eps",basename,i);
		if(evts[i].ev == Enter)
			RSB_BITVECTOR_SET(bos,n,evts[i].smi);
		else
			RSB_BITVECTOR_UNSET(bos,n,evts[i].smi);
		errval |= rsb__dump_postscript_recursion_from_mtx_t(NULL,filename,mtxAp,br,bc,width,height,rflags,want_blocks,z_dump,want_nonzeros,pv,otv,bos);
	}
	return errval;
}

static rsb_err_t rsb__dump_postscript_operands(FILE*fd, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_marf_t rflags, rsb_bool_t want_blocks, rsb_bool_t z_dump, rsb_bool_t want_nonzeros, const rsb_submatrix_idx_t *pv, const struct rsb_optrace_t * otv, const rsb_bitmap_data_t * smb)
{
	const rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_submatrix_idx_t smi;
	const rsb_bool_t is_sym = RSB_DO_FLAG_HAS_INTERSECTION(mtxAp->flags,(RSB_FLAG_ANY_SYMMETRY));
	const char * dstr = "\
/NLRM1 { NLR 1 sub } bind def %% \n\
/RN { NLR } bind def\n\
/LL { L0A length 1 sub } bind def\n\
\n\
%% in plot coords:\n\
\
/MY { BY 2 div } bind def %% half (middle) height\n\
/FS { BY 20 div } bind def %% font size\n\
/NRHS { 1 } bind def %% number of right hand sides (you change this)\n\
/SRW { 10 } bind def %% single right hand side width\n\
/HSRW { SRW 2 div } bind def %% half SRW\n\
/OPW { SRW NRHS mul } bind def %% operand width\n\
/XB { 1 } bind def %% X border\n\
/YB { 1 } bind def %% Y border\n\
/LX0 { EX -1.0 mul OPW add XB add } bind def %% LHS leftmost X\n\
/LX1 { LX0 OPW sub } bind def %% LHS\n\
/LXM { LX0 2 div FS 0.5 mul sub } bind def %% LHS operand X\n\
/RX0 { BX EX +1.0 mul add OPW sub } bind def %% RHS\n\
/RX1 { RX0 OPW add XB sub } bind def %% RHS\n\
/RXM { BX RX0 BX sub FS 2 div sub 2 div add } bind def %% RHS operand X\n\
\n\
%% colors:\n\
/LHSSC { 1.0 0.0 1.0 setrgbcolor } bind def\n\
/XHSSC { 0 0 0 setrgbcolor } bind def\n\
/RHSSC { 0.5 1.0 0.5 setrgbcolor } bind def\n\
%% colors for order lines:\n\
/CLDOL { 0.1 0.1 1.0 setrgbcolor } bind def %%       C Layout Dashed Order Line\n\
/CLDOLL{ 0.7 0.7 1.0 setrgbcolor } bind def %%       C Layout Dashed Order Line Light\n\
/FLDOL { CLDOL  } bind def %% Fortran Layout Dashed Order Line\n\
/FLDOLL{ CLDOLL } bind def %% Fortran Layout Dashed Order Line\n\
%%/FLDOL { 0.0 0.6 0.0 setrgbcolor } bind def %% Fortran Layout Dashed Order Line\n\
%%/FLDOLL{ 0.6 0.8 0.6 setrgbcolor } bind def %% Fortran Layout Dashed Order Line\n\
";
	const char * mstr = "\
\n\
/SCATRA { %% matrix to plot Y coords\n\
1 -1 BY NLR div mul scale\n\
0 NLR -1 mul translate\n\
} bind def\n\
\n\
0 0 moveto %% bottom lower left\n\
\n\
0 1 LL {\n\
	0 index /LIDX exch def\n\
	gsave\n\
\n\
	/X0 { LX0 } bind def %% \n\
	/X1 { LX1 } bind def %% \n\
	/Y0 { L0A LIDX get } bind def %% \n\
	/Y1 { L1A LIDX get } bind def %% \n\
\n\
	SCATRA\n\
\n\
	newpath\n\
	X0 Y0 moveto\n\
	X1 Y0 lineto\n\
	X1 Y1 lineto\n\
	X0 Y1 lineto\n\
	closepath\n\
	LHSSC\n\
	1 setlinewidth\n\
	fill\n\
\n\
	/X0 { RX0 } bind def %% \n\
	/X1 { RX1 } bind def %% \n\
	/Y0 { R0A LIDX get } bind def %% \n\
	/Y1 { R1A LIDX get } bind def %% \n\
\n\
	newpath\n\
	X0 Y0 moveto\n\
	X1 Y0 lineto\n\
	X1 Y1 lineto\n\
	X0 Y1 lineto\n\
	closepath\n\
	RHSSC\n\
	1 setlinewidth\n\
	fill\n\
\n\
	grestore\n\
} for\n\
\n\
gsave\n\
SCATRA\n\
\n\
/IOPL { %% Intra-OPerand Lines\n\
	/FOCE { BY SRW div } bind def %% Fictional Operand Columns Elements\n\
	/DI { NLR FOCE div } bind def %% dash interval (matrix row height units)\n\
	/HDI { DI 2 div } bind def %% half DI\n\
\n\
	0 1 eq { %% set equal to enable Fortran dense order\n\
	0 1 NRHS 1 sub {\n\
		0 index /RHSI exch def %% pop RHS index RHSI\n\
		/X1 { TX1 RHSI SRW mul add HSRW add } bind def %% TX1 is this operand.s X1\n\
		/Y0 { HDI } bind def %%\n\
		/Y1 { NLRM1 HDI sub } bind def %%\n\
		FLDOL\n\
		newpath\n\
		X1 Y0 moveto\n\
		X1 Y1 lineto\n\
		1 setlinewidth\n\
		stroke\n\
		%%\n\
		RHSI NRHS 1 sub lt { %% if not last\n\
			[DI DI ] DI setdash\n\
			[] 0 setdash %% solid line\n\
			FLDOLL\n\
			newpath\n\
			X1 Y1 moveto\n\
			X1 SRW add Y0 lineto\n\
			stroke\n\
		} {} ifelse\n\
	%%\n\
	} for\n\
	} {} ifelse \n\
\n\
	0 1 eq { %% set equal to enable C dense order\n\
	0 1 FOCE 1 sub {\n\
		0 index /FVI exch def %% pop fictional vertical index RHSI\n\
		/X0 { TX1 HSRW add } bind def %% TX1 is this operand.s X1\n\
		/X1 { X0 NRHS 1 sub SRW mul add } bind def %%\n\
		/Y0 { HDI FVI DI mul add } bind def %%\n\
\n\
		FVI FOCE 2 sub lt { %% if not last\n\
			CLDOLL\n\
			1 setlinewidth\n\
			%%\n\
			[HSRW HSRW] HSRW 2 div setdash\n\
			newpath\n\
			X1 Y0 moveto\n\
			X0 Y0 DI add lineto\n\
			stroke\n\
		} {} ifelse\n\
		%%\n\
		CLDOL\n\
		[] 0 setdash %% solid line\n\
		newpath\n\
		X0 Y0 moveto\n\
		X1 Y0 lineto\n\
		stroke\n\
	} for\n\
	} {} ifelse \n\
\n\
	%% NRHS separator lines\n\
	1 1 NRHS 1 sub {\n\
		0 index /RHSI exch def %% pop RHS index RHSI\n\
		/X1 { TX1 RHSI SRW mul add } bind def %% TX1 is this operand's X1\n\
		/Y0 { 0 } bind def %%\n\
		/Y1 { NLRM1 } bind def %%\n\
		newpath\n\
		X1 Y0 moveto\n\
		X1 Y1 lineto\n\
		XHSSC\n\
		1 setlinewidth\n\
		stroke\n\
		[] 0 setdash %% solid line\n\
	} for\n\
} bind def %% Intra-OPerand Lines\n\
%%\n\
/TX1 { LX1 } bind def %% this X1\n\
IOPL\n\
/TX1 { RX0 } bind def %% this X1\n\
IOPL\n\
\n\
\n\
/X0 { LX0 } bind def %% \n\
/X1 { LX1 } bind def %% \n\
/Y0 { 0 } bind def %% \n\
/Y1 { NLRM1 } bind def %% \n\
newpath\n\
X0 Y0 moveto\n\
X1 Y0 lineto\n\
X1 Y1 lineto\n\
X0 Y1 lineto\n\
closepath\n\
XHSSC\n\
1 setlinewidth\n\
stroke\n\
\n\
/X0 { 0 } bind def %% \n\
/X1 { BX } bind def %% \n\
/Y0 { 0 } bind def %% \n\
/Y1 { NLRM1 } bind def %% \n\
newpath\n\
X0 Y0 moveto\n\
X1 Y0 lineto\n\
X1 Y1 lineto\n\
X0 Y1 lineto\n\
closepath\n\
XHSSC\n\
1 setlinewidth\n\
stroke\n\
\n\
/X0 { RX0 } bind def %% \n\
/X1 { RX1 } bind def %% \n\
/Y0 { 0 } bind def %% \n\
/Y1 { NLRM1 } bind def %% \n\
newpath\n\
X0 Y0 moveto\n\
X1 Y0 lineto\n\
X1 Y1 lineto\n\
X0 Y1 lineto\n\
closepath\n\
XHSSC\n\
1 setlinewidth\n\
stroke\n\
grestore\n\
\n\
LXM MY moveto\n\
/Courier-Bold findfont FS scalefont setfont (+=) 0 0 0 setrgbcolor show\n\
RXM MY moveto\n\
/Courier-Bold findfont FS scalefont setfont (*) 0 0 0 setrgbcolor show\n\
";

	RSB_FPRINTF(fd,"%% in plot coords:\n");
	RSB_FPRINTF(fd,"/BX { %d } bind def %% Matrix X width\n",(int)width);
	RSB_FPRINTF(fd,"/EX { %d } bind def %% Extra X on each side (negative on left, positive on right)\n",(int)(width/RSB_OPS_RENDERING_EXTRA_FRAC));
	RSB_FPRINTF(fd,"/BY { %d } bind def %% Matrix X height\n",(int)height);
	RSB_FPRINTF(fd,"\n");

	if(is_sym)
		RSB_FPRINTF(fd,"%% Note: matrix considered as symmetric.\n");
	RSB_FPRINTF(fd,"%% in matrix coords:\n");
	RSB_FPRINTF(fd,"/L0A [");
	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
	if(!smb || (RSB_BITVECTOR_GET(smb,mtxAp->all_leaf_matrices_n,smi)))
	{
		RSB_FPRINTF(fd,"%d ",(int)(submatrix->broff));
		if(is_sym)
			RSB_FPRINTF(fd,"%d ",(int)(submatrix->bcoff));
	}
	RSB_FPRINTF(fd,"] def %% LHS \n");

	RSB_FPRINTF(fd,"/L1A [");
	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
	if(!smb || (RSB_BITVECTOR_GET(smb,mtxAp->all_leaf_matrices_n,smi)))
	{
		RSB_FPRINTF(fd,"%d ",(int)(submatrix->broff+RSB_MTX_EFF_R(submatrix)));
		if(is_sym)
			RSB_FPRINTF(fd,"%d ",(int)(submatrix->bcoff+RSB_MTX_EFF_C(submatrix)));
	}
	RSB_FPRINTF(fd,"] def %% LHS \n");

	RSB_FPRINTF(fd,"/R0A [");
	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
	if(!smb || (RSB_BITVECTOR_GET(smb,mtxAp->all_leaf_matrices_n,smi)))
	{
		RSB_FPRINTF(fd,"%d ",(int)(submatrix->bcoff));
		if(is_sym)
			RSB_FPRINTF(fd,"%d ",(int)(submatrix->broff));
	}
	RSB_FPRINTF(fd,"] def %% RHS \n");

	RSB_FPRINTF(fd,"/R1A [");
	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
	if(!smb || (RSB_BITVECTOR_GET(smb,mtxAp->all_leaf_matrices_n,smi)))
	{
		RSB_FPRINTF(fd,"%d ",(int)(submatrix->bcoff+RSB_MTX_EFF_C(submatrix)));
		if(is_sym)
			RSB_FPRINTF(fd,"%d ",(int)(submatrix->broff+RSB_MTX_EFF_R(submatrix)));
	}
	RSB_FPRINTF(fd,"] def %% RHS \n");

	RSB_FPRINTF(fd,"/NLR {%d} bind def %% LHS,RHS\n",(int)(RSB_MAX(mtxAp->nr,mtxAp->nc)));
	RSB_FPRINTF(fd,"%s%s\n",dstr,mstr);

	return errval;
}

static rsb_err_t rsb__dump_postscript_recursion_blocks(FILE*fd, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_marf_t rflags, rsb_bool_t want_blocks, rsb_bool_t z_dump, rsb_bool_t want_nonzeros, const rsb_submatrix_idx_t *pv, const struct rsb_optrace_t * otv, const rsb_bitmap_data_t * smb)
{
	const rsb_coo_idx_t nrA = mtxAp->nr, ncA = mtxAp->nc;
	const rsb_rf_t ys = ((rsb_rf_t)height)/nrA, xs = ((rsb_rf_t)width)/ncA;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * submatrix = NULL;
	rsb_submatrix_idx_t smi;
	double mnnz = 0, annz = 0;
	const rsb_coo_idx_t hoo = -RSB_MIN(RSB_MIN(mtxAp->nr,mtxAp->nc)/1000,10); /* how much out; this shall turn some coloured lines inwards the box; on smaller matrices this shall be limited */
	const rsb_bool_t want_cs = !RSB_DO_FLAG_HAS(rflags,RSB_MARF_EPS_O);
	const rsb_bool_t want_bb = RSB_BOOL_TRUE;

	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
		mnnz = RSB_MAX(mnnz,submatrix->nnz),
		annz += submatrix->nnz;
	annz /= mtxAp->all_leaf_matrices_n;

	RSB_FPRINTF(fd,"%% colored boxes dump\n");
	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
	if(!smb || (RSB_BITVECTOR_GET(smb,mtxAp->all_leaf_matrices_n,smi)))
	{
		rsb_rf_t fx,fy;
		double rv,gv,bv, iv;
		double r0,dr,c0,dc;

		rv = gv = bv = 0.7;

		if( want_bb )
		{
			r0 = submatrix->broff;
			c0 = submatrix->bcoff;
			dr = RSB_MTX_EFF_R(submatrix);
			dc = RSB_MTX_EFF_C(submatrix);
		}
		else
		{
			r0 = submatrix->roff;
			c0 = submatrix->coff;
			dr = submatrix->nr;
			dc = submatrix->nc;
		}

		RSB_EPS_TRSL(fx,fy,c0,r0,xs,ys,mtxAp);
		RSB_FPRINTF(fd,RSB_EPS_NEWPATH" %g %g "RSB_EPS_MOVETO" ",fx,fy);
		RSB_FPRINTF(fd,"%g %g "RSB_EPS_RLINETO" ", xs*dc,0.0);
		RSB_FPRINTF(fd,"%g %g "RSB_EPS_RLINETO" ",0.0,-ys*dr);
		RSB_FPRINTF(fd,"%g %g "RSB_EPS_RLINETO" ",-xs*dc,0.0);

		if( want_cs )
		{

			if(submatrix->nnz > annz)
				iv = 0.3 * ( ( - annz + submatrix->nnz ) / submatrix->nnz ),
				rv+=iv;
			else
				iv = 0.3 * ( ( + annz - submatrix->nnz ) / annz ),
				gv+=iv;
		}

		RSB_FPRINTF(fd,RSB_EPS_CLOSEPATH" %g %g %g "RSB_EPS_SETRGB" 1 "RSB_EPS_SETLINEWIDTH" fill ",rv,gv,bv);
		RSB_FPRINTF(fd,"%% submatrix %d square\n",smi);

		if(!want_bb)
		{
			rsb_rf_t fx,fy,tx,ty;

			RSB_EPS_TRSL(fx,fy,-hoo+c0,-hoo+r0      ,xs,ys,mtxAp);
			RSB_EPS_TRSL(tx,ty,-hoo+c0,-hoo+r0+dr,xs,ys,mtxAp);
			RSB_FPRINTF(fd,RSB_EPS_NEWPATH" %g %g "RSB_EPS_MOVETO" %g %g "RSB_EPS_LINETO" ",fx,fy,tx,ty);
			RSB_FPRINTF(fd,RSB_EPS_CLOSEPATH RSB_EPS_LHS_COLOR RSB_EPS_SETRGB" 1 "RSB_EPS_SETLINEWIDTH" 1 stroke ");
			RSB_FPRINTF(fd,"%% submatrix %d lhs\n",smi);

			RSB_EPS_TRSL(fx,fy,-hoo+c0,             -hoo+submatrix->broff,xs,ys,mtxAp);
			RSB_EPS_TRSL(tx,ty,-hoo+c0+dc,-hoo+r0,xs,ys,mtxAp);
			RSB_FPRINTF(fd,RSB_EPS_NEWPATH" %g %g "RSB_EPS_MOVETO" %g %g "RSB_EPS_LINETO" ",fx,fy,tx,ty);
			RSB_FPRINTF(fd,RSB_EPS_CLOSEPATH RSB_EPS_RHS_COLOR RSB_EPS_SETRGB" 1 "RSB_EPS_SETLINEWIDTH" 1 stroke ");
			RSB_FPRINTF(fd,"%% submatrix %d rhs\n",smi);
		}
	}

	RSB_FPRINTF(fd,"%% lhs dump\n");
	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
	if(want_bb)
	if(!smb || (RSB_BITVECTOR_GET(smb,mtxAp->all_leaf_matrices_n,smi)))
	{
		rsb_rf_t fx,fy,tx,ty;
	/*
		RSB_FPRINTF(fd,RSB_EPS_NEWPATH" ");
		RSB_EPS_TRSL(fx,fy,submatrix->coff ,submatrix->roff,xs,ys,mtxAp);
		RSB_EPS_TRSL(tx,ty,0               ,submatrix->roff,xs,ys,mtxAp);
		RSB_FPRINTF(fd,"%g %g "RSB_EPS_MOVETO" ",fx,fy);
		RSB_FPRINTF(fd,"%g %g lineto ",tx,ty);
		px=tx,py=ty;
		RSB_EPS_TRSL(fx,fy,submatrix->coff ,submatrix->roff+submatrix->nr,xs,ys,mtxAp);
		RSB_EPS_TRSL(tx,ty,0               ,submatrix->roff+submatrix->nr,xs,ys,mtxAp);
		RSB_FPRINTF(fd,"%g %g "RSB_EPS_MOVETO" ",fx,fy);
		RSB_FPRINTF(fd,"%g %g lineto ",tx,ty);

		RSB_FPRINTF(fd,"closepath "); RSB_FPRINTF(fd,"1 0 1 "RSB_EPS_SETRGB" "); RSB_FPRINTF(fd,"1 setlinewidth "); RSB_FPRINTF(fd,"1 stroke ");
		RSB_FPRINTF(fd,"%% submatrix %d to-lhs\n",smi);

		RSB_FPRINTF(fd,RSB_EPS_NEWPATH" ");
		RSB_FPRINTF(fd,"%g %g "RSB_EPS_MOVETO" ",px,py);
		RSB_FPRINTF(fd,"%g %g lineto ",tx,ty);
		RSB_FPRINTF(fd,"closepath "); RSB_FPRINTF(fd,"1 0 1 "RSB_EPS_SETRGB" "); RSB_FPRINTF(fd,"5 setlinewidth "); RSB_FPRINTF(fd,"1 stroke ");
		RSB_FPRINTF(fd,"%% submatrix %d lhs\n",smi);
		*/

		RSB_EPS_TRSL(fx,fy,-hoo+submatrix->bcoff,-hoo+submatrix->broff      ,xs,ys,mtxAp);
		RSB_EPS_TRSL(tx,ty,-hoo+submatrix->bcoff,-hoo+RSB_MTX_LAR(submatrix),xs,ys,mtxAp);
		RSB_FPRINTF(fd,RSB_EPS_NEWPATH" %g %g "RSB_EPS_MOVETO" %g %g "RSB_EPS_LINETO" ",fx,fy,tx,ty);
		RSB_FPRINTF(fd,RSB_EPS_CLOSEPATH RSB_EPS_LHS_COLOR RSB_EPS_SETRGB" 1 "RSB_EPS_SETLINEWIDTH" 1 stroke ");
		RSB_FPRINTF(fd,"%% submatrix %d lhs\n",smi);
	}

	RSB_FPRINTF(fd,"%% rhs dump\n");
	RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
	if(want_bb)
	if(!smb || (RSB_BITVECTOR_GET(smb,mtxAp->all_leaf_matrices_n,smi)))
	{
		rsb_rf_t fx,fy,tx,ty;
		//rsb_rf_t ih = (RSB_DO_FLAG_HAS(submatrix->flags,RSB_FLAG_USE_HALFWORD_INDICES)) ? 1.0 : 0.0;
	/*
		RSB_FPRINTF(fd,"%% submatrix %d\n",smi);
		RSB_FPRINTF(fd,RSB_EPS_NEWPATH" ");

		RSB_EPS_TRSL(fx,fy,submatrix->coff ,submatrix->roff,xs,ys,mtxAp);
		RSB_EPS_TRSL(tx,ty,mtxAp->nr       ,submatrix->coff,xs,ys,mtxAp);
		RSB_FPRINTF(fd,"%g %g moveto ",fx,fy);
		RSB_FPRINTF(fd,"%g %g lineto ",tx,ty);
		px=tx,py=ty;
		RSB_EPS_TRSL(fx,fy,submatrix->coff+submatrix->nc,submatrix->roff,              xs,ys,mtxAp);
		RSB_EPS_TRSL(tx,ty,mtxAp->nr                    ,submatrix->coff+submatrix->nc,xs,ys,mtxAp);
		RSB_FPRINTF(fd,"%g %g moveto ",fx,fy);
		RSB_FPRINTF(fd,"%g %g lineto ",tx,ty);
		RSB_FPRINTF(fd,"closepath "); RSB_FPRINTF(fd,"0.5 1 0.5 setrgbcolor "); RSB_FPRINTF(fd,"1 setlinewidth "); RSB_FPRINTF(fd,"1 stroke ");
		RSB_FPRINTF(fd,"%% submatrix %d to-rhs\n",smi);

		RSB_FPRINTF(fd,RSB_EPS_NEWPATH" ");
		RSB_FPRINTF(fd,"%g %g moveto ",px,py);
		RSB_FPRINTF(fd,"%g %g lineto ",tx,ty);
		RSB_FPRINTF(fd,"closepath "); RSB_FPRINTF(fd,"0.5 1 0.5 setrgbcolor "); RSB_FPRINTF(fd,"5 setlinewidth "); RSB_FPRINTF(fd,"1 stroke ");
		RSB_FPRINTF(fd,"%% submatrix %d rhs\n",smi);
	*/
		RSB_EPS_TRSL(fx,fy,-hoo+submatrix->bcoff,             -hoo+submatrix->broff,xs,ys,mtxAp);
		RSB_EPS_TRSL(tx,ty,-hoo+submatrix->coff+submatrix->bk,-hoo+submatrix->broff,xs,ys,mtxAp);
		RSB_FPRINTF(fd,RSB_EPS_NEWPATH" %g %g "RSB_EPS_MOVETO" %g %g "RSB_EPS_LINETO" ",fx,fy,tx,ty);
		RSB_FPRINTF(fd,RSB_EPS_CLOSEPATH RSB_EPS_RHS_COLOR RSB_EPS_SETRGB" 1 "RSB_EPS_SETLINEWIDTH" 1 stroke ");
		RSB_FPRINTF(fd,"%% submatrix %d rhs\n",smi);
	}

	if( !RSB_DO_FLAG_HAS(rflags,RSB_MARF_EPS_NO_TEXT) )
	if( !RSB_DO_FLAG_HAS(rflags,RSB_MARF_EPS_O) )
	{
		RSB_FPRINTF(fd,"%% node content labels\n");
		RSB_SUBMATRIX_FOREACH_LEAF(mtxAp,submatrix,smi)
		{
			char fstr[RSB_MAX_STRERRLEN];
			sprintf(fstr," %d/%d %s%s %0.1e",(int)(1+smi),(int)(mtxAp->all_leaf_matrices_n),(RSB_DO_FLAG_HAS(submatrix->flags,RSB_FLAG_USE_HALFWORD_INDICES))?"H":"",(submatrix->matrix_storage == RSB_MATRIX_STORAGE_BCOR)?"COO":"CSR",(double)(submatrix->nnz));
			rsb__render_eps_centered_text(fd, xs, ys, submatrix, mtxAp, fstr, 1, 2);
		}
	}
	return errval;
}

static rsb_err_t rsb__dump_latex_matrix_recursion_internal(FILE*fd, const struct rsb_mtx_t*mtxAp, rsb_submatrix_idx_t lvl, rsb_submatrix_idx_t * smip)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(rsb__is_terminal_recursive_matrix(mtxAp))
	{
		const int smi = (int)(++*smip);

		switch(3)
		{
			case(0):
			RSB_FPRINTF(fd,"M_{%d}^{%s%s}",smi,(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))?"H":"",(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR)?"COO":"CSR");
			break;
			case(1):
			RSB_FPRINTF(fd,"\\%s{M}_{%d}^{%s}",(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR)?"dot":"overline",smi,(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))?"\\prime":"");
			break;
			case(2):
			RSB_FPRINTF(fd,"\\%s{%.1le}_{%d}^{%s}",(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR)?"dot":"overline",(double)mtxAp->nnz,smi,(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))?"\\prime":"");
			case(3):
			RSB_FPRINTF(fd,"\\%s{\\left[\\frac{%.1le nz}{%.1le pr}\\right]}_{}^{%s}",(mtxAp->matrix_storage == RSB_MATRIX_STORAGE_BCOR)?"dot":"overline",(double)mtxAp->nnz,((double)mtxAp->nnz)/RSB_MTX_EFF_R(mtxAp),(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_USE_HALFWORD_INDICES))?"\\prime":"");
			break;
		}
	}
	else
	{
		const char * ca[] = {"&","\\\\","&","\\\\"};
		const char mc = 'p'; // p: round  b: compact  v: lines
		rsb_submatrix_idx_t i;

		++lvl;
		RSB_FPRINTF(fd,"\\begin{%cmatrix}\n",mc);

		for(i=0;i<4;++i)
		{
			if(mtxAp->sm[i])
				rsb__dump_latex_matrix_recursion_internal(fd, mtxAp->sm[i], lvl, smip);
			else
				RSB_FPRINTF(fd,"0");
			RSB_FPRINTF(fd,"%s",ca[i]);
		}

		RSB_FPRINTF(fd,"\\end{%cmatrix}\n",mc);
	}
	return errval;
}

static rsb_err_t rsb__dump_latex_matrix_recursion(FILE*fd, const struct rsb_mtx_t*mtxAp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const rsb_submatrix_idx_t lvl = 0;
	rsb_submatrix_idx_t smi = 0;

	RSB_FPRINTF(fd,"\\documentclass[9pt]{article}\n");
	RSB_FPRINTF(fd,"\\usepackage[paperwidth=40cm, paperheight=20cm, top=0cm, bottom=0cm, outer=0cm, inner=0cm]{geometry}\n");
	RSB_FPRINTF(fd,"\\usepackage{amsmath}\n");
	RSB_FPRINTF(fd,"\\begin{document}\n");
	RSB_FPRINTF(fd,"\\begin{tiny}\n");
	RSB_FPRINTF(fd,"\\begin{math}\\begin{aligned}\n");
	errval = rsb__dump_latex_matrix_recursion_internal(fd, mtxAp, lvl, &smi);
	RSB_FPRINTF(fd,"\\end{aligned}\\end{math}\n");
	RSB_FPRINTF(fd,"\\end{tiny}\n");
	RSB_FPRINTF(fd,"\\end{document}\n");

	return errval;
}

rsb_err_t rsb__dump_postscript_recursion_from_mtx_t(FILE*fd, const char * filename, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_marf_t rflags, rsb_bool_t want_blocks, rsb_bool_t z_dump, rsb_bool_t want_nonzeros, const rsb_submatrix_idx_t *pv, const struct rsb_optrace_t * otv, const rsb_bitmap_data_t * smb)
{
	/*
	 * ( rflags == RSB_FLAG_NOFLAGS ) is allowed and implies defaults.
	 * */
	const rsb_coo_idx_t nrA = mtxAp->nr, ncA = mtxAp->nc;
	const int want_structure_comments_dump = 1;
	const rsb_rf_t ys = ((rsb_rf_t)height)/nrA, xs = ((rsb_rf_t)width)/ncA;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( (fd && filename) || ( !fd && !filename) )
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	if(filename)
	{
		fd = rsb__util_fopen(filename,"w");
		if( fd == NULL )
		{
			errval=RSB_ERR_GENERIC_ERROR;
			goto err;
		}
	}

	errval = rsb__do_print_postscript_header(fd, rflags, width, height, xs, ys );

	if( RSB_DO_FLAG_HAS(rflags,RSB_MARF_EPS_O) )
		errval = rsb__dump_postscript_operands(fd, mtxAp, br, bc, width, height, rflags, want_blocks, z_dump, want_nonzeros, pv, otv, smb);

	if ( rflags == RSB_MARF_EPS_T )
	{
		errval = rsb__dump_operation_trace(fd, mtxAp, br, bc, pv, otv);
		goto done;
	}

	if(want_structure_comments_dump)
	{
		RSB_POSTSCRIPT_DUMP_COMMENT(fd,RSB_PRINTF_MTX_SUMMARY_ARGS(mtxAp));
		RSB_POSTSCRIPT_DUMP_COMMENT(fd,RSB_PRINTF_MATRIX_TIME_ARGS(mtxAp));
		RSB_FPRINTF(fd,"%%%% ");
		rsb__fprint_matrix_implementation_code(mtxAp, "", mtxAp->flags, fd);
		RSB_FPRINTF(fd,"\n");
	}

	RSB_POSTSCRIPT_DUMP_COMMENT(fd,"sparse blocks dump");
	errval = rsb__dump_block_rectangles(fd,mtxAp,0,0,xs,ys,mtxAp->nr,0);

	if( RSB_DO_FLAG_HAS(rflags,RSB_MARF_EPS_L) )
		errval = rsb__dump_postscript_recursion_blocks(fd, mtxAp, br, bc, width, height, rflags, want_blocks, z_dump, want_nonzeros, pv, otv, smb);

	if(z_dump)
	{
		rsb_submatrix_idx_t p = 0;
		RSB_POSTSCRIPT_DUMP_COMMENT(fd,"z dump\nnewpath");

		if(z_dump==1)
			errval = rsb_dump_postscript_z_curve(fd,rflags,mtxAp, 0,0,xs,ys,mtxAp->nr,0,&p,pv,otv);
		else
			errval = rsb__dump_postscript_ussv_order_curve(mtxAp,(rsb_rf_t)height,(rsb_rf_t)width,&p);
		RSB_FPRINTF(fd,
			"%d %d  %d setrgbcolor\n"
			"1 setlinewidth\n"
			"stroke\n\n"
			,
			0,0,1
			);
	}
	if(want_blocks)
		;/* dead code removed from here */
done:
	if(filename)
	{
		fclose(fd);
	}
err:
	return errval;
}

static rsb_err_t rsb_dump_postscript_from_coo(FILE*fd, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, void *VA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_bool_t all_nnz, rsb_type_t typecode)
{
	/**
	 \ingroup gr_internals
	 Need better error handling.
	 This function is experimentally used to render the sparse matrix.
	 */
#if RSB_ALLOW_STDOUT
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t n=0;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_rf_t ys,xs;
	rsb_rf_t csh,csw,dx=0.0,dy=0.0,rd=1.0;

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE) ;
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT) ;

if(1)
{
	       	rsb_coo_idx_t /* ri=0,ci=0,*/nr=m,nc=k;
	       	rsb_nnz_idx_t nzi=0;
	       	rsb_nnz_idx_t onnz=nnz;
		// RSB_STDERR("%s","FIXME: this is functioning code to render PostScript spy plots; it just needs to be called the right way ..\n");
		rsb_aligned_t max[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
		rsb_aligned_t min[RSB_CONST_ENOUGH_ALIGNED_FOR_ANY_TYPE];
		const int want_points_plot = 1;

		rsb__fill_with_ones(VA,typecode,nnz,1);
#if 1
		while( (nnz>100000 || nr>RSB_DEFAULT_MATRIX_RENDERING_ROWS || nc>RSB_DEFAULT_MATRIX_RENDERING_COLS ) && (nr>2 && nc>2))
		{
			/*
				May be better to write a function *resize_to_nnz*.
				This code is quite poor but does the job.
			 */
		       	rsb_coo_idx_t nnr=nr/2, nnc=nc/2;
			rsb_flags_t flags = RSB_FLAG_NOFLAGS;
			// RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT) ;
			// RSB_STDERR("will rescale %d %d (%d nz) to %d %d...\n",nr,nc,nnz,nnr,nnc);
			errval = rsb__mtx_as_pixmap_resize(VA,IA,JA,nnz,&nnz,nr,nc,nnr,nnc,typecode,flags);
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT) ;
			// nnz = rsb__weed_out_duplicates(IA,JA,VA,nnz,typecode,RSB_FLAG_DUPLICATES_SUM/*|RSB_FLAG_SORTED_INPUT*/);
			if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
			nc=nnc; nr=nnr;
			nnc/=2; nnr/=2;
		}
#endif
#if 1
		// if( (nnz>100000 || nr>height || nc>width ) && (nr>2 && nc>2))
		if(1)
		{
		       	rsb_coo_idx_t nnr=height, nnc=width;
			rsb_flags_t flags = RSB_FLAG_NOFLAGS;
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT) ;
			// RSB_STDERR("will rescale further %d %d (%d nz) to %d %d...\n",nr,nc,nnz,nnr,nnc);
			if(nnr>nr)
				rd=((rsb_rf_t)nnr)/((rsb_rf_t)nr);
			errval = rsb__mtx_as_pixmap_resize(VA,IA,JA,nnz,&nnz,nr,nc,nnr,nnc,typecode,flags);
			if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
			nc=nnc; nr=nnr;
		}
#endif
		errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,nr,nc, typecode,RSB_FLAG_NOFLAGS);
		if(RSB_SOME_ERROR(errval)) 
		{
		       	RSB_ERROR(RSB_ERRM_ES); /* RSB_PERR_GOTO(err,RSB_ERRM_ES)*/ /*not a critical error; however shall emit a warning : */
		}
		nnz = rsb__weed_out_duplicates(IA,JA,VA,nnz,typecode,RSB_FLAG_DUPLICATES_SUM|RSB_FLAG_SORTED_INPUT);

		RSB_FPRINTF(fd,""
			"%%!PS-Adobe-3.0 EPSF-3.0\n"
			"%%%%Creator: "RSB_PACKAGE_STRING"\n"
			"%%%%Title: matrix spy plot (originally %ld x %ld / %ld, here %ld x %ld/%ld)\n"
			"%%%%CreationDate: \n"
			"%%%%DocumentData: Clean7Bit\n"
			"%%%%Origin: 0 0\n"
			"%%%%BoundingBox: 0 0 %ld %ld\n"
			"%%%%LanguageLevel: 2\n"
			"%%%%Pages: 1\n"
			"%%%%Page: 1 1\n"
			/* "0.5 0.5 0.5 setrgbcolor\n" */
			,(long int)m,(long int)k,(long int)onnz,(long int)nr,(long int)nc,(long int)nnz,(long int)nc,(long int)nr);
		RSB_FPRINTF(fd,"save /$LIBRSB_DICT 3 dict def $LIBRSB_DICT begin /M {moveto} bind def /Z {gsave currentpoint lineto %g setlinewidth 1 setlinecap stroke grestore} bind def /D {M Z} bind def /K {0.5 0.5 setrgbcolor} bind def\n",rd);
		RSB_FPRINTF(fd,"/R {rlineto} bind def\n");
		RSB_FPRINTF(fd,"/N {newpath} bind def\n");
		RSB_FPRINTF(fd,"/L {lineto} bind def\n");
		RSB_FPRINTF(fd,"/C {closepath} bind def\n");
		RSB_FPRINTF(fd,"/SLW {setlinewidth} bind def\n");
		RSB_FPRINTF(fd,"/SRGB {setrgbcolor} bind def\n");
		RSB_FPRINTF(fd,"/SCF {scalefont} bind def\n");
		RSB_FPRINTF(fd,"/SF {setfont} bind def\n");
		rsb__util_find_max(&max[0],VA,typecode,nnz,1);
		rsb__util_find_min(&min[0],VA,typecode,nnz,1);
		// RSB_STDERR("%lf %lf\n", *(double*)(&min[0]), *(double*)(&max[0]));
		rsb__util_do_negate(&min[0],typecode,1);
		rsb__util_vector_add(&max[0],&min[0],typecode,1);
		// RSB_STDERR("%lf %lf\n", *(double*)(&min[0]), *(double*)(&max[0]));
		if(RSB_IS_ELEMENT_NONZERO(&max[0],typecode))
			rsb__vector_scale_inv(VA,&max[0],typecode,nnz); /* scale elements in [0,1] */
		else
			rsb__fill_with_ones(VA,typecode,nnz,1);

		if(want_points_plot)
		{
			RSB_FPRINTF(fd,"%% dots plot\n");
			for(nzi=0;nzi<nnz;++nzi)
			{
				// RSB_FPRINTF(fd,"%d %d D\n",IA[nzi],JA[nzi]);
				// RSB_FPRINTF(fd,"%d %d D ",nr-1-IA[nzi],JA[nzi]);
				rsb_rf_t cv=0.0;
				RSB_NUMERICAL_TYPE_CAST_TO_ANY_P(rsb_rf_t,cv,typecode,VA,nzi);
				// cv=1.0f-cv;
				cv=0.5+cv/2.0; /* stronger */
				//cv=0.0+cv/2;
				// cv=0.5;
				// gray ... red
				// RSB_FPRINTF(fd,"%0.2f %0.2f %0.2f setrgbcolor\n",cv,0.5,0.5);
				RSB_FPRINTF(fd,"%.2f K ",cv);
				//RSB_FPRINTF(fd,"%d %d D ",nr-1-IA[nzi],JA[nzi]);
				RSB_FPRINTF(fd,"%ld %ld D ",(long int)JA[nzi],(long int)(nr-1-IA[nzi]));
				if(nzi%32==0)
					RSB_FPRINTF(fd,"\n");
			}
		}
		// RSB_FPRINTF(fd,"gsave grestore showpage\n");
		RSB_FPRINTF(fd,"stroke\n");
		goto err;
	}

	if(!all_nnz)
	{
		/* rsb__mtx_as_pixmap_resize is optional */
		if( RSB_SOME_ERROR(errval = rsb__mtx_as_pixmap_resize(VA, IA, JA, nnz, &nnz, m, k, height, width, typecode, flags)))
			goto err;
	}
	/*	if(m<=height)
			ys=1.0;
		else
			ys=((rsb_rf_t)height)/m;
		if(k<=width)
			xs=1.0;
		else*/
		ys=((rsb_rf_t)height)/m;
		xs=((rsb_rf_t)width)/k;
		csw=ys;
		csh=xs;
	
//	{
//		ys=((rsb_rf_t)height)/m;
//		xs=((rsb_rf_t)width)/k;
//	}
	if(width>k)
		dx=.1*csw, csw*=.8;
	else
		xs=csw=1.0;
	if(height>m)
		dy=.1*csh, csh*=.8;
	else
		ys=csh=1.0;
/*
	if(height>m)
		yps=ys*.8;
	else
		yps=ys;
	if(width>m)
		xps=xs*.8;
	else
		xps=xs;

	if(!all_nnz)
	{
		m=height;
		k=width;
	}*/

	rsb__do_print_postscript_header(fd, RSB_MARF_EPS, width, height, csw, csh);

	RSB_FPRINTF(fd,"%%%% nnz dump\n");
	RSB_FPRINTF(fd,"%%%% scales : %g %g\n",xs,ys);


	if(xs>1.0) xs=1.0;
	if(ys>1.0)ys=1.0;

	for(n=0;n<nnz;++n)
	{
		RSB_FPRINTF(fd,"%%%% at : %d %d\n",(int)IA[n],(int)JA[n]);
		RSB_FPRINTF(fd,
			"%g %g translate\n"
			".85 .85 .85 csquare\n"
			"-%g -%g translate\n"
			, dx+((rsb_rf_t) (JA[n]))*xs, -dy+((rsb_rf_t)height)-((rsb_rf_t) (IA[n]))*ys
			, dx+((rsb_rf_t) (JA[n]))*xs, -dy+((rsb_rf_t)height)-((rsb_rf_t) (IA[n]))*ys);
	}

	//RSB_FPRINTF(fd, "%%%%EOF\n");

err:
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE);
#endif /* RSB_ALLOW_STDOUT */
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_err_t rsb__dump_postscript_from_matrix(const char * filename, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_bool_t all_nnz)
{
	/**
	 \ingroup gr_internals
	 This function is experimentally used to render the sparse matrix.
	 Needs better error handling.
	 */
#if RSB_ALLOW_STDOUT
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
	rsb_time_t t=0;

	RSB_DO_FLAG_ADD(flags,RSB_FLAG_WANT_BCSS_STORAGE) ;
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORTED_INPUT) ;

	if(!filename )
		return RSB_ERR_BADARGS;

	t = - rsb_time();
	if( RSB_SOME_ERROR(errval = rsb__util_mm_load_matrix_f(filename, &IA, &JA,&VA , &m, &k, &nnz , typecode, flags, NULL, NULL)) )
		goto err;
	t += rsb_time();

#if 1
	errval = rsb_dump_postscript_from_coo(/*fd*/RSB_DEFAULT_FD, IA, JA, VA, m, k, nnz, br, bc, width, height, all_nnz, typecode);
#else
#if 0
	{
		RSB_STDERR("%s","FIXME: this is functioning code to render PostScript raster spy plots; it just needs the right place to be employed ..\n");
		FILE*fd = RSB_DEFAULT_FD;
	       	rsb_coo_idx_t ri=0,ci=0,nr=m,nc=k;
	       	const rsb_coo_idx_t nnr=16/*0*2*/;
		const rsb_coo_idx_t nnc=nc/(nr/nnr);
	       	rsb_nnz_idx_t nzi=0;
		errval = rsb__mtx_as_pixmap_resize(VA,IA,JA,nnz,&nnz,nr,nc,nnr,nnc,typecode,RSB_FLAG_NOFLAGS);
		if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
		nc=nnc;
		nr=nnr;
		errval = rsb__util_sort_row_major_inner(VA,IA,JA,nnz,nr,nc, typecode,RSB_FLAG_NOFLAGS);
		if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
		errval = rsb__util_compress_to_row_pointers_array(NULL,nnz,nr,RSB_FLAG_NOFLAGS,RSB_FLAG_NOFLAGS,IA);
		if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
		// errval = rsb__do_switch_rsb_mtx_to_csr_sorted(mtxAp, &VA, &IA, &JA, RSB_FLAG_NOFLAGS);
		if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
		// RSB_POSTSCRIPT_DUMP_COMMENT(fd,"raster dump\n");
		RSB_FPRINTF(fd,""
			"%%!PS-Adobe-3.0 EPSF-3.0\n"
			"%%%%Creator: "RSB_PACKAGE_STRING"\n"
			"%%%%Title: matrix spy plot\n"
			"%%%%CreationDate: \n"
			"%%%%DocumentData: Clean7Bit\n"
			"%%%%Origin: 0 0\n"
			"%%%%BoundingBox: 0 0 %d %d\n"
			"%%%%LanguageLevel: 2\n"
			"%%%%Pages: 1\n"
			"%%%%Page: 1 1\n"
		,nc,nr);
		RSB_FPRINTF(fd,"gsave\n""0 %d translate\n""%d %d scale\n""%d %d 8 [%d 0 0 -%d 0 0]\n"" {<",nr,nc,nr,nc,nr,nc,nr);
		for(ri=0;ri<nr;++ri)
		{
	       		rsb_coo_idx_t fc=0,lc=0;
	       		rsb_coo_idx_t crp=IA[ri],nrp=IA[ri+1];
			if(nrp==crp)
				lc=nc-1;
			else
				lc=JA[crp]-1;
			for(ci=fc;ci<=lc;++ci)
				RSB_FPRINTF(fd,"FF");
			for(nzi=crp;nzi<nrp;++nzi)
			{
				RSB_FPRINTF(fd,"00");
				fc=JA[nzi]+1;
				lc=fc-1;
				if(JA[nzi]==nc-1)
					lc=nc-1;
				else
				{
					if(nzi+1 < nrp)
						lc=JA[nzi+1]-1;
					else
						lc=nc-1;
				}
				for(ci=fc;ci<=lc;++ci)
					RSB_FPRINTF(fd,"FF");
			}
			RSB_FPRINTF(fd,"\n");
		}
		RSB_FPRINTF(fd,">}\n""image\n""grestore\n""showpage\n");
	}
#else
		goto err;
#endif
#endif
err:
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(VA);
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE);
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__dump_postscript_recursion_from_matrix(const char * filename, rsb_marf_t rflags, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_flags_t flags, rsb_bool_t want_blocks, rsb_bool_t z_dump , rsb_bool_t want_nonzeros, rsb_bool_t want_recursion, rsb_type_t typecode)
{
	/**
	 \ingroup gr_internals
	 */
#if RSB_ALLOW_STDOUT
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t *IA=NULL, *JA=NULL;
	void *VA=NULL;
	rsb_coo_idx_t m=0,k=0;
       	rsb_nnz_idx_t nnz=0;
	struct rsb_mtx_t * mtxAp=NULL;
	FILE*fd = RSB_DEFAULT_FD;

	if(!filename )
	{
		RSB_ERROR(RSB_ERRM_ES); 
		return RSB_ERR_BADARGS;
	}

	errval = rsb__util_mm_load_matrix_f(filename,&IA,&JA,&VA,&m,&k,&nnz,typecode,flags,NULL,NULL);
	if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,m,k,br,bc,flags,&errval);
	if(!mtxAp || RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err,RSB_ERRM_ES); 
	}

	if( rflags == RSB_MARF_LATEX_RECURSION )
	{
		errval = rsb__dump_latex_matrix_recursion(fd,mtxAp);
		goto ret;
	}

	if(want_recursion)
	{
		errval = rsb__dump_postscript_recursion_from_mtx_t(fd, NULL, mtxAp, br, bc, width, height, rflags|RSB_MARF_EPS_L /* FIXME */, want_blocks, z_dump , 0 /*want_nonzeros*/, NULL, NULL, NULL );
		if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}

	if( want_nonzeros )
	{
		// RSB_POSTSCRIPT_DUMP_COMMENT(fd,"nonzeros structure dump");
		errval = rsb_dump_postscript_from_coo(fd, IA, JA, VA, m, k, nnz, br, bc, width, height, want_nonzeros, typecode);
		want_nonzeros = 0;
		if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}

	RSB_CONDITIONAL_FREE(IA);
       	RSB_CONDITIONAL_FREE(JA);
       	RSB_CONDITIONAL_FREE(VA);

	RSB_MTX_FREE(mtxAp);/* we don't need it anymore here */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
	if( want_nonzeros )
	{
		RSB_POSTSCRIPT_DUMP_COMMENT(fd,"nonzeros structure dump");
		errval = rsb__dump_postscript_from_matrix(filename, br, bc, width, height, 1);
		if(RSB_SOME_ERROR(errval)) { RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

	errval = RSB_ERR_NO_ERROR;
	goto ret;
err:
	errval = RSB_ERR_GENERIC_ERROR;
ret:
	RSB_CONDITIONAL_FREE(IA); RSB_CONDITIONAL_FREE(JA); RSB_CONDITIONAL_FREE(VA);
	RSB_MTX_FREE(mtxAp);
	return errval;
#else  /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

static rsb_err_t rsb_dump_postscript_from_mtx_t(FILE*fd, const struct rsb_mtx_t*mtxAp, rsb_blk_idx_t br, rsb_blk_idx_t bc, int width, int height, rsb_bool_t all_nnz)
{
	struct rsb_coo_mtx_t coo;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!fd)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}

	RSB_INIT_COO_FROM_MTX(&coo,mtxAp);

	if(rsb__allocate_coo_matrix_t(&coo)!=&coo)
       	{
	       	RSB_PERR_GOTO(err,RSB_ERRM_ES); 
	}
	errval = rsb__do_get_coo(mtxAp,(rsb_byte_t**)(&coo.VA),&coo.IA,&coo.JA,RSB_FLAG_NOFLAGS);
	if(!RSB_SOME_ERROR(errval))
		errval = rsb_dump_postscript_from_coo(fd, coo.IA, coo.JA, coo.VA, coo.nr, coo.nc, coo.nnz, br, bc, width, height, all_nnz, mtxAp->typecode);
	if(RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(merr,RSB_ERRM_ES);
	}
merr:
	rsb__destroy_coo_matrix_t(&coo);
err:
	return errval;
}

rsb_err_t rsb__do_mtx_render(const char * filename, const struct rsb_mtx_t*mtxAp, rsb_coo_idx_t pmWidth, rsb_coo_idx_t pmHeight, rsb_marf_t rflags)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!mtxAp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP);
	}

	switch(rflags)
	{
		case(RSB_MARF_RGB):
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_UNIMPLEMENTED_YET);
		break;
		case(RSB_MARF_EPS_L):
		case(RSB_MARF_EPS_B):
		case(RSB_MARF_EPS_S):
		case(RSB_MARF_EPS):
		{
			FILE * fd = NULL;
			rsb_time_t dt = rsb_time();
			/* filename = filename ? filename : RSB_DEFAULT_DUMPFILENAME; */
			if( ! filename )
			{
		       		fd = RSB_DEFAULT_FD;
			}
			else
			{
		       		fd = rsb__util_fopen(filename,"w");
				if (!fd)
				{
					errval = RSB_ERR_BADARGS;
					RSB_PERR_GOTO(err,RSB_ERRM_ES);
				}
			}

			if( rflags == RSB_MARF_EPS || rflags == RSB_MARF_EPS_S || rflags == RSB_MARF_EPS_L )
				RSB_DO_ERROR_CUMULATE(errval,rsb_dump_postscript_from_mtx_t(fd,mtxAp,1,1,pmWidth,pmHeight,1));
			if( rflags == RSB_MARF_EPS || rflags == RSB_MARF_EPS_B || rflags == RSB_MARF_EPS_L )
				RSB_DO_ERROR_CUMULATE(errval,rsb__dump_postscript_recursion_from_mtx_t(fd,NULL,mtxAp,1,1,pmWidth,pmHeight,rflags,0,1,0,NULL,NULL,NULL));
			if( fd )
			{
				dt = rsb_time() - dt;
				RSB_FPRINTF(fd,"%% rendering time ~ %lg s\n",dt);
			}
			
			if( filename )
				fclose(fd);
		}
		break;
		default: {errval = RSB_ERR_UNIMPLEMENTED_YET; goto err;}
	}
err:
	return errval;
}

/* @endcond */
