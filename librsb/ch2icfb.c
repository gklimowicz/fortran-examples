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

/**
 * @file
 * @author Michele Martone
 * @brief
 * This standalone C 99 program produces ISO C BINDING oriented Fortran code of the rsb.h C header.
 * FIXME/TODO: BUFLEN of 1024 is too tight, and provokes segfaults at -O3 !!
 *
 * It is easy to break this program; e.g.:
 * #define a 1 comment_begin ..
 *  ... comment_end  
 * \internal
 *
 * Missing features:
 *  - handling of 'extern' 
 *  - proper preprocessor macros expansion
 *  - skipping function definitions
 *  - struct definitions
 *  - skipping RSB_RESTRICT, RSB_INLINE, restrict, double complex, RSB_EMPTY_FILE_FILLER , RSB_INNER_NRHS_SPMV_ARGS, RSB_OUTER_NRHS_SPMV_ARGS
 *  - skipping some enum's
 *  - rsb_time_t/void as return type
 *  - argv[], typedef char a[  as return types]
 * */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> // isblank, ...

#define BS 1
#define VERBOSE 0
#define WANT_BINDINGS 1
#define WANT_REAL_STAR 1
#define WANT_C_STRINGS 1 /* Use CHARACTER(C_CHAR) instead of TYPE(C_PTR) */
#define WANT_PURE 0
#define WANT_POINTERS_TO_INTEGERS 1 /* if 1, will have to invoke C_LOC(IA), ... */
#define PRINT_FUNC_LIST 0
#define WANT_CHAR_FOR_TYPECODE 1
#define WANT_PTR_COMMENT 1 /* FIXME */

#if WANT_REAL_STAR
#define C2F_DOUBLE "REAL(C_DOUBLE)"
#define C2F_FLOAT  "REAL(C_FLOAT)"
#else
/*
 integer, parameter :: rsb_dpk_=KIND(1.d0)
 integer, parameter :: rsb_spk_=KIND(1.s0)
 integer, parameter :: rsb_dpk_=SELECTED_REAL_KIND(15,300)
 integer, parameter :: rsb_spk_=SELECTED_REAL_KIND(6,37)
*/
#define C2F_DOUBLE "REAL(rsb_dpk_)"
#define C2F_FLOAT  "REAL(rsb_spk_)"
#endif

#if VERBOSE
#define INSPECT fprintf(stderr,"at 0:%c, at %d:%c\n",s[0],n,s[n])
#if 0
//#define WARN printf("in %s:%d\n",__FUNCTION__,__LINE__)
#else /* */
#define WARN fprintf(stderr,"in %s:%d (%d chars)\n",__FUNCTION__,__LINE__,n)
#endif /* */
#define DEBUG printf
#define INFO fprintf
#else /* VERBOSE */
#define WARN
#define INSPECT
#define DEBUG( ... )
#define INFO( ... )
#endif /* VERBOSE */

#define PRINTF  printf 
#if 0
#define IPRINTF printf
#else
#define IPRINTF( ... )
#endif
#define WANT_UNTYPED_VA 1
	const char * g_fortran_append = "";
	const char * g_fortran_prepend = "";

enum syntypes { NULLST, POPST, FUNCSIG, RETTYPE, FUNCSYM, ARGTYPES, ARGLABELS, INITIALIZE, FINALIZE }; 
struct ccs_t{ enum syntypes st; const char *s; size_t sl ; struct ccs_t * ss; } ; 
typedef struct ccs_t cc_t;
typedef int (*ccfp_t)(cc_t *);
struct pts_t{ ccfp_t fp; cc_t*ca; size_t cn,cc; int isstatic; } ; 
typedef struct pts_t pt_t;
#if 0
#define C2IFS (enum syntypes st, char *s, size_t sl)
#else
#define C2IFS (cc_t *cc)
#endif
typedef int (*c2ifp_t) C2IFS;
#define CCTMAX 100
#define TKLMAX 300
#define BUFLEN 1024*16 /* e.g. 1024 is too little ! */
#define ARGMAX 100
#define CCFPAPPLY(PT,ST,SS,SL) if((PT)&&((PT)->fp)){cc_t cc = {(ST),(SS),(SL),NULL};((PT)->fp)(&cc);/**/}
#define CCPRINT(CC) {if((CC)){PRINTF("%c..%c\n",(CC)->s[0],(CC)->s[(CC)->sl?(CC)->sl-1:0]);}}
#define CH2_MIN(X,Y) ((X)<(Y)?(X):(Y))		/*!< quick macro for minimum */

#if PRINT_FUNC_LIST
#define g_fnlmax 10000
 size_t g_fnlc = 0;
 char g_fns[g_fnlmax][TKLMAX];
 int  g_fnl[g_fnlmax];
#define FNADD(FNAME,FNL) {strncpy(&g_fns[g_fnlc][0],FNAME,FNL);g_fnl[g_fnlc]=FNL;/* printf("%s",g_fns[g_fnlc]);*/g_fnlc++;}
#define FDUMP() {int g_fni;printf("! ");for (g_fni=0;g_fni<g_fnlc;++g_fni){  printf("%s ",g_fns[g_fni]); }printf("\n"); }
#else
#define FNADD(FNAME,FNL) 
#define FDUMP() 
#endif /* PRINT_FUNC_LIST */

#define DOXY_SEE_STR "!> ISO C BINDING interface to ::"

size_t parse_c_substr_gen(const char*s, size_t l, const char *kw, int aat)
{
	size_t kl;

	kl = strlen(kw);
	if( l>=kl && strncmp(s,kw,kl)==0
		&& ( aat || !isalnum(s[kl]) ) )
		return kl;
	return 0;
}

size_t parse_c_type_kw(const char*s, size_t l, const char *kw)
{
	return parse_c_substr_gen(s, l, kw, 0);
}

size_t parse_c_substr(const char*s, size_t l, const char *kw)
{
	return parse_c_substr_gen(s, l, kw, 1);
}

const char * get_ptr_to_type_with_spec(const char*s, size_t l);

int a_primitive_type(const char *s, size_t l)
{
	if(parse_c_type_kw(s,l,"double")) return 1;
	if(parse_c_type_kw(s,l,"float")) return 1;
	if(parse_c_type_kw(s,l,"complex")) return 1;
	if(parse_c_type_kw(s,l,"double complex")) return 1;
	if(parse_c_type_kw(s,l,"float complex")) return 1;
#if 0
	if(parse_c_type_kw(s,l,"enum rsb_mif_t")) return 1;
	if(parse_c_type_kw(s,l,"enum rsb_elopf_t")) return 1;
#endif
	if(parse_c_type_kw(s,l,"int")) return 1;
	return 0;
}


int should_use_pointer_as_err_ref_in_type_string(const char *ls, size_t ll,const char *s, size_t l, char*ds)
{
	if(parse_c_type_kw(ls,ll,"errvalp"))
	{
#if WANT_PTR_COMMENT
		if(ds)
		       	strcpy(ds, " ! INTEGER(C_INT)");
#endif
		return 1;
	}
#if 0
	if(parse_c_type_kw(ls,ll,"rnz")) return 1;
#endif
	return 0;
}

int should_use_pointer_as_value_ref_in_type_string(const char *ls, size_t ll,const char *s, size_t l, char *ds)
{
	if(parse_c_type_kw(ls,ll,"VA")) goto isnt;
	if(parse_c_type_kw(ls,ll,"VAp")) goto isnt;
	if(parse_c_type_kw(ls,ll,"alphap")) goto isnt;
	if(parse_c_type_kw(ls,ll,"betap")) goto isnt;
	if(parse_c_type_kw(ls,ll,"alpha")) goto isnt;
	if(parse_c_type_kw(ls,ll,"beta")) goto isnt;
#if 0
	if(parse_c_type_kw(ls,ll,"infinity_norm")) return 1;
#endif
	if(parse_c_type_kw(ls,ll,"Np")==ll) goto isnt;
	if(parse_c_type_kw(ls,ll,"Dp")==ll) goto isnt;
	if(parse_c_type_kw(ls,ll,"Xp")==ll) goto isnt;
	if(parse_c_type_kw(ls,ll,"Yp")==ll) goto isnt;
	if(parse_c_type_kw(ls,ll,"iop")) goto isat;
	goto boh;
isnt:
#if WANT_PTR_COMMENT
	if(ds)
	       	strcpy(ds, " ! A single variable of the same numerical type of the matrix.");
#endif
	return 1;
isat:
#if WANT_PTR_COMMENT
	if(ds)
	       	strcpy(ds, " ! C_NULL_PTR is a safe value. Please consult the rsb.h documentation for other options.");
#endif
	return 1;
boh:
	return 0;
}

int should_use_pointer_as_index_array_in_type_string(const char *ls, size_t ll,const char *s, size_t l)
{
	if(parse_c_type_kw(ls,ll,"JA")) goto iia;
	if(parse_c_type_kw(ls,ll,"JAp")) goto iia;
	if(parse_c_type_kw(ls,ll,"IA")) goto iia;
	if(parse_c_type_kw(ls,ll,"IAp")) goto iia;
	if(parse_c_type_kw(ls,ll,"JAc")) goto iia;
	if(parse_c_type_kw(ls,ll,"IAc")) goto iia;
	return 0;
iia:
	return 1;
}

int should_use_pointer_as_ref_in_type_string(const char *ls, size_t ll,const char *s, size_t l, char *ds)
{
	return	
#if WANT_POINTERS_TO_INTEGERS
		should_use_pointer_as_index_array_in_type_string(ls,ll,s,l)+
#endif
		should_use_pointer_as_err_ref_in_type_string(ls,ll,s,l,ds)+
		should_use_pointer_as_value_ref_in_type_string(ls,ll,s,l,ds);
}

int should_use_pointer_as_value_array_in_type_string(const char *ls, size_t ll,const char *s, size_t l, char *ds)
{
#if 0
	if(parse_c_type_kw(ls,ll,"VA")==ll) return 1;
#endif
	s = get_ptr_to_type_with_spec(s,l);
	if(memchr(s,'*',l) && a_primitive_type(s,l))
		goto upava;
	return 0;
upava:
	return 1;
}

int should_use_pointer_as_char_c_string(const char *ls, size_t ll,const char *s, size_t l, char*ds)
{
#if WANT_C_STRINGS 
	if(parse_c_type_kw(ls,ll,"filename")) goto iss;
	if(parse_c_type_kw(ls,ll,"opvp")) goto iss;
	if(parse_c_type_kw(ls,ll,"opnp")) goto iss;
	if(parse_c_type_kw(ls,ll,"buf")) goto iss;
	if(parse_c_type_kw(ls,ll,"mis")) goto iss;
#endif
	return 0;
#if WANT_C_STRINGS 
iss:
#endif
#if WANT_PTR_COMMENT
	if(ds)
		strcpy(ds, " ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .");
#endif
	return 1;
}

int should_use_pointer_as_array_in_type_string(const char *ls, size_t ll,const char *s, size_t l, char*ds)
{
	return 
#if !WANT_POINTERS_TO_INTEGERS
		should_use_pointer_as_index_array_in_type_string(ls,ll,s,l)+
#endif
		should_use_pointer_as_value_array_in_type_string(ls,ll,s,l,ds);
}

size_t parse_c_blanks_and_comments(const char*s, size_t l);
size_t parse_substr(const char*s, size_t l, char*ss);

const char * c2f_rettype(const char * s, size_t l, char*ds)
{
	/* The following is very important.
	 * If we were capable of interpreting typedef's it would be much better.
	 * */
	s = get_ptr_to_type_with_spec(s,l);
	if(parse_c_type_kw(s,l,"char")) return "CHARACTER(C_CHAR)";
#if WANT_C_STRINGS 
	if(parse_c_type_kw(s,l,"rsb_char_t")) return "CHARACTER(C_CHAR), DIMENSION(*)";
#endif
	if(parse_c_type_kw(s,l,"double")) return C2F_DOUBLE;
	if(parse_c_type_kw(s,l,"float")) return C2F_FLOAT;
	if(parse_c_type_kw(s,l,"complex")) return "REAL(C_COMPLEX)";
	if(parse_c_type_kw(s,l,"double complex")) return "REAL(C_DOUBLE_COMPLEX)";
	if(memchr(s,'*',l))
		goto ip;
	if(parse_c_type_kw(s,l,"rsb_time_t")) return C2F_DOUBLE;
	if(parse_c_type_kw(s,l,"rsb_err_t")) goto ii;
	if(parse_c_type_kw(s,l,"int")) goto ii;
	if(parse_c_type_kw(s,l,"rsb_coo_idx_t")) goto li;
	if(parse_c_type_kw(s,l,"rsb_blk_idx_t")) goto ii;
	if(parse_c_type_kw(s,l,"rsb_nnz_idx_t")) goto li;
	if(parse_c_type_kw(s,l,"rsb_opt_t")) goto ii;
#if WANT_CHAR_FOR_TYPECODE
	if(parse_c_type_kw(s,l,"rsb_type_t")) return "INTEGER(C_SIGNED_CHAR)";
#else
	if(parse_c_type_kw(s,l,"rsb_type_t")) goto ii;
#endif
	if(parse_c_type_kw(s,l,"rsb_flags_t")) { sprintf(ds," " DOXY_SEE_STR "rsb_flags_t");goto ii; }
	/* if(parse_c_type_kw(s,l,"rsb_prec_flags_t")) goto ii; */
	if(parse_c_type_kw(s,l,"rsb_extff_t")) { sprintf(ds," " DOXY_SEE_STR "rsb_extff_t");goto ii; }
	if(parse_c_type_kw(s,l,"rsb_trans_t")) { /*sprintf(ds," " DOXY_SEE_STR" rsb_trans_t");*/goto ii; /* FIXME */ }
	if(parse_c_type_kw(s,l,"rsb_marf_t"))  { sprintf(ds," " DOXY_SEE_STR "rsb_marf_t");goto ii; }
	if(parse_c_type_kw(s,l,"rsb_bool_t"))  { sprintf(ds," " DOXY_SEE_STR "rsb_bool_t");goto ii; }
	if(parse_c_type_kw(s,l,"rsb_int_t")) goto ii;
	if(parse_c_type_kw(s,l,"blas_sparse_matrix")) goto ii;
#if 0
	if(parse_c_type_kw(s,l,"enum rsb_mif_t")) return "INTEGER(C_ENUM)";
	if(parse_c_type_kw(s,l,"enum rsb_elopf_t")) return "INTEGER(C_ENUM)";
	if(parse_c_type_kw(s,l,"rsb_mif_t")) return "INTEGER(C_ENUM)"; // FIXME: skipped "enum" !
	if(parse_c_type_kw(s,l,"rsb_elopf_t")) return "INTEGER(C_ENUM)"; // FIXME: skipped "enum" !
#endif
	if(parse_c_type_kw(s,l,"rsb_mif_t")) goto ii; // FIXME: skipped "enum" !
	if(parse_c_type_kw(s,l,"rsb_elopf_t")) goto ii; // FIXME: skipped "enum" !
	if(parse_c_type_kw(s,l,"rsb_precf_t")) goto ii;
	if(parse_c_type_kw(s,l,"blas_base_type")) goto ii;
	if(parse_c_type_kw(s,l,"blas_conj_type")) goto ii;
	if(parse_c_type_kw(s,l,"size_t")) return "INTEGER(C_SIZE_T)";
	return "(UNKNOWN TYPE)";
li:
#ifdef RSB_WANT_LONG_IDX_TYPE 
	return "INTEGER(C_RSB_INT_KND_)";
#else /* RSB_WANT_LONG_IDX_TYPE */
	return "INTEGER(C_RSB_INT_KND_)";
#endif /* RSB_WANT_LONG_IDX_TYPE */
ii:
	return "INTEGER(C_INT)";
ip:
#if WANT_PTR_COMMENT
	if(ds)
	{
		/* FIXME: very fragile and cheap test: breaks with substrings */
		if(parse_substr(s,l,"rsb_coo_idx_t")) sprintf(ds," ! INTEGER(C_INT)");
		if(parse_substr(s,l,"rsb_nnz_idx_t")) sprintf(ds," ! INTEGER(C_INT)");
		if(parse_substr(s,l,"rsb_int_t")) sprintf(ds," ! INTEGER(C_INT)");
		if(parse_substr(s,l,"rsb_flags_t")) sprintf(ds," ! INTEGER(C_INT)");
		if(parse_substr(s,l,"char")) sprintf(ds," ! CHARACTER(C_CHAR)");
		if(parse_substr(s,l,"void")) sprintf(ds," ! A numerical type");
		if(parse_substr(s,l,"rsb_mtx_t")) sprintf(ds," ! A matrix pointer variable: (TYPE(C_PTR),TARGET)");
		if(parse_substr(s,l,"rsb_real_t")) sprintf(ds," ! REAL*8");
	}
#endif
	return "TYPE(C_PTR)";
}

const char * c2f_argtype(const char * ls, size_t ll, const char * s, size_t l, char * ds)
{
	/* static const char*etl = "INTEGER(C_INT)"; */
#if WANT_UNTYPED_VA
	static const char*rtl = "TYPE(C_PTR),VALUE";
#else
	static const char*rtl = C2F_DOUBLE;
#endif
#if WANT_POINTERS_TO_INTEGERS
#else
	static const char*itl = "INTEGER(C_INT)";
#endif

	if(should_use_pointer_as_index_array_in_type_string(ls,ll,s,l))
	{
#if WANT_POINTERS_TO_INTEGERS
#if WANT_PTR_COMMENT
		if(ds)
		       	strcpy(ds, " ! INTEGER(C_INT)");
#endif
		return rtl;
#else
		return itl;
#endif
	}
	if(should_use_pointer_as_value_array_in_type_string(ls,ll,s,l,NULL))
	{
#if WANT_PTR_COMMENT
		if(ds)
		       	strcpy(ds, " ! An array of numerical type");
#endif
		return rtl;
	}
	if(should_use_pointer_as_err_ref_in_type_string(ls,ll,s,l,ds))
	{
		return rtl;
	}
	if(should_use_pointer_as_value_ref_in_type_string(ls,ll,s,l,ds))
	{
		return rtl;
	}
	return c2f_rettype(s,l,ds);
}

const char * c2f_varlabel(const char * s, size_t l)
{
	static char buf[BUFLEN];

	strncpy(buf,s,l);
	buf[l] = '\0';
	return buf;
}

const char * c2f_funname(const char * s, size_t l)
{
	static char buf[BUFLEN];

	strcpy(buf,"");
	strncat(buf,s,l);
#if 0
	strcat(buf,"_c2f");
#endif
	strcat(buf,g_fortran_append);
	strcat(buf,g_fortran_prepend);
	return buf;
}

int c2i C2IFS
{
	static enum syntypes s0 = NULLST;
	enum syntypes st = cc->st;
#if 0
	static cc_t ccts[CCTMAX];
	static pt_t pt = { c2i,ccts,0,CCTMAX,1 };
#endif
	static const char*funname = ""; 
	static size_t funnamel = 0; 
	static const char*rettype = ""; 
	static size_t rettypel = 0;
	static const char*arglabels[ARGMAX];
	static const char*argtypes[ARGMAX];
	static size_t arglabelsl[ARGMAX];
	static size_t argtypesl[ARGMAX];
	static size_t argn = 0;
#if WANT_UNTYPED_VA
	const char*modname = "rsb";
#else
	const char*modname = "rsb_d_mod";
#endif

#if 0
	IPRINTF("STATIC %d, ARG %d\n",s0,st);
#endif

	if(st==INITIALIZE)
	{
		/* PRINTF("MODULE %s\n   USE ISO_C_BINDING, ONLY: C_INT,C_SIZE_T,C_DOUBLE,C_PTR,C_NULL_PTR,C_CHAR\n\n",modname); */
		PRINTF("MODULE %s\n   USE ISO_C_BINDING, ONLY: C_INT,C_INT64_T,C_PTR,C_NULL_PTR"
#if WANT_CHAR_FOR_TYPECODE
				",C_SIGNED_CHAR"
#endif
				"\n\n",modname);
		/* PRINTF("MODULE %s\n   USE ISO_C_BINDING\n\n",modname); */
#if 0
		PRINTF("! MODULE constants:\n");
#endif

        	PRINTF("#ifdef RSB_WANT_LONG_IDX_TYPE\n");
        	PRINTF("INTEGER,PARAMETER :: RSB_IDX_KIND=8\n");
        	PRINTF("#define C_RSB_INT_KND_ C_INT64_T\n");
        	PRINTF("#else\n");
        	PRINTF("INTEGER,PARAMETER :: RSB_IDX_KIND=4\n");
        	PRINTF("#define C_RSB_INT_KND_ C_INT\n");
        	PRINTF("#endif\n");
	}
	else
	if(st==FINALIZE)
	{
		PRINTF("END MODULE %s\n",modname);
	}
	else
	if(s0==NULLST && st==FUNCSIG)
	{
#if 0
		IPRINTF("in %d\n",s0);
#endif
		IPRINTF("func begin:");
		s0 = st;
		argn = 0;
#if 0
		IPRINTF("in %d\n",s0);
#endif
	}
	else
	if(s0==FUNCSIG && st==FUNCSIG)
	{
#if 0
		/int ti,tc = 4;
		for(ti = 0;ti<tc;++ti)
#endif
	{
		size_t i;
		size_t nl=0,ol=0,mll=50;
		char buf[BUFLEN],*linbrk = "&\n  &";

		s0 = NULLST;
		IPRINTF("func end (%d args):",argn);
		strcpy(buf,"");
		strcat(buf,DOXY_SEE_STR);
		strcat(buf,c2f_funname(funname,funnamel));
		strcat(buf,".\n");
		strcat(buf,"INTERFACE\n");
		strcat(buf," ");
#if WANT_PURE
		if(0)
		{
			int ras = 0, rar = 0, rav =0;
			for(i = 0;i<argn;++i)
			{
				ras += should_use_pointer_as_array_in_type_string(arglabels[i],arglabelsl[i],argtypes[i],argtypesl[i]),NULL;
				rar += should_use_pointer_as_ref_in_type_string(arglabels[i],arglabelsl[i],argtypes[i],argtypesl[i]);
				rav += should_use_pointer_as_char_c_string(arglabels[i],arglabelsl[i],argtypes[i],argtypesl[i],NULL);
			}
			if( ras + rar + rav == 0)
				strcat(buf," PURE ");
		}
		else
		if(parse_c_type_kw(funname,funnamel,"rsb_strerror_r")) /* need to replace this with is_rsb_pure(funname,funnamel,arglabelsl)  */
			strcat(buf," PURE ");
#endif
		strcat(buf,c2f_rettype(rettype,rettypel,NULL));
		strcat(buf," FUNCTION ");
		strcat(buf,linbrk);
		strcat(buf,c2f_funname(funname,funnamel));
		strcat(buf,linbrk);
		strcat(buf,"(");
		ol=strlen(buf);
		for(i = 0;i<argn;++i)
		{
			if(i)strcat(buf,",");
			strncat(buf,arglabels[i],arglabelsl[i]);
			if(strlen(buf+ol)-nl*mll>mll)
				strcat(buf,"&\n&"),
				++nl;
		}
		strcat(buf,")");
		strcat(buf,linbrk);
		strcat(buf,"BIND(c,NAME = '");
		strncat(buf,funname,funnamel);
		strcat(buf,"')\n");
		strcat(buf," USE ISO_C_BINDING\n");
		for(i = 0;i<argn;++i)
		{
#if WANT_PTR_COMMENT
			char ds[BUFLEN];
#else
			char * ds = NULL;
#endif
			int as, ar, av;

#if WANT_PTR_COMMENT
			ds[0] = '\0';
#endif
			as = should_use_pointer_as_array_in_type_string(arglabels[i],arglabelsl[i],argtypes[i],argtypesl[i],ds);
			ar = should_use_pointer_as_ref_in_type_string(arglabels[i],arglabelsl[i],argtypes[i],argtypesl[i],ds);
			av = should_use_pointer_as_char_c_string(arglabels[i],arglabelsl[i],argtypes[i],argtypesl[i],ds);

			strcat(buf," ");
			strcat(buf,c2f_argtype(arglabels[i],arglabelsl[i],argtypes[i],argtypesl[i],ds));
			if(!as && !ar && !av)strcat(buf,", VALUE ");
			strcat(buf," :: ");
			strcat(buf,c2f_varlabel(arglabels[i],arglabelsl[i]));
			if(as)strcat(buf,"(*)");
#if 0
			strcat(buf," ! ");
			strncat(buf,argtypes[i],argtypesl[i]);
			strcat(buf," ");
			strncat(buf,arglabels[i],arglabelsl[i]);
#endif
#if WANT_PTR_COMMENT
#else
			if(ds)
#endif
				strcat(buf,ds);
			strcat(buf,"\n");
		}
		strcat(buf," END FUNCTION ");
		strcat(buf,c2f_funname(funname,funnamel));
		FNADD(c2f_funname(funname,funnamel),funnamel);
		strcat(buf,"\n");
		strcat(buf,"END INTERFACE");
		/* FIXME: if buf matches 'DOUBLE' shall reallocate it with different types and replicate it */
#if WANT_BINDINGS
		PRINTF("\n%s\n",buf);
#endif /* WANT_BINDINGS */
#if 0
		dump ...
		free ...
		IPRINTF("in %d\n",s0);
#endif
	}
	}
	else
	if(s0==FUNCSIG && st==RETTYPE)
	{
		/* ... */
		IPRINTF("rettype:");
		rettype = cc->s;
		rettypel = cc->sl;
	}
	else
	if(s0==FUNCSIG && st==FUNCSYM)
	{
		/* ... */
		IPRINTF("funcsym:");
		funname = cc->s;
		funnamel = cc->sl;
	}
	else
	if(s0==FUNCSIG && st==ARGTYPES)
	{
		/* append type ... */
		IPRINTF("argtypes %d:",argn);
		argtypes[argn] = cc->s;
		argtypesl[argn] = cc->sl;
	}
	else
	if(s0==FUNCSIG && st==ARGLABELS)
	{
		/* append type ... */
		IPRINTF("arglabels %d:",argn);
		arglabels[argn] = cc->s;
		arglabelsl[argn] = cc->sl;
	       	argn++;
	}
	else
	if(s0==FUNCSIG && st==POPST)
	{
		IPRINTF("pop %d:\n",argn);
		s0 = NULLST;
	       	argn = 0;
	}
	else
	{
		/*
		IPRINTF("unprocessed transition %d -> %d:%s\n",s0,st,cc->s);
		if(s0!=NULLST && st==POPST)
		if(st==POPST)
		if(s0==FUNCSIG && st==POPST)
			s0 = NULLST;
		*/
		goto no;
	}
	//CCPRINT(cc);
	goto ok;
ok:
	return 0;
no:
	return -1;
}

pt_t*realloc_pt(pt_t*pt)
{
	/* FIXME if isstatic, could not realloc ! */
	(pt)->ca = realloc((pt)->ca,sizeof(cc_t)*(pt)->cc);
	return pt;
}

pt_t*append_pc(pt_t*pt,cc_t *cc)
{
	if((pt)->cc-(pt)->cn==0)
		(pt)->cn += 16,
		realloc_pt(pt);
	IPRINTF("APPENDING %d\n",(pt)->ca[(pt)->cn].st);
	(pt)->ca[(pt)->cn++] = *cc;
	return pt;
}

pt_t*append_pt_l(pt_t*pt,enum syntypes st, char *s, size_t sl)
{
	cc_t cc = {st,s,sl,NULL};

	return append_pc(pt,&cc);
}

size_t parse_substr(const char*s, size_t l, char*ss)
{
	size_t sl = strlen(ss),n = 0;

	if(sl>l)
		return 0;
	for( ; n <= l-sl ; ++n )
		if(strncmp(s+n,ss,sl)==0)
		{
			n += sl;
			/* WARN; */
			return n;
		}
	return 0;
}

size_t parse_c_multi_line_comment(const char*s, size_t l)
{
	size_t n;

	if(l>=4 && s[0]=='/' && s[1]=='*')
	{
		n = parse_substr(s+2,l-2,"*/");
		if(n>=2)
		{
			n += 2;
			/* WARN; */
			return n;
		}
	}
	return 0;
}

size_t parse_c_line(const char*s, size_t l)
{
	size_t n = 0;

again:
	while(n<l && s[n]!='\n')
		++n;
	if(n>0 && s[n-1]=='\\' && s[n]=='\n')
	{n++;goto again;}
	if(s[n]=='\n')
		++n;
	/* DEBUG("%d\n",n); */
	return n;
}



size_t parse_c_one_line_comment(const char*s, size_t l)
{
	if(l>=2 && s[0]=='/' && s[1]=='/')
	{
		return 2+parse_c_line(s+2,l-2);
	}
	return 0;
}

size_t parse_c_comment(const char*s, size_t l)
{
	size_t n = parse_c_one_line_comment(s,l);

	if(n)
		goto ok;
	n = parse_c_multi_line_comment(s,l);
	if(n)
		goto ok;
	return 0;
ok:
	if(n)WARN;	
	return n;
}



int is_id_fchar(const char c)
{
	return isalpha(c) || c=='_';
}

int is_id_lchar(const char c)
{
	return is_id_fchar(c) || isdigit(c);
}

size_t parse_c_blanks(const char*s, size_t l)
{
	size_t n = 0;

	while(n<l && ( isblank(s[n]) || s[n]=='\n'  || s[n]=='\r' )) /* FIXME: does not handle DOS style newlines */
	{
		++n;
	}
	if(n)WARN;
	return n;
}

size_t parse_c_type_pointer_spec(const char*s, size_t l)
{
	size_t n = 0;

	if((n = parse_c_substr(s,l,"***"))) return n;
	if((n = parse_c_substr(s,l,"**"))) return n;
	if((n = parse_c_substr(s,l,"*"))) return n;
	return 0;
}
size_t parse_c_type_spec(const char*s, size_t l)
{
	size_t n = 0;

	if((n = parse_c_type_kw(s,l,"static"))) return n;
	if((n = parse_c_type_kw(s,l,"const"))) return n;
	if((n = parse_c_type_kw(s,l,"enum"))) return n;
	if((n = parse_c_type_kw(s,l,"struct"))) return n;
#if 0
	if((n = parse_c_type_kw(s,l,"void"))) return n;
#endif
	return parse_c_type_pointer_spec(s,l);
}

size_t parse_c_blanks_and_comments(const char*s, size_t l)
{
	size_t n = 0,nc = 0,nb = 0;

	do
	{
		nc = parse_c_comment(s+n,l-n);
		n += nc;
		nb = parse_c_blanks(s+n,l-n);
		n += nb;
	}
	while(nb+nc);
	return n;
}

size_t parse_c_int(const char*s, size_t l)
{
	size_t n = 0;

	while(isxdigit(s[n]) && l-n>0)
		++n;
	return n;
}

size_t parse_c_hexa(const char*s, size_t l)
{
	size_t n = 0;

	n = parse_substr(s+n,l-n,"0x");
	if(n==0)
		return 0;
	while(isxdigit(s[n]) && l-n>0)
		++n;
	return n;
}

size_t parse_c_type_specs(const char*s, size_t l)
{
	size_t n = 0,nn = parse_c_type_spec(s+n,l-n);

	n += nn;
	if(nn)
	do
	{
		size_t nc,ni;
		nc = parse_c_blanks_and_comments(s+n,l-n);
		if(!nc)goto ok;
		ni = parse_c_type_spec(s+n+nc,l-n-nc);
		if(!ni)goto ok;
		nn = ni+nc;
		n += nn;
	}
	while(nn>0 && n-l);
ok:
	if(n)WARN;
	return n;
}

size_t parse_c_identifier(const char*s, size_t l)
{
	if(l>0 && is_id_fchar(s[0]))
	{
		size_t n = 1;
		for(;is_id_lchar(s[n]);++n)
				;
		return n;
	}
	return 0;
}
size_t parse_c_tidentifier(const char*s, size_t l)
{
	return parse_c_identifier(s,l);
}

size_t parse_c_type_with_spec(const char*s, size_t l)
{
	size_t n = 0,nc = 0,nn = parse_c_type_specs(s+n,l-n);

	/* PRINTF("SPECS %d..%d:%s\n",0,nn,s); */
	if(nn)
	{
		size_t nc = parse_c_blanks_and_comments(s+nn,l-nn);
		if(!nc) goto ok;
		nn += nc;
	}
	n = parse_c_tidentifier(s+nn,l-nn);
	if(n)n += nn;
	nc = parse_c_blanks_and_comments(s+n,l-n);
	nn = parse_c_type_pointer_spec(s+n+nc,l-n-nc);
	if(nn)n += nc+nn;
#if 0
	n += parse_c_blanks_and_comments(s+n,l-n);
	PRINTF("STYPE %d..%d:%s\n",0,n,s+n);
#endif
	if(n)WARN;
ok:
	return n;
}

const char * get_ptr_to_type_with_spec(const char*s, size_t l)
{
	size_t n = 0;

	n += parse_c_type_specs(s+n,l-n);
	n += parse_c_blanks_and_comments(s+n,l-n);
#if 0
	n += parse_c_tidentifier(s+n,l-n);//type name
	n += parse_c_blanks_and_comments(s+n,l-n);
	n += parse_c_type_pointer_spec(s+n,l-n);
#endif
	return s+n;
}

size_t parse_c_void_func_args(const char*s, size_t l)
{
	size_t nc = 0,nv = 0;

	nc += parse_c_blanks_and_comments(s,l);
	nv = parse_c_type_kw(s+nc,l-nc,"void");
	nc += parse_c_blanks_and_comments(s+nv+nc,l-nc-nv);
	if(nv)
		return nc+nv;
	return 0;
}

size_t parse_c_func_args(const char*s, size_t l, pt_t*pt)
{
	size_t n = 0,nn = parse_c_void_func_args(s+n,l-n);

	if(nn==4) /* FIXME */
	{
		n = nn+1;
		goto ok;
	}
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	if(nn<l)
	do
	{
		size_t nc = 0,nv = 0,ni = 0;

		nn = 0;
		nv = parse_c_type_with_spec(s+n,l-n);
		if(!nv)goto ok;
		nc = parse_c_blanks_and_comments(s+n+nv,l-n-nv);
		/* if(!nc)goto ok; */
       		ni = parse_c_identifier(s+n+nc+nv,l-n-nc-nv);
		if(!ni)goto ok;
		CCFPAPPLY(pt,ARGTYPES,s+n,nv);
		CCFPAPPLY(pt,ARGLABELS,s+n+nc+nv,ni);
		nc += parse_c_blanks_and_comments(s+n+nv+nc+ni,l-n-nv-nc+ni);
		n += ni+nv+nc;
		if(n>=l)
			goto no;
		if(s[n]!=',')
		{
			++n;
		       	goto ok;/* FIXME: hack */
		}
		++n;
		nc = parse_c_blanks_and_comments(s+n,l-n);
		n += nc;
		nn = nc+ni+nv;
		/*
		PRINTF("OK ->%s\n",s+n);
		... parse label ...
		... parse comments ...
		... parse comma ...
		FIXME: UNFINISHED.CONTINUE HERE
		*/
	}
	while(nn>0 && n-l);
ok:
	if(n)WARN;
	return n;
no:
	return 0;
}

size_t parse_c_typename(const char*s, size_t l)
{
	size_t n = parse_c_identifier(s,l),nn = n;

	while(nn>0 && l-n>0)
	{
		size_t nc = parse_c_blanks_and_comments(s+n,l-n), ni = 0;
		if(!nc) goto ok;
		ni = parse_c_identifier(s+n+nc,l-n-nc);
		if(!ni) goto ok;
		nn = nc+ni;
		n += nn;
	}
ok:
	if(n)WARN;
	return n;
}

size_t parse_c_typedec(const char*s, size_t l)
{
	/* FIXME: UNFINISHED */
	return parse_c_identifier(s,l);
}

size_t parse_c_var_decl(const char*s, size_t l)
{
	size_t nn = 0,n = parse_c_typename(s,l);

	while(l-n>0 && (s[n]=='*' || s[n]==' '))++n;
	/* size_t nn = 0,n = parse_c_type_with_spec(s,l); */
	if(!n)
		return 0;
       	nn = parse_c_identifier(s+n,l-n);
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	/* nn = parse_substr(s+n,l-n,";"); */
	nn = parse_c_substr(s+n,l-n,";");
	if(!nn)
		return 0;
	/* FIXME: QUICK HACK */
	n += nn;
	WARN;
	return n;
}

size_t parse_c_var_mdef(const char*s, size_t l)
{
	size_t nn = 0, n = parse_c_typename(s,l);

	if(!n)
		return 0;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	nn = parse_c_substr(s+n,l-n,"=");
	if(nn==0) return 0;
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	nn = parse_c_hexa(s+n,l-n);
	if(nn==0)
	{
		nn = parse_c_int(s+n,l-n);
		if(nn==0)
			return 0;
	}
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	/* FIXME: QUICK HACK */
	n += nn;
	WARN;
	return n;
}

size_t parse_c_enum_decls(const char*s, size_t l)
{
	size_t n = 0, ll = 0;
	
	if(l-n>0)
	do
	{
		ll = 0;
		ll += parse_c_var_mdef(s+n+ll,l-n-ll);
		if(ll==0)return 0;
		/* WARN; */
		n += ll;
		/* INSPECT; */
		ll = parse_c_substr(s+n,l-n,",");
		n += ll;
		n += parse_c_blanks_and_comments(s+n,l-n);
	}
	while(l-n>0 && ll>0);
	if(n)WARN;
	return n;
}

size_t parse_c_var_decls(const char*s, size_t l)
{
	size_t n = 0,ll = 0;
	
	if(l-n>0)
	do
	{
		ll = 0;
		ll += parse_c_blanks_and_comments(s+n+ll,l-n-ll);
		ll += parse_c_var_decl(s+n+ll,l-n-ll);
		/* WARN; */
		n += ll;
		/* INSPECT; */
	}
	while(l-n>0 && ll>0);
	if(n)WARN;
	return n;
}

size_t parse_c_func_body(const char*s, size_t l)
{
	/* FIXME: UNFINISHED */
	return 0;
}

size_t parse_c_func_signature(const char*s, size_t l, pt_t*pt)
{
	size_t n = 0,nn = parse_c_type_with_spec(s,l);

	if(!nn)goto no;
	CCFPAPPLY(pt,FUNCSIG,s,0);
	CCFPAPPLY(pt,RETTYPE,s,nn);
	/* PRINTF("PIPPO %d..%d:%s\n",0,nn,s); */
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	if(!nn && (s[n-1]!='*'))goto no;
	n += nn;
       	nn = parse_c_identifier(s+n,l-n);
	if(!nn)goto no;
	CCFPAPPLY(pt,FUNCSYM,s+n,nn);
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	if(l-n==0 || s[n]!='(')
		goto no;
	++n;
#if 1
	nn = parse_c_func_args(s+n,l-n,pt);
#else
	/* FIXME: a quick hack. need a cycle, here. */
	nn = parse_substr(s+n,l-n,")");
#endif
	if(s[n+nn-1]!=')')goto no; // (void) will go, with this hack
	n += nn;
	goto yes;
yes:
	CCFPAPPLY(pt,FUNCSIG,s,n);
	WARN;
	/* *pt = append_pt_l(pt,FUNCSIG,s,n); */
	return n;
no:
	CCFPAPPLY(pt,POPST,s,0);
	return 0;
}

size_t parse_c_func_def(const char*s, size_t l, pt_t*pt)
{
	/* FIXME: UNFINISHED */
	size_t n = parse_c_func_signature(s,l,pt),nn = 0;

	if(!n)
		return 0;
	nn = parse_c_func_body(s+n,l-n);
	if(!nn)
		return 0;
	/* FIXME: HACK */
	return n+nn;
}

size_t parse_c_func_decl(const char*s, size_t l, pt_t*pt)
{
	size_t n = parse_c_func_signature(s,l,pt);

	if(!n)
		return 0;
	n += parse_c_blanks_and_comments(s+n,l-n);
	if(n>=l || s[n]!=';')
		return 0;
	n++;
	if(n)WARN;
	return n;
}

size_t parse_c_struct_def(const char*s, size_t l)
{
	/* size_t n = parse_substr(s,l,"struct"),nn = 0; */
	size_t n = parse_c_type_kw(s,l,"struct"),nn = 0;

	if(!n) return 0;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	if(!nn) return 0;
	n += nn;
       	nn = parse_c_identifier(s+n,l-n);
	if(!nn) return 0;
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	nn = parse_c_substr(s+n,l-n,"{");
	/* nn = parse_substr(s+n,l-n,"{"); */
	if(!nn) return 0;
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	/* FIXME: identifier declarations */
	nn = parse_c_var_decls(s+n,l-n);
	if(!nn) return 0;
	n += nn;
	nn = parse_c_substr(s+n,l-n,"}");
	/* nn = parse_substr(s+n,l-n,"}"); */
	if(!nn) return 0;
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	nn = parse_substr(s+n,l-n,";");
	if(!nn) return 0;
	n += nn;
	WARN;
	return n;

}

size_t parse_c_preproc_line(const char*s, size_t l)
{
	size_t n = 0;
	size_t nn = parse_c_blanks_and_comments(s+n,l-n);

	n = parse_c_substr(s+nn,l-nn,"#");
	WARN;
	if(!n)
		return 0;
	n += nn;
	n += parse_c_line(s+n,l-n);
	return n;
}

size_t parse_c_enum_def(const char*s, size_t l)
{
	/* size_t n = parse_substr(s,l,"enum"),nn = 0; */
	size_t n = parse_c_type_kw(s,l,"enum"),nn = 0;

	if(!n) return 0;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	if(!nn) return 0;
	n += nn;
       	nn = parse_c_identifier(s+n,l-n);
	if(!nn) return 0;
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	nn = parse_substr(s+n,l-n,"{");
	if(!nn) return 0;
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	/* FIXME: identifier declarations */
	/* nn = parse_c_var_decls(s+n,l-n); */
	nn = parse_c_enum_decls(s+n,l-n);
	if(!nn) return 0;
	n += nn;
	nn = parse_c_substr(s+n,l-n,"};");
	if(!nn) return 0;
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	WARN;
	return n;
}

size_t parse_c_special(const char*s, size_t l)
{
	/* FIXME: UNFINISHED */
	const char *el = "extern \"C\" {";
	size_t n = 0,ell = strlen(el);

	if(l>=ell && strncmp(s,el,ell)==0)
		n = ell;
	if(n)
	{
		n += parse_c_line(s+n,l-n);
		goto ok;
	}
	if(l>0 && s[0]=='}') /* FIXME: this is as a complement to 'extern "C"'  */
	{
		n = 1;
		goto ok;
	}
ok:
	if(n)WARN;
	return n;
}

size_t parse_c_typedef(const char*s, size_t l)
{
	const char *ts = "typedef";
	size_t n = 0,nn = 0,sl = strlen(ts);

	if(l>=sl && strncmp(s,ts,sl)==0)
		n = sl;
	if(!n)
 	      return 0;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	if(!nn)
		return 0;
	n += nn;
	nn = parse_c_typename(s+n,l-n);
	if(!nn)
		return 0;
	/* printf("nn at %d:%c%c\n",nn,s[n],s[n+1]); */
	n += nn;
	nn = parse_c_blanks_and_comments(s+n,l-n);
	n += nn;
	/* WARN; INSPECT; */
	if(l-n<1 || s[n]!=';')
		return 0;
	/* FIXME: simplifications */
	n += 1;
	WARN; INSPECT;
	return n;
}

size_t parse_c_multiline_lines(const char*s, size_t l)
{
	size_t n = parse_c_line(s,l),nn = n;
	
	while(n>1 && ( s[nn-2]=='\\' || s[nn-1]=='\\' ) && l-nn>0)
	{
		n = parse_c_line(s+nn,l-nn);
		nn += n;
	}
	if(nn)n = nn;
	/* INSPECT; */
	if(n)WARN;
	return n;
}

size_t parse_c_prepcode(const char*s, size_t l)
{
	size_t n = 0;

	if(l>0 && s[0]=='#')
	{
		n++;
		n += parse_c_multiline_lines(s+1,l-n);
	}
	if(!n)
		goto no;
	WARN;
no:
	return n;
}

size_t parse_c_header(const char*s, size_t l, pt_t*pt, const int pmbe)
{
	size_t tl = 0,ll = 0;
#if 0
	while((ll = parse_c_line(s+tl,l-tl))>0)
		tl += ll;
	//DEBUG("%d\n",tl);
#else
	if(pmbe)
	CCFPAPPLY(pt,INITIALIZE,NULL,0);

	if(l-tl>0)
	do
	{
		ll = 0;
		ll += parse_c_blanks(s+tl+ll,l-tl-ll);
	 	ll += parse_c_prepcode(s+tl+ll,l-tl-ll);
		ll += parse_c_special(s+tl+ll,l-tl-ll);
		ll += parse_c_comment(s+tl+ll,l-tl-ll);
		ll += parse_c_typedef(s+tl+ll,l-tl-ll);
		ll += parse_c_struct_def(s+tl+ll,l-tl-ll);
		ll += parse_c_enum_def(s+tl+ll,l-tl-ll);
		ll += parse_c_preproc_line(s+tl+ll,l-tl-ll);
		ll += parse_c_func_decl(s+tl+ll,l-tl-ll,pt);
		ll += parse_c_var_decl(s+tl+ll,l-tl-ll);
	/*	
		ll += parse_c_func_def(s+tl+ll,l-tl-ll,pt);
		*/
		tl += ll;
		if(ll==0 && tl<l)
		{
			size_t bl=CH2_MIN(100,l-tl);
			char buf[bl];
			strncpy(buf,s,bl-1);
			buf[bl-1]='\0';
			INFO(stderr,"terminating prematurely and char %d / %d:\n...\n\"%s\"\n...\n!\n",tl,l,buf);
		}
	}
	while(l-tl>0 && ll>0);
	if(pmbe)
	CCFPAPPLY(pt,FINALIZE,NULL,0);
#endif
	return tl;
}

void * do_realloc(void *p, size_t n)
{
	return realloc(p,n);
}

int main(void)
{
	char *s = NULL;
	size_t rb = 0,sb = 0,cs = BS,pc = 0;
	ssize_t rn = -1;
	int ret = 0;
	pt_t pt = { c2i,NULL,0,CCTMAX,1 };
	const int pmbe = !getenv("CH2ICFB_NH");

	/* bzero(&pt,sizeof(pt)); */

	if((s = do_realloc(s,cs)))
		rb += cs,cs *= 2;

	/* chomp */
	while(rn)
	{
		while( rb>0 && (rn = read(0,s+sb,rb))>0 )
			sb += rn,rb -= rn;

		/* printf("REALLO %d+%d %d: %p(%d)\n",sb,rb,rn,s,cs); */
		if(rn && rb==0)
			if((s = do_realloc(s,sb+cs)))
				rb += cs, cs *= 2;
	}

	if(pmbe)
	{
		PRINTF("!> @file.\n");
		PRINTF("!! @brief Header file automatically generated from <rsb.h>, offering ISO-C-BINDING interfaces to <rsb.h>'s functions.\n");
		PRINTF("!! Defines \\c MODULE \\c rsb.\n");
		PRINTF("!! For examples of usage, see Fortran examples in \\ref rsb_doc_examples.\n");
		PRINTF("!! The official documentation is that of <rsb.h>.\n");
		PRINTF("!! Make sure you are using a modern Fortran compiler.\n\n");
		/* PRINTF("!> @cond INNERDOC\n"); */
		PRINTF("!DEC$IF .NOT. DEFINED (RSB_FORTRAN_HEADER)\n!DEC$DEFINE RSB_FORTRAN_HEADER\n\n");
	}
	if((pc = parse_c_header(s,sb,&pt,pmbe))==sb)
	{ INFO(stderr,"header file parsed (%d chars parsed).\n",pc);ret = 0; }
	else
	{ INFO(stderr,"header file NOT parsed (%d chars parsed out of %d).\n",pc,sb);ret = 1; }
	FDUMP()
	if(pmbe)
		PRINTF("\n!DEC$ENDIF\n\n");
	/* PRINTF("!> @endcond\n"); */
	goto err;
err:
	if(s)
		free(s);
	return ret;
}
/* @endcond */
