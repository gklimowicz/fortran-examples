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
 * @brief Perfomance measuring/reporting code (this is mostly devel stuff).
 * @author Michele Martone
 * */

#include <strings.h>		/* bzero */
#include "rsb_internals.h"
#include "rsb-config.h"
#include <stdint.h> /* int64_t / uint64_t */

static struct rsb_global_reference_performance_info_t rsb_gpi;
static struct rsb_mbw_cm_t rsb_gmpi;		/**< a memory bandwidth benchmark record for each level+extra_level */
RSB_INTERNALS_COMMON_HEAD_DECLS

#if RSB_WANT_PERFORMANCE_FILE
static rsb_err_t rsb__read_global_reference_performance_info(struct rsb_global_reference_performance_info_t *gpip)
{
	/*!
	  \ingroup gr_internals
	   Reads in the system specific performance information binary file, which
	   has been created by exactly this library instance on this system.
	   TODO: Error handling is insufficient, e.g. in closing files.
	 */
	FILE *fp = rsb__util_fopen(rsb_global_session_handle.performance_binary_dump_file,"r+b");
	char sigbuf[RSB_PERFORMANCE_BINARY_DUMP_FILE_SIGNATURE_MAX_CHARS];
	rsb_int sl = rsb__strlen(RSB_PERFORMANCE_BINARY_DUMP_FILE_SIGNATURE);
	size_t s o =0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(!fp)
	{
		/* for a while, this will not be displayed so loudly */
		/*RSB_ERROR("error opening performance file (no such file?).\n");*/
		errval = RSB_ERR_NO_USER_CONFIGURATION;
		goto err;
	}
	if( fread(sigbuf,sl,1,fp)!=1 )
	{
		RSB_ERROR("problems reading performance file signature .\n");
		errval = RSB_ERR_GENERIC_ERROR;
		goto err;
	}
	sigbuf[sizeof(sigbuf)-1]='\0';
	if( strncmp(sigbuf,RSB_PERFORMANCE_BINARY_DUMP_FILE_SIGNATURE,sl) != 0 )
	{
		/* Warning: This could print out unterminated junk */
		RSB_ERROR("read an unknown performance file signature: %s.\n",sigbuf);
		errval = RSB_ERR_GENERIC_ERROR;
		goto err;
	}
	if( fread(&so,sizeof(size_t),1,fp)!=1 )
	{
		RSB_ERROR("problems reading performance file size.\n");
		errval = RSB_ERR_GENERIC_ERROR;
		goto err;
	}
	if( so != sizeof(*gpip) )
	{
		RSB_STDERR("perfomance file size (%zd) should be %zd!\n",so,sizeof(*gpip));
		errval = RSB_ERR_GENERIC_ERROR;
		goto err;
	}
	if(
		fread(gpip,sizeof(*gpip),1,fp)!=1 
	)
	{
		RSB_ERROR(RSB_ERRM_EQRPF);
		errval = RSB_ERR_GENERIC_ERROR;
		fclose(fp);
		goto err;
	}

	if(RSB_SOME_ERROR(rsb_load_bw_info(&rsb_gmpi, fp)))
	{
		errval = RSB_ERR_GENERIC_ERROR;
		goto err;
/*	if(rsb__print_mem_hier_timings(&rsb_gmpi))
  		goto err;*/
	}

	if( fclose(fp) != 0)
		RSB_DO_ERROR_CUMULATE(errval,RSB_ERR_GENERIC_ERROR);
err:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_WANT_PERFORMANCE_FILE */

static rsb_err_t rsb__load_performance_info(rsb_bool_t force_benchmark)
{
	/**
	  \ingroup gr_internals
	   loads performance info from file.
	 */

#if RSB_ALLOW_STDOUT
	/* NOTE : just a shortcut, not a real check */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(rsb_gmpi.mb!=NULL)
		return RSB_ERR_NO_ERROR;/* assume already loaded */

#if RSB_WANT_PERFORMANCE_FILE
	if(RSB_ERR_NO_ERROR != rsb__read_global_reference_performance_info(&rsb_gpi))
#else /* RSB_WANT_PERFORMANCE_FILE */
#endif /* RSB_WANT_PERFORMANCE_FILE */
	{
	
		if(!force_benchmark)
			return RSB_ERR_NO_ERROR;/* assume we'll load it later */

		RSB_STDOUT("there is no reference performance information available.\n");
		RSB_STDOUT("running reference performance benchmark ...\n");
		errval = rsb__do_referencebenchmark();
		RSB_STDOUT("..done.\n");
		if(RSB_ERR_NO_ERROR!=errval)
			goto err;/* NEW */
#if RSB_WANT_PERFORMANCE_FILE
		if(RSB_ERR_NO_ERROR != rsb__read_global_reference_performance_info(&rsb_gpi))
		{
			RSB_ERROR("A problem occurred reading global reference performance info.\n");
			/* filesystem may be full : it is not necessarily an internal error */
			return RSB_ERR_GENERIC_ERROR;
		/*	return RSB_ERR_INTERNAL_ERROR; */
		}
#else /* RSB_WANT_PERFORMANCE_FILE */
#endif /* RSB_WANT_PERFORMANCE_FILE */
	}
	rsb_gpi.initialized=1;
err:
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif /* RSB_ALLOW_STDOUT */
}

rsb_err_t rsb__perf_init(void)
{
	/**
	  \ingroup gr_internals

	  Blanks performance info structures.
	  Loads some performance info.
	*/
	RSB_BZERO_P(&rsb_gmpi);
	RSB_BZERO_P(&rsb_gpi);
	return rsb__load_performance_info(0);
}

rsb_err_t rsb__perf_exit(void)
{
	/**
	  \ingroup gr_internals
	  Frees performance info structures.
	 */
	if(rsb_gmpi.mb!=NULL)
	{
		RSB_CONDITIONAL_FREE(rsb_gmpi.mb); /*  we'll rather need a destructor function */
	}

	return RSB_ERR_NO_ERROR;
}

rsb_err_t rsb__print_mop_reference_performance_info_header(void)
{
#if RSB_ALLOW_STDOUT
	RSB_INFO("#type\top\t");
	RSB_INFO("rows\tcols\tbr\tbc\tmflops p.s.\te_mflops p.s.\tnnz\tfillin\n");
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

#if RSB_WANT_PERFORMANCE_FILE
static rsb_err_t rsb_print_mop_maxmins(const struct rsb_mop_reference_performance_info_t *pi)
{
	/**
	  \ingroup gr_internals
	   Analyzing this data further we can bound the optimization gains.
 	 */
#if RSB_ALLOW_STDOUT
	rsb_int si=0,ci=0,ri=0;
	rsb_int rua[] = RSB_ROWS_UNROLL_ARRAY;
	rsb_int cua[] = RSB_COLUMNS_UNROLL_ARRAY;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_int mri=0,mci=0,Mri=0,Mci=0;

	for(si=0;si<RSB_FITTING_SAMPLES;++si)
	{
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
		{
			if( pi->pipfs[si].m_flops[ri][ci] > pi->pipfs[si].m_flops[Mri][Mci] )
				Mri = ri,Mci=ci;
		}
		
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
		{
			if( pi->pipfs[si].m_flops[ri][ci] < pi->pipfs[si].m_flops[mri][mci] )
				mri = ri,mci=ci;
		}
		RSB_STDOUT("sample %zd : %zd %zd : max:%lg\n",
		(rsb_printf_int_t)si,(rsb_printf_int_t)rua[Mri],(rsb_printf_int_t)cua[Mci],pi->pipfs[si].m_flops[Mri][Mci]);
		RSB_STDOUT("sample %zd : %zd %zd : min:%lg\n",
		(rsb_printf_int_t)si,(rsb_printf_int_t)rua[mri],(rsb_printf_int_t)cua[mci],pi->pipfs[si].m_flops[mri][mci]);
	}
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_WANT_PERFORMANCE_FILE */

rsb_err_t rsb__print_mop_reference_performance_info(const struct rsb_mop_reference_performance_info_t *pi, char *s)
{
	/**
	  \ingroup gr_internals
	 * NEW
	 * but OBSOLETE
	 */
#if RSB_ALLOW_STDOUT
	rsb_int si=0,ci=0,ri=0;
	rsb_int rua[] = RSB_ROWS_UNROLL_ARRAY;
	rsb_int cua[] = RSB_COLUMNS_UNROLL_ARRAY;
	/*char *s_="";
	if(s)s_=s;*/
	for(si=0;si<RSB_FITTING_SAMPLES;++si)
	for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		RSB_STDOUT("%s"
			"%zd\t%zd\t"
			"%zd\t%zd\t"
			"%lg\t%lg\t"
			"%zd\t%lg\n",s,
			(rsb_printf_int_t)pi->pipfs[si].rows,
			(rsb_printf_int_t)pi->pipfs[si].cols,
			(rsb_printf_int_t)rua[ri],
			(rsb_printf_int_t)cua[ci],
			pi->pipfs[si].m_flops[ri][ci] /pi->pipfs[si].seconds[ri][ci],
			pi->pipfs[si].e_mflops[ri][ci]/pi->pipfs[si].seconds[ri][ci],
			(rsb_printf_int_t)pi->pipfs[si].nnz,
			pi->pipfs[si].fillin[ri][ci]
		);
	}
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
static rsb_int rsb_dump_mops_performance_info(const struct rsb_mops_performance_info_t *mpi)
{
	/**
	  \ingroup gr_internals
	   Dumps the whole struct in C format.
	   but OBSOLETE
	 */
#if RSB_ALLOW_STDOUT
	rsb_int oi;
	RSB_STDOUT("{\n");
	RSB_STDOUT(".pipmo={\n");
	for(oi=0;oi<RSB_IMPLEMENTED_META_MOPS;++oi)
	{
		rsb__dump_performance_info(mpi->pipmo+oi,"");
		RSB_STDOUT(",\n");
	}
	RSB_STDOUT("}\n");
	RSB_STDOUT("}\n");
	return 0;
#else /* RSB_ALLOW_STDOUT */
	return -1;
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_global_performance_info(const struct rsb_global_performance_info_t *gpip)
{
	/**
	  \ingroup gr_internals
	 * Dumps the whole struct in C format.
	 * but OBSOLETE
	 */
#if RSB_ALLOW_STDOUT
	rsb_int ti;
	RSB_STDOUT("#include \"rsb_krnl.h\"\n");
	RSB_STDOUT("struct rsb_global_performance_info_t gpi=\n");
	RSB_STDOUT("{\n");
	RSB_STDOUT(".gpi={\n");
	for(ti=0;ti<RSB_IMPLEMENTED_TYPES;++ti)
	{
		rsb_dump_mops_performance_info(gpip->gpi+ti);
		RSB_STDOUT(",\n");
	}
	RSB_STDOUT("}\n");
	RSB_STDOUT("}\n");
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#if RSB_WANT_PERFORMANCE_FILE
static rsb_int rsb_dump_reference_mop_performance_info(const struct rsb_mop_reference_performance_info_t *mpi)
{
	/**
	  \ingroup gr_internals
	   Dumps the whole struct in C format.
	   but OBSOLETE
	 */
#if RSB_ALLOW_STDOUT
	rsb_int oi;
	RSB_STDOUT("{ /* struct rsb_mop_reference_performance_info_t  */ \n");
	RSB_STDOUT(".pipfs={\n");
	for(oi=0;oi<RSB_FITTING_SAMPLES;++oi)
	{
		rsb__dump_performance_info(mpi->pipfs+oi,"");
		RSB_STDOUT(",\n");
	}
	RSB_STDOUT("},\n");

	RSB_STDOUT(".blocks_per_row=\n{");
	for(oi=0;oi<RSB_FITTING_SAMPLES;++oi)
		RSB_STDOUT("%lg,",mpi->blocks_per_row[oi]);
	RSB_STDOUT("},\n");

	/** alpha, beta, gamma parameterization as in the accels experimental setup*/
	rsb__dump_performance_array("alpha" ,(const double*)mpi->alpha);
	rsb__dump_performance_array("beta"  ,(const double*)mpi->beta);
	rsb__dump_performance_array("gamma" ,(const double*)mpi->gamma);

	RSB_STDOUT("}\n");
	return 0;
#else /* RSB_ALLOW_STDOUT */
	return -1;
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_WANT_PERFORMANCE_FILE */

#if RSB_WANT_PERFORMANCE_FILE
static rsb_int rsb_dump_reference_mops_performance_info(const struct rsb_mops_reference_performance_info_t *mpi)
{
	/**
	  \ingroup gr_internals
	   Dumps the whole struct in C format.
	   but OBSOLETE
	 */
#if RSB_ALLOW_STDOUT
	rsb_int oi;
	const char * mops[] = RSB_MATRIX_OPS_ARRAY;

	RSB_STDOUT("{ /* struct rsb_mops_reference_performance_info_t  */ \n");
	RSB_STDOUT(".pipmo={\n");
	for(oi=0;oi<RSB_IMPLEMENTED_META_MOPS;++oi)
	{
		RSB_STDOUT("/* mop is %s */\n",mops[oi]);
		rsb_dump_reference_mop_performance_info(mpi->pipmo+oi);
		RSB_STDOUT(",\n");
	}
	RSB_STDOUT("}\n");
	RSB_STDOUT("}\n");
	return 0;
#else /* RSB_ALLOW_STDOUT */
	return -1;
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_WANT_PERFORMANCE_FILE */

FILE *rsb__util_fopen(const char *path, const char *mode)
{
	/**
	  \ingroup gr_internals
	  A fopen wrapper.
	 */
/*	struct stat stat_s;
	if(-1==stat(path,&stat_s))return NULL;
	if( S_IFREG(stat_s.st_mode))return NULL;*/
	return fopen(path,mode);/* Flawfinder: ignore */
}

#if RSB_WANT_PERFORMANCE_FILE
static rsb_err_t rsb__save_bw_info(const struct rsb_mbw_cm_t *mi, FILE *fp)
{
	/*!
	  \ingroup gr_internals
	  Saves performance information on memory bandwidth.
	 */
	size_t sw=0;

	if(!fp || !mi)
	{
		RSB_ERROR("error.\n");
		return RSB_ERR_BADARGS;
	}

	sw = sizeof(*(mi->mb))*(mi->cln+mi->extra_level);

	if( fwrite(mi,sizeof(*mi),1,fp)!=1 )
		goto err;

	if(!mi->mb)
		goto err;/* should give an internal error ? */

	if( fwrite(mi->mb,sw,1,fp)!=1 )
		goto err;

	return RSB_ERR_NO_ERROR;
err:
	RSB_ERROR("error writing memory performance file.\n");
	return RSB_ERR_GENERIC_ERROR;
}
#endif /* RSB_WANT_PERFORMANCE_FILE */

#if RSB_WANT_PERFORMANCE_FILE
rsb_err_t rsb__save_global_reference_performance_info(const struct rsb_global_reference_performance_info_t *gpip)
{
	/*!
	  \ingroup gr_internals
	 * Writes out the system specific performance information binary file, which
	 * should be read by exactly this library instance on this system.
	 */
	FILE *fp = rsb__util_fopen(rsb_global_session_handle.performance_binary_dump_file,"w+b");
	size_t so=sizeof(*gpip);
	rsb_int sl=rsb__strlen(RSB_PERFORMANCE_BINARY_DUMP_FILE_SIGNATURE);

	if(!fp)
	{
		RSB_ERROR("error opening performance file %s.\n",rsb_global_session_handle.performance_binary_dump_file);
		return RSB_ERR_GENERIC_ERROR;
	}

	if( fwrite(RSB_PERFORMANCE_BINARY_DUMP_FILE_SIGNATURE,sl,1,fp)!=1 )
		goto err;

	if( fwrite(&so,sizeof(size_t),1,fp)!=1 )
		goto err;
	if( fwrite(gpip,sizeof(*gpip),1,fp)!=1 )
		goto err;

	if(rsb__mem_hier_timings(&rsb_gmpi))
		goto err;

/*	if(rsb__print_mem_hier_timings(&rsb_gmpi))
  		goto err;*/

	if(RSB_SOME_ERROR(rsb__save_bw_info(&rsb_gmpi, fp)))
		goto err;

	if( fclose(fp) == 0)
		return RSB_ERR_NO_ERROR;
	else
		return RSB_ERR_GENERIC_ERROR;

err:
	RSB_ERROR("error writing performance file.\n");
	fclose(fp);
	return RSB_ERR_GENERIC_ERROR;
}
#endif /* RSB_WANT_PERFORMANCE_FILE */

#if RSB_WANT_PERFORMANCE_FILE
static rsb_err_t rsb_load_bw_info(struct rsb_mbw_cm_t *mi, FILE *fp)
{
	/*!
	  \ingroup gr_internals
	  Loads performance information on memory bandwidth.
	  TODO: Error handling is insufficient.
	 */
	size_t sr=0;

	if(!fp || !mi)
	{
		RSB_ERROR(RSB_ERRM_ERROR);
		return RSB_ERR_BADARGS;
	}

	if( fread(mi,sizeof(*mi),1,fp)!=1 )
		goto err;

	sr = sizeof(*(mi->mb))*(mi->cln+mi->extra_level);
	mi->mb = rsb__calloc(sr);

	if(!mi->mb)
		goto ferr;
	
	if( fread(mi->mb,sr,1,fp)!=1 )
		goto ferr;

	return RSB_ERR_NO_ERROR;
ferr:
	RSB_CONDITIONAL_FREE(mi->mb);
err:
	RSB_BZERO_P(mi);
	RSB_ERROR(RSB_ERRM_ELMPF);
	return RSB_ERR_GENERIC_ERROR;
}
#endif /* RSB_WANT_PERFORMANCE_FILE */

#if RSB_WANT_PERFORMANCE_FILE
static rsb_err_t rsb__dump_global_reference_performance_info(const struct rsb_global_reference_performance_info_t *gpip)
{
	/**
	  \ingroup gr_internals
	   Dumps the whole struct in C format.
	   but OBSOLETE
	 */
#if RSB_ALLOW_STDOUT
	rsb_int ti;
	const char * types[] = RSB_MATRIX_TYPES_ARRAY;
	RSB_STDOUT("{ /* struct rsb_global_reference_performance_info_t */ \n");
	RSB_STDOUT(".initialized=%d,\n",gpip->initialized);
	RSB_STDOUT(".gpi={ \n");
	for(ti=0;ti<RSB_IMPLEMENTED_TYPES;++ti)
	{
		RSB_STDOUT("/* type is %s */\n",types[ti]);
		rsb_dump_reference_mops_performance_info(gpip->gpi+ti);
		RSB_STDOUT(",\n");
	}
	RSB_STDOUT("}\n");
	RSB_STDOUT("}\n");
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_WANT_PERFORMANCE_FILE */


rsb_err_t rsb__dump_current_global_reference_performance_info(void)
{
	/**
	  \ingroup gr_internals
	   but OBSOLETE
	 */
#ifdef 	RSB_WITH_FEEDBACK
	/* Warning: this is a dirty hack */
	return rsb__dump_global_reference_performance_info(&rsb_gpi);
#else /* RSB_WITH_FEEDBACK */
#if RSB_WANT_PERFORMANCE_FILE
	RSB_BZERO_P(&rsb_gpi);
	return rsb__dump_global_reference_performance_info(&rsb_gpi);
#else
	return RSB_ERR_NO_ERROR;
#endif
#endif /* RSB_WITH_FEEDBACK */
}

#if RSB_WANT_EXPERIMENTS_CODE
rsb_err_t rsb__dump_performance_info_line(const struct rsb_mop_performance_info_t * pi)
{
	/**
	  \ingroup gr_internals
	   NEW
	 */
#if RSB_ALLOW_STDOUT
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	rsb_blk_idx_t ri,ci;	/* row index, columns index */
	rsb_blk_idx_t rua[] = RSB_ROWS_UNROLL_ARRAY;
	rsb_blk_idx_t cua[] = RSB_COLUMNS_UNROLL_ARRAY;

	RSB_STDOUT(
		"#m\tk\t"
		"br\tbc\t"
		"nnz\t"
		"fillin\tm_flops\t"
		"e_mflops\t"
		"m.p.s.\t"
		"seconds\t"
		"\n"
		);

	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		{

			double perf;
			if( !pi->seconds[ri][ci] )continue;/* this for handling special cases of quasi empty records */

			perf= ( pi->m_flops[ri][ci]/pi->seconds[ri][ci])/pi->fillin[ri][ci] ;

			RSB_STDOUT(
				"%zd\t%zd\t"
				"%zd\t%zd\t"
				"%zd\t"
				"%lg\t"
				"%lg\t"
				"%lg\t"
				"%lg\t"
				"%lg"
				"\n"
				//"\n"
				,
				(rsb_printf_int_t)pi->rows,(rsb_printf_int_t)pi->cols,
				(rsb_printf_int_t)rua[ri],(rsb_printf_int_t)cua[ci],
				(rsb_printf_int_t)pi->nnz,
				pi->fillin[ri][ci],
				pi->m_flops[ri][ci],
				pi->e_mflops[ri][ci],/* effective mflops */
				perf,
				pi->seconds[ri][ci]
				);
		}
	}
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif /* RSB_ALLOW_STDOUT */
}
#endif /* RSB_WANT_EXPERIMENTS_CODE */

rsb_err_t rsb__dump_performance_info(const struct rsb_mop_performance_info_t * pi, const char * pid)
{
	/*!
	  \ingroup gr_internals
	   Another benchmark info dumping function.
          
           FIXME : UNFINISHED
          
	 * \return RSB_ERR_NO_ERROR on correct operation, an error code (see \ref errors_section) otherwise.
	 */
#if RSB_ALLOW_STDOUT
	if(!pi)
		return RSB_ERR_BADARGS;
	if(!pid)
		pid="pi";
	
	if(0)
	RSB_STDOUT("\n"
	"#define RSB_ROWS_UNROLL_ARRAY_LENGTH 4\n"
	"#define RSB_COLUMNS_UNROLL_ARRAY_LENGTH 4\n");

	RSB_STDOUT("{\n");
	RSB_STDOUT("/* rsb_mop_performance_info_t */\n");
	RSB_STDOUT(".rows=%zd,.cols=%zd,.nnz=%zd, /** some matrix info : size_t rows,cols,nnz; */\n",pi->rows,pi->cols,pi->nnz);
	rsb__dump_performance_array("m_flops" ,(const double*)pi->m_flops);
	rsb__dump_performance_array("e_mflops",(const double*)pi->e_mflops);
	rsb__dump_performance_array("fillin"  ,(const double*)pi->fillin);
	rsb__dump_performance_array("seconds" ,(const double*)pi->seconds);
	RSB_STDOUT("}\n");
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

#if 0
rsb_err_t rsb_print_all_system_info(void)
{
	/*
	  \ingroup gr_internals
	   NEW
	   but OBSOLETE
	 */

	rsb__sys_info();

	/* temporary dumpout */
	if(rsb_gpi.initialized)	/* FIXME : only partially */
		return rsb__dump_global_reference_performance_info(&rsb_gpi);
	else
		RSB_INFO("There is no hardcoded machine performance information in this build.\n");
	return RSB_ERR_NO_ERROR;
}
#endif

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_nnz_idx_t rsb__fillin_estimation_nnz_count(
	const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, 
	const  rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_flags_t flags, rsb_int nprobes
)
{
	rsb_nnz_idx_t pnnz=0;/* probing non zeros */
	// size_t el_size=0;
	rsb_int fraction=100;
//	rsb_err_t errval = RSB_ERR_NO_ERROR;

	/* NEW */
	/**
	  \ingroup gr_internals
		FIXME : unfinished
		\return 0 in case of error

		TODO : should follow some cache blocking/throughput based criteria
	*/
	const rsb_nnz_idx_t minnnz=10000;

	if(!IA || !JA || nnz<1)
	{
		RSB_ERROR("bad args to rsb__fillin_estimation_nnz_count()!");
		return 0;
	}

	// el_size = RSB_SIZEOF(typecode);

	/* FIXME : should be a FRACTION ! :) */
	pnnz=(nnz>minnnz)?minnnz:nnz;

	if(pnnz<minnnz)
		goto ok;

	if((nnz/fraction)*nprobes<=nnz)/* this is one among many probes */
	{
		pnnz=nnz/fraction;
	}

	if(pnnz<minnnz)
		pnnz=minnnz;
		
ok:
	return pnnz;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__estimate_expected_raw_performance_for_blocking(
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	const  rsb_nnz_idx_t nnz, rsb_type_t typecode, 
	rsb_flags_t flags,
	double efillin,
	double*eperf)
{
	/**
		FIXME : only BCSR !	
		FIXME : unfinished
	*/
#ifdef RSB_OPTYPE_INDEX_SPMV_UAUZ
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_blk_idx_t ri=0,ci=0;	/* row index, columns index */
	rsb_blk_idx_t rua[] = RSB_ROWS_UNROLL_ARRAY;
	rsb_blk_idx_t cua[] = RSB_COLUMNS_UNROLL_ARRAY;
	rsb_int si=0,oi = RSB_NUMERICAL_OP_INDEX_FROM_CODE(RSB_OPTYPE_INDEX_SPMV_UAUZ);
	rsb_int ti = RSB_NUMERICAL_TYPE_INDEX_FROM_CODE(typecode);
	/* FIXME */

	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		if(rua[ri]==mB)
			goto okr;
	}
	goto failr;
okr:
	for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
	{
		if(cua[ci]==kB)
			goto okc;
	}
	goto failc;
okc:

	#if 1
	*eperf = rsb_gpi.gpi[ti].pipmo[oi].pipfs[si].m_flops[ri][ci];
	#else
	/* mB per s */
	*eperf= (rsb_gmpi.mb[rsb_gmpi.cln].nr[RSB_MB_READ].t/
		 rsb_gmpi.mb[rsb_gmpi.cln].times) 
		* rsb_gmpi.mb[rsb_gmpi.cln].sz ;
	/* spmv per s  */
	*eperf/=((double)rsb_spmv_memory_accessed_bytes_(
		mB, kB,
		m,k,
		efillin*nnz,
		(efillin*nnz)*mB*kB,
		m/mB,
		RSB_SIZEOF(typecode)
		));

	/* mflops */
	*e/perf*=((2*nnz)*efillin)*(1.e-6);

	//RSB_STDERR("cache levels : %d + %d\n",rsb_gmpi.cln,rsb_gmpi.extra_level);
	rsb_int i=0;
	//for(i=0;i<rsb_gmpi.cln; ++i)
	//	RSB_STDERR("time at %d : %lg\n",i,rsb_gmpi.mb[i].nr[RSB_MB_READ].t);
	#endif
	goto err;
failr:
failc:
	*eperf=0;
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_OPTYPE_INDEX_SPMV_UAUZ */
	errval = RSB_ERR_UNSUPPORTED_OPERATION;
#endif /* RSB_OPTYPE_INDEX_SPMV_UAUZ */
err:
	RSB_DO_ERR_RETURN(errval)
}

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__estimate_expected_fillin_for_blocking(
	const void * VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, 
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	const  rsb_nnz_idx_t nnz, rsb_type_t typecode, 
	rsb_flags_t flags,
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	double *efillinp)
{
	/**
	  	\ingroup gr_internals

		Should estimate the fillin for a single matrix blocking.
		Assumes the input arrays to be already sorted for this blocking.

		FIXME : unfinished
		TODO  : should use const array arguments.
		FIXME : should not re-sort if sorted flag on.
		FIXME : only BCSR !
		FIXME : to work safely, should support Z ordering (which should perform nice for this purpose)
			AND copy the probing area in some temporary array..
		FIXME : IT IS SLOW SLOW SLOW
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	rsb_time_t mt;
//	rsb_int verbose=1;
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_nnz_idx_t pnnz=0;

	void *new_VA  = NULL;
	rsb_coo_idx_t *new_IA = NULL, *new_JA = NULL;

	if(!VA || !IA || !JA || nnz<1 || !efillinp)
	{
		RSB_ERROR(RSB_ERRM_ES);
		errval = RSB_ERR_BADARGS;
		goto err;
	}
				
	if(
	 ( flags & RSB_FLAG_SORTED_INPUT ) || 
	 ( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR ) 
	)
	{
		/* FIXME : disabled, because we are not prepared for Z sorted input ... */
		errval = RSB_ERR_UNIMPLEMENTED_YET;
		RSB_ERROR(RSB_ERRM_ES);
		goto err;
	}

#if 0
	if((!(flags & RSB_FLAG_SORTED_INPUT)) || (flags & RSB_FLAG_SORT_INPUT) )
	{
		RSB_ERROR(RSB_ERRM_ES);
		/* we want sorted input */
		errval = RSB_ERR_BADARGS;
		goto err;
	}				
#endif
	RSB_DO_FLAG_ADD(flags,RSB_FLAG_SORT_INPUT);

	new_VA = rsb__clone_area( VA , RSB_SIZEOF(typecode)    * nnz );
	new_IA = rsb__clone_area( IA , sizeof(rsb_coo_idx_t) * nnz );
	new_JA = rsb__clone_area( JA , sizeof(rsb_coo_idx_t) * nnz );

	if( !new_VA || !new_IA || !new_JA )
	{
		RSB_ERROR(RSB_ERRM_ES);
		errval = RSB_ERR_ENOMEM;goto err;
	}
	
	if(RSB_SOME_ERROR(errval))
        {
                RSB_ERROR(RSB_ERRM_SLIINS);
                goto err;/* NOTE : this jump will cause memory leaks */
        }

	pnnz = rsb__fillin_estimation_nnz_count( new_IA, new_JA, nnz, typecode, flags, 1 );

	if(pnnz<1)
	{
		RSB_WARN("rsb__fillin_estimation_nnz_count() gave pnnz<1!\n");
		errval = RSB_ERR_BADARGS;
		goto err;
	}				
	
#if 1
	if( flags & RSB_FLAG_QUAD_PARTITIONING)
		RSB_WARN("ignoring RSB_FLAG_QUAD_PARTITIONING in %s\n",__func__);
	if( flags & RSB_FLAG_AUTO_BLOCKING)
		RSB_WARN("ignoring RSB_FLAG_AUTO_BLOCKING in %s\n",__func__);
	if( flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR)
		RSB_WARN("ignoring RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR in %s\n",__func__);
	RSB_DO_FLAG_DEL(flags,RSB_FLAG_AUTO_BLOCKING);
        RSB_DO_FLAG_DEL(flags,RSB_FLAG_QUAD_PARTITIONING);   /* problems otherwise */
        RSB_DO_FLAG_DEL(flags,RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR);   /* problems otherwise */
#else
	if( flags & RSB_FLAG_AUTO_BLOCKING && 0) /* if 1, segfault */
	{
		*efillinp = RSB_REAL_ZERO;
		//RSB_DO_FLAG_DEL(flags,RSB_FLAG_AUTO_BLOCKING); // this flags causes trouble here (FIXME)
		goto ok;
	}
#endif
		
	RSB_DO_FLAG_DEL(flags,RSB_FLAG_OWN_PARTITIONING_ARRAYS);

	mt = - rsb_time();
	mtxAp = rsb__do_mtx_alloc_from_coo_const(new_VA,new_IA,new_JA,nnz,typecode,m,k,mB,kB,flags,&errval);
	mt += rsb_time();

	if(!mtxAp || (RSB_SOME_ERROR(errval)))
	{
		RSB_ERROR(RSB_ERRM_MBE);
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

	*efillinp = rsb__do_get_matrix_fillin(mtxAp);


	if(mtxAp)
		rsb__do_mtx_free(mtxAp);
err:
	rsb__do_perror(NULL,errval);
	RSB_CONDITIONAL_FREE(new_VA);
	RSB_CONDITIONAL_FREE(new_IA);
	RSB_CONDITIONAL_FREE(new_JA);
	goto ok;
ok:
	RSB_DO_ERR_RETURN(errval)
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
void rsb__pinfo_init(struct rsb_mtx_partitioning_info_t * pinfop,
	rsb_blk_idx_t M_b, rsb_blk_idx_t K_b,
	rsb_coo_idx_t *rpntr,rsb_coo_idx_t *cpntr,
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_blk_idx_t br, rsb_blk_idx_t bc)
{
	/**
	  \ingroup gr_internals
           Using this function as an initializer serves as a reminder
	   when changing this datatype definition.
	 */
	if(!pinfop)return;

	pinfop->rpntr = rpntr;
	pinfop->cpntr=cpntr;
	pinfop->nr=m;
	pinfop->nc=k;
	pinfop->M_b=M_b;
	pinfop->K_b=K_b;
	pinfop->br=br;
	pinfop->bc=bc;
}
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

rsb_err_t rsb__dump_system_performance_summary(void)
{
	/* TODO: find a better placement for this. */
#if RSB_ALLOW_STDOUT
#if RSB_WANT_PERFORMANCE_FILE
	rsb_int oi, ti;
	const char * types[] = RSB_MATRIX_TYPES_ARRAY;
	const char * mops[] = RSB_MATRIX_OPS_ARRAY;
#endif /* RSB_WANT_PERFORMANCE_FILE */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if RSB_WANT_PERFORMANCE_FILE
	errval = rsb__read_global_reference_performance_info(&rsb_gpi);
	if(RSB_SOME_ERROR(errval))
		goto err;

	errval = rsb__dump_global_reference_performance_info(&rsb_gpi);
	if(RSB_SOME_ERROR(errval))
		goto err;

	for(oi=0;oi<RSB_IMPLEMENTED_META_MOPS;++oi)
	for(ti=0;ti<RSB_IMPLEMENTED_TYPES;++ti)
	{
		RSB_STDOUT("%s %s:\n",types[ti],mops[oi]);
		errval = rsb_print_mop_maxmins(&(rsb_gpi.gpi[ti].pipmo[oi]));
		if(RSB_SOME_ERROR(errval))
			goto err;
	}
#else /* RSB_WANT_PERFORMANCE_FILE */
	errval = RSB_ERR_NO_USER_CONFIGURATION;
	if(errval == RSB_ERR_NO_USER_CONFIGURATION)
	{
		/* not a critical error; we can restore no error condition */
		errval = RSB_ERR_NO_ERROR;
		goto err;
	}
#endif /* RSB_WANT_PERFORMANCE_FILE */
err:
	RSB_DO_ERR_RETURN(errval)
#else /* RSB_ALLOW_STDOUT */
	RSB_DO_ERR_RETURN(RSB_ERR_UNSUPPORTED_FEATURE)
#endif /* RSB_ALLOW_STDOUT */
}

size_t rsb_spmv_memory_accessed_bytes_max(const struct rsb_mtx_t * mtxAp)
{
	/* upper bound is pessimistic access pattern */
	rsb_blk_idx_t columns;
	rsb_blk_idx_t rows;

	RSB_DEBUG_ASSERT(mtxAp);

	if(!rsb__is_bcsr_matrix(mtxAp))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	rsb__get_blocking_size(mtxAp, &rows, &columns);

	if(rows < 0 || columns < 0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	return 	
		mtxAp->el_size *
			(
			 mtxAp->element_count	/* 1 time each element */
			+(mtxAp->element_count /*/ rows*/)	/* the rhs, multiple times (one time per block width)  */
			+mtxAp->Mdim		/* the out vector, one time */
			)
		+
#if RSB_WANT_DBC
		sizeof(rsb_nnz_idx_t) * ( mtxAp->Mdim	/* bpntr */)*mtxAp->block_count
		+
		sizeof(rsb_nnz_idx_t) * ( mtxAp->block_count	/* bindx */)*mtxAp->block_count
#else
		sizeof(rsb_nnz_idx_t) * ( mtxAp->Mdim	/* bpntr */)*mtxAp->nnz
		+
		sizeof(rsb_nnz_idx_t) * ( mtxAp->nnz	/* bindx */)*mtxAp->nnz
#endif
		;
err:
	return 0;
}

size_t rsb_spmv_memory_accessed_bytes_min(const struct rsb_mtx_t * mtxAp)
{
	/* lower bound is reading matrix once */
	rsb_blk_idx_t columns;
	rsb_blk_idx_t rows;

	RSB_DEBUG_ASSERT(mtxAp);

	if(!rsb__is_bcsr_matrix(mtxAp))
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	rsb__get_blocking_size(mtxAp, &rows, &columns);

	if(rows < 0 || columns < 0)
		RSB_PERR_GOTO(err,RSB_ERRM_ES);

	return 	
		mtxAp->el_size *
			(
			 mtxAp->element_count	/* 1 time each element */
			+mtxAp->Mdim		/* the out vector, one time */
			+mtxAp->mdim		/* the rhs vector, one time */
			)
		+
		sizeof(rsb_nnz_idx_t) * ( mtxAp->Mdim	/* bpntr */)
		+
#if RSB_WANT_DBC
		sizeof(rsb_nnz_idx_t) * ( mtxAp->block_count	/* bindx */)
#else
		sizeof(rsb_nnz_idx_t) * ( mtxAp->nnz	/* bindx */)
#endif
		;
err:
	return 0;
}

size_t rsb_spmv_memory_accessed_bytes_(
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_nnz_idx_t element_count,
	rsb_nnz_idx_t block_count,
	rsb_blk_idx_t Mdim,
	size_t el_size
)
{
	/* (quasi) pessimistic  */

	if(mB < 0 || kB < 0)
	{
		RSB_ERROR("no blocking info supplied : can't estimate memory footprint.");
		return 0; /* error */
	}

	return 	
		el_size *
			(
			 element_count	/* 1 time each element */
			+(element_count / mB)	/* the rhs, multiple times (one time per block width)  */
			+m		/* the out vector, one time */
			)
		+
		sizeof(rsb_nnz_idx_t) * ( Mdim	/* bpntr */)
		+
		sizeof(rsb_nnz_idx_t) * ( block_count	/* bindx */)
		;
}

static size_t rsb_spmv_memory_accessed_bytes_leaf(const struct rsb_mtx_t * mtxAp)
{
	rsb_blk_idx_t bcolumns;
	rsb_blk_idx_t brows;
	
	RSB_DEBUG_ASSERT(mtxAp);

	if(!rsb__is_bcsr_matrix(mtxAp))
		return 0;

	rsb__get_blocking_size(mtxAp, &brows, &bcolumns);

	if(brows < 0 || bcolumns < 0)
		return 0; /* error */

	/* (quasi) pessimistic , in the sense that there is a lot of accesses */
	return   /* FIXME : possible overflow */
		mtxAp->el_size *
			(
			 mtxAp->element_count	/* 1 time each element */
			//+(mtxAp->element_count / brows)	/* the rhs, multiple times (one time per block width)  */
			+(mtxAp->element_count)	/* the rhs, multiple times (one time per block width)  */
			+mtxAp->nr		/* the out vector, one time */
			)
		+
		sizeof(rsb_nnz_idx_t) * ( mtxAp->Mdim	/* bpntr */)
		+
#if RSB_WANT_DBC
		sizeof(rsb_nnz_idx_t) * ( mtxAp->block_count	/* bindx */)
#else
		sizeof(rsb_nnz_idx_t) * ( mtxAp->nnz	/* bindx */)
#endif
		;
}

size_t rsb_spmv_memory_accessed_bytes(const struct rsb_mtx_t * mtxAp)
{
	rsb_submatrix_idx_t i,j;
	const struct rsb_mtx_t * submatrix;
	size_t sum = 0;

	RSB_DEBUG_ASSERT(mtxAp);

	if(rsb__is_recursive_matrix(mtxAp->flags))
	{
		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
			sum += rsb_spmv_memory_accessed_bytes(submatrix);
	}
	else
		sum = rsb_spmv_memory_accessed_bytes_leaf(mtxAp);
	
	return sum;
}

double rsb_spmv_memory_accessed_bytes_wr_ratio(const struct rsb_mtx_t * mtxAp)
{
	rsb_blk_idx_t columns;
	rsb_blk_idx_t rows;
	double rb,wb;

	RSB_DEBUG_ASSERT(mtxAp);

	if(!rsb__is_bcsr_matrix(mtxAp))
		return 0;

	rsb__get_blocking_size(mtxAp, &rows, &columns);

	if(rows < 0 || columns < 0)
		return 0; /* error */

	/* (quasi) pessimistic , in the sense that there is a lot of accesses */
	rb=(double)
		mtxAp->el_size *
			(
			 mtxAp->element_count	/* 1 time each element */
			+(mtxAp->element_count / rows)	/* the rhs, multiple times (one time per block width)  */
			)
		+
		sizeof(rsb_nnz_idx_t) * ( mtxAp->Mdim	/* bpntr */)
		+
#if RSB_WANT_DBC
		sizeof(rsb_nnz_idx_t) * ( mtxAp->block_count	/* bindx */)
#else
		sizeof(rsb_nnz_idx_t) * ( mtxAp->nnz	/* bindx */)
#endif
		;
	wb=(double)mtxAp->el_size*mtxAp->Mdim;		/* the out vector, one time */

	return wb/rb;
}

rsb_err_t rsb__dump_performance_record(const char * s, const struct rsb_mtx_t * mtxAp, rsb_real_t rsb_NMflops_ps, rsb_real_t rsb_RMflops_ps, const char *op, rsb_flags_t inflags)
{
	/**
		\ingroup gr_internals
		writes on stdout a line suitable for plotting and later analysis
		s and matrix are optional
	*/
#if RSB_ALLOW_STDOUT
	/* FIXME : buffer overflow risk */
	char buf[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];/* Flawfinder: ignore */
	/* rsb_NMflops_psm rsb_RMflops_ps :
	   algorithmic millions of ops per second */

	/* single line output, ideal for benchmark data to be processed later */
	RSB_STDOUT ("%-20s	%s",s,rsb__sprint_matrix_implementation_code2(mtxAp,buf,inflags));
	RSB_STDOUT ("	%.3lf	%lg",rsb_RMflops_ps,rsb_NMflops_ps);
	{
		RSB_STDOUT ("	");
		rsb__fprint_matrix_implementation_code(mtxAp,op,inflags,stdout);
	}
	RSB_STDOUT ("\n");

	return RSB_ERR_NO_ERROR ;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}

/* @endcond */
