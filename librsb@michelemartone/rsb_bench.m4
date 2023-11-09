dnl
dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
/*!
 @file
 @brief

 Former performance info gathering code; now obsoleted and used as test.
 */
dnl
include(`rsb_misc.m4')dnl
dnl
RSB_M4_HEADER_MESSAGE()dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_BENCH_H_INCLUDED
#define RSB_BENCH_H_INCLUDED
')
dnl
include(`do_unroll.m4')dnl
include(`rsb_krnl_vb_macros.m4')dnl
include(`rsb_krnl_macros.m4')dnl
dnl
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
dnl
#include "rsb_internals.h"
dnl
#ifdef RSB_HAVE_CBLAS_H
#include <cblas.h>
#endif /* RSB_HAVE_CBLAS_H */
#ifdef RSB_HAVE_CLAPACK_H
#include <clapack.h>
#endif /* RSB_HAVE_CLAPACK_H */
#include <math.h>
dnl
dnl
dnl
dnl	RSB_M4_HYPERBOLIC_FITTING_FUNCTION_ARGS()
dnl	---------------------------------
dnl
define(`RSB_M4_HYPERBOLIC_FITTING_FUNCTION_ARGS',`dnl
dnl
`(double x[], double y[], size_t nb_loop, double * a, double * b, double *c, double c_s)'dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_HYPERBOLIC_FITTING_FUNCTION_IDENTIFIER()
dnl	---------------------------------
dnl
define(`RSB_M4_HYPERBOLIC_FITTING_FUNCTION_IDENTIFIER',`dnl
dnl
`rsb__fit_hyp'dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_HYPERBOLIC_FITTING_FUNCTION()
dnl	---------------------------------
dnl
define(`RSB_M4_HYPERBOLIC_FITTING_FUNCTION',`dnl
dnl
rsb_err_t RSB_M4_HYPERBOLIC_FITTING_FUNCTION_IDENTIFIER()`'dnl
RSB_M4_HYPERBOLIC_FITTING_FUNCTION_ARGS()`'dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
{
#if !(RSB_HAVE_CLAPACK && RSB_HAVE_CBLAS)
	return RSB_ERR_UNSUPPORTED_OPERATION;
#else
	/**
	 * \ingroup gr_bench
         * Note : 
	 * 
	 * This function will compute a performance predictor based on
         * nonzero per row ratio, by fitting the two input x (non zeros per row)
         * and y (megaflops) vectors (both with n = RSB_FITTING_SAMPLES points) to
         * the following formula :
         *
         *           `megaflops (nnz_per_row) a + b / ( c + nnz_per_row )'
         *
         * The c_s and nb_loop arguments will be documented some day.
         *
	 * This model is discussed in the following article :

@article{ButtEijkLang:spmvp,
  title = {Performance Optimization and Modeling of Blocked Sparse Kernels},
  author = {Buttari, Alfredo and Eijkhout, Victor and Langou, Julien and Filippone, Salvatore},
  pages = {467--484},
  year = 2007,
  journal = {IJHPCA},
  volume = 21,
  url = {\url{{http://www.tacc.utexas.edu/~eijkhout/Articles/2007-buttari-spmvp.pdf}}}
}
         *
         */

	rsb_int nparms=3;
	rsb_int n = RSB_FITTING_SAMPLES;
	/* Fortran arrays */
#define RSB_FORTRAN_ARRAY(AI,ROWS,COLS) AI[(ROWS)*(COLS)]

	rsb_int nj = 3;
	rsb_int i,j;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	double RSB_FORTRAN_ARRAY(G ,n,3);
	double RSB_FORTRAN_ARRAY(G1,n,3);
	double RSB_FORTRAN_ARRAY(GG,3,3);
	double RSB_FORTRAN_ARRAY(z ,n,1);
	double RSB_FORTRAN_ARRAY(z0,n,1);
	double RSB_FORTRAN_ARRAY(dy,n,1);
	double RSB_FORTRAN_ARRAY(ddy,3,1);
	double RSB_FORTRAN_ARRAY(xj ,nj,1);
	double RSB_FORTRAN_ARRAY(yj ,nj,1);
	double RSB_FORTRAN_ARRAY(zj ,nj,1);

	double xcpy[n];
	double a_t,b_t,sum1,sum2,sum3,sum4,error,tmp_a,tmp_b,tmp_c, min_err,max,min,avg,intl;
  	int /*i,*/info,ipivot[3],/*nj,j,*/k,cnt;
	rsb__memcpy(xcpy,x,sizeof(xcpy));	/* not a bit more .. and please note that sizeof(x)=sizeof(double*) != sizeof(x[n])*/


	RSB_INFO("starting analysis...\n");
	RSB_STDOUT("\n");
	RSB_STDOUT("performance data:\n");
	for(i=0;i<n;++i)
	{
		RSB_STDOUT("%lg %lg\n",xcpy[i],y[i]);
	}

	sum1=0;
	sum2=0;
	sum3=0;
	sum4=0;


  	*a=y[n-1];
	
	rsb__memcpy(xj,x,sizeof(xj));	/* not a bit more */
	rsb__memcpy(yj,y,sizeof(yj));	/* not a bit more */

	for(i=0;i<nj;++i)
  	{
		zj[i]=yj[i]-*a;
  		zj[i]=1/zj[i];
	}

	for(i=0;i<nj;++i)
	{
		sum1=sum1 + xj[i]*zj[i];
		sum2=sum2 + xj[i];
		sum3=sum3 + zj[i];
		sum4=sum4 + xj[i]*xj[i];
	}

	a_t= (sum3*sum4-sum2*sum1)/(nj*sum4-sum2*sum2);
	b_t=(nj*sum1 - sum2*sum3) / (nj*sum4 - sum2*sum2);

  	*b=1/b_t;
	*c=a_t* *b;

	for(i=0;i<n;++i)
		z0[i]= *a +*b/(x[i]+*c);

	error = 0;
	for(j=0;j<n;++j)
		error = error + (fabs( z0[j] - y[j] ) / y[j] );

	error = error / n * 100;

	min_err=error;

	tmp_a=*a;
	tmp_b=*b;
	tmp_c=*c;

	for(i=0;i<nb_loop;++i)
	{
		for(j=0;j<n;++j)
			dy[j] = z0[j]-y[j];

		for(j=0;j<n;++j)
		{
			G[j+0*n]=1;
			G[j+1*n]=1/(x[j]+tmp_c);
			G[j+2*n]=-tmp_b/( (x[j]+tmp_c)*(x[j]+tmp_c) );

			G1[j+0*n]= G[j+0*n];
			G1[j+1*n]= G[j+1*n];
			G1[j+2*n]= G[j+2*n];
		}

#if 
		cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,3,3,n,1.0,G,n,G1,n,0.0,GG,3);
		errval =  clapack_dgetrf(CblasColMajor,3,3,GG,3,ipivot);
		if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
		cblas_dgemv(CblasColMajor,CblasTrans,n,3,1.0,G,n,dy,1,0.0,ddy,1);
		errval =  clapack_dgetrs(CblasColMajor,CblasNoTrans,3,1,GG,3,ipivot,ddy,3);
		if(RSB_SOME_ERROR(errval))
			RSB_PERR_GOTO(err,RSB_ERRM_ES);
#else /* (RSB_HAVE_CLAPACK && RSB_HAVE_CBLAS) */
#endif /* (RSB_HAVE_CLAPACK && RSB_HAVE_CBLAS) */
	
		tmp_a = tmp_a-ddy[1-1];
		tmp_b = tmp_b-ddy[2-1];
		tmp_c = tmp_c-ddy[3-1];

		for(j=0;j<n;++j)
			z0[j]= tmp_a +tmp_b/(x[j]+tmp_c);

		error = 0;
		for(j=0;j<n;++j)
	       		error = error + (fabs( z0[j] - y[j] ) / y[j] );

		error = error / n * 100;
		if(error < min_err)
		{
		        *a=tmp_a;
		        *b=tmp_b;
		        *c=tmp_c;
		}
	}

	if((*c< 0) && (*c  < c_s))
	{
		*c=10000;
		*b=10000;
		avg=0;
		max=y[0];
		min=y[0];
		for(i=0;i<n;++i)
		{
		        if (y[i] > max) max=y[i];
		        if (y[i] < min) min=y[i];
		        avg=avg+y[i];
		}
		avg=avg/(double)(n);
		*a=avg;
		intl=max-min;
		avg=0;
		cnt=0;
		for(/*i=0*/;i<n;++i)
		//for(i=0;i<n;++i)
		{
        		if (fabs(y[i]-avg) < (0.3*intl))
			{
				avg = avg + y[i];
				cnt=cnt+1;
			}
		}
     		if(cnt > 0) *a=avg/(double)cnt;
	}
	else
  	if (*b >= 0)
	{
		*c=10000;
		*b=10000;
		avg=0;
		max=y[0];
		min=y[0];
		for(i=0;i<n;++i)
		{
			if (y[i] > max) max=y[i];
			if (y[i] < min) min=y[i];
			avg=avg+y[i];
		}
		avg=avg/(double)n;
		intl=max-min;
		avg=0;
		cnt=0;
		//for(i=0;i<n;++i)
		for(/*i=0*/;i<n;++i)
		{
		        if (fabs(y[i]-avg) < (0.3*intl))
			{
				avg = avg + y[i];
				cnt=cnt+1;
			}
		}
		if(cnt > 0) *a=avg/ (double) cnt;
	}


	RSB_STDOUT("\n");
	RSB_STDOUT("alpha:%lg beta:%lg gamma:%lg\n",*a,*b,*c);

	RSB_STDOUT("\nfitting:\n");
	for(i=0;i<n;++i)
	{
		RSB_STDOUT("%lg %lg\n", xcpy[i], *a+*b/(xcpy[i]+*c));
	}

	return RSB_ERR_NO_ERROR;
err:
	RSB_ERROR(RSB_ERRM_ES);
	RSB_DO_ERR_RETURN(errval)
#endif /* RSB_HAVE_CLAPACK && RSB_HAVE_CBLAS */
}
')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_REFERENCEBENCHMARK_FUNCTION_ARGS()
dnl	------------------------------------------------------------------
dnl
define(`RSB_M4_REFERENCEBENCHMARK_FUNCTION_ARGS',`dnl
dnl
`(void)'dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_REFERENCEBENCHMARK_FUNCTION_IDENTIFIER()
dnl	--------------------------------------------
dnl
define(`RSB_M4_REFERENCEBENCHMARK_FUNCTION_IDENTIFIER',`dnl
dnl
`rsb__do_referencebenchmark'dnl
dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_REFERENCEBENCHMARK_FUNCTION_NAME()
dnl	--------------------------------------
dnl
define(`RSB_M4_REFERENCEBENCHMARK_FUNCTION_NAME',`dnl
dnl
rsb_err_t RSB_M4_REFERENCEBENCHMARK_FUNCTION_IDENTIFIER`'dnl
dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_REFERENCEBENCHMARK_FUNCTION()
dnl	---------------------------------
dnl
define(`RSB_M4_REFERENCEBENCHMARK_FUNCTION',`dnl
dnl
RSB_M4_REFERENCEBENCHMARK_FUNCTION_NAME`'dnl
RSB_M4_REFERENCEBENCHMARK_FUNCTION_ARGS`'dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
{
	/*!
	 * \ingroup gr_bench
	 * Benchmark/test all supported matrix operations over all supported types.
         *
	 * \return \rsb_errval_inp_param_msg
	 */
	struct rsb_global_reference_performance_info_t grpi;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_blk_idx_t ri,ci;	/* row index, columns index */
	rsb_coo_idx_t order=20000;
	rsb_coo_idx_t rows=order,cols=order;	/* FIXME : TEMPORARY */
	rsb_blk_idx_t rua[] = RSB_ROWS_UNROLL_ARRAY;
	rsb_blk_idx_t cua[] = RSB_COLUMNS_UNROLL_ARRAY;
	double tot_secs=0.0,pred_secs=1.0;
	rsb_trans_t transA = RSB_DEFAULT_TRANSPOSITION;
	size_t kernels_n = RSB_ROWS_UNROLL_ARRAY_LENGTH*RSB_COLUMNS_UNROLL_ARRAY_LENGTH*RSB_IMPLEMENTED_MOPS*RSB_IMPLEMENTED_TYPES;
	rsb_int ti=0;	/* type index */
	int fbw,bwi;
	const rsb_time_t mrbt = rsb__getenv_real_t("RSB_BENCHMARK_MIN_SECONDS", RSB_BENCHMARK_MIN_SECONDS);
	RSB_BZERO_P(&grpi);

	/* if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))) RSB_PERR_GOTO(err,RSB_ERRM_ES); we skip this to enable calling this from within our library */

	if(RSB_FITTING_SAMPLES<2)
	{	
		fbw=(RSB_FIRST_FITTING_SAMPLE_BW_MAX + RSB_FIRST_FITTING_SAMPLE_BW_MIN)/2;
		bwi=fbw;
	}
	else
	{
		fbw = RSB_FIRST_FITTING_SAMPLE_BW_MIN;
		bwi=(RSB_FIRST_FITTING_SAMPLE_BW_MAX - RSB_FIRST_FITTING_SAMPLE_BW_MIN)/(RSB_FITTING_SAMPLES-1);
	}
	
	tot_secs = -rsb_time();
	pred_secs *= RSB_ROWS_UNROLL_ARRAY_LENGTH * RSB_COLUMNS_UNROLL_ARRAY_LENGTH * RSB_FITTING_SAMPLES * RSB_IMPLEMENTED_META_MOPS *  RSB_IMPLEMENTED_TYPES * mrbt;
	RSB_STDERR("#reference benchmarking of %zd kernels (no transposed, no symmetric, and so on) should take at least %lg seconds..\n",kernels_n,pred_secs);

foreach(`mtype',RSB_M4_MATRIX_TYPES,`dnl
	/* mtype type benchmarking */
/*	RSB_INFO("#mtype type benchmarking\n");*/
	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		{
			rsb_blk_idx_t br = rua[ri];
			rsb_blk_idx_t bc = cua[ci];
			rsb_coo_idx_t bw,mbw=(cols/bc);
			rsb_int si=0;	/* sample index */
			mbw=(cols-bc)/bc;	/* tune here to fill further our matrix */
			/* FIXME : there is the danger of empty samples! */
			for(bw=fbw;bw<=mbw && si< RSB_FITTING_SAMPLES ;bw+=bwi)	/* this parameter should be tunable, too */
			{
				//RSB_INFO("bw = %d\n",bw);
				rsb_int moi=0;	/* matrix operation index */
				double time,*timep=&time;
				struct rsb_mtx_t * mtxAp =
					rsb__generate_blocked_banded(br,bc,rows,cols,bw,timep,RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype),RSB_BOOL_TRUE );	/* FIXME : generating triangular factors always ! */
				if(!mtxAp)
				{
					RSB_STDERR(RSB_ERRM_IE);
					{errval = RSB_ERR_GENERIC_ERROR; RSB_PERR_GOTO(err,RSB_ERRM_ES); }
				}
dnl				struct rsb_options_t * o = mtxAp->options;

foreach(`mop',RSB_M4_MATRIX_META_OPS,`dnl
				{
/*					RSB_INFO("#mtype type, ");*/
/*					RSB_INFO("mop operation benchmarking\n");*/
					/* mop operation benchmarking */
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
					mtype *out=NULL,*rhs=NULL;
')dnl
ifelse(RSB_M4_IS_SPXX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
			double alpha=1.0;/* FIXME */
			double * alphap = &alpha;
')dnl
ifelse(RSB_M4_IS_SPXX_SCALING_KERNEL_MOP(mop),1,`dnl
			double beta =1.0;/* FIXME */
			double * betap  = &beta ;
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
					mtype * row_sums;
')dnl

					
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
					row_sums = rsb__malloc(mtxAp->el_size*(rows+br));
					if(!row_sums) {RSB_ERROR(RSB_ERRM_ES);errval = RSB_ERR_ENOMEM;goto erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop;}
					if(rsb__fill_with_ones(row_sums,mtxAp->typecode,cols,1))     {RSB_ERROR(RSB_ERRM_ES);errval = RSB_ERR_ENOMEM;goto erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop;}
')dnl
ifelse(mop,`scale',`dnl
					mtype * scale_factors = rsb__malloc(mtxAp->el_size*(rows+br));
					if(!scale_factors) {RSB_ERROR(RSB_ERRM_ES);errval = RSB_ERR_ENOMEM;goto erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop;}
					if(rsb__fill_with_ones(scale_factors,mtxAp->typecode,rows,1))     {RSB_ERROR(RSB_ERRM_ES);errval = RSB_ERR_ENOMEM;goto erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop;}
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
')dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
					rsb_coo_idx_t nrhs=4;
					rsb_coo_idx_t bstride = cols+bc;
					rsb_coo_idx_t cstride = rows+br;
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
					rsb_coo_idx_t incx=1,incy=1;
					incx=1,incy=1;	/* this is just a pacifier for "unused variable"-like warnings */
',`dnl
dnl					rsb_coo_idx_t incx=1,incy=1;
')dnl
					rhs = rsb__malloc(mtxAp->el_size*(bstride)*nrhs);
					out = rsb__malloc(mtxAp->el_size*(cstride)*nrhs);
					if(!out || rsb__fill_with_ones(out,mtxAp->typecode,cstride*nrhs,1)){RSB_ERROR(RSB_ERRM_ES);errval = RSB_ERR_ENOMEM;goto erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop;}
					if(!rhs || rsb__fill_with_ones(rhs,mtxAp->typecode,bstride*nrhs,1)){RSB_ERROR(RSB_ERRM_ES);errval = RSB_ERR_ENOMEM;goto erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop;}
')dnl
ifelse(mop,`negation',`dnl
					int please_fix_RSB_M4_ARGS_TO_ACTUAL_ARGS=-1;/* here to fix negation */
')dnl

					grpi.gpi[ti].pipmo[moi].blocks_per_row[si]=bw*bc; /* FIXME : TEMPORARY !!  */

					/* we benchmark our mtype library implementation for operation mop */
					grpi.gpi[ti].pipmo[moi].pipfs[si].seconds[ri][ci] = mrbt; /* min seconds */
					grpi.gpi[ti].pipmo[moi].pipfs[si].m_flops[ri][ci] = (double)RSB_BENCHMARK_MIN_RUNS; /* min runs */

					errval = dnl
ifelse(RSB_M4_MATRIX_OP_IS_META_OP(mop),`1',dnl
`0;/* meta-op : we already measured matrix creation time  */
grpi.gpi[ti].pipmo[moi].pipfs[si].seconds[ri][ci]=time;
grpi.gpi[ti].pipmo[moi].pipfs[si].m_flops[ri][ci]=((double)rsb__do_get_matrix_nnz(mtxAp))/1000000;
/* FIXME : this is experimental and unfinished code */
',`
RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_IDENTIFIER(mop,mtype)(dnl 
&(grpi.gpi[ti].pipmo[moi].pipfs[si].seconds[ri][ci]),dnl
&(grpi.gpi[ti].pipmo[moi].pipfs[si].m_flops[ri][ci]),dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ACTUAL_ARGS(mop,mtype));')

					grpi.gpi[ti].pipmo[moi].pipfs[si].fillin[ri][ci]  = rsb__do_get_matrix_fillin(mtxAp); 
					grpi.gpi[ti].pipmo[moi].pipfs[si].rows = rows;
					grpi.gpi[ti].pipmo[moi].pipfs[si].cols = cols;
					grpi.gpi[ti].pipmo[moi].pipfs[si].nnz  = rsb__do_get_matrix_nnz(mtxAp) ;
					grpi.gpi[ti].pipmo[moi].pipfs[si].flags= mtxAp->flags ;
					grpi.gpi[ti].pipmo[moi].pipfs[si].storage= mtxAp->matrix_storage ;
					grpi.gpi[ti].pipmo[moi].pipfs[si].typecode= mtxAp->typecode ;
					grpi.gpi[ti].pipmo[moi].pipfs[si].element_count= mtxAp->element_count;

					grpi.gpi[ti].pipmo[moi].pipfs[si].e_mflops[ri][ci] = 
						grpi.gpi[ti].pipmo[moi].pipfs[si].m_flops[ri][ci] /
						grpi.gpi[ti].pipmo[moi].pipfs[si].fillin[ri][ci];

					if(RSB_SOME_ERROR(errval)){RSB_ERROR(RSB_ERRM_ES);goto erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop;}
					++moi;

					erri_`'RSB_M4_CHOPSPACES(mtype)`'`_'`'mop:
					if(RSB_SOME_ERROR(errval))
						RSB_PERR_GOTO(err,RSB_ERRM_ES);

					RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
					RSB_CONDITIONAL_FREE(row_sums);
')dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
					RSB_CONDITIONAL_FREE(out);
					RSB_CONDITIONAL_FREE(rhs);
')dnl
ifelse(mop,`scale',`dnl
					RSB_CONDITIONAL_FREE(scale_factors);
')dnl
				}
')dnl
				RSB_MTX_FREE(mtxAp);
				++si;
			}	
		}
	}
	{
		rsb_int moi;
		rsb_char_t * mops[] = RSB_M4_MATRIX_META_OPS_ARRAY;
		rsb_char_t * types[] = RSB_M4_MATRIX_TYPES_ARRAY;
		rsb_char_t s[RSB_M4_BUFLEN];
		rsb__print_mop_reference_performance_info_header();
		for(moi=0;moi<RSB_IMPLEMENTED_META_MOPS;++moi)
		{	
/*			rsb_int si;*/
			/* informational printout */
			sprintf(s,"%s\t%s\t",types[ti], mops[moi]);
			rsb__print_mop_reference_performance_info(&(grpi.gpi[ti].pipmo[moi]),s);
/*			for(si=0;si<RSB_FITTING_SAMPLES;++si)*/
/*				rsb__dump_performance_info(&(grpi.gpi[ti].pipmo[moi].pipfs[si]), NULL);*/
		}
	}
	++ti;
')dnl
	tot_secs += rsb_time();
	RSB_STDERR("#reference benchmarking took %lg seconds (predicted %lg :)....\n",tot_secs,pred_secs);

	grpi.initialized=1;	/* FIXME : only partially */
	//rsb__dump_global_reference_performance_info(&grpi);
#if RSB_WANT_PERFORMANCE_FILE
	rsb__save_global_reference_performance_info(&grpi);
#endif /* RSB_WANT_PERFORMANCE_FILE */

	ti=0;	/* type index */
	for(ti=0;ti<RSB_IMPLEMENTED_TYPES	;++ti)
	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		{
			rsb_blk_idx_t bc = cua[ci];
			rsb_int moi=0;	/* matrix operation index */
			for(moi=0;moi<RSB_IMPLEMENTED_META_MOPS ;++moi)
			{
				rsb_int si=0;	/* sample index */

				double y[RSB_FITTING_SAMPLES];
				double * x = grpi.gpi[ti].pipmo[moi].blocks_per_row;

				for(si=0;si< RSB_FITTING_SAMPLES ;++si)
				{
					/* we tune our mtype library implementation for operation mop */
						y[si] = 
							grpi.gpi[ti].pipmo[moi].pipfs[si].m_flops[ri][ci]/
							grpi.gpi[ti].pipmo[moi].pipfs[si].seconds[ri][ci];
				}

				/*
				 * FIXME : make this fitting analysis offline respect our benchmark!
				 */
				errval = RSB_M4_HYPERBOLIC_FITTING_FUNCTION_IDENTIFIER()(
						x, y, 3, 
						&(grpi.gpi[ti].pipmo[moi].alpha[ri][ci]),
						&(grpi.gpi[ti].pipmo[moi].beta [ri][ci]),
						&(grpi.gpi[ti].pipmo[moi].gamma[ri][ci]), (double)bc
						/* FIXME : is this right ?*/
					);
				if(RSB_SOME_ERROR(errval))
				{
					if(errval==RSB_ERR_UNSUPPORTED_OPERATION)
						; /* not a problem: this model is obsolete */
						/* RSB_ERROR(RSB_ERRM_UNSUPPORTED_OPERATION); */
					else
						RSB_PERR_GOTO(err,RSB_ERRM_ES);
				}

			}
		}
	}

	errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
	if( RSB_SOME_ERROR(errval) )
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err,RSB_ERRM_ES);
	}
err:
	RSB_DO_ERR_RETURN(errval)
}
')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
RSB_M4_HYPERBOLIC_FITTING_FUNCTION()
RSB_M4_REFERENCEBENCHMARK_FUNCTION()
dnl
dnl
dnl
dnl
dnl
dnl
dnl
#ifdef __cplusplus
}
#endif  /* __cplusplus */
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_BENCH_H_INCLUDED */
')
/* @endcond */
dnl
