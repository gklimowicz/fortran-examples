dnl
dnl
dnl	@author: Michele Martone
dnl
dnl	The generated code will expand from here
dnl	Take care of compiling this code without loop unrolling optimizations (-fno-unroll-loops, or -ON with N<=2 on gcc)
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
dnl
dnl
dnl
dnl	FIXME : THERE ARE NO TRANSPOSED KERNELS
dnl
dnl
dnl
dnl	RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME(mtype,matrix_storage,transposition,symmetry,b_rows,b_columns,unrolling,mop) dnl	-----------------------------------------------------------------------------------------------
dnl
define(`RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`transposition',$3)dnl
pushdef(`symmetry',$4)dnl
pushdef(`b_rows',$5)dnl		block rows
pushdef(`b_columns',$6)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t')dnl integer type (for indices)
pushdef(`unrolling',$7)dnl	
pushdef(`mop',$8)dnl	
pushdef(`citype',$9)dnl	
dnl
dnl pushdef(`diagonal',ifelse($10,`',`e',diagonal,uplo))dnl	FIXME: for now this prameter is optional
pushdef(`diagonal',$10)dnl	FIXME: new
pushdef(`uplo',$11)dnl
dnl
RSB_M4_PREFIX`'matrix_storage`_'mop`_'RSB_M4_TYPE_CODE(mtype)`_'dnl
RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE(citype)`_'dnl
ifelse(matrix_storage,`fixed_block',`_r'b_rows`_c'b_columns)`'dnl
`_t'touppercase(RSB_M4_TRANSPOSITION_CODE(transposition))`'dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`_r'b_rows`_c'b_columns)`'dnl
ifelse(RSB_M4_IS_FORMAT_BCOO(matrix_storage),1,`_r'b_rows`_c'b_columns)`'dnl
dnl ifelse(matrix_storage,`VBR',`_r'b_rows`_c'b_columns)`'dnl
dnl ifelse(matrix_storage,`VBC',`_r'b_rows`_c'b_columns)`'dnl
`_u'unrolling`_s'touppercase(symmetry)`'dnl
`_d'touppercase(diagonal)`'dnl
`_u'touppercase(uplo)`'dnl
popdef(`uplo')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`symmetry')dnl
popdef(`transposition')dnl
popdef(`b_rows')dnl
popdef(`b_columns')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
popdef(`diagonal')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_ARGS(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
dnl	-----------------------------------------------------------------------------------------------
dnl
define(`RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_ARGS',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`b_rows',$3)dnl		block rows
pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t')dnl integer type (for indices)
pushdef(`unrolling',$5)dnl	
pushdef(`mop',$6)dnl	
pushdef(`matrix_structs',`const itype Mdim, const itype mdim, const rsb_nnz_idx_t * RSB_M4_RESTRICT bindx, const rsb_nnz_idx_t * RSB_M4_RESTRICT bpntr, const rsb_nnz_idx_t *RSB_M4_RESTRICT indptr, const rsb_coo_idx_t * RSB_M4_RESTRICT rpntr, const rsb_coo_idx_t * RSB_M4_RESTRICT cpntr, const rsb_coo_idx_t dummy_br, const rsb_coo_idx_t dummy_bc')dnl	
(
ifelse(mop,`scale',`dnl
	mtype * VA, 
	matrix_structs, 
	const mtype *scale_factors
')dnl
ifelse(RSB_M4_MEMBER(mop,`spmv_uauz',`spmv_uaua',`spmv_unua'),1,`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs
')dnl
ifelse(mop,`spmm_az',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs,
	const itype bstride, const itype cstride, const itype nrhs
')dnl
ifelse(mop,`spmv_uxux',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs,
	const mtype * alphap, const mtype * betap
')dnl
ifelse(mop,`infty_norm',`dnl
ifelse(matrix_storage,`_fixed_block',`dnl
	const mtype * VA, mtype * local_row_sums,
	matrix_structs
',`dnl
	const mtype * VA, mtype * global_row_sums, 
	matrix_structs
')dnl
')dnl
ifelse(mop,`negation',`dnl
ifelse(matrix_storage,`_fixed_block',`dnl
	mtype * VA, 
	matrix_structs
',`dnl
	mtype * VA, 
	matrix_structs
')dnl
')dnl
)dnl
popdef(`matrix_structs')dnl	
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`b_rows')dnl
popdef(`b_columns')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_PROTOTYPE(mtype,matrix_storage,transposition,b_rows,b_columns,unrolling,mop)
dnl	----------------------------------------------------------------------------------------------------
dnl
define(`RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_PROTOTYPE',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`transposition',$3)dnl
pushdef(`b_rows',$4)dnl		block rows
pushdef(`b_columns',$5)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t')dnl integer type (for indices)
pushdef(`unrolling',$6)dnl	
pushdef(`mop',$7)dnl	
rsb_err_t RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME(mtype,matrix_storage,transposition,symmetry,b_rows,b_columns,unrolling,mop)dnl
RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_ARGS(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)dnl
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`transposition')dnl
popdef(`b_rows')dnl
popdef(`b_columns')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_BODY(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
dnl	-----------------------------------------------------------------------------------------------
dnl
define(`RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_BODY',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`b_rows',$3)dnl		block rows
pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t')dnl integer type (for indices)
pushdef(`unrolling',$5)dnl
pushdef(`mop',$6)dnl
{
ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
pushdef(`out_dim',rpntr[Mdim])dnl
pushdef(`mi',`i')dnl
pushdef(`Mi',`j')dnl
')dnl
ifelse(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),1,`dnl
pushdef(`out_dim',rpntr[mdim])dnl
pushdef(`mi',`j')dnl
pushdef(`Mi',`i')dnl
')dnl
dnl

	/**
	 * \ingroup rsb_doc_kernels
ifelse(matrix_storage,`VBR',`dnl
	 * This code is intended for a pure VBR partitioned matrix.
	 * It does not dispatch the kernel function for each block,
	 * but employ explicitly inlined kernels.
')dnl
	 *
	 * \return \rsb_errval_inp_param_msg
	 */
	register rsb_coo_idx_t i,j;
	register rsb_nnz_idx_t k;
	rsb_coo_idx_t columns;
	rsb_coo_idx_t rows   ;

ifelse(mop,`spmv_uxux',`dnl
	/* this is slow, however this is column based scanning.  FIXME : optimize spmv_uxux */
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type), out_dim, betap, out, 1);/* we scale the destination vector */
')dnl
ifelse(mop,`spmv_uauz',`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type), out_dim, NULL, out, 1);// FIXME: stride is bad
')dnl

	for(Mi=0;Mi<Mdim;++Mi)
	{
		/* Mi is the working block row */
		if(bpntr[Mi]==bpntr[Mi+1])continue;/* empty block row (or column) */
		for(k=bpntr[Mi];k<bpntr[Mi+1];++k)	/* k is the index of the block */
		{
			/* mi is the working block column (or row) */
			mi=bindx[k];
ifelse(matrix_storage,`VBR',`dnl
			columns=cpntr[j+1]-cpntr[j];
			rows   =rpntr[i+1]-rpntr[i];
')dnl
ifelse(1,0,`dnl
ifelse(RSB_M4_MEMBER(matrix_storage,`VBR',`VBC'),1,`dnl
/* matrix_storage ! mop  */
			DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_IDENTIFIER(mop,mtype,unrolling)(rows,columns)
ifelse(mop,`scale',`dnl
			(
					VA+indptr[(k)],
					scale_factors+rpntr[i],
					rows, columns
					);
')dnl
ifelse(mop,`spmv_uauz',`dnl
			(
					VA+indptr[(k)],
					rhs+cpntr[j],
					out+rpntr[i],
					rows, columns
					);
')dnl
ifelse(mop,`spmv_unua',`dnl
			(
					VA+indptr[(k)],
					rhs+cpntr[j],
					out+rpntr[i],
					rows, columns
					);
')dnl
ifelse(mop,`spmv_uaua',`dnl
			(
					VA+indptr[(k)],
					rhs+cpntr[j],
					out+rpntr[i],
					rows, columns
					);
')dnl
ifelse(mop,`spmm_az',`dnl
			(
					VA+indptr[(k)],
					rhs+cpntr[j],
					out+rpntr[i],
					rows, columns,
					bstride, cstride, nrhs
					);
')dnl
ifelse(mop,`spmv_uxux',`dnl
			(
					VA+indptr[(k)],
					rhs+cpntr[j],
					out+rpntr[i],
					rows, columns,
					alphap, betap
					);
')dnl
ifelse(mop,`infty_norm',`dnl
			(
					VA+indptr[(k)],
					global_row_sums+cpntr[j],
					rows, columns
					);
')dnl
ifelse(mop,`negation',`dnl
			(
					VA+indptr[(k)],
					rows, columns
					);
')dnl
')dnl
')dnl
ifelse(RSB_M4_MEMBER(matrix_storage,`VBR',`VBC'),1,`dnl
ifelse(mop,`scale',`dnl
			mtype *a = VA+indptr[(k)];
			const mtype *d = scale_factors+rpntr[i];
')dnl
ifelse(RSB_M4_MEMBER(mop,`spmv_uaua',`spmv_unua',`spmv_uauz',`spmm_az',`spmv_uxux'),1,`dnl
			const mtype *a = VA+indptr[(k)];
			const mtype *b = rhs+cpntr[j];
			mtype *c = out+rpntr[i];
')dnl
ifelse(mop,`infty_norm',`dnl
			const mtype *a = VA+indptr[(k)];
			mtype *local_row_sums = global_row_sums+rpntr[i];
')dnl
ifelse(mop,`negation',`dnl
			mtype *a = VA+indptr[(k)];
')dnl
			columns=cpntr[j+1]-cpntr[j];
			rows   =rpntr[i+1]-rpntr[i];
			/* we jump to the right point */
dnl
ifelse(1,0,`
#if 0
/* it will not happen */
#ifdef RSB_WANT_BLOCK_TRAILING_STRUCT	/* EXPERIMENTAL */
ifelse(mop,`scale',`dnl
			a = (mtype*) (((char*)a) + (RSB_BLOCK_EXTRA_BYTES)*(k+1)) ;
')dnl
ifelse(mop,`spmv_uauz',`dnl
			a = (const mtype*) (((const char*)a) + (RSB_BLOCK_EXTRA_BYTES)*(k+1)) ;
')dnl
ifelse(mop,`spmm_az',`dnl
			a = (const mtype*) (((const char*)a) + (RSB_BLOCK_EXTRA_BYTES)*(k+1)) ;
')dnl
ifelse(mop,`spmv_uxux',`dnl
			a = (const mtype*) (((const char*)a) + (RSB_BLOCK_EXTRA_BYTES)*(k+1)) ;
')dnl
ifelse(mop,`infty_norm',`dnl
			a = (const mtype*) (((const char*)a) + (RSB_BLOCK_EXTRA_BYTES)*(k+1)) ;
')dnl
ifelse(mop,`negation',`dnl
			a = (mtype*) (((char*)a) + (RSB_BLOCK_EXTRA_BYTES)*(k+1)) ;
')dnl
dnl
#endif /* RSB_WANT_BLOCK_TRAILING_STRUCT */
#endif /* 0 */
')dnl

switch(rows)
{
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
	case rowsu:
	switch(columns)
	{
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
		case colsu:
dnl			goto RSB_M4_KERNEL_FUNCTION_NAME(mtype,rowsu,colsu,`',mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE);
{RSB_M4_KERNEL_FUNCTION_BODY(`row',`rows',rowsu,`column',`columns',colsu,mtype,,mop,`')}
			break;
')dnl
	default: goto RSB_M4_KERNEL_FUNCTION_NAME(mtype,RSB_M4_ROWS_FALLBACK_UNROLL,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE);
	}
	break;
')dnl
	default: goto RSB_M4_KERNEL_FUNCTION_NAME(mtype,RSB_M4_ROWS_FALLBACK_UNROLL,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE);
	RSB_M4_KERNEL_FUNCTION_NAME(mtype,RSB_M4_ROWS_FALLBACK_UNROLL,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE):
{RSB_M4_KERNEL_FUNCTION_BODY(`row',`rows',RSB_M4_ROWS_FALLBACK_UNROLL,`column',`columns',RSB_M4_COLUMNS_FALLBACK_UNROLL,mtype,,mop,`l')}
};
')dnl end VBR
		}
	}
popdef(`mi')dnl
popdef(`Mi')dnl
dnl
dnl
dnl
popdef(`out_dim')dnl
	return RSB_ERR_NO_ERROR;
}
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`b_rows')dnl
popdef(`b_columns')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION(mtype,matrix_storage,transposition,b_rows,b_columns,unrolling,mop)
dnl	------------------------------------------------------------------------------------------
dnl
dnl	These functions will perform their operations on fixed block matrices.
dnl
define(`RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`transposition',$3)dnl
pushdef(`b_rows',$4)dnl		block rows
pushdef(`b_columns',$5)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t')dnl integer type (for indices)
pushdef(`unrolling',$6)dnl	
pushdef(`mop',$7)dnl	
ifelse(matrix_storage,`fixed_block',dnl
`/* This code is intended for a purely blocked matrix. */',dnl
`/* This code is intended for a pure VBR partitioned matrix. */')
RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_PROTOTYPE(mtype,matrix_storage,transposition,b_rows,b_columns,unrolling,mop)dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_BODY(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
popdef(`transposition')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
