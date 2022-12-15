dnl
dnl	@author: Michele Martone
dnl
dnl
dnl	RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_ARGS(mtype,matrix_storage,unrolling,mop)
dnl	--------------------------------------------------------------------------------
dnl
dnl	FIXME : THERE ARE NO TRANSPOSED KERNELS
dnl	FIXME : THERE ARE NO SYMMETRY HANDLING KERNELS
dnl
define(`RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_ARGS',`dnl
dnl
include(`do_unroll.m4')dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl	
pushdef(`unrolling',$3)dnl	
pushdef(`mop',$4)dnl	
pushdef(`matrix_structs',`const itype Mdim, const itype mdim, const rsb_nnz_idx_t * RSB_M4_RESTRICT bindx, const rsb_nnz_idx_t * RSB_M4_RESTRICT bpntr, const rsb_nnz_idx_t *RSB_M4_RESTRICT indptr, const rsb_coo_idx_t * RSB_M4_RESTRICT rpntr, const rsb_coo_idx_t * RSB_M4_RESTRICT cpntr, const rsb_coo_idx_t dummy_br, const rsb_coo_idx_t dummy_bc')dnl	
(
ifelse(mop,`scale',`dnl
	mtype * VA, 
	matrix_structs,
	const mtype *scale_factors
')dnl
ifelse(mop,`spmv_uauz',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs
')dnl
ifelse(mop,`spmv_uaua',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs
')dnl
ifelse(mop,`spmv_unua',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs
')dnl
ifelse(mop,`spmv_uxux',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs,
	const mtype * alphap, const mtype * betap
')dnl
ifelse(mop,`spmm_az',`dnl
	const mtype * VA, const mtype * mrhs, mtype * mout,
	matrix_structs,
	const itype bstride, 
	const itype cstride,
	const itype nrhs
')dnl
ifelse(mop,`infty_norm',`dnl
	const mtype * VA, mtype * global_row_sums,/* ! */
	matrix_structs
')dnl
ifelse(mop,`negation',`dnl
	mtype * VA, 
	matrix_structs
')dnl
)dnl
popdef(`matrix_structs')dnl
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl	RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_PROTOTYPE
dnl	-------------------------------------------------
dnl
define(`RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_PROTOTYPE',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`unrolling',$3)dnl
pushdef(`mop',$4)dnl
pushdef(`citype',$5)dnl
pushdef(`diagonal',RSB_M4_DEFAULT_DIAGONAL_TYPE)dnl
rsb_err_t RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,symmetry,unrolling,mop,citype,diagonal,uplo)dnl
RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_ARGS(mtype,matrix_storage,unrolling,mop,citype)dnl
popdef(`diagonal')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_BODY(mtype,matrix_storage,unrolling,mop)
dnl	--------------------------------------------------------------------------------
dnl
define(`RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_BODY',`
dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`unrolling',$3)dnl
pushdef(`mop',$4)dnl
pushdef(`citype',`rsb_coo_idx_t')dnl
{
RSB_M4_DEBUGINFO(``$0'')dnl
	/*!
	 * This function is experimental.
	 */

ifelse(RSB_M4_IS_FORMAT_LINKED_LIST(matrix_storage),1,`dnl
dnl	FIXME: the following call uses non declared arguments
pushdef(`args',`RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo))')dnl
	/* this is sample code for scanning a whole linked list format matrix */
	struct rsb_block_tail_t * bt;
	const char *data = (const char*) VA;


ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
pushdef(`out_dim',rpntr[Mdim])dnl
',`dnl
pushdef(`out_dim',rpntr[mdim])dnl
')dnl
ifelse(mop,`spmv_uxux',`dnl
	/* this is slow, however this is column based scanning.  FIXME : optimize spmv_uxux */
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type), out_dim, betap, out, 1);/* we scale the destination vector */
')dnl
ifelse(mop,`spmv_uauz',`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),out_dim,NULL,out,incy);
')dnl

	bt = (struct rsb_block_tail_t*)data;
	while(RSB_LIKELY(bt->block_rows))
	{
		data += sizeof(struct rsb_block_tail_t);
		/* do stuff */
/*		RSB_STDERR("%d %d ",bt->block_columns , bt->block_rows);
		RSB_STDERR("%d %d ",bt->block_column , bt->block_row);
		RSB_STDERR("%d %d ",bt->base_column , bt->base_row);
		RSB_STDERR(" : %d\n",*(int*)data);*/

		//int j = bt->block_column;
		//int i = bt->block_row;
		register rsb_coo_idx_t columns = bt->block_columns;
		register rsb_coo_idx_t rows = bt->block_rows;
		//int colsu = bt->block_columns;
		//int rowsu = bt->block_rows;
		type * a = (type*) data;

ifelse(mop,`scale',`dnl
			const mtype *d = scale_factors+bt->base_row;
')dnl
ifelse(mop,`spmv_uauz',`dnl
			const mtype *b = rhs+bt->base_column;
dnl			/*mtype *c = out+(rowsu*i);*/ /* experimentally commented out and put up */
			mtype *c = out+bt->base_row;
')dnl
ifelse(mop,`spmv_uaua',`dnl
			const mtype *b = rhs+bt->base_column;
dnl			/*mtype *c = out+(rowsu*i);*/ /* experimentally commented out and put up */
			mtype *c = out+bt->base_row;
')dnl
ifelse(mop,`spmv_unua',`dnl
			const mtype *b = rhs+bt->base_column;
dnl			/*mtype *c = out+(rowsu*i);*/ /* experimentally commented out and put up */
			mtype *c = out+bt->base_row;
')dnl
ifelse(mop,`spmv_uxux',`dnl
			/*const mtype *b = rhs+(colsu*j);*/
			const mtype *b = rhs+bt->base_column;
			/*mtype *c = out+(rowsu*i);*/
			mtype *c = out+bt->base_row;
')dnl
ifelse(mop,`spmm_az',`dnl
			/*const mtype *b = mrhs+(colsu*j);
			mtype *c = mout+(rowsu*i);*/
			const mtype *b = mrhs+bt->base_column;
			mtype *c = mout+bt->base_row;
')dnl
ifelse(mop,`infty_norm',`dnl
			/*mtype *row_sums = global_row_sums+(rowsu*i);*/
			mtype *row_sums = global_row_sums+bt->base_row;
')dnl
ifelse(mop,`negation',`dnl
')dnl
			/* `mop' is mop */
dnl {RSB_M4_KERNEL_FUNCTION_BODY(`row',`rows',b_rows,`column',`columns',b_columns,mtype,,mop,unrolling)}

	DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_IDENTIFIER(mop,mtype,unrolling)(rows,columns)
ifelse(mop,`scale',`dnl
			(
					a,
					scale_factors+bt->base_row,
					rows, columns
					);
')dnl
ifelse(mop,`spmv_uauz',`dnl
			(
					a,b,c,
					rows, columns
					);
')dnl
ifelse(mop,`spmv_unua',`dnl
			(
					a,b,c,
					rows, columns
					);
')dnl
ifelse(mop,`spmv_uaua',`dnl
			(
					a,b,c,
					rows, columns
					);
')dnl
ifelse(mop,`spmm_az',`dnl
			(
					a,
					mrhs+bt->base_column,
					mout+bt->base_row,
					rows, columns,
					bstride, cstride, nrhs
					);
')dnl
ifelse(mop,`spmv_uxux',`dnl
			(
					a,
					rhs+bt->base_column,
					out+bt->base_row,
					rows, columns,
					alphap, betap
					);
')dnl
ifelse(mop,`infty_norm',`dnl
			(
					a,
					global_row_sums+bt->base_row,
					rows, columns
					);
')dnl
ifelse(mop,`negation',`dnl
			(
					a,
					rows, columns
					);
')dnl


		data += /*(size_t)*/(int)(bt->block_columns * bt->block_rows)*(int)sizeof(type);
		bt = (struct rsb_block_tail_t*)data;
	}

dnl

popdef(`args')dnl
')dnl
	return RSB_ERR_NO_ERROR;
popdef(`out_dim')dnl
}
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
dnl	These functions dispatch on the column size, calling the
dnl	proper kernels.
dnl
dnl	They assume type dispatching has just been performed.
dnl
dnl
dnl	RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION(mtype,matrix_storage,unrolling,mop)
dnl	---------------------------------------------------------------------------
dnl
define(`RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION',`
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`unrolling',$3)dnl	
dnl pushdef(`b_rows',$3)dnl		block rows
dnl pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t')dnl integer type (for indices)
pushdef(`mop',`$4')dnl
`/* This code is intended for a block compressed sparse stripe matrix. */'
RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_PROTOTYPE(mtype,matrix_storage,unrolling,mop)dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`dnl
RSB_M4_LL_KERNEL_SIZE_DISPATCH_FUNCTION_BODY(mtype,matrix_storage,unrolling,mop)
')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl	RSB_M4_LL_KERNEL_FUNCTION_NAME(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
dnl	-----------------------------------------------------------------------------------
dnl
define(`RSB_M4_LL_KERNEL_FUNCTION_NAME',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`b_rows',$3)dnl		block rows
pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`int')dnl integer type (for indices)
pushdef(`unrolling',$5)dnl	
pushdef(`mop',$6)dnl	
RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME(mtype,matrix_storage,transposition,symmetry,b_rows,b_columns,unrolling,mop) dnl
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
dnl	RSB_M4_LL_KERNEL_FUNCTION_ARGS(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
dnl	-----------------------------------------------------------------------------------
dnl
define(`RSB_M4_LL_KERNEL_FUNCTION_ARGS',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`b_rows',$3)dnl		block rows
pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`int')dnl integer type (for indices)
pushdef(`unrolling',$5)dnl	
pushdef(`mop',$6)dnl	
pushdef(`matrix_structs',`const itype Mdim, const itype mdim, const rsb_nnz_idx_t * bindx, const rsb_nnz_idx_t * bpntr, const rsb_nnz_idx_t *indptr, const rsb_coo_idx_t * rpntr, const rsb_coo_idx_t * cpntr')dnl	
(
ifelse(mop,`scale',`dnl
	mtype * VA, 
	matrix_structs,
	const mtype *scale_factors
')dnl
ifelse(mop,`spmv_uauz',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs
')dnl
ifelse(mop,`spmv_uxux',`dnl
	const mtype * VA, const mtype * rhs, mtype * out,
	matrix_structs,
	const mtype * alphap, const mtype * betap
')dnl
ifelse(mop,`spmm_az',`dnl
	const mtype * VA, const mtype * mrhs, mtype * mout,
	matrix_structs,
	const itype bstride, 
	const itype cstride,
	const itype nrhs 
')dnl
ifelse(mop,`infty_norm',`dnl
	const mtype * VA, mtype * global_row_sums, 
	matrix_structs
')dnl
ifelse(mop,`negation',`dnl
	mtype * VA, 
	matrix_structs
')dnl
dnl
dnl
dnl
dnl
)
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
dnl	RSB_M4_LL_KERNEL_FUNCTION_PROTOTYPE(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
dnl	----------------------------------------------------------------------------------------
dnl
define(`RSB_M4_LL_KERNEL_FUNCTION_PROTOTYPE',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`b_rows',$3)dnl		block rows
pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`int')dnl integer type (for indices)
pushdef(`unrolling',$5)dnl	
pushdef(`mop',$6)dnl	
ifelse(RSB_M4_IS_FORMAT_LINKED_LIST(matrix_storage),1,`dnl
rsb_err_t RSB_M4_LL_KERNEL_FUNCTION_NAME(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)dnl
RSB_M4_LL_KERNEL_FUNCTION_ARGS(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)dnl
')dnl
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
dnl	RSB_M4_LL_KERNEL_FUNCTION_BODY(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
dnl	-----------------------------------------------------------------------------------
dnl
define(`RSB_M4_LL_KERNEL_FUNCTION_BODY',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl
pushdef(`b_rows',$3)dnl		block rows
pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`int')dnl integer type (for indices)
pushdef(`unrolling',$5)dnl
pushdef(`mop',$6)dnl
{
RSB_M4_DEBUGINFO(``$0'')dnl
	/* FIXME : STUB */
dnl
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
dnl	RSB_M4_LL_KERNEL_FUNCTION(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)
dnl	------------------------------------------------------------------------------
dnl
dnl	These functions will perform their operations on fixed block matrices.
dnl
define(`RSB_M4_LL_KERNEL_FUNCTION',`dnl
dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl	
pushdef(`b_rows',$3)dnl		block rows
pushdef(`b_columns',$4)dnl	block columns
pushdef(`itype',`int')dnl integer type (for indices)
pushdef(`unrolling',$5)dnl	
pushdef(`mop',$6)dnl	
RSB_M4_LL_KERNEL_FUNCTION_PROTOTYPE(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`dnl
RSB_M4_LL_KERNEL_FUNCTION_BODY(mtype,matrix_storage,b_rows,b_columns,unrolling,mop)dnl
')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl	
dnl	
dnl
