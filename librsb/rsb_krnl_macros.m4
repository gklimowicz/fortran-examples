dnl
dnl
dnl	@author: Michele Martone
dnl
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
include(`rsb_krnl_vb_macros.m4')dnl
include(`rsb_krnl_bcss_macros.m4')dnl
include(`rsb_krnl_bcoo_macros.m4')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_BCXX_KERNEL_FUNCTION',`dnl
dnl	---------------------------------
dnl	Except for ID, expands to either RSB_M4_BCOO_KERNEL_FUNCTION or RSB_M4_BCSS_KERNEL_FUNCTION.
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl	
pushdef(`transposition',$4)dnl	
pushdef(`k_symmetry',$5)dnl	
pushdef(`b_rows',$6)dnl		block rows
pushdef(`b_columns',$7)dnl	block columns
pushdef(`unrolling',$8)dnl	
pushdef(`mop',$9)dnl	
pushdef(`citype',$10)dnl	
pushdef(`k_diagonal',$11)dnl	
pushdef(`uplo',$12)dnl	
dnl
ifelse(want_what,`ID',`dnl
dnl	FIXME: seemingly broken otherwise.
pushdef(`fid',RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,b_rows,b_columns,unrolling,mop,citype,k_diagonal,uplo))dnl
fid`'dnl
popdef(`fid')dnl
',`dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCOO(matrix_storage),1,`dnl
RSB_M4_BCOO_KERNEL_FUNCTION(want_what,mtype,matrix_storage,transposition,k_symmetry,b_rows,b_columns,unrolling,mop,citype,k_diagonal,uplo)`'dnl
')`'dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
RSB_M4_BCSS_KERNEL_FUNCTION(want_what,mtype,matrix_storage,transposition,k_symmetry,b_rows,b_columns,unrolling,mop,citype,k_diagonal,uplo)`'dnl
')`'dnl
dnl
')dnl
dnl
popdef(`uplo')dnl
popdef(`k_diagonal')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`unrolling')dnl
popdef(`b_columns')dnl
popdef(`b_rows')dnl
popdef(`k_symmetry')dnl
popdef(`transposition')dnl
popdef(`matrix_storage')dnl
popdef(`mtype')dnl
popdef(`want_what')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(mop)
dnl	------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS',`dnl
pushdef(`mop',$1)dnl
dnl
`('RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(,mop,`function_args')`)'dnl
dnl (const struct rsb_mtx_t * mtxAp, const struct rsb_options_t *o, const void * rhs, void * out)dnl
dnl
dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(mop,transposition)
dnl	------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER',`dnl
pushdef(`mop',$1)dnl
pushdef(`transposition',`')dnl
dnl
RSB_M4_PREFIX`'do_`'mop`'`'dnl
dnl
popdef(`transposition')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_NAME(mop)
dnl	------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_NAME',`dnl
pushdef(`mop',$1)dnl
pushdef(`transposition',`')dnl
dnl
rsb_err_t RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(mop,transposition)dnl
dnl
popdef(`mop')dnl
popdef(`transposition')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_GET_NEXT_BLOCK_POINTER_MACRO(matrix_storage)
dnl	---------------------------------------------------
dnl
define(`RSB_M4_GET_NEXT_BLOCK_POINTER_MACRO',`dnl
pushdef(`matrix_storage',$1)dnl
dnl
ifelse(matrix_storage,`VBR',`RSB_GET_NEXT_BLOCK_POINTER')`'dnl
ifelse(matrix_storage,`BCSR',`RSB_BCSR_GET_NEXT_BLOCK_POINTER')`'dnl
dnl	else  ? should give error! fixme :)
dnl
popdef(`matrix_storage')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_GET_FIRST_BLOCK_POINTER_MACRO(matrix_storage)
dnl	----------------------------------------------------
dnl
define(`RSB_M4_GET_FIRST_BLOCK_POINTER_MACRO',`dnl
pushdef(`matrix_storage',$1)dnl
dnl
ifelse(matrix_storage,`VBR',`RSB_GET_FIRST_BLOCK_POINTER')`'dnl
ifelse(matrix_storage,`BCSR',`RSB_BCSR_GET_FIRST_BLOCK_POINTER')`'dnl
dnl	else  ? should give error! fixme :)
dnl
popdef(`matrix_storage')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_GOT_LAST_BLOCK_POINTER_MACRO(matrix_storage)
dnl	---------------------------------------------------
dnl
define(`RSB_M4_GOT_LAST_BLOCK_POINTER_MACRO',`dnl
pushdef(`matrix_storage',$1)dnl
dnl
ifelse(matrix_storage,`VBR',`RSB_GOT_LAST_BLOCK_POINTER')`'dnl
ifelse(matrix_storage,`BCSR',`RSB_BCSR_GOT_LAST_BLOCK_POINTER')`'dnl
dnl	else  ? should give error! fixme :)
dnl
popdef(`matrix_storage')dnl
')dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION(types,mop)
dnl	-------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION',`dnl
pushdef(`types',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`unrolling',`u')dnl
dnl pushdef(`transposition',`')dnl
dnl pushdef(`transposition',$3)dnl
dnl
dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_NAME(mop,`')dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(mop)dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
RSB_M4_DEBUGINFO(``$0'')dnl
ifelse(RSB_M4_IS_IMPLEMENTED_MOP(mop),0,`dnl
{
	return RSB_ERR_UNSUPPORTED_OPERATION;
}
',`dnl
{
	/*!
	 * \ingroup rsb_doc_kernels
	 * Run-time "`'mop`'" kernels dispatching function for a single sparse matrix block.
	 *
	 * \return \rsb_errval_inp_param_msg
	 */
	rsb_err_t errval = RSB_ERR_INTERNAL_ERROR;
dnl	rsb_err_t errval = RSB_ERR_NO_ERROR;
dnl	register rsb_coo_idx_t baserow,basecolumn,rows,columns;
dnl	register rsb_coo_idx_t blockrow,blockcolumn;
dnl	register char *bp=NULL;
dnl	rsb_flags_t symmetry, diagonal;
`#ifdef' RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(`rsb_half_idx_t')
	const rsb_int_t half_storage = rsb__do_is_candidate_size_for_halfword(mtxAp->Mdim,mtxAp->mdim,/*nnz*/0,mtxAp->flags)?`'dnl
RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(`rsb_half_idx_t'):`'dnl
RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(`rsb_coo_idx_t');
#else /* RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(`rsb_half_idx_t') */
	const rsb_int_t half_storage=`'dnl
RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(`rsb_coo_idx_t');
#endif /* RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(`rsb_half_idx_t') */
	const rsb_flags_t symmetry = rsb__get_symmetry_type_flag(mtxAp);
	const rsb_flags_t diagonal = rsb__get_diagonal_type_flag(mtxAp);
	const rsb_type_t typecode = mtxAp->`typecode';

dnl	if(!mtxAp /*|| !mtxAp->options */)
dnl		return RSB_ERR_BADARGS;
	if(RSB_MTX_NOT_OK(mtxAp))
		return RSB_ERR_BADARGS;
dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
dnl	if(RSB_DO_FLAG_HAS_INTERSECTION(mtxAp->flags,RSB_FLAG_ANY_SYMMETRY))
	if(RSB_DO_FLAG_HAS_INTERSECTION(symmetry,RSB_FLAG_ANY_SYMMETRY))
		goto ssoerr;
')dnl
dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
pushdef(`symmetries_here',(RSB_M4_SYMBOL_UNSYMMETRIC))dnl
',`dnl
pushdef(`symmetries_here',`RSB_M4_MATRIX_SYMMETRY')dnl
')dnl

dnl ifelse(mop,`spmv_uxux',`dnl
dnl	if(RSB_IS_ELEMENT_ZERO(betap,mtxAp->`typecode') && RSB_IS_ELEMENT_ONE(alphap,mtxAp->`typecode') )
dnl		return rsb_spmv(...);
dnl')dnl
	switch(diagonal)
	{
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
	case(RSB_M4_MATRIX_DIAGONAL_PREPROCESSOR_SYMBOL(k_diagonal)):
	switch(half_storage)
	{
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
	case(RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(citype)):
	switch(transA)
	{
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
	case(RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL(transposition)):
dnl //	switch(mtxAp->`flags' | RSB_M4_MATRIX_SYMMETRY_PREPROCESSOR_SYMBOL(RSB_M4_SYMBOL_SYMMETRIC))
	switch(symmetry)
	{
foreach(`k_symmetry',symmetries_here,`dnl
	case(RSB_M4_MATRIX_SYMMETRY_PREPROCESSOR_SYMBOL(k_symmetry)):
	switch(mtxAp->`matrix_storage')
	{
foreach(`matrix_storage',RSB_M4_MATRIX_STORAGE,`dnl
	case(RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(matrix_storage)):
dnl		/* return RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(matrix_storage,unrolling,mop,identifier)(...); */
	switch(`typecode')
	{
foreach(`mtype',types,`dnl
	case(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)):`'dnl
dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),1,`dnl
	goto ssoerr;
',`dnl
	if(rsb__is_lower_triangle(mtxAp->flags))
		errval = RSB_M4_BCXX_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,k_symmetry,1,1,unrolling,mop,citype,k_diagonal,`l')( RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCXX_KERNEL_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,1,1,unrolling,mop,citype,k_diagonal,`l'))) );
dnl		errval = RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,unrolling,mop,citype,k_diagonal,`l')`'dnl
dnl (RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCXX_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,`l'))));
	else
		errval = RSB_M4_BCXX_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,k_symmetry,1,1,unrolling,mop,citype,k_diagonal,`u')( RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCXX_KERNEL_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,1,1,unrolling,mop,citype,k_diagonal,`u'))) );
dnl		errval = RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,unrolling,mop,citype,k_diagonal,`u')`'dnl
dnl (RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCXX_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,`u'))));
')dnl	RSB_M4_IS_SPSX_KERNEL_MOP
dnl
',`dnl
dnl
ifelse(RSB_M4_AND(RSB_M4_IS_FORMAT_BCXX(matrix_storage),RSB_M4_IS_SPXX_KERNEL_MOP(mop)),1,`dnl skip_register_block_dispatcher
dnl pushdef(`args',`RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCXX_KERNEL_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,1,1,unrolling,mop,citype,k_diagonal,`g')))')dnl
pushdef(`args',`RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCXX_KERNEL_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,1,1,unrolling,mop,citype,k_diagonal,`g')))')dnl
		errval = RSB_M4_BCXX_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,k_symmetry,1,1,unrolling,mop,citype,k_diagonal,`g')( args );
popdef(`args')dnl
',`dnl skip_register_block_dispatcher
		errval = RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,unrolling,mop,citype,k_diagonal,`g')`'dnl
(RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCXX_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,`g'))));
')dnl skip_register_block_dispatcher
dnl
')dnl
dnl
	break;
')dnl
	default:	goto typerr;
dnl	default:
dnl		RSB_ERROR("Sorry, data type \"%c\" currently not supported.\n",mtxAp->typecode);
dnl		errval = RSB_ERR_UNSUPPORTED_TYPE	;
	}
	break;
')dnl
	default:	goto fmterr;
dnl	default:
dnl	{
dnl		RSB_ERROR("Sorry, matrix storage \"%c\" currently not supported.\n",mtxAp->`matrix_storage');
dnl		FIXME : SOMEWHERE SOMEONE FORGETS TO POPDEF(`matrix_storage') ...
dnl		errval = RSB_ERR_UNSUPPORTED_FORMAT;
dnl	}
	}
	break;
')dnl
	default:	goto symerr;
dnl	default:
dnl	{
dnl		RSB_ERROR("Sorry, this symmetry case (0x%x) is not supported.\n",(rsb_int)symmetry);
dnl		errval = RSB_ERR_UNSUPPORTED_TYPE	;
dnl	}
	}
	break;
')dnl
	default:	goto traerr;
dnl	default:
dnl	{
dnl		RSB_ERROR("Sorry, this transposition case (0x%x) is not supported.\n",(rsb_int)transA);
dnl		errval = RSB_ERR_UNSUPPORTED_TYPE	;
dnl	}
	}
	break;
')dnl
	default:	goto idxerr;
dnl	default:
dnl	{
dnl		RSB_ERROR("Sorry, this coordinate index (0x%x) is not supported.\n",(rsb_int)half_storage);
dnl		errval = RSB_ERR_UNSUPPORTED_FEATURE;
dnl	}
	}
	break;
')dnl
	default:	goto diaerr;
dnl	default:
dnl	{
dnl		RSB_ERROR("Sorry, this diagonal type (0x%x) is not supported.\n",(rsb_int)diagonal);
dnl			errval = RSB_ERR_UNSUPPORTED_FEATURE;
dnl	}
	}
dnl
popdef(`symmetries_here')dnl
dnl
	goto ret;
typerr: errval = RSB_ERR_UNSUPPORTED_TYPE;
	RSB_ERROR("Sorry, data type \"%c\" currently not supported.\n",mtxAp->typecode);
	goto ret;
fmterr:	errval = RSB_ERR_UNSUPPORTED_FORMAT;
	RSB_ERROR("Sorry, matrix storage \"%c\" currently not supported.\n",mtxAp->`matrix_storage');
	goto ret;
symerr:	errval = RSB__ERR_UNSUPPORTED_SYMM;
	RSB_ERROR("Sorry, this symmetry case (0x%x) is not supported.\n",(rsb_int)symmetry);
	goto ret;
traerr: errval = RSB__ERR_UNSUPPORTED_TRANSA;
	RSB_ERROR("Sorry, this transposition case (0x%x) is not supported.\n",(rsb_int)transA);
	goto ret;
idxerr:	errval = RSB__ERR_UNSUPPORTED_IDX_TYPE;
	RSB_ERROR("Sorry, this coordinate index (0x%x) is not supported.\n",(rsb_int)half_storage);
	goto ret;
diaerr:	errval = RSB__ERR_UNSUPPORTED_DIAG;
	RSB_ERROR("Sorry, this diagonal type (0x%x) is not supported.\n",(rsb_int)diagonal);
dnl	goto ret;
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
	goto ret;
ssoerr:	errval = RSB__ERR_NO_SYM_SPSV;
	RSB_ERROR("Sorry, triangular solve of a non-unsymmetric matrix is nonsense.\n");
')dnl
dnl
ret:	return errval;
dnl	return RSB_ERR_INTERNAL_ERROR;	
}
')dnl
')dnl
dnl
popdef(`unrolling')dnl
popdef(`types')dnl
popdef(`mop')dnl
dnl popdef(`transposition')dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ACTUAL_ARGS(mop)
dnl	--------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ACTUAL_ARGS',`dnl
pushdef(`mop',$1)dnl
dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(,mop,`function_args')))dnl
dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ARGS(mop)
dnl	-------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ARGS',`dnl
pushdef(`mop',$1)dnl
dnl
`('double * elapsed_time, RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(,mop,`function_args')`)'dnl
dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_IDENTIFIER(mop)
dnl	-------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_IDENTIFIER',`dnl
pushdef(`mop',$1)dnl
dnl
rsb_do_time_`'mop`'dnl
dnl
dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_NAME(mop)
dnl	-------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_NAME',`dnl
pushdef(`mop',$1)dnl
dnl
rsb_err_t RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_IDENTIFIER($1)`'dnl
dnl
dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION(types,mop)
dnl	--------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION',`dnl
pushdef(`types',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`transposition',RSB_M4_TRANS_N)dnl FIXNE
dnl
dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_NAME(mop)dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ARGS(mop)dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
RSB_M4_DEBUGINFO(``$0'')dnl
{
	/*!
	 * \ingroup gr_bench
	 * This wrapper function will perform the "mop" operation, 
	 * measuring the time elapsed in seconds and writing it in a
	 * user set variable.
         * 
	 * Note that this dispatch function is matrix type indipendent.
         *
	 * \return \rsb_errval_inp_param_msg
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( ! elapsed_time ) return RSB_ERR_BADARGS;

	*elapsed_time = - rsb_time();
	errval = RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(mop,transposition)dnl
	(RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ACTUAL_ARGS(mop));
	
	*elapsed_time += rsb_time(); dnl 	FIXME!

	return errval;
}
')dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_ACTUAL_ARGS(mop)
dnl	-----------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_ACTUAL_ARGS',`dnl
pushdef(`mop',$1)dnl
dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(,mop,`function_args')))dnl
dnl
popdef(`mop')dnl
popdef(`transposition')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_ARGS(mop,mtype)
dnl	----------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_ARGS',`dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
double * total_elapsed_time, double * m_flops, RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(,mop,`function_args')`'dnl
dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_IDENTIFIER(mop,mtype)
dnl	----------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_IDENTIFIER',`dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
RSB_M4_PREFIX`'do_benchmark_`'RSB_M4_CHOPSPACES(mtype)`_'mop`'dnl
dnl
dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_NAME(mop,mtype)
dnl	----------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_NAME',`dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
rsb_err_t RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_IDENTIFIER(mop,mtype)`'dnl
dnl
dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION(types,mop)
dnl	-----------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION',`dnl
pushdef(`types',$1)dnl
pushdef(`mop',$2)dnl
dnl
dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_NAME(mop,mtype)dnl
(RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_ARGS(mop,mtype))dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
{
	/*!
	 * \ingroup gr_bench
	 * This wrapper function will benchmark the "mop" operation
	 * a number of times, measuring the elapsed time in seconds
	 * and writing it in a user set location for a specified matrix.
	 *
	 * It will also add  the performed millions of floating point
	 * operation count in another user specified location.
	 *
	 * \param total_elapsed_time if > 0 on input, will benchmark at least total_elapsed_time seconds
	 * \param m_flops if m_flops > 0 on input, will benchmark at least m_flops times
	 *
	 * If neither of the two input arguments will be set on input,
	 * the benchmark will cease after RSB_BENCHMARK_MIN_RUNS runs or RSB_BENCHMARK_MIN_SECONDS seconds.
	 *
	 * Assuming time_limit = *total_elapsed_time :
	 *
	 * if(time_limit <= 0) will benchmark at least min_runs times
	 * if(time_limit >  0) will benchmark at least min_runs times and for time_limit seconds
	 *
	 * \return \rsb_errval_inp_param_msg
	 */

	double time_limit, elapsed_time;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	int runs = 0, min_runs;

        if( ! total_elapsed_time || ! m_flops)
		return RSB_ERR_BADARGS;

	time_limit = *total_elapsed_time;	/* we read input (FIXME) */
	min_runs   = (int)*m_flops;			/* we read input (FIXME) */

	*total_elapsed_time = RSB_TIME_ZERO;
	*m_flops = RSB_TIME_ZERO;

	if(time_limit <= 0 )
		time_limit = rsb__getenv_real_t("RSB_BENCHMARK_MIN_SECONDS", RSB_BENCHMARK_MIN_SECONDS);

	if(min_runs   <= 0 )
		min_runs = RSB_BENCHMARK_MIN_RUNS;

dnl	//RSB_INFO("will perform min  %d runs, for %lg seconds\n",min_runs, time_limit);
dnl
dnl	// FIXME : seems like this affects performance ...
dnl	// *total_elapsed_time = - rsb_time();
	*total_elapsed_time =0;

	while( ( time_limit? ( *total_elapsed_time < time_limit):0 ) || ( min_runs ? ( runs < min_runs ) : 0 ) )
	{
dnl		//elapsed_time = RSB_TIME_ZERO;
		errval |= dnl
dnl		/* FIXME : use an even more general function here (the following is vbr-only!) */
RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_IDENTIFIER(mop)dnl
`'(&elapsed_time,RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_ACTUAL_ARGS(mop));
dnl		errval = RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(mop)dnl
dnl (RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ACTUAL_ARGS(mop));
dnl
dnl		//*total_elapsed_time += rsb_time();
dnl/*		RSB_INFO("tl : %lg\n",time_limit );*/
dnl/*		RSB_INFO("ss : %lg\n",*total_elapsed_time );*/
dnl/*		RSB_INFO("sse : %lg\n",elapsed_time );*/

		*total_elapsed_time  +=  elapsed_time;
		*m_flops += RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_IDENTIFIER(mop)(mtxAp);
		if(RSB_SOME_ERROR(errval)) {RSB_ERROR(""); return errval;}
		++runs;
	}
dnl	/* FIXME : get rid of this line */
	{
		rsb__fprint_matrix_implementation_code(mtxAp,"mop",RSB_FLAG_NOFLAGS,stderr);
		RSB_STDERR(" : ");
		RSB_STDERR("performed %d runs, %lg/%lg seconds (mop,mtype) \n",runs, *total_elapsed_time,time_limit);
	}

dnl	/*
dnl         * FIXME : this is a candidate location for a conditional performance data printout
dnl         */
dnl
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION_ACTUAL_ARGS(mop,mtype)
dnl	--------------------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION_ACTUAL_ARGS',`dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
dnl RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_ARGS(mop)))`'dnl
dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION(mop,mtype)
dnl	--------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION',`dnl
dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`want_what',$3)dnl
pushdef(`args',`$1,$2')dnl
dnl
ifelse(want_what,`function_identifier',`dnl
RSB_M4_PREFIX`do_fullrangebenchmark_'RSB_M4_CHOPSPACES(mtype)`_'mop`'dnl
')dnl
ifelse(want_what,`function_declaration',`dnl
rsb_err_t $0(args,`function_identifier')dnl
($0(args,`function_args'));
')dnl
ifelse(want_what,`function_args',`dnl
void *VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t rows, rsb_coo_idx_t cols, struct rsb_mop_performance_info_t * mpi, rsb_flags_t flags`'dnl
')dnl
ifelse(want_what,`function_definition',`dnl
rsb_err_t $0(args,`function_identifier')dnl
($0(args,`function_args'))
RSB_M4_DEBUGINFO(``$0'')dnl
{
	/*!
	 * \ingroup gr_bench
	 * Benchmark "mtype" type kernels of operation "mop" 
	 * for a single sparse matrix block.
         * 
dnl         * Therefore, the VBR features of this library will be NOT used here.
dnl	 *
dnl	 * The performance information will be written in a user supplied structure.
dnl         *
	 * \return \rsb_errval_inp_param_msg
	 */
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	int ri=0,ci=0;
	rsb_blk_idx_t br=0,bc=0;
	//rsb_blk_idx_t M_b,K_b;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	rsb_trans_t transAa [] = RSB_ROWS_TRANSPOSITIONS_ARRAY;
	int transi;
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	mtype *out=NULL,*rhs=NULL;
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
	mtype * row_sums=NULL;
')dnl
	rsb_blk_idx_t rua[]=RSB_ROWS_UNROLL_ARRAY;
	rsb_blk_idx_t cua[]=RSB_COLUMNS_UNROLL_ARRAY;
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`1',`dnl
	int incyi, incxi;
	rsb_coo_idx_t incxa[]={1,2}, incya[]={1,2};
',`dnl
dnl	
')dnl

	if(!VA || !IA || !JA || !mpi)
		return RSB_ERR_BADARGS;

	RSB_BZERO_P(mpi);
	mpi->rows = rows;
	mpi->cols=cols;
	mpi->nnz=nnz;

ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`1',`dnl
	for(incyi=0;incyi<sizeof(incya)/sizeof(incya[0]);++incyi)
	for(incxi=0;incxi<sizeof(incxa)/sizeof(incxa[0]);++incxi)
')dnl
	for (transi=0;transi<RSB_TRANSPOSITIONS_ARRAY_LENGTH;++transi)
	{
	rsb_trans_t transA = transAa[transi];
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`1',`dnl
	rsb_coo_idx_t incx = incxa[incxi], incy = incya[incyi];
',`dnl
dnl	incx is sometimes needed for scaling a vector, even if the op itself is strided 1 (and so incyi, incxi above)
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
	rsb_coo_idx_t incx=1,incy=1;
')dnl
')dnl
	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		{
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
			rsb_coo_idx_t bstride = 0;
			rsb_coo_idx_t cstride = 0;
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
ifelse(RSB_M4_IS_SPXM_KERNEL_MOP(mop),`1',`dnl
			rsb_coo_idx_t nrhs=4;
',`dnl
			rsb_coo_idx_t nrhs=1;
')dnl
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
			double alpha=1.0;/* FIXME */
			double * alphap = &alpha;
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_SCALING_KERNEL_MOP(mop),1,`dnl
			double beta=1.0;/* FIXME */
			double * betap=&beta ;
')dnl
dnl
ifelse(mop,`scale',`dnl
			mtype * scale_factors = NULL;
')dnl
dnl
			br = rua[ri];
			bc = cua[ci];
dnl			mtxAp = rsb_allocate_bcsr_sparse_matrix(VA, IA, JA, nnz, typecode, rows, cols, br,bc, flags,&errval);
			mtxAp = rsb__do_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,rows,cols,br,bc,flags,&errval);
			if(!mtxAp||RSB_SOME_ERROR(errval)) {goto erri;}

ifelse(RSB_M4_MAXIMAL_CONFIGURED_BLOCK_SIZE,`1',`',`dnl
			if( ( flags & RSB_FLAG_AUTO_BLOCKING ) != 0)
			{

				/* no need for further benchmarks (FIXME : a temporary, horrible hack! ) */
				ri=ci=-1;
				for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
					for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
						if( rua[ri] == mtxAp->rpntr[1] - mtxAp->rpntr[0] )
							if( cua[ci] == mtxAp->cpntr[1] - mtxAp->cpntr[0] )
								goto ok; /* lol */
				errval = RSB_ERR_INTERNAL_ERROR;
				goto erri;
			}

			ok:
				br = rua[ri];
				bc = cua[ci];
				/* autoblocking found a blocking among the supported ones.
				 * we fill in performance info and quit.
				 */

')dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
			bstride=cols+bc;
			cstride = rows+br;
			rhs = rsb__malloc(mtxAp->el_size*(bstride)*nrhs*incx);
			out = rsb__malloc(mtxAp->el_size*(cstride)*nrhs*incy);
			if(!out || rsb__fill_with_ones(out,mtxAp->typecode,cstride*nrhs,incy)){errval = RSB_ERR_ENOMEM;goto erri;}
			if(!rhs || rsb__fill_with_ones(rhs,mtxAp->typecode,bstride*nrhs,incx)){errval = RSB_ERR_ENOMEM;goto erri;}
			if(!out || !rhs) {errval = RSB_ERR_ENOMEM;goto erri;}
			if(rsb__fill_with_ones(rhs,mtxAp->typecode,(cols)*nrhs,incx))     {errval = RSB_ERR_ENOMEM;goto erri;}
			if(rsb__cblas_Xscal(mtxAp->typecode,(rows+br)*nrhs,NULL,out,incy)) {errval = RSB_ERR_ENOMEM;goto erri;}
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
			row_sums = rsb__malloc(mtxAp->el_size*(rows+br));
			if(!row_sums) {errval = RSB_ERR_ENOMEM;goto erri;}
			if(rsb__fill_with_ones(row_sums,mtxAp->typecode,cols,1))     {errval = RSB_ERR_ENOMEM;goto erri;}
')dnl
ifelse(mop,`scale',`dnl
			scale_factors = rsb__malloc(mtxAp->el_size*(rows+br));
			if(!scale_factors) {errval = RSB_ERR_ENOMEM;goto erri;}
			if(rsb__fill_with_ones(scale_factors,mtxAp->typecode,rows,1))     {errval = RSB_ERR_ENOMEM;goto erri;}
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
')dnl
ifelse(mop,`negation',`dnl
			int please_fix_RSB_M4_ARGS_TO_ACTUAL_ARGS=-1;
')dnl
			
			mpi->seconds[ri][ci] = rsb__getenv_real_t("RSB_BENCHMARK_MIN_SECONDS", RSB_BENCHMARK_MIN_SECONDS); /* min seconds */
			mpi->m_flops[ri][ci] = (double)RSB_BENCHMARK_MIN_RUNS; /* min runs */

dnl			struct rsb_options_t * o = mtxAp->options;
			errval = `'dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_BENCHMARK_FUNCTION_IDENTIFIER(mop,mtype)dnl
( &(mpi->seconds[ri][ci]), &(mpi->m_flops[ri][ci]), RSB_M4_DIRECT_KERNEL_DISPATCH_TIMING_FUNCTION_ACTUAL_ARGS(mop,mtype));
dnl			if(RSB_SOME_ERROR(errval))
dnl				{RSB_PERR_GOTO(erri,"");}
			mpi->fillin[ri][ci] = rsb__do_get_matrix_fillin(mtxAp);
			mpi->e_mflops[ri][ci] =	mpi->m_flops[ri][ci] / mpi->fillin[ri][ci] ;/* new */
	erri:
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
			RSB_MTX_FREE(mtxAp);
			if(RSB_SOME_ERROR(errval))
			{rsb__do_perror(NULL,errval);RSB_PERR_GOTO(err,"");}

			if( RSB_DO_FLAG_HAS( flags , RSB_FLAG_AUTO_BLOCKING ) )
				goto err;
dnl	/* no need for further benchmarks (TODO: solve this better.. */
		}
	}
	}
err:
	return errval;
}
')dnl
dnl
popdef(`args')dnl
popdef(`want_what')dnl
popdef(`mtype')dnl
popdef(`mop')dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ACTUAL_ARGS(mop,mtype)
dnl	-----------------------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ACTUAL_ARGS',`dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
dnl	FIXME
dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ARGS(mtype)
dnl	------------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ARGS',`dnl
pushdef(`mtype',$1)dnl
dnl
(const char * filename, struct rsb_mops_performance_info_t * mspi)dnl
dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_IDENTIFIER(mtype)
dnl	------------------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_IDENTIFIER',`dnl
pushdef(`mtype',$1)dnl
dnl
`rsb_do_completetypebenchmark_'RSB_M4_CHOPSPACES(mtype)`'dnl
dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_NAME(mtype)
dnl	------------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_NAME',`dnl
pushdef(`mtype',$1)dnl
dnl
rsb_err_t RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_IDENTIFIER(mtype)`'dnl
dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION(mtype)
dnl	-------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION',`dnl
dnl
pushdef(`mtype',$1)dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
',`
static RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_NAME(mtype)dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ARGS(mtype)dnl
RSB_M4_DEBUGINFO(``$0'')dnl
{
        /*!
	 * \ingroup gr_bench
	 * Will benchmark all supported matrix operations over the "mtype" type.
	 * over all supported matrix partitionings for a fixed block size.
         *
	 * \return \rsb_errval_inp_param_msg
	 */

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t * IA=NULL,*JA=NULL;
	rsb_coo_idx_t rows=0,cols=0;
	rsb_nnz_idx_t nnz=0;
	void *VA=NULL;

	struct rsb_mop_performance_info_t * mpi = &(mspi->pipmo[0]);
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	rsb_flags_t flagsa [] = { RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS, RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS, RSB_FLAG_USE_HALFWORD_INDICES_CSR,  RSB_FLAG_USE_HALFWORD_INDICES_COO };
	rsb_bool_t is_symmetric = RSB_BOOL_FALSE,is_hermitian = RSB_BOOL_FALSE;
	rsb_bool_t is_lower = RSB_BOOL_FALSE,is_upper= RSB_BOOL_FALSE;
	int flagsi;
	int diagi;

	RSB_BZERO(mspi,sizeof(*mspi));

	if( RSB_SOME_ERROR(rsb__util_mm_info_matrix_f(filename,NULL,NULL,NULL,NULL,&is_symmetric,&is_hermitian,NULL,NULL,NULL,NULL))
	 || RSB_SOME_ERROR(rsb__util_mm_load_matrix_f(filename,&IA,&JA,&VA,&rows,&cols,&nnz,typecode,RSB_FLAG_NOFLAGS,&is_lower,&is_upper)) )
	{
		RSB_PERR_GOTO(err,RSB_ERRMSG_NOTMTXMKT " : %s ..\n",filename);
	}
	
for (diagi=0;diagi<2;++diagi)
for (flagsi=0;flagsi<(sizeof(flagsa)/sizeof(flagsa[0]));++flagsi)
{
	rsb_flags_t flags = flagsa[flagsi];
	if ( diagi)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT);

	if ( is_symmetric )
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_SYMMETRIC);
	else
	{
		if ( is_hermitian )
			RSB_DO_FLAG_ADD(flags,RSB_FLAG_HERMITIAN);
		else
		{
			if ( is_lower )
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
			if ( is_upper )
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
			if ( !is_lower &&  is_upper )
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
			if (  is_lower && !is_upper )
				RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
		}
	}
pushdef(`mopcode',0)dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
pushdef(`mopcode',incr(mopcode))dnl

dnl	/* we benchmark our mtype library implementation for operation mop */
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
	if(!(flags & RSB_FLAG_ANY_SYMMETRY))
	if( (flags & RSB_FLAG_TRIANGULAR ))
	{
')dnl
	errval = dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION(mop,mtype,`function_identifier')dnl
(RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION(mop,mtype,`function_args'))));
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
	}
')dnl
	++mpi;
	if(RSB_SOME_ERROR(errval))goto err;
')dnl
	mpi-=mopcode;
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
popdef(`mopcode')dnl
')dnl
popdef(`mopcode')dnl
dnl
	
dnl		performance info dumpout

pushdef(`mopcode',0)dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
pushdef(`mopcode',incr(mopcode))dnl
dnl	may dump extra performance info here
	errval = rsb__dump_performance_info(mpi,"RSB_M4_DUMP_PERFOMANCE_INFO_RECORD_IDENTIFIER(mtype,mop)");
	if(RSB_SOME_ERROR(errval))goto err;
	++mpi;
')dnl
	mpi-=mopcode;
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
popdef(`mopcode')dnl
')dnl
popdef(`mopcode')dnl
dnl
}

	err:
	RSB_CONDITIONAL_FREE(IA);
	RSB_CONDITIONAL_FREE(JA);
	RSB_CONDITIONAL_FREE(VA);
	return errval;
}
popdef(`mtype')dnl
')dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_DUMP_PERFORMANCE_INFO_ARRAY_FUNCTION_IDENTIFIER',`dnl
dnl
RSB_M4_PREFIX`dump_performance_array'dnl
dnl
')dnl
dnl
define(`RSB_M4_DUMP_PERFORMANCE_INFO_ARRAY_FUNCTION_ARGS',`dnl
dnl
`(const char * an, const double*array)'dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_DUMP_PERFORMANCE_INFO_ARRAY_FUNCTION
dnl	-------------------------------------------
dnl	FIXME : this should go in some other file
dnl
define(`RSB_M4_DUMP_PERFORMANCE_INFO_ARRAY_FUNCTION',`dnl
dnl
rsb_err_t` 'dnl
RSB_M4_DUMP_PERFORMANCE_INFO_ARRAY_FUNCTION_IDENTIFIER()dnl
RSB_M4_DUMP_PERFORMANCE_INFO_ARRAY_FUNCTION_ARGS()dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
RSB_M4_DEBUGINFO(``$0'')dnl
{
	/*!
	 * \ingroup gr_bench
	 * A benchmark info dumping function.
         *
	 * \return \rsb_errval_inp_param_msg
         *
	 * FIXME : UNFINISHED
	 */
#if RSB_ALLOW_STDOUT
	int ri,ci;
	rsb_blk_idx_t rua[]=RSB_ROWS_UNROLL_ARRAY;
	rsb_blk_idx_t cua[]=RSB_COLUMNS_UNROLL_ARRAY;
	if(!array || !an)
		return RSB_ERR_BADARGS;

/*	RSB_STDOUT("const double %s [RSB_ROWS_UNROLL_ARRAY_LENGTH][RSB_COLUMNS_UNROLL_ARRAY_LENGTH] = \n",an);*/
	RSB_STDOUT(".%s = \n",an);
	RSB_STDOUT("{");
	RSB_STDOUT("\t/*");
	for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci) RSB_STDOUT("%d, ",cua[ci]);
	RSB_STDOUT("columns per block */\n");
		
	for(ri=0;ri<RSB_ROWS_UNROLL_ARRAY_LENGTH;++ri)
	{
		RSB_STDOUT("\t{");
		for(ci=0;ci<RSB_COLUMNS_UNROLL_ARRAY_LENGTH;++ci)
		{
			if(ci)RSB_STDOUT(",");
			RSB_STDOUT(" %lg",array[ri*RSB_ROWS_UNROLL_ARRAY_LENGTH+ci]);
		}
		RSB_STDOUT(" }, /* %d rows per block */\n",rua[ri]);
	}
	RSB_STDOUT("},\n");
	return RSB_ERR_NO_ERROR;
#else /* RSB_ALLOW_STDOUT */
	return RSB_ERR_UNSUPPORTED_FEATURE;
#endif /* RSB_ALLOW_STDOUT */
}
dnl
')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DUMP_PERFOMANCE_INFO_RECORD_IDENTIFIER(MTYPE,MOP)
dnl	--------------------------------------------------------
dnl
dnl
define(`RSB_M4_DUMP_PERFOMANCE_INFO_RECORD_IDENTIFIER',`dnl
pushdef(`mtype',$1)dnl
pushdef(`mop',$2)dnl
dnl
`pi_'RSB_M4_CHOPSPACES(mtype)`_'mop`'dnl
dnl
popdef(`mop')dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ACTUAL_ARGS(mtype)
dnl	-------------------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ACTUAL_ARGS',`dnl
pushdef(`mtype',$1)dnl
dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ARGS(mtype))dnl
dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_ARGS()
dnl	---------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_ARGS',`dnl
dnl
(const int argc, char *const argv[])dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_IDENTIFIER(mop)
dnl	-------------------------------------------------------
dnl
define(`RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_IDENTIFIER',`dnl
pushdef(`mop',$1)dnl
dnl
RSB_M4_PREFIX`estimate_mflops_per_op_'mop`'dnl
dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl	RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_ARGS(mop)
dnl	-------------------------------------------------
dnl
define(`RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_ARGS',`dnl
pushdef(`mop',$1)dnl
dnl
`(const struct rsb_mtx_t * mtxAp)'dnl
dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION(mop)
dnl	--------------------------------------------
dnl
define(`RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION',`dnl
pushdef(`mop',$1)dnl
dnl
`double 'dnl
RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_IDENTIFIER(mop)dnl
RSB_M4_ESTIMATE_MFLOPS_PER_MOP_FUNCTION_ARGS(mop)dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
RSB_M4_DEBUGINFO(``$0'')dnl
{
	/*!
	 * \ingroup gr_internals
	 * Return the canonical count of floating point operations
	 * needed to perform the "mop" matrix operation.
	 * In the case of a complex type, the number of operations is increased by six:
	 *  (a +bi)*(c+di) = (ac-bd)+(ad+bc)i
	 * accounting for an extra of four real multiplications and two real additions.
	 * In the case of symmetric/hermitian, this is further multiplied by two.
	 * Note that this count is very rough: e.g. ignores diagonal implicit or diagonal with symmetry.
	 */

	const double M_  = 1000000.0;
	const double Ec = ((double)mtxAp->element_count);
	double Me = Ec;
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
	if(RSB_IS_MATRIX_TYPE_COMPLEX(mtxAp->typecode)) { Me=8*Ec; } else { Me=2*Ec; }
')dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),`1',`dnl
	if(rsb__is_not_unsymmetric(mtxAp)){ Me*=2; }
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
')dnl
ifelse(mop,`negation',`dnl
')dnl
ifelse(mop,`scale',`dnl
')dnl
	Me /= M_;
	return Me;
}
dnl
')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_IDENTIFIER()
dnl	---------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_IDENTIFIER',`dnl
dnl
`rsb__do_completebenchmark'dnl
dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_NAME()
dnl	---------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_NAME',`dnl
dnl
rsb_err_t`' RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_IDENTIFIER`'dnl
dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION()
dnl	----------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION',`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_NAME`'dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_ARGS`'dnl
;
',`
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_NAME`'dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETEBENCHMARK_FUNCTION_ARGS`'dnl
RSB_M4_DEBUGINFO(``$0'')dnl
{
	/*!
	 * \ingroup gr_bench
	 * A complete benchmark program.
	 * Benchmark all supported matrix operations over all supported types
	 * over all supported matrix partitionings for a fixed block size.
         *
	 * \return \rsb_errval_inp_param_msg
dnl	 * Originaly this info was meant to be processed and dumped in a header file.
	 */
	struct rsb_global_performance_info_t mspis;
	struct rsb_mops_performance_info_t * mspi = &(mspis.gpi[0]);

	rsb_option_t options[] = {
	    {"matrix-filename",	required_argument, NULL, 0x66},  /* f */
	    {0,0,0,0}
	};
	const char * filename=NULL;
	int c=0;
	int opt_index=0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if(rsb_lib_init(RSB_NULL_INIT_OPTIONS))goto err;

	for (;;)
	{
		c = rsb__getopt_long(argc, argv, "f:" , options, &opt_index);/* Flawfinder: ignore */
		if (c == -1)break;
		switch (c)
		{
			case 0x66:/* f */
			filename = optarg;
			break;
	    	}
	}

foreach(`mtype',types,`dnl

	errval=dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_IDENTIFIER(mtype)dnl
(RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ACTUAL_ARGS(mtype));
	if(RSB_SOME_ERROR(errval)) return errval;
	++mspi;
')dnl

	if( rsb_lib_exit(RSB_NULL_EXIT_OPTIONS) )
		return RSB_ERR_INTERNAL_ERROR;
	return RSB_ERR_NO_ERROR;
	err:
	return RSB_ERR_INTERNAL_ERROR;
}
')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION(types,mop,want_what)
dnl	---------------------------------------------------------------------
dnl
dnl	EDIT THIS MACRO TO SPECIFY ARGS TO NEW KERNELS
dnl
define(`RSB_M4_MULTI_BLOCK_KERNEL_TYPE_DISPATCH_FUNCTION',`dnl
pushdef(`types',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`want_what',$3)dnl
pushdef(`args',`$1,$2')dnl
dnl
dnl
ifelse(want_what,`function_identifier',`dnl
RSB_M4_PREFIX`do_'mop`_with_macros_vbr'dnl
')dnl
dnl
ifelse(want_what,`function_declaration',`dnl
rsb_err_t $0(args,`function_identifier')dnl
($0(args,`function_args'));dnl
')dnl
dnl
ifelse(want_what,`function_args',`dnl
dnl
dnl
ifelse(RSB_M4_IS_READONLY_KERNEL_MOP(mop),1,`dnl
const struct rsb_mtx_t * mtxAp`'dnl
',`dnl
struct rsb_mtx_t * mtxAp`'dnl
')dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
,const void * RSB_M4_RESTRICT rhs, void * RSB_M4_RESTRICT out`'dnl
')dnl
ifelse(RSB_M4_IS_OP_SCALING_KERNEL_MOP(mop),`1',`dnl
,const void * alphap`'dnl
')dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),`1',`dnl
,const void * betap`'dnl
')dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`1',`dnl
,rsb_coo_idx_t incx, rsb_coo_idx_t incy`'dnl
')dnl
dnl
,const rsb_trans_t transA`'dnl
ifelse(mop,`scale',`dnl
,const void * scale_factors`'dnl
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
,void * row_sums`'dnl
')dnl
ifelse(mop,`negation',`dnl
,int please_fix_RSB_M4_ARGS_TO_ACTUAL_ARGS`'dnl
')dnl
dnl
dnl
dnl
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
rsb_err_t $0(args,`function_identifier')dnl
($0(args,`function_args'))
{
RSB_M4_DEBUGINFO(``$0'')dnl
	/*!
	 * \ingroup rsb_doc_kernels
	 * Kernel function dispatching will be performed inline, after type dispatching, in a separate function.
         *
	 * \return \rsb_errval_inp_param_msg
	 */
dnl
dnl	removed old junk in revision 625 ...
dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
}
')dnl body
dnl
popdef(`args')dnl
popdef(`want_what')dnl
popdef(`mop')dnl
popdef(`types')dnl
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_TESTING_FUNCTION(types,mop)
dnl	---------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_TESTING_FUNCTION',`dnl
dnl
pushdef(`types',$1)dnl
pushdef(`mop',$2)dnl
dnl
dnl	
#if defined(RSB_WANT_KERNELS_DEBUG) && (RSB_WANT_KERNELS_DEBUG>0)
rsb_err_t RSB_M4_PREFIX`'mop`_testing'dnl
RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(mop)dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
{
RSB_M4_DEBUGINFO(``$0'')dnl
	/*!
	 * \ingroup gr_debug
	 * This is a trivial reference implementation of the "mop" kernel, and 
	 * its numerical results will be used to shed some evidence if bugs 
	 * should be introduced in performance computational kernels.
	 * 
	 * It should be used for debugging or comparing with performance optimized
	 * functions.
         *
	 * \return \rsb_errval_inp_param_msg
	 */
	
	register rsb_coo_idx_t baserow = RSB_INI,basecolumn = RSB_INI,rows = RSB_INI,columns = RSB_INI;
	register rsb_coo_idx_t blockrow = RSB_INI,blockcolumn = RSB_INI;
dnl	register char *bp=0;
	register rsb_byte_t *bp=0;
ifelse(RSB_M4_NOT(RSB_M4_IS_STRIDED_KERNEL_MOP(mop)),1,`dnl
dnl	rsb_coo_idx_t incx=2,incy=2;
	rsb_coo_idx_t incxa[]={1,2}, incya[]={1,2};
dnl	incx=2,incy=2;
dnl	to avoid "unused variable"-like warnings
')dnl

	if(!mtxAp /*|| !mtxAp->options*/ )return RSB_ERR_BADARGS;
	{
	RSB_GET_FIRST_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
foreach(`mtype',types,`dnl
	if(mtxAp->`typecode' == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))
	{

ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),1,`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype), betap, mtxAp->nr, out, 1);/* we scale the destination vector */
')dnl
ifelse(RSB_M4_IS_ZEROING_KERNEL_MOP(mop),1,`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype),mtxAp->nr,NULL,out,incy);
')dnl
ifelse(RSB_M4_IS_SPXM_KERNEL_MOP(mop),`1',`dnl
dnl	if(mtxAp && mout) rsb__cblas_Xscal(mtxAp->typecode,nrhs*mtxAp->nr,NULL,mout,incy);// NEW
	if(mtxAp && out) rsb__cblas_Xscal(mtxAp->typecode,nrhs*mtxAp->nr,NULL,out,incy);// NEW
')dnl
	
	while(!RSB_GOT_LAST_BLOCK_POINTER(mtxAp))
	{
ifelse(mop,`scale',`dnl
		mtype* a = (mtype*)bp;
		rsb_coo_idx_t i,j;
		for(i=0;i<rows;++i)
			for(j=0;j<columns;++j)
				a[i*columns+j]*=((const mtype*)scale_factors)[i];
')dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
		const mtype* a = (const mtype*)bp;
		const mtype* b = ((const mtype*)rhs)+mtxAp->cpntr[blockcolumn];
		mtype* c = ((mtype*)out)+mtxAp->rpntr[blockrow];
		rsb_coo_idx_t i,j;
		c=c;/* here just to prevent from compiler warning */

#if 0
		for(i=0;i<rows;++i)
			for(j=0;j<columns;++j)
				c[i]+=a[i*columns+j]*b[j];
#else
		/*
		 * this code will emulate the same kernel order!
		 * it should generate the same numerical roundoff errors the current kernel would
		 * */
		for(i=0;i<rows;++i)
		{
			mtype rs=0;
			for(j=0;j<columns;++j)
				rs+=a[i*columns+j]*b[j];
ifelse(mop,`spmv_uaua',`dnl
			c[i]+=rs;
');
ifelse(mop,`spmv_unua',`dnl
			c[i]-=rs;
');
		}
#endif /* 0 */

')dnl
ifelse(RSB_M4_IS_SPXM_KERNEL_MOP(mop),`1',`dnl
		const mtype* a = (const mtype*)bp;
dnl		const mtype* b = ((const mtype*)mrhs)+mtxAp->cpntr[blockcolumn];
		const mtype* b = ((const mtype*)rhs)+mtxAp->cpntr[blockcolumn];
		mtype* c = ((mtype*)out)+mtxAp->rpntr[blockrow];
dnl		mtype* c = ((mtype*)mout)+mtxAp->rpntr[blockrow];
		rsb_coo_idx_t i,j,k;
		for(i=0;i<rows;++i)
			for(j=0;j<columns;++j)
				for(k=0;k<nrhs;++k)
					c[k*cstride+i]+=a[i*columns+j]*b[j+k*bstride];

')dnl
ifelse(RSB_M4_MEMBER(mop,`spsv_uxua',`spsv_sxsx'),1,`dnl
/*	FIXME : UNFINISHED */
')dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	dnl	rsb_coo_idx_t incx=2,incy=2;
')dnl
ifelse(RSB_M4_MEMBER(mop,`spmv_sxsx',`spsv_sa',`spmv_uxux'),1,`dnl
/*	FIXME : UNFINISHED */
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
		const mtype* a = (const mtype*)bp;
		mtype* row_sums_=row_sums;
		rsb_coo_idx_t i,j;
		for(i=0;i<rows;++i)
			for(j=0;j<columns;++j)
ifelse(mop,`infty_norm',`dnl
				row_sums_[mtxAp->rpntr[blockrow]+i]+=RSB_M4_ABS(mtype,a[i*columns+j]);
')dnl
ifelse(mop,`rowssums',`dnl
				row_sums_[mtxAp->rpntr[blockrow]+i]+=a[i*columns+j];
')dnl
')dnl
	
ifelse(mop,`negation',`dnl
		mtype* a = (mtype*)bp;
		rsb_coo_idx_t i,j;
		for(i=0;i<rows;++i)
			for(j=0;j<columns;++j)
				a[i*columns+j]=-a[i*columns+j];
')dnl
		RSB_GET_NEXT_BLOCK_POINTER(bp,mtxAp,baserow,basecolumn,rows,columns,blockrow,blockcolumn);
	}
	}
	else
')dnl
	{
		RSB_ERROR("Sorry, data type \"%c\" currently not supported.\n",mtxAp->typecode);
		return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	}
	return RSB_ERR_NO_ERROR;	
}
dnl
')dnl
dnl
#endif /* RSB_WANT_KERNELS_DEBUG */
dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ACTUAL_ARGS(mtype)
dnl	-------------------------------------------------------------------------------
dnl
define(`RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ACTUAL_ARGS',`dnl
pushdef(`mtype',$1)dnl
dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_DIRECT_KERNEL_DISPATCH_COMPLETETYPEBENCHMARK_FUNCTION_ARGS(mtype))dnl
dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
