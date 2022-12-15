dnl
dnl	@author: Michele Martone
dnl
dnl
dnl
define(`RSB_M4_ARE_KERNEL_GENERATION_PARMS_ALLOWED',`dnl
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl
pushdef(`transposition',$4)dnl
pushdef(`k_symmetry',$5)dnl
pushdef(`unrolling',$6)dnl	
pushdef(`b_rows',$7)dnl		block rows
pushdef(`b_columns',$8)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t ')dnl integer type (for indices)
pushdef(`mop',`$9')dnl
pushdef(`citype',`$10')dnl
pushdef(`k_diagonal',`$11')dnl
pushdef(`uplo',$12)dnl
dnl
RSB_M4_AND(dnl
RSB_M4_IMPLY(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_NOT(RSB_M4_SAME(uplo,`g'))),dnl
RSB_M4_IMPLY(RSB_M4_NOT(RSB_M4_IS_SPSX_KERNEL_MOP(mop)),RSB_M4_SAME(uplo,`g')),dnl
1)`'dnl
dnl
dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
popdef(`b_rows')dnl
popdef(`b_columns')dnl
popdef(`transposition')dnl
popdef(`k_symmetry')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
popdef(`k_diagonal')dnl
popdef(`want_what')dnl
popdef(`uplo')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	These functions dispatch on the column size, calling the
dnl	proper kernels.
dnl
dnl	They assume type dispatching has just been performed.
dnl
dnl
dnl	RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(want_what,mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
dnl	-----------------------------------------------------------------------------------------------------------------------------------
dnl
define(`RSB_M4_BCXX_KERNEL_SIZE_DISPATCH_FUNCTION',`dnl
dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCOO(matrix_storage),`1',`dnl
dnl
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION($@)`'dnl
dnl
')dnl
dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),`1',`dnl
dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION($@)`'dnl
dnl
')dnl
dnl
dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION',`dnl
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl
pushdef(`transposition',$4)dnl
pushdef(`k_symmetry',$5)dnl
pushdef(`unrolling',$6)dnl	
dnl pushdef(`b_rows',$7)dnl		block rows
dnl pushdef(`b_columns',$8)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t ')dnl integer type (for indices)
pushdef(`mop',`$9')dnl
pushdef(`citype',`$10')dnl
pushdef(`k_diagonal',`$11')dnl
pushdef(`uplo',$12)dnl
dnl
ifelse(RSB_M4_ARE_KERNEL_GENERATION_PARMS_ALLOWED(want_what,mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo),`1',`dnl
dnl
ifelse(want_what,`DOC',`dnl
	/*  TODO */
')dnl
dnl
ifelse(want_what,`all',`dnl
dnl `/* This code is intended for a block compressed sparse stripe matrix. */'
ifdef(`ONLY_WANT_HEADERS',`dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`function_declaration',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
',`dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`function_definition',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
')dnl
dnl
dnl
dnl
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
rsb_err_t RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,unrolling,mop,citype,k_diagonal,uplo)dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`BODY',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
')dnl
dnl
ifelse(want_what,`function_declaration',`dnl
rsb_err_t RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,unrolling,mop,citype,k_diagonal,uplo)dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo);dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
dnl
dnl
pushdef(`matrix_structs',`const itype Mdim,const itype mdim,const citype * RSB_M4_RESTRICT bindx,const rsb_nnz_idx_t * RSB_M4_RESTRICT bpntr,const rsb_nnz_idx_t *RSB_M4_RESTRICT indptr,const rsb_coo_idx_t * RSB_M4_RESTRICT rpntr,const rsb_coo_idx_t * RSB_M4_RESTRICT cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags')dnl
(`'dnl
ifelse(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop),1,`dnl
dnl
dnl	no restrict on aliasing ops
dnl
ifelse(RSB_M4_IS_ALLOWING_ALIASING_KERNEL_MOP(mop),1,`dnl
const mtype * RSB_M4_RESTRICT VA, const mtype * rhs, mtype * out, matrix_structs`'dnl
',`dnl
const mtype * RSB_M4_RESTRICT VA, const mtype * RSB_M4_RESTRICT rhs, mtype * RSB_M4_RESTRICT out, matrix_structs`'dnl
')dnl
')dnl
ifelse(RSB_M4_IS_OP_SCALING_KERNEL_MOP(mop),`1',`dnl
,const mtype * RSB_M4_RESTRICT alphap`'dnl
')dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),`1',`dnl
,const mtype * RSB_M4_RESTRICT betap`'dnl
')dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`1',`dnl
,rsb_coo_idx_t incx, rsb_coo_idx_t incy`'dnl
')dnl
ifelse(mop,`spmm_az',`dnl
dnl
dnl	FIXME
dnl
const itype bstride, const itype cstride, const itype nrhs`'dnl
')dnl
ifelse(mop,`scale',`dnl
mtype * VA, matrix_structs, const mtype *scale_factors`'dnl
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
const mtype * VA, mtype * row_sums, matrix_structs`'dnl
')dnl
ifelse(mop,`negation',`dnl
mtype * VA, matrix_structs`'dnl
')dnl
)dnl
dnl
')dnl
dnl
dnl
ifelse(want_what,`BODY',`dnl
dnl
dnl
{
	RSB_M4_DEBUGINFO(``$0'')dnl
dnl	/*!  \ingroup rsb_doc_kernels
	/*
	 * This function will dispatch the specialized looped kernel function for 
	 * performing the desired matrix operation ("mop") for the current fixed
	 * block size.
	 *
	 * \return \rsb_errval_inp_param_msg
ifelse(RSB_M4_IS_KERNEL_REGISTER_BLOCKED(b_rows,b_columns),1,`dnl
	 *
	 * Since this is strictly blocked code, you should allow the rhs and the out
	 * vector to accept a small overflow not bigger, respectively, than
	 *       mod(blockrows-mod(matrixrows,blockrows),blockrows)
	 * and
	 *       mod(blockcols-mod(matrixcols,blockcols),blockcols)
')dnl
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
pushdef(`args',`RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo))')dnl

	register rsb_coo_idx_t columns,rows;
#if !RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS
	if(cpntr && rpntr)
	{
		columns=cpntr[1]-cpntr[0];
		rows   =rpntr[1]-rpntr[0];
		RSB_ASSERT(columns==1 && rows==1);
	}
	else
#endif /* RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS */
dnl #if RSB_EXPERIMENTAL_WANT_PURE_BCSS
ifelse(RSB_M4_WANT_20110206_BOUNDED_BOX_PATCH,1,`dnl
dnl 20110206	set the following 
		columns = rows=1;
',`dnl
dnl 20110206	and commented the following 
		columns=bc,rows=br;
')dnl
dnl #else
dnl 		columns = rows=1;
dnl #endif

switch(rows)
{
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
	case rowsu:
	{switch(columns)
	{
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
		case colsu:/* rowsu colsu matrix_storage */
		errval = RSB_M4_BCSS_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)( args );
		break;
')dnl

#ifdef RSB_WANT_LOOPING_KERNELS 
default: 	errval = RSB_M4_BCSS_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,k_symmetry,rowsu,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,citype,k_diagonal,uplo)( args );
#else /* RSB_WANT_LOOPING_KERNELS */
default:	errval = RSB_ERR_UNSUPPORTED_OPERATION;
#endif /* RSB_WANT_LOOPING_KERNELS */
	}}
	break;
')dnl

#ifdef RSB_WANT_LOOPING_KERNELS 
default:	errval = RSB_M4_BCSS_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,k_symmetry,RSB_M4_ROWS_FALLBACK_UNROLL,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,citype,k_diagonal,uplo)( args );
#else /* RSB_WANT_LOOPING_KERNELS */
default:	errval = RSB_ERR_UNSUPPORTED_OPERATION;
#endif /* RSB_WANT_LOOPING_KERNELS */
};
popdef(`args')dnl
')dnl
	return errval;
}
dnl
')dnl
dnl
')dnl
dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
dnl popdef(`b_rows')dnl
dnl popdef(`b_columns')dnl
popdef(`transposition')dnl
popdef(`k_symmetry')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
popdef(`k_diagonal')dnl
popdef(`want_what')dnl
popdef(`uplo')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	These functions will perform their operations on fixed block matrices.
dnl
define(`RSB_M4_BXXX_KERNEL_FUNCTION_HAS_IMPLEMENTATION',`dnl
dnl
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl	
pushdef(`transposition',$4)dnl	
pushdef(`k_symmetry',$5)dnl	
pushdef(`b_rows',$6)dnl		block rows
pushdef(`b_columns',$7)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t ')dnl integer type (for indices)
pushdef(`unrolling',$8)dnl	
pushdef(`mop',$9)dnl	
pushdef(`citype',$10)dnl	
pushdef(`k_diagonal',$11)dnl	
pushdef(`uplo',$12)dnl	
dnl
ifelse(dnl
dnl	Emit a non-empty character (non-`') if no code is to be generated:
dnl ifelse(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),RSB_M4_NOT(transposed)),1,`no',`')`'dnl	CSC SPSV
dnl ifelse(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),transposed),1,`no'`')dnl	CSR transposed SPSV
ifelse(RSB_M4_IMPLY(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_KERNEL_REGISTER_UNBLOCKED(b_rows,b_columns)),`1',`',`no')`'dnl SPSV for non 1x1 blockings gets blocked
ifelse(RSB_M4_OR(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),RSB_M4_IS_KERNEL_REGISTER_UNBLOCKED(b_rows,b_columns)),`1',`',`no')`'dnl any symmetric kernel for non 1x1 blockings
dnl	TODO : should modify RSB_M4_EXTRA_SYMMETRIC_DIAGONAL_FIXING_KERNEL to support k_symmetry and blocking
ifelse(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry)),1,`no',`')`'dnl any SPSV symmetric
dnl
,`',`1',`0')dnl
dnl
popdef(`uplo')dnl
popdef(`want_what')dnl
popdef(`k_diagonal')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
popdef(`k_symmetry')dnl
popdef(`transposition')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_IMPOSSIBLE_IMPLEMENTATION_ERROR',`dnl
dnl	RSB_ERR_UNIMPLEMENTED_YET
dnl	RSB_ERR_UNSUPPORTED_FEATURE
dnl	RSB_ERR_UNSUPPORTED_OPERATION
dnl	RSB_ERR_BADARGS|RSB_ERR_UNSUPPORTED_OPERATION
RSB_ERR_BADARGS`'dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_BXXX_KERNEL_FUNCTION_HELP',`dnl
dnl
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl	
pushdef(`transposition',$4)dnl	
pushdef(`k_symmetry',$5)dnl	
pushdef(`b_rows',$6)dnl		block rows
pushdef(`b_columns',$7)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t ')dnl integer type (for indices)
pushdef(`unrolling',$8)dnl	
pushdef(`mop',$9)dnl	
pushdef(`citype',$10)dnl	
pushdef(`k_diagonal',$11)dnl	
pushdef(`uplo',$12)dnl	
dnl
	/**
	 * \ingroup rsb_doc_kernels
ifelse(RSB_M4_MEMBER(mop,`spsv_uxua'),1,`dnl
	 * Computes \f$y \leftarrow RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A')^{-1} \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`spmv_unua',`dnl
	 * Computes \f$y \leftarrow y - RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`spmv_uaua',`dnl
	 * Computes \f$y \leftarrow y + RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`spmv_sxsa',`dnl
	 * Computes \f$y \leftarrow \beta \cdot y + \alpha \cdot RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
	 * with incx and incy as x and y vector strides
')dnl
ifelse(mop,`spmv_sxsx',`dnl
	 * Computes \f$y \leftarrow \beta \cdot y + \alpha \cdot RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
	 * with incx and incy as x and y vector strides
')dnl
ifelse(mop,`spmv_sasa',`dnl
	 * Computes \f$y \leftarrow y + RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`spmv_uxua',`dnl
	 * Computes \f$y \leftarrow y + \alpha \cdot RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`spmv_uxux',`dnl
	 * Computes \f$y \leftarrow \beta \cdot y + \alpha \cdot RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`spmm_az',`dnl
	 * Computes \f$y \leftarrow RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`infty_norm',`dnl
	 * Computes \f$ \|A\|_{\infty} \f$ (or rather, \f$ row\_sums_i \leftarrow \sum_{j=0}^{mdim} A_{ij} ), where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A').\f$
')dnl
ifelse(mop,`rowssums',`dnl
	 * Computes \f$ \|A\|_{1} \f$ (or rather, \f$ row\_sums_i \leftarrow \sum_{i=0}^{Mdim} A^{T}_{ij} ), where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A').\f$
')dnl
ifelse(mop,`spmv_uauz',`dnl
	 * Computes \f$y \leftarrow RSB_M4_TRANSPOSITION_OP_EFFECT(transposition,`A') \cdot x, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A'). \f$
')dnl
ifelse(mop,`scale',`dnl
	 * Computes \f$A \leftarrow A\cdot P, P_{ii}=s_{i}, where RSB_M4_SYMMETRY_EFFECT(k_symmetry,`A').\f$
')dnl
ifelse(mop,`negation',`dnl
	 * Computes \f$A \leftarrow - A \f$
')dnl
         * A blocked b_rows x b_columns, stored in matrix_storage format, RSB_M4_MATRIX_DIAGONAL_DENOMINATION(k_diagonal), of `type' mtype, with citype column indices.
dnl
ifelse(RSB_M4_BXXX_KERNEL_FUNCTION_HAS_IMPLEMENTATION($@),`1',`dnl
	 * \return \rsb_errval_inp_param_msg
	 */
',`dnl
dnl
	 * \return RSB_M4_IMPOSSIBLE_IMPLEMENTATION_ERROR (this function is not implemented).
	 */
dnl
')dnl
dnl
popdef(`uplo')dnl
popdef(`want_what')dnl
popdef(`k_diagonal')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
popdef(`k_symmetry')dnl
popdef(`transposition')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl	These functions will perform their operations on fixed block matrices.
dnl
define(`RSB_M4_BCSS_KERNEL_FUNCTION',`dnl
dnl
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl	
pushdef(`transposition',$4)dnl	
pushdef(`k_symmetry',$5)dnl	
pushdef(`b_rows',$6)dnl		block rows
pushdef(`b_columns',$7)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t ')dnl integer type (for indices)
pushdef(`unrolling',$8)dnl	
pushdef(`mop',$9)dnl	
pushdef(`citype',$10)dnl	
pushdef(`k_diagonal',$11)dnl	
pushdef(`uplo',$12)dnl	
dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
pushdef(`fid',RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,b_rows,b_columns,unrolling,mop,citype,k_diagonal,uplo))dnl
dnl
ifelse(RSB_M4_ARE_KERNEL_GENERATION_PARMS_ALLOWED(want_what,mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo),`1',`dnl
dnl
ifelse(want_what,`all',`dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
rsb_err_t fid`'dnl
RSB_M4_BCSS_KERNEL_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,b_rows,b_columns,unrolling,mop,citype,k_diagonal,uplo)dnl
')dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
dnl /* begin of fid function */
RSB_M4_BCSS_KERNEL_FUNCTION(`BODY',mtype,matrix_storage,transposition,k_symmetry,b_rows,b_columns,unrolling,mop,citype,k_diagonal,uplo)dnl
')dnl
')dnl
dnl
ifelse(want_what,`ID',`dnl
fid`'dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)`'dnl
')dnl
dnl
ifelse(want_what,`BODY',`dnl
dnl
{
dnl
dnl	The body of a CSR/CSC computational kernel.
dnl
dnl	RSB_M4_DEBUGINFO(``$0'')dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
dnl
pushdef(`total_block_columns',ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`Mdim',`mdim'))dnl
pushdef(`total_block_rows',ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`mdim',`Mdim'))dnl
pushdef(`total_rows',ifelse(unrolling,`l',rpntr[total_block_rows],total_block_rows*b_rows))dnl
pushdef(`total_columns',ifelse(unrolling,`l',cpntr[total_block_columns],total_block_columns*b_columns))dnl
dnl
ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
pushdef(`mi',`i')dnl
pushdef(`Mi',`j')dnl
')dnl
ifelse(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),1,`dnl
pushdef(`mi',`j')dnl
pushdef(`Mi',`i')dnl
')dnl
dnl
dnl	FIXME : out_dim should depend on the operation!
dnl
pushdef(`out_dim',ifelse(transposition,RSB_M4_TRANS_N,total_rows,total_columns))dnl
dnl
pushdef(`is_zero_acc_spsv_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_OR(RSB_M4_AND(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),RSB_M4_SAME(transposition,RSB_M4_TRANS_N)),RSB_M4_AND(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)))))')dnl
dnl pushdef(`is_zero_acc_spsv_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_OR(RSB_M4_AND(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),RSB_M4_SAME(transposition,RSB_M4_TRANS_N))))')dnl
dnl
pushdef(`is_diag_d_spsv_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_NOT(RSB_M4_OR(RSB_M4_AND(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),RSB_M4_SAME(transposition,RSB_M4_TRANS_N)),RSB_M4_AND(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N))))))')dnl
dnl
dnl pushdef(`is_an_externally_backward_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_XOR(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),RSB_M4_SAME(transposition,RSB_M4_TRANS_N)))')dnl
dnl pushdef(`is_an_externally_backward_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)))')dnl
pushdef(`is_an_externally_backward_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_XOR(RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)),RSB_M4_SAME(uplo,`u')))')dnl
dnl
pushdef(`is_a_backward_kernel',is_an_externally_backward_kernel)dnl
dnl pushdef(`is_a_backward_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)))')dnl
dnl
pushdef(`block_backward',`ifelse(is_a_backward_kernel,1,`a += rows*columns',`a -= rows*columns')')dnl
pushdef(`block_forward',`ifelse(is_a_backward_kernel,1,`a -= rows*columns',`a += rows*columns')')dnl
dnl
dnl
dnl	FIXME : and so the stride x/y association
dnl
dnl pushdef(`extra_xstride',ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`incx',`0'))dnl
dnl pushdef(`extra_ystride',ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`incy',`0'))dnl
pushdef(`xstride',ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`(incx)',`1'))dnl
pushdef(`ystride',ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`(incy)',`1'))dnl
pushdef(`extra_xstride',ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`xstride',`1'))dnl
pushdef(`extra_ystride',ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`ystride',`1'))dnl
dnl
dnl	NEW:
dnl
pushdef(`transposed',ifelse(transposition,RSB_M4_TRANS_N,0,1))dnl
dnl pushdef(`transposed',dnl
dnl ifelse(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),1,eval(transposed),eval(1-transposed))dnl
dnl )dnl
dnl
dnl
pushdef(`brin',`(i*extra_ystride)')dnl
pushdef(`bcin',`(j*extra_xstride)')dnl
dnl
ifelse(transposed,`1',`dnl
dnl
dnl	block row index, block column index
dnl
pushdef(`bci',`(i*extra_xstride)')dnl
pushdef(`bri',`(j*extra_ystride)')dnl
pushdef(`bcit',`(j*extra_xstride)')dnl
pushdef(`brit',`(i*extra_ystride)')dnl
',`dnl
pushdef(`bri',`(i*extra_ystride)')dnl
pushdef(`bci',`(j*extra_xstride)')dnl
pushdef(`brit',`(j*extra_ystride)')dnl
pushdef(`bcit',`(i*extra_xstride)')dnl
')dnl
dnl
pushdef(`should_init_out_vector_before_outer_loop',`dnl
RSB_M4_OR(RSB_M4_IS_SCALING_KERNEL_MOP(mop),dnl
RSB_M4_AND(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),RSB_M4_NOT(eval(transposed))),dnl
RSB_M4_AND(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),eval(transposed)),dnl
RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry))dnl
')dnl
dnl
dnl
dnl
pushdef(`has_implementation',`dnl
RSB_M4_BXXX_KERNEL_FUNCTION_HAS_IMPLEMENTATION($@)`'dnl
')dnl
dnl
')
RSB_M4_BXXX_KERNEL_FUNCTION_HELP($@)
ifelse(RSB_M4_AND(RSB_M4_IS_SPMX_KERNEL_MOP(mop),RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal)),1,`dnl
	RSB_M4_FAKE_DIAG_IMPLICIT_MSG
')dnl
ifelse(has_implementation,`1',`dnl
',`dnl
dnl	/* or RSB_ERR_UNSUPPORTED_FEATURE ? */
	return RSB_M4_IMPOSSIBLE_IMPLEMENTATION_ERROR;`'
')dnl
dnl
ifelse(has_implementation,`1',`dnl
dnl	Comments
dnl
ifelse(RSB_M4_AND(RSB_M4_IS_SPMX_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry)),1,`dnl
	/*
ifelse(RSB_M4_want_verbose_comments,`1',`dnl
		This function assumes the matrix symmetric, and will therefore update the
		output vector in the 0,Mdim and -roff+coff,-roff+coff+Mdim range.
		If you are using several threads and kernels, you should synchronize them
		to update disjoing ranges.
')dnl
ifelse(RSB_M4_AND(RSB_M4_IS_SPMX_SCALING_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry)),1,`dnl
		The output vector zero-ing is impacted, too, so if you are using this kernel with
		recursive storage, you should care about the proper zeroing of the whole output vector.
')dnl
	*/
')dnl
dnl
dnl
dnl
dnl
ifelse(RSB_M4_AND(RSB_M4_NOT(RSB_M4_IS_COMPLEX_TYPE(type)),RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),RSB_M4_NOT(transposition,RSB_M4_TRANS_N)),1,`dnl
dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_C),1,`dnl
	/* `For non complex types, hermitian defaults to plain transposition.' */
	return RSB_M4_BCSS_KERNEL_FUNCTION(`ID',type,matrix_storage,RSB_M4_H2T_TRANSPOSITION(transposition),k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)dnl
(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCSS_KERNEL_FUNCTION(`ARGS',type,matrix_storage,RSB_M4_H2T_TRANSPOSITION(transposition),k_symmetry,rowsu,colsu,unrolling,mop,citype)));
')dnl
dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_T),1,`dnl
	/* `This kernel performs the same as its transposed', transposition -> RSB_M4_TRANSPOSE_TRANSPOSITION(transposition). */
	return RSB_M4_BCSS_KERNEL_FUNCTION(`ID',type,matrix_storage,RSB_M4_TRANSPOSE_TRANSPOSITION(transposition),k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)dnl
(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCSS_KERNEL_FUNCTION(`ARGS',type,matrix_storage,RSB_M4_TRANSPOSE_TRANSPOSITION(transposition),k_symmetry,rowsu,colsu,unrolling,mop,citype)));
')dnl
dnl
',`dnl
ifelse(RSB_M4_OR(RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(type),RSB_M4_SAME(k_symmetry,`hNEVEROCCURINGFIXME'),RSB_M4_SAME(transposition,RSB_M4_TRANS_C)),RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(type),RSB_M4_SAME(k_symmetry,`s'),RSB_M4_SAME(transposition,RSB_M4_TRANS_T))),1,`dnl
dnl
	/* `This kernel performs the same as its transposed', transposition -> RSB_M4_TRANSPOSE_TRANSPOSITION(transposition). */
	return RSB_M4_BCSS_KERNEL_FUNCTION(`ID',type,matrix_storage,RSB_M4_TRANSPOSE_TRANSPOSITION(transposition),k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)dnl
(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCSS_KERNEL_FUNCTION(`ARGS',type,matrix_storage,RSB_M4_TRANSPOSE_TRANSPOSITION(transposition),k_symmetry,rowsu,colsu,unrolling,mop,citype)));
dnl
',`dnl
dnl
ifelse(RSB_M4_AND(RSB_M4_NOT(RSB_M4_IS_COMPLEX_TYPE(type)),RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N))),1,`dnl
dnl
	/* Symmetric `transposed' reverts to symmetric `not transposed' */
	return RSB_M4_BCSS_KERNEL_FUNCTION(`ID',type,matrix_storage,RSB_M4_TRANS_N,k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)dnl
(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCSS_KERNEL_FUNCTION(`ARGS',type,matrix_storage,RSB_M4_TRANS_N,k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)));
dnl
',`dnl
dnl
dnl
ifelse(unrolling,`l',/* FIXME : l-unrolled functions are broken */)dnl
dnl
dnl	BEGIN VARIABLES DECLARATIONS
dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),`1',`dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),1,`dnl
	register rsb_coo_idx_t i=0,j=0;
',`dnl
	register rsb_coo_idx_t i=0;
')dnl
',`dnl
ifelse(is_a_backward_kernel,1,`
	register rsb_coo_idx_t i=0;
',`
	register rsb_coo_idx_t i=0,j=0;
')dnl
')dnl
	register rsb_nnz_idx_t k=0;
dnl
ifelse(RSB_M4_NOT(RSB_M4_IS_SPMV_KERNEL_MOP(mop)),`1',`dnl
ifelse(unrolling,`l',`dnl
	const register rsb_coo_idx_t columns=cpntr[1]-cpntr[0];	/* we assume that block_count >= 1 */
	const register rsb_coo_idx_t rows   =rpntr[1]-rpntr[0];	/* we assume that block_count >= 1 */
',`dnl
	const register rsb_coo_idx_t columns=b_columns,rows=b_rows;
')dnl
')dnl
')dnl
dnl
ifelse(RSB_M4_IS_READONLY_KERNEL_MOP(mop),1,`dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),`0',`dnl
	const mtype *a=VA;
')dnl
')dnl
ifelse(RSB_M4_IS_WRITEONLY_KERNEL_MOP(mop),1,`dnl
	mtype *a=VA;
')dnl
dnl
ifelse(RSB_M4_IS_RC_BIASED_KERNEL_MOP(mop),`0',`dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`0',`dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),`0',`dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),`0',`dnl
	const rsb_coo_idx_t incx=1,incy=1;`'
')dnl
')dnl
')dnl
')dnl
dnl
ifelse(RSB_M4_IS_OP_SCALING_KERNEL_MOP(mop),`1',`dnl
	const mtype alpha=*alphap;`'dnl
')dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),`1',`dnl
	const mtype beta=*betap;`'dnl
')dnl
dnl
ifelse(RSB_M4_is_transposed_spmv,1,`dnl
	const mtype *trhs = rhs+xstride*(roff-coff);`'
	mtype *tout=out+ystride*(coff-roff);`'

')dnl
ifelse(RSB_M4_IS_SPXX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
')dnl
dnl
dnl
dnl	BEGIN CONDITIONAL DEBUG SECTION
dnl
ifelse(RSB_M4_DEBUG,`1',`	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in fid\n");
',`')dnl
dnl
dnl	END CONDITIONAL DEBUG SECTION
dnl
dnl	BEGIN CONDITIONAL VECTOR SCALING
dnl
ifelse(should_init_out_vector_before_outer_loop,1,`dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),1,`dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	if(beta!=1)rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),out_dim,&beta,out,ystride);
',`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type), out_dim,&beta, out, 1);
') /* we scale the destination vector */
')dnl
ifelse(RSB_M4_IS_ZEROING_KERNEL_MOP(mop),1,`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),out_dim,NULL,out,ystride);
')dnl
')dnl
dnl
dnl	END CONDITIONAL VECTOR SCALING
dnl
dnl	BEGIN COMMON EXTERNAL LOOP BEGINNING
dnl
ifelse(RSB_M4_want_verbose_comments,`1',` /*	Outer loop. Occurs on the major dimension.	*/ ')dnl
dnl
ifelse(is_an_externally_backward_kernel,1,`
dnl	/*trick for unsigned indices: */
dnl	//RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N))
	for(Mi=Mdim-1; RSB_LIKELY((Mi+1)>0);--Mi)
	{
',`dnl
dnl
ifelse(RSB_M4_AND(RSB_M4_WANT_20110206_BOUNDED_BOX_PATCH,RSB_M4_NOT(RSB_M4_IS_SCALING_OR_ZEROING_KERNEL_MOP(mop))),1,`dnl
dnl	really, the above condition should also check for transposition! but in this way it does no wrong.
	for(Mi=br;RSB_LIKELY(Mi<bc);++Mi)
',`dnl
	for(Mi=0;RSB_LIKELY(Mi<Mdim);++Mi)
')dnl
dnl
	{
')dnl
dnl
ifelse(RSB_M4_want_verbose_comments,`1',`dnl
		/* logically,  i is the working block row, j is the working block column */
		/* physically, Mi is the working block row, mi is the working block column */
')dnl
dnl
pushdef(`colsu',ifelse(unrolling,`l',columns,colsu))dnl
pushdef(`rowsu',ifelse(unrolling,`l',rows,rowsu))dnl
pushdef(`tcolsu',ifelse(transposition,RSB_M4_TRANS_T,rowsu,colsu))dnl
pushdef(`trowsu',ifelse(transposition,RSB_M4_TRANS_T,colsu,rowsu))dnl
dnl
ifelse(RSB_M4_IS_SPXX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
pushdef(`postalphamult',`(alpha)*')dnl
',`dnl
dnl
ifelse(RSB_M4_IS_SPMX_OP_NEGATING_KERNEL_MOP(mop),1,`dnl
pushdef(`postalphamult',`(-1)*')dnl
',`dnl
pushdef(`postalphamult',`')dnl
')dnl
dnl
')dnl
dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),`1',`dnl
ifelse(transposed,`0',`dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),0,`dnl
		const mtype *a=VA;
')dnl
')dnl
ifelse(RSB_M4_OR(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),RSB_M4_AND(RSB_M4_IS_UNSYMMETRIC(k_symmetry),RSB_M4_NOT(transposed))),1,`dnl
		register mtype cacc = RSB_M4_ZERO(mtype);
dnl		mtype *outi=out+(trowsu*i*ystride);
',`dnl
')dnl
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_KERNEL_MOP(mop),`1',`dnl
ifelse(RSB_M4_is_transposed_spmv,1,`dnl
		const mtype bt=postalphamult`'trhs[(tcolsu*xstride*(Mi))];
')dnl
')dnl
		const rsb_nnz_idx_t fk=bpntr[Mi],lk=bpntr[Mi+1];
dnl
dnl
dnl	END COMMON EXTERNAL LOOP BEGINNING
dnl
dnl	BEGIN EXTERNAL LOOP VECTOR SCALING
dnl
ifelse(RSB_NOT(RSB_M4_IS_ALLOWING_ALIASING_KERNEL_MOP(mop)),1,`
ifelse(should_init_out_vector_before_outer_loop,0,`dnl
ifelse(unrolling,`l',`
ifelse(mop,`spmv_uxux',`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type), b_rows,&beta, out+rows*i, 1);/* we scale the destination vector */
')dnl
ifelse(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),1,`dnl
ifelse(mop,`spmv_uauz',`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),b_rows,NULL,out+rows*i,ystride);
')dnl
')dnl
dnl
',`dnl
dnl
ifelse(RSB_M4_IS_ZEROING_KERNEL_MOP(mop),1,`dnl
		forloop(`row',0,decr(trowsu),`out[trowsu*bri+row]=0;')
')dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),1,`dnl
		forloop(`row',0,decr(trowsu),`out[trowsu*bri+row]*=beta;')
')dnl
')dnl
')dnl
')dnl
dnl
dnl
ifelse(should_init_out_vector_before_outer_loop,0,`dnl
ifelse(RSB_M4_IS_ZEROING_KERNEL_MOP(mop),1,`dnl
		forloop(`row',0,decr(trowsu),`out[trowsu*bri+row]=0;')
')dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),1,`dnl
		forloop(`row',0,decr(trowsu),`out[trowsu*bri+row]*=beta;')
')dnl
')dnl
dnl
dnl	END EXTERNAL LOOP VECTOR SCALING
dnl
ifelse(RSB_M4_want_verbose_comments,`1',` /*		Inner loop. Occurs on the minor dimension.	*/ ')dnl
dnl
dnl	BEGIN KERNELS DEFINITION
dnl
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),`1',`dnl
dnl
dnl	BEGIN SPMV KERNEL DEF
dnl		/* SPMV KERNEL BEGINS HERE */
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),1----,`dnl
ifelse(RSB_M4_IS_SPMX_KERNEL_MOP(mop),1,`dnl
ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),1,`',`dnl
ifelse(RSB_M4_want_verbose_comments,`1',`dnl
/* Symmetric kernels should process the first block separately, if it contains `diagonal' elements. */
dnl	In register blocked code this is not the case.
')dnl
		k=fk;
		if(RSB_UNLIKELY(lk==k)) continue;/* nothing to do here */
		mi=bindx[k];
		if(mi==Mi && ((lk-k)>1) && roff==coff)	/* a `diagonal' element, and not the only one, on a diagonally positioned matrix */
		{
			const mtype *b = rhs+(tcolsu*bci);
			mtype *c=out+(trowsu*bri);
dnl			const mtype *b = rhs+(trowsu*bri);
dnl			mtype *c=out+(tcolsu*bci);
dnl
dnl	/* FIXME : THIS IS AN EXAMPLE : SHOULD INTRODUCE DIAGONAL-SUBTRACTION CODELET */
dnl
{RSB_M4_EXTRA_SYMMETRIC_DIAGONAL_FIXING_KERNEL(`row',`rows',b_rows,`column',`columns',b_columns,mtype,,mop,unrolling,transposition,RSB_M4_SYMMETRY_SWITCH(k_symmetry))}
		}
')dnl
')dnl
')dnl
dnl
ifelse(RSB_M4_AND(RSB_M4_IS_SPMX_KERNEL_MOP(mop),RSB_M4_SAME(transposed,1)),1,`dnl
ifelse(RSB_M4_want_verbose_comments,`1',`dnl
dnl		/* `Since this is a transposed kernel, we apply a correction to the output vector locations.' */
')dnl
dnl		rhs=(rhs-coff*(xstride))+roff*(xstride); out=(out-roff*(ystride))+coff*(ystride);
')dnl
dnl
dnl
ifelse(RSB_M4_IS_UNSYMMETRIC(k_symmetry),1,`dnl
ifelse(transposed,`0',`dnl
dnl
dnl	RSB_M4_EARLY_EVICT_INSTRUCTION((a+k,bindx+k))`'dnl
dnl
dnl RSB_M4_SIMPLE_LOOP_UNROLL_2S_J..
RSB_M4_SIMPLE_LOOP_UNROLL_5S(`k',`LI',`fk',`lk',`dnl
',`dnl
dnl
			`const rsb_coo_idx_t' `j_'``''LI`'=bindx[k+LI];
			`const mtype b_'``''LI`'=rhs[tcolsu*(`j_'``''LI`')*xstride];
			`const mtype a_'``''LI`'=a[k+LI];
dnl
',`dnl
',`dnl
dnl			cacc+=a[k+LI]*b_``''LI;
dnl			cacc+=a_``''LI*b_``''LI;
			``cacc+=a_''``''LI``*b_''``''LI;
',`dnl RSB_M4_EARLY_EVICT_INSTRUCTION((a+k,bindx+k))`'dnl
',RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_SMALL)
dnl
dnl	RSB_M4_EARLY_EVICT_INSTRUCTION((a+k,bindx+k))`'dnl
dnl	RSB_M4_EARLY_EVICT_INSTRUCTION((outi+k-12))`'dnl
dnl
')dnl
')dnl
dnl
dnl
ifelse(RSB_M4_IS_UNSYMMETRIC(k_symmetry),1,`dnl
ifelse(transposed,`1',`dnl
dnl
RSB_M4_SIMPLE_LOOP_UNROLL_2S_J(`k',`LI',`fk',`lk',`dnl
dnl
			`const rsb_coo_idx_t' `j_'``''LI`'=bindx[k+LI];
			`const mtype a_'``''LI`'=RSB_M4_CONJ(VA[k+LI],mtype,transposition,k_symmetry);
			`mtype c_'``''LI`'=a_``''LI*bt;
dnl
',`dnl
			tout[(tcolsu)*(`j_'``''LI`')*ystride]+=`c_'``''LI`';
',RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_SMALL)
dnl
dnl
')dnl
')dnl
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),1,`dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_HERMITIAN,`dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_C),1,`dnl
pushdef(`ntransposition',transposition)dnl
pushdef(`ttransposition',RSB_M4_TRANSPOSE_TRANSPOSITION(transposition))dnl
')dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_T),1,`dnl
pushdef(`ntransposition',transposition)dnl
pushdef(`ttransposition',RSB_M4_TRANS_C)dnl
')dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_N),1,`dnl
pushdef(`ntransposition',RSB_M4_TRANS_C)dnl
pushdef(`ttransposition',transposition)dnl
')dnl
',`dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_C),1,`dnl
pushdef(`ntransposition',transposition)dnl
pushdef(`ttransposition',transposition)dnl
',`dnl
pushdef(`ntransposition',RSB_M4_TRANSPOSE_TRANSPOSITION(transposition))dnl
pushdef(`ttransposition',RSB_M4_TRANSPOSE_TRANSPOSITION(transposition))dnl
')dnl
')dnl
dnl			// nt stay for either ntransposition or ttransposition
			k=fk;
			if(k==lk)continue;
			j=bindx[k];
			cacc += RSB_M4_CONJ(VA[k],mtype,ntransposition,k_symmetry)*rhs[tcolsu*j*xstride];
			if(roff!=coff || (j!=i))
				tout[(tcolsu)*(j)*ystride]+=RSB_M4_CONJ(VA[k],mtype,ttransposition,k_symmetry)*bt;
			++k;
dnl RSB_M4_SIMPLE_LOOP_UNROLL_2S..
RSB_M4_SIMPLE_LOOP_UNROLL_2S_J(`k',`LI',`fk+1',`lk-1',`dnl
dnl
			`const rsb_coo_idx_t' `j_'``''LI`'=bindx[k+LI];
			`const mtype b_'``''LI`'=rhs[tcolsu*(`j_'``''LI`')*xstride];
			`const mtype a_'``''LI`'=VA[k+LI];
			`mtype c_'``''LI`'=RSB_M4_CONJ_SYM(mtype,ttransposition,k_symmetry)( `a_'``''LI)*bt;
dnl			`mtype c_'``''LI`'=RSB_M4_CONJ(( `a_'``''LI),mtype,transposition,k_symmetry) *bt ;
dnl
',`dnl
			cacc += RSB_M4_CONJ_SYM(mtype,ntransposition,k_symmetry)(`a_'``''LI)*b_``''LI;
			tout[(tcolsu)*(`j_'``''LI`')*ystride]+=`c_'``''LI`';
',RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_SMALL)
			if(k<lk)
			{
				j=bindx[k];
				cacc += RSB_M4_CONJ(VA[k],mtype,ntransposition,k_symmetry)*rhs[trowsu*j*xstride];
				if(roff!=coff || (j!=i))
					tout[(tcolsu)*(j)*ystride]+=RSB_M4_CONJ(VA[k],mtype,ttransposition,k_symmetry)*bt;
				++k;
			}
popdef(`ntransposition')dnl
popdef(`ttransposition')dnl
dnl
')dnl
dnl
ifelse(RSB_M4_should_merge_value_after_inner_loop,`1',`dnl
dnl			outi[0]+=postalphamult`cacc';
			out[(trowsu*i*ystride)]+=postalphamult`cacc';
')dnl
dnl
dnl		}
dnl
dnl		/* SPMV KERNEL ENDS HERE */
popdef(`postalphamult')dnl
dnl	END SPMV KERNEL DEF
')dnl
dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),`1',`dnl
dnl	BEGIN SPSV KERNEL DEF
dnl	/* SPSV KERNEL BEGINS HERE */
dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),1,`dnl
dnl		const mtype bb_0=rhs[(trowsu*bri)];
ifelse(is_diag_d_spsv_kernel,1,`',`dnl
ifelse(RSB_M4_OR(RSB_M4_IS_SPSX_OP_SCALING_KERNEL_MOP(mop),RSB_M4_IS_SPSX_OP_SETTING_KERNEL_MOP(mop)),1,`dnl
		const mtype bb_0=rhs[(trowsu*Mi*extra_xstride)];
')dnl
')dnl
		mtype ax_0;
dnl
ifelse(is_diag_d_spsv_kernel,1,`dnl
dnl	
dnl	FIXME: missing incx, incy support here!
dnl
ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),1,`dnl
		const mtype aa=1;
',`dnl
		const mtype aa=VA[ifelse(uplo,`u',`fk',`lk-1')];
ifelse(RSB_M4_WANT_SPSM_DIAG_CHECK(),1,`dnl
		if(aa == RSB_M4_ZERO(mtype))return RSB_ERR_INVALID_NUMERICAL_DATA;
')dnl
')dnl
dnl

ifelse(RSB_M4_IS_SPSX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
dnl
dnl		out[tcolsu*bci]/=RSB_M4_CONJ(VA[bpntr[Mi+1]-1],mtype,transposition,k_symmetry);
dnl
',`dnl
dnl		out[tcolsu*bci]/=RSB_M4_CONJ(VA[bpntr[Mi+1]-1],mtype,transposition,k_symmetry);
dnl
')dnl
dnl		
		out[tcolsu*brit]/=aa;
dnl
')dnl
dnl
ifelse(is_zero_acc_spsv_kernel,1,`dnl
		ax_0=0;
',`dnl
		ax_0=out[tcolsu*bci];
')dnl
dnl
dnl
')dnl
dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),`1',`dnl
pushdef(`skip_head_row_elements',ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),`1',`0',ifelse(uplo,`u',`1',`0')))dnl
pushdef(`skip_tail_row_elements',ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),`1',`0',ifelse(uplo,`u',`0',`1')))dnl
',`dnl
pushdef(`skip_head_row_elements',0)dnl
pushdef(`skip_tail_row_elements',0)dnl
')dnl
dnl
ifelse(is_a_backward_kernel,1,`
dnl
dnl	FIXME : backward kernels are only used for SPSV, and they start with one element less
dnl
		for(k=lk-1-skip_tail_row_elements`'dnl
,a=VA+k;k+1>=fk+1+skip_head_row_elements;--k,block_forward)
dnl	/* k is the index of the block */
',`dnl
		ifelse(skip_head_row_elements,1,block_forward;)
		for(k=fk+skip_head_row_elements,mi=bindx[k];k<lk-skip_tail_row_elements  ;++k,block_forward,mi=bindx[k])
dnl	/* k is the index of the block */
')dnl
		{
ifelse(is_a_backward_kernel,1,`			const rsb_coo_idx_t mi=bindx[k];')
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_N),1,`dnl
			const mtype *b=out + (tcolsu*bci);
			mtype *c=&ax_0;
')dnl
dnl
dnl	Fixed for Hermitian k_symmetry.
dnl
ifelse(is_diag_d_spsv_kernel,1,`dnl
		out[trowsu*bri]-=RSB_M4_CONJ(*a,mtype,transposition,k_symmetry)*ax_0;
',`dnl
{RSB_M4_KERNEL_FUNCTION_BODY(`row',`rows',b_rows,`column',`columns',b_columns,mtype,,mop,unrolling,RSB_M4_SYMBOL_UNSYMMETRIC)}
')dnl
dnl
		}
dnl
ifelse(is_diag_d_spsv_kernel,1,`dnl
ifelse(RSB_M4_IS_SPSX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
		out[tcolsu*brit]*=alpha;
')dnl
')dnl
dnl
ifelse(is_diag_d_spsv_kernel,1,`',`dnl
ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),1,`',`dnl
		if(lk-fk>0)
dnl	/* if this row block was not empty */
')dnl
		{
			/* `the last element (which for a lower triangular solve is on the diagonal')*/
dnl			block_backward;
			/* Lx=y ; x_0=y_0/L_1_1  */
			mtype *c_0=out+(trowsu*bri);
ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),1,`dnl
			const mtype aa=1;
',`dnl
dnl			elements on the diagonal are real, and no conjugation is needed 
			const mtype aa=VA[ifelse(uplo,`u',`fk',`lk-1')];
ifelse(RSB_M4_WANT_SPSM_DIAG_CHECK(),1,`dnl
		if(aa == RSB_M4_ZERO(mtype))return RSB_ERR_INVALID_NUMERICAL_DATA;
')dnl
')dnl
dnl
dnl
ifelse(RSB_M4_IS_SPSX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
			*c_0 =(alpha*bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
')dnl
ifelse(RSB_M4_IS_SPSX_OP_SETTING_KERNEL_MOP(mop),1,`dnl
			*c_0=(bb_0 - ax_0)/aa;	/* ax_0 + *a * *c_0=bb_0 -> (*c_0)=(bb_0 - ax_0 )/(*a) */
')dnl
dnl
ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),1,`',`dnl
			block_forward;
')dnl
		}
')dnl
dnl
popdef(`skip_head_row_elements')dnl
popdef(`skip_tail_row_elements')dnl
dnl
dnl
dnl		/* SPSV KERNEL ENDS HERE */
dnl	END SPSV KERNEL DEF
')dnl
dnl
ifelse(RSB_M4_NOT(RSB_M4_IS_SPXX_KERNEL_MOP(mop)),`1',`dnl
dnl	BEGIN MISC KERNEL DEF
dnl
 		/* touppercase(mop) KERNEL HERE */
dnl		for(k=fk,mi=bindx[k];k<lk;++k,block_forward,mi=bindx[k]) 20120915 /*buggy loop */
		for(k=fk;k<lk;++k,block_forward)
		{
		mi=bindx[k];
		{
ifelse(mop,`scale',`dnl
			/*a=VA+indptr[(k)];*/
			const mtype *d=scale_factors+(trowsu*bri);
')dnl
ifelse(mop,`negation',`dnl
			/*a=VA+indptr[k];*/
')dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
			/*a=VA+indptr[k];*/
			mtype *local_row_sums = row_sums+(trowsu*bri);
')dnl
dnl {RSB_M4_KERNEL_FUNCTION_BODY(`row',`rows',b_rows,`column',`columns',b_columns,mtype,,mop,unrolling,RSB_M4_SYMBOL_UNSYMMETRIC)}
{RSB_M4_KERNEL_FUNCTION_BODY(`row',`rows',b_rows,`column',`columns',b_columns,mtype,,mop,unrolling,k_symmetry)}
		}
		}
dnl
dnl	END MISC KERNEL DEF
')dnl
dnl
dnl	END KERNELS DEFINITION
dnl
dnl	BEGIN COMMON EXTERNAL LOOP CLOSING
	}
dnl	END COMMON EXTERNAL LOOP CLOSING
dnl
dnl ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`1',dnl
dnl	`incx--;incy--;/* we are interested in the increment off 1 */
dnl ')dnl
dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
dnl	this check would be good for no-looped functions only!
dnl	if(columns != b_columns || rows != b_rows)return RSB_ERR_BADARGS; /* a non comprehensive check of course*/

dnl	FIXME : ONLY EXPERIMENTAL OPENMP SUPPORT
dnl
dnl
ifelse(RSB_M4_WANT_OMP_IN_KERNELS,`1',`dnl
	size_t tn;
	size_t nt;
`#'dnl
       pragma omp parallel num_threads(rsb_global_session_handle.rsb_g_threads) private(mi,Mi,k,tn,nt) 
	{
	tn = omp_get_thread_num();
	nt = omp_get_num_threads();
	/*RSB_INFO("working on %d / %d threads\n",tn,nt);*/
	//for(Mi=tn;Mi<Mdim;Mi+=nt)
	size_t ui=((Mdim/nt)*(tn+1));
	size_t li=(Mdim/nt)*tn;
	if(ui>Mdim)ui=Mdim;
dnl	#pragma omp for schedule(static,1)		/* shared L1 cache */
#pragma omp for schedule(static,(Mdim+1)/2)		/* separate L1 caches */
	for(Mi=li;RSB_LIKELY(Mi<ui);++Mi)
	{
	//RSB_INFO("row %d working on %d / %d threads\n",mi,tn,nt);
',`dnl
dnl
')dnl
dnl ifelse(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),1,`dnl
dnl 		/* should zero output block here (for efficiency) instead of function top */
dnl ')dnl
dnl
dnl		FIXME: the following is NEW, and useful also for SYMMETRIC
dnl		/* transpose.. is transposed */
dnl		/* useless for storage matrix_storage */
dnl		/*if(bpntr[Mi]==bpntr[Mi+1])continue;*/ /* empty  */
ifelse(mop,`spmv_uauz',`dnl
dnl		mtype *c=out+(rowsu*mi); /* declaration of c put here for experimental purposes */
')dnl
dnl
dnl
popdef(`is_diag_d_spsv_kernel')dnl
popdef(`tcolsu')dnl
popdef(`trowsu')dnl
popdef(`colsu')dnl
popdef(`rowsu')dnl
popdef(`transposed')dnl 1/2
dnl popdef(`transposed')dnl 2/2
popdef(`should_init_out_vector_before_outer_loop')dnl
popdef(`total_block_columns')dnl
popdef(`total_block_rows')dnl
popdef(`total_rows')dnl
popdef(`total_columns')dnl
dnl
dnl
ifelse(RSB_M4_WANT_OMP_IN_KERNELS,`1',`dnl
	}
')dnl
popdef(`mi')dnl
popdef(`Mi')dnl
popdef(`brit')dnl
popdef(`bcit')dnl
popdef(`brin')dnl
popdef(`bcin')dnl
popdef(`bri')dnl
popdef(`bci')dnl
')dnl
dnl
	return RSB_ERR_NO_ERROR;
dnl
')dnl
')dnl
dnl
')')dnl
dnl
dnl
popdef(`skip_implementation')dnl
popdef(`out_dim')dnl
popdef(`is_a_backward_kernel')dnl
popdef(`is_an_externally_backward_kernel')dnl
popdef(`is_zero_acc_spsv_kernel')dnl
popdef(`block_forward')dnl
popdef(`block_backward')dnl
popdef(`extra_xstride')dnl
popdef(`extra_ystride')dnl
dnl } /* end of fid function */
}
dnl
')dnl
dnl
')dnl
dnl
popdef(`fid')dnl
dnl
')dnl
popdef(`uplo')dnl
popdef(`want_what')dnl
popdef(`k_diagonal')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`matrix_storage')dnl
popdef(`k_symmetry')dnl
popdef(`transposition')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
')dnl
dnl
dnl
define(`RSB_M4_BCSS_MISC_KERNELS',`dnl
dnl
pushdef(`unrollings',$1)dnl
dnl
dnl	FIXED BLOCK SIZE KERNELS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
foreach(`unrolling',unrollings,`dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop)RSB_M4_IS_SPMV_KERNEL_MOP(mop),00,`dnl
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
RSB_M4_BCSS_KERNEL_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
dnl	FIXED BLOCK SIZE DISPATCHERS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop)RSB_M4_IS_SPMV_KERNEL_MOP(mop),00,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
foreach(`unrolling',unrollings,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
dnl
popdef(`unrollings')dnl
dnl	
')dnl
dnl	
dnl	
dnl
dnl
define(`RSB_M4_BCSS_SPMV_KERNELS',`dnl
dnl
pushdef(`unrollings',$1)dnl
dnl
dnl	FIXED BLOCK SIZE KERNELS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
foreach(`unrolling',unrollings,`dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
RSB_M4_BCSS_KERNEL_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
dnl	FIXED BLOCK SIZE DISPATCHERS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
ifelse(RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_FORMAT_BCXX(matrix_storage),RSB_M4_IS_SPXX_KERNEL_MOP(mop))),1,`dnl skip_register_block_dispatcher
foreach(`unrolling',unrollings,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl skip_register_block_dispatcher
')dnl
')dnl
')dnl
')dnl
dnl
dnl
popdef(`unrollings')dnl
dnl	
')dnl
dnl	
dnl	
dnl
dnl
define(`RSB_M4_BCSS_SPSV_KERNELS',`dnl
dnl
pushdef(`unrollings',$1)dnl
dnl
dnl	FIXED BLOCK SIZE KERNELS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
foreach(`unrolling',unrollings,`dnl
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
ifelse(RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry))),1,`dnl
RSB_M4_BCSS_KERNEL_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)
',`dnl
dnl	no_symm_spsx
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
dnl	FIXED BLOCK SIZE DISPATCHERS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
ifelse(RSB_M4_IS_SPSV_KERNEL_MOP(mop),1,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
ifelse(RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_FORMAT_BCXX(matrix_storage),RSB_M4_IS_SPXX_KERNEL_MOP(mop))),1,`dnl skip_register_block_dispatcher
foreach(`unrolling',unrollings,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
ifelse(RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry))),1,`dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
',`dnl
dnl	no_symm_spsx
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl skip_register_block_dispatcher
')dnl
')dnl
')dnl
')dnl
dnl
dnl
popdef(`unrollings')dnl
dnl	
')dnl
dnl	
dnl	
dnl
dnl
dnl
define(`RSB_M4_BCSS_KERNELS',`dnl
dnl
pushdef(`unrollings',$1)dnl
dnl
dnl	FIXED BLOCK SIZE KERNELS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`unrolling',unrollings,`dnl
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
RSB_M4_BCSS_KERNEL_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,rowsu,colsu,unrolling,mop,citype,k_diagonal,uplo)
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
dnl	FIXED BLOCK SIZE DISPATCHERS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`k_diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`matrix_storage',RSB_M4_BCSS_FORMATS,`dnl
foreach(`unrolling',unrollings,`dnl
foreach(`k_symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
RSB_M4_BCSS_KERNEL_SIZE_DISPATCH_FUNCTION(`all',type,matrix_storage,transposition,k_symmetry,unrolling,,,mop,citype,k_diagonal,uplo)
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
dnl
dnl
popdef(`unrollings')dnl
dnl	
')dnl
dnl	
dnl	
dnl
