dnl
dnl	@author: Michele Martone
dnl
/*!
 @file
 @brief
 Performance kernels dispatching code, for each type, submatrix size, operation.
 For block coordinates format.
 Kernels unrolled, with no loops, for only user-specified blockings.
 */
dnl
include(`rsb_misc.m4')dnl
RSB_M4_HEADER_MESSAGE()dnl
RSB_M4_HEADER_EXTRA_DECLARATIONS()dnl
include(`rsb_krnl_bcss_macros.m4')dnl
include(`rsb_krnl_vb_macros.m4')dnl FIXME : RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME
dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_BCOO_SPMV_KERNELS',`dnl
dnl
pushdef(`unrollings',$1)dnl
dnl
dnl	FIXED BLOCK SIZE KERNELS :
dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`matrix_storage',RSB_M4_BCOO_FORMATS,`dnl
foreach(`unrolling',unrollings,`dnl
dnl ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),
ifelse(1,1,`dnl
foreach(`diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
foreach(`symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
ifelse(RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(symmetry))),1,`dnl
RSB_M4_BCOO_KERNEL_FUNCTION(`all',type,matrix_storage,transposition,symmetry,rowsu,colsu,unrolling,mop,citype,diagonal,uplo)
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
dnl ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,..
ifelse(1,1,`dnl
foreach(`matrix_storage',RSB_M4_BCOO_FORMATS,`dnl
foreach(`unrolling',unrollings,`dnl
foreach(`symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
foreach(`diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
foreach(`uplo',RSB_M4_MATRIX_UPLO_TYPES,`dnl
ifelse(RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_FORMAT_BCXX(matrix_storage),RSB_M4_IS_SPXX_KERNEL_MOP(mop))),1,`dnl skip_register_block_dispatcher
ifelse(RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(symmetry))),1,`dnl
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`all',type,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo)
',`dnl
dnl	no_symm_spsx
')dnl
')dnl skip_register_block_dispatcher
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
define(`RSB_M4_BCOO_KERNEL_FUNCTION',`dnl
dnl
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl	
pushdef(`transposition',$4)dnl	
pushdef(`symmetry',$5)dnl	
pushdef(`b_rows',$6)dnl		block rows
pushdef(`b_columns',$7)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t ')dnl integer type (for indices)
pushdef(`unrolling',$8)dnl	
pushdef(`mop',$9)dnl	
pushdef(`citype',$10)dnl	
pushdef(`diagonal',$11)dnl	
pushdef(`uplo',$12)dnl	
dnl
pushdef(`total_columns',ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`Mdim',`mdim'))dnl
pushdef(`total_rows',ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`mdim',`Mdim'))dnl
pushdef(`out_dim',ifelse(transposition,RSB_M4_TRANS_N,total_rows,total_columns))dnl
pushdef(`fid',RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME(mtype,matrix_storage,transposition,symmetry,b_rows,b_columns,unrolling,mop,citype,diagonal,uplo))dnl
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
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_N),1,`dnl
pushdef(`tmi',mi)dnl
pushdef(`tMi',Mi)dnl
')dnl
ifelse(RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)),1,`dnl
pushdef(`tmi',Mi)dnl
pushdef(`tMi',mi)dnl
')dnl
dnl
ifelse(RSB_M4_IS_SPXX_OP_SCALING_KERNEL_MOP(mop),1,`dnl
pushdef(`postmult',`(alpha)*')dnl
',`dnl
dnl
ifelse(RSB_M4_IS_SPMX_OP_NEGATING_KERNEL_MOP(mop),1,`dnl
pushdef(`postmult',`(-1)*')dnl
',`dnl
pushdef(`postmult',`')dnl
')dnl
dnl
')dnl
dnl
pushdef(`ttransposition',`RSB_M4_TRANSPOSE_TRANSPOSITION(transposition)')dnl
pushdef(`htransposition',`ifelse(symmetry,RSB_M4_SYMBOL_HERMITIAN,`RSB_M4_H2T_TRANSPOSITION(transposition)',transposition)')dnl
dnl
pushdef(`tsymmetry',`ifelse(symmetry,RSB_M4_SYMBOL_HERMITIAN,`RSB_M4_TRANSPOSE_SYMMETRY(symmetry)',symmetry)')dnl
dnl
pushdef(`toskipbecauseofsymmetry',`RSB_M4_AND(RSB_M4_IS_SPMX_KERNEL_MOP(mop),RSB_M4_NOT(RSB_M4_IS_COMPLEX_TYPE(mtype)),RSB_M4_IS_NOT_UNSYMMETRIC(symmetry),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)))')dnl
pushdef(`no_spsx_symm',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_IS_NOT_UNSYMMETRIC(symmetry))')`'dnl
pushdef(`gen_no_body',`RSB_M4_OR(no_spsx_symm,0)')`'dnl
dnl
dnl
ifelse(RSB_M4_ARE_KERNEL_GENERATION_PARMS_ALLOWED(want_what,mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo),`1',`dnl
dnl
ifelse(want_what,`DOC',`dnl
	/*  TODO */
')dnl
ifelse(want_what,`all',`dnl
dnl
ifelse(RSB_M4_IS_FORMAT_BCOO(matrix_storage),1,`dnl
rsb_err_t fid`'dnl
RSB_M4_BCOO_KERNEL_FUNCTION(`ARGS',mtype,matrix_storage,transposition,symmetry,b_rows,b_columns,unrolling,mop,citype,diagonal,uplo)dnl
')dnl
ifdef(`ONLY_WANT_HEADERS',`;
',`
dnl /* begin of fid function */
RSB_M4_BCOO_KERNEL_FUNCTION(`BODY',mtype,matrix_storage,transposition,symmetry,b_rows,b_columns,unrolling,mop,citype,diagonal,uplo)dnl
')dnl
')dnl
dnl
ifelse(want_what,`ID',`dnl
fid`'dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo)`'dnl
')dnl
dnl
dnl
ifelse(want_what,`BODY',`dnl
dnl
{
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`0',`dnl
pushdef(`incx',`1')dnl
pushdef(`incy',`1')dnl
')dnl
RSB_M4_BXXX_KERNEL_FUNCTION_HELP($@)
dnl
ifelse(RSB_M4_AND(RSB_M4_IS_SPMX_KERNEL_MOP(mop),RSB_M4_IS_DIAGONAL_IMPLICIT(diagonal)),1,`dnl
	RSB_M4_FAKE_DIAG_IMPLICIT_MSG
')dnl
dnl
ifelse(toskipbecauseofsymmetry,1,`dnl
dnl
	/* Symmetric `transposed' reverts to symmetric `not transposed' */
	return RSB_M4_BCOO_KERNEL_FUNCTION(`ID',mtype,matrix_storage,RSB_M4_TRANS_N,symmetry,b_rows,b_columns,unrolling,mop,citype,diagonal,uplo)dnl
(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCOO_KERNEL_FUNCTION(`ARGS',mtype,matrix_storage,RSB_M4_TRANS_N,symmetry,b_rows,b_columns,unrolling,mop,citype,diagonal,uplo)));
dnl
')dnl
dnl
dnl	Common Variables Declaration Section BEGIN
dnl	
ifelse(gen_no_body,`0',`dnl 
ifelse(toskipbecauseofsymmetry,0,`dnl
dnl
dnl	The i,j type has to be the same as the arrays one.
dnl	If not, mismatch on the copied bytes will occur.
ifelse(RSB_M4_AND(RSB_M4_NOT(RSB_M4_IS_RC_BIASED_KERNEL_MOP(mop)),RSB_M4_NOT(RSB_M4_AND(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),RSB_M4_NOT(RSB_M4_IS_NOT_UNSYMMETRIC(symmetry))))),`1',`dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	register rsb_coo_idx_t i=0,j=0;
',`dnl
	register citype i=0,j=0;
dnl 20110227 if declaring short indices, we should care about proper conversion
')dnl
	const citype *IA=(const citype*)bpntr, *JA=(const citype*)bindx;
dnl
',`dnl
dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_N),`0',`dnl
	const citype *JA=(const citype*)bindx;
	register citype j=0;
',`dnl
	const citype *IA=(const citype*)bpntr;
	register citype i=0;	
')dnl
')dnl
dnl ifelse(mop,`scale',`',`dnl
dnl ')dnl	20121005 shall change this condition when enabling transpose scale as well
	register rsb_nnz_idx_t n=0;
ifelse(RSB_M4_IS_OP_SCALING_KERNEL_MOP(mop),`1',`dnl
	const mtype alpha=*alphap;`'dnl
')dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),`1',`dnl
	const mtype beta=*betap;`'dnl
')dnl
dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`0',`dnl
	dnl const rsb_coo_idx_t incx=1,incy=1;`'
')dnl
dnl
dnl
dnl
dnl	Common Variables Declaration Section END
dnl	
dnl	SPMV Section BEGIN
dnl	
ifelse(RSB_M4_AND(RSB_M4_IS_SPMX_KERNEL_MOP(mop)),1,`dnl
dnl

dnl
dnl
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(symmetry),1,`dnl
	const mtype *trhs = rhs+incx*(roff-coff);`'// symmetry
	mtype *tout=out+incy*(coff-roff);`'

')dnl
dnl
ifelse(RSB_M4_IS_ZEROING_KERNEL_MOP(mop),1,`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype),out_dim,NULL,out,incy);
')dnl
dnl
ifelse(RSB_M4_IS_SCALING_KERNEL_MOP(mop),1,`dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),1,`dnl
	if(beta!=1)rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype),out_dim,&beta,out,ystride);
',`dnl
	rsb__cblas_Xscal(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype), out_dim,&beta, out, 1);
')dnl
')dnl
dnl
ifelse(transposition,RSB_M4_TRANS_N,`dnl
',`dnl
ifelse(RSB_M4_IS_UNSYMMETRIC(symmetry),1,`dnl
	rhs=(rhs-coff*(incx))+roff*(incx);
	out=(out-roff*(incy))+coff*(incy);
')dnl
')dnl
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(symmetry),1,`dnl
	if(roff==coff)
')dnl
dnl
dnl
ifelse(RSB_M4_IS_UNSYMMETRIC(symmetry),1,`dnl
dnl
ifelse(1,1,`dnl
dnl
dnl	RSB_M4_SIMPLE_LOOP_UNROLL..
	RSB_M4_SIMPLE_LOOP_UNROLL_5S(`n',`LI',`0',`nnz',`dnl
',`dnl
	i=IA[n+LI]; j=JA[n+LI];
	out[tMi*incy]+=`'postmult`'RSB_M4_CONJ(VA[n+LI],mtype,transposition,RSB_M4_SYMBOL_UNSYMMETRIC)*rhs[tmi*incx];
dnl
',`',`',`RSB_M4_EARLY_EVICT_INSTRUCTION((IA+n,JA+n,VA+n))`'dnl
',RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_SMALL)
dnl
',`dnl
dnl
RSB_M4_SIMPLE_LOOP_UNROLL_5S(`n',`LI',`0',`nnz',`dnl
dnl
dnl
',`dnl
dnl
			`const rsb_coo_idx_t' `i_'``''LI`'=IA[n+LI];
			`const rsb_coo_idx_t' `j_'``''LI`'=JA[n+LI];
			`const mtype b_'``''LI`'=rhs[tmi``_''LI`'*incx];
			`const mtype a_'``''LI`'=VA[n+LI];
dnl
',`dnl
			if(tMi``_''0`'== tMi``_''eval(RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_MEDIUM-1)`')
			{
				mtype cacc = RSB_M4_ZERO(mtype);
forloop(`_LI_',0,decr(RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_MEDIUM),`dnl
				cacc+=`'postmult`'RSB_M4_CONJ(`a_'``''_LI_,mtype,transposition,RSB_M4_SYMBOL_UNSYMMETRIC)`*b_'``''_LI_;
')dnl
			out[tMi``_''0`'*incy]+=cacc;
`'dnl
			}
			else
			{
',`dnl
				out[tMi``_''LI`'*incy]+=`'postmult`RSB_M4_CONJ(a``_''``''LI`',mtype,transposition,RSB_M4_SYMBOL_UNSYMMETRIC)'`*b_'``''LI;
',`dnl
			}
',RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_MEDIUM)
dnl
')dnl
dnl
',`dnl
dnl
	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
		i=IA[n];
		j=JA[n];
dnl		assert(i< Mdim);
dnl		assert(j< mdim);
		out[tMi*incy]+=`'postmult`'RSB_M4_UIM_CONJ(VA[n],mtype,transposition,symmetry)*rhs[tmi*incx];
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(symmetry),1,`dnl
		if(RSB_LIKELY(tMi!=tmi))
			out[tmi*incy]+=`'postmult`'RSB_M4_CIM_CONJ(VA[n],mtype,transposition,symmetry)*rhs[tMi*incx];
')dnl
dnl
	}
dnl
')dnl
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(symmetry),1,`dnl
	if(roff!=coff)
	RSB_M4_SIMPLE_LOOP_UNROLL(`n',`LI',`0',`nnz',`dnl
		i=IA[n+LI];
		j=JA[n+LI];
dnl		assert(i< Mdim);
dnl		assert(j< mdim);
ifelse(transposition,RSB_M4_TRANS_N,`dnl
		out[Mi*incy]+=`'postmult`'RSB_M4_UIM_CONJ(VA[n+LI],mtype,transposition,symmetry)*rhs[mi*incx];
		tout[mi*incy]+=`'postmult`'RSB_M4_CIM_CONJ(VA[n+LI],mtype,transposition,symmetry)*trhs[Mi*incx];
',`dnl
		tout[tMi*incy]+=`'postmult`'RSB_M4_UIM_CONJ(VA[n+LI],mtype,transposition,symmetry)*trhs[tmi*incx];
		out[tmi*incy]+=`'postmult`'RSB_M4_CIM_CONJ(VA[n+LI],mtype,transposition,symmetry)*rhs[tMi*incx];
')dnl
dnl
	',RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_SMALL)
')dnl
dnl
dnl	BEGIN CONDITIONAL DEBUG SECTION
dnl
ifelse(RSB_M4_DEBUG,`1',`	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in fid\n");
',`')dnl
dnl
dnl	END CONDITIONAL DEBUG SECTION
dnl
	return RSB_ERR_NO_ERROR;
')dnl
dnl
dnl	SPMV Section END
dnl
')dnl
')dnl gen_no_body
dnl
dnl	SPSV Section BEGIN
dnl
ifelse(RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop)),1,`dnl
ifelse(no_spsx_symm,`1',`dnl no_spsx_symm
dnl	FIXME: shall trim code generator earlier on.
	return RSB_M4_IMPOSSIBLE_IMPLEMENTATION_ERROR;
',`dnl no_spsx_symm
dnl
dnl
dnl	FIXME: and roff and coff ?
dnl
dnl
pushdef(`is_an_externally_backward_kernel',`RSB_M4_AND(RSB_M4_IS_SPSX_KERNEL_MOP(mop),RSB_M4_XOR(RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)),RSB_M4_SAME(uplo,`u')))')dnl
pushdef(`is_vector_updating_spsv',RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_N)))dnl
dnl
	rsb_coo_idx_t ii;
ifelse(is_an_externally_backward_kernel,1,`
	for(n=nnz-1,ii=Mdim-1;RSB_LIKELY(ii+1>0) ;--ii)
',`dnl
	for(n=0,ii=0;RSB_LIKELY(ii<Mdim);++ii)
')dnl
	{
		mtype ax;
ifelse(is_vector_updating_spsv,1,`dnl
ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(diagonal),1,`dnl
dnl	..
',`dnl
dnl		const mtype aa;
		mtype aa;
ifelse(RSB_M4_WANT_SPSM_DIAG_CHECK,1,`dnl
		if(n>=nnz)return RSB_ERR_INVALID_NUMERICAL_DATA;
')dnl
		aa=VA[n];
ifelse(RSB_M4_WANT_SPSM_DIAG_CHECK,1,`dnl
		if(VA[n]==RSB_M4_ZERO(mtype))return RSB_ERR_INVALID_NUMERICAL_DATA;
')dnl
ifelse(is_an_externally_backward_kernel,1,`
		n--;
',`dnl
		n++;
')dnl
		out[ii*incy]/=aa;
')dnl
		ax=out[ii*incy];
',`dnl
		ax=0;
')dnl
ifelse(is_an_externally_backward_kernel,1,`
		for(;RSB_LIKELY(n+1>0);--n)
',`dnl
		for(;RSB_LIKELY(n<nnz);++n)
')dnl
		{
			i=IA[n];
			j=JA[n];
ifelse(is_vector_updating_spsv,1,`dnl
			if(RSB_UNLIKELY(!(i==ii )))
',`dnl
			if(RSB_UNLIKELY(!(i==ii && j!=i)))
')dnl
				break;
ifelse(is_vector_updating_spsv,1,`dnl
			out[j*incy]-=RSB_M4_CONJ(VA[n],mtype,transposition,symmetry)*ax;
',`dnl
			ax += RSB_M4_CONJ(VA[n],mtype,transposition,symmetry)*out[j*incy];
')dnl
		}

ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(diagonal),1,`dnl
ifelse(is_vector_updating_spsv,1,`dnl
		out[ii*incy]=(`'postmult`'out[ii*incy]);
',`dnl
		out[ii*incy]=(`'postmult`'out[ii*incy]-ax);
')dnl
',`dnl
dnl
dnl	FIXME: goto err is illegal for nnz=0 ...
dnl
dnl		if(!(i==ii && i==j))
dnl			goto err;
ifelse(is_vector_updating_spsv,1,`dnl
		out[ii*incy]=(`'postmult`'out[ii*incy]);
',`dnl
ifelse(RSB_M4_WANT_SPSM_DIAG_CHECK,1,`dnl
		if(n==nnz || VA[n]==RSB_M4_ZERO(mtype))return RSB_ERR_INVALID_NUMERICAL_DATA;
')dnl
		out[ii*incy]=(`'postmult`'out[ii*incy]-ax)/VA[n];
ifelse(is_an_externally_backward_kernel,1,`dnl
		--n;
',`dnl
		++n;
')dnl
')dnl
')dnl
	}
dnl
dnl	BEGIN CONDITIONAL DEBUG SECTION
dnl
ifelse(RSB_M4_DEBUG,`1',`	if(rsb__getenv_int_t("RSB_VERBOSE_KERNELS",0))RSB_STDOUT("in fid\n");
',`')dnl
dnl
dnl	END CONDITIONAL DEBUG SECTION
dnl
	return RSB_ERR_NO_ERROR;
dnl
popdef(`is_an_externally_backward_kernel')dnl
popdef(`is_vector_updating_spsv')dnl
dnl
')`'dnl no_spsx_symm
')dnl
dnl
dnl	SPSV Section END
dnl
dnl
dnl ifelse(RSB_M4_NOT(RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP(mop)),1,`dnl
dnl 	return RSB_ERR_UNIMPLEMENTED_YET;
dnl ')dnl
dnl
ifelse(mop,`scale',`dnl
	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
	dnl
dnl	FIXME: what about hermitian ?
dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_N),1,`dnl
		i=IA[n];
		VA[n]*=scale_factors[i];
',`dnl
		j=JA[n];
dnl		i=IA[n];
dnl		VA[n]*=scale_factors[i];
		VA[n]*=scale_factors[j];
')dnl
dnl
	}
	return RSB_ERR_NO_ERROR;
')dnl
dnl
ifelse(RSB_M4_IS_ACC_WRITING_KERNEL_MOP(mop),`1',`dnl
	dnl
	dnl	TODO: do we need vector blank ?
	dnl
	for(n=0;RSB_LIKELY(n<nnz);++n)
	{
dnl
ifelse(RSB_M4_IS_UNSYMMETRIC(symmetry),1,`dnl
dnl
ifelse(RSB_M4_SAME(transposition,RSB_M4_TRANS_N),1,`dnl
		i=IA[n];
ifelse(mop,`infty_norm',`dnl
		row_sums[roff+i]+=RSB_M4_ABS(mtype,VA[n]);
')dnl
ifelse(mop,`rowssums',`dnl
		row_sums[roff+i]+=VA[n];
')dnl
',`dnl
		j=JA[n];
ifelse(mop,`infty_norm',`dnl
		row_sums[coff+j]+=RSB_M4_ABS(mtype,VA[n]);
')dnl
ifelse(mop,`rowssums',`dnl
		row_sums[coff+j]+=VA[n];
')dnl
')dnl
')dnl
dnl
dnl
ifelse(RSB_M4_IS_UNSYMMETRIC(symmetry),0,`dnl
dnl
		i=IA[n];
		j=JA[n];
dnl
ifelse(mop,`infty_norm',`dnl
		row_sums[roff+i]+=RSB_M4_ABS(mtype,VA[n]);
')dnl
ifelse(mop,`rowssums',`dnl
		row_sums[roff+i]+=VA[n];
')dnl
		if( roff+i != coff+j )
ifelse(mop,`infty_norm',`dnl
			row_sums[coff+j]+=RSB_M4_ABS(mtype,VA[n]);
')dnl
ifelse(mop,`rowssums',`dnl
			row_sums[coff+j]+=VA[n];
')dnl
')dnl
dnl
	}
	return RSB_ERR_NO_ERROR;
')dnl
dnl
ifelse(RSB_M4_IS_STRIDED_KERNEL_MOP(mop),`0',`dnl
popdef(`incx')dnl
popdef(`incy')dnl
')dnl
dnl
}
dnl } /* end of fid function */
dnl
')dnl
dnl
')dnl
dnl
dnl
popdef(`gen_no_body')dnl
popdef(`no_spsx_symm')dnl
popdef(`toskipbecauseofsymmetry')dnl
popdef(`htransposition')dnl
popdef(`ttransposition')dnl
popdef(`tsymmetry')dnl
popdef(`postmult')dnl
popdef(`tmi')dnl
popdef(`tMi')dnl
popdef(`mi')dnl
popdef(`Mi')dnl
popdef(`total_columns')dnl
popdef(`total_rows')dnl
popdef(`out_dim')dnl
popdef(`fid')dnl
dnl
popdef(`uplo')dnl
popdef(`diagonal')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`unrolling')dnl
popdef(`itype')dnl
popdef(`b_columns')dnl
popdef(`b_rows')dnl
popdef(`symmetry')dnl
popdef(`transposition')dnl
popdef(`matrix_storage')dnl
popdef(`mtype')dnl
popdef(`want_what')dnl
')dnl
dnl
dnl
define(`RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION',`dnl
dnl
pushdef(`want_what',$1)dnl
pushdef(`mtype',$2)dnl
pushdef(`matrix_storage',$3)dnl
pushdef(`transposition',$4)dnl
pushdef(`symmetry',$5)dnl
pushdef(`unrolling',$6)dnl	
dnl pushdef(`b_rows',$7)dnl		block rows
dnl pushdef(`b_columns',$8)dnl	block columns
pushdef(`itype',`rsb_coo_idx_t ')dnl integer type (for indices)
pushdef(`mop',`$9')dnl
pushdef(`citype',`$10')dnl
pushdef(`diagonal',`$11')dnl
pushdef(`uplo',$12)dnl
dnl
dnl
dnl
ifelse(RSB_M4_ARE_KERNEL_GENERATION_PARMS_ALLOWED(want_what,mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo),`1',`dnl
dnl
ifelse(want_what,`DOC',`dnl
	/*  TODO */
')dnl
dnl
ifelse(want_what,`all',`dnl
dnl `/* This code is intended for a block compressed sparse stripe matrix. */'
ifdef(`ONLY_WANT_HEADERS',`dnl
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`function_declaration',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo)
',`dnl
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`function_definition',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo)
')dnl
dnl
dnl
dnl
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
rsb_err_t RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,symmetry,unrolling,mop,citype,diagonal,uplo)dnl
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo)
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`BODY',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo)
')dnl
dnl
ifelse(want_what,`function_declaration',`dnl
rsb_err_t RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,symmetry,unrolling,mop,citype,diagonal,uplo)dnl
RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo);dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
dnl
dnl
pushdef(`matrix_structs',`const itype Mdim,const itype mdim,const citype * RSB_M4_RESTRICT bindx,const rsb_nnz_idx_t * RSB_M4_RESTRICT bpntr,const rsb_nnz_idx_t *RSB_M4_RESTRICT indptr,const rsb_coo_idx_t * RSB_M4_RESTRICT rpntr,const rsb_coo_idx_t * RSB_M4_RESTRICT cpntr,const rsb_coo_idx_t br,const rsb_coo_idx_t bc,const rsb_coo_idx_t roff,const rsb_coo_idx_t coff,const rsb_flags_t flags,const rsb_nnz_idx_t nnz')dnl
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
	 * Select kernel function for operation "mop".
ifelse(RSB_M4_MAXIMAL_CONFIGURED_BLOCK_SIZE,`1',`',`dnl
	 *
	 * Since this is strictly blocked code, you should allow the rhs and the out
	 * vector to accept a small overflow not bigger, respectively, than
	 *       mod(blockrows-mod(matrixrows,blockrows),blockrows)
	 * and
	 *       mod(blockcols-mod(matrixcols,blockcols),blockcols)
')dnl
	 *
	 * \return \rsb_errval_inp_param_msg
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
ifelse(RSB_M4_MAXIMAL_CONFIGURED_BLOCK_SIZE,`1',`dnl
dnl 		non-unitary blocking error shall get caught earlier than this, so no check needed
',`dnl
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
dnl #if RSB_EXPERIMENTAL_WANT_PURE_BCOO
ifelse(RSB_M4_WANT_20110206_BOUNDED_BOX_PATCH,1,`dnl
dnl 20110206	set the following 
		columns = rows=1;
',`dnl
dnl 20110206	and commented the following 
		columns=bc,rows=br;
')dnl
		RSB_ASSERT(columns==1 && rows==1);
')dnl
dnl #else
dnl 		columns = rows=1;
dnl #endif

dnl
ifelse(RSB_M4_IS_FORMAT_BCOO(matrix_storage),1,`dnl
pushdef(`args',`RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_BCOO_KERNEL_SIZE_DISPATCH_FUNCTION(`ARGS',mtype,matrix_storage,transposition,symmetry,unrolling,,,mop,citype,diagonal,uplo))')dnl
ifelse(RSB_M4_MAXIMAL_CONFIGURED_BLOCK_SIZE,`1',`dnl
	errval = RSB_M4_BCOO_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,symmetry,1,1,unrolling,mop,citype,diagonal,uplo)( args );
',`dnl
switch(rows)
{
foreach(`rowsu',RSB_M4_ROWS_UNROLL,`dnl
	case rowsu:
	{switch(columns)
	{
foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`dnl
		case colsu:/* rowsu colsu matrix_storage */
		errval = RSB_M4_BCOO_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,symmetry,rowsu,colsu,unrolling,mop,citype,diagonal,uplo)( args );
		break;
')dnl
	default: 
#ifdef RSB_WANT_LOOPING_KERNELS 
		errval = RSB_M4_BCOO_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,symmetry,rowsu,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,citype,diagonal,uplo)( args );
#else /* RSB_WANT_LOOPING_KERNELS  */
	errval = RSB_ERR_UNSUPPORTED_OPERATION;
#endif /* RSB_WANT_LOOPING_KERNELS  */
	}}
	break;
')dnl
	default:
#ifdef RSB_WANT_LOOPING_KERNELS 
		errval = RSB_M4_BCOO_KERNEL_FUNCTION(`ID',mtype,matrix_storage,transposition,symmetry,RSB_M4_ROWS_FALLBACK_UNROLL,RSB_M4_COLUMNS_FALLBACK_UNROLL,`l',mop,citype,diagonal,uplo)( args );
#else /* RSB_WANT_LOOPING_KERNELS */
	errval = RSB_ERR_UNSUPPORTED_OPERATION;
#endif /* RSB_WANT_LOOPING_KERNELS */
};
')dnl
popdef(`args')dnl
')dnl
	dnl errval = RSB_ERR_UNSUPPORTED_TYPE;
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
popdef(`symmetry')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`unrolling')dnl
popdef(`diagonal')dnl
popdef(`want_what')dnl
popdef(`uplo')dnl
')dnl
dnl
dnl
dnl
