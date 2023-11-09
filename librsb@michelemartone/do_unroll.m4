dnl
dnl
dnl	@author: Michele Martone
dnl
dnl	This file contains loop unrolled kernels for many operations.
dnl
dnl	TODO  : eliminate all traces of "negation" kernel, as it is not legal anymore.
dnl	FIXME : the only kernel working with transposition is spmv.
dnl
include(`rsb_misc.m4')dnl
include(`wisdom.m4')dnl
dnl 
dnl 
dnl  Follows an m4 documentation macro
dnl 
dnl
dnl 
dnl	---------------------------------------------------------------------------	dnl
dnl				Function name and declaration macros
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl	RSB_M4_KERNEL_FUNCTION_NAME(TYPE,ROWS_UNROLL,COLS_UNROLL,UNROLLING,MOP,TRANSPOSITION)
dnl	-----------------------------------------------------------------------
dnl	Expands to the function name of the kernel specified by the macro arguments.
dnl
dnl
define(`RSB_M4_KERNEL_FUNCTION_NAME',dnl
pushdef(`type',$1)dnl
pushdef(`rows_unroll',$2)dnl
pushdef(`cols_unroll',$3)dnl
pushdef(`unrolling',$4)dnl
pushdef(`mop',$5)dnl
pushdef(`transposition',$6)dnl
pushdef(`citype',$7)dnl
dnl
dnl	FIXME : using symbolic names here is troublesome !
dnl	
dnl pushdef(`unrolling',ifelse($2,`l',`l',`'))dnl FIXME : this is a temporary fix (setting to `' instead of `u')
``rsb_m'`$5'`_'RSB_M4_TYPE_CODE(`type')`_'dnl
RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE($7)`_'dnl
`_r'`$2'`_c'`$3'ifelse(`$4',l,_l,)`'RSB_M4_TRANSPOSITION_CODE(transposition)'dnl
dnl popdef(`unrolling')
popdef(`citype')dnl
popdef(`transposition')dnl
popdef(`rows_unroll')dnl
popdef(`cols_unroll')dnl
popdef(`unrolling')dnl
popdef(`type')dnl
popdef(`mop')dnl
)dnl
dnl
dnl
dnl	RSB_M4_KERNEL_FUNCTION_ARGS(TYPE,UNROLLING,MOP)
dnl	-----------------------------------------------
dnl	Expands to the function arguments of the kernel specified by the macro arguments.
dnl
define(`RSB_M4_KERNEL_FUNCTION_ARGS',`dnl
dnl
dnl
pushdef(`mtype',`$1')dnl
pushdef(`itype',`rsb_coo_idx_t')dnl
pushdef(`unrolling',ifelse($2,`l',`l',`'))dnl FIXME : this is a temporary fix (setting to `' instead of `u')
pushdef(`optype',$3)dnl
dnl	not fully unrolled :
ifelse(unrolling`'optype, `spmm_az',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns, const itype bstride, const itype cstride, const itype nrhs)')dnl
ifelse(unrolling`'optype, `l'`spmm_az',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns, const itype bstride, const itype cstride, const itype nrhs)')dnl
dnl	fully unrolled :
ifelse(unrolling`'optype, `infty_norm',dnl
`(const mtype *a, mtype *local_row_sums, const itype dummy_rows, const itype dummy_columns)')dnl
ifelse(unrolling`'optype,`l'`infty_norm',dnl
`(const mtype *a, mtype *local_row_sums, const itype rows, const itype columns)')dnl
ifelse(unrolling`'optype, `rowssums',dnl
`(const mtype *a, mtype *local_row_sums, const itype dummy_rows, const itype dummy_columns)')dnl
ifelse(unrolling`'optype,`l'`rowssums',dnl
`(const mtype *a, mtype *local_row_sums, const itype rows, const itype columns)')dnl
ifelse(unrolling`'optype, `negation',dnl
`(mtype *a, const itype dummy_rows, const itype dummy_columns)')dnl
ifelse(unrolling`'optype,`l'`negation',dnl
`(mtype *a, const itype rows, const itype columns)')dnl
ifelse(unrolling`'optype, `spmv_uauz',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype dummy_rows, const itype dummy_columns)')dnl
ifelse(unrolling`'optype,`l',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns)')dnl
ifelse(unrolling`'optype, `spsv_uxua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype dummy_rows, const itype dummy_columns)')dnl
ifelse(unrolling`'optype,`l'`spsv_uxua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns)')dnl
ifelse(unrolling`'optype, `spmv_unua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype dummy_rows, const itype dummy_columns)')dnl
ifelse(unrolling`'optype,`l'`spmv_unua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns)')dnl
ifelse(unrolling`'optype, `spmv_uaua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype dummy_rows, const itype dummy_columns)')dnl
ifelse(unrolling`'optype,`l'`spmv_uaua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns)')dnl
ifelse(unrolling`'optype, `spmv_uxua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype dummy_rows, const itype dummy_columns, const mtype *alphap)')dnl
ifelse(unrolling`'optype,`l'`spmv_uxua',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns, const mtype * alphap)')dnl
ifelse(unrolling`'optype, `spmv_uxux',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype dummy_rows, const itype dummy_columns, const mtype *alphap, const mtype *betap)')dnl
ifelse(unrolling`'optype,`l'`spmv_uxux',dnl
`(const mtype *a, const mtype *b, mtype *c, const itype rows, const itype columns, const mtype * alphap, const mtype *betap)')dnl
ifelse(unrolling`'optype,`scale',dnl
`(mtype *a, const mtype *d, const itype rows, const itype columns /* dummy rows and columns */)')dnl
dnl `(mtype *a, const mtype *d)')dnl
ifelse(unrolling`'optype,`l'`scale',dnl
`(mtype *a, const mtype *d, const itype rows, const itype columns)')dnl
popdef(`optype')dnl
popdef(`unrolling')dnl
popdef(`itype')dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl	RSB_M4_KERNEL_FUNCTION_PROTOTYPE(TYPE,ROWS_UNROLL,COLS_UNROLL,UNROLLING,MOP,TRANSPOSITION)
dnl	----------------------------------------------------------------------------
dnl	Expands to the function prototype of the kernel specified by the macro arguments.
dnl
define(`RSB_M4_KERNEL_FUNCTION_PROTOTYPE',`dnl
void RSB_M4_KERNEL_FUNCTION_NAME($1,$2,$3,$4,$5,$6,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE) RSB_M4_KERNEL_FUNCTION_ARGS($1,$4,$5)'dnl
)dnl
dnl
dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl				Function body macros
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_NEGATE_FUNCTION_BODY_UNROLLED()
dnl	--------------------------------------
dnl	Negates the whole matrix;
dnl	Unrolled.
dnl
define(`RSB_M4_NEGATE_FUNCTION_BODY_UNROLLED',`dnl
dnl
	/* NOTE : should better use some intrinsic here. */
RSB_M4_DEBUGINFO(``$0'')dnl
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/* major_increment times unrolling on the major dimension, for each unroll on the minor dimension */
')dnl
	forloop(`major_unrolling',0,decr(major_increment),`forloop(`minor_unrolling',0,decr(minor_increment),
	`a[(minor_unrolling*major_increment)+major_unrolling]=-a[(minor_unrolling*major_increment)+major_unrolling];
	')'
	)
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl	RSB_M4_UNROLL_L_NEGATE_FUNCTION_BODY()
dnl	--------------------------------------
dnl	Negates the whole matrix; 
dnl
define(`RSB_M4_UNROLL_L_NEGATE_FUNCTION_BODY',`dnl
dnl
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/*!
	 * loop fully unrolled minor_increment times on the minor dimension, major_increment times on the major one
	 */
')dnl
RSB_M4_DEBUGINFO(``$0'')dnl
	register itype minor_index,major_index;
	for(minor_index=0;minor_index+eval(minor_increment-1)<minor_maximum;minor_index+=minor_increment)
	{
		for(major_index=0;major_index+eval(major_increment-1)<major_maximum;major_index+=major_increment)
		{
			forloop(`minor_unrolling',0,decr(minor_increment),
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
			`/* major_increment times unrolling on the major dimension,  */
')dnl
			forloop(`major_unrolling',0,decr(major_increment),
			`a[major_maximum*(minor_index+minor_unrolling)+major_index+major_unrolling]=-a[major_maximum*(minor_index+minor_unrolling)+major_index+major_unrolling];
			')
			')
		}
		ifelse(major_increment,1,`',
		`/* we handle the last (columns mod major_increment) columns */
		for(;major_index<major_maximum;++major_index)
		{forloop(`minor_unrolling',0,decr(minor_increment),`
		a[major_maximum*(minor_index+minor_unrolling)+major_index]=-a[major_maximum*(minor_index+minor_unrolling)+major_index];
		')
		}
		')
	}
	ifelse(minor_increment,1,`',
	`/* we handle the last (rows mod minor_increment) rows entirely */
	for(;minor_index<minor_maximum;++minor_index) for(major_index=0;major_index<major_maximum;++major_index)a[major_maximum*(minor_index)+major_index]=-a[major_maximum*(minor_index)+major_index];
	')
')dnl
dnl
dnl
dnl	RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED()
dnl	------------------------------------------
dnl	Expands to the fully unrolled infinity norm kernel.
dnl
define(`RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED',`dnl
dnl
ifelse(RSB_M4_IS_COMPLEX_TYPE(mtype),1,`/* FIXME : THE FOLLOWING CODE IS NOT CORRECT */')dnl

`#define CABS(X)' RSB_M4_ABS(mtype,X)
ifelse(transposition,RSB_M4_TRANS_T,`dnl
pushdef(`transposed_major_increment',minor_increment)`'dnl
pushdef(`transposed_minor_increment',major_increment)`'dnl
pushdef(`transposed_minor_unrolling',major_unrolling)`'dnl
',`dnl
pushdef(`transposed_minor_increment',minor_increment)`'dnl
pushdef(`transposed_major_increment',major_increment)`'dnl
pushdef(`transposed_minor_unrolling',minor_unrolling)`'dnl
')dnl

	/* NOTE : should better use some intrinsic here. */
RSB_M4_DEBUGINFO(``$0'')dnl
	forloop(`minor_unrolling',0,decr(transposed_minor_increment),`register mtype `sum_'minor_unrolling=0;
	')
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/* major_increment times unrolling on the major dimension, for each unroll on the minor dimension */
')dnl
	forloop(`major_unrolling',0,decr(major_increment),`forloop(`minor_unrolling',0,decr(minor_increment),
	``sum_'transposed_minor_unrolling+=CABS(RSB_M4_CONJ(a[(minor_unrolling*major_increment)+major_unrolling],mtype,transposition))<0?-RSB_M4_CONJ(a[(minor_unrolling*major_increment)+major_unrolling],mtype,transposition,k_symmetry):RSB_M4_CONJ(a[(minor_unrolling*major_increment)+major_unrolling],mtype,transposition,k_symmetry);
	')'
	)
	forloop(`minor_unrolling',0,decr(transposed_minor_increment),
	`local_row_sums[minor_unrolling]+=`sum_'minor_unrolling;
	')
`#undef CABS'
dnl
popdef(`transposed_major_increment')`'dnl
popdef(`transposed_minor_increment')`'dnl

')dnl
dnl
dnl	end RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl
dnl	RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED()
dnl	------------------------------------------
dnl	Expands to the fully unrolled infinity norm kernel.
dnl
dnl	FIXME: should document this is used also for rowssums
dnl
define(`RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED',`dnl
dnl
ifelse(RSB_M4_IS_COMPLEX_TYPE(mtype),1,`/* FIXME : THE FOLLOWING CODE IS NOT CORRECT */')dnl

ifelse(mop,`rowssums',`dnl
pushdef(`abs_or_not',`0')`'dnl
')dnl
dnl
ifelse(mop,`infty_norm',`dnl
pushdef(`abs_or_not',`1')`'dnl
dnl`#define CABS(X)' RSB_M4_ABS(mtype,X)
dnl pushdef(`RSB_M4_CABS_',`CABS')`'dnl
dnl pushdef(`RSB_M4_CABS_',`RSB_M4_ABS')`'dnl
')dnl
dnl
ifelse(transposition,RSB_M4_TRANS_T,`dnl
pushdef(`transposed_major_increment',minor_increment)`'dnl
pushdef(`transposed_minor_increment',major_increment)`'dnl
pushdef(`transposed_minor_unrolling',major_unrolling)`'dnl
',`dnl
pushdef(`transposed_minor_increment',minor_increment)`'dnl
pushdef(`transposed_major_increment',major_increment)`'dnl
pushdef(`transposed_minor_unrolling',minor_unrolling)`'dnl
')dnl
dnl
ifelse(transposition,RSB_M4_TRANS_N,`dnl
pushdef(`transposed_row_sums_off',roff)`'dnl
pushdef(`retransposed_row_sums_off',coff)`'dnl
',`dnl
pushdef(`transposed_row_sums_off',coff)`'dnl
pushdef(`retransposed_row_sums_off',roff)`'dnl
')dnl

	/* NOTE : should better use some intrinsic here. */
RSB_M4_DEBUGINFO(``$0'')dnl
	forloop(`minor_unrolling',0,decr(transposed_minor_increment),`register mtype `sum_'minor_unrolling=0;
	')
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/* major_increment times unrolling on the major dimension, for each unroll on the minor dimension */
')dnl
dnl
	forloop(`major_unrolling',0,decr(major_increment),`forloop(`minor_unrolling',0,decr(minor_increment),
	``sum_'transposed_minor_unrolling += RSB_M4_ABS_IF_1(mtype,RSB_M4_CONJ(a[(minor_unrolling*major_increment)+major_unrolling],mtype,transposition),abs_or_not);
	')'
	)
	forloop(`minor_unrolling',0,decr(transposed_minor_increment),
	`local_row_sums[transposed_row_sums_off+minor_unrolling]+=`sum_'minor_unrolling;
	')
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),0,`dnl
dnl	unsymmetric case: ok
',`dnl
dnl	symmetric case: 
	if(roff!=coff || i!=j)
	forloop(`minor_unrolling',0,decr(transposed_minor_increment),
	`	row_sums[retransposed_row_sums_off+minor_unrolling+bci]+=`sum_'minor_unrolling;
	')
')dnl
ifelse(mop,`infty_norm',`dnl
dnl `#undef CABS'
')dnl
dnl
popdef(`retransposed_row_sums_off')`'dnl
popdef(`transposed_row_sums_off')`'dnl
popdef(`abs_or_not')`'dnl
popdef(`transposed_major_increment')`'dnl
popdef(`transposed_minor_increment')`'dnl

')dnl
dnl
dnl	end RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl	RSB_M4_UNROLL_L_INFTY_NORM_FUNCTION_BODY()
dnl	------------------------------------------
dnl	Expands to the unrolled (with loops) infinity norm kernel.
dnl
define(`RSB_M4_UNROLL_L_INFTY_NORM_FUNCTION_BODY',`dnl
dnl
ifelse(RSB_M4_IS_COMPLEX_TYPE(mtype),1,`/* FIXME : THE FOLLOWING CODE IS NOT CORRECT */')dnl
dnl
ifelse(RSB_M4_IS_BCSR(minor_increment,major_increment),`1',`dnl
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/*!
	 * loop fully unrolled minor_increment times on the minor dimension, major_increment times on the major one
	 */
')dnl
')dnl
dnl
RSB_M4_DEBUGINFO(``$0'')dnl

`#define CABS(X)' RSB_M4_ABS(mtype,X)

	register itype minor_index,major_index;
	for(minor_index=0;minor_index+eval(minor_increment-1)<minor_maximum;minor_index+=minor_increment)
	{
		for(major_index=0;major_index+eval(major_increment-1)<major_maximum;major_index+=major_increment)
		{
			forloop(`minor_unrolling',0,decr(minor_increment),
			`/* major_increment times unrolling on the major dimension,  */
			forloop(`major_unrolling',0,decr(major_increment),
			`local_row_sums[minor_index+minor_unrolling]+=CABS(RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index+major_unrolling],mtype,transposition))<0? - RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index+major_unrolling],mtype,transposition,k_symmetry):RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index+major_unrolling],mtype,transposition,k_symmetry);
			')
			')
		}
		ifelse(major_increment,1,`',
		`/* we handle the last (columns mod major_increment) columns */
		for(;major_index<major_maximum;++major_index)
		{forloop(`minor_unrolling',0,decr(minor_increment),`
		local_row_sums[minor_index+minor_unrolling]+=RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index],mtype,transposition,k_symmetry);
		')
		}
		')
	}
`#undef CABS'
	ifelse(minor_increment,1,`',
	`/* we handle the last (rows mod minor_increment) rows entirely */
	for(;minor_index<minor_maximum;++minor_index) for(major_index=0;major_index<major_maximum;++major_index)local_row_sums[minor_index]+=RSB_M4_CONJ(a[major_maximum*(minor_index)+major_index],mtype,transposition,k_symmetry);
	')
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl	RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED()
dnl	-----------------------------------------
dnl	Expands to the fully unrolled row scaling kernel.
dnl
define(`RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED',`dnl
dnl
dnl	* Apply a scaling to the rows of the matrix, or apply a scalar scaling to all the coefficients.
dnl	* 
dnl	* Equivalent to multiply against a k x m matrix (the d vector is sized m) where column i values 
dnl	* are all d[i].
dnl
dnl	TODO : C ORDER ONLY
dnl
RSB_M4_DEBUGINFO(``$0'')dnl
ifelse(transposition,RSB_M4_TRANS_T,`
dnl	Fortran order
	forloop(`major_unrolling',0,decr(major_increment),`forloop(`minor_unrolling',0,decr(minor_increment),
	`a[(minor_unrolling*major_increment)+major_unrolling]*=d[major_unrolling];
	')'
	)
',`
dnl	C order
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/* major_increment times unrolling on the major dimension, for each unroll on the minor dimension */
')dnl
	forloop(`major_unrolling',0,decr(major_increment),`forloop(`minor_unrolling',0,decr(minor_increment),
	`a[(minor_unrolling*major_increment)+major_unrolling]*=d[minor_unrolling];
	')'
	)
')dnl
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl	RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED_L()
dnl	-------------------------------------------
dnl	Expands to the unrolled (with loops) row scaling kernel.
dnl
define(`RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED_L',`dnl
dnl
RSB_M4_DEBUGINFO(``$0'')dnl
	/*
	* Apply a scaling to the rows of the matrix, or apply a scalar scaling to all the coefficients.
	* 
	* Equivalent to multiply against a k x m matrix (the d vector is sized m) where column i values 
	* are all d[i].
	*/

		register itype minor_index,major_index;
		for(minor_index=0;minor_index+eval(minor_increment-1)<minor_maximum;minor_index+=minor_increment)
		{
			for(major_index=0;major_index+eval(major_increment-1)<major_maximum;major_index+=major_increment)
			{
				forloop(`minor_unrolling',0,decr(minor_increment),
				`/* major_increment times unrolling on the major dimension,  */
				forloop(`major_unrolling',0,decr(major_increment),
				`a[(minor_index+minor_unrolling)*columns+major_index+major_unrolling]*=d[minor_index+minor_unrolling];
				')
				')
			}
			ifelse(major_increment,1,`',
			`/* we handle the last (columns mod major_increment) columns */
			for(;major_index<major_maximum;++major_index)
			{forloop(`minor_unrolling',0,decr(minor_increment),`
				a[(minor_index+minor_unrolling)*columns+major_index]*=d[minor_index+minor_unrolling];')
			}
			')
		}
	ifelse(minor_increment,1,`',
		`/* we handle the last (rows mod minor_increment) rows entirely */
		for(;minor_index<minor_maximum;++minor_index) for(major_index=0;major_index<major_maximum;++major_index)
		a[minor_index * major_maximum + major_index ]*=d[minor_index];
		')
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl	RSB_M4_MV_FUNCTION_BODY_UNROLLED()
dnl	----------------------------------
dnl	Expands to the fully unrolled matrix vector multiplication kernel.
dnl
dnl	TODO : This macro is unused when register blocking is on. Should unify both.
dnl
define(`RSB_M4_MV_FUNCTION_BODY_UNROLLED',`dnl
dnl
RSB_M4_DEBUGINFO(``$0'')dnl
	forloop(`minor_unrolling',0,decr(minor_increment),`register mtype `c_'minor_unrolling=0;
	')
dnl	register mtype c_=0;
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/* major_increment times unrolling on the major dimension, for each unroll on the minor dimension */
')dnl
	forloop(`major_unrolling',0,decr(major_increment),`forloop(`minor_unrolling',0,decr(minor_increment),
	``c_'minor_unrolling += RSB_M4_CONJ(a[(minor_unrolling*major_increment)+major_unrolling],mtype,transposition,k_symmetry)*b[major_unrolling];
	')'
	)
	forloop(`minor_unrolling',0,decr(minor_increment),
	`c[minor_unrolling]+=`c_'minor_unrolling;
	/*c[minor_unrolling]+= alpha * `c_'minor_unrolling + beta * c[minor_unrolling];*/
	')
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl	RSB_M4_MV_FUNCTION_BODY_UNROLLED_REGISTER_BLOCKED()
dnl	---------------------------------------------------
dnl	Expands to the fully unrolled, register blocked matrix vector multiplication kernel.
dnl
define(`RSB_M4_MV_FUNCTION_BODY_UNROLLED_REGISTER_BLOCKED',`dnl
dnl
dnl	FIXME : missing documentation
dnl
dnl	20101112 spsv.* subkernel is supported here	
dnl
dnl	/* experimental RSB_M4_REGISTERS registers tuned blocking technique. */
dnl
pushdef(`should_merge_value_after_inner_loop',`RSB_M4_should_merge_value_after_inner_loop_inner')dnl
dnl
forloop(`RB',0,eval(minor_increment/RSB_M4_REGISTERS),`dnl
ifelse(eval(RB*RSB_M4_REGISTERS<minor_increment),1,`dnl
ifelse(eval(RB*RSB_M4_REGISTERS<minor_increment),1,`dnl
	{
dnl		/* using alpha and beta seems like a 10% hit for some 8x8 setups */
dnl		/*mtype *alphap=1, *betap=0;*/
dnl		/* beginning of a register block */',`
')dnl
dnl
dnl	FIXME : new
dnl
dnl
dnl ifelse(transposition,RSB_M4_TRANS_T,`
dnl /* this is a transposed kernel */
dnl ',`')dnl
dnl
pushdef(`c_unrolling_rb',`ifelse(transposition,RSB_M4_TRANS_T,`major_unrolling',`minor_unrolling_rb')')dnl
pushdef(`b_unrolling_rb',`ifelse(transposition,RSB_M4_TRANS_T,`minor_unrolling_rb',`major_unrolling')')dnl
pushdef(`c_increment',`ifelse(transposition,RSB_M4_TRANS_T,`major_increment',`minor_increment')')dnl
pushdef(`b_increment',`ifelse(transposition,RSB_M4_TRANS_T,`minor_increment',`major_increment')')dnl
dnl pushdef(`c_dest',ifelse(should_merge_value_after_inner_loop,`1',`cacc',`c[c_unrolling_rb]'))dnl
pushdef(`c_dest',ifelse(should_merge_value_after_inner_loop,`1',`cacc',`c[c_unrolling_rb]'))dnl
dnl
forloop(`c_unrolling_rb',eval((RB*RSB_M4_REGISTERS)),ifelse(eval((RB+1)*RSB_M4_REGISTERS>=c_increment),1,decr(c_increment),decr(eval((RB+1)*RSB_M4_REGISTERS))),`
dnl	register block variables declaration
ifelse(should_merge_value_after_inner_loop,`1',`',`dnl
		register mtype `c_'c_unrolling_rb = RSB_M4_ZERO(mtype);
')dnl
dnl		format(`register %s `c_'%02d=0;',mtype,c_unrolling_rb)dnl
		')',`')dnl
dnl	register block variables declaration
		ifelse(eval(RB*RSB_M4_REGISTERS<minor_increment),1,`
dnl	register block variables declaration
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
		/* major_increment times unrolling on the major dimension, for each unroll on the minor dimension */
')dnl
dnl	FIXME : THE FOLLOWING LOOP UNROLLINGS ORDER HAS INFLUENCE ON PERFORMANCE AND SHOULD BE STUDIED AT LOW LEVEL !
forloop(`minor_unrolling_rb',eval((RB*RSB_M4_REGISTERS)),ifelse(eval((RB+1)*RSB_M4_REGISTERS>=minor_increment),1,decr(minor_increment),decr(eval((RB+1)*RSB_M4_REGISTERS))),dnl
`
forloop(`major_unrolling',0,decr(major_increment),`dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
ifelse(should_merge_value_after_inner_loop,`1',`dnl
`		cacc'+=RSB_M4_CONJ(a[(minor_unrolling_rb*major_increment)+major_unrolling],mtype,transposition,k_symmetry)*bn[b_unrolling_rb];
',`dnl
`		c_'c_unrolling_rb += RSB_M4_CONJ(a[(minor_unrolling_rb*major_increment)+major_unrolling],mtype,transposition,k_symmetry)*bn[b_unrolling_rb];
')dnl
',`dnl
'ifelse(should_merge_value_after_inner_loop,`1',`dnl
`		cacc'+=RSB_M4_CONJ(a[(minor_unrolling_rb*major_increment)+major_unrolling],mtype,transposition,k_symmetry)*b[b_unrolling_rb];
',`dnl
`		c_'c_unrolling_rb += RSB_M4_CONJ(a[(minor_unrolling_rb*major_increment)+major_unrolling],mtype,transposition,k_symmetry)*b[b_unrolling_rb];
')dnl
)dnl
')dnl
')dnl
dnl
forloop(`c_unrolling_rb',eval((RB*RSB_M4_REGISTERS)),ifelse(eval((RB+1)*RSB_M4_REGISTERS>=c_increment),1,decr(c_increment),decr(eval((RB+1)*RSB_M4_REGISTERS))),`dnl
ifelse(should_merge_value_after_inner_loop,`1',`',`dnl
ifelse(RSB_M4_IS_SPMX_OP_NEGATING_KERNEL_MOP(mop),`1',`dnl
			c[c_unrolling_rb]-=`c_'c_unrolling_rb;
')dnl
ifelse(RSB_M4_IS_SPMX_OP_ADDING_KERNEL_MOP(mop),`1',`dnl
			c[c_unrolling_rb]+=`c_'c_unrolling_rb;
')dnl
ifelse(RSB_M4_IS_SPMX_OP_SCALING_KERNEL_MOP(mop),`1',`dnl
ifelse(RSB_M4_is_transposed_spmv,1,`dnl
			c[c_unrolling_rb]+= `c_'c_unrolling_rb ;
',`dnl
			c[c_unrolling_rb]+= alpha * `c_'c_unrolling_rb ;
')dnl
')dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),`1',`dnl
dnl
dnl	FIXME: an optimized implementation places this out.
dnl
			c[c_unrolling_rb]+= `c_'c_unrolling_rb ;
')dnl
')dnl
')dnl
dnl
	}',`')dnl
	')
dnl
popdef(`should_merge_value_after_inner_loop')dnl
popdef(`c_dest')dnl
popdef(`b_unrolling_rb')dnl
popdef(`c_unrolling_rb')dnl
popdef(`c_increment')dnl
popdef(`b_increment')dnl
dnl
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl	RSB_M4_UNROLL_L_MV_FUNCTION_BODY()
dnl	----------------------------------
dnl	Expands to the unrolled (but with loops), matrix vector multiplication kernel.
dnl
define(`RSB_M4_UNROLL_L_MV_FUNCTION_BODY',`dnl
dnl
RSB_M4_DEBUGINFO(``$0'')dnl
	/*!
	 * Loop fully unrolled minor_increment times on the minor dimension, major_increment times on the major one.
	 * Assumes a,b,c, are mtype arrays : a is minor_maximum x major_maximum, b is major_maximum x 1, c is minor_maximum x 1.
	 * Assumes matrices are in column major (C) order.
	 */
	register itype minor_index,major_index;
dnl
pushdef(`spmv_uxux_alpha',`ifelse(RSB_M4_IS_SPMX_OP_SCALING_KERNEL_MOP(mop),`1',`alpha*',`')')dnl
dnl
	for(minor_index=0;minor_index+eval(minor_increment-1)<minor_maximum;minor_index+=minor_increment)
	{
		for(major_index=0;major_index+eval(major_increment-1)<major_maximum;major_index+=major_increment)
		{
			forloop(`minor_unrolling',0,decr(minor_increment),
			`/* major_increment times unrolling on the major dimension,  */
			forloop(`major_unrolling',0,decr(major_increment),
			`c[minor_index+minor_unrolling]+=spmv_uxux_alpha`'RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index+major_unrolling],mtype,transposition,k_symmetry)*b[major_index+major_unrolling];
			')
			')
		}
		ifelse(major_increment,1,`',
		`/* we handle the last (columns mod major_increment) columns */
		for(;major_index<major_maximum;++major_index)
		{forloop(`minor_unrolling',0,decr(minor_increment),
ifelse(RSB_M4_IS_SPMX_OP_ADDING_KERNEL_MOP(mop),1,`dnl
		c[minor_index+minor_unrolling]+=RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index],mtype,transposition,k_symmetry)*b[major_index];
		')dnl
ifelse(RSB_M4_IS_SPMX_OP_NEGATING_KERNEL_MOP(mop),`1',`dnl
		c[minor_index+minor_unrolling]-=RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index],mtype,transposition,k_symmetry)*b[major_index];
		')dnl
ifelse(RSB_M4_IS_SPMX_OP_SCALING_KERNEL_MOP(mop),`1',`dnl
		c[minor_index+minor_unrolling]+=alpha*RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index],mtype,transposition,k_symmetry)*b[major_index];
	')dnl
		)
		}
		')
	}
	ifelse(minor_increment,1,`',
	`/* we handle the last (rows mod minor_increment) rows entirely */
	for(;minor_index<minor_maximum;++minor_index) for(major_index=0;major_index<major_maximum;++major_index)c[minor_index]+=spmv_uxux_alpha`'RSB_M4_CONJ(a[major_maximum*(minor_index)+major_index],mtype,transposition,k_symmetry)*b[major_index];
	')
dnl
popdef(`spmv_uxux_alpha')dnl
dnl
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl	RSB_M4_UNROLL_MM_FUNCTION_BODY()
dnl	--------------------------------
dnl	Expands to the matrix matrix multiplication kernel.
dnl
define(`RSB_M4_UNROLL_MM_FUNCTION_BODY',`dnl
RSB_M4_DEBUGINFO(``$0'')dnl
pushdef(`r',`rh')
	/*!
	 * `RSB_M4_UNROLL_MM_FUNCTION_BODY kernel'
	 * Loop fully unrolled minor_increment times on the minor dimension, major_increment times on the major one.
	 * Assumes a,b,c, are mtype arrays : a is minor_increment x major_increment, b is major_increment x nrhs, c is minor_increment x nrhs.
	 * Assumes matrices are in column major (C) order.
	 */
	/* WARNING : THIS IS JUST A SKETCH OF THE CODE : IT IS NOT SUPPOSED TO BE CORRECT NOR FAST */
	register itype r;
	for(r=0;r<nrhs;++r)/* right hand side columns could have a big stride, though */
	{
		/* loop fully unrolled minor_increment times on the minor dimension, major_increment times on the major one */
		/* this function will touch a block eval(minor_increment*major_increment)*sizeof(mtype) matrix bytes */
forloop(`minor_unrolling',0,decr(minor_increment),`dnl
		/* major_increment times unrolling on the major dimension,  */
forloop(`major_unrolling',0,decr(major_increment),`dnl
		c[r*cstride+minor_unrolling]+=RSB_M4_CONJ(a[minor_unrolling*major_increment+major_unrolling],mtype,transposition,k_symmetry)*b[r*bstride+major_unrolling];
')dnl
')dnl
	}
popdef(`r')
')dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl	RSB_M4_UNROLL_L_MM_FUNCTION_BODY()
dnl	----------------------------------
dnl	Expands to the matrix matrix multiplication kernel.
dnl
define(`RSB_M4_UNROLL_L_MM_FUNCTION_BODY',`dnl
RSB_M4_DEBUGINFO(``$0'')dnl
pushdef(`r',`rh')
	/*!
	 * Loop unrolled minor_increment times on the minor dimension, major_increment times on the major one.
	 * Assumes a,b,c, are mtype arrays : a is minor_maximum x major_maximum, b is MM x nrhs, c is mM x nrhs.
	 * Assumes matrices are in column major (C) order.
	 */
	/* WARNING : THIS IS JUST A SKETCH OF THE CODE : IT IS NOT SUPPOSED TO BE CORRECT NOR FAST */
	register itype minor_index=0,major_index=0,r=0;
	for(r=0;r<nrhs;++r)/* right hand side columns could have a big stride, though */
	{
		for(minor_index=0;minor_index+eval(minor_increment-1)<minor_maximum;minor_index+=minor_increment)
		{
			for(major_index=0;major_index+eval(major_increment-1)<major_maximum;major_index+=major_increment)
			{
forloop(`minor_unrolling',0,decr(minor_increment),dnl
				`/* major_increment times unrolling on the major dimension,  */
forloop(`major_unrolling',0,decr(major_increment),dnl
				`c[r*cstride+minor_index+minor_unrolling]+=RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index+major_unrolling],mtype,transposition,k_symmetry)*b[r*bstride+major_index+major_unrolling];
')dnl
')dnl
			}
			ifelse(major_increment,1,`',
			`/* we handle the last (columns mod major_increment) columns */
			for(;major_index<major_maximum;++major_index)
			{forloop(`minor_unrolling',0,decr(minor_increment),`
				c[r*cstride+minor_index+minor_unrolling]+=RSB_M4_CONJ(a[major_maximum*(minor_index+minor_unrolling)+major_index],mtype,transposition,k_symmetry)*b[r*bstride+major_index];')
			}
			')
		}
		
dnl	in the following, we omit the generation of a whole double loop if minor_increment is 1
		ifelse(minor_increment,1,`',
		`/* we handle the last (rows mod minor_increment) rows entirely */
		for(;minor_index<minor_maximum;++minor_index)
			for(major_index=0;major_index<major_maximum;++major_index)
				c[r*cstride+minor_index]+=RSB_M4_CONJ(a[major_maximum*(minor_index)+major_index],mtype,transposition,k_symmetry)*b[r*bstride+major_index];
			')
	}
popdef(`r')
')dnl
dnl	---------------------------------------------------------------------------	dnl
dnl				Function body dispatcher
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl	RSB_M4_KERNEL_FUNCTION_BODY(	r,minor_maximum,minor_increment,major_index_basename,major_maximum,major_increment,
dnl					mtype,want_header,mop,unrolling)
dnl	-------------------------------------------------------------------------------------------------------------------
dnl	Expands to the function body of a particular computational kernel.
dnl
define(`RSB_M4_KERNEL_FUNCTION_BODY',`dnl
dnl
pushdef(`r',$1)dnl
pushdef(`minor_index',`$1_$3')dnl minor dimension unroll variable identifier
pushdef(`minor_increment',$3)dnl minor dimension increment
pushdef(`minor_maximum',$2)dnl minor dimension maximum
pushdef(`major_index',$4_$6)dnl major dimension unroll variable identifier
pushdef(`major_increment',$6)dnl minor dimension maximum
pushdef(`major_maximum',$5)dnl major dimension maximum
pushdef(`mtype',$7)dnl
pushdef(`itype',`rsb_coo_idx_t')dnl
pushdef(`want_header',$8)dnl
pushdef(`mop',$9)dnl
pushdef(`unrolling',`ifelse($10,`l',`l',)')dnl
pushdef(`k_symmetry',$11)dnl
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/*!
	 * Loop fully unrolled minor_increment times on the minor dimension, major_increment times on the major one.
	 * Assumes a,b,c, are mtype arrays : a is minor_increment x major_increment, b is major_increment x 1, c is minor_increment x 1.
	 * Assumes matrices are in column major (C) order.
	 */
')dnl
ifelse(RSB_M4_IS_SPSX_KERNEL_MOP(mop),`1',`dnl
dnl
dnl	FIXME : THIS IS A TEMPORARY HACK
dnl
RSB_M4_MV_FUNCTION_BODY_UNROLLED_REGISTER_BLOCKED(mop)dnl
dnl RSB_M4_KERNEL_FUNCTION_BODY($1,$2,$3,$4,$5,$6,$7,$8,`spmv_uauz',$10,$11)dnl
')dnl
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
ifelse(unrolling,`l',`dnl
dnl	Matrix vector unroll with loop
RSB_M4_UNROLL_L_MV_FUNCTION_BODY(mop)
',`dnl unrolling else
dnl	Matrix vector unroll without any loop
dnl
	ifelse(eval(RSB_M4_WANT_BLOCKING>=1),`1',`dnl
RSB_M4_MV_FUNCTION_BODY_UNROLLED_REGISTER_BLOCKED(mop)',`
RSB_M4_MV_FUNCTION_BODY_UNROLLED')
')dnl end unrolling ifelse
')dnl 
dnl
dnl same as spmv_uauz:
dnl
ifelse(mop,`spmm_az',`dnl
dnl`spmm_az',`dnl mop else
ifelse(unrolling,`l',`dnl
dnl	Matrix - matrix unroll with loop
RSB_M4_UNROLL_L_MM_FUNCTION_BODY
',`dnl unrolling else
dnl	Matrix - matrix unroll without loop
/* WRONG */
RSB_M4_UNROLL_MM_FUNCTION_BODY
/* WRONG */
')
dnl,`unknown `unrolling' : urolling ?')dnl end unrolling ifelse
popdef(`spmm_az')
dnl ',`unknown `mop' : mop ?'
')dnl end mop ifelse
ifelse(mop,`scale',`dnl
ifelse(unrolling,`l',`dnl
RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED_L
',`dnl unrolling else
RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED
')
')dnl end mop ifelse
dnl
ifelse(mop,`rowssums',`dnl
ifelse(unrolling,`l',`dnl
RSB_M4_ERROR_UNIMPLEMENTED
',`dnl unrolling else
RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED(mop)
')dnl
')dnl end mop ifelse
dnl
ifelse(mop,`infty_norm',`dnl
ifelse(unrolling,`l',`dnl
RSB_M4_UNROLL_L_INFTY_NORM_FUNCTION_BODY(mop)
',`dnl unrolling else
RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED(mop)
')
')dnl end mop ifelse
dnl
ifelse(mop,`negation',`dnl
ifelse(unrolling,`l',`dnl
RSB_M4_UNROLL_L_NEGATE_FUNCTION_BODY(mop)
',`dnl unrolling else
RSB_M4_NEGATE_FUNCTION_BODY_UNROLLED(mop)
')
')dnl end mop ifelse
dnl
popdef(`k_symmetry')dnl
popdef(`unrolling')dnl
popdef(`mop')dnl
popdef(`want_header')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`major_maximum')dnl
popdef(`major_increment')dnl
popdef(`major_index')dnl
popdef(`minor_maximum')dnl
popdef(`minor_increment')dnl
popdef(`spmm_az')dnl
popdef(`r')dnl
')dnl end RSB_M4_KERNEL_FUNCTION_BODY
dnl
dnl
dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl				Function definitions
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
define(`RSB_M4_EXTRA_SYMMETRIC_DIAGONAL_FIXING_KERNEL',`dnl
dnl
dnl	UNFINISHED
dnl
pushdef(`r',$1)dnl
pushdef(`minor_index',`$1_$3')dnl minor dimension unroll variable identifier
pushdef(`minor_increment',$3)dnl minor dimension increment
pushdef(`minor_maximum',$2)dnl minor dimension maximum
pushdef(`major_index',$4_$6)dnl major dimension unroll variable identifier
pushdef(`major_increment',$6)dnl minor dimension maximum
pushdef(`major_maximum',$5)dnl major dimension maximum
pushdef(`mtype',$7)dnl
pushdef(`itype',`rsb_coo_idx_t')dnl
pushdef(`want_header',$8)dnl
pushdef(`mop',$9)`'dnl
pushdef(`unrolling',`ifelse($10,`l',`l',`')')dnl
pushdef(`transposition',$11)`'dnl
pushdef(`k_symmetry',$12)`'dnl
RSB_M4_DEBUGINFO(``$0'')
ifelse(RSB_M4_want_verbose_comments,`1',`dnl
dnl /* : UNFINISHED : FIXME */
dnl
/*
	Should determine the offset (in terms of elements) to first diagonal element.
	Should determine the intersection length.
	For each diagonal intersecting element, subtract the outcome of operation mop.
	`transposition' : "transposition", `symmetry' : "k_symmetry"
*/
')dnl
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),1,`pushdef(`transposition',RSB_M4_TRANS_N)')dnl
dnl
ifelse(RSB_M4_IS_SPMX_OP_SCALING_KERNEL_MOP(mop),`1',`dnl
			*c -= RSB_M4_CONJ(*a,mtype,transposition,k_symmetry)**b*alpha; /* no matrix pointer advance is needed, as this is a corrective term */
')dnl
ifelse(RSB_M4_IS_SPMX_OP_NEGATING_KERNEL_MOP(mop),`1',`dnl
			*c += RSB_M4_CONJ(*a,mtype,transposition,k_symmetry)**b; /* no matrix pointer advance is needed, as this is a corrective term */
')dnl
ifelse(RSB_M4_IS_SPMX_OP_ADDING_KERNEL_MOP(mop),`1',`dnl
			*c -= RSB_M4_CONJ(*a,mtype,transposition,k_symmetry)**b; /* no matrix pointer advance is needed, as this is a corrective term */
')dnl
dnl
dnl	FIXME : symmetry is not among the macro arguments, and this is DANGEROUS
dnl
ifelse(RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),1,`popdef(`transposition')')dnl
dnl
dnl	UNFINISHED
dnl
dnl
dnl	UNFINISHED
dnl
popdef(`k_symmetry')`'dnl
popdef(`unrolling')dnl
popdef(`mop')dnl
popdef(`want_header')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`major_maximum')dnl
popdef(`major_increment')dnl
popdef(`major_index')dnl
popdef(`minor_maximum')dnl
popdef(`minor_increment')dnl
popdef(`spmm_az')dnl
popdef(`r')dnl
popdef(`transposition')dnl
')dnl end RSB_M4_EXTRA_SYMMETRIC_DIAGONAL_FIXING_KERNEL
dnl	
dnl
dnl
dnl
dnl
dnl
dnl	---------------------------------------------------------------------------	dnl
dnl				Function definitions
dnl	---------------------------------------------------------------------------	dnl
dnl
dnl
dnl	RSB_M4_UNROLL_KERNEL(	r,minor_maximum,minor_increment,major_index_basename,major_maximum,major_increment,
dnl				mtype,want_header,mop,unrolling)
dnl	-----------------------------------------------------------------------------------------------------------
dnl	A general macro for matrix-matrix and matrix-vector multiplication unrolled kernels.
dnl	FIXME : THIS MACRO (WILL BE) DEPRECATED.
dnl
define(`RSB_M4_UNROLL_KERNEL',`dnl
dnl
dnl
pushdef(`r',$1)dnl
pushdef(`minor_index',`$1_$3')dnl minor dimension unroll variable identifier
pushdef(`minor_increment',$3)dnl minor dimension increment
pushdef(`minor_maximum',$2)dnl minor dimension maximum
pushdef(`major_index',$4_$6)dnl major dimension unroll variable identifier
pushdef(`major_increment',$6)dnl minor dimension maximum
pushdef(`major_maximum',$5)dnl major dimension maximum
pushdef(`mtype',$7)dnl
pushdef(`itype',`rsb_coo_idx_t')dnl
pushdef(`want_header',$8)dnl
pushdef(`mop',$9)`'dnl
pushdef(`unrolling',`ifelse($10,`l',`l',`')')dnl
pushdef(`transposition',$11)`'dnl
RSB_M4_DEBUGINFO(``$0'')dnl
dnl
ifelse(want_header,`h',`dnl
dnl
dnl Only the function header is expanded.
dnl
RSB_M4_KERNEL_FUNCTION_PROTOTYPE(mtype,minor_increment,major_increment,unrolling,mop,transposition);
',`dnl
dnl
dnl The entire function declaration is expanded.
dnl
RSB_M4_KERNEL_FUNCTION_PROTOTYPE(mtype,minor_increment,major_increment,unrolling,mop,transposition)
{
ifelse(RSB_M4_IS_SPMV_KERNEL_MOP(mop),1,`dnl
dnl
dnl	Matrix Vector product
dnl
ifelse(unrolling,`l',`dnl
dnl
dnl	Matrix Vector unroll with loop
dnl
RSB_M4_UNROLL_L_MV_FUNCTION_BODY
',`dnl unrolling else
dnl
dnl	Matrix vector unroll without any loop
dnl
ifelse(RSB_M4_there_is_real_blocking,`1',`dnl
	/*!
	 * Loop fully unrolled minor_increment times on the minor dimension, major_increment times on the major one.
	 * Assumes a,b,c, are mtype arrays : a is minor_increment x major_increment, b is major_increment x 1, c is minor_increment x 1.
	 * Assumes matrices are in column major (C) order.
	 */
')dnl
	ifelse(eval(RSB_M4_WANT_BLOCKING>=1),`1',`dnl
RSB_M4_MV_FUNCTION_BODY_UNROLLED_REGISTER_BLOCKED',`
RSB_M4_MV_FUNCTION_BODY_UNROLLED')
')dnl end unrolling ifelse
')
dnl
ifelse(mop,`spmm_az',`dnl
dnl
dnl	Matrix Matrix product
dnl
pushdef(`spmm_az',`$1_$3')dnl minor dimension unroll variable identifier
ifelse(unrolling,`l',`dnl
dnl
dnl	Matrix - matrix unroll with loop
dnl
RSB_M4_UNROLL_L_MM_FUNCTION_BODY
',`dnl unrolling else
dnl
dnl	Matrix - matrix unroll without loop
dnl
RSB_M4_UNROLL_MM_FUNCTION_BODY
')
popdef(`spmm_az')
')dnl
dnl
dnl
ifelse(mop,`scale',`dnl
dnl
dnl	Matrix Scaling
dnl
ifelse(unrolling,`l',`dnl
/* WARNING : THIS IS NOT LOOPED ! */
RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED_L
',`dnl unrolling else
RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED
')
')dnl
dnl
ifelse(mop,`spmv_uxux',`dnl
dnl
dnl	Matrix Scaling
dnl
ifelse(unrolling,`l',`dnl
/* WARNING : THIS IS NOT LOOPED ! */
dnl RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED_L(mop)
RSB_M4_UNROLL_L_MV_FUNCTION_BODY(mop)
',`dnl unrolling else
dnl RSB_M4_ROW_SCALE_FUNCTION_BODY_UNROLLED(mop)
RSB_M4_MV_FUNCTION_BODY_UNROLLED_REGISTER_BLOCKED(mop)
')
')dnl
dnl
ifelse(mop,`rowssums',`dnl
ifelse(unrolling,`l',`dnl
RSB_M4_ERROR_UNIMPLEMENTED
',`dnl unrolling else
RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED(mop)
')dnl
')dnl end mop ifelse
dnldnl
ifelse(mop,`infty_norm',`dnl
dnl
dnl	Infinity Norm
dnl
ifelse(unrolling,`l',`dnl
RSB_M4_UNROLL_L_INFTY_NORM_FUNCTION_BODY(mop)
',`dnl unrolling else
RSB_M4_INFTY_NORM_FUNCTION_BODY_UNROLLED(mop)
')
')dnl end mop ifelse
dnl
ifelse(mop,`negation',`dnl
ifelse(unrolling,`l',`dnl
RSB_M4_UNROLL_L_NEGATE_FUNCTION_BODY(mop)
',`dnl unrolling else
RSB_M4_NEGATE_FUNCTION_BODY_UNROLLED(mop)
')
')dnl end mop ifelse
dnl
dnl
	return;
}
')dnl end header ifelse
popdef(`unrolling')dnl
popdef(`mop')dnl
popdef(`want_header')dnl
popdef(`mtype')dnl
popdef(`itype')dnl
popdef(`major_maximum')dnl
popdef(`major_increment')dnl
popdef(`major_index')dnl
popdef(`minor_maximum')dnl
popdef(`minor_increment')dnl
popdef(`spmm_az')dnl
popdef(`r')dnl
popdef(`transposition')dnl
')dnl end RSB_M4_UNROLL_KERNEL
dnl	
dnl	
dnl	DOUBLE_LINEAR_KERNEL_SEARCH(MOP,TYPE,...)
dnl	-----------------------------------------
dnl	
dnl	The following M4 macro generates C macro code for selecting the right completely unrolled function.
dnl	It is not optimal, but it exists.
dnl	
define(`DOUBLE_LINEAR_KERNEL_SEARCH',`dnl
pushdef(`mop',$1)dnl
pushdef(`type',$2)dnl
pushdef(`transposition',$3)dnl
ifelse($#,4,`$4\'
,`dnl
( (R)==($4) && (C)==($5)?  pushdef(`rowsu',$4)pushdef(`colsu',$5)dnl
RSB_M4_KERNEL_FUNCTION_NAME(type,rowsu,colsu,looped,mop,transposition,RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE) \
popdef(`rowsu')popdef(`colsu') :  (DOUBLE_LINEAR_KERNEL_SEARCH(mop,type,transposition,shift(shift(shift(shift(shift($@)))))) )) dnl
')dnl
popdef(`transposition')dnl
popdef(`type')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl	DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_IDENTIFIER(MOP,TYPE,UNROLLING)
dnl	----------------------------------------------------------------
dnl
define(`DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_IDENTIFIER',`dnl
pushdef(`mop',$1)dnl
pushdef(`type',$2)dnl
pushdef(`unrolling',$3)dnl
`RSB_'RSB_M4_TYPE_CODE(type)`_kernel_'mop`'dnl
popdef(`unrolling')dnl
popdef(`type')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl	DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_(mop,type,unrolling,...)
dnl	----------------------------------------------------------
dnl	FIXME : someday there will be a binary search macro
dnl
define(`DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_',`dnl
pushdef(`mop',$1)dnl
pushdef(`type',$2)dnl
pushdef(`unrolling',$3)dnl
pushdef(`unrolls',shift(shift(shift($@))))dnl
`#define' DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_IDENTIFIER(mop,type,unrolling)`'(R,C) DOUBLE_LINEAR_KERNEL_SEARCH(mop,type,unrolls)
popdef(`unrolls')dnl
popdef(`unrolling')dnl
popdef(`mop')dnl
popdef(`type')dnl
')dnl
dnl
dnl
dnl	UNLOOP_R_C_PAIRS()
dnl	------------------
dnl
dnl	generates a list of (R,C) pairs; the cartesian product of rowsu and colsu lists
dnl
define(`UNLOOP_R_C_PAIRS',`foreach(`rowsu',RSB_M4_ROWS_UNROLL,`foreach(`colsu',RSB_M4_COLUMNS_UNROLL,`rowsu,colsu,')')')dnl
dnl
dnl 
dnl	RSB_M4_TYPES : ...
dnl 
define(`RSB_M4_TYPES',(WANT_TYPES))dnl
dnl 
dnl 
dnl 
dnl	RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)
dnl	-------------------------------------------------------
dnl	FIXME
dnl
dnl #define RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)
dnl
define(`RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING',`dnl
pushdef(`mtype',$1)`'dnl
RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(mtype)`'`_PRINTF_STRING'dnl
popdef(`mtype')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl 
dnl	RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG(mtype)
dnl	-------------------------------------------------------
dnl	FIXME
dnl
dnl #define RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG(mtype)
dnl
define(`RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG',`dnl
pushdef(`mtype',$1)`'dnl
pushdef(`arg',$2)`'dnl
ifelse(mtype,`long double complex',`creall(arg),cimagl(arg)',`dnl
ifelse(mtype,`double complex',`creal(arg),cimag(arg)',`dnl
ifelse(mtype,`float complex',`crealf(arg),cimagf(arg)',`dnl
ifelse(mtype,`complex',`creal(arg),cimag(arg),creal(arg),cimag(arg)',`dnl
arg`'dnl
')dnl
')dnl
')dnl
')dnl
popdef(`arg')`'dnl
popdef(`mtype')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl 
dnl	RSB_M4_MATRIX_STORAGE_PREPROCESSOR_STRING(matrix_storage)
dnl	---------------------------------------------------------
dnl	FIXME
dnl
dnl #define RSB_M4_MATRIX_STORAGE_PREPROCESSOR_STRING(matrix_storage)
dnl
define(`RSB_M4_MATRIX_STORAGE_PREPROCESSOR_STRING',`dnl
pushdef(`matrix_storage',$1)`'dnl
RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(matrix_storage)`'`_STRING'dnl
popdef(`matrix_storage')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl	RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL()
dnl	-----------------------------------------------
dnl
dnl #define RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL() 
dnl
define(`RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL',`dnl
pushdef(`citype',$1)`'dnl
`RSB_COORDINATE_TYPE_'RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE(citype)`'dnl
popdef(`citype')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl	RSB_M4_MATRIX_DIAGONAL_PREPROCESSOR_SYMBOL()
dnl	-----------------------------------------------
dnl
dnl #define RSB_M4_MATRIX_DIAGONAL_PREPROCESSOR_SYMBOL() 
dnl
define(`RSB_M4_MATRIX_DIAGONAL_PREPROCESSOR_SYMBOL',`dnl
pushdef(`k_diagonal',$1)`'dnl
`RSB_DIAGONAL_'RSB_M4_MATRIX_DIAGONAL_CHAR(k_diagonal)`'dnl
popdef(`k_diagonal')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl	RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL()
dnl	-----------------------------------------------
dnl
dnl #define RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL() 
dnl
define(`RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL',`dnl
pushdef(`transposition',$1)`'dnl
`RSB_TRANSPOSITION_'RSB_M4_MATRIX_TRANSPOSITION_CHAR(transposition)`'dnl
popdef(`transposition')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl	RSB_M4_MATRIX_SYMMETRY_PREPROCESSOR_SYMBOL(type)
dnl	-----------------------------------------------
dnl
dnl #define RSB_M4_MATRIX_SYMMETRY_PREPROCESSOR_SYMBOL(type) 
dnl
define(`RSB_M4_MATRIX_SYMMETRY_PREPROCESSOR_SYMBOL',`dnl
pushdef(`k_symmetry',$1)`'dnl
`RSB_SYMMETRY_'touppercase(RSB_M4_CHOPSPACES(k_symmetry))`'dnl
popdef(`k_symmetry')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl	RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(type)
dnl	-----------------------------------------------
dnl	FIXME
dnl
dnl #define RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(type) 
dnl
define(`RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL',`dnl
pushdef(`matrix_storage',$1)`'dnl
`RSB_MATRIX_STORAGE_'touppercase(RSB_M4_CHOPSPACES(matrix_storage))`'dnl
popdef(`matrix_storage')`'dnl
')dnl
dnl 
dnl 
dnl 
dnl	RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)
dnl	-----------------------------------------------
dnl	Converts a type name in a preprocessor symbol used to indicate the type availability.
dnl
dnl #define RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type) 
dnl
define(`RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL',`dnl
pushdef(`type',$1)`'dnl
`RSB_NUMERICAL_TYPE_'touppercase( RSB_M4_CHOPSPACES(type) )dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_MATRIX_DIAGONAL_CHAR()
dnl	------------------------------------
dnl	FIXME
dnl
define(`RSB_M4_MATRIX_DIAGONAL_CHAR',`dnl
pushdef(`k_diagonal',$1)`'dnl
`'touppercase(RSB_M4_CHOPSPACES(k_diagonal))`'dnl
popdef(`k_diagonal')`'dnl
')dnl
dnl
dnl
dnl	RSB_M4_MATRIX_TRANSPOSITION_CHAR()
dnl	------------------------------------
dnl	FIXME
dnl
define(`RSB_M4_MATRIX_TRANSPOSITION_CHAR',`dnl
pushdef(`transposition',$1)`'dnl
`'touppercase(RSB_M4_CHOPSPACES(transposition))`'dnl
popdef(`transposition')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_DEFAULT_COORDINATE_INDEX_TYPE',`rsb_coo_idx_t')dnl
dnl
dnl	RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE()
dnl	------------------------------------
dnl	FIXME !
dnl
define(`RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE_',`dnl
pushdef(`citype',$1)`'dnl
ifelse(citype,`rsb_coo_idx_t',`0x01')`'dnl
ifelse(citype,`rsb_half_idx_t',`0x02')`'dnl
popdef(`citype')`'dnl
')dnl
dnl
dnl
dnl	RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE()
dnl	------------------------------------
dnl	FIXME
dnl
define(`RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE',`dnl
pushdef(`citype',$1)`'dnl
ifelse(citype,`rsb_coo_idx_t',`C')`'dnl
ifelse(citype,`rsb_half_idx_t',`H')`'dnl
popdef(`citype')`'dnl
')dnl
dnl
dnl
dnl	RSB_M4_MATRIX_TRANSPOSITION_CHARCODE()
dnl	------------------------------------
dnl	FIXME
dnl
define(`RSB_M4_MATRIX_TRANSPOSITION_CHARCODE',`dnl
pushdef(`transposition',$1)`'dnl
dnl ifelse(transposition,RSB_M4_TRANS_T,`0x01 /*!< Transposed flag value, valid for \ref rsb_trans_t valued variables. */')`'dnl
dnl ifelse(transposition,RSB_M4_TRANS_N,`0x00 /*!< Non transposed flag, valid for \ref rsb_trans_t typed variables. */')`'dnl
dnl ifelse(transposition,RSB_M4_TRANS_C,`0x02 /*!< Conjugated transpose flag, valid for \ref rsb_trans_t typed variables. */')`'dnl
ifelse(transposition,RSB_M4_TRANS_T,`0x54 /*!< \brief T: Transposed flag value, valid for \ref rsb_trans_t valued variables. */')`'dnl
ifelse(transposition,RSB_M4_TRANS_N,`0x4E /*!< \brief N: Non transposed flag, valid for \ref rsb_trans_t typed variables. */')`'dnl
ifelse(transposition,RSB_M4_TRANS_C,`0x43 /*!< \brief C: Conjugated transpose flag, valid for \ref rsb_trans_t typed variables. */')`'dnl
ifelse(transposition,RSB_M4_TRANS_INVALID,`0x3F /*!< \brief ?: Transposition type flag value guaranteed to be invalid. Useful for tests. Valid as char. */')`'dnl
popdef(`transposition')`'dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_MATRIX_DIAGONAL_CHARCODE(k_diagonal)
dnl	------------------------------------
dnl	FIXME
dnl
define(`RSB_M4_MATRIX_DIAGONAL_CHARCODE',`dnl
pushdef(`k_diagonal',$1)`'dnl
ifelse(k_diagonal,`e',`0x01 /* Explicit diagonal (default, implicit) */')`'dnl
ifelse(k_diagonal,`i',`0x02 /* Implicit diagonal */')`'dnl FIXME : new
popdef(`k_diagonal')`'dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_MATRIX_SYMMETRY_CHARCODE(k_symmetry)
dnl	------------------------------------
dnl	FIXME
dnl
define(`RSB_M4_MATRIX_SYMMETRY_CHARCODE',`dnl
pushdef(`k_symmetry',$1)`'dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_SYMMETRIC,`RSB_FLAG_SYMMETRIC /*  */')`'dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_HERMITIAN,`RSB_FLAG_HERMITIAN /*  */')`'dnl FIXME : new
ifelse(k_symmetry,RSB_M4_SYMBOL_UNSYMMETRIC,`0x00 /*  */')`'dnl
popdef(`k_symmetry')`'dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_MATRIX_STORAGE_CHARCODE(type)
dnl	------------------------------------
dnl	FIXME
dnl
define(`RSB_M4_MATRIX_STORAGE_CHARCODE',`dnl
pushdef(`matrix_storage',$1)`'dnl
ifelse(matrix_storage,`BCSR',`0x01 /*  */')`'dnl
ifelse(matrix_storage,`BCSC',`0x02 /* */ ')`'dnl
ifelse(matrix_storage,`VBR',`0x04 /* */')`'dnl
ifelse(matrix_storage,`VBC',`0x08 /* */')`'dnl
ifelse(matrix_storage,`LR',`0x10 /* */')`'dnl
ifelse(matrix_storage,`LC',`0x20 /* */')`'dnl
ifelse(matrix_storage,`BCOR',`0x40 /* */')`'dnl
ifelse(matrix_storage,`BCOC',`0x80 /* */')`'dnl
dnl ifelse(matrix_storage,`LVBR',`0x10 /* */')`'dnl
dnl ifelse(matrix_storage,`LVBC',`0x20 /* */')`'dnl
dnl ifelse(matrix_storage,`BCOO',`0x10 /* */')`'dnl
popdef(`matrix_storage')`'dnl
')dnl
dnl 
dnl
dnl
dnl
dnl	RSB_M4_TYPE_CHARCODE_ASCII_VALUE(type)
dnl	--------------------------------------
dnl
define(`RSB_M4_TYPE_CHARCODE_ASCII_VALUE',`dnl
pushdef(`type',$1)`'dnl
ifelse(type,`long double',`76')`'dnl /* L (? compliance)*/
ifelse(type,`double',`68')`'dnl /* D (BLAS compliance)*/
ifelse(type,`float',`83')`'dnl /*F (BLAS compliance)*/ 
ifelse(type,`int',`73')`'dnl /*I*/
ifelse(type,`unsigned int',`85')`'dnl /*U*/
ifelse(type,`char',`h')`'dnl /*h*/
ifelse(type,`unsigned char',`72')`'dnl /*H*/
dnl ifelse(type,`complex',`99')`'dnl /*c*/
ifelse(type,`float complex',`67')`'dnl /*c (BLAS compliance)*/
ifelse(type,`double complex',`90')`'dnl /*z (BLAS compliance)*/
ifelse(type,`long double complex',`81')`'dnl /*Q (? compliance)*/
ifelse(type,RSB_M4_INVALID_TYPE,`?')`'dnl /* invalid type */
popdef(`type')`'dnl
')dnl
dnl 
dnl	
dnl	dnl
dnl
dnl
dnl	RSB_M4_TYPE_CHARCODE(type)
dnl	--------------------------
dnl	Converts a type name to a char type value (in hex, for quoting reasons).
dnl
define(`RSB_M4_TYPE_CHARCODE',`dnl
pushdef(`type',$1)`'dnl
ifelse(type,`long double',`L')`'dnl /* L (? compliance)*/
ifelse(type,`double',`D')`'dnl /* D (BLAS compliance)*/
ifelse(type,`float',`S')`'dnl /* F (BLAS compliance)*/ 
ifelse(type,`int',`I')`'dnl /* I */
ifelse(type,`unsigned int',`U')`'dnl /*U*/
ifelse(type,`char',`h')`'dnl /*h*/
ifelse(type,`unsigned char',`H')`'dnl /*H*/
ifelse(type,`float complex',`C')`'dnl /* c (BLAS compliance)*/
ifelse(type,`double complex',`Z')`'dnl /* z (BLAS compliance)*/
ifelse(type,`long double complex',`Q')`'dnl /* Q (? compliance)*/
ifelse(type,RSB_M4_INVALID_TYPE,`?')`'dnl /* invalid type */
popdef(`type')`'dnl
')dnl
dnl 
dnl	
dnl	
dnl	DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH(...)
dnl	-------------------------------------------
dnl	FIXME: REDOCUMENT
dnl	
dnl	The following code generates macro for selecting the right completely unrolled function.
dnl	It is not optimal, but it exists.
dnl	
define(`DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH',`dnl
pushdef(`mop',$1)dnl
pushdef(`type',$2)dnl
pushdef(`transposition',$3)dnl
ifelse($#,4,`$4\'
,`dnl
( (R)==($4) && (C)==($5)?  pushdef(`rowsu',$4)pushdef(`colsu',$5)dnl
RSB_M4_KERNEL_DIRECT_DISPATCHER_FUNCTION_NAME(type,`BCSR',transposition,k_symmetry,rowsu,colsu,unrolling,mop,,) \
popdef(`rowsu')popdef(`colsu') :  (DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH(mop,type,transposition,shift(shift(shift(shift(shift($@)))))) )) dnl
')dnl
popdef(`type')dnl
popdef(`mop')dnl
popdef(`transposition')dnl
')dnl
dnl
dnl
dnl	DOUBLE_LINEAR_KERNEL_SEARCH_MACRO_IDENTIFIER(...)
dnl	-------------------------------------------------
dnl	FIXME: REDOCUMENT
dnl
define(`DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH_MACRO_IDENTIFIER',`dnl
pushdef(`mop',$1)dnl
pushdef(`type',$2)dnl
pushdef(`unrolling',$3)dnl
pushdef(`matrix_storage',$4)dnl
pushdef(`transposition',$5)dnl
`RSB_'RSB_M4_TYPE_CODE(type)`_kernel_dispatcher_'matrix_storage`_'mop`_'unrolling`_'RSB_M4_TRANSPOSITION_CODE(transposition)`'dnl
popdef(`matrix_storage')dnl
popdef(`transposition')dnl
popdef(`unrolling')dnl
popdef(`type')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl	DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH_MACRO_(mop,type,unrolling,...)
dnl	---------------------------------------------------------------------
dnl	FIXME: REDOCUMENT
dnl
define(`DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH_MACRO_',`dnl
pushdef(`mop',$1)dnl
pushdef(`type',$2)dnl
pushdef(`unrolling',$3)dnl
pushdef(`matrix_storage',$4)dnl
pushdef(`transposition',$5)dnl
pushdef(`k_symmetry',$5)dnl
pushdef(`unrolls',shift(shift(shift(shift(shift(shift($@)))))))dnl
/* a macro is faster than a switch construct */
`#define' DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH_MACRO_IDENTIFIER(mop,type,unrolling,matrix_storage,transposition)`'(R,C) DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH(mop,type,transposition,unrolls)
popdef(`k_symmetry')dnl
popdef(`matrix_storage')dnl
popdef(`transposition')dnl
popdef(`unrolls')dnl
popdef(`unrolling')dnl
popdef(`mop')dnl
popdef(`type')dnl
')dnl
dnl
dnl
dnl	KERNEL_TYPE_DISPATCHER_SEARCH_MACRO_IDENTIFIER(MOP,TYPE,UNROLLING,..)
dnl	------------------------------------------------------------------
dnl	FIXME: REDOCUMENT
dnl
define(`KERNEL_TYPE_DISPATCHER_SEARCH_MACRO_IDENTIFIER',`dnl
pushdef(`mop',$1)dnl
pushdef(`type',$2)dnl
pushdef(`unrolling',$3)dnl
pushdef(`matrix_storage',$4)dnl
pushdef(`transposition',$5)dnl
`RSB_type_kernel_dispatcher_'matrix_storage`_'mop`_'unrolling`_'RSB_M4_TRANSPOSITION_CODE(transposition)`'dnl
popdef(`transposition')dnl
popdef(`matrix_storage')dnl
popdef(`unrolling')dnl
popdef(`type')dnl
popdef(`mop')dnl
')dnl
dnl
dnl
dnl	DOUBLE_LINEAR_KERNEL_DISPATCHER_TYPE_SEARCH_MACRO_(MOP,UNROLLING,..)
dnl	-----------------------------------------------------------------
dnl	FIXME: REDOCUMENT
dnl
define(`DOUBLE_LINEAR_KERNEL_DISPATCHER_TYPE_SEARCH_MACRO_',`dnl
pushdef(`mop',$1)dnl
pushdef(`unrolling',$2)dnl
pushdef(`matrix_storage',$3)dnl
pushdef(`transposition',$3)dnl
/* a macro is faster than a switch construct */
`#define' KERNEL_TYPE_DISPATCHER_SEARCH_MACRO_IDENTIFIER(mop,type,unrolling,matrix_storage,transposition)`'(TYPE,R,C) \
(dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
  (TYPE)==RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type) ? (void*)DOUBLE_LINEAR_KERNEL_DISPATCHER_SEARCH_MACRO_IDENTIFIER(mop,type,unrolling,matrix_storage,transposition)(R,C) : \
dnl
')dnl
NULL ) dnl
dnl
popdef(`transposition')dnl
popdef(`matrix_storage')dnl
popdef(`unrolling')dnl
popdef(`type')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME(mtype,matrix_storage,transposition,k_symmetry,unrolling,mop)
dnl	-----------------------------------------------------------------------------
dnl
define(`RSB_M4_KERNEL_SIZE_DISPATCH_FUNCTION_NAME',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`matrix_storage',$2)dnl	
pushdef(`transposition',$3)dnl	
pushdef(`k_symmetry',$4)dnl	
pushdef(`unrolling',$5)dnl	
pushdef(`mop',$6)dnl	
pushdef(`citype',$7)dnl	
pushdef(`k_diagonal',$8)dnl	
pushdef(`uplo',$9)dnl
dnl
dnl	FIXME : should handle citype
dnl
RSB_M4_PREFIX`'matrix_storage`_'mop`_'RSB_M4_TYPE_CODE(mtype)`_'RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE(citype)`_'unrolling`'dnl
dnl
dnl ifelse(RSB_M4_MEMBER(matrix_storage,`VBR',`VBC'),1,`dnl
dnl RSB_M4_PREFIX`'matrix_storage`_'mop`_'RSB_M4_TYPE_CODE(mtype)`_'RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE(citype)`_'unrolling`'dnl
dnl ')dnl
dnl ifelse(RSB_M4_MEMBER(matrix_storage,`BCSR',`BCSC'),1,`dnl
dnl RSB_M4_PREFIX`'matrix_storage`_'mop`_'RSB_M4_TYPE_CODE(mtype)`_'RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE(citype)`_'unrolling`'dnl
dnl ')dnl
dnl ifelse(RSB_M4_MEMBER(matrix_storage,`LR',`LC'),1,`dnl
dnl dnl NEW
dnl RSB_M4_PREFIX`'matrix_storage`_'mop`_'RSB_M4_TYPE_CODE(mtype)`_'RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE(citype)`_'unrolling`'dnl
dnl ')dnl
dnl ifelse(transposition,RSB_M4_TRANS_T,`_T',`_N')`'dnl
`_t'RSB_M4_MATRIX_TRANSPOSITION_CHAR(transposition)`'dnl
`_s'touppercase(k_symmetry)`'dnl
`_d'touppercase(k_diagonal)`'dnl
`_u'touppercase(uplo)`'dnl
dnl else should give error : fixme
popdef(`uplo')dnl
popdef(`k_diagonal')dnl
popdef(`citype')dnl
popdef(`mop')dnl
popdef(`mtype')dnl
popdef(`matrix_storage')dnl
popdef(`k_symmetry')dnl
popdef(`transposition')dnl
popdef(`unrolling')dnl
')dnl
dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR',WANT_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR)dnl
dnl define(`RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR',`16')dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_MEDIUM',`8')dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR_SMALL',`4')dnl
dnl
dnl	RSB_M4_SIMPLE_LOOP_UNROLL_2S_WITH_JUMP()
dnl	---------------------------
dnl	A quick nd dirty way to unroll simple loops.
dnl	Could cause infinite recursion on identifiers clash.
dnl
dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL_2S_WITH_JUMP',`dnl
pushdef(`_ii',$1)dnl
pushdef(`_LI',`$2')dnl
pushdef(`_li',$3)dnl
pushdef(`_ui',$4)dnl
pushdef(`_st1',$5)dnl
pushdef(`_st2',$6)dnl
pushdef(`_uff',`eval(ifelse($7,,RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR,$7))')dnl
pushdef(`_uf',`eval(4*_uff)')dnl
pushdef(`_ff',`eval(1*_uff)')dnl
pushdef(`_hf',`eval(_uff/2)')dnl
{
	_ii=_li;
	if((_ui-_ii)<(_ff))goto anolu;
dnl	if((_ui-_ii)<(_uf+_ff+_hf))goto nolu;
dnl	switch(_ui%eval(_uf/4)){case }
dnl
for(;_ii+decr(_uf)<_ui;_ii+=_uf){`'pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',0,decr(_uf),`_st1`'')`'forloop(`_LI_',0,decr(_uf),`_st2`'')dnl
}`'popdef(_LI)dnl
dnl
for(;_ii+decr(_ff)<_ui;_ii+=_ff){`'pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',0,decr(_ff),`_st1`'')`'forloop(`_LI_',0,decr(_ff),`_st2`'')dnl
}`'popdef(_LI)dnl
dnl
anolu:
dnl
ifelse(eval(_hf),0,`',`dnl
for(;_ii+decr(_hf)<_ui;_ii+=_hf){`'pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',0,decr(_hf),`_st1`'')`'forloop(`_LI_',0,decr(_hf),`_st2`'')dnl
}`'popdef(_LI)dnl
')dnl
dnl
pushdef(`_LI_',`0')dnl
for(     ;_ii<_ui;++_ii){`'_st1`'_st2`'}
}`'popdef(`_LI_')dnl
dnl
popdef(`_ii')dnl
popdef(`_LI')dnl
popdef(`_li')dnl
popdef(`_ui')dnl
popdef(`_st2')dnl
popdef(`_st1')dnl
popdef(`_uf')dnl
popdef(`_ff')dnl
popdef(`_hf')dnl
popdef(`_uff')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_SIMPLE_LOOP_UNROLL_5S()
dnl	---------------------------
dnl	A quick nd dirty way to unroll simple loops.
dnl	Could cause infinite recursion on identifiers clash.
dnl
dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL_5S',`dnl
pushdef(`_ii',$1)dnl
pushdef(`_LI',`$2')dnl
pushdef(`_li',$3)dnl
pushdef(`_ui',$4)dnl
pushdef(`_st0',$5)dnl
pushdef(`_st1',$6)dnl
pushdef(`_st2',$8)dnl
pushdef(`_st3',$7)dnl
pushdef(`_st4',$9)dnl
pushdef(`_uf',`ifelse($10,,RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR,$10)')dnl
{
_st0`'dnl
for(_ii=_li;_ii+decr(_uf)<_ui;_ii+=_uf){
pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',0,decr(_uf),`dnl
_st1`'dnl
')dnl
popdef(_LI)dnl
_st3`'dnl
pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',0,decr(_uf),`dnl
_st2`'dnl
')dnl
popdef(_LI)dnl
_st4`'dnl
}
pushdef(`_LI_',`0')dnl
for(     ;_ii<_ui;++_ii){`'_st1`'_st2`'}
popdef(`_LI_')dnl
}
popdef(`_ii')dnl
popdef(`_LI')dnl
popdef(`_li')dnl
popdef(`_ui')dnl
popdef(`_st4')dnl
popdef(`_st3')dnl
popdef(`_st2')dnl
popdef(`_st1')dnl
popdef(`_st0')dnl
popdef(`_uf')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_SIMPLE_LOOP_UNROLL_2S()
dnl	---------------------------
dnl	A quick nd dirty way to unroll simple loops.
dnl	Could cause infinite recursion on identifiers clash.
dnl
dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL_2S',`dnl
RSB_M4_SIMPLE_LOOP_UNROLL_5S(`$1',`$2',`$3',`$4',`',`$5',`',`$6',`',`$7')`'dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_SIMPLE_LOOP_UNROLL()
dnl	---------------------------
dnl	A quick nd dirty way to unroll simple loops.
dnl	Could cause infinite recursion on identifiers clash.
dnl
dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL',`dnl
pushdef(`_ii',$1)dnl
pushdef(`_LI',`$2')dnl
pushdef(`_li',$3)dnl
pushdef(`_ui',$4)dnl
pushdef(`_st',$5)dnl
pushdef(`_uf',`ifelse($6,,RSB_M4_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR,$6)')dnl
ifelse(_uf,1,`dnl
for(_ii=_li;_ii+decr(_uf)<_ui;_ii+=_uf)
{
pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',0,decr(_uf),`dnl
_st`'dnl
')dnl
popdef(_LI)dnl
}
',`dnl
{
for(_ii=_li;_ii+decr(_uf)<_ui;_ii+=_uf){
pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',0,decr(_uf),`dnl
_st`'dnl
')dnl
popdef(_LI)dnl
}
pushdef(`_LI_',`0')dnl
for(     ;_ii<_ui;++_ii){ _st }
popdef(`_LI_')dnl
}')
popdef(`_ii')dnl
popdef(`_LI')dnl
popdef(`_li')dnl
popdef(`_ui')dnl
popdef(`_st')dnl
popdef(`_uf')dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_SIMPLE_LOOP_UNROLL_2S_J',`dnl
pushdef(`RSB_DO_WANT_PATCH_20101213',`0')dnl
ifelse(RSB_DO_WANT_PATCH_20101213,`1',`RSB_M4_SIMPLE_LOOP_UNROLL_2S_WITH_JUMP($@)',`RSB_M4_SIMPLE_LOOP_UNROLL_2S($@)')dnl
dnl
popdef(`RSB_DO_WANT_PATCH_20101213')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_SIMPLE_UNROLL()
dnl	----------------------
dnl	A quick and dirty way to unroll simple loops.
dnl	Could cause infinite recursion on identifiers clash.
dnl	It is buggy : doesn't support nested unrollings. (FIXME)
dnl
dnl
define(`RSB_M4_SIMPLE_UNROLL',`dnl
pushdef(`_LI',`$1')dnl
pushdef(`_li',$2)dnl
pushdef(`_st',$4)dnl
pushdef(`_uf',`ifelse($3,,16,$3)')dnl
{
pushdef(_LI,`_LI_ ')dnl
forloop(`_LI_',_li,decr(_uf),`dnl
_st`'dnl
')dnl
popdef(_LI)dnl
}
popdef(`_uf')dnl
popdef(`_st')dnl
popdef(`_li')dnl
popdef(`_LI')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_RESTRICT
dnl	---------------
dnl	The 'restrict' keyword of C99.
dnl	It forces strict aliasing in parameter arguments, leaving room for more compiler optimizations.
dnl
define(`RSB_M4_RESTRICT',`ifelse(RSB_M4_USE_RESTRICT,`1',`restrict',`')')`'dnl
dnl
dnl
dnl
