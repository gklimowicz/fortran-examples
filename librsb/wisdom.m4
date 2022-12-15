dnl
dnl
include(`rsb_misc.m4')dnl
dnl
dnl	Matrix modifying kernel operation codes.
define(`RSB_M4_MATRIX_WRITEONLY_KERNEL_MOPS',``scale',`negation'')dnl
dnl
dnl	SPMV-like kernel operation codes.
define(`RSB_M4_SPMV_KERNEL_MOPS',``spmv_sasa',`spmv_uauz',`spmv_uxux',`spmv_uaua',`spmv_uxua',`spmv_unua',`spmv_sxsx',`spmv_sxsa'')dnl
dnl
dnl	SPMM-like kernel operation codes.
define(`RSB_M4_SPMM_KERNEL_MOPS',``spmm_az'')dnl
dnl
dnl	SPMV/SPMM-like kernel operation codes.
define(`RSB_M4_SPMX_KERNEL_MOPS',`RSB_M4_SPMV_KERNEL_MOPS,RSB_M4_SPMM_KERNEL_MOPS')dnl
dnl
dnl	kernel operation codes scaling the output vector.
define(`RSB_M4_SPSX_ZEROING_KERNEL_MOPS',``spsv_uxua'')dnl
define(`RSB_M4_SPMX_ZEROING_KERNEL_MOPS',``spmv_uauz',`spmm_az'')dnl
dnl define(`RSB_M4_SPMX_ZEROING_KERNEL_MOPS',``spmv_uauz',`spmm_az',`spsv_uxua'')dnl
define(`RSB_M4_SPXX_ZEROING_KERNEL_MOPS',`RSB_M4_SPMX_ZEROING_KERNEL_MOPS')dnl
dnl define(`RSB_M4_SPXX_ZEROING_KERNEL_MOPS',`RSB_M4_SPSX_ZEROING_KERNEL_MOPS,RSB_M4_SPMX_ZEROING_KERNEL_MOPS')dnl
define(`RSB_M4_SPXX_ALLOW_ALIASING_KERNEL_MOPS',``spsv_uxua'')dnl
dnl
dnl	SP**-like kernel operation codes scaling the output vector.
define(`RSB_M4_SPMX_SCALING_KERNEL_MOPS',``spmv_uxux',`spmv_sxsx'')dnl
dnl
dnl	SP**-like kernel operation codes scaling,scaling, or adding only to the computation result
define(`RSB_M4_SPSX_OP_ADDING_KERNEL_MOPS',``spsv_uxua'')dnl
define(`RSB_M4_SPMX_OP_ADDING_KERNEL_MOPS',``spmv_uauz',`spmv_uaua',`spmv_sasa'')dnl
define(`RSB_M4_SPMX_OP_NEGATING_KERNEL_MOPS',``spmv_unua'')dnl
define(`RSB_M4_SPMX_OP_SCALING_KERNEL_MOPS',``spmv_uxua',`spmv_uxux',`spmv_sxsx',`spmv_sxsa'')dnl
define(`RSB_M4_SPSX_OP_SCALING_KERNEL_MOPS',``spsv_sxsx'')dnl
define(`RSB_M4_SPXX_OP_ACC_WRITING_KERNEL_MOPS',``infty_norm',`rowssums'')dnl
dnl
define(`RSB_M4_SPXX_RC_BIASED_KERNEL_MOPS',``scale'')dnl
dnl
dnl	Matrix readonly kernel operation codes.
define(`RSB_M4_MATRIX_READONLY_KERNEL_MOPS',`RSB_M4_SPXX_OP_ACC_WRITING_KERNEL_MOPS,RSB_M4_SPSX_KERNEL_MOPS,RSB_M4_SPMX_KERNEL_MOPS')dnl
dnl
dnl	SPSV-like kernel operation codes.
define(`RSB_M4_SPSV_KERNEL_MOPS',``spsv_uxua',`spsv_sxsx'')dnl
define(`RSB_M4_SPSX_KERNEL_MOPS',`RSB_M4_SPSV_KERNEL_MOPS')dnl
dnl
dnl	SPSV-like kernel operation codes scaling the output vector.
dnl define(`RSB_M4_SPSX_SCALING_KERNEL_MOPS',`spsv_sxsx')dnl
define(`RSB_M4_SPSX_SCALING_KERNEL_MOPS',`')dnl
define(`RSB_M4_SPSX_OP_SETTING_KERNEL_MOPS',`spsv_uxua')dnl
dnl
dnl	kernel operation codes scaling the output vector.
define(`RSB_M4_ZEROING_KERNEL_MOPS',`RSB_M4_SPXX_ZEROING_KERNEL_MOPS')dnl
define(`RSB_M4_SCALING_KERNEL_MOPS',`RSB_M4_SPMX_SCALING_KERNEL_MOPS,RSB_M4_SPSX_SCALING_KERNEL_MOPS')dnl
dnl
define(`RSB_M4_STRIDED_KERNEL_MOPS',``spmv_sxsx',`spmv_sxsa',`spsv_sxsx',`spmv_sasa'')dnl
dnl
define(`RSB_M4_OP_SCALING_KERNEL_MOPS',`RSB_M4_SPMX_OP_SCALING_KERNEL_MOPS,RSB_M4_SPSX_OP_SCALING_KERNEL_MOPS')dnl
define(`RSB_M4_OP_NEGATING_KERNEL_MOPS',`RSB_M4_SPMX_OP_NEGATING_KERNEL_MOPS')dnl
define(`RSB_M4_OP_ADDING_KERNEL_MOPS',`RSB_M4_SPMX_OP_ADDING_KERNEL_MOPS,RSB_M4_SPSX_OP_ADDING_KERNEL_MOPS')dnl
define(`RSB_M4_SPSX_NEGATING_KERNEL_MOPS',`')dnl
dnl
define(`RSB_M4_MATRIX_ALL_COMPLEX_TYPES',(`complex',`long complex',`float complex',`long double complex',`double complex'))dnl
define(`RSB_M4_MATRIX_INT_TYPES',(`int'))dnl
dnl
define(`RSB_M4_HAVE_COMPLEX_TYPE',`(RSB_M4_INTERSECTION(RSB_M4_MATRIX_ALL_COMPLEX_TYPES,(WANT_TYPES)))')dnl
define(`RSB_M4_HAVE_INT_TYPE',`(RSB_M4_INTERSECTION(RSB_M4_MATRIX_ALL_INT_TYPES,(WANT_TYPES)))')dnl
dnl
define(`RSB_M4_IS_READONLY_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_MATRIX_READONLY_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_WRITEONLY_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_MATRIX_WRITEONLY_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPMV_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMV_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPSV_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSV_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPMX_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMX_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPXX_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMX_KERNEL_MOPS,RSB_M4_SPSX_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPXV_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMV_KERNEL_MOPS,RSB_M4_SPSV_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPXM_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMM_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPMX_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMX_SCALING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPMX_OP_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMX_OP_SCALING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPMX_OP_NEGATING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMX_OP_NEGATING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPMX_OP_ADDING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMX_OP_ADDING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPSX_OP_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSX_OP_SCALING_KERNEL_MOPS)')dnl
dnl define(`RSB_M4_IS_SPSX_OP_NEGATING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSX_OP_NEGATING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPSX_OP_SETTING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSX_OP_SETTING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPSX_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSX_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPSX_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSX_SCALING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPXX_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSX_SCALING_KERNEL_MOPS,RSB_M4_SPMX_SCALING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPXX_OP_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPSX_OP_SCALING_KERNEL_MOPS,RSB_M4_SPMX_OP_SCALING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPXX_TWO_VECTORS_OPERATING_KERNEL_MOP',`RSB_M4_IS_SPXX_KERNEL_MOP($1)')dnl
define(`RSB_M4_IS_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SCALING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SCALING_OR_ZEROING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SCALING_KERNEL_MOPS,RSB_M4_ZEROING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_OP_SCALING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_OP_SCALING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_OP_NEGATING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_OP_NEGATING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_OP_ADDING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_OP_ADDING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_ZEROING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_ZEROING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_SPMX_ZEROING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPMX_ZEROING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_ACC_WRITING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPXX_OP_ACC_WRITING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_ALLOWING_ALIASING_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPXX_ALLOW_ALIASING_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_STRIDED_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_STRIDED_KERNEL_MOPS)')dnl
define(`RSB_M4_IS_RC_BIASED_KERNEL_MOP',`RSB_M4_MEMBER($1,RSB_M4_SPXX_RC_BIASED_KERNEL_MOPS)')dnl
define(`RSB_M4_MAXIMAL_CONFIGURED_BLOCK_SIZE',`RSB_M4_MAX2(RSB_M4_MAXN(WANT_COLUMN_UNLOOP_FACTORS),RSB_M4_MAXN(WANT_ROW_UNLOOP_FACTORS))')dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_there_is_real_blocking',`ifelse(`'minor_increment`'major_increment`',`11',`0',`1')')dnl
define(`RSB_M4_want_verbose_comments',`0')dnl
define(`RSB_M4_want_old_fortran_float_types',`0')dnl
dnl
define(`RSB_M4_there_is_no_real_blocking',`ifelse(`'minor_increment`'major_increment`',`11',`1',`0')')dnl
dnl
dnl
define(`RSB_M4_should_merge_value_after_inner_loop',`dnl
dnl
pushdef(`transposed',ifelse(tolowercase(transposition),RSB_M4_TRANS_N,0,1))dnl
RSB_M4_AND(dnl
RSB_M4_IS_SPMX_KERNEL_MOP(mop),dnl
RSB_M4_OR(dnl
RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry),dnl
RSB_M4_AND(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),transposed),dnl
RSB_M4_AND(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),RSB_M4_NOT(transposed))))dnl
popdef(`transposed')dnl
dnl
')dnl
dnl
define(`RSB_M4_should_merge_value_after_inner_loop_inner',`dnl
dnl
pushdef(`transposed',ifelse(tolowercase(transposition),RSB_M4_TRANS_N,0,1))dnl
RSB_M4_AND(dnl
RSB_M4_IS_SPMX_KERNEL_MOP(mop),dnl
RSB_M4_IS_UNSYMMETRIC(k_symmetry),dnl
RSB_M4_OR(dnl
RSB_M4_AND(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),transposed),dnl
RSB_M4_AND(RSB_M4_IS_FORMAT_ROW_MAJOR(matrix_storage),RSB_M4_NOT(transposed))))dnl
popdef(`transposed')dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_is_transposed_spmv',`dnl
dnl
pushdef(`transposed',ifelse(tolowercase(transposition),RSB_M4_TRANS_N,0,1))dnl
RSB_M4_AND(RSB_M4_OR(transposed,RSB_M4_IS_NOT_UNSYMMETRIC(k_symmetry)),RSB_M4_IS_SPMX_KERNEL_MOP(mop))dnl
popdef(`transposed')dnl
dnl
')dnl
dnl
dnl
