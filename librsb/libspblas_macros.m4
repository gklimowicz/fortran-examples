dnl
dnl
dnl	Sparse BLAS interface code generating macros.
dnl	Preliminary code.
dnl
define(`RSB_M4_SPBLAS_DOC_CONFIGUREDOUT_MSG',`dnl
	/*!
	  \ingroup rsb_doc_sparse_blas
	  \warning \rsb_warn_configuredout_msg
	*/
')dnl
dnl
define(`RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_MSG',`dnl
	/*!
	  \ingroup rsb_doc_sparse_blas
	  \warning \rsb_warn_unimplemented_msg
	*/
')dnl
dnl
dnl
define(`RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_SBL1_MSG',`dnl
	/*!
	  \ingroup rsb_doc_sparse_blas
	  RSB_M4_SPBLAS_HELP_INFO(mop)
	  \warning \rsb_spblasl1_msg
	*/
ifelse(lang,`lang_c',`dnl
	RSB_SPB_INTERFACE_PREAMBLE
')dnl
dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_DOC_COMMENT',`dnl
         /*!
           \ingroup rsb_doc_sparse_blas
           RSB_M4_SPBLAS_HELP_INFO(mop)
         */
ifelse(lang,`lang_c',`dnl
	RSB_SPB_INTERFACE_PREAMBLE
')dnl
dnl
')dnl
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
dnl
dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_LANGUAGES',(`lang_c',`f90'))dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_TYPES',(`float',`double',`float complex',`double complex'))dnl
define(`RSB_M4_SPBLAS_SYMMETRY_UL_CHARCODE',(`u'))dnl	FIXME
dnl
dnl
define(`RSB_SPBLAS_FUNCTION_IDENTIFIER',`dnl
pushdef(`type',$2)dnl
pushdef(`mop',$1)dnl
pushdef(`lang',$3)dnl
ifelse(lang,`f90',RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C,`')`'dnl
ifelse(lang,`f90',`dnl
`blas_'dnl
',`dnl
`BLAS_'dnl
')dnl
dnl
RSB_M4_SPBLAS_TYPE_CHARCODE(type)`'dnl
dnl
dnl	"US" stands for "Unstructured Sparse"
dnl
`us'dnl
mop`'dnl
dnl
dnl	FIXME: the trailing underscore should be removed.
dnl
ifelse(lang,`f90',RSB_M4_FORTRAN_SYMBOL_ADD_TO_C,`')`'dnl
dnl
popdef(`lang')dnl
popdef(`mop')dnl
popdef(`type')dnl
')dnl
dnl
dnl
dnl
define(`RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE',`dnl
pushdef(`type',$1)dnl
ifelse(RSB_M4_MEMBER(type,`double complex',`float complex'),`1',`',`dnl
ifelse(RSB_M4_MEMBER(type,`double complex',`float complex'),`1',`',`&')`'')dnl
popdef(`type')dnl
')dnl
dnl
dnl define(`RSB_SPBLAS_OVER_TYPE_POINTER_DEREFERENCE',`dnl
dnl pushdef(`type',$1)dnl
dnl ifelse(RSB_M4_MEMBER(type,`double complex',`float complex'),`1',`*',`dnl
dnl ifelse(RSB_M4_MEMBER(type,`double complex *',`float complex *'),`1',`*',` ')`'')dnl
dnl popdef(`type')dnl
dnl ')dnl
dnl
define(`RSB_M4_SPBLAS_FUNCTION',`dnl
pushdef(`type',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`tri',$3)dnl
pushdef(`want_what',$4)dnl
pushdef(`over',$5)dnl
pushdef(`lang',$6)dnl
dnl
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,want_what,over,lang)`'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,want_what,over,lang)`'dnl
RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(type,mop,tri,want_what,over,lang)`'dnl
RSB_M4_SPBLAS_EXTRA_FUNCTION(type,mop,tri,want_what,over,lang)`'dnl
dnl
popdef(`lang')dnl
popdef(`over')dnl
popdef(`want_what')dnl
popdef(`tri')dnl
popdef(`mop')dnl
popdef(`type')dnl
dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_SPBLAS_TO_RSB_FIX_ARGS',`dnl
pushdef(`type',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`tri',$3)dnl
pushdef(`want_what',$4)dnl
pushdef(`over',$5)dnl
pushdef(`lang',$6)dnl
dnl
dnl
ifelse(`'RSB_M4_MEMBER(mop,`cr_insert_entry')`'RSB_M4_IS_COMPLEX_TYPE(type)`',`10',`A,RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)val,i,j',`dnl
ifelse(`'RSB_M4_MEMBER(mop,`axpy')`'RSB_M4_MEMBER(lang,`lang_c')`',`11',`dnl
ifelse(`'RSB_M4_IS_COMPLEX_TYPE(type)`',`1',`/* FIXME: this is an exception; shall use a formal substitution technique, rather */nnz,alpha,x,indx,y,incy,index_base ',`dnl
nnz,&alpha,x,indx,y,incy,index_base')`'dnl
',`dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_SPBLAS_FUNCTION(type,mop,tri,`ARGS',over,lang)))`'dnl
')`'dnl
')`'dnl
dnl
dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`,RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)')`'dnl
dnl
popdef(`lang')dnl
popdef(`want_what')dnl
popdef(`over')dnl
popdef(`tri')dnl
popdef(`mop')dnl
popdef(`type')dnl
')dnl
dnl
define(`RSB_SPBLAS_OVER_TYPE',`dnl
pushdef(`type',$1)dnl
pushdef(`over',$2)dnl
ifelse(over,`1',`dnl
ifelse(RSB_M4_MEMBER(type,`double complex',`float complex'),`1',`const void *',`dnl
ifelse(RSB_M4_MEMBER(type,`double complex *',`float complex *'),`1',`void *',type )`'')dnl
',`dnl
type `'dnl
')dnl
dnl ifelse(RSB_M4_MEMBER(type,`double complex',`float complex'),`1',`const 'type` *',`dnl
dnl ifelse(RSB_M4_MEMBER(type,`double complex *',`float complex *'),`1',`'type`',type )`'')dnl
popdef(`over')dnl
popdef(`type')dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_TYPE_CHARCODE',`dnl
pushdef(`type',$1)`'dnl
dnl
tolowercase(`RSB_M4_TYPE_CHARCODE(type)')`'dnl
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_SPBLAS_MOP_CODE_TRANSLATE',`dnl
pushdef(`mop',$1)dnl
dnl
dnl	FIXME
dnl
ifelse(mop,`mv',`spmv_uxux')`'dnl
ifelse(mop,`mm',`spmm_xx')`'dnl
ifelse(mop,`sv',`spsv_uxua')`'dnl
ifelse(mop,`sm',`spsm_xxl')`'dnl
popdef(`mop')dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_SPBLAS_HELP_INFO',`dnl
pushdef(`mop',$1)dnl
dnl
dnl	FIXME: for some reason we are not using these!
dnl
ifelse(lang,`f90',`dnl
')dnl
ifelse(RSB_M4_MEMBER(mop,`dot'),`1',`\rsb_spblasl1_dot_msg')`'dnl
ifelse(RSB_M4_MEMBER(mop,`axpy'),`1',`\rsb_spblasl1_axpy_msg')`'dnl
ifelse(RSB_M4_MEMBER(mop,`ga'),`1',`\rsb_spblasl1_ga_msg')`'dnl
ifelse(RSB_M4_MEMBER(mop,`gz'),`1',`\rsb_spblasl1_gz_msg')`'dnl
ifelse(RSB_M4_MEMBER(mop,`sc'),`1',`\rsb_spblasl1_sc_msg')`'dnl
dnl
dnl
ifelse(mop,`mv',`\rsb_spblasl2_mv_msg')`'dnl
ifelse(mop,`mm',`\rsb_spblasl2_mm_msg')`'dnl
ifelse(mop,`sv',`\rsb_spblasl2_sv_msg')`'dnl
ifelse(mop,`sm',`\rsb_spblasl2_sm_msg')`'dnl
ifelse(mop,`cr_begin',`\rsb_spblasl2_cr_begin_msg')`'dnl
ifelse(mop,`cr_block_begin',`\rsb_spblasl2_cr_block_msg')`'dnl
ifelse(mop,`cr_variable_block_begin',`\rsb_spblasl2_cr_vbr_msg')`'dnl
ifelse(mop,`cr_insert_entry',`\rsb_spblasl2_cr_insert_entry_msg')`'dnl
ifelse(mop,`cr_insert_entries',`\rsb_spblasl2_cr_insert_entries_msg')`'dnl
ifelse(mop,`cr_insert_col',`\rsb_spblasl2_cr_insert_col_msg')`'dnl
ifelse(mop,`cr_insert_row',`\rsb_spblasl2_cr_insert_row_msg')`'dnl
ifelse(mop,`cr_insert_block',`\rsb_spblasl2_cr_insert_block_msg')`'dnl
ifelse(mop,`cr_insert_clique',`\rsb_spblasl2_cr_insert_clique_msg')`'dnl
ifelse(mop,`cr_end',`\rsb_spblasl2_cr_end_msg')`'dnl
ifelse(mop,`ds',`\rsb_spblasl2_ds_msg')`'dnl
dnl
ifelse(mop,`sp',`\rsb_spblasl2_sp_msg')`'dnl
ifelse(mop,`gp',`\rsb_spblasl2_gp_msg')`'dnl
dnl
dnl	Extra operations:
dnl
ifelse(mop,`get_diag',`\rsb_spblasl2e_usget_diag_msg')`'dnl
ifelse(mop,`get_element',`\rsb_spblasl2e_usget_element_norm_msg')`'dnl
ifelse(mop,`get_matrix_nnz',`\rsb_spblasl2e_usget_matrix_nnz_msg')`'dnl
ifelse(mop,`get_infinity_norm',`\rsb_spblasl2e_usget_infinity_norm_msg')`'dnl
ifelse(mop,`get_rows_nnz',`\rsb_spblasl2e_usget_rows_nnz_msg')`'dnl
ifelse(mop,`get_rows_sparse',`\rsb_spblasl2e_usget_rows_sparse_msg')`'dnl
ifelse(mop,`rows_scale',`\rsb_spblasl2e_usrows_scale_msg')`'dnl
ifelse(mop,`set_elements',`\rsb_spblasl2e_usset_elements_norm_msg.')`'dnl
ifelse(mop,`set_element',`\rsb_spblasl2e_usset_element_norm_msg')`'dnl
dnl
dnl
dnl
dnl
dnl
dnl
dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`dnl
ifelse(lang,`f90',`\rsb_spblasl2_Ap_msg \rsb_spblas_istat_msg \rsb_spblas_set_mtx_msg',`\rsb_spblas_return_mtx_msg')`'dnl
',`dnl
ifelse(lang,`f90',`\rsb_spblas_istat_msg',`\rsb_spblas_return_msg')`'dnl
')dnl
popdef(`mop')dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_LIST_MEMBER',`dnl
pushdef(`E',$1)dnl
pushdef(`L',$2)dnl
dnl
pushdef(`M',`F')dnl
dnl
foreach(`X',L,`dnl
ifelse(X,E,`ifelse(M,`T',`',`pushdef(`M',`T')')')dnl
')dnl
ifelse(M,`T',`1'popdef(`M'),`0')dnl
popdef(`M')dnl
dnl
popdef(`E')dnl
popdef(`L')dnl
dnl
')dnl
dnl
define(`RSB_M4_DIFFERENCE',`dnl
pushdef(`L1',$1)dnl
pushdef(`L2',$2)dnl
pushdef(`fel',`T')dnl
dnl
dnl
foreach(`E1',L1,`dnl
dnl E1 L2 : RSB_M4_LIST_MEMBER(E1,L2)
ifelse(RSB_M4_LIST_MEMBER(E1,L2),`1',`',`ifelse(fel,`T',`pushdef(`fel',`F')E1',`,E1')')`'dnl
')dnl
ifelse(fel,`T',popdef(`fel'),popdef(`fel')popdef(`fel'))dnl
popdef(`L1')dnl
popdef(`L2')dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_INTERSECTION',`dnl
pushdef(`L1',$1)dnl
pushdef(`L2',$2)dnl
pushdef(`fel',`T')dnl
dnl
foreach(`E1',L1,`dnl
foreach(`E2',L2,`dnl
ifelse(E1,E2,`dnl
ifelse(fel,`T',`E1`'pushdef(`fel',`F')',``,'E1`'')dnl
')dnl
')')dnl
dnl
ifelse(fel,`T',popdef(`fel'),popdef(`fel')popdef(`fel'))dnl
popdef(`L1')dnl
popdef(`L2')dnl
')dnl
dnl
define(`RSB_M4_SBLAS_MATRIX_SUPPORTED_TYPES',`dnl
dnl
RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_PSBLAS_MATRIX_SUPPORTED_TYPES',`dnl
dnl
(RSB_M4_INTERSECTION(RSB_M4_SPBLAS_MATRIX_ALL_TYPES,(WANT_TYPES)))dnl
dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST',`dnl
dnl
RSB_M4_INTERSECTION(RSB_M4_SPBLAS_MATRIX_ALL_TYPES,(WANT_TYPES))dnl
dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST_LENGTH',`dnl
dnl	FIXME: this is broken.
dnl
ifelse(RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST,`',0,RSB_M4_LIST_LENGTH(RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST))dnl
dnl
')dnl
dnl
ifelse(RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST,`',`dnl
define(`RSB_M4_DEFAULT_POSSIBLY_BLAS_TYPE',RSB_M4_INVALID_TYPE)dnl
define(`RSB_M4_DEFAULT_POSSIBLY_BLAS_TYPE_OR_DEFAULT',RSB_M4_DEFAULT_TYPE)dnl
',`
define(`RSB_M4_DEFAULT_POSSIBLY_BLAS_TYPE',`RSB_M4_FIRST(RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST,RSB_M4_INVALID_TYPE)')dnl
define(`RSB_M4_DEFAULT_POSSIBLY_BLAS_TYPE_OR_DEFAULT',`RSB_M4_FIRST(RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST,RSB_M4_INVALID_TYPE)')dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES',`dnl
dnl
(RSB_M4_INTERSECTION(RSB_M4_SPBLAS_MATRIX_ALL_TYPES,(WANT_TYPES)))dnl
dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_UNSUPPORTED_TYPES',`dnl
dnl
(RSB_M4_DIFFERENCE(RSB_M4_SPBLAS_MATRIX_ALL_TYPES,(WANT_TYPES)))dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS',``cr_begin',`cr_block_begin',`cr_variable_block_begin'')dnl
define(`RSB_M4_SPBLAS_MATRIX_END_MOPS',``cr_end'')dnl
define(`RSB_M4_SPBLAS_MATRIX_INSERTION_MOPS',``cr_insert_entry',`cr_insert_entries',`cr_insert_col',`cr_insert_row',`cr_insert_clique',`cr_insert_block'')dnl
define(`RSB_M4_SPBLAS_MATRIX_CREATION_MOPS_LIST',`RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS,RSB_M4_SPBLAS_MATRIX_END_MOPS,RSB_M4_SPBLAS_MATRIX_INSERTION_MOPS')dnl
define(`RSB_M4_SPBLAS_MATRIX_CREATION_MOPS',`(RSB_M4_SPBLAS_MATRIX_CREATION_MOPS_LIST)')dnl
dnl
dnl
dnl	FIXME : level 1 sparse blas is not implemented
define(`RSB_M4_SPBLAS_MATRIX_ALL_L1_MOPS_LIST',``dot',`axpy',`ga',`gz',`sc'')dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L1_MOPS',(RSB_M4_SPBLAS_MATRIX_ALL_L1_MOPS_LIST))dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L2_MOPS_LIST',``mv',`sv'')dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L2_MOPS',(RSB_M4_SPBLAS_MATRIX_ALL_L2_MOPS_LIST))dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L3_MOPS_LIST',``mm',`sm'')dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L3_MOPS',(RSB_M4_SPBLAS_MATRIX_ALL_L3_MOPS_LIST))dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L23_MOPS_LIST',``mv',`sv',`mm',`sm'')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_IMPLEMENTED_CODE_FOR_BLAS_CALL',`dnl
pushdef(`type',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`tri',$3)dnl
dnl
dnl
RSB_M4_LIST_MEMBER(type,RSB_M4_MATRIX_TYPES)dnl
RSB_M4_OR(`dnl
RSB_M4_LIST_MEMBER(RSB_M4_SPBLAS_MOP_CODE_TRANSLATE(mop),RSB_M4_MATRIX_OPS)dnl
dnl
dnl	FIXME : this is fake, always-positive shortcut to get "11" !
,1')dnl
dnl
dnl
popdef(`tri')dnl
popdef(`mop')dnl
popdef(`type')dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION',`dnl
pushdef(`type',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`tri',$3)dnl
pushdef(`want_what',$4)dnl
pushdef(`over',$5)dnl
pushdef(`lang',$6)dnl
pushdef(`args',`$1,$2,$3')dnl
dnl
dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_ALL_L1_MOPS_LIST),`1',`dnl
dnl
ifelse(want_what,`function_declaration',`dnl
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`TYPE',1,lang)` 'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`ID',1,lang)dnl
(RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`ARGS',1,lang));
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`TYPE',1,lang)` 'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`ID',1,lang)dnl
(RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`ARGS',1,lang))
RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(type,mop,tri,`BODY',1,lang)dnl
')dnl
dnl
ifelse(want_what,`BODY',`dnl
{
dnl
ifelse(RSB_M4_IMPLEMENTED_CODE_FOR_BLAS_CALL(type,mop,tri),`11',`dnl
dnl

ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_ALL_L1_MOPS_LIST),`1',`dnl
ifelse($0(type,mop,tri,`TYPE',1,lang),`void',`dnl
RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_SBL1_MSG
	int istatv = $0(type,mop,tri,`ID',1,`lang_c')`'(RSB_M4_FORTRAN_ADDRS_TO_C_VALUES(($0(type,mop,tri,`ARGS',1,`lang_c'))) );
	RSB_SET_IF_NOT_NULL(istat,istatv);
	return;
',`dnl C:
RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_SBL1_MSG
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_SPB_INTERFACE_RETURN(`rsb__BLAS_'X`us'`'mop`'(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),RSB_SPBLAS_TO_RSB_FIX_ARGS(type,mop,tri,`ID',over,lang)))
#else  /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
#endif /* RSB_WANT_SPARSE_BLAS_LEVEL_1 */
')dnl
')`'dnl
dnl
',`dnl
dnl RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_MSG
dnl	return RSB_BLAS_ERROR;
ifelse($0(type,mop,tri,`TYPE',1,lang),`void',`dnl
RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_SBL1_MSG
RSB_M4_SPBLAS_DOC_CONFIGUREDOUT_MSG
	RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
	return;
',`dnl
RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_SBL1_MSG
RSB_M4_SPBLAS_DOC_CONFIGUREDOUT_MSG
	RSB_SPB_INTERFACE_RETURN(RSB_BLAS_ERROR);
')dnl
')dnl
dnl
}
')dnl
dnl
ifelse(want_what,`TYPE',`dnl
ifelse(lang,`f90',`dnl
void`'dnl
',`dnl
int`'dnl
')dnl
')dnl
dnl
ifelse(lang,`f90',`dnl
ifelse(want_what,`ARGS',`dnl
RSB_M4_C_VALUES_TO_FORTRAN_ADDRS(($0(type,mop,tri,`ARGS',1,`lang_c')`'dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`,blas_sparse_matrix A')`'dnl
`,int istat'`'dnl
))`'dnl
')dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
ifelse(lang,`lang_c',`dnl
ifelse(RSB_M4_MEMBER(mop,`dot'),`1',`const enum blas_conj_type conj, const rsb_blas_int_t nnz, const RSB_SPBLAS_OVER_TYPE(type `*',over)x,
		const rsb_blas_int_t *indx, const RSB_SPBLAS_OVER_TYPE(type `*',over)y, const int incy, RSB_SPBLAS_OVER_TYPE(type `*',over)r,
		const enum blas_base_type index_base')`'dnl
ifelse(RSB_M4_MEMBER(mop,`axpy'),`1',`const rsb_blas_int_t nnz, RSB_SPBLAS_OVER_TYPE(type,over) alpha, const RSB_SPBLAS_OVER_TYPE(type `*',over)x, const rsb_blas_int_t *indx,
                 RSB_SPBLAS_OVER_TYPE(type `*',over)y, const int incy, const enum blas_base_type index_base')`'dnl
ifelse(RSB_M4_MEMBER(mop,`ga'),`1',`const rsb_blas_int_t nnz, const RSB_SPBLAS_OVER_TYPE(type `*',over)y, const int incy, RSB_SPBLAS_OVER_TYPE(type `*',over)x, const rsb_blas_int_t *indx,
              const enum blas_base_type index_base')`'dnl
ifelse(RSB_M4_MEMBER(mop,`gz'),`1',`const rsb_blas_int_t nnz, RSB_SPBLAS_OVER_TYPE(type `*',over)y, const int incy, RSB_SPBLAS_OVER_TYPE(type `*',over)x, const rsb_blas_int_t *indx,
              const enum blas_base_type index_base')`'dnl
ifelse(RSB_M4_MEMBER(mop,`sc'),`1',`const rsb_blas_int_t nnz, const RSB_SPBLAS_OVER_TYPE(type `*',over)x, RSB_SPBLAS_OVER_TYPE(type `*',over)y, const int incy, const rsb_blas_int_t *indx,
              const enum blas_base_type index_base')`'dnl
')dnl
')dnl
dnl
dnl
ifelse(want_what,`ID',`dnl
RSB_SPBLAS_FUNCTION_IDENTIFIER(mop,type,lang)`'dnl
')dnl
dnl
')dnl
dnl
popdef(`args')dnl
popdef(`lang')dnl
popdef(`want_what')dnl
popdef(`over')dnl
popdef(`tri')dnl
popdef(`mop')dnl
popdef(`type')dnl
')dnl
dnl
define(`RSB_M4_C_ARG_TYPE',`dnl
patsubst(patsubst($1,`\(.*\)\( \|\*\)\([a-zA-Z_0-9]+\)$',`\1\2'),`\([^ ]*\) *$',`\1')`'dnl
')dnl
dnl
define(`RSB_M4_C_ARG_ID',`dnl
patsubst(patsubst($1,`\(.*\)\( \|\*\)\([a-zA-Z_0-9]+\)$',`\3'),`\([^ ]*\) *$',`\1')`'dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_FORTRAN_ADDR_TO_C_VAL',`dnl
dnl
ifelse(patsubst(RSB_M4_C_ARG_TYPE($1),`[a-zA-Z_0-9 ]',`'),`*',`',`*')RSB_M4_C_ARG_ID($1)`'dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_C_VALUE_REFS_TO_ADDR',`dnl
dnl
RSB_M4_C_ARG_TYPE($1)`'dnl
ifelse(patsubst(RSB_M4_C_ARG_TYPE($1),`[a-zA-Z_0-9 ]',`'),`*',`',`*')RSB_M4_C_ARG_ID($1)`'dnl
')dnl
dnl
dnl
define(`RSB_M4_C_VALUES_TO_FORTRAN_ADDRS',`dnl
dnl
dnl	WARNING : this is THIN ICE :)
pushdef(`firstarg',`0')dnl
foreach(`arg',`$1',`ifelse(firstarg,`0',`pushdef(`firstarg',1)',`,')`'RSB_M4_C_VALUE_REFS_TO_ADDR(arg)')`'dnl
ifelse(firstarg,`1',`popdef(`firstarg')')dnl
popdef(`firstarg')dnl
')dnl
dnl
define(`RSB_M4_FORTRAN_ADDRS_TO_C_VALUES',`dnl
dnl
dnl	WARNING : this is THIN ICE :)
pushdef(`firstarg',`0')dnl
foreach(`arg',`$1',`ifelse(firstarg,`0',`pushdef(`firstarg',1)',`,')`'RSB_M4_FORTRAN_ADDR_TO_C_VAL(arg)')`'dnl
ifelse(firstarg,`1',`popdef(`firstarg')')dnl
popdef(`firstarg')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION',`dnl
dnl
pushdef(`type',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`tri',$3)dnl
pushdef(`want_what',$4)dnl
pushdef(`over',$5)dnl
pushdef(`lang',$6)dnl
pushdef(`args',`$1,$2,$3')dnl
dnl
dnl	FIXME
dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_ALL_L23_MOPS_LIST),`1',`dnl
dnl
dnl
ifelse(want_what,`function_declaration',`dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`TYPE',1,lang)` 'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`ID',1,lang)dnl
(RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`ARGS',1,lang));
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`TYPE',1,lang)` 'dnl
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`ID',1,lang)dnl
(RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`ARGS',1,lang))
RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`BODY',1,lang)dnl
')dnl
dnl
ifelse(want_what,`BODY',`dnl
{
dnl
ifelse(lang,`lang_c',`dnl
dnl
ifelse(RSB_M4_IMPLEMENTED_CODE_FOR_BLAS_CALL(type,mop,tri),`11',`dnl
dnl
dnl
ifelse(RSB_M4_MEMBER(mop,`mv'),`1',`dnl
RSB_M4_SPBLAS_DOC_COMMENT
{
	const type beta = RSB_M4_ONE(type);
dnl	const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(A);
dnl	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb_spmv(rsb__blas_trans_to_rsb_trans(transA),RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)alpha,mtxAp,x,incx,&beta,y,incy)))
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmv(transA,RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)alpha,A,x,incx,&beta,y,incy))
}
	')`'dnl
dnl
dnl	FIXME : no & operator should be used when type is complex !
dnl
ifelse(RSB_M4_MEMBER(mop,`sv'),`1',`dnl
RSB_M4_SPBLAS_DOC_COMMENT
dnl
dnl	FIXME : no & operator should be used when type is complex !
dnl
{
	const struct rsb_mtx_t *mtxAp = rsb__BLAS_inner_matrix_retrieve(T);
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsv(rsb__blas_trans_to_rsb_trans(transT),RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)alpha,mtxAp,x,incx,x,incx)))
}
	')`'dnl
dnl
ifelse(RSB_M4_MEMBER(mop,`mm'),`1',`dnl
RSB_M4_SPBLAS_DOC_COMMENT
{
	const type beta = RSB_M4_ONE(type);
dnl	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spmm(rsb__blas_trans_to_rsb_trans(transA),RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)alpha,rsb__BLAS_inner_matrix_retrieve(A),nrhs,rsb__blas_order_to_rsb_order(order),b,ldb,&beta,c,ldc,RSB_OP_FLAG_DEFAULT)))
dnl	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spmm(rsb__BLAS_inner_matrix_retrieve(A),b,c,ldb,ldc,nrhs,rsb__blas_trans_to_rsb_trans(transA),RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)alpha,&beta,rsb__blas_order_to_rsb_order(order),RSB_OP_FLAG_DEFAULT)))
	RSB_SPB_INTERFACE_RETURN(rsb__BLAS_Xusmm(transA,RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)alpha,A,b,ldb,&beta,c,ldc,nrhs,order))
}
	')`'dnl
dnl
ifelse(RSB_M4_MEMBER(mop,`sm'),`1',`dnl
RSB_M4_SPBLAS_DOC_COMMENT
{
	const type beta = RSB_M4_ZERO(type);
	RSB_SPB_INTERFACE_RETURN(RSB_ERROR_TO_BLAS_ERROR(rsb__do_spsm(rsb__blas_trans_to_rsb_trans(transT),RSB_SPBLAS_OVER_TYPE_ARGVAR_REFERENCE(type)alpha,rsb__BLAS_inner_matrix_retrieve(T),nrhs,rsb__blas_order_to_rsb_order(order),&beta,b,ldb,b,ldb)))
}
	')`'dnl
dnl
',`dnl
dnl	/* FIXME : missing implementation */
dnl	RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_MSG
	/*!
	 RSB_M4_SPBLAS_HELP_INFO(mop)
	*/
RSB_M4_SPBLAS_DOC_CONFIGUREDOUT_MSG
	return RSB_BLAS_ERROR;
')dnl
dnl
dnl
',`dnl
dnl
RSB_M4_SPBLAS_DOC_COMMENT
	int istatv = $0(type,mop,tri,`ID',1,`lang_c')`'(RSB_M4_FORTRAN_ADDRS_TO_C_VALUES((RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(type,mop,tri,`ARGS',1,`lang_c'))));
	RSB_SET_IF_NOT_NULL(istat,istatv);
dnl
')dnl
dnl
dnl
}
')dnl
dnl
ifelse(want_what,`TYPE',`dnl
ifelse(lang,`f90',`dnl
void`'dnl
',`dnl
int`'dnl
')dnl
')dnl
dnl
ifelse(lang,`f90',`dnl
ifelse(want_what,`ARGS',`dnl
dnl $0(type,mop,tri,`ARGS',1,`lang_c')`'dnl
RSB_M4_C_VALUES_TO_FORTRAN_ADDRS(($0(type,mop,tri,`ARGS',1,`lang_c')`'dnl
ifelse(`RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS)',`1',`,blas_sparse_matrix A')`'dnl
`,int istat'`'dnl
))`'dnl
')dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
dnl
ifelse(lang,`lang_c',`dnl
dnl
ifelse(RSB_M4_MEMBER(mop,`mv'),`1',`const enum blas_trans_type transA, RSB_SPBLAS_OVER_TYPE(type,over)alpha,
    const blas_sparse_matrix A, const RSB_SPBLAS_OVER_TYPE(type `*',over)x, const int incx, RSB_SPBLAS_OVER_TYPE(type `*',over)y, const int incy')`'dnl
ifelse(RSB_M4_MEMBER(mop,`sv'),`1',`enum blas_trans_type transT, RSB_SPBLAS_OVER_TYPE(type,over)alpha,
    const blas_sparse_matrix T, RSB_SPBLAS_OVER_TYPE(type `*',over)x, const int incx')`'dnl
ifelse(RSB_M4_MEMBER(mop,`mm'),`1',`const enum blas_order_type order, const enum blas_trans_type transA,
   const int nrhs, RSB_SPBLAS_OVER_TYPE(type,over)alpha, const blas_sparse_matrix A, const RSB_SPBLAS_OVER_TYPE(type `*',over)b, const int ldb,
       RSB_SPBLAS_OVER_TYPE(type `*',over) c, const int ldc')`'dnl
ifelse(RSB_M4_MEMBER(mop,`sm'),`1',`const enum blas_order_type order, const enum blas_trans_type transT,
               const int nrhs, RSB_SPBLAS_OVER_TYPE(type,over)alpha, const blas_sparse_matrix T, RSB_SPBLAS_OVER_TYPE(type `*',over)b, const int ldb')`'dnl
')dnl
')dnl
dnl
ifelse(want_what,`ID',`dnl
RSB_SPBLAS_FUNCTION_IDENTIFIER(mop,type,lang)`'dnl
')dnl
dnl
dnl
')dnl
dnl
popdef(`args')dnl
popdef(`want_what')dnl
popdef(`over')dnl
popdef(`tri')dnl
popdef(`mop')dnl
popdef(`type')dnl
popdef(`lang')dnl
dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS',`dnl
pushdef(`type',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`tri',$3)dnl
pushdef(`want_what',$4)dnl
pushdef(`over',$5)dnl
pushdef(`lang',$6)dnl
pushdef(`args',`$1,$2,$3')dnl
dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_CREATION_MOPS_LIST,`ds',`sp',`gp',`cr_end'),`1',`dnl
dnl
dnl
ifelse(want_what,`function_declaration',`dnl
$0(type,mop,tri,`TYPE',1,lang)` 'dnl
$0(type,mop,tri,`ID',1,lang)dnl
($0(type,mop,tri,`ARGS',1,lang));
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
$0(type,mop,tri,`TYPE',1,lang)` 'dnl
$0(type,mop,tri,`ID',1,lang)dnl
( $0(type,mop,tri,`ARGS',1,lang) )
$0(type,mop,tri,`BODY',1,lang)dnl
')dnl
dnl
ifelse(want_what,`BODY',`dnl
{
dnl
dnl
ifelse(RSB_M4_OR(RSB_M4_LIST_MEMBER(type,RSB_M4_MATRIX_TYPES),RSB_M4_AND(RSB_M4_SAME(type,`'),RSB_M4_LIST_MEMBER(mop,(`cr_end',`ds',`sp',`gp')))),`1',`dnl
dnl
dnl
RSB_M4_SPBLAS_DOC_COMMENT
dnl
ifelse(lang,`f90',`dnl
	int istatv = $0(type,mop,tri,`ID',1,`lang_c')`'dnl
(RSB_M4_FORTRAN_ADDRS_TO_C_VALUES(($0(type,mop,tri,`ARGS',1,`lang_c'))) );
	RSB_SET_IF_NOT_NULL(istat,istatv);
dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`dnl
	RSB_SET_IF_NOT_NULL(A,istatv);
	if(*A && (*A != RSB_BLAS_INVALID_VAL))
	{
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_NO_ERROR);
		rsb__BLAS_ussp(*A,blas_one_base);
	}
	else
		RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_ERROR);
')`'dnl
dnl
',`dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`dnl
	RSB_SPB_INTERFACE_RETURN_HDL(`rsb__BLAS_'X`us'`'mop`'(RSB_SPBLAS_TO_RSB_FIX_ARGS(type,mop,tri,`ID',over,lang)))
',`dnl
	RSB_SPB_INTERFACE_RETURN(`rsb__BLAS_'X`us'`'mop`'(RSB_SPBLAS_TO_RSB_FIX_ARGS(type,mop,tri,`ID',over,lang)))
')dnl
')dnl
dnl
',`dnl
dnl
dnl	RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_MSG
	RSB_M4_SPBLAS_DOC_CONFIGUREDOUT_MSG
	/*!
          RSB_M4_SPBLAS_HELP_INFO(mop)
	 */
dnl	/* FIXME : missing implementation */
dnl	return RSB_BLAS_ERROR;
ifelse(lang,`f90',`dnl
	RSB_SET_IF_NOT_NULL(istat,RSB_BLAS_INVALID_VAL);
',`dnl
	return RSB_BLAS_INVALID_VAL;
')dnl
dnl
')dnl
}
')dnl
dnl
ifelse(want_what,`TYPE',`dnl
ifelse(lang,`f90',`dnl
void`'dnl
',`dnl
ifelse(RSB_M4_MEMBER(mop,`cr_begin',`cr_block_begin',`cr_variable_block_begin'),`1',`blas_sparse_matrix',`int')`'dnl
')dnl
')dnl
dnl
ifelse(lang,`f90',`dnl
ifelse(want_what,`ARGS',`dnl
RSB_M4_C_VALUES_TO_FORTRAN_ADDRS(($0(type,mop,tri,`ARGS',1,`lang_c')`'dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`,blas_sparse_matrix A')`'dnl
`,int istat'`'dnl
))`'dnl
')dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
ifelse(lang,`lang_c',`dnl
dnl
dnl
ifelse(RSB_M4_MEMBER(mop,`cr_begin'),`1',`rsb_blas_int_t m, rsb_blas_int_t n')`'dnl
ifelse(RSB_M4_MEMBER(mop,`cr_block_begin'),`1',`rsb_blas_int_t Mb, rsb_blas_int_t Nb, rsb_blas_int_t k, rsb_blas_int_t l')`'dnl
ifelse(RSB_M4_MEMBER(mop,`cr_variable_block_begin'),`1',`rsb_blas_int_t Mb, rsb_blas_int_t Nb,
		const rsb_blas_int_t *K, const rsb_blas_int_t *L')`'dnl
dnl
ifelse(RSB_M4_MEMBER(mop,`cr_insert_entry'),`1',`blas_sparse_matrix A, RSB_SPBLAS_OVER_TYPE(type,over) val, rsb_blas_int_t i, rsb_blas_int_t j')`'dnl
dnl	FIXME : complex cr_insert_entry originally has not const pointers !?
ifelse(RSB_M4_MEMBER(mop,`cr_insert_entries'),`1',`blas_sparse_matrix A, rsb_blas_int_t nnz, const RSB_SPBLAS_OVER_TYPE(type `*',over)val,
                            const rsb_blas_int_t *indx, const rsb_blas_int_t *jndx')`'dnl
ifelse(RSB_M4_MEMBER(mop,`cr_insert_col'),`1',`blas_sparse_matrix A, rsb_blas_int_t j, rsb_blas_int_t nnz,
                           const RSB_SPBLAS_OVER_TYPE(type `*',over)val, const rsb_blas_int_t *indx')`'dnl
ifelse(RSB_M4_MEMBER(mop,`cr_insert_row'),`1',`blas_sparse_matrix A, rsb_blas_int_t i, rsb_blas_int_t nnz,
                           const RSB_SPBLAS_OVER_TYPE(type `*',over)val, const rsb_blas_int_t *indx')`'dnl
ifelse(RSB_M4_MEMBER(mop,`cr_insert_clique'),`1',`blas_sparse_matrix A, const rsb_blas_int_t k, const rsb_blas_int_t l,
                       const RSB_SPBLAS_OVER_TYPE(type `*',over)val, const rsb_blas_int_t row_stride,
                       const rsb_blas_int_t col_stride, const rsb_blas_int_t *indx,
                       const rsb_blas_int_t *jndx')`'dnl
ifelse(RSB_M4_MEMBER(mop,`cr_insert_block'),`1',`blas_sparse_matrix A, const RSB_SPBLAS_OVER_TYPE(type `*',over)val,
                        rsb_blas_int_t row_stride, rsb_blas_int_t col_stride, rsb_blas_int_t i, rsb_blas_int_t j')`'dnl
ifelse(RSB_M4_MEMBER(mop,`cr_end'),`1',`blas_sparse_matrix A')`'dnl
ifelse(RSB_M4_MEMBER(mop,`ds'),`1',`blas_sparse_matrix A')`'dnl
ifelse(RSB_M4_MEMBER(mop,`sp'),`1',`blas_sparse_matrix A, rsb_blas_int_t pname')`'dnl
ifelse(RSB_M4_MEMBER(mop,`gp'),`1',`blas_sparse_matrix A, rsb_blas_int_t pname')`'dnl
')dnl
dnl
')dnl
dnl
ifelse(want_what,`ID',`dnl
RSB_SPBLAS_FUNCTION_IDENTIFIER(mop,type,lang)`'dnl
')dnl
dnl
')dnl
dnl
popdef(`args')dnl
popdef(`over')dnl
popdef(`want_what')dnl
popdef(`tri')dnl
popdef(`mop')dnl
popdef(`type')dnl
popdef(`lang')dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_SPBLAS_EXTRA_FUNCTION',`dnl
pushdef(`type',$1)dnl
pushdef(`mop',$2)dnl
pushdef(`tri',$3)dnl
pushdef(`want_what',$4)dnl
pushdef(`over',$5)dnl
pushdef(`lang',$6)dnl
pushdef(`args',`$1,$2,$3')dnl
dnl
ifelse(RSB_M4_LIST_MEMBER(mop,RSB_M4_SBLAS_EXTRA_INTERFACE_OPS),`1',`dnl
dnl
dnl
ifelse(want_what,`function_declaration',`dnl
$0(type,mop,tri,`TYPE',1,lang)` 'dnl
$0(type,mop,tri,`ID',1,lang)dnl
($0(type,mop,tri,`ARGS',1,lang));
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
$0(type,mop,tri,`TYPE',1,lang)` 'dnl
$0(type,mop,tri,`ID',1,lang)dnl
( $0(type,mop,tri,`ARGS',1,lang) )
$0(type,mop,tri,`BODY',1,lang)dnl
')dnl
dnl
ifelse(want_what,`BODY',`dnl
{
dnl
dnl
ifelse(RSB_M4_LIST_MEMBER(mop,RSB_M4_SBLAS_EXTRA_INTERFACE_OPS),`1',`dnl
dnl
dnl
RSB_M4_SPBLAS_DOC_COMMENT
dnl
ifelse(lang,`f90',`dnl
dnl	RSB_M4_SPBLAS_DOC_CONFIGUREDOUT_MSG
	int istatv = $0(type,mop,tri,`ID',1,`lang_c')`'dnl
(RSB_M4_FORTRAN_ADDRS_TO_C_VALUES(($0(type,mop,tri,`ARGS',1,`lang_c'))) );
	RSB_SET_IF_NOT_NULL(istat,istatv);
dnl
dnl
',`dnl
dnl	RSB_M4_SPBLAS_DOC_CONFIGUREDOUT_MSG
	RSB_SPB_INTERFACE_RETURN(`rsb__BLAS_'X`us'`'mop`'(RSB_SPBLAS_TO_RSB_FIX_ARGS(type,mop,tri,`ID',over,lang)))
')dnl
dnl
',`dnl
dnl
dnl	/* FIXME : missing implementation */
RSB_M4_SPBLAS_DOC_UNIMPLEMENTED_MSG
	return RSB_BLAS_ERROR;
dnl
')dnl
}
')dnl
dnl
ifelse(want_what,`TYPE',`dnl
ifelse(lang,`f90',`dnl
void`'dnl
',`dnl
int`'dnl
')dnl
')dnl
dnl
ifelse(lang,`f90',`dnl
ifelse(want_what,`ARGS',`dnl
RSB_M4_C_VALUES_TO_FORTRAN_ADDRS(($0(type,mop,tri,`ARGS',1,`lang_c')`'dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`,blas_sparse_matrix A')`'dnl
`,int istat'`'dnl
))`'dnl
')dnl
')dnl
dnl
ifelse(want_what,`ARGS',`dnl
ifelse(lang,`lang_c',`dnl
dnl
dnl
`blas_sparse_matrix A'dnl
ifelse(RSB_M4_MEMBER(mop,`rows_scale'),`1',`,const RSB_SPBLAS_OVER_TYPE(type *,over) d, const enum blas_trans_type trans')`'dnl
ifelse(RSB_M4_MEMBER(mop,`get_diag'),`1',`,RSB_SPBLAS_OVER_TYPE(type *,over) d')`'dnl
ifelse(RSB_M4_MEMBER(mop,`get_rows_sparse'),`1',`, RSB_SPBLAS_OVER_TYPE(type *,over) VA, rsb_blas_int_t * IA, rsb_blas_int_t * JA, rsb_blas_int_t * nnz, const rsb_blas_int_t fr, const rsb_blas_int_t lr')`'dnl
ifelse(RSB_M4_MEMBER(mop,`get_rows_nnz'),`1',`, const rsb_blas_int_t fr, const rsb_blas_int_t lr, rsb_blas_int_t * nnzp')`'dnl
ifelse(RSB_M4_MEMBER(mop,`get_matrix_nnz'),`1',`,rsb_blas_int_t * nnz')`'dnl
ifelse(RSB_M4_MEMBER(mop,`get_infinity_norm'),`1',`,RSB_SPBLAS_OVER_TYPE(type *, over)in, const enum blas_trans_type trans')`'dnl
ifelse(RSB_M4_MEMBER(mop,`set_elements'),`1',`,const rsb_blas_int_t * ia, const rsb_blas_int_t *ja, const RSB_SPBLAS_OVER_TYPE(type *,over) va, const rsb_blas_int_t nnz')`'dnl
ifelse(RSB_M4_MEMBER(mop,`get_element'),`1',`,const rsb_blas_int_t i, const rsb_blas_int_t j, RSB_SPBLAS_OVER_TYPE(type *,over) v')`'dnl
ifelse(RSB_M4_MEMBER(mop,`set_element'),`1',`,const rsb_blas_int_t i, const rsb_blas_int_t j, RSB_SPBLAS_OVER_TYPE(type *,over) v')`'dnl
dnl
')dnl
dnl
')dnl
dnl
ifelse(want_what,`ID',`dnl
RSB_SPBLAS_FUNCTION_IDENTIFIER(mop,type,lang)`'dnl
')dnl
dnl
')dnl
dnl
popdef(`args')dnl
popdef(`over')dnl
popdef(`want_what')dnl
popdef(`tri')dnl
popdef(`mop')dnl
popdef(`type')dnl
popdef(`lang')dnl
')dnl
dnl
