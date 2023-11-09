dnl
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
include(`libspblas_macros.m4')dnl
dnl
dnl
define(`RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_RSB_INTERFACE',`dnl
dnl RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_PSB_INTERFACE($@)`'dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS($@)`'dnl
')dnl
dnl
define(`RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_SB_INTERFACE',`dnl
dnl RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_PSB_INTERFACE($@)`'dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS($@)`'dnl
')dnl
dnl
define(`RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_PSB_INTERFACE',`dnl
dnl
dnl	FIXME
dnl
patsubst(dnl
`patsubst(dnl
`patsubst(dnl
`patsubst(dnl
`patsubst(dnl
`patsubst(dnl
`patsubst(dnl
`RSB_M4_ARGS_TO_ACTUAL_ARGS($@)',dnl
`\<a\>',`a%infoa(psb_const_infoa_rsb_)'dnl
)',dnl
`\<append\>',`appendi'dnl
)',dnl
`\<has_iren\>',`has_ireni'dnl
)',dnl
`\<has_gtl\>',`has_gtli'dnl
)',dnl
`\<do_rebuild\>',`do_rebuildi'dnl
)',dnl
`\<real_in\>',`in'dnl
)',dnl
`\<extra\>',`extra,typecode'dnl
)dnl
dnl
')dnl
dnl
define(`RSB_M4_PSBLAS_IARRAY_TYPE',`dnl
INTEGER :: dnl
pushdef(`firstarg',`0')dnl
foreach(`arg',`($@)',`ifelse(firstarg,`0',`pushdef(`firstarg',1)arg`(:)'',`,arg`(:)'')`'')`'
dnl	DO NOT REMOVE THE FOLLOWING LINE
ifelse(firstarg,`1',`popdef(`firstarg')')dnl
popdef(`firstarg')dnl
')dnl
dnl
define(`RSB_M4_C2F_NORM_TYPE',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(RSB_M4_want_old_fortran_float_types,`1',`dnl
ifelse(type,`float',`REAL*4')`'dnl
ifelse(type,`double',`REAL*8')`'dnl
ifelse(type,`float complex',`REAL*4')`'dnl
ifelse(type,`double complex',`REAL*8')`'dnl
',`dnl
dnl ifelse(type,`float',`REAL(rsb_spk_)')`'dnl
dnl ifelse(type,`double',`REAL(rsb_dpk_)')`'dnl
dnl ifelse(type,`float complex',`REAL(rsb_spk_)')`'dnl
dnl ifelse(type,`double complex',`REAL(rsb_dpk_)')`'dnl
ifelse(type,`float',`REAL(KIND(1.e0))')`'dnl
ifelse(type,`double',`REAL(KIND(1.d0))')`'dnl
ifelse(type,`float complex',`REAL(KIND(1.e0))')`'dnl
ifelse(type,`double complex',`REAL(KIND(1.d0))')`'dnl
')`'dnl
dnl
ifelse(type,`int',`INTEGER(KIND=RSB_BLAS_IDX_KIND)')`'dnl lda,ldb,nrhs,...
dnl FIXME : and, for other, non canonical types ? FIXME 
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_C2F_TYPE',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(RSB_M4_want_old_fortran_float_types,`1',`dnl
ifelse(type,`float',`REAL*4')`'dnl
ifelse(type,`double',`REAL*8')`'dnl
ifelse(type,`float complex',`COMPLEX*8')`'dnl
ifelse(type,`double complex',`COMPLEX*16')`'dnl
',`dnl
dnl ifelse(type,`float',`REAL(rsb_spk_)')`'dnl
dnl ifelse(type,`double',`REAL(rsb_dpk_)')`'dnl
dnl ifelse(type,`float complex',`COMPLEX(rsb_spk_)')`'dnl
dnl ifelse(type,`double complex',`COMPLEX(rsb_dpk_)')`'dnl
ifelse(type,`float',`REAL(KIND(1.e0))')`'dnl
ifelse(type,`double',`REAL(KIND(1.d0))')`'dnl
ifelse(type,`float complex',`COMPLEX(KIND(1.e0))')`'dnl
ifelse(type,`double complex',`COMPLEX(KIND(1.d0))')`'dnl
')`'dnl
dnl
ifelse(type,`int',`INTEGER')`'dnl
dnl FIXME : and, for other, non canonical types ? FIXME 
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_PSBLAS_INTERFACE_RADIX',`psb_rsb')dnl
define(`RSB_M4_RSBLAS_INTERFACE_RADIX',`rsb')dnl
define(`RSB_M4_SBLAS_INTERFACE_RADIX',`us')dnl unstructured sparse
dnl
define(`RSB_M4_FORTRAN_INTERFACE_RADIX',`rsb__')dnl
dnl
dnl define(`RSB_M4_SBLAS_INTERFACE_OPS',`RSB_M4_PSBLAS_INTERFACE_OPS')dnl
dnl
define(`RSB_M4_SBLAS_INTERFACE_OPS',`(dnl
RSB_M4_COMMA_LIST(RSB_M4_SPBLAS_MATRIX_CREATION_MOPS)`'dnl
,RSB_M4_COMMA_LIST(RSB_M4_SPBLAS_MATRIX_ALL_L2_MOPS)`'dnl
,RSB_M4_COMMA_LIST(RSB_M4_SPBLAS_MATRIX_ALL_L3_MOPS)`'dnl
)')dnl
dnl
define(`RSB_M4_SBLAS_GENERIC_OPS',`(dnl
RSB_M4_COMMA_LIST((RSB_M4_SPBLAS_MATRIX_INSERTION_MOPS))`'dnl
,RSB_M4_COMMA_LIST(RSB_M4_SPBLAS_MATRIX_ALL_L2_MOPS)`'dnl
,RSB_M4_COMMA_LIST(RSB_M4_SPBLAS_MATRIX_ALL_L3_MOPS)`'dnl
)')dnl
dnl dnl
define(`RSB_M4_RSBLAS_INTERFACE_OPS',`RSB_M4_PSBLAS_INTERFACE_OPS')dnl
dnl
define(`RSB_M4_PSBLAS_INTERFACE_OPS',`(scale,getdiag,get_rows_nnz,get_rows_sparse,destroy_sparse_matrix,allocate_sparse_matrix,get_matrix_nnz,infinity_norm,usmm,ussm,usmv,set_elements,set_element,get_element)')dnl
dnl
dnl
dnl	FIXME: new stuff
dnl
define(`RSB_M4_SBLAS_EXTRA_INTERFACE_OPS_LIST',``rows_scale',`get_diag',`get_rows_nnz',`get_rows_sparse',`get_matrix_nnz',`get_infinity_norm',`set_elements',`set_element',`get_element'')dnl
dnl define(`RSB_M4_SBLAS_EXTRA_INTERFACE_OPS_LIST',`')dnl
define(`RSB_M4_SBLAS_EXTRA_INTERFACE_OPS',`(RSB_M4_SBLAS_EXTRA_INTERFACE_OPS_LIST)')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_RSBLAS_INTERFACE_IDENTIFIER',`dnl
pushdef(`pmop',$1)dnl
dnl
RSB_M4_RSBLAS_INTERFACE_RADIX`_'pmop`'dnl
dnl
popdef(`pmop')dnl
')dnl
dnl
define(`RSB_M4_PSBLAS_INTERFACE_IDENTIFIER',`dnl
pushdef(`pmop',$1)dnl
dnl
RSB_M4_PSBLAS_INTERFACE_RADIX`_'pmop`'dnl
dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
define(`RSB_M4_SBLAS_INTERFACE_IDENTIFIER',`dnl
pushdef(`pmop',$1)dnl
dnl
RSB_M4_SBLAS_INTERFACE_RADIX`'pmop`'dnl
dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
define(`RSB_M4_RSB_TYPE_CHARCODE',`dnl
RSB_M4_PSB_TYPE_CHARCODE($@)`'dnl
')dnl
dnl
define(`RSB_M4_SB_TYPE_CHARCODE',`dnl
RSB_M4_PSB_TYPE_CHARCODE($@)`'dnl
')dnl
dnl
dnl
define(`RSB_M4_PSB_TYPE_CHARCODE',`dnl
pushdef(`type',$1)`'dnl
tolowercase(RSB_M4_TYPE_CHARCODE(type))`'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_PSB_MTYPE_CHARCODE',`dnl
pushdef(`type',$1)`'dnl
`psb_'RSB_M4_PSB_TYPE_CHARCODE(type)`'`spmat_type'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_RSB_MTYPE_CHARCODE',`dnl
dnl
RSB_M4_SB_MTYPE_CHARCODE($@)`'dnl
dnl
')dnl
dnl
define(`RSB_M4_SB_MTYPE_CHARCODE',`dnl
dnl
RSB_M4_SB_MTYPE_CHARCODE($@)`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_RSBLAS2VBR_SUBROUTINE_RADIX',`dnl
dnl
dnl
RSB_M4_FORTRAN_INTERFACE_RADIX`do_fortran_rsb_'`'dnl
dnl RSB_M4_FORTRAN_INTERFACE_RADIX`do_fortran_'RSB_M4_PSB_TYPE_CHARCODE(mtype)`_'pmop`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_SBLAS2VBR_SUBROUTINE_RADIX',`dnl
pushdef(`mop',$1)`'dnl
pushdef(`mtype',$2)`'dnl
pushdef(`lang',$3)`'dnl
dnl
dnl RSB_M4_FORTRAN_INTERFACE_RADIX`do_fortran_rsb_blas_'`'dnl
dnl RSB_M4_FORTRAN_INTERFACE_RADIX`do_fortran_'RSB_M4_PSB_TYPE_CHARCODE(mtype)`_'pmop`'dnl
ifelse(lang,`f90',`dnl
dnl `rsb_'`'f90_blas_`'RSB_M4_PSB_TYPE_CHARCODE(mtype)`us'`'dnl
`'blas_`'RSB_M4_PSB_TYPE_CHARCODE(mtype)`us'`'dnl
',`dnl
BLAS_`'RSB_M4_PSB_TYPE_CHARCODE(mtype)`us'`'dnl
')`'dnl
dnl
popdef(`mtype')`'dnl
popdef(`lang')`'dnl
popdef(`mop')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_PSBLAS2VBR_SUBROUTINE_RADIX',`dnl
dnl
dnl
RSB_M4_FORTRAN_INTERFACE_RADIX`do_fortran_'`'dnl
dnl RSB_M4_FORTRAN_INTERFACE_RADIX`do_fortran_'RSB_M4_PSB_TYPE_CHARCODE(mtype)`_'pmop`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_RSBLAS2VBR_SUBROUTINE_IDENTIFIER',`dnl
dnl
RSB_M4_SBLAS2VBR_SUBROUTINE_RADIX`'pmop`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_SBLAS2VBR_SUBROUTINE_IDENTIFIER',`dnl
pushdef(`mop',$1)`'dnl
pushdef(`mtype',$2)`'dnl
pushdef(`lang',$3)`'dnl
dnl
dnl
dnl
RSB_M4_SBLAS2VBR_SUBROUTINE_RADIX(mop,mtype,lang)`'mop`'dnl
ifelse(lang,`f90',RSB_M4_FORTRAN_SYMBOL_ADD_TO_F,`')`'dnl
dnl
popdef(`lang')`'dnl
popdef(`mtype')`'dnl
popdef(`mop')`'dnl
dnl
')dnl
dnl
define(`RSB_M4_PSBLAS2VBR_SUBROUTINE_IDENTIFIER',`dnl
pushdef(`pmop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
dnl
RSB_M4_PSBLAS2VBR_SUBROUTINE_RADIX`'pmop`'dnl
dnl RSB_M4_FORTRAN_INTERFACE_RADIX`do_fortran_'RSB_M4_PSB_TYPE_CHARCODE(mtype)`_'pmop`'dnl
dnl
popdef(`mtype')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
define(`RSB_M4_SBLAS_SUBROUTINE_IDENTIFIER',`dnl
pushdef(`pmop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
RSB_M4_SB_TYPE_CHARCODE(mtype)`'dnl
RSB_M4_SBLAS_INTERFACE_RADIX`'dnl
pmop`'dnl
dnl
popdef(`mtype')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_RSBLAS_SUBROUTINE_IDENTIFIER',`dnl
pushdef(`pmop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
RSB_M4_RSBLAS_INTERFACE_RADIX`'RSB_M4_RSB_TYPE_CHARCODE(mtype)`'pmop`'dnl
dnl
popdef(`mtype')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_PSBLAS_SUBROUTINE_IDENTIFIER',`dnl
pushdef(`pmop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
RSB_M4_PSBLAS_INTERFACE_RADIX`_'RSB_M4_PSB_TYPE_CHARCODE(mtype)`_'pmop`'dnl
dnl RSB_M4_PSBLAS_INTERFACE_RADIX`_'RSB_M4_C2F_TYPE(mtype)`_'pmop`'dnl
dnl
popdef(`mtype')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_RSBLAS_SUBROUTINE_INFO_DECLARATION',`dnl
RSB_M4_PSBLAS_SUBROUTINE_INFO_DECLARATION($@)`'dnl
dnl
')dnl
dnl
define(`RSB_M4_SBLAS_SUBROUTINE_INFO_DECLARATION',`dnl
pushdef(`id',$1)`'dnl
dnl RSB_M4_PSBLAS_SUBROUTINE_INFO_DECLARATION($@)`'dnl
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::id
dnl
popdef(`id')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_PSBLAS_SUBROUTINE_INFO_DECLARATION',`dnl
dnl
           INTEGER, INTENT(OUT) ::info
dnl
')dnl
dnl
dnl
define(`RSB_M4_RSBLAS_SUBROUTINE_ARGS_DECLARATION',`dnl
RSB_M4_SBLAS_SUBROUTINE_ARGS_DECLARATION()`'dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_C_POINTER_TO_FORTRAN_ARRAY',`dnl
pushdef(`type',$1)`'dnl
pushdef(`id',$2)`'dnl
ifelse(patsubst(id,`[bc]',`-'),`-',`(:,:)',`dnl
dnl patsubst(`patsubst(`patsubst(`$1',`const',`')',` *[a-zA-Z_0-9]+ *[a-zA-Z_0-9]+ *[a-zA-Z_0-9]+ *',`')',` *\* *',`(:)')`'dnl
ifelse(patsubst(patsubst($1,`[a-zA-Z_0-9]*',`'),` *',`'),`',`',`(:)')`'dnl
')dnl
dnl patsubst(`patsubst(`patsubst(`$1',`const',`')',` *[a-zA-Z_0-9]+ *[a-zA-Z_0-9]+ *[a-zA-Z_0-9]+ *',`')',` *\* *',`(:)')`'dnl
popdef(`id')`'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_C_TYPE_TO_FORTRAN_TYPE',`dnl
pushdef(`type',$1)`'dnl
ifelse(type,`blas_sparse_matrix',`INTEGER',`dnl
ifelse(type,`rsb_blas_int_t',`INTEGER(KIND=RSB_BLAS_IDX_KIND)',`dnl
ifelse(type,`rsb_blas_int_t*',`INTEGER(KIND=RSB_BLAS_IDX_KIND)',`dnl
ifelse(type,`enum blas_trans_type',`INTEGER',`dnl
ifelse(type,`enum`'blas_trans_type',`INTEGER',`dnl
ifelse(type,`enum blas_trans_type ',`INTEGER',`dnl
ifelse(type,`enum blas_order_type',`INTEGER',`dnl
ifelse(type,`enum`'blas_order_type',`INTEGER',`dnl FIXME
ifelse(type,`enum blas_order_type*',`INTEGER',`dnl
RSB_M4_C2F_NORM_TYPE(patsubst(type,` *\* *',`'))`'dnl
')`'dnl
')`'dnl
')`'dnl
')`'dnl
')`'dnl
')`'dnl
')`'dnl
')`'dnl
')`'dnl
popdef(`type')`'dnl
')dnl
dnl
define(`RSB_M4_SPBLAS_FIXTYPE',`dnl
pushdef(`type',$1)`'dnl
patsubst(`patsubst(`$2',`\n',`')',`void',type)`'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
dnl
dnl define(`RSB_M4_C_ARG_TO_FORTRAN_ARG',`a')dnl
dnl define(`RSB_M4_C_ARG_TO_FORTRAN_ARG',`$@ ')dnl
define(`RSB_M4_C_ARG_TO_FORTRAN_ARG',`patsubst($1,`\(.*\) *\(\<[a-zA-Z_0-9]+$\)',`dnl
RSB_M4_C_TYPE_TO_FORTRAN_TYPE(`patsubst(`patsubst(\1,` * ',`')',`const ',`')') `::' \2 pushdef(`kiki',\1)`'RSB_M4_C_POINTER_TO_FORTRAN_ARRAY(kiki,\2)`'popdef(`kiki')')')dnl
dnl
dnl
define(`RSB_M4_C_ARGS_TO_FORTRAN_ARGS',`dnl
dnl
dnl	WARNING : this is THIN ICE :)
pushdef(`firstarg',`0')dnl
foreach(`arg',`$1',`ifelse(firstarg,`0',`pushdef(`firstarg',1)',`,')`'RSB_M4_C_ARG_TO_FORTRAN_ARG(arg)')`'dnl
ifelse(firstarg,`1',`popdef(`firstarg')')dnl
popdef(`firstarg')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_SBLAS_SUBROUTINE_ARGS_DECLARATION',`dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
dnl        INTEGER :: a
RSB_M4_SBLAS_SUBROUTINE_DECL(mop,mtype)`'dnl
dnl	RSB_M4_PSBLAS_SUBROUTINE_ARGS_DECLARATION($@)`'dnl
dnl	dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
define(`RSB_M4_PSBLAS_SUBROUTINE_ARGS_DECLARATION',`dnl
pushdef(`pmop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
          type(RSB_M4_PSB_MTYPE_CHARCODE(mtype)) :: a ! psblas matrix type (mtype)
dnl
dnl
ifelse(pmop,`scale',`dnl
          RSB_M4_C2F_TYPE(mtype) :: d(:)
')dnl
ifelse(pmop,`getdiag',`dnl
          RSB_M4_C2F_TYPE(mtype) :: d(:)
')dnl
ifelse(pmop,`get_element',`dnl
          INTEGER :: i,j
          RSB_M4_C2F_TYPE(mtype) :: v
')dnl
ifelse(pmop,`set_elements',`dnl
          INTEGER :: ia(:),ja(:)
          INTEGER :: imin,imax,jmin,jmax,nnz
          INTEGER :: has_gtli,gtl(:)
dnl ,do_rebuildi
          logical :: has_gtl
dnl ,do_rebuild=.FALSE.
          RSB_M4_C2F_TYPE(mtype) :: val(:)
')dnl
ifelse(pmop,`set_element',`dnl
          INTEGER :: i,j
          RSB_M4_C2F_TYPE(mtype) :: v
')dnl
ifelse(pmop,`get_rows_nnz',`dnl
          INTEGER :: fr,lr
          INTEGER :: nnz
')`'dnl
ifelse(pmop,`get_rows_sparse',`dnl
          RSB_M4_C2F_TYPE(mtype) :: val(:)
          INTEGER :: fr,lr
          logical :: append
          logical :: has_iren
          INTEGER :: appendi=0
dnl          INTEGER :: nzin ! FIXME, still unused !
          INTEGER :: nzin,has_ireni=0,iren(:) ! FIXME, still unused !
          RSB_M4_PSBLAS_IARRAY_TYPE(ia,ja)dnl
          INTEGER :: nnz
')`'dnl
ifelse(pmop,`destroy_sparse_matrix',`dnl
')`'dnl
ifelse(pmop,`allocate_sparse_matrix',`dnl
          RSB_M4_C2F_TYPE(mtype) :: val(:)
          RSB_M4_PSBLAS_IARRAY_TYPE(ia,ja)dnl
          INTEGER :: nnz,m,k,br,bc
dnl          character(len=16),parameter :: extra="a"
          character(len=*) :: extra
          INTEGER,parameter :: typecode = RSB_M4_TYPE_CHARCODE_ASCII_VALUE(mtype)
dnl          !extra=""
')`'dnl
ifelse(pmop,`get_matrix_nnz',`dnl
          INTEGER :: ires
')`'dnl
ifelse(pmop,`infinity_norm',`dnl
          RSB_M4_C2F_NORM_TYPE(mtype) :: real_in
          RSB_M4_C2F_TYPE(mtype) :: in
          character :: trans ! 
dnl          INTEGER :: itrans ! FIXME: this value is ignored !
')`'dnl
ifelse(pmop,`ussm',`dnl
          RSB_M4_C2F_TYPE(mtype) :: b(:,:),c(:,:)
dnl          INTEGER :: lb,lc,nrhs ! FIXME : new, still unused! (FIXME: nrhs is PSBLAS s iwsz)
          INTEGER :: lb,lc,nrhs ! FIXME : new, still unused! (FIXME: nrhs is iwsz)
          RSB_M4_C2F_TYPE(mtype) :: alpha
          character :: trans ! 
dnl          INTEGER :: itrans ! FIXME: this value is ignored !
')`'dnl
ifelse(pmop,`usmv',`dnl
          RSB_M4_C2F_TYPE(mtype) :: b(:),c(:)
dnl          INTEGER :: lb,lc,nrhs ! FIXME : new, still unused! (FIXME: nrhs is PSBLAS s iwsz)
          RSB_M4_C2F_TYPE(mtype) :: alpha,beta
          character :: trans ! 
dnl          INTEGER :: itrans ! FIXME: this value is ignored !
')`'dnl
dnl
ifelse(pmop,`usmm',`dnl
          RSB_M4_C2F_TYPE(mtype) :: b(:,:),c(:,:)
          INTEGER :: lb,lc,nrhs ! FIXME : new, still unused! (FIXME: nrhs is iwsz)
dnl          INTEGER :: lb,lc,nrhs ! FIXME : new, still unused! (FIXME: nrhs is PSBLAS s iwsz)
          RSB_M4_C2F_TYPE(mtype) :: alpha,beta
          character :: trans ! 
dnl          INTEGER :: itrans ! FIXME: this value is ignored !
')`'dnl
dnl
popdef(`mtype')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
define(`RSB_M4_RSBLAS_SUBROUTINE_HELP_COMMENT',`dnl
RSB_M4_PSBLAS_SUBROUTINE_HELP_COMMENT($@)`'dnl
')dnl
define(`RSB_M4_SBLAS_SUBROUTINE_HELP_COMMENT',`dnl
RSB_M4_PSBLAS_SUBROUTINE_HELP_COMMENT($@)`'dnl
')dnl
dnl
define(`RSB_M4_SBLAS_SUBROUTINE_EXTRA_FORTRAN_HELP_COMMENT',`dnl
pushdef(`pmop',$1)dnl
\rsb_spblas_f_istat_msg`'dnl
ifelse(RSB_M4_MEMBER(pmop,`cr_begin',`cr_block_begin',`cr_variable_block_begin'),`1',`dnl
\rsb_spblasl2_A_msg_ftn
')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
define(`RSB_M4_PSBLAS_SUBROUTINE_HELP_COMMENT',`dnl
pushdef(`pmop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
ifelse(pmop,`getdiag',`gets the diagonal of A into array D')`'dnl
ifelse(pmop,`set_element',`sets the matrix at element i,j (if in the nonzero pattern)')`'dnl
ifelse(pmop,`set_elements',`sets the matrix at elements ia,ja (if in the nonzero pattern)')`'dnl
ifelse(pmop,`get_element',`gets the matrix at element i,j (if in the nonzero pattern)')`'dnl
ifelse(pmop,`scale',`scales each row of A by multiplying it to a value of D')`'dnl
ifelse(pmop,`get_rows_sparse',`writes in ia,ja,va the row index, column index, and value of nonzeros from row fr to lr')`'dnl
ifelse(pmop,`get_rows_nnz',`gets the number of nonzeros in the specified rows interval')`'dnl
ifelse(pmop,`allocate_sparse_matrix',`allocates a sparse matrix A')`'dnl
ifelse(pmop,`destroy_sparse_matrix',`frees all allocated resources to the descriptor of matrix A')`'dnl
ifelse(pmop,`get_matrix_nnz',`gets the nonzeros count of matrix A')`'dnl
ifelse(pmop,`infinity_norm',`gets the infinity norm (the maximal sum of rows elements) of matrix A')`'dnl
ifelse(pmop,`sm',`triangular solve: b <- alpha A^-1 b')`'dnl
ifelse(pmop,`sv',`triangular solve: b <- alpha A^-1 b')`'dnl
ifelse(pmop,`mm',`multiplication  : c <- beta c + alpha A b')`'dnl
ifelse(pmop,`mv',`multiplication  : c <- beta c + alpha A b')`'dnl
dnl
ifelse(pmop,`cr',`matrix creation')`'dnl
ifelse(pmop,`cr_insert_row',`inserts a sparse row')`'dnl
ifelse(pmop,`cr_insert_col',`inserts a sparse column')`'dnl
ifelse(pmop,`cr_insert_block',`inserts a dense block')`'dnl
ifelse(pmop,`cr_insert_clique',`inserts a clique')`'dnl
ifelse(pmop,`cr_insert_entry',`inserts a single entry')`'dnl
ifelse(pmop,`cr_insert_entries',`inserts multiple entries')`'dnl
dnl
popdef(`mtype')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_RSBLAS_SUBROUTINE_ARGS',`dnl
RSB_M4_PSBLAS_SUBROUTINE_ARGS($@)`'dnl
')dnl
dnl
dnl
define(`RSB_M4_SBLAS_SUBROUTINE_DECL',`dnl
pushdef(`mop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
dnl
foreach(`d',`(dnl
RSB_M4_C_ARGS_TO_FORTRAN_ARGS((`RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(mtype,mop,`',`ARGS',0,`lang_c')'))dnl
RSB_M4_C_ARGS_TO_FORTRAN_ARGS((`RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(mtype,mop,`?tri?',`ARGS',0,`lang_c')'))dnl
RSB_M4_C_ARGS_TO_FORTRAN_ARGS((`RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,mop,`?tri?',`ARGS',0,`lang_c')'))dnl
)',`dnl
          patsubst(patsubst(d,` *:: \(trans\|order\)',`INTEGER :: \1'),`^ *::',RSB_M4_C2F_TYPE(mtype) ::)
')`'dnl
dnl
ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`dnl
          INTEGER,INTENT(OUT) :: A
')`'dnl
dnl
dnl RSB_M4_C_ARGS_TO_FORTRAN_ARGS(RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(mtype,mop,`?tri?',`ARGS'))`'dnl
dnl RSB_M4_C_ARGS_TO_FORTRAN_ARGS(RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,mop,`?tri?',`ARGS'))`'dnl
dnl
dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
define(`RSB_M4_SBLAS_SUBROUTINE_ARGS',`dnl
pushdef(`mop',$1)`'dnl
pushdef(`mtype',$2)`'dnl
pushdef(`lang',$3)`'dnl
dnl RSB_M4_PSBLAS_SUBROUTINE_ARGS($@)`'dnl
dnl
dnl
(RSB_M4_ARGS_TO_ACTUAL_ARGS((dnl
RSB_M4_SPBLAS_FUNCTION(mtype,mop,`',`ARGS',0,lang)`'dnl
dnl RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_SPBLAS_MATRIX_CREATION_FUNCS(mtype,mop,`',`ARGS',0,lang)))`'dnl
dnl RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_SPBLAS_MATRIX_ALL_L2_FUNCTION(mtype,mop,`?tri?',`ARGS',0,lang)))`'dnl
dnl RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,mop,`?tri?',`ARGS',0,lang)))`'dnl
dnl ifelse(RSB_M4_MEMBER(mop,RSB_M4_SPBLAS_MATRIX_BEGIN_MOPS),`1',`,A')`'dnl
)))`'
dnl
popdef(`lang')dnl
popdef(`mtype')dnl
popdef(`mop')dnl
')dnl
dnl
define(`RSB_M4_PSBLAS_SUBROUTINE_ARGS',`dnl
pushdef(`pmop',$1)dnl
pushdef(`mtype',$2)dnl
dnl
(a`'dnl
ifelse(pmop,`getdiag',`,d')`'dnl
ifelse(pmop,`set_elements',`,val,ia,ja,nnz,imin,imax,jmin,jmax,has_gtl,gtl')`'dnl
dnl ifelse(pmop,`set_elements',`,val,ia,ja,nnz,imin,imax,jmin,jmax,has_gtl,gtl,do_rebuild')`'dnl
ifelse(pmop,`set_element',`,v,i,j')`'dnl
ifelse(pmop,`get_element',`,v,i,j')`'dnl
ifelse(pmop,`scale',`,d')`'dnl
dnl ifelse(pmop,`get_rows_sparse',`,val,fr,lr,ia,ja,nnz,nzin,append')`'dnl
ifelse(pmop,`get_rows_sparse',`,val,fr,lr,ia,ja,nnz,nzin,append,has_iren,iren')`'dnl
ifelse(pmop,`get_rows_nnz',`,fr,lr,nnz')`'dnl
ifelse(pmop,`allocate_sparse_matrix',`,val,ia,ja,nnz,m,k,br,bc,extra')`'dnl
ifelse(pmop,`destroy_sparse_matrix',`')`'dnl
ifelse(pmop,`get_matrix_nnz',`,ires')`'dnl
ifelse(pmop,`infinity_norm',`,real_in,trans')`'dnl
ifelse(pmop,`ussm',`,b,c,lb,lc,nrhs,alpha,trans')`'dnl
ifelse(pmop,`usmm',`,b,c,lb,lc,nrhs,alpha,beta,trans')`'dnl
ifelse(pmop,`usmv',`,b,c,alpha,beta,trans')`'dnl
,info)`'dnl
dnl
popdef(`mtype')dnl
popdef(`pmop')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_INTERFACE_LIST',`dnl
pushdef(`firstarg',`0')dnl
foreach(`arg',`($@)',`ifelse(firstarg,`0',`pushdef(`firstarg',1)arg &',`
        &, arg &')`'')`'
dnl	DO NOT REMOVE THE FOLLOWING LINE
        & ;
ifelse(firstarg,`1',`popdef(`firstarg')')dnl
popdef(`firstarg')dnl
')dnl
dnl
dnl
dnl
dnl
