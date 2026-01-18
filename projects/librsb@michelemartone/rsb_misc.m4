dnl
dnl	@author: Michele Martone
dnl
dnl	forloop(TOKEN,LOWERI,UPPERI,ACTION)
dnl	-----------------------------------
dnl	Expands every occurrence of TOKEN in ACTION to each numerical
dnl	value in the [LOWERI,UPPERI] interval.
dnl	Therefore, ACTION is expanded [UPPERI-LOWERI+1] times, with a
dnl	varying index value being substituted to TOKEN.
dnl
divert(`-1')dnl
# forloop(var, from, to, stmt) - simple version
dnl define(`forloop', `pushdef(`$1', `$2')_forloop($@)popdef(`$1')')dnl
define(`forloop',`pushdef(`$1',`$2')_forloop($@)popdef(`$1')')dnl
dnl define(`_forloop',dnl
dnl        `$4`'ifelse($1, `$3', `', `define(`$1', incr($1))$0($@)')')dnl
define(`_forloop',dnl
`$4`'ifelse($1,`$3',`',`define(`$1',incr($1))$0($@)')')dnl
divert`'dnl
dnl	this `foreach' macro is the one in the ./examples  directory distributed with the M4 package
dnl
dnl
dnl
include(`rsb_config.m4')dnl		we include essential directives there
dnl
define(`RSB_M4_HEADER_MESSAGE',`
include(`rsb_license_header.inc')dnl
/*
 The code in this file was generated automatically by an M4 script. 
 It is not meant to be used as an API (Application Programming Interface).
 p.s.: right now, only row major matrix access is considered.

dnl File generated syscmd(`date')
 */
')dnl
dnl
dnl
dnl
dnl
dnl	foreach(TOKEN,VALUES,ACTION)
dnl	----------------------------
dnl	Expands each occurrence of TOKEN in ACTION to each value in the
dnl	VALUES list.
dnl	Therefore, ACTION is expanded a number of times equal to the
dnl	number of items in the VALUES list.
dnl
divert(`-1')dnl
# foreach(x, (item_1, item_2, ..., item_n), stmt)
#   parenthesized list, simple version
define(`foreach',`pushdef(`$1')_foreach($@)popdef(`$1')')dnl
define(`_arg1',`$1')dnl
define(`_foreach',`ifelse(`$2',`()',`',dnl
`define(`$1',_arg1$2)$3`'$0(`$1',(shift$2),`$3')')')dnl
divert`'dnl
dnl
dnl 
dnl	tolowercase(TOKEN)
dnl	------------------
dnl	Expands TOKEN to lowercase.
dnl 
define(`tolowercase',`translit($1,`A-Z',`a-z')')dnl
dnl 
dnl 
dnl	touppercase(TOKEN)
dnl	------------------
dnl	Expands TOKEN to uppercase.
dnl 
define(`touppercase',`translit($1,`a-z',`A-Z')')dnl
dnl 
dnl 
dnl	touppercase_(TOKEN)
dnl	------------------
dnl	Expands TOKEN to uppercase.
dnl
changequote([,])dnl
define([touppercase_],dnl
[dnl comment
translit([$1], [a-z], [A-Z])
])dnl
dnl 
dnl
dnl	singlequote(TOKEN)
dnl	------------------
dnl	Expands to TOKEN surrounded by single quotes.
dnl	That is, 'TOKEN'.
define([singlequote],
[dnl comment
'$1'])dnl
dnl 
dnl 
dnl	tolowercase_(TOKEN)
dnl	------------------
dnl	Expands TOKEN to lowercase
dnl	FIXME : remove this
dnl
define([tolowercase_],
[dnl comment
translit([$1], [A-Z], [a-z])
])dnl
dnl 
dnl 
changequote(`,')dnl
dnl 
dnl	RSB_M4_CHOPTRAILINGSPACES(STRING)
dnl	-------------------------
dnl	FIXME : document
dnl
define(`RSB_M4_CHOPTRAILINGSPACES',`dnl
pushdef(`type',$1)dnl
patsubst($1,` *$',)`'dnl
popdef(`type')dnl
')dnl
dnl 
dnl	RSB_M4_CHOPSPACES(STRING)
dnl	-------------------------
dnl	Expands to the input STRING with underscores ('_') substituted
dnl	to spaces (' ').
dnl
define(`RSB_M4_CHOPSPACES',`dnl
pushdef(`type',$1)dnl
patsubst($1,` ',_)`'dnl
popdef(`type')dnl
')dnl
dnl
dnl	
dnl	RSB_M4_TRANSPOSITION_CODE(transposition)
dnl	----------------------
dnl	DOCUMENT ME
dnl
define(`RSB_M4_TRANSPOSITION_CODE',`dnl
pushdef(`transposition',$1)dnl
`'touppercase(transposition)`'dnl
popdef(`transposition')dnl
')dnl
dnl
dnl	
dnl	RSB_M4_TYPE_CODE(TYPE)
dnl	----------------------
dnl	Expands to a code assigned to the specified numerical type.
dnl
define(`RSB_M4_TYPE_CODE',`dnl
pushdef(`type',$1)dnl
RSB_M4_CHOPSPACES(type)dnl
popdef(`type')dnl
')dnl
dnl 
dnl	RSB_M4_HAVE_TYPE(TYPE)
dnl	----------------------
dnl	...
dnl
define(`RSB_M4_HAVE_TYPE',`dnl
pushdef(`type',$1)dnl
RSB_M4_MEMBER(mtype,WANT_TYPES)dnl
popdef(`type')dnl
')dnl
dnl 
dnl
dnl
dnl	RSB_M4_HAVE_TYPE_PREPROCESSOR_SYMBOL(TYPE)
dnl	------------------------------------------
dnl	Converts a matrix type code in a preprocessor symbol used to
dnl	indicate the type availability.
dnl
dnl #define RSB_M4_HAVE_TYPE_PREPROCESSOR_SYMBOL(mop)
dnl
define(`RSB_M4_HAVE_TYPE_PREPROCESSOR_SYMBOL',`dnl
pushdef(`type',$1)`'dnl
`RSB_HAVE_TYPE_'touppercase( RSB_M4_CHOPSPACES(type) )dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_TYPE_INDEX_PREPROCESSOR_SYMBOL(type)
dnl	--------------------------------------------
dnl	Converts a matrix type code in a preprocessor symbol used to
dnl	index various mop-related arrays.
dnl
define(`RSB_M4_TYPE_INDEX_PREPROCESSOR_SYMBOL',`dnl
pushdef(`type',$1)`'dnl
`RSB_TYPE_INDEX_'touppercase( RSB_M4_CHOPSPACES(type) )dnl
popdef(`type')`'dnl
')dnl
dnl
dnl 
dnl
dnl	RSB_M4_OPTYPE_INDEX_PREPROCESSOR_SYMBOL(MOP)
dnl	--------------------------------------------
dnl	Converts a matrix operation code in a preprocessor symbol used to
dnl	index various mop-related arrays.
dnl
define(`RSB_M4_OPTYPE_INDEX_PREPROCESSOR_SYMBOL',`dnl
pushdef(`mop',$1)`'dnl
`RSB_OPTYPE_INDEX_'touppercase( RSB_M4_CHOPSPACES(mop) )dnl
popdef(`mop')`'dnl
')dnl
dnl
dnl 
dnl
dnl	RSB_M4_HAVE_OPTYPE_PREPROCESSOR_SYMBOL(MOP)
dnl	-------------------------------------------
dnl	Converts a matrix operation code in a preprocessor symbol used to
dnl	indicate the operation availability.
dnl
dnl #define RSB_M4_HAVE_OPTYPE_PREPROCESSOR_SYMBOL(mop)
dnl
define(`RSB_M4_HAVE_OPTYPE_PREPROCESSOR_SYMBOL',`dnl
pushdef(`mop',$1)`'dnl
`RSB_HAVE_OPTYPE_'touppercase( RSB_M4_CHOPSPACES(mop) )dnl
popdef(`mop')`'dnl
')dnl
dnl
dnl
dnl	RSB_M4_DEBUGINFO(``MACRO_SYMBOL'')dnl
dnl	-------------------------------------
dnl	Will expand to a C comment stating debug info about the given macro.
dnl
dnl	e.g.:
dnl	RSB_M4_DEBUGINFO(``RSB_M4_UNROLL_KERNEL'')dnl
dnl
dnl
define(`RSB_M4_DEBUGINFO',dnl
ifelse(`RSB_M4_DEBUG',`1',`/* generated by the $1 macro */',`')dnl
)dnl
dnl
dnl
dnl
dnl	RSB_M4_SPACED_LIST
dnl	-----------------
dnl	Expands the given list inserting spaces between each consecutive element.
dnl
define(`RSB_M4_SPACED_LIST',`dnl
patsubst(`'dnl
patsubst(`'dnl
foreach(`listel',$1,listel )`'dnl
,` $',`')`'dnl
,` ',` ')`'dnl
')dnl
dnl
dnl
dnl	RSB_M4_COMMA_LIST
dnl	-----------------
dnl	Expands the given list inserting commas between each consecutive element.
dnl
define(`RSB_M4_COMMA_LIST',`dnl
patsubst(`'dnl
patsubst(`'dnl
foreach(`listel',$1,listel )`'dnl
,` $',`')`'dnl
,` ',`,')`'dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_QUOTED_COMMA_LIST
dnl	------------------------
dnl
define(`RSB_M4_QUOTED_COMMA_LIST',`dnl
dnl
patsubst(`'dnl
patsubst(`'dnl
foreach(`listel',$1,"listel" )`'dnl
,` $',`')`'dnl
,` ',`,')`'dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_ARG_TO_ACTUAL_ARG',`patsubst($1,`.*\(\<[a-zA-Z_0-9]+$\)',`\1')')dnl
dnl
dnl	RSB_M4_ARGS_TO_ACTUAL_ARGS
dnl	--------------------------
dnl	Takes a C prototype string and cleans it from the type and array declarators,
dnl	making from it an argument list.
dnl	e.g.:
dnl	RSB_M4_ARGS_TO_ACTUAL_ARGS(`const struct rsb_mtx_t * m, const struct rsb_options_t *o, const void * rhs, void * out')
dnl	=> m, o, rhs, out
dnl
dnl was:
dnl `'patsubst( patsubst( foreach(`arg',$1,`patsubst(patsubst(arg`',`.+[ *]',` '),`\[.*\]',`')')`', `^ ', `'),` ',`,')`'dnl
dnl
define(`RSB_M4_ARGS_TO_ACTUAL_ARGS',`dnl
dnl
dnl	WARNING : this is THIN ICE :)
pushdef(`firstarg',`0')dnl
foreach(`arg',`$1',`ifelse(firstarg,`0',`pushdef(`firstarg',1)',`,')`'RSB_M4_ARG_TO_ACTUAL_ARG(arg)')`'dnl
ifelse(firstarg,`1',`popdef(`firstarg')')dnl
popdef(`firstarg')dnl
')dnl
dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP
dnl	-----------------------------------
dnl	Takes a string extracted from applying RSB_M4_ARGS_TO_ACTUAL_ARGS, and then
dnl	prepending patterns looking as rsb_mtx_t members with a "mtxAp->".
dnl
dnl	e.g.:
dnl	RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP(`const struct rsb_mtx_t * m, const struct rsb_options_t *o, const void * rhs, void * out, const int rpntr')
dnl	=> m, o, rhs, out, mtxAp->rpntr
dnl
define(`RSB_M4_ACTUAL_ARGS_APPLY_MEMBERSHIP',`dnl
dnl
dnl	WARNING : this is THIN ICE :)
dnl	patsubst(`$@',`\<flags\>\|\<rpntr\>\|\<cpntr\>\|\<bindx\>\|\<bpntr\>\|\<VA\>\|\<indptr\>\|\<Mdim\>\|\<mdim\>\|\<br\>\|\<bc\>\|\<roff\>\|\<coff\>\|\<nnz\>',`mtxAp->\&')`'dnl was \0, but gnu m4 told me to use \&
ifelse(RSB_M4_WANT_20110206_BOUNDED_BOX_PATCH,1,`dnl
	patsubst(`patsubst(`patsubst(`patsubst(`$@',`\<flags\>\|\<rpntr\>\|\<cpntr\>\|\<bpntr\>\|\<VA\>\|\<indptr\>\|\<Mdim\>\|\<mdim\>\|\<br\>\|\<bc\>\|\<roff\>\|\<coff\>\|\<nnz\>',`mtxAp->\&')',`\<bindx\>',(citype*)`mtxAp->\&')',`\<br\>',`broff-mtxAp->roff')',`\<bc\>',`bm')`'dnl was \0, but gnu m4 told me to use \&
',`dnl
	patsubst(`patsubst(`$@',`\<flags\>\|\<rpntr\>\|\<cpntr\>\|\<bpntr\>\|\<VA\>\|\<indptr\>\|\<Mdim\>\|\<mdim\>\|\<br\>\|\<bc\>\|\<roff\>\|\<coff\>\|\<nnz\>',`mtxAp->\&')',`\<bindx\>',(citype*)`mtxAp->\&')`'dnl was \0, but gnu m4 told me to use \&
')dnl
')dnl
dnl
dnl
dnl	RSB_M4_LIST_LENGTH(LIST)
dnl	------------------------
dnl	This macro returns the input list length.
dnl
define(`RSB_M4_LIST_LENGTH',`$#')dnl
dnl
dnl
dnl	The following values specify the row and column unroll factors
dnl	to be chosen when a specialized submatrix size function is not found.
dnl
dnl	FIXME : these default unrollings should be in!
define(`RSB_M4_ITH_LIST_ELEMENT',`ifelse($#,1,$2,`ifelse($1,`1',$2,`RSB_M4_ITH_LIST_ELEMENT(decr($1),shift(shift($@)))')')')dnl
dnl
define(`RSB_M4_MIDDLE_LIST_ELEMENT',`RSB_M4_ITH_LIST_ELEMENT(eval(($#+1)/2),$@)')dnl
dnl
define(`RSB_M4_LAST_LIST_ELEMENT',`RSB_M4_ITH_LIST_ELEMENT($#,$@)')dnl
dnl
define(`RSB_M4_FIRST_LIST_ELEMENT',`RSB_M4_ITH_LIST_ELEMENT(1,$@)')dnl
dnl
dnl	RSB_M4_EMPTY_LIST
dnl	-----------------
dnl	Gives 1 if the list is empty, otherwise 0
dnl
define(`RSB_M4_EMPTY_LIST',`ifelse($#,1,0)')dnl
dnl
dnl
dnl
dnl	RSB_M4_FIRST
dnl	------------
dnl	This macro returns the first input argument.
define(`RSB_M4_FIRST',`$1')dnl
dnl
dnl
dnl
dnl	RSB_M4_MAX2
dnl	-----------
dnl	This macro returns the maximum input argument among 2 arguments.
dnl
define(`RSB_M4_MAX2',`ifelse(eval(`$1>$2'),1,`$1',`$2')')dnl
dnl
dnl
dnl
dnl	RSB_M4_MIN2
dnl	-----------
dnl	This macro returns the maximum input argument among 2 arguments.
dnl
define(`RSB_M4_MIN2',`ifelse($2,,$1,`ifelse(eval(`$1<$2'),1,`$1',`$2')')')dnl
dnl
dnl
dnl
dnl	RSB_M4_MAXN
dnl	-----------
dnl	This macro returns the maximum input argument.
dnl
define(`RSB_M4_MAXN',`ifelse($#,1,$1,`ifelse($#,2,`RSB_M4_MAX2($1,$2)',`RSB_M4_MAX2($1,RSB_M4_MAXN(shift($@)))')')')dnl
dnl
dnl
dnl
dnl	RSB_M4_MINN
dnl	-----------
dnl	This macro returns the maximum input argument.
dnl
define(`RSB_M4_MINN',`ifelse($2,,$1,`ifelse($#,2,`RSB_M4_MIN2($1,$2)',`RSB_M4_MIN2($1,RSB_M4_MINN(shift($@)))')')')dnl
dnl
dnl
dnl
dnl	RSB_M4_EXCEPT
dnl	-------------
dnl	This macro processes the input arguments ($@) removing all occurrences of the first element ($1).
dnl
define(`RSB_M4_EXCEPT',`pushdef(`GOT',`0')ifelse($#,0,,`foreach(`exel',($@),`ifelse(exel,$1,,`ifelse(GOT,`1',`,')`'exel`'ifelse(GOT,`0',`popdef(`GOT')pushdef(`GOT',`1')')')')')popdef(`GOT')')dnl
dnl
dnl
dnl
dnl	RSB_M4_SORT
dnl	-----------
dnl	This macro sorts the input arguments.
dnl
define(`RSB_M4_SORT',`ifelse($#,1,$1,`pushdef(`ALL',`$@')pushdef(`REST',RSB_M4_MINN(ALL))REST,RSB_M4_SORT(RSB_M4_EXCEPT(REST,ALL))popdef(`ALL')popdef(`REST')')')dnl
dnl
dnl
dnl
dnl	RSB_M4_SAME
dnl	-----------
dnl
define(`RSB_M4_SAME',`ifelse($1,$2,`1',`0')')dnl
dnl
dnl
dnl	RSB_M4_OR
dnl	---------
dnl
define(`RSB_M4_OR',`pushdef(`GOT',`0')foreach(`exel',($@),`ifelse(exel,`1',`ifelse(GOT,`0',`popdef(`GOT')pushdef(`GOT',`1')')')')GOT`'popdef(`GOT')dnl
')dnl
dnl
dnl
dnl	RSB_M4_AND
dnl	----------
dnl
define(`RSB_M4_AND',`pushdef(`GOT',`1')foreach(`exel',($@),`ifelse(exel,`0',`ifelse(GOT,`1',`popdef(`GOT')pushdef(`GOT',`0')')')')GOT`'popdef(`GOT')dnl
')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_XOR
dnl	----------
dnl
define(`RSB_M4_XOR',`RSB_M4_NOT(RSB_M4_OR(RSB_M4_AND($1,$2),RSB_M4_AND(RSB_M4_NOT($1),RSB_M4_NOT($2))))')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_IMPLY
dnl	------------
dnl
define(`RSB_M4_IMPLY',`RSB_M4_OR(RSB_M4_NOT($1),$2)')dnl
dnl
dnl
dnl
dnl	RSB_M4_NOT
dnl	----------
dnl
define(`RSB_M4_NOT',`ifelse($1,`0',`1',`0')`'foreach(`exel',(shift($@)),`ifelse(exel,`0',`,1',`,0')')')dnl
dnl
dnl
dnl
dnl	RSB_M4_MEMBER
dnl	-------------
dnl	if $1 is among (shift($@)), returns 1, otherwise 0
dnl
define(`RSB_M4_MEMBER',`pushdef(`GOT',`0')foreach(`exel',(shift($@)),`ifelse(exel,$1,`ifelse(GOT,`0',`popdef(`GOT')pushdef(`GOT',`1')')')')GOT`'popdef(`GOT')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_ALL_LIST_ARGS
dnl	--------------------
dnl	All macro arguments : used as a syntactical trick.
define(`RSB_M4_ALL_LIST_ARGS',`$@')dnl
dnl
dnl
dnl
dnl
dnl	RSB_M4_LIST_PUSH_BACK
dnl	---------------------
dnl	Will give back the second argument list and the second argument.
dnl	Will not work with an empty list.
dnl
define(`RSB_M4_LIST_PUSH_BACK',`RSB_M4_ALL_LIST_ARGS$1,$2')dnl
dnl
dnl
dnl
dnl
dnl	ROW AND COLUMN UNROLL FACTORS
dnl	-----------------------------
dnl	Unroll count / type parameterizations
dnl
define(`RSB_M4_ROWS_UNROLL',(RSB_M4_SORT(WANT_ROW_UNLOOP_FACTORS)))dnl
define(`RSB_M4_COLUMNS_UNROLL',(RSB_M4_SORT(WANT_COLUMN_UNLOOP_FACTORS)))dnl
define(`RSB_M4_MATRIX_TYPES',(WANT_TYPES))dnl
define(`RSB_M4_ALL_MATRIX_TYPES',(WANT_MATRIX_ALL_TYPES))dnl
dnl
define(`RSB_M4_INVALID_TYPE',`invalid_type')dnl
define(`RSB_M4_DEFAULT_TYPE',RSB_M4_FIRST(WANT_TYPES))dnl
define(`RSB_M4_DEFAULT_POSSIBLY_INTEGER_TYPE',ifelse(RSB_M4_MEMBER(`int',WANT_TYPES),`1',`int',RSB_M4_FIRST(WANT_TYPES)))dnl
dnl
define(`RSB_M4_DEFAULT_SYMMETRY',RSB_M4_FIRST(RSB_M4_WANT_MATRIX_SYMMETRY))dnl
define(`RSB_M4_DEFAULT_TRANSPOSITION',RSB_M4_FIRST(RSB_M4_WANT_MATRIX_TRANSPOSITIONS))dnl
dnl
define(`RSB_M4_MATRIX_OPS',(WANT_MATRIX_OPS))dnl
dnl
define(`RSB_M4_TRANS_N',`n')dnl
define(`RSB_M4_TRANS_T',`t')dnl
define(`RSB_M4_TRANS_C',`c')dnl
define(`RSB_M4_TRANS_INVALID',`?')dnl
dnl
define(`RSB_M4_WANT_MATRIX_TRANSPOSITIONS',`RSB_M4_TRANS_N,RSB_M4_TRANS_T,RSB_M4_TRANS_C')dnl
dnl define(`RSB_M4_WANT_MATRIX_TRANSPOSITIONS',``n',`t',`h'')dnl
dnl define(`RSB_M4_WANT_MATRIX_TRANSPOSITIONS',``n',`t'')dnl
dnl define(`RSB_M4_WANT_MATRIX_TRANSPOSITIONS',``n'')dnl
define(`RSB_M4_SYMBOL_SYMMETRIC',`S')dnl
define(`RSB_M4_SYMBOL_UNSYMMETRIC',`U')dnl
define(`RSB_M4_SYMBOL_HERMITIAN',`H')dnl
define(`RSB_M4_WANT_MATRIX_SYMMETRY',`RSB_M4_SYMBOL_UNSYMMETRIC,RSB_M4_SYMBOL_SYMMETRIC,RSB_M4_SYMBOL_HERMITIAN')dnl
dnl define(`RSB_M4_WANT_MATRIX_SYMMETRY',`RSB_M4_SYMBOL_UNSYMMETRIC')dnl
dnl define(`RSB_M4_WANT_MATRIX_SYMMETRY',`RSB_M4_SYMBOL_SYMMETRIC')dnl
dnl
define(`RSB_M4_IS_NOT_UNSYMMETRIC',`pushdef(`k_symmetry',$1)RSB_M4_NOT(RSB_M4_SAME(k_symmetry,RSB_M4_SYMBOL_UNSYMMETRIC))popdef(`k_symmetry')')dnl
define(`RSB_M4_IS_UNSYMMETRIC',`pushdef(`k_symmetry',$1)RSB_M4_SAME(k_symmetry,RSB_M4_SYMBOL_UNSYMMETRIC)popdef(`k_symmetry')')dnl
define(`RSB_M4_IS_SYMMETRIC',`pushdef(`k_symmetry',$1)RSB_M4_SAME(k_symmetry,RSB_M4_SYMBOL_SYMMETRIC)popdef(`k_symmetry')')dnl
dnl
dnl
define(`RSB_M4_MATRIX_TRANSPOSITIONS',(RSB_M4_WANT_MATRIX_TRANSPOSITIONS))dnl
dnl
define(`RSB_M4_MATRIX_UPLO_TYPES',(`u',`l',`g'))dnl
dnl define(`RSB_M4_MATRIX_UPLO_TYPES',(`g'))dnl
dnl
dnl
define(`RSB_M4_MATRIX_DIAGONAL_DENOMINATION',``diagonal 'ifelse(RSB_M4_IS_DIAGONAL_IMPLICIT(k_diagonal),1,`implicit',`explicit')')dnl
define(`RSB_M4_MATRIX_DIAGONAL_TYPES',(`e',`i'))dnl
define(`RSB_M4_DEFAULT_DIAGONAL_TYPE',`e')dnl
define(`RSB_M4_IS_DIAGONAL_IMPLICIT',`pushdef(`k_diagonal',$1)RSB_M4_SAME(k_diagonal,`i')popdef(`k_diagonal')')dnl
dnl
dnl define(`RSB_M4_MATRIX_COORDINATE_TYPES',(`rsb_coo_idx_t',`rsb_half_idx_t'))dnl	FIXME : new
dnl define(`RSB_M4_MATRIX_COORDINATE_TYPES',(`rsb_half_idx_t'))dnl	FIXME : new
define(`RSB_M4_MATRIX_COORDINATE_TYPES',(`rsb_coo_idx_t'ifelse(WANT_HALFWORD_INDICES,`yes',`,rsb_half_idx_t',`')))dnl	FIXME : new
define(`RSB_M4_WANT_SPSM_DIAG_CHECK',ifelse(WANT_SPSM_DIAG_CHECK,`yes',`1',`0'))dnl	FIXME : new
dnl define(`RSB_M4_WANT_SPSM_DIAG_CHECK',1)dnl	FIXME : new
define(RSB_M4_TRANSPOSITION_OP_EFFECT,`dnl
dnl
pushdef(`transposition',$1)`'dnl
pushdef(`operand',$2)`'dnl
dnl
`{'dnl
ifelse(transposition,RSB_M4_TRANS_T,`'operand`^T',`dnl
ifelse(transposition,RSB_M4_TRANS_C,`'operand`^H',`dnl
ifelse(transposition,RSB_M4_TRANS_N,`'operand`',dnl
`'operand`')')')`'dnl
`}'dnl
dnl
dnl
popdef(`operand')`'dnl
popdef(`transposition')`'dnl
')dnl
dnl
define(`RSB_M4_THRESHOLD_VALUE',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(type,`long double complex',`1e-6')`'dnl
ifelse(type,`double complex',`1e-6')`'dnl
ifelse(type,`float complex',`1e-4')`'dnl
ifelse(type,`long double',`1e-6')`'dnl
ifelse(type,`double',`1e-6')`'dnl
ifelse(type,`float',`1e-5')`'dnl
ifelse(type,`int',`0')`'dnl
dnl FIXME : and, for other, non canonical types ? FIXME
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_REALT',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(type,`long double complex',`long double')`'dnl
ifelse(type,`double complex',`double')`'dnl
ifelse(type,`float complex',`float')`'dnl
ifelse(type,`long double',`long double')`'dnl
ifelse(type,`double',`double')`'dnl
ifelse(type,`float',`float')`'dnl
ifelse(type,`int',`int')`'dnl
dnl FIXME : and, for other, non canonical types ? FIXME 
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_CREAL',`dnl
pushdef(`type',$1)`'dnl
pushdef(`value',$2)`'dnl
dnl
ifelse(type,`long double complex',`creall(value)')`'dnl
ifelse(type,`double complex',`creal(value)')`'dnl
ifelse(type,`float complex',`crealf(value)')`'dnl
ifelse(type,`long double',`(value)')`'dnl
ifelse(type,`double',`(value)')`'dnl
ifelse(type,`float',`(value)')`'dnl
ifelse(type,`int',`(value)')`'dnl
dnl FIXME : and, for other, non canonical types ? FIXME 
dnl
popdef(`value')`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_CIMAG',`dnl
pushdef(`type',$1)`'dnl
pushdef(`value',$2)`'dnl
dnl
ifelse(type,`long double complex',`cimagl(value)')`'dnl
ifelse(type,`double complex',`cimag(value)')`'dnl
ifelse(type,`float complex',`cimagf(value)')`'dnl
ifelse(type,`long double',`(value)')`'dnl
ifelse(type,`double',`(value)')`'dnl
ifelse(type,`float',`(value)')`'dnl
ifelse(type,`int',`(value)')`'dnl
dnl FIXME : and, for other, non canonical types ? FIXME 
dnl
popdef(`value')`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_BACKUP_SIZEOF',`dnl
pushdef(`type',$1)`'dnl
dnl
ifelse(type,`double complex',`16')`'dnl
ifelse(type,`float complex',`8')`'dnl
ifelse(type,`double',`8')`'dnl
ifelse(type,`float',`4')`'dnl
ifelse(type,`int',`4')`'dnl
dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_ASSIGN',`dnl
pushdef(`dtype',$1)`'dnl
pushdef(`stype',$2)`'dnl
pushdef(`lval',$3)`'dnl
pushdef(`rval',$4)`'dnl
dnl
dnl
ifelse(RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(dtype),RSB_M4_NOT(RSB_M4_IS_COMPLEX_TYPE(stype))),1,`lval = rval + 0*I;',`')dnl
ifelse(RSB_M4_AND(RSB_M4_NOT(RSB_M4_IS_COMPLEX_TYPE(dtype)),RSB_M4_IS_COMPLEX_TYPE(stype)),1,`lval = RSB_M4_CREAL(stype,rval);',`')dnl
ifelse(RSB_M4_NOT(RSB_M4_XOR(RSB_M4_IS_COMPLEX_TYPE(dtype),RSB_M4_IS_COMPLEX_TYPE(stype))),1,`lval = rval;',`')dnl
dnl
dnl FIXME : and, for other, non canonical types ? FIXME 
dnl
popdef(`rval')`'dnl
popdef(`lval')`'dnl
popdef(`stype')`'dnl
popdef(`dtype')`'dnl
dnl
')dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_ABS',`dnl
pushdef(`type',$1)`'dnl
pushdef(`value',$2)`'dnl
dnl
ifelse(type,`long double complex',`cabsl(value)')`'dnl
ifelse(type,`double complex',`cabs(value)')`'dnl
ifelse(type,`float complex',`cabsf(value)')`'dnl
ifelse(type,`long double',`fabsl(value)')`'dnl
ifelse(type,`double',`fabs(value)')`'dnl
ifelse(type,`float',`fabsf(value)')`'dnl
ifelse(type,`int',`abs(value)')`'dnl
dnl FIXME : and, for other, non canonical types ? FIXME 
dnl
popdef(`value')`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_ABS_IF_1',`dnl
pushdef(`type',$1)`'dnl
pushdef(`value',$2)`'dnl
pushdef(`condition',$3)`'dnl
dnl
ifelse(condition,`1',`RSB_M4_ABS(type,value)',`value')`'dnl
dnl
popdef(`condition')`'dnl
popdef(`value')`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_POW',`dnl
pushdef(`type',$1)`'dnl
pushdef(`value',$2)`'dnl
pushdef(`exp',$3)`'dnl
dnl
dnl	FIXME
dnl
ifelse(type,`long double complex',`cpowl(value,exp)')`'dnl
ifelse(type,`double complex',`cpow(value,exp)')`'dnl
ifelse(type,`float complex',`cpowf(value,exp)')`'dnl
dnl
ifelse(type,`long double',`powl(value,exp)')`'dnl
ifelse(type,`double',`pow(value,exp)')`'dnl
ifelse(type,`float',`powf(value,exp)')`'dnl
ifelse(type,`int',`(int)pow((int)(value),(int)(exp))')`'dnl	yeah, it is dumb.
dnl
popdef(`exp')`'dnl
popdef(`value')`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
dnl
dnl
define(`RSB_M4_SQRT',`dnl
pushdef(`type',$1)`'dnl
pushdef(`value',$2)`'dnl
dnl
dnl	FIXME
dnl
ifelse(type,`long double complex',`csqrtl(value)')`'dnl
ifelse(type,`double complex',`csqrt(value)')`'dnl
ifelse(type,`float complex',`csqrtf(value)')`'dnl
dnl
ifelse(type,`long double',`sqrtl(value)')`'dnl
ifelse(type,`double',`sqrt(value)')`'dnl
ifelse(type,`float',`sqrtf(value)')`'dnl
ifelse(type,`int',`(int)sqrt((int)(value))')`'dnl
dnl
popdef(`value')`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_ONE',`dnl
pushdef(`type',$1)`'dnl
`'((type)(1.0))`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_ZERO',`dnl
pushdef(`type',$1)`'dnl
`'((type)(0))`'dnl
popdef(`type')`'dnl
dnl
')dnl
dnl
dnl
define(`RSB_M4_CONJ_SYM',`dnl
pushdef(`type',$1)`'dnl
pushdef(`transposition',$2)`'dnl
pushdef(`k_symmetry',$3)`'dnl
dnl
ifelse(dnl
RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(type),dnl
RSB_M4_XOR(RSB_M4_SAME(transposition,RSB_M4_TRANS_C),RSB_M4_SAME(k_symmetry,RSB_M4_SYMBOL_HERMITIAN))),1,`dnl
ifelse(type,`long double complex',`conjl')`'dnl FIXME : long double complex is not supported
ifelse(type,`double complex',`conj')`'dnl
ifelse(type,`float complex',`conjf')`'dnl
',`dnl
`'dnl
')`'dnl
dnl
popdef(`k_symmetry')`'dnl
popdef(`transposition')`'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_CONJ',`dnl
pushdef(`exp',$1)`'dnl
pushdef(`type',$2)`'dnl
pushdef(`transposition',$3)`'dnl
pushdef(`k_symmetry',$4)`'dnl
dnl
ifelse(dnl
RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(type),dnl
RSB_M4_XOR(RSB_M4_SAME(transposition,RSB_M4_TRANS_C),RSB_M4_SAME(k_symmetry,RSB_M4_SYMBOL_HERMITIAN))),1,`dnl
ifelse(type,`long double complex',`conjl(exp)')`'dnl FIXME : long double complex is not supported
ifelse(type,`double complex',`conj(exp)')`'dnl
ifelse(type,`float complex',`conjf(exp)')`'dnl
',`dnl
exp`'dnl
')`'dnl
dnl
popdef(`k_symmetry')`'dnl
popdef(`transposition')`'dnl
popdef(`type')`'dnl
popdef(`exp')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_MATRIX_TYPE_AS_CXX',`dnl
pushdef(`type',$1)`'dnl
ifelse(type,`long double complex',`std::complex<long double>')`'dnl
ifelse(type,`double complex',`std::complex<double>')`'dnl
ifelse(type,`float complex',`std::complex<float>')`'dnl
ifelse(type,`float',`type')`'dnl
ifelse(type,`long double',`type')`'dnl
ifelse(type,`double',`type')`'dnl
ifelse(type,`int',`type')`'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_MATRIX_TYPES_LIST_CXX',`dnl
pushdef(`lastel',RSB_M4_LAST_LIST_ELEMENT(WANT_TYPES))`'dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`RSB_M4_MATRIX_TYPE_AS_CXX(type)'`ifelse(type,lastel,`',`,')')dnl
popdef(`lastel')`'dnl
dnl RSB_M4_MATRIX_TYPE_AS_CXX(RSB_M4_FIRST_LIST_ELEMENT(WANT_TYPES))
dnl
')dnl
dnl
dnl
define(`RSB_M4_H2T_TRANSPOSITION',`dnl
pushdef(`transposition',$1)`'dnl
dnl
ifelse(transposition,RSB_M4_TRANS_T,RSB_M4_TRANS_C)`'dnl
ifelse(transposition,RSB_M4_TRANS_C,RSB_M4_TRANS_T)`'dnl
ifelse(transposition,RSB_M4_TRANS_N,RSB_M4_TRANS_N)`'dnl
dnl
popdef(`transposition')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_TRANSPOSE_SYMMETRY',`dnl
pushdef(`k_symmetry',$1)`'dnl
dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_HERMITIAN,RSB_M4_SYMBOL_UNSYMMETRIC)`'dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_UNSYMMETRIC,RSB_M4_SYMBOL_HERMITIAN)`'dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_SYMMETRIC,RSB_M4_SYMBOL_HERMITIAN)`'dnl
dnl
popdef(`k_symmetry')`'dnl
')dnl
dnl
define(`RSB_M4_TRANSPOSE_TRANSPOSITION',`dnl
pushdef(`transposition',$1)`'dnl
dnl
ifelse(transposition,RSB_M4_TRANS_T,RSB_M4_TRANS_N)`'dnl
ifelse(transposition,RSB_M4_TRANS_C,RSB_M4_TRANS_N)`'dnl
ifelse(transposition,RSB_M4_TRANS_N,RSB_M4_TRANS_T)`'dnl
dnl
popdef(`transposition')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_MATRIX_SYMMETRY',(RSB_M4_WANT_MATRIX_SYMMETRY))dnl
dnl
dnl
define(RSB_M4_SYMMETRY_SWITCH,`dnl
dnl
pushdef(`k_symmetry',$1)`'dnl
dnl
dnl	FIXME
dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_UNSYMMETRIC,RSB_M4_SYMBOL_UNSYMMETRIC,dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_SYMMETRIC,RSB_M4_SYMBOL_UNSYMMETRIC,dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_HERMITIAN,RSB_M4_SYMBOL_UNSYMMETRIC,dnl
`')))`'dnl
dnl
popdef(`k_symmetry')`'dnl
')dnl
dnl
dnl
define(RSB_M4_SYMMETRY_EFFECT,`dnl
dnl
dnl
pushdef(`k_symmetry',$1)`'dnl
pushdef(`operand',$2)`'dnl
dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_UNSYMMETRIC,`'operand` \neq 'operand`^T',dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_SYMMETRIC,`'operand` == 'operand`^T',dnl
ifelse(k_symmetry,RSB_M4_SYMBOL_HERMITIAN,`'operand` == 'operand`^H',dnl
`')))`'dnl
dnl
dnl
popdef(`operand')`'dnl
popdef(`k_symmetry')`'dnl
')dnl
dnl
dnl
define(`RSB_M4_IS_COMPLEX_TYPE',`dnl
pushdef(`type',$1)dnl
RSB_M4_MEMBER(type,`double complex',`float complex',`long double complex')`'dnl
popdef(`type')dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_UIM_CONJ',`dnl
dnl Unconditional iteration [multiplication] conjugation
pushdef(`exp',$1)`'dnl
pushdef(`type',$2)`'dnl
pushdef(`transposition',$3)`'dnl
pushdef(`symmetry',$4)`'dnl
dnl
ifelse(RSB_M4_OR(dnl
`RSB_M4_AND(RSB_M4_SAME(symmetry,RSB_M4_SYMBOL_SYMMETRIC),RSB_M4_SAME(transposition,RSB_M4_TRANS_C))',dnl
`RSB_M4_AND(RSB_M4_SAME(symmetry,RSB_M4_SYMBOL_HERMITIAN),RSB_M4_SAME(transposition,RSB_M4_TRANS_C))'dnl
),1,RSB_M4_CONJ(exp,type,RSB_M4_TRANS_N,RSB_M4_SYMBOL_HERMITIAN),exp)dnl
dnl
popdef(`exp')`'dnl
popdef(`symmetry')`'dnl
popdef(`transposition')`'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_CIM_CONJ',`dnl
dnl Conditional iteration [multiplication] conjugation
pushdef(`exp',$1)`'dnl
pushdef(`type',$2)`'dnl
pushdef(`transposition',$3)`'dnl
pushdef(`symmetry',$4)`'dnl
dnl
ifelse(RSB_M4_OR(dnl
`RSB_M4_AND(RSB_M4_SAME(symmetry,RSB_M4_SYMBOL_SYMMETRIC),RSB_M4_SAME(transposition,RSB_M4_TRANS_C))',dnl
`RSB_M4_AND(RSB_M4_SAME(symmetry,RSB_M4_SYMBOL_HERMITIAN),RSB_M4_NOT(RSB_M4_SAME(transposition,RSB_M4_TRANS_C)))'dnl
),1,RSB_M4_CONJ(exp,type,RSB_M4_TRANS_N,RSB_M4_SYMBOL_HERMITIAN),exp)dnl
dnl
popdef(`exp')`'dnl
popdef(`symmetry')`'dnl
popdef(`transposition')`'dnl
popdef(`type')`'dnl
')dnl
dnl
dnl
dnl
define(`RSB_M4_MATRIX_ALL_OPS',(WANT_MATRIX_ALL_OPS))dnl
define(`RSB_M4_MATRIX_STORAGE',(WANT_MATRIX_STORAGE))dnl
define(`RSB_M4_LONG_IDX',WANT_LONG_IDX)dnl
define(`RSB_M4_WANT_IHI',WANT_IHI)dnl
define(`RSB_M4_ROWS_FALLBACK_UNROLL',RSB_M4_MIDDLE_LIST_ELEMENT(RSB_M4_SORT(WANT_ROW_UNLOOP_FACTORS)))dnl
define(`RSB_M4_COLUMNS_FALLBACK_UNROLL',RSB_M4_MIDDLE_LIST_ELEMENT(RSB_M4_SORT(WANT_COLUMN_UNLOOP_FACTORS)))dnl
dnl
define(`RSB_M4_BCOO_FORMATS',(WANT_MATRIX_BCOO_STORAGE))dnl
define(`RSB_M4_BCSS_FORMATS',(WANT_MATRIX_BCSS_STORAGE))dnl
define(`RSB_M4_VB_FORMATS',(WANT_MATRIX_VB_STORAGE))dnl
dnl
define(`RSB_M4_IS_IMPLEMENTED_MOP',`RSB_M4_MEMBER($1,WANT_MATRIX_OPS)')dnl
define(`RSB_M4_IS_FORMAT_BCSS',`RSB_M4_MEMBER($1,`BCSR',`BCSC')')dnl
define(`RSB_M4_IS_FORMAT_BCOO',`RSB_M4_MEMBER($1,`BCOR',`BCOC')')dnl
define(`RSB_M4_IS_FORMAT_BCXX',`RSB_M4_OR(RSB_M4_IS_FORMAT_BCSS($1),RSB_M4_IS_FORMAT_BCOO($1))')dnl
define(`RSB_M4_IS_FORMAT_ROW_MAJOR',`RSB_M4_MEMBER($1,`BCSR',`VBR',`LR',`BCOR')')dnl
define(`RSB_M4_IS_FORMAT_COLUMN_MAJOR',`RSB_M4_MEMBER($1,`BCSC',`VBC',`LC',`BCOC')')dnl
define(`RSB_M4_IS_FORMAT_LINKED_LIST',`RSB_M4_MEMBER($1,`LR',`LC')')dnl NEW
define(`RSB_M4_IS_KERNEL_REGISTER_UNBLOCKED',`RSB_M4_AND(RSB_M4_SAME($1,1),RSB_M4_SAME($2,1))')dnl use as RSB_M4_IS_KERNEL_REGISTER_UNBLOCKED(b_rows,b_columns)
define(`RSB_M4_IS_KERNEL_REGISTER_BLOCKED',`RSB_M4_NOT(RSB_M4_IS_KERNEL_REGISTER_UNBLOCKED($1,$2))')dnl use as RSB_M4_IS_KERNEL_REGISTER_BLOCKED(b_rows,b_columns)
dnl
dnl define(`RSB_M4_PREFIX',`rsb_')dnl
define(`RSB_M4_PREFIX',`rsb__')dnl
dnl
dnl	RSB_M4_MATRIX_META_OPS
dnl	----------------------
dnl
define(`RSB_M4_MATRIX_META_OPS',dnl
(RSB_M4_LIST_PUSH_BACK(dnl
RSB_M4_MATRIX_OPS`'dnl
dnl (RSB_M4_LIST_PUSH_BACK(RSB_M4_MATRIX_OPS,`sort_only'))dnl
,`mat_stats'))dnl
)dnl
dnl
dnl	RSB_M4_MATRIX_META_OPS_REDUCED
dnl	----------------------
dnl
define(`RSB_M4_MATRIX_META_OPS_REDUCED',dnl
(RSB_M4_LIST_PUSH_BACK(dnl
dnl (`spmv_uaua',`spsv_uxua')`'dnl
(`spmv_sxsa',`spsv_sxsx')`'dnl The most general kernels.
,`mat_stats'))dnl
)dnl
dnl
dnl
dnl	RSB_M4_MATRIX_OP_IS_META_OP(OP)
dnl	-------------------------------
dnl	TODO: should differentiate from RSB_M4_IS_IMPLEMENTED_MOP
dnl
define(`RSB_M4_MATRIX_OP_IS_META_OP',`dnl
pushdef(`mop',$1)`'dnl
ifelse(RSB_M4_MEMBER(mop,WANT_MATRIX_OPS),`1',`0',`1')dnl
popdef(`mop')`'dnl
')dnl
dnl
dnl	RSB_M4_INTERVAL_LIST(LOWER_INDEX,UPPER_INDEX,[INCREMENT])
dnl	---------------------------------------------------------
dnl
dnl
define(`RSB_M4_INTERVAL_LIST',`dnl
ifelse($3,,`dnl
forloop(`i',$1,decr($2),i`,')$2`'dnl
',`dnl
forloop(`i',0,decr(eval(($2-$1)/$3)),`eval($1+i*$3),')eval($1+eval(($2-$1)/$3)*$3)`'dnl
')dnl
')dnl
dnl
dnl
define(`RSB_M4_ERROR_UNIMPLEMENTED',`#error "missing implementation! Contact the author!"')dnl
dnl
define(`RSB_M4_MATRIX_TYPES_ARRAY',`RSB_MATRIX_TYPES_ARRAY')dnl
define(`RSB_M4_MATRIX_META_OPS_ARRAY',`RSB_MATRIX_OPS_ARRAY')dnl dnl
dnl
define(`RSB_M4_ZEROS_ARRAY',`dnl
`{'`0'forloop(`__dummy',0,decr(eval($1-1)),``,0'')`}'dnl
')dnl dnl
dnl
dnl	RSB_M4_MAKE_FUNCTION_POINTER_TABLE()
dnl	-----------------------------------------------------------------------
dnl	The resulting table should be easily addressable by C code.
dnl
dnl	UNFINISHED : YOU COULD DELETE THIS CODE NOW AND NO ONE WOULD NOTICE
dnl
define(`RSB_M4_MAKE_FUNCTION_POINTER_TABLE',`dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
foreach(`mtype',RSB_M4_MATRIX_TYPES,`dnl
dnl RSB_M4_DIRECT_KERNEL_DISPATCH_FULLRANGEBENCHMARK_FUNCTION_ARGS(mop,mtype)dnl
')dnl
')dnl
')dnl
dnl
dnl
define(`RSB_M4_IS_BCSR',dnl
pushdef(`rows_unroll',$1)dnl
pushdef(`cols_unroll',$2)dnl
dnl
ifelse(`ifelse(rows_unroll,1,1,0),ifelse(cols_unroll,1,1,0)',11,1,0)`'dnl
dnl
popdef(`rows_unroll')dnl
popdef(`cols_unroll')dnl
)dnl
dnl
define(`RSB_M4_FAKE_DIAG_IMPLICIT_MSG',`/* NOTE: Diagonal implicit is not really handled here: look at caller level. */')dnl
dnl
dnl	------------------------------------------------- 20110206
dnl
define(`RSB_M4_EARLY_EVICT_INSTRUCTION',`dnl
dnl	foreach(`location',$1,_mm_prefetch(location+24,_MM_HINT_NTA);
dnl	)`'dnl
')dnl dnl
dnl
define(`RSB_M4_HEADER_EXTRA_DECLARATIONS',`dnl
dnl	#include <xmmintrin.h>
')dnl dnl
dnl
dnl	------------------------------------------------- 20110206
dnl
dnl	define(`RSB_M4_FORTRAN_SYMBOL_ADD_TO_C',`_f_')dnl
ifelse(RSB_M4_FORTRAN_CONVENTION,`xlf',`dnl
define(`RSB_M4_FORTRAN_SYMBOL_ADD_TO_C',`')dnl
',`dnl
define(`RSB_M4_FORTRAN_SYMBOL_ADD_TO_C',`_')dnl
')dnl
define(`RSB_M4_FORTRAN_SYMBOL_ADD_TO_F',`')dnl
define(`RSB_M4_FORTRAN_SYMBOL_PREPEND_TO_C',`')dnl
dnl
dnl
dnl
define(`RSB_M4_INCLUDE_HEADERS',`dnl
#include "rsb.h"
#include "rsb_common.h"
#include "rsb_internals.h"
')dnl
dnl
dnl
