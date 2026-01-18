dnl
dnl
dnl	@author: Michele Martone
dnl
ifelse(LIBMMVBR_INCLUDED_TYPES_M4,1,`',`
define(`LIBMMVBR_INCLUDED_TYPES_M4',`1')dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
include(`libspblas_macros.m4')dnl
dnl
dnl 
dnl 
dnl	FIXME : this should go out of here:
dnl
define(`RSB_M4_MATRIX_OPS',(WANT_MATRIX_OPS))dnl
dnl
dnl 
dnl 
dnl 
/** @file
    @brief
    Macros and constants, which are type specific and used by rsb.h.
    \n
    Here reside declarations related to supported matrix numerical types, and other declarations
    according to the build time options.
    \n
    If you wish to use this library with different matrix numerical types, you shall regenerate
     the library source code accordingly; see the README file how to do this.
    \n
    Only a those declarations which are documented are meant to be part of the API.
    \n
    The rest is meant as internals.
    \n
    See also \ref matrix_type_symbols_section.
  */
RSB_M4_HEADER_MESSAGE()dnl
dnl 
dnl 	FIXME: move RSB_TYPES to RSB_CONFIG
dnl 
ifdef(`ONLY_WANT_HEADERS',`dnl
#ifndef RSB_TYPES_H_INCLUDED
#define RSB_TYPES_H_INCLUDED
')dnl
dnl 
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

dnl	***********************************************************************
dnl
/* 
   Each of the following symbols corresponds to a type opted in or out at code generation time.
   Other types may be enabled by regenerating the whole library code.
 */

/* Miscellaneous version strings.
dnl  Adopting a naming scheme similar to that of png.h.
 */
`#define 'RSB_LIBRSB_VER_STRING`		'"RSB_M4_WANT_LIBRSB_VER_MAJOR.RSB_M4_WANT_LIBRSB_VER_MINOR.RSB_M4_WANT_LIBRSB_VER_PATCH`'"`'	/*!< \brief Library version string. */
`#define 'RSB_HEADER_VERSION_STRING`		'"librsb version RSB_M4_WANT_LIBRSB_VER_MAJOR.RSB_M4_WANT_LIBRSB_VER_MINOR.RSB_M4_WANT_LIBRSB_VER_PATCH`'RSB_M4_WANT_LIBRSB_VER_PRERS` - 'RSB_M4_WANT_LIBRSB_VER_DATE"`'	/*!< \brief Library header version string. */
`#define 'RSB_LIBRSB_VER_MAJOR`		'RSB_M4_WANT_LIBRSB_VER_MAJOR`'	/*!< \brief Major version. */
`#define 'RSB_LIBRSB_VER_MINOR`		'RSB_M4_WANT_LIBRSB_VER_MINOR`'	/*!< \brief Minor version. */
`#define 'RSB_LIBRSB_VER_PATCH`		'RSB_M4_WANT_LIBRSB_VER_PATCH`'	/*!< \brief Patch version. */
`#define 'RSB_LIBRSB_VER`		'RSB_M4_WANT_LIBRSB_LIBRSB_VER`'	/*!< \brief Version number. */
`#define 'RSB_LIBRSB_VER_DATE`		'RSB_M4_WANT_LIBRSB_VER_DATE`'	/*!< \brief Version release date. */

dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
dnl
`#define' RSB_M4_HAVE_TYPE_PREPROCESSOR_SYMBOL(type) 1 /*!< \brief Type type is supported, so RSB_M4_HAVE_TYPE_PREPROCESSOR_SYMBOL(type) is defined .*/
dnl
')dnl
`#define' RSB_DEFAULT_TYPE RSB_M4_DEFAULT_TYPE	/*!< \brief The default numerical matrix type (can be used for declarations), used in the example programs. */
`#define' RSB_DEFAULT_POSSIBLY_INTEGER_TYPE RSB_M4_DEFAULT_POSSIBLY_INTEGER_TYPE/*!< \brief The default, integer if possible , numerical type (can be used for declarations). */
`#define' RSB_DEFAULT_POSSIBLY_FIRST_BLAS_TYPE RSB_M4_FIRST(RSB_M4_DEFAULT_POSSIBLY_BLAS_TYPE_OR_DEFAULT) ` '/*!< \brief The default, blas if possible , numerical type (can be used for declarations). */
`#define' RSB_DEFAULT_TYPE_STRING RSB_M4_QUOTED_COMMA_LIST((RSB_M4_DEFAULT_TYPE))	/*!< \brief A string specifying the name of the default type. */
`#define' RSB_DEFAULT_POSSIBLY_INTEGER_TYPE_STRING RSB_M4_QUOTED_COMMA_LIST((RSB_M4_DEFAULT_POSSIBLY_INTEGER_TYPE)) /*!< \brief A string specifying the name of the default possibly integer type.*/
`#define' RSB_DEFAULT_SYMMETRY RSB_M4_MATRIX_SYMMETRY_PREPROCESSOR_SYMBOL(RSB_M4_DEFAULT_SYMMETRY)	/*!< \brief The default symmetry flag. */
`#define' RSB_DEFAULT_TRANSPOSITION RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL(RSB_M4_DEFAULT_TRANSPOSITION)	/*!< \brief The default transposition flag (no transposition). */
dnl
dnl
ifelse(RSB_M4_HAVE_COMPLEX_TYPE,`()',`',`dnl
#define RSB_HAVE_ANY_COMPLEX_TYPE 1	/*!< \brief Defined to 1 if any complex type has been configured in. Currently: 'RSB_M4_HAVE_COMPLEX_TYPE`. */
')dnl
dnl
`#define RSB_ROWS_TRANSPOSITIONS_ARRAY	{'dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL(transposition), ') RSB_TRANSPOSITION_INVALID} /*!< \brief An array with transposition constants. */
`#define RSB_TRANSPOSITIONS_ARRAY_LENGTH 3 /* valid transpositions in RSB_ROWS_TRANSPOSITIONS_ARRAY */'

dnl
pushdef(`counter',`0')dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
/*!  This preprocessor index can be used to address the type-related arrays.  */
`#define' RSB_M4_TYPE_INDEX_PREPROCESSOR_SYMBOL(type) counter
pushdef(`counter',eval(counter+1))dnl
')dnl
foreach(`type',RSB_M4_MATRIX_TYPE,`popdef(`counter')')dnl
popdef(`counter')dnl

dnl
dnl
dnl	***********************************************************************
/* @cond INNERDOC  */
dnl
/*
   Each of the following symbols corresponds to an operation opted in or out at code generation time.
   \n
   Other operations may be enabled by regenerating the whole library code.
 */
dnl
foreach(`mop',RSB_M4_MATRIX_OPS,`dnl
dnl
`#define' RSB_M4_HAVE_OPTYPE_PREPROCESSOR_SYMBOL(mop) 1
dnl
')dnl
dnl

/*!
 * These preprocessor indices can be used to address various mop-related arrays.
 */
dnl
pushdef(`counter',`0')dnl
foreach(`mop',RSB_M4_MATRIX_META_OPS,`dnl
`#define' RSB_M4_OPTYPE_INDEX_PREPROCESSOR_SYMBOL(mop) counter
pushdef(`counter',eval(counter+1))dnl
')dnl
foreach(`mop',RSB_M4_MATRIX_META_OPS,`popdef(`counter')')dnl
popdef(`counter')dnl
dnl
dnl

dnl
/**
 \name Values for valid matrix coordinate index types flags.
 */
dnl
foreach(`citype',RSB_M4_MATRIX_COORDINATE_TYPES,`dnl
`#define ' RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_PREPROCESSOR_SYMBOL(citype) RSB_M4_MATRIX_INDEX_COORDINATE_TYPE_CHARCODE_(citype) /*!< \brief Character code for type citype.*/
')dnl
dnl
/* @endcond */
dnl
/**
 \name Values for matrix transposition flags (rsb_trans_t).
 \anchor matrix_transposition_flags_section
 Note that for non complex types, the Hermitian flag will act as simple transposed.
 */
dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
`#define ' RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL(transposition) RSB_M4_MATRIX_TRANSPOSITION_CHARCODE(transposition)
')dnl
`#define ' RSB_TRANSPOSITION_INVALID RSB_M4_MATRIX_TRANSPOSITION_CHARCODE(`RSB_M4_TRANS_INVALID')
dnl
/* @cond INNERDOC  */
dnl
/**
 \name Values for matrix symmetry flags (rsb_flags_t).
 \anchor matrix_symmetry_flags_section
 */
dnl
foreach(`symmetry',RSB_M4_MATRIX_SYMMETRY,`dnl
`#define ' RSB_M4_MATRIX_SYMMETRY_PREPROCESSOR_SYMBOL(symmetry) RSB_M4_MATRIX_SYMMETRY_CHARCODE(symmetry)
')dnl
dnl
/* @endcond */
dnl
/**
dnl \name Values for valid matrix symmetry flags.
\name Values for diagonal specification flags (rsb_flags_t).
 \anchor matrix_diagonal_flags_section
 */
dnl
/* @cond INNERDOC  */
dnl
foreach(`diagonal',RSB_M4_MATRIX_DIAGONAL_TYPES,`dnl
`#define ' RSB_M4_MATRIX_DIAGONAL_PREPROCESSOR_SYMBOL(diagonal) RSB_M4_MATRIX_DIAGONAL_CHARCODE(diagonal) /*!< \brief */
')dnl
dnl
/* @endcond */
dnl
/* @cond INNERDOC  */
dnl
/**
 \name Values for valid matrix storage formats.
 \anchor matrix_storage_flags_section
 */
dnl
foreach(`matrix_storage',RSB_M4_MATRIX_STORAGE,`dnl
`#define ' RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(matrix_storage) RSB_M4_MATRIX_STORAGE_CHARCODE(matrix_storage)
')dnl
dnl
/**
 \name Values for valid matrix storage formats strings.
 \anchor matrix_storage_strings_section
 */
dnl
foreach(`matrix_storage',RSB_M4_MATRIX_STORAGE,`dnl
`#define ' RSB_M4_MATRIX_STORAGE_PREPROCESSOR_STRING(matrix_storage) "touppercase(RSB_M4_CHOPSPACES(matrix_storage))"
')dnl
dnl
/* @endcond */
dnl

dnl
/**
 \name Supported matrix numerical types.
 \anchor matrix_supported_numerical_types_section
 */
`#define 'RSB_MATRIX_TYPES_LIST_CXX`		'RSB_M4_MATRIX_TYPES_LIST_CXX /*!< \brief list of C++ types configured in this \librsb build, usable in \c <rsb.hpp> (since RSB_LIBRSB_VER>=10300) */

dnl
/**
 \name Valid symbol values for matrix numerical type specification (type codes).
 \anchor matrix_type_symbols_section
 */
dnl
foreach(`citype',RSB_M4_MATRIX_TYPES,`dnl
`#define' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(SAME_TYPE) 1 /*!< \brief a bogus type flag for specifying no type conversion */
`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(citype) singlequote(RSB_M4_TYPE_CHARCODE(citype)) /*!< \brief Character code for type citype. See \ref rsb_type_t . */
')dnl

`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(`FORTRAN_'`SAME_TYPE') 1 /*!< \brief a bogus type flag for specifying no type conversion */
foreach(`citype',RSB_M4_ALL_MATRIX_TYPES,`dnl
`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(`FORTRAN_'citype) ICHAR(singlequote(RSB_M4_TYPE_CHARCODE(citype))) /*!< \brief Character code for type citype, to be used (only) from Fortran. */
')dnl

`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(`DEFAULT')  RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(RSB_M4_DEFAULT_TYPE)  /*!< \brief A default numerical matrix type. */
`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(`DEFAULT_INTEGER')  RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(RSB_M4_DEFAULT_POSSIBLY_INTEGER_TYPE)  /*!< \brief A default numerical matrix type; if possible, an integer one. */
`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(RSB_M4_INVALID_TYPE) singlequote(RSB_M4_TYPE_CHARCODE(RSB_M4_INVALID_TYPE)) /*!< \brief By definition, an invalid type code. */
`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(`FIRST_BLAS')  RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(RSB_M4_DEFAULT_POSSIBLY_BLAS_TYPE)  /*!< \brief A default numerical matrix type; if possible, not integer one. If no such type is configured in, then the invalid type. */

`#define ' RSB_CHAR_AS_TRANSPOSITION(TRANSC)	\
(														\
foreach(`transA',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
		(TRANSC) == (touppercase(singlequote(transA))) ? (RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL(transA)) : 		\
		(TRANSC) == (tolowercase(singlequote(transA))) ? (RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL(transA)) : 		\
')dnl
		singlequote(RSB_M4_TRANS_INVALID)	\
) /*!< \brief Get the right transposition flag out of either n, c, t chars. */


/**
 \name Miscellaneous constants.
 */
dnl
#define RSB_CONST_MAX_TUNING_ROUNDS 16 /*!< \brief Maximal count of tuning rounds in one invocation of (rsb_tune_spmm/rsb_tune_spsm). */

dnl
/* @cond INNERDOC  */
dnl
/**
 \name Values for other numerical type related macros.
*/
`#define ' RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS "foreach(`type',RSB_M4_MATRIX_TYPES,`RSB_M4_TYPE_CHARCODE(type) ')"
`#define ' RSB_BLAS_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS "foreach(`type',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`RSB_M4_TYPE_CHARCODE(type) ')"

/* a bogus type for pattern input (TODO : should also implement ANY, just for matrix input) */
`#define' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(PATTERN) 0
dnl
/* @endcond */
dnl
/* @cond INNERDOC */
dnl
dnl

dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
`#define ' RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(type) dnl
dnl
ifelse(type,`int',"%d")dnl
ifelse(type,`char',"%c")dnl
dnl
ifelse(type,`float',"%.9g")dnl
ifelse(type,`float complex',"%.9g %.9g")dnl
ifelse(type,`double',"%.17g")dnl
ifelse(type,`double complex',"%.17g %.17g")dnl
ifelse(type,`long double',"%.17Lg")dnl
ifelse(type,`long double complex',"%.17Lg %.17Lg")dnl
dnl
dnl ifelse(type,`long double',"%Lg")dnl
dnl ifelse(type,`double',"%lg")dnl
dnl ifelse(type,`float',"%g")dnl
dnl ifelse(type,`complex',"%g %g")dnl
dnl ifelse(type,`long double complex',"%Lg %Lg")dnl
dnl ifelse(type,`double complex',"%lg %lg")dnl
dnl ifelse(type,`float complex',"%g %g")dnl

')dnl
dnl
dnl
dnl	#define RSB_SIZEOF(type) rsb__do_sizeof(type)
dnl
dnl
dnl	UNLOOP_PAIRS
dnl	------------
dnl
define(`UNLOOP_PAIRS',`foreach(`type',RSB_M4_MATRIX_TYPES,`RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),sizeof(type),')')
dnl
dnl define(`UNLOOP_PAIRS_ALL',`foreach(`type',RSB_M4_ALL_MATRIX_TYPES,`RSB_M4_TYPE_CHARCODE_ASCII_VALUE(type),sizeof(type),')') dnl depends on complex.h
define(`UNLOOP_PAIRS_ALL_GUESSED',`foreach(`type',RSB_M4_ALL_MATRIX_TYPES,`RSB_M4_TYPE_CHARCODE_ASCII_VALUE(type),RSB_M4_BACKUP_SIZEOF(type),')')
dnl
define(`REALT__PAIRS',`foreach(`type',RSB_M4_MATRIX_TYPES,`RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(RSB_M4_REALT(type)),')')
#if 1
dnl
dnl
dnl	SINGLE_LINEAR_SEARCH
dnl	--------------------
dnl	Expands to a linear search implemented in the C macro preprocessor language.
dnl
define(`SINGLE_LINEAR_SEARCH',`ifelse($#,1,$1,`( (TYPE)==($1) ?  pushdef(`type',$1)$2 popdef(`type'): \
	(SINGLE_LINEAR_SEARCH(shift(shift($@))) ) ) ')') 
dnl 
dnl
dnl
`#define RSB_ROWS_TRANSPOSITIONS_ARRAY_AS_CHAR	{'dnl
foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`singlequote(transposition), ')singlequote(RSB_M4_TRANS_INVALID) }


`#define ' RSB_TRANSPOSITIONS_PREPROCESSOR_SYMBOLS "foreach(`transposition',RSB_M4_MATRIX_TRANSPOSITIONS,`transposition ')"

#define RSB_TRANSPOSITION_AS_CHAR(TRANSA) 										\
(														\
foreach(`transA',RSB_M4_MATRIX_TRANSPOSITIONS,`dnl
		(TRANSA) == (RSB_M4_MATRIX_TRANSPOSITION_PREPROCESSOR_SYMBOL(transA)) ? (touppercase(singlequote(transA))) : 		\
')dnl
		singlequote(RSB_M4_TRANS_INVALID)	\
)

`#define RSB__ORDER_AS_SINGLE_CHAR(ORDER) ( ((ORDER) == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER) ?  'singlequote(C)` : 'singlequote(R)` )'
`#define RSB__ORDER_AS_STRING(ORDER) ( ((ORDER) == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER) ?  "cols" : "rows" )'
`#define RSB__ORDER_AS_LANG_CHAR(ORDER) ( ((ORDER) == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER) ?  'singlequote(F)` : 'singlequote(C)` )'

#define RSB_NUMERICAL_TYPE_STRING(CSP,TYPE) \
		{ \
		switch(TYPE) \
		{ \
			/* supported RSB_M4_MATRIX_TYPES */ \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:CSP="RSB_M4_CHOPSPACES(type)";break; 	\
')dnl
			/* unsupported type */ \
			default : CSP="RSB_M4_TRANS_INVALID"; \
		} \
		}



#define RSB_NUMERICAL_TYPE_SIZE(TYPE) \
	SINGLE_LINEAR_SEARCH( UNLOOP_PAIRS 0 )

#define RSB_SIZEOF_BACKUP(TYPE) /* This is for rsb__pr_load. Please feed in upper case char codes (toupper(...)). */ \
    	SINGLE_LINEAR_SEARCH( UNLOOP_PAIRS_ALL_GUESSED 0 )

#define RSB_NUMERICAL_TYPE_REAL_TYPE(TYPE) \
	SINGLE_LINEAR_SEARCH( REALT__PAIRS 0 )
dnl

#define RSB_NUMERICAL_TYPE_CAST_TO_ANY_P(CTYPE,CVAR,TYPE,TP,TOFF) \
		{ \
		switch(TYPE) \
		{ \
			/* supported RSB_M4_MATRIX_TYPES */ \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:\
				(CVAR)=(CTYPE)((type*)TP)[TOFF] ; break; 	\
')dnl
			/* unsupported type */ \
			default : ; \
		} \
		}

/* *A += abs(*B) */
#define RSB_NUMERICAL_TYPE_ABS_SUM_AND_STORE_ELEMENTS(A,B,TYPE) \
		{ \
		switch(TYPE) \
		{ \
			/* supported RSB_M4_MATRIX_TYPES */ \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:	*(type*)(A)+= (	\
				*(type*)(B) < (type)(0) ? - *(type*)(B) : *(type*)(B) ); break; 	\
')dnl
			/* unsupported type */ \
			default : ; \
		} \
		}

/* *A += *B */
#define RSB_NUMERICAL_TYPE_SUM_AND_STORE_ELEMENTS(A,B,TYPE) \
		{ \
		switch(TYPE) \
		{ \
			/* supported RSB_M4_MATRIX_TYPES */ \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:	*(type*)(A)+=*(type*)(B); break; \
')dnl
			/* unsupported type */ \
			default : ; \
		} \
		}

#define RSB_NUMERICAL_TYPE_SET_ELEMENT(DST,SRC,TYPE) \
		{ \
		switch(TYPE) \
		{ \
			/* supported RSB_M4_MATRIX_TYPES */ \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:	*(type*)(DST)=*(type*)(SRC); break; \
')dnl
			/* unsupported type */ \
			default : ; \
		} \
		}

#define RSB_NUMERICAL_TYPE_SET_ELEMENT_REAL(DST,SRC,TYPE) \
		{ \
		switch(TYPE) \
		{ \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:	*(RSB_M4_REALT(type)*)(DST)=RSB_M4_CREAL(type,*(type*)(SRC)); break; \
')dnl
			/* unsupported type */ \
			default : ; \
		} \
		}

#define RSB_NUMERICAL_TYPE_SET_ELEMENT_FROM_DOUBLE(DST,DSRC,TYPE) \
		{ \
		switch(TYPE) \
		{ \
			/* supported RSB_M4_MATRIX_TYPES */ \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:	*(type*)(DST)=(type)(DSRC); break; \
')dnl
			/* unsupported type */ \
			default : ; \
		} \
		}

/* CODE NOT DEBUGGED */
#define RSB_VECTOR_FIND_MAXIMAL_ELEMENT(INDEX,ARRAY,ELEMENTS,TYPE) 								\
		{ 													\
		int _index;												\
		switch(TYPE) 												\
		{ 													\
			/* supported RSB_M4_MATRIX_TYPES */ 									\
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:						\
			{												\
				type * _array = (type*)(ARRAY);								\
				type _maxel=(type)(0);									\
				int  _maxindex=0;									\
				_maxel=_maxel-_maxel;	/* could this be evil ? */					\
				for(_index=0;_index<(ELEMENTS);++_index)						\
					if(RSB_M4_ABS(type,_maxel)<RSB_M4_ABS(type,_array[_index])){_maxel=_array[_index];_maxindex=_index;}	\
					(INDEX)=_maxindex;								\
			}												\
			break;			\
')dnl
			/* unsupported type */ \
			default :  (INDEX)=-1; \
		} \
		}

dnl
dnl	***********************************************************************
dnl
#define RSB_NUMERICAL_OP_INDEX_FROM_CODE(CODE) 								\
pushdef(`code',`0')dnl
foreach(`op',RSB_M4_MATRIX_META_OPS,`dnl
( ((CODE)==RSB_M4_OPTYPE_INDEX_PREPROCESSOR_SYMBOL(op))?(code):			\
pushdef(`code',eval(code+1))dnl
')dnl
-1 dnl
foreach(`op',RSB_M4_MATRIX_META_OPS,`dnl
popdef(`code')dnl
) \
')dnl
popdef(`code')dnl
dnl
/* uhm. does it seem redundant ? */
dnl
dnl	***********************************************************************
dnl
#define RSB_NUMERICAL_TYPE_INDEX_FROM_CODE(CODE) 								\
pushdef(`code',`0')dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
( ((CODE)==RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type))?(code):			\
pushdef(`code',eval(code+1))dnl
')dnl
-1 dnl
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
popdef(`code')dnl
) \
')dnl
popdef(`code')dnl
dnl
/* uhm. seems redundant ? */
dnl
dnl	***********************************************************************
dnl


#define RSB_IS_ELEMENT_MINUS_ONE(SRC,TYPE) 										\
(														\
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
		(TYPE) == (RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)) ? (*(type*)(SRC)==(type)(-1)) : 		\
')dnl
		0												\
)

#define RSB_IS_ELEMENT_ONE(SRC,TYPE) 										\
(														\
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
		(TYPE) == (RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)) ? (*(type*)(SRC)==(type)1) : 		\
')dnl
		0												\
)

#define RSB_IS_ELEMENT_ZERO(SRC,TYPE) 										\
(														\
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
		(TYPE) == (RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)) ? (*(type*)(SRC)==(type)0) : 		\
')dnl
		0												\
)

#define RSB_IS_ELEMENT_NONZERO(SRC,TYPE) 		(!(RSB_IS_ELEMENT_ZERO(SRC,TYPE)))

#define RSB_MATRIX_UNSUPPORTED_TYPE(TYPE) ( \
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			(TYPE)!=RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type) && \
')dnl
			1 )

#define RSB_IS_MATRIX_TYPE_COMPLEX(TYPE) 										\
(														\
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
		(TYPE) == (RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)) ? RSB_M4_IS_COMPLEX_TYPE(type) : 		\
')dnl
		0												\
)

#define RSB_IS_MATRIX_TYPE_BLAS_TYPE(TYPE) 										\
(/*RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST*/														\
foreach(`type',(RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST),`dnl
		(TYPE) == (RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)) ? 1 : 		\
')dnl
		0												\
)

#define RSB_IS_ELEMENT_LESS_THAN(SRC,CMPSRC,TYPE) \
( foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			( (TYPE)==RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type) && RSB_M4_CREAL(type,*(type*)(SRC))<RSB_M4_CREAL(type,*(type*)(CMPSRC)) ) || \
')dnl
			0 )


dnl
dnl
/** use RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE to oversize your arrays safely */
`#define RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE	' RSB_M4_MAX2(RSB_M4_MAXN(WANT_COLUMN_UNLOOP_FACTORS),RSB_M4_MAXN(WANT_ROW_UNLOOP_FACTORS)) 
dnl
dnl
dnl
/** use RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE_EXTRA to oversize your arrays safely */
`#define RSB_MAXIMAL_CONFIGURED_BLOCK_SIZE_EXTRA	' (RSB_M4_MAXIMAL_CONFIGURED_BLOCK_SIZE-1)
dnl
dnl
#define RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH (2*1024)	/** chars to reserve for a matrix implementation code */

/* Section dedicated to implemented operations on matrices. */



dnl
dnl
`#define RSB_ROWS_UNROLL_ARRAY		{' RSB_M4_COMMA_LIST(RSB_M4_COLUMNS_UNROLL) }
dnl
`#define RSB_COLUMNS_UNROLL_ARRAY	{' RSB_M4_COMMA_LIST(RSB_M4_ROWS_UNROLL) }
dnl


`#define RSB_ROWS_UNROLL_ARRAY_LENGTH		'RSB_M4_LIST_LENGTH(WANT_ROW_UNLOOP_FACTORS)
`#define RSB_COLUMNS_UNROLL_ARRAY_LENGTH		'RSB_M4_LIST_LENGTH(WANT_COLUMN_UNLOOP_FACTORS)
`#define RSB_IMPLEMENTED_META_MOPS		'RSB_M4_LIST_LENGTH(RSB_M4_QUOTED_COMMA_LIST(RSB_M4_MATRIX_META_OPS))
`#define RSB_IMPLEMENTED_MOPS		'RSB_M4_LIST_LENGTH(RSB_M4_QUOTED_COMMA_LIST(RSB_M4_MATRIX_OPS))
`#define RSB_IMPLEMENTED_TYPES		'RSB_M4_LIST_LENGTH(WANT_TYPES)
dnl `#define RSB_IMPLEMENTED_BLAS_TYPES		'RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST_LENGTH
`#define RSB_IMPLEMENTED_SOME_BLAS_TYPES		'RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST_LENGTH

`#define 'RSB_M4_MATRIX_META_OPS_ARRAY`	{' RSB_M4_QUOTED_COMMA_LIST(RSB_M4_MATRIX_META_OPS) }
dnl
dnl NOTE: the following maps "double complex" to "double,complex"
dnl `#define 'RSB_M4_MATRIX_TYPES_ARRAY`	{' RSB_M4_QUOTED_COMMA_LIST(RSB_M4_MATRIX_TYPES) }
dnl
`#define 'RSB_M4_MATRIX_TYPES_ARRAY`	{' foreach(`type',RSB_M4_MATRIX_TYPES,`"type",') }
`#define 'RSB_MATRIX_TYPE_CODES_ARRAY`	{' foreach(`type',RSB_M4_MATRIX_TYPES,`RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),') }
`#define 'RSB_MATRIX_SPBLAS_TYPE_CODES_ARRAY`	{' foreach(`type',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type),') }

/* Trick to exercise lightweight Sparse BLAS Fortran wrappers */
`#define 'RSB__BLAS_CALL_TF_(TYPECODE,CALL,ISTAT,...)`	( \'dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`
	( TYPECODE==RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type) ) ? RSB__BLAS_CALL_F(tolowercase(RSB_M4_TYPE_CHARCODE(type)),CALL,__VA_ARGS__,&ISTAT),ISTAT : ( \')
	RSB_BLAS_ERROR \
	foreach(`type',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`)') `)'

`#define RSB__FORTRAN_APPEND_UNDERSCORE 'ifelse(RSB_M4_FORTRAN_CONVENTION,`xlf',`0',`1')

`#define 'RSB_M4_MATRIX_META_OPS_STRING`	'"WANT_MATRIX_OPS"
`#define 'RSB_M4_MATRIX_TYPES_STRING`		'"WANT_TYPES"
`#define 'RSB_M4_WANT_COLUMN_UNLOOP_FACTORS_STRING`		'"RSB_M4_SPACED_LIST((WANT_COLUMN_UNLOOP_FACTORS))"
`#define 'RSB_M4_WANT_ROW_UNLOOP_FACTORS_STRING`		'"RSB_M4_SPACED_LIST((WANT_ROW_UNLOOP_FACTORS))"

/**
 \name Macro to check matrix storage flags correctness
 */
dnl
`#define ' RSB_IS_MATRIX_STORAGE_ALLOWED_FOR_LEAF(MATRIX_STORAGE)	(( \
foreach(`matrix_storage',RSB_M4_MATRIX_STORAGE,`dnl
	((MATRIX_STORAGE)==RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(matrix_storage)) || \
')dnl
	0 ) ? RSB_BOOL_TRUE:RSB_BOOL_FALSE )
dnl

/**
 \name Short internal M4 macros check
 */
dnl
`#define ' RSB__M4_CHECK()              \
	\
        assert(RSB_M4_XOR(0,0)==0);      \
        assert(RSB_M4_XOR(0,1)==1);      \
        assert(RSB_M4_XOR(1,0)==1);      \
        assert(RSB_M4_XOR(1,1)==0);      \
	\
        assert(RSB_M4_OR(0,0)==0);      \
        assert(RSB_M4_OR(0,1)==1);      \
        assert(RSB_M4_OR(1,0)==1);      \
        assert(RSB_M4_OR(1,1)==1);      \
	\
        assert(RSB_M4_AND(0,0)==0);      \
        assert(RSB_M4_AND(0,1)==0);      \
        assert(RSB_M4_AND(1,0)==0);      \
        assert(RSB_M4_AND(1,1)==1);      \
	\
        assert(RSB_M4_IMPLY(0,0)==1);      \
        assert(RSB_M4_IMPLY(0,1)==1);      \
        assert(RSB_M4_IMPLY(1,0)==0);      \
        assert(RSB_M4_IMPLY(1,1)==1);      \
	\
        assert(RSB_M4_NOT(0)==1);      \
        assert(RSB_M4_NOT(1)==0);      \
	\
        assert(RSB_M4_MAX2(0,1)==1);      \
        assert(RSB_M4_MAX2(1,0)==1);      \
        assert(RSB_M4_MIN2(0,1)==0);      \
        assert(RSB_M4_MIN2(1,0)==0);      \
	\
        assert(RSB_M4_MAXN(0,1)==1);      \
        assert(RSB_M4_MAXN(1,0)==1);      \
        assert(RSB_M4_MINN(0,1)==0);      \
        assert(RSB_M4_MINN(1,0)==0);      \
	\
        assert(RSB_M4_MAXN(0,1,2)==2);      \
        assert(RSB_M4_MAXN(1,2,0)==2);      \
        assert(RSB_M4_MAXN(2,0,1)==2);      \
	\
        assert(RSB_M4_MINN(0,1,2)==0);      \
        assert(RSB_M4_MINN(1,2,0)==0);      \
        assert(RSB_M4_MINN(2,0,1)==0);      \
        ;   dnl TODO: many more: RSB_M4_MEMBER, RSB_M4_FIRST, RSB_M4_SORT
dnl

ifelse(RSB_M4_LONG_IDX,`0',`',`dnl
#ifndef RSB_WANT_LONG_IDX_TYPE
#define RSB_WANT_LONG_IDX_TYPE int64_t
#else
#error RSB_WANT_LONG_IDX_TYPE shall be defined by rsb_types.h, not elsewhere!
#endif
')dnl

`#define RSB_HAVE_IHI' RSB_M4_WANT_IHI /* is rsb-librsb-internals.h installed ? */

#ifdef __cplusplus
}
#endif  /* __cplusplus */

dnl
#endif
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl
dnl	FIXME: WHY THE HECK UNCOMMENTING THE FOLLOWING TRIGGERS AN ERROR ?
dnl
dnl #endif
#endif /* RSB_TYPES_H_INCLUDED */
')dnl
dnl
/* @endcond */
dnl
',`dnl
static int foo(){return 0;}
')dnl ONLY_WANT_HEADERS
dnl
dnl
')dnl the whole file
dnl
dnl
