dnl
dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
dnl
ifelse(LIBMMVBR_INCLUDED_UTIL_M4,1,`',`
define(`LIBMMVBR_INCLUDED_TYPES_M4',`1')dnl
dnl `PACK' format will be experimented with in the future :)
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
include(`libspblas_macros.m4')dnl RSB_M4_SPBLAS...
/**
 * @file
 * @brief
 * Auxiliary functions.
 */
RSB_M4_HEADER_MESSAGE()dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_UTIL_H_INCLUDED
#define RSB_UTIL_H_INCLUDED
')
dnl
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

dnl
#include "rsb_common.h"
dnl #include "rsb_types.h"
dnl 
dnl
dnl
dnl	FIXME : COMMENT THIS FILE
dnl	-------------------------
dnl
dnl
dnl
/* non blas-like functions */
dnl

dnl
rsb_err_t rsb__util_m4_sanity_check(void)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
		There are bugs in the m4 macros or a bad m4 implementation which will trigger this test to fail.
		We are interested in catching them, as we should rely on a sane m4 environment.
	*/
	RSB_M4_DEBUGINFO(``$0'')
	if(
		RSB_M4_XOR(0,0)!=0 ||
		RSB_M4_XOR(1,0)!=1 || 
		RSB_M4_XOR(0,1)!=1 || 
		RSB_M4_XOR(1,1)!=0 ||
		RSB_M4_AND(0,0)!=0 ||
		RSB_M4_AND(1,0)!=0 ||
		RSB_M4_AND(0,1)!=0 ||
		RSB_M4_AND(1,1)!=1 ||
		RSB_M4_OR(0,0)!=0 ||
		RSB_M4_OR(1,0)!=1 ||
		RSB_M4_OR(0,1)!=1 ||
		RSB_M4_OR(1,1)!=1 ||
		0
		)
		goto err;
	return RSB_ERR_NO_ERROR;
err:
	return RSB_ERR_INTERNAL_ERROR;
}
')dnl
dnl

dnl
const void * rsb__util_increase_by_one(void *p, rsb_nnz_idx_t n, rsb_type_t typecode)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ) {(((mtype*)p)[n])+=1;return p;}
	else 
#endif
')dnl
	return NULL;
}
')dnl
dnl

dnl
void rsb__util_set_area_to_fraction_of_integer(void *p, const int alphai, rsb_type_t typecode)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*
		alpha NULL will imply 1
	*/
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ) {*(mtype*)p = 1;*(mtype*)p/=alphai;}
	else 
#endif
')dnl
	return;
}
')dnl
dnl

dnl
void rsb__util_set_area_to_negated_fraction(void *p, const void *alpha, rsb_type_t typecode)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*
		alpha NULL will imply 1
	*/
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ) {*(mtype*)p = -1;if(alpha)*(mtype*)p/=(*(mtype*)alpha);}
	else 
#endif
')dnl
	return;
}
')dnl
dnl

dnl
void rsb__util_set_area_to_converted_integer(void *p, rsb_type_t typecode, const rsb_int n)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ) {*(mtype*)p = (mtype)n;}
	else 
#endif
')dnl
	return;
}
')dnl
dnl

dnl
rsb_coo_idx_t * rsb__util_get_partitioning_array( size_t bs, size_t X , rsb_blk_idx_t * X_b, rsb_flags_t flags)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * Given a block size (be it rows or columns), an element size X in bytes,
	 * and a dimension (rows or columns), returns an array containing the 
	 * indices of the elements in each block.
	 *
	 * Therefore, the allocated arrays 
	 *
	 * \param bs	the block size
	 * \param X	the rows or columns count
	 * \param X_b	on output, the allocated array elements count : (X+bs-1)/bs
	 * \return NULL on error;  a valid array pointer on success
	 *
	 * FIXME : why not size_t ? or maybe rsb_size_t ?
	 * */
	size_t i;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t * p_x = NULL;

	*X_b = (X+bs-1)/bs;

	/* WARNING : 1 is the extreme limit before overflow :) */
	if( ( ((size_t)(*X_b)) < ((size_t)((X+bs-1)/bs))) || (RSB_BLK_ADD_OVERFLOW(*X_b,1)) )
	{
		/* overflow. should print some message. */
		errval = RSB_ERR_LIMITS;goto err;
	}

	p_x = rsb__malloc(sizeof(rsb_coo_idx_t)*(*X_b+1));
	if(! p_x) goto err;
	/* note: should use some perrno some day */

	/* note the last block size : it is the same, regardless congruences */
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`*X_b',` p_x[i+LI] = (i+LI)*bs;
	')
dnl	for(i = 0;i<*X_b;++i)p_x[i] = i*bs;

	/* FIXME : this point should be remarked and documented way better ! */
	if(flags&(RSB_FLAG_WANT_BCSS_STORAGE|RSB_FLAG_WANT_FIXED_BLOCKING_VBR))
		p_x[*X_b] = *X_b*bs;	/* the last element of p_x is the index of the last matrix row/column    + 1  */
	else
		p_x[*X_b] = X;	/* the last element of p_x is the index of the last matrix row/column    + 1  */
	
	return p_x;
err:
	RSB_CONDITIONAL_FREE(p_x);
	rsb__do_perror(NULL,errval);
	return NULL;
}
')dnl
dnl

dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
dnl
rsb_err_t rsb__vector_diff(void * c, const void * a, const void * b, rsb_type_t type, size_t n)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * c <- a-b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see daxpy,dcopy in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype*ta = a,*tb = b;mtype *tc = c;
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		tc[i+LI] = ta[i+LI]-tb[i+LI];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')
dnl
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_vector_norm_square(void * c, const void * a, rsb_type_t type, size_t n)
{
	/*!
	 * c <- a^T*a
         *
	 * \param a	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype*ta = a;mtype *tc = c;
		tc[0] = RSB_M4_ZERO(mtype);
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		tc[0]+=ta[i+LI]*ta[i+LI];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_vector_norm(void * c, const void * a, rsb_type_t type, size_t n)
{
	/*!
	 * c <- sqrt(a^T*a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	rsb_err_t errval;
	if(!c)
		return RSB_ERR_BADARGS;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		mtype*cp = (mtype*)c;
		errval = rsb_vector_norm_square(cp,a,type,n);
		*cp = RSB_M4_SQRT(mtype,*cp);
	}
	else 
#endif
')dnl
		errval = RSB_ERR_UNSUPPORTED_TYPE;
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl

dnl
ifelse(1,0,`dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_vector_norm_square_strided(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
{
	/*!
	 * c <- a^T*a
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(inc==1)
		return rsb_vector_norm_square(c,a,type,n);
dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype*ta = a;mtype *tc = c;
		tc[0] = RSB_M4_ZERO(mtype);
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		tc[0]+=ta[(i+LI)*inc]*ta[(i+LI)*inc];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
dnl
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl
')dnl
dnl

dnl
rsb_err_t rsb__vector_norm_strided(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * c <- sqrt(a^T*a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	rsb_err_t errval;
	if(!c)
		return RSB_ERR_BADARGS;
	if(inc==1)
		return rsb_vector_norm(c,a,type,n);
ifelse(1,0,`dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		mtype*cp = (mtype*)c;
		errval = rsb_vector_norm_square_strided(cp,a,type,n,inc);
		*cp = RSB_M4_SQRT(mtype,*cp);
	}
	else 
#endif
')dnl
		errval = RSB_ERR_UNSUPPORTED_TYPE;
dnl
',`
	errval = RSB_ERR_INTERNAL_ERROR;
')dnl
dnl
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl

dnl
rsb_err_t rsb__util_vector_sum_strided(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * c <- sum(a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	/* See also rsb__cblas_Xnrm2 */

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		register mtype acc = RSB_M4_ZERO(mtype); const mtype*ta = a; mtype*tc = c;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		acc+=ta[(i+LI)*inc];
	'); 
		tc[0] = acc;
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_vector_sum(void * c, const void * a, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * c <- sum(a)
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	const mtype*ta = a; mtype*tc = c; tc[0] = RSB_M4_ZERO(mtype);
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	tc[0]+=ta[i+LI];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
dnl
dnl	rsb_err_t rsb_blas_Xdot(void * c, const void * a, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
dnl	{
dnl		cblas_ddot(n,a,1,a,1)
dnl	}
dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb__vector_mult_sum(const void * a, const void * b, void * c, rsb_type_t type, size_t n, const int inca, const int incb)
{
	/*!
	 * c <- sum(a*b)
	 * It is allowed to give c == a or c == b or a==b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see ddot in BLAS
	 *
	 * \return \rsberrcodemsg
	 * 
	 * p.s.: this routine is, numerically speaking, a crime!
	 * 
	 * */
	size_t i;
	if(a==b && inca==incb)
		return rsb_vector_norm_square_strided(c,a,type,n,inca);
	if(inca == 1 && incb == 1)
	{
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype*tb = b; const mtype*ta = a; mtype*tc = c,cacc = RSB_M4_ZERO(mtype);
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		cacc+=ta[i+LI]*tb[i+LI];
	'); 
		*tc = cacc;
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	else
	{
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype*tb = b; const mtype*ta = a; mtype*tc = c,cacc = RSB_M4_ZERO(mtype);
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		cacc+=ta[inca*(i+LI)]*tb[incb*(i+LI)];
	'); 
		*tc = cacc;
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_fill_with_zeros_nostride(void * array, rsb_type_t type, size_t n)
{
	/*!
	 * \ingroup gr_vec
	 * Will zero the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ){
	mtype*ta = array;
RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`ta[i+LI] = RSB_M4_ZERO(mtype);')}
	else 
#endif
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_fill_with_zeros(void * array, rsb_type_t type, size_t n, size_t incx)
{
	/*!
	 * \ingroup gr_vec
	 * Will zero the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(incx==1)
		return rsb_fill_with_zeros_nostride(array,type,n);

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ){
	mtype*ta = array;
RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`ta[(i+LI)*incx] = RSB_M4_ZERO(mtype);')}
	else 
#endif
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_vector_scale(void * a, const void * alphap, rsb_type_t type, size_t n)
{
	/*!
	 * a <- a * alpha
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param alphap scaling value (if NULL assumed to be zero)
	 * \param n	the input array length
	 * \note see dscal in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(alphap==NULL || RSB_IS_ELEMENT_ZERO(alphap,type))
		return rsb_fill_with_zeros(a,type,n,1);
		
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	const mtype alpha = *(mtype*)alphap; mtype*ta = a;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI]*=alpha;
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_strided_vector_scale(void * a, const void * alphap, rsb_type_t type, size_t n, size_t stride)
{
	/*!
	 * a <- a * alpha
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see dscal in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(stride==1)
		return rsb_vector_scale(a,alphap,type,n);
	if(alphap==NULL || RSB_IS_ELEMENT_ZERO(alphap,type))
		return rsb_fill_with_zeros(a,type,n,stride);

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype alpha = *(mtype*)alphap; mtype*ta = a;
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		ta[stride*(i+LI)]*=alpha;
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_vector_add(void * a, const void * alphap, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * a <- a + alpha
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
		
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype alpha = *(mtype*)alphap; mtype*ta = a;
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		ta[i+LI]+=alpha;
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_vector_div(void * a, const void * alphap, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * this is a benchmark-oriented function only..
	 * \return \rsberrcodemsg
	 * */
	size_t i;
		
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	const mtype alpha = *(mtype*)alphap; mtype*ta = a;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI]/=alpha;
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__vector_increase_by_one(void * a, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
		
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{ mtype*ta = a;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI]+=RSB_M4_ONE(mtype);
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_vector_pow(void * a, rsb_type_t type, const void *y, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(!a || !y)
		goto err;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		mtype ty = *(mtype*)y,*ta = a;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI] = RSB_M4_POW(mtype,ta[i+LI],ty);
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
err:
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_vector_sqrt(void * a, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(!a)goto err;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{mtype*ta = a;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI] = RSB_M4_SQRT(mtype,ta[i+LI]);
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
err:
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__vector_scale_inv(void * a, const void * alphap, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * a <- 1/a * alpha
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see dscal in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	if(!alphap)
		return RSB_ERR_BADARGS;
		
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		mtype alphai = RSB_M4_ONE(mtype)/(*(mtype*)alphap);
		return rsb_vector_scale(a,&alphai,type,n);
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__vector_sum_of_abs_diffs(void * c, const void * a, const void * b, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype*ap = a,*bp = b;
		mtype ac = RSB_M4_ZERO(mtype);
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		ac += RSB_M4_ABS(mtype,ap[i+LI]-bp[i+LI]);
		'); 
		*((mtype*)(c)) = ac;
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__vector_sum_of_abs(void * c, const void * a, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype*ap = a;
		mtype ac = RSB_M4_ZERO(mtype);
		RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
		ac += RSB_M4_ABS(mtype,ap[i+LI]);
		'); 
		*((mtype*)(c)) = ac;
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__vector_to_abs(void * a, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{mtype*ta = a;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI] = RSB_M4_ABS(mtype,ta[i+LI]);
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_alpha_sum(void * a, const void * b, const void * alphap, rsb_type_t type, size_t n)
{
	/*!
	 * a <- a + alpha * b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see daxpy in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	const mtype alpha = alphap ? *(mtype*)alphap : RSB_M4_ONE(mtype);
	mtype*ta = a; const mtype*tb = b;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI]+=alpha*tb[i+LI];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifelse(1,0,`dnl
/* redundant code (see rsb__cblas_Xaxpy) */
rsb_err_t rsb_vectors_sum(void * a, const void * b, rsb_type_t typecode, const void *alphap, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * a <- a + alpha * b
         *
	 * \param array	an array pointer
	 * \param typecode	a valid type code
	 * \param n	the input array length
	 * \note see daxpy in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if( !alphap || RSB_IS_ELEMENT_ONE(alphap,typecode))
	{
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{mtype*ta = a; const mtype*tb = b;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[i+LI]+=tb[i+LI];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	}
	else
		return rsb_alpha_sum(a,b,alphap,typecode,n);
dnl 	{
dnl foreach(`mtype',RSB_M4_TYPES,`dnl
dnl `#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
dnl 	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
dnl 	{const mtype alpha = *((const mtype*)alphap);
dnl 	mtype*ta = a; const mtype*tb = b;
dnl 	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
dnl 	ta[i+LI]+=alpha*tb[i+LI];
dnl 	'); 
dnl 	}
dnl 	else 
dnl #endif
dnl ')dnl
dnl 	return RSB_ERR_UNSUPPORTED_TYPE	;
dnl 	}

	return RSB_ERR_NO_ERROR;
}
')dnl
dnl
')dnl
dnl

dnl
rsb_err_t rsb__util_set_array_to_converted_integer(void *p, rsb_type_t typecode, const rsb_nnz_idx_t n, const rsb_nnz_idx_t incp, const rsb_int v)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * */
	size_t i;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	mtype*tp = p; const mtype tv = (mtype)v;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	tp[((i+LI)*incp)] = tv;
	'); 
	}
	else
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__vectors_left_sum_reduce_and_zero(void * d, void * s, const rsb_type_t typecode, const size_t n, const size_t incd, const size_t off)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * d[off:off+n-1] <- d[off:off+n-1] + s[off:off+n-1] 
	 * s[off:off+n-1] <- 0
         *
	 * \param array	an array pointer
	 * \param typecode	a valid type code
	 * \param incd	the stride of d
	 * \param off offset in the vectors
	 * \return \rsberrcodemsg
	 * */
	size_t i;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	mtype*td = d,*ts = s;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	td[(off+i+LI)*incd]+=ts[(off+i+LI)];
	ts[(off+i+LI)] = RSB_M4_ZERO(mtype);
	'); 
	}
	else
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
ifelse(1,0,`dnl
/* redundant code -- see rsb__cblas_Xaxpy */
rsb_err_t rsb_vectors_sum_scale_strided(void * a, const void * b, rsb_type_t typecode, const void *alphap, size_t n, size_t inca, size_t incb)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * a <- a + alpha * b
         *
	 * \param array	an array pointer
	 * \param typecode	a valid type code
	 * \param n	the input array length
	 * \note see daxpy in BLAS
	 * TODO: declare alpha as a const local variable, so the compiler will not contempt aliasing.
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(inca==1 && incb==1 /*&& ( !alphap || RSB_IS_ELEMENT_ONE(alphap,typecode))*/ )
		return rsb_vectors_sum(a,b,typecode,alphap,n);

	if( !alphap || RSB_IS_ELEMENT_ONE(alphap,typecode))
	{
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	mtype*ta = a; const mtype*tb = b;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[(i+LI)*inca]+=tb[(i+LI)*incb];
	'); 
	}
#endif
')dnl
	}
	else
	{
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{const mtype alpha = *((const mtype*)alphap);
	mtype*ta = a; const mtype*tb = b;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[(i+LI)*inca]+=alpha*tb[(i+LI)*incb];
	'); 
	}
#endif
')dnl
	}
dnl	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb_alpha_sum_strided(void * a, const void * b, const void * alphap, rsb_type_t type, size_t n, int inca, int incb)
{
	/*!
	 * a <- a + alpha * b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see daxpy in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(inca == 1 && incb == 1)
		return rsb_alpha_sum(a,b,alphap,type,n);
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	const mtype alpha = alphap ? *(mtype*)alphap : RSB_M4_ONE(mtype);
	mtype*ta = a; const mtype*tb = b;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ta[inca*(i+LI)]+=alpha*tb[incb*(i+LI)];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__cblas_Xaxpy(rsb_type_t type, size_t n, const void * alphap, const void * x, const int incx, void * y, const int incy)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * y <- y + alpha * x
         */
	return rsb_alpha_sum_strided(y,x,alphap,type,n,incy,incx);
}
')dnl
dnl

dnl
rsb_err_t rsb__vector_mult(const void * a, const void * b, void * c, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * c <- a*b
	 * It is allowed to give c == a or c == b or a == b
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * 
	 * FIXME : useless ?
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	const mtype*ta = a; const mtype*tb = b; mtype*tc = c;
dnl	//const mtype omega = *(mtype*)omegap;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	tc[i+LI] = ta[i+LI]*tb[i+LI];
	'); 
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__xcopy(void * a, const void * b, rsb_nnz_idx_t toi, rsb_nnz_idx_t foi, rsb_nnz_idx_t n,size_t el_size)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * a[toi:toi+n] <- b[foi:foi+n] 
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 * */
	rsb__memcpy(((rsb_byte_t*)a)+el_size*toi,((const rsb_byte_t*)b)+el_size*foi,el_size*n);
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_are_similar_parametric(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy, int extra_decimals)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 *
	 * For cases like 1+0I differing from 1-0I .
	 * */
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		rsb_nnz_idx_t i;
		const mtype *a = ap;
		const mtype *b = bp;
		const RSB_M4_REALT(mtype) threshold = RSB_M4_THRESHOLD_VALUE(mtype) * RSB_M4_POW(mtype,10*RSB_M4_ONE(mtype),(mtype)extra_decimals);

		for(i=0;i<n;++i)
		{
			const mtype av = a[incx*(i)];
			const mtype bv = b[incy*(i)];
			if( av - bv )
			{
				const mtype aav = RSB_M4_ABS(mtype,av);
				const mtype abv = RSB_M4_ABS(mtype,bv);
				if( av && RSB_M4_ABS(mtype,(aav - abv)) / RSB_M4_ABS(mtype,av) > threshold )
					goto differing;
				if( bv && RSB_M4_ABS(mtype,(aav - abv)) / RSB_M4_ABS(mtype,bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE;
differing:
	return RSB_ERR_GENERIC_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_are_similar(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 *
	 * For cases like 1+0I differing from 1-0I .
	 * */
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		rsb_nnz_idx_t i;
		const mtype *a = ap;
		const mtype *b = bp;
		const RSB_M4_REALT(mtype) threshold = RSB_M4_THRESHOLD_VALUE(mtype);

		for(i=0;i<n;++i)
		{
			const mtype av = a[incx*(i)];
			const mtype bv = b[incy*(i)];
			if( av - bv )
			{
				const mtype aav = RSB_M4_ABS(mtype,av);
				const mtype abv = RSB_M4_ABS(mtype,bv);
				if( av && RSB_M4_ABS(mtype,(aav - abv)) / RSB_M4_ABS(mtype,av) > threshold )
					goto differing;
				if( bv && RSB_M4_ABS(mtype,(aav - abv)) / RSB_M4_ABS(mtype,bv) > threshold )
					goto differing;
			}
		}
		return RSB_ERR_NO_ERROR;
	}
	else
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE;
differing:
	return RSB_ERR_GENERIC_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_are_same(const void * ap, const void * bp, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	return rsb__do_are_similar_parametric(ap, bp, n, typecode, incx, incy, 0);
}
')dnl
dnl

dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
static rsb_err_t rsb__xcopy_strided_typed(void * a, const void * b, rsb_nnz_idx_t toi, rsb_nnz_idx_t foi, rsb_nnz_idx_t n,rsb_type_t typecode, rsb_nnz_idx_t incx, rsb_nnz_idx_t incy)
{
	/*!
	 * a[toi:toi+n] <- b[foi:foi+n] 
         *
	 * \param array	an array pointer
	 * \param type	a valid type code
	 *
	 * \return \rsberrcodemsg
	 * */
	if(incx==1 && incy==1)
		return rsb__xcopy(a,b,toi,foi,n,RSB_SIZEOF(typecode));
	/* else */
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	rsb_nnz_idx_t i;
	mtype *ap = a; const mtype *bp = b;
	ap+=toi;
	bp+=foi;
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	ap[(i+LI)*incx] = bp[(i+LI)*incy];
	'); 
		return RSB_ERR_NO_ERROR;
	}
	else 
#endif
')dnl
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__sqrt_of_sum_of_fabs_diffs(const void * a, const void * b, void *err, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * Will compute the square root of the sum of the squares of the vectors elements differences.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * 
	 * FIXME
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
ifelse(mtype,`int',`dnl
	/* UHM ...  */
	{
	double acc = RSB_M4_ZERO(double);
	const mtype*ta = a, *tb = b;
	*((double*)err) = RSB_M4_ZERO(double);
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	acc+=(ta[i+LI]-tb[i+LI])*(ta[i+LI]-tb[i+LI]);
	'); 
	*((double*)err) = sqrt(acc);
	}
',`dnl
	{
	const mtype*ta = a; const mtype*tb = b;
	*((mtype*)err) = RSB_M4_ZERO(mtype);
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	*((mtype*)(err))+=(ta[i+LI]-tb[i+LI])*(ta[i+LI]-tb[i+LI]);
	'); 
	*((mtype*)err) = RSB_M4_SQRT(mtype,(*((mtype*)err)));
	}
')dnl
	else 
#endif
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__fill_with_increasing_values(void * array, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 * FIXME : document me
	 * starts with one.
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{ 
	mtype*ta = array;
RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`ta[i+LI] = (const mtype)(i+LI+1);')
	}
	else 
#endif
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_do_conjugate(void * array, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
ifelse(RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(mtype)),1,`dnl
	{
		mtype*ta = array;
RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`ta[i+LI] = RSB_M4_CONJ_SYM(mtype,`n',RSB_M4_SYMBOL_HERMITIAN)(ta[i+LI]);')
	}
',`dnl
		return RSB_ERR_NO_ERROR;
')dnl
	else 
`#endif'
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_do_negate(void * array, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 * Will negate the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{ 
	mtype*ta = array;
RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`ta[i+LI] = -ta[i+LI];')}
	else 
#endif
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_find_min(void * minp, const void * array, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(n<1)return RSB_ERR_BADARGS;
	if(inc<1)return RSB_ERR_BADARGS;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{const mtype * ap = array;mtype *mp = minp;
	*mp = *ap;for(i = 1;i<n;++i){if(RSB_M4_ABS(mtype,ap[i*inc])<RSB_M4_ABS(mtype,*mp) )*mp = ap[i*inc];
	}}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_find_max(void * maxp, const void * array, rsb_type_t type, size_t n, rsb_nnz_idx_t inc)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
	if(n<1)return RSB_ERR_BADARGS;
	if(inc<1)return RSB_ERR_BADARGS;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{const mtype * ap = array;mtype *mp = maxp;
	*mp = *ap;for(i=1;i<n;++i){if(RSB_M4_ABS(mtype,ap[i*inc])>RSB_M4_ABS(mtype,*mp))*mp = ap[i*inc];
	}}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_drop_to_zero_if_above_threshold(void * array, rsb_type_t type, size_t n, const void * threshold)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{const mtype th = (*(const mtype*)(threshold)); mtype*ta = array;
	for(i = 0;i<n;++i)
	{if(RSB_M4_ABS(mtype,th)<RSB_M4_ABS(mtype,ta[i]))ta[i] = RSB_M4_ZERO(mtype);}}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_nnz_idx_t rsb__util_count_positive(void * array, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i, c = 0;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{	mtype*ta = array;
		 for(i=0;i<n;++i)
			c+=(RSB_M4_CREAL(mtype,ta[i])>(RSB_M4_REALT(mtype))0);
	}else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return c;
}
')dnl
dnl

dnl
rsb_nnz_idx_t rsb__util_count_negative(void * array, rsb_type_t type, size_t n)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i, c = 0;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{	mtype*ta = array;
		 for(i=0;i<n;++i)
			c+=(RSB_M4_CREAL(mtype,ta[i])<(RSB_M4_REALT(mtype))0);
	}else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return c;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_drop_to_zero_if_under_threshold(void * array, rsb_type_t type, size_t n, const void * threshold)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ) {
	const mtype th = (*(mtype*)(threshold)); mtype*ta = ((mtype*)(array));
	for(i=0;i<n;++i){if(RSB_M4_ABS(mtype,th)>RSB_M4_ABS(mtype,ta[i]))ta[i] = RSB_M4_ZERO(mtype);}}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__fill_with_ones(void * array, rsb_type_t type, size_t n, size_t incx)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * \ingroup gr_vec
	 * Will set to one the input n elements long array of type type.
	 * \param array	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 *
	 * \return \rsberrcodemsg
	 * TODO:RENAME: rsb__fill_with_ones -> rsb__val_fill_with_ones.
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) ){
	mtype*ta = ((mtype*)(array));
 for(i=0;i<n;++i) {ta[i*incx] = RSB_M4_ONE(mtype);}}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__debug_print_vectors_diff_fd(const void * v1, const void * v2, size_t n, rsb_type_t type, size_t incx, size_t incy, int onlyfirst, FILE*fd)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*! 
	 * A debug function for printing the difference of two vectors of a specified type, in parallel.
	 * FIXME : It should take into account thresholds specific to each numerical type.
	 **/
#if RSB_ALLOW_STDOUT
	size_t i, differing = 0;
	if(!v1 || !v2)return RSB_ERR_BADARGS;

	/*RSB_STDERR("\t vectors diff :\n"); */
	RSB_FPRINTF(fd,"\t vectors diff :\n");
	
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype *v1p = v1,*v2p = v2;
		const RSB_M4_REALT(mtype) th = RSB_M4_THRESHOLD_VALUE(mtype);
		for(i=0;i<n ;++i) 
		ifelse(mtype,`long double complex',`if(creall(v1p[i*incx])-creall(v2p[i*incy])>th)/*FIXME : incomplete check*/',`dnl
		ifelse(mtype,`double complex',`if(creal(v1p[i*incx])-creal(v2p[i*incy])>th)/*FIXME : incomplete check*/',`dnl
		ifelse(mtype,`float complex',`if(crealf(v1p[i*incx])-crealf(v2p[i*incy])>th)/*FIXME : incomplete check*/',`dnl
		ifelse(mtype,`complex',       `if(creal(v1p[i*incx])-creal(v2p[i*incy])>th)/*FIXME : incomplete check*/',`dnl
		ifelse(mtype,`int',       `if(v1p[i*incx]-v2p[i*incy])',`dnl
		ifelse(mtype,`char',       `if(v1p[i*incx]-v2p[i*incy])',`dnl
if(fabs((double)(v1p[i*incx]-v2p[i*incy]))>th)/*FIXME : incomplete check*/
')dnl
')dnl
')dnl
')dnl
')dnl
')dnl
{		differing++;
		if((onlyfirst==0)||(onlyfirst>differing))
		RSB_FPRINTF(fd,"%zd : "RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)" "RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)"\n",(rsb_printf_int_t)i,dnl
		ifelse(mtype,`long double complex',`creall(v1p[i*incx]),cimagl(v1p[i*incx]),creall(v2p[i*incy]),cimagl(v2p[i*incy])',`dnl
		ifelse(mtype,`double complex',`creal(v1p[i*incx]),cimag(v1p[i*incx]),creal(v2p[i*incy]),cimag(v2p[i*incy])',`dnl
		ifelse(mtype,`float complex',`crealf(v1p[i*incx]),cimagf(v1p[i*incx]),crealf(v2p[i*incy]),cimagf(v2p[i*incy])',`dnl
		ifelse(mtype,`complex',`creal(v1p[i*incx]),cimag(v1p[i*incx]),creal(v2p[i*incy]),cimag(v2p[i*incy])',`dnl
v1p[i*incx],v2p[i*incy]`'dnl
')dnl
')dnl
')dnl
')dnl
		);
}
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	if(differing>onlyfirst)RSB_FPRINTF(fd,"...(for a total of %zd differing entries)...\n",(rsb_printf_int_t)(differing-onlyfirst));
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}
')dnl
dnl

dnl
rsb_err_t rsb__debug_print_vectors_diff(const void * v1, const void * v2, size_t n, rsb_type_t type, size_t incx, size_t incy, int onlyfirst)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	return rsb__debug_print_vectors_diff_fd(v1, v2, n, type, incx, incy, onlyfirst, RSB_DEFAULT_STREAM);
}
')dnl
dnl

dnl
rsb_err_t rsb__debug_print_value(const void * v, rsb_type_t type)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*! 
	 **/
#if RSB_ALLOW_STDOUT
	if(!v)return RSB_ERR_BADARGS;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype *v1p = v;
		RSB_STDOUT(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype),dnl
		ifelse(mtype,`long double complex',`creall(v1p[0]),cimagl(v1p[0])',`dnl
		ifelse(mtype,`double complex',`creal(v1p[0]),cimag(v1p[0])',`dnl
		ifelse(mtype,`float complex',`crealf(v1p[0]),cimagf(v1p[0])',`dnl
		ifelse(mtype,`complex',`creal(v1p[0]),cimag(v1p[0])',`dnl
v1p[0]`'dnl
')dnl
')dnl
')dnl
')dnl
		);
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}
')dnl
dnl

dnl
rsb_err_t rsb__debug_print_vector_extra(const void * v1, size_t n, rsb_type_t type, size_t inc, int style, FILE*stream)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*! 
	 * A debug function for printing two vectors of a specified type, in parallel.
	 **/
#if RSB_ALLOW_STDOUT
	rsb_nnz_idx_t i;
	int want_header = ( style == 0x1 );
	const char * ts = RSB_IS_MATRIX_TYPE_COMPLEX(type)?"complex":"real";
	const char * ss = RSB_SYMMETRY_STRING(RSB_FLAG_NOFLAGS);
	
	if( n <= 0 )
		goto errb;

	if(!v1 || !stream)
		goto errb;

	/*if(!want_header)
		RSB_STDERR("\t vectors  :\n");*/
	
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype *v1p = v1;
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix array %s %s\n%zd %zd\n",ts,ss,(rsb_printf_int_t)n,(rsb_printf_int_t)1);
		for(i=0;i<n;++i) 
		RSB_FPRINTF(stream,RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype) "\n",dnl
		ifelse(mtype,`long double complex',`creall(v1p[i*inc]),cimagl(v1p[i*inc])',`dnl
		ifelse(mtype,`double complex',`creal(v1p[i*inc]),cimag(v1p[i*inc])',`dnl
		ifelse(mtype,`float complex',`crealf(v1p[i*inc]),cimagf(v1p[i*inc])',`dnl
		ifelse(mtype,`complex',`creal(v1p[i*inc]),cimag(v1p[i*inc])',`dnl
v1p[i*inc]`'dnl
')dnl
')dnl
')dnl
')dnl
		);
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
errb:
	return RSB_ERR_BADARGS;
}
')dnl
dnl

dnl
rsb_err_t rsb__debug_print_vector(const void * v1, size_t n, rsb_type_t type, size_t inc)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	return rsb__debug_print_vector_extra(v1, n, type, inc, 0x0, stdout);
}
')dnl
dnl

dnl
rsb_err_t rsb__debug_print_vectors(const void * v1, const void * v2, size_t n, size_t incx, size_t incy, rsb_type_t type)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*! 
	 * A debug function for printing two vectors of a specified type, in parallel.
	 **/
#if RSB_ALLOW_STDOUT
	size_t i;
	if(!v1 || !v2)return RSB_ERR_BADARGS;

	RSB_STDERR("\t vectors  :\n");
	
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		const mtype *v1p = v1,*v2p = v2;
		for(i=0;i<n;++i) 
		RSB_STDOUT(RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)" "RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)"\n",dnl
RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG(mtype,`v1p[(i)*incx]'),dnl
RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG(mtype,`v2p[(i)*incy]')dnl
);
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}
')dnl
dnl
dnl
')
dnl
dnl
dnl
dnl ifdef(`ONLY_WANT_HEADERS',`dnl
dnl #ifndef RSB_UTIL_H_INCLUDED
dnl #define RSB_UTIL_H_INCLUDED
dnl ')
dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
dnl
dnl
')dnl
dnl
dnl
dnl

#ifdef RSB_WANT_OSKI_BENCHMARKING 
rsb_err_t rsb__do_account_sorted_optimized_css(
	 const rsb_coo_idx_t * MIndx, const rsb_coo_idx_t * mIndx,
	 const rsb_coo_idx_t Mdim, const rsb_coo_idx_t mdim,
	 const rsb_nnz_idx_t nnz, rsb_nnz_idx_t * elements_per_block_row, rsb_nnz_idx_t * blocks_per_block_row
)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	 	\ingroup gr_internals

		elements_per_block_row and blocks_per_block_row arrays should be blanked.
		FIXME : missing error handling.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t n = 0;

	if(blocks_per_block_row)
	for(n=0;n<nnz;++n)
	{
		RSB_DEBUG_ASSERT(MIndx[n]<Mdim);
		RSB_DEBUG_ASSERT(MIndx[n]>=0);
		elements_per_block_row[MIndx[n]]++;
		blocks_per_block_row  [MIndx[n]]++;
	}
	else
	for(n=0;n<nnz;++n)
	{
		RSB_DEBUG_ASSERT(MIndx[n]<Mdim);
		RSB_DEBUG_ASSERT(MIndx[n]>=0);
		elements_per_block_row[MIndx[n]]++;
	}
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl

#endif /* RSB_WANT_OSKI_BENCHMARKING */

dnl
#if RSB_OBSOLETE_QUARANTINE_UNUSED

rsb_err_t rsb__do_account_sorted_optimized(
	 struct rsb_mtx_t * mtxAp,
	 const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA,
	 const rsb_coo_idx_t Idim, const rsb_coo_idx_t Jdim,
	 const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop,
rsb_nnz_idx_t * elements_per_block_row, 
rsb_nnz_idx_t * blocks_per_block_row
)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	 *	\ingroup gr_internals
	 * 	FIXME : document this
	 */
	rsb_coo_idx_t blockrows = 0;
	rsb_coo_idx_t blockcolumns = 0;
	rsb_coo_idx_t baserow = 0;
	rsb_coo_idx_t basecolumn = 0;
dnl	rsb_nnz_idx_t block_count = 0;
dnl	rsb_nnz_idx_t *indptr = mtxAp->indptr;
dnl	rsb_nnz_idx_t *bpntr = mtxAp->bpntr;
dnl	rsb_nnz_idx_t *bindx = mtxAp->bindx;
	const rsb_coo_idx_t *Mpntr = NULL;
	const rsb_coo_idx_t *mpntr = NULL;
	const rsb_coo_idx_t *MIndx = NULL;
	const rsb_coo_idx_t *mIndx = NULL;
	rsb_blk_idx_t mI = 0, MI = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;	/* will index a nnz sized array */
	int K = 0;
	
dnl	if(0)
dnl	//if( flags & RSB_FLAG_SHOULD_DEBUG )
dnl		errval = rsb__do_account_sorted( mtxAp, IA, JA, nnz, pinfop, elements_per_block_row, blocks_per_block_row);

	if(nnz==0)
	{
		/* FIXME: new case, incomplete (useful for diagonal implicit matrices) */
		return RSB_ERR_NO_ERROR;
	}

#if RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS
	if(!pinfop)
	{
		/* a performance fix */
		if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			return rsb__do_account_sorted_optimized_css(JA,IA,Jdim,Idim,nnz,elements_per_block_row,blocks_per_block_row);
		else
			return rsb__do_account_sorted_optimized_css(IA,JA,Idim,Jdim,nnz,elements_per_block_row,blocks_per_block_row);
	}
#endif
	
	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	{
		mpntr = pinfop->rpntr;
		Mpntr = pinfop->cpntr;
		mIndx = IA;
		MIndx = JA;
	}
	else
	{
		Mpntr = pinfop->rpntr;
		mpntr = pinfop->cpntr;
		MIndx = IA;
		mIndx = JA;
	}

foreach(`matrix_storage',RSB_M4_MATRIX_STORAGE,`dnl
	/*	storage matrix_storage	*/
	if( mtxAp->`matrix_storage'==RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(matrix_storage) )
{
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
pushdef(`RSB_M4_Mpntr',`(blockrows   *($1))')dnl
pushdef(`RSB_M4_mpntr',`(blockcolumns*($1))')dnl
',`dnl
pushdef(`RSB_M4_Mpntr',`Mpntr[$1]')dnl
pushdef(`RSB_M4_mpntr',`mpntr[$1]')dnl
')dnl


	k = mI = MI = K=0;
	while( MIndx[k] >= RSB_M4_Mpntr(MI+1) )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= RSB_M4_mpntr(mI+1) )++mI;	/* skipping preceding block columns .. */
	blockrows    = RSB_M4_Mpntr(MI+1) - RSB_M4_Mpntr(MI);
	blockcolumns = RSB_M4_mpntr(mI+1) - RSB_M4_mpntr(mI);
	baserow = RSB_M4_Mpntr(MI);
	basecolumn = RSB_M4_mpntr(mI);
	*elements_per_block_row = 0;
	*blocks_per_block_row   = 0;	
	elements_per_block_row[MI*0] += blockrows * blockcolumns;
	blocks_per_block_row[MI]   +=1;

	while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k,(rsb_printf_int_t) (MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
			mI = mIndx[k]/blockcolumns;
',`
			while( mIndx[k] >= RSB_M4_mpntr(mI+1) )++mI;
			blockcolumns = RSB_M4_mpntr(mI+1) - RSB_M4_mpntr(mI);
')dnl
			basecolumn = RSB_M4_mpntr(mI);

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
				MI = MIndx[k]/blockrows;
',`
				while( MIndx[k] >= RSB_M4_Mpntr(MI+1) )++MI;
				blockrows    = RSB_M4_Mpntr(MI+1) - RSB_M4_Mpntr(MI);
')dnl
				baserow = RSB_M4_Mpntr(MI);
			}
			else
			{
				/* same block row  */
			}
			*elements_per_block_row += blockrows * blockcolumns;
			blocks_per_block_row[MI]   +=1;
			++K;
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
			MI = MIndx[k]/blockrows;
',`
			while( MIndx[k] >= RSB_M4_Mpntr(MI+1) )++MI;
			blockrows    = RSB_M4_Mpntr(MI+1) - RSB_M4_Mpntr(MI);
')dnl
			baserow = RSB_M4_Mpntr(MI);

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
				mI = mIndx[k]/blockcolumns;
',`
				mI = 0;
				while( mIndx[k] >= RSB_M4_mpntr(mI+1) )++mI;
				blockcolumns = RSB_M4_mpntr(mI+1) - RSB_M4_mpntr(mI);
')dnl
				basecolumn = RSB_M4_mpntr(mI);
			}
			else
			{
				/* new row block, same column  */
			}
			/* get rid of this var : elements_per_block_row */
			*elements_per_block_row += blockrows * blockcolumns;
			blocks_per_block_row[MI]   +=1;
			++K;
		}
		else
		{
			/* same block row for sure */
		}
		++k;
	}
	errval = RSB_ERR_NO_ERROR;goto ret;
	}
popdef(`RSB_M4_Mpntr')dnl
popdef(`RSB_M4_mpntr')dnl
')dnl
dnl
	errval = RSB_ERR_INTERNAL_ERROR;
ret:	return errval;
}
dnl
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
dnl
static rsb_err_t rsb__do_insert_sorted_optimized_css( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * MIndx, const rsb_coo_idx_t * mIndx, const rsb_nnz_idx_t nnz)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	 	\ingroup gr_internals

		elements_per_block_row and blocks_per_block_row arrays should be blanked.
		FIXME : missing error handling.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t n = 0;

	/* in case of RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR, they are equal */
	if(mtxAp->VA != VA)
		rsb__memcpy(mtxAp->VA  ,VA  ,mtxAp->el_size*nnz);

	for(n=0;n<nnz+1;++n)
		mtxAp->indptr[n] = n;

	for(n=0;n<mtxAp->nnz;++n)
		mtxAp->bindx [n] = mIndx[n];
	mtxAp->bindx [nnz] = 0;

	// should also set bindx, indptr, 
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
dnl

dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
dnl
rsb_err_t rsb__do_insert_sorted_optimized( struct rsb_mtx_t * mtxAp, const void *VA, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const struct rsb_mtx_partitioning_info_t * pinfop)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*
	 *	FIXME ! UNFINISHED 
	 * 	and please note that linked format is incomplete, so it does not support well block column major
	 */
	rsb_coo_idx_t blockrows = 0;
	rsb_coo_idx_t blockcolumns = 0;
	rsb_coo_idx_t baserow = 0;
	rsb_coo_idx_t basecolumn = 0;
dnl	rsb_nnz_idx_t block_count = 0;
	rsb_nnz_idx_t *indptr = mtxAp->indptr;
dnl	rsb_nnz_idx_t *bpntr = mtxAp->bpntr;
	rsb_nnz_idx_t *bindx = mtxAp->bindx;
	const rsb_coo_idx_t *Mpntr = NULL;
	const rsb_coo_idx_t *mpntr = NULL;
	const rsb_coo_idx_t *MIndx = NULL;
	const rsb_coo_idx_t *mIndx = NULL;
	rsb_blk_idx_t mI = 0, MI = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;	/* will index a nnz sized array */
	rsb_nnz_idx_t K = 0;

	if(nnz==0)
	{
		/* FIXME: new case, incomplete (useful for diagonal implicit matrices) */
		K = 0;		/* if nnz == 0 then K == 0 */
		bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
		return RSB_ERR_NO_ERROR;
	}

dnl	if(0)
dnl		return rsb__do_insert_sorted( mtxAp, VA, IA, JA, nnz, pinfop);

#if RSB_WANT_EXPERIMENTAL_NO_EXTRA_CSR_ALLOCATIONS
	if(!pinfop)
	{
		/* a performance fix */
		if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			return rsb__do_insert_sorted_optimized_css( mtxAp, VA, JA, IA, nnz );
		else
			return rsb__do_insert_sorted_optimized_css( mtxAp, VA, IA, JA, nnz );
	}
#endif

	if(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
	{
		mpntr = pinfop->rpntr;
		Mpntr = pinfop->cpntr;
		mIndx = IA;
		MIndx = JA;
	}
	else
	{
		Mpntr = pinfop->rpntr;
		mpntr = pinfop->cpntr;
		MIndx = IA;
		mIndx = JA;
	}


foreach(`matrix_storage',RSB_M4_MATRIX_STORAGE,`dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
	/*	type mtype, storage matrix_storage	*/
	if( mtxAp->typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	if( mtxAp->`matrix_storage'==RSB_M4_MATRIX_STORAGE_PREPROCESSOR_SYMBOL(matrix_storage) )
{
	mtype * dst = mtxAp->VA;
	k = mI = MI = 0;K = 0;
#if RSB_EXPERIMENTAL_USE_PURE_BCSS_FOR_CONSTRUCTOR
	rsb__get_physical_blocking_size(mtxAp, &blockrows, &blockcolumns);
	RSB_ASSERT( blockrows && blockcolumns);
#else
	blockrows    = Mpntr[MI+1] - Mpntr[MI];
	blockcolumns = mpntr[mI+1] - mpntr[mI];
#endif
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
pushdef(`RSB_M4_Mpntr',`(blockrows   *($1))')dnl
pushdef(`RSB_M4_mpntr',`(blockcolumns*($1))')dnl
',`dnl
pushdef(`RSB_M4_Mpntr',`Mpntr[$1]')dnl
pushdef(`RSB_M4_mpntr',`mpntr[$1]')dnl
')dnl
ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
pushdef(`RSB_M4_IBO',`(MIndx[k]-baserow)+(mIndx[k]-basecolumn)*blockrows')dnl
',`dnl
pushdef(`RSB_M4_IBO',`(MIndx[k]-baserow)*blockcolumns+(mIndx[k]-basecolumn)')dnl
')dnl
pushdef(`RSB_M4_BLOCK_OFFSET',`
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
	K * blockrows * blockcolumns
',`dnl
ifelse(RSB_M4_IS_FORMAT_LINKED_LIST(matrix_storage),1,`dnl
	RSB_BLOCK_OFFSET(mtxAp,K)
dnl	indptr[ K ] seems not adequate
',`dnl
	indptr[ K ]
	/*K * blockrows * blockcolumns*/
	/*RSB_BLOCK_OFFSET(mtxAp,K)/mtxAp->el_size*/ /* FIXME : unfinished ! */ 
')dnl
')dnl
dnl
')dnl

	while( MIndx[k] >= RSB_M4_Mpntr(MI+1) )++MI;	/* skipping preceding block rows .. */
	while( mIndx[k] >= RSB_M4_mpntr(mI+1) )++mI;	/* skipping preceding block columns .. */
dnl	blockrows    = RSB_M4_Mpntr(MI+1) - RSB_M4_Mpntr(MI);
dnl	blockcolumns = RSB_M4_mpntr(mI+1) - RSB_M4_mpntr(mI);
	baserow = RSB_M4_Mpntr(MI);
	basecolumn = RSB_M4_mpntr(mI);
	bindx [ K ] = mI;			/* a new block */
	indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;	/* FIXME : DUPLICATION ?! see later */

dnl 	DELETE THIS
ifelse(RSB_M4_IS_FORMAT_LINKED_LIST(matrix_storage),1,`dnl
	{
		if(RSB_WANT_VERBOSE_MESSAGES)
			RSB_INFO("initializing linked lists stuff.\n");
ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),mI,MI,blockcolumns,blockrows,basecolumn,baserow)
',`dnl
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),MI,mI,blockrows,blockcolumns,baserow,basecolumn)
')dnl
	}
')dnl

dnl	dst =  mtxAp->VA;
dnl	dst += RSB_M4_IBO;
dnl	dst += RSB_M4_BLOCK_OFFSET;
dnl	{rsb_blk_idx_t ibo = 0;/* FIXME */
dnl ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
dnl 		ibo = RSB_GET_INTRA_BLOCK_OFFSET(mIndx[k],MIndx[k],mI,MI,mtxAp) ;
dnl ',`dnl
dnl 		ibo = RSB_GET_INTRA_BLOCK_OFFSET(MIndx[k],mIndx[k],MI,mI,mtxAp) ;
dnl ')dnl
dnl 		dst += ibo;
dnl	}

	if( (mtxAp->flags & RSB_FLAG_SORTED_INPUT ) != 0 && 1 /* ONLY FOR 1 X 1 BLOCKED */)
	{
		//RSB_STDERR("rsb__do_insert_sorted_optimized : TODO : please specialize for specific blockings ! \n");
	}

while(RSB_LIKELY(k<nnz))
	{
#ifdef DEBUG
		if( MIndx[k] < baserow  )
		{
			RSB_ERROR("k=%zd : (%zd %zd) is not ok\n",k, (rsb_printf_int_t)(MIndx[k]+1),(rsb_printf_int_t)(mIndx[k]+1));
			RSB_STDERR("(minor dim. index %zd < base row %zd)\n",(rsb_printf_int_t)MIndx[k] , (rsb_printf_int_t)baserow);
			errval = RSB_ERR_INTERNAL_ERROR;
			goto err;/* NOTE : this jump could be evil */
		}
#endif

		if( mIndx[k] >= basecolumn+blockcolumns  )
		{
			/* new block column, for sure */
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
			mI = mIndx[k]/blockcolumns;
',`
			while( mIndx[k] >= RSB_M4_mpntr(mI+1) )++mI;
			blockcolumns = RSB_M4_mpntr(mI+1) - RSB_M4_mpntr(mI);
')dnl
			basecolumn = RSB_M4_mpntr(mI);

			if( MIndx[k] >= baserow+blockrows  )
			{
				/* new block row  */
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
				MI = MIndx[k]/blockrows;
',`
				while( MIndx[k] >= RSB_M4_Mpntr(MI+1) )++MI;
				blockrows    = RSB_M4_Mpntr(MI+1) - RSB_M4_Mpntr(MI);
')dnl
				baserow = RSB_M4_Mpntr(MI);
			}
			else
			{
				/* same block row  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
ifelse(RSB_M4_IS_FORMAT_LINKED_LIST(matrix_storage),1,`dnl
ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),mI,MI,blockcolumns,blockrows,basecolumn,baserow)
',`dnl
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),MI,mI,blockrows,blockcolumns,baserow,basecolumn)
')dnl
')dnl
		}
		else
		if( MIndx[k] >= baserow+blockrows  )
		{
			/* new row block, for sure */
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
				MI = MIndx[k]/blockrows;
',`
				while( MIndx[k] >= RSB_M4_Mpntr(MI+1) )++MI;
				blockrows    = RSB_M4_Mpntr(MI+1) - RSB_M4_Mpntr(MI);
')dnl
				baserow = RSB_M4_Mpntr(MI);

			if( mIndx[k] < basecolumn  )
			{
				/* new row block, new block column  */
				mI = 0;
ifelse(RSB_M4_IS_FORMAT_BCSS(matrix_storage),1,`dnl
				mI = mIndx[k]/blockcolumns;
',`
				while( mIndx[k] >= RSB_M4_mpntr(mI+1) )++mI;
				blockcolumns = RSB_M4_mpntr(mI+1) - RSB_M4_mpntr(mI);
')dnl
				basecolumn = RSB_M4_mpntr(mI);
			}
			else
			{
				/* new row block, same column  */
			}
			++K;
			bindx [ K ] = mI;			/* a new block */
			indptr[ K+1 ] = indptr[ K  ] + blockrows * blockcolumns;
ifelse(RSB_M4_IS_FORMAT_LINKED_LIST(matrix_storage),1,`dnl
ifelse(RSB_M4_IS_FORMAT_COLUMN_MAJOR(matrix_storage),1,`dnl
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),mI,MI,blockcolumns,blockrows,basecolumn,baserow)
',`
	RSB_BLOCK_TRAILING_STRUCT_SET(RSB_BLOCK_TRAILING_STRUCT_GET(mtxAp,K),MI,mI,blockrows,blockcolumns,baserow,basecolumn)
')dnl
')dnl
		}
		else
		{
			/* same block row for sure */
		}
		dst =  mtxAp->VA;
                RSB_DEBUG_ASSERT(mI>=0);
                RSB_DEBUG_ASSERT(MI>=0);

ifelse(RSB_M4_IS_FORMAT_LINKED_LIST(matrix_storage),1,`dnl
		/* :( */
		dst = (mtype*) ( ((rsb_byte_t*)dst)+RSB_BLOCK_OFFSET(mtxAp,K) );
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
dnl		dst = (mtype*)((rsb_byte_t*)dst)+(K+1)*RSB_BLOCK_EXTRA_BYTES;
		rsb_nnz_idx_t ibo = 0;
		if(RSB_UNLIKELY(mtxAp->flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
			ibo = RSB_GET_INTRA_BLOCK_OFFSET(mIndx[k],MIndx[k],mI,MI,mtxAp) ;
		else
			ibo = RSB_GET_INTRA_BLOCK_OFFSET(MIndx[k],mIndx[k],MI,mI,mtxAp) ;
		dst = (mtype*) (((rsb_byte_t*)dst)+ibo);
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
',`
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += RSB_M4_BLOCK_OFFSET;
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst += RSB_M4_IBO;
')dnl
		RSB_DEBUG_ASSERT(((rsb_byte_t*)dst)>=((rsb_byte_t*)mtxAp->VA));
		dst[0] = ((const mtype*)VA)[k];
		++k;
	}
	if(nnz)++K;	/* if nnz == 0 then K = 0 */
	bindx[K] = 0;	// the first element off the working bindx should be set to a safe value
	return RSB_ERR_NO_ERROR;	/* FIXME ! */
}
popdef(`RSB_M4_Mpntr')dnl
popdef(`RSB_M4_mpntr')dnl
popdef(`RSB_M4_IBO')dnl
popdef(`RSB_M4_BLOCK_OFFSET')dnl
dnl
')dnl
')dnl
	errval = RSB_ERR_INTERNAL_ERROR;
	return errval;
}
dnl
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_block(rsb_type_t type, const void * VA, rsb_blk_idx_t roff, rsb_blk_idx_t coff, rsb_blk_idx_t rows, rsb_blk_idx_t cols )
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * Will dump to stdout a dense matrix.
	 * Used for debugging purposes.
	 *
	 * FIXME : should be integrated with the macro subsystem in util.m4, and support column major order, and debugged.
	 */
#if RSB_ALLOW_STDOUT
	register rsb_coo_idx_t i, j;

	if(RSB_BLK_MUL_OVERFLOW(rows,cols))
		return RSB_ERR_LIMITS;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if(type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))
	{
		for(i=0;i<rows;++i)for(j=0;j<cols;++j)
		if(((mtype*)VA)[cols*i+j]!=RSB_M4_ZERO(mtype) )
		{ RSB_STDOUT(""
dnl :( were %10 %10 % 20 ...
		"%zd"/* FIXME : this could be any index type! */
		"\t"
		"%zd"
		"\t"
		RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)
		"\n",(rsb_printf_int_t)(roff+i+1),(rsb_printf_int_t)(coff+j+1),
RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG(mtype,`((mtype*)VA)[cols*i+j]'));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
dnl

dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__dump_blocks(const struct rsb_mtx_t *mtxAp)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	return RSB_ERR_UNIMPLEMENTED_YET;
#if 0
	/*! 
	 * \ingroup gr_internals
	 * A debug function for printing out the matrix structure.
	 *
	 * FIXME : UNFINISHED
	 * Note : it is extremely slow.
	 **/
	rsb_blk_idx_t i,j;
	if(!mtxAp)return RSB_ERR_BADARGS;
	if(!mtxAp->options)return RSB_ERR_BADARGS;

	RSB_STDERR("\t block structure :\n");
	
	/* this prints out the matrix blocks nnz structure */
	for(i=0;i<mtxAp->M_b;++i)
	{
		for(j=0;j<mtxAp->K_b;++j)
		if((RSB_BITMAP_GET(mtxAp->options->bitmap,mtxAp->M_b,mtxAp->K_b,i,j)))
		{
			RSB_STDERR("1");
		}
		else
		{
			RSB_STDERR("0");
		}
		RSB_STDERR("\n");
	}
	return RSB_ERR_NO_ERROR;
#endif
}
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
dnl

dnl
rsb_err_t rsb__test_print_csr(rsb_type_t type, rsb_flags_t flags, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t want_header, FILE*stream)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	 * \ingroup gr_internals
	 * Dumps out a whole matrix, from its CSR representation.
	 * 
	 * Warning : the nonzeros should be sorted on input.
	 */
#if RSB_ALLOW_STDOUT
	rsb_coo_idx_t k;
	if( !stream )goto err;
	if( !IA )goto err;
	if( ( !JA || !VA ) && nnz>0  )goto err;

	RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t)rows);
	/* RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t) nnz); */
	for(k=0;k<rows+1;++k) { RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t)(IA[k]+1)); }
	for(k=0;k<nnz   ;++k) { RSB_FPRINTF(stream,"%zd\n",(rsb_printf_int_t)(JA[k]+1)); }
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if(type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))
	{
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)
				"\n"
				,RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG(mtype,`((mtype*)VA)[k]'));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
')dnl
err:
	return RSB_ERR_GENERIC_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}
')dnl
dnl

dnl
rsb_err_t rsb__test_print_coo_mm(rsb_type_t type, rsb_flags_t flags, const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const void * VA, rsb_coo_idx_t rows, rsb_coo_idx_t cols, rsb_nnz_idx_t nnz, rsb_bool_t want_header, FILE*stream)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	 * \ingroup gr_internals
	 * Dumps out a whole matrix, from its coordinates, in matrix market format.
	 * 
	 * Warning : the nonzeros should be sorted on input.
	 */
#if RSB_ALLOW_STDOUT
	rsb_coo_idx_t k;
	const char * ts = RSB_IS_MATRIX_TYPE_COMPLEX(type)?"complex":"real";
	const char * ss = RSB_SYMMETRY_STRING(flags);
	
	if( !stream )
	{
		goto err;
	}

	if( ( !IA || !JA || !VA ) && nnz > 0 )
		goto err;
	if( rows < 0 || cols < 0 || nnz < 0 )
		goto err;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if(type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))
	{
		if(want_header)RSB_FPRINTF(stream,"%%%%MatrixMarket matrix coordinate %s %s\n%zd %zd %zd\n",ts,ss,(rsb_printf_int_t)rows,(rsb_printf_int_t)cols,(rsb_printf_int_t)nnz);
/*		for(k=0;k<nnz;++k) { RSB_FPRINTF(stream,"%6zd %6zd %20g\n",(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),((float*)VA)[k]); }*/
		for(k=0;k<nnz;++k)
		{
			RSB_FPRINTF(stream,
				"%zd"
				"\t"
				"%zd"
				"\t"
				RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_STRING(mtype)
				"\n"
				,(rsb_printf_int_t)(IA[k]+1),(rsb_printf_int_t)(JA[k]+1),dnl
RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_PRINTF_ARG(mtype,`((mtype*)VA)[k]'));
		}
		return RSB_ERR_NO_ERROR;
	}
#endif
')dnl
err:
	return RSB_ERR_GENERIC_ERROR;
#else
	return RSB_ERR_UNSUPPORTED_FEATURE; 
#endif
}
')dnl
dnl

dnl
/*static*/ /*inline*/ size_t rsb__do_sizeof(rsb_type_t type)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
	{
		/*
		 * FIXME : UNUSED ?
		 */
		size_t so = 0;
		switch(type)
		{
			/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
			case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:
				so = sizeof(type);
			break;
')dnl
			/* unsupported type */
			default :
			RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS 
		}
		return so;
	}
')dnl
dnl

dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
dnl
ifdef(`1',`0',`dnl
rsb_err_t rsb__do_coo_sum( struct rsb_coo_mtx_t*coocp, const void *alphap, const struct rsb_coo_mtx_t*cooap, const void *betap,  const struct rsb_coo_mtx_t*coobp)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	struct rsb_coo_mtx_t cooa = *cooap, coob = *coobp, cooc = *coocp;
	rsb_nnz_idx_t /*rnz = 0,*/an, bn, cn;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if(cooa.typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))
	{
	mtype alpha = alphap?*(mtype*)alphap:RSB_M4_ONE(mtype);
	mtype beta  = betap ?*(mtype*)betap :RSB_M4_ONE(mtype);
	for(cn = 0, an = 0, bn = 0;an<cooa.nnz || bn<coob.nnz;)
	{
		rsb_nnz_idx_t ap = an, bp = bn;
		if(cooa.IA[an]==coob.IA[bn] && cooa.JA[an]==coob.JA[bn])
			cooc.IA[cn] = cooa.IA[an],cooc.JA[cn] = cooa.JA[an],
			((mtype*)cooc.VA)[cn] = alpha * ((mtype*)cooa.VA)[an] + beta * ((mtype*)coob.VA)[bn],
			ap = an, bp = bn, ++cn, ++an, ++bn;

		for(;an<cooa.nnz && cooa.IA[an]==cooa.IA[ap] && cooa.JA[an]==cooa.JA[ap] ;++an)
			//RSB_STDOUT("x> %d %d\n",cooa.IA[an],cooa.JA[an])
			((mtype*)cooc.VA)[cn] += alpha * ((mtype*)cooa.VA)[an];

		for(;bn<coob.nnz && coob.IA[bn]==coob.IA[bp] && coob.JA[bn]==coob.JA[bp] ;++bn)
			//RSB_STDOUT("x> %d %d\n",coob.IA[bn],coob.JA[bn])
			((mtype*)cooc.VA)[cn] += beta  * ((mtype*)coob.VA)[bn];

		if( bn<coob.nnz )
		for(;an<cooa.nnz && (cooa.IA[an]<coob.IA[bn] ||
			       	(cooa.IA[an] <= coob.IA[bn] && cooa.JA[an]<coob.JA[bn]))
			       	;++an)
				//RSB_STDOUT("-> %d %d\n",cooa.IA[an],cooa.JA[an]),
			cooc.IA[cn] = cooa.IA[an], cooc.JA[cn] = cooa.JA[an],
			((mtype*)cooc.VA)[cn] = alpha * ((mtype*)cooa.VA)[an],
			++cn;

		if( an<cooa.nnz )
		for(;bn<coob.nnz && (cooa.IA[an]>coob.IA[bn] ||
			       	(cooa.IA[an]>=coob.IA[bn] && cooa.JA[an]>coob.JA[bn]))
			       	;++bn)
			//	RSB_STDOUT("-> %d %d\n",coob.IA[bn],coob.JA[bn]),
			cooc.IA[cn] = coob.IA[bn],cooc.JA[cn] = coob.JA[bn],
			((mtype*)cooc.VA)[cn] = beta * ((mtype*)coob.VA)[bn],
			++cn;
		//RSB_STDOUT("? %d %d\n",an,bn);
	}
	}
	else 
#endif
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
dnl

dnl
dnl
rsb_err_t rsb__cor_merge_dups(rsb_type_t typecode, void* RSB_RESTRICT VA, rsb_coo_idx_t * RSB_RESTRICT IA, rsb_coo_idx_t * RSB_RESTRICT JA, rsb_nnz_idx_t offB, rsb_nnz_idx_t nnzB, rsb_nnz_idx_t nnzC, const int wv, int wp, rsb_nnz_idx_t *onzp, struct rsb_coo_mtx_t*RSB_RESTRICT coop)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
		See rsb__cor_merge.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VB = NULL, *VC = NULL, *VT = NULL;
	rsb_coo_idx_t * IB = NULL, *JB = NULL;
	rsb_coo_idx_t * IC = NULL, *JC = NULL;
	rsb_coo_idx_t * IT = NULL, *JT = NULL;
	rsb_nnz_idx_t bi = 0, ci = 0, ti = 0;
	rsb_nnz_idx_t b0 = 0, c0 = 0;
	rsb_nnz_idx_t onz = 0;
	struct rsb_coo_mtx_t coo;
	size_t es = RSB_SIZEOF(typecode);

	if( nnzB == 0 || nnzC == 0 )
	{
		goto ret;
	}

	b0 = offB;
	c0 = offB + nnzB;
	VB = RSB_TYPED_OFF_PTR(typecode,VA,b0);
	VC = RSB_TYPED_OFF_PTR(typecode,VA,c0);
	IB = IA + b0;
	IC = IA + c0;
	JB = JA + b0;
	JC = JA + c0;

	RSB_BZERO_P(&coo);
	coo.nnz = nnzB + nnzC;
	coo.typecode = typecode;

	if( coop && coop->nnz)
	{
		coo = *coop;
		coo.nnz = nnzB + nnzC; /* necessary */
	}
	else
	{
		if( NULL == rsb__allocate_coo_matrix_t(&coo) )
			goto err;
	}

	IT = coo.IA;
	JT = coo.JA;
	VT = coo.VA;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if(typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))
	{
	mtype * vT = VT;
	mtype * vB = VB;
	mtype * vC = VC;

again`_'RSB_M4_CHOPSPACES(mtype):

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_LT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

	/* FIXME: this works as RSB_FLAG_DUPLICATES_SUM but should support either merge, last, first, ...  */
       	if   ( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi] + vC[ci];
		++bi,++ci,++ti;
		++onz;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_EQ(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi] + vC[ci];
		++bi;
		++ci;
		++ti;
		++onz;
	}

       	if   ( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( bi<nnzB && ci<nnzC && RSB_COO_GT(IB[bi],JB[bi],IC[ci],JC[ci]) && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC && bi < nnzB )
		goto again`_'RSB_M4_CHOPSPACES(mtype);

again`_once_'RSB_M4_CHOPSPACES(mtype):

       	if   ( bi<nnzB && ci==nnzC )
	{
		IT[ti] = IB[bi];
		JT[ti] = JB[bi];
		vT[ti] = vB[bi];
		++bi,++ti;
	}

       	while( bi<nnzB && ci==nnzC && ti > 0 && RSB_COO_EQ(IB[bi],JB[bi],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti] += vB[bi];
		++bi;
		++ti;
		++onz;
	}

       	if   ( ci<nnzC && bi==nnzB )
	{
		IT[ti] = IC[ci];
		JT[ti] = JC[ci];
		vT[ti] = vC[ci];
		++ci,++ti;
	}

       	while( ci<nnzC && bi==nnzB && ti > 0 && RSB_COO_EQ(IC[ci],JC[ci],IT[ti-1],JT[ti-1]) )
	{
		--ti;
		vT[ti]+= vC[ci];
		++ci;
		++ti;
		++onz;
	}

	if( ci < nnzC || bi < nnzB )
		goto again`_once_'RSB_M4_CHOPSPACES(mtype);
	}
	else 
#endif
')dnl
		errval = RSB_ERR_INTERNAL_ERROR;

	coo.nnz -= onz;
	RSB_COA_MEMCPY(IA,IT,offB,0,(coo.nnz));
	RSB_COA_MEMCPY(JA,JT,offB,0,(coo.nnz));
	if(wp)
	{
		RSB_A_MEMCPY_parallel(  VA,VT,offB,0,(coo.nnz),es);
	}
	else
	{
		RSB_A_MEMCPY(  VA,VT,offB,0,(coo.nnz),es);
	}
	RSB_ASSERT(rsb__util_is_coo_array_sorted_up_partial_order(IA,coo.nnz));
	goto done;
err:
	errval = RSB_ERR_ENOMEM;
done:
	if( coop && coop->nnz)
		;
	else
		rsb__destroy_coo_matrix_t(&coo);
	RSB_ASSIGN_IF(onzp,onz);
ret:
	return errval;
}
')dnl
dnl

dnl
rsb_err_t rsb__do_copy_converted_scaled(const void *RSB_RESTRICT  src, void *RSB_RESTRICT dst, const void *RSB_RESTRICT  alphap, rsb_type_t stype,rsb_type_t dtype, size_t nnz, rsb_trans_t transA)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * Copies scaled and conj-transposed.
	 * alpha according to src code type.
	 * \return \rsberrcodemsg
	 * */
	rsb_nnz_idx_t nzi;

	if((!dst) || (!src))
		return RSB_ERR_BADARGS;

foreach(`mtypea',RSB_M4_TYPES,`dnl
foreach(`mtypeb',RSB_M4_TYPES,`dnl
	if( stype == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtypea) && dtype == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtypeb) )
	{
		const mtypea alpha = alphap?*(mtypea*)alphap:RSB_M4_ONE(mtypea);
		const mtypea*tsrc = src;
		mtypeb*tdst = dst;
ifelse(RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(mtypeb)),1,`dnl
		if(RSB_DOES_CONJUGATE(transA))
			for(nzi=0;nzi<nnz;++nzi) RSB_M4_ASSIGN(mtypeb,mtypea,`tdst[nzi]',`(mtypeb)(alpha*RSB_M4_CONJ(`tsrc[nzi]',mtypeb,RSB_M4_TRANS_C,RSB_M4_SYMBOL_UNSYMMETRIC))')
		else
')dnl
			for(nzi=0;nzi<nnz;++nzi) RSB_M4_ASSIGN(mtypeb,mtypea,`tdst[nzi]',`(mtypeb)(alpha*tsrc[nzi])')
	}
	else 
')dnl
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_csc2csr(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t*flagsp)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * */
	rsb_nnz_idx_t nzi = 0, nzo;
	rsb_coo_idx_t nr, nc;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_bool_t islowtri = RSB_BOOL_TRUE, isupptri = RSB_BOOL_TRUE;
	rsb_nnz_idx_t lowtrin = 0, upptrin = 0;

	RSB_BZERO(oIA, sizeof(*oIA)*(m+1));
	oIA[0] = offo;
	for(nzi=0;nzi<nnz;++nzi)
		oIA[IA[nzi]-offi+1]++;
	for(nr=0;nr<m;++nr)
		oIA[nr+1]+=oIA[nr];
foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	for(nc=0;nc<k;++nc)
	for(nzi = JA[nc]-offi;nzi<JA[nc+1]-offi;++nzi)
	{
		nzo = oIA[IA[nzi]-offi]++;
		oJA[nzo] = nc+offo;
		((mtype*)oVA)[nzo] = ((const mtype*)VA)[nzi];
	}
	else 
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE;
	for(nc=0;nc<k;++nc)
	for(nzi=JA[nc]-offi;nzi<JA[nc+1]-offi;++nzi)
	{
		oIA[IA[nzi]-offi]--;
		if(IA[nzi]-offi>nc)
			lowtrin++;
		else
			if(IA[nzi]-offi<nc)
				upptrin++;
	}
	if(upptrin)
		islowtri = RSB_BOOL_FALSE;
	if(lowtrin)
		isupptri = RSB_BOOL_FALSE;
	if(isupptri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
	if(islowtri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
	if( RSB_XOR(upptrin,lowtrin) )
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
	if(*flagsp) RSB_DO_FLAG_ADD(*flagsp,flags);
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_coo_copy_and_stats(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, rsb_coo_idx_t*m, rsb_coo_idx_t*k, const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t iflags, rsb_flags_t*flagsp)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
         * FIXME: unfinished! shall support also typecode-based zeros removal
	 * */
	rsb_nnz_idx_t nzi = 0;
	rsb_coo_idx_t maxi = 0,maxj = 0;
	rsb_bool_t islowtri = RSB_BOOL_TRUE,isupptri = RSB_BOOL_TRUE;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	rsb_nnz_idx_t lowtrin = 0,upptrin = 0;

	if(nnz<1)
		goto done;

foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{
	rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
	maxi = i, maxj = j;
	((mtype*)oVA)[nzi] = ((mtype*)VA)[nzi];
	oIA[nzi] = i-offi+offo;
	oJA[nzi] = j-offi+offo;
	lowtrin |= (i>j), upptrin |= (i<j);
	for(nzi=1;RSB_LIKELY(nzi<nnz);++nzi)
dnl	RSB_M4_SIMPLE_LOOP_UNROLL(`nzi',`LI',`1',`nnz',`
dnl	/* if ( is non zero ... ) */
	{
dnl		if(IA[nzi+LI]>maxi) maxi = IA[nzi+LI];
dnl		if(JA[nzi+LI]>maxj) maxj = JA[nzi+LI];
dnl		const rsb_coo_idx_t i = IA[nzi+LI],j = JA[nzi+LI];
		rsb_coo_idx_t i = IA[nzi],j = JA[nzi];
		maxi = RSB_MAX(maxi, i);
		maxj = RSB_MAX(maxj, j);
		((mtype*)oVA)[nzi] = ((mtype*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
dnl		((mtype*)oVA)[nzi+LI] = ((mtype*)VA)[nzi+LI];
dnl		oIA[nzi+LI] = i-offi+offo;
dnl		oJA[nzi+LI] = j-offi+offo;
dnl		if(IA[nzi+LI]>JA[nzi+LI])isupptri = RSB_BOOL_FALSE;
dnl		else if(IA[nzi+LI]<JA[nzi+LI])islowtri = RSB_BOOL_FALSE;
		lowtrin |= (i>j);
		upptrin |= (i<j);
	}
dnl	')
}
	else
')dnl
		return RSB_ERR_UNSUPPORTED_TYPE;
	if(upptrin)
		islowtri = RSB_BOOL_FALSE;
	if(lowtrin)
		isupptri = RSB_BOOL_FALSE;
	if(isupptri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_UPPER);
	if(islowtri)
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_LOWER);
	if( RSB_XOR(upptrin,lowtrin) )
		RSB_DO_FLAG_ADD(flags,RSB_FLAG_TRIANGULAR);
	if(flagsp) RSB_DO_FLAG_ADD(*flagsp,flags);
	if(flagsp)
		if(RSB_DO_FLAG_HAS_INTERSECTION(*flagsp,RSB_FLAG_ANY_SYMMETRY))
			maxi = maxj = RSB_MAX(maxi,maxj);
	if(m) *m = maxi+1;
	if(k) *k = maxj+1;
done:
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
rsb_err_t rsb__util_coo_copy(const void *RSB_RESTRICT VA, const rsb_coo_idx_t * RSB_RESTRICT IA, const rsb_coo_idx_t * RSB_RESTRICT JA, void *RSB_RESTRICT oVA, rsb_coo_idx_t * RSB_RESTRICT oIA, rsb_coo_idx_t * RSB_RESTRICT oJA, const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
         * FIXME: unfinished! shall support also typecode-based zeros removal
	 * */
	rsb_nnz_idx_t nzi = 0;
dnl	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{
	for(nzi=0;RSB_LIKELY(nzi<nnz);++nzi)
dnl	RSB_M4_SIMPLE_LOOP_UNROLL(`nzi',`LI',`1',`nnz',`
dnl	/* if ( is non zero ... ) */
	{
dnl		const rsb_coo_idx_t i = IA[nzi+LI],j = JA[nzi+LI];
		rsb_coo_idx_t i = IA[nzi], j = JA[nzi];
		((mtype*)oVA)[nzi] = ((mtype*)VA)[nzi];
		oIA[nzi] = i-offi+offo;
		oJA[nzi] = j-offi+offo;
dnl		((mtype*)oVA)[nzi+LI] = ((mtype*)VA)[nzi+LI];
dnl		oIA[nzi+LI] = i-offi+offo;
dnl		oJA[nzi+LI] = j-offi+offo;
	}
dnl	')
}
	else
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
dnl done:
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

/* sparse blas level 1 equivalent functions */

dnl
int rsb__BLAS_Xusdot(const rsb_type_t typecode, const enum blas_conj_type conj_arg, const rsb_blas_int_t nz, const void*x, const rsb_blas_int_t*indx, const void*y, const rsb_blas_int_t incy, void*r, const enum blas_base_type index_base)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
		\rsb_spblasl1_dot_msg
		\rsb_warn_untested_msg
	*/
foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{
	mtype*xa = (mtype*)x;
	mtype*ya = (mtype*)y;
	mtype*rp = (mtype*)r;
	mtype ac = RSB_M4_ZERO(mtype);
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
ifelse(RSB_M4_AND(RSB_M4_IS_COMPLEX_TYPE(mtype)),1,`dnl
	if( conj_arg == blas_conj )
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += RSB_M4_CONJ(xa[nzi],mtype,RSB_M4_TRANS_C,RSB_M4_SYMBOL_UNSYMMETRIC) * ya[xi*incy];
	}
	else
')dnl
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		ac += xa[nzi] * ya[xi*incy];
	}
	RSB_SET_IF_NOT_NULL(rp,ac);
}
	else
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
int rsb__BLAS_Xusaxpy(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*alpha, const void*x, const rsb_blas_int_t*indx, const void*y, const rsb_blas_int_t incy, const enum blas_base_type index_base)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
		\rsb_spblasl1_axpy_msg
		\rsb_warn_untested_msg
	*/
foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{
	const mtype*xa = (const mtype*)x;
	mtype*ya = (mtype*)y;
	const mtype alphav = *(mtype*)alpha;
	rsb_blas_int_t nzi, xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[nzi*incy] += alphav*xa[xi];
	}
}
	else
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
int rsb__BLAS_Xusga(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*y, const rsb_blas_int_t incy, void*x, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
		\rsb_spblasl1_ga_msg
		\rsb_warn_untested_msg
	*/
foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{
	mtype*xa = (mtype*)x;
	const mtype*ya = (const mtype*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
	}
}
	else
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
int rsb__BLAS_Xusgz(const rsb_type_t typecode, const rsb_blas_int_t nz, void*y, const rsb_blas_int_t incy, void*x, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
		\rsb_spblasl1_gz_msg
		\rsb_warn_untested_msg
	*/
foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{
	mtype*xa = (mtype*)x;
	mtype*ya = (mtype*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
    		xa[nzi] = ya[xi*incy];
		ya[xi*incy] = RSB_M4_ZERO(mtype);
	}
}
	else
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
int rsb__BLAS_Xussc(const rsb_type_t typecode, const rsb_blas_int_t nz, const void*x, void*y, const rsb_blas_int_t incy, const rsb_blas_int_t*indx, const enum blas_base_type index_base)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
		\rsb_spblasl1_sc_msg
		\rsb_warn_untested_msg
	*/
foreach(`mtype',RSB_M4_TYPES,`dnl
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
{
	const mtype*xa = (const mtype*)x;
	mtype*ya = (mtype*)y;
	rsb_blas_int_t nzi,xi;
	if( index_base == blas_one_base )
		ya-=incy;
	for(nzi=0;RSB_LIKELY(nzi<nz);++nzi)
	{
		xi = indx[nzi];
		ya[xi*incy] = xa[nzi];
	}
}
	else
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
/* blas level 1 equivalent functions */
dnl

dnl
rsb_err_t rsb__cblas_Xcopy(rsb_type_t typecode, rsb_nnz_idx_t n, const void * x, rsb_nnz_idx_t incx, void * y, rsb_nnz_idx_t incy)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	return rsb__xcopy_strided_typed(y,x,0,0,n,typecode,incy,incx);
}
')dnl
dnl

dnl
rsb_err_t rsb__cblas_Xnrm2(rsb_type_t type, size_t n, const void * a, rsb_nnz_idx_t incA, void * c)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * c <- sqrt(sum(|a_i|^2))
         *
	 * \param a	an array pointer
	 * \param type	a valid type code
	 * \param n	the input array length
	 * \note see dznrm2 in BLAS
	 *
	 * \return \rsberrcodemsg
	 * */
	size_t i;
foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( type == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
	const mtype*ta = a;RSB_M4_REALT(mtype) *tc = c,acc = RSB_M4_ZERO(mtype),tmp = RSB_M4_ZERO(mtype);
	RSB_M4_SIMPLE_LOOP_UNROLL(`i',`LI',`0',`n',`dnl
	acc = RSB_M4_ABS(mtype,ta[(i+LI)*incA]);tmp += acc*acc;
	'); 
	tc[0] = RSB_M4_CREAL(mtype,RSB_M4_SQRT(mtype,tmp));
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
dnl
rsb_err_t rsb__cblas_Xdotu_sub(rsb_type_t type, size_t n, const void * x, rsb_nnz_idx_t incx, const void * y, rsb_nnz_idx_t incy, void *dotu)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * */
	return rsb__vector_mult_sum(x,y,dotu,type,n,incx,incy);
}
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */

dnl
rsb_err_t rsb__cblas_Xscal(rsb_type_t type, size_t n, const void * alphap, void * a, size_t stride)dnl
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*!
	 * a <- a * alpha
	 * */
	return rsb_strided_vector_scale(a,alphap,type,n,stride);
}
')dnl
dnl

dnl
dnl
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__coo_insertion_sort(rsb_type_t typecode, void* VB, rsb_coo_idx_t * IB, rsb_coo_idx_t * JB, rsb_nnz_idx_t offA, rsb_nnz_idx_t nnzA)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/* only for *small* arrays, where allocation of a temporary array is not justified */
	rsb_coo_idx_t * IA = NULL, *JA = NULL;
	rsb_nnz_idx_t i, j;

	IA = IB + offA;
	JA = JB + offA;

foreach(`mtype',RSB_M4_TYPES,`dnl
`#ifdef 'RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype)
	if( typecode == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) )
	{
		mtype * VA = (mtype*) RSB_TYPED_OFF_PTR(typecode,VB,offA);
		for(i=1;i<nnzA;++i)
		for(j=i;j>0 && RSB_COO_LT(IA[j],JA[j],IA[j-1],JA[j-1]);--j)
		{
			RSB_SWAP(rsb_coo_idx_t,IA[j],IA[j-1]);
			RSB_SWAP(rsb_coo_idx_t,JA[j],JA[j-1]);
			RSB_SWAP(mtype        ,VA[j],VA[j-1]);
		}
	}
	else 
#endif /* RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype) */
')dnl
	return RSB_ERR_UNSUPPORTED_TYPE	;
	return RSB_ERR_NO_ERROR;
}
')dnl
dnl

#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
dnl

void rsb__coo_to_lr( void * RSB_RESTRICT VBu, rsb_coo_idx_t*RSB_RESTRICT IB, rsb_coo_idx_t*RSB_RESTRICT JB, void * RSB_RESTRICT VAu, rsb_coo_idx_t*RSB_RESTRICT IA, rsb_coo_idx_t*RSB_RESTRICT JA, rsb_coo_idx_t mj, rsb_nnz_idx_t nnzA, rsb_nnz_idx_t nzoffB, rsb_nnz_idx_t nzoffA, rsb_nnz_idx_t*RSB_RESTRICT nzlp, rsb_nnz_idx_t*RSB_RESTRICT nzrp, rsb_coo_idx_t iadd, rsb_coo_idx_t jadd, rsb_type_t typecode)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/*
	 * Given COO arrays matrices A and (temporary) B, store the coefficients left of the mj-th column before the one coming after it, respecting their original ordering.
	 * A serial function.
	 * Each of arrays VAu, IA, JA fit nzoffA+nnzA elements.
	 * Each of arrays VBu, IB, JB fit nzoffB+nnzA elements.
	 * *nzlp will be set to the count of the elements left  of mj.
	 * *nzrp will be set to the count of the elements right of mj.
	 * iadd will be added to the each of the nnzA IA element.
	 * jadd will be added to the each of the nnzA JA element which is >= mj.
	 * */
	rsb_nnz_idx_t nzl = 0, nzr = 0, nzi = 0;

	RSB_DEBUG_ASSERT(IA!=IB);
	RSB_DEBUG_ASSERT(JA!=JB);

	IA += nzoffA;
	JA += nzoffA;

	IB += nzoffB;
	JB += nzoffB;
	
switch(typecode)
{
			/* supported RSB_M4_MATRIX_TYPES */
foreach(`type',RSB_M4_MATRIX_TYPES,`dnl
case RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(type)	:
{
	type * VA = ((type *)VAu) + nzoffA; 
	type * VB = ((type *)VBu) + nzoffB; 

	RSB_DEBUG_ASSERT(VA!=VB);

	/* visit IA, JA, VA */
	for(nzi=0;nzi<nnzA;++nzi)
	{
		if( JA[nzi] < mj )
		{
			/* elements left of mj go to the beginning of IB,JB,VB */
			IB[nzl] = IA[nzi] + iadd;
			JB[nzl] = JA[nzi];
			VB[nzl] = VA[nzi];
			nzl++;
		}
		else
		{
			/* elements right of mj go to the end of IB,JB,VB, reversed */
			nzr++;
			IB[nnzA-nzr] = IA[nzi] + iadd;
			JB[nnzA-nzr] = JA[nzi] + jadd;
			VB[nnzA-nzr] = VA[nzi];
		}
	}

	/* first nzl elements lay left of mj; the last nzr are right of mj, reversed */
	RSB_DEBUG_ASSERT( nzl+nzr == nnzA );

	/* copy left elements back to A */
	for(nzi=0;nzi<nzl ;++nzi)
	{
		IA[nzi] = IB[nzi];
		JA[nzi] = JB[nzi];
		VA[nzi] = VB[nzi];
	}
	
	/* copy right elements back to A, reversed to original relative order */
	for(     ;nzi<nnzA;++nzi)
	{
		IA[nzi] = IB[nnzA-(1+nzi-nzl)];
		JA[nzi] = JB[nnzA-(1+nzi-nzl)];
		VA[nzi] = VB[nnzA-(1+nzi-nzl)];
	}
}
	break;
')dnl
	/* unsupported type */
	default :
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS 
}

	*nzlp = nzl;
	*nzrp = nzr;
}')dnl
dnl

dnl
rsb_err_t rsb__util_testall(void)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
dnl
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_blk_idx_t inc;
dnl
	RSB_STDOUT("MTX PRINT TEST BEGIN\n");
foreach(`mtype',RSB_M4_TYPES,`dnl
{
	const mtype iVA [] = {1,2};
	const rsb_coo_idx_t iIA[] = {0,1};
	const rsb_coo_idx_t iJA[] = {0,1};
	mtype oVA [] = {-1,-2};
	rsb_coo_idx_t oJA[] = {-1,-2};
	rsb_coo_idx_t oIA[] = {-1,-2};
	rsb_coo_idx_t m=0,k=0;
	const rsb_nnz_idx_t nnz=2;
	const rsb_type_t typecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const rsb_coo_idx_t offi=0, offo=0;
	mtype exp = RSB_M4_ZERO(mtype);
	const mtype alpha = 2;
	mtype res = 1;
 	rsb_flags_t iflags=RSB_FLAG_NOFLAGS, *flagsp=&iflags;

	inc = 1;

	// {-1,-2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, &m, &k, nnz, typecode, offi, offo, iflags, flagsp);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+1,+2}
	errval = rsb__util_coo_copy_and_stats(iVA, iIA, iJA, oVA, oIA, oJA, NULL, NULL, nnz, typecode, offi, offo, iflags, NULL);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+2} -> {+2,+4}
	errval = rsb__cblas_Xscal(typecode, nnz, &alpha, oVA, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=6) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+2,+4} -> {+4,+8}
	errval = rsb_strided_vector_scale(oVA,&alpha,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=12) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+8} -> {+2,+1}
	oVA[0] /= 2; // avoid pow roundoff errors
	oVA[1] /= 8; // avoid pow roundoff errors
	errval = rsb__util_vector_pow(oVA, typecode, &alpha, nnz);
	// {+2,+1} -> {+4,+1}
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=5) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+4,+1} -> {+1,+1}
	errval = rsb__util_vector_pow(oVA, typecode, &exp, nnz);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	// {+1,+1} -> {+2,+1}
	errval = rsb_strided_vector_scale(oVA, &alpha, typecode, nnz*0+1, inc*2); // test for inc==2
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,oVA,typecode,nnz,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__test_print_coo_mm(typecode, iflags, oIA, oJA, oVA, m, k, nnz, RSB_BOOL_TRUE, RSB_DEFAULT_STREAM);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')
dnl
	RSB_STDOUT("MTX PRINT TEST END\n");

foreach(`mtype',RSB_M4_TYPES,`dnl
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const mtype VA [] = {1,0,2,0,3,0};
	const mtype VAL[17] = {1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1};
	mtype nrm2[2] = {0,-1};
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);

	// TODO: maybe get rid of redundance here.
	errval = rsb__util_vector_sum_strided(&nrm2[0], VA, type, n, inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, n, VA, inc, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(nrm2[0]!=nrm2[0]) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__cblas_Xnrm2(type, 17, VAL, 1, &nrm2[1]);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(roundf(nrm2[1]) != 4 ) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}

for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	mtype VA [] = {+1,-1,+3,-3,+4,-4};
	const mtype VAU[] = {+0,+0,+3,-3,+4,-4};
	const mtype VAA[] = {+0,+0,+3,-3,+0,+0};
	mtype soad = RSB_M4_ZERO(mtype);
	const rsb_nnz_idx_t n = sizeof(VA)/(sizeof(VA[0])*inc);
	#if defined(__INTEL_COMPILER) /* meant to circumvent e.g. icc -O2 -fp-model fast=2 (observed on 19.1.2.254 or 2021.3.0) */
	        const mtype lthr = 3 * 0.95, hthr = 3 * 1.05;
	#else
	        const mtype lthr = 3, hthr = 3;
	#endif

	errval = rsb__util_drop_to_zero_if_under_threshold(VA,type,n,&lthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAU,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_drop_to_zero_if_above_threshold(VA,type,n,&hthr);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__vector_sum_of_abs_diffs(&soad,VA,VAA,type,n);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(soad) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')dnl

dnl
RSB_STDOUT("DIFF PRINT TEST BEGIN\n");
foreach(`mtype',RSB_M4_TYPES,`dnl
{
	const rsb_type_t type = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const mtype VAU[] = {+0,+0,+3,-3,+4,-4};
	const mtype VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0]));

	// TODO: beef this test up; add e.g. quiet option to rsb__debug_print_vectors_diff
	errval = rsb__debug_print_vectors_diff(VAU,VAA,n,type,1,1,0);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')dnl
dnl

dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
for(inc=1;inc<=2;++inc)
{
	const rsb_type_t type = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const mtype VAU[] = {+0,+0,+3,-3,+4,-4};
	const mtype VAA[] = {+0,+0,+3,-3,+0,+0};
	const rsb_nnz_idx_t n = sizeof(VAU)/(sizeof(VAU[0])*inc);

	errval = rsb__debug_print_vectors(VAU,VAA,n,1,inc,type);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__debug_print_vector(VAU,n,type,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')dnl
RSB_STDOUT("DIFF PRINT TEST END\n");
dnl

dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
{
	const int inc = 1;
	const rsb_type_t type = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	mtype v1[] = {+1,+0,+5,-2,+4,-6};
	mtype v2[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t off = 1;
	const rsb_nnz_idx_t n = sizeof(v1)/(sizeof(v1[0])*inc)-off;
	const mtype zero = RSB_M4_ZERO(mtype), sum = 2;
	mtype res = zero;
	rsb_nnz_idx_t cnt;

	cnt = rsb__util_count_positive(v1,type,n+off);
	if(cnt!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	cnt = rsb__util_count_negative(v1,type,n+off);
	if(cnt!=2) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	
	errval = rsb__vectors_left_sum_reduce_and_zero(v1,v2,type,n,inc,off);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v1+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=sum) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__util_vector_sum_strided(&res,v2+off,type,n,inc);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=zero) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')dnl
dnl

dnl
#if RSB_WITH_SPARSE_BLAS_INTERFACE
foreach(`mtype',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
{
	/* BLAS types and interfaces only here  */
	const int incy = 1;
	const rsb_type_t type = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const mtype x[] = {+1,+0,+5,-2,+4,-6};
	mtype y[] = {+2,+0,+3,-2,+0,+0};
	mtype w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const mtype zero = RSB_M4_ZERO(mtype);
	mtype alpha = 2;
	mtype dotr = zero;
	mtype res = zero;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	int istat;
dnl	rsb_blas_int_t istat;
#endif
	const enum blas_base_type bzb = blas_zero_base;
	const enum blas_conj_type bc = blas_conj;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`dot',`u',`ID',`0',`f90')(&bc, &nz, x, &(indx[0]), y, &incy, &dotr, &bzb, &istat);
dnl	istat=RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`dot',`u',`ID',`0',`lang_c')(bc, nz, x, indx, y, incy, &dotr, bzb);
#else
	errval = rsb__BLAS_Xusdot(type, bj, nz, x, indx, y, incy, &dotr, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
ifelse(RSB_M4_IS_COMPLEX_TYPE(mtype),`1',`dnl
dnl	istat=RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`axpy',`u',`ID',`0',`lang_c')(nz,&alpha, x, indx, y, incy, bzb);
	RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`axpy',`u',`ID',`0',`f90')(&nz,&alpha, x, &(indx[0]), y, &incy, &bzb, &istat);
',`dnl
dnl	istat=RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`axpy',`u',`ID',`0',`lang_c')(nz, alpha, x, indx, y, incy, bzb);
	RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`axpy',`u',`ID',`0',`f90')(&nz,&alpha, x, &(indx[0]), y, &incy, &bzb, &istat);
')dnl
#else
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, bzb);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
dnl	istat=RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`ga',`u',`ID',`0',`lang_c')(nz, x, incy, y, indx, bzb);
	RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`ga',`u',`ID',`0',`f90')(&nz, x, &incy, y, &(indx[0]), &bzb, &istat);
#else
	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, bzb);
#endif
	errval=RSB_BLAS_ERROR_TO_RSB_ERROR(istat);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
#if RSB_WANT_SPARSE_BLAS_LEVEL_1
	RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`gz',`u',`ID',`0',`f90')(&nz, y, &incy, w, &(indx[0]), &bzb, &istat);
dnl	istat=RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`gz',`u',`ID',`0',`lang_c')(nz, y, incy, w, indx, bzb);
#else
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, bzb);
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

#if RSB_WANT_SPARSE_BLAS_LEVEL_1
dnl	istat=RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`sc',`u',`ID',`0',`lang_c')(nz,x,y,incy,indx,bzb);
	RSB_M4_SPBLAS_MATRIX_ALL_L1_FUNCTION(mtype,`sc',`u',`ID',`0',`f90')(&nz,x,y,&incy,&(indx[0]),&bzb,&istat);
#else
	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,bzb); 
#endif
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')dnl
#endif /* RSB_WITH_SPARSE_BLAS_INTERFACE */
dnl

dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
{
	const int incy = 1;
	const rsb_type_t type = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype);
	const mtype x[] = {+1,+0,+5,-2,+4,-6};
	mtype y[] = {+2,+0,+3,-2,+0,+0};
	mtype w[] = {+0,+0,+0,+0,+0,+0};
	const rsb_blas_int_t indx[] = {0,1,2,3,3,3};
	const rsb_nnz_idx_t nz = sizeof(indx)/(sizeof(indx[0]));
	const rsb_nnz_idx_t ny = sizeof(y)/(sizeof(y[0]));
	const mtype zero = RSB_M4_ZERO(mtype);
	mtype alpha = 2;
	mtype dotr = zero;
	mtype res = zero;

	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;

	errval = rsb__BLAS_Xusdot(type, blas_conj, nz, x, indx, y, incy, &dotr, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(dotr!=25) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	alpha=-alpha;
	errval = rsb__BLAS_Xusaxpy(type, nz, &alpha, x, indx, y, incy, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} <= {+4,+0,+13,-6,-4,-4} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xusga(type, nz, x, incy, y, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// {+2,+0,+3,-2,+0,+0} => {+1,+0,+5,-2,-2,-2} 
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=0) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	y[3]=-3;
	errval = rsb__BLAS_Xusgz(type, nz, y, incy, w, indx, blas_zero_base);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// w: {0,0,0,0,0,0} => {+1,+0,+5,-3,+0,+0} // the two zeroes becaue of usgz zeroing after first copy
	// y: {+1,+0,+5,-2,-2,-2} => {+0,+0,+0,+0,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,w,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=3) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }

	errval = rsb__BLAS_Xussc(type,nz,x,y,incy,indx,blas_zero_base); 
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	// y: {+0,+0,+0,+0,-2,-2} => {+1,+0,+5,-6,-2,-2}
	errval = rsb__util_vector_sum_strided(&res,y,type,ny,incy);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=-4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')dnl
dnl

dnl
foreach(`stype',RSB_M4_TYPES,`dnl
foreach(`dtype',RSB_M4_TYPES,`dnl
{
	const rsb_type_t stypecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(stype);
	const rsb_type_t dtypecode = RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(dtype);
	const stype x[] = {+1,+0,+5,-2,+4,-6};
	dtype y[] = {+2,+0,+3,-2,+0,+0};
	const rsb_nnz_idx_t nnz = sizeof(y)/(sizeof(y[0]));
	stype alpha = 2;
	dtype res = RSB_M4_ZERO(dtype);

	errval = rsb__do_copy_converted_scaled(x, y, &alpha, stypecode, dtypecode, nnz, RSB_TRANSPOSITION_N);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	errval = rsb__util_vector_sum_strided(&res,y,dtypecode,nnz,1);
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
	if(res!=4) errval = RSB_ERR_INTERNAL_ERROR;
	if(RSB_SOME_ERROR(errval)){ RSB_PERR_GOTO(err,RSB_ERRM_ES) }
}
')
')
dnl

err:
	return errval;;
}
')
dnl
ifdef(`ONLY_WANT_HEADERS',`dnl
`#define 'RSB__MAX_SIZEOF`	( \'dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`
	RSB_MAX(sizeof(type),\')
	`0'foreach(`type',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`)'))
',`')dnl
dnl
dnl

dnl
#ifdef __cplusplus
}
#endif  /* __cplusplus */
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_UTIL_H_INCLUDED */
')
dnl
/* @endcond */
dnl
