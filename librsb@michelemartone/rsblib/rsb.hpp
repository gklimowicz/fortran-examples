/*
Copyright (C) 2017-2022 Michele Martone

This file is part of librsb.

librsb is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

librsb is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with librsb; see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

*/
/*!
 *  \file
 *  \brief  Classes \c RsbLib and \c RsbMatrix provide  native C++ access to \librsb.
 *
 *  \details
 *  
 *  
 *  Most of the \librsb functionality is available via C++ classes \c RsbLib and \c RsbMatrix.
 *  \n
 *  These classes are defined in header file \c <rsb.hpp>, which wraps functionality of \librsb's C interface \c <rsb.h>.
 *  \n
 *  The \c RsbMatrix class can manipulate sparse matrices of several numerical types (same ones as \librsb: \ref matrix_supported_numerical_types_section).
 *  \n
 *  Before using \c RsbMatrix, the library must be initialized by having a \c RsbLib object.
 *  \n
 *  To avoid problems when including this header, don't define preprocessor macros prefixed with  \c RSB_ or \c RSBP_.
 *
 *  For a quick start, check out \ref examples/example.cpp or other examples in its directory.
 *  \include examples/example.cpp
 *  
 *  \author Michele Martone
 * */
/*!
 \if STANDALONEDOC
 \mainpage
 \endif

  Header-only (\c rsb.hpp) wrapper classes for librsb (http://librsb.sf.net).
*/

#ifndef RSB_RSB_HPP_INCLUDED
#define RSB_RSB_HPP_INCLUDED

#include <iostream>	// cout
#include <rsb.h>	// rsb_*
#include <cstdlib>	// std::exit
#include <vector>	// std::vector
#include <tuple>	// std::tuple
#include <limits>	// std::limits
#include <exception>	// std::exception
#include <memory>	// std::unique_ptr
#if defined(RSB_LIBRSB_VER) && ( RSB_LIBRSB_VER >= 10300 )
#include <complex>	/* if std::complex supported, <complex> included by rsb.h */
#endif

#if __cplusplus >= 201709L
#define RSBP_WANT_CPP20 1
#endif

#if RSBP_WANT_CPP20
#include <span>	// std::span
#endif /* RSBP_WANT_CPP20 */

/* ****** librsb C++ interface: INTERNALS SECTION BEGIN (SCROLL DOWN FOR PUBLIC SECTION) ****** */
/* @cond INNERDOC  */
/* option and debug flags */

#define RSBP_WANT_DBG_LWCP 0 /* want librsb wrappers call printout */
#ifndef RSBP_WANT_CXX11
# if defined(__cplusplus) && (__cplusplus>=201103L)
#  define RSBP_WANT_CXX11 1
# else
#  define RSBP_WANT_CXX11 0
# endif 
#endif /* RSBP_WANT_CXX11 */

/* macros */

/* Disable exceptions by defining RSBP_NOTHROW */
// #define RSBP_NOTHROW 1

#define RSBP_ALWAYS_THROW(E) throw(std::exception{}) /* e.g. constructors */
#ifndef RSBP_NOTHROW
#define RSBP_THROW(E) throw(std::exception{})
#else
#define RSBP_THROW(E) 
#endif /* RSBP_NOTHROW */
#define RSBP_PERRORA(E) rsb_perror(stdout, E)
#define RSBP_PERROR(ERRVAL) if ( (ERRVAL) != RSB_ERR_NO_ERROR ) { rsb_perror(stdout,(ERRVAL));RSBP_THROW(E);}
#define RSBP_PERROR_ALWAYS_THROW(ERRVAL) if ( (ERRVAL) != RSB_ERR_NO_ERROR ) { rsb_perror(stdout,(ERRVAL));RSBP_ALWAYS_THROW(E);}

#if RSBP_WANT_CXX11
#define RSBP_NULL nullptr
#define RSBP_CONSTEXPR constexpr
#else /* RSBP_WANT_CXX11 */
#define RSBP_NULL NULL
#define RSBP_CONSTEXPR /* nothing */
#endif /* RSBP_WANT_CXX11 */

#define RSBP_STOPP(EC) {std::exit(EC);}	/* stop program */
#define RSBP_CRERR() RSBP_STOPP(-1)	/* critical error */
#define RSBP_MFRV(RETVAL) RSBP_PERROR(errval_); return (RETVAL); /* member function return value */
#define RSBP_PLOC(OS) OS << __func__ << " in " <<__FILE__ << ":" << __LINE__  /* print location */
#define RSBP_PFLOC(OS) OS << "at " <<__FILE__ << "@" << __LINE__<< ":" << __func__  /* print location */
#define RSBP_LWCP() if (RSBP_WANT_DBG_LWCP) RSBP_PLOC(RSBP_EOS) << "\n" /* librsb wrappers call printout */
#define RSBP_EOS std::cout /* error output stream */
#define RSBP_WANT_VFMC() getenv("RSBP_VERBOSE_CALLS")  /* want verbose function member call */
#define RSBP_VWARNM(MSG) if ( RSBP_WANT_VFMC() ) RSBP_PLOC(RSBP_EOS) << "\n" /* verbose warning message */
#define RSBP_ERRORM(MSG) RSBP_PLOC(RSBP_EOS) << ": error: " << MSG << "\n" /* print error message */

#define RSBP_EFE(ERRVAL,EXP) ERRVAL = EXP; RSBP_PERROR(ERRVAL); /* error value of function evaluate, followed by error handling */
#define RSBP_NFE(EXP) RSBP_EFE(errval_,EXP) /* non const function evaluate */

#define RSBP_MFR(ERRVAL) return return_output<Err_t>(ERRVAL); /* non-void member function return */
#define RSBP_DFE(EXP) RSBP_EFE(errval,EXP) /*     const derived function evaluate*/
#define RSBP_CFE(EXP) RSBP_EFE(errval,EXP) /*     const function evaluate*/

#define RSBP_MDR() RSBP_MFR(errval_); /* default member function return */
#define RSBP_STCHK() RSBP_VWARNM(""); if( this->state_ != Valid ) { RSBP_ERRORM("Invalid matrix state (did you initialize it in the proper sequence ?) !\n"); RSBP_ALWAYS_THROW(); }  /* status check */

// #define RSBP_MSC() RSBP_VWARNM(""); if( this->errval_ != RSB_ERR_NO_ERROR ) { RSBP_ERRORM("Invalid library state !\n"); }   /* matrix state check */ // NOTE: may use this in RSBP_CWMFHO, RSBP_NWMFHO, RSBP_CDMFHO ...
#define RSBP_LFMH() RSBP_STCHK() /* librsb function member header */
#define RSBP_CFMH(DC) rsb_err_t errval = RSB_ERR_NO_ERROR; if (DC) { RSBP_STCHK() ;}   /* librsb const function member header */
#define RSBP_CNMH(DC)                                      if (DC) { RSBP_STCHK() ;}   /* librsb const no-error function member header */

#define RSBP_NWMFHO() RSBP_LWCP(); RSBP_LFMH( ) /* non constant wrapper member function heading opening */
#define RSBP_CWMFHO() RSBP_LWCP(); RSBP_CFMH(1) /*     constant wrapper member function heading opening */
#define RSBP_CDMFHO() RSBP_LWCP(); RSBP_CFMH(1) /*     constant derived member function heading opening */
#define RSBP_CNMFHO()              RSBP_CNMH(1) /*     constant no-err. member function heading opening */
#define RSBP_KWMFHO() RSBP_LWCP(); RSBP_CNMH(0) /* construction wrapper member function heading opening */
#define RSBP_LWMFHO() RSBP_LWCP();  		/* construction wrapper member function heading opening */ \
	rsb_err_t errval = RSB_ERR_NO_ERROR; \
	RSBP_VWARNM(""); if( this->state_ != Begin ) { errval = RSB_ERR_BADARGS; RSBP_ERRORM("Invalid matrix state (did you initialize it in the proper sequence ?) !\n"); RSBP_THROW(errval); RSBP_MFR(errval); }  /* status check */

#define RSBP_CWMFRE() RSBP_MFR(errval) /*     constant wrapper member function return error value */

#define RSBP_CFMO() RSBP_RMCI(); RSBP_VWARNM(""); {  /* constructor function member heading opening */
#define RSBP_CFMC() } /* constructor function member closing */

#define RSBP_WANT_IMPLICIT_ORDER 0 // pertains spmm overloads with implicit order
#define RSBP_WANT_INLINE_INIT 0 // Putting this to 1 is for debug mode only.
#define RSBP_ERRGOTO(ERRVAL,ERRLABEL) { if ( (ERRVAL) != RSB_ERR_NO_ERROR ) { goto ERRLABEL; } };
#define RSBP_CANNOT_HONOUR_CSR(NNZ,ROWS,FLAGSA) ( NNZ < 2*((ROWS)+1) && ! RSB_DO_FLAG_HAS(FLAGSA,RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS) )
#ifdef RSB_TRANSPOSITION_INVALID
#define RSBP_INVALID_TRANS_CHAR RSB_TRANSPOSITION_INVALID /* since librsb-1.2.0.10 */
#else /* RSB_TRANSPOSITION_INVALID */
#define RSBP_INVALID_TRANS_CHAR '?'
#endif /* RSB_TRANSPOSITION_INVALID */
/* ****** librsb C++ interface: INTERNALS SECTION END (SCROLL DOWN FOR PUBLIC SECTION) ****** */

#ifndef RSBP_NO_NAMESPACE
namespace rsb {
#endif /* RSBP_NO_NAMESPACE */

using rsb_string_t = std::string; // string type

#ifdef RSBP_WANT_CPP20
template<typename S>
concept RSBP_Scalar_t =
#ifdef RSB_NUMERICAL_TYPE_FLOAT
     std::is_same_v<S, float>
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
  || std::is_same_v<S, long double>
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
  || std::is_same_v<S, double>
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
  || std::is_same_v<S, std::complex<float>>
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
  || std::is_same_v<S, std::complex<double>>
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
  || std::is_same_v<S, std::complex<long double>>
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
  || std::is_same_v<S, int>
#endif
   ;
#else
#define RSBP_Scalar_t typename
#endif

#ifdef RSBP_WANT_CPP20
template<typename out>
concept RSBP_Err_t =
  std::is_same_v<out, void> ||
  std::is_same_v<out, rsb_err_t>;
#else
#define RSBP_Err_t typename
#endif

#ifdef RSBP_WANT_CPP20
template<RSBP_Err_t Err_t>
Err_t return_output(rsb_err_t errval) {
  if constexpr (std::is_same_v<Err_t, void>){
  } else {
    return errval;
  }
}
#else
template<typename Err_t>
Err_t return_output(rsb_err_t ) {
    RSBP_PERROR(RSB_ERR_BADARGS);
}
template<>
void return_output<void>(rsb_err_t errval) {
    RSBP_PERROR(errval);
}
template<>
rsb_err_t return_output<rsb_err_t>(rsb_err_t errval) {
        return errval;
}
#endif
/* @endcond */

/*! \brief Internal attribute specifier for deprecated member functions. */
#if __cplusplus >= 201709
#define RSBP_DEPRECATED [[deprecated]]
#else
#define RSBP_DEPRECATED
#endif

/*! \brief Internal attribute. */
#if __cplusplus >= 201709
#define RSBP_NODISCARD  [[nodiscard]]
#else
#define RSBP_NODISCARD
#endif

#ifndef RSBP_WANT_REV
/*! \brief If this is defined to 1 before including <rsb.hpp>, \ref rsb_err_t is the default return type. Otherwise the default is void. */
#define RSBP_WANT_REV 0
#endif

#if RSBP_WANT_REV
/*! \brief The librsb error flag type. */
#define RSBP_RVT template <typename Err_t=rsb_err_t> RSBP_NODISCARD
#else /* RSBP_WANT_REV */
/*! \brief No return type. */
#define RSBP_RVT template <typename Err_t=void>
#endif /* RSBP_WANT_REV */
/*! \brief Minimal supported librsb version (value of RSB_LIBRSB_VER, defined via rsb.h) */
#define RSBP_MSLVRV  10201

/*!
\brief	Class initializing/finalizing \librsb state.

\details
Before creating any \c RsbMatrix objects, make sure you have initialized \librsb  by creating one single \c RsbLib object.
Similarly, \c RsbMatrix object shall leave scope after the \c RsbLib object has been deallocated.

Several \librsb options can be queried or changed via e.g.
 \ref get_opt(),
 \ref set_opt(),
 \ref get_num_threads(),
 \ref set_num_threads().
*/
class RsbLib
{
private:
	bool rsb_cpp_initialized_ {false};

	void initialize(struct rsb_initopts * iop = RSBP_NULL)
	{
		if( ! rsb_cpp_initialized_ )
		{
			rsb_err_t errval = RSB_ERR_NO_ERROR;
#ifdef RSB_LIBRSB_VER
			if( RSB_LIBRSB_VER < RSBP_MSLVRV )
			{
				errval = RSB_ERR_UNSUPPORTED_FEATURE;
				RSBP_ERRORM("This librsb version (see RSB_LIBRSB_VER) is " << RSB_LIBRSB_VER << " -- and is too old: need at least " << RSBP_MSLVRV );
				RSBP_ALWAYS_THROW();
			}
#endif /* RSB_LIBRSB_VER */
			errval = rsb_lib_init( iop ? iop : RSB_NULL_INIT_OPTIONS );
			if( errval != RSB_ERR_NO_ERROR )
			{
				RSBP_ALWAYS_THROW();
			}
			rsb_cpp_initialized_ = true;
		}
	}
public:

	RSBP_RVT Err_t set_opt_str(const rsb_char_t* opnp, const rsb_char_t* opvp)
	{
		/*!
			Interface to rsb_lib_set_opt_str.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = rsb_lib_set_opt_str(opnp,opvp);
		RSBP_MFR(errval);
	}

	RSBP_RVT Err_t set_opt(enum rsb_opt_t iof, const void*iop)
	{
		/*!
			Interface to rsb_lib_set_opt.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = rsb_lib_set_opt(iof,iop);
		RSBP_MFR(errval);
	}

	RSBP_RVT Err_t get_opt(enum rsb_opt_t iof, void*iop) const
	{
		/*!
			Interface to rsb_lib_get_opt.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = rsb_lib_get_opt(iof, iop);
		RSBP_MFR(errval);
	}

	rsb_string_t get_opt(enum rsb_opt_t iof) const
	{
		/*!
			Interface to rsb_lib_get_opt.
			\warning	Only RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING is supported for now.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		const char * sp {nullptr};
		if(iof != RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING)
		{
			errval = RSB_ERR_UNIMPLEMENTED_YET;
			RSBP_PERROR(errval);
		}
		else
		{
			RSB_REINIT_SINGLE_VALUE_GET(RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING,&sp,errval);
		}
		RSBP_PERROR(errval);
		return rsb_string_t(sp?sp:"");
	}

	RSBP_RVT Err_t set_num_threads(rsb_int_t nt)
	{
		/*!
			Indirect interface to rsb_lib_set_opt.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = set_opt<rsb_err_t>(RSB_IO_WANT_EXECUTING_THREADS, &nt);
		RSBP_MFR(errval);
	}

	rsb_int_t get_num_threads(void) const
	{
		/*!
			Indirect interface to rsb_lib_get_opt.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		rsb_int_t nt{};
		errval = get_opt<rsb_err_t>(RSB_IO_WANT_EXECUTING_THREADS, &nt);
		RSBP_PERROR(errval);
		return nt;
	}

	RSBP_RVT Err_t reinit(struct rsb_initopts * iop)
	{
		/*!
			Interface to \ref rsb_lib_reinit().
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = rsb_lib_reinit(iop);
		RSBP_MFR(errval);
	}

	RSBP_DEPRECATED RsbLib ( bool verbose ):
		rsb_cpp_initialized_(false)
	{
		/*!
			Initialize librsb.  \sa \ref rsb_lib_reinit().
		*/
		struct rsb_initopts tio;
		auto of { stdout };
		std::pair<enum rsb_opt_t,void*> vio { RSB_IO_WANT_VERBOSE_INIT, &of };
		tio.keys = & vio.first;
		tio.values = & vio.second;
		tio.n_pairs = 1;
		tio.action = RSB_IO_SPECIFIER_SET;
		initialize(verbose ? &tio : RSBP_NULL);
	}

	RsbLib (void):
		rsb_cpp_initialized_(false)
	{
		/*!
			Initialize librsb.  \sa \ref rsb_lib_reinit().
		*/
		initialize();
	}

	size_t meminfo(void)
	{
		/*!
			Provide memory debug info from librsb and return usage amount.
			Only effective if librsb configured accordingly.

			On error, throw an exception.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		size_t totmem = 0, totall = 0;
		RSB_REINIT_SINGLE_VALUE_GET(RSB_IO_WANT_MEM_ALLOC_CNT,&totall,errval);
		if(errval==RSB_ERR_UNSUPPORTED_FEATURE)
			return 0;
		RSBP_PERROR(errval);
		RSB_REINIT_SINGLE_VALUE_GET(RSB_IO_WANT_MEM_ALLOC_TOT,&totmem,errval);
		RSBP_PERROR(errval);
		if(totmem || totall)
			RSBP_EOS << "librsb has allocated " << totmem << " bytes in " << totall << " chunks.\n";
		return totmem;
	}

	~RsbLib (void)
	{
		/*!
			Destructor: finalize \librsb.

			This is being invoked at the end of the scope of a \c RsbLib object:
			typically at application's end.
		*/
		rsb_err_t errval = RSB_ERR_NO_ERROR;

		auto of { stdout }; /* can't assume is an lvalue and therefore usable in '&stdout' */
		RSB_REINIT_SINGLE_VALUE_SET(RSB_IO_WANT_VERBOSE_EXIT,&of,errval);
		errval |= rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
		if( errval != RSB_ERR_NO_ERROR )
			RSBP_CRERR();
	}
}; /* RsbLib */

/* @cond INNERDOC  */
#if RSBP_WANT_INLINE_INIT
static RsbLib rsb_lib;
#define RSBP_RMCI() rsb_lib.initialize();
#else /* RSBP_WANT_INLINE_INIT */
#define RSBP_RMCI() 
#endif /* RSBP_WANT_INLINE_INIT */

namespace rsb_internal
{
template<typename NT> constexpr rsb_type_t rsb_type_t_for(void) { return RSB_NUMERICAL_TYPE_INVALID_TYPE; }
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
template<> constexpr rsb_type_t rsb_type_t_for<long double> (void) { return RSB_NUMERICAL_TYPE_LONG_DOUBLE; }
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
template<> constexpr rsb_type_t rsb_type_t_for<double> (void) { return RSB_NUMERICAL_TYPE_DOUBLE; }
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
template<> constexpr rsb_type_t rsb_type_t_for<std::complex<long double>> (void) { return RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX; }
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
template<> constexpr rsb_type_t rsb_type_t_for<std::complex<double>> (void) { return RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX; }
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
template<> constexpr rsb_type_t rsb_type_t_for<std::complex<float>> (void) { return RSB_NUMERICAL_TYPE_FLOAT_COMPLEX; }
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
template<> constexpr rsb_type_t rsb_type_t_for<float> (void) { return RSB_NUMERICAL_TYPE_FLOAT; }
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
template<> constexpr rsb_type_t rsb_type_t_for<int> (void) { return RSB_NUMERICAL_TYPE_INT; }
#endif
}
/* @endcond */

/*!
\brief
Represent a sparse matrix in RSB format by means of \librsb.

\details
Manage construction, destruction, and std::move of numerical matrices.
\n
Most of the member functions here translate directly to a single function call to \librsb (\c rsb.h),
and pass the parameters as they are, so the error checking is done by \librsb.
\n
While most of \librsb \c C  functions use \c void* pointers instead of numerical data, \c RsbMatrix is templated by a type parameter.
This introduces type safety at compile time.

Users of member functions can choose among several overloads.
So in additional to the more direct overloads passing e.g. \f$ \alpha \f$ and \f$ \beta \f$ by reference, here a user can pass them by value.

\param	NT	the numerical type, at least for the four canonical ones (\c float, \c double, \c std::complex<float>, \c std::complex<double>); see \ref matrix_supported_numerical_types_section and #rsb_type_t for more.

\note
Default error propagation is by exception throw for all constructors and most member functions.
\n
Functions declared to return Err_t can be specialized in rsb_err_t so not to throw exceptions, but to return an error code instead.
\n
Exceptions thrown by member functions (not constructors) can be deactivated at build time by defining \ref RSBP_NOTHROW before including \c <rsb.hpp>.

\note
One may turn on return error value as default at build time by defining \ref RSBP_WANT_REV=1 .

\warning
The error model is work in progress and subject to change.

\todo
While the \c rsb.h  interface is stable, 
the \c rsb.hpp  interface is neither stable, nor complete: early users' feedback is very welcome.
\n
Shall all \f$ \alpha \f$ and \f$ \beta \f$ be passed by values only?
This is natural and fits C++.
But rsb_tune_spmm() / rsb_tune_spsm() have more parameters with a nullptr-like `default' value: how to deal with them in RsbMatrix::tune_spmm() / RsbMatrix::tune_spsm() without introducing too many special cases ?
\n
Similarly for the consistency of
 RsbMatrix::get_flags_t(),
 RsbMatrix::get_type_t(),
 RsbMatrix::get_info_coo_t(),
and similar: shall they throw an exception if given a flag not matching the return type ?
\n
Or maybe keeping only RsbMatrix::get_info(), with per-type reference overloads ?
\n
Furthermore: if working under C++20, shall \rsblib avoid pointer-based interfaces completely (using \c std::span only)?
*/
template<RSBP_Scalar_t NT>
class RsbMatrix
{
private:
#ifdef RSBP_TESTING_ONLY
	friend class MatrixConstructors_Test;
	FRIEND_TEST(MatrixConstructors_Test,MoveAssignment);
	FRIEND_TEST(MatrixConstructors_Test,Move);
	FRIEND_TEST(MatrixConstructors_Test,InternalsMoveBadType);
	FRIEND_TEST(MatrixConstructors_Test,InternalsCtorBadType);
	FRIEND_TEST(MatrixConstructors_Test,InternalsCtorTypeOk);
	FRIEND_TEST(MatrixConstructors_Test,InternalsCtorTypeBad);
	FRIEND_TEST(MatrixConstructors_Test,clone);
	FRIEND_TEST(LowerTest,compare_same);
	FRIEND_TEST(LowerTest,compare_different);
#endif /* RSBP_TESTING_ONLY */
	rsb_mtx_t * mtxAp_ {RSBP_NULL};
private:
	rsb_err_t errval_ { RSB_ERR_NO_ERROR };
	enum State { Invalid, Begin, Valid };
	enum State state_ { Invalid };
	const static rsb_type_t typecode_ = rsb_internal::rsb_type_t_for<NT>();
	RSBP_CONSTEXPR static const NT one = 1.0;
	RSBP_CONSTEXPR static const NT zero = 0.0;
	RSBP_CONSTEXPR static const NT defmultbeta = 0.0;
	RSBP_CONSTEXPR static const NT defmultalpha = 1.0;
	RSBP_CONSTEXPR static const NT RSBP_NUMT_INVALID = std::numeric_limits<const NT>::max(); // aka DBL_MAX

public:

	/*! Matrix structure: either general, symmetric, hermitian, or triangular (also see #rsb_flags_t). */
	enum RsbSym { IsGen = RSB_FLAG_NOFLAGS /*!< General matrix, no triangle structure or symmetry assumed. */, IsHer = RSB_FLAG_HERMITIAN  /*!< Hermitian (\f$ A == A^H \f$). Please pass only lower/upper triangle. */, IsSym =  RSB_FLAG_SYMMETRIC /*!< Symmetric (\f$ A == A^T \f$). Please pass only lower/upper triangle. */, IsTri = RSB_FLAG_TRIANGULAR /*!< Triangular (required for \ref spsv/\ref spsm). */ };

private:
	rsb_flags_t _adjustSym( const RsbSym sym = IsGen ) const
	{
		rsb_flags_t flags = RSB_FLAG_NOFLAGS;

		if ( sym == IsGen )
			RSB_DO_FLAG_ADD( flags, RSB_FLAG_NOFLAGS);
		if ( sym == IsHer )
			RSB_DO_FLAG_ADD( flags, RSB_FLAG_HERMITIAN);
		if ( sym == IsSym )
			RSB_DO_FLAG_ADD( flags, RSB_FLAG_SYMMETRIC);
		if ( sym == IsTri )
			RSB_DO_FLAG_ADD( flags, RSB_FLAG_TRIANGULAR);
		return flags;
	}

	bool _is_integer(void) const
	{
		return std::numeric_limits<const NT>::is_integer;
	}

	RSBP_RVT Err_t perror (const rsb_err_t errval) const
	{
		RSBP_MFR(errval);
	}

	RSBP_RVT Err_t perror(void) const
	{
		this->perror(errval_);
		return errval_;
	}

public:

	RsbMatrix(rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, const RsbSym sym = IsGen ):
		mtxAp_ (RSBP_NULL), 
		errval_ (RSB_ERR_NO_ERROR),
		state_ ( Invalid )
       	{
		/*!
			Begin assembling a sparse matrix of given dimensions and type.
			\n
			Then you can use \ref set_val() or \ref set_vals() repeatedly to populate the matrix.
			\n
			After populating the matrix, use \ref close() to terminate its assembly.
			\n

			\rsblib_asm_snip

			\rsblib_based_on_pred rsb_mtx_alloc_from_coo_begin() \rsblib_based_on_post
		*/
		RSBP_CFMO();
		const rsb_nnz_idx_t nnzA = 0;

		this->mtxAp_ = rsb_mtx_alloc_from_coo_begin(nnzA, typecode_, nrA, ncA, _adjustSym( sym ), &errval_);
		if( errval_ == RSB_ERR_NO_ERROR )
			state_ = Begin;
		RSBP_PERROR_ALWAYS_THROW(errval_);
		RSBP_CFMC();
	}

	RsbMatrix( rsb_coo_idx_t nrA, const rsb_coo_idx_t * RP, const rsb_coo_idx_t * JA, const NT * VA, const RsbSym sym = IsGen ):
		mtxAp_ (RSBP_NULL), 
		errval_ (RSB_ERR_NO_ERROR),
		state_ ( Invalid )
       	{
		/*!
			Assemble a sparse matrix given CSR input.
			\n
			\rsblib_based_on_pred rsb_mtx_alloc_from_csr_const() \rsblib_based_on_post
		*/
		RSBP_CFMO();
		const rsb_coo_idx_t ncA = 0;
		const rsb_blk_idx_t brA = 0, bcA = 0;

		if ( RP )
			this->mtxAp_ = rsb_mtx_alloc_from_csr_const(VA, RP, JA, RP[nrA], typecode_, nrA, ncA, brA, bcA, _adjustSym( sym ), &errval_);
		else
			this->mtxAp_ = RSBP_NULL,
			errval_ = RSB_ERR_GENERIC_ERROR;

		if( errval_ == RSB_ERR_NO_ERROR )
			state_ = Valid;

		RSBP_PERROR_ALWAYS_THROW(errval_);
		RSBP_CFMC();
	}

#if RSBP_WANT_CPP20
	RsbMatrix( const std::span<const rsb_coo_idx_t> IA, const std::span<const rsb_coo_idx_t> JA, const std::span<const NT> VA, rsb_nnz_idx_t nnzA = RSB_INVALID_NNZ_IDX_VAL, const rsb_flags_t flagsA = RSB_FLAG_NOFLAGS ):
		mtxAp_ (RSBP_NULL), 
		errval_ (RSB_ERR_NO_ERROR),
		state_ ( Invalid )
	{
		/*!
			Assemble a sparse matrix given COO input.
			\n
			\rsblib_based_on_pred rsb_mtx_alloc_from_coo_const() \rsblib_based_on_post

			Example snip from examples/span.cpp: \snippet examples/span.cpp snip__span_RsbMatrix
		*/
		RSBP_CFMO();
		const rsb_coo_idx_t nrA = 0, ncA = 0;
		const rsb_blk_idx_t brA = 0, bcA = 0;
		if ( nnzA == RSB_INVALID_NNZ_IDX_VAL )
			nnzA = VA.size();
		this->mtxAp_ = rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode_, nrA, ncA, brA, bcA, flagsA, &errval_);

		if( errval_ == RSB_ERR_NO_ERROR )
			state_ = Valid;

		RSBP_PERROR_ALWAYS_THROW(errval_);
		RSBP_CFMC();
	}
#endif /* RSBP_WANT_CPP20 */

	RsbMatrix( const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const NT * VA, rsb_nnz_idx_t nnzA, const rsb_flags_t flagsA = RSB_FLAG_NOFLAGS ):
		mtxAp_ (RSBP_NULL), 
		errval_ (RSB_ERR_NO_ERROR),
		state_ ( Invalid )
       	{
		/*!
			Assemble a sparse matrix given COO input.
			\n
			\rsblib__asm_coo_snip
			\n
			\rsblib_based_on_pred rsb_mtx_alloc_from_coo_const() \rsblib_based_on_post
		*/
		RSBP_CFMO();
		const rsb_coo_idx_t nrA = 0, ncA = 0;
		const rsb_blk_idx_t brA = 0, bcA = 0;

		if ( ( IA && JA && VA ) || nnzA == 0 )
			this->mtxAp_ = rsb_mtx_alloc_from_coo_const(VA, IA, JA, nnzA, typecode_, nrA, ncA, brA, bcA, flagsA, &errval_);
		else
			this->mtxAp_ = RSBP_NULL,
			errval_ = RSB_ERR_GENERIC_ERROR;

		if( errval_ == RSB_ERR_NO_ERROR )
			state_ = Valid;

		RSBP_PERROR_ALWAYS_THROW(errval_);
		RSBP_CFMC();
	}

	RsbMatrix( const rsb_char_t * filename, const RsbSym sym = IsGen ):
		mtxAp_ (RSBP_NULL), 
		errval_ (RSB_ERR_NO_ERROR),
		state_ ( Invalid )
       	{
		/*!
			Assemble a sparse matrix given filename input.
			\n
			\rsblib_based_on_pred rsb_file_mtx_load() \rsblib_based_on_post
		*/
		RSBP_CFMO();

		this->mtxAp_ = rsb_file_mtx_load(filename, _adjustSym( sym ), typecode_, &errval_);
		if( errval_ == RSB_ERR_NO_ERROR )
			state_ = Valid;
		RSBP_PERROR_ALWAYS_THROW(errval_);
		RSBP_CFMC();
	}

	RsbMatrix( const RsbMatrix & A_Rsb, bool do_trans = false, rsb_flags_t flagsA = RSB_FLAG_NOFLAGS ):
		mtxAp_ (RSBP_NULL), 
		errval_ (RSB_ERR_NO_ERROR),
		state_ ( Invalid )
       	{
		/*!
			Copy a sparse matrix given example input.
			\n
			Can either clone it, or transpose it or change flags (structure) in the process.
		*/
		RSBP_CFMO();
		const rsb_trans_t transA = do_trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;

		if ( RSBP_CANNOT_HONOUR_CSR(A_Rsb.nnz(), A_Rsb.rows(), flagsA) )
			RSB_DO_FLAG_ADD( flagsA, RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS ); // Note: this is parachute code.
		errval_ = A_Rsb.clone<rsb_err_t>(&this->mtxAp_/*, rsb_type_t typecode*/, transA, RSBP_NULL, flagsA);

		if( errval_ == RSB_ERR_NO_ERROR )
			state_ = Valid;
		RSBP_PERROR_ALWAYS_THROW(errval_);
		RSBP_CFMC();
	}

#if defined(RSBP_TESTING_ONLY) || defined(RSBT_TESTING_ONLY) /* it's an users's option to define this */
	RsbMatrix( rsb_mtx_t*mtxAp ):
		mtxAp_ (mtxAp), 
		errval_ (RSB_ERR_NO_ERROR),
		state_ ( mtxAp ? Valid : Invalid )
       	{
		/*!
			Create a sparse matrix given internal pointer input (owning).
			This is an internal function: don't use it.
		*/
		RSBP_CFMO();
		errval_ = type_check<rsb_err_t>();
		RSBP_PERROR_ALWAYS_THROW(errval_);
		RSBP_CFMC();
	}
#endif /* */

private:
	RsbMatrix(void) = default; // ugly, isn't it ? only here for the move ctor
public:
	RsbMatrix(RsbMatrix && other)
	{
		/*!
			Move constructor.
			\n
			The moved matrix object will be invalid afterwards.
			\n
			Example snip from examples/misc.cpp: \snippet examples/misc.cpp snip__move_RsbMatrix
		*/
		{
			RsbMatrix swappable;
			this->swap(swappable);
		}
		this->swap(other);
	}

	~RsbMatrix(void)
	{
		/*!
			Destructor.
			\n
			Frees matrix object memory.
			\n
			\rsblib_based_on_pred rsb_mtx_free() \rsblib_based_on_post
		*/
		this->mtxAp_ = rsb_mtx_free (this->mtxAp_);
		state_ = Invalid;
	}

	RSBP_RVT RSBP_DEPRECATED Err_t _add(rsb_coo_idx_t i, rsb_coo_idx_t j, NT val)
	{
		/*!
			\deprecated Use \ref set_val() and \ref set_vals() instead.
		*/
		RSBP_KWMFHO();
		RSBP_NFE(set_val<rsb_err_t> (val, i, j, RSB_FLAG_NOFLAGS));
		RSBP_MDR();
	}
	      
	RSBP_RVT Err_t close(void)
	{
		/*!
			Terminate assembly of a previously started and populated matrix.
			\n
			Shall be called once.

			\rsblib_asm_snip

			\sa RsbMatrix::RsbMatrix(rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, const RsbSym sym = IsGen );
			\rsblib_based_on_pred rsb_mtx_alloc_from_coo_end() \rsblib_based_on_post
		*/
		RSBP_LWMFHO();
		RSBP_NFE(rsb_mtx_alloc_from_coo_end (&this->mtxAp_));
		if( errval_ == RSB_ERR_NO_ERROR )
			state_ = Valid;
		RSBP_MDR();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t _close(void)
	{
		/*!
			\deprecated Use \ref close() instead.
		*/
		RSBP_LWMFHO();
		RSBP_NFE(close<rsb_err_t>());
		RSBP_MDR();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t spmv(rsb_trans_t transA, const NT *alphap, const NT * Xp, rsb_coo_idx_t incX, const NT * betap, NT * Yp, rsb_coo_idx_t incY) const
	{
		/*!
		\rsblib_based_on_pred rsb_spmv() \rsblib_based_on_post
		\rl_deprecate
		*/
		RSBP_CWMFHO();
		RSBP_CFE( rsb_spmv(transA, alphap, this->mtxAp_, Xp, incX, betap, Yp, incY) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spmv(rsb_trans_t transA, const NT alpha, const NT * Xp, rsb_coo_idx_t incX, const NT beta, NT * Yp, rsb_coo_idx_t incY) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmv() \rsblib_based_on_post

		\rsblib_spmv_snip
		*/
		RSBP_CDMFHO();
		RSBP_DFE( spmv<rsb_err_t>(transA, &alpha, Xp, incX, &beta, Yp, incY) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spmv(rsb_trans_t transA, const NT alpha, const NT * Xp, const NT beta, NT * Yp) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmv() \rsblib_based_on_post

		\rsblib_spmv_snip
		*/
		RSBP_CDMFHO();
		const rsb_coo_idx_t incX = 1, incY = 1;
		RSBP_DFE( spmv<rsb_err_t>(transA, &alpha, Xp, incX, &beta, Yp, incY) );
		RSBP_CWMFRE();
	}

private:
	std::pair<rsb_nnz_idx_t,rsb_nnz_idx_t> _get_ldX_ldY(rsb_trans_t transA, rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, rsb_coo_idx_t nrhs = 1) const
	{
		if (order != RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
			return {nrhs,nrhs};
		const auto ldX = ( transA == RSB_TRANSPOSITION_N ) ? this->cols() : this->rows();
		const auto ldY = ( transA == RSB_TRANSPOSITION_N ) ? this->rows() : this->cols();
		return {ldX,ldY};
	}
public:

#if RSBP_WANT_CPP20
	RSBP_RVT Err_t spmv(rsb_trans_t transA, const NT alpha, const std::span<const NT> X, const NT beta, const std::span<NT> Y) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmv() \rsblib_based_on_post

		\rsblib_spmv_snip
		*/
		RSBP_CWMFHO();
		const rsb_coo_idx_t incX = 1, incY = 1;
		const auto [ldX, ldY] = _get_ldX_ldY(transA);

		if ( static_cast<rsb_nnz_idx_t>(X.size()) < ldX || static_cast<rsb_nnz_idx_t>(Y.size()) < ldY )
		{
			errval = RSB_ERR_BADARGS;
			RSBP_PERROR(errval);
			RSBP_MFR(errval);
		}
		RSBP_CFE( rsb_spmv(transA, &alpha, this->mtxAp_, X.data(), incX, &beta, Y.data(), incY) );
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_CPP20 */

	RSBP_RVT Err_t spmv(NT * y, const NT * x, bool do_trans = false) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmv() \rsblib_based_on_post

		\rsblib_spmv_snip
		*/
		RSBP_CDMFHO();
		const rsb_trans_t transA = do_trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;
		const NT alpha = defmultalpha;
		const NT beta = defmultbeta;
		const rsb_coo_idx_t incX = 1, incY = 1;
		RSBP_DFE( spmv<rsb_err_t>(transA, &alpha, x, incX, &beta, y, incY) );
		RSBP_CWMFRE();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t spmm(rsb_trans_t transA, const NT * alphap, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * Bp, rsb_nnz_idx_t ldB, const NT * betap, NT * Cp, rsb_nnz_idx_t ldC) const
	{
		/*!
		\rsblib_based_on_pred rsb_spmm() \rsblib_based_on_post
		\rl_deprecate
		*/
		RSBP_CWMFHO();
		RSBP_CFE ( rsb_spmm (transA,alphap,this->mtxAp_,nrhs,order,Bp,ldB,betap,Cp,ldC) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spmm(rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * Bp, rsb_nnz_idx_t ldB, const NT beta, NT * Cp, rsb_nnz_idx_t ldC) const
	{
		/*!
		\rsblib_spmm_snip

		\rsblib_based_on_prei rsb_spmm() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		RSBP_DFE ( spmm<rsb_err_t> (transA,&alpha,nrhs,order,Bp,ldB,&beta,Cp,ldC) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spmm(rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * Bp, const NT beta, NT * Cp) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmm() \rsblib_based_on_post
		*/
		const rsb_nnz_idx_t ldC {}, ldB {};
		RSBP_CDMFHO();
		RSBP_DFE ( spmm<rsb_err_t> (transA,&alpha,nrhs,order,Bp,ldB,&beta,Cp,ldC) );
		RSBP_CWMFRE();
	}

#if RSBP_WANT_IMPLICIT_ORDER
	RSBP_RVT RSBP_DEPRECATED Err_t spmm(rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, const NT * Bp, const NT beta, NT * Cp) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmm() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		bool do_trans = ( transA == RSB_TRANSPOSITION_N ? false : true );
		const rsb_coo_idx_t ldC = do_trans ? cols():rows(), ldB = do_trans ? rows():cols();
		const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
		RSBP_DFE ( spmm<rsb_err_t> (transA,&alpha,nrhs,order,Bp,ldB,&beta,Cp,ldC) );
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_IMPLICIT_ORDER */

#if RSBP_WANT_IMPLICIT_ORDER
	RSBP_RVT RSBP_DEPRECATED Err_t spmm(rsb_trans_t transA, const NT * alphap, rsb_coo_idx_t nrhs, const NT * Bp, const NT * betap, NT * Cp) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmm() \rsblib_based_on_post

		\rl_deprecate
		*/
		RSBP_CDMFHO();
		RSBP_DFE ( spmm<rsb_err_t> (transA,*alphap,nrhs,Bp,*betap,Cp) );
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_IMPLICIT_ORDER */

#if RSBP_WANT_IMPLICIT_ORDER
	RSBP_RVT RSBP_DEPRECATED Err_t spmm(NT * y, const NT * x, rsb_coo_idx_t nrhs, bool do_trans = false) const
	{
		/*!
		\rsblib_based_on_prei rsb_spmm() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		const rsb_trans_t transA = do_trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;
		const NT alpha = defmultalpha;
		const NT beta = defmultbeta;
		const rsb_coo_idx_t ldC = do_trans ? cols():rows(), ldB = do_trans ? rows():cols();
		const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
		RSBP_DFE( spmm<rsb_err_t>(transA, &alpha, nrhs, order, x, ldB, &beta, y, ldC) );
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_IMPLICIT_ORDER */

#if RSBP_WANT_CPP20
	RSBP_RVT Err_t spmm(rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const std::span<const NT> x, const NT beta, const std::span<NT> y) const
	{
		/*!
		\rsblib_based_on_pred rsb_spmm() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		const auto [ldB, ldC] = _get_ldX_ldY(transA,order,nrhs);
		RSBP_CFE ( rsb_spmm (transA,&alpha,this->mtxAp_,nrhs,order,x.data(),ldB,&beta,y.data(),ldC) );
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_CPP20 */

	RSBP_RVT RSBP_DEPRECATED Err_t spsm(rsb_trans_t transT, const NT * alphap, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * betap, const NT * Bp, rsb_nnz_idx_t ldB, NT * Cp, rsb_nnz_idx_t ldC) const
	{
		/*!
		\rsblib_based_on_pred rsb_spsm() \rsblib_based_on_post
		\rl_deprecate
		*/
		RSBP_CWMFHO();
		RSBP_CFE( rsb_spsm (transT,alphap,this->mtxAp_,nrhs,order,betap,Bp,ldB,Cp,ldC) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spsm(rsb_trans_t transT, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT beta, const NT * Bp, rsb_nnz_idx_t ldB, NT * Cp, rsb_nnz_idx_t ldC) const
	{
		/*!
		\rsblib_based_on_prei rsb_spsm() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		RSBP_DFE( spsm<rsb_err_t> (transT,&alpha,nrhs,order,&beta,Bp,ldB,Cp,ldC) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spsm(rsb_trans_t transT, const NT alpha, rsb_coo_idx_t nrhs, const NT * Bp, NT * Cp) const
	{
		/*!
		\rsblib_based_on_prei rsb_spsm() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		const rsb_coo_idx_t ldC = rows(), ldB = ldC;
		const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
		const NT beta = defmultbeta;
		RSBP_DFE( spsm<rsb_err_t> (transT,&alpha,nrhs,order,&beta,Bp,ldB,Cp,ldC) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spsm(NT * y, const NT * x, rsb_coo_idx_t nrhs, bool do_trans = false) const
	{
		/*!
		\rsblib_based_on_prei rsb_spsm() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		const rsb_trans_t transA = do_trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;
		const NT alpha = one;
		const NT beta = zero;
		const rsb_coo_idx_t ldC = rows(), ldB = ldC;
		RSBP_DFE( spsm<rsb_err_t> (transA,&alpha,nrhs,RSB_FLAG_WANT_COLUMN_MAJOR_ORDER,&beta,x,ldB,y,ldC));
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spsm(NT * y, rsb_coo_idx_t nrhs, bool do_trans = false) const
	{
		/*!
		\rsblib_based_on_prei rsb_spsm() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		RSBP_DFE( spsm<rsb_err_t> (y,y,nrhs,do_trans ));
		RSBP_CWMFRE();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t spsv(rsb_trans_t transT, const NT * alphap, const NT * Xp, rsb_coo_idx_t incX, NT * Yp, rsb_coo_idx_t incY) const
	{
		/*!
		\rsblib_based_on_pred rsb_spsv() \rsblib_based_on_post
		\rl_deprecate
		*/
		RSBP_CWMFHO();
		RSBP_CFE( rsb_spsv(transT,alphap,this->mtxAp_,Xp,incX,Yp,incY) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spsv(rsb_trans_t transT, const NT alpha, const NT * Xp, NT * Yp) const
	{
		/*!
		\rsblib_based_on_prei rsb_spsv() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		const rsb_coo_idx_t incX = 1, incY = 1;
		RSBP_CFE( spsv<rsb_err_t>(transT,&alpha,Xp,incX,Yp,incY) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spsv(NT * y, const NT * x, bool do_trans = false) const
	{
		/*!
		\rsblib_based_on_prei rsb_spsv() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		const rsb_trans_t transA = do_trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;
		const NT alpha = one;
		const rsb_coo_idx_t incX = 1, incY = 1;
		RSBP_DFE( spsv<rsb_err_t>(transA,&alpha,x,incX,y,incY) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t spsv(NT * y, bool do_trans = false) const
	{
		/*!
		\rsblib_based_on_prei rsb_spsv() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		const rsb_trans_t transA = do_trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;
		const NT alpha = one;
		const rsb_coo_idx_t incY = 1;
		RSBP_DFE( spsv<rsb_err_t>(transA,&alpha,y,incY,y,incY) );
		RSBP_CWMFRE();
	}

	size_t get_info_size_t(enum rsb_mif_t mif) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		size_t val = 0;
		RSBP_DFE( rsb_mtx_get_info(this->mtxAp_,mif,&val));
		RSBP_MFRV(val);
	}

	rsb_flags_t get_info_rsb_flags_t(enum rsb_mif_t mif) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		rsb_flags_t val = 0;
		RSBP_DFE( rsb_mtx_get_info(this->mtxAp_,mif,&val));
		RSBP_MFRV(val);
	}

	rsb_blk_idx_t get_info_blk_t(enum rsb_mif_t mif) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		rsb_blk_idx_t val = 0;
		RSBP_DFE( rsb_mtx_get_info(this->mtxAp_,mif,&val));
		RSBP_MFRV(val);
	}

	rsb_nnz_idx_t get_info_nnz_t(enum rsb_mif_t mif) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		rsb_nnz_idx_t val = 0;
		RSBP_DFE( rsb_mtx_get_info(this->mtxAp_,mif,&val));
		RSBP_MFRV(val);
	}

	rsb_flags_t get_flags_t(enum rsb_mif_t mif) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		rsb_flags_t val = 0;
		RSBP_DFE( rsb_mtx_get_info(this->mtxAp_,mif,&val));
		RSBP_MFRV(val);
	}

	rsb_type_t get_type_t(enum rsb_mif_t mif) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		rsb_type_t val = 0;
		RSBP_DFE( rsb_mtx_get_info(this->mtxAp_,mif,&val));
		RSBP_MFRV(val);
	}

	rsb_coo_idx_t get_info_coo_t(enum rsb_mif_t mif) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		rsb_coo_idx_t val = 0;
		RSBP_DFE( rsb_mtx_get_info(this->mtxAp_,mif,&val));
		RSBP_MFRV(val);
	}

	size_t _get_index_storage_bytes(void) const
	{
		/*!
		\rl_name_to_be_changed
		*/
		return get_info_size_t(RSB_MIF_INDEX_STORAGE_IN_BYTES__TO__SIZE_T);
	}

	size_t _get_storage_bytes(void) const
	{
		/*!
		\rl_name_to_be_changed
		*/
		return _get_index_storage_bytes() + 2 * sizeof(NT) * nnz();
	}

private:
	void swap(RsbMatrix & other)
	{
		std::swap(this->mtxAp_,other.mtxAp_);
		std::swap(this->state_,other.state_);
	}

	RSBP_RVT Err_t type_check(void)
	{
		RSBP_KWMFHO();
		if ( rsbtype() != get_type_t(RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T) )
		{
			errval_ = RSB_ERR_BADARGS;
			RSBP_PERROR(RSB_ERR_BADARGS);
		}
		RSBP_MDR();
	}

	RSBP_RVT Err_t use_mtx_ptr(struct rsb_mtx_t * mtxBp)
	{
		RSBP_KWMFHO();
		if( mtxBp != RSBP_NULL && mtxBp != this->mtxAp_ )
		{
   			rsb_mtx_free (this->mtxAp_);
			this->mtxAp_ = mtxBp;
		}
		errval_ = type_check<rsb_err_t>();
		RSBP_MDR();
	}
public:

	rsb_nnz_idx_t nnz(void) const
	{
		/*!
		*/
		RSBP_CNMFHO();
		RSBP_MFRV(get_info_nnz_t(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T));
	}

	rsb_blk_idx_t blocks(void) const
	{
		/*!
		*/
		RSBP_CNMFHO();
		RSBP_MFRV(get_info_blk_t(RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T));
	}

	rsb_coo_idx_t rows(void) const
	{
		/*!
		\rsblib_dims_snip
		*/
		RSBP_CNMFHO();
		RSBP_MFRV(get_info_coo_t(RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T));
	}

	rsb_coo_idx_t cols(void) const
	{
		/*!
		\rsblib_dims_snip
		*/
		RSBP_CNMFHO();
		RSBP_MFRV(get_info_coo_t(RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T));
	}

	RSBP_RVT Err_t get_vals(NT* VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_vals() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		RSBP_CFE(rsb_mtx_get_vals(this->mtxAp_, VA, IA, JA, nnz, flags));
		RSBP_MDR();
	}

	NT get_val(const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_flags_t flags = RSB_FLAG_NOFLAGS) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_vals() \rsblib_based_on_post
		*/
		NT val {};
		rsb_mtx_get_vals(this->mtxAp_, &val, &i, &j, 1, flags);
		return val;
	}

	RSBP_RVT Err_t set_val(const NT val, const rsb_coo_idx_t i, const rsb_coo_idx_t j, rsb_flags_t flags = RSB_FLAG_NOFLAGS)
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_set_vals() \rsblib_based_on_post

		\rsblib_asm_snip

		\sa \ref set_vals().
		*/
		RSBP_KWMFHO();
		RSBP_NFE(rsb_mtx_set_vals (this->mtxAp_, &val, &i, &j, 1, flags));
		RSBP_MDR();
	}

	RSBP_RVT Err_t set_vals(const NT * VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_flags_t flags)
	{
		/*!
			Add a single entry during the assembly of a matrix created empty.
			\n
			Use close() to terminate matrix assembly.

			\rsblib_asm_snip

			\rsblib_based_on_pred rsb_mtx_set_vals() \rsblib_based_on_post
		*/
		RSBP_KWMFHO();
		RSBP_NFE(rsb_mtx_set_vals (this->mtxAp_, VA, IA, JA, nnz, flags));
		RSBP_MDR();
	}

	RSBP_RVT Err_t get_vec(NT * Dp, enum rsb_extff_t flags) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_vec() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		RSBP_CFE(rsb_mtx_get_vec(this->mtxAp_, Dp, flags));
		RSBP_MDR();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t get_coo(rsb_trans_t transA, NT * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_flags_t flags ) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_coo() \rsblib_based_on_post

		\rl_no_transa
		\rl_deprecate
		*/
		RSBP_CWMFHO();

		if(transA != RSB_TRANSPOSITION_N)
			RSBP_PERROR(RSB_ERR_UNIMPLEMENTED_YET);
		RSBP_CFE(rsb_mtx_get_coo(this->mtxAp_, VA, IA, JA, flags));
		RSBP_CWMFRE();
	}

#if RSBP_WANT_CPP20
	RSBP_RVT Err_t get_coo(rsb_trans_t transA, const std::span<NT> VA, const std::span<rsb_coo_idx_t> IA, const std::span<rsb_coo_idx_t> JA, rsb_flags_t flags ) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_coo() \rsblib_based_on_post

		\rl_no_transa
		*/
		RSBP_CWMFHO();

		if(transA != RSB_TRANSPOSITION_N)
			RSBP_PERROR(RSB_ERR_UNIMPLEMENTED_YET);
		RSBP_CFE(rsb_mtx_get_coo(this->mtxAp_, VA.data(), IA.data(), JA.data(), flags));
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_CPP20 */

	RSBP_RVT RSBP_DEPRECATED Err_t get_csr(rsb_trans_t transA, NT * VA, rsb_coo_idx_t * RP, rsb_coo_idx_t * JA, rsb_flags_t flags ) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_csr() \rsblib_based_on_post

		\rl_no_transa
		\rl_deprecate
		*/
		RSBP_CWMFHO();

		if(transA != RSB_TRANSPOSITION_N)
		{
			errval = RSB_ERR_UNIMPLEMENTED_YET;
			RSBP_PERROR(errval);
		}
		else
		{
			RSBP_CFE(rsb_mtx_get_csr(RSB_NUMERICAL_TYPE_SAME_TYPE, this->mtxAp_, VA, RP, JA, flags));
		}
		RSBP_CWMFRE();
	}

#if RSBP_WANT_CPP20
	RSBP_RVT RSBP_DEPRECATED Err_t get_csr(rsb_trans_t transA, const std::span<NT> VA, const std::span<rsb_coo_idx_t> RP, const std::span<rsb_coo_idx_t> JA, rsb_flags_t flags ) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_csr() \rsblib_based_on_post

		\rl_no_transa
		*/
		RSBP_CWMFHO();

		if(transA != RSB_TRANSPOSITION_N)
		{
			errval = RSB_ERR_UNIMPLEMENTED_YET;
			RSBP_PERROR(errval);
		}
		else
		{
			RSBP_CFE(rsb_mtx_get_csr(RSB_NUMERICAL_TYPE_SAME_TYPE, this->mtxAp_, VA.data(), RP.data(), JA.data(), flags));
		}
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_CPP20 */

	RSBP_RVT RSBP_DEPRECATED Err_t get_rows_sparse(rsb_trans_t transA, const NT * alphap, NT * VA, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, rsb_coo_idx_t frA, rsb_coo_idx_t lrA, rsb_nnz_idx_t *rnzp, rsb_flags_t flags ) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_rows_sparse() \rsblib_based_on_post

		\rl_deprecate
		*/
		RSBP_CWMFHO();
		RSBP_CFE(rsb_mtx_get_rows_sparse(transA, alphap, this->mtxAp_, VA, IA, JA, frA, lrA, rnzp, flags));
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t upd_vals(enum rsb_elopf_t elop_flags, const NT & omega)
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_upd_vals() \rsblib_based_on_post

		\rsblib_upd_vals_snip
		*/
		RSBP_NWMFHO();
		RSBP_NFE(rsb_mtx_upd_vals(this->mtxAp_, elop_flags, &omega));
		RSBP_MDR();
	}

	RSBP_RVT Err_t upd_vals(enum rsb_elopf_t elop_flags, const NT * omegap)
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_upd_vals() \rsblib_based_on_post

		\rl_deprecate
		*/
		RSBP_NWMFHO();
		RSBP_NFE(rsb_mtx_upd_vals(this->mtxAp_, elop_flags, omegap));
		RSBP_MDR();
	}

	RSBP_RVT Err_t get_nrm(NT * Np, enum rsb_extff_t flags) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_nrm() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		RSBP_CFE(rsb_mtx_get_nrm(this->mtxAp_, Np, flags));
		RSBP_MDR();
	}

	rsb_type_t rsbtype(void) const
	{
		/*!
		\rsblib_rsbtype_snip
		*/
		RSBP_CNMFHO();
		RSBP_MFRV(typecode_);
	}

	rsb_flags_t rsbflags(void) const
       	{
		/*!
		*/
		RSBP_CNMFHO();
		RSBP_MFRV(get_flags_t(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T));
	}

	rsb_string_t get_info_str(const char * key) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info_str() \rsblib_based_on_post
		*/
		RSBP_CDMFHO();
		const auto RSBP_INFOBUF { 256 };
		char ss[RSBP_INFOBUF];
		RSBP_CFE(rsb_mtx_get_info_str(this->mtxAp_,key,ss,RSBP_INFOBUF));
		RSBP_MFRV(rsb_string_t(ss));
	}

	RSBP_RVT Err_t get_info(enum rsb_mif_t miflags, void* minfop)const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_get_info() \rsblib_based_on_post
		*/
		RSBP_CWMFHO();
		RSBP_CFE(rsb_mtx_get_info(this->mtxAp_, miflags, minfop));
		RSBP_MDR();
	}

	rsb_string_t _info(void) const
	{
		/*!
		\rl_name_to_be_changed
		*/
		RSBP_CNMFHO();
		RSBP_MFRV(get_info_str("RSB_MIF_MATRIX_INFO__TO__CHAR_P"));
	}

	RSBP_RVT RSBP_DEPRECATED Err_t tune_spsm_threads(/*RsbMatrix * mtxp,*/ rsb_real_t *sfp=RSBP_NULL, rsb_int_t *tnp=RSBP_NULL, rsb_int_t maxr=0, rsb_time_t maxt=0, rsb_trans_t transA=RSB_TRANSPOSITION_N, const NT * alphap=RSBP_NULL, rsb_coo_idx_t nrhs=1, rsb_flags_t order=RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, const NT * Bp=RSBP_NULL, rsb_nnz_idx_t ldB=0, const NT * betap=RSBP_NULL, NT * Cp=RSBP_NULL, rsb_nnz_idx_t ldC=0) const
	{
		/*!
		\rsblib_based_on_pred rsb_tune_spsm() \rsblib_based_on_post

		\rl_deprecate
		*/
		RSBP_CWMFHO();

		RSBP_CFE(rsb_tune_spsm(RSBP_NULL, sfp, tnp, maxr, maxt, transA, alphap, this->mtxAp_, nrhs, order, Bp, ldB, betap, Cp, ldC));
		RSBP_MDR();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t tune_spmm_threads(/*RsbMatrix * mtxp,*/ rsb_real_t *sfp=RSBP_NULL, rsb_int_t *tnp=RSBP_NULL, rsb_int_t maxr=0, rsb_time_t maxt=0, rsb_trans_t transA=RSB_TRANSPOSITION_N, const NT * alphap=RSBP_NULL, rsb_coo_idx_t nrhs=1, rsb_flags_t order=RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, const NT * Bp=RSBP_NULL, rsb_nnz_idx_t ldB=0, const NT * betap=RSBP_NULL, NT * Cp=RSBP_NULL, rsb_nnz_idx_t ldC=0) const
	{
		/*!
		\rsblib_based_on_pred rsb_tune_spmm() \rsblib_based_on_post

		\rsblib_autotune_snip
		\rl_deprecate
		*/
		RSBP_CWMFHO();

		RSBP_CFE(rsb_tune_spmm(RSBP_NULL, sfp, tnp, maxr, maxt, transA, alphap, this->mtxAp_, nrhs, order, Bp, ldB, betap, Cp, ldC));
		RSBP_MDR();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t tune_spmm(/*RsbMatrix * mtxp,*/ rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const NT * alphap, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * Bp, rsb_nnz_idx_t ldB, const NT * betap, NT * Cp, rsb_nnz_idx_t ldC)
	{
		/*!
		\rsblib_based_on_pred rsb_tune_spmm() \rsblib_based_on_post

		\rsblib_autotune_snip

		\rl_deprecate
		*/
		RSBP_NWMFHO();

		if(true)
		{
			struct rsb_mtx_t * mtxBp {};
			RSBP_NFE(rsb_tune_spmm(&mtxBp, sfp, tnp, maxr, maxt, transA, alphap, this->mtxAp_, nrhs, order, Bp, ldB, betap, Cp, ldC));
			errval_ = use_mtx_ptr<rsb_err_t>(mtxBp);
		}
		else
		{
			// Note: case not ready yet; may want a new member function for this.
			RSBP_NFE(rsb_tune_spmm(&(this->mtxAp_), sfp, tnp, maxr, maxt, transA, alphap, RSBP_NULL, nrhs, order, Bp, ldB, betap, Cp, ldC));
		}
		RSBP_MDR();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t tune_spmm(/*RsbMatrix * mtxp,*/ rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * Bp, rsb_nnz_idx_t ldB, const NT beta, NT * Cp, rsb_nnz_idx_t ldC)
	{
		/*!
		\rsblib_based_on_prei rsb_tune_spmm() \rsblib_based_on_post

		\rsblib_autotune_snip
		*/

		RSBP_CWMFHO();
		RSBP_CFE( tune_spmm<rsb_err_t>(/*RsbMatrix * mtxp,*/ sfp, tnp, maxr, maxt, transA, &alpha, nrhs, order, Bp, ldB, &beta, Cp, ldC) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t tune_spmm(rsb_real_t & sf, rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * Bp, const NT beta, NT * Cp)
	{
		/*!
		\rsblib_based_on_prei rsb_tune_spmm() \rsblib_based_on_post

		\rsblib_autotune_snip
		*/
		rsb_int_t *tnp {nullptr};
		const rsb_int_t maxr {0};
		const rsb_time_t maxt {0.0};
		const std::pair<rsb_nnz_idx_t,rsb_nnz_idx_t> ldX = _get_ldX_ldY(transA,order,nrhs);
# if defined(__cplusplus) && (__cplusplus>=201709L)
		const auto [ldB, ldC] = ldX;
#else
		const rsb_nnz_idx_t ldB = ldX.first, ldC = ldX.second;
#endif

		RSBP_CWMFHO();
		RSBP_CFE( tune_spmm<rsb_err_t>(/*RsbMatrix * mtxp,*/ &sf, tnp, maxr, maxt, transA, &alpha, nrhs, order, Bp, ldB, &beta, Cp, ldC) );
		RSBP_CWMFRE();
	}

#if RSBP_WANT_CPP20
	RSBP_RVT RSBP_DEPRECATED Err_t tune_spmm(rsb_real_t & sf, rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const std::span<const NT> B, rsb_nnz_idx_t ldB, const NT beta, const std::span<NT> C, rsb_nnz_idx_t ldC)
	{
		/*!
		\rsblib_based_on_prei rsb_tune_spmm() \rsblib_based_on_post

		\rsblib_autotune_snip
		*/
		rsb_int_t *tnp {nullptr};
		const rsb_int_t maxr {0};
		const rsb_time_t maxt {0.0};
		RSBP_CWMFHO();
		RSBP_CFE( tune_spmm<rsb_err_t>(/*RsbMatrix * mtxp,*/ &sf, tnp, maxr, maxt, transA, &alpha, nrhs, order, B.data(), ldB, &beta, C.data(), ldC) );
		RSBP_CWMFRE();
	}

	RSBP_RVT Err_t tune_spmm(rsb_real_t & sf, rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const std::span<const NT> B, const NT beta, const std::span<NT> C)
	{
		/*!
		\rsblib_based_on_pred rsb_tune_spmm() \rsblib_based_on_post

		\rsblib_autotune_snip
		*/
		rsb_int_t *tnp {nullptr};
		const rsb_int_t maxr {0};
		const rsb_time_t maxt {0.0};
		const auto [ldB, ldC] = _get_ldX_ldY(transA,order,nrhs);
		RSBP_NWMFHO();
		struct rsb_mtx_t * mtxBp {};
		RSBP_NFE(rsb_tune_spmm(&mtxBp, &sf, tnp, maxr, maxt, transA, &alpha, this->mtxAp_, nrhs, order, B.data(), ldB, &beta, C.data(), ldC));
		errval_ = use_mtx_ptr<rsb_err_t>(mtxBp);
		RSBP_MDR();
	}

	RSBP_RVT RSBP_DEPRECATED Err_t tune_spmm(/*RsbMatrix * mtxp,*/ rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const std::span<const NT> B, rsb_nnz_idx_t ldB, const NT beta, const std::span<NT> C, rsb_nnz_idx_t ldC)
	{
		/*!
		\rsblib_based_on_prei rsb_tune_spmm() \rsblib_based_on_post

		\rsblib_autotune_snip
		*/

		RSBP_CWMFHO();
		RSBP_CFE( tune_spmm<rsb_err_t>(/*RsbMatrix * mtxp,*/ sfp, tnp, maxr, maxt, transA, &alpha, nrhs, order, B.data(), ldB, &beta, C.data(), ldC) );
		RSBP_CWMFRE();
	}
#endif /* RSBP_WANT_CPP20 */

	RSBP_RVT RSBP_DEPRECATED Err_t tune_spsm(/*RsbMatrix * mtxp,*/ rsb_real_t *sfp, rsb_int_t *tnp, rsb_int_t maxr, rsb_time_t maxt, rsb_trans_t transA, const NT * alphap, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT * Bp, rsb_nnz_idx_t ldB, const NT * betap, NT * Cp, rsb_nnz_idx_t ldC)
	{
		/*!
		\rsblib_based_on_pred rsb_tune_spsm() \rsblib_based_on_post

		\rl_deprecate
		*/
		RSBP_NWMFHO();

		if(true)
		{
			struct rsb_mtx_t * mtxBp {};
			RSBP_NFE(rsb_tune_spsm(&mtxBp, sfp, tnp, maxr, maxt, transA, alphap, this->mtxAp_, nrhs, order, Bp, ldB, betap, Cp, ldC));
			use_mtx_ptr(mtxBp);
		}
		else
		{
			// Note: case not ready yet; may want a new member function for this.
			RSBP_NFE(rsb_tune_spsm(&(this->mtxAp_), sfp, tnp, maxr, maxt, transA, alphap, RSBP_NULL, nrhs, order, Bp, ldB, betap, Cp, ldC));
		}
		RSBP_MDR();
	}

	RSBP_RVT Err_t file_save(const rsb_char_t * filename=RSBP_NULL) const
	{
		/*!
		\rsblib_based_on_pred rsb_file_mtx_save() \rsblib_based_on_post

		Example snip from examples/misc.cpp: \snippet examples/misc.cpp snip__RsbMatrix_file_save_stdout
		*/
		RSBP_CWMFHO();
		RSBP_CFE(rsb_file_mtx_save(this->mtxAp_, filename));
		RSBP_MDR();
	}

private:
	RSBP_RVT RSBP_DEPRECATED Err_t clone(struct rsb_mtx_t ** mtxBpp/*, rsb_type_t typecode*/, rsb_trans_t transA = RSB_TRANSPOSITION_N, const NT *alphap = RSBP_NULL, rsb_flags_t flags = RSB_FLAG_NOFLAGS) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_clone() \rsblib_based_on_post

		\rl_deprecate
		*/
		// Note: ignoring typecode argument for now.
		RSBP_CWMFHO();
		if(!mtxBpp)
			RSBP_PERROR(RSB_ERR_BADARGS);
		RSBP_CFE(rsb_mtx_clone(mtxBpp, this->typecode_, transA, alphap, this->mtxAp_, flags));
		RSBP_MDR();
	}
public:

	RsbMatrix & operator=(const RsbMatrix & A_Rsb)
	{
		/*!
		A copy constructor.
		Will clone the input matrix contents.
		*/
		RSBP_NWMFHO();
		errval_ = A_Rsb.clone<rsb_err_t>(&this->mtxAp_/*, rsb_type_t typecode*/);
		RSBP_MFRV(*this);
	}

	bool _is_complex(void) const
	{
		/*!
		\rl_name_to_be_changed
		*/
		return !std::is_scalar<NT>();
	}

private:
	bool _equal_to(const RsbMatrix & B_Rsb) const
	{
		/*!
		\rl_name_to_be_changed
		*/
		// Note: space-inefficient.
		// An iterator or a submatrices-descent based method would be more space efficient.
		RSBP_CNMFHO();
		bool retval = false;
		const rsb_coo_idx_t invalid_coo_ = RSB_INVALID_COO_IDX_VAL;
		const rsb_nnz_idx_t nnzA = nnz();

		if( nnzA != B_Rsb.nnz() )
			goto ret;
		if( this->rows() != B_Rsb.rows() )
			goto ret;
		if( this->cols() != B_Rsb.cols() )
			goto ret;

		if( ! _is_integer() )
		{
			if( this->normOne() != B_Rsb.normOne() )
				goto ret;
			if( this->normInf() != B_Rsb.normInf() )
				goto ret;
		}

		{
			std::vector<rsb_coo_idx_t > AIA(nnzA,invalid_coo_),AJA(nnzA,invalid_coo_);
			std::vector<rsb_coo_idx_t > BIA(nnzA,invalid_coo_),BJA(nnzA,invalid_coo_);
			std::vector<NT > AVA(nnzA,invalid_coo_);
			std::vector<NT > BVA(nnzA,invalid_coo_);
			rsb_err_t errval;

			errval = rsb_mtx_get_coo(this->mtxAp_, &AVA[0], &AIA[0], &AJA[0], RSB_FLAG_SORTED_INPUT);
			RSBP_ERRGOTO(errval,ret)
			errval = rsb_mtx_get_coo(B_Rsb.mtxAp_, &BVA[0], &BIA[0], &BJA[0], RSB_FLAG_SORTED_INPUT);
			RSBP_ERRGOTO(errval,ret)
			retval = (
				std::equal(AIA.begin(),AIA.end(),BIA.begin()) &&
				std::equal(AJA.begin(),AJA.end(),BJA.begin()) &&
				std::equal(AVA.begin(),AVA.end(),BVA.begin()) );
		}
ret:
		RSBP_MFRV(retval);
	} /* _equal_to */

public:
	bool operator==(const RsbMatrix & B_Rsb) const
	{
		/*!
		Deep comparison: compare if the two matrices have same dimensions, nonzeroes count, nonzeroes pattern and value.
		Meant for very sporadic use.
		Inefficient: it can involve matrices copying.

		\rsblib__cmp_snip
		*/
		RSBP_CNMFHO();
		const bool retval = _equal_to(B_Rsb);
		RSBP_MFRV(retval);
	}

	bool operator!=(const RsbMatrix & B_Rsb) const
	{
		/*!
			\rsblib__cmp_snip

			\see RsbMatrix::operator==(const RsbMatrix & B_Rsb) const;
		*/
		return ! ( (*this) == B_Rsb );
	}

	NT normOne(void) const
	{
		/*!
		*/
		NT nrm = RSBP_NUMT_INVALID;
		this->get_nrm(&nrm, RSB_EXTF_NORM_ONE);
		return nrm;
	}

	NT normInf(void) const
	{
		/*!
		*/
		NT nrm = RSBP_NUMT_INVALID;
		this->get_nrm(&nrm, RSB_EXTF_NORM_INF);
		return nrm;
	}

	RSBP_RVT Err_t rndr(const rsb_char_t * filename=RSBP_NULL, rsb_coo_idx_t pmWidth=512, rsb_coo_idx_t pmHeight=512, rsb_marf_t rflags=RSB_MARF_EPS) const
	{
		/*!
		\rsblib_based_on_pred rsb_mtx_rndr() \rsblib_based_on_post

		\rsblib_rndr_snip
		*/
		RSBP_CWMFHO();
		RSBP_CFE(rsb_mtx_rndr(filename, mtxAp_, pmWidth, pmHeight, rflags));
		RSBP_MDR();
	}

}; /* class RsbMatrix */
#ifndef RSBP_NO_NAMESPACE
} // namespace rsb
#endif /* RSBP_NO_NAMESPACE */

#endif	/* RSB_RSB_HPP_INCLUDED */
