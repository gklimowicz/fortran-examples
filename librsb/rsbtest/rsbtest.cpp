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
#include "config.h"
#if defined(__cplusplus) && (__cplusplus>= 201703L) 
#define USE_CXX17 1
#else
#define USE_CXX17 0
#endif
#ifdef HAVE_RSB_H
#include <rsb.h>
#include <rsb-config.h>
#if 0
#if defined(RSB_HAVE_IHI) && (RSB_HAVE_IHI)
#include <rsb-librsb-internals.h>
#endif /* RSB_HAVE_IHI */
#endif /* 0 */
#else /* HAVE_RSB_H */
#error
#endif /* HAVE_RSB_H */
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#endif /* _OPENMP */
#if HAVE_GETOPT_H
#include <getopt.h>
#endif
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <algorithm> // min_element, max_element
#include <cmath>
#include <functional> // std::function
#include <unistd.h>
#include <random>
#include <type_traits>
#include <thread>
#include <regex>
#include <cassert>
#ifdef HAVE_MUTEX_H
#include <mutex>
#endif /* HAVE_MUTEX_H */
#include <fstream>
#if HAVE_FILESYSTEM && USE_CXX17
#include <filesystem>	// c++17 only
#endif /* HAVE_FILESYSTEM && USE_CXX17 */
#if defined(HAVE_RSB_HPP) && USE_CXX17
#define RSBT_WANT_RSB_HPP 1
#endif
#ifdef RSBT_WANT_RSB_HPP
#define RSBT_TESTING_ONLY 1 /* activates test features */
#define RSBP_WANT_REV 1 /* Necessary to activate rsb_err_t as default return type. */
#include <rsb.hpp>
#endif /* RSBT_WANT_RSB_HPP */
#if defined(WANT_MKL)
#if defined(HAVE_MKL_GET_VERSION) && defined(HAVE_MKL_MKL_H)
#define RSBT_USE_MKL
#include <mkl/mkl.h>
#endif /* */
#if defined(HAVE_MKL_GET_VERSION) && defined(HAVE_MKL_MKL_SPBLAS_H)
#define RSBT_USE_MKL_SPBLAS
#include <mkl/mkl_spblas.h> // want_mkl
#endif /* HAVE_MKL_MKL_SPBLAS_H */
#endif /* WANT_MKL */
#ifdef HAVE_RSBPP_HPP
#if RSBT_WANT_MINIMAL_LIBRSBPP
#undef HAVE_RSBPP_HPP
#endif
#define RSBPP_HAS_RSB_H 1 /*force <rsbpp.hpp> to include <rsb.h> (and avoid half redeclaration)*/
#include <rsbpp.hpp>
#endif /* HAVE_RSBPP_HPP */
#ifdef HAVE_SYS_UTSNAME_H
#include <sys/utsname.h>	/* uname */
#endif /* HAVE_SYS_UTSNAME_H */
#include <cstdlib> // EXIT_SUCCESS
#define RSB_XSTRINGIFY(X) RSB_STRINGIFY(X)
#define RSB_STRINGIFY(X) #X
#define RSBT_INFOBUF 256
#define RSBT_ETOL 1e-4
#define RSBT_ERR_MISMATCH RSB_ERR_CAST(0x15BADBAD)
#define RSBT_FAIL_REC -1 /* recover failed case */
#define RSBT_SAVE_VEC(VFN,NT,VP,VL) if(want_dump_failed){ rsb_type_t typecode=rsb_type_t_for<NT>(); std::string fn = VFN; rsb_file_vec_save(fn.c_str(),typecode,VP,VL); std::cout << "Saved to vector " << fn << std::endl; }
#define RSBT_SAVE_MAT(MFN,MTX) if(want_dump_failed){ std::string fn = MFN; rsb_file_mtx_save(MTX,fn.c_str()); std::cout << "Saved matrix to " << fn << std::endl; }
#define RSBT_LOAD_VEC(VFN, TC, V) { rsb_type_t typecode=rsb_type_t_for<NT>(); rsb_coo_idx_t yvl; rsb_file_vec_load(VFN, typecode, nullptr, &yvl); V.resize(yvl); rsb_file_vec_load(VFN, typecode, V.data(), &yvl); }

#define report_at_line(X) { std::cout << "@" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << std::endl; }
#define return_at_line(X) { const rsb_err_t ral_errval {X}; if (ral_errval!=RSB_ERR_NO_ERROR ) report_at_line(X); return ral_errval; }

#ifdef HAVE_RSBPP_HPP
//#define RSBPP_MATRIX_T_ORDERING RsbPP_Matrix_T<NT>::OrderingRSB // not ready yet
#define RSBPP_MATRIX_T_ORDERING RsbPP_Matrix_T<NT>::OrderingCOR
#endif /* HAVE_RSBPP_HPP */
const std::string ltx_preamble = R"LaTeXDelimiter(
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[xcolor=dvipsnames]{beamer}
%rm *.pdf *.aux *.nav *.out *.dvi *.log *.snm *.toc -f && ls && latex slides.tex && pdflatex slides.tex && pdflatex slides.tex   && evince slides.pdf
\usepackage{color}
\usepackage{sverb} % for verbinput
\usepackage{hyperref} % for url
\usepackage{epsfig} 
\hypersetup{pdfpagemode=fullscreen}
\setbeamertemplate{footline}[page number] % no footline
\setbeamertemplate{footline}[text line{...}] % no footline
\setbeamertemplate{navigation symbols}{}
\usepackage[]{verbatim}
\usepackage[]{ulem}
\usepackage{graphicx}
\usepackage{ucs}
\usepackage[utf8x]{inputenc}
\usepackage{epstopdf}
\usepackage[]{comment}
\usepackage[]{color}
%\logo{\includegraphics[scale=0.20]{logo.png}}
\def\mytitle{librsb autotuning benchmark report\\
EXTRA-TITLE}%
\title{\mytitle}
%\author[{michele DOT martone AT lrz DOT de}]{Michele MARTONE}
%\institute{LRZ\\
%Garching bei Muenchen, Germany
%}
%
%\date[]{Place\\Date}%
%%%%%
\def\plotwidth{0.6}
\begin{document}
\section{Title Page}
\frame { \titlepage }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)LaTeXDelimiter";

const std::string ltx_pgtmplt = R"LaTeXDelimiter(
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{frame}[fragile]
\frametitle{FRAME-TITLE}
\begin{figure} \begin{center}
\begin{tabular}{ c c  } {\resizebox{\plotwidth\columnwidth}{!}{\includegraphics{EPS-PLOT-FILE}}} \end{tabular}
\begin{small} \caption{CAPTION} \end{small}
\end{center} \end{figure}%
\end{frame}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%)LaTeXDelimiter";

const std::string ltx_fgtmplt = R"LaTeXDelimiter(
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{frame}[fragile]
\frametitle{FRAME-TITLE}
FRAME-TEXT
\end{frame}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%)LaTeXDelimiter";

const std::string ltx_epilogue = R"LaTeXDelimiter(
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\end{document}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%)LaTeXDelimiter";

static int re2ce(rsb_err_t errval)
{
	return errval == RSB_ERR_NO_ERROR ? EXIT_SUCCESS : EXIT_FAILURE;
}

std::string library_version_string(void)
{
	std::ostringstream lvs;
	lvs << "Linked to " << RSB_PACKAGE << "-" RSB_PACKAGE_VERSION << std::endl;
	lvs << " compiled by " << RSB_XSTRINGIFY(RSB_CC) << std::endl;
	lvs << " with flags " << RSB_XSTRINGIFY(RSB_CFLAGS) << std::endl;
	lvs << " supporting up to " << RSB_XSTRINGIFY(RSB_CONST_MAX_SUPPORTED_THREADS) << " threads" << std::endl;
	lvs << " and hardcoded cache parameters " << RSB_XSTRINGIFY(RSB_DETECTED_MEM_HIERARCHY_INFO) << "" << std::endl;
	if(RSB_WANT_OMP_RECURSIVE_KERNELS)
		lvs << " and OpenMP enabled.";
	else
		lvs << " and OpenMP disabled.";
	lvs  << std::endl;
	lvs  << " Supported types: " << RSB_XSTRINGIFY(RSB_M4_MATRIX_TYPES_STRING) << std::endl;
#ifdef HAVE_RSBPP_HPP
	lvs <<     " Using <rsbpp.hpp> C++ RSB implementation." << std::endl;
#else /* HAVE_RSBPP_HPP */
	lvs << " Not using <rsbpp.hpp> C++ RSB implementation." << std::endl;
#endif /* HAVE_RSBPP_HPP */
#ifdef RSBT_WANT_RSB_HPP
	lvs <<     " Using <rsb.hpp> C++ wrapper to librsb." << std::endl;
#else /* RSBT_WANT_RSB_HPP */
	lvs << " Not using <rsb.hpp> C++ wrapper to librsb." << std::endl;
#endif /* RSBT_WANT_RSB_HPP */
#ifdef RSBT_USE_MKL
	lvs  << " Linked to the Intel MKL: ";
#ifdef mkl_get_version
	{
		MKLVersion mv;
		mkl_get_version(&mv);
		lvs << mv.MajorVersion << "." << mv.MinorVersion << "-" << mv.UpdateVersion << ", ";
		lvs << mv.ProductStatus << ", " << mv.Build << ", " << mv.Processor << ", " << mv.Platform;
	}
#endif /* RSBT_USE_MKL_SPBLAS */
#else /* RSBT_USE_MKL */
	lvs  << " Not linked to the Intel MKL.";
#endif /* RSBT_USE_MKL */
	lvs << std::endl;
#ifdef __cplusplus
	lvs  << " __cplusplus: " << __cplusplus << std::endl;
#endif
#ifdef __GNUC__
	lvs  << " __GNUC__: " << __GNUC__ << std::endl;
#endif
#ifdef __GNUC_MINOR__
	lvs  << " __GNUC_MINOR__: " << __GNUC_MINOR__ << std::endl;
#endif
	return lvs.str();
}

#ifdef RSBT_USE_MKL_SPBLAS
template <class NT> sparse_status_t rsb_mkl_sparse_X_create_coo( sparse_matrix_t        *A, sparse_index_base_t    indexing, MKL_INT    rows, MKL_INT    cols, MKL_INT    nnz, MKL_INT    *row_indx, MKL_INT    *col_indx, NT     *values ) { return SPARSE_STATUS_NOT_SUPPORTED; }
sparse_status_t rsb_mkl_sparse_X_create_coo( sparse_matrix_t        *A, sparse_index_base_t    indexing, MKL_INT    rows, MKL_INT    cols, MKL_INT    nnz, MKL_INT    *row_indx, MKL_INT    *col_indx, double *values ) { return mkl_sparse_d_create_coo( A, indexing, rows, cols, nnz, row_indx, col_indx, values ); }
sparse_status_t rsb_mkl_sparse_X_create_coo( sparse_matrix_t        *A, sparse_index_base_t    indexing, MKL_INT    rows, MKL_INT    cols, MKL_INT    nnz, MKL_INT    *row_indx, MKL_INT    *col_indx, float *values ) { return mkl_sparse_s_create_coo( A, indexing, rows, cols, nnz, row_indx, col_indx, values ); }
sparse_status_t rsb_mkl_sparse_X_create_coo( sparse_matrix_t        *A, sparse_index_base_t    indexing, MKL_INT    rows, MKL_INT    cols, MKL_INT    nnz, MKL_INT    *row_indx, MKL_INT    *col_indx, std::complex<double> *values ) { return mkl_sparse_z_create_coo( A, indexing, rows, cols, nnz, row_indx, col_indx, reinterpret_cast<MKL_Complex16*>(values) ); }
sparse_status_t rsb_mkl_sparse_X_create_coo( sparse_matrix_t        *A, sparse_index_base_t    indexing, MKL_INT    rows, MKL_INT    cols, MKL_INT    nnz, MKL_INT    *row_indx, MKL_INT    *col_indx, std::complex<float> *values ) { return mkl_sparse_c_create_coo( A, indexing, rows, cols, nnz, row_indx, col_indx, reinterpret_cast<MKL_Complex8*>(values) ); }

template <class NT> sparse_status_t rsb_mkl_sparse_X_mv ( sparse_operation_t operation, NT alpha, const  sparse_matrix_t A, struct matrix_descr descr, const  NT *x, NT beta, NT *y ) { return SPARSE_STATUS_NOT_SUPPORTED; }
sparse_status_t rsb_mkl_sparse_X_mv ( sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const  double *x, double beta, double*y ) { return mkl_sparse_d_mv(operation, alpha, A, descr, x, beta, y); }
sparse_status_t rsb_mkl_sparse_X_mv ( sparse_operation_t operation, float alpha, const sparse_matrix_t A, struct matrix_descr descr, const  float *x, float beta, float*y ) { return mkl_sparse_s_mv(operation, alpha, A, descr, x, beta, y); }
sparse_status_t rsb_mkl_sparse_X_mv ( sparse_operation_t operation, std::complex<float> alpha, const sparse_matrix_t A, struct matrix_descr descr, const  std::complex<float> *x, std::complex<float> beta, std::complex<float>*y ) { MKL_Complex8 malpha {alpha.real(),alpha.imag()},mbeta{beta.real(),beta.imag()}; return mkl_sparse_c_mv(operation, malpha, A, descr, reinterpret_cast<const MKL_Complex8*>(x), mbeta, reinterpret_cast<MKL_Complex8*>(y)); }
sparse_status_t rsb_mkl_sparse_X_mv ( sparse_operation_t operation, std::complex<double> alpha, const sparse_matrix_t A, struct matrix_descr descr, const  std::complex<double> *x, std::complex<double> beta, std::complex<double>*y ) { MKL_Complex16 malpha{alpha.real(),alpha.imag()},mbeta{beta.real(),beta.imag()}; return mkl_sparse_z_mv(operation, malpha, A, descr, reinterpret_cast<const MKL_Complex16*>(x), mbeta, reinterpret_cast<MKL_Complex16*>(y)); }
#else /* RSBT_USE_MKL_SPBLAS */
#endif /* RSBT_USE_MKL_SPBLAS */

template<typename NT>
class rsb_mtx_lw_t
{
	private:
	template<typename RT>
	RT get_mif(enum rsb_mif_t miflags)const
	{
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		RT qty;
		errval = this->get_info(miflags, &qty);
		assert(errval == RSB_ERR_NO_ERROR);
		return qty;
	}
	public:
	struct rsb_mtx_t*mtxAp_;
	rsb_mtx_lw_t(struct rsb_mtx_t *mtxAp):mtxAp_(mtxAp){}
	auto blocks(void)const
	{
		return get_mif<rsb_blk_idx_t>(RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T);
	}
	auto nnz(void)const
	{
		return get_mif<rsb_nnz_idx_t>(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T);
	}
	auto rows(void)const
	{
		return get_mif<rsb_coo_idx_t>(RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T);
	}
	auto cols(void)const
	{
		return get_mif<rsb_coo_idx_t>(RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T);
	}
	auto rsbflags(void)const
	{
		return get_mif<rsb_flags_t>(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T);
	}
	rsb_err_t spmv(rsb_trans_t transA, const NT *alphap, const NT * Xp, rsb_coo_idx_t incX, const NT * betap, NT * Yp, rsb_coo_idx_t incY)const
	{
		return rsb_spmv(transA, alphap, mtxAp_, Xp, incX, betap, Yp, incY);
	}
	rsb_err_t spmm(rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, const rsb_flags_t order, const NT * Bp, const NT beta, NT * Cp)const
	{
		return spmm(transA, &alpha, nrhs, order, Bp, &beta, Cp);
	}
	rsb_err_t spmm(rsb_trans_t transA, const NT * alphap, rsb_coo_idx_t nrhs, const rsb_flags_t order, const NT * Bp, const NT * betap, NT * Cp)const
	{
		rsb_nnz_idx_t ldB = 0;
		rsb_nnz_idx_t ldC = 0;

		if(RSB_LIBRSB_VER>=10201)
			; // auto since 1.2.0.9
		else
		{
			// not yet auto
			rsb_coo_idx_t nr, nc;
			rsb_mtx_get_info(this->mtxAp_,RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T,&nr);
			rsb_mtx_get_info(this->mtxAp_,RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T,&nc);
			ldB = transA == RSB_TRANSPOSITION_N ? nr : nc;
			ldC = transA != RSB_TRANSPOSITION_N ? nr : nc;
		}
		return rsb_spmm(transA, alphap, mtxAp_, nrhs, order, Bp, ldB, betap, Cp, ldC);
	}
	rsb_err_t get_info(enum rsb_mif_t miflags, void* minfop)const
	{
		return rsb_mtx_get_info(mtxAp_, miflags, minfop);
	}
	std::string _info(void) const
	{
		char ss[RSBT_INFOBUF];
		rsb_mtx_get_info_str(this->mtxAp_,"RSB_MIF_MATRIX_INFO__TO__CHAR_P",ss,RSBT_INFOBUF);
		return std::string(ss);
	}

	~rsb_mtx_lw_t(void)
	{
		rsb_mtx_free(mtxAp_);
		mtxAp_ = nullptr;
	}
}; /* rsb_mtx_lw_t */
#ifdef RSBT_WANT_RSB_HPP
using RsbLib = rsb::RsbLib;
template<typename NT>
using RsbMatrix = rsb::RsbMatrix<NT>;
template<typename NT>
using rsb_mtx_w_t = RsbMatrix<NT>;
#else /* RSBT_WANT_RSB_HPP */
template<typename NT>
using rsb_mtx_w_t = rsb_mtx_lw_t<NT>;
#endif /* RSBT_WANT_RSB_HPP */

using times_t = int; // TODO: rename this type
using rsb_trans_c = rsb_trans_t;
//using rsb_trans_c = char; // just because rsb_trans_t is not printable

template<typename T>
static
typename std::enable_if< std::is_scalar<T>::value,T>::type rsbt_conj(T v) {return v;}

template<typename T>
static
typename std::enable_if<!std::is_scalar<T>::value,T>::type rsbt_conj(T v) {return std::conj(v);}

template<typename NT> static rsb_type_t rsb_type_t_for(void) { return RSB_NUMERICAL_TYPE_INVALID_TYPE; }
template<typename NT> static std::string rsb_type_s_for(void) { return "?"; }
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
template<> rsb_type_t rsb_type_t_for<long double> (void) { return RSB_NUMERICAL_TYPE_LONG_DOUBLE; }
template<> std::string rsb_type_s_for<long double> (void) { return "long double"; }
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
template<> rsb_type_t rsb_type_t_for<double> (void) { return RSB_NUMERICAL_TYPE_DOUBLE; }
template<> std::string rsb_type_s_for<double> (void) { return "double"; }
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
template<> rsb_type_t rsb_type_t_for<std::complex<double>> (void) { return RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX; }
template<> std::string rsb_type_s_for<std::complex<double>> (void) { return "std::complex<double>"; }
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
template<> rsb_type_t rsb_type_t_for<std::complex<long double>> (void) { return RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX; }
template<> std::string rsb_type_s_for<std::complex<long double>> (void) { return "std::complex<long double>"; }
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
template<> rsb_type_t rsb_type_t_for<std::complex<float>> (void) { return RSB_NUMERICAL_TYPE_FLOAT_COMPLEX; }
template<> std::string rsb_type_s_for<std::complex<float>> (void) { return "std::complex<float>"; }
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
template<> rsb_type_t rsb_type_t_for<float> (void) { return RSB_NUMERICAL_TYPE_FLOAT; }
template<> std::string rsb_type_s_for<float> (void) { return "float"; }
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
template<> rsb_type_t rsb_type_t_for<int> (void) { return RSB_NUMERICAL_TYPE_INT; }
#else
template<> rsb_type_t rsb_type_t_for<int> (void) { return 'I'; }
#endif
template<> std::string rsb_type_s_for<int> (void) { return "int"; }
#if RSB_WANT_LONG_IDX_TYPE
/* no rsb_type_t_for */
template<> std::string rsb_type_s_for<RSB_WANT_LONG_IDX_TYPE> (void) { return std::string("rsb_coo_idx_t"); }
#endif

#define RSB_INVALID_TRANS_ '?'

bool want_dump_failed = true;
const bool want_dump_vecs = false;
static int cnt {0};
static int scnt {0};
static int sntc {0};
static int sent {0};
static int sern {0};
static int ontc {0};
static int otcn {0};
static bool want_quiet {false};
static bool want_no_timings {false}; // for reproducibility
const rsb_time_t strt_t {rsb_time()};
rsb_time_t max_tt {std::numeric_limits<rsb_time_t>::max()};

#define RSBT_ENOUGH_TEST_CASES ( ( ontc > 0 && cnt >= ontc ) || (otcn != 0 && cnt > otcn) || ( rsb_time()-strt_t > max_tt ) )

static char char_for_flag(rsb_flags_t flags)
{
	return 
			 RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN)?'H':
			(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC)?'S':
			( RSB_DO_FLAG_HAS(flags,RSB_FLAG_TRIANGULAR)?'T':'G') );
}

static char trichar_for_flag(rsb_flags_t flags)
{
	return 
			 RSB_DO_FLAG_HAS(flags,RSB_FLAG_DIAGONAL)?'D':
			 (RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER)?'U':
			  (RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER)?'L':' '));
}

template<typename NT>
static
bool is_complex_type(void){return !std::is_scalar<NT>();}

template <typename NT=RSB_DEFAULT_TYPE, typename MWT=rsb_mtx_lw_t<NT>>
static rsb_err_t report(rsb_err_t errval, rsb_nnz_idx_t nr=0, rsb_nnz_idx_t nc=0, rsb_trans_t transA=RSB_INVALID_TRANS_, rsb_nnz_idx_t nnz=0, NT alpha=0, NT beta=0, rsb_blk_idx_t ns=0, const rsb_blk_idx_t nrhs=1, const rsb_blk_idx_t incX=1, const rsb_blk_idx_t incY=1, rsb_flags_t flags=RSB_FLAG_NOFLAGS, rsb_flags_t order=RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, NT mar=0, rsb_time_t dt=0.0)
{
	char c = rsb_type_t_for<NT>();
	const NT zero = 0;
	const NT one = 1;

	++cnt;

	if(want_quiet)
		return errval;

	std::cout << " #" << cnt << ": ";
	if(errval==RSB_ERR_NO_ERROR)
		std::cout << "[success]";
	else
		std::cout << "[failure]";
	std::cout << " sy:" << char_for_flag(flags) << trichar_for_flag(flags);
	std::cout << " ty:" << c << "";
	if(nr>0)
		std::cout << " nr:" << nr ;
	if(nc>0)
		std::cout << " nc:" << nc ;
	if(ns>0)
		std::cout << " ns:" << ns ;
	if(nnz>0)
		std::cout << " nnz:" << nnz ;
	if( transA != RSB_INVALID_TRANS_ )
		std::cout << " tr:" << static_cast<char>(transA);
	if(alpha!=zero)
		std::cout << " al:" << alpha ;
	if(nrhs>1)
		std::cout << " nrhs:" << nrhs ;
	if(incX>1)
		std::cout << " iX:" << incX;
	if(incY>1)
		std::cout << " iY:" << incY;
	if(beta!=one)
		std::cout << " be:" << beta;
	if(order!=RSB_FLAG_WANT_COLUMN_MAJOR_ORDER)
		std::cout << " order:rows";
	else
		std::cout << " order:cols";
	if(mar!=zero)
		std::cout << " err:" << mar;
	if(dt!=0.0)
	if(!want_no_timings)
	{
		const bool is_sym = flags&(RSB_SYMMETRY_S|RSB_SYMMETRY_H) ? true : false;
		const auto flops = (is_complex_type<NT>()?8:2)*(is_sym?2:1)*(nrhs)*(nnz/dt);
		std::cout << " dt:" << std::scientific << std::setprecision(2) << dt << std::defaultfloat;
		std::cout << " flops:" << std::scientific << std::setprecision(2) << flops << std::defaultfloat;
	}
#ifdef RSBT_WANT_RSB_HPP
	if(typeid(MWT) == typeid(RsbMatrix<NT>))
		std::cout << " [rsb.hpp:rsblib]";
#endif /* RSBT_WANT_RSB_HPP */
	if(typeid(MWT) == typeid(rsb_mtx_lw_t<NT>))
		std::cout << " [rsb.h:librsb]";
#ifdef HAVE_RSBPP_HPP
	if(typeid(MWT) == typeid(RsbPP_Matrix_T<NT>))
		std::cout << " [rsbpp.hpp:librsbpp]";
#endif /* HAVE_RSBPP_HPP */
	std::cout << std::endl;

	return errval;
} /* report */

template<typename NT, typename IT>
static 
typename std::enable_if< std::is_scalar<NT>::value,IT>::type
read_triplets(std::vector<NT>&VA, std::vector<IT>&IA, std::vector<IT>&JA, std::istream & ins, bool cmplx)
{
	// non-scalar output
	const static IT ione {1};
	IT rnz{};
	IT i,j;

	if(cmplx)
	{
		NT rp,ip;
		while( ins >> i && ins >> j && ins >> rp && ins >> ip )
			VA[rnz]={rp}, // we ignore ip
			IA[rnz]={i-ione}, JA[rnz]={j-ione},
			++rnz;
	}
	else
	{
		NT v;
		while( ins >> i && ins >> j && ins >> v )
			VA[rnz]={v},
			IA[rnz]={i-ione}, JA[rnz]={j-ione},
			++rnz;
	}
	return rnz;
} /* read_triplets */

template<typename NT, typename IT>
static 
typename std::enable_if<!std::is_scalar<NT>::value,IT>::type
read_triplets(std::vector<NT>&VA, std::vector<IT>&IA, std::vector<IT>&JA, std::istream & ins, bool cmplx)
{
	const static IT ione {1};
	IT rnz{};
	IT i,j;
	typename NT::value_type rp{},ip{};

	if(cmplx)
		while( ins >> i && ins >> j && ins >> rp && ins >> ip )
			VA[rnz]={rp,ip},
			IA[rnz]={i-ione}, JA[rnz]={j-ione},
			++rnz;
	else
		while( ins >> i && ins >> j && ins >> rp              )
			VA[rnz]={rp,ip}, // imaginary part set to zero
			IA[rnz]={i-ione}, JA[rnz]={j-ione},
			++rnz;
	return rnz;
} /* read_triplets */

template<class NT, class IT>
static 
int read_mm_as_coo(const std::string & fn, std::vector<NT>&VA, std::vector<IT>&IA, std::vector<IT>&JA)
{
	std::string tok;
	bool cmplx{false};
	IT nr{0},nc{0},nnz{0};
	std::fstream ins(fn,std::ios::in);
	
	if( (!(ins >> tok)) || tok!="%%MatrixMarket" )
		throw;
	if( (!(ins >> tok)) || tok!="matrix" )
		throw;
	if( (!(ins >> tok)) || tok!="coordinate" )
		throw;
	while(ins.peek()!='\n' && ins >> tok && !tok.empty() && std::isalpha(tok[0]) )
	//while(ins.peek()!='%' && std::isalpha(ins.peek()) ins >> tok && !tok.empty() && std::isalpha(tok[0]) )
	{
		if( tok=="complex" )
			cmplx=true;
		if( tok=="symmetric" )
			{}
		if( tok=="hermitian" )
			{}
		if( tok=="general" )
			{}
	}
	if(ins.peek()!='\n')
	       throw;
	ins.get();
	while(ins.peek()=='%')
		std::getline(ins,tok);
	ins >> nr >> nc >> nnz;
	std::vector<IT> my_IA(nnz,0),my_JA(nnz,0);
	std::vector<NT> my_VA(nnz,0);

	if(IT rnz = read_triplets(my_VA,my_IA,my_JA,ins,cmplx) != nnz)
	{
		std::cout << "read only " << rnz << " out of " << nnz << " values\n";
		throw;
	}
	//else std::cout << "read all " << nnz << " values\n";
	IA=std::move(my_IA);
	JA=std::move(my_JA);
	VA=std::move(my_VA);
	assert(*std::max_element(IA.begin(),IA.end())<nr);
	assert(*std::max_element(JA.begin(),JA.end())<nc);
	assert(*std::min_element(IA.begin(),IA.end())>=0);
	assert(*std::min_element(JA.begin(),JA.end())>=0);
	ins.close();
	return 0;
} /* read_mm_as_coo */

template<typename NT>
static 
void add(std::vector<NT> &y, const NT alpha)
{
	for(auto & v: y)
		v += alpha;
}

template<typename NT>
static 
void scale_by(std::vector<NT> &y, const NT beta)
{
	const NT one {1};
	if(beta!=one)
		for(auto & v: y)
			v *= beta;
}

using rsb_cv_t = std::vector<rsb_coo_idx_t>;
template<typename NT>
class coo_mtx_t
{
	public:// TODO: temporary
	using rsb_vt_t = std::vector<NT>;
	rsb_vt_t VA;
	rsb_cv_t IA;
	rsb_cv_t JA;
	public:
	typedef NT value_type;
	coo_mtx_t(const rsb_vt_t &VA, const rsb_cv_t &IA, const rsb_cv_t &JA):
		VA(VA),
		IA(IA),
		JA(JA)
	{
	}
	coo_mtx_t(rsb_vt_t&&VA, rsb_cv_t&&IA, rsb_cv_t&&JA):
		VA(VA),
		IA(IA),
		JA(JA)
	{
	}
	coo_mtx_t(std::string fn): VA{}, IA{}, JA{}
	{
		read_mm_as_coo<NT,rsb_coo_idx_t>(fn,VA,IA,JA);
	}
#ifdef RSBT_WANT_RSB_HPP
	coo_mtx_t(const RsbMatrix<NT> & mtx_A):
		VA(mtx_A.get_info_nnz_t(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T)),
		IA(mtx_A.get_info_nnz_t(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T)),
		JA(mtx_A.get_info_nnz_t(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T))
	{
		rsb_nnz_idx_t nnz{0};
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = mtx_A.get_rows_sparse(RSB_TRANSPOSITION_N,nullptr,VA.data(),IA.data(),JA.data(),0,mtx_A.rows()-1,&nnz,RSB_FLAG_NOFLAGS);
		assert(nnz==mtx_A.get_info_nnz_t(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T) && errval == RSB_ERR_NO_ERROR);
	}
#endif /* RSBT_WANT_RSB_HPP */
	rsb_nnz_idx_t get_nnz(void)const{return VA.size();}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_gn(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		const rsb_nnz_idx_t nnz=get_nnz();
		scale_by(y,beta);
		for( auto k = 0 ; k < nnz; ++k )
			y[IA[k]*incY]+=alpha*VA[k]*x[JA[k]*incX];
		return RSB_ERR_NO_ERROR;
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_gt(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		const rsb_nnz_idx_t nnz=get_nnz();
		scale_by(y,beta);
		for( auto k = 0 ; k < nnz; ++k )
			y[JA[k]*incY]+=alpha*VA[k]*x[IA[k]*incX];
		return RSB_ERR_NO_ERROR;
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_gc(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		const rsb_nnz_idx_t nnz=get_nnz();
		scale_by(y,beta);
		for( auto k = 0 ; k < nnz; ++k )
			y[JA[k]*incY]+=alpha*rsbt_conj(VA[k])*x[IA[k]*incX];
		return RSB_ERR_NO_ERROR;
	}

	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv(rsb_flags_t flags, rsb_trans_t transA, const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		rsb_err_t errval = RSB_ERR_GENERIC_ERROR;

		if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN))
		if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC))
		switch(transA)
		{
			case RSB_TRANSPOSITION_N:
				errval = spmv_sn(x, y, alpha, beta, incX, incY );
			break;
			case RSB_TRANSPOSITION_T:
				errval = spmv_st(x, y, alpha, beta, incX, incY );
			break;
			case RSB_TRANSPOSITION_C:
				errval = spmv_sc(x, y, alpha, beta, incX, incY );
			break;
		}

		if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC))
		if( RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN))
		switch(transA)
		{
			case RSB_TRANSPOSITION_N:
				errval = spmv_hn(x, y, alpha, beta, incX, incY );
			break;
			case RSB_TRANSPOSITION_T:
				errval = spmv_ht(x, y, alpha, beta, incX, incY );
			break;
			case RSB_TRANSPOSITION_C:
				errval = spmv_hc(x, y, alpha, beta, incX, incY );
			break;
		}

		if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN))
		if(!RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC))
		switch(transA)
		{
			case RSB_TRANSPOSITION_N:
				errval = spmv_gn(x, y, alpha, beta, incX, incY );
			break;
			case RSB_TRANSPOSITION_T:
				errval = spmv_gt(x, y, alpha, beta, incX, incY );
			break;
			case RSB_TRANSPOSITION_C:
				errval = spmv_gc(x, y, alpha, beta, incX, incY );
			break;
		}

		return errval;
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmm(rsb_flags_t flags, rsb_trans_t transA, rsb_coo_idx_t nrhs, const rsb_flags_t order, const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta )const
	{
		const rsb_nnz_idx_t ldX = x.size()/nrhs;
		const rsb_nnz_idx_t ldY = y.size()/nrhs;
		rsb_err_t errval = RSB_ERR_GENERIC_ERROR;

		rsb_vt_t tx(ldX);
		rsb_vt_t ty(ldY);

		if( order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER )
		for ( auto nrhsi = 0; nrhsi < nrhs; ++ nrhsi )
		{
			std::copy(x.begin()+ldX*nrhsi, x.begin()+ldX*(nrhsi+1), tx.begin());
			std::copy(y.begin()+ldY*nrhsi, y.begin()+ldY*(nrhsi+1), ty.begin());
			spmv(flags, transA, tx, ty, alpha, beta );
			std::copy(ty.begin(),ty.end(),y.begin()+ldY*nrhsi);
			errval = RSB_ERR_NO_ERROR;
		}
		else
		if( order == RSB_FLAG_WANT_ROW_MAJOR_ORDER )
		for ( auto nrhsi = 0; nrhsi < nrhs; ++ nrhsi )
		{
			for ( auto i = 0; i < ldX ; ++i)
				tx[i] = x[nrhsi+nrhs*i];
			for ( auto i = 0; i < ldY ; ++i)
				ty[i] = y[nrhsi+nrhs*i];
			spmv(flags, transA, tx, ty, alpha, beta );
			for ( auto i = 0; i < ldY ; ++i)
				y[nrhsi+nrhs*i] = ty[i];
			errval = RSB_ERR_NO_ERROR;
		}

		return errval;
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_sn(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		const rsb_nnz_idx_t nnz=get_nnz();
		scale_by(y,beta);
		for( auto k = 0 ; k < nnz; ++k )
		{
			y[JA[k]*incY]+=alpha*VA[k]*x[IA[k]*incX];
			if(IA[k]!=JA[k])
				y[IA[k]*incY]+=alpha*VA[k]*x[JA[k]*incX];
		}
		return RSB_ERR_NO_ERROR;
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_sc(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		const rsb_nnz_idx_t nnz=get_nnz();
		scale_by(y,beta);
		for( auto k = 0 ; k < nnz; ++k )
		{
			y[JA[k]*incY]+=alpha*rsbt_conj(VA[k])*x[IA[k]*incX];
			if(IA[k]!=JA[k])
				y[IA[k]*incY]+=alpha*rsbt_conj(VA[k])*x[JA[k]*incX];
		}
		return RSB_ERR_NO_ERROR;
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_st(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		return spmv_sn(x, y, alpha, beta, incX, incY );
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_hn(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		const rsb_nnz_idx_t nnz=get_nnz();
		scale_by(y,beta);
		for( auto k = 0 ; k < nnz; ++k )
		{
			y[IA[k]*incY]+=alpha*VA[k]*x[JA[k]*incX];
			if(IA[k]!=JA[k])
				y[JA[k]*incY]+=alpha*rsbt_conj(VA[k])*x[IA[k]*incX];
		}
		return RSB_ERR_NO_ERROR;
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_hc(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		return spmv_hn(x, y, alpha, beta, incX, incY );
	}
	template<typename IT=rsb_coo_idx_t>
	rsb_err_t spmv_ht(const rsb_vt_t &x, rsb_vt_t &y, NT alpha, NT beta, IT incX=1, IT incY=1 )const
	{
		const rsb_nnz_idx_t nnz=get_nnz();
		scale_by(y,beta);
		for( auto k = 0 ; k < nnz; ++k )
		{
			y[IA[k]*incY]+=alpha*rsbt_conj(VA[k])*x[JA[k]*incX];
			if(IA[k]!=JA[k])
				y[JA[k]*incY]+=alpha*VA[k]*x[IA[k]*incX];
		}
		return RSB_ERR_NO_ERROR;
	}
}; /* coo_mtx_t */

template<typename NT>
static
bool diovth(const NT &l, const NT &r)
{
	// differ over threshold
	const NT zero = 0;

	if(l == r)
		return false;
	else
	if ( std::abs(l) - std::abs(r) > RSBT_ETOL || std::abs(r) - std::abs(l) > RSBT_ETOL )
	{
		if(r != zero)
			return std::abs((l-r)/r) > RSBT_ETOL;
		else
			return std::abs((l-r)/l) > RSBT_ETOL;
	}
	return false;
}

template<typename NT, typename IT=rsb_coo_idx_t>
static 
bool any_diff(const std::vector<NT> &l, const std::vector<NT> &r, IT incY=1)
{
	const rsb_nnz_idx_t n = l.size() / incY;
	// assert(l.size()==r.size());
	for(auto i = 0; i<n; ++i)
		if(diovth(l[i*incY],r[i*incY]))
			return true;
	return false;
}

template<typename NT, typename IT=rsb_coo_idx_t>
static 
NT print_if_diff(const std::vector<NT> &l, const std::vector<NT> &r, IT incY=1, const rsb_nnz_idx_t dm=4)
{
	const rsb_nnz_idx_t n = l.size() / incY;
	rsb_nnz_idx_t dd=0;
	auto mar { std::abs(NT{}) };
	// assert(l.size()==r.size());
	for(auto oi = 0; oi<n; ++oi)
	{
		const auto i = oi * incY;

		if(l[i]!=r[i])
		{
			const NT zero = 0;
			if( l[i] != zero )
				mar = std::max(mar,std::abs((l[i]-r[i])/l[i]));

			if(diovth(l[i],r[i]))
			if(++dd<=dm)
			{
				std::cout << " diff @ " << i << " : " << l[i] << " vs " << r[i] << " (res vs ref)";
				std::cout << ", diff.: " << l[i]-r[i] << "";
				if( r[i] != zero )
				{
					std::cout << ", abs. rel. diff.: " << std::abs((l[i]-r[i])/l[i]) << " ";
				}
		       		std::cout << std::endl;
			}
		}
	}
	if(dd>dm)
		std::cout << " ... < a total of " << dd << " entries differ > ..."<< std::endl;

	return mar;
}

static 
void print_cxx_cmtln(const std::string & cmt)
{
	std::cout << "/* " << cmt << " */" << std::endl;
}

template<typename NT>
static 
void print_cxx_vec(std::string lbl, std::vector<NT> vec)
{
	std::cout << "std::vector<" << rsb_type_s_for<NT>() << "> " << lbl << " = {";
	for (auto x : vec)
		std::cout << " " << x << ",";
	std::cout << "};" << std::endl;
}

template<typename NT>
static 
void print_cxx_var(std::string lbl, NT x)
{
	std::cout << rsb_type_s_for<NT>() << " " << lbl << " = {";
	std::cout << " " << x << ",";
	std::cout << "};" << std::endl;
}

static 
rsb_coo_idx_t get_rci(rsb_coo_idx_t max)
{
	static std::random_device rd;
	static std::mt19937 g(rd());
	if(max>0)
		return g() % max;
	return 0;
}

static
bool shall_skip_case(void)
{
	/* uses globals */
	if((cnt+1 != otcn && otcn != 0) || cnt < sntc || ( sent > 0 && (cnt % sent) != 0 ) || RSBT_ENOUGH_TEST_CASES || ( sern > 0 && (cnt % (1+get_rci(sern)) ) != 0 ) )
		return true;
	else
		return false;
}

template <typename NT>
static 
void save_case(const rsb_trans_t transA, const rsb_flags_t order, const std::vector<NT> &cx, const std::vector<NT> &rx, const std::vector<NT> &ay, const std::vector<NT> &cy, const std::vector<NT> &ry, const NT &alpha, const NT &beta, const rsb_blk_idx_t nrhs, const rsb_coo_idx_t incX, const rsb_coo_idx_t incY)
{
		RSBT_SAVE_VEC("ry.mtx",NT,ry.data(),ry.size());
		RSBT_SAVE_VEC("cy.mtx",NT,cy.data(),cy.size());
		RSBT_SAVE_VEC("rx.mtx",NT,rx.data(),rx.size());
		RSBT_SAVE_VEC("cx.mtx",NT,cx.data(),cx.size());
		RSBT_SAVE_VEC("ay.mtx",NT,ay.data(),ay.size());
		if( want_dump_vecs )
		{
			print_cxx_cmtln("reference lhs before op:");
			print_cxx_vec("cy",cy);
			print_cxx_cmtln("rhs after reference op:");
			print_cxx_vec("cx",cx);
			print_cxx_cmtln("rhs after RSB op:");
			print_cxx_vec("rx",rx);
			print_cxx_cmtln("lhs after reference op:");
			print_cxx_vec("ry",ry);
			print_cxx_cmtln("lhs after RSB op:");
			print_cxx_vec("ay",ay);
			print_cxx_var("incX",incY);
			print_cxx_var("incY",incX);
			print_cxx_var("transA",transA);
			print_cxx_var("order",order);
			print_cxx_var("nrhs",nrhs);
		}
		if(std::equal(ay.begin(),ay.end(),cy.begin()))
			print_cxx_cmtln("problem: above ry does not match cy: SEEMS LIKE A NO OP!");
		else
			print_cxx_cmtln("problem: above ry does not match cy");
		print_cxx_var("alpha",alpha);
		print_cxx_var("beta",beta);
}

template <typename NT, typename MWT>
static 
rsb_err_t spmv_comparison(const rsb_trans_t transA, const MWT & mtxAw, const coo_mtx_t<NT> &cm, const rsb_flags_t order, const std::vector<NT> &cx, const std::vector<NT> &cy, const NT &alpha, const NT &beta, const rsb_blk_idx_t nrhs, const rsb_coo_idx_t incX, const rsb_coo_idx_t incY, int mr=0)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const std::vector<NT> rx(cx);
	std::vector<NT> ry(cy);
	std::vector<NT> ay(cy);
	const auto nr = mtxAw.rows();
	const auto nc = mtxAw.cols();
	//const auto nnzA = mtxAw.nnz();
	const auto flagsA = mtxAw.rsbflags();
	const auto ns = mtxAw.blocks();
	rsb_time_t dt;

	assert(nr>0);
	assert(nc>0);
	if(mr!=RSBT_FAIL_REC)
	{
		assert(!(nr!=nc && RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_LOWER)));
		assert(!(nr!=nc && RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_UPPER)));
	}
	assert(incX==1 || nrhs==1);
	assert(incY==1 || nrhs==1);

	if(nrhs == 1)
	{
		errval = cm.spmv(flagsA,transA,cx,ry,alpha,beta,incX,incY);
		if(errval != RSB_ERR_NO_ERROR)
			return_at_line(report<NT>(errval,nr,nc,transA));
		dt = -rsb_time();
#if defined(RSBP_WANT_CPP20) || defined(RSBPP_WANT_CPP20)
		if(incX==1 && incY==1)
			errval = mtxAw.spmv(transA, alpha, rx, beta, ay );
		else
#endif
			errval = mtxAw.spmv(transA, &alpha, rx.data(), incX, &beta, ay.data(), incY);
		dt += rsb_time();
		if(errval != RSB_ERR_NO_ERROR)
			return_at_line(report<NT>(errval,nr,nc,transA));
	}
	else
	{
		errval = cm.spmm(flagsA,transA,nrhs,order,cx,ry,alpha,beta);
		if(errval != RSB_ERR_NO_ERROR)
			return_at_line(report<NT>(errval,nr,nc,transA));
		dt = -rsb_time();
#if defined(RSBP_WANT_CPP20) || defined(RSBPP_WANT_CPP20)
		errval = mtxAw.spmm(transA, alpha, nrhs, order, rx, beta, ay);
#else
		errval = mtxAw.spmm(transA, alpha, nrhs, order, rx.data(), beta, ay.data());
#endif
		dt += rsb_time();
		if(errval != RSB_ERR_NO_ERROR)
			return_at_line(report<NT>(errval,nr,nc,transA));
	}
#ifdef HAVE_RSBPP_HPP
	if(errval == RSB_ERR_UNSUPPORTED_OPERATION)
	{
		std::cout << "WARNING: operation skipped" << RSB_XSTRINGIFY(RSB_ERR_UNSUPPORTED_OPERATION) << std::endl;
		errval=RSB_ERR_NO_ERROR;
		return_at_line(report<NT>(errval,nr,nc));
	}
	if(errval != RSB_ERR_NO_ERROR)
	{
		rsb_perror(stdout, errval);
		return_at_line(report<NT>(errval,nr,nc));
	}
#endif /* HAVE_RSBPP_HPP */

	const auto mar = print_if_diff(ay,ry,incY);
	if(any_diff(ay,ry,incY))
	{
		save_case(transA, order, cx, rx, ay, cy, ry, alpha, beta, nrhs, incX, incY);
		errval = RSBT_ERR_MISMATCH;
	}
	return_at_line((report<NT,decltype(mtxAw)>(errval,nr,nc,transA,cm.get_nnz(),alpha,beta,ns,nrhs,incX,incY,flagsA,order,mar,dt)));
} /* spmv_comparison */

template<typename NT>
static 
void fill_rhs(std::vector<NT> & v, int bv, int mr)
{
	NT iv = 0;
	auto n = v.size();
	decltype(n) i{};

	if(is_complex_type<NT>())
		iv = std::sqrt(NT{-1});

	if(mr)
		for(i = 0; i<n; i++)
			v[i]= bv + get_rci(mr),
			v[i]+=iv*static_cast<NT>(get_rci(mr));
	else
		for(i = 0; i<n; i++)
			v[i]=bv,
			v[i]+=iv;
}

struct rsb_lib_opts final
{
	rsb_int_t wvt {0};
	FILE*wvi{stdout};
	FILE*wve{stdout};
	FILE*wos{stdout};
	FILE*wvr{stdout};
	rsb_int_t evi{0};
	rsb_int_t llm{0};

	times_t maxtimes {1000*1000};
	times_t mintimes {3};
	rsb_time_t min_t {0.5};
	rsb_time_t max_t {std::numeric_limits<rsb_time_t>::max()};
	rsb_flags_t af {~RSB_FLAG_NOFLAGS}; // skip if any of those flags not intersecting
	rsb_flags_t sf { RSB_FLAG_NOFLAGS}; // skip if any of those flags     intersecting
	rsb_nnz_idx_t min_nnz {0};
	rsb_nnz_idx_t max_nnz {std::numeric_limits<rsb_nnz_idx_t>::max()};

template <typename T>
void set_opt(enum rsb_opt_t iof, const T*vp)const
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	
	if( ( errval = rsb_lib_set_opt(iof, vp)) != RSB_ERR_NO_ERROR )
		if( errval != RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT)
			throw;
}
template <typename T>
void get_opt(enum rsb_opt_t iof, T*vp)const
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	
	if( ( errval = rsb_lib_get_opt(iof, vp)) != RSB_ERR_NO_ERROR )
		if( errval != RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT)
			throw;
}
void set(void)const
{
	// TODO: need check to enforce execution after rsb_lib_init() and before rsb_lib_exit()

	set_opt(RSB_IO_WANT_VERBOSE_TUNING,&wvt);
	set_opt(RSB_IO_WANT_VERBOSE_INIT,&wvi);
	set_opt(RSB_IO_WANT_VERBOSE_EXIT,&wve);
	set_opt(RSB_IO_WANT_OUTPUT_STREAM,&wos);
	set_opt(RSB_IO_WANT_VERBOSE_ERRORS,&wvr);
	set_opt(RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE,&evi);
	set_opt(RSB_IO_WANT_LEAF_LEVEL_MULTIVEC ,&llm); // TODO: shall quantify benefit
#if 0
	// the following wait for a systematic approach:
	RSB_IO_WANT_VERBOSE_INIT =0x000001	/* (FILE*) */
	RSB_IO_WANT_VERBOSE_EXIT =0x000002	/* (FILE*) */
	RSB_IO_WANT_OUTPUT_STREAM =0x000003	/* (FILE*) */
	RSB_IO_WANT_SORT_METHOD =0x000004	/* (rsb_int_t) */
	RSB_IO_WANT_CACHE_BLOCKING_METHOD =0x000005	/* (rsb_int_t) */
	RSB_IO_WANT_SUBDIVISION_MULTIPLIER =0x000006	/* (rsb_real_t) */
	RSB_IO_WANT_VERBOSE_ERRORS =0x000007	/* (FILE*) */
	RSB_IO_WANT_BOUNDED_BOX_COMPUTATION =0x000008	/* (rsb_int_t) */
	RSB_IO_WANT_EXECUTING_THREADS =0x000009	/* (rsb_int_t) */
	RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE =0x000010	/* (rsb_int_t) */
	RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING =0x000011	/* (const rsb_char_t*) */
	RSB_IO_WANT_IS_INITIALIZED_MARKER =0x000012	/* (const rsb_bool_t*) */
	RSB_IO_WANT_MEM_ALLOC_CNT =0x000013	/* (const size_t*) */
	RSB_IO_WANT_MEM_ALLOC_TOT =0x000014	/* (const size_t*) */
	RSB_IO_WANT_LEAF_LEVEL_MULTIVEC =0x000015	/* (rsb_int_t) */
	RSB_IO_WANT_MAX_MEMORY_ALLOCATIONS =0x000016	/* (size_t) */
	RSB_IO_WANT_MAX_MEMORY_ALLOCATED =0x000017	/* (size_t) */
	RSB_IO_WANT_LIBRSB_ETIME =0x000018	/* (rsb_time_t) */
	RSB_IO_WANT_VERBOSE_TUNING =0x000019	/* (rsb_int_t) */
#endif
} /* set */
}; /* rsb_lib_opts */

struct rsbt_bopts final
{
	std::vector<std::string> fna {};
	std::vector<rsb_int_t> tna={0};
	std::vector<rsb_coo_idx_t> nrhsa={1,2};
	std::vector<rsb_coo_idx_t> incxa={1,2};
	std::vector<rsb_coo_idx_t> incya={1,2};
	std::vector<int> alphaia={1,-1,4,-4};
	const std::vector<int> kappaia={1,-1};
	std::vector<int> dnst_pct_a={50,100};
	std::vector<int> betaia={0,1,4};
	std::vector<rsb_coo_idx_t> csf_fa={1}; // coordinates stretch factor
	std::vector<rsb_coo_idx_t> bsf_fa={0}; // band stretch factor (<0: randomized additive, >0: additive)
	std::vector<rsb_flags_t> symfa = {RSB_FLAG_NOFLAGS,RSB_FLAG_UPPER_SYMMETRIC,RSB_FLAG_UPPER_HERMITIAN,RSB_FLAG_LOWER_SYMMETRIC,RSB_FLAG_LOWER_HERMITIAN}; // symmetry flags array
	std::string rfn{""}; // report file name
	std::string ts{""}; // types string
	std::vector<rsb_trans_c> transAa = {RSB_TRANSPOSITION_N,RSB_TRANSPOSITION_T,RSB_TRANSPOSITION_C};
	std::vector<rsb_flags_t> ordera = { RSB_FLAG_WANT_ROW_MAJOR_ORDER, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER };
	rsb_int_t maxr{RSB_CONST_MAX_TUNING_ROUNDS};
	rsb_time_t maxt{0.0};
	enum RenderOpts { No, Once, All };
#ifdef RSBT_WANT_RSB_HPP
	RenderOpts ropts { No}; // used with --report --render-only --render hpp_test hpp_tests 
#endif
	const bool wrl{true}; // want rsblib
	bool want_basenameplot{true};
	bool want_autotune{true};
	bool want_mkl{false};
	bool want_tmm{false}; // tolerate mismatches
	bool want_any_herm{false}; // any type as hermitian is ok
	rsb_int_t mor{10}; // margins of randomization; on if >0
	rsb_flags_t eflags {RSB_FLAG_NOFLAGS}; /* extra flags */
	bool want_rpp {false}; /* want rsbpp.hpp */
	rsb_lib_opts dro;
}; /* rsbt_bopts */

template <typename NT, typename IT>
std::vector<NT> i2tv(const std::vector<IT>& iv, bool add_imag_scaled=false)
{
	std::vector<NT> tv(iv.size());
	std::copy(iv.begin(), iv.end(), tv.begin());

	if(add_imag_scaled)
	{
		const NT mo = -1;
		const NT iv = std::sqrt(mo);
		const std::vector<NT> alpha_c { tv };
		for (const auto ac : alpha_c)
			tv.push_back(iv*ac);
	}

	return tv;
}

std::vector<rsb_coo_idx_t> if_no_trans_then_nc_else_nr(const std::vector<rsb_trans_t> & trans_a, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA)
{
	std::vector<rsb_coo_idx_t> ld_a(trans_a.size());

	std::transform(trans_a.begin(), trans_a.end(), ld_a.begin(),
		[&](rsb_trans_t t) -> rsb_coo_idx_t { return t == RSB_TRANSPOSITION_N ? ncA : nrA; });
	return ld_a;
}

static 
int half_square(const int n)
{
	return (n * (n+1))/2;
}

struct rsbt_btrs final // benchmarking timing results structure
{
	rsb_time_t avt{0.0},lot{0.0},hit{0.0};
	times_t times{0};
};

template<typename T>
static
T pcnt_of(T what, T of)
{
	return 100.0*((what)/of);
}

template<typename T>
static
T pcnt_diff(T what, T of)
{
	return (pcnt_of(what-of,of));
}

#ifdef RSBT_WANT_RSB_HPP
static
std::ostream& operator<<(std::ostream &os, const rsbt_btrs & rs)
{
	std::ostringstream rsos;
	rsos.precision(2);
	if(rs.times==1)
		rsos	<<  rs.avt;
	else
		rsos	<<  rs.times << " samples for " << rs.avt << " +" << pcnt_diff(rs.hit,rs.avt) << "%/" << pcnt_diff(rs.lot,rs.avt) << "%";
		/*rsos 
			<<  "min: " << rs.lot 
			<< " max: " << rs.hit
			<< " avg: " << rs.avt;*/
	os << rsos.str();
	return os;
}
#endif /* RSBT_WANT_RSB_HPP */

struct rsbt_bress final
{
	std::ostringstream rs; // report stream
	struct lrpi_t // latex report page info:
	{ std::string title{"Sample Title"},pfn{""},caption{"Sample caption"},text{""}; }; // title, plot name
	std::vector<lrpi_t> lrpv;
}; /* rsbt_bress */

struct rsbt_topts final
{
	rsb_coo_idx_t n{RSB_INVALID_NNZ_IDX_VAL};
	rsb_flags_t fflags {RSB_FLAG_NOFLAGS}; /* format flags */
}; /* rsbt_topts */

struct rsb_lib final
{
private:
#ifdef HAVE_MUTEX_H
	static std::mutex mutex;
#endif /* HAVE_MUTEX_H */
public:
	rsb_lib(bool vinit=true, bool vexit=true, bool verrs=true, bool verbt=true, bool vintf=false)
	{
		{
#ifdef HAVE_MUTEX_H
		std::lock_guard< std::mutex > lock( rsb_lib::mutex );
#endif /* HAVE_MUTEX_H */
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS);
		if(errval != RSB_ERR_NO_ERROR)
		{
			report_at_line(errval);
			throw;
		}
		}
		auto sp { stdout }; /* can't assume is an lvalue and therefore usable in '&stdout' */
		rsb_int_t vil = 1;

		if(vinit)
			set_opt(RSB_IO_WANT_VERBOSE_INIT, &sp);
		if(vexit)
			set_opt(RSB_IO_WANT_VERBOSE_EXIT, &sp);
		if(verrs)
			set_opt(RSB_IO_WANT_VERBOSE_ERRORS, &sp);
		if(verbt)
			set_opt(RSB_IO_WANT_VERBOSE_TUNING, &vil);
		if(vintf)
			set_opt(RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE, &vil);
	}
template <typename T>
	void set_opt(enum rsb_opt_t iof, const T*vp)const
	{
#ifdef HAVE_MUTEX_H
		std::lock_guard< std::mutex > lock( rsb_lib::mutex );
#endif /* HAVE_MUTEX_H */
		if(rsb_lib_set_opt(iof, vp)!=RSB_ERR_NO_ERROR)
			throw;
	}
	~rsb_lib(void)
	{
#ifdef HAVE_MUTEX_H
		std::lock_guard< std::mutex > lock( rsb_lib::mutex );
#endif /* HAVE_MUTEX_H */
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
		if(errval != RSB_ERR_NO_ERROR)
		{
			report_at_line(errval);
		}
	}
};

#ifdef HAVE_MUTEX_H
	std::mutex rsb_lib::mutex;
#endif /* HAVE_MUTEX_H */

template<typename NT>
static 
void set_diag_coo(std::vector<NT>& VA,const std::vector<rsb_coo_idx_t>& IA,const std::vector<rsb_coo_idx_t>& JA, const NT dv)
{
	const auto nnzA = VA.size();

	for( auto k {nnzA-nnzA}; k<nnzA; ++k )
		if (IA[k] == JA[k])
			VA[k] = dv;
}

template<typename NT>
static 
rsb_err_t gen_rectangular(std::vector<NT>& VA,std::vector<rsb_coo_idx_t>& IA,std::vector<rsb_coo_idx_t>& JA, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const rsb_flags_t sflags, const int mr, int dnst_pct)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;
	const NT one = 1;
	const rsb_nnz_idx_t nnzA=nr*nc;
	NT val = one;

	assert(dnst_pct <= 100);
	assert(dnst_pct >    0);

	if(is_complex_type<NT>())
	{
		const NT mo = -1;
		const NT i = std::sqrt(mo);
		val += i;
	}

	for(auto i = 0; i<nr; ++i)
		for(auto j = 0; j<nc; ++j)
		{
			if ( dnst_pct < 100 && get_rci(100) < ( 100 - dnst_pct ) )
				VA[k]=0;
			else
			{
				if (mr > 0)
					VA[k]=val*NT(get_rci(mr));
				VA[k]+=get_rci(mr);
			}
			JA[k]=j;
		       	IA[k++]=i;
		}
	if(k!=nnzA)
		return_at_line(errval=RSB_ERR_GENERIC_ERROR);

	if( RSB_DO_FLAG_HAS(sflags,RSB_FLAG_UNIT_DIAG_IMPLICIT) )
		set_diag_coo(VA,IA,JA,one);

	if( RSB_DO_FLAG_HAS(sflags,RSB_FLAG_FORTRAN_INDICES_INTERFACE) )
		add(IA,rsb_coo_idx_t{1}), add(JA,rsb_coo_idx_t{1});

	return errval;
}

template<typename NT>
static 
rsb_err_t gen_banded(std::vector<NT>& VA,std::vector<rsb_coo_idx_t>& IA,std::vector<rsb_coo_idx_t>& JA, const rsb_coo_idx_t b, const rsb_coo_idx_t n, const rsb_flags_t sflags, const int mr, int dnst_pct = 100)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t k = 0;
	const NT one = 1;
	const rsb_nnz_idx_t nnzA=half_square(n)-half_square(n-b-1);
	NT val = one;

	assert(dnst_pct <= 100);
	assert(dnst_pct >    0);

	if(is_complex_type<NT>())
	{
		const NT mo = -1;
		const NT i = std::sqrt(mo);
		val += i;
	}

	for(auto i = k; i<n; ++i)
		for(auto j = i; j<std::min(i+b+1,n); ++j)
		{
			if(i==j && is_complex_type<NT>() && RSB_DO_FLAG_HAS(sflags,RSB_FLAG_HERMITIAN ))
				VA[k]=one;
			else
			{
				if ( dnst_pct < 100 && get_rci(100) < ( 100 - dnst_pct ) )
					VA[k]=0;
				else
					VA[k]=val*NT(get_rci(mr)),
					VA[k]+=get_rci(mr);
			}

			if(RSB_DO_FLAG_HAS(sflags,RSB_FLAG_UPPER))
				IA[k]=i, JA[k++]=j;
			else
				JA[k]=i, IA[k++]=j;
		}
	if(k!=nnzA)
		return_at_line(errval=RSB_ERR_GENERIC_ERROR);

	if( RSB_DO_FLAG_HAS(sflags,RSB_FLAG_UNIT_DIAG_IMPLICIT) )
		set_diag_coo(VA,IA,JA,one);

	if( RSB_DO_FLAG_HAS(sflags,RSB_FLAG_FORTRAN_INDICES_INTERFACE) )
		add(IA,rsb_coo_idx_t{1}), add(JA,rsb_coo_idx_t{1});
	//std::cout << " " << n << " " << nnzA << " " << k << std::endl;
	return errval;
}

#define timed_exec_macro(callable, mintimes, maxtimes, min_t, max_t, dtv)	\
{	\
	times_t times=0;	\
	if( maxtimes > 1 ) callable; /* warmup call, unless maxtimes <= 1 */ \
	\
	rsb_time_t bdt = rsb_time();	\
	rsb_time_t min_bdt = std::numeric_limits<rsb_time_t >::max(), max_bdt = 0;	\
	\
	do	/* at least once */	\
	{			\
		rsb_time_t dt=-rsb_time();	\
		callable;	\
		dt+=rsb_time();	\
		++times;	\
		min_bdt=std::min(min_bdt,dt); \
		max_bdt=std::max(max_bdt,dt); \
	}			\
	while( ( rsb_time() - bdt < min_t || times < mintimes ) && ( times < maxtimes && rsb_time() - bdt < max_t )  ); \
	bdt = rsb_time() - bdt;	\
	bdt/=times;	\
	dtv.avt=bdt;	\
	dtv.lot=min_bdt;	\
	dtv.hit=max_bdt;	\
	dtv.times=times;	\
}
template<class C>
static
void timed_exec_func(std::function<C> & callable, const rsb_lib_opts ro, rsbt_btrs &dtv)
{
	/* it seems like serial C++ implementation is slower here in the lambda --- straight macro invocation is faster. */
	timed_exec_macro(callable(), ro.mintimes, ro.maxtimes, ro.min_t, ro.max_t, dtv);
}
#define timed_exec_m(expr, opts, dtv)                                      timed_exec_macro(expr,    opts.mintimes, opts.maxtimes, opts.min_t, opts.max_t, dtv)
#define timed_exec_f(expr, opts, dtv) {std::function callable = [&]()expr; timed_exec_func(callable, opts                                                , dtv);}

#ifdef RSBT_USE_MKL_SPBLAS
void rsb_mkl_err_handle(sparse_status_t status)
{
	if(SPARSE_STATUS_SUCCESS!=status)
	{
		std::cout << "MKL gave error code: " << status << std::endl;
		throw;
	}
}
#endif /* RSBT_USE_MKL_SPBLAS */

#ifdef RSBT_USE_MKL_SPBLAS
template<class IT, class NT>
static
int test_mkl(rsbt_btrs &mbtr, const NT alpha, const NT beta, std::vector<NT>&VA, std::vector<IT>&IA, std::vector<IT>&JA, rsb_flags_t flagsA, IT rows, IT cols, /*IT trows, IT tcols,*/ rsb_trans_t transA, const rsb_lib_opts ro, const std::vector<NT> & x, std::vector<NT> & y)
{
	struct matrix_descr descrA;
	sparse_matrix_t mklAcoo, mklA;
	std::map<rsb_trans_t,sparse_operation_t> rsbttomklt{
		{RSB_TRANSPOSITION_N,SPARSE_OPERATION_NON_TRANSPOSE},
		{RSB_TRANSPOSITION_T,SPARSE_OPERATION_TRANSPOSE},
		{RSB_TRANSPOSITION_C,SPARSE_OPERATION_CONJUGATE_TRANSPOSE},
	};
	std::map<rsb_flags_t,sparse_matrix_type_t> rsbstomkls{
		{RSB_SYMMETRY_U,SPARSE_MATRIX_TYPE_GENERAL},
		{RSB_SYMMETRY_S,SPARSE_MATRIX_TYPE_SYMMETRIC},
		{RSB_SYMMETRY_H,SPARSE_MATRIX_TYPE_HERMITIAN},
	};
	sparse_status_t status=SPARSE_STATUS_SUCCESS;
	//const int is_sym = flagsA&(RSB_SYMMETRY_S|RSB_SYMMETRY_H) ? 1 : 0;
	//const auto op_nnz=VA.size()*2*is_sym ;

	if(sizeof(MKL_INT)!=sizeof(IT))
		throw; // IA,JA shall be MKL_INT
	status = rsb_mkl_sparse_X_create_coo ( &mklAcoo, SPARSE_INDEX_BASE_ZERO, rows,  cols,  VA.size(), IA.data(), JA.data(), VA.data());
	//mkl_sparse_d_create_coo(&mklA,SPARSE_INDEX_BASE_ZERO,rows,cols,VA.size(),IA.data(),  JA.data(),VA.data());
	//mkl_sparse_d_create_csr(&mklA,SPARSE_INDEX_BASE_ZERO,rows,cols,PA.data(),PA.data()+1,JA.data(),VA.data());
	rsb_mkl_err_handle(status);
	status = mkl_sparse_convert_csr ( mklAcoo, SPARSE_OPERATION_NON_TRANSPOSE, &mklA);
	rsb_mkl_err_handle(status);
	status = mkl_sparse_destroy ( mklAcoo );
	rsb_mkl_err_handle(status);
	descrA.type = rsbstomkls[flagsA&(RSB_SYMMETRY_U|RSB_SYMMETRY_S|RSB_SYMMETRY_H)];

	if(RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_UPPER))
		descrA.mode = SPARSE_FILL_MODE_UPPER;
	if(RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_LOWER))
		descrA.mode = SPARSE_FILL_MODE_LOWER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	status = mkl_sparse_optimize ( mklA );
	rsb_mkl_err_handle(status);
	timed_exec_m({ status = rsb_mkl_sparse_X_mv(rsbttomklt[transA],alpha,mklA,descrA,x.data(),beta,y.data()); },ro,mbtr);
	rsb_mkl_err_handle(status);
	status = mkl_sparse_destroy ( mklA );
	rsb_mkl_err_handle(status);
	return 0;
}
#endif /* RSBT_USE_MKL_SPBLAS */

template <typename NT, typename MWT, typename EAT>
static
rsb_err_t spmv_comparison_rectangular(const rsbt_bopts & bopts, rsb_nnz_idx_t nrA, rsb_nnz_idx_t ncA, MWT & mtxAw, const coo_mtx_t<NT> &cm, int mr, EAT soe)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const std::vector<NT> beta_a { i2tv<NT>(bopts.betaia) };
	const std::vector<NT> kappa_a { i2tv<NT>(bopts.kappaia) };
	const std::vector<rsb_trans_t> trans_a{bopts.transAa};
	const std::vector<rsb_coo_idx_t> ldX_a=if_no_trans_then_nc_else_nr(trans_a,nrA,ncA);
	const std::vector<rsb_coo_idx_t> ldY_a=if_no_trans_then_nc_else_nr(trans_a,ncA,nrA);
	const int tn = trans_a.size();
	const bool add_imag_scaled = is_complex_type<NT>() && kappa_a.size() > 1;
	const std::vector<NT> alpha_a { i2tv<NT>(bopts.alphaia,add_imag_scaled) };
	const auto nrhs_a=bopts.nrhsa;
	rsb_blk_idx_t ns=0;

	errval = mtxAw.get_info(RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T, &ns);
	if(errval != RSB_ERR_NO_ERROR)
		return_at_line(report<NT>(errval,nrA,ncA));

	for (const auto nrhs : nrhs_a)
	for (int bv = 0; bv < 1 + (nrA == ncA ? 1 : 0); ++bv)
	for (const NT alpha : alpha_a)
	for (const NT beta  : beta_a )
	for (int ti =0; ti< tn; ++ti)
	for (const auto incX: bopts.incxa)
	for (const auto incY: bopts.incya)
	for(const rsb_flags_t order : bopts.ordera)
	if( incX==1 || nrhs==1 )
	if( incY==1 || nrhs==1 )
	if( ! ( nrhs == 1 && bopts.ordera.size() == 2 && order == bopts.ordera[1] ) )
	if( ! RSBT_ENOUGH_TEST_CASES )
	{
	if( shall_skip_case() )
	{
		++cnt;
		++scnt;
		continue;
	}
	else
	{
		auto transA = trans_a[ti];
		auto ldX = ldX_a[ti];
		auto ldY = ldY_a[ti];
		std::vector<NT> cx(ldX*nrhs*incX,1);
		std::vector<NT> cy(ldY*nrhs*incY,1);

		fill_rhs(cx,bv,mr);
		fill_rhs(cy,bv,mr);
		if(mr==RSBT_FAIL_REC)
		{
			std::cout << "Recover y and x from cy.mtx and cx.mtx" << std::endl;
			RSBT_LOAD_VEC("cx.mtx", NT, cx);
			RSBT_LOAD_VEC("cy.mtx", NT, cy);
			std::cout << "Recovered y sized " << cy.size() << " and x sized " << cx.size() << std::endl;
		}
		errval = spmv_comparison(transA, mtxAw, cm, order, cx, cy, alpha, beta, nrhs, incX, incY, mr);
		if(errval != RSB_ERR_NO_ERROR)
		{
			soe();
			if ( bopts.want_tmm && ( errval == RSBT_ERR_MISMATCH ) )
				std::cout << "A numeric mismatch occurred." << std::endl;
			else
				return_at_line(errval);
		}

#ifdef RSBT_USE_MKL_SPBLAS
		if(bopts.want_mkl)
		if(nrhs==1)
		if(incX==1)
		if(incY==1)
		{
			coo_mtx_t<NT> ccm{cm};
			rsbt_btrs mbtr;
			rsb_lib_opts ro{};
			const auto flagsA { mtxAw.rsbflags() };
			ro.min_t = 0.0;
			auto ay {cy}, ry {cy};
			test_mkl(mbtr,alpha,beta,ccm.VA,ccm.IA,ccm.JA,flagsA,nrA,ncA,/*ldY,ldX,*/transA,ro,cx,ay);

			rsb_err_t errval = cm.spmv(flagsA,transA,cx,ry,alpha,beta,incX,incY);
			if(errval != RSB_ERR_NO_ERROR)
				return_at_line(report<NT>(errval,nrA,ncA,transA));
			const auto mar = print_if_diff(ay,ry,incY);
			if(any_diff(ay,ry,incY))
			{
				save_case(transA, order, cx, cx, ay, cy, ry, alpha, beta, nrhs, incX, incY);
				errval = RSBT_ERR_MISMATCH;
			}
			std::cout << "[MKL:]" << std::endl;
			return_at_line((report<NT,decltype(mtxAw)>(errval,nrA,ncA,transA,cm.get_nnz(),alpha,beta,ns,nrhs,incX,incY,flagsA,order,mar)));
			//oss << "MKL mkl_sparse_*_mv [csr] took " << mbtr << " s, " << 1.0/(mbtr.avt/op_nnz) << " nnz/s (macro)" << std::endl;
		}
#endif /* RSBT_USE_MKL_SPBLAS */
	}
	}
	return errval;
} /* spmv_comparison_rectangular */

static bool unfit_for_strict_csr(rsb_flags_t flagsA, rsb_coo_idx_t nrA, rsb_coo_idx_t nnzA )
{
	// for --csr-tester
	return ( (!RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_QUAD_PARTITIONING)) && RSB_DO_FLAG_HAS_INTERSECTION(flagsA,(RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS|RSB_FLAG_USE_CSR_RESERVED)) && nnzA < 2*(nrA+1) );
}

template<typename NT>
static 
rsb_err_t rectangular_tester(const rsbt_bopts & bopts, const rsbt_topts topts, const rsb_coo_idx_t b, rsb_flags_t sflags, int mr)
{
	const auto n{topts.n};
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_lib rli{};

	bopts.dro.set();

	rsb_nnz_idx_t nnzA=n*b; // actually, max nnzA
	const rsb_type_t typecode=rsb_type_t_for<NT>();
	const NT one = 1;
	std::vector<NT> VA(nnzA,one);
	std::vector<rsb_coo_idx_t> IA(nnzA,0);
	std::vector<rsb_coo_idx_t> JA(nnzA,0);
	const rsb_flags_t flagsA=topts.fflags | sflags | RSB_FLAG_DISCARD_ZEROS;
	struct rsb_mtx_t * mtxAp = nullptr;

	//mtxAp = rsb_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA, flagsA, &errval);
	assert(b>0);

	assert(!RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_LOWER));
	assert(!RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_UPPER));

	for(auto dnst_pct : bopts.dnst_pct_a)
	for(const auto csf_f : bopts.csf_fa)
	for(const bool want_wide : {true,false})
	if (! ( dnst_pct < 100 && n == 1 && b == 1 )) // avoid [0] matrix
	if (! ( want_wide && n == 1 && b == 1 )) // avoid repetition
	{
		const rsb_coo_idx_t nrA = (want_wide?b:n)*csf_f;
		const rsb_coo_idx_t ncA = (want_wide?n:b)*csf_f;

		errval = gen_rectangular(VA,IA,JA,nrA/csf_f,ncA/csf_f,flagsA,mr,dnst_pct);
		if(errval != RSB_ERR_NO_ERROR)
			return_at_line(errval);

		scale_by(IA,csf_f);
		scale_by(JA,csf_f);

		//errval = rsb_mtx_set_vals(mtxAp, VA.data(), IA.data(), JA.data(), nnzA, RSB_FLAG_NOFLAGS);
		//if(errval != RSB_ERR_NO_ERROR)
		//	return_at_line(report<NT>(errval,nrA,ncA));

		if( unfit_for_strict_csr(flagsA, nrA, nnzA) )
			continue;

		if(mr==RSBT_FAIL_REC)
		{
			std::cout << "Recover matrix from fail.mtx" << std::endl;
			mtxAp = rsb_file_mtx_load("fail.mtx", flagsA, typecode, &errval);
			if(errval == RSB_ERR_NO_ERROR)
			{
				errval = rsb_mtx_get_info(mtxAp, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &nnzA);
				IA.resize(nnzA);
				JA.resize(nnzA);
				VA.resize(nnzA);
				errval |= rsb_mtx_get_coo(mtxAp, VA.data(), IA.data(), JA.data(), RSB_FLAG_NOFLAGS);
			}
		}
		else
			mtxAp = rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, 1, 1, flagsA, &errval);

		if(errval != RSB_ERR_NO_ERROR)
			return_at_line(report<NT>(errval,nrA,ncA));

		if( RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_FORTRAN_INDICES_INTERFACE) )
			add(IA,rsb_coo_idx_t{-1}), add(JA,rsb_coo_idx_t{-1});

		if( bopts.wrl || bopts.want_rpp )
		{
			const coo_mtx_t<NT> cm(VA,IA,JA);
			const rsb_mtx_w_t<NT> mtxAw(mtxAp);

			auto soe = [&](void)
			{
				std::cout << mtxAw._info() << std::endl << std::endl; // FIXME: works only if RsbMatrix
				if( want_dump_vecs )
				{
					print_cxx_vec("IA",IA);
					print_cxx_vec("JA",JA);
					print_cxx_vec("VA",VA);
					for(auto k = 0; k<nnzA; ++k)
						std::cout << IA[k] << " " << JA[k] << " " << VA[k] << std::endl;
				}
				RSBT_SAVE_MAT("fail.mtx",mtxAp);
				return_at_line(errval);
			};

			if( bopts.wrl )
				errval = spmv_comparison_rectangular<NT>(bopts,nrA,ncA,mtxAw,cm,mr,soe);

#ifdef HAVE_RSBPP_HPP
			if( bopts.want_rpp )
			{
				if ( false )
				{
					std::cout << "WARNING: Skipping RsbPP_Matrix_T<NT> because of symmetry." << std::endl;
				}
				else
				{
					/* RsbPP_Matrix_T<NT>::spmv is still unfinished */
					const RsbPP_Matrix_T<NT> mtxAppp(VA,IA,JA,flagsA,RSBPP_MATRIX_T_ORDERING);
					errval = spmv_comparison_rectangular<NT>(bopts,nrA,ncA,mtxAppp,cm,mr,soe);
					if(errval != RSB_ERR_NO_ERROR)
					{
						mtxAppp.print(std::cout);
					}
				}
			}
#endif /* HAVE_RSBPP_HPP */
		}
	}
	return errval;
} /* rectangular_tester */

template<typename NT>
static 
rsb_err_t band_tester(const rsbt_bopts & bopts, const rsbt_topts topts, const rsb_coo_idx_t b, rsb_flags_t sflags, int mr, int dnst_pct, rsb_coo_idx_t bsf_f, rsb_coo_idx_t csf_f)
{
	const auto n{topts.n};
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_lib rli{};

	bopts.dro.set();

	if(b>n-1)
		return errval;

	rsb_nnz_idx_t nnzA=half_square(n)-half_square(n-b-1);
	const rsb_type_t typecode=rsb_type_t_for<NT>();
	const rsb_coo_idx_t nrA=n*csf_f;
	const rsb_coo_idx_t ncA=n*csf_f;
	const NT one = 1;
	std::vector<NT> VA(nnzA,one);
	std::vector<rsb_coo_idx_t> IA(nnzA,0);
	std::vector<rsb_coo_idx_t> JA(nnzA,0);
	const rsb_flags_t flagsA=topts.fflags | sflags | RSB_FLAG_DISCARD_ZEROS | ( bsf_f ? RSB_FLAG_DUPLICATES_SUM : RSB_FLAG_NOFLAGS );
	struct rsb_mtx_t * mtxAp = nullptr;
	//mtxAp = rsb_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA, flagsA, &errval);

	errval = gen_banded(VA,IA,JA,b,nrA/csf_f,flagsA,mr,dnst_pct);
	if(errval != RSB_ERR_NO_ERROR)
		return_at_line(errval);

	scale_by(IA,csf_f);
	scale_by(JA,csf_f);

	if(bsf_f > 0)
		for(auto k = 0; k<nnzA; ++k)
			if ( IA[k] != JA[k] )
				JA[k] = ( ncA + JA[k] + (JA[k] - IA[k]) + bsf_f ) % ncA;
	if(bsf_f < 0)
		for(auto k = 0; k<nnzA; ++k)
			if ( IA[k] != JA[k] )
				JA[k] = ( ncA + JA[k] + (JA[k] - IA[k]) + get_rci(ncA) ) % ncA;

	//errval = rsb_mtx_set_vals(mtxAp, VA.data(), IA.data(), JA.data(), nnzA, RSB_FLAG_NOFLAGS);
	//if(errval != RSB_ERR_NO_ERROR)
	//	return_at_line(report<NT>(errval,nrA));

	if( unfit_for_strict_csr(flagsA, nrA, nnzA) )
		return RSB_ERR_NO_ERROR;

	if(mr==RSBT_FAIL_REC)
	{
		std::cout << "Recover matrix from fail.mtx" << std::endl;
		mtxAp = rsb_file_mtx_load("fail.mtx", flagsA, typecode, &errval);
		if(errval == RSB_ERR_NO_ERROR)
		{
			errval = rsb_mtx_get_info(mtxAp, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &nnzA);
			IA.resize(nnzA);
			JA.resize(nnzA);
			VA.resize(nnzA);
			errval |= rsb_mtx_get_coo(mtxAp, VA.data(), IA.data(), JA.data(), RSB_FLAG_NOFLAGS);
		}
	}
	else
		mtxAp = rsb_mtx_alloc_from_coo_const(VA.data(), IA.data(), JA.data(), nnzA, typecode, nrA, ncA, 1, 1, flagsA, &errval);

	if(errval != RSB_ERR_NO_ERROR)
		return_at_line(report<NT>(errval,nrA));

	if( RSB_DO_FLAG_HAS(flagsA,RSB_FLAG_FORTRAN_INDICES_INTERFACE) )
		add(IA,rsb_coo_idx_t{-1}), add(JA,rsb_coo_idx_t{-1});

	if(bopts.wrl || bopts.want_rpp )
	{
		const coo_mtx_t<NT> cm(VA,IA,JA);
		const rsb_mtx_w_t<NT> mtxAw(mtxAp);

		auto soe = [&](void)
		{
			std::cout << mtxAw._info() << std::endl << std::endl; // FIXME: works only if RsbMatrix
			if( want_dump_vecs )
			{
				print_cxx_vec("IA",IA);
				print_cxx_vec("JA",JA);
				print_cxx_vec("VA",VA);
				for(auto k = 0; k<nnzA; ++k)
					std::cout << IA[k] << " " << JA[k] << " " << VA[k] << std::endl;
			}
			RSBT_SAVE_MAT("fail.mtx",mtxAp);
			return_at_line(errval);
		};

		if(bopts.wrl)
			errval = spmv_comparison_rectangular<NT>(bopts,nrA,nrA,mtxAw,cm,mr,soe);

#ifdef HAVE_RSBPP_HPP
		if( ! bopts.want_rpp )
			return errval;
		//if ( flagsA & ( RSB_FLAG_SYMMETRIC|RSB_FLAG_HERMITIAN ) )
		//if ( flagsA & ( RSB_FLAG_HERMITIAN ) )
		if ( false )
		{
			std::cout << "WARNING: Skipping RsbPP_Matrix_T<NT> because of symmetry." << std::endl;
		}
		else
		{
			/* RsbPP_Matrix_T<NT>::spmv is still unfinished */
			RsbPP_Matrix_T<NT> mtxAppp(VA,IA,JA,flagsA,RSBPP_MATRIX_T_ORDERING);
			errval = spmv_comparison_rectangular<NT>(bopts,n,n,mtxAppp,cm,mr,soe);
			if(errval != RSB_ERR_NO_ERROR)
			{
				mtxAppp.print(std::cout);
			}
		}
#endif /* HAVE_RSBPP_HPP */
	}
	
	return errval;
} /* band_tester */

static 
rsb_err_t band_testers(const rsbt_bopts & bopts, const std::string & ops, const rsbt_topts topts)
{
	const auto ts{bopts.ts};
	const auto n{topts.n};
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	//const std::vector<rsb_coo_idx_t> bs={1};
	//const std::vector<rsb_coo_idx_t> bs={0,1,n/2,n-1};
	//const std::vector<rsb_coo_idx_t> bs={0,1};
	const std::vector<rsb_coo_idx_t> bs={0,1,2,3,10};
	rsb_coo_idx_t mr = bopts.mor;

	if(ops.find('r')==std::string::npos)
	{
		if(ops.find('f')!=std::string::npos)
			mr = RSBT_FAIL_REC;
		else
			mr = 0; // no randomization
	}

	for(const auto csf_f : bopts.csf_fa)
	for(const auto dnst_pct : bopts.dnst_pct_a)
	if (! ( dnst_pct < 100 && n == 1 )) // avoid [0] matrix
	for(const rsb_coo_idx_t b : bs)
	for(const char          c : ts)
	for(const rsb_flags_t flags: bopts.symfa )
	for(const auto bsf_f : bopts.bsf_fa)
	if( ! ( bsf_f && RSB_DO_FLAG_HAS_INTERSECTION(flags,RSB_FLAG_SYMMETRIC|RSB_FLAG_HERMITIAN) ) ) /* need to respect given triangle boundaries.. */
	if( ! RSBT_ENOUGH_TEST_CASES )
	{
		if(b==0 && flags!=RSB_FLAG_NOFLAGS)
			continue;
		if ( !bopts.want_any_herm && !RSB_IS_MATRIX_TYPE_COMPLEX(std::toupper(c)) && RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN) )
			continue;
		switch(std::toupper(c))
		{
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
			case 'L':
				errval |= band_tester<long double>(bopts,topts,b,flags,mr,dnst_pct,bsf_f,csf_f );
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
			case 'D':
				errval |= band_tester<double>(bopts,topts,b,flags,mr,dnst_pct,bsf_f,csf_f );
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
			case 'S':
				errval |= band_tester<float>(bopts,topts,b,flags,mr,dnst_pct,bsf_f,csf_f);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
			case 'Q':
				errval |= band_tester<std::complex<long double>>(bopts,topts,b,flags,mr,dnst_pct,bsf_f,csf_f);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
			case 'Z':
				errval |= band_tester<std::complex<double>>(bopts,topts,b,flags,mr,dnst_pct,bsf_f,csf_f);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
			case 'C':
				errval |= band_tester<std::complex<float>>(bopts,topts,b,flags,mr,dnst_pct,bsf_f,csf_f);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
			case 'I':
				errval |= band_tester<int>(bopts,topts,b,flags,mr,dnst_pct,bsf_f,csf_f);
			break;
#endif
			default:
			errval = RSB_ERR_GENERIC_ERROR;
		}
		if(errval != RSB_ERR_NO_ERROR)
			return errval;
	}

	return errval;
} /* band_testers */

static 
rsb_err_t rectangular_testers(const rsbt_bopts & bopts, const std::string & ops, const rsbt_topts topts)
{
	const auto ts{bopts.ts};
	const auto no_square = true;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const std::vector<rsb_coo_idx_t> bs={1,2,3,10};
	rsb_coo_idx_t mr = bopts.mor;

	if(ops.find('r')==std::string::npos)
	{
		if(ops.find('f')!=std::string::npos)
			mr = RSBT_FAIL_REC;
		else
			mr = 0; // no randomization
	}

	for(const rsb_coo_idx_t b : bs)
	if( topts.n != b || !no_square )
	for(const char          c : ts)
	for(const rsb_flags_t flags: {RSB_FLAG_NOFLAGS} )
	if( ! RSBT_ENOUGH_TEST_CASES )
	if(errval == RSB_ERR_NO_ERROR)
	{
		errval = RSB_ERR_GENERIC_ERROR;
		switch(std::toupper(c))
		{
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
			case 'L':
				errval = rectangular_tester<long double>(bopts,topts,b,flags,mr);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
			case 'D':
				errval = rectangular_tester<double>(bopts,topts,b,flags,mr);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
			case 'S':
				errval = rectangular_tester<float>(bopts,topts,b,flags,mr);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
			case 'Q':
				errval = rectangular_tester<std::complex<long double>>(bopts,topts,b,flags,mr);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
			case 'Z':
				errval = rectangular_tester<std::complex<double>>(bopts,topts,b,flags,mr);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
			case 'C':
				errval = rectangular_tester<std::complex<float>>(bopts,topts,b,flags,mr);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
			case 'I':
				errval = rectangular_tester<int>(bopts,topts,b,flags,mr);
			break;
#endif
		}
	}

	return errval;
} /* rectangular_testers */

#ifdef HAVE_GETOPT_H
struct rsbt_options_t{
  const char *name{nullptr};
  int has_arg;
  int *flag{nullptr};
  int val;
  const char *desc{nullptr}; /* this is ctb specific */
};

struct rsbt_options_t rsbt_options[] = {
	{"help", no_argument, nullptr, 'h'},
	{"version", no_argument, nullptr, 'V'},
	{"serial", no_argument, nullptr, 's'},
	{"verbose", no_argument, nullptr, 'v'},
	{"quiet", no_argument, nullptr, 'q', "quiet progress of test cases"},
	{"no-timings", no_argument, nullptr, 0x70726e74, "print no timings"}, /* prnt */
	{"rand", no_argument, nullptr, 0x77720a },
	{"cplusplus", no_argument, nullptr, '+' },
	{"types", required_argument, nullptr, 't' },
	{"inc", required_argument, nullptr, 0x696e630a },
	{"mul", required_argument, nullptr, 0x6d756c0a },
	{"min", required_argument, nullptr, 0x6d696e0a },
	{"max", required_argument, nullptr, 0x6d61780a },
	{"mkl", no_argument, nullptr, 0x6d6b6c0a },
	{"nrhs", required_argument, nullptr, 0x6e726873},
	{"incx", required_argument, nullptr, 0x696e6378},
	{"incy", required_argument, nullptr, 0x696e6379},
	{"density", required_argument, nullptr, 0x646e7374}, /* dnst */
	{"alpha", required_argument, nullptr, 0x616c6661},
	{"beta", required_argument, nullptr, 0x62657461},
	{"coo-stretch-factor", required_argument, nullptr, 0x636f7366}, /*cosf*/
	{"band-stretch-factor", required_argument, nullptr, 0x62617368}, /*basf*/
	{"nthreads", required_argument, nullptr, 'n'},
	{"maxtimes ", required_argument, nullptr, 0x6d617854},
	{"recover-failed", no_argument, nullptr, 0x7265636f},
	{"rectangular", no_argument, nullptr, 0x72656374},
	{"no-rectangular", no_argument, nullptr, 0x6e726374}, /* nrct */
	{"render-only", no_argument, nullptr, 0x72656e64},
	{"render", no_argument, nullptr, 0x72656e74},
	{"no-basename-render", no_argument, nullptr, 0x6e626e72},/* nbnr */
	{"mintimes ", required_argument, nullptr, 0x6d696e54},
	{"min_t", required_argument, nullptr, 0x6d696e74},
	{"rsb-tester", no_argument, nullptr, 0x7273620a },
	{"coo-tester", no_argument, nullptr, 0x636f6f0a },
	{"csr-tester", no_argument, nullptr, 0x6373720a },
	{"hcoo-tester", no_argument, nullptr, 0x68636f6f},
	{"hcsr-tester", no_argument, nullptr, 0x68637372},
	{"half-indices-tester", no_argument, nullptr, 0x68616c66},
	{"recursive-tester", no_argument, nullptr, 0x7265630a},
	{"diagonal-implicit", no_argument, nullptr, 0x696d6469}, /* imdi */
	{"fortran-indices", no_argument, nullptr, 0x666f6969}, /* foii */
	{"only-test-case-n", required_argument, nullptr, 0x6f74636e},/*otcn*/
	{"only-n-test-cases", required_argument, nullptr, 0x6f6e7463},/*ontc*/
	{"skip-n-test-cases", required_argument, nullptr, 0x736e7463},/*sntc*/
	{"skip-except-every-n-test-cases", required_argument, nullptr, 0x73656e74},/*sent*/
	{"skip-except-every-random-n-test-cases", required_argument, nullptr, 0x7365726e},/*sern*/
	{"skip-loading-if-more-nnz-matrices", required_argument, nullptr, 0x736c6d6},/*slmn*/
	{"skip-loading-if-less-nnz-matrices", required_argument, nullptr, 0x736c6e6e},/* slnn */
	{"max_t", required_argument, nullptr, 0x746d6178},/*tmax*/
	{"max-test-time", required_argument, nullptr, 0x6d617474},/*matt*/
	{"self-test", no_argument, nullptr, 0x73657465},/*sete*/
	{"no-rand", no_argument, nullptr, 0x6e77720a },
	{"max-random-val", required_argument, nullptr, 0x72616d61}, /* rama */
	{"extra-verbose-interface", no_argument, nullptr, 0x6576690a},
	{"skip-loading-symmetric-matrices",  no_argument, nullptr, 0x736c736d},/* slsm */
	{"skip-loading-unsymmetric-matrices",no_argument, nullptr, 0x736c756d},/* slum */
	{"skip-loading-hermitian-matrices",no_argument, nullptr, 0x736c686d},/* slhm */
	{"skip-loading-not-unsymmetric-matrices",no_argument, nullptr, 0x736c6e75},/* slnu */
	{"symmetric",  no_argument, nullptr, 0x73796d6d},/* symm */
	{"hermitian",  no_argument, nullptr, 0x6865726d},/* herm */
	{"unsymmetric",  no_argument, nullptr, 0x756e7379},/* unsy */
	//{   "leaf-level-multivec", no_argument, nullptr, 0x6c6c6d0a},
	{"no-leaf-multivec", no_argument, nullptr, 0x6e6c6d6d}, // nlmm
	{"verbose-tuning", no_argument, nullptr, 0x76657475},
	{"report", required_argument, nullptr, 0x72707274},/*rprt*/
	{"transA", required_argument, nullptr, 0x7472616e},/*tran*/
	{"no-trans", no_argument, nullptr, 0x6e6f7472},/*notr*/
	{"no-tune", no_argument, nullptr, 0x74756e65},/*tune*/
	{"tolerate-mismatch", no_argument, nullptr, 0x746f6d69},/*want_tmm*/
	{"any-hermitian-ok", no_argument, nullptr, 0x61686f6b},/*ahok*/
	{"tune-maxr", required_argument, nullptr, 0x6d617872},/*maxr*/
	{"tune-maxt", required_argument, nullptr, 0x6d617874},/*maxt*/
	{0,no_argument,0,0}
};
const int options_count=sizeof(rsbt_options)/sizeof(rsbt_options_t);
struct option options[options_count];
#else /* HAVE_GETOPT_H */
#endif /* HAVE_GETOPT_H */

template <typename NT=double, typename TT=double>
static 
rsb_err_t serial_test(const bool wv)
{
	// TODO: eliminate this test -- it's obsolete.
	// TODO: a good deal of duplicated code here.
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<TT> elapsed_seconds;
	rsb_lib rli{};
	
	if(wv) std::cout << "Will perform serial test.\n";

	const rsb_flags_t flagsA=RSB_FLAG_NOFLAGS;
	//const rsb_flags_t flags =RSB_FLAG_NOFLAGS;
	const rsb_nnz_idx_t n=100;
	const rsb_nnz_idx_t nrA=n,ncA=n;
	const rsb_nnz_idx_t nnzA=half_square(n);
	const std::vector<NT> cx(nrA,1);
	std::vector<NT> VA(nnzA,1);
	std::vector<rsb_coo_idx_t> IA(nnzA,0);
	std::vector<rsb_coo_idx_t> JA(nnzA,0);
	struct rsb_mtx_t * mtxAp = nullptr;
	const rsb_type_t typecode=rsb_type_t_for<NT>();
	const rsb_coo_idx_t b=0;
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;
	int nt = 1;

	mtxAp = rsb_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA, flagsA, &errval);

	NT val = 1;

	if(is_complex_type<NT>())
	{
		const NT mo = -1;
		const NT i = std::sqrt(mo);
		val += i;
	}

	rsb_nnz_idx_t k = 0;
	for(auto i = k; i<n; ++i)
		for(auto j = i; j<std::min(i+b+1,n); ++j)
		{
			VA[k]=val,
			VA[k]+=get_rci(0);
//			if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER))
//				IA[k]=i, JA[k++]=j;
//			else
				JA[k]=i, IA[k++]=j;
		}

	errval = rsb_mtx_set_vals(mtxAp, VA.data(), IA.data(), JA.data(), nnzA, RSB_FLAG_NOFLAGS);
	if(errval != RSB_ERR_NO_ERROR)
		return_at_line(report<NT>(errval,n));

	errval = rsb_mtx_alloc_from_coo_end(&mtxAp);
	if(errval != RSB_ERR_NO_ERROR)
		return_at_line(report<NT>(errval,n));

	start = std::chrono::system_clock::now();
#pragma omp parallel
#pragma omp master
{
	nt = omp_get_num_threads();
}
	if(wv) std::cout << "Will run  " << nt << " threads loose." << std::endl;
	if(wv) std::cout << "Using a " << nrA << " x " << ncA << " matrix with " << nnzA << " nonzeroes." << std::endl;
//#pragma omp parallel
//#pragma omp master
{
	for(int i=0; i<omp_get_num_threads(); i++)
	{
//#pragma omp task
{
		//errval = band_testers(.. ,ops,10);
		//if(errval != RSB_ERR_NO_ERROR)
		//	return errval;
		const rsb_coo_idx_t incX=1,incY=1;
		const int nit = 100;
		std::vector<NT> alpha_a={-4,-1,1,4};
		std::vector<NT> beta_a={1};
		std::vector<NT> rx(nrA,1);
		std::vector<NT> ry(ncA,1);

		if(wv) std::cout << "Thread " << i << " executing " << nit << " iterations..." << std::endl;
		for(int r=0;r<nit;++r)
			for (NT alpha : alpha_a)
				for (NT beta  : beta_a )
					rsb_spmv(transA, &alpha, mtxAp, rx.data(), incX, &beta, ry.data(), incY);
		if(wv) std::cout << "Thread " << i << " done." << std::endl;
}
	}
}
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	TT st = elapsed_seconds.count();
	if(wv) std::cout << "Done in " << st << " s." << std::endl;
	std::cout << "Ran with " << nt << " threads using a " << nrA << " x " << ncA << " matrix with " << nnzA << " nonzeroes in " << st << " s." << std::endl;
	rsb_mtx_free(mtxAp);
	return errval;
	// might use a tuple to encode information about what each task shall execute.
} /* serial_test */

struct rsbt_mtx_ptn_t
{
       	std::vector<rsb_coo_idx_t> IA{},JA{};
	rsb_flags_t flags{RSB_FLAG_NOFLAGS};
	rsb_nnz_idx_t nnz{0}; 
	void reserve(size_t nnz_){IA.reserve(nnz_); JA.reserve(nnz_); nnz=nnz_;}
};

static
bool has_extension( const std::string s, const std::string e )
{
	const auto es = e.size();
	return ( s.size() > es && 0==s.compare(s.size()-es,es,e));
}

#ifdef RSBT_WANT_RSB_HPP
static
std::string strip_mm_ext ( const std::string s )
{
	if (has_extension(s,".mtx"))
		return s.substr(0,s.size()-4);
	if (has_extension(s,".mtx.gz"))
		return s.substr(0,s.size()-7);
	return s;
}

template <typename NT, typename ST>
static
void render_matrix_file(std::string rn, const RsbMatrix<NT> & mtx_A_base, ST &os)
{
	rsb_err_t errval;
	//rsb_marf_t rflags{RSB_MARF_EPS};
	rsb_marf_t rflags{RSB_MARF_EPS_L};
	rsb_coo_idx_t pmWidth{512}, pmHeight{512};
	os << "Rendering matrix into " << rn << std::endl;
#if 1
	errval=mtx_A_base.rndr(rn.c_str(), pmWidth, pmHeight, rflags);
#else
	errval=rsb_mtx_rndr(rn.c_str(), mtx_A_base.mtxAp_, pmWidth, pmHeight, rflags);
#endif
	if(errval!=RSB_ERR_NO_ERROR)
	{
		rsb_perror(stdout, errval);
		throw;
	}
}
#endif /* RSBT_WANT_RSB_HPP */

static
std::string str_to_tex(const std::string & s)
{
	return std::regex_replace(s, std::regex("_"), "\\_");
}

#ifdef RSBT_WANT_RSB_HPP
template<typename NT>
static
rsb_err_t hpp_test(rsbt_bress & bress, const rsbt_bopts & bopts, const std::string & ops, const rsbt_mtx_ptn_t & mp, const std::string & fn )
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RsbLib rsb_lib;
	std::ostringstream ors; // report stream
	const rsb_coo_idx_t mr = bopts.mor;

	bopts.dro.set();

	try
	{
	for (int bv= 0; bv< 1; ++bv)
	for(rsb_int_t tn: bopts.tna)// TODO: can be moved further  inside
	{
		//bopts.dro.get_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn);// TODO: need rsblib own wrapper
		bopts.dro.set_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn);// TODO: need rsblib own wrapper
		bopts.dro.get_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn);// TODO: need rsblib own wrapper
		if(tn==0)
		{
			if(RSB_LIBRSB_VER>=10300) // 
				throw;
			else
				{}//oss << "librsb using as many threads as in OpenMP\n";
#ifdef _OPENMP
			tn = omp_get_max_threads();
#else /* _OPENMP */
			tn = -1;
#endif /* _OPENMP */
		}
#ifdef RSBT_USE_MKL
		if(tn>0)
			mkl_set_num_threads(tn);
#endif /* RSBT_USE_MKL */
		// TODO: need cross-type conversion, pattern-based
		//std::cout << fn << std::endl;
		//TODO: spmv time and results now; then librsbppp benchmark, transposition
		//std::cout << "loading " << fn << " using C++ interface to librsb" << std::endl;
		//RsbMatrix<NT> mtx_A_base(fn.c_str());
		//RsbMatrix<NT> mtx_A_base("fn.c_str()");
		//TODO: conversion/copy time
		std::string sn{strip_mm_ext(fn)};// stripped name
		std::string mn{sn.find_last_of('/')==std::string::npos ? sn : sn.substr(1+sn.find_last_of('/'))};//
		if(bopts.want_basenameplot)
			sn=mn;// basename by default 
		std::ostringstream oss;
		std::ostringstream vcptn; // vanilla caption
		std::vector<NT>VA(mp.nnz);
		fill_rhs(VA,bv,mr);
		rsb_time_t dt = rsb_time();
		const RsbMatrix<NT> mtx_A_base(mp.IA.data(),mp.JA.data(),VA.data(),mp.nnz,mp.flags);
		dt = rsb_time() - dt;
		oss.precision(3);
		oss << "Matrix converted from pattern in " << dt << " s (" << 1.0/(dt/mp.nnz)<< " nnz/s)" << std::endl;
		std::cout << oss.str();
		ors << oss.str();
		oss = std::ostringstream{};

		vcptn << mn << ": " << mtx_A_base.rows() << " rows, " << mtx_A_base.cols() << " columns, " << mtx_A_base.nnz() << " nnz";
		//oss << mtx_A_base << std::endl;
		oss << mtx_A_base._info() << std::endl << std::endl;
#ifdef HAVE_RSBPP_HPP
		rsb_flags_t flagsA=RSB_FLAG_NOFLAGS;
		mtx_A_base.get_info(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &flagsA);
#endif /* HAVE_RSBPP_HPP */

		std::string tm,ts="----------------------";

		for(const NT alpha: bopts.alphaia)
		for(const NT beta:  bopts.betaia )
		for(const auto incX: bopts.incxa)
		for(const auto incY: bopts.incya)
		for(const rsb_coo_idx_t nrhs: bopts.nrhsa)
		for(const rsb_trans_t transA: bopts.transAa)
		for(const rsb_flags_t order : bopts.ordera)
		if( ! ( nrhs == 1 && bopts.ordera.size() == 2 && order == bopts.ordera[1] ) )
		{
			if(!(nrhs==1 || (incX*incY==1)))
				{ oss << "Skipping: incX>1 && (incX*incY>=1)!" << std::endl ; continue;}
			std::ostringstream dscs{};
			dscs<<" t:"; dscs<<rsb_type_t_for<NT>();
			dscs<<" T:"; dscs<<RSB_TRANSPOSITION_AS_CHAR(transA);
			dscs<<" s:"; dscs<<char_for_flag(mp.flags);
			dscs<<" a:"; dscs<<alpha;
			dscs<<" b:"; dscs<<beta;
			dscs<<" n:"; dscs<<nrhs;
			dscs<<" x:"; dscs<<incX;
			dscs<<" y:"; dscs<<incY;
			dscs<<" N:"; dscs<<tn;
			std::string dsc{(str_to_tex(mn) + ": ") + dscs.str()};
			std::string tag{(sn + "--") + dscs.str()};
			tag=std::regex_replace(tag, std::regex(":"), "");
			tag=std::regex_replace(tag, std::regex(" "), "");

			dt = rsb_time();
			RsbMatrix<NT> mtx_A{mtx_A_base};
			//RsbMatrix<NT> mtx_A(mp.IA.data(),mp.JA.data(),VA.data(),mp.nnz,mp.flags);
			dt = rsb_time() - dt;
			oss.precision(3);

			tm = "Testing C++ interface to librsb";
			oss << ts << std::endl;
			oss << tm << "... : BEGIN" << std::endl;

			//oss << "Matrix converted from COO to RSB [" << mtx_A.blocks()<< " subm] in " << dt << " s (" << 1.0/(dt/mp.nnz)<< " nnz/s )" << std::endl;
			//oss << "Matrix cloned from RSB to RSB in " << dt << " s (" << 1.0/(dt/mp.nnz)<< " nnz/s )" << std::endl;
			oss << "Matrix cloned from RSB [" << mtx_A_base.blocks()<< " subm] to RSB [" << mtx_A.blocks()<< " subm] ("
						<< mtx_A_base._info() << ")" << " in " << dt << " s (" << 1.0/(dt/mp.nnz)<< " nnz/s )"
						<< std::endl;
			oss << ts << std::endl << "Testing specifying "; 
			if(tn>0)
				oss << tn << " threads,";
			else
				oss << "auto" << " threads,";
			oss     <<  " alpha=" << alpha
				<< ", beta=" << beta
				<< ", type=" << rsb_type_t_for<NT>()
				<< ", nrhs=" << nrhs 
				<< ", incX=" << incX 
				<< ", incY=" << incY 
				<< ", transA=" << RSB_TRANSPOSITION_AS_CHAR(transA)
				<< std::endl;

			if(bopts.ropts==rsbt_bopts::RenderOpts::All)
			{
				auto pfn = tag+"-untuned.eps";
				render_matrix_file(pfn,mtx_A,oss),
				bress.lrpv.push_back(rsbt_bress::lrpi_t{ dsc+" before tuning",pfn,str_to_tex(vcptn.str())});
			}

			// rename to hpp_bench
			const auto trows = ( transA == RSB_TRANSPOSITION_N ) ? mtx_A.rows() : mtx_A.cols();
			const auto tcols = ( transA != RSB_TRANSPOSITION_N ) ? mtx_A.rows() : mtx_A.cols();
			const rsb_nnz_idx_t ldB = ( order == RSB_FLAG_WANT_ROW_MAJOR_ORDER ) ? nrhs : tcols;
			const rsb_nnz_idx_t ldC = ( order == RSB_FLAG_WANT_ROW_MAJOR_ORDER ) ? nrhs : trows;
			std::vector<NT> x(tcols*incX*nrhs);
		       	const NT * Xp{x.data()};
			std::vector<NT> y(trows*incY*nrhs);
		       	NT * Yp{y.data()};
	       		const NT * Bp=Xp;
			NT * Cp=Yp;

			const rsb_blk_idx_t nb = mtx_A.blocks();
			rsb_real_t tsf;
			rsbt_btrs cmtr,cftr,attr; // 
			rsbt_btrs pmtr,pftr,atts; // 

			assert(nrhs==1 || (incX*incY==1)); // TODO: FIXME

			fill_rhs(x,bv,mr);
			fill_rhs(y,bv,mr);

			std::cout << oss.str();
			ors << oss.str();
			oss = std::ostringstream{};

			if(nrhs==1)
			{
				timed_exec_m({ errval |= mtx_A.spmv(transA, alpha, Xp, incX, beta, Yp, incY);},bopts.dro,cmtr);
				timed_exec_f({ errval |= mtx_A.spmv(transA, alpha, Xp, incX, beta, Yp, incY);},bopts.dro,cftr);
			}
			else
			{
				timed_exec_m({ errval |= mtx_A.spmm(transA, alpha, nrhs, order, Bp, ldB, beta, Cp, ldC);},bopts.dro,cmtr);
				timed_exec_f({ errval |= mtx_A.spmm(transA, alpha, nrhs, order, Bp, ldB, beta, Cp, ldC);},bopts.dro,cftr);
			}

			if( bopts.want_autotune )
			{
				// warning: the following is memory hungry
				timed_exec_macro({ 
					rsb_int_t *tnp=nullptr; /* nullptr -> current threads count */
					// TODO: is the following taking care of incX, incY ?
					errval |= mtx_A.tune_spmm(&tsf, tnp, bopts.maxr, bopts.maxt, transA, &alpha, nrhs, order, Bp, ldB, &beta, Cp, ldC);// FIXME: Bp and Cp are declared void in tune_spmm ?!
				},1,1,0.0,0.0,atts);
				// TODO: might rely on tune_spmm solely.
				if(nrhs==1)
				{
					timed_exec_m({ errval |= mtx_A.spmv(transA, alpha, Xp, incX, beta, Yp, incY);},bopts.dro,attr);
				}
				else
				{
					timed_exec_m({ errval |= mtx_A.spmm(transA, alpha, nrhs, order, Bp, ldB, beta, Cp, ldC);},bopts.dro,attr);
				}
			} /* bopts.want_autotune */
			if(errval!=RSB_ERR_NO_ERROR)
				throw;

			std::string opid;
			if(nrhs==1)
				opid="spmv";
			else
				opid="spmm";
			std::string os{"RsbMatrix<NT> " + opid};

			// TODO: per-op/per-type flops count.
			const auto op_nnz=mtx_A.nnz()*nrhs;// TODO: symmetry? *2 ?!
			const auto tol {1.001};
			std::ostringstream cptn; // vanilla caption
			cptn=std::ostringstream{};
			cptn.precision(3);
			if(attr.avt*tol<cmtr.avt)
				cptn
					<< cmtr.avt << " s/" << opid << " (" << 1.0/(cmtr.avt/op_nnz) << " nnz/s)" 
					<< nb << " blocks"
					<< " $\\rightarrow$ "
					<< attr.avt << " s/" << opid << " (" << 1.0/(attr.avt/op_nnz) << " nnz/s)"
				 	<< ", " << mtx_A.blocks() << " blocks"
				 	<< ": " << cmtr.avt/attr.avt << "x speedup"
				 	<< " (" << tsf << "x reported by rsb tuner)"
					;
			else
				cptn
					<< cmtr.avt << " s/" << opid << " (" << 1.0/(cmtr.avt/op_nnz) << " nnz/s)" 
					<< ", tuning time (" << atts.avt << " s) seems lost"; // TODO: duplication; TODO: equivalent as how many operations ?

			if(bopts.ropts==rsbt_bopts::RenderOpts::All)
			{
				const auto pfn = tag+"-tuned.eps";
				render_matrix_file(pfn,mtx_A,oss);
				bress.lrpv.push_back(rsbt_bress::lrpi_t{ dsc+" after tuning",pfn,cptn.str()});
			}

			oss.precision(3);
			oss << os << " took " << cmtr << " s, " << 1.0/(cmtr.avt/op_nnz) << " nnz/s (macro)" << std::endl;
			oss << os << " took " << cftr << " s, " << 1.0/(cftr.avt/op_nnz) << " nnz/s (lambda)" << std::endl;
			if( bopts.want_autotune )
			{
				oss << os << " tuning took "<< atts.avt << " s, as " << atts.avt/cmtr.avt                      << " ops" << std::endl;
				oss << os << " took " << attr << " s, " << 1.0/(attr.avt/op_nnz) << " nnz/s (tuned)" << std::endl; 
				oss.precision(4);
				oss << os << ": lambda -> macro speedup is: " << cftr.avt/cmtr.avt << " x " << std::endl;
				oss << os << ": default -> tuned speedup is: " << cmtr.avt/attr.avt << " x, by changing " 
					<< nb << " -> " << mtx_A.blocks() << " blocks" << std::endl;
				if(attr.avt*tol<cmtr.avt)
					oss << "Tuning time (" << atts.avt << " s) can be amortized in " << atts.avt/(cmtr.avt-attr.avt) << " ops; afterwards, will save " << (cmtr.avt-attr.avt) << "s per op " << std::endl;
				else
					oss << "Tuning time (" << atts.avt << " s) seems lost" << std::endl;
			}
			oss << tm << " : END" << std::endl << ts << std::endl;
			std::cout << oss.str();
			ors << oss.str();
			oss = std::ostringstream{};

			if(ops.find('=')!=std::string::npos)
			{
				coo_mtx_t<NT> coo_A{mtx_A}; // via librsb: can handle *.gz

				oss << "Checking if results match " << "...";
				errval = spmv_comparison(transA, mtx_A, coo_A, order, x, y, alpha, beta, nrhs, incX, incY);
				if(errval != RSB_ERR_NO_ERROR)
				{
					oss << "a numeric mismatch occurred!" << std::endl;
					std::exit(-1);
					return {};
				}
				oss << "OK" << std::endl;
			}

			if(bopts.want_mkl)
			if(nrhs==1)
			if(incX==1)
			if(incY==1)
			{
#ifdef RSBT_USE_MKL_SPBLAS
				//tm = "Testing Intel MKL implementation";
				//oss << tm << "...";
				coo_mtx_t<NT> coo_A{mtx_A}; // via librsb: can handle *.gz
				rsbt_btrs mbtr;
				test_mkl(mbtr,alpha,beta,coo_A.VA,coo_A.IA,coo_A.JA,mp.flags,mtx_A.rows(),mtx_A.cols(),/*trows,tcols,*/transA,bopts.dro,x,y);
				oss << "MKL mkl_sparse_*_mv [csr] took " << mbtr << " s, " << 1.0/(mbtr.avt/op_nnz) << " nnz/s (macro)" << std::endl;
				//oss << std::endl;
#endif
			}

			if( bopts.want_rpp )
			{
#ifdef HAVE_RSBPP_HPP
				tm = "Testing C++ RSB implementation";
				oss << tm << "...";

				if(nrhs!=1){ oss << "skipping: nrhs>1!" << std::endl ; continue;}
				if(incX!=1){ oss << "skipping: incX>1!" << std::endl ; continue;}
				if(incY!=1){ oss << "skipping: incY>1!" << std::endl ; continue;}
				assert(incX==1); // TODO: FIXME
				assert(incY==1); // TODO: FIXME
				assert(nrhs==1); // TODO: FIXME

				//coo_mtx_t<NT> coo_A{fn}; // own: cannot handle *.gz
				coo_mtx_t<NT> coo_A{mtx_A}; // via librsb: can handle *.gz
				// TODO: shall use move semantics to empty mtx_A.
				//oss << coo_A << std::endl;

				if(coo_A.get_nnz() != mtx_A.nnz())
				{
					oss << coo_A.get_nnz() << " != " << mtx_A.nnz() << "\n";
					throw;
				}
				oss << " : BEGIN" << std::endl;

				std::cout << oss.str();
				ors << oss.str();
				oss = std::ostringstream{};
				RsbPP_Matrix_T<NT> mtxAppp(coo_A.VA,coo_A.IA,coo_A.JA,flagsA,RSBPP_MATRIX_T_ORDERING);

				timed_exec_m({mtxAppp.spmv(transA, x, y, alpha, beta);},bopts.dro,pmtr);
				timed_exec_f({mtxAppp.spmv(transA, x, y, alpha, beta);},bopts.dro,pftr);

				oss.precision(3);
				// TODO: per-op/per-type flops count.
				oss << "RsbPP_Matrix_T<NT> spmv took " << pmtr << " s, " << 1.0/(pmtr.avt/op_nnz) << " nnz/s (macro)" << std::endl;
				oss << "RsbPP_Matrix_T<NT> spmv took " << pftr.avt << " s, " << 1.0/(pftr.avt/op_nnz) << " nnz/s (lambda)" << std::endl;
				oss << "RsbPP_Matrix_T<NT> spmv: lambda -> macro speedup is: " << pftr.avt/pmtr.avt << " x " << std::endl;
				// TODO: need sentinel check for 'speed of light' spmv and error in case

				std::cout << oss.str();
				ors << oss.str();
				oss = std::ostringstream{};

				oss.precision(3);
				oss << " librsbpp -> librsb speedup is: " << pmtr.avt/cmtr.avt << " x (macro)" << std::endl;
				oss << " librsbpp -> librsb speedup is: " << pftr.avt/cftr.avt << " x (lambda) " << std::endl;
				//std::cout << " p.s.:" << y[0] << std::endl;
				oss << oss.str();
				oss << tm << " : END" << std::endl << ts << std::endl;
#else /* HAVE_RSBPP_HPP */
				oss << "will skip librsbpp testing: not built in" << std::endl;
#endif /* HAVE_RSBPP_HPP */
				std::cout << oss.str();
				ors << oss.str();
				oss = std::ostringstream{};
			}
			ors << oss.str();
		}
	}
	}
	catch(std::exception & e)
       	{
		std::cout << "Caught an exception!\n" << std::endl;
	}
	bress.rs << ors.str();
	return {};
} /* hpp_test */
#endif /* RSBT_WANT_RSB_HPP */

template <typename T>
static
T suffixed_val_read(std::istringstream & iss)
{
	T val;
	if(!isdigit(iss.peek()))
	{
		if(!iss.eof() && iss.peek()=='-' && iss.get() == '-' && isdigit(iss.peek()))
			iss.unget();// TODO: this shall be only allowed if type is signed
		else
			throw std::exception{};
	}
	iss >> val;
	while( !iss.eof() && isalpha(iss.peek()) )
	{
		switch(const char c=iss.get())
		//switch(const char c=tolower(iss.get()))
		{
			case 'k':val*=1000;break;
			case 'm':val*=1000*1000;break;
			case 'g':val*=1000*1000*1000;break;
			// yes I know is out of any standard but..
			case 'K':val*=1024;break;
			case 'M':val*=1024*1024;break;
			case 'G':val*=1024*1024*1024;break;
			default: throw std::exception{};
		}
	}
	return val;
}

template <typename T>
static
void suffixed_val_read(std::string s, T & v)
{
	std::istringstream iss(s);
	v = suffixed_val_read<T>(iss);
}

template <typename T>
static
std::vector<T> vals_array(std::string optarg)
{
	std::vector<T> vec{};
	T val;
	std::istringstream iss(optarg);
	while( isprint(iss.peek()) && ! iss.eof() )
	{
		try { val = suffixed_val_read<T>(iss); }
		catch ( ... ) { throw; }
		if( iss.peek() == ',' )
			iss.get();
		if( iss.peek() == ':' )
		{
			T max,inc=0,mul=1;
			char pc {'+'};
			iss.get();
			try { max = suffixed_val_read<T>(iss); }
			catch ( ... ) { throw; }
			if( iss.peek() == ':' )
			{
				iss.get();
				pc = iss.peek();
				if(!isdigit(pc))
				{
					if(pc!='+' && pc!='*')
						break;
					iss.get();
				}
				else 
					pc='+';
				try { mul = suffixed_val_read<T>(iss); } catch ( ... ) { break; }
			}
			switch( pc )
			{
				case( '*' ): mul = mul, inc = 0; break;
				case( '+' ): inc = mul, mul = 1; break;
			} 
			//if( std::abs(mul + inc) == 0 || int(std::abs(mul + inc)) == 0)
			//	break;
					
			for(;std::abs(val)<=std::abs(max);val=val*mul+inc)
				// std::cout << "adding " << val << " = " << " val*" << mul << " + " << inc << std::endl,
				vec.push_back(val);
		}
		else
			vec.push_back(val);
	}
	if( ! iss.eof() )
	{
		std::cout << "Warning: last part of the string (\"" << std::string(iss.str()).substr(iss.tellg()) << "\") not recognized and therefore ignored!" << std::endl;
		throw std::exception{};
	}
	return vec;
}

#if HAVE_FILESYSTEM && USE_CXX17
static
bool is_mtx_market_filename( std::string s )
{
	return ( has_extension(s,".mtx") || has_extension(s,".mtx.gz") );
}
#endif /* HAVE_FILESYSTEM && USE_CXX17 */

static
bool rsbt_passes_filter( const std::string & fn, const rsb_lib_opts ro)
{
	// TODO: might print once at the end and collect reasons along the way
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_coo_idx_t nr,nc;
	rsb_nnz_idx_t nnz;
	rsb_flags_t flags;

	errval = rsb_file_mtx_get_dims(fn.c_str(), &nr, &nc, &nnz, &flags);
	if( errval != RSB_ERR_NO_ERROR )
	{
		std::cout << "skip " << fn << " : problems opening ?!" << std::endl;
		return false;
	}
	if(RSB_DO_FLAG_HAS_INTERSECTION(flags,ro.sf))
	{
		std::cout << "skip " << fn << " : bad flags !" << std::endl;
		return false;
	}
	if(!RSB_DO_FLAG_HAS_INTERSECTION(flags,ro.af))
	{
		std::cout << "skip " << fn << " : missing required flags !" << std::endl;
		return false;
	}
	if( nnz > ro.max_nnz )
	{
		std::cout << "skip " << fn << " : exceeding " << ro.max_nnz << " nnz!" << std::endl;
		return false;
	}
	if( nnz < ro.min_nnz )
	{
		std::cout << "skip " << fn << " : not reaching " << ro.min_nnz << " nnz!" << std::endl;
		return false;
	}
	std::cout << "adding " << fn << " (" << nr << " x " << nc << ", " << nnz << " nnz, " << char_for_flag(flags) << trichar_for_flag(flags) << ")" <<  std::endl;

	return true;
} /* rsbt_passes_filter */

template<typename PT>
std::ostream& operator<<(std::ostream &os,const std::vector<PT> & tv)
{
	for (const auto e : tv)
		os << e << " ";
	return os;
}

#ifdef HAVE_SYS_UTSNAME_H 
std::string get_utsname(void)
{
	struct utsname un;
	std::string utsname;
	if(uname(&un)==0)
		utsname = un.nodename;
	return utsname;
}
#endif /* HAVE_SYS_UTSNAME_H */

template <typename NT>
static 
int vec_self_test(void)
{
	const std::vector<NT> a(3,0), b(3,1), c{1,0,1};

	int failed { 0 };

	failed += ( ! any_diff(a,b) ) ;
	failed += (   any_diff(b,b) ) ;
	failed += ( ! any_diff(b,c,1) ) ;
	failed += (   any_diff(b,c,2) ) ;

	print_cxx_vec("test",a);
	print_cxx_var("beta",a[0]);
;
	failed += ( std::abs(print_if_diff(b,b)) );
	failed += ( ! std::abs(print_if_diff(b,a)) );
	failed += ( ! std::abs(print_if_diff(b,a,1)) );
;
	failed += ( vals_array<NT>("1:1") != std::vector<NT>{1} );
	failed += ( vals_array<NT>("1:3") != std::vector<NT>{1,2,3} );
	failed += ( vals_array<NT>("1:6:+4") != std::vector<NT>{1,5} );
	failed += ( vals_array<NT>("1:6:*2") != std::vector<NT>{1,2,4} );
	failed += ( vals_array<NT>("1:6:2") != std::vector<NT>{1,3,5} );
	failed += ( vals_array<NT>("1k:1001") != std::vector<NT>{1000,1001} );
	failed += ( vals_array<NT>("1K:1024") != std::vector<NT>{1024} );
	failed += ( vals_array<NT>("1m:1m") != std::vector<NT>{1000*1000} );
	failed += ( vals_array<NT>("2M:2M") != std::vector<NT>{2*1024*1024} );

	// notice the following two stay specialized:
	failed += ( vals_array<rsb_int_t>("1g:1g") != std::vector<rsb_int_t>{1000*1000*1000} );
	failed += ( vals_array<rsb_int_t>("1G:1G") != std::vector<rsb_int_t>{1024*1024*1024} );

	try { std::istringstream iss("1z");suffixed_val_read<NT>(iss); failed ++; } catch ( ... ) { }
	try { std::istringstream iss("z1");suffixed_val_read<NT>(iss); failed ++; } catch ( ... ) { }
	try { vals_array<NT>("1z"); failed ++; } catch ( ... ) { }
	try { vals_array<NT>("1g:1z"); failed ++; } catch ( ... ) { }
	try { vals_array<NT>("1k:1k:/4"); failed ++; } catch ( ... ) { }

	failed += ( rsb_type_s_for<NT>().size() < 3 );
	failed += ( report(RSB_ERR_GENERIC_ERROR) != RSB_ERR_GENERIC_ERROR );
	failed += ( report(RSB_ERR_NO_ERROR) != RSB_ERR_NO_ERROR);

	failed += rsb_type_t_for<NT>() == rsb_type_t_for<void>();

	std::cout << "self-test " << rsb_type_s_for<NT>() << ( failed ? " FAILED" : " PASSED" ) << "\n";

	return failed;
}

#ifdef RSBT_WANT_RSB_HPP
void hpp_tests(const rsbt_bopts & bopts, rsbt_bress & bress, const std::string & ops)
{
	rsb_err_t retval = RSB_ERR_NO_ERROR;
		for(const auto & fn : bopts.fna)
		{
			rsbt_mtx_ptn_t mp{};
			{
				RsbLib rsb_lib;
				bopts.dro.set();
				std::cout << "loading " << fn << " using C++ interface to librsb" << std::endl;
				rsb_time_t dt = rsb_time();
				RsbMatrix<RSB_DEFAULT_TYPE> mtx_A_base(fn.c_str());//FIXME
				dt = rsb_time() - dt;
				std::cout << fn << " loaded in " << dt << "s" << std::endl;
				if(bopts.ropts==rsbt_bopts::RenderOpts::Once)
				{
					render_matrix_file(strip_mm_ext(fn)+".eps",mtx_A_base,std::cout);
					break;
				}
				mp.reserve(mtx_A_base.get_info_nnz_t(RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T));
				std::vector<RSB_DEFAULT_TYPE>VA(mp.nnz,1);
				retval |= mtx_A_base.get_rows_sparse(RSB_TRANSPOSITION_N,nullptr,VA.data(),mp.IA.data(),mp.JA.data(),0,mtx_A_base.rows()-1,&mp.nnz,/*RSB_FLAG_DEFAULT_STORAGE_FLAGS*/RSB_FLAG_NOFLAGS);
				retval |= mtx_A_base.get_info(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T,&mp.flags);
			}

			for(const char c: bopts.ts)
			switch(std::toupper(c))
			{
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
			case 'D':
				retval |= hpp_test<double>(bress,bopts,ops,mp,fn);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
			case 'S':
				retval |= hpp_test<float>(bress,bopts,ops,mp,fn);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
			case 'Z':
				retval |= hpp_test<std::complex<double>>(bress,bopts,ops,mp,fn);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
			case 'C':
				retval |= hpp_test<std::complex<float>>(bress,bopts,ops,mp,fn);
			break;
#endif
#ifdef RSB_NUMERICAL_TYPE_INT
			case 'I':
				retval |= hpp_test<int>(bress,bopts,ops,mp,fn);
			break;
#endif
			default:
			retval = RSB_ERR_GENERIC_ERROR;
			}
		}
}
#endif /* RSBT_WANT_RSB_HPP */

int self_test(void)
{
	return 
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE
		vec_self_test<long double>()|
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
		vec_self_test<double>()|
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT
		vec_self_test<float>()|
#endif
#ifdef RSB_NUMERICAL_TYPE_LONG_DOUBLE_COMPLEX
		vec_self_test<std::complex<long double>>()|
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
		vec_self_test<std::complex<double>>()|
#endif
#ifdef RSB_NUMERICAL_TYPE_FLOAT_COMPLEX
		vec_self_test<std::complex<float>>()|
#endif
		vec_self_test<int>();
}

void eps_report(const std::string & rfn, const rsbt_bress & bress)
{
	std::ofstream ostrm(rfn);
	auto rstr = bress.rs.str();

	if (has_extension(rfn,".tex"))
	{
		//rstr = std::regex_replace(rstr, std::regex("^"), "%$&", std::regex::ECMAScript|std::regex::multiline);
		std::string pfx = "% ";
		rstr = "%%" + std::regex_replace(rstr, std::regex("\n"), "$&" + pfx);
		std::string pcs = "% for f in ";
		for ( const auto lrpi : bress.lrpv )
			pcs += str_to_tex(lrpi.pfn) + " ";
		pcs += "; do epstopdf $f; done;";
		rstr = ( pcs + std::regex_replace(ltx_preamble, std::regex("EXTRA-TITLE"), str_to_tex(library_version_string()))) + rstr;
		for ( const auto lrpi : bress.lrpv )
		{
			if(lrpi.pfn.size())
			{
				std::string pg{ltx_pgtmplt};
				pg = std::regex_replace(pg, std::regex("FRAME-TITLE"), lrpi.title);
				pg = std::regex_replace(pg, std::regex("CAPTION"), lrpi.caption);
				pg = std::regex_replace(pg, std::regex("EPS-PLOT-FILE"), lrpi.pfn);
				rstr += pg;
			}
			else
			if(lrpi.text.size())
			{
				std::string pg{ltx_fgtmplt};
				pg = std::regex_replace(pg, std::regex("FRAME-TITLE"), lrpi.title);
				pg = std::regex_replace(pg, std::regex("FRAME-TEXT"), lrpi.text);
				rstr += pg;
			}
		}
		rstr += ltx_epilogue;
	}

	ostrm << rstr << std::endl;
	std::cout << "Report written to " << rfn << std::endl;
}

int main(int argc,char *argv[])
{
	rsb_err_t retval = RSB_ERR_NO_ERROR;
	rsbt_bopts bopts;
	rsbt_topts topts;
	rsbt_bress bress; 
	std::ostringstream ssos; // setup string

	std::string ops{""}; // operations string
	bool want_comp=true;
	bool want_rand=true;
	bool want_recf=false;
	bool want_band=true;
	bool want_rect=true;
	bool want_verb=false;
	rsb_nnz_idx_t inc = 10, mul = 3, minn = 1, maxn = 1000000;
	std::vector<rsb_flags_t> symfa = {};
	
#ifdef HAVE_SYS_UTSNAME_H 
	std::cout << "running on " << get_utsname() << std::endl;
#endif /*HAVE_SYS_UTSNAME_H */ 
#if HAVE_GETOPT_H
	for(auto coi=0;coi<options_count;++coi)
	{
		options[coi].name=rsbt_options[coi].name;
		options[coi].has_arg=rsbt_options[coi].has_arg;
		options[coi].flag=rsbt_options[coi].flag;
		options[coi].val=rsbt_options[coi].val;
	}

    	for (;;) {
	int opt_index = 0;
	const int c = getopt_long(argc, argv, "hn:qst:Vv+"
				    ,
				options, &opt_index);
		if (c == -1)
		    break;
		switch (c) {
			case '+': // --cplusplus
				bopts.want_rpp = true;
			break;
			case 't': // --types
				if( optarg && ( *optarg==':' || std::string(optarg) == "all" ) )
					bopts.ts=RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS;
				else
				{
					if( optarg && std::string(optarg) == "blas" )
						bopts.ts=RSB_BLAS_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS ;
					else
						bopts.ts=optarg;
				}
				if(optarg && *optarg=='?')
				{
					std::cout << "librsb configured to support types " << RSB_XSTRINGIFY(RSB_M4_MATRIX_TYPES_STRING) << " -- " << RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS << std::endl;
					++optarg;
				}
			break;
			case 'n': // --nthreads  // TODO: shall affect ALL implementations
				bopts.tna = vals_array<rsb_int_t>(optarg);
			break;
			case 0x73657465: // --self-test
				std::exit(self_test());
			break;
			case 0x736c736d: // --skip-loading-symmetric-matrices
				RSB_DO_FLAG_ADD(bopts.dro.sf,RSB_FLAG_SYMMETRIC);
			break;
			case 0x736c686d: // --skip-loading-hermitian-matrices
				RSB_DO_FLAG_ADD(bopts.dro.sf,RSB_FLAG_HERMITIAN);
			break;
			case 0x736c6e75: // --skip-loading-not-unsymmetric-matrices
				RSB_DO_FLAG_ADD(bopts.dro.sf,RSB_FLAG_HERMITIAN);
				RSB_DO_FLAG_ADD(bopts.dro.sf,RSB_FLAG_SYMMETRIC);
			break;
			case 0x736c756d: // --skip-loading-unsymmetric-matrices
				bopts.dro.af = RSB_FLAG_HERMITIAN | RSB_FLAG_SYMMETRIC;
			break;
			case 0x73796d6d: // --symmetric
				symfa.push_back(RSB_FLAG_UPPER_SYMMETRIC);
				symfa.push_back(RSB_FLAG_LOWER_SYMMETRIC);
			break;
			case 0x6865726d: // --hermitian
				symfa.push_back(RSB_FLAG_UPPER_HERMITIAN);
				symfa.push_back(RSB_FLAG_LOWER_HERMITIAN);
			break;
			case 0x756e7379: // --unsymmetric
				symfa.push_back(RSB_FLAG_NOFLAGS);
			break;
			case 0x6e726873: // --nrhs
				bopts.nrhsa = vals_array<rsb_coo_idx_t>(optarg);
			break;
			case 0x696e6378: // --incx
				bopts.incxa = vals_array<rsb_coo_idx_t>(optarg);
			break;
			case 0x696e6379: // --incy
				bopts.incya = vals_array<rsb_coo_idx_t>(optarg);
			break;
			case 0x646e7374: // --density
				bopts.dnst_pct_a = vals_array<int>(optarg);
			break;
			case 0x616c6661: // --alpha
				bopts.alphaia = vals_array<int>(optarg);
			break;
			case 0x62657461: // --beta
				bopts.betaia = vals_array<int>(optarg);
			break;
			case 0x636f7366: // --coo-stretch-factor
				bopts.csf_fa = vals_array<rsb_coo_idx_t>(optarg);
			break;
			case 0x62617368: // --band-stretch-factor
				bopts.bsf_fa = vals_array<rsb_coo_idx_t>(optarg);
			break;
			case 's': // --serial
				ops+="s";
			break;
			case 0x70726e74: // --no-timings
				want_no_timings=true;
				std::cout << "won't print timings" << std::endl;
			break;
			case 0x77720a: // --rand
				want_rand=true;
				std::cout << "enabling randomization" << std::endl;
			break;
			case 0x6e77720a: // --no-rand
				want_rand=false;
				std::cout << "disabling randomization" << std::endl;
			break;
			case 0x72616d61: // --max-random-val
				// TODO: decide whether to merge with --rand
				std::istringstream(optarg) >> bopts.mor;
				std::cout << "if enabled, randomize values in [1," << bopts.mor << std::endl;
				assert(bopts.mor>=1);
			break;
			case 'h': // --help
			std::cout << argv[0] << " [OPTIONS] [FILENAMES/DIRS] " << std::endl;
			std::cout << "# with OPTIONS from:" << std::endl;
			for (auto o : rsbt_options)
			if( o.name )
			{
				std::cout << " --" << o.name;
				if(isascii(o.val))
					std::cout << " / -" << static_cast<char>( o.val );
				if(o.has_arg == required_argument)
					std::cout << " {arg}";
				//if(o.has_arg == optional_argument)
				//	std::cout << " [arg]";
				if(o.desc && *o.desc )
					std::cout << " // " << o.desc;
				std::cout << std::endl;
			}
			std::exit(0);
			break;
			case 'V': // --version
				std::cout << library_version_string();
				std::exit(0);
			break;
			case 'q': // --quiet
			want_quiet=true; // TODO: want_quiet vs want_verb ?
			break;
			case 'v': // --verbose
			want_verb=true;
			break;
			case 0x696e630a: // --inc
				std::istringstream(optarg) >> inc;
			break;
			case 0x6d756c0a: // --mul
				std::istringstream(optarg) >> mul;
			break;
			case 0x6d61780a : // --max
				std::istringstream(optarg) >> maxn;
			break;
			case 0x6d696e0a: // --min
				std::istringstream(optarg) >> minn;
			break;
			case 0x6d6b6c0a: // --mkl
#ifdef RSBT_USE_MKL_SPBLAS
			bopts.want_mkl=true;
			std::cout << "Will perform MKL tests, too."<< std::endl;
#else /* RSBT_USE_MKL_SPBLAS */
			std::cout << "Built without the MKL."<< std::endl;
			return re2ce(RSB_ERR_GENERIC_ERROR);
#endif /* RSBT_USE_MKL_SPBLAS */
			break;
			case 0x7265636f: // --recover-failed
				std::cout << "Will recover saved failed case (EXPERIMENTAL,UNFINISHED)." << std::endl;
				want_recf=true;
				want_rand=false;
				want_dump_failed=false;
			break;
			case 0x6e726374: // --no-rectangular
				std::cout << "Will NOT test rectangular matrices." << std::endl;
				want_rect=false;
			break;
			case 0x72656374: // --rectangular
				std::cout << "Will test rectangular matrices." << std::endl;
				want_rect=true;
				want_band=false;
			break;
			case 0x72656e64: // --render-only
#ifdef RSBT_WANT_RSB_HPP
				bopts.ropts=rsbt_bopts::RenderOpts::Once; 
#else
				std::cout << "Built without render support!"<< std::endl;
				return re2ce(RSB_ERR_GENERIC_ERROR);
#endif
			break;
			case 0x72656e74: // --render
#ifdef RSBT_WANT_RSB_HPP
				bopts.ropts=rsbt_bopts::RenderOpts::All; 
#else
				std::cout << "Built without render support!"<< std::endl;
				return re2ce(RSB_ERR_GENERIC_ERROR);
#endif
			break;
			case 0x6e626e72:/* --no-basename-render */
				bopts.want_basenameplot=false;
				std::cout << "Any rending will be produced in current directory"<< std::endl;
			break;
			case 0x6d617854: // --maxtimes
				suffixed_val_read(optarg,bopts.dro.maxtimes);
				std::cout << "Benchmark will iterate at most " << bopts.dro.maxtimes << " times" << std::endl;
			break;
			case 0x6d696e54: // --mintimes
				suffixed_val_read(optarg,bopts.dro.mintimes);
				std::cout << "Benchmark will iterate at least " << bopts.dro.mintimes << " times" << std::endl;
			break;
			case 0x6d696e74: // --min_t
				suffixed_val_read(optarg,bopts.dro.min_t);
				std::cout << "Benchmark will sample for at least " << bopts.dro.min_t << " s" << std::endl;
			break;
			case 0x746d6178: // --max_t
				suffixed_val_read(optarg,bopts.dro.max_t);
				std::cout << "Benchmark will sample for at most " << bopts.dro.max_t << " s" << std::endl;
			break;
			case 0x6d617474: // --max-test-time
				suffixed_val_read(optarg,max_tt);
				std::cout << "Test will run for at most " << max_tt << " s" << std::endl;
			break;
			case 0x6f74636e: // --only-test-case-n
				suffixed_val_read(optarg,otcn);
				std::cout << "Only test case " << otcn << " will be executed." << std::endl;
			break;
			case 0x736e7463: // --skip-n-test-cases
				suffixed_val_read(optarg,sntc);
				std::cout << "First " << sntc << " test cases will be skipped." << std::endl;
			break; 
			case 0x73656e74: // --skip-except-every-n-test-cases
				suffixed_val_read(optarg,sent);
				std::cout << "Execute once every " << sent << " test cases." << std::endl;
			break; 
			case 0x7365726e: // --skip-except-every-random-n-test-cases
				suffixed_val_read(optarg,sern);
				std::cout << "Execute once every " << sern << " randomized test cases." << std::endl;
			break; 
			case 0x6f6e7463: // --only-n-test-cases
				suffixed_val_read(optarg,ontc);
				std::cout << "Only " << ontc << " test cases will be executed." << std::endl;
			break; 
			case 0x736c6d6: // --skip-loading-if-more-nnz-matrices
				suffixed_val_read(optarg,bopts.dro.max_nnz);
				std::cout << "Benchmark will skip loading matrices with more than " << bopts.dro.max_nnz << " nnz" << std::endl;
			break;
			case 0x736c6e6e: // --skip-loading-if-less-nnz-matrices
				suffixed_val_read(optarg,bopts.dro.min_nnz);
				std::cout << "Benchmark will skip loading matrices with less than " << bopts.dro.min_nnz << " nnz" << std::endl;
			break;
			//default: 1;
			case 0x6576690a: // --extra-verbose-interface
				bopts.dro.evi++;
			break;
	/*		case 0x6c6c6d0a: // --leaf-level-multivec
				bopts.dro.llm=0; // this (l.l.m.) is default anyway
			break;*/
			case 0x6e6c6d6d: // --no-leaf-multivec
				bopts.dro.llm=1;
				std::cout << "Will pass " << bopts.dro.llm << " to RSB_IO_WANT_LEAF_LEVEL_MULTIVEC" << std::endl;
			break;
			case 0x76657475: // --verbose-tuning
				bopts.dro.wvt=1;
			break;
			case 0x72707274: // --report
				bopts.rfn=optarg;
				std::cout << "Before termination will dump a benchmark report to " << bopts.rfn << std::endl;
				if (has_extension(bopts.rfn,".tex"))
				{
					std::cout << "LaTeX dump with renderings." << std::endl;
#ifdef RSBT_WANT_RSB_HPP
					bopts.ropts=rsbt_bopts::RenderOpts::All; 
#endif
				}
			break;
			case 0x7273620a : // --rsb-tester
				RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_DEFAULT_MATRIX_FLAGS);
			break;
			case 0x636f6f0a : // --coo-tester
				RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS);
				//RSB_FLAG_WANT_COO_STORAGE
			break;
			case 0x68636f6f : // --hcoo-tester
				RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_USE_HALFWORD_INDICES_COO);
			break;
			case 0x6373720a : // --csr-tester
				RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS);
				//RSB_FLAG_WANT_BCSS_STORAGE
			break;
			case 0x68637372 : // --hcsr-tester
				//RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS|RSB_FLAG_USE_HALFWORD_INDICES);
				RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_USE_HALFWORD_INDICES_CSR);
			break;
			case  0x68616c66 : // --half-indices-tester
				RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_USE_HALFWORD_INDICES);
			break;
			case  0x7265630a: // --recursive-tester
				RSB_DO_FLAG_ADD(topts.fflags,RSB_FLAG_QUAD_PARTITIONING);
			break;
			case  0x696d6469: // --diagonal-implicit
				RSB_DO_FLAG_ADD(bopts.eflags,RSB_FLAG_UNIT_DIAG_IMPLICIT);
			break;
			case  0x666f6969: // --fortran-indices
				std::cout << "Use Fortran indices interface." << std::endl;
				RSB_DO_FLAG_ADD(bopts.eflags,RSB_FLAG_FORTRAN_INDICES_INTERFACE);
			break;
			case  0x746f6d69: // --tolerate-mismatch
				std::cout << "Will tolerate mismatch errors." << std::endl;
				bopts.want_tmm = true;
			break;
			case  0x61686f6b: // --any-hermitian-ok
				std::cout << "Will not skip hermitian non-complex." << std::endl;
				bopts.want_any_herm = true; // any type as hermitian is ok
			break;
			case 0x6e6f7472: // --no-trans
				std::cout << "Will use no transposition." << std::endl;
				bopts.transAa = {RSB_TRANSPOSITION_N};
			break;
			case 0x74756e65: // --no-tune
				std::cout << "Will not invoke autotuning routine." << std::endl;
				bopts.want_autotune = false;
			break;
			case 0x6d617872: // --tune-maxr
				{ std::istringstream iss(optarg); iss >> bopts.maxr; }
				std::cout << "Setting maxr autotuning parameter to " << bopts.maxr << std::endl;
			break;
			case 0x6d617874: // --tune-maxt
				{ std::istringstream iss(optarg); iss >> bopts.maxt; }
				std::cout << "Setting maxt autotuning parameter to " << bopts.maxt << std::endl;
			break;
			case 0x7472616e: // --transA
			{
				const auto btransAa = bopts.transAa; // backup
				if(optarg)
					bopts.transAa = {};
				std::cout << "Passed transA string " << optarg  << " using among possible ones: (" << RSB_TRANSPOSITIONS_PREPROCESSOR_SYMBOLS << ") so:";
				for(auto c: std::string{optarg})
					if(RSB_CHAR_AS_TRANSPOSITION(c)!='?') // RSB_INVALID_TRANS
					{
						bopts.transAa.push_back(RSB_CHAR_AS_TRANSPOSITION(c));
						std::cout << " " << char(std::toupper(c)) << "/" << std::toupper(c);
					}
				std::cout << std::endl;
				if( bopts.transAa.empty() )
				{
					std::cout << "String " << std::string{optarg} << " had no valid transposition char!" << std::endl;
					std::exit(-1);
				}
			}
			break;
			case '?':
				std::cout << "unrecognized option, aborting." << std::endl;
				return re2ce(RSB_ERR_GENERIC_ERROR);
			break;
		} /* end options parsing */
	}
#endif /* HAVE_GETOPT_H */
	if(!symfa.empty())
		bopts.symfa = symfa;
	if ( topts.fflags == RSB_FLAG_NOFLAGS )
		topts.fflags = RSB_FLAG_DEFAULT_MATRIX_FLAGS;
	topts.fflags |= bopts.eflags;

	if ( bopts.ts.empty() )
	{
		bopts.ts = RSB_MATRIX_SPBLAS_TYPE_CODES_ARRAY; // default type string; vs. all types in RSB_MATRIX_TYPE_CODES_ARRAY
		std::cout << "no types string specified; ";
	}
	std::transform(bopts.ts.begin(), bopts.ts.end(), bopts.ts.begin(),
		[ ](char c) -> char { return std::toupper(c); });
	bopts.ts.erase(std::remove_if(bopts.ts.begin(), bopts.ts.end(), [&](auto x){return RSB_MATRIX_UNSUPPORTED_TYPE(x);}), bopts.ts.end());
#ifdef RSB_NUMERICAL_TYPE_INT
	if ( !bopts.ts.empty() )
		if(bopts.ts[0] == 'I')
		{
			std::cout << "found you have I as first type: rotating types to avoid this" << std::endl;
			std::rotate( bopts.ts.begin(), bopts.ts.begin()+1, bopts.ts.end() );
		}
#endif
	std::cout << "using types in " << bopts.ts << std::endl;

	if(want_rand)
		ops+="r";
	if(want_comp)
		ops+="=";
	if(want_recf)
		ops+="f";

#if 0
	if(want_recf)
	if(optind < argc)
	{
		//for(const auto & fn : bopts.fna)
		for (int i = optind; i < argc; i++)
		{
			if(i-optind==0)
				std::cout << "Got matrix: " << argv[i] << std::endl;
			if(i-optind==1)
				std::cout << "Got rhs: " << argv[i] << std::endl;
		}
		// TODO: unfinished here.
		// ...
	}
#endif


	if(optind < argc)
	{
		for (int i = optind; i < argc; i++)
		{
#if HAVE_FILESYSTEM && USE_CXX17
			namespace fs = std::filesystem;
			if(fs::is_directory(argv[i]))
			{
				std::string s;
				std::cout << "dir: " << argv[i] << std::endl;
				for(auto& p: fs::recursive_directory_iterator(argv[i]))
				if(fs::is_regular_file(p))
				if (is_mtx_market_filename( s = p.path().string() ))
				if (rsbt_passes_filter(s, bopts.dro))
					bopts.fna.push_back(s);
			}
			else
#endif /* HAVE_FILESYSTEM && USE_CXX17 */
				if (rsbt_passes_filter(argv[i], bopts.dro))
					bopts.fna.push_back(argv[i]);
		}
		
		if(!bopts.rfn.empty())
		{
			ssos	<< "Experiment setup using: " << std::endl
				<< " types in: " << bopts.ts << std::endl
				<< " files in: " << bopts.fna << std::endl
				<< " thread spec in: " << bopts.tna << std::endl
				<< " transpositions in: " << bopts.transAa << std::endl
				<< " nrhs in: " << bopts.nrhsa << std::endl
				<< " incx in: " << bopts.incxa << std::endl
				<< " incy in: " << bopts.incya << std::endl
				<< " density in: " << bopts.dnst_pct_a << std::endl
				<< " alpha in: " << bopts.alphaia << std::endl
				<< " beta in: " << bopts.betaia << std::endl
				<< " tune-maxr in: " << bopts.maxr << std::endl
				<< " tune-maxt in: " << bopts.maxt << " s" << std::endl
				<< " bench-mintimes in: " << bopts.dro.mintimes << std::endl
				<< " bench-maxtimes in: " << bopts.dro.maxtimes << std::endl
				<< " bench-min_t in: " << bopts.dro.min_t << " s" << std::endl
				<< " bench-max_t in: " << bopts.dro.max_t << " s" << std::endl
			;
			bress.lrpv.push_back(rsbt_bress::lrpi_t{ "Experiment Setup","","", str_to_tex( std::regex_replace(ssos.str(), std::regex("\n"), "$&\\\\") ) });
			std::cout << ssos.str() << std::endl;
		}
		if( bopts.fna.empty() )
		{
			std::cout << "no matrix file loaded !" << std::endl;
			retval = RSB_ERR_GENERIC_ERROR;
		}

#ifdef RSBT_WANT_RSB_HPP
		hpp_tests(bopts,bress,ops);
#else /* RSBT_WANT_RSB_HPP */
		std::cout << "C++ interface to librsb not found: terminating!" << std::endl;
#endif /* RSBT_WANT_RSB_HPP */
		eps_report(bopts.rfn, bress);
		std::cout << "All done." << std::endl;
		return re2ce(retval);
	}
	else
		std::cout << "no matrix files specified. will perform a general test." << std::endl;

	if(ops.find('s')!=std::string::npos)
	{
		// NOTE: this ignores --types
#ifdef RSB_NUMERICAL_TYPE_DOUBLE
		retval = serial_test<double>(want_verb);
#endif
#ifdef RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
		retval = serial_test<std::complex<double>>(want_verb);
#endif
		return re2ce(retval);
	}

	std::cout << "generating matrices sized " << minn << " + " << inc << "  * " << mul << " ... " << maxn << "." << std::endl;
	for(topts.n=minn; topts.n<=maxn; topts.n=topts.n*mul+inc)
	{
		if(want_band)
			retval = band_testers(bopts,ops,topts);
		if(retval != RSB_ERR_NO_ERROR)
			return re2ce(retval);

		if(want_rect)
			retval = rectangular_testers(bopts,ops,topts);
		if(retval != RSB_ERR_NO_ERROR)
			return re2ce(retval);
	}
	std::cout << cnt-scnt << " comparison / correctness tests executed successfully";
	if(scnt)
		std::cout << " (" << scnt << " cases skipped)";
	std::cout << "." << std::endl;

	return re2ce(retval);
} /* main */
