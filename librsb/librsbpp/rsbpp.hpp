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
#ifndef RSBPP_HPP_INCLUDED
#define RSBPP_HPP_INCLUDED

#include <algorithm>
#include <vector>
#include <cstdlib>	// malloc, free
#include <sstream>
#include <cmath>
#include <cfloat>
#include <typeinfo>
#include <map>
#include <set>
#include <stdexcept>	// TODO: throw std::runtime_error{"error msg.."}; .. try {..} catch (..) {}
#include <iostream>
#include <numeric> // std::accumulate
#include <chrono>
#include <string>
#include <fstream>
#include <cassert>
#include <memory>
#include <limits>
#include <cstdint>	// int64_t
#include <complex>	// std::complex<...>
#include <ccomplex>	// complex aka _Complex
#if defined(RSBPP_HAVE_STD_THREAD) && defined(RSBPP_HAVE_STD_MUTEX )
#define RSBPP_WANT_CPP_THREADS 1
#include <thread>
#include <mutex>
#else
/* no C++ threading */
#endif
#include <valarray>
#include <array>
#include <type_traits>
#include <ctype.h>	// toupper
#if defined(_OPENMP)
#include <omp.h>	// 
#else /* _OPENMP */
#define omp_get_max_threads() 1
#define omp_get_thread_num() 1
#endif /* _OPENMP */
#include <cstddef>    // ptrdiff_t
#include "openmp_allocator.hpp"	// OpenMP_Allocator
#if __cplusplus >= 201709
#define RSBPP_WANT_CPP20 1
#include <span>
#endif
#ifdef RSBPP_HAS_RSB_H
#include <rsb.h>	// flags etc; TODO: what about constants only ?
#else /* RSBPP_HAS_RSB_H */
/* FIXME: temporary from rsb.h: */
/*! \ingroup rsb_doc_misc rsb_doc_rsb
\brief Flags for getting matrix information via #rsb_mtx_get_info()/#rsb_mtx_get_info_str().
*/
enum rsb_mif_t
{
  RSB_MIF_INDEX_STORAGE_IN_BYTES__TO__SIZE_T		=0x00000001	/*!< Index storage occupation, in bytes. (size_t) */
, RSB_MIF_INDEX_STORAGE_IN_BYTES_PER_NNZ__TO__RSB_REAL_T	=0x00000002	/*!< Index storage occupation per nnz, in bytes. (#rsb_real_t) */
, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T		=0x00000004	/*!< Rows count(#rsb_coo_idx_t) */
, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T		=0x00000008	/*!< Columns count (#rsb_coo_idx_t) */
, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T		=0x00000010	/*!< Nonzeroes count (#rsb_nnz_idx_t) */
, RSB_MIF_TOTAL_SIZE__TO__SIZE_T			=0x00000020	/*!< Total size, in bytes (size_t) */
, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T			=0x00000040	/*!< Matrix flags (#rsb_flags_t) */
, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T		=0x00000080	/*!< Matrix type code (#rsb_type_t) */
, RSB_MIF_MATRIX_INFO__TO__CHAR_P			=0x00000100	/*!< Matrix info string, only for #rsb_mtx_get_info_str() (#rsb_char_t*) */
, RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T		=0x00000200	/*!< Leaf submatrices count (#rsb_blk_idx_t) */
};typedef signed int rsb_flags_t;
typedef rsb_flags_t rsb_trans_t;
typedef char rsb_type_t;
#define RSB_BOOL_FALSE	0 /*!< A "false" value for #rsb_bool_t. */
#define RSB_BOOL_TRUE	1	/*!< A "true"  value for #rsb_bool_t. */
#define RSB_DO_FLAG_HAS(V,F)	((((V)&(F))==(F))?RSB_BOOL_TRUE:RSB_BOOL_FALSE)	 /*!< Presence check for flag \c F. */
#define  RSB_TRANSPOSITION_N 0x4E /*!< N: Non transposed flag, valid for \ref rsb_trans_t typed variables. */
#define  RSB_TRANSPOSITION_T 0x54 /*!< T: Transposed flag value, valid for \ref rsb_trans_t valued variables. */
#define  RSB_TRANSPOSITION_C 0x43 /*!< C: Conjugated transpose flag, valid for \ref rsb_trans_t typed variables. */
/*! Default storage flags. */
#define RSB_FLAG_DEFAULT_STORAGE_FLAGS		 	(RSB_FLAG_WANT_BCSS_STORAGE|RSB_FLAG_WANT_COO_STORAGE)

/*! A flag combination specifying a pure COO matrix.  */
#define RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS		 	RSB_FLAG_WANT_COO_STORAGE 

/*! A flag combination specifying a pure CSR matrix.  */
#define RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS		 	RSB_FLAG_WANT_BCSS_STORAGE

/*! A flag combination specifying a pure RSB matrix.  */
#define RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS (RSB_FLAG_QUAD_PARTITIONING|RSB_FLAG_USE_HALFWORD_INDICES|RSB_FLAG_WANT_COO_STORAGE|RSB_FLAG_WANT_BCSS_STORAGE)

/*! A flag combination specifying a matrix in a default, supported format.  */
#define RSB_FLAG_DEFAULT_MATRIX_FLAGS			RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS

/*! The null (empty) flag. */
#define RSB_FLAG_NOFLAGS		 		0x000000

/*! The identical flag (used in cloning function #rsb_mtx_clone). */
#define RSB_FLAG_IDENTICAL_FLAGS RSB_FLAG_NOFLAGS

/*! If set, the input/output coordinate indices will be assumed to be 1 based. */
#define RSB_FLAG_FORTRAN_INDICES_INTERFACE		0x000001

/*! If set, the input/output coordinate indices will be assumed to be 0 based (default). */
#define RSB_FLAG_C_INDICES_INTERFACE		0x000000

/*! If set, the matrix will internally use a half word (16 bit) type for indices. */
#define RSB_FLAG_USE_HALFWORD_INDICES	0x000002

/*! Used to specify multi-vector (dense matrix) operations. */
#define RSB_FLAG_WANT_ROW_MAJOR_ORDER 			0x000000

/*! Used to specify multi-vector (dense matrix) operations. */
#define RSB_FLAG_WANT_COLUMN_MAJOR_ORDER 		0x4000000

/*! If set, the code will assume the input nonzeroes as sorted.	*/
#define RSB_FLAG_SORTED_INPUT				0x000004

/*! If set, the matrix is considered as triangular. \see #RSB_FLAG_LOWER,#RSB_FLAG_UPPER. */ 
#define RSB_FLAG_TRIANGULAR 				0x000008

/*! If set, the matrix will be stored in as lower (triangular or symmetric). \see #RSB_FLAG_TRIANGULAR,#RSB_FLAG_SYMMETRIC,#RSB_FLAG_UPPER. */
#define RSB_FLAG_LOWER		 			0x000010

/*! If set, the matrix will be stored in as upper (triangular or symmetric). \see #RSB_FLAG_LOWER*/
#define RSB_FLAG_UPPER		 			0x000020

/*! If set, the (whole super-)matrix will not store the diagonal, which will be assumed to be unitary. */
#define RSB_FLAG_UNIT_DIAG_IMPLICIT			0x000040

/* ghost flag ( moved in a non-public header, and reserved): 0x000080	*/
/* ghost flag ( moved in a non-public header, and reserved): 0x80000000	*/

/*! If set, the matrix will use COO storage, where necessary. */
#define RSB_FLAG_WANT_COO_STORAGE		0x000100

/*! Keep the last nonzero duplicate, at matrix assembly time. */ 
#define RSB_FLAG_DUPLICATES_KEEP_LAST				0x000000

/*! The default nonzeroes duplicates handling.  */ 
#define RSB_FLAG_DUPLICATES_DEFAULT_HANDLE			0x000000

/*! Compute and keep the sum of nonzero duplicates, at matrix assembly time.  */ 
#define RSB_FLAG_DUPLICATES_SUM				0x000200

/*! If set, explicit zeros will not be inserted	\warning: this flag is active by default	*/
#define RSB_FLAG_DISCARD_ZEROS				0x000400

/* ghost flag ( moved in a non-public header, and reserved): 0x000800 */
/* ghost flag ( moved in a non-public header, and reserved): 0x001000 */

/*! If set, matrix will be organized as a quad tree of submatrices. */
#define RSB_FLAG_QUAD_PARTITIONING 			0x002000

/*! If set, the block partitioning will be fixed (BCSS: BCSR or BCSC, but no VBR).	*/
#define RSB_FLAG_WANT_BCSS_STORAGE 			0x004000

/* ghost flag ( moved in a non-public header, and reserved): 0x008000 */
/* ghost flag ( moved in a non-public header, and reserved): 0x010000 */
/* ghost flag ( moved in a non-public header, and reserved): 0x020000 */

/*! If set, matrices will be fit in the three input coo arrays, after conversion. */ 
#define RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS		0x040000

/*! \internal \todo: should remove this. */ 
#define RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT	0x080000

/* ghost flag ( moved in a non-public header, and reserved): 0x100000*/
/* ghost flag (temporarily reserved): 0x200000*/

/*! If set, the input matrix will be treated as symmetric (stored as a lower triangular one by default). \see #RSB_FLAG_LOWER,#RSB_FLAG_LOWER. */
#define RSB_FLAG_SYMMETRIC 			0x400000

/*! If set, the input matrix will be treated as symmetric hermitian (stored as a lower triangular one). \see #RSB_FLAG_LOWER,#RSB_FLAG_LOWER. */
#define RSB_FLAG_HERMITIAN 			0x800000

/*! If set, recursion on small matrices will last at least the number of active threads. */
#define RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS	0x1000000

/* ghost flag ( moved in a non-public header, and reserved): 0x2000000	*/

/*! Combined flags for a lower hermitian matrix. */
#define RSB_FLAG_LOWER_HERMITIAN			(RSB_FLAG_HERMITIAN | RSB_FLAG_LOWER)

/*! Combined flags for an upper hermitian matrix. */
#define RSB_FLAG_UPPER_HERMITIAN			(RSB_FLAG_HERMITIAN | RSB_FLAG_UPPER)

/*! Combined flags for a lower triangular matrix. */
#define RSB_FLAG_LOWER_TRIANGULAR 			(RSB_FLAG_TRIANGULAR | RSB_FLAG_LOWER)

/*! Combined flags for an upper triangular matrix. */
#define RSB_FLAG_UPPER_TRIANGULAR 			(RSB_FLAG_TRIANGULAR | RSB_FLAG_UPPER)

/*! Combined flags for a symmetric, lower-stored matrix. */

#define RSB_FLAG_LOWER_SYMMETRIC 			(RSB_FLAG_SYMMETRIC | RSB_FLAG_LOWER)

/*! Combined flags for a diagonal matrix. */
#define RSB_FLAG_DIAGONAL 				(RSB_FLAG_UPPER | RSB_FLAG_LOWER)

/*! Combined flags for a symmetric, upper-stored matrix. */
#define RSB_FLAG_UPPER_SYMMETRIC 			(RSB_FLAG_SYMMETRIC | RSB_FLAG_UPPER)

/*! If set, the matrix will be subdivided at a finer grain on diagonal blocks. */
#define RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG 	0x8000000

/*! If set, the input COO arrays to the assembly functions will not be freed at matrix destruction time.
  \warning Please do NOT use this flag, for the default memory allocation handling is still not specified. Instead, use the in place allocation functions: #rsb_mtx_alloc_from_csr_inplace() and #rsb_mtx_alloc_from_coo_inplace().
 */
#define RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS 		0x40000000

/* Reserved, undocumented flags. Not for use. */
#define RSB_FLAG_USE_CSR_RESERVED	0x200000

/*! \internal Combined flags for half word CSR. */
#define RSB_FLAG_USE_HALFWORD_INDICES_CSR	(RSB_FLAG_USE_HALFWORD_INDICES|RSB_FLAG_USE_CSR_RESERVED)

/*! Combined flags for half word COO. */
#define RSB_FLAG_USE_HALFWORD_INDICES_COO	(RSB_FLAG_USE_HALFWORD_INDICES|RSB_FLAG_WANT_COO_STORAGE)

/*! A combination of flags which is forbidden (so don't use it). */
#define RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES	(RSB_FLAG_USE_HALFWORD_INDICES_COO|RSB_FLAG_USE_HALFWORD_INDICES_CSR)
/*!@}*/

typedef signed int rsb_err_t;
typedef signed int rsb_nnz_idx_t;
typedef signed int rsb_coo_idx_t;
#define RSB_ERR_CAST(E) (-(E))
#define RSB_ERR_NO_ERROR		RSB_ERR_CAST(0x000)
#define RSB_ERR_UNSUPPORTED_OPERATION	RSB_ERR_CAST(0x002)
#define RSB_ERR_UNSUPPORTED_TYPE	RSB_ERR_CAST(0x004)
#if 0
#define RSB_ERR_UNSUPPORTED_FORMAT	RSB_ERR_CAST(0x008)
#endif
#define RSB_ERR_INTERNAL_ERROR		RSB_ERR_CAST(0x010)
#define RSB_ERR_BADARGS			RSB_ERR_CAST(0x020)
#if 0
#define RSB_ERR_ENOMEM			RSB_ERR_CAST(0x040)
#endif
#define RSB_ERR_UNIMPLEMENTED_YET	RSB_ERR_CAST(0x100)
#if 0
#define RSB_ERR_LIMITS			RSB_ERR_CAST(0x200)
#define RSB_ERR_FORTRAN_ERROR		RSB_ERR_GENERIC_ERROR
#define RSB_ERR_UNSUPPORTED_FEATURE	RSB_ERR_CAST(0x400)
#define RSB_ERR_NO_USER_CONFIGURATION	RSB_ERR_CAST(0x800)
#define RSB_ERR_CORRUPT_INPUT_DATA	RSB_ERR_CAST(0x1000)
#define RSB_ERR_FAILED_MEMHIER_DETECTION	RSB_ERR_CAST(0x2000)
#define RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS	RSB_ERR_CAST(0x4000)
#define RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT	RSB_ERR_CAST(0x8000)
#define RSB_ERR_INVALID_NUMERICAL_DATA	RSB_ERR_CAST(0x10000)
#define RSB_ERR_MEMORY_LEAK	RSB_ERR_CAST(0x20000)
#endif

/*! Program success error code (int). */
#ifdef EXIT_SUCCESS
#define RSB_PROGRAM_SUCCESS	(EXIT_SUCCESS)
#else
#define RSB_PROGRAM_SUCCESS		(0)
#endif


/* Note these are not part of rsb.h : */
#define  RSB_NUMERICAL_TYPE_DOUBLE  'D' /*!< Character code for type double. */
#define  RSB_NUMERICAL_TYPE_FLOAT  'S' /*!< Character code for type float. */
#define  RSB_NUMERICAL_TYPE_FLOAT_COMPLEX  'C' /*!< Character code for type float complex. */
#define  RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX  'Z' /*!< Character code for type double complex. */
#endif /* RSBPP_HAS_RSB_H */
typedef unsigned short int rsb_half_idx_t; /* this is not part of rsb.h */



#define RSBPP_ERR_CAST(E) (-(E))
#define RSBPP_ERR_UNSUPPORTED_OPERATION	RSBPP_ERR_CAST(0x002)
#ifndef RSBPP_PRESET_SPMM_BLOCKINGS
#define RSBPP_PRESET_SPMM_BLOCKINGS 4 /* 1:no blocking; 4: otherwise */
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */

#define STRINGIFY(X) #X
#define STRINGIFY1(X) #X << " = " << X 
#define HALFSQUARE(N) ((N*(N+1))/2)
template <typename X> constexpr X tone = (X) 1;
template <typename X> constexpr X tzero= (X) 0.0;
template <typename X> inline X VALFROMIJ(int & I, int & J) {  return ((I+1)*100+(J+1)); }
#define RSB_MAX_THREADS 64
#define RSB_MIN_THREADS 1
#define RSB_ERRM_FIFL "Failed in " << __FILE__ << ":" << __LINE__
#if defined(__cplusplus) && (__cplusplus>= 201703L) 
#define USE_CXX17 1
#else
#define USE_CXX17 0
#endif

// #define NDEBUG 1 // uncomment this to deactivate assert() statements

//#define VALFROMIJ(I,J) (I+1)*1000+(J+1)

#if 0
// FIXME: UNFINISHED/BROKEN
template <typename IT=int>
class CooAllocator final/*:public std::allocator <IT>*/
{
	public:
	using value_type = IT;
	using pointer = IT*;
	using const_pointer = const IT*;
	using reference = IT& ;
	// using const_reference = const IT&;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t ;
	pointer allocate(const size_t n, const void* = 0) {return static_cast <pointer> ( std::malloc(sizeof(value_type)*n));}
	void deallocate(pointer p, size_type n) { std::free(p); }
	void construct(pointer p, const IT& v) { /* no initialization */ }
	void destroy(pointer p) { /* no destruction */ }
};

template <class T1, class T2>
bool operator== (const CooAllocator<T1>&, const CooAllocator<T2>&) throw() { return true; }
template <class T1, class T2>
bool operator!= (const CooAllocator<T1>&, const CooAllocator<T2>&) throw() { return false; }
#endif
#include <librsbpp.h>

#ifndef RSBPP_WANT_MINIMAL_LIBRSBPP
#define RSBPP_WANT_MINIMAL_LIBRSBPP 0
#endif /* RSBPP_WANT_MINIMAL_LIBRSBPP */
#define RSBPP_WANT_ALL !RSBPP_WANT_MINIMAL_LIBRSBPP

template<typename T>
static
typename std::enable_if< std::is_scalar<T>::value,T>::type rsbpp_conj(T v) {return v;}

template<typename T>
static
typename std::enable_if<!std::is_scalar<T>::value,T>::type rsbpp_conj(T v) {return std::conj(v);}

#if RSBPP_WANT_ALL
static int g_verbosity = 0;
static size_t g_cbs = (1024*1024/4); // default, one fourth of a MiB
static std::vector<char> g_types_a {'d'};
static std::vector<int> g_nrhs_a {1,2,3,4,5,6,7,8};
static std::vector<rsb_trans_t> g_trans_a {RSB_TRANSPOSITION_N, RSB_TRANSPOSITION_T, RSB_TRANSPOSITION_C};
using mytime_t = double ; // TODO: move

/* class BitVector begin */
using BitVector = std::vector<bool>;
/* class BitVector end */

/* class QuadrantT begin */
template <typename IT=int>
class QuadrantT final
{
	enum {Five=5, Four=4};
	std::array<IT,Five> off_; /* offset */
	public:
	explicit QuadrantT(IT n0, IT n1, IT n2, IT n3, IT n4):off_({n0,n1,n2,n3,n4}) {}
	IT noff(void)const{return off_[0];        }
	IT&o_at(IT i){chck(i);return off_[i];        }
	const IT o_at(IT i)const{return off_[i];        }
	const IT n_at(IT i)const{chck(i);chck(i-1);return (o_at(i)-o_at(i-1));}
	IT n_ul(void)const{return off_[1]-off_[0];}
	IT n_ur(void)const{return off_[2]-off_[1];}
	IT n_ll(void)const{return off_[3]-off_[2];}
	IT n_lr(void)const{return off_[4]-off_[3];}
	IT ntot(void)const{return off_[4]-off_[0];}
	IT diff(IT n)const{chck(n);chck(n-1);return off_[n]-off_[n-1];}
	IT sum(void)const { IT tot{0}; tot+=std::accumulate(off_.begin(),off_.end(),tot); return tot;}
	void add(IT v) { std::for_each(off_.begin(),off_.end(),[&v](IT & x){ x+=v; });}
	IT & operator[](IT i)  {return o_at(i);}
	const IT & operator[](IT i)const{return off_[i];}
	void chck(IT n)const{assert(n>=0);assert(n<=Four);}
}; /* QuadrantT class */

template <typename IT>
inline std::ostream & operator << (std::ostream & os, const QuadrantT<IT> & q)
{
	//os <<  "[" << "@ " << q.noff() << " : " << q.n_ul() << " " << q.n_ur() << " " << q.n_ll() << " " << q.n_lr() << "]" ;
	os <<  "[" << "@ " << q.o_at(0) << " :";
	for (auto i = 1 ; i<=4; ++i )
       		os << " " << q.o_at(i);
	os << "]" ;
	return os;
}
/* class QuadrantT end */

/* class LocalSubmatrixT begin */
template <typename IT>
struct LocalSubmatrixT final
{
	public:
	IT fnz,nnz,fr,nr,fc,nc;
	//LocalSubmatrixT(IT nnz_, IT nr_, IT nc_): nnz(nnz_),nr(nr_),nc(nc_){}
	
	bool intersects_on_rows(const LocalSubmatrixT & lsm) const
	{
		return std::max(fr,lsm.fr) < std::min(fr+nr,lsm.fr+lsm.nr);
	}
	bool intersects_on_cols(const LocalSubmatrixT & lsm) const
	{
		return std::max(fc,lsm.fc) < std::min(fc+nc,lsm.fc+lsm.nc);
	}
	bool operator < (const LocalSubmatrixT & lsm) const { return ( nnz < lsm.nnz || nr < lsm.nr || nc < lsm.nc ); }
	bool fits_csr(void) const
	{
		return ( nnz > nr );
	}
	bool fits_short_indices(void) const
	{
		constexpr auto msi = std::numeric_limits<short int>::max();
		return (nr < msi && nc < msi );
	}
};

template <typename IT>
inline std::ostream & operator << (std::ostream & os, const LocalSubmatrixT<IT> & lsm)
{
	os << " nz:" << lsm.fnz << ".." << lsm.fnz+lsm.nnz;
	os <<  " r:" << lsm.fr  << ".." << lsm.fr+lsm.nr;
	os <<  " c:" << lsm.fc  << ".." << lsm.fc+lsm.nc;
	return os;
}
/* class LocalSubmatrixT end */

/* class NodeT begin */

template <typename IT>
//using NodeT = std::array<QuadrantT<IT>,2>;
using NodeT = std::tuple<LocalSubmatrixT<IT>,QuadrantT<IT>>;
//using NodeT = std::tuple<QuadrantT<IT>,QuadrantT<IT>>;
//using NodeT = std::pair<QuadrantT<IT>,QuadrantT<IT>>;

template <typename IT>
bool isLeaf(const NodeT<IT> & node)
{
	// TODO: unfinished
	return std::get<1>(node).ntot()==0;
}

template <typename IT>
inline std::ostream & operator << (std::ostream & os, const NodeT<IT> & qp)
{
	const auto & lsm = std::get<0>(qp);
	const auto & smq = std::get<1>(qp);

	if(isLeaf(qp))
		//os << "(" << "l:" << " nnz:" << qnz[0] << ".." << qnz[4] << " " << "subm." << smq[0] << ")";
		os << "(" << "l:" << lsm << " " << "subm." << smq[0] << ")";
		//os << "(" << "l:" << " nnz:" << lsm.fnz << "..+" << lsm.nnz << " " << "subm." << smq[0] << ")";
	else
		//os << "(" << "n:" << " nnz:" << lsm.nnz << " " << "subm:" << smq << ")";
		os << "(" << "n:" << lsm << " " << "subm." << smq    << ")";
		//os << "(" << "n:" << " nnz:" << qnz << " " << "subm:" << smq << ")";
	return os;
}
/* class NodeT end */

inline int rsbpp_getenv_int_t(const char * evid, int defval)
{
	const auto evval { getenv(evid) };
	int val = evval ? atoi( evval ) : defval;
	return val;
}

template <typename IT, typename NT, typename CC, typename XC, typename YC, const int UF=4, const int nrhs=1>
rsb_err_t spmm_cooi_partial_unrolled_by_cols(const CC & coo, rsb_flags_t flags_, rsb_trans_t transA, const IT ldX, XC & x, const IT ldY, YC & y, const IT frhsi, const NT alpha, const IT bidx, const IT eidx_)
{
	/* TODO: Slated for removal. */
	const int uf = UF;
	auto nzi = bidx;
	const IT eidx=eidx_;
	const auto roff = coo.roff();
	const auto coff = coo.coff();
	auto const * __restrict const xp  = & x[(frhsi*ldX)];
	auto       * __restrict const yp  = & y[(frhsi*ldY)];
	const rsb_flags_t flags{flags_};

	if (( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) ||
	      RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) ) &&
		( coo.nr() != coo.nc() && roff == coff ) )
		return RSB_ERR_BADARGS;

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_N )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
		}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_T )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
		}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_C )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
		}

	if ( transA == RSB_TRANSPOSITION_N )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		constexpr bool assume_diagonal_alignment = true; // would benefit in hermitian loops as well

		if(assume_diagonal_alignment)
		{
			auto yn = & yp[roff];
			auto yt = & yp[coff];
			const auto xn = & xp[coff];
			const auto xt = & xp[roff];

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);

				if(i==j)
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i+nrhsi*ldY]+=alpha*v*xn[j+nrhsi*ldX];
				else
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i+nrhsi*ldY]+=alpha*v*xn[j+nrhsi*ldX],
						yt[j+nrhsi*ldY]+=alpha*v*xt[i+nrhsi*ldX];
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i+nrhsi*ldY]+=alpha*v*xn[j+nrhsi*ldX],
						yt[j+nrhsi*ldY]+=alpha*v*xt[i+nrhsi*ldX];
			}
		}
		else // simpler but a bit less efficient
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX],
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
		}
	}

	if ( transA == RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX],
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
		}
	}

	if ( transA == RSB_TRANSPOSITION_C )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX],
					yp[roff+i+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX];
		}
	}

	if ( transA == RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX];
			else
			{
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX],
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
			}
		}
	}

	if ( transA != RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
			else
			{
				const NT c = rsbpp_conj(v);
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX],
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
			}
		}
	}

	if( nzi == bidx && eidx - bidx > uf )
	       return RSB_ERR_UNSUPPORTED_OPERATION;

	assert( ! ( nzi < eidx && uf == 1 ) );
	if( nzi < eidx )
		spmm_cooi_partial_unrolled_by_cols<IT,NT,CC,XC,YC,1,nrhs>(coo,flags_,transA,ldX,x,ldY,y,frhsi,alpha,nzi,eidx); /*the remainder */
	return RSB_ERR_NO_ERROR;
}
#endif /* RSBPP_WANT_ALL */

template <typename IT, typename NT, typename CC, typename XC, typename YC, const int UF=4, const int nrhs=1>
rsb_err_t spmm_coo_partial_unrolled_by_cols(const CC & coo, rsb_flags_t flags_, rsb_trans_t transA, const IT ldX, XC & x, const IT ldY, YC & y, const IT frhsi, const NT alpha, const IT bidx, const IT eidx_)
{
	const int uf = UF;
	auto nzi = bidx;
	const IT eidx=eidx_;
	const auto roff = coo.roff();
	const auto coff = coo.coff();
	auto const * __restrict const xp  = & x[(frhsi*ldX)];
	auto       * __restrict const yp  = & y[(frhsi*ldY)];
	const rsb_flags_t flags{flags_};

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_N )
	{
		const int use_lol = 4; // loop optimization level
		auto const * __restrict const va = coo.va_;
		auto const * __restrict const ia = coo.ia_;
		auto const * __restrict const ja = coo.ja_;

		if( use_lol == 0 )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
		}

		if( use_lol == 1 )
		{
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=ia[nzi+ui], j=ja[nzi+ui];
				const NT v=va[nzi+ui];
				const NT av=v*alpha;

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=av*xp[coff+j+nrhsi*ldX];
			}
		}

		if( use_lol == 2 )
		{
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=ia[nzi+ui], j=ja[nzi+ui];
				const NT v=va[nzi+ui];
				const NT av=v*alpha;

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=av*xp[coff+j+nrhsi*ldX];
			}
		}

		if( use_lol == 3 )
		{
			auto const * __restrict const xc = xp + coff;
			auto       * __restrict const yr = yp + roff;

			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=ia[nzi+ui], j=ja[nzi+ui];
				const NT v=va[nzi+ui];
				const NT av=v*alpha;

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i+nrhsi*ldY]+=av*xc[j+nrhsi*ldX];
			}
		}

		if( use_lol == 4 )
		{
			auto const * __restrict const xc = xp + coff;
			auto       * __restrict const yr = yp + roff;

			std::array<NT,uf> au;
			std::array<IT,uf> iu;
			std::array<IT,uf> ju;

			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					iu[ui]=ia[nzi+ui];
					ju[ui]=ja[nzi+ui];
					au[ui]=alpha*va[nzi+ui];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				for (auto ui = 0 ; ui < uf ; ++ui )
					yr[iu[ui]+nrhsi*ldY]+=au[ui]*xc[ju[ui]+nrhsi*ldX];
			}
		}
	}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_T )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
		}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_C )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
		}

	if ( transA == RSB_TRANSPOSITION_N )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		constexpr bool assume_diagonal_alignment = true; // would benefit in hermitian loops as well

		if(assume_diagonal_alignment)
		{
			auto yn = & yp[roff];
			auto yt = & yp[coff];
			const auto xn = & xp[coff];
			const auto xt = & xp[roff];

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);

				if(i==j)
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i+nrhsi*ldY]+=alpha*v*xn[j+nrhsi*ldX];
				else
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i+nrhsi*ldY]+=alpha*v*xn[j+nrhsi*ldX],
						yt[j+nrhsi*ldY]+=alpha*v*xt[i+nrhsi*ldX];
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i+nrhsi*ldY]+=alpha*v*xn[j+nrhsi*ldX],
						yt[j+nrhsi*ldY]+=alpha*v*xt[i+nrhsi*ldX];
			}
		}
		else // simpler but a bit less efficient
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX],
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
		}
	}

	if ( transA == RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX],
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
		}
	}

	if ( transA == RSB_TRANSPOSITION_C )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX],
					yp[roff+i+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX];
		}
	}

	if ( transA == RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX];
			else
			{
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX],
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
			}
		}
	}

	if ( transA != RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
			else
			{
				const NT c = rsbpp_conj(v);
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX],
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
			}
		}
	}

	if( nzi == bidx && eidx - bidx >= uf )
		return RSB_ERR_BADARGS;

	if( nzi < eidx )
		return spmm_coo_partial_unrolled_by_cols<IT,NT,CC,XC,YC,1,nrhs>(coo,flags_,transA,ldX,x,ldY,y,frhsi,alpha,nzi,eidx); /*the remainder */
	return RSB_ERR_NO_ERROR;
}

template <typename IT, typename CT, typename NT, typename CC, typename XC, typename YC, const int UF=4, const int nrhs=1>
rsb_err_t spmm_coo_partial_unrolled_by_rows(const CC & coo, rsb_flags_t flags_, rsb_trans_t transA, const IT ldX, XC & x, const IT ldY, YC & y, const IT frhsi, const NT alpha, const IT bidx, const IT eidx_)
{
	const int uf = UF;
	auto nzi = bidx;
	const IT eidx=eidx_;
	const auto roff = coo.roff();
	const auto coff = coo.coff();
	auto const * __restrict const xp  = & x[frhsi];
	auto       * __restrict const yp  = & y[frhsi];
	const rsb_flags_t flags{flags_};
	constexpr NT one {1};
	constexpr NT zero {};

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_N )
	{
		const int use_lol = 4; // loop optimization level
		auto const * __restrict const va = coo.va_;
		auto const * __restrict const ia = coo.ia_;
		auto const * __restrict const ja = coo.ja_;

		if( use_lol == 0 )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
		}

		if( use_lol == 1 )
		{
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=ia[nzi+ui], j=ja[nzi+ui];
				const NT v=va[nzi+ui];
				const NT av=v*alpha;

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=av*xp[(coff+j)*ldX+nrhsi];
			}
		}

		if( use_lol == 2 )
		{
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=ia[nzi+ui], j=ja[nzi+ui];
				const NT v=va[nzi+ui];
				const NT av=v*alpha;

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=av*xp[(coff+j)*ldX+nrhsi];
			}
		}

		if( use_lol == 3 )
		{
			auto const * __restrict const xc = xp + coff*ldY;
			auto       * __restrict const yr = yp + roff*ldX;

			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=ia[nzi+ui], j=ja[nzi+ui];
				const NT v=va[nzi+ui];
				const NT av=v*alpha;

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i*ldY+nrhsi]+=av*xc[j*ldX+nrhsi];
			}
		}

		if( use_lol == 4 )
		{
			auto const * __restrict const xc = xp + coff*ldY;
			auto       * __restrict const yr = yp + roff*ldY;

			std::array<NT,uf> au;
			std::array<IT,uf> iu;
			std::array<IT,uf> ju;

			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					iu[ui]=ia[nzi+ui]*ldY;
					ju[ui]=ja[nzi+ui]*ldX;
					au[ui]=alpha*va[nzi+ui];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				for (auto ui = 0 ; ui < uf ; ++ui )
					yr[iu[ui]+nrhsi]+=au[ui]*xc[ju[ui]+nrhsi];
			}
		}
	}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_T )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[(coff+j)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
		}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_C )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[(coff+j)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
		}

	if ( transA != RSB_TRANSPOSITION_C )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		const int use_lol = 5; // loop optimization level

		if ( use_lol == -1 )
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi],
					yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
		}

		constexpr bool assume_diagonal_alignment = true; // would benefit in hermitian loops as well

		if ( use_lol == 0 )
		{
		if(assume_diagonal_alignment)
		{
			auto yn = & yp[roff*ldY];
			auto yt = & yp[coff*ldY];
			const auto xn = & xp[coff*ldX];
			const auto xt = & xp[roff*ldX];

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);

				if(i==j)
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i*ldY+nrhsi]+=alpha*v*xn[j*ldX+nrhsi];
				else
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i*ldY+nrhsi]+=alpha*v*xn[j*ldX+nrhsi],
						yt[j*ldY+nrhsi]+=alpha*v*xt[i*ldX+nrhsi];
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i*ldY+nrhsi]+=alpha*v*xn[j*ldX+nrhsi],
						yt[j*ldY+nrhsi]+=alpha*v*xt[i*ldX+nrhsi];
			}
		}
		else // simpler but a bit less efficient
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi],
					yp[(coff+j)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
		}
		}

		if ( use_lol == 1 )
		{
		if(assume_diagonal_alignment)
		{
			auto yn = & yp[roff*ldY];
			auto yt = & yp[coff*ldY];
			const auto xn = & xp[coff*ldX];
			const auto xt = & xp[roff*ldX];

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);
				const NT av=v*alpha;

				if(i==j)
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i*ldY+nrhsi]+=av*xn[j*ldX+nrhsi];
				else
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i*ldY+nrhsi]+=av*xn[j*ldX+nrhsi],
						yt[j*ldY+nrhsi]+=av*xt[i*ldX+nrhsi];
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			for (auto ui = 0 ; ui < uf ; ++ui )
			{
				const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
				const NT v=coo.va(nzi+ui);
				const NT av=v*alpha;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[i*ldY+nrhsi]+=av*xn[j*ldX+nrhsi],
						yt[j*ldY+nrhsi]+=av*xt[i*ldX+nrhsi];
			}
		}
		else // simpler but a bit less efficient
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT av=v*alpha;

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=av*xp[(coff+j)*ldX+nrhsi];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=av*xp[(coff+j)*ldX+nrhsi],
					yp[(coff+j)*ldY+nrhsi]+=av*xp[(roff+i)*ldX+nrhsi];
		}
		}

		if ( use_lol == 2 )
		{
			auto yn = & yp[roff*ldY];
			auto yt = & yp[coff*ldY];
			const auto xn = & xp[coff*ldX];
			const auto xt = & xp[roff*ldX];
			auto const * __restrict const va = coo.va_;
			auto const * __restrict const ia = coo.ia_;
			auto const * __restrict const ja = coo.ja_;

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> vu;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					vu[ui]=va[nzi+ui],
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					const NT av = vu[ui]*alpha;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=av*xn[ju[ui]*ldX+nrhsi],
						yt[ju[ui]*ldY+nrhsi]+=av*xt[iu[ui]*ldX+nrhsi] * (iu[ui] == ju[ui] ?one:zero);
				}
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> vu;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					vu[ui]=va[nzi+ui],
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					const NT av = vu[ui]*alpha;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=av*xn[ju[ui]*ldX+nrhsi],
						yt[ju[ui]*ldY+nrhsi]+=av*xt[iu[ui]*ldX+nrhsi];
				}
			}
		}

		if ( use_lol == 3 )
		{
			auto yn = & yp[roff*ldY];
			auto yt = & yp[coff*ldY];
			const auto xn = & xp[coff*ldX];
			const auto xt = & xp[roff*ldX];
			auto const * __restrict const va = coo.va_;
			auto const * __restrict const ia = coo.ia_;
			auto const * __restrict const ja = coo.ja_;

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> vu;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					vu[ui]=va[nzi+ui],
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					const NT av = vu[ui]*alpha;
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=av*xn[ju[ui]*ldX+nrhsi];
				}

				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					const NT av = vu[ui]*alpha;
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yt[ju[ui]*ldY+nrhsi]+=av*xt[iu[ui]*ldX+nrhsi] * (iu[ui] != ju[ui] ?one:zero);
				}
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> vu;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					vu[ui]=va[nzi+ui],
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					const NT av = vu[ui]*alpha;
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=av*xn[ju[ui]*ldX+nrhsi];
				}

				for (auto ui = 0 ; ui < uf ; ++ui )
				{
					const NT av = vu[ui]*alpha;
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yt[ju[ui]*ldY+nrhsi]+=av*xt[iu[ui]*ldX+nrhsi];
				}
			}
		}

		if ( use_lol == 4 )
		{
			auto yn = & yp[roff*ldY];
			auto yt = & yp[coff*ldY];
			const auto xn = & xp[coff*ldX];
			const auto xt = & xp[roff*ldX];
			auto const * __restrict const va = coo.va_;
			auto const * __restrict const ia = coo.ia_;
			auto const * __restrict const ja = coo.ja_;

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> au;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					au[ui]=va[nzi+ui]*alpha,
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=au[ui]*xn[ju[ui]*ldX+nrhsi];

				for (auto ui = 0 ; ui < uf ; ++ui )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yt[ju[ui]*ldY+nrhsi]+=au[ui]*xt[iu[ui]*ldX+nrhsi] * (iu[ui] != ju[ui] ?one:zero);
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> au;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					au[ui]=va[nzi+ui]*alpha,
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=au[ui]*xn[ju[ui]*ldX+nrhsi];

				for (auto ui = 0 ; ui < uf ; ++ui )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yt[ju[ui]*ldY+nrhsi]+=au[ui]*xt[iu[ui]*ldX+nrhsi];
			}
		}

		if ( use_lol == 5 )
		{
			auto yn = & yp[roff*ldY];
			auto yt = & yp[coff*ldY];
			const auto xn = & xp[coff*ldX];
			const auto xt = & xp[roff*ldX];
			auto const * __restrict const va = coo.va_;
			auto const * __restrict const ia = coo.ia_;
			auto const * __restrict const ja = coo.ja_;

			if(roff==coff)
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> au;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					au[ui]=va[nzi+ui]*alpha,
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=au[ui]*xn[ju[ui]*ldX+nrhsi];

				for (auto ui = 0 ; ui < uf ; ++ui )
					if( iu[ui] != ju[ui] )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yt[ju[ui]*ldY+nrhsi]+=au[ui]*xt[iu[ui]*ldX+nrhsi];
			}
			else
			for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
			{
				std::array<NT,uf> au;
				std::array<IT,uf> iu;
				std::array<IT,uf> ju;

				for (auto ui = 0 ; ui < uf ; ++ui )
					au[ui]=va[nzi+ui]*alpha,
					iu[ui]=ia[nzi+ui],
					ju[ui]=ja[nzi+ui];

				for (auto ui = 0 ; ui < uf ; ++ui )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yn[iu[ui]*ldY+nrhsi]+=au[ui]*xn[ju[ui]*ldX+nrhsi];

				for (auto ui = 0 ; ui < uf ; ++ui )
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yt[ju[ui]*ldY+nrhsi]+=au[ui]*xt[iu[ui]*ldX+nrhsi];
			}
		}
	}

	if ( transA == RSB_TRANSPOSITION_C )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi],
					yp[(roff+i)*ldY+nrhsi]+=alpha*c*xp[(coff+j)*ldX+nrhsi];
		}
	}

	if ( transA != RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
			else
			{
				const NT c = rsbpp_conj(v);
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi],
					yp[(coff+j)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
			}
		}
	}

	if ( transA == RSB_TRANSPOSITION_T )
	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN ) )
	{
		for ( ; nzi+(uf-1) < eidx ; nzi+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=coo.ia(nzi+ui), j=coo.ja(nzi+ui);
			const NT v=coo.va(nzi+ui);
			const NT c=rsbpp_conj(v);

			if(i+roff==j+coff)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=alpha*c*xp[(coff+j)*ldX+nrhsi];
			else
			{
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=alpha*c*xp[(coff+j)*ldX+nrhsi],
					yp[(coff+j)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
			}
		}
	}

	if( nzi == bidx && eidx - bidx >= uf )
		return RSB_ERR_BADARGS;

	if( nzi < eidx )
		return spmm_coo_partial_unrolled_by_rows<IT,CT,NT,CC,XC,YC,1,nrhs>(coo,flags_,transA,ldX,x,ldY,y,frhsi,alpha,nzi,eidx); /*the remainder */
	return RSB_ERR_NO_ERROR;
}

template <typename IT, typename CT, typename NT, typename CC, typename XC, typename YC, const int UF=4, const int nrhs=1>
rsb_err_t spmm_csr_partial_unrolled_by_rows(const CC & csr, rsb_flags_t flags_, rsb_trans_t transA, const IT ldX, XC & x, const IT ldY, YC & y, const IT frhsi, const NT alpha, const IT bidx, const IT eidx_)
{
	const int uf = UF;
	auto nri = bidx;
	const IT eidx = eidx_;
	const auto roff = csr.roff();
	const auto coff = csr.coff();
	auto const * __restrict const xp  = & x[(frhsi)];
	auto       * __restrict const yp  = & y[(frhsi)];
	const rsb_flags_t flags{flags_};
	constexpr NT none {1};

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if( transA == RSB_TRANSPOSITION_N )
	{
		const int use_lol = 2; // loop optimization level

		if ( use_lol == -1)
		{
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				auto lnz { csr.ip(nri+1)};

				for (                       ; nzi < lnz; nzi+=1 )
				{
					const IT j = csr.ja(nzi);
					const NT v = csr.va(nzi);

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
				}
			}
			return RSB_ERR_NO_ERROR;
		}

		if ( use_lol == 2 )
		{
			auto const * __restrict const xc = xp + coff * ldX;
			auto       * __restrict const yr = yp + roff * ldY;
			auto const * __restrict const va = csr.va_;
			auto const * __restrict const ja = csr.ja_;
			auto const * __restrict const ip = csr.ip_;

			if ( alpha == none )
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi+uf-1 < lnz; nzi+=uf )
				{
					std::array<NT,uf> vu;
					std::array<CT,uf> ju;

					for ( auto ui = 0 ; ui < uf ; ++ui )
						vu[ui]=va[nzi+ui];
					for ( auto ui = 0 ; ui < uf ; ++ui )
						ju[ui]=ja[nzi+ui];

					for ( auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					for ( auto ui = 0 ; ui < uf ; ++ui )
						racc[nrhsi]+=      vu[ui]*xc[ju[ui]*ldX+nrhsi];
				}

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = ja[nzi];
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=      v*xc[j*ldX+nrhsi];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i*ldY+nrhsi]+=racc[nrhsi];
			}
			else
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi+uf-1 < lnz; nzi+=uf )
				{
					std::array<NT,uf> vu;
					std::array<CT,uf> ju;

					for ( auto ui = 0 ; ui < uf ; ++ui )
						vu[ui]=va[nzi+ui];
					for ( auto ui = 0 ; ui < uf ; ++ui )
						ju[ui]=ja[nzi+ui];

					for ( auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					for ( auto ui = 0 ; ui < uf ; ++ui )
						racc[nrhsi]+=      vu[ui]*xc[ju[ui]*ldX+nrhsi];
				}

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = ja[nzi];
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=      v*xc[j*ldX+nrhsi];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i*ldY+nrhsi]+=alpha*racc[nrhsi];
			}
			return RSB_ERR_NO_ERROR;
		}

		if ( use_lol == 3 )
		{
			auto const * __restrict const xc = xp + coff * ldX;
			auto       * __restrict const yr = yp + roff * ldY;
			auto const * __restrict const va = csr.va_;
			auto const * __restrict const ja = csr.ja_;
			auto const * __restrict const ip = csr.ip_;

			if ( alpha == none )
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri*ldY};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi+uf-1 < lnz; nzi+=uf )
				{
					std::array<NT,uf> vu;
					std::array<IT,uf> ju;

					for ( auto ui = 0 ; ui < uf ; ++ui )
						vu[ui]=va[nzi+ui];
					for ( auto ui = 0 ; ui < uf ; ++ui )
						ju[ui]=ja[nzi+ui]*ldX;

					for ( auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					for ( auto ui = 0 ; ui < uf ; ++ui )
						racc[nrhsi]+=      vu[ui]*xc[ju[ui]+nrhsi];
				}

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const IT j = ja[nzi]*ldX;
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=      v*xc[j+nrhsi];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i+nrhsi]+=racc[nrhsi];
			}
			else
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri*ldY};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi+uf-1 < lnz; nzi+=uf )
				{
					std::array<NT,uf> vu;
					std::array<IT,uf> ju;

					for ( auto ui = 0 ; ui < uf ; ++ui )
						vu[ui]=va[nzi+ui];
					for ( auto ui = 0 ; ui < uf ; ++ui )
						ju[ui]=ja[nzi+ui]*ldX;

					for ( auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					for ( auto ui = 0 ; ui < uf ; ++ui )
						racc[nrhsi]+=      vu[ui]*xc[ju[ui]+nrhsi];
				}

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const IT j = ja[nzi]*ldX;
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=      v*xc[j+nrhsi];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i+nrhsi]+=alpha*racc[nrhsi];
			}
			return RSB_ERR_NO_ERROR;
		}
		return RSB_ERR_INTERNAL_ERROR;
	}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_T )
	{
		for ( ; nri < eidx ; ++nri )
		{
			const auto i {nri};
			auto nzi { csr.ip(nri)};
			auto lnz { csr.ip(nri+1)};

			for (                       ; nzi+uf-1 < lnz; nzi+=uf )
			for ( auto ui = 0 ; ui < uf ; ++ui )
			{
				const CT j = csr.ja(nzi+ui);
				const NT v = csr.va(nzi+ui);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
			}

			for (                       ; nzi      < lnz; nzi+= 1 )
			{
				const CT j = csr.ja(nzi);
				const NT v = csr.va(nzi);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
			}
		}
		return RSB_ERR_NO_ERROR;
	}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_C )
	{
		for ( ; nri < eidx ; ++nri )
		{
			const auto i {nri};
			auto nzi { csr.ip(nri)};
			auto lnz { csr.ip(nri+1)};

			for (                       ; nzi+uf-1 < lnz; nzi+=uf )
			for ( auto ui = 0 ; ui < uf ; ++ui )
			{
				const CT j = csr.ja(nzi+ui);
				const NT v = csr.va(nzi+ui);
				const NT c = rsbpp_conj(v);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
			}

			for (                       ; nzi      < lnz; nzi+= 1 )
			{
				const CT j = csr.ja(nzi);
				const NT v = csr.va(nzi);
				const NT c = rsbpp_conj(v);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(coff+j)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
			}
		}
		return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if ( transA != RSB_TRANSPOSITION_C )
	{
		const int use_lol = 2; // loop optimization level

		if ( use_lol == -1 )
		for ( ; nri < eidx ; ++nri )
		{
			const auto i {nri};
			auto nzi { csr.ip(nri)};
			const auto lnz { csr.ip(nri+1)};

			for (                       ; nzi      < lnz; nzi+= 1 )
			{
				const IT j = csr.ja(nzi);
				const NT v = csr.va(nzi);

				if(true)
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
				if(roff+i != coff+j)
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[(j+coff)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
			}
		}

		if ( use_lol == 2 )
		{
			{
				auto const * __restrict const xr = xp + roff * ldY;
				auto const * __restrict const xc = xp + coff * ldX;
				auto       * __restrict const yr = yp + roff * ldY;
				auto       * __restrict const yc = yp + coff * ldX;
				auto const * __restrict const va = csr.va_;
				auto const * __restrict const ja = csr.ja_;
				auto const * __restrict const ip = csr.ip_;

				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { ip[nri]};
					const auto lnz { ip[nri+1]};
					std::array<NT,nrhs> racc{{}};
					std::array<NT,nrhs> alxt;

					if(nzi==lnz)
						continue;

					{
						const CT j = ja[nzi];
						const NT v = va[nzi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j*ldX+nrhsi];
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							alxt[nrhsi] = alpha*xr[i*ldX+nrhsi];

						if(i+roff!=j+coff)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yc[j*ldY+nrhsi]+=v*alxt[nrhsi];
						nzi++;
					}

					for (                       ; nzi+1      < lnz; nzi+= 1 )
					{
						const CT j = ja[nzi];
						const NT v = va[nzi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j*ldX+nrhsi],
							yc[j*ldY+nrhsi]+=v*alxt[nrhsi];
					}

					if( nzi < lnz )
					{
						const CT j = ja[nzi];
						const NT v = va[nzi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j*ldX+nrhsi];
						if(i+roff!=j+coff)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yc[j*ldY+nrhsi]+=v*alxt[nrhsi];
						nzi++;
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i*ldY+nrhsi]+=alpha*racc[nrhsi];
				}
			}
		}

		if ( use_lol == 3 ) // seems inferior to 2
		{
			{
				auto const * __restrict const xr = xp + roff * ldY;
				auto const * __restrict const xc = xp + coff * ldX;
				auto       * __restrict const yr = yp + roff * ldY;
				auto       * __restrict const yc = yp + coff * ldX;
				auto const * __restrict const va = csr.va_;
				auto const * __restrict const ja = csr.ja_;
				auto const * __restrict const ip = csr.ip_;

				for ( ; nri < eidx ; ++nri )
				{
					auto nzi { ip[nri]};
					const auto lnz { ip[nri+1]};

					if(nzi==lnz)
						continue;

					const auto i {nri};
					std::array<NT,nrhs> racc{{}};
					std::array<NT,nrhs> alxt{{}};

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xr[i*ldX+nrhsi];

					{
						const CT j = ja[nzi];
						const NT v = va[nzi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j*ldX+nrhsi];

						if(i+roff!=j+coff)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yc[j*ldY+nrhsi]+=v*alxt[nrhsi];
						nzi++;
					}

					for (                       ; nzi+uf-1+1 < lnz; nzi+=uf )
					{
						std::array<NT,uf> vu;
						std::array<IT,uf> ju;

						for ( auto ui = 0 ; ui < uf ; ++ui )
							vu[ui]=va[nzi+ui];
						for ( auto ui = 0 ; ui < uf ; ++ui )
							ju[ui]=ja[nzi+ui];

						for ( auto ui = 0 ; ui < uf ; ++ui )
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								racc[nrhsi]+=vu[ui]*xc[ju[ui]*ldX+nrhsi];
						for ( auto ui = 0 ; ui < uf ; ++ui )
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yc[ju[ui]*ldY+nrhsi]+=vu[ui]*alxt[nrhsi];
					}

					for (                       ; nzi      < lnz; nzi+= 1 )
					{
						const CT j = ja[nzi];
						const NT v = va[nzi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j*ldX+nrhsi];

						if(i+roff!=j+coff)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yc[j*ldY+nrhsi]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i*ldY+nrhsi]+=alpha*racc[nrhsi];
				}
			}
		}
		return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if ( transA == RSB_TRANSPOSITION_C )
	{
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
					const NT c = rsbpp_conj(v);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[(j+coff)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
					if(roff+i!=coff+j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[(i+roff)*ldY+nrhsi]+=alpha*c*xp[(coff+j)*ldX+nrhsi];
				}
			}
			return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA != RSB_TRANSPOSITION_T )
	{
		const int use_lol = 1; // loop optimization level

			if ( use_lol == -1)
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const IT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
					const NT c = rsbpp_conj(v);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[(roff+i)*ldY+nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
					if(roff+i != coff+j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[(j+coff)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
				}
			}

			if ( use_lol != -1)
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
					const NT c = rsbpp_conj(v);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=alpha*v*xp[(coff+j)*ldX+nrhsi];
					if(roff+i != coff+j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[(j+coff)*ldY+nrhsi]+=alpha*c*xp[(roff+i)*ldX+nrhsi];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[(roff+i)*ldY+nrhsi]+=racc[nrhsi];
			}
			return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_T )
	{
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
					const NT c = rsbpp_conj(v);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[(roff+i)*ldY+nrhsi]+=alpha*c*xp[(coff+j)*ldX+nrhsi];
					if(roff+i != coff+j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[(j+coff)*ldY+nrhsi]+=alpha*v*xp[(roff+i)*ldX+nrhsi];
				}
			}
			return RSB_ERR_NO_ERROR;
	}

	return RSB_ERR_BADARGS;
}

template <typename IT, typename CT, typename NT, typename CC, typename XC, typename YC, const int UF=4, const int nrhs=1>
rsb_err_t spmm_csr_partial_unrolled_by_cols(const CC & csr, rsb_flags_t flags_, rsb_trans_t transA, const IT ldX, XC & x, const IT ldY, YC & y, const IT frhsi, const NT alpha, const IT bidx, const IT eidx_)
{
	const int uf = UF;
	auto nri = bidx;
	const IT eidx = eidx_;
	const auto roff = csr.roff();
	const auto coff = csr.coff();
	auto const * __restrict const xp  = & x[(frhsi*ldX)];
	auto       * __restrict const yp  = & y[(frhsi*ldY)];
	const rsb_flags_t flags{flags_};
	constexpr NT none {1};

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_N )
	{
		const int use_lol = 2; // loop optimization level

		if ( use_lol == -1)
		{
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				auto lnz { csr.ip(nri+1)};

				for (                       ; nzi < lnz; nzi+=1 )
				{
					const IT j = csr.ja(nzi);
					const NT v = csr.va(nzi);

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
				}
			}
			return RSB_ERR_NO_ERROR;
		}

		if ( use_lol == 0)
		{
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				auto lnz { csr.ip(nri+1)};
				std::array<NT,nrhs> racc{{}};
	
				for (                       ; nzi+uf-1 < lnz; nzi+=uf )
				for ( auto ui = 0 ; ui < uf ; ++ui )
				{
					const CT j = csr.ja(nzi+ui);
					const NT v = csr.va(nzi+ui);
	
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=alpha*v*xp[coff+j+nrhsi*ldX];
				}
	
				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
	
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=alpha*v*xp[coff+j+nrhsi*ldX];
				}
	
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[roff+i+nrhsi*ldY]+=racc[nrhsi];
			}
			return RSB_ERR_NO_ERROR;
		}

		if ( use_lol == 1 )
		{
			auto const * __restrict const xc = xp + coff;
			auto       * __restrict const yr = yp + roff;
			auto const * __restrict const va = csr.va_;
			auto const * __restrict const ja = csr.ja_;
			auto const * __restrict const ip = csr.ip_;

			if ( alpha == none )
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = ja[nzi];
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=v*xc[j+nrhsi*ldX];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i+nrhsi*ldY]+=      racc[nrhsi];
			}
			else
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = ja[nzi];
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=v*xc[j+nrhsi*ldX];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
			}
			return RSB_ERR_NO_ERROR;
		}

		if ( use_lol == 2 )
		{
			auto const * __restrict const xc = xp + coff;
			auto       * __restrict const yr = yp + roff;
			auto const * __restrict const va = csr.va_;
			auto const * __restrict const ja = csr.ja_;
			auto const * __restrict const ip = csr.ip_;

			if ( alpha == none )
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi+uf-1 < lnz; nzi+=uf )
				{
					std::array<NT,uf> vu;
					std::array<IT,uf> ju;

					for ( auto ui = 0 ; ui < uf ; ++ui )
						vu[ui]=va[nzi+ui];
					for ( auto ui = 0 ; ui < uf ; ++ui )
						ju[ui]=ja[nzi+ui];

					for ( auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					for ( auto ui = 0 ; ui < uf ; ++ui )
						racc[nrhsi]+=      vu[ui]*xc[ju[ui]+nrhsi*ldX];
				}

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const IT j = ja[nzi];
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=      v*xc[j+nrhsi*ldX];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i+nrhsi*ldY]+=racc[nrhsi];
			}
			else
			for ( ; nri < eidx ; ++nri )
			{
				const auto i {nri};
				auto nzi { ip[nri]};
				const auto lnz { ip[nri+1]};
				std::array<NT,nrhs> racc{{}};

				for (                       ; nzi+uf-1 < lnz; nzi+=uf )
				{
					std::array<NT,uf> vu;
					std::array<CT,uf> ju;

					for ( auto ui = 0 ; ui < uf ; ++ui )
						vu[ui]=va[nzi+ui];
					for ( auto ui = 0 ; ui < uf ; ++ui )
						ju[ui]=ja[nzi+ui];

					for ( auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					for ( auto ui = 0 ; ui < uf ; ++ui )
						racc[nrhsi]+=      vu[ui]*xc[ju[ui]+nrhsi*ldX];
				}

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const IT j = ja[nzi];
					const NT v = va[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=      v*xc[j+nrhsi*ldX];
				}

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
			}
			return RSB_ERR_NO_ERROR;
		}
		return RSB_ERR_INTERNAL_ERROR;
	}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_T )
	{
		for ( ; nri < eidx ; ++nri )
		{
			const auto i {nri};
			auto nzi { csr.ip(nri)};
			auto lnz { csr.ip(nri+1)};

			for (                       ; nzi+uf-1 < lnz; nzi+=uf )
			for ( auto ui = 0 ; ui < uf ; ++ui )
			{
				const CT j = csr.ja(nzi+ui);
				const NT v = csr.va(nzi+ui);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
			}

			for (                       ; nzi      < lnz; nzi+= 1 )
			{
				const CT j = csr.ja(nzi);
				const NT v = csr.va(nzi);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
			}
		}
		return RSB_ERR_NO_ERROR;
	}

	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if (!RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_C )
	{
		for ( ; nri < eidx ; ++nri )
		{
			const auto i {nri};
			auto nzi { csr.ip(nri)};
			auto lnz { csr.ip(nri+1)};

			for (                       ; nzi+uf-1 < lnz; nzi+=uf )
			for ( auto ui = 0 ; ui < uf ; ++ui )
			{
				const CT j = csr.ja(nzi+ui);
				const NT v = csr.va(nzi+ui);
				const NT c = rsbpp_conj(v);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
			}

			for (                       ; nzi      < lnz; nzi+= 1 )
			{
				const CT j = csr.ja(nzi);
				const NT v = csr.va(nzi);
				const NT c = rsbpp_conj(v);

				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[coff+j+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
			}
		}
		return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if ( transA != RSB_TRANSPOSITION_C )
	{
		if ( csr.roff() == csr.coff() )
		{
			if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_LOWER) )
			{
				// TODO: need to honour uf here as well
				const int use_lol = 1; // loop optimization level

				if ( use_lol == -1 )
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};

					for (                       ; nzi      < lnz; nzi+= 1 )
					{
						const IT j = csr.ja(nzi);
						const NT v = csr.va(nzi);

						if(true)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
						if(i != j)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
					}
				}

				if ( use_lol == 0 )
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};
					std::array<NT,nrhs> racc{{}};

					if(nzi==lnz)
						continue;

					for (                       ; nzi+1      < lnz; nzi+= 1 )
					{
						const CT j = csr.ja(nzi);
						const NT v = csr.va(nzi);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=alpha*v*xp[coff+j+nrhsi*ldX],
							yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
					}

					const CT j = csr.ja(nzi);
					const NT v = csr.va(nzi);

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						racc[nrhsi]+=alpha*v*xp[coff+j+nrhsi*ldX];
					if(i!=j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[roff+i+nrhsi*ldY]+=racc[nrhsi];
				}

				if ( use_lol == 1 )
				{
					auto const * __restrict const xr = xp + roff;
					auto const * __restrict const xc = xp + coff;
					auto       * __restrict const yr = yp + roff;
					auto       * __restrict const yc = yp + coff;
					auto const * __restrict const va = csr.va_;
					auto const * __restrict const ja = csr.ja_;
					auto const * __restrict const ip = csr.ip_;

					for ( ; nri < eidx ; ++nri )
					{
						const auto i {nri};
						auto nzi { ip[nri]};
						const auto lnz { ip[nri+1]};
						std::array<NT,nrhs> racc{{}};
						std::array<NT,nrhs> alxt;

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							alxt[nrhsi] = alpha*xr[i+nrhsi*ldX];

						if(nzi==lnz)
							continue;

						for (                       ; nzi+1      < lnz; nzi+= 1 )
						{
							const CT j = ja[nzi];
							const NT v = va[nzi];

							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								racc[nrhsi]+=v*xr[j+nrhsi*ldX],
								yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
						}

						const CT j = ja[nzi];
						const NT v = va[nzi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX];
						if(i!=j)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yc[j+nrhsi*ldY]+=v*alxt[nrhsi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
					}
				}
			}
			else
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const IT j = csr.ja(nzi);
					const NT v = csr.va(nzi);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
					if(i!=j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
				}
			}
			return RSB_ERR_NO_ERROR;
		}
		else
		{
			const int use_lol = 7; // loop optimization level
			const bool with_peel_loop = (uf-1>0);

			if ( use_lol == -1 )
			{
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};

					for (                       ; nzi      < lnz; nzi+= 1 )
					{
						const IT j = csr.ja(nzi);
						const NT v = csr.va(nzi);

						if(true)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];

						if(roff+i != coff+j)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
					}
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 0 )
			{
				if ( alpha != none )
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};
					std::array<NT,nrhs> racc{{}};

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					for ( auto ui = 0 ; ui < uf ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=alpha*v*xp[coff+j+nrhsi*ldX];
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
					}

					if(with_peel_loop)
					for (                       ; nzi      < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1  ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=alpha*v*xp[coff+j+nrhsi*ldX];
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
					}
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[roff+i+nrhsi*ldY]+=racc[nrhsi];
				}
				else
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};
					std::array<NT,nrhs> racc{{}};

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					for ( auto ui = 0 ; ui < uf ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=      v*xp[coff+j+nrhsi*ldX];
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=      v*xp[roff+i+nrhsi*ldX];
					}

					if(with_peel_loop)
					for (                       ; nzi      < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1  ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=      v*xp[coff+j+nrhsi*ldX];
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=      v*xp[roff+i+nrhsi*ldX];
					}
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[roff+i+nrhsi*ldY]+=racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 1 )
			{
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};
					std::array<NT,nrhs> racc{{}};
					std::array<NT,nrhs> alxt;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xp[roff+i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					for ( auto ui = 0 ; ui < uf ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xp[coff+j+nrhsi*ldX],
							yp[j+coff+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi      < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1 ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xp[coff+j+nrhsi*ldX],
							yp[j+coff+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[roff+i+nrhsi*ldY]+=alpha*racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 2 )
			{
				auto const * __restrict const xr = xp + roff;
				auto const * __restrict const xc = xp + coff;

				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};
					std::array<NT,nrhs> racc{{}};
					std::array<NT,nrhs> alxt;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xr[i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					for ( auto ui = 0 ; ui < uf ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yp[j+coff+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi      < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yp[j+coff+nrhsi*ldY]+=v*alxt[nrhsi];
					}
					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yp[roff+i+nrhsi*ldY]+=alpha*racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 3 )
			{
				auto const * __restrict const xr = xp + roff;
				auto const * __restrict const xc = xp + coff;
				auto       * __restrict const yr = yp + roff;
				auto       * __restrict const yc = yp + coff;

				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};
					std::array<NT,nrhs> racc{{}};
					std::array<NT,nrhs> alxt;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xr[i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					for ( auto ui = 0 ; ui < uf ; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1; ++ui )
					{
						const CT j = csr.ja(nzi+ui);
						const NT v = csr.va(nzi+ui);

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 4 )
			{
				auto const * __restrict const xr = xp + roff;
				auto const * __restrict const xc = xp + coff;
				auto       * __restrict const yr = yp + roff;
				auto       * __restrict const yc = yp + coff;
				auto       * __restrict const va = csr.va_;
				auto       * __restrict const ja = csr.ja_;

				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { csr.ip(nri)};
					const auto lnz { csr.ip(nri+1)};
					std::array<NT,nrhs> racc{{}};
					std::array<NT,nrhs> alxt;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xr[i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					for ( auto ui = 0 ; ui < uf ; ++ui )
					{
						const CT j = ja[nzi+ui];
						const NT v = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi      < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1; ++ui )
					{
						const CT j = ja[nzi+ui];
						const NT v = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 5 )
			{
				auto const * __restrict const xr = xp + roff;
				auto const * __restrict const xc = xp + coff;
				auto       * __restrict const yr = yp + roff;
				auto       * __restrict const yc = yp + coff;
				auto const * __restrict const va = csr.va_;
				auto const * __restrict const ja = csr.ja_;
				auto const * __restrict const ip = csr.ip_;

				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					auto nzi { ip[nri]};
					const auto lnz { ip[nri+1]};
					std::array<NT,nrhs> racc{{}};
					std::array<NT,nrhs> alxt;

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xr[i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					for ( auto ui = 0 ; ui < uf ; ++ui )
					{
						const CT j = ja[nzi+ui];
						const NT v = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1; ++ui )
					{
						const CT j = ja[nzi+ui];
						const NT v = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 6 )
			{
				auto const * __restrict const xr = xp + roff;
				auto const * __restrict const xc = xp + coff;
				auto       * __restrict const yr = yp + roff;
				auto       * __restrict const yc = yp + coff;
				auto const * __restrict const va = csr.va_;
				auto const * __restrict const ja = csr.ja_;
				auto const * __restrict const ip = csr.ip_;
				std::array<NT,nrhs> alxt;
				auto nzi { ip[nri]};

				if ( alpha != none )
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					const auto lnz { ip[nri+1]};
					std::array<NT,nrhs> racc{{}};

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xr[i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					{
						std::array<IT,uf> jr;
						std::array<NT,uf> vr;

						for ( auto ui = 0 ; ui < uf ; ++ui )
							jr[ui] = ja[nzi+ui],
							vr[ui] = va[nzi+ui];

						for ( auto ui = 0 ; ui < uf ; ++ui )
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								racc[nrhsi]+=vr[ui]*xc[jr[ui]+nrhsi*ldX],
								yc[jr[ui]+nrhsi*ldY]+=vr[ui]*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1; ++ui )
					{
						const CT j = ja[nzi+ui];
						const NT v = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
				}
				else
				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					const auto lnz { ip[nri+1]};
					std::array<NT,nrhs> racc{{}};

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] =       xr[i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz; nzi+=uf )
					{
						std::array<IT,uf> jr;
						std::array<NT,uf> vr;

						for ( auto ui = 0 ; ui < uf ; ++ui )
							jr[ui] = ja[nzi+ui],
							vr[ui] = va[nzi+ui];

						for ( auto ui = 0 ; ui < uf ; ++ui )
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								racc[nrhsi]+=vr[ui]*xc[jr[ui]+nrhsi*ldX],
								yc[jr[ui]+nrhsi*ldY]+=vr[ui]*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi < lnz; nzi++ )
					for ( auto ui = 0 ; ui < 1; ++ui )
					{
						const CT j = ja[nzi+ui];
						const NT v = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i+nrhsi*ldY]+=      racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}

			if ( use_lol == 7 )
			{
				auto const * __restrict const xr = xp + roff;
				auto const * __restrict const xc = xp + coff;
				auto       * __restrict const yr = yp + roff;
				auto       * __restrict const yc = yp + coff;
				auto const * __restrict const va = csr.va_;
				auto const * __restrict const ja = csr.ja_;
				auto const * __restrict const ip = csr.ip_;
				std::array<NT,nrhs> alxt;

				for ( ; nri < eidx ; ++nri )
				{
					const auto i {nri};
					const auto lnz { ip[nri+1]};
					auto nzi { ip[nri]};

					if (nzi==lnz)
						continue;

					const auto axv {alpha * va[nzi]};
					const auto j = ja[nzi];

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i+nrhsi*ldY]+=axv*xc[j+nrhsi*ldX];
					if (roff!=coff || j != i)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yc[j+nrhsi*ldY]+=axv*xr[i+nrhsi*ldX];

					++nzi;
					std::array<NT,nrhs> racc{{}};

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						alxt[nrhsi] = alpha*xr[i+nrhsi*ldX];

					for (                       ; nzi+uf-1 < lnz - 1 ; nzi+=uf )
					{
						std::array<IT,uf> jr;
						for ( auto ui = 0 ; ui < uf ; ++ui )
							jr[ui] = ja[nzi+ui];

						std::array<NT,uf> xu[nrhs];
						for ( auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							for ( auto ui = 0 ; ui < uf ; ++ui )
								xu[nrhsi][ui] = xc[jr[ui]+nrhsi*ldX];

						std::array<NT,uf> vr;
						for ( auto ui = 0 ; ui < uf ; ++ui )
							vr[ui] = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							for ( auto ui = 0 ; ui < uf ; ++ui )
								racc[nrhsi]+=vr[ui]*xu[nrhsi][ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							for ( auto ui = 0 ; ui < uf ; ++ui )
								yc[jr[ui]+nrhsi*ldY]+=vr[ui]*alxt[nrhsi];
					}

					if(with_peel_loop)
					for (                       ; nzi < lnz - 1; nzi++ )
					for ( auto ui = 0 ; ui < 1; ++ui )
					{
						const IT j = ja[nzi+ui];
						const NT v = va[nzi+ui];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX],
							yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					if( nzi < lnz )
					{
						const auto j = ja[nzi];
						const NT v = va[nzi];

						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							racc[nrhsi]+=v*xc[j+nrhsi*ldX];
						if (roff!=coff || j != i)
							for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
								yc[j+nrhsi*ldY]+=v*alxt[nrhsi];
					}

					for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
						yr[i+nrhsi*ldY]+=alpha*racc[nrhsi];
				}
				return RSB_ERR_NO_ERROR;
			}
			return RSB_ERR_INTERNAL_ERROR;
		}
		return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) )
	if ( transA == RSB_TRANSPOSITION_C )
	{
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
					const NT c = rsbpp_conj(v);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
					if(roff+i!=coff+j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[i+roff+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX];
				}
			}
			return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA != RSB_TRANSPOSITION_T )
	{
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const IT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
					const NT c = rsbpp_conj(v);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[roff+i+nrhsi*ldY]+=alpha*v*xp[coff+j+nrhsi*ldX];
					if(roff+i != coff+j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=alpha*c*xp[roff+i+nrhsi*ldX];
				}
			}
			return RSB_ERR_NO_ERROR;
	}

	if ( RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) )
	if ( transA == RSB_TRANSPOSITION_T )
	{
			for ( ; nri < eidx ; ++nri )
			{
				// TODO: need to honour uf here as well
				const auto i {nri};
				auto nzi { csr.ip(nri)};
				const auto lnz { csr.ip(nri+1)};

				for (                       ; nzi      < lnz; nzi+= 1 )
				{
					const CT j = csr.ja(nzi);
					const NT v = csr.va(nzi);
					const NT c = rsbpp_conj(v);

					if(true)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[roff+i+nrhsi*ldY]+=alpha*c*xp[coff+j+nrhsi*ldX];
					if(roff+i != coff+j)
						for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
							yp[j+coff+nrhsi*ldY]+=alpha*v*xp[roff+i+nrhsi*ldX];
				}
			}
			return RSB_ERR_NO_ERROR;
	}
	return RSB_ERR_BADARGS;
}

#if RSBPP_WANT_ALL
template <typename IT, typename NT>
class Coo final
{
	// TODO: Coo -> CooT
	enum {Five=5, Four=4};
	public:
	enum SpmvFlavour : char {SpmvRangeFor='r',SpmvPartialUnrolled='u',SpmvForIterator='i',SpmvForIndex='x',SpmvForIndexU4='U',SpmvForIteratorU='I',SpmvForPtrU4='P',SpmvParForPtr='p',SpmvDefault=SpmvRangeFor,SpmvAll='*',SpmvForIteratorSym='S'};
	enum Ordering : char {OrderingRSB='r',OrderingCOR='u',OrderingCOC='i',OrderingAny='x',OrderingDefault=OrderingRSB};
	static constexpr std::array<enum SpmvFlavour,9> AllSpmvFlavours = {SpmvRangeFor,SpmvPartialUnrolled,SpmvForIterator,SpmvForIndex,SpmvForIndexU4,SpmvForIteratorU,SpmvForPtrU4,SpmvParForPtr,SpmvDefault};
	static constexpr std::array<enum Ordering,Five> AllOrderings = {OrderingRSB,OrderingCOR,OrderingCOC,OrderingAny,OrderingDefault};
	using ij_t = std::pair < IT, IT >;
	using triple_pp_t = std::pair < ij_t, NT >;
	using triple_tp_t = std::tuple < IT, IT, NT >;
	using Quadrant = QuadrantT<IT>;
	using Node = NodeT<IT>;
	using LocalSubmatrix = LocalSubmatrixT<IT>;
	using QuadTree = std::vector<Node>; // TODO: a linearized tree; need more info here; maybe via tuple ..
	static constexpr IT End = std::numeric_limits<IT>::max(); // TODO: maximum ...
	static constexpr IT Begin = 0;

	// FIXME: following class does not work as triple_t
	class triple_tpd_t: public triple_tp_t { public:
		constexpr triple_tpd_t(const IT & i, const IT & j, const NT & v): triple_tp_t ({   i,   j,v}){}
		constexpr triple_tpd_t(const ij_t & ij,            const NT & v): triple_tp_t {ij.first,ij.second,v} {}
		constexpr const IT & i(void)const{return std::get<0>(*this);}
		constexpr const IT & j(void)const{return std::get<1>(*this);}
		constexpr const NT & v(void)const{return std::get<2>(*this);}
		IT & i(void){return std::get<0>(*this); } 
      		IT & j(void){return std::get<1>(*this); } 
		NT & v(void){return std::get<2>(*this); }
	};

	class triple_ppd_t: public triple_pp_t { public:
		constexpr triple_ppd_t(const ij_t & ij, const NT & v): triple_pp_t (ij,v){}
		constexpr triple_ppd_t(void) = default;
		triple_ppd_t(const triple_ppd_t &t) = default;
		constexpr const NT & v(void)const{return this->second      ;}
		constexpr const IT & i(void)const{return this->first.first ;}
		constexpr const IT & j(void)const{return this->first.second;}
		NT & v(void){return this->second      ; }
      		IT & i(void){return this->first.first ; }
		IT & j(void){return this->first.second; }
		std::ostream & print (std::ostream & os) const { os << "i:"<< i() << " j:" << j() << " v:" << v(); return os; }
	};
	//using triple_t = triple_tpd_t; // ok
	using triple_t = triple_ppd_t; // ok
	using trivec_t = std::vector<triple_t,OpenMP_Allocator<triple_t>> ;
	//using trivec_t = std::vector< triple_t, std::allocator<triple_t> > ; // FIXME: broken ?
	//using trivec_t = std::vector< triple_t, CooAllocator <triple_t> > ; // FIXME: broken
	private:
	enum SpmvFlavour sf_{SpmvDefault};
	enum Ordering io_{OrderingDefault}; // TODO: this might transition into leaf nodes
	//enum Ordering lo_{OrderingCOC}; // TODO: new
	//enum Ordering lo_{OrderingCOR}; // TODO: new
	enum Ordering lo_{OrderingAny}; // TODO: new
	trivec_t coo_;
	IT nr_{},nc_{},nnz_{};
	QuadTree qt_{}; // FIXME: temporarily here
	rsb_flags_t flags_{RSB_FLAG_NOFLAGS}; // TODO: move out of this struct

	const static IT base {0}; /* 0 or 1 */
	const static IT ione {1}; /* 0 or 1 */
	const static IT imb {ione - base}; /* 1 or 0 */

#if USE_CXX17
	constexpr static const          NT            none {1};
	constexpr static const          NT            nzero{0.0};
#else
	const          static NT            none ;
	const          static NT            nzero;
#endif

	void sanity_check(void)
	{
		using ct_t = typename trivec_t::size_type;

		assert(nr_ > 0);
		assert(nc_ > 0);
		assert(nnz_>=0);
		assert(coo_.size() >= static_cast<ct_t>(nnz_));
		assert(coo_.size() == static_cast<ct_t>(nnz_));
		assert(checkSpmvFlavour(sf_)==true);
		assert(checkOrdering(io_)==true);
	}
	void erase_qt(void)
	{
		qt_.erase(qt_.begin(),qt_.end());
	}
	void erase_coo(void)
	{
		coo_.erase(coo_.begin(),coo_.end());
	}
	void make_empty(void)
	{
		nnz_=0;
		nr_=0;
		nc_=0;
		erase_coo();
		erase_qt();
	}

	public:
	rsb_coo_idx_t get_ldB(rsb_trans_t transA)const
	{
		return rsb_coo_idx_t (transA == RSB_TRANSPOSITION_N) ? nc() : nr();
	}

	rsb_coo_idx_t get_ldC(rsb_trans_t transA)const
	{
		return rsb_coo_idx_t (transA != RSB_TRANSPOSITION_N) ? nc() : nr();
	}

	public:
	/* explicit */ Coo(Coo && other): // move constructor
		coo_(std::move(other.coo_)), nr_(other.nr_), nc_(other.nc_), nnz_(other.nnz_),qt_(other.qt_), flags_(other.flags_)
	{
		other.make_empty();
		sanity_check();
	}
	/* explicit */ Coo(const Coo & other): // copy constructor (to be redefined is move constructor present)
		coo_(other.coo_), nr_(other.nr_), nc_(other.nc_), nnz_(other.nnz_),qt_(other.qt_), flags_(other.flags_)
	{
		// TODO: make it unnecessary
		sanity_check();
	}
	IT nr(void)const{return nr_;}
	IT roff(void)const{return 0;}
	IT coff(void)const{return 0;}
	IT rows(void)const{return nr();} // TODO: redundant
	IT nc(void)const{return nc_;}
	IT cols(void)const{return nc();} // TODO: redundant
	rsb_flags_t rsbflags(void)const{return flags_;} // TODO: redundant
	IT blocks(void)const{return qt_.size();} // TODO: unfinished
	IT mid(IT of)const{return (of+1)/2;} // middle 
	IT midr(void)const{return mid(nr());} // middle row
	IT midc(void)const{return mid(nc());} // middle column
	IT nnz(void)const{return nnz_;}
	IT ia(IT at)const{ return coo_[at].i();  }
	IT ja(IT at)const{ return coo_[at].j(); }
	NT va(IT at)const{ return coo_[at].v();       }
	bool is_empty(void)const{return nnz()==0;}

	template <typename XC, typename YC, const int UF=4, const int nrhs=1>
	void spmm_partial_unrolled(rsb_trans_t transA, const IT ldX, XC & x, const IT ldY, YC & y, const IT frhsi=0, const NT alpha = none, const IT bidx=Begin, const IT eidx_=End) const
	{
#if 1
		const IT eidx=(eidx_==End)?nnz():eidx_;
		spmm_cooi_partial_unrolled_by_cols<IT, NT, decltype(*this), XC, YC, UF, nrhs>(*this,flags_, transA, ldX, x, ldY, y, frhsi, alpha, bidx, eidx);
#else
		const int uf = UF; /* customizable */
		auto const * __restrict cip = & coo_[bidx];
		const IT eidx=(eidx_==End)?nnz():eidx_;
		auto const * __restrict const clp = & coo_[eidx];
		auto const * __restrict const xp  = & x[(frhsi*ldX)];
		auto       * __restrict const yp  = & y[(frhsi*ldY)];

		if (!RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC  ) )
		if (!RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN  ) )
		if ( transA == RSB_TRANSPOSITION_N )
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=cip[ui].v();

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[i+nrhsi*ldY]+=alpha*v*xp[j+nrhsi*ldX];
		}

		if (!RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC  ) )
		if (!RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN  ) )
		if ( transA == RSB_TRANSPOSITION_T )
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=cip[ui].v();

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[j+nrhsi*ldY]+=alpha*v*xp[i+nrhsi*ldX];
		}

		if (!RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC  ) )
		if (!RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN  ) )
		if ( transA == RSB_TRANSPOSITION_C )
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=rsbpp_conj(cip[ui].v());

			for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
				yp[j+nrhsi*ldY]+=alpha*v*xp[i+nrhsi*ldX];
		}

		if ( transA != RSB_TRANSPOSITION_C )
		if ( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC  ) )
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=cip[ui].v();

			if(i==j)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[i+nrhsi*ldY]+=alpha*v*xp[j+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[i+nrhsi*ldY]+=alpha*v*xp[j+nrhsi*ldX],
					yp[j+nrhsi*ldY]+=alpha*v*xp[i+nrhsi*ldX];
		}

		if ( transA == RSB_TRANSPOSITION_C )
		if ( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC  ) )
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=cip[ui].v();
			const NT c=rsbpp_conj(v);

			if(i==j)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[j+nrhsi*ldY]+=alpha*c*xp[i+nrhsi*ldX];
			else
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[i+nrhsi*ldY]+=alpha*c*xp[j+nrhsi*ldX],
					yp[j+nrhsi*ldY]+=alpha*c*xp[i+nrhsi*ldX];
		}

		if ( transA == RSB_TRANSPOSITION_T )
		if ( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN ) )
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=cip[ui].v();
			const NT c=rsbpp_conj(v);

			if(i==j)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[i+nrhsi*ldY]+=alpha*c*xp[j+nrhsi*ldX];
			else
			{
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[i+nrhsi*ldY]+=alpha*c*xp[j+nrhsi*ldX],
					yp[j+nrhsi*ldY]+=alpha*v*xp[i+nrhsi*ldX];
			}
		}

		if ( transA != RSB_TRANSPOSITION_T )
		if ( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN ) )
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=cip[ui].v();

			if(i==j)
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[i+nrhsi*ldY]+=alpha*v*xp[j+nrhsi*ldX];
			else
			{
				const NT c = rsbpp_conj(v);
				for (auto nrhsi = 0 ; nrhsi < nrhs ; ++nrhsi )
					yp[i+nrhsi*ldY]+=alpha*v*xp[j+nrhsi*ldX],
					yp[j+nrhsi*ldY]+=alpha*c*xp[i+nrhsi*ldX];
			}
		}

		IT ridx = cip-&coo_[0];
		assert( ! ( ridx < eidx && uf == 1 ) );
		if( ridx < eidx )
			spmm_partial_unrolled<XC,YC,1,nrhs>(transA,ldX,x,ldY,y,frhsi,alpha,ridx,eidx); /*the remainder */
#endif
	}

	template <typename C, int UF=4>
	void spmv_partial_unrolled(C const & x, C & y, const NT alpha = none, const IT bidx=Begin, const IT eidx_=End) const
	{
		int const uf = UF; /* customizable */
		auto const * __restrict cip = & coo_[bidx];
		const IT eidx=(eidx_==End)?nnz():eidx_;
		auto const * __restrict const clp = & coo_[eidx];
		auto const * __restrict const xp  = & x[0];
		auto       * __restrict const yp  = & y[0];

		for ( ; cip+(uf-1) < clp ; cip+=uf )
		for (auto ui = 0 ; ui < uf ; ++ui )
		{
			const IT i=cip[ui].i(), j=cip[ui].j();
			const NT v=cip[ui].v(); 
			yp[i]+=alpha*v*xp[j];
		}
		IT ridx = cip-&coo_[0];
		assert( ! ( ridx < eidx && uf == 1 ) );
		if( ridx < eidx )
			spmv_partial_unrolled<C,1>(x,y,alpha,ridx,eidx); /*the remainder */
	}

	template <typename C>
	void spmv_range_for(C const & x, C & y, const NT alpha = none ) const
	{
		for (const auto & coi : coo_)
		{
#if USE_CXX17
			auto [ij,v] = coi; auto [i,j] = ij; //c++17 (structured bindings)
#else /* USE_CXX17 */
			auto i=coi.first.first, j=coi.first.second;
			auto v=coi.second;
#endif /* USE_CXX17 */
			y[i]+=alpha*v*x[j];
		}
	}

	template <typename C>
	void spmv_for_iterator(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin()+bidx, cie = coo_.cbegin()+eidx;
      		for (auto ci = cib ; ci < cie; ++ci)
		{
#if USE_CXX17
			auto [i,j] = ci->first; //c++17 (structured bindings)
#else /* USE_CXX17 */
			const IT i=ci->first.first, j=ci->first.second;
#endif /* USE_CXX17 */
		        const NT v=ci->second;
			y[i]+=alpha*v*x[j];
		}
	}

	template <typename C>
	void spmv_for_iterator_gen_t(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin()+bidx, cie = coo_.cbegin()+eidx;
      		for (auto ci = cib ; ci < cie; ++ci)
		{
#if USE_CXX17
			auto [i,j] = ci->first; //c++17 (structured bindings)
#else /* USE_CXX17 */
			const IT i=ci->first.first, j=ci->first.second;
#endif /* USE_CXX17 */
		        const NT v=ci->second;
			y[j]+=alpha*v*x[i];
		}
	}

	template <typename C>
	void spmv_for_iterator_gen_c(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin()+bidx, cie = coo_.cbegin()+eidx;
      		for (auto ci = cib ; ci < cie; ++ci)
		{
#if USE_CXX17
			auto [i,j] = ci->first; //c++17 (structured bindings)
#else /* USE_CXX17 */
			const IT i=ci->first.first, j=ci->first.second;
#endif /* USE_CXX17 */
		        const NT v=ci->second;
			y[j]+=alpha*rsbpp_conj(v)*x[i];
		}
	}

	template <typename C/*, int symmetry*/>
	void spmv_for_iterator_her(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin()+bidx, cie = coo_.cbegin()+eidx;
      		for (auto ci = cib ; ci < cie; ++ci)
		{
#if USE_CXX17
			auto [i,j] = ci->first; //c++17 (structured bindings)
#else /* USE_CXX17 */
			const IT i=ci->first.first, j=ci->first.second;
#endif /* USE_CXX17 */
		        const NT v=ci->second;
			y[i]+=alpha*v*x[j];
			if(i!=j)
				y[j]+=alpha*rsbpp_conj(v)*x[i];
		}
	}

	template <typename C/*, int symmetry*/>
	void spmv_for_iterator_her_t(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin()+bidx, cie = coo_.cbegin()+eidx;
      		for (auto ci = cib ; ci < cie; ++ci)
		{
#if USE_CXX17
			auto [i,j] = ci->first; //c++17 (structured bindings)
#else /* USE_CXX17 */
			const IT i=ci->first.first, j=ci->first.second;
#endif /* USE_CXX17 */
		        const NT v=ci->second;
			y[i]+=alpha*rsbpp_conj(v)*x[j];
			if(i!=j)
				y[j]+=alpha*v*x[i];
		}
	}

	template <typename C/*, int symmetry*/>
	void spmv_for_iterator_sym_n(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin()+bidx, cie = coo_.cbegin()+eidx;
      		for (auto ci = cib ; ci < cie; ++ci)
		{
#if USE_CXX17
			auto [j,i] = ci->first; //c++17 (structured bindings)
#else /* USE_CXX17 */
			const IT j=ci->first.first, i=ci->first.second;
#endif /* USE_CXX17 */
		        const NT v=ci->second;
			y[i]+=alpha*v*x[j];
			if(i!=j) // lower symmetry
				y[j]+=alpha*v*x[i];
		}
	}

	template <typename C/*, int symmetry*/>
	void spmv_for_iterator_sym_c(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin()+bidx, cie = coo_.cbegin()+eidx;
      		for (auto ci = cib ; ci < cie; ++ci)
		{
#if USE_CXX17
			auto [i,j] = ci->first; //c++17 (structured bindings)
#else /* USE_CXX17 */
			const IT i=ci->first.first, j=ci->first.second;
#endif /* USE_CXX17 */
		        const NT v=ci->second;
			y[j]+=alpha*rsbpp_conj(v)*x[i];
			if(i!=j)
				y[i]+=alpha*rsbpp_conj(v)*x[j];
		}
	}

	template <typename C, IT UF=4>
	void spmv_for_iterator_u(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin(), cie = cib+eidx;
      		auto ci = cib+bidx;

      		for ( ; ci+(UF-1) < cie; ci+=UF)
		{
			for (auto ui = 0 ; ui < UF ; ++ui )
			{
				IT i=ci[ui].first.first, j=ci[ui].first.second; NT v=ci[ui].second;
			       	y[i]+=alpha*v*x[j];
			}
		}
		spmv_for_iterator(x, y, alpha, ci-cib, cie-cib ); /* remainder ... */
	}

	template <typename C>
	void spmv_for_ptr_u4(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		int const uf = 4; /* fixed */
		auto const * __restrict cip = & coo_[bidx];
		auto const * __restrict const clp = & coo_[eidx];
		auto const * __restrict const xp  = & x[0];
		auto       * __restrict const yp  = & y[0];
		for ( ; cip+(uf-1) < clp ; cip+=uf )
		{
			const IT i0=cip[0].first.first, j0=cip[0].first.second; const NT v0=cip[0].second; 
			const IT i1=cip[1].first.first, j1=cip[1].first.second; const NT v1=cip[1].second; 
			const IT i2=cip[2].first.first, j2=cip[2].first.second; const NT v2=cip[2].second; 
			const IT i3=cip[3].first.first, j3=cip[3].first.second; const NT v3=cip[3].second; 
			yp[i0]+=alpha*v0*xp[j0];
			yp[i1]+=alpha*v1*xp[j1];
			yp[i2]+=alpha*v2*xp[j2];
			yp[i3]+=alpha*v3*xp[j3];
		}
		spmv_partial_unrolled<C,1>(x,y,alpha,cip-&(coo_[0]),eidx); /*the remainder */
	}

#if RSBPP_WANT_CPP_THREADS
	template <typename C>
	//void par_spmv_for_ptr(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	void par_spmv_for_ptr(C const & x, C & y, const NT alpha = none) const
	{
		// TODO: use bidx eidx
		// FIXME: this is probably buggy.
		IT nzm;// = mid(nnz_);
		//IT nrm = coo_[nzm].i() + 1;
		// TODO: this mechanism only makes sense if midr() really partitions the matrix.
		IT nrm = midr();
		nzm = std::find_if(std::begin(coo_)/*+nzm*/,std::end(coo_),[&](const triple_t & a){ return a.i() >= nrm; }) - std::begin(coo_);
		IT nz0=0,nz1=nzm,nz2=nnz_;
		std::thread t1 ( [&] { spmv_partial_unrolled<C,4>(x,y,alpha,nz0,nz1); });
		std::thread t2 ( [&] { spmv_partial_unrolled<C,4>(x,y,alpha,nz1,nz2); });
		t1.join();
		t2.join();
	}
#endif

	template <typename C>
	void spmv_for_index_u4(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
		const auto cib = coo_.cbegin(), cie = cib+eidx;
		const IT uf = 4; /* fixed */
      		auto ci = cib+bidx;
      		for ( ; ci+(uf-1) < cie; ci+=uf)
		{
			IT i0=ci[0].first.first, j0=ci[0].first.second; NT v0=ci[0].second; y[i0]+=alpha*v0*x[j0];
			IT i1=ci[1].first.first, j1=ci[1].first.second; NT v1=ci[1].second; y[i1]+=alpha*v1*x[j1];
			IT i2=ci[2].first.first, j2=ci[2].first.second; NT v2=ci[2].second; y[i2]+=alpha*v2*x[j2];
			IT i3=ci[3].first.first, j3=ci[3].first.second; NT v3=ci[3].second; y[i3]+=alpha*v3*x[j3];
		}
		spmv_for_index(x, y, alpha, ci-cib, cie-cib); /* remainder ... */
	}

	template <typename C>
	void spmv_for_index(C const & x, C & y, const NT alpha = none, IT bidx=Begin, IT eidx=End ) const
	{
		if(eidx==End)
			eidx=nnz();
      		for (auto ci = bidx ; ci < eidx; ++ci)
		{
			const IT i=ia(ci), j=ja(ci);
			const NT v=va(ci); 
			y[i]+=alpha*v*x[j];
		}
	}

	template <typename C>
	void beta_scale(C & y, const NT beta = none) const
       	{
		const NT one = none;
		if(beta != one)
      			for (auto & yel : y)
				yel*=beta;
	}

	int get_info(enum rsb_mif_t miflags, void* minfop)const
	//int get_info(int, void*)const
	{
		//rsb_blk_idx_t*nsp=minfop;
		//rsb_flags_t*flagsAp=minfop;
		// ...
		// if(miflags==RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T)
		// 	*(static_cast<rsb_blk_idx_t*>(minfop))=...;
		if(miflags==RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T)
			*(static_cast<rsb_flags_t*>(minfop))=flags_;
		return 0; // FIXME
	}

	template <typename C>
	rsb_err_t spmv(rsb_trans_t transA, const C & x, C & y, const NT alpha = none, const NT beta = none, const IT bidx=Begin, /*const*/ IT eidx=End) const
       	{
		rsb_err_t errval = RSB_ERR_NO_ERROR;
		bool total = ( bidx == Begin && (eidx==End || eidx==nnz()) );
		SpmvFlavour sf=sf_;

		// TODO: can this be invoked on intermediate nodes ?

		if(eidx==End)
			eidx=nnz();

		beta_scale(y, beta); // beta considered to be one, now on

		if(io_==OrderingRSB)
		if(total && qt_.size()>1)
		{
			// FIXME: this branch is experimental (btw, build it with c++17).
			if(g_verbosity>0)
				std::cout << "Recursing SPMV (" << qt_.size() << ")" << std::endl;
#if RSBPP_WANT_CPP_THREADS
			if(omp_get_max_threads()>1)
				rsb_spmx_parallel(transA, qt_, [&](IT qidx) {
					rsb_spmv_node(transA, qt_, qidx, x, y, alpha);
				});
			else
#endif
			{
				rsb_spmv_serial(transA, qt_, x, y, alpha); // FIXME: test this extensively
			}
			return errval;
		}

		if ( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN ) )
		{
			switch(transA)
			{
				case(RSB_TRANSPOSITION_N):
					spmv_for_iterator_her(x, y, alpha, bidx, eidx);
				break;
				case(RSB_TRANSPOSITION_T):
					spmv_for_iterator_her_t(x, y, alpha, bidx, eidx);
				break;
				case(RSB_TRANSPOSITION_C):
					spmv_for_iterator_her(x, y, alpha, bidx, eidx);
				break;
			}
			return errval;
		}

		if ( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC ) )
		{
			switch(transA)
			{
				case(RSB_TRANSPOSITION_N):
				case(RSB_TRANSPOSITION_T):
				       	spmv_for_iterator_sym_n(x, y, alpha, bidx, eidx);
				break;
				case(RSB_TRANSPOSITION_C):
				       	spmv_for_iterator_sym_c(x, y, alpha, bidx, eidx);
				break;
			}
			return errval;
		}

		if(transA != RSB_TRANSPOSITION_N)
			switch(transA)
			{
				case(RSB_TRANSPOSITION_T):
					spmv_for_iterator_gen_t(x, y, alpha, bidx, eidx);
					return errval;
				break;
				case(RSB_TRANSPOSITION_C):
					spmv_for_iterator_gen_c(x, y, alpha, bidx, eidx);
					return errval;
				break;
				// default:
				// need error handling
			}

		if(!total)
		switch(sf)
		{
			case(SpmvRangeFor):
				sf=SpmvForIterator;
			break;
			case(SpmvParForPtr):
				sf=SpmvPartialUnrolled;
			break;
			default:
			;
			// other flavours work on ranges
		};

		if(g_verbosity>2)
			std::cout << "Range SPMV on " << bidx << ".." << eidx << std::endl;
		assert( bidx < eidx );

		switch(sf)
		{
			case(SpmvRangeFor):
			spmv_range_for(x, y, alpha);
			break;
			case(SpmvPartialUnrolled):
			spmv_partial_unrolled(x, y, alpha, bidx, eidx);
			break;
			case(SpmvForIterator):
			spmv_for_iterator(x, y, alpha, bidx, eidx);
			break;
			case(SpmvForIndex):
			spmv_for_index(x, y, alpha, bidx, eidx);
			break;
			case(SpmvForIndexU4):
			spmv_for_index_u4(x, y, alpha, bidx, eidx);
			break;
			case(SpmvForIteratorU):
			spmv_for_iterator_u(x, y, alpha, bidx, eidx);
			break;
			case(SpmvForPtrU4):
			spmv_for_ptr_u4(x, y, alpha, bidx, eidx);
			break;
#if RSBPP_WANT_CPP_THREADS
			case(SpmvParForPtr):
			par_spmv_for_ptr(x, y, alpha);
			break;
#endif
			default:
			std::cerr << "Unknown Spmv type \"" << char(sf_) << "\"/" << int(sf_) << std::endl;
			throw;
		}
		return errval;
	}

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	[[deprecated("only transitory")]]
	rsb_err_t spmv(rsb_trans_t transA, const NT * x, NT * y, const NT*alphap, const NT*betap, const IT bidx=Begin, /*const*/ IT eidx=End) const
       	{
		//TODO: temporary: just for testing (via rsbtest) purposes.
		//FIXME: need no-copy implementation
		const rsb_coo_idx_t ldX = get_ldB(transA);
		const rsb_coo_idx_t ldY = get_ldC(transA);
		rsb_err_t errval = RSB_ERR_NO_ERROR;

		std::vector<NT> xv(ldX),yv(ldY);
		std::copy(y,y+ldY,yv.begin());
		std::copy(x,x+ldX,xv.begin());
		spmv(transA, xv, yv, *alphap, *betap, bidx, eidx);
		std::copy(yv.begin(),yv.end(),y);

		return errval;
	}

	rsb_err_t spmv(rsb_trans_t transA, const NT*alphap, const NT*x, IT incX, const NT*betap, NT*y, IT incY)const
	{
		//TODO: temporary: just for testing purposes.
		//TODO: adapt spmv and beta_scale to handle pointers (or iterable ranges..), then you can get rid of this.
		rsb_err_t errval = RSB_ERR_NO_ERROR;

		if(incX == 1 && incY == 1)
			errval = spmv(transA, x, y, alphap, betap);
		else
		{
			std::cerr << RSB_ERRM_FIFL << " UNSUPPORTED: " <<
				" transA:" << static_cast<char>(transA) <<
				" alphap:" << *alphap <<
				" incX:" << incX <<
				" incY:" << incY <<
			       	std::endl ; 
			// throw;
			errval = RSBPP_ERR_UNSUPPORTED_OPERATION;
		}
		return errval;
	}

#if RSBPP_WANT_CPP20
	rsb_err_t spmv(rsb_trans_t transA, const NT alpha, const std::span<const NT> X, const NT beta, std::span<NT> Y) const
	{
		return spmv(transA, &alpha, X.data(), 1, &beta, Y.data(), 1);
	}
#endif /* RSBPP_WANT_CPP20 */

	rsb_err_t spmm(rsb_trans_t transA, const NT* const alphap, rsb_coo_idx_t nrhs, const NT* const Bp, const NT* const betap, NT* Cp, const IT bidx=Begin, const IT eidx=End)const
	{
		// NOTE: not all layout / symmetry / transposition combinations are parallel
		const rsb_coo_idx_t ldB = get_ldB(transA);
		const rsb_coo_idx_t ldC = get_ldC(transA);
		const bool total = ( bidx == Begin && (eidx==End || eidx==nnz()) );
		constexpr int uf = 4;
		const NT alpha = alphap ? *alphap : none;
		const NT beta = betap ? *betap : none;
		const NT one = none;
		const bool want_spmvs = false;
		rsb_err_t errval = RSB_ERR_NO_ERROR;

		assert(sf_ == SpmvRangeFor);

		if(beta != one)
		{
			const auto vl = ldC * nrhs;
			#pragma omp parallel for if ( ( sizeof(NT) * vl ) / omp_get_num_threads() >= 1024 )
			for (auto i = 0; i < vl; ++i)
				Cp [i] *= beta;
		}

#if RSBPP_WANT_CPP_THREADS
		if(io_==OrderingRSB)
		if(total && qt_.size()>1)
		if(omp_get_max_threads()>1)
		{
			if(g_verbosity>0)
				std::cout << "Recursing SPMM (" << qt_.size() << ")" << std::endl;
			rsb_spmx_parallel(transA, qt_, [&](IT qidx) {
				if(!isLeaf(qt_[qidx]))
					return;
				if(g_verbosity>2)
					std::cout << "Leaf " << qidx << " SPMM on " << qt_[qidx] << std::endl;
				const auto & lsm = std::get<0>(qt_[qidx]);
				spmm(transA, alphap, nrhs, Bp, &one, Cp, lsm.fnz,lsm.fnz+lsm.nnz);
			});
			return errval;
		}
#endif

		if ( ! want_spmvs )
		{
			const IT nt = omp_get_max_threads();

			if(io_!=OrderingRSB)
			if( !RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC  ) )
			if( !RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN  ) )
			if( transA == RSB_TRANSPOSITION_N )
			if( total )
			if(omp_get_max_threads()>1)
			switch(nt)
			{
			case(3):
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
			if(nnz()>=3)
			if(nr()>=3)
			{
				const IT cnz = (nnz()+2)/3;
				const IT fti = coo_[1*cnz].i() >  0  ? coo_[1*cnz].i() : 1;
				const IT sti = coo_[2*cnz].i() > fti ? coo_[2*cnz].i() : fti + 1;
				const auto sidx = std::lower_bound(std::begin(coo_),std::end(coo_),triple_t({fti,0},0) ,[&](const triple_t & a, const triple_t & b){ return a.i() < b.i(); }) - std::begin(coo_);
				const auto tidx = std::lower_bound(std::begin(coo_),std::end(coo_),triple_t({sti,0},0) ,[&](const triple_t & a, const triple_t & b){ return a.i() < b.i(); }) - std::begin(coo_);

				if(sidx > 0 && sidx<nnz())
					assert(coo_[sidx-1].i() < coo_[sidx+0].i());
				if(tidx > 0 && tidx<nnz())
					assert(coo_[tidx-1].i() < coo_[tidx+0].i());
#pragma omp parallel
#pragma omp sections
{
#pragma omp section
				if(sidx> 0)
					spmm(transA, alphap, nrhs, Bp, &one, Cp, 0  , sidx);
#pragma omp section
				if(tidx>sidx)
					spmm(transA, alphap, nrhs, Bp, &one, Cp, sidx, tidx);
#pragma omp section
				if(nnz()>tidx)
					spmm(transA, alphap, nrhs, Bp, &one, Cp, tidx, nnz());
}
				return errval;
			}
			case(2):
			if(nnz()>=2)
			if(nr()>=2)
			{
				const IT mnz = (nnz()+1)/2;
				const IT mid = coo_[mnz].i();
				const auto lidx = std::lower_bound(std::begin(coo_),std::end(coo_),triple_t({mid,0},0)
						,[&](const triple_t & a, const triple_t & b){ return a.i() < b.i(); }
						) - std::begin(coo_);
#pragma omp parallel
#pragma omp sections
{
#pragma omp section
				if(lidx>0)
					spmm(transA, alphap, nrhs, Bp, &one, Cp, 0   , lidx );
#pragma omp section
				if(lidx<nnz())
					spmm(transA, alphap, nrhs, Bp, &one, Cp, lidx, nnz());
}
				return errval;
			}
			default:
			if(nnz()>=nt)
			if(nr()>=nt)
			{
				const IT nt = omp_get_max_threads();
				std::vector<IT> nidxv(nt+1), ridxv(nt+1);

				nidxv[0] = 0, ridxv[0] = 0;
				for(IT n=0;n<nt;++n)
				{
					// approximate next chunk nnz
					const IT cnz = (nnz()-nidxv[n]+(nt-n-1))/(nt-n);
					assert(cnz>=0);
					assert(cnz<=nnz());
					// find row index of balanced next chunk
					ridxv[n+1] = nidxv[n]+cnz < nnz() ? coo_[nidxv[n]+cnz].i() : nr();
					if( ridxv[n+1] <= ridxv[n] )
						ridxv[n+1] = ridxv[n] + 1; // eg a very long row;
					// find first i >= ridxv[n+1] for next rows block
					nidxv[n+1] = std::lower_bound(std::begin(coo_)+nidxv[n],std::end(coo_),triple_t({ridxv[n+1],0},0) ,[&](const triple_t & a, const triple_t & b){ return a.i() < b.i(); }) - std::begin(coo_);

					// std::cout << " i: " << coo_[nidxv[n]].i() << "--" << coo_[nidxv[n+1]].i() ;
					// std::cout << " n: " << nidxv[n] << "--" << nidxv[n+1] ;
					// std::cout << std::endl;

					assert( nidxv[n+1] == nnz() || nidxv[n] < nidxv[n+1] );
					if( nidxv[n+1] != nnz() && ( nidxv[n+1] > 0 ) )
						assert( coo_[nidxv[n+1]-1].i() < coo_[nidxv[n+1]].i() );
				}
#pragma omp parallel
{
				const IT th_id = omp_get_thread_num();
				if( nidxv[th_id] < nidxv[th_id+1] )
					spmm(transA, alphap, nrhs, Bp, &one, Cp, nidxv[th_id], nidxv[th_id+1]);
}
#pragma omp barrier
				return errval;
			}
			}

			for(rsb_coo_idx_t nrhsi=0;nrhsi<nrhs;)
			{
				using ThisMT = Coo<IT,NT>;
				const auto spmm_kf = std::make_tuple(
					&ThisMT::spmm_partial_unrolled<decltype(Bp),decltype(Cp),uf,0>,
					&ThisMT::spmm_partial_unrolled<decltype(Bp),decltype(Cp),uf,1>,
#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
					&ThisMT::spmm_partial_unrolled<decltype(Bp),decltype(Cp),uf,2>,
					&ThisMT::spmm_partial_unrolled<decltype(Bp),decltype(Cp),uf,3>,
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */
					&ThisMT::spmm_partial_unrolled<decltype(Bp),decltype(Cp),uf,RSBPP_PRESET_SPMM_BLOCKINGS>
				);

				switch(nrhs-nrhsi)
				{
#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
					case(1):
					{
						constexpr auto nrhsb=1;
						(this->*std::get<nrhsb>(spmm_kf))(transA,ldB,Bp,ldC,Cp,nrhsi,alpha,bidx,eidx);
						nrhsi+=nrhsb;
						break;
					}
					case(2):
					{
						constexpr auto nrhsb=2;
						(this->*std::get<nrhsb>(spmm_kf))(transA,ldB,Bp,ldC,Cp,nrhsi,alpha,bidx,eidx);
						nrhsi+=nrhsb;
						break;
					}
					case(3):
					{
						constexpr auto nrhsb=3;
						(this->*std::get<nrhsb>(spmm_kf))(transA,ldB,Bp,ldC,Cp,nrhsi,alpha,bidx,eidx);
						nrhsi+=nrhsb;
						break;
					}
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */
					default:
					{
						constexpr auto nrhsb=RSBPP_PRESET_SPMM_BLOCKINGS;
						(this->*std::get<nrhsb>(spmm_kf))(transA,ldB,Bp,ldC,Cp,nrhsi,alpha,bidx,eidx);
						nrhsi+=nrhsb;
						break;
					}
				}
			}
			return errval;
		}

		assert( want_spmvs );

#if 0
		for(rsb_coo_idx_t nrhsi = 0; nrhsi < nrhs; ++nrhsi)
			spmv(transA, Bp+nrhsi*ldB, Cp+nrhsi*ldC, alphap, &one, bidx, eidx);
#else
		for(rsb_coo_idx_t nrhsi = 0; nrhsi < nrhs; ++nrhsi)
			spmm_partial_unrolled<decltype(Bp),decltype(Cp),uf,1>(transA,ldB,Bp,ldC,Cp,nrhsi,alpha,bidx,eidx);
#endif
		return errval;
	}

	rsb_err_t spmm(rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, rsb_flags_t order, const NT* const Bp, const NT beta, NT* Cp, const IT bidx=Begin, const IT eidx=End)const
	{
		assert(order == RSB_FLAG_WANT_COLUMN_MAJOR_ORDER);
		return spmm(transA, &alpha, nrhs, /* order,*/ Bp, &beta, Cp, bidx, eidx);
	}

	rsb_err_t spmm(const rsb_trans_t transA, const NT alpha, rsb_coo_idx_t nrhs, const NT* const Bp, const NT beta, NT* Cp)const
	{
		return spmm(transA, &alpha, nrhs, Bp, &beta, Cp);
	}

#if RSBPP_WANT_CPP20
	rsb_err_t spmm(rsb_trans_t transA, const NT alpha, IT nrhs, rsb_flags_t order, const std::span<const NT> x, const NT beta, std::span<NT> y) const
	{
		return spmm(transA, alpha, nrhs, order, x.data(), beta, y.data());
	}
#endif /* RSBPP_WANT_CPP20 */


	template <typename C>
	void spmm(const rsb_trans_t transA, const C & x, C & y, rsb_coo_idx_t nrhs, const NT alpha = none, const NT beta = none) const
	{
		spmm(transA, &alpha, nrhs, x.data(), &beta, y.data());
	}

	template< class Function, class... Args >
	void rsb_postorder(const QuadTree & qt, IT qidx, Function&& f, Args&&... args /*,...*/) const
	{
		// TODO: unfinished
		if(isLeaf(qt[qidx]))
		{
			if(g_verbosity>2)
				std::cout << "L:" << qidx << ":" << qt[qidx] << std::endl;
			f(qidx);
		}
		else
		{
			for (const auto midx : {0,1,2,3}) // ... : qindices
			{
				assert(std::get<1>(qt[qidx])[midx] != qidx);
				if( std::get<1>(qt[qidx])[midx +1] > std::get<1>(qt[qidx])[midx] )
					rsb_postorder(qt, std::get<1>(qt[qidx])[midx], f, args ... );
			}
			//std::cout << "N:" << qidx << ":" << qt[qidx] << std::endl;
			f(qidx);
		}
	}

	int rsb_num_threads(void)const
	{
		int ent=rsbpp_getenv_int_t("RSB_NUM_THREADS", 1);
		ent=std::min(std::max(ent,RSB_MIN_THREADS),RSB_MAX_THREADS);
		return ent;
	}

#if RSBPP_WANT_CPP_THREADS
	template <class Function, class... Args >
	void rsb_spmx_parallel(const rsb_trans_t transA, const QuadTree & qt, Function && spmx_f, Args&&... args /*,...*/) const
	{
		// TODO: unfinished
		// TODO: trace infrastructure
		const int nq=qt.size();
		const int ent=rsb_num_threads();
		const int nt=std::min(ent,nq); // number of threads // TODO: shall not exceed leaves count
		BitVector psb(nq,false); // processed submatrices indices
		BitVector asi(nq,false); // // active submatrices indices
		std::mutex ss_mutex; // shared structures mutex
		const bool vrbs=(g_verbosity>3);
		const bool any_sym =
			( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_SYMMETRIC  ) ) ||
			( RSB_DO_FLAG_HAS( flags_, RSB_FLAG_HERMITIAN  ) );
		//IT qidx=0;
		//IT tp=0;
		auto rtl = ( [&] 
		(IT tidx)
		{
#if 0
			IT lidx=0;
			for(lidx=0;lidx<qt.size();++lidx)
				rsb_spmv_node(transA, qt, lidx, x, y, alpha);
#else
			IT lidx{0};
			bool done = false;
			//std::string mb("RSB thread "+('0'+tidx));
			std::string mb("RSB thread ");
			if(vrbs)
				std::cout << mb << " BEGIN" << std::endl;
		while( !done )
		{
			IT dc {0};
			{
				std::lock_guard<std::mutex> guard(ss_mutex);
				// ss_mutex.lock();
				for(lidx=0;lidx<nq;++lidx)
				{
					if(psb[lidx]==false) // not processed yet
					if(asi[lidx]==false) // no thread processing it
					{
						psb[lidx]=true; // no candidate anymore
						if(isLeaf(qt[lidx]))
						{
							bool is_free = true;
							for(auto kidx=0;kidx<nq;++kidx)
							if(asi[kidx]==true) // if busy
							if(kidx!=lidx)
							if(isLeaf(qt[kidx]))
							{
								if( any_sym || ( transA == RSB_TRANSPOSITION_N ) )
								if(std::get<0>(qt[lidx]).intersects_on_rows(std::get<0>(qt[kidx])))
								{
									is_free = false;
									break;
								}
								if( any_sym || ( transA != RSB_TRANSPOSITION_N ) )
								if(std::get<0>(qt[lidx]).intersects_on_cols(std::get<0>(qt[kidx])))
								{
									is_free = false;
									break;
								}
							}
							if(is_free)
							{
								psb[lidx]=true; // no candidate anymore
								asi[lidx]=true; // declare lidx busy
								break;
							}
							else
								psb[lidx]=false; // for later
						}
					}
				}
				if((dc=std::count(psb.begin(),psb.end(),true))==nq)
					done=true;
				// for(lidx=0;lidx<qt.size();++lidx) std::cout << psb[lidx]; std::cout << std::endl;
				// if find a free matrix whose range is free then lock
				// ss_mutex.unlock();
			}

			if( lidx<nq )
			{
				spmx_f(lidx,args...);
				{
					std::lock_guard<std::mutex> guard(ss_mutex);
					asi[lidx]=false;
				}
			}

			if(vrbs) if(  lidx<nq )
				std::cout << mb << tidx<< " PROCESS " << lidx << " " << dc << "/" << nq <<  std::endl;
			if(vrbs) if(!(lidx<nq ))
				std::cout << mb << " WAIT    " << lidx << " " << dc << "/" << nq <<  std::endl;
		}
			if(vrbs)
				std::cout << mb << "END" << std::endl;
#endif
		});
		if(vrbs)
			std::cout << "RSB SPMV BEGIN (" << nt << " threads)" << std::endl;
		//std::thread rt1 (rtl,1);              rt1.join();
		//std::thread rt1 (rtl,1), rt2 (rtl,2); rt1.join(); rt2.join();
		std::vector<std::thread*> tv;
		for(auto ti=0;ti<nt;++ti) tv.push_back(new std::thread(rtl,ti));
		for(auto ti=0;ti<nt;++ti) if(tv[ti]) tv[ti]->join();
		for(auto ti=0;ti<nt;++ti) if(tv[ti]) delete tv[ti];
		if(g_verbosity>2)
			std::cout << "RSB SPMV END" << std::endl;
	}
#endif

	void rsb_postorder_print(const QuadTree & qt) const
	{
		rsb_postorder(qt, 0, [&](IT qidx)
		{ 
			if(isLeaf(qt[qidx])) std::cout << "L"; else std::cout << "N";
			std::cout << ":" << qidx << ":" << qt[qidx] << std::endl;return;
		} /*,...*/);
	}

	template <typename C>
	void rsb_spmv_node(const rsb_trans_t transA, const QuadTree & qt, IT qidx, const C & x, C & y, const NT alpha) const
       	{
		if(!isLeaf(qt[qidx]))
			return;

		if(g_verbosity>2)
			std::cout << "Leaf " << qidx << " SPMV on " << qt[qidx] << std::endl;

		const auto & lsm = std::get<0>(qt[qidx]);
		spmv(transA,x,y,alpha,none,lsm.fnz,lsm.fnz+lsm.nnz);
	}

	template <typename C>
	void rsb_spmv_serial(const rsb_trans_t transA, const QuadTree & qt, const C & x, C & y, const NT alpha = none, const NT beta = none) const
       	{
		beta_scale(y, beta);

		rsb_postorder(qt, 0, [&](IT qidx) { rsb_spmv_node(transA, qt, qidx, x, y, alpha); });
	}
#if 0
	bool minitest()
	{
		triple_t t {{1,2},3.0};
		std::cout << t.first.first << " " << t.first.second << " " << t.second << std::endl;
		//trivec_t = std::vector< triple_t, CooAllocator <triple_t> > ;
	}
#endif

	std::ostream & print_qt (std::ostream & os)const
	{
		os << "QuadTree [" << qt_.size() << "]:\n" ;
		for (auto & q : qt_)
			os << " " << ((&q)-&qt_[0]) << " [" << q << "]" << std::endl;
		//os << std::endl;
		return os;
	}

	std::ostream & print (std::ostream & os)const
	{
		const auto cib = coo_.cbegin();

		os	<< STRINGIFY1 (nr_) << "  " 
			<< STRINGIFY1 (nc_) << "  " 
			<< STRINGIFY1 (nnz_) << std::endl;

      		if (coo_.cbegin() == coo_.cend() )
			os << "<empty>";
      		for (auto coi = coo_.cbegin(); coi < coo_.cend() ; ++coi )
			//os << (coi-cib) << "  " << "(*coi)" << std::endl;
			//os << (coi-cib) << "  " << (*coi) << std::endl;
			os << (coi-cib) << "  ", coi->print(os) << std::endl;
		return os;
	}
	explicit Coo(const rsb_flags_t flags=RSB_FLAG_NOFLAGS):flags_(flags){}// FIXME: what about index ordering and spmv flavour ? FIXME: move flags_ out of here.
	explicit Coo(std::string filename)
	{
		std::ifstream infile(filename.c_str());
		std::string line;
		std::getline(infile,line);
		while(line.size() > 0 && line[0] == '%')
			// std::cout << " " << line << std::endl,
			std::getline(infile,line);
		std::istringstream (line) >> this->nr_ >> this->nc_ >> this->nnz_;
		coo_.resize(nnz_);
		if ( ! infile.is_open())
			return; // TODO
		for ( IT nzi = 0; nzi < nnz_ && !infile.eof() ; ++nzi )
		{
			IT i,j; NT v;
			infile >> i >> j >> v;
			nr_=std::max(nr_,i);
			nc_=std::max(nc_,j);
			triple_t triple = {{i-ione,j-ione},v};
			coo_[nzi]=triple;
		}
	}
	Coo(const std::vector<NT>&VA,const std::vector<IT>&IA,const std::vector<IT>&JA,rsb_flags_t flags = RSB_FLAG_NOFLAGS,const Ordering ordering = OrderingDefault)
 	{
		// TODO: overall state consistency check!
		init_coo(VA,IA,JA,flags,ordering);
	}

	private:
	bool checkSpmvFlavour(enum SpmvFlavour sf)const
	{
		const auto af = AllSpmvFlavours;
		if(std::find(af.begin(),af.end(),sf)!=af.end())
			return true;
		return false;
	}
	public:
	void setSpmvFlavour(enum SpmvFlavour sf)
	{
		if(!checkSpmvFlavour(sf))
			throw;
		sf_=sf;
		return;
	}
	private:
	bool checkOrdering(enum Ordering io)const
	{
		const auto ao = AllOrderings;
		if(std::find(ao.begin(),ao.end(),io)!=ao.end())
			return true;
		return false;
	}
	public:
	void setOrdering(enum Ordering io)
	{
		// TODO: soon back private..
		if(!checkOrdering(io))
			throw;
		io_=io;
		return;
	}
	void init_coo(trivec_t && coo, rsb_flags_t sym_flags, enum Ordering ordering)
	{
		coo_=coo;
		nr_=nc_=0;
		for(auto t : coo_)
			nr_ = std::max(nr_,t.i()),
			nc_ = std::max(nc_,t.j());
		nr_++;
		nc_++;
		nnz_=coo_.size();
		flags_=sym_flags;
		setOrdering(ordering);
	}
	void init_coo(const std::vector<NT>&VA,const std::vector<IT>&IA,const std::vector<IT>&JA,rsb_flags_t flags = RSB_FLAG_NOFLAGS,const Ordering ordering = OrderingDefault)
	{
		// TODO: overall state consistency check!
		nnz_=VA.size();
		coo_.clear();
		coo_.reserve(nnz_);
		flags_=flags;
		nr_=nc_=0;
		IT i,j;
		NT v;

		for(auto k=0; k < nnz_; ++k)
			v = VA[k],
			i = IA[k],
			j = JA[k],
			coo_.emplace_back( triple_t {{i,j},v} ),
			nr_ = std::max(nr_,i),
			nc_ = std::max(nc_,j);
		nr_++;
		nc_++;
		if(ordering != OrderingAny)
		{
			io_ = ordering;
			this->sort();
		}
		assert(coo_.begin()+nnz_==coo_.end());
	}
	void push_back(IT i, IT j, NT v)
	{
#if 1
		triple_ppd_t ijv {{i,j},v};
		coo_.emplace_back(std::move(ijv)); // implements perfect forwarding (in different situations spares duplicated construction)
#else
		coo_.push_back({{i,j},v});
#endif
		i+=imb;
		nr_=std::max(nr_,i);
		j+=imb;
		nc_=std::max(nc_,j);
		++nnz_;
	}
	void sort(void)
	{
		// TODO: say in future we will use sorting-insertors, this will be useless
		switch(io_)
		{
			case(OrderingRSB):
				zort();
			break;
			case(OrderingCOR):
				sort_cor();
			break;
			case(OrderingCOC):
				sort_coc();
			break;
			case(OrderingAny):
				;// bogus
			break;
			default:
				throw;
		}
	}
	void sort_cor(IT bidx, IT eidx)
	{
		if(eidx==End)
			eidx=nnz();

		std::sort(coo_.begin()+bidx, coo_.begin()+eidx,
			[](triple_t a, triple_t b){ return ( ( a.i() < b.i() || ( a.i() == b.i() && a.j() < b.j()) ) ? true : false ); }
			);
	}
	void sort_cor()
	{
		io_ = OrderingCOR;

		sort_cor(Begin,End);
		erase_qt();
	}
	bool sort_cor_check(void) const
	{
		return std::is_sorted(coo_.begin(), coo_.end(),
			[](triple_t a, triple_t b) { return b.i() > a.i() || ( b.i() == a.i() && b.j() > a.j()) ; }
			);
	}
	bool sort_coc_check(void) const
	{
		return std::is_sorted(coo_.begin(), coo_.end(),
			[](triple_t a, triple_t b) { return b.j() > a.j() || ( b.j() == a.j() && b.i() > a.i()) ;   }
			);
	}
	void sort_coc(IT bidx, IT eidx)
	{
		if(eidx==End)
			eidx=nnz();

		std::sort(coo_.begin()+bidx, coo_.begin()+eidx,
			[](triple_t a, triple_t b){ return ( ( a.j() < b.j() || ( a.j() == b.j() && a.i() < b.i()) ) ? true : false ); }
		);
	}
	void sort_coc()
	{
		io_ = OrderingCOC;

		sort_coc(Begin,End);
		erase_qt();
	}
	void zort(typename trivec_t::iterator bi, const IT lfnz, const IT lnnz, const IT lfr, const IT lnr, const IT lfc, const IT lnc, QuadTree & qt)
	{
		// TODO: rename to zpart.
		// TODO: make return type template 
		// TODO: shall not take trivec_t iterators but only offsets, and then generic triple structure
		// need a max heap structure here.
		const auto midrv = mid(lnr) + lfr;
		const auto midcv = mid(lnc) + lfc;
		const auto qidx=static_cast<IT>(qt.size());
	       	const Quadrant smq {qidx,qidx,qidx,qidx,qidx};// TODO: static_const // FIXME: big type
		const size_t cbs = g_cbs;
		const IT minnr{Four},minnc{Four};
		//const LocalSubmatrix cbk {lnnz,lnr,lnc};
		const LocalSubmatrix lbk {0,static_cast<IT>(cbs/sizeof(triple_t)),0,minnr,0,minnc};
		//const LocalSubmatrix lbk {     4,    4,    4};
		const LocalSubmatrix lsm {lfnz,lnnz,lfr,lnr,lfc,lnc};

		if(g_verbosity>0)
		{
			std::cout << "lfnz:" << lfnz << " " ;
			std::cout << "lnnz:" << lnnz << " " ;
			std::cout << "lfr:" << lfr << " " << "lfc:" << lfc << " ";
			std::cout << "lnr:" << lnr << " " << "lnc:" << lnc << " ";
			std::cout << "midrv:" << midrv << " " << "midcv:" << midcv << std::endl;
		}

		qt.push_back(Node {lsm,smq});

		const auto fi = bi+lfnz;
		const auto ei = bi+lfnz+lnnz;

		if( lsm < lbk )
		{
			// accept quadrant as it is
			if(g_verbosity>0)
				std::cout << "Small Leaf Quadrant " << qidx << " " << qt.back() << "\n";

			if(g_verbosity>1)
			{
				std::cout << "Best formats: ";
				if( lsm.fits_short_indices() )
					std::cout << " H";
				if( lsm.fits_csr() )
					std::cout << "-CSR";
				else
					std::cout << "-COO";
				std::cout << "\n";
			}

			if(lo_!=OrderingAny)
			{
				if(lo_==OrderingCOR)
					sort_cor(fi-bi,ei-bi);
				if(lo_==OrderingCOC)
					sort_coc(fi-bi,ei-bi);
			}
			return;
		}
		
		assert(lnr  >0);
		assert(lnc  >0);
		assert(lnnz >0);
		assert(lfnz>=0);
		assert(lfr >=0);
		assert(lfc >=0);
		assert(midrv > lfr);
		assert(midcv > lfc);

		// TODO: partitioning shall be able to e.g. encompass paddings or alignments
		const auto li = std::stable_partition(bi+lfnz, ei, [&](triple_t a){return a.i() < midrv ;});
		const auto uj = std::stable_partition(bi+lfnz, li, [&](triple_t a){return a.j() < midcv ;});
		const auto lj = std::stable_partition(li     , ei, [&](triple_t a){return a.j() < midcv ;});
		const IT nrl=(lfr+lnr)-midrv,nru=midrv-lfr;
		const IT ncl=midcv-lfc,ncr=(lfc+lnc)-midcv;
		assert(nrl >0);
		assert(nru >0);
		assert(ncl >0);
		assert(ncr >0);
		assert(static_cast<IT>(uj>=bi));
		assert(static_cast<IT>(li>=bi));
		assert(static_cast<IT>(lj>=bi));
		assert(static_cast<IT>(ei> bi));
		const Quadrant qnz {lfnz,static_cast<IT>(uj-bi),static_cast<IT>(li-bi),static_cast<IT>(lj-bi),static_cast<IT>(ei-bi)};
		// the static_cast's here are to avoid warnings stemming from promotion of e.g. short int to int during sum.
		const Quadrant qro {lfr,lfr,lfr,static_cast<IT>(lfr+nru),static_cast<IT>(lfr+nru)}; // quadrant's rows offset
		const Quadrant qco {lfc,lfc,static_cast<IT>(lfc+ncl),lfc,static_cast<IT>(lfc+ncl)}; // quadrant's cols offset
		const Quadrant qnr {lnr,nru,nru,nrl,nrl}; // quadrant's number of rows
		const Quadrant qnc {lnc,ncl,ncr,ncl,ncr}; // quadrant's number of cols
		assert(ei>bi);
		std::get<1>(qt[qidx])[0]=qt.size(); // only leaves can point to themselves
		// TODO: make parallel tasks out of these ?!
		// TODO: need further decision e.g. revert local recursion when only one leaf is enough
		for (const auto midx : {1,2,3,4}) // ... : qindices
		{
			if(qnz.n_at(midx) > 0)
				std::get<1>(qt[qidx])[midx-1]=qt.size(),
				zort(bi,qnz[midx-1],qnz.n_at(midx),qro[midx],qnr[midx],qco[midx],qnc[midx],qt);
			std::get<1>(qt[qidx])[midx]=qt.size();
		}
		for (const auto midx : {0,1,2,3,4}) // ... : qindices
			assert(std::get<1>(qt[qidx])[midx] != qidx);
		if(g_verbosity>3)
			std::cout << "     Quadrant " << qidx << " " << qt[qidx] << "\n";
		return;
	} /* zort */
	Ordering getOrdering(void)const
	{
		return io_;
	}
	void zort(void)
	{
		erase_qt();
		zort(coo_.begin(),0,nnz(),0,nr(),0,nc(),qt_);
		if(g_verbosity>2)
			std::cout << "Zorted. Resulting quadtree:" << std::endl,
			print_qt(std::cout),
			std::cout << "." << std::endl;
		io_ = OrderingRSB;
	}
#define LIN_TRANS(X,TIMES,PLUS) (((X)*TIMES) + PLUS)
	void lin_trans(NT plus, NT times)
	{
		// TODO: need test case
		auto lin_trans_lambda = [&](triple_t & triple){ triple_t rt{{triple.i(),triple.j()},LIN_TRANS(triple.v(),times,plus) }; triple = rt; };
		std::for_each(coo_.begin(),coo_.end(),lin_trans_lambda );
	}

	void lin_trans_idx(IT plus, IT times)
	{
		// TODO: need test case
		auto lin_trans_lambda = [&](triple_t & triple){ 
			triple_t rt{{LIN_TRANS(triple.i(),times,plus),LIN_TRANS(triple.j(),times,plus)},triple.v()}; triple = rt; };
		std::for_each(coo_.begin(),coo_.end(),lin_trans_lambda );
	}

#if 0
	triple_t & at(IT idx)
	{
		return coo_.at(idx); // at has bounds checking 
	}
	const triple_t & at(IT idx)const
	{
		return coo_.at(idx); // at has bounds checking 
	}
	triple_t & operator[](IT idx)
	{
		return coo_[idx];
	}
	const triple_t & operator[](IT idx)const 
	{
		return coo_[idx];
	}
#else
	NT & at(IT idx)
	{
		return coo_.at(idx).v(); // at has bounds checking 
	}
	const NT & at(IT idx)const
	{
		return coo_.at(idx).v(); // at has bounds checking 
	}
	NT & operator[](IT idx)
	{
		return coo_[idx].v();
	}
	const NT & operator[](IT idx)const 
	{
		return coo_[idx].v();
	}
#endif

	void selftests(IT expected_nnz = -1, enum SpmvFlavour sf=SpmvAll) 
	{
#if USE_CXX17
		static_assert(!std::is_aggregate<IT>::value,"");
#endif /* USE_CXX17 */
		static_assert(std::is_pod<IT>::value,"!"); // no e.g. std::complex
		static_assert(std::is_arithmetic<IT>::value,"!");
		if( sf == SpmvAll )
		{
			const auto af = AllSpmvFlavours;
			for ( auto sfa : af)
			{
				setSpmvFlavour(sfa);
				std::cout << "Now will use spmv algorithm " << char(sfa) << std::endl;
				const_selftest(expected_nnz);
				selftest();
			}
		}
		else
		{
			std::cout << "Will use only spmv algorithm " << char(sf) << std::endl;
			setSpmvFlavour(sf);
			const_selftest(expected_nnz);
			selftest();
		}
	}

	void selftests(IT expected_nnz, char sfc)
	{
		const auto sf=(enum Coo<IT,NT>::SpmvFlavour)sfc;// FIXME
		selftests(expected_nnz, sf);
	}

	void selftest(void) 
	{
		if(!is_empty())
		{
			auto & A = *this;
			auto fel = A[0];

			A[0]*=2;

			if( (fel+fel) != A[0] )
			       	throw;
		}
	}

	void spmv_selftest(void) const
	{
		const auto & A = *this;
		rsb_trans_t transA{RSB_TRANSPOSITION_N};

		if ( true /* this test can make sense if all values are near integers ... */ )
		{
			constexpr NT zero{};
			constexpr NT none {1};
			const auto azero=std::abs(zero);
			const auto eps=std::numeric_limits<decltype(azero)>::epsilon();
			const auto tol=(eps*nr())*100;

			// TODO: valarray, vector, ...
			std::valarray<NT> x(nc());
			std::valarray<NT> y(nr());

			x=none;
			// std::cout << "Type: " << '?' << "  epsilon: " << eps << std::endl;
			if(! (x.sum() == NT(nc())) )
				throw;
			y=0;
			if(y.sum() != zero)
				std::cout << "check sum is not zero: " << y.sum() << std::endl,
				throw;
			A.spmv(transA,x,y,+2,+1); // y' = 1 * y  + 2 * A * x = 2 * A * x;
			if(y.sum() == zero)
				std::cout << "check sum is zero where it should not! " << std::endl,
				throw;
			A.spmv(transA,x,y,+3,+2); // y" = 2 * y' + 3 * A * x;
			// y" = 4 * A * x + 3 * A * x = 7 * A * x;
			A.spmv(transA,x,y,+7,-1); // y shall become 0
			
			if(y.sum() != zero)
			{
				std::cout << "These ideally shall be zero but are not:" << std::endl;
				for(int i=0;i<nr();++i)
					//if(std::abs(y[i])!=zero)
					if(std::abs(y[i]) > std::abs(tol))
						std::cout << i << " : " << y[i] << std::endl;
				std::cout << std::endl;
			}
			//if(y.sum() != zero)
			if(std::abs(y.sum()) > std::abs(tol))
				std::cout << "check sum is too much: " << y.sum() << " wrt tolerance " << tol << std::endl,
				//std::cout << "check sum is not zero: " << y.sum() << std::endl,
				throw;
		}
	}

	void const_selftest(IT expected_nnz = -1) const
	{
		const auto & A = *this;

		if ( expected_nnz >  0 )
			if(A.is_empty())
			       	throw;
		if ( expected_nnz >= 0 )
			if(A.nnz()!=expected_nnz )
			       	throw;
		if(!is_empty())
		{
			auto C = *this; // copy
			if( this->is_empty())
			       	throw;
			if( C.is_empty())
			       	throw;
			auto B = std::move(C);
			if(!C.is_empty())
			       	throw;
			if( B.is_empty())
			       	throw;
			spmv_selftest();
		}
	}
	Coo & operator += (const Coo & other)
	{
		trivec_t c_coo(coo_);
		triple_t t{other.coo_[0]};
		c_coo.resize(c_coo.size()+other.coo_.size(),t);// TODO: need triple_t ctor
		std::copy_n(other.coo_.begin(),other.coo_.size(),c_coo.begin()+nnz_);
		assert(nr_==other.nr_);
		assert(nc_==other.nc_);
		//nr_=std::max(nr_,other.nr_);
		//nc_=std::max(nc_,other.nc_);
		// TODO: flags check ...
		// ----- 
		coo_=std::move(c_coo);
		nnz_=coo_.size();
		this->sort(); // shall go from member function to library function
		// TODO: incomplete: need duplicate/zeros removal/collapse.
		return *this;
	}
	Coo & operator = (const Coo & other)
	{
		coo_=other.coo_; nr_=other.nr_; nc_=other.nc_; nnz_=other.nnz_; qt_=other.qt_; flags_=other.flags_;
		return *this;
	}
}; /* Coo */

#if USE_CXX17
#else
template <typename IT, typename NT>
const          NT            Coo<IT,NT>::none {1} ;
template <typename IT, typename NT>
const          NT            Coo<IT,NT>::nzero {0.0} ;
#endif

template <typename NT>
using RsbPP_Matrix_T = Coo<rsb_coo_idx_t,NT>;

template <typename IT,typename NT>
using Csr = Coo<IT,NT>;

#endif /* RSBPP_WANT_ALL */

template<typename NT, typename IT, typename CT = IT>
rsb_err_t rsbpp_coo_spmx(rsb_flags_t flags, const IT nnz, const IT nr, const IT nc, const NT* VA, const CT* IA, const CT* JA, const IT nrhs, const IT ldX, const NT* rhs, const IT ldY, NT* out, const NT* alphap, IT incx, IT incy, const rsb_trans_t transA, IT roff = 0, IT coff = 0, bool by_rows = false)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const NT alpha = alphap ? *alphap : 1;

	if (( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) ||
	      RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) ) &&
		( nr != nc && roff == coff ) )
		return RSB_ERR_BADARGS;

	struct CooP {
		const IT nnz_;
		const IT nr_;
		const IT nc_;
		const NT* va_;
		const CT* ia_, * ja_;
		const IT roff_, coff_;

		inline IT roff(void)const{ return roff_; }
		inline IT coff(void)const{ return coff_; }
		inline IT nr(void)const{ return nr_; }
		inline IT nc(void)const{ return nc_; }
		inline CT ia(IT at)const{ return ia_[at]; }
		inline CT ja(IT at)const{ return ja_[at]; }
		inline NT va(IT at)const{ return va_[at]; }
	};

	const CooP coop {nnz,nr,nc,VA,IA,JA,roff,coff};
	const IT bidx{0};
	const IT eidx{nnz};

	if( by_rows )
	{
		const auto ldY = incy;
		const auto ldX = incx;
		constexpr int uf = 4;
		using RT = decltype(rhs);
		using OT = decltype(out);

		rsb_coo_idx_t nrhsi=0;

#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
		for(;nrhsi+32-1<nrhs;)
				errval |= spmm_coo_partial_unrolled_by_rows<IT,CT,NT,decltype(coop),RT,OT,uf,32>
					(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 32;

		for(;nrhsi+16-1<nrhs;)
				errval |= spmm_coo_partial_unrolled_by_rows<IT,CT,NT,decltype(coop),RT,OT,uf,16>
					(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 16;

		for(;nrhsi+8-1<nrhs;)
				errval |= spmm_coo_partial_unrolled_by_rows<IT,CT,NT,decltype(coop),RT,OT,uf,8>
					(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 8;

		for(;nrhsi+4-1<nrhs;)
				errval |= spmm_coo_partial_unrolled_by_rows<IT,CT,NT,decltype(coop),RT,OT,uf,4>
					(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 4;

		for(;nrhsi+2-1<nrhs;)
				errval |= spmm_coo_partial_unrolled_by_rows<IT,CT,NT,decltype(coop),RT,OT,uf,2>
					(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 2;
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */

		for(;nrhsi+1-1<nrhs;)
				errval |= spmm_coo_partial_unrolled_by_rows<IT,CT,NT,decltype(coop),RT,OT,uf,1>
					(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 1;
		return errval;
	}


	for(rsb_coo_idx_t nrhsi=0;nrhsi<nrhs;)
	{
		constexpr int uf = 4;
		using RT = decltype(rhs);
		using OT = decltype(out);
		const auto spmm_kf = std::make_tuple(
			&spmm_coo_partial_unrolled_by_cols<IT,NT,decltype(coop),RT,OT,uf,0>,
			&spmm_coo_partial_unrolled_by_cols<IT,NT,decltype(coop),RT,OT,uf,1>,
#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
			&spmm_coo_partial_unrolled_by_cols<IT,NT,decltype(coop),RT,OT,uf,2>,
			&spmm_coo_partial_unrolled_by_cols<IT,NT,decltype(coop),RT,OT,uf,3>,
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */
			&spmm_coo_partial_unrolled_by_cols<IT,NT,decltype(coop),RT,OT,uf,4>
		);

		switch(nrhs-nrhsi)
		{
#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
			case(1):
			{
				constexpr auto nrhsb=1;
				errval |= (std::get<nrhsb>(spmm_kf))(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
			case(2):
			{
				constexpr auto nrhsb=2;
				errval |= (std::get<nrhsb>(spmm_kf))(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
			case(3):
			{
				constexpr auto nrhsb=3;
				errval |= (std::get<nrhsb>(spmm_kf))(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */
			default:
			{
				constexpr auto nrhsb=RSBPP_PRESET_SPMM_BLOCKINGS;
				errval |= (std::get<nrhsb>(spmm_kf))(coop,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
		}
	}
	return errval;
}

template<typename NT, typename IT, typename CT = IT>
	struct rsbpp_CsrP {
		const IT nnz_;
		const IT nr_;
		const IT nc_;
		const NT* va_;
		const IT* ip_;
		const CT* ja_;
		const IT roff_, coff_;

		inline IT roff(void)const{ return roff_; }
		inline IT coff(void)const{ return coff_; }
		inline IT nr(void)const{ return nr_; }
		inline IT nc(void)const{ return nc_; }
		inline IT ip(IT at)const{ return ip_[at]; }
		inline CT ja(IT at)const{ return ja_[at]; }
		inline NT va(IT at)const{ return va_[at]; }
	};

template<typename NT, typename IT, typename CT = IT>
rsb_err_t rsbpp_csr_spmx(rsb_flags_t flags, const IT nnz, const IT nr, const IT nc, const NT* VA, const IT* IP, const CT* JA, const IT nrhs, const IT ldX, const NT* rhs, const IT ldY, NT* out, const NT* alphap, IT incx, IT incy, const rsb_trans_t transA, IT roff = 0, IT coff = 0, bool by_rows = false)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	const NT alpha = alphap ? *alphap : 1;

	if (( RSB_DO_FLAG_HAS( flags, RSB_FLAG_SYMMETRIC  ) ||
	      RSB_DO_FLAG_HAS( flags, RSB_FLAG_HERMITIAN  ) ) &&
		( nr != nc && roff == coff ) )
		return RSB_ERR_BADARGS;

	const rsbpp_CsrP<NT,IT,CT> csrp {nnz,nr,nc,VA,IP,JA,roff,coff};
	const IT bidx = std::lower_bound(IP+0,IP+nr+1,IT{  1})-IP-1;
	const IT eidx = std::lower_bound(IP+0,IP+nr+1,IT{nnz})-IP;

	if( by_rows )
	{
		const auto ldY = incy;
		const auto ldX = incx;
		constexpr int uf = 2;
		using RT = decltype(rhs);
		using OT = decltype(out);

		rsb_coo_idx_t nrhsi=0;

#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
		for(;nrhsi+100-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,100>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 100;

		for(;nrhsi+50-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,50>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 50;

		for(;nrhsi+32-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,32>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 32;

		for(;nrhsi+16-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,16>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 16;

		for(;nrhsi+8-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,8>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 8;

		for(;nrhsi+4-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,4>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 4;

		for(;nrhsi+3-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,3>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 3;

		for(;nrhsi+2-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,2>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 2;

#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */
		for(;nrhsi+1-1<nrhs;)
				errval |= spmm_csr_partial_unrolled_by_rows<IT,CT,NT,decltype(csrp),RT,OT,uf,1>
					(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx),
				nrhsi += 1;
		return errval;
	}

	for(rsb_coo_idx_t nrhsi=0;nrhsi<nrhs;)
	{
		constexpr int uf = 4;
		using RT = decltype(rhs);
		using OT = decltype(out);
		const auto spmm_kf = std::make_tuple(
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,1>,
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,1>,
#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,2>,
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,3>,
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,4>,
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,50>,
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,100>
#else /* RSBPP_PRESET_SPMM_BLOCKINGS */
			&spmm_csr_partial_unrolled_by_cols<IT,CT,NT,decltype(csrp),RT,OT,uf,1>
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */
		);

		switch(nrhs-nrhsi)
		{
#if (RSBPP_PRESET_SPMM_BLOCKINGS>1)
			case(1):
			{
				constexpr auto nrhsb=1;
				errval |= (std::get<nrhsb>(spmm_kf))(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
			case(2):
			{
				constexpr auto nrhsb=2;
				errval |= (std::get<nrhsb>(spmm_kf))(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
			case(3):
			{
				constexpr auto nrhsb=3;
				errval |= (std::get<nrhsb>(spmm_kf))(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
			case(50):
			{
				constexpr auto nrhsb=50;
				errval |= (std::get<5>(spmm_kf))(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
			case(100):
			{
				constexpr auto nrhsb=100;
				errval |= (std::get<6>(spmm_kf))(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
#endif /* RSBPP_PRESET_SPMM_BLOCKINGS */
			default:
			{
				constexpr auto nrhsb=RSBPP_PRESET_SPMM_BLOCKINGS;
				errval |= (std::get<nrhsb>(spmm_kf))(csrp,flags,transA,ldX,rhs,ldY,out,nrhsi,alpha,bidx,eidx);
				nrhsi+=nrhsb;
				break;
			}
		}
	}
	return errval;
}

inline bool is_valid_trans(rsb_trans_t transA)
{
	const std::array<rsb_trans_t,3> vt = { RSB_TRANSPOSITION_N, RSB_TRANSPOSITION_T, RSB_TRANSPOSITION_C};

	return std::find(vt.begin(),vt.end(),transA) != vt.end();
}

#endif /* RSBPP_HPP_INCLUDED */
