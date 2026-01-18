/*

Copyright (C) 2021-2022 Michele Martone

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
 \ingroup rsb_doc_examples
 @file
 @author Michele Martone

 @brief C++ example based on <rsb.hpp> converting a Matrix Market file into a custom format using RsbMatrix.get_coo().

 \include mtx2bin.cpp
 */
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
#include <fstream>
#include <cstdlib> // EXIT_FAILURE
#include <complex>
#include <array>
#include <vector>
#include <iomanip>
#include <iostream>
#if HAVE_FILESYSTEM
#include <filesystem>
#endif
#include <string>
#include <rsb.hpp>

using namespace ::rsb;

template <typename T>
size_t vecdump(std::ofstream & ofs, const std::vector<T> & v) {
	const size_t wb {v.size()*sizeof(T)};
	ofs.write( reinterpret_cast<const char*>(v.data()), wb);
	return wb;
}

template <typename IT, typename NT>
size_t coodump(std::ofstream & ofs, std::vector<IT> & ia, const std::vector<IT> & ja, const std::vector<NT> & va) {
	size_t wb { 0 };
	wb += vecdump(ofs, ia);
	wb += vecdump(ofs, ja);
	wb += vecdump(ofs, va);
	return wb;
}

template <typename IT, typename NT>
size_t mtxdump(std::ofstream & ofs, const std::string & hdr, std::vector<IT> & ia, const std::vector<IT> & ja, const std::vector<NT> & va) {
	size_t wb { 0 };
	ofs << hdr;
	wb = hdr.size() + coodump(ofs, ia, ja, va);
	return wb;
}

template <typename IT, typename NT>
size_t mtxdump(std::ofstream & ofs, IT m, IT k, IT nnz, rsb_flags_t flagsA, std::vector<IT> & ia, const std::vector<IT> & ja, const std::vector<NT> & va) {
	const std::string rsb_bin_tag { "BINARY LIBRSB MATRIX MARKET EXTENSION COO FORMAT" };
	const char rsb_bin_cookie { '\0' };
	const std::string typstr{ std::is_scalar<NT>() ? "real" : "complex" };
	const std::string symstr{ ( flagsA & RSB_FLAG_SYMMETRIC ) ? "symmetric" : ( (flagsA & RSB_FLAG_HERMITIAN) ? "hermitian" : "general" ) };
	std::ostringstream oss;
	oss << "%%MatrixMarket matrix coordinate " << typstr << " " << symstr << "\n";
	oss << m << " " << k << " " << nnz << "\n";
	oss << rsb_bin_tag << "\n";
	oss << rsb_internal::rsb_type_t_for<NT>() << rsb_bin_cookie;
	auto hdr { oss.str() };
	return mtxdump(ofs, hdr, ia, ja, va);
}

template <typename nt_t>
void mtx2bin(const std::string ifilename, const std::string ofilename) {
  	RsbLib rsblib;
  	RsbMatrix<nt_t> mtx(ifilename.c_str());
	const rsb_trans_t transA = RSB_TRANSPOSITION_N;

	std::vector<nt_t> VA(mtx.nnz());
	std::vector<rsb_coo_idx_t> IA(mtx.nnz());
	std::vector<rsb_coo_idx_t> JA(mtx.nnz());

	mtx.get_coo(transA,VA.data(),IA.data(),JA.data(),RSB_FLAG_FORTRAN_INDICES_INTERFACE);

	auto ofs = std::ofstream(ofilename,std::ios::binary);
	mtxdump(ofs, mtx.rows(), mtx.cols(), mtx.nnz(), mtx.rsbflags(), IA, JA, VA);
}

#ifdef RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS
const char dtypechar {RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS[0]};
#else
const char dtypechar {'S'};
#endif

int err_msg(const std::string argv0) {
		std::cout << "usage: " << argv0 << " matrix-input-file [matrix-output-file [type]]" << std::endl;
#if defined(RSB_BLAS_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS)
		std::cout << "with [type] among " << RSB_BLAS_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS << " ; default " << dtypechar << std::endl;
#endif
	return EXIT_FAILURE;
}

auto main(const int argc, char * argv[]) -> int
{
	const std::string ifilename{ argc > 1 ? argv[1] : "../A.mtx"};
	const std::string ofilename{ argc > 2 ? argv[2] : ifilename + ".bin"};
	const char typechar = toupper( argc > 3 ? *argv[3] : dtypechar );

#if defined(HAVE_FILESYSTEM) && (__cplusplus >= 201703)
        if ( std::filesystem::exists(ifilename) )
#else
        if( std::ifstream(ifilename) )
#endif
	switch (typechar)
	{
		case 'D': mtx2bin<double>(ifilename,ofilename); break;
		case 'S': mtx2bin<float>(ifilename,ofilename); break;
		case 'C': mtx2bin<std::complex<float>>(ifilename,ofilename); break;
		case 'Z': mtx2bin<std::complex<double>>(ifilename,ofilename); break;
		default:
			return err_msg(argv[0]);
	}
	else
	{
		return err_msg(argv[0]);
	}
}
