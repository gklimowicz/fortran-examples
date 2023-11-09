/*
Copyright (C) 2015-2022 Michele Martone

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
// Goal of this effort is to implement the RSB (Recursive Sparse Blocks) format in modern C++.
#include "config.h"
#if HAVE_FILESYSTEM
#include <filesystem>
#endif /* HAVE_FILESYSTEM */
#include <iomanip> // quoted
#include <cstdlib> // EXIT_SUCCESS
#include "rsbpp.hpp"
#include "rsbck.hpp"

#if RSBPP_WANT_ALL
template <typename IT, typename NT>
const // const here forbids assignment to return values, namely (A+A)=A
Coo<IT,NT> operator +(const Coo<IT,NT> &a, const Coo<IT,NT> &b)
{
	Coo<IT,NT> rv{a};
	rv+=b;
	return rv;
}

#if 0
//template <typename IT=int, typename NT=double>
template <typename IT, typename NT>
std::ostream & operator << (std::ostream & os, 
	//const typename Coo<int,double>::triple_t & coo)
	const typename Coo<IT,NT>::triple_t & coo)
	//std::pair < std::pair < IT, IT >, NT > &coo)
	//const std::pair < std::pair < IT, IT >, NT > & coo)
{
	os << "i:"<< coo.i() << " j:" << coo.j() << " v:" << coo.v();
	return os;
	//return coo.print(os);
}
#endif

static char g_ordering = '?';

template <typename IT, typename NT>
std::ostream & operator << (std::ostream & os, const Coo<IT,NT> & cm)
{
	return cm.print(os);
}

template <typename NT>
std::ostream & operator << (std::ostream & os, const std::vector<NT> & vec)
{
	os << "[" << vec.size() << "]: ";
     	for (const auto vel : vec)
		os << vel << " ";
	os << std::endl;
	return os;
}

class MatrixFactory final
{
	public:

	enum Type {Diagonal, Lower, Upper, Full, Symmetric};
	
	static
	enum Type charToType(const char c)
	{
		switch(toupper(c))
		{
			case('D'):
				return Diagonal;
			case('L'):
				return Lower;
			case('U'):
				return Upper;
			case('S'):
				return Symmetric;
			case('F'):
			default:
			   return Full;
		};
	}

	template <typename IT=int>
	static
	IT expect_nnz(enum Type type, IT n)
	{
		switch(type)
		{
			case(Diagonal):
				return n;
			case(Lower):
				return HALFSQUARE(n);
			case(Upper):
				return HALFSQUARE(n);
			case(Full):
				return n*n;
			default:
			   return -1;
		};
	}

	template <typename IT, typename NT>
	static
	void test_shaped(IT m, IT M, IT I, char symchar, char pc, char sfc)
	{
		for (auto n = m; n <= M ; n+=I)
		{
			enum Coo<IT,NT>::Ordering ordering = Coo<IT,NT>::OrderingDefault;
			auto A = build<IT,NT>(charToType(symchar),n);
			std::cout << n << " x " << n << " " << pc << ":" ;
			std::cout << " typeid(" << typeid(A).name() << ") :"<< std::endl;
			if(g_verbosity>3)
				std::cout << A << std::endl;
			A.selftests(expect_nnz(charToType(pc),n),sfc);
			if(g_ordering == 'R')
				ordering = Coo<IT,NT>::OrderingRSB;
			if(g_ordering == 'r')
				ordering = Coo<IT,NT>::OrderingCOR;
			A.setOrdering(ordering);
			A.sort();
			spmv_benchmark(A,100,10.0);
		}
	}

	template <typename IT>
	static
	bool test(const std::string ss = "**vqm1M100sDsLsUsFsS")
	{
		std::cout << "Tester will be driven by string \"" << ss << "\"" << std::endl;
		// TODO: indicate on which type.
#if 1
		std::istringstream is(ss);
		IT m=1,M=3,I=1;
		char pc;
		char sfc{'r'}; // FIXME

		g_verbosity = 0;

		while( (pc=is.peek()) != -1 && isprint(pc) )
		switch (pc)
		{
			case('*'):
			{
				is.ignore(1);
	       			is >> sfc;
			}
			break;
			case('q'):
				is.ignore(1);
				g_verbosity--;
			break;
			case('v'):
				is.ignore(1);
				g_verbosity++;
			break;
			case('m'):
				is.ignore(1);
			       	is >> m;
			break;
			case('M'):
				is.ignore(1);
				is >> M;
			break;
			case('I'):
				is.ignore(1);
				is >> I;
			break;
			case('C'):
				is.ignore(1);
				is >> g_cbs;
				std::cout << "Cache size option: " << g_cbs << "\n";
			break;
			case('r'):
			{
				g_nrhs_a = {};
				do {
					is.ignore(1);
					IT nrhs;
					is >> nrhs;
					if(nrhs)
					{
						std::cout << "Adding " << nrhs << " to nrhs\n";
						g_nrhs_a.push_back(nrhs);
					}
					else
					{
						std::cout << "Use empty nrhs (no spmm)\n";
						g_nrhs_a.clear();
					}
				} while( (pc=is.peek()) == ',');
			}
			break;
			case('s'):
			{
				is.ignore(1);
				char symchar;
				is >> symchar;
				switch (symchar)
				{
					case('S'):
					case('D'):
					case('L'):
					case('U'):
					case('F'):
						break;
					default:
						std::cout << symchar << " is a wrong symmetry char!\n";
						throw;
				}
				std::cout << "Will test matrix shaped " << pc <<" on sizes " << m << " to " << M << " with increase " << I << "." << std::endl;

				for(auto tc : g_types_a)
				switch(tc)
				{
					case('d'): test_shaped<IT,double>(m,M,I,symchar,pc,sfc); break;
					case('s'): test_shaped<IT,float>(m,M,I,symchar,pc,sfc); break;
					case('c'): test_shaped<IT,std::complex<float>>(m,M,I,symchar,pc,sfc); break;
					case('z'): test_shaped<IT,std::complex<double>>(m,M,I,symchar,pc,sfc); break;
					default: throw;
				}
				is.ignore(1);
			}
			break;
			case('o'):
			{
				is.ignore(1);
				is >> g_ordering;
				std::cout << "Ordering spec: " << g_ordering << "\n";
				assert(g_ordering=='R' || g_ordering=='r');
			}
			break;
			case('t'):
			{
				g_trans_a = {};
				do {
					is.ignore(1);
					std::map<char,rsb_trans_t> c2t = {
						{'N',RSB_TRANSPOSITION_N},
						{'C',RSB_TRANSPOSITION_C},
						{'T',RSB_TRANSPOSITION_T},
					};
					is >> pc;
					pc = std::toupper(pc);
					std::cout << "Adding " << pc << " to transA\n";
					g_trans_a.push_back(c2t[pc]);
				} while( (pc=is.peek()) == ',');
			}
			break;
			case('T'):
			{
				g_types_a = {};
				do {
					is.ignore(1);
					is >> pc;
					pc = std::tolower(pc);
					std::cout << "Adding " << pc << " to g_types_a\n";
					g_types_a.push_back(pc);
				} while( (pc=is.peek()) == ',');
			}
			break;
			case('0'): case('1'): case('2'): case('3'): case('4'): 
			case('5'): case('6'): case('7'): case('8'): case('9'):
				is >> M, m=M;
			break;
			default:
				std::cout << "Ignoring specifier '" << pc << "'" << std::endl;
				is.ignore(1);
		}
#else
		for (auto n = 1; n < 3 ; ++n) 
		{
			auto A = build<IT,NT>(Diagonal,n);
			std::cout << n << " x " << n << " Diagonal:" << std::endl << A << std::endl;
			A.selftests(n);
		}

		for (auto n = 1; n <= 4 ; ++n) 
		{
			auto A = build<IT,NT>(Lower,n);
			std::cout << n << " x " << n << " Lower:" << std::endl << A << std::endl;
			A.selftests(HALFSQUARE(n));
		}

		for (auto n = 1; n <= 4 ; ++n) 
		{
			auto A = build<IT,NT>(Upper,n);
			std::cout << n << " x " << n << " Upper:" << std::endl << A << std::endl;
			A.selftests(HALFSQUARE(n));
		}

		for (auto n = 1; n <= 4 ; ++n) 
		{
			auto A = build<IT,NT>(Full,n);
			std::cout << n << " x " << n << " Full:" << std::endl << A << std::endl;
			A.selftests(n*n);
		}
#endif
		return true;
	}

	template <typename IT, typename NT>
	static
	Coo<IT,NT>
	build_from_coo(typename Coo<IT,NT>::trivec_t && coo, rsb_flags_t sym_flags, enum Coo<IT,NT>::Ordering ordering)
	{
		auto A = Coo<IT,NT>{};
		A.init_coo(std::forward<typename Coo<IT,NT>::trivec_t>(coo), sym_flags, ordering);
		A.sort();
		return A;
	}

	template <typename IT, typename NT>
	[[deprecated("only transitory")]]
	static
	Coo<IT,NT>
	build(Type type, const IT n, bool also_check=true)
	{
		const rsb_flags_t flags{type == Symmetric ? RSB_FLAG_SYMMETRIC : RSB_FLAG_NOFLAGS};
		auto A = Coo<IT,NT>{flags};

		switch (type)
		{
			case(Diagonal):
			for (auto i = 0; i< n; ++i)
				A.push_back(i,i,(VALFROMIJ<NT>(i,i)));
			break;
			case(Full):
			for (auto i = 0; i< n; ++i)
				for (auto j = 0; j< n; ++j)
					A.push_back(i,j,(VALFROMIJ<NT>(i,j)));
			break;
			case(Symmetric):
			case(Lower):
			for (auto i = 0; i< n; ++i)
				for (auto j = 0; j<=i; ++j)
					A.push_back(i,j,(VALFROMIJ<NT>(i,j)));
			break;
			case(Upper):
			for (auto i = 0; i< n; ++i)
				for (auto j = i; j< n; ++j)
					A.push_back(i,j,(VALFROMIJ<NT>(i,j)));
			break;
			default:
		       	throw;
		}

		A.sort();
		// TODO: check here, per-type

		if(also_check==false)
			return A;

		// continue with checks.
		// FIXME: finish this.
		const NT alpha{3},beta{1},xV{1},yV{1};
		rsb_trans_t transA {RSB_TRANSPOSITION_N};// FIXME
		std::valarray<NT> x(A.nc()), y(A.nr()), yR(A.nr());
		x=xV;
		y=yR=yV;
      		//for (auto p = x.begin(); p < x.end() ; ++p ) *p = 2.0;
      		//for (auto p = y.begin(); p < y.end() ; ++p ) *p = 1.0;
		//NT ysum0=std::accumulate(y.begin(),y.end(),NT());
		A.spmv(transA,x,y,alpha,beta);
		
		switch (type)
		{
			case(Diagonal):
			for (auto i = 0; i< n; ++i)
			{
				NT rR = (alpha*VALFROMIJ<NT>(i,i)*x[i]) + yR[i]*beta;
		       		if( rR != y[i] )
					std::cout << RSB_ERRM_FIFL << " @ " << i << " : " << rR << " != " << y[i] << std::endl,
		       			throw;
			}
			break;
			case(Full):
			for (auto i = 0; i< n; ++i)
			{
				NT rR = yR[i]*beta;
				for (auto j = 0; j< n; ++j)
					rR += (alpha*VALFROMIJ<NT>(i,j)*x[j]);
		       		if( rR != y[i] )
					std::cout << RSB_ERRM_FIFL << " @ " << i << " : " << rR << " != " << y[i] << std::endl,
		       			throw;
			}
			break;
			case(Symmetric):
			for (auto i = 0; i< n; ++i)
			{
				NT rR = yR[i]*beta;
				for (auto j = 0; j<=i; ++j)
					rR += (alpha*VALFROMIJ<NT>(i,j)*x[j]) + tzero<NT>*yR[i]*beta;
				for (auto j = i+1; j<n; ++j)
					rR += (alpha*VALFROMIJ<NT>(j,i)*x[i]) + tzero<NT>*yR[i]*beta;
		       		if( rR != y[i] )
					std::cout << RSB_ERRM_FIFL << " @ " << i << " : " << rR << " != " << y[i] << std::endl,
		       			throw;
			}
			break;
			case(Lower):
			for (auto i = 0; i< n; ++i)
			{
				NT rR = yR[i]*beta;
				for (auto j = 0; j<=i; ++j)
					rR += (alpha*VALFROMIJ<NT>(i,j)*x[j]);
		       		if( rR != y[i] )
					std::cout << RSB_ERRM_FIFL << " @ " << i << " : " << rR << " != " << y[i] << std::endl,
		       			throw;
			}
			break;
			case(Upper):
			for (auto i = 0; i< n; ++i)
			{
				NT rR = yR[i]*beta;
				for (auto j = i; j< n; ++j)
					rR += (alpha*VALFROMIJ<NT>(i,j)*x[j]);
		       		if( rR != y[i] )
					std::cout << RSB_ERRM_FIFL << " @ " << i << " : " << rR << " != " << y[i] << std::endl,
		       			throw;
			}
			break;
			default:
		       	throw;
		}

		return A;
	} /* build */
}; /* MatrixFactory */

template <typename IT, typename TT=mytime_t>
void test(/*const int nit, */std::string ts="")
{
	if(ts!="")
		MatrixFactory::test<IT>(ts);
	else
		MatrixFactory::test<IT>();
}

template <typename IT, typename NT, typename TT=mytime_t>
void spmv_benchmark(const Coo<IT,NT> & A, const int maxit = 100, TT maxt=10.0/*, const std::vector<rsb_trans_t> & trans_a_ = {}*/)
{
	const bool sym {RSB_DO_FLAG_HAS(A.rsbflags(),RSB_FLAG_SYMMETRIC)?true:false};
	const IT cxf { std::is_scalar<NT>::value ? 1 : 3 };
	std::chrono::time_point<std::chrono::system_clock> start, end;
	const IT mps=1; /* max printable size */
	std::chrono::duration<TT> elapsed_seconds;
	int nit = 0;
	TT dt=TT();
	TT st=TT();
	TT zt=TT();
	//NT zero{};
	const NT zero{};
	const NT alpha{2};
	const NT beta{1}; // Note: checksum only effective as long as beta==1
	const NT one{1};
	//const std::string irs = "";
	const auto irs = "   "; // intra-record separator (not quite the same as endl)
	//const auto irs = "\n"; // intra-record separator (not quite the same as endl)

	assert(maxit>0);
	assert(maxt>0);

	std::cout << "Benchmarking " << A.nc() << "x" << A.nr() << ":" << A.nnz() << " with numerical type " << typeid(NT).name();
	if ( sym )
		 std::cout << " symmetric ";
	std::cout << " and ordering " << (char)A.getOrdering();
	std::cout << " and " << A.blocks() << " blocks.";
	std::cout << std::endl;

#if 0
	// TODO: these modify A: move them elsewhere.
	start = std::chrono::system_clock::now();
	A.zort();
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	zt = elapsed_seconds.count();
	if (A.nr() < mps)
		std::cout << "after zort:" << std::endl,
		std::cout << A;
	A.sort_coc();
	assert( A.sort_coc_check() );
	A.sort_cor();
	if (A.nr() < mps)
		std::cout << "after sort_cor:" << std::endl,
		std::cout << A;
	assert( A.sort_cor_check() );
#endif

	{
		auto B{A};
		start = std::chrono::system_clock::now();
		B.sort();
		end = std::chrono::system_clock::now();
		elapsed_seconds = end-start;
		st = elapsed_seconds.count();
	}

	for( const rsb_trans_t transA : g_trans_a )
{
	std::cout << std::endl;
	std::cout << " use transA = " << (char)(transA) << std::endl;
#if 1
	std::vector<NT> x(A.get_ldB(transA)), y(A.get_ldC(transA));
	if (A.nc() < mps)
		std::cout << "x:" << x;
	if (A.nr() < mps)
		std::cout << "y:" << y;
	//std::valarray<NT> x(A.nc()), y(A.nr());
	std::uninitialized_fill(x.begin(),x.end(),2.0);
	std::uninitialized_fill(y.begin(),y.end(),1.0);

	const NT ysum0=std::accumulate(y.begin(),y.end(),NT());
	const auto mocc = (sizeof(NT)+2*sizeof(IT))*A.nnz(); // matrix occupation
	const auto lws = sizeof(NT)*(sym?2:1)*1*(A.nr()+A.nc()) + mocc; // assume all data read once
	const auto hws = sizeof(NT)*(sym?2:1) * A.nnz() * 2 * 1 + mocc; // matrix once, vectors elements once per nnz

	A.spmv(transA,x,y,+alpha,beta);
	NT ysum=std::accumulate(y.begin(),y.end(),NT());
	if(alpha==zero)
		std::cout << "WARNING: Using null alpha!" << std::endl;

	if(ysum == ysum0 && alpha!=zero)
	{
		std::cout << "WARNING: check sum is same as " << ysum0 << " : " << ysum << "!!" << std::endl,
		std::cout << "Maybe this SPMV case is UNFINISHED!! " << std::endl;
		throw;
	}
	A.spmv(transA,x,y,-alpha,beta);
	ysum=std::accumulate(y.begin(),y.end(),NT());
	if(beta==one)
	if(ysum != ysum0)
	{
		std::cout << "WARNING: first check sum shall be " << ysum0 << " but is  " << ysum << " (differs by " << (ysum-ysum0) << " )" << std::endl,
		std::cout << "Maybe this SPMV case is UNFINISHED!! " << std::endl;
		throw;
	}
	start = std::chrono::system_clock::now();
	for(nit = 0, dt = 0; nit <maxit && dt < maxt; nit+=2)
	{
		A.spmv(transA,x,y,+alpha,beta);
		A.spmv(transA,x,y,-alpha,beta);
		end = std::chrono::system_clock::now();
		elapsed_seconds = end-start;
		dt = elapsed_seconds.count();
		// TODO: cache flushing is missing; via lambda ?
	}

	ysum=std::accumulate(y.begin(),y.end(),NT());
	if(beta==one)
	if(ysum != ysum0)
	{
		std::cout << "WARNING: last check sum shall be " << ysum0 << " but is  " << ysum << " (differs by " << (ysum-ysum0) << " )" << std::endl,
		std::cout << "SPMV is UNFINISHED!! " << std::endl;
		throw;
	}
	//std::cout << "check sum is : " << ysum << std::endl;
/*
	if (A.nr() < mps)
		std::cout << "y:" << y;
		*/

	std::cout << " mflops: " << cxf*(sym?2:1)*((2.0 * A.nnz()) / (dt/nit))/1.0e6 << irs;
	std::cout << " elapsed   time: " << dt << "s" << irs;
	std::cout << " iterations: " << nit << irs;
	if(sym)
	       	std::cout << " symmetric";
	std::cout << " spmv time: " << dt/nit << "s" << irs;
	std::cout << " min memory bw: " << (lws/1e9)/(dt/nit) << "GB/s" << irs;
	std::cout << " max memory bw: " << (hws/1e9)/(dt/nit) << "GB/s" << irs;
	std::cout << " Z-sort to spmv time: " << zt/(dt/nit) << "x" << irs;
	std::cout << " R-sort to spmv time: " << st/(dt/nit) << "x" << std::endl;

#endif

	for( const IT nrhs : g_nrhs_a )
	{
		const IT ldB = A.get_ldB(transA);
		const IT ldC = A.get_ldC(transA);
		const std::vector<NT> x(ldB*nrhs,2.0);

		std::cout << std::endl;

		std::vector<NT> y(ldC*nrhs);

		std::uninitialized_fill(y.begin(),y.end(),1.0);
		A.spmm(transA,+alpha,nrhs,x.data(),beta,y.data());
		const NT ysum0=std::accumulate(y.begin(),y.end(),NT());
		const std::vector<NT> y0(y);

		start = std::chrono::system_clock::now();
		for(nit = 0, dt = 0; nit <maxit && dt < maxt; nit+=2)
		{
			A.spmm(transA,+alpha,nrhs,x.data(),beta,y.data());
			A.spmm(transA,-alpha,nrhs,x.data(),beta,y.data());
			end = std::chrono::system_clock::now();
			elapsed_seconds = end-start;
			dt = elapsed_seconds.count();
			// TODO: cache flushing is missing; via lambda ?
		}
		NT ysum=std::accumulate(y.begin(),y.end(),NT());
		if(beta==one)
		if(ysum != ysum0)
		{
			std::cout << "WARNING: check sum is not " << ysum0 << " but  " << ysum << " (differs by " << (ysum0-ysum) << " )" << std::endl,
			std::cout << "SPMM is UNFINISHED!! " << std::endl;

			for (decltype(y.size()) i=0;i<y.size();++i)
				if(y[i]!=y0[i])
					std::cout << " differs at " << i << "; has " << y[i] << " vs " << y0[i] << std::endl;

			throw;
		}

		const auto lws = sizeof(NT)*(sym?2:1)*nrhs*(A.nr()+A.nc()) + mocc; // assume all data read once
		const auto hws = sizeof(NT)*(sym?2:1) * A.nnz() * 2 * nrhs + mocc; // matrix once, vectors elements once per nnz
		std::cout << " mflops: " << cxf*(sym?2:1)*nrhs*((2.0 * A.nnz()) / (dt/nit))/1.0e6 << irs;
		std::cout << " elapsed   time: " << dt << "s" << irs;
		std::cout << " iterations: " << nit << irs;
		if(sym)
		       	std::cout << " symmetric";
		std::cout << " spmm-" << nrhs << " time: " << dt/nit << "s" << irs;
		std::cout << " min memory bw: " << (lws/1e9)/(dt/nit) << "GB/s" << irs;
		std::cout << " max memory bw: " << (hws/1e9)/(dt/nit) << "GB/s" << irs;
		std::cout << " Z-sort to spmm-" << nrhs << " time: " << st/(dt/nit) << "x" << std::endl;

		ysum=std::accumulate(y.begin(),y.end(),NT());
	}
}
	std::cout << "." << std::endl;

	//(A+A)=A; // illegal :-)
	auto nnzA=A.nnz();
	//std::cout << "   A .nnz(): " << A.nnz() << std::endl;
	//std::cout << "(A+A).nnz(): " << (A+A).nnz() << std::endl;
	{
		auto B{A};
		B=(B+B); // TODO: waiting for an implementation
		//std::cout << "(A+A).nnz(): " << A.nnz() << std::endl;
		assert( B.nnz() == 2*nnzA );// TODO: valid e.g. as long as no coo sorting/compression exists..
	}
} /* benchmark */

template <typename IT, typename NT, typename TT=mytime_t>
void benchmark(const int nit, std::string filename="")
{
	Coo<IT,NT> A;

	if(filename == "")
	{
		const IT n=20500;
		for (auto i=0;i<n;++i)
		       	A.push_back(i,i,i*1.0);
	}
	else
		std::cout << "benchmarking " << filename << std::endl,
		A=Coo<IT,NT>(filename);

	A.sort();
	spmv_benchmark(A, nit);
} /* benchmark */

// using myint = int64_t;
using myint = int;
// using myint = size_t ;
// using myint = long long int ;
// using myint = short unsigned int ;
// using myint = short int ;
// using myint = uint16_t ;
// using myint = uint8_t ;

// using myfloat = float ;
using myfloat = double ;
//  using myfloat = long double ;
//  using myfloat = uint8_t ;
//  using myfloat = int ;
//  using myfloat = uint64_t ;
//  using myfloat = std::complex<float> ;
//  using myfloat = std::complex<double> ;
//  using myfloat = std::complex<long double> ; // FIXME

template <typename IT, typename TT=mytime_t>
int main_all(int argc, char **argv)
{
	//std::cout << "TEST BEGIN: " << STRINGIFY(IT) << std::endl; 

	if(argc > 1)
		test<IT>(argv[1]);//,
	else
		test<IT>();//,
	//std::cout << "TEST DONE: " << STRINGIFY(IT) << std::endl; 
	return 0;
}

template <typename TT=mytime_t>
TT rsbpp_time(void)
{
	static std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
	std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
	std::chrono::duration<TT> elapsed_seconds = end-start;
	return elapsed_seconds.count();
}

#define WANT_MICROBENCHMARK 0
#if WANT_MICROBENCHMARK
template <int NS, typename AT >
auto arithmetic_loop(const AT & a)
{
	const int nb = NS;
	const auto b = a.cbegin(), e = a.cend();

	assert( a.size() >= nb );
	assert( a.size() >  0 );
	assert( nb >= 1 );
	assert( nb <= 64 );
	assert( a.size() % 4 == 0 );

	auto d = a[0]-a[0];

#if 0
      	for (auto ci = b+nb ; ci < e; ++ci)
	{
      		for (auto bi = 0 ; bi < nb; ++bi)
			d += ci[-bi];
	}
#endif
#if 0
      	for (auto ci = b+nb ; ci < e; ++ci)
	{
		auto e = ci[0];
      		for (auto bi = 0 ; bi < nb; ++bi)
			e*=2,
			d += e;
	}
#endif
#if 0
      	for (auto ci = b+nb ; ci < e; ++ci)
	{
		auto e = ci[0];
      		for (auto bi = 1 ; bi < nb; ++bi)
			e* = ci[-bi];
		d += e;
	}
#endif
#if 0
	auto f = b[0];
	assert(f*f==f);
      	for (auto ci = b+nb ; ci < e; ++ci)
	{
		f = ci[0];
      		for (auto bi = 1 ; bi < nb; ++bi)
			f = f*f;
		d += f;
	}
#endif
#if 0
	//AT::value_type f[nb];
	//decltype(d) f[nb];
	std::array<decltype(d),nb>  f;

	for(auto i=0;i<nb;++i)
		f[i]=b[i];

      	for (auto ci = b ; ci+nb-1 < e; ++ci)
	{
      		for (auto bi = 0 ; bi < nb; ++bi)
			d += f[bi]*ci[bi];
	}
#endif
#if 1
	std::array<decltype(d),nb>  f;
	const auto maxi = (e-b)-nb+1;

	for(auto i=0;i<nb;++i)
		f[i]=b[i];
#pragma omp parallel for schedule(static,1) default(shared)
      	for (int i = 0 ; i < maxi; ++i)
      	//for (auto ci = b ; ci < (e-nb+1); ++ci)
	{
      		for (auto bi = 0 ; bi < nb; ++bi)
			d += f[bi]*b[i+bi];
		//for (auto bi = 0 ; bi < nb; ++bi)
		//	d += f[bi]*ci[bi];
	}
#endif
#if 0
      	for (auto ci = b+nb ; ci+4-1 < e; ci+=4)
	{
      		for (auto bi = 0 ; bi < nb; ++bi)
			d += ci[0-bi];
      		for (auto bi = 0 ; bi < nb; ++bi)
			d += ci[1-bi];
      		for (auto bi = 0 ; bi < nb; ++bi)
			d += ci[2-bi];
      		for (auto bi = 0 ; bi < nb; ++bi)
			d += ci[3-bi];
	}
#endif
	return d;
}

template <int NS, typename AT >
auto timed_arithmetic_loop(const AT & a, const mytime_t mt=1.0)
{
	mytime_t t0=rsbpp_time(), t1=t0;
	int nit;
	auto r = a[0]-a[0];

	for(nit=0,t0=rsbpp_time();(t1=rsbpp_time())-t0<mt;++nit)
       		r += arithmetic_loop<NS,decltype(a)>(a);

	std::cout << "in " << (t1-t0) << "s, " 
		<< nit << " iterations, "
	       	<< NS << " ops per iteration, " 
	       	<< omp_get_max_threads() << " threads, " 
		<< (((1.e-6*a.size())*nit)*NS)/(t1-t0)  << " MOPS"
	       	<< std::endl;

	return r;
}

template <typename AT = int>
int roofline_benchmark(void)
{
	using rt_t = AT;
	//using rt_t = int;
	rt_t iv=1;
	std::vector<rt_t> a(16*1024*1024,iv);
	rt_t r = a[0]-a[0];
	const mytime_t mt{.5};

	r+=timed_arithmetic_loop<1,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<2,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<4,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<8,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<10,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<12,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<13,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<14,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<15,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<16,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<17,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<18,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<19,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<20,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<24,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<28,decltype(a)>(a,mt);
	r+=timed_arithmetic_loop<32,decltype(a)>(a,mt);
	//r+=timed_arithmetic_loop<64,decltype(a)>(a,mt);

	//std::cout << ((r%99991==0)?" ":"  ") << std::endl;
	//std::cout << ((r%3==0)?" ":"  ") << std::endl;
	std::cout << ((r>0)?" ":"  ") << std::endl;

	return r;
}
#endif /* WANT_MICROBENCHMARK */

#if 0
template<class IT, class NT>
typename Coo<IT,NT>::trivec_t read_mm_triplets_c(std::istream & is, IT nnz)
{
	// catch-all version
	typename Coo<IT,NT>::trivec_t coo;
		IT i,j;
		NT rp,ip;
		NT v;
		NT ci=0;
		if(!std::is_scalar<NT>())
			ci=std::sqrt(-1);
		while(nnz-- && is >> i && is >> j && is >> rp && is >> ip )
		{
			v=rp+ci*ip;
			coo.push_back({{i-1,j-1},{v}});
		}
	return coo; 
}
#endif

template<class IT, class NT>
typename Coo<IT,
	typename std::enable_if<! std::is_scalar<NT>::value,NT>::type
	 >::trivec_t read_mm_triplets_c(std::istream & is, IT nnz)
{
	typename Coo<IT,NT>::trivec_t coo;
	typename NT::value_type rp,ip;
	IT i,j;
	while(nnz-- && is >> i && is >> j && is >> rp && is >> ip )
		coo.push_back({{i-1,j-1},{rp,ip}});
	return coo; 
}

template<class IT, class NT>
typename Coo<IT,
	typename std::enable_if<std::is_scalar<NT>::value,NT>::type
	 >::trivec_t read_mm_triplets_c(std::istream & is, IT nnz)
{
	typename Coo<IT,NT>::trivec_t coo;
	IT i,j;
	NT rp,ip;
	while(nnz-- && is >> i && is >> j && is >> rp && is >> ip ) // ignore ip
		coo.push_back({{i-1,j-1},{rp}});
	return coo; 
}

template<class IT, class NT>
typename Coo<IT,NT>::trivec_t read_mm_triplets_r(std::istream & is, IT nnz)
{
	typename Coo<IT,NT>::trivec_t coo;

	IT i,j;
	NT v;
	while(nnz-- && is >> i && is >> j && is >> v )
		coo.push_back({{i-1,j-1},{v}});
	return coo; 
}

template<class IT, class NT, bool cmplxMtx>
typename Coo<IT,NT>::trivec_t read_mm_triplets(std::istream & is, IT nnz)
{
	if(cmplxMtx)
		return read_mm_triplets_c<IT,NT>(is,nnz);
	else
		return read_mm_triplets_r<IT,NT>(is,nnz);
}

template<class IT, class NT>
typename Coo<IT,NT>::trivec_t read_mm_as_coo(std::istream & is, rsb_flags_t & sym_flags)
{
	using trivec_t = typename Coo<IT,NT>::trivec_t;
	trivec_t coo;
	using ct_t = typename trivec_t::size_type;
	std::string tok;
	bool cmplx{false};
	IT nr{0},nc{0},nnz{0};
	
	if( (!(is >> tok)) || tok!="%%MatrixMarket" )
		throw;
	if( (!(is >> tok)) || tok!="matrix" )
		throw;
	if( (!(is >> tok)) || tok!="coordinate" )
		throw;
	while(is.peek()!='\n' && is >> tok && !tok.empty() && std::isalpha(tok[0]) )
	{
		if( tok=="complex" )
			cmplx=true;
		if( tok=="symmetric" )
			sym_flags = RSB_FLAG_SYMMETRIC;
		if( tok=="hermitian" )
			sym_flags = RSB_FLAG_HERMITIAN;
		if( tok=="general" )
			sym_flags = RSB_FLAG_NOFLAGS;
	}
	if(is.peek()!='\n')
	       throw;
	is.get();
	while(is.peek()=='%')
		std::getline(is,tok);
	is >> nr >> nc >> nnz;
	std::vector<IT> IA(nnz,0),JA(nnz,0);
	coo.reserve(nnz);
	if(cmplx)
		coo = read_mm_triplets<IT,NT,true>(is,nnz);
	else
		coo = read_mm_triplets<IT,NT,false>(is,nnz);

	if(coo.size() != static_cast<ct_t>(nnz))
		throw;

	return coo;
}


template<class IT, class NT>
void bench_mtx_file(const std::string & fn)
{
	rsb_flags_t sym_flags {}; // FIXME: still unused
	auto is {std::ifstream(fn,std::ios::in)};

	std::cout << "Reading matrix " << std::quoted(fn) << " ..." << std::endl;

	auto coo = read_mm_as_coo<IT,NT>(is,sym_flags);
	mytime_t dt = rsbpp_time();
	dt = rsbpp_time() - dt;
	const bool to_ones = true;
#if HAVE_FILESYSTEM && __cplusplus >= 201703L
	const auto fs = std::filesystem::file_size(fn);
	std::cout << "Matrix I/O took " << dt << "s at " << fs/(1024*1024*dt) << " MB/s" << std::endl;
#endif /* HAVE_FILESYSTEM */
	if(to_ones)
		for(decltype(coo.size()) i=0;i<coo.size();++i)
			coo[i].v()={1};
	enum Coo<IT,NT>::Ordering ordering = Coo<IT,NT>::OrderingDefault;
	if(g_ordering == 'R')
		ordering = Coo<IT,NT>::OrderingRSB;
	if(g_ordering == 'r')
		ordering = Coo<IT,NT>::OrderingCOR;
	auto A = MatrixFactory::build_from_coo<IT,NT>(std::move(coo), sym_flags, ordering);
	A.selftests();
	A.sort();
	spmv_benchmark(A,100,10.0);
}

int main(int argc, char **argv)
{
	int res {EXIT_SUCCESS};

	for(int argi = 1; argi<argc; ++argi)
	{
		std::string fn = argv[argi];
		const bool mmf = ( fn.size()>4 && fn.rfind(".mtx")==fn.size()-4 );
		const bool iss = ( fn == "-" );

		if(iss)
			fn="/dev/stdin";

		if( mmf || iss )
		{
			if(std::find(g_types_a.begin(),g_types_a.end(),'d')!=g_types_a.end())
				bench_mtx_file<int,double>(fn);
			if(std::find(g_types_a.begin(),g_types_a.end(),'z')!=g_types_a.end())
				bench_mtx_file<int,std::complex<double>>(fn);
			if(std::find(g_types_a.begin(),g_types_a.end(),'s')!=g_types_a.end())
				bench_mtx_file<int,float>(fn);
			if(std::find(g_types_a.begin(),g_types_a.end(),'c')!=g_types_a.end())
				bench_mtx_file<int,std::complex<float>>(fn);
		}
		else
		{
			std::cout << "TEST BEGIN ALL" << std::endl;
#if WANT_MICROBENCHMARK
			roofline_benchmark<double>();
			roofline_benchmark<float>();
			roofline_benchmark<int>();
			roofline_benchmark<char>();
			return res;
#endif /* WANT_MICROBENCHMARK */

#if 0
			// FIXME: not all yet ready !
			res += main_all<int,float>(argc,argv);
			res += main_all<int,double>(argc,argv);
			res += main_all<int,long double>(argc,argv);
			res += main_all<int,std::complex<double> >(argc,argv);
			res += main_all<int,std::complex<float> >(argc,argv);
			res += main_all<int,_Complex>(argc,argv);	// aka complex, in <ccomplex>

			res += main_all<int16_t,float>(argc,argv);
			res += main_all<int16_t,double>(argc,argv);
			res += main_all<int16_t,long double>(argc,argv);
			res += main_all<int16_t,std::complex<double> >(argc,argv);
			res += main_all<int16_t,std::complex<float> >(argc,argv);

			res += main_all<myint,myfloat>(argc,argv);
#else
			res += main_all<myint>(argc,argv);
			//res += main_all<myint,myfloat>(argc,argv);
			//res += main_all<int16_t,double>(argc,argv);
#endif
			std::cout << "TEST DONE ALL" << std::endl;
		}
	}
	return res;
}
#else /* RSBPP_WANT_ALL */
int main(void) {}
#endif /* RSBPP_WANT_ALL */
/*
 vim:number:
*/
