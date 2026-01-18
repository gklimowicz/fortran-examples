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
/* rsb.cpp */

#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <rsb.hpp>
#include <numeric>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_GETOPT_H
#include <getopt.h> /* we want GNU's getopt_long */
#endif /* HAVE_GETOPT_H */
#include <stdlib.h> /* atoi */

#ifndef RSBEP_NO_STUB
#define RSBEP_NO_STUB 0
#endif /* RSBEP_NO_STUB */

using namespace ::rsb;

#define RSBEP_NUMT double
typedef RsbMatrix<RSBEP_NUMT> Rsb_Matrix;

#if RSBEP_NO_STUB
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else /* HAVE_MPI */
#include "Epetra_SerialComm.h"
#endif /* HAVE_MPI */
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
//#include "MyEpetra_CrsMatrix.h"
#endif /* RSBEP_NO_STUB */

#define RSBEP_ABS(X) (((X) < (0)) ? -(X) : (X))
#define RSBEP_MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define RSBEP_MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define RSBEP_FREE(P) free(P)
#define RSBEP_CALLOC(N,S) calloc((N),(S))
#define RSBEP_ERR_OK RSB_ERR_NO_ERROR /* 0 */
#define RSBEP_ERR_GENERIC_ERROR RSB_ERR_GENERIC_ERROR /* -1 */
#define RSBPE_RTE(TEC) return (TEC);
#define RSBEP_ERRREACT(ERRVAL) { if (ERRVAL != RSB_ERR_NO_ERROR) return -1;/*goto ERRLABEL;*/ };
#define RSBEP_EERRREACT(ERRVAL) { if (ERRVAL != RSBEP_ERR_OK ) { RSBEP_TEST_MSG("!");return -1;}/*goto ERRLABEL;*/ };
#define RSBEP_ERRGOTO(ERRVAL,ERRLABEL) { if ( (ERRVAL) != RSB_ERR_NO_ERROR ) { RSBEP_ERROR(); goto ERRLABEL; } };
#define RSBEP_ERROR() RSBP_ERRORM("error") 
//#define RSBEP_ERROR() RSBP_ERRORM("unspecified error") 
#define RSBEP_PERR_GOTO(LABEL,MSG) {RSBEP_ERROR();goto LABEL;} /* FIXME */
#define RSBEP_IDXT rsb_coo_idx_t
#define RSBEP_ERRT int
#define RSBEP_BOOL_MAYBE true
#define RSBEP_OS std::cout
#define RSBEP_PRERR_GOTO(LABEL,MSG) {RSBEP_ERROR(); RSBEP_OS << __func__ << " : ERROR: \"" << MSG << "\"\n"; goto LABEL;} /* FIXME */
#define RSBEP_INCOMPLETE_BODY { RSBEP_VI() return RSBEP_ERR_GENERIC_ERROR; }
#define RSBEP_TEST_MSG(OKNOK) RSBEP_OS << OKNOK << "\n"
#define RSBP_PRINTF printf

/* debug macros */
#define RSBEP_GETENVI(SYM, DEF) ( getenv(SYM) ? atoi(getenv(SYM)) : (DEF) )
#define RSBEP_GETENVB(SYM, DEF) ( getenv(SYM) ? (atoi(getenv(SYM))?true:false) : (DEF) )
#define RSBEP_GETENVS(SYM, DEF) ( getenv(SYM) ?     (getenv(SYM)) : (DEF) )
//#define RSBEP_MAX_PRINT 10
//#define RSBEP_MAX_PRINTV 10
#define RSBEP_MAX_PRINT  RSBEP_GETENVI("RSBEP_MAX_PRINT" ,10)
#define RSBEP_MAX_PRINTV RSBEP_GETENVI("RSBEP_MAX_PRINTV",10)

#ifndef RSBP_WANT_CXX11
#define RSBEP_IFCXX11_BRACKETS(X) {X}
#else /* RSBP_WANT_CXX11 */
#define RSBEP_IFCXX11_BRACKETS(X) = X
#endif /* RSBP_WANT_CXX11 */

#define RSBEP_INVALID_TIMEV -1 // FIXME !!
#define RSBEP_INVALID_TRANS -1 // FIXME: need RSB_INVALID_TRANS!
#define RSBEP_ITC(TC) ( toupper(TC)=='N' || toupper(TC)=='T' )
#define RSBEP_CTT(TC) ((rsb_trans_t)toupper(TC)) // FIXME !!

#define RSBEP_VI() if(RSBEP_GETENVB("RSBP_VERBOSE_INTERFACE",false)==true) { RSBEP_OS << "RSBP_VERBOSE_INTERFACE  : ", RSBP_PFLOC(RSBEP_OS), RSBEP_OS << "\n"; } /* class Epetra_RsbMatrix verbose interface */

class Epetra_RsbMultiVector : public std::vector<RSBEP_NUMT> 
{
	int nrhs_;
public:
	Epetra_RsbMultiVector(RSBEP_IDXT n, RSBEP_NUMT v, RSBEP_IDXT nrhs=1): std::vector<RSBEP_NUMT> (n*nrhs,v), nrhs_(nrhs) {}
	RSBEP_IDXT NumVectors(void) const { return nrhs_; } 
	RSBEP_IDXT MyLength(void) const { return (this->end() - this->begin()); } 
	RSBEP_IDXT ldX(void) const { return (MyLength())/NumVectors(); } 
	bool conform(const Epetra_RsbMultiVector & x) const
	{
		// TODO: launch exception here.
		return ! ( this->NumVectors() != x.NumVectors() || this->MyLength() != x.MyLength() );
	}
	RSBEP_NUMT sum(void) const { return std::accumulate(this->begin(),this->end(),0.0); } 
	RSBEP_NUMT diff(const Epetra_RsbMultiVector & x) const { RSBEP_NUMT ss = sum(), xs = x.sum(); return ss>xs ? ss-xs:xs-ss; };
	RSBEP_NUMT errnorm(const Epetra_RsbMultiVector & x) const
       	{
	       	RSBEP_NUMT enorm = 0.0;
		for (RSBEP_IDXT i=0;i<MyLength();++i)
			enorm = RSBEP_MAX(enorm, RSBEP_ABS((*this)[i]-x[i]));
	       	return enorm;
	}
	bool diverging (const Epetra_RsbMultiVector & x) const { return errnorm(x) > 1e-9; }
	Epetra_RsbMultiVector & operator=(const RSBEP_NUMT x)
	{
		for (RSBEP_IDXT i=0;i<MyLength();++i)
		{
			(*this)[i]=x;
		}
		return *this;
	}
	Epetra_RsbMultiVector & operator=(const Epetra_RsbMultiVector & x)
	{
		if ( ! conform(x) )
		{
			RSBEP_PRERR_GOTO(ret,"Vector assignment is non conform!\n");
			goto ret;
		}

		for (RSBEP_IDXT i=0;i<MyLength();++i)
		{
			//std::cout << (*this)[i] << " <- " << x[i] << "\n";
			// std::cout << "COPYING\n" << x[i];
			//(*this)[i]=x[i];
		}
		std::copy(x.begin(),x.end(),this->begin());
		//std::copy(this->begin(),this->end(),x.begin());
		//std::cout << "COPYING\n";
ret:
		return *this;// FIXME
	}
#if RSBEP_NO_STUB
	Epetra_RsbMultiVector & operator=(const Epetra_Vector & x)
	{
		// FIXME: no conformity check !
		for (RSBEP_IDXT i=0;i<MyLength();++i)
		{
			(*this)[i]=x[i];
		}
		return *this;
	}
	Epetra_RsbMultiVector & operator=(const Epetra_MultiVector & x)
	{
		// FIXME: no conformity check !
		std::vector<RSBEP_NUMT>::const_iterator xi = this->begin();
		int maxvi = this->ldX();

		for( RSBEP_IDXT vi = 0; vi < maxvi;        vi++ )
		{
			for( RSBEP_IDXT ri = 0; ri < this->NumVectors(); ri++ )
				(*this)[ri*this->ldX()+vi] = x[ri][vi];
		}
		return *this;
	}
	void assign_to(Epetra_MultiVector & x) const
	{
		// FIXME: no conformity check !
		std::vector<RSBEP_NUMT>::const_iterator xi = this->begin();
		int maxvi = this->ldX();

		for( RSBEP_IDXT vi = 0; vi < maxvi;        vi++ )
		{
			for( RSBEP_IDXT ri = 0; ri < this->NumVectors(); ri++ )
					x[ri][vi] = (*this)[ri*this->ldX()+vi];
		}
	}
	bool operator==(const Epetra_Vector & x) const
	{
		// FIXME: no conformity check !
		for (RSBEP_IDXT i=0;i<MyLength();++i)
			if( (*this)[i] != x[i] )
				return false;
		return true;
	}
	void assign_to(Epetra_Vector & x) const
	{
		// FIXME: no conformity check !
		for (RSBEP_IDXT i=0;i<MyLength();++i)
		{
			x[i]=(*this)[i];
		}
	}
#endif /* RSBEP_NO_STUB */
	bool operator==(const Epetra_RsbMultiVector & x) const
	{
		if ( ! conform(x) )
			return false;
		// problem with denormals
		// return std::equal(this->begin(),this->end(),x.begin());
		// return ( 0 == bcmp( &((*this)[0]), &x[0], sizeof(RSBEP_NUMT)*MyLength() ) );
		// better
		return ( this->sum() - x.sum() ) == 0.0;
	}
	bool operator==(const RSBEP_NUMT v) const
	{
		for (RSBEP_IDXT i=0;i<MyLength();++i)
			if( (*this)[i] != v )
				return false;
		return true;
	}
	bool operator!=(const RSBEP_NUMT v) const
	{
		return !(*this == v);
	}
}; /* Epetra_RsbMultiVector  */
std::ostream& operator<<(std::ostream &os, const Epetra_RsbMultiVector & x)
{
	std::vector<RSBEP_NUMT>::const_iterator xi = x.begin();
	RSBEP_IDXT maxvi = RSBEP_MIN( x.ldX(), RSBEP_MAX_PRINTV);

	RSBEP_IDXT extra = x.ldX() - RSBEP_MAX_PRINTV;

	for( RSBEP_IDXT vi = 0; vi < maxvi + (extra == 1 ? 1 : 0);        vi++ )
	{
		for( RSBEP_IDXT ri = 0; ri < x.NumVectors(); ri++ )
			os << xi[ri*x.ldX()+vi] << " ";
		os << "\n";
	}
	if( extra > 1 )
		os << "...<" << extra << " rows follow>...\n";
	return os;
}
template <class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> & X)
{
	os << "[" << ( X.end() - X.begin() ) << "]: ";
	for( auto x : X )
		os << x << " ";
	// os << "\n";
	return os;
}
#if RSBEP_NO_STUB
	/* */
#else /* RSBEP_NO_STUB */
	/* Use stub declarations */
//#define Epetra_Vector std::vector<RSBEP_NUMT> 
#define Epetra_Vector Epetra_RsbMultiVector 
class Epetra_CrsMatrix {
		RSBEP_IDXT nrA_; RSBEP_IDXT *RP_; RSBEP_IDXT *JA_; RSBEP_NUMT*VA_;
	public:
		Epetra_CrsMatrix () {} 
		Epetra_CrsMatrix (RSBEP_IDXT nrA, RSBEP_IDXT *RP, RSBEP_IDXT *JA, RSBEP_NUMT *VA):
			nrA_(nrA), RP_(RP), JA_(JA), VA_(VA)
		{}
		RSBEP_IDXT NumMyRows() const { return nrA_; }
	       	RSBEP_IDXT *ExpertExtractIndexOffset() const { return RP_; }
		RSBEP_IDXT *ExpertExtractIndices() const { return JA_; }
		RSBEP_NUMT * ExpertExtractValues() const { return VA_; }
};

#endif /* RSBEP_NO_STUB */

#if 0
void pv(const std::vector<RSBEP_NUMT> & x, RSBEP_IDXT num = 0)
{
	if(num == 0)
		num = x.end() - x.begin();
	else
	if(num > 0)
		num = RSBEP_MIN(num,x.end() - x.begin());

	if( num > RSBEP_MAX_PRINTV)
		return;

	for( std::vector<RSBEP_NUMT>::const_iterator i = x.begin() ; i < x.begin()+num ; ++i )
		RSBEP_OS << *i << "\n";
	RSBEP_OS << "\n";
}

#if RSBEP_NO_STUB
void pv(const Epetra_Vector & x, RSBEP_IDXT num = 0)
{
	if(num == 0)
		num = x.MyLength();
	else
	if(num > 0)
		num = RSBEP_MIN(num,x.MyLength());

	if( num > RSBEP_MAX_PRINTV)
		return;

	for( RSBEP_IDXT i = 0; i < num ; ++i )
		RSBEP_OS << x[i] << "\n";

	RSBEP_OS << "\n";
}
#endif /* RSBEP_NO_STUB */
#endif /* 0 */

class Epetra_RsbMatrix: public Rsb_Matrix 
			,Epetra_CrsMatrix
{
	private:
	Rsb_Matrix * Rsb_p(void) { return this; } /* cheap style delegation */
	const Rsb_Matrix * Rsb_p(void) const { return this; } /* cheap style delegation */

	public:

#if !RSBEP_NO_STUB
	Epetra_RsbMatrix(rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, const RsbSym sym = IsGen ):
		Rsb_Matrix(nrA, ncA, sym )
       	{
       		/* FIXME: need some status */
	}
#endif /* RSBEP_NO_STUB */

#if !RSBEP_NO_STUB
	Epetra_RsbMatrix(const std::string filename):
		Rsb_Matrix(filename.c_str())
       	{
	}
#endif /* RSBEP_NO_STUB */

       	/* FIXME */
#if RSBEP_NO_STUB
	Epetra_RsbMatrix(/* const*/ Epetra_CrsMatrix& Matrix, /*const Teuchos::ParameterList& List,*/ const RsbSym sym = IsGen ):
		Rsb_Matrix(Matrix.NumMyRows(), Matrix.ExpertExtractIndexOffset().Values(), Matrix.ExpertExtractIndices().Values(),Matrix.ExpertExtractValues(),sym) ,Epetra_CrsMatrix (Matrix)
		//Rsb_Matrix(Source.NumMyRows(), RSBP_NULL, RSBP_NULL,RSBP_NULL,sym) ,Epetra_CrsMatrix (Source)
#else
	Epetra_RsbMatrix(/* const*/ Epetra_CrsMatrix& Source, /*const Teuchos::ParameterList& List,*/ const RsbSym sym = IsGen ):
		Rsb_Matrix(Source.NumMyRows(), Source.ExpertExtractIndexOffset(), Source.ExpertExtractIndices(),Source.ExpertExtractValues(),sym)
		,Epetra_CrsMatrix (Source)
#endif /* RSBEP_NO_STUB */
	{
	}

	RSBEP_ERRT SumIntoMyValues(RSBEP_IDXT MyRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT* Values, const RSBEP_IDXT* Indices)
	{
		RSBEP_VI()
		// Version 1.2-rc3 implements this; versions precedinf it do not.
		const RSBEP_NUMT * VA = Values;
		rsb_coo_idx_t *IA = (rsb_coo_idx_t *) RSBEP_CALLOC(NumEntries,sizeof(RSBEP_IDXT));
		const rsb_coo_idx_t *JA = Indices;
	       	rsb_nnz_idx_t nnz = NumEntries;
	       	rsb_flags_t flags = RSB_FLAG_DUPLICATES_SUM;

		if (!IA)
			RSBEP_PERR_GOTO(err,"");

		for ( rsb_nnz_idx_t i = 0; i < nnz; ++i )
			IA[i] = MyRow;
		Rsb_p()->set_vals(VA, IA, JA, nnz, flags); // FIXME: RSB_FLAG_DUPLICATES_SUM still unsupported on librsb side
		if(IA) RSBEP_FREE ( IA );

		return RSBEP_ERR_OK;
err:
		return RSBEP_ERR_GENERIC_ERROR;
	}

	RSBEP_ERRT ReplaceMyValues(RSBEP_IDXT MyRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT* Values, const RSBEP_IDXT* Indices)
	{
		RSBEP_VI()
		const RSBEP_NUMT * VA = Values;
		rsb_coo_idx_t *IA = (rsb_coo_idx_t *) RSBEP_CALLOC(NumEntries,sizeof(RSBEP_IDXT));
		const rsb_coo_idx_t *JA = Indices;
	       	rsb_nnz_idx_t nnz = NumEntries;
	       	rsb_flags_t flags = RSB_FLAG_DUPLICATES_KEEP_LAST;

		if (!IA)
			RSBEP_PERR_GOTO(err,"");

		for ( rsb_nnz_idx_t i = 0; i < nnz; ++i )
			IA[i] = MyRow;
		Rsb_p()->set_vals(VA, IA, JA, nnz, flags);
		if(IA) RSBEP_FREE ( IA );

		return RSBEP_ERR_OK;
err:
		return RSBEP_ERR_GENERIC_ERROR;
	}

	RSBEP_ERRT ExtractDiagonalCopy(Epetra_Vector& Diagonal) const
       	{
		RSBEP_VI()
		Rsb_p()->get_vec(&Diagonal[0], RSB_EXTF_DIAG);
		return RSBEP_ERR_OK;
	}

  	RSBEP_ERRT ReplaceDiagonalValues(const Epetra_Vector& Diagonal)
	{
		RSBEP_VI()
	       	rsb_flags_t flags = RSB_FLAG_NOFLAGS;

		for ( rsb_nnz_idx_t ii = 0; ii < rows(); ++ii )
			Rsb_p()->set_vals(&Diagonal[ii], &ii, &ii, 1, flags);
		// TODO: error checking is missing

		return RSBEP_ERR_OK;
	}
	/*
RSBEP_IDXT 	ExtractGlobalRowView (RSBEP_IDXT GlobalRow, RSBEP_IDXT &NumEntries, RSBEP_NUMT *&Values, RSBEP_IDXT *&Indices) const
RSBEP_IDXT 	ExtractGlobalRowView (long long GlobalRow, RSBEP_IDXT &NumEntries, RSBEP_NUMT *&Values, long long *&Indices) const
RSBEP_IDXT 	ExtractMyRowView (RSBEP_IDXT MyRow, RSBEP_IDXT &NumEntries, RSBEP_NUMT *&Values, RSBEP_IDXT *&Indices) const
RSBEP_IDXT 	ExtractGlobalRowView (RSBEP_IDXT GlobalRow, RSBEP_IDXT &NumEntries, RSBEP_NUMT *&Values) const
RSBEP_IDXT 	ExtractGlobalRowView (long long GlobalRow, RSBEP_IDXT &NumEntries, RSBEP_NUMT *&Values) const
RSBEP_IDXT 	ExtractMyRowView (RSBEP_IDXT MyRow, RSBEP_IDXT &NumEntries, RSBEP_NUMT *&Values) const
RSBEP_IDXT MakeDataContiguous ()
	*/

	RSBEP_ERRT OptimizeStorage (void)
	{
		RSBEP_VI()
		Rsb_p()->tune_spmm(/*RSBP_NULL,*/ RSBP_NULL, RSBP_NULL, 0, 0.0, RSB_TRANSPOSITION_N, RSBP_NULL, 1, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, RSBP_NULL, 0, RSBP_NULL, RSBP_NULL, 0);
		return RSBEP_ERR_OK; // FIXME
	}

#if RSBEP_NO_STUB
	RSBEP_ERRT ExtractGlobalRowCopy (long long GlobalRow, RSBEP_IDXT Length, RSBEP_IDXT &NumEntries, RSBEP_NUMT *Values, long long *Indices) const
	{
		// TODO: 
		RSBEP_VI()
		return RSBEP_ERR_GENERIC_ERROR;
	}

	RSBEP_ERRT ExtractGlobalRowCopy (long long GlobalRow, RSBEP_IDXT Length, RSBEP_IDXT &NumEntries, RSBEP_NUMT *Values) const
	{
		// TODO: 
		RSBEP_VI()
		return RSBEP_ERR_GENERIC_ERROR;
	}
#endif

	RSBEP_ERRT ExtractMyRowCopy (RSBEP_IDXT MyRow, RSBEP_IDXT Length, RSBEP_IDXT &NumEntries, RSBEP_NUMT *Values, RSBEP_IDXT *Indices) const
	{
		RSBEP_VI()
		rsb_trans_t transA = RSB_TRANSPOSITION_N;
		const RSBEP_NUMT * alphap = RSBP_NULL;
		RSBEP_NUMT * VA = Values;
		rsb_coo_idx_t * IA = RSBP_NULL;
		rsb_coo_idx_t * JA = Indices;
		rsb_coo_idx_t frA = MyRow, lrA = MyRow;
		rsb_nnz_idx_t *rnzp = &NumEntries;
		rsb_flags_t flags = RSB_FLAG_NOFLAGS;

		if(!rnzp)
		{
			// TODO: error checking is missing
			RSBEP_PERR_GOTO(err,"");
		}
		*rnzp = 0;
		Rsb_p()->get_rows_sparse(transA, alphap, RSBP_NULL, RSBP_NULL, RSBP_NULL, frA, lrA, rnzp, flags );
		if( *rnzp > Length )
		{
			// TODO: error checking is missing
			*rnzp = 0;
			RSBEP_PERR_GOTO(err,"");
		}
		Rsb_p()->get_rows_sparse(transA, alphap, VA, IA, JA, frA, lrA, rnzp, flags );
		return RSBEP_ERR_OK;
err:
		return RSBEP_ERR_GENERIC_ERROR;
	}

	RSBEP_ERRT PutScalar (RSBEP_NUMT ScalarConstant)
	{
		RSBEP_VI()
		const RSBEP_NUMT * omegap = &ScalarConstant;
		RSBEP_NUMT dzero = 0.0;

		this->upd_vals(RSB_ELOPF_POW, &dzero);
		this->upd_vals(RSB_ELOPF_MUL, omegap);
		return RSBEP_ERR_OK;
	}

	RSBEP_ERRT Scale (RSBEP_NUMT ScalarConstant)
	{
		RSBEP_VI()
		const RSBEP_NUMT * omegap = &ScalarConstant;

		this->upd_vals(RSB_ELOPF_MUL, omegap);
		return RSBEP_ERR_OK;
	}

	RSBEP_ERRT ExtractMyRowCopy (RSBEP_IDXT MyRow, RSBEP_IDXT Length, RSBEP_IDXT &NumEntries, RSBEP_NUMT *Values) const
	{
		RSBEP_VI()
		return ExtractMyRowCopy (MyRow, Length, NumEntries, Values, RSBP_NULL);
	}

	RSBEP_ERRT InsertMyValues(RSBEP_IDXT Row, RSBEP_IDXT NumEntries, const RSBEP_NUMT* values, const RSBEP_IDXT* Indices)
	{
		RSBEP_VI()
		// TODO: difference with ReplaceMyValues ?
		return ReplaceMyValues(Row, NumEntries, values, Indices);
	}

	RSBEP_ERRT InsertMyValues(RSBEP_IDXT Row, RSBEP_IDXT NumEntries, RSBEP_NUMT* values, RSBEP_IDXT* Indices)
	{
		RSBEP_VI()
		return InsertMyValues(Row, NumEntries, values, Indices);
	}

     	RSBEP_ERRT NumMyRowEntries(RSBEP_IDXT MyRow, RSBEP_IDXT & NumEntries) const
	{
		RSBEP_VI()
		rsb_flags_t flags = RSB_FLAG_NOFLAGS;
		rsb_trans_t transA = RSB_TRANSPOSITION_N;
	
		Rsb_p()->get_rows_sparse(transA, RSBP_NULL, RSBP_NULL, RSBP_NULL, RSBP_NULL, MyRow, MyRow, &NumEntries, flags );
		return RSBEP_ERR_OK;
	}

#if RSBEP_NO_STUB
     	RSBEP_ERRT MaxNumEntries(RSBEP_IDXT MyRow, RSBEP_IDXT & NumEntries) const
	{
		RSBEP_VI()
		rsb_nnz_idx_t rnz = 0;
		rsb_flags_t flags = RSB_FLAG_NOFLAGS;
		rsb_trans_t transA = RSB_TRANSPOSITION_N;

		NumEntries = 0;
		for ( rsb_nnz_idx_t i = 0; i < Rsb_p()->rows(); ++i )
			Rsb_p()->get_rows_sparse(transA, RSBP_NULL, RSBP_NULL, RSBP_NULL, RSBP_NULL, i, i, &rnz, flags ),
			NumEntries = rnz > NumEntries ? rnz : NumEntries;
		return RSBEP_ERR_OK;
	}
#endif

     	RSBEP_ERRT LeftScale(const Epetra_Vector& x)
	{
		RSBEP_VI()
		this->upd_vals(RSB_ELOPF_SCALE_COLS, &x[0]);
		return RSBEP_ERR_OK;
	}

     	RSBEP_ERRT RightScale(const Epetra_Vector& x)
	{
		RSBEP_VI()
		this->upd_vals(RSB_ELOPF_SCALE_ROWS, &x[0]);
		return RSBEP_ERR_OK;
	}

	RSBEP_IDXT NumMyRows() const
       	{
		RSBEP_VI()
     		// NumGlobalRows NumGlobalRows64
	       	return rows();
       	}
	RSBEP_IDXT NumMyCols() const
       	{
		RSBEP_VI()
     		// NumGlobalCols NumGlobalCols64
	       	return cols();
       	}
	RSBEP_IDXT NumMyNonzeros() const
       	{
		RSBEP_VI()
     		// NumGlobalNonzeros NumGlobalNonzeros64 
	       	return nnz(); 
	}
     	bool HasNormInf() 
	{
		RSBEP_VI()
	       	return true;
       	}

  	bool LowerTriangular() const
       	{
		RSBEP_VI()
	       	return RSB_DO_FLAG_HAS(Rsb_p()->rsbflags(),RSB_FLAG_LOWER_TRIANGULAR)?true:false;
       	}
  	bool UpperTriangular() const
       	{
		RSBEP_VI()
	       	return RSB_DO_FLAG_HAS(Rsb_p()->rsbflags(),RSB_FLAG_UPPER_TRIANGULAR)?true:false;
       	}
	RSBEP_IDXT NumMyDiagonals() const
	{
		RSBEP_VI()
     		// TODO: NumGlobalDiagonals NumGlobalDiagonals64
		// Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
		RSBEP_NUMT eii = 0.0;
		RSBEP_IDXT nmd = 0;
		rsb_flags_t flags = RSB_FLAG_NOFLAGS;

		for ( rsb_nnz_idx_t ii = 0; ii < Rsb_p()->rows(); ++ii )
		{
			eii = 0.0;
			this->get_vals(&eii, &ii, &ii, 1, flags);
			if(eii)
				++nmd;
		}
		return nmd;
	}

#if RSBEP_NO_STUB
     	RSBEP_ERRT SetUseTranspose(bool UseTranspose)
		RSBEP_INCOMPLETE_BODY
#endif

     	bool UseTranspose() const
	{
		RSBEP_VI()
		// Returns the current UseTranspose setting. 
		// FIXME: write me 
		return false;
	}

	RSBEP_ERRT FillComplete (bool OptimizeDataStorage=true)
	{
		RSBEP_VI()
		// Signal that data entry is complete. Perform transformations to local index space. 
		Rsb_p()->close();
		if( OptimizeDataStorage )
			return this->OptimizeStorage();
		return RSBEP_ERR_OK;
	}

/*
	RSBEP_IDXT FillComplete (const Epetra_Map &DomainMap, const Epetra_Map &RangeMap, bool OptimizeDataStorage=true)
	{
		RSBEP_VI()
	}
*/
private:
	virtual RSBEP_ERRT _InsertGlobalValues (RSBEP_IDXT GlobalRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT *Values, const RSBEP_IDXT *Indices)
	{
  		RSBEP_ERRT ierr = 0;
		for(RSBEP_IDXT nzi=0;nzi<NumEntries;++nzi)
			Rsb_p()->_add(GlobalRow,Indices[nzi],Values[nzi]);
		ierr = RSBEP_ERR_OK; // FIXME
		return ierr;
	}
public:

	virtual RSBEP_ERRT InsertGlobalValues (RSBEP_IDXT GlobalRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT *Values, const RSBEP_IDXT *Indices)
	{
		RSBEP_VI()
		return _InsertGlobalValues (GlobalRow, NumEntries, Values, Indices);
	}

#if RSBEP_NO_STUB
	virtual RSBEP_ERRT InsertGlobalValues (long long GlobalRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT *Values, const long long *Indices)
		RSBEP_INCOMPLETE_BODY
	virtual RSBEP_ERRT InsertGlobalValues (RSBEP_IDXT GlobalRow, RSBEP_IDXT NumEntries, RSBEP_NUMT *Values, RSBEP_IDXT *Indices)
	{
		RSBEP_VI()
		return _InsertGlobalValues (GlobalRow, NumEntries, Values, Indices);
	}
	virtual RSBEP_ERRT InsertGlobalValues (long long GlobalRow, RSBEP_IDXT NumEntries, RSBEP_NUMT *Values, long long *Indices)
		RSBEP_INCOMPLETE_BODY
	virtual RSBEP_ERRT ReplaceGlobalValues (RSBEP_IDXT GlobalRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT *Values, const RSBEP_IDXT *Indices)
		RSBEP_INCOMPLETE_BODY
	virtual RSBEP_ERRT ReplaceGlobalValues (long long GlobalRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT *Values, const long long *Indices)
		RSBEP_INCOMPLETE_BODY
	virtual RSBEP_ERRT SumIntoGlobalValues (RSBEP_IDXT GlobalRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT *Values, const RSBEP_IDXT *Indices)
		RSBEP_INCOMPLETE_BODY
	virtual RSBEP_ERRT SumIntoGlobalValues (long long GlobalRow, RSBEP_IDXT NumEntries, const RSBEP_NUMT *Values, const long long *Indices)
		RSBEP_INCOMPLETE_BODY
  	
	virtual RSBEP_ERRT Solve1(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const
	{
		RSBEP_VI()
		// FIXME: translate and return error code
		RSBEP_ERRT errval = spsv<rsb_err_t>(&y[0], &x[0], Trans);
		return errval;
	}

	virtual RSBEP_ERRT Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const
	{
		RSBEP_VI()
		// FIXME: translate and return error code
		RSBEP_ERRT errval = spmv<rsb_err_t>(&y[0], &x[0], TransA);
		return errval;
	}
#endif

#if RSBEP_NO_STUB
  	virtual RSBEP_ERRT Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
	{
		RSBEP_VI()
		// FIXME: unfinished
		// FIXME: untested
		// FIXME: error handling missing
  		RSBEP_ERRT ierr = 0;
		RSBEP_IDXT nrhs = X.NumVectors();
		const rsb_trans_t transA = TransA?RSB_TRANSPOSITION_T:RSB_TRANSPOSITION_N;

		if( ! ( X.ConstantStride() && Y.ConstantStride()) )
		{
			for (RSBEP_IDXT ri = 0; ri < nrhs ; ++ri)
				ierr |= spmm<rsb_err_t>(transA,alpha,nrhs,order,&X[ri][0],beta,&Y[ri][0]);
			goto ret;
		}
		else
		{
			RSBEP_IDXT incX = 0, incY = 0;
			incX = X.Stride();
			incY = Y.Stride();
 			// FIXME: ignoring stride !, ..
			ierr |= spmm<rsb_err_t>(transA,alpha,nrhs,order,&X[0][0],beta,&Y[0][0]);
		}

		ierr = RSBEP_ERR_OK; // FIXME
		goto ret;
err:
		ierr = RSBEP_ERR_GENERIC_ERROR;
ret: 		RSBPE_RTE(ierr);
	}
#endif /* RSBEP_NO_STUB */

#if RSBEP_NO_STUB
	// unlike Epetra_CrsMatrix member functions, these are const, and as well is x.
	void GeneralMTV(const RSBEP_NUMT * x, RSBEP_NUMT * y) const
	{
		RSBEP_VI()
		// Note: new, untested
		spmv<rsb_err_t>(y, x, true);
	}
	void GeneralMTM(const RSBEP_NUMT ** X, RSBEP_IDXT LDX, RSBEP_NUMT ** Y, RSBEP_IDXT LDY, RSBEP_IDXT NumVectors) const
	{
		RSBEP_VI()
		// Note: new, untested
		const RSBEP_NUMT alpha = 1, beta = 0;
		const RSBEP_IDXT incX = 1, incY = 1;
		rsb_trans_t transA = RSB_TRANSPOSITION_T;
		if (LDX!=0 && LDY!=0)
			Rsb_p()->spmm<rsb_err_t>(transA, &alpha, NumVectors, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, X[0], LDX, &beta, Y[0], LDY);
		else
			for (RSBEP_IDXT k=0; k < NumVectors; k++)
				Rsb_p()->spmv<rsb_err_t>(transA, &alpha, X[k], incX, &beta, Y[k], incY);
	}
	void GeneralMM(const RSBEP_NUMT ** X, RSBEP_IDXT LDX, RSBEP_NUMT ** Y, RSBEP_IDXT LDY, RSBEP_IDXT NumVectors) const
	{
		RSBEP_VI()
		// Note: new, untested
		const RSBEP_NUMT alpha = 1, beta = 0;
		const RSBEP_IDXT incX = 1, incY = 1;
		rsb_trans_t transA = RSB_TRANSPOSITION_N;
		if (LDX!=0 && LDY!=0)
			Rsb_p()->spmm<rsb_err_t>(transA, &alpha, NumVectors, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, X[0], LDX, &beta, Y[0], LDY);
		else
			for (RSBEP_IDXT k=0; k < NumVectors; k++)
				Rsb_p()->spmv<rsb_err_t>(transA, &alpha, X[k], incX, &beta, Y[k], incY);
		return;
	}
	void GeneralMV(const RSBEP_NUMT * x, RSBEP_NUMT * y) const
	{
		RSBEP_VI()
		// Note: new, untested
		spmv<rsb_err_t>(y, x, false);
	}
	void GeneralSV(bool Upper, bool Trans, bool UnitDiagonal, RSBEP_NUMT * xp, RSBEP_NUMT * yp)  const
	{
		RSBEP_VI()
    		const RSBEP_NUMT alpha = 1;
		const RSBEP_IDXT incX = 1, incY = 1;
		Rsb_p()->spsv<rsb_err_t>(Trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N, &alpha, xp, incX, yp, incY);
	}
	void GeneralSM(bool Upper, bool Trans, bool UnitDiagonal, RSBEP_NUMT ** Xp, RSBEP_IDXT LDX, RSBEP_NUMT ** Yp, RSBEP_IDXT LDY, RSBEP_IDXT NumVectors) const
	{
		RSBEP_VI()
    		/* FIXME: UnitDiagonal and Upper args ignored */
		const RSBEP_NUMT alpha = 1, beta = 0.0;
		const RSBEP_IDXT incX = 1, incY = 1;
		rsb_trans_t transA = Trans ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;
		if (LDX!=0 && LDY!=0)
			Rsb_p()->spsm<rsb_err_t>(transA, &alpha, NumVectors, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, &beta, Xp[0], LDX, Yp[0], LDY);
		else
			for (RSBEP_IDXT k=0; k < NumVectors; k++)
				Rsb_p()->spsv<rsb_err_t>(transA, &alpha, Xp[k], incX, Yp[k], incY);
	}
#endif /* RSBEP_NO_STUB */

	RSBEP_ERRT _print(std::ostream &os) const
	{
		// TODO: shall maybe produce a string instead.
  		RSBEP_ERRT ierr = 0;
		os << "nr:" << rows();
		os << " nc:" << cols();
		os << " nnz:" << nnz();
		os << " normOne:" << normOne();
		os << " normInf:" << normInf();
		os << "\n";
		if( nnz() < RSBEP_MAX_PRINT )
			file_save(RSBP_NULL);
		// os << "\n";
		RSBPE_RTE(ierr);
	}

	bool isSym(void) const
	{
		RSBEP_VI()
		return (RSB_DO_FLAG_HAS(Rsb_p()->get_info_rsb_flags_t(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T),RSB_FLAG_SYMMETRIC));
	}

	RSBEP_NUMT NormInf(void) const
	{
		RSBEP_VI()
		return Rsb_p()->normInf();
	}

	RSBEP_NUMT NormOne(void) const
	{
		RSBEP_VI()
		return Rsb_p()->normOne();
	}
#if 0
  	virtual RSBEP_ERRT Solve1(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
		RSBEP_INCOMPLETE_BODY
  	virtual RSBEP_ERRT Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
		RSBEP_INCOMPLETE_BODY
#endif
}; /* Epetra_RsbMatrix */

std::ostream& operator<<(std::ostream &os,const Epetra_RsbMatrix & s)
{
	s._print(os);
	return os;
}

std::ostream& operator<<(std::ostream &os,const Rsb_Matrix & s)
{
	os << s._info() << std::endl;
	os << "normOne:" << s.normOne();
	os << " normInf:" << s.normInf();
	os << "\n";
	if( s.nnz() < RSBEP_MAX_PRINT )
		s.file_save(RSBP_NULL);
	return os;
}

#if !RSBP_WANT_INLINE_INIT
	RsbLib rsb_lib;
#endif /* RSBP_WANT_INLINE_INIT */

int Rsb_Matrix_test_ms( const Rsb_Matrix & Rsb_A, RSBEP_IDXT nrhs, rsb_trans_t transA, const RSBEP_NUMT alphav, Epetra_RsbMultiVector * Rsb_Yp, std::string mcs, bool wq)
{
	// relies on defmultbeta = 0.0
	// only makes sense on triangular matrices
	// TODO: a triangularity + diagonal non singularity check
	int perrval = EXIT_FAILURE;
	RSBEP_IDXT nr(Rsb_A.rows());
	RSBEP_IDXT nc(Rsb_A.cols());
	RSBEP_IDXT tnr = ( transA == RSB_TRANSPOSITION_N ) ? nr : nc;
	RSBEP_IDXT tnc = ( transA == RSB_TRANSPOSITION_N ) ? nc : nr;
	const Epetra_RsbMultiVector x(tnc,2*1.0,nrhs); /* e.g. 1.1 would be unstable here */
	Epetra_RsbMultiVector y(tnr,-1*1.0,nrhs);
	Epetra_RsbMultiVector z(tnc,-1*0.0,nrhs);
	bool wv = RSBEP_GETENVB("RSBP_VERBOSE_TEST",false);
	std::string msg,lmsg;
	bool islowtri = RSB_DO_FLAG_HAS(Rsb_A.rsbflags(),RSB_FLAG_LOWER_TRIANGULAR)?true:false;
	bool btransA = ( transA == RSB_TRANSPOSITION_N ) ? false : true;
	bool issym = (RSB_DO_FLAG_HAS(Rsb_A.get_info_rsb_flags_t(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T),RSB_FLAG_SYMMETRIC) );
	const RSBEP_NUMT beta = 0.0; // FIXME: shall try different alpha's
	const RSBEP_NUMT alphai = 1.0/alphav;
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
	const RSBEP_NUMT alpha = 1.0;

	islowtri = islowtri && (Rsb_A.rows()==Rsb_A.cols()); // FIXME: temporary

       	msg += __func__;
      	msg += ": ";
      	msg += mcs;
      	msg += ": ";
	msg += std::to_string(Rsb_A.rows()) + ("x" + std::to_string(Rsb_A.cols()));
	msg += " nnz:" + std::to_string(Rsb_A.nnz());
	msg += " sym:" + std::to_string(issym);
	msg += " tri:" + std::to_string(islowtri);
	msg += " nrhs:" + std::to_string(nrhs);
	msg += " trans:" + std::to_string(btransA);
	msg += " alpha:" + std::to_string(alphav);
	/*
	msg += std::to_string(Rsb_A.normInf());
	msg += " normInf  ";
	msg += std::to_string(Rsb_A.normOne());
	msg += " normOne  "; 
	*/
	msg += ( std::string(" onrm:") + std::to_string(Rsb_A.normOne()));
	msg += ( std::string(" inrm:") + std::to_string(Rsb_A.normInf()));
	msg += " blk:" + std::to_string(Rsb_A.blocks());

	lmsg = msg;
	lmsg += Rsb_A._info() + "\n";

	if (wv) RSBEP_OS << __func__ << ":\n";
	if (wv) RSBEP_OS << " x :\n" << x;
	if (wv) RSBEP_OS << " y :\n" << y;

	if(alphav!=1.0)
	{
		if(nrhs==1)
		{
			RSBEP_ERRGOTO(Rsb_A.spmv<rsb_err_t>(transA,alphav,&x[0],beta,&y[0]),err);
			if( Rsb_Yp )
			       	*Rsb_Yp = y;
			if( islowtri )
			{
				RSBEP_ERRGOTO(Rsb_A.spsv<rsb_err_t>(transA,alphai,&y[0],&y[0]),err);
				z=y;
				RSBEP_ERRGOTO(Rsb_A.spmv<rsb_err_t>(transA,alphav,&z[0],beta,&y[0]),err);
				RSBEP_ERRGOTO(Rsb_A.spsv<rsb_err_t>(transA,alphai,&y[0],&y[0]),err);
			}
		}
		else
		{

			RSBEP_ERRGOTO(Rsb_A.spmm<rsb_err_t>(transA, alphav, nrhs, order, &x[0], beta, &y[0]), err);
			if( Rsb_Yp )
			       	*Rsb_Yp = y;
			if(islowtri )
			{
				RSBEP_ERRGOTO(Rsb_A.spsm<rsb_err_t>(transA,alphai,nrhs,&y[0],&y[0]),err);
				z=y;
				RSBEP_ERRGOTO(Rsb_A.spmm<rsb_err_t>(transA,alphav, nrhs, order, &z[0], beta, &y[0]), err);
				RSBEP_ERRGOTO(Rsb_A.spsm<rsb_err_t>(transA,alphai,nrhs,&y[0],&y[0]),err);
			}
		}
		goto postop;
	}

	if(nrhs==1)
	{
		RSBEP_ERRGOTO(Rsb_A.spmv<rsb_err_t>(&y[0],&x[0],btransA),err);
		if( Rsb_Yp )
		       	*Rsb_Yp = y;
		if( islowtri )
		{
			RSBEP_ERRGOTO(Rsb_A.spsv<rsb_err_t>(&y[0],btransA),err);
			z=y;
			RSBEP_ERRGOTO(Rsb_A.spmv<rsb_err_t>(&y[0],&z[0],btransA),err);
			RSBEP_ERRGOTO(Rsb_A.spsv<rsb_err_t>(&y[0],btransA),err);
		}
	}
	else
	{
		RSBEP_ERRGOTO(Rsb_A.spmm<rsb_err_t>(transA,alpha,nrhs,order,&x[0],beta,&y[0]),err);
		if( Rsb_Yp )
		       	*Rsb_Yp = y;
		if(islowtri )
		{
			RSBEP_ERRGOTO(Rsb_A.spsm<rsb_err_t>(&y[0],nrhs,btransA),err);
			z=y;
			RSBEP_ERRGOTO(Rsb_A.spmm<rsb_err_t>(transA,alpha,nrhs,order,&z[0],beta,&y[0]),err);
			RSBEP_ERRGOTO(Rsb_A.spsm<rsb_err_t>(&y[0],nrhs,btransA),err);
		}
	}
postop:
	if(!islowtri)
	{
		perrval = EXIT_SUCCESS;
		goto err;
	}
	if (wv) RSBEP_OS << " z :\n" << z;
	if (wv) RSBEP_OS << " y :\n" << y;
	if (wv) RSBEP_OS << "Rsb_A:\n" << Rsb_A;

	if( ! x.diverging(y)  )
	{
		perrval = EXIT_SUCCESS;
		if (wv) RSBEP_OS << "OK ! Seems like multiply and solve is working (diff of sums:" << x.diff(y) << ")!\n";
	}
	else
	{
		RSBEP_OS << "NO ! Seems like multiply and solve is NOT working !\n";
		if (wv) RSBEP_OS << " expected approximated equality here:\n";
		if (wv) RSBEP_OS << " x :\n" << x;
		if (wv) RSBEP_OS << " y :\n" << y;
		RSBEP_OS << "sum(x)-sum(y): " << x.diff(y) << "\n";
		RSBEP_OS << "errnorm(x,y): " << x.errnorm(y) << "\n";
	}
err:
	if ( perrval == EXIT_FAILURE )
		RSBEP_TEST_MSG("ERROR: " + lmsg);
	else
		if(!wq)
			RSBEP_TEST_MSG(("OK: " + msg) + ( islowtri ? ( " errnorm:" + std::to_string(x.errnorm(y)) ) : std::string("") ) );

	return perrval;
} /* Rsb_Matrix_test_ms */

int Rsb_Matrix_test_ms_mnrhs( const Rsb_Matrix & Rsb_A, RSBEP_IDXT nrhs, rsb_trans_t transA, std::string mcs )
{
	bool wq = ( RSBEP_GETENVB("RSBP_QUIET",false) );
	int perrval = EXIT_FAILURE;
	rsb_flags_t flags = RSB_FLAG_NOFLAGS;
	bool issym = (RSB_DO_FLAG_HAS(Rsb_A.get_info_rsb_flags_t(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T),RSB_FLAG_SYMMETRIC) );

	if (issym) 
		flags |= RSB_FLAG_LOWER_SYMMETRIC;
	if (RSB_DO_FLAG_HAS(Rsb_A.get_info_rsb_flags_t(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T),RSB_FLAG_LOWER_TRIANGULAR) )
		flags |= RSB_FLAG_LOWER_TRIANGULAR;
	Rsb_Matrix Coo_A( Rsb_A, false, flags|RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS ); // a COO clone
	Rsb_Matrix Csr_A( Rsb_A, false, flags|RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS ); // a CRS clone

	if( Rsb_A != Coo_A )
		RSBEP_PRERR_GOTO(ret,"");
	if( Rsb_A != Csr_A )
		RSBEP_PRERR_GOTO(ret,"");

	if( Coo_A.blocks() != 1 )
		RSBEP_PRERR_GOTO(ret,"");

	if( Coo_A.blocks() != Csr_A.blocks() )
		RSBEP_PRERR_GOTO(ret,"");
//#if (RSB_LIBRSB_VER >= 10300)
	/* FIXME: this will easily break :-) */
	/*
RSBEP_OS << " " << Csr_A.normOne() << "\n";
RSBEP_OS << " " << Csr_A.normInf() << "\n";
RSBEP_OS << " " << Rsb_A.normOne() << "\n";
RSBEP_OS << " " << Rsb_A.normInf() << "\n";
RSBEP_OS << " " << Coo_A.normOne() << "\n";
RSBEP_OS << " " << Coo_A.normInf() << "\n";
*/
	if (issym) 
	{
		if( Coo_A.normInf() != Coo_A.normOne() )
			RSBEP_PRERR_GOTO(ret,"");
		if( Csr_A.normInf() != Csr_A.normOne() )
			RSBEP_PRERR_GOTO(ret,"");
		if( Rsb_A.normInf() != Rsb_A.normOne() )
			RSBEP_PRERR_GOTO(ret,"");
	}

	// TODO: in general, shall limit these checks to when Rsb_A._is_integer() or when no roundoff can occur.
	if( Rsb_A.normOne() != Coo_A.normOne() )
		RSBEP_PRERR_GOTO(ret,"");
	if( Rsb_A.normOne() != Csr_A.normOne() )
		RSBEP_PRERR_GOTO(ret,"");
	if( Rsb_A.normInf() != Coo_A.normInf() )
		RSBEP_PRERR_GOTO(ret,"");
	if( Rsb_A.normInf() != Csr_A.normInf() )
		RSBEP_PRERR_GOTO(ret,"");
//#endif /* (RSB_LIBRSB_VER >= 10300) */

	for (RSBEP_NUMT alpham =  1 ; alpham <= 2 ; ++alpham )
	for (RSBEP_NUMT alphas = -1 ; alphas <= 1 ; ++alphas )
	for (RSBEP_IDXT nrhsi = 1 ; nrhsi <= nrhs ; ++nrhsi)
	if(alphas)
	{
		RSBEP_NUMT alpha = alpham*alphas;
		RSBEP_NUMT ival = -1.0;
		RSBEP_IDXT nr(Rsb_A.rows());
		RSBEP_IDXT nc(Rsb_A.cols());
		RSBEP_IDXT tnr = ( transA == RSB_TRANSPOSITION_N ) ? nr : nc;
		Epetra_RsbMultiVector Rsb_Y(tnr,ival,nrhsi);
		Epetra_RsbMultiVector Coo_Y(tnr,ival,nrhsi);
		Epetra_RsbMultiVector Csr_Y(tnr,ival,nrhsi);

		if( (perrval = Rsb_Matrix_test_ms(Rsb_A,nrhsi,transA,alpha,&Rsb_Y,mcs,wq)) != EXIT_SUCCESS)
			RSBEP_PRERR_GOTO(ret,"");
		perrval = EXIT_FAILURE; /* ready for further tests */
		if( Rsb_Y == ival )
			RSBEP_PRERR_GOTO(ret,"");
		if( (perrval = Rsb_Matrix_test_ms(Csr_A,nrhsi,transA,alpha,&Csr_Y,mcs,wq)) != EXIT_SUCCESS)
			RSBEP_PRERR_GOTO(ret,"");
		perrval = EXIT_FAILURE; /* ready for further tests */
		if( Csr_Y == ival )
			RSBEP_PRERR_GOTO(ret,"");
		if( Csr_Y.diverging(Rsb_Y) )
			RSBEP_PRERR_GOTO(ret,"");
		if( (perrval = Rsb_Matrix_test_ms(Coo_A,nrhsi,transA,alpha,&Coo_Y,mcs,wq)) != EXIT_SUCCESS)
			RSBEP_PRERR_GOTO(ret,"");
		perrval = EXIT_FAILURE; /* ready for further tests */
		if( Coo_Y == ival )
			RSBEP_PRERR_GOTO(ret,"");
		if( Coo_Y.diverging(Rsb_Y) )
			RSBEP_PRERR_GOTO(ret,"");
		perrval = EXIT_FAILURE; /* ready for further tests */
	}
	perrval = EXIT_SUCCESS;
ret:
	return perrval;

} /* Rsb_Matrix_test_ms_mnrhs */

int Rsb_Matrix_test_multimatrix_ms_mnrhs( RSBEP_IDXT nrhs, const char * const transs )
{
	/* A test matrix generator. */
	int perrval = EXIT_FAILURE;
	RSBEP_IDXT dim, ldim = 1, udim = RSBEP_GETENVI("RSBP_DMAX", 20);
	//bool wv = RSBEP_GETENVB("RSBP_VERBOSE",false);
	rsb_trans_t transA = RSBEP_INVALID_TRANS;
	const char * transp = RSBP_NULL;
bool dnsm = false; // debug non square matrices // FIXME: temporary
	// if(wv)
	       	RSBEP_OS << __func__ << "\n";

	if( RSBEP_GETENVB("RSBP_TRYDIAG",true) )
	for ( transp = transs, transA = RSBEP_CTT(*transp); RSBEP_ITC(*transp) ; ++transp, transA = RSBEP_CTT(*transp) )
	for (dim = ldim; dim < udim; ++dim )
	{
		// pure diagonal
		std::vector<RSBEP_NUMT> VA(dim  ,0.0);
		std::vector<RSBEP_IDXT> RP(dim+1,0);
		std::vector<RSBEP_IDXT> JA(dim  ,0);
		RSBEP_IDXT nzi = 0;
		RSBEP_IDXT ri = 0;

		RP[ri]=0;
		for ( ri = 0; ri < dim ; ++ri )
		{
			RSBEP_NUMT rsum = (1+ri) *  1.0;
			VA[nzi]=rsum;
			JA[nzi]=ri;
			nzi++;
			RP[ri+1]=nzi;
		}
		Rsb_Matrix Rsb_T(dim,&RP[0],&JA[0],&VA[0],Epetra_RsbMatrix::IsTri);
		if ( ( perrval = Rsb_Matrix_test_ms_mnrhs( Rsb_T, nrhs, transA, "PD" ) ) != EXIT_SUCCESS )
			RSBEP_PRERR_GOTO(err,"")
	}

	if( RSBEP_GETENVB("RSBP_TRYVERT",true) )
	for ( transp = transs, transA = RSBEP_CTT(*transp); RSBEP_ITC(*transp) ; ++transp, transA = RSBEP_CTT(*transp) )
	for (dim = ldim; dim < udim; ++dim )
	for (RSBEP_IDXT cci = 0; cci < dim; ++cci ) // constant column index
	{
		// pure vertical
		// please note that CSR constructor here uses maximal column index to extrapolate width
		std::vector<RSBEP_NUMT> VA(dim  ,0.0);
		std::vector<RSBEP_IDXT> RP(dim+1,0);
		std::vector<RSBEP_IDXT> JA(dim  ,0);
		RSBEP_IDXT nzi = 0;
		RSBEP_IDXT ri = 0, ci = 0;

		RP[ri]=0;
		for ( ri = 0; ri < dim ; ++ri )
			RP[ri+1]=RP[ri]+1;
		for ( ci = 0; ci < dim ; ++ci )
		{
			VA[nzi]= ci+1;
			JA[nzi]=cci+0;
			nzi++;
		}
		Rsb_Matrix Rsb_T(dim,&RP[0],&JA[0],&VA[0]);
		if ( ( perrval = Rsb_Matrix_test_ms_mnrhs( Rsb_T, nrhs, transA, "PV" ) ) != EXIT_SUCCESS )
			RSBEP_PRERR_GOTO(err,"")
	}

	if( RSBEP_GETENVB("RSBP_TRYHORZ",true) )
	for ( transp = transs, transA = RSBEP_CTT(*transp); RSBEP_ITC(*transp) ; ++transp, transA = RSBEP_CTT(*transp) )
	for (dim = ldim; dim < udim; ++dim )
	for (RSBEP_IDXT cri = 0; cri < dim; ++cri ) // constant row index
	{
		// pure horizontal
		std::vector<RSBEP_NUMT> VA(dim  ,0.0);
		std::vector<RSBEP_IDXT> RP(dim+1,0);
		std::vector<RSBEP_IDXT> JA(dim  ,0);
		RSBEP_IDXT nzi = 0;
		RSBEP_IDXT ri = 0, ci = 0;

		RP[ri]=0;
		for ( ri = 0; ri < cri ; ++ri )
			RP[ri+1]=0;
		for ( ci = 0; ci < dim ; ++ci )
		{
			VA[nzi]=ci+1;
			JA[nzi]=ci+0;
			nzi++;
		}
		for ( ri = cri; ri < dim ; ++ri )
		{
			RP[ri+1]=nzi;
		}
		Rsb_Matrix Rsb_T(dim,&RP[0],&JA[0],&VA[0]);
		if ( ( perrval = Rsb_Matrix_test_ms_mnrhs( Rsb_T, nrhs, transA, "PH" ) ) != EXIT_SUCCESS )
			RSBEP_PRERR_GOTO(err,"")
	}

	if( RSBEP_GETENVB("RSBP_TRYTRI",true) )
	for ( transp = transs, transA = RSBEP_CTT(*transp); RSBEP_ITC(*transp) ; ++transp, transA = RSBEP_CTT(*transp) )
	for ( int sym = 0 ; sym <= 3 ; ++ sym )
	for (dim = ldim; dim < udim; ++dim )
	{
		// lower triangle
		// RSBEP_OS << "Testing for " << dim << " equations.\n" ;
		const RSBEP_IDXT ff = 2; // form factor
		RSBEP_IDXT rmf = (sym == 2) ? ff : 1;
		RSBEP_IDXT cmf = (sym == 3) ? ff : 1;
		RSBEP_IDXT nr = dim * rmf;
		RSBEP_IDXT nc = dim * cmf;
		RSBEP_IDXT nnzT = rmf * ((dim*(dim+1))/2);
		std::vector<RSBEP_NUMT> VA(nnzT,0.0);
		std::vector<RSBEP_IDXT> RP(nr+1,0);
		std::vector<RSBEP_IDXT> JA(nnzT,0);
		RSBEP_IDXT nzi = 0;
		RSBEP_IDXT ri = 0, ci = 0;

		RP[ri]=0;
		for ( ri = 0; ri < nr ; ++ri )
		{
			RSBEP_NUMT rsum = - 1.0, val = 0.0;
//std::cout << "a "<< (cmf) << "\n";
//std::cout << "a "<< (nc) << "\n";
//std::cout << "a "<< (cmf*(ri+1))%nc << "\n";
			// for ( ci = cmf-1; ci <= cmf-1+cmf*ri ; ci+= cmf )
			for ( ci = 0; ci <= ri%nc ; ci++ )
			{
				val = (ci+1)+(ri+1);
				VA[nzi] = val;
				val *= (ci%2?-1:1);
				rsum += val;
				//JA[nzi]=ci;
				JA[nzi]=(cmf-1+cmf*ci)%nc;
//std::cout << "ci "<< (ci) << ":" << JA[nzi]<< "\n";
				nzi++;
			}
			assert(nzi>0);
			VA[nzi-1]=-rsum+val;
			RP[ri+1]=nzi;
		}
			assert(nzi==nnzT);
if(dnsm)
{
std::cout << "nz:" << nnzT << "\n";
std::cout << "nzi:" << nzi << "\n";
std::cout << "nr:" << nr << "\n";
std::cout << "nc:" << nc << "\n";
std::cout << "RP:" << RP << "\n";
std::cout << "JA:" << JA << "\n";
std::cout << "VA:" << (VA) << "\n";
}
		Rsb_Matrix Rsb_T(nr,&RP[0],&JA[0],&VA[0], (sym == 1 ? Epetra_RsbMatrix::IsSym : ( sym == 0 ? Epetra_RsbMatrix::IsTri : Epetra_RsbMatrix::IsGen) ));
if(dnsm)
std::cout << Rsb_T << "\n";
		if ( ( perrval = Rsb_Matrix_test_ms_mnrhs( Rsb_T, nrhs, transA, "LT" ) ) != EXIT_SUCCESS )
			RSBEP_PRERR_GOTO(err,"")
	}

	perrval = EXIT_SUCCESS;
err:
	return perrval;
} /* Rsb_Matrix_test_multimatrix_ms_mnrhs */

int Rsb_Matrix_file_demo(const rsb_char_t * filename)
{
	int perrval = EXIT_FAILURE;
	RSBEP_OS << "BEGIN\n";
	RSBEP_IDXT nrhs = RSBEP_GETENVI("RSBP_NRHS", 9);
	//const rsb_char_t * filename = RSBEP_GETENVS("RSBP_FILENAME","A.mtx");
	bool isTri = RSBEP_GETENVB("RSBP_ISTRI",true); // attempts building a triangular matrix
	rsb_lib.meminfo();
	Rsb_Matrix Rsb_A( filename, isTri ? Epetra_RsbMatrix::IsTri : Epetra_RsbMatrix::IsGen );
	bool wv = RSBEP_GETENVB("RSBP_VERBOSE",false);
	RSBEP_IDXT nr(Rsb_A.rows());
	RSBEP_IDXT nc(Rsb_A.cols());
	Epetra_RsbMultiVector x(nc,1*1.1,nrhs);
	Epetra_RsbMultiVector y(nr,0*2.2,nrhs);
	Epetra_RsbMultiVector z(nr,0*2.2     );
	Rsb_Matrix Rsb_B( Rsb_A, false ); // a clone// TODO: also test optimized instances !
	rsb_real_t sf = RSBEP_INVALID_TIMEV;
	std::string sfb RSBEP_IFCXX11_BRACKETS("Tuned with speedup factor of ");
	bool btransA = RSBEP_GETENVB("RSBP_TRANS", false);
	rsb_trans_t transA = RSBEP_INVALID_TRANS;
	const char * transs = RSBEP_GETENVS("RSBP_TRANSS", "NT");
	const char * transp = RSBP_NULL;
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
	const RSBEP_NUMT alpha = 1.0;
	const RSBEP_NUMT beta = 0.0;

	if(!*transs)
		transs = "nt";
	for ( transp = transs, transA = RSBEP_CTT(*transp); (*transp) ; ++transp, transA = RSBEP_CTT(*transp) )
	if(transA != RSB_TRANSPOSITION_N && transA != RSB_TRANSPOSITION_T)
	{
		RSBEP_OS << "RSBP_TRANSS string must be of the form {NT}*!\n";
		goto ret;
	}

	if (RSB_DO_FLAG_HAS(Rsb_A.get_info_rsb_flags_t(RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T),RSB_FLAG_SYMMETRIC) ) // because the matrix was e.g. symmetric 
		isTri = false;

	transA = btransA ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N; /* after transs stuff */

	if(isTri)
	{
		if( (perrval = Rsb_Matrix_test_ms_mnrhs(Rsb_A,nrhs,transA,filename)) != EXIT_SUCCESS)
			RSBEP_PRERR_GOTO(ret,"")
		else
			perrval = EXIT_FAILURE; /* ready for further tests */
	}

	if(wv) RSBEP_OS << "Rsb_A:\n" << Rsb_A;
	rsb_lib.meminfo();

	if(wv) RSBEP_OS << "Rsb_A:\n" << Rsb_A;
	if(wv) RSBEP_OS << "x:\n" << x;

	Rsb_B = Rsb_A;
	Rsb_B.spmv(&z[0], &x[0],btransA);
	if(wv) RSBEP_OS << "z:\n" << z;
	RSBEP_OS << "" << Rsb_B._info() << std::endl;
	Rsb_B.tune_spmm(/*RSBP_NULL,*/ &sf, RSBP_NULL, 0, 0.0, transA, RSBP_NULL, 1, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, RSBP_NULL, 0, RSBP_NULL, RSBP_NULL, 0);
	RSBEP_OS << sfb << sf << ":\n" << Rsb_B._info() << std::endl;

	if(isTri)
	{
		Rsb_B = Rsb_A;
		Rsb_B.spsv<rsb_err_t>(&z[0],&z[0],btransA);
		if(wv) RSBEP_OS << "z:\n" << z;
		RSBEP_OS << "" << Rsb_B._info() << std::endl;
		Rsb_B.tune_spsm(/*RSBP_NULL,*/ &sf, RSBP_NULL, 0, 0.0, transA, RSBP_NULL, 1, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, RSBP_NULL, 0, RSBP_NULL, RSBP_NULL, 0);
		RSBEP_OS << sfb << sf << ":\n" << Rsb_B._info() << std::endl;
	}

	Rsb_B = Rsb_A;
	if(wv) RSBEP_OS << "y:\n" << y;
	Rsb_B.spmm<rsb_err_t>(transA,alpha,nrhs,order,&x[0],beta,&y[0]);
	if(wv) RSBEP_OS << "y:\n" << y;
	RSBEP_OS << "" << Rsb_B._info() << std::endl;
	Rsb_B.tune_spmm(/*RSBP_NULL,*/ &sf, RSBP_NULL, 0, 0.0, transA, RSBP_NULL, nrhs, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, RSBP_NULL, 0, RSBP_NULL, RSBP_NULL, 0);
	RSBEP_OS << sfb << sf << ":\n" << Rsb_B._info() << std::endl;
	if(isTri)
	{
		Rsb_B = Rsb_A;
		Rsb_B.spsm<rsb_err_t>(&y[0],&y[0],nrhs,btransA);
		if(wv) RSBEP_OS << "y:\n" << y;
		RSBEP_OS << "" << Rsb_B._info() << std::endl;
		Rsb_B.tune_spsm(/*RSBP_NULL,*/ &sf, RSBP_NULL, 0, 0.0, transA, RSBP_NULL, nrhs, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER, RSBP_NULL, 0, RSBP_NULL, RSBP_NULL, 0);
		RSBEP_OS << sfb << sf << ":\n" << Rsb_B._info() << std::endl;
	}

	RSBEP_OS << "END\n";
	perrval = EXIT_SUCCESS;
ret:
	return perrval;
} /* Rsb_Matrix_file_demo */

int Rsb_Matrix_demo(void)
{
	int perrval = EXIT_FAILURE;
	RSBEP_OS << "BEGIN\n";
	RSBEP_IDXT nrhs = RSBEP_GETENVI("RSBP_NRHS", 9);
	rsb_lib.meminfo();
	rsb_trans_t transA = RSBEP_INVALID_TRANS;
	const char * transs = RSBEP_GETENVS("RSBP_TRANSS", "NT");
	const char * transp;

	if(!*transs)
		transs = "nt";
	for ( transp = transs, transA = RSBEP_CTT(*transp); (*transp) ; ++transp, transA = RSBEP_CTT(*transp) )
	if(transA != RSB_TRANSPOSITION_N && transA != RSB_TRANSPOSITION_T)
	{
		RSBEP_OS << "RSBP_TRANSS string must be of the form {NT}*!\n";
		RSBEP_PRERR_GOTO(ret,"")
	}

	if( (perrval = Rsb_Matrix_test_multimatrix_ms_mnrhs(nrhs, transs)) != EXIT_SUCCESS)
		RSBEP_PRERR_GOTO(ret,"")
	else
		perrval = EXIT_FAILURE; /* ready for further tests */

	perrval = EXIT_SUCCESS;
ret:
	return perrval;
} /* Rsb_Matrix_demo */

#if defined(RSBEP_NO_STUB) && defined(EPETRA_MPI)
int Epetra_RsbMatrix_demo(int argc, char *argv[])
#else
int Epetra_RsbMatrix_demo()
#endif
{
	std::ostream & os = RSBEP_OS;
	int perrval = EXIT_FAILURE;
	const RSBEP_IDXT n = RSBEP_GETENVI("RSBP_N", 3);
	const RSBEP_IDXT nrhs = RSBEP_GETENVI("RSBP_NRHS", 2);
	const bool tune = RSBEP_GETENVI("RSBP_TUNE", false);
	RSBEP_IDXT MyRow = n - 1;
	RSBEP_IDXT NumEntries = n;
	RSBEP_NUMT xiv = 1.1, yiv = 99.;
	const rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
	const RSBEP_NUMT alpha = 1.0;
	const RSBEP_NUMT beta = 0.0;
#if RSBEP_NO_STUB
#ifdef EPETRA_MPI
  	MPI_Init(&argc,&argv);
#ifdef HAVE_MPI
  	Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else /* HAVE_MPI */
  	Epetra_SerialComm Comm;
#endif /* HAVE_MPI */
#else /* EPETRA_MPI */
  	Epetra_SerialComm Comm;
#endif /* EPETRA_MPI */
   	const int MyPID = Comm.MyPID();
   	const int NumProc = Comm.NumProc(); 
   	const bool verbose = (MyPID==0);
   	const RSBEP_IDXT NumGlobalElements = n;
   	const Epetra_Map Map(NumGlobalElements, 0, Comm);
   	const RSBEP_IDXT NumMyElements = Map.NumMyElements();
   	std::vector<RSBEP_IDXT> NumNz(NumMyElements);
	bool egal = true;
	bool isUpper = egal;
	bool unitDiagonal = egal;

	if ( NumProc > 1 )
	{
		os << "Only one MPI rank is supported !\n";
		return perrval;
	}

   	for (RSBEP_IDXT i=0; i<NumMyElements; i++)
     		NumNz[i]=i+1;
#endif /* RSBEP_NO_STUB */
	bool btransA = RSBEP_GETENVB("RSBP_TRANS", false);
	rsb_trans_t transA = btransA ? RSB_TRANSPOSITION_T : RSB_TRANSPOSITION_N;

#if RSBEP_NO_STUB
	Epetra_CrsMatrix Crs_T(Copy, Map, &NumNz[0]);
	Epetra_CrsMatrix Crs_A(Copy, Map, &NumNz[0]);
	//MyEpetra_CrsMatrix MyCrs_A(Copy, Map, &NumNz[0]);
#else /* RSBEP_NO_STUB */
	Epetra_RsbMatrix Rsb_T = Epetra_RsbMatrix(n,n,Epetra_RsbMatrix::IsTri);
	Epetra_RsbMatrix Rsb_A = Epetra_RsbMatrix(n,n,Epetra_RsbMatrix::IsGen);
	//Epetra_RsbMatrix = Epetra_RsbMatrix(n,n,Epetra_RsbMatrix::IsSym);
#endif /* RSBEP_NO_STUB */
	Epetra_RsbMultiVector x(n,  1.1);
	Epetra_RsbMultiVector y(n,0*2.2);
	Epetra_RsbMultiVector z(n,0*2.2);
	Epetra_RsbMultiVector B(n,  1.1,nrhs);
	Epetra_RsbMultiVector C(n,0*2.2,nrhs);

	for ( RSBEP_IDXT i = 0; i <  n; ++i )
	for ( RSBEP_IDXT j = 0; j <= i; ++j )
	{
		RSBEP_NUMT val;
		val = (1+i)*1.0+(1+j)*0.1;
#if RSBEP_NO_STUB
		RSBEP_EERRREACT( Crs_A.InsertGlobalValues( i ,1,&val,&j) );
		RSBEP_EERRREACT( Crs_T.InsertGlobalValues( i ,1,&val,&j) );
#else /* RSBEP_NO_STUB */
		RSBEP_EERRREACT( Rsb_A.InsertGlobalValues( i ,1,&val,&j) );
		RSBEP_EERRREACT( Rsb_T.InsertGlobalValues( i ,1,&val,&j) );
#endif /* RSBEP_NO_STUB */
	}
#if RSBEP_NO_STUB
	RSBEP_EERRREACT( Crs_A.FillComplete() );
	RSBEP_EERRREACT( Crs_T.FillComplete() );

	RSBEP_EERRREACT( Crs_A.OptimizeStorage());
	RSBEP_EERRREACT( Crs_T.OptimizeStorage());

	Epetra_RsbMatrix Rsb_A = Epetra_RsbMatrix(Crs_A,Rsb_Matrix::IsGen);
	Epetra_RsbMatrix Rsb_T = Epetra_RsbMatrix(Crs_T,Rsb_Matrix::IsTri);
#else /* RSBEP_NO_STUB */
	RSBEP_EERRREACT( Rsb_A.FillComplete() );
	RSBEP_EERRREACT( Rsb_T.FillComplete() );
#endif /* RSBEP_NO_STUB */
#if RSBEP_NO_STUB
	Epetra_Vector ex(View, Crs_A.Map(), &x[0]); /* internally x is reused */
	Epetra_Vector ey(View, Crs_A.Map(), &y[0]); /* internally y is reused */
	Epetra_Vector ez(View, Crs_A.Map(), &z[0]); /* internally z is reused */
	Epetra_MultiVector eB(Crs_A.Map(), nrhs);
	Epetra_MultiVector eC(Crs_A.Map(), nrhs);
	Epetra_MultiVector eZ(Crs_A.Map(), nrhs);
#endif /* RSBEP_NO_STUB */
	os << "RSB constructed\n";

	os << "A:" << Rsb_A;
	os << "T:" << Rsb_T;

	os << "SPMM:\n";
	os << "B:\n" << B;
	os << "C:\n" << C;
	RSBEP_ERRREACT(Rsb_A.spmm<rsb_err_t>(transA,alpha,nrhs,order,&B[0],beta,&C[0]));
	os << "before tuning for SPMV:\n";
	os << Rsb_A._info() << std::endl;
	if (tune)
	{
		RSBEP_EERRREACT(Rsb_A.OptimizeStorage());
		RSBEP_EERRREACT(Rsb_T.OptimizeStorage());
		os << "after tuning for SPMV:\n";
		os << Rsb_A._info() << std::endl;
	}

	os << "**\n";
	x=xiv;
	y=yiv;
	os << "x:\n" << x;
	os << "y:\n" << y;
	//RSBEP_ERRREACT(Rsb_A.spmv<rsb_err_t>( &y[0],&x[0]));
#if RSBEP_NO_STUB
	//x.assign_to(ex);
	//y.assign_to(ey);
	os << "SPMV:\n";
	Crs_T.Multiply(btransA,ex,ey);
	//y=ey;
	os << "y <- Crs_T    * x:\n" << y;
	ez=ey;
	os << "SPSV:\n";
	Crs_T.Solve(/*Upper*/Crs_T.UpperTriangular(),btransA,/*UnitDiagonal*/false/*Crs_T.NoDiagonal()*/,ez,ey);
	//y=ey;
	os << "y <- Crs_T    \\ y:\n" << y;

	if( x.diverging(y) )
		os << "Crs_T MV/SV test failed! " << x.diff(y) << " \n";
#endif /* RSBEP_NO_STUB */

	os << "SPMV:\n";
	x=xiv;
	y=yiv;
	RSBEP_ERRREACT(Rsb_A.spmv<rsb_err_t>(transA,1.0,&x[0],1,0.0,&y[0],1));
	os << "y <- Rsb_A    * x:\n" << y;
	os << "SPSV:\n";
	RSBEP_ERRREACT(Rsb_T.spsv<rsb_err_t>(&y[0],&y[0],btransA));
	os << "y <- Rsb_A    \\ y:\n" << y;
	if( x.diverging(y) )
		os << "Rsb_T MV/SV test failed! " << x.diff(y) << " \n";

	os << "**\n";
	B=xiv;
	C=yiv;
	os << "x:\n" << B;
	os << "y:\n" << C;
#if RSBEP_NO_STUB
	B.assign_to(eB);
	C.assign_to(eC);
	os << "SPMM:\n";
	Crs_T.Multiply(btransA,eB,eC);
	C=eC;
	os << "y <- Crs_T    * x:\n" << C;
	eZ=eC;
	os << "SPSM:\n";
	Crs_T.Solve(/*Upper*/Crs_T.UpperTriangular(),btransA,/*UnitDiagonal*/false/*Crs_T.NoDiagonal()*/,eZ,eC);
	C=eC;
	os << "y <- Crs_T    \\ y:\n" << C;
	if( B.diverging(C) )
		os << "Crs_T MM/SM test failed! " << B.diff(C) << " \n";
#endif /* RSBEP_NO_STUB */

	B=xiv;
	C=yiv;
	os << "SPMM:\n";
	RSBEP_ERRREACT(Rsb_A.spmm<rsb_err_t>(transA,alpha,nrhs,order,&B[0],beta,&C[0]));
	os << "y <- Rsb_A    * x:\n" << C;
	os << "SPSM:\n";
	RSBEP_ERRREACT(Rsb_T.spsm<rsb_err_t>(&C[0],&C[0],nrhs,btransA));
	os << "y <- Rsb_A    \\ y:\n" << C;
	if( B.diverging(C) )
		os << "Rsb_T MM/SM test failed! " << B.diff(C) << " \n";

	Epetra_RsbMultiVector Values(n,0.0);
	std::vector<RSBEP_IDXT> Indices(n,0);

	for (RSBEP_IDXT i = 0 ; i < n ; ++i)
		Indices[i] = i,
		Values[i] = 10*(n)+(i+1);

	os << "SumIntoMyValues:\n";
	os << "A:" << Rsb_A;

	RSBEP_EERRREACT(Rsb_A.SumIntoMyValues(MyRow, NumEntries, &Values[0], &Indices[0]));
	os << "Matrix after SumIntoMyValues:\n" << Values << "\n";
	os << "A:" << Rsb_A;

	os << "ReplaceMyValues:\n";
	for ( RSBEP_IDXT i = 0; i < NumEntries; ++i )
		Values[i]*=100.0;
	os << "Adjusted Values:\n" << Values;
	RSBEP_EERRREACT(Rsb_A.ReplaceMyValues(MyRow, NumEntries, &Values[0], &Indices[0]));

	os << "A:" << Rsb_A;

#if RSBEP_NO_STUB
	Epetra_Vector Diagonal(Crs_A.RowMap());
	Epetra_Vector Result(Crs_A.RowMap());
#else /* RSBEP_NO_STUB */
	Epetra_RsbMultiVector Diagonal(n,0.0);
	Epetra_RsbMultiVector Result(n,0.0);
#endif /* RSBEP_NO_STUB */

	RSBEP_EERRREACT(Rsb_A.ExtractDiagonalCopy(Diagonal));
	os << "Diagonal Values Before:\n" << Diagonal;

	for ( RSBEP_IDXT i = 0; i < n; ++i )
		Diagonal[i]*=-1.0;
	RSBEP_EERRREACT(Rsb_A.ReplaceDiagonalValues(Diagonal));
	for ( RSBEP_IDXT i = 0; i < n; ++i )
		Diagonal[i]=0.0;
	RSBEP_EERRREACT(Rsb_A.ExtractDiagonalCopy(Diagonal));
	os << "Diagonal Values After :\n" << Diagonal;

#if 0
{
	RSBEP_IDXT MyRow = n-1, Length = n, NumEntries = 0;
	std::vector<RSBEP_NUMT> Values(Length,0.0);
	std::vector<RSBEP_IDXT> Indices(Length,0);

	os << "Values Before ExtractMyRowCopy (=uninitialized):\n" << Values<< "\n";
	Rsb_A.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0], &Indices[0]);
	os << "Values After ExtractMyRowCopy :\n" << Values;

	Rsb_A.PutScalar (55.0);
	Rsb_A.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0], &Indices[0]);
	os << "Values After PutScalar:\n" << Values;

	Rsb_A.Scale (10.0);
	Rsb_A.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0], &Indices[0]);
	os << "Values After ExtractMyRowCopy:\n" << Values;
}
#endif

#if RSBEP_NO_STUB
{
	const RSBEP_IDXT nnz = ((n+1)*n)/2;
	RSBEP_IDXT Length = n;
	std::vector<RSBEP_NUMT> VA(nnz,0.0);
	std::vector<RSBEP_IDXT> JA(nnz,0);
	std::vector<RSBEP_IDXT> RP(n+1,0);
   	std::vector<RSBEP_IDXT> NumNz(NumMyElements,0);

	RP[0]=0;
	for ( RSBEP_IDXT i = 1; i <= n; ++i )
		NumNz[i-1]=i,
		RP[i]=RP[i-1]+NumNz[i-1];
		// os << "numnz/rp   " << NumNz[i-1] << " / " << RP[i] << "\n";

	for ( RSBEP_IDXT i = 0; i < n; ++i )
	for ( RSBEP_IDXT p = RP[i]; p < RP[i+1] ; ++p )
	{
		RSBEP_IDXT j = 0 + p-RP[i];

		JA[p]=j;
		VA[p]=(1+i)*100000+j*100;
		//os << JA[p] << " " << VA[p];
	}

#if RSBEP_NO_STUB
	Epetra_CrsMatrix ecm(Copy, Map, &NumNz[0]);
#if 0
	// element per element
	for ( RSBEP_IDXT i = 0; i < n; ++i )
	for ( RSBEP_IDXT p = RP[i]; p < RP[i+1] ; ++p )
		RSBEP_EERRREACT( ecm.InsertGlobalValues( i ,1,&VA[p],&JA[p]) );
#endif
#if 1
	// lower triangle, element per elemnt
	for ( RSBEP_IDXT i = 0; i <  n; ++i )
	for ( RSBEP_IDXT j = 0; j <= i; ++j )
	{
		RSBEP_NUMT val = (1+i)*1.0+(1+j)*0.1;
		RSBEP_EERRREACT( ecm.InsertGlobalValues( i ,1,&val,&j) );
	}
#endif
#if 0
	// row per row
	for ( RSBEP_IDXT i = 0; i < n; ++i )
	{
		// RSBEP_EERRREACT( ecm.InsertGlobalValues( i ,1,&VA[RP[i]],&JA[RP[i]]));
		// RSBEP_EERRREACT( ecm.InsertGlobalValues( i ,i+1,&VA[0],&JA[0]) );
		RSBEP_EERRREACT( ecm.InsertGlobalValues( i ,NumNz[i],&VA[RP[i]],&JA[RP[i]]) );
	}
#endif
	RSBEP_EERRREACT( ecm.FillComplete() );
	ecm.OptimizeStorage();
#else /* RSBEP_NO_STUB */
	Epetra_CrsMatrix ecm (n, RP, JA, VA); /* FIXME: non existent interface ! */
#endif /* RSBEP_NO_STUB */
	Epetra_RsbMatrix R(ecm);
	os << "R:\n" << R;
	os << "Diagonal Before Multiply1:\n" << Diagonal;
	RSBEP_EERRREACT(Rsb_T.Multiply1(btransA, Diagonal, Result));
	os << "Diagonal After  Multiply1:\n" << Diagonal;
	RSBEP_EERRREACT(Rsb_T.Solve1(isUpper,btransA,unitDiagonal, Diagonal, Diagonal));
	os << "Diagonal After  Solve1:   \n" << Diagonal;
}
#endif /* RSBEP_NO_STUB */
	//Rsb_A.spsm(&y[0],&x[0],nrhs);
	// os << y << "\n";
  	// RSBEP_ERRT Solve1(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  	// RSBEP_ERRT Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  	// RSBEP_ERRT Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
	perrval = EXIT_SUCCESS;
	os << " terminating run with RSBEP_NO_STUB=" << RSBEP_NO_STUB << ", exit code=" << perrval << "\n";
	return perrval;
} /* Epetra_RsbMatrix_demo */

#ifdef HAVE_GETOPT_H
struct rsbp_option { struct option o;const char*dv;const char*od;};

void rsbp_help(char * argv[], const struct option*options, const struct rsbp_option*ioptions, int nopts)
{
	int c;
	RSBP_PRINTF("Usage: {mpi-launcher} %s OPTIONS\n",argv[0]);
	for(c=0;c<nopts;++c)
	{
		int ia=isprint(options[c].val);
		RSBP_PRINTF(" [ { --%s%s%c}%s%s%s ]%s%s\n",
			options[c].name,
			ia?" | -":" ",
			ia?options[c].val:' ',
			options[c].has_arg==required_argument?(ioptions[c].dv?" {val=":" {val"):
			(options[c].has_arg==optional_argument?"[{val=":""),
			options[c].has_arg==required_argument?(ioptions[c].dv?ioptions[c].dv:""):
			(((options[c].has_arg==optional_argument)&&ioptions[c].dv)?ioptions[c].dv:""),
			options[c].has_arg==required_argument?"}":
			(options[c].has_arg==optional_argument?"}]":""),
			(ioptions[c].od?"\t#":""),
			(ioptions[c].od?ioptions[c].od:"")
			);
	}
} /* rsbp_help */
#endif /* HAVE_GETOPT_H */

int main(int argc, char *argv[])
{
	int perrval = EXIT_SUCCESS;
	bool wantrsbtest = false;
	bool wanterctest = false;
	bool wantrfbtest = false;

	const rsb_char_t * filename = RSBEP_GETENVS("RSBP_FILENAME",RSBP_NULL);
#ifdef HAVE_GETOPT_H
	const struct rsbp_option ioptions[] = {
    		{"epetra",no_argument,RSBP_NULL,'e',RSBP_NULL,RSBP_NULL},
    		{"rsb"   ,no_argument,RSBP_NULL,'r',RSBP_NULL,RSBP_NULL},
    		{"file-based"   ,no_argument,RSBP_NULL,'f',RSBP_NULL,RSBP_NULL},
    		{"help"  ,no_argument,RSBP_NULL,'h',RSBP_NULL,RSBP_NULL},
	};
	const int ion=sizeof(ioptions)/sizeof(struct rsbp_option);
	struct option options[ion];

	for(int ioi=0;ioi<ion;++ioi)options[ioi]=ioptions[ioi].o;/* translation */

	for (;;)
	{
		int c,opt_index;

		c = getopt_long(argc, argv, "efhr", options, &opt_index);
		if (c == -1)
		break;
		switch (c)
		{
			case 'e':
				wanterctest = true;
			break;
			case 'f':
				wantrfbtest = true;
			break;
			case 'r':
				wantrsbtest = true;
			break;
			case 'h':
			default:
				rsbp_help(argv,options,ioptions,sizeof(options)/sizeof(options[0]));
				goto ret;
		}
	}

	for (int i = optind; i < argc; i++)
		filename = argv[i];
#endif /* HAVE_GETOPT_H */

	if( ! ( wantrsbtest || wanterctest || wantrfbtest ) )
		wantrsbtest =  wanterctest =  wantrfbtest = true;

	rsb_lib.meminfo();
#if !RSBP_WANT_INLINE_INIT
	if ( RSBEP_GETENVI("RSBP_WANT_VERBOSE_TUNING", 0) )
		rsb_lib_set_opt_str("RSB_IO_WANT_VERBOSE_TUNING","1");
#endif /* RSBP_WANT_INLINE_INIT */

	if( wanterctest )
	{
		rsb_lib.meminfo();
#if defined(RSBEP_NO_STUB) && defined(EPETRA_MPI)
		perrval |= Epetra_RsbMatrix_demo(argc,argv);
#else
		perrval |= Epetra_RsbMatrix_demo();
#endif
		if( perrval != EXIT_SUCCESS)
			RSBEP_PRERR_GOTO(pret,"")
	}
	
	if( wantrsbtest )
	{
		rsb_lib.meminfo();
		perrval |= Rsb_Matrix_demo();
		if( perrval != EXIT_SUCCESS)
			RSBEP_PRERR_GOTO(pret,"")
	}

	if( wantrfbtest && filename )
	{
		rsb_lib.meminfo();
		perrval |= Rsb_Matrix_file_demo(filename);
		if( perrval != EXIT_SUCCESS)
			RSBEP_PRERR_GOTO(pret,"")
	}
pret:
	if( 0 == rsb_lib.meminfo() )
	{
		if( perrval == EXIT_SUCCESS )
			RSBEP_OS << "OK: terminating with no allocations registered in librsb\n";
	}
#ifdef HAVE_GETOPT_H
ret:
#endif /* HAVE_GETOPT_H */
	return perrval;
} /* main */

/* rsb.cpp */
