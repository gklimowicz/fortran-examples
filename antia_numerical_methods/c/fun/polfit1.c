/*	Least squares polynomial fit using orthogonal polynomials in 1 dimension
	Modified version of POLFIT to fit multiple sets of function values
	This function is meant to be used for fit in multiple dimensions

	N : (input) Number of data points
	M : (input) Required degree of polynomial
	NUM : (input) Number of different RHS (function values) to be fitted
		Each set must be defined over the same abscissas and
		with same weights.
	X : (input) Array of length N containing the abscissas
	F : (input) Array of length N*NUM containing the function values
		F[J][I] is the value at X[I] in Jth data set
		The second dimension of F is assumed to be exactly equal
		to N to minimise storage requirement.
	SIG : (input) Array of length N containing the errors associated
		with each point. Errors are same for all data sets.
		SIG[I] is the error for F[J][I].
	A : (output) Array of length (M+1)*NUM containing the coefficients
		for the fit for each RHS. The second dimension of A
		is assumed to be M+1 to minimise storage requirements
		A[j][i] is the ith coefficient for jth RHS.
	ALP, BETA : (output) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
	GAM : (output) Array of length M+1, containing the quantities
		\gamma_i for the orthogonal polynomials
		
	Error status is returned by the value of the function POLFIT1.
		0 value implies successful execution
		601 implies that N<M+1 or M<0 or N < 1
		621 implies that GAM[I] vanishes at some I and calculations
			are abandoned
	
	The fitted polynomial can be calculated at any value of x using POLEVL

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int polfit1(int n, int m, int num, double x[], double *f, double sig[],
	double *a, double alp[], double beta[], double gam[])

{
	int i,j,i1,i2,it,j1,m1;
	double gam0,s;
	double *wk;

	if(m>=n || m<0 || n<1) return 601;

/*	Initialisation */
	i1=0; i2=1; m1=m+1;
	gam0=0.0;
	wk=(double *) calloc((size_t) (2*n), sizeof(double));
	for(i=0; i<n; ++i) {
		wk[i]=0.0; wk[i+n]=1.0;
		gam0=gam0+1/(sig[i]*sig[i]);
	}
	gam[0]=gam0;
	beta[0]=0.0;

/*	Loop over the degree of polynomial */
	for(j=0; j<=m; ++j) {
		if(gam[j] <=0.0) {free(wk); return 621;}

		for(j1=0; j1<num; ++j1) {
			s=0.0;
			for(i=0; i<n; ++i) s=s+f[i+j1*n]*wk[i+i2*n]/(sig[i]*sig[i]);

/*	The coefficient a_j */
			a[j+j1*m1]=s/gam[j];
		}
		if(j == m) {free(wk); return 0;}

		s=0.0;
		for(i=0; i<n; ++i) s=s+x[i]*wk[i+i2*n]*wk[i+i2*n]/(sig[i]*sig[i]);
/*	The coefficient \alpha_{j+1} */
		alp[j]=s/gam[j];
		gam0=0.0;
		for(i=0; i<n; ++i) {
			wk[i+i1*n]=(x[i]-alp[j])*wk[i+i2*n]-beta[j]*wk[i+i1*n];
			gam0=gam0+wk[i+i1*n]*wk[i+i1*n]/(sig[i]*sig[i]);
		}
/*	The coefficient \beta_{j+1} */
		beta[j+1]=gam0/gam[j];
/*	The coefficient \gamma_{j+1} */
		gam[j+1]=gam0;
/*	Interchange indices I1, I2 so that only last two columns of WK are stored */
		it=i1; i1=i2; i2=it;
	}
	free(wk);
	return 0;
}
