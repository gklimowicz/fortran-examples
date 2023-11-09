/*	Least squares polynomial fit using orthogonal polynomials in 1 dimension

	N : (input) Number of data points
	M : (input) Required degree of polynomial
	X : (input) Array of length N containing the abscissas
	F : (input) Array of length N containing the function values
		F[I] is the value at X[I]
        SIG : (input) Real array of length N containing the estimated
	        error in the function values, SIG(I) is the error in F(I)
	A : (output) Array of length M+1 containing the coefficients for the fit
	ALP, BETA : (output) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
	Y : (output) Array of length N containing the fitted values at X[I]
	H : (output) Array of length M+1 containing the \chi^2 values
		for residuals using polynomial fit of degrees 0,...,M
	GAM : (output) Array of length M+1, containing the quantities
		\gamma_i for the orthogonal polynomials
		
	Error status is returned by the value of the function POLFIT.
		0 value implies successful execution
		601 implies that N<M+1 or M<0 or N < 1
		621 implies that GAM[I] vanishes at some I and calculations
			are abandoned
	
	The fitted polynomial can be calculated at any value of x using POLEVL

	Required functions : None

	THE ARGUMENTS OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
	VERSION AS THE 5TH ARGUMENT IS NOW ERROR INSTEAD OF WEIGHT

*/

#include <math.h>
#include <stdlib.h>

int polfit(int n, int m, double x[], double f[], double sig[], double a[],
	double alp[], double beta[], double y[], double h[], double gam[])

{
	int i,j,i1,i2,it;
	double gam0,h0,s;
	double *wk;

	if(m>=n || m<0 || n<1) return 601;

/*	Initialisation */
	i1=0; i2=1;
	gam0=0.0; h0=0.0;
	wk=(double *) calloc((size_t) (2*n), sizeof(double));
	for(i=0; i<n; ++i) {
		wk[i]=0.0; wk[i+n]=1.0;
		y[i]=0.0;
		gam0=gam0+1/(sig[i]*sig[i]);
		h0=h0+f[i]*f[i]/(sig[i]*sig[i]);
	}
	gam[0]=gam0;
	beta[0]=0.0;

/*	Loop over the degree of polynomial */
	for(j=0; j<=m; ++j) {
		s=0.0;
		for(i=0; i<n; ++i) s=s+f[i]*wk[i+i2*n]/(sig[i]*sig[i]);
		if(gam[j] <=0.0) {free(wk); return 621;}

/*	The coefficient a_j */
		a[j]=s/gam[j];
		h0=h0-a[j]*a[j]*gam[j];
		h[j]=h0;
		for(i=0; i<n; ++i) y[i]=y[i]+a[j]*wk[i+i2*n];
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
