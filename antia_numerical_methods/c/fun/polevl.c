/*	Evaluating the fitted polynomial and its derivatives at any value
	of x using known coefficients of orthogonal polynomials
	Should be used to evaluate the polynomial using coefficients calculated
	by POLFIT.

	M : (input) Degree of polynomial
	A : (input) Array of length M+1 containing the coefficients
		of the fit
	ALP, BETA : (input) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
		A, ALP, BETA could be calculated using POLFIT
	X : (input) Value of x at which polynomial needs to be evaluated
	F : (output) Calculated value of the polynomial at X
	DF : (output) First derivative of F(X) at X
	DDF : (output) Second derivative of F(X) at X

	The returned value of POLEVL is always zero
	
	Required functions : None
*/

#include <math.h>

int polevl(int m, double a[], double alp[], double beta[], double x,
	double *f, double *df, double *ddf)

{
	int j;
	double f1,df1,ddf1,ff,d,dd;

	*f=a[m-1]+(x-alp[m-1])*a[m];
	f1=a[m];
	*df=a[m]; df1=0.0;
	*ddf=0.0; ddf1=0.0;
	if(m<=1) return 0;

/*	Clenshaw's recurrence for F, DF and DDF */
	for(j=m-2; j>=0; --j) {
		dd=2.*(*df)+(x-alp[j])*(*ddf)-beta[j+1]*ddf1;
		d=(*f)+(x-alp[j])*(*df)-beta[j+1]*df1;
		ff=a[j]+(x-alp[j])*(*f)-beta[j+1]*f1;
		f1=(*f); *f=ff;
		df1=(*df); (*df)=d;
		ddf1=(*ddf); *ddf=dd;
	}
	return 0;
}
