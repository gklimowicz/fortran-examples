/*	Evaluating the orthogonal polynomial basis functions at any value
	of x using known coefficients
	Should be used to evaluate the basis using coefficients calculated
	by POLFIT or POLFIT1.

	M : (input) Degree of polynomial
	ALP, BETA : (input) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
		ALP, BETA could be calculated using POLFIT
	X : (input) Value of x at which polynomials needs to be evaluated
	F : (output) Array of length M+1 containing the value of
		orthogonal polynomials at X
	DF : (output) Array of length M+1 containing first derivative of F at X
	DDF : (output) Array of length M+1 containing second derivative of F at X

	Returned value POLORT is always zero.
	
	Required functions : None
*/	

#include <math.h>

int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[])

{
	int j;
	
	f[0]=1.0; f[1]=x-alp[0];
	df[0]=0.0; df[1]=1.0;
	ddf[0]=0.0; ddf[1]=0.0;

 
/*	The recurrence relations */
	for(j=2; j<=m; ++j) {
		ddf[j]=2.*df[j-1]+(x-alp[j-1])*ddf[j-1]-beta[j-1]*ddf[j-2];
		df[j]=f[j-1]+(x-alp[j-1])*df[j-1]-beta[j-1]*df[j-2];
		f[j]=(x-alp[j-1])*f[j-1]-beta[j-1]*f[j-2];
	}
	return 0;
}
