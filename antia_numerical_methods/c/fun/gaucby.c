/*	To integrate a function over finite interval using Gauss-Chebyshev formulas
	Calculates the integral of FUN(X)/SQRT((X-A)*(B-X))

	RINT : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	NPT : (output) Number of function evaluations used
	FUN : (input) Name of the function to calculate the
		integrand multiplied by SQRT((X-A)*(B-X))
		
	Error status is returned by the value of the function GAUCBY.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy

		Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>

#define PI 3.14159265358979324e0

int gaucby(double *rint, double a, double b, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double))

{
	int i,j,n, nmax=13;
	double r1,a1,a0,dx,da,r2;

	n=1;
	*rint=0.0;
	*npt=0;
	a0=(b+a)/2.;
	da=(b-a)/2.;

	for(i=0; i<nmax; ++i) {
		n=n*2;
		dx=PI/(2.*n);
 
/*	Apply N-point formula after transforming the range to (-1,1) */
		r1=0.0;
		for(j=1; j<=n; ++j) {
			a1=(2*j-1)*dx;
			r1=r1+fun(a0+da*cos(a1));
		}
		r1=r1*dx*2.;

		*dif=r1-(*rint);
		*rint=r1;
		*npt=(*npt)+n;
		r2=reps*fabs(*rint); if(aeps>r2) r2=aeps;
		if(i>2 && fabs(*dif) < r2) return 0;
	}

	return 30;
}
