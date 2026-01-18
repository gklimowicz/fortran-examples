/*	To integrate a function over finite interval using Simpson's rule

	RI : (output) Calculated value of the integral
	XL : (input) The lower limit
	XU : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) estimated (absolute) error achieved by the function
	N : (output) Number of function evaluations used by SIMSON
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function SIMSON.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy

	Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>


int simson(double *ri, double xl, double xu, double reps, double aeps,
	double *dif, int *n, double (*fun)(double ))

{
	int i,j,n2, nmin=3, nmax=13;
	double r1,r2,fend,odd,even,x,x1,h2,h;

	fend=fun(xl)+fun(xu);
	even=0.0; odd=0.0;
	*ri=0.0; *dif=0.0;

	*n=2;
/*	starting with 2+1 points, subdivide the intervals into 2 until convergence */
	h=(xu-xl);
	if(h==0) return 0;

	for(i=1; i<= nmax; ++i) {
		h=h/2.;
		even=even+odd;
		odd=0.0;
		x1=xl+h;
		n2=(*n)/2;
		h2=2.*h;

		for(j=0; j<n2; ++j) {
			x=x1+h2*j;
			odd=odd+fun(x);
		}

/*	Estimate for the integral */
		r1=(fend+4.*odd+2.*even)*h/3.;
		*dif=r1-*ri;
		*ri=r1;
/*	To avoid spurious convergence in first few trials skip the convergence test */
		if(i>nmin) {
			r2=reps*fabs(r1); if(aeps>r2) r2=aeps;
			if(fabs(*dif) < r2) return 0;
		}
		*n = (*n)*2;
	}

	*n = (*n)/2;
	return 30;
}
