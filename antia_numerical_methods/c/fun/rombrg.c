/*	To integrate a function over finite interval using Romberg integration

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	GI : (input/output) Array of length NMAX (=13), containing
		the expected values of exponents \gamma_i in error expansion
		If GI[I]<=0 it will be set to 2I+2, the correct value for
		a smooth function
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) Estimated (absolute) error achieved by the function
	N : (input/output) On input it should contain the number of function evaluations
		to be used for first estimate. If N<2 or N>NPT it is set to 2.
		After execution it will contain the number of function
		evaluations actually used by ROMBRG.
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function ROMBRG.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy
		33 implies that N>NPT (=100) in which case it is set to 2

	Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/		

#include <math.h>

int rombrg(double *ri, double a, double b, double gi[], double reps,
	double aeps, double *dif, int *n, double (*fun)(double ))

{
	int i,j,nd,ier, nmin=3, nmax=13, npt=100;
	double dif1,fj,r1,s1,h, t[13][13];

	for(i=0; i<nmax; ++i) {
		if(gi[i] <= 0.0) gi[i]=2*i+2;
	}

	ier=0;
	if((*n)<=1) *n=2;
	if((*n)>npt) {ier=33; *n=2;}

/*	Contribution from the end points */
	s1=0.5*(fun(a)+fun(b));
/*	First time use all points */
	nd=1; *dif=0.0;

	for(i=0; i<nmax; ++i) {
		h=(b-a)/(*n-1);

/*	Add new points to the sum */
		for(j=1; j<=(*n)-2; j=j+nd) s1=s1+fun(a+j*h);
/*	The trapezoidal rule approximation */
		t[i][0]=s1*h;

/*	The Richardson's extrapolation */
		for(j=0; j<i-1; ++j) {
			fj=pow(2.,gi[j]);
			t[i][j+1]=t[i][j]+(t[i][j]-t[i-1][j])/(fj-1.0);
			dif1=fabs(t[i][j]-t[i-1][j]);
/*	Find the minimum difference between the last two rows of T-table */
			if(dif1 < (*dif) || j == 0) {*dif=dif1; *ri=t[i][j];}
		}

/*	On second and subsequent pass add only new points to the sum */
		nd=2;
		if(i>=nmin) {
			r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
			if(*dif < r1) return ier;
		}
		*n=2*(*n)-1;
	}

	*n=((*n)+1)/2;
	return 30;
}

