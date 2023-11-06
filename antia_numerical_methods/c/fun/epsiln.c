/*	To integrate a function over finite interval using Epsilon algorithm

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) estimated (absolute) error achieved by the function
	N : (input/output) On input it should contain the number of function
		evaluations to be used for first estimate. If N<2 or N>NPT it
		is set to 2. After execution it will contain the number of
		function evaluations actually used by the function
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function EPSILN.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy
		33 implies that N>NPT (=100) in which case it is set to 2
		34 implies that at some stage denominator vanished while
			calculating epsilon table. This value is ignored.
		35 implies that roundoff error appears to be dominating
		
	Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>

int epsiln(double *ri, double a, double b, double reps, double aeps,
	double *dif, int *n, double (*fun)(double ))

{
	int i,j,nd,ier, nmin=3, nmax=13, npt=100;
	double h,ri1,dif0,dif1,s1,den,r1,t[13][13];

	ier=0;
	if((*n)<=1) *n=2;
	if((*n)>npt) {ier=33; *n=2;}

/*	Sum of the end points for trapezoidal rule */
	s1=0.5*(fun(a)+fun(b));
	nd=1; *ri=0.0; t[0][0]=0.0;

	for(i=0; i<nmax-1; ++i) {
		h=(b-a)/(*n-1);

		for(j=1; j<=(*n)-2; j=j+nd) s1=s1+fun(a+j*h);
/*	The trapezoidal rule approximation */
		t[i][1]=s1*h;
		t[i+1][0]=0.0;
		ri1=*ri;
		if(i>0) {*dif=fabs(t[i][1]-t[i-1][1]); *ri=t[i][1];}

/*	Construct the Epsilon table */
		for(j=2; j<=i+1; ++j) {
			den=t[i-j+2][j-1]-t[i-j+1][j-1];

/*	If denominator is zero set the error flag */
			if(den != 0.0) t[i-j+1][j]=t[i-j+2][j-2]+1./den;
			else {ier=34; t[i-j+1][j]=t[i-j+2][j-2];}
		}

		if(i>3) {

/*	DIF is the minimum difference between two rows of epsilon table */
			for(j=3; j<=i-1; j=j+2) {
				dif1=fabs(t[i-j+1][j]-t[i-j][j]);
				if(dif1 < (*dif)) {*dif=dif1; *ri=t[i-j+1][j];}
			}
		}

		nd=2;
		if(i>5 && *dif>dif0) {
/*	Roundoff error appears to be dominating, retain the previous value of RI */
			*ri=ri1;
			return 35;
		}
		if(i>nmin) {
			dif0=*dif;
			r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
			if(*dif < r1) return ier;
		}
		*n=2*(*n)-1;
	}

	*n=(*n+1)/2;
	return 30;
}
