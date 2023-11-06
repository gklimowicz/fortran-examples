/*	Multiple integration over a hyper-rectangle in n-dimensions
	using Monte-Carlo technique

	A : (input) Array of length N containing the lower limit along each dimension
	B : (input) Array of length N containing the upper limit along each dimension
	N : (input) The number of dimensions
	NPT : (input) Maximum number of function evaluations to be used
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RI : (output) The calculated value of the integral
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	ERR : (output) estimated (absolute) error achieved by the function
		It is 2.576 times the estimated standard deviation
		divided by SQRT(NP).
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function MCARLO.
		0 value implies successful execution
		39 implies specified accuracy was not achieved in
			which case ERR will contain the estimated accuracy
		311 implies N<1, and no calculations are done

	Function F(N,X) should be supplied by the user
	
	Required functions : RANF, F
*/

#include <math.h>
#include <stdlib.h>

double ranf(int *iseed);
double ran(double *seed);


int mcarlo(double a[], double b[], int n, int npt, double (*f) (int , double * ),
	double *ri, double reps, double aeps, double *err, int *np)

{
	int i,j,iseed,npt1, nchk=100;
	double ri1,var,var1,hh,f1,r1;
	double *h, *xa;

	if(n<1) return 311;
	*ri=0.0; *np=0;
	h=(double *) calloc((size_t) n, sizeof(double));
	xa=(double *) calloc((size_t) n, sizeof(double));

	hh=1.0;
	for(i=0; i<n; ++i) {h[i]=b[i]-a[i]; hh=hh*h[i];}

	ri1=0.0; var1=0.0;

/*	Seed for random number generator, should be changed if another routine is used */
	iseed=-12345;
	npt1=nchk;

	for(i=1; i<=npt; ++i) {
/*	Generating the abscissas */
		for(j=0; j<n; ++j) xa[j]=a[j]+h[j]*ranf(&iseed);
		f1=f(n,xa);
		ri1=ri1+f1;
		var1=var1+f1*f1;

		if(i == npt1) {
/*	Compute intermediate sums to check for convergence */
			*ri=ri1/i;
			var=var1/i-(*ri)*(*ri);
			if(var<0.0) var=0.0;
			*err=2.576*hh*sqrt(var/npt1);
			*ri=(*ri)*hh;
			*np=i;
			r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
			if(*err < r1) {free(xa); free(h); return 0;}
			npt1=2*npt1;
		}
	}

/*	Integral fails to converge */
	*ri=ri1/npt;
	var=var1/npt-(*ri)*(*ri);
	*err=2.576*hh*sqrt(var/npt);
	*ri=(*ri)*hh;
	*np=npt;
	r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
	if(*err > r1) {free(xa); free(h); return 39;}
	free(xa); free(h);
	return 0;
}
