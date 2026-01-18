/*	Integrate  f(x)*sqrt((x-a)/(b-x))  over (a,b) */

#include <stdio.h>
#include <math.h>

double fun(double x);
int gaucb1(double *rint, double a, double b, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double));

double b;
main()
{
	int i,i1,j,nuse, id, iflg, ier,np;
	double xl, xu, rint, reps, aeps, df, fb[20];

/*	Exercise 6.18 */

	aeps=1.e-19; reps=1.e-14;

	for(i1=0; i1<99; ++i1) {
		printf("type  xl = lower limit,    xu = upper limit \n");
		printf("                    (quits when xl = xu)\n");
		scanf(" %le %le", &xl, &xu);
		if(xl==xu) return 0;

		b=xu;
		i=gaucb1(&rint,xl,xu,reps,aeps,&df,&nuse,fun);
		printf(" ier = %d  no. of function evaluations = %d  interval = %e , %e \n", i,nuse,xl,xu);
		printf(" integral = %e   estimated error = %e  \n", rint, df);
	}
	return;
}

/*   Th integrand for caucb1 */

double fun(double x)

{
	return exp(x);
}

 
/*	To integrate a function over finite interval using Gauss-Chebyshev formulas
	Calculates the integral of FUN(X)*SQRT(X-A)/SQRT(B-X)

	RINT : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function

	FUN : (input) Name of the function to calculate the
		integrand multiplied by SQRT((B-X)/(X-A))
		
	Error status is returned by the value of the function GAUCB1.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy

		Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>

#define PI 3.14159265358979324e0

int gaucb1(double *rint, double a, double b, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double))

{
	int i,j,n, nmax=13;
	double r1,a1,a0,dx,da,r2;

	n=1;
	*rint=0.0;
	*npt=0;
	da=(b-a);

	for(i=0; i<nmax; ++i) {
		n=n*2;
		dx=PI/(4.*n+2);
 
/*	Apply N-point formula after transforming the range to (0,1) */
		r1=0.0;
		for(j=1; j<=n; ++j) {
			a1=cos((2*j-1)*dx); a1=a1*a1;
			r1=r1+a1*fun(a+da*a1);
		}
		r1=r1*dx*4.*da;

		*dif=r1-(*rint);
		*rint=r1;
		*npt=(*npt)+n;
		r2=reps*fabs(*rint); if(aeps>r2) r2=aeps;
		if(i>2 && fabs(*dif) < r2) return 0;
	}

	return 30;
}
