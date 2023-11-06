/*	To integrate a function over semi-infinite interval using a
	combination of Gauss-Legendre and Gauss-Laguerre formulas

	RINT : (output) Calculated value of the integral
	A : (input) The lower limit
	A1 : (input/output) The point at which integral has to be broken
		A1 will be adjusted by the function.
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	F2 : (input) Name of the function to calculate F(X)*EXP(X)
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function GAULAG.
		0 value implies successful execution
		30 implies specified accuracy was not achieved by GAUSS
		37 implies specified accuracy was not achieved by LAGURE
		38 implies specified accuracy was not achieved by both
			GAUSS and LAGURE
		In all cases DIF will contain the estimated accuracy

	Functions F(X) and F2(X) must be supplied by the user.

	Required functions : GAUSS, LAGURE, F, F2
*/

#include <math.h>

int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ));
int lagure(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt);

int gaulag(double *rint, double a, double *a1, double reps, double aeps,
	double *dif, double (*f) (double ), double (*f2) (double ), int *np)

{
	int npt,npt1,ier,ier1;
	double r1,r2,t1,dif1,a2, amax=50.0;

	r1=0.0; r2=0.0; *dif=0.0;
	*np=0; npt=0; ier1=0;
	if( (*a1)< a) *a1=a;

/*	To calculate integral over [A1,Infinity) */
	do {
		ier=lagure(&r1,*a1,aeps,reps,&dif1,f2,&npt1);
		*np=(*np)+npt1;
		if(ier == 0) break;

/*	If LAGURE fails then increase A1 */
		t1= (*a1);
		a2=(*a1)*2.; if(*a1+2. > a2) a2=(*a1)+2.;
		*a1=a2;
	} while (*a1 < amax);

	if(ier >0) {ier1=37; *a1=t1;}

/*	To calculate integral over [A,A1] */
	npt1=16;
	ier=0;
	if((*a1)-a > aeps) ier=gauss(&r2,a,*a1,&npt1,reps,aeps,dif,&npt,f);
	*rint=r1+r2;
	*dif=fabs(*dif)+fabs(dif1);
	*np=(*np)+npt;
	ier=ier+ier1;
	if(ier>ier1 && ier1>0) ier=38;
	return ier;
}
