/*     To integrate a function with logarithmic singularity over (0,A]
     using a combination of Gaussian formulas

	RINT : (output) Calculated value of the integral
	A : (input) The upper limit
	A1 : (input/output) The point at which integral has to be broken
		A1 will be adjusted by the function.
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	F2 : (input) Name of the function to calculate F(X)/LOG(X)
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function GAULG2.
		0 value implies successful execution
		31 implies specified accuracy was not achieved by GAUSS over [A1,A]
		32 implies specified accuracy was not achieved by GAULOG
		34 implies specified accuracy was not achieved by GAUSS over (0,A1]
		In case of multiple failures second digit will be sum
			of these values.
		In all cases DIF will contain the estimated accuracy

     	FUNCTION F(X) and F2(X) must be supplied by the user.

	Required functions : GAUSS, GAULOG, F, F2
*/

#include <math.h>

int gaulog(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt);
int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ));


int gaulg2(double *rint, double a, double *a1, double reps, double aeps,
	double *dif, double (*f) (double ), double (*f2) (double ), int *np)

{
	int npt,npt1,npt2,ier,ier1,ier2;
	double r1,r2,r3,a2,dif1,dif2,t1, amn=0.01;

	r1=0.0; r2=0.0; r3=0.0; a2=0.0; *dif=0.0; dif2=0.0;
	*np=0; npt=0; npt2=0; ier1=0; ier2=0;
	if( *a1 > a) *a1=a;
	if( *a1<= 0.0) *a1=a;
 
/*     Evaluate the integral over (0,A1] */
	do {
		ier=gaulog(&r1,*a1,aeps,reps,&dif1,f2,&npt1);
		*np=(*np)+npt1;
		if(ier == 0) break;

/*     If GAULOG fails decrease A1 */
		t1= (*a1);
		*a1=(*a1)/2.;
	} while (*a1 > amn);
 
	if(ier>0) {ier1=2; *a1=t1;}

/*     Evaluate the integral over [A1,A] */
	npt1=16;
	ier=0;
	if(a-(*a1) > aeps) ier=gauss(&r2,*a1,a,&npt1,reps,aeps,dif,&npt,f);
	if(ier>0) ier=1;

/*     Evaluate the regular part over [0,A1] */
	if((*a1) != 1) ier2=gauss(&r3,a2,*a1,&npt1,reps,aeps,&dif2,&npt2,f2);
	if(ier2>0) ier2=4;
	*rint=r1+r2-r3*log(*a1);
	*dif=fabs(*dif)+fabs(dif1)+fabs(dif2);
	*np = (*np)+npt+npt2;
	ier=ier+ier1+ier2;
	if(ier>0) ier=ier+30;
	return ier;
}
