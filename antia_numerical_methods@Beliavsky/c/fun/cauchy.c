/*	To evaluate the Cauchy principal value of an integral over finite interval

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit (B > A)
	C : (input) Location of the singularity (A < C < B)
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	FUNP : (input) Name of the function to calculate F(C+x)+F(C-x)
	NPT : (output) Number of function evaluations used
		
	Error status is returned by the value of the function CAUCHY.
		0 value implies successful execution
		304 implies A>B, A>C or C>B in which case no calculations
			are done
		Other values may be set by function ADPINT which is called
		twice. The returned value is set to IER+IER1*2, where IER
		and IER1 are the values returned by ADPINT after the two calls.
		In this case DIF will contain the estimated accuracy

	Function F(X) must be supplied by the user.
	Function FUNP(X) to calculate F(C+X)+F(C-X) should also be
		supplied by the user. The value of C is available through
		global variable CAU_CC. Simplest version for FUNP may be

	double funp(double x)

	{
		return f(CAU_CC+x)+f(CAU_CC-x);
	}

		If F(C+X)+F(C-X) can be combined roundoff error may be reduced.
		There is no provision to pass the name F to FUNP, so it
		will have to put explicitly.

	Required functions : ADPINT, KRONRD, F, FUNP
*/

#include <math.h>

/*	To pass the value of C to function F or FUNP */
double CAU_CC;

int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax);

int cauchy(double *ri, double a, double b, double c, double reps, double aeps,
	double *dif, double (*f) (double ), double (*funp) (double ), int *npt)

{
	int npt1,npt2,nmax,ier,ier1;
	double r,aa,dif1,dif2,ri1,ri2,bb;

	if(a>b || a>c || c>b) return 304; 
 
/*     FIRST EVALUATE THE SINGULAR PART */

	CAU_CC=c;
	r=c-a; if(b-c < r) r=b-c;
	aa=0.0; ri1=0.0; dif1=0.0;
	nmax=0; npt1=0;
	ier=adpint(&ri1,aa,r,reps,aeps,&dif1,funp,&npt1,&nmax);
 
/*     EVALUATE THE REMAINING PORTION */

	if(c-a > b-c) {aa=a; bb=c-r;}
	else {aa=c+r; bb=b;}

	npt2=0; dif2=0.0; ri2=0.0; ier1=0;
	if(fabs(bb-aa)>aeps) ier1=adpint(&ri2,aa,bb,reps,aeps,&dif2,f,&npt2,&nmax);

	*ri=ri1+ri2;
 
/*     FUNCTION FUNP REQUIRES TWO EVALUATIONS OF FUN */
	*npt=2*npt1+npt2;
	*dif=fabs(dif1)+fabs(dif2);
	return ier+2*ier1;
}

/*
	---------------------------
	double funp(double x)

	{
		return f(CAU_CC+x)+f(CAU_CC-x);
	}
*/

