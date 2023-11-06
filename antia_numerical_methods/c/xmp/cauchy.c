/*	To evaluate the Cauchy principal value of an integral */

#include <stdio.h>
#include <math.h>
int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));
int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax);

int cauchy(double *ri, double a, double b, double c, double reps, double aeps,
	double *dif, double (*f) (double ), double (*funp) (double ), int *npt);

double fun(double x);
double f2(double x);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np;
	double xl, xu, c, rint, reps, aeps, df;

/*	Example 6.26 */

	aeps=1.e-15; reps=1.e-13;
	for(i1=0; i1<99; ++i1) {
		printf("type  xl=lower limit    xu=upper limit    c=singularity\n");
		printf("                          (quits when xl > xu)\n");
		scanf(" %le %le %le", &xl, &xu, &c);
		if(xl>xu) return 0;

		i=cauchy(&rint,xl,xu,c,reps,aeps,&df,fun,f2,&nuse);
		printf(" ier = %d  no. of function evaluations = %d  interval = %e , %e \n", i,nuse,xl,xu);
		printf(" c = %e   integral = %e   error estimate = %e  \n", c,rint, df);
	}
	return;
}

double fun(double x)

{
	return exp(x)/x;
}




	double f2(double x)

{
/*	return fun(CAU_CC+x)+fun(CAU_CC-x);  */
	return 2*sinh(x)/x; 
}



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



/*	To integrate a function over finite interval using adaptive control
	of step size

	RINT : (output) Calculated value of the integral
	XL : (input) The lower limit
	XU : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	NPT : (output) Number of function evaluations
	NMAX : (input/output) Maximum number of function evaluations to be tried
		If NMAX<=0 it is set to MAXPT (=100000)

	Error status is returned by the value of the function ADPINT.
		0 value implies successful execution
		31 implies specified accuracy was not achieved on
			at least one subinterval
		32 implies that this failure occurred more than IFMAX (=5) times
		325 implies that function failed to attain required
			accuracy using NMAX function evaluations
		In all cases DIF will contain the estimated accuracy

		Function F(X) must be supplied by the user.

	Required functions : KRONRD (or GAUS16), F
*/

#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));
int gaus16(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));


int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax)

{
	int i,j,k,ier,ifail,iu,np,q, ipmax=100, ifmax=5, maxpt=100000;
	double aepsl,rl,ru,rm,fint,dif0,r1,xu1[100];
 
	*rint=0.0; *dif=0.0;
	if(xl == xu) return 0;
	ifail=0; ier=0;
	if(*nmax <= 0) *nmax = maxpt;
	aepsl= aeps;
	*npt=0;
	rl=xl; ru=xu; iu=0;

	do {
/*	To evaluate the integral over [rl,ru] */
		i=kronrd(&fint,rl,ru,&dif0,&np,f);
/*		i=gaus16(&fint,rl,ru,&dif0,&np,f);  */
		*npt = (*npt)+np;
		rm=0.5*(rl+ru);
/*	q=.TRUE. if the interval cannot be divided further */
		q = ((iu>=ipmax) || (rm == rl) || (rm == ru));
		r1=fabs(fint)*reps; if(aepsl>r1) r1=aepsl;

		if( (dif0 < r1) || q) {
/*	Accept the value of fint if adequate convergence or if the interval
	cannot be subdivided further */
			*rint=(*rint)+fint;
			*dif=(*dif)+dif0;
			r1=fabs(*rint)*reps; if(aepsl>r1) r1=aepsl;

			if(q && (dif0>r1)) {
/*	Integration fails to converge on this subinterval. Go to the next subinterval */
				ier=31; ifail=ifail+1;
/*	If failure is frequent then adjust the convergence criterion. */
				if(ifail > ifmax) {ier=32; aepsl=(*dif)*0.5;}
			}

/*	If all subintervals are exhausted then return */
			if(iu <= 0) return ier;

/*	otherwise try next subinterval */
			rl=ru; ru=xu1[iu]; iu=iu-1;
		}
		else {

/*	Subdivide the current interval and try again */
			iu=iu+1;
			xu1[iu]=ru;
			ru=rm;
		}
	} while(*npt < (*nmax));


/*	If the number of function evaluations has exceeded the limit then
	try a last call to estimate the integral over the remaining interval */
	ru=xu;
	i=kronrd(&fint,rl,ru,&dif0,&np,f);
/*	i=gaus16(&fint,rl,ru,&dif0,&np,f); */
	*npt=(*npt)+np;
	*rint=(*rint)+fint;
	*dif=(*dif)+dif0;
	return 325;
}


/*	To integrate a function over a finite interval using Gauss-Kronrod formula
	For use with ADPINT

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	DIF : (output) estimated (absolute) error achieved by the function
	N : (output) Number of function evaluations used
	F : (input) Name of the function to calculate the integrand

	Function F(X) must be supplied by the user

	The returned value is always zero.

	Required functions : F
*/

#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ))

{
	int k;
	double at,bt,fbt,r1,f1,f2;

/*	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
	WK7 are the weights for these points in Kronrod formula
	WK15 and AK15 are the weights and abscissas for the remaining points
	in Kronrod formula.
	Because of symmetry only half the points are given.
*/

	double w7[4] = {0.12948496616886969327e0, 0.27970539148927666790e0,
                        0.38183005050511894495e0, 0.41795918367346938775e0};
	double a7[4] = {0.94910791234275852452e0, 0.74153118559939443986e0,
                        0.40584515137739716690e0, 0.0};
	double wk7[4] = {0.06309209262997855329e0, 0.14065325971552591874e0,
                         0.19035057806478540991e0, 0.20948214108472782801e0};
	double wk15[4] = {0.02293532201052922496e0, 0.10479001032225018383e0,
                          0.16900472663926790282e0, 0.20443294007529889241e0};
	double ak15[4] = {0.99145537112081263920e0, 0.86486442335976907278e0,
                          0.58608723546769113029e0, 0.20778495500789846760e0};

	at=(b-a)/2.;
	bt=(b+a)/2.;
	fbt=f(bt);
	r1=w7[3]*fbt;
	*ri=wk7[3]*fbt;
	for(k=0; k<3; ++k) {
		f1=f(bt+at*a7[k]);
		f2=f(bt-at*a7[k]);
/*	7-point Gauss-Legendre formula */
		r1=r1+w7[k]*(f1+f2);
/*	15-point Kronrod formula */
		*ri=(*ri)+wk7[k]*(f1+f2);
	}

	for(k=0; k<4; ++k) *ri=(*ri)+wk15[k]*(f(bt+at*ak15[k]) + f(bt-at*ak15[k]));

	*ri=(*ri)*at;
	r1=r1*at;
	*dif=fabs(*ri-r1);
	*n=15;
	return 0;
}
