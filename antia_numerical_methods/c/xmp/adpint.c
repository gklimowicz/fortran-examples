/*	Adaptive integration over a finite interval */
/*	This program uses GAUS16 instead of KRONRD as in quad.f
	and hence will be less efficient */

#include <stdio.h>
#include <math.h>

double fun(double x);
int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double xl, xu, rint, reps, aeps, dif;

/*	Example 6.1 */

	aeps=1.e-18; reps=1.e-13;
	nmax=0;
	for(i1=0; i1<99; ++i1) {
		printf("type xl = lower limit, xu = upper limit\n");
		printf("             quits when xl=xu)\n");
		scanf(" %le %le", &xl, &xu);
		if(xl==xu) return 0;

		i=adpint(&rint,xl,xu,reps,aeps,&dif,fun,&nuse,&nmax);
		printf("ier = %d  lower limit = %e  upper limit =  %e \n", i,xl,xu);
		printf("No. of function evaluations = %d   integral = %e   estimated error = %e\n",nuse,rint,dif);
	}
	return;
}

/*  Specify the integrand */
double fun(double x)

{
	return sqrt(x); 
}

 

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
/*		i=kronrd(&fint,rl,ru,&dif0,&np,f); */
		i=gaus16(&fint,rl,ru,&dif0,&np,f);  
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
/*	i=kronrd(&fint,rl,ru,&dif0,&np,f); */
	i=gaus16(&fint,rl,ru,&dif0,&np,f); 
	*npt=(*npt)+np;
	*rint=(*rint)+fint;
	*dif=(*dif)+dif0;
	return 325;
}




/*	To integrate a function over a finite interval using 16 point
	Gauss-Legendre formula, for use with ADPINT

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

int gaus16(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ))

{
	int k;
	double at,bt,fbt,r1,f1,f2;

/*	W8 and A8 are the weights and abscissas for the 8-point Gauss formula
	W16 and A16 are the weights and abscissas for the 16-point Gauss formula
	Because of symmetry only half the points are given.
*/

	double w8[4] = {0.10122853629037625915e0, 0.22238103445337447054e0,
                        0.31370664587788728734e0, 0.36268378337836198297e0};
	double a8[4] = {0.96028985649753623168e0, 0.79666647741362673959e0,
                        0.52553240991632898582e0, 0.18343464249564980494e0};
	double w16[8] ={0.02715245941175409485e0, 0.06225352393864789286e0,
                        0.09515851168249278481e0, 0.12462897125553387205e0,
                        0.14959598881657673208e0, 0.16915651939500253819e0,
                        0.18260341504492358887e0, 0.18945061045506849629e0};
	double a16[8] ={0.98940093499164993260e0, 0.94457502307323257608e0,
                        0.86563120238783174388e0, 0.75540440835500303390e0,
                        0.61787624440264374845e0, 0.45801677765722738634e0,
                        0.28160355077925891323e0, 0.09501250983763744019e0};

	at=(b-a)/2.;
	bt=(b+a)/2.;
	r1=0.0;
/*	8-point Gauss-Legendre formula */
	for(k=0; k<4; ++k) r1=r1+w8[k]*(f(bt+at*a8[k])+f(bt-at*a8[k]));

	*ri=0.0;
/*	16-point Gauss-Legendre formula */
	for(k=0; k<8; ++k) *ri=(*ri)+w16[k]*(f(bt+at*a16[k])+f(bt-at*a16[k]));

	*ri=(*ri)*at;
	r1=r1*at;
	*dif=fabs(*ri-r1);
	*n=24;
	return 0;
}
