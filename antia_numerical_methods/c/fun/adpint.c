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
