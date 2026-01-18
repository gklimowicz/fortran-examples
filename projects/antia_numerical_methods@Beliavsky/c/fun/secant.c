/*	Real zero of a given function using secant iteration

	X0 : (input) Initial guess for the zero
	XL : (input) Lower limit of interval where zero is expected
	XU : (input) Upper limit of interval where zero is expected
	X : (output) Computed value of the zero
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))
	FUN : (input) Name of the function routine to calculate the function
		
	Error status is returned by the value of the function SECANT.
		0 value implies successful execution
		40 implies that function value is equal at two points
			and it is not possible to continue the iteration
		402 implies XL>X0 or XU<X0, in which case no calculations are done
		422 implies that iteration goes outside the specified limits
		423 implies that iteration failed to converge to specified accuracy

	Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>

int secant(double x0, double xl, double xu, double *x, double reps,
	double aeps, double (*fun) (double ))

{
	int i, l,nis=75;
	double dx,f1,f,dx1,r1;

	if(xl>x0 || xu<x0) return 402;

	*x=x0;
/*	Select the increment for the next point X+DX */
	dx=(xu-x0)/100.0;
	r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
	if(fabs(dx)<5.0*r1) dx=(xl-x0)/5.;
	r1=fabs(*x); if(r1<100.0*aeps) r1=100.0*aeps;
	if(fabs(dx)>0.1*r1) {
		if(dx>=0.0) dx=r1;
		else dx=-r1;
	}
	f1=0.0;

	for(l=1; l<=nis; ++l) {
		f=fun(*x);
		dx1=dx;

		if(f1-f == 0.0) {
			if(f == 0.0) return 0;
/*	If F1=F and F!=0, then quit */
			else return 40;
		}

/*	The secant iteration */
		if(l>1) dx=dx1*f/(f1-f);
		*x=(*x)+dx;
		f1=f;

		r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
		if(fabs(dx)<r1 && l>2) return 0;
		if(*x<xl || (*x)>xu) return 422;
	}

	return 423;
}
