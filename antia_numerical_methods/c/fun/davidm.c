/*	To minimise a function in one dimension using Hermite cubic interpolation

	X1,X2 : (input/output) Two starting values for iteration
		These will be updated during execution
		with X2 giving the best estimate for minimiser.
	F2 : (output) The function value at X2, which should be the minimum
	D2F : (output) Estimated value of second derivative at X2, to detect
		if the extremum is a minimum or maximum. If D2F>0 then
		X2 should be minimum, otherwise it will be maximum
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The Minimiser will be calculated with accuracy max(AEPS, REPS*fabs(X2))
	F : (input) Name of the function to calculate the function
		which is to be minimised
		
	Error status is returned by the value of the function DAVIDM.
		0 value implies successful execution
		52 implies that iteration has converged to a maximum
		502 implies that initial values X1=X2 and no calculations are done
		524 implies that iteration cannot be continued 
		525 implies that iteration failed to converge to specified accuracy

	Function F(X,DX) to calculate the required function, must be supplied
		by the user. Here DX is the first derivative of F(X) at X

	Required functions : F
*/

#include <math.h>

int davidm(double *x1, double *x2, double *f2, double *d2f, double reps,
	double aeps, double (*f) (double , double *))

{
	int i, nit=75;
	double f1,dx,dx1,df12,r,r1,r2,df1,df2;

	if(x1==x2) return 502;

	f1=f(*x1,&df1);
	*f2=f(*x2,&df2);
	dx1=(*x2)-(*x1);

/*	Loop for iteration */
	for(i=1; i<=nit; ++i) {
		df12=(*f2-f1)/dx1;
		r=2.*df2+df1-3.*df12;
		r1=3.*df12-df2-df1; r1=r1*r1-df1*df2;
		if(r1>0.0) {
/*	Perform Hermite cubic interpolation */
			if(r>=0.0) r1=sqrt(r1);
			else r1= -sqrt(r1);
			r=r+r1;
			dx=-df2*dx1/r;
		}
		else {
/*	Perform parabolic interpolation */
			r2=2.*(df12-df2);
			if(r2 != 0.0) dx=dx1*df2/r2;
			else return 524;
		}

		*x1=(*x2);
		*x2=(*x2)+dx;
		*d2f=(df2-df1)/dx1;
		dx1=dx;
		f1=(*f2);
		df1=df2;
		*f2=f(*x2,&df2);

/*	Convergence check */
		r1=reps*fabs(*x2); if(aeps>r1) r1=aeps;
		if(x1==x2 || (fabs(dx)<r1 && i>2)) {
/*	If it is maximum set error flag */
			if(*d2f<=0.0) return 52;
			else return 0;
		}
	}

	return 525;
}
