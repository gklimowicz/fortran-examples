/*      Real zero of a given function using secant iteration
	Function is calculated as FX*2**JF
	This function uses reverse communication to calculate function
	values. If IER<0 the function should be evaluated and SECANI should
	be called back with new function value. Calculation of function
	value should not change any other variables in the call statement.

	X0 : (input) Initial guess for the zero
	XL : (input) Lower limit of interval where zero is expected
	XU : (input) Upper limit of interval where zero is expected
	X : (output) Value of x at which the function evaluation is required.
		If IER=0 then it will contain the final value of zero computed
		by the function.
	F : (input) Calculated value of the function at X.
		If function exits with IER<0, then the calling function should
		calculate the function value at X and call SECANI with this value
		stored in F and JF. Other variables should not be changed.
	JF : (input) The exponent of function value, the function value
		should be F*2**JF
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
     		The estimated error should be less than max(AEPS, REPS*fabs(X))
	IER : (input/output) Error parameter, IER=0 implies successful execution
		Before the first call IER should be set to zero
		IER<0 implies that execution is not over and the function needs
			a new function evaluation at X. After calculating the
			function value SECANI should be called back.
     		IER=40 implies that function value is equal at two points
     			and it is not possible to continue the iteration
     		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
     		IER=422 implies that iteration goes outside the specified limits
     		IER=423 implies that iteration failed to converge to specified accuracy


	Required functions : None (Function value is calculated by the calling program)
*/

#include <math.h>


void secani(double x0, double xl, double xu, double *x, double f, int jf,
	double reps,double aeps, int *ier)

{
	int nis=75;
	static int l,jf1;
	static double f1,r1,dx,dx1;

	if(*ier==-1) goto func;

	if(xl>x0 || xu<x0) {*ier=402; return;}

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
	f1=0.0; jf1=0; l=0;

loop:	l=l+1; *ier=-1;
/*	To evaluate the function at x */
		return;

func:	dx1=dx;
		f1=f1*pow(2.0, (double) (jf1-jf));

		if(f1-f == 0.0) {
			if(f == 0.0) {*ier = 0; return ;}
/*	If F1=F and F!=0, then quit */
			else {*ier=40; return;}
		}

/*	The secant iteration */
		if(l>1) dx=dx1*f/(f1-f);
		*x=(*x)+dx;
		f1=f; jf1=jf;

		r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
		if(fabs(dx)<r1 && l>2) {*ier=0; return;}
		if(*x<xl || (*x)>xu) {*ier=422; return;}
	if(l<nis) goto loop;

	*ier=423; return;
}
