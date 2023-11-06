/*	Real roots of a nonlinear equation using Newton-Raphson method */

#include <stdio.h>
#include <math.h>

double fun(double x, double *dx);
int newrap(double x0, double xl, double xu, double *x, double reps,
	double aeps, double (*fun) (double , double * ));

main()
{
	int i,i1,j, id, iflg, ier,np;
	double xl, xu, root, reps, aeps, x0;

/*	Example 7.4 */

	aeps=1.e-9; reps=1.e-7;
	for(i1=0; i1<99; ++i1) {
		printf("type  xl = lower limit  xu = upper limit   x0 = initial guess\n");
		printf("                          (quits when xl=xu)\n");
		scanf("%le %le %le", &xl, &xu, &x0);
		if(xl==xu) return 0;

		i=newrap(x0,xl,xu,&root,reps,aeps,fun);
		printf(" ier = %d   initial guess = %e    root = %e   \n", i,x0,root);

	}
	return;
}

double fun(double x, double *df)

{
	*df=cos(x);
	*df=1.-1./((*df)*(*df));
	return x-tan(x);
}



/*	Real zero of a given function using Newton-Raphson iteration

	X0 : (input) Initial guess for the zero
	XL : (input) Lower limit of interval where zero is expected
	XU : (input) Upper limit of interval where zero is expected
	X : (output) Computed value of the zero
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))
	FUN : (input) Name of the function routine to calculate the function
		
	Error status is returned by the value of the function NEWRAP.
		0 value implies successful execution
		-k implies that the zero is multiple and iteration
			is modified to account for multiplicity k.
		126-100*k implies that the derivative is zero and the zero is
			detected to be multiple. Iteration is terminated.
		403 implies XL>X0 or XU<X0, in which case no calculations are done
		424 implies that iteration goes outside the specified limits
		425 implies that iteration failed to converge to specified accuracy
		426 implies that the derivative is zero and iteration is terminated

	Function FUN(X,DX) must be supplied by the user.
		Here DX is the first derivative of FUN at X.

	Required functions : FUN
*/

#include <math.h>

int newrap(double x0, double xl, double xu, double *x, double reps,
	double aeps, double (*fun) (double , double * ))

{
	int l,l0,mr,mr1,mr2,ier,nit=75;
	double dx,dxr,dx1,dxr1,f,df,r1;

	if(xl>x0 || xu<x0) return 403;

	mr=1; l0=1; ier=0;
	*x=x0; dx=1.0; dxr=1.0;

	for(l=1; l<=nit; ++l) {
		f=fun(*x,&df);
		dx1=dx;
		dxr1=dxr;

		if(df==0.0) {
/*	If the derivative is zero, then quit */
			if(f==0.0) return ier;
			else {
				ier=426; if(mr>1) ier=126-100*mr;
				return ier;
			}
		}

/*	Newton-Raphson iteration for root with multiplicity MR */
		dx=-mr*f/df;
		dxr=dx/dx1;
		*x=(*x)+dx;

		r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
		if(fabs(dx) < r1 && l>2) return ier;
/*	If the iteration goes outside the specified limits, then quit */
		if(*x<xl || (*x)>xu) return 424;

		if(l-l0>3) {
/*	If 3 iterations have been done with same multiplicity, then
	get a new estimate for multiplicity MR */
			mr1=mr;
			mr2=mr;
			if(dxr<0.99) mr1=mr/(1.0-dxr)+0.5;
			if(dxr1<0.99) mr2=mr/(1.0-dxr1)+0.5;
/*	Accept the new value of MR only if both estimates match */
			if(mr1==mr2 && mr1 != mr) {
				ier=-mr1; l0=l;
				mr=mr1;
				if(mr==1) ier=0;
			}
		}
	}

	return 425;
}
