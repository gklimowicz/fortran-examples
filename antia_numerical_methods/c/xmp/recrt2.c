/*	Recursive solution of a system of 2 nonlinear equations 
	It scans for zeros along the first variable by looking for sign
	changes at specified search step */

#include <stdio.h>
#include <math.h>

double fun(double x);
double fy(double y);
int brent(double *a, double *b, double *x, double reps, double aeps,
	double (*f) (double ));

double XX,YY,FY;

main()
{
	int i,i1,j, id, iflg, ier,np;
	double a, b, xl, xu, root, reps, aeps, x0, x1, x2, dx, f1, fx;

/*	Exercise 7.49  (I) */

	aeps=1.e-8; reps=1.e-7;
	for(i1=0; i1<99; ++i1) {
		printf("type dx=search step,  xl=lower limit,  xu=upperlimit\n");
		printf("                            (quits when xl=xu)\n");
		scanf(" %le %le %le", &dx, &xl, &xu);
		if(xl==xu) return 0;
		printf("search step = %e     lower limit = %e    upper limit = %e\n",dx,xl,xu);

/*	locate zeros by looking for sign changes */
		f1=fun(xl);
		for(x1=xl+dx; x1<=xu; x1=x1+dx) {
			fx=fun(x1);
			if(fx>0 != f1>0) {

/*	Find accurate value of root using brent */
				a=x1-dx;
				b=x1;
				i=brent(&a,&b,&root,reps,aeps,fun);
				printf(" ier = %d  root = %e  %e    function values = %e , %e \n", i,root,YY,fun(root),FY);
			}
			f1=fx;
		}
	}
	return;
}




/*	Real zero of a given function using Brent's method

	A : (input/output) Lower limit of interval where zero is located
	B : (input/output) Upper limit of interval where zero is located
		The limits A, B are updated by the function
	X : (output) Computed value of the zero
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))
	F : (input) Name of the function routine to calculate the function

	Error status is returned by the value of the function BRENT.
		0 value implies successful execution
		427 implies that function has same sign at x=a and b
			in which case no calculations are done
		428 implies that iteration failed to converge to specified
			accuracy

		Function F(X) must be supplied by the user.

	Required functions : F
*/

#include <math.h>

int brent(double *a, double *b, double *x, double reps, double aeps,
	double (*f) (double ))

{

	int i,nit=75;
	double fa,fb,c,fc,d,e,dx,r1,r2,r3,p,p1,eps;

	fa=f(*a); fb=f(*b);
/*	If the function has the same sign at both the end points, then quit */
	if( (fa>0.0) == (fb>0.0)) return 427;

	c=(*a); fc=fa;
	d=(*b)-(*a); e=d;

/*	The iteration loop */
	for(i=1; i<=nit; ++i) {
		if(fabs(fc) < fabs(fb)) {
/*	If f(c) < f(b), then interchange b and c */
			*a=(*b); *b=c; c=(*a);
			fa=fb; fb=fc; fc=fa;
		}

		eps=reps*fabs(*b); if(aeps>eps) eps=aeps;
		dx=0.5*(c-(*b));
		if(fabs(dx)<eps || fb==0.0) {*x=(*b); return 0;}

/*	If the previous step is too small or if f(a)<=f(b), perform bisection */
		if(fabs(e)<eps || fabs(fa)<=fabs(fb)) {d=dx; e=dx;}
		else {
			r3=fb/fa;
			if((*a)==c) {
/*	Try linear interpolation */
				p=2.*dx*r3;
				p1=1.-r3;
			}
			else {
/*	Try inverse quadratic interpolation */
				r1=fa/fc; r2=fb/fc;
				p=r3*(2.*dx*r1*(r1-r2)-(*b-(*a))*(r2-1.0));
				p1=(r1-1.0)*(r2-1.0)*(r3-1.0);
			}

			if(p>0.0) p1=-p1;
			else p=-p;

			e=d;
			if(2.0*p < (3.*dx*p1-fabs(eps*p1)) && p<fabs(0.5*e*p1)) d=p/p1;
			else {d=dx; e=dx;}
		}

		*a=(*b); fa=fb;
		if(fabs(d)>eps) *b=(*b)+d;
		else {
			if(dx>=0.0) *b=(*b)+eps;
			else *b=(*b)-eps;
		}
		fb=f(*b);
		if((fb>0.0) == (fc>0.0)) {
/*	If f(b) and f(c) have the same sign, then replace c by a */
			c=(*a); fc=fa;
			d=(*b)-(*a); e=d;
		}
	}

	*x=(*b);
	return 428;
}





double fun(double x)

{
	double reps, aeps, a, b, y;
	int i,j;

	reps=1.e-7; aeps=1.e-8;

/*	limits on y  */
	a=-1.57;  b=1.57;

/*	Store the value of x in global variable XX for use by fy */
	XX=x;

/*	For a given value of x solve the second eq. fy=0 for y */

	i=brent(&a,&b,&y,reps,aeps,fy);
	if(i>0) exit(1);

/*	Store the value of y and fy(y) in global variables YY and FY */
	YY=y;
	FY=fy(y);

/*	The first equation */
	return sin(4/(x*x+y*y+0.1));
}


/*	The second equation */

double fy(double y)

{
	return XX-tan(y);
}
