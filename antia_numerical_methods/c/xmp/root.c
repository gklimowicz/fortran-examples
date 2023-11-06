/*	PROGRAM FOR FINDING REAL ROOTS OF A NONLINEAR EQUATION IN A
	SPECIFIED INTERVAL.  THE ROOTS ARE LOCATED BY LOOKING FOR SIGN
	CHANGES IN THE FUNCTION VALUE AT THE SPECIFIED SEARCH STEP. 
	AFTER THAT SECANT ITERATION OR BRENT'S METHOD IS USED TO CALCULATE
	THE ROOT.
*/


#include <stdio.h>
#include <math.h>

double f(double x);
int brent(double *a, double *b, double *x, double reps, double aeps,
	double (*f) (double ));
int secant(double x0, double xl, double xu, double *x, double reps,
	double aeps, double (*fun) (double ));


main()
{
    int i,it,ier;
    double dx,xl,xu,x0,x1,x,f1,fx,a,b,aeps,reps;

    reps=1e-7; aeps=1e-8;

    for(i=0;i<99;++i){
	printf("Type dx=Searc step, xl=lower limit, xu=upper limit\n");
	printf("		(Quits when dx(xl-xu)>0)\n");
	scanf("%le %le %le",&dx,&xl,&xu);
	if(dx*(xl-xu)>0) return;
	printf("Type it=1/2    for   secant/brent\n");
	scanf("%d",&it);
	printf("Search step = %e   for roots in (%e , %e)  it = %d\n",dx,xl,xu,it);

/*	LOCATE THE ROOTS BY LOOKING FOR SIGN CHANGES */

	f1=f(xl);
	for(x1=xl+dx; x1<=xu; x1=x1+dx){
		fx=f(x1);
/*	to eliminate the singularity check the sign of tan(x) also */
		if(fx*f1<0 && tan(x1)*tan(x1-dx)>0) {

/*	Find accurate value of the root  */
			if(dx>0) { a=x1-dx; b=x1;}
			else {b=x1-dx; a=x1;}
			x0=(a+b)/2;
			if(it==1) {
				ier=secant(x0,a,b,&x,reps,aeps,f);
				printf("ier = %d   root = %e   f(x) = %e\n",ier,x,f(x));}
			else {
				ier=brent(&a,&b,&x,reps,aeps,f);
				printf("ier = %d   root = %e   f(x) = %e\n",ier,x,f(x));}
		}
		f1=fx;
	}
    }
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


/*	THE REQUIRED FUNCTION WHOSE ZERO IS REQUIRED */

double f(double x)

{
	return x-tan(x);
}


