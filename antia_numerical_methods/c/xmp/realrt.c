/*	Real roots of a nonlinear equation using bisection, Brent's method
	or secant iteration*/

/*	Bisect and Brent may converge to a singularity if sign change
	is across a singularity, while secant will not do so
*/

#include <stdio.h>
#include <math.h>

double fun(double x);
int bisect(double *x, double *xl, double *xu, int nb, double (*f) (double ));
int secant(double x0, double xl, double xu, double *x, double reps,
	double aeps, double (*fun) (double ));
int brent(double *a, double *b, double *x, double reps, double aeps,
	double (*f) (double ));

main()
{
	int i,i1,j, id, iflg, ier,np;
	double xl, xu, root, reps, aeps, x0, x1, x2;

/*	Example 7.2 */

	aeps=1.e-8; reps=1.e-7;
	for(i1=0; i1<99; ++i1) {
		printf("type np=No. of bisections,  xl=lower limit,  xu=upperlimit, x0=initial guess\n");
		printf("                            (quits when xl=xu)\n");
		scanf(" %d %le %le %le", &np, &xl, &xu, &x0);
		if(xl==xu) return 0;
		x1=xl;  x2=xu;

		printf("  No. of bisections = %d   lower limit = %e   upper limit = %e\n",np,xl,xu);
		i=bisect(&root,&xl,&xu,np,fun);
		printf(" bisect: ier = %d  root = %e  final interval = %e , %e \n", i,root,xl,xu);

/*	restore the limits */
		xl=x1;  xu=x2;
		i=secant(x0,xl,xu,&root,reps,aeps,fun);
		printf(" secant: ier = %d   initial guess = %e   root = %e   \n", i,x0,root);

		i=brent(&xl,&xu,&root,reps,aeps,fun);
		printf(" brent: ier = %d  root = %e  final interval = %e , %e \n", i,root,xl,xu);
	}
	return;
}

double fun(double x)

{
	return x-tan(x);
}



/*	Real zero of a continuous function using bisection

	X : (output) Computed value of the zero using interpolation
	XL : (input/output) Lower limit of interval containing the zero
	XU : (input/output) Upper limit of interval containing the zero
		These limits will be replaced by refined limits after execution
	NB : (input) Number of bisections to be performed
	F : (input) Name of the function routine to calculate the function

	Error status is returned by the value of the function BISECT.
		0 value implies successful execution
		-1 implies that function value is zero at some trial point
		401 implies NB<=0 in which case no calculations are done
		421 implies that function has same sign at both limits
		and it is not possible to apply bisection

	Function F(X) must be supplied by the user.

	Required functions : F
*/

#include <math.h>

int bisect(double *x, double *xl, double *xu, int nb, double (*f) (double ))

{
	int k;
	double fl,fu,fx;

	if(nb<=0) return 401;

	fl=f(*xl); fu=f(*xu);
/*	If the function value is zero at either point then quit */
	if(fl==0.0) {*x=(*xl); return -1;}
	else if(fu==0) {*x=(*xu); return -1;}

/*	If the function has the same sign at both end points then quit */
	if((fl>0.0) == (fu>0.0)) return 421;

/*	Loop for bisection */
	for(k=1; k<=nb; ++k) {
		*x=(*xl+(*xu))/2.;
		fx=f(*x);
/*	If function is zero then quit */
		if(fx == 0.0) return -1;
		if((fx>0.0) == (fu>0.0)) {*xu=(*x); fu=fx;}
		else {*xl=(*x); fl=fx;}
	}

/*	linear interpolation between XL and XU */
	*x=((*xl)*fu-(*xu)*fl)/(fu-fl);
	return 0;
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
