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


