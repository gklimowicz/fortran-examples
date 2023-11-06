/*	To minimise a function in one dimension using Brent's method

	A,B,X : (input/output) Triplet which brackets the minimum.
		After execution these will contain the final triplet
		with X giving the best estimate for minimiser.
		X must be between A and B; and F(X)<min(F(A),F(B))
	FX : (output) The function value at X.
	reps : (input) Required relative accuracy
	aeps : (input) Required absolute accuracy
		The bracketing interval is subdivided until
		fabs(B-A) < MAX(AEPS, REPS*fabs(X))
	F : (input) Name of the function to calculate the function
		which is to be minimised

	Error status is returned by the value of the function BRENTM.
		0 value implies successful execution
		51 implies that function failed to reduce the
			bracketing interval to required level
		523 implies that initial values of A,X,B do not bracket
			the minimum

	Function F(X) to calculate the required function, must be supplied
		by the user.

	Required functions : F
*/

#include <math.h>

int brentm(double *a, double *b, double *x, double *fx, double reps,
	double aeps, double (*f) (double ))

{
	int i, nit=75;
	double fa,fb,t,e,v,w,fv,fw,xm,eps,eps2,p,r,r1,d,u,fu;
	double gr=0.381966e0;

	if(*a>(*b)) {t=(*a); *a=(*b); *b=t;}
	fa=f(*a); fb=f(*b); *fx=f(*x);
	if(fa<(*fx) || fb<(*fx) || (*x)<=(*a) || (*x)>=(*b) ) return 523;

	e=0.0; v=(*x); w=(*x);
	fv=(*fx); fw=(*fx);

/*	Loop for iteration */
	for(i=1; i<=nit; ++i) {
		xm=0.5*(*a+(*b));
		eps2=reps*fabs(*x); if(aeps>eps2) eps2=aeps;
		eps=0.5*eps2;

/*	The convergence test */
		if(fabs(*x-xm)<eps2-0.5*(*b-(*a))) return 0;

		p=0.0; t=0.0; r=0.0;
		if(fabs(e)>eps) {
/*	Parabolic interpolation */
			r=(*x-w)*(*fx-fv);
			t=(*x-v)*(*fx-fw);
			p=(*x-v)*t-(*x-w)*r;
			t=2.*(t-r);
			if(t>0) p=-p;
			t=fabs(t);
			r=e;
			e=d;
		}

		if(fabs(p)<fabs(0.5*t*r) && p>t*(*a-(*x)) && p<t*(*b-(*x))) {
/*	accept the interpolated point */
			d=p/t; u=(*x)+d;
			if(u-(*a)<eps2 || (*b)-u<eps2) {
/*	If it is too close to end points shift it by eps at least */
				d=eps;
				if((*x)>=xm) d=-eps;
			}
		}
		else {
/*	Perform golden section search */
			e=(*b)-(*x);
			if((*x)>=xm) e=(*a)-(*x);
			d=gr*e;
		}
		if(fabs(d)>=eps) u=(*x)+d;
		else {
/*	Shift the point by at least eps */
			if(d>=0) u=(*x)+eps;
			else u=(*x)-eps;
		}
		fu=f(u);

/*	Updating the bracketing triplet */
		if(fu<=(*fx)) {
/*	(a, u, x) is the triplet */
			if(u<(*x)) *b=(*x);
/*	(x, u, b) is the triplet */
			else (*a)=(*x);
			v=w; fv=fw;
			w=(*x); fw=(*fx);
			*x=u; *fx=fu;
		}
		else {
/*	(u, x, b) is the triplet */
			if(u<(*x)) *a=u;
/*	(a, x, u) is the triplet */
			else *b=u;
			if(fu<=fw || w==(*x)) {
				v=w; fv=fw;
				w=u; fw=fu;
			}
			else if(fu<=fv || v==(*x) || v==w) {v=u; fv=fw;}
		}
	}

	return 51;
}
