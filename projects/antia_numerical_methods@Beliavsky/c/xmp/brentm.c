/*	Minimisation in one dimension using Brent's method 
	It first finds a bracketing triplet to locate a minimum and then
	uses Brent's method to find the accurate value */

#include <stdio.h>
#include <math.h>

double fun(double x);
int brackm(double *a, double *b, double *x, double (*f) (double ));
int brentm(double *a, double *b, double *x, double *fx, double reps,
	double aeps, double (*f) (double ));

main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double xl, xu, x, reps, aeps, x0, fx;

/*	Example 8.2 */

	aeps=1.e-7; reps=1.e-6;
	for(i1=0; i1<99; ++i1) {
		printf("type xl, xu = initial pts for bracketing    (quits when xl=xu)\n");
		scanf(" %le %le", &xl, &xu);
		if(xl==xu) return 0;
		printf("initial choice for bracketing = %e , %e\n",xl,xu);

		i=brackm(&xl,&xu,&x,fun);
		printf(" ier = %d    bracketing triplet = %e %e %e   \n", i,xl,x,xu);

/*	Use the bracketing triplet xl,x,xu for Brent's method */
		i=brentm(&xl,&xu,&x,&fx,reps,aeps,fun);
		printf(" ier = %d  minimiser = %e   minimum = %e\n", i,x,fx);
	}
	return;
}

double fun(double x)

{
	return (x*x-0.01)*exp(-10.0*x);
}



/*	To bracket a minimum in one dimension

	A,B : (input/output) Starting values for the search. After
		execution they will contain the limits of bracketing triplet
	X : (output) A variable which will contain the point X between
		A and B, such that F(X) < MIN(F(A),F(B))
	F : (input) Name of the function routine to calculate the function
		which is to be minimised

	Error status is returned by the value of the function BRACKM.
		0 value implies successful execution
		501 implies that A=B, in which case no calculations are done
		521 implies that function failed to find a bracketing triplet

	Function F(X) to calculate the required function, must be supplied
		by the user.

	Required functions : F
*/

#include <math.h>

int brackm(double *a, double *b, double *x, double (*f) (double ))

{
	int i, ntry=200;
	double t,fa,fb,fx,r,p,xp,fp;
	double rmax=10.0, gr=1.618034;

	if(*a == (*b)) return 501;

	fa=f(*a); fb=f(*b);
	if(fb > fa) {
/*	Exchange a and b */
		t=(*a); *a=(*b); *b=t;
		t=fa; fa=fb; fb=t;
	}

/*	First guess for x */
	*x=(*b)+gr*(*b-(*a));
	fx=f(*x);

	for(i=1; i<=ntry; ++i) {
		if(fb<fx) {
/*	The minimum is bracketed by the triplet (a,b,x), exchange b and x */
			t=(*b); *b=(*x); *x=t;
			return 0;
		}

/*	Next guess by parabolic interpolation */
		r=(*b-(*a))*(fb-fx);
		t=(*b-(*x))*(fb-fa);
		p=(*b-(*x))*t-(*b-(*a))*r;
		t=2.*(t-r);
		if(t>0) p=-p;
		t=fabs(t);

		if((t>0) && ( (p>0) == ( (x-b)>0))) {
			if(fabs(p) < fabs(t*(*x-(*b))) ) {
/*	Interpolated point is between b and x */
				xp=(*b)+p/t;
				fp=f(xp);
				if(fp<fx) {
/*	Minimum is bracketed by (b,xp,x) */
					*a=(*b); *b=(*x); *x=xp;
					return 0;
				}
				else if(fp>fb) {
/*	Minimum is bracketed by (a,b,xp) */
					*x = (*b); *b=xp;
					return 0;
				}

/*	If the minimum is not bracketed, then reject the interpolated point */
				xp=(*x)+gr*(*x-(*b));
				fp=f(xp);
			}
			
			else if(fabs(p) < fabs(t*(*x-(*b))*rmax)) {
/*	Interpolated point is between x and the allowed limit */
				xp=(*b)+p/t;
				fp=f(xp);
				if(fp<fx) {
/*	If the minimum is not bracketed, then reject the interpolated point */
					xp=(*x)+gr*(*x-(*b));
					fp=f(xp);
				}
			}
			else {
/*	Reject the interpolated point if it is beyond the allowed limit */
				xp=(*b)+rmax*(*x-(*b));
				fp=f(xp);
			}
		}
		else {
/*	Reject the interpolated point if it is on the other side of b */
			xp=(*x)+gr*(*x-(*b));
			fp=f(xp);
		}

		*a=(*b); *b=(*x); *x=xp;
		fa=fb; fb=fx; fx=fp;
	}

	return 521;
}


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
