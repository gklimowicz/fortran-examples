/*	Minimisation in one dimension using the golden section search 
	It first finds a bracketing triplet to locate a minimum and then
	uses golden section search to find the accurate value */

#include <stdio.h>
#include <math.h>

double fun(double x);
int golden(double *a, double *b, double *x, double *fx, double reps,
	double aeps, double (*f) (double ));
int brackm(double *a, double *b, double *x, double (*f) (double ));

main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double xl, xu, x, reps, aeps, x0, fx;

/*	Example 8.1 */

	aeps=1.e-7; reps=1.e-6;
	for(i1=0; i1<99; ++i1) {
		printf("type xl, xu = initial pts for bracketing    (quits when xl=xu)\n");
		scanf(" %le %le", &xl, &xu);
		if(xl==xu) return 0;
		printf("initial choice for bracketing = %e , %e\n",xl,xu);

		i=brackm(&xl,&xu,&x,fun);
		printf(" ier = %d    bracketing triplet = %e %e %e   \n", i,xl,x,xu);

/*	Use the bracketing triplet xl,x,xu for golden section search */
		i=golden(&xl,&xu,&x,&fx,reps,aeps,fun);
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


/*	To minimise a function in one dimension using Golden section search

	A,B,X : (input/output) Triplet which brackets the minimum.
		After execution these will contain the final triplet
		with X giving the best estimate for minimiser.
		X must be between A and B; and F(X)<MIN(F(A),F(B))
	FX : (output) The function value at X.
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The bracketing interval is subdivided until
		ABS(B-A) < max(AEPS, REPS*fabs(X))
	F : (input) Name of the function routine to calculate the function
		which is to be minimised
		
	Error status is returned by the value of the function GOLDEN.
		0 value implies successful execution
		50 implies that function values at A, X, B are not
			distinct, most likely due to roundoff errors
		51 implies that iteration failed to reduce the
			bracketing interval to required level
		522 implies that initial values of A,X,B do not bracket the minimum

	Function F(X) to calculate the required function, must be supplied
		by the user.

	Required functions : F
*/

#include <math.h>

int golden(double *a, double *b, double *x, double *fx, double reps,
	double aeps, double (*f) (double ))

{
	int i,nit=100;
	double fa,fb,xp,fp,r1, gr=0.381966e0;

	fa=f(*a); fb=f(*b); *fx=f(*x);
	r1=fa; if(fb<r1) r1=fb;
/*	(A, X, B) does not bracket a minimum */
	if((*fx)>r1 || ( (*x-(*a)>0) == (*x-(*b)>0)) ) return 522;

	for(i=1; i<=nit; ++i) {
/*	Choose the next point in the larger section */
		if(fabs(*x-(*a))>fabs(*x-(*b))) {
			xp=(*x)+gr*(*a-(*x));
			fp=f(xp);
			if(fp<(*fx)) {
/*	(A, XP, X) bracket the minimum */
				*b=(*x); *x=xp;
				fb=(*fx); *fx=fp;
			}
			else {
/*	(XP, X, B) bracket the minimum */
				*a=xp; fa=fp;
			}
		}
		else {
			xp=(*x)+gr*(*b-(*x));
			fp=f(xp);
			if(fp<(*fx)) {
/*	(X, XP, B) bracket the minimum */
				*a=(*x); *x=xp;
				fa=(*fx); *fx=fp;
			}
			else {
/*	(A, X, XP) bracket the minimum */
				*b=xp; fb=fp;
			}
		}

		r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
		if(fabs(*b-(*a)) < r1) return 0;
		
/*	Roundoff errors may be dominating, then quit */
		if(((*fx)>=fa) && ((*fx)>=fb)) return 50;
	}

	return 51;
}
