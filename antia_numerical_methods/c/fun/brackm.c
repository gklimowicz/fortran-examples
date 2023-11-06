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
