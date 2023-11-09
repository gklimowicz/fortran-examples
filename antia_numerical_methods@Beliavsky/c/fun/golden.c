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
