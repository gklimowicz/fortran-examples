/*	To perform a line search for minimum of several variables as required
	by direction set method
	This function should not be used for any other purpose

	X0 : (input/output) Starting value for the line search. After execution
		it should contain the distance to minimiser along the line
	X1 : (input/output) Initial estimate for the minimum along the line.
		 This value will be modified by the function.
	F0 : (input/output) The function value at X1, this value must be supplied
	REPS : (input) Required relative accuracy 
	AEPS : (input) Required absolute accuracy,
		This criterion is only used to terminate line search under
		certain conditions and not generally applicable to the line
		search. These values should be same as what is used in function NMINF
	F : (input) Name of the function routine to calculate the function value
	V : (input/output) Array of length 2N. First N element specify
		the direction in which minimisation is required. The next
		N elements will contain the coordinates of the minimiser
		found by LINMNF.
	XI : (input) Array of length N containing the coordinates of
		starting point for line search
	N : (input) Number of variables in the function to be minimised
	NUM : (output) Integer variable to keep count of the number of
		function evaluations used so far
		
	Error status is returned by the value of the function LINMNF.
		0 value implies successful execution
		56 implies that the function failed to find acceptable point
			in this case the starting point itself is accepted

	Function F(N,X,FX) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, FX is the
		function value at X. X is an array of length N.

	Required functions :  FLN, F
*/

#include <math.h>

double fln(void fcn(int , double * , double * ), double x, double v[],
	double x0[], int n, int *num);

int linmnf(double *x0, double *x1, double *f0, double reps, double aeps,
	void (*f) (int , double * , double * ), double v[], double xi[],
	int n, int *num)

{
	int i,nit=25;
	double f1,x2,f2,r,r1,r2,t,p,xp,fp;
	double t1=9.0, gr=1.618034, gc=0.381966;

	f1=fln(f,*x1,v,xi,n,num);
	if(f1<(*f0)) x2=(*x0)+gr*(*x1-(*x0));
	else x2=(*x0)-gr*(*x1-(*x0));
	f2=fln(f,x2,v,xi,n,num);

	for(i=1; i<=nit; ++i) {
/*	Parabolic interpolation */
		r=(x2-(*x1))*(f2-(*f0));
		t=(x2-(*x0))*(f2-f1);
		p=(x2-(*x0))*t-(x2-(*x1))*r;
		t=2.*(t-r);
		if(t>0.0) p=-p;
		t=fabs(t);
		xp=(*x0);
		if(fabs(p)<fabs(t1*t*(x2-(*x1)))) {
/*	Try the interpolated value */
			xp=x2+p/t;
			fp=fln(f,xp,v,xi,n,num);
			if(fp<(*f0)) {
/*	accept the point */
				*x0=xp; *f0=fp;
				return 0;
			}
		}

		if(f1<(*f0)) {*x0=(*x1); *f0=f1; return 0;}
		else if(f2<(*f0)) {*x0=x2; *f0=f2; return 0;}
		else if(xp != (*x0)) {
/*	subdivide the interval */
			if( (xp-(*x0) > 0.0) == (*x1-(*x0) > 0.0) ) {
				if(fabs(xp-(*x0))<0.5*fabs(*x1-(*x0))) {*x1=xp; f1=fp;}
				else {
/*	use golden section */
					*x1=(*x0)+gc*(*x1-(*x0));
					f1=fln(f,*x1,v,xi,n,num);
				}
			} else {
				if(fabs(xp-(*x0))<0.5*fabs(x2-(*x0))) {x2=xp; f2=fp;}
				else {
/*	use golden section */
					x2=(*x0)+gc*(x2-(*x0));
					f2=fln(f,x2,v,xi,n,num);
				}
			}

/*	If interpolated point is not acceptable use golden section */
		} else if(fabs(x2-(*x0))>fabs(*x1-(*x0))) {
			x2=(*x0)+gc*(x2-(*x0));
			f2=fln(f,x2,v,xi,n,num);
		} else {
			*x1=(*x0)+gc*(*x1-(*x0));
			f1=fln(f,*x1,v,xi,n,num);
		}

/*	If the change in function value is too small, then quit */
		r1=reps*fabs(*f0); if(aeps>r1) r1=aeps;
		r2=f1-(*f0); if(f2-(*f0) < r2) r2=f2-(*f0);
		if(r2<r1) return 0;
	}

	return 56;
}
