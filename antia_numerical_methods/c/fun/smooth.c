/*	To draw a smooth curve passing through a set of data points
	    using cubic spline interpolation

	NTAB : (input) Number of points in the table
	X : (input) Array of length NTAB containing X values
	F : (input) Array of length NTAB containing function values at X[I]
	C : (output) Array of length 3*NTAB which will contain the spline coefficients
	NP : (input) Number of points at which interpolation is to be calculated
	XP : (output) Array of length NP containing the x values at
	           NP uniformly spaced points for use in plotting
	FP : (output) Array of length NP containing interpolated
		function values at XP[I]
		
	Error status is returned by the value of the function SMOOTH.
		0 value implies successful execution
		202 implies NP<=1
		other values may be set by SPLINE

	Arrays XP and FP can be used to draw a smooth curve through the
	tabulated points.

	Required functions : SPLINE, SPLEVL
*/

#include <math.h>

int spline(double x[], double f[], int n, double c[][3]);
double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier);

int smooth(int ntab, double x[], double f[], double c[][3], int np,
	double xp[], double fp[])

{
	int i,ier;
	double dx,dfb,ddfb;

	i=spline(x,f,ntab,c);
	if(i>100) return i;

	if(np <= 1) return 202;

	dx=(x[ntab-1]-x[0])/(np-1);
	for(i=0; i<np; ++i) {
		xp[i]=x[0]+dx*i;
		fp[i]=splevl(xp[i], ntab, x, f, c, &dfb, &ddfb, &ier);
	}
	return 0;
}
