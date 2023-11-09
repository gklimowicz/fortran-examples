/*	To evaluate the cubic spline interpolant at a specified point

	XB : (input) point at which interpolation is required
	N : (input) Number of points in the table
	X : (input) Array of length N, containing the abscissas
	F : (input) Array of length N, containing the function values at X[I]
	C : (input) Array of length 3*N containing the spline coefficients
		which should have been calculated using SPLINE
	DFB : (output) First derivative of spline at x=XB
	DDFB : (output) Second derivative of spline at x=XB
	IER : (output) error parameter, IER=0 if execution is successful
		IER=24 implies XB is outside the range of table on higher side
		IER=25 implies XB is outside the range of table on lower side
		IER=201 implies N<2
	SPLEVL will be the interpolated value at x=XB

	Required functions : None
*/

#include <math.h>

double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier)

{
	int i,j,igh,nigh, mid;
	double dx, r1;
	static int low=-1;
	
	if(n<2) {*ier=201; return 0.0;}

	*ier=0;
/*	If the previous value of LOW is inadmissible, set the range to (0,N-1) */
	if(low<0 || low>=n-1) {low=0; igh=n-1;}
	else igh=low+1;

	while((xb<x[low] && xb<x[igh]) || (xb>x[low] && xb>x[igh])) {
/*	Extend the range */
		if((xb>x[low]) == (x[n-1]>x[0])) {
/*	Extend the range on higher side */
			if(igh >= n-1) {*ier=24; low=n-2; break;}
			else {
				nigh=igh+2*(igh-low); if(n-1 < nigh) nigh=n-1;
				low=igh; igh=nigh;
			}
		}

		else {
/*	Extend the range on lower side */
			if(low <= 0) {*ier=25; igh=low+1; break;}
			else {
				nigh=low;
				low=low-2*(igh-low); if(low<0) low=0;
				igh=nigh;
			}
		}
	}


/*	Once the point is bracketed between two tabular points locate it by bisection */
	while((igh-low > 1) && (xb != x[low])) {
		mid=(low+igh)/2;
		if((xb<= x[mid]) == (xb<= x[low])) low=mid;
		else igh=mid;
	}

	dx=xb-x[low];
	r1=((c[low][2]*dx+c[low][1])*dx+c[low][0])*dx+f[low];
	*dfb=(3.0*c[low][2]*dx+2.*c[low][1])*dx+c[low][0];
	*ddfb=6.*c[low][2]*dx+2.*c[low][1];
	return r1;
}
