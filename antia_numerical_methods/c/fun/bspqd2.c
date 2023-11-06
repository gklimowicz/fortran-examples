/*	To calculate the integral of a function defined by B-spline expansion
	over 2 dimensions

	NX : (input) Number of knots for B-spline expansion along x direction
	NY : (input) Number of knots for B-spline expansion along y direction
	X : (input) Array of length NX containing the knots along x direction
		The knots must be distinct and in ascending order.
	Y : (input) Array of length NX containing the knots along y direction
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Array of length IW*(NY+K-2) containing coefficients
		of B-spline expansion. 
	IW : (input) The second dimension of WT as declared in calling function
	XL : (input) Lower limit of integration along x
	XU : (input) Upper limit of integration along x
	YL : (input) Lower limit of integration along y
	YU : (input) Upper limit of integration along y
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values may be set by BSPQD which is called

	BSPQD2 = Integral of \sum WT[J][I]\phi_I(x)\psi_J(y) over [XL,XU] x [YL,YU]

	Required functions : BSPLIN, BSPQD
*/

#include <math.h>
#include <stdlib.h>

double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);

double bspqd2(int nx, int ny, double *x, double *y, int k, double *wt,
	int iw, double xl, double xu, double yl, double yu, int *ier)

{
	int i;
	double ri;
	double *wk;

	wk=(double *) calloc((size_t) (ny+k), sizeof(double));

	for(i=0; i<ny+k-2; ++i) {
/*	Calculate integral along x */
		wk[i]=bspqd(nx,x,k,&wt[i*iw],xl,xu,ier);
		if(*ier>100) {free(wk); return 0.0;}
	}
		
/*	Calculate integral along y */
	ri=bspqd(ny,y,k,wk,yl,yu,ier);
	free(wk);
	return ri;
}
