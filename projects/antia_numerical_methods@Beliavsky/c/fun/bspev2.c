/*	To evaluate a B-spline expansion in 2 dimensions

	NX : (input) Number of knots to define B-splines along 1st dimension
	NY : (input) Number of knots to define B-splines along 2nd dimension
	X : (input) Array of length NX containing the knots.
		The knots must be distinct and in ascending order.
	Y : (input) Array of length NY containing the knots.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	NDERIV : (input) Number of derivatives required
		For NDERIV<=0 only function value is calculated
		For NDERIV=1 first derivative is also calculated
		For NDERIV>1 both first and second derivatives are calculated
	WT : (input) Array of length IW*(NY+K-2) containing coefficients
		of B-spline expansion,
	IW : (input) Second dimension of array WT as defined in the calling function
	X0,Y0 : (input) The point at which expansion has to be evaluated
	DFX : (output) First derivative of function w.r.t. X at X0, Y0
	DFY : (output) First derivative of function w.r.t. Y at X0, Y0
	DFXX : (output) Second derivative of function w.r.t X,X at X0, Y0
	DFXY : (output) Second derivative of function w.r.t X,Y at X0, Y0
	DFYY : (output) Second derivative of function w.r.t Y,Y at X0, Y0
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values of IER may be set by BSPLIN which is called

	BSPEV2 =SUM_{i=1}^{NX+K-2} SUM_{j=1}^{NY+K-2} WT[j][i]\phi_i(X0)\psi_j(Y0)
	where \phi_i(x) are B-spline basis functions on knots X
	and \psi_j(y) are B-spline basis functions on knots Y

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

double bspev2(int nx, int ny, double *x, double *y, int k, int nderiv,
	double *wt, int iw, double x0, double y0, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy, int *ier)

{
	int i,j,lx,ly;
	double f;
	double *wx, *wx1, *wx2, *wy, *wy1, *wy2;

	wx=(double *) calloc((size_t) (nx+k),sizeof(double));
	wy=(double *) calloc((size_t) (ny+k),sizeof(double));
	wx1=(double *) calloc((size_t) (nx+k),sizeof(double));
	wy1=(double *) calloc((size_t) (ny+k),sizeof(double));
	wx2=(double *) calloc((size_t) (nx+k),sizeof(double));
	wy2=(double *) calloc((size_t) (ny+k),sizeof(double));
	*ier=bsplin(x,nx,k,x0,nderiv,wx,wx1,wx2,&lx);
	if(*ier > 100) {
		free(wy2); free(wx2); free(wy1); free(wx1); free(wy); free(wx);
		return 0.0;
	}
	*ier=bsplin(y,ny,k,y0,nderiv,wy,wy1,wy2,&ly);
	if(*ier > 100) {
		free(wy2); free(wx2); free(wy1); free(wx1); free(wy); free(wx);
		return 0.0;
	}

	f=0.0; *dfx=0.0; *dfy=0.0;
	*dfxx=0.0; *dfyy=0.0; *dfxy=0.0;

	for(i=lx; i<=lx+k-1; ++i) {
		for(j=ly; j<=ly+k-1; ++j) {
			f=f+wt[i+j*iw]*wx[i]*wy[j];
			*dfx = *dfx + wt[i+j*iw]*wx1[i]*wy[j];
			*dfy = *dfy + wt[i+j*iw]*wx[i]*wy1[j];
			*dfxx = *dfxx + wt[i+j*iw]*wx2[i]*wy[j];
			*dfyy = *dfyy + wt[i+j*iw]*wx[i]*wy2[j];
			*dfxy = *dfxy + wt[i+j*iw]*wx1[i]*wy1[j];
		}
	}
	free(wy2); free(wx2); free(wy1); free(wx1); free(wy); free(wx);
	return f;
}
