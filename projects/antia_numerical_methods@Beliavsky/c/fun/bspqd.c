/*	To calculate the integral of a function defined by B-spline expansion

	N : (input) Number of knots for B-spline expansion
	X : (input) Array of length N containing the knots.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Array of length N+K-2 containing the coefficients
		of B-spline expansion
	XL : (input) Lower limit of integration
	XU : (input) Upper limit of integration
	IER : (output) Error parameter, IER=0 implies successful execution
		IER=31 implies that XL is outside the range of table
		IER=32 implies that XU is outside the range of table

	BSPQD = Integral of \sum WT[I]\phi_I(x) over [XL,XU]

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
                double db[], double ddb[], int *left);

double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier)

{
	int i,j,k1,jk,lf1,lf2,nderiv;
	double f,f1,f2;
	double *wk1, *wk2, *wk;

	*ier=0;
	if(xl == xu) return 0.0;
 
/*	Calculate the B-spline basis of order K+1 at XL and XU */
	k1=k+1;
	nderiv=0;
	wk=(double *)calloc((size_t) (n+3*k), sizeof(double));
	wk1=(double *)calloc((size_t) (n+k), sizeof(double));
	wk2=(double *)calloc((size_t) (n+k), sizeof(double));

	*ier=bsplin(x,n,k1,xl,nderiv,wk1,wk,wk,&lf1);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	*ier=bsplin(x,n,k1,xu,nderiv,wk2,wk,wk,&lf2);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	if(xl<x[0] || xl>x[n-1]) *ier=31;
	if(xu<x[0] || xu>x[n-1]) *ier=32;

	for(i=0; i<n; ++i) wk[k+i]=x[i];
	for(i=1; i<=k; ++i) {
		wk[k-i]=x[0];
		wk[n+i-1+k]=x[n-1];
	}
 
/*	The sum for x=XL */
	f1=0.0;
	for(i=lf1; i<=lf1+k; ++i) {
		f=0.0;
		for(j=1; j<=i; ++j) f=f+wt[j-1]*(wk[k+j]-wk[j]);
		f1=f1+f*wk1[i]/k;
	}

/*	The sum for x=XU */
	f2=0.0;
	for(i=lf2; i<=lf2+k; ++i) {
		f=0.0;
		for(j=1; j<=i; ++j) f=f+wt[j-1]*(wk[k+j]-wk[j]);
		f2=f2+f*wk2[i]/k;
	}

	free(wk2); free(wk1); free(wk);
	return f2-f1;
}
