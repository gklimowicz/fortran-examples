/*	To calculate the integral of a function defined by B-spline expansion
	over hyper-rectangular region in N dimensions

	N : (input) Number of dimensions
	NK : (input) Integer array of length N containing the number of knots
		for B-spline representation along each direction
	X : (input) Array of length NXD*N containing the knots along each direction
		X[j][i] is the ith knot along jth dimension
	NXD : (input) Second dimension of array X as declared in the calling function
		NXD must be greater than or equal to maximum of NK[I]
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Array of length
		(NK[0]+K-2)(NK[1]+K-2) ... (NK[N-1]+K-2)
		containing coefficients of B-spline expansion. 
		The coefficients are assumed to be stored in natural
		Fortran order with no gaps in data. In the calling
		function the array should have dimension
		WT[NK[N-1]+K-2]...[NK[1]+K-2][NK(0)+K-2]
	XL : (input) Array of length N containing the lower limits of
		integration along each direction
	XU : (input) Array of length N containing the upper limits of
		integration along each direction
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values may be set by BSPQD which is called

	BSPQDN = Integral over[XL[0],XU[0]] x ... x [XL[N-1],XU[N-1]] of
	\sum WT[i_N-1]...[i_1][i_0]\phi_i0(x_0)\phi_i1(x_1) ... \phi_in-1(x_n-1) 

	Required functions : BSPLIN, BSPQD
*/

#include <math.h>
#include <stdlib.h>

double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);

double bspqdn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double xl[], double xu[], int *ier)
 
{
	int i,j,n1,nk1,n0,na,nb,nt;
	double ri;
	double *wk;

	n1=1;
	for(i=1; i<n; ++i) n1=n1*(nk[i]+k-2);
	wk=(double *) calloc((size_t) (2*n1+6), sizeof(double));
	n0=n1+3;
	
/*	Integration along the first dimension */
	nk1=nk[0]+k-2;
	for(j=0; j<n1; ++j) {
		wk[j]=bspqd(nk[0],&x[0],k,&wt[j*nk1],xl[0],xu[0],ier);
		if(*ier >100) {free(wk); return 0.0;}
	}
 
	na=0; nb=n0;
/*	Integrate along the remaining dimensions */
	for(i=1; i<n; ++i) {
		nk1=nk[i]+k-2;
		n1=n1/nk1;
		for(j=0; j<n1; ++j) {
			wk[nb+j]=bspqd(nk[i],&x[i*nxd],k,&wk[na+j*nk1],xl[i],xu[i],ier);
			if(*ier >100) {free(wk); return 0.0;}
		}
		nt=na; na=nb; nb=nt;
	}
	ri=wk[na];
	free(wk);
	return ri;
}
