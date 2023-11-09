/*	To calculate coefficients for B-spline interpolation in n-dimensions
 
	N : (input) Number of dimensions
	NK : (input) Integer array of length N giving the number of
		tabular points in each dimension
	X : (input) Array of length NXD*N containing the abscissas
		X[j][i] is the ith abscissa in jth dimension
	NXD : (input) Second dimension of arrays X, XF, INTX as specified
		in the calling function. NXD must be greater than or equal
		to the maximum of NK[I]
	F : (input) Array of length NK[0]*NK[1]* ... *NK[N-1] containing
		the table of values. The dimensions of F in the calling
		function must exactly match the size of table so that
		there are no gaps in memory allocation.
	K : (input) Order of B-spline required. K=4 gives cubic B-splines
	AX : (output) Array of length NXD*3*K*N containing the
		triangular decomposition of equation matrix in band form
		for each dimension.
	C : (output) Array of length NK[0]*NK[1]*... *NK[N-1]
		containing the coefficients of expansion.
	XF : (output) Array of size NXD*N, containing
		the knots used for B-spline calculations.
		XF[j][i] is the ith knot along jth dimension.
	MK : (output) Integer array of length N containing number of
		knots for B-splines in each dimension.
	INTX : (output) Integer array of length NXD*N containing information
		about pivoting during solution of system of linear equations
		for each dimension.
		
	Error status is returned by the value of the function BSPINTN.
		0 value implies successful execution
		Nonzero values may be set by BSPINT, BSPLIN or GAUBND

	Required functions : BSPINT, BSPLIN, GAUBND
*/

#include <math.h>
#include <stdlib.h>

int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);

int bspintn(int n, int nk[], double *x, int nxd, double *f, int k, double *ax,
	double *c, double *xf, int mk[], int *intx)

{
	int i,i1,i2,j,lj,iflg1,kb,nx,ny,no,ju,nx1,ny1,lj1,j1,idet;
	double det;
	double *wk;

/*	Calculate the triangular decomposition of matrices for interpolation
	along each dimension */
	for(i=0; i<n; ++i) {
		iflg1=1;
		lj=nk[i];
		j=bspint(lj,(x+i*nxd),f,k,(ax+i*nxd*3*k),lj,c,(xf+i*nxd),&mk[i],
			&iflg1,(intx+i*nxd));
		if(j>100) return j;
	}

	kb=k-1;
	nx=1;
	for(i=0; i<n-1; ++i) nx=nx*nk[i];
/*	no=nx*nk[n-1]+1 */
	ny=1;
	ju=n-1;
	wk=(double *) calloc((size_t) (nx*nk[n-1]), sizeof(double));
 
/*	If N is odd interpolate along last dimension outside the loop */
	if((n - 2*(n/2)) == 1) {
		lj=nk[n-1];
/*	Set up the RHS */
		for(i=0; i<nx; ++i) {
			for(j=0; j<nk[n-1]; ++j) wk[j+i*lj]=f[i+j*nx];
		}
		iflg1=2;
		i=gaubnd(nk[n-1],kb,nx,(ax+(n-1)*nxd*3*k),wk,&det,&idet,(intx+(n-1)*nxd),
			lj,&iflg1);
		if(i>100) {free(wk); return i;}

/*	Set up the RHS for next interpolation along N-2 th dimension */
		nx1=nx/nk[n-2];
		ny=nk[n-1];
		lj1=nk[n-2];
		for(i1=0; i1<nk[n-2]; ++i1) {
			for(i=0; i<nx1; ++i) {
				for(j=0; j<ny; ++j)
					c[i1+i*lj1+j*nx1*lj1] = wk[j+i*lj+i1*nx1*lj];
			}
		}
		nx=nx1;
		ju=n-2;
	}
	else {
 
/*	Set up the RHS for interpolation along N-1 th dimension */
		lj=nk[n-1];
		for(i=0; i<nx; ++i) {
			for(j=0; j<nk[n-1]; ++j) c[j+i*lj]=f[i+j*nx];
		}
	}
 
/*	Loop for interpolation in each dimension, each pass
	interpolates along 2 dimensions */
	for(j1=ju; j1>=0; j1=j1-2) {
		iflg1=2;
		lj=nk[j1];
		i=gaubnd(nk[j1],kb,nx*ny,(ax+j1*nxd*3*k),c,&det,&idet,(intx+j1*nxd),
			lj,&iflg1);
		if(i>100) {free(wk); return i;}
		
/*	Set up the RHS for interpolation along the next dimension */
		nx1=nx/nk[j1-1];
		ny1=ny*nk[j1];
		lj1=nk[j1-1];
		for(i1=0; i1<ny; ++i1) {
			for(i2=0; i2<nk[j1]; ++i2) {
				for(i=0; i<nk[j1-1]; ++i) {
					for(j=0; j<nx1; ++j)
						wk[i+j*lj1+i2*nx+i1*nx*nk[j1]]=c[i2+j*lj+i*nx1*lj+i1*nx*lj];
				}
			}
		}
		nx=nx1; ny=ny1;
		iflg1=2;
		lj=nk[j1-1];
		i=gaubnd(nk[j1-1],kb,nx*ny,(ax+(j1-1)*nxd*3*k),wk,&det,&idet,
			(intx+(j1-1)*nxd),lj,&iflg1);
		if(i>100) {free(wk); return i;}

		if(j1==1) {
/*	Store the coefficients in array */
			for(i=0; i<nk[0]; ++i) {
				for(j=0; j<ny; ++j) c[i+j*nk[0]]=wk[i+j*lj];
			}
		}
		else {

/*	Set up the RHS for interpolation along the next dimension */
			lj1=nk[j1-2];
			nx1=nx/nk[j1-2];
			ny1=ny*nk[j1-1];
			for(i1=0; i1<nk[j1-1]; ++i1) {
				for(i2=0; i2<nk[j1-2]; ++i2) {
					for(i=0; i<nx1; ++i) {
						for(j=0; j<ny; ++j)
							c[i2+i*lj1+i1*nx+j*nx*nk[j1-1]]=wk[i1+i*lj+i2*lj*nx1+j*lj*nx];
					}
				}
			}
			nx=nx1; ny=ny1;
		}
	}
	free(wk);
	return 0;
}
