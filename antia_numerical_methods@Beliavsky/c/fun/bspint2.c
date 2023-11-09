/*	To calculate coefficients for B-spline interpolation in 2 dimensions

	NX, NY : (input) Number of entries in the table along X, Y directions
	X, Y : (input) Array of length NX, NY containing the abscissas
	F : (input) Array of length LA*NY containing the function values
	K : (input) Order of B-spline required. K=4 gives cubic B-splines
	AX : (input/output) Array of length LA*3K containing the
		triangular decomposition of equation matrix in band form
		for interpolation along X.
		For IFLG>=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	AY : (input/output) Array of length LA*3K containing the
		triangular decomposition of equation matrix in band form
		for interpolation along Y.
		For IFLG>=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	LA : (input) the second dimension of arrays AX, AY, F, C as specified
			 in calling function, LA >= max(NX,NY)
	C : (output) Array of length LA*NY containing coefficients
		of expansion, which will be calculated by the function
	XF : (input/output) Array of size MX, containing
		the knots used for B-spline calculations along x.
		The knots must be distinct and in ascending order.
		For IFLG>1, this array must be supplied, for other values
		of IFLG it is calculated by the function
	YF : (input/output) Array of size MY, containing
		the knots used for B-spline calculations along y.
		The knots must be distinct and in ascending order.
		For IFLG>1, this array must be supplied, for other values
		of IFLG it is calculated by the function
	MX : (input/output) Number of knots for B-splines along X
		For IFLG>1, this value must be supplied, for other values
		of IFLG it is calculated by the function
	MY : (input/output) Number of knots for B-splines along Y
		For IFLG>1, this value must be supplied, for other values
		of IFLG it is calculated by the function
	IFLG : (input/output) Integer specifying the type of calculation required
		IFLG<2 The matrix will be calculated and solved for coefficients
		IFLG>1 The triangular decomposition of matrix is assumed
			to be available in AX, AY and coefficients C are calculated
	INTX, INTY : (input/output) Integer arrays of length NX, NY
		containing information about
		pivoting during solution of system of linear equations
		For IFLG>1, these arrays must be supplied, for other values
		of IFLG they are calculated by the function
		
	Error status is returned by the value of the function BSPINT2.
		0 value implies successful execution
		nonzero values may be set by BSPLIN, BSPINT or GAUBND

	Required functions : BSPINT, BSPLIN, GAUBND
*/

#include <math.h>
#include <stdlib.h>

int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);

int bspint2(int nx, int ny, double x[], double y[], double *f, int k, double *ax,
	double *ay, int la, double *c, double *xf, double *yf, int *mx, int *my,
	int *iflg, int intx[], int inty[])

{
	int i,j,kb,iflg1,idet;
	double det;
	double *wk;

	if(*iflg < 2) {
		iflg1=1;
		i=bspint(nx,x,f,k,ax,la,c,xf,mx,&iflg1,intx);
		if(i>100) return i;
		iflg1=1;
		i=bspint(ny,y,f,k,ay,la,c,yf,my,&iflg1,inty);
		if(i>100) return i;
		*iflg=2;
	}

	wk=(double *) calloc((size_t) (la*nx), sizeof(double));
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) wk[j+i*la]=f[i+j*la];
	}
	iflg1=2;
	kb=k-1;
	i=gaubnd(ny,kb,nx,ay,wk,&det,&idet,inty,la,&iflg1);
	if(i>100) {free(wk); return i;}
 
	for(j=0; j<ny; ++j) {
		for(i=0; i<nx; ++i) c[i+j*la]=wk[j+i*la];
	}
	iflg1=2;
	i=gaubnd(nx,kb,ny,ax,c,&det,&idet,intx,la,&iflg1);
	free(wk);
	return i;
}
