/*	To calculate linear least squares fit to B-spline basis functions in 2 dimension
	Weights are assumed to be unity

	NX : (input) Number of data points along x-axis
	NY : (input) Number of data points along y-axis
	X,Y : (input) Arrays of length NX,NY containing the coordinates
		of points at which function values are available
	F : (input) Array of length LA*NY containing the function values
		F[J][I] should be function value at X[I],Y[J]
	K : (input) Order of B-splines required, K=4 gives cubic B-splines
	AX : (output) Array of length IV*2*NX containing the matrix
		U of SVD of the design matrix for fit along x-axis
	AY : (output) Array of length IV*2*NY containing the matrix
		U of SVD of the design matrix for fit along y-axis
	LA : (input) Second dimension of arrays F, C, FY as declared
		in the calling function (LA >= 2*MAX(NX,NY))
	VX : (output) Array of length IV*(MX+K-2) containing the matrix
		V of SVD of the design matrix for fit along x-axis
	VY : (output) Array of length IV*(MY+K-2) containing the matrix
		V of SVD of the design matrix for fit along y-axis
	IV : (input) Second dimension of AX, AY, VX, VY in the calling function
		IV >= max(MX,MY)+K-2
	SIGMAX : (output) Array of length MX+K-2 containing the singular
		values of the design matrix for fit along x-axis
	SIGMAY : (output) Array of length MY+K-2 containing the singular
		values of the design matrix for fit along y-axis
	C : (output) Array of length LA*NY containing the fitted coefficients
		Note that although the number of coefficients is (MX+K-2)*(MY+K-2)
		the rest of array is used as scratch space
	XF : (input) Array of size MX, containing the knots
		along x-axis used for defining B-spline basis functions.
		The knots must be distinct and in ascending order.
	YF : (input) Array of size MY, containing the knots
		along y-axis used for defining B-spline basis functions.
		The knots must be distinct and in ascending order.
	MX : (input) Number of knots for B-splines along x-axis,
		the number of basis functions would be MX+K-2
	MY : (input) Number of knots for B-splines along y-axis,
		the number of basis functions would be MY+K-2
	FY : (output) Array of length LA*NY containing the values of fitted
		function at each of the tabular points
	REPS : (input) Required accuracy for solution of equations using SVD
		singular values less than REPS times maximum will be set to zero
	RLM : (input) Parameter lambda for smoothing. If RLM<=0 no smoothing
		is applied
	IDE : (input) Order of derivative to be used for smoothing
		This is used only when RLM>0. IDE=1 for first derivative
		and IDE=2 for second derivative smoothing
	CHISQ : (output) The value of Chi square at minimum

	Error status is returned by the value of the function BSPFIT2.
		0 value implies successful execution
		608 implies that MX+K-2>NX, MY+K-2>NY or K<2
		609 implies that RLM>0 and IDE is not acceptable
		No calculations are done in all these cases
		Other values may be set by SVD or BSPLIN

	Required functions : BSPFIT, BSPLIN, BSPEVL, BSPEV2, SVD, SVDEVL
*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int bspfit(int n, double x[], double f[], double ef[], int k, double *a,
	int la, double *v, int iv, double sigma[], double c[], double xf[],
	int no, double y[], int *iflg, double reps, double rlm, int ide,
	double *chisq, double *cov);
double bspev2(int nx, int ny, double *x, double *y, int k, int nderiv,
	double *wt, int iw, double x0, double y0, double *dfx, double *dfy,
	double *dfxx, double *dfyy, double *dfxy, int *ier);

int bspfit2(int nx, int ny, double x[], double y[], double *f, int k,
	double *ax, double *ay, int la, double *vx, double *vy, int iv,
	double sigmax[], double sigmay[], double *c, double xf[], double yf[],
	int mx, int my, double *fy, double reps, double rlm, int ide, double *chisq)

{
	int i,j,iflg,ier,n1,m,m1,nderiv;
	double r1,dfx,dfy,dfxx,dfxy,dfyy;
	double *wk, *cov;

	wk=(double *) calloc((size_t) (la*nx), sizeof(double));
/*	Set the weights to 1 */
	n1=nx; if(ny>nx) n1=ny;
	cov=(double *) calloc((size_t) (n1*n1), sizeof(double));
	for(i=0; i<=n1; ++i) wk[i]=1.0;

/*	Calculate the SVD of matrix for fit along x-axis */
	iflg=1;
	ier=bspfit(nx,x,f,wk,k,ax,iv,vx,iv,sigmax,c,xf,mx,fy,&iflg,reps,rlm,ide,chisq,cov);
	if(ier>100) {free(wk); return ier;}

/*	Calculate the SVD of matrix for fit along y-axis */
	iflg=1;
	ier=bspfit(ny,y,f,wk,k,ay,iv,vy,iv,sigmay,c,yf,my,fy,&iflg,reps,rlm,ide,chisq,cov);
	if(ier>100) {free(wk); return ier;}

/*	Set up the RHS for calculating the fits along y-axis */
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) {
			wk[j+i*la]=f[i+j*la];
			if(rlm>0.0) wk[j+ny+i*la]=0.0;
		}
	}

/*	N1 is the number of equations for fit along y-axis */
	n1=ny;
	if(rlm>0.0) n1=2*ny;
	m=my+k-2;
	for(i=0; i<nx; ++i) ier=svdevl(m,n1,ay,vy,sigmay,iv,iv,&wk[i*la],reps);

/*	Set up the RHS for calculating the fits along x-axis */
	for(j=0; j<m; ++j) {
		for(i=0; i<nx; ++i) {
			c[i+j*la]=wk[j+i*la];
			if(rlm>0.0) c[i+nx+j*la]=0.0;
		}
	}
	free(wk);
	m1=mx+k-2;
	n1=nx; if(rlm>0.0) n1=2*nx;
	for(i=0; i<m; ++i) ier=svdevl(m1,n1,ax,vx,sigmax,iv,iv,&c[i*la],reps);
 
/*	Calculate the CHI square */
	*chisq=0.0;
	nderiv=0;
	for(i=0; i<ny; ++i) {
		for(j=0; j<nx; ++j) {
			fy[j+i*la]=bspev2(mx,my,xf,yf,k,nderiv,c,la,x[j],y[i],&dfx,
					&dfy,&dfxx,&dfxy,&dfyy,&ier);
			r1=(f[j+i*la]-fy[j+i*la]);
			*chisq=(*chisq)+r1*r1;
		}
	}
	return 0;
}
