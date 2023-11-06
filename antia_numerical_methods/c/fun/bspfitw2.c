/*	To calculate linear least squares fit to B-spline basis functions in 2 dimension
	version for general weights, but much slower than BSPFIT2

	NX : (input) Number of data points along x-axis
	NY : (input) Number of data points along y-axis
	X,Y : (input) Array of length NX,NY containing the coordinates
		of points at which function values is available
	F : (input) Array of length IC*NY containing the function values
		F[j][i] should be function value at X[i], Y[j]
	EF : (input) Array of length IC*NY containing the error estimates in F[J][I]
	K : (input) Order of B-splines required, K=4 gives cubic B-splines
	A : (output) Array of length LA*3*NX*NY containing
		the matrix U of SVD of the design matrix
	LA : (input) Second dimension of array A as declared
		in the calling function  LA >= (MX+K-2)*(MY+K-2)
	V : (output) Array of length IV*(MX+K-2)*(MY+K-2) containing
		the matrix V of SVD of the design matrix
	IV : (input) Second dimension of V in the calling function
		IV >= (MX+K-2)*(MY+K-2)
	SIGMA : (output) Array of length (MX+K-2)*(MY+K-2)
		containing the singular values of the design matrix
	C : (output) Array of length IC*(MY+K-2) containing the fitted coefficients
	IC : (input) Second dimension of arrays C, F, EF, FY as declared
		in the calling function (IC >= NX)
	XF : (input) Array of size MX, containing the knots
		along x-axis used for defining B-spline basis functions.
	YF : (input) Array of size MY, containing the knots
		along y-axis used for defining B-spline basis functions.
	MX : (input) Number of knots for B-splines along x-axis,
		the number of basis functions would be MX+K-2
	MY : (input) Number of knots for B-splines along y-axis,
		the number of basis functions would be MY+K-2
	FY : (output) Array of length IC*NY containing the values of fitted
		function at each of the tabular points
	REPS : (input) Required accuracy for solution of equations using SVD
		singular values less than REPS times maximum will be set to zero
	RLM : (input) Parameter lambda for smoothing. If RLM<=0 no smoothing
		is applied
	IDE : (input) Order of derivative to be used for smoothing
		This is used only when RLM>0. IDE=1 for first derivative
		and IDE=2 for second derivative smoothing
	CHISQ : (output) The value of Chi square at minimum

	Error status is returned by the value of the function BSPFITW2.
		0 value implies successful execution
		608 implies that MX+K-2>NX, MY+K-2>NY or K<2
		609 implies that RLM>0 and IDE is not acceptable
		610 implies that EF[I][J]<=0 for some I,J
		No calculations are done in all these cases
		Other values may be set by SVD or BSPLIN

	Required functions : BSPLIN, BSPEV2, SVD, SVDEVL
*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
double bspev2(int nx, int ny, double *x, double *y, int k, int nderiv,
	double *wt, int iw, double x0, double y0, double *dfx, double *dfy,
	double *dfxx, double *dfyy, double *dfxy, int *ier);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier);

int bspfitw2(int nx, int ny, double x[], double y[], double *f, double *ef,
	int k, double *a, int la, double *v, int iv, double sigma[], double *c,
	int ic, double xf[], double yf[], int mx, int my, double *fy,
	double reps, double rlm, int ide, double *chisq)

{
	int i,j,n1,m1,m2,m,ndb,nderiv,nb,left,ier,k1,k2,ij,ij1;
	double xb,yb,r1,dfx,dfy,dfxx,dfxy,dfyy;
	double *wk, *wk1;

	if(nx<k || ny<k || k<2) return 608;

	if(rlm>0.0 && (ide<1 || ide>2)) return 609;

/*     N1 is the number of equations to be solved */
	n1=nx*ny;
	if(rlm>0.0) n1=3*n1;
 
/*     Set up the matrix equation and obtain SVD 
       M is the number of coefficients to be determined */
	m1=mx+k-2; m2=my+k-2; m=m1*m2;
	ndb=mx+4; if(my>mx) ndb=my+4;
	nderiv=0;
	if(rlm>0) {
		nderiv=ide;
		nb=ndb; if(ide==2) nb=2*ndb;
	}
 
	wk=(double *) calloc((size_t) (3*ndb), sizeof(double));
	wk1=(double *) calloc((size_t) (3*ndb), sizeof(double));
/*     Set up the matrix for equations */
	for(i=0; i<nx; ++i) {
		xb=x[i];
		ier=bsplin(xf,mx,k,xb,nderiv,wk,&wk[ndb],&wk[2*ndb],&left);
		if(ier>100) {free(wk1); free(wk); return ier;}
		for(j=0; j<ny; ++j) {
			yb=y[j];
			ier=bsplin(yf,my,k,yb,nderiv,wk1,&wk1[ndb],&wk1[2*ndb],&left);
			if(ier>100) {free(wk1); free(wk); return ier;}
			if(ef[i+j*ic]<= 0.0) {free(wk1); free(wk); return 610;}
			ij=i+j*nx;
			for(k1=0; k1<m1; ++k1) {
				for(k2=0; k2<m2; ++k2) {
					ij1=k1+k2*m1;
					a[ij1+ij*la]=wk[k1]*wk1[k2]/ef[i+j*ic];
					if(rlm>0.0) {
						a[ij1+(ij+nx*ny)*la]=rlm*wk[nb+k1]*wk1[k2];
						a[ij1+(ij+2*nx*ny)*la]=rlm*wk[k1]*wk1[nb+k2];
					}
				}
			}
		}
	}
	free(wk1); free(wk);

	ier=svd(m,n1,a,v,sigma,la,iv);
	if(ier>100)  return ier;

	wk=(double *) calloc((size_t) n1, sizeof(double)); 
/*     Setup the RHS and solve the equations */
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) {
			ij=i+j*nx;
			wk[ij]=f[i+j*ic]/ef[i+j*ic];
			if(rlm>0.0) {
				wk[ij+nx*ny]=0.0;
				wk[ij+2*nx*ny]=0.0;
			}
		}
	}

	ier=svdevl(m,n1,a,v,sigma,la,iv,wk,reps);

	for(i=0; i<m1; ++i) {
		for(j=0; j<m2; ++j) c[i+j*ic]=wk[i+m1*j];
	}
	free(wk);
 
/*     Calculate the \chi^2 */
	*chisq=0.0;
	nderiv=0;
	for(i=0; i<ny; ++i) {
		for(j=0; j<nx; ++j) {
			fy[j+i*ic]=bspev2(mx,my,xf,yf,k,nderiv,c,ic,x[j],y[i],&dfx,
					&dfy,&dfxx,&dfxy,&dfyy,&ier);
			r1=(f[j+i*ic]-fy[j+i*ic])/ef[j+i*ic];
			*chisq=(*chisq)+r1*r1;
		}
	}
	return 0;
}
