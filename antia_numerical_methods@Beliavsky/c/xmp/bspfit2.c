/*	Least squares fit using B-spline basis functions in 2 dimensions */

#include <stdio.h>
#include <math.h>

double fun(double x, double y);
double rangau(double *seed);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspev2(int nx, int ny, double *x, double *y, int k, int nderiv,
	double *wt, int iw, double x0, double y0, double *dfx, double *dfy,
	double *dfxx, double *dfyy, double *dfxy, int *ier);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier);
int bspfit(int n, double x[], double f[], double ef[], int k, double *a,
	int la, double *v, int iv, double sigma[], double c[], double xf[],
	int no, double y[], int *iflg, double reps, double rlm, int ide,
	double *chisq, double *cov);
double bspevn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double x0[], int *ier);

int bspfit2(int nx, int ny, double x[], double y[], double *f, int k,
	double *ax, double *ay, int la, double *vx, double *vy, int iv,
	double sigmax[], double sigmay[], double *c, double xf[], double yf[],
	int mx, int my, double *fy, double reps, double rlm, int ide, double *chisq);
int bspfitw2(int nx, int ny, double x[], double y[], double *f, double *ef,
	int k, double *a, int la, double *v, int iv, double sigma[], double *c,
	int ic, double xf[], double yf[], int mx, int my, double *fy,
	double reps, double rlm, int ide, double *chisq);

double seed;

main()
{
	int i,i1,j,nx,ny,mx,my,ide,id,k, ia,iv, ier,np,nmax;
	double hh, x[100], f[100][100], w[100][100], ax[650000],ay[150000],
	dfxx,dfyy,dfxy,dfx,dfy, fx, xx,yy, chi, a[100][100],y[200],fy[100][100];
	double reps,rlm,vx[190000],vy[5000],sx[390],sy[50],xf[150],yf[50],c[100][100];

	seed=2; np=100; id=60;

/*	For simplicity assume that the no. of points and knots is the same
	along both directions */

	printf("type n=no. of points,   k=order of B-spline,  m=no. of knots,\n");
	printf("    ide=order of derivative for regularisation,  rlm=regularisation parameter\n");
	scanf(" %d %d %d %d %le", &nx,&k,&mx,&ide,&rlm);

/*	Set up the table of function values */
	hh=1.0/(nx-1.0); ny=nx; reps=1.e-7;
	for(i=0; i<nx; ++i) {
		x[i]=i*hh;
		y[i]=i*hh;
	}

	my=mx;
	hh=1./(mx-1.0);
	for(i=0; i<mx; ++i) {
		xf[i]=i*hh;
		yf[i]=i*hh;
	}
	for(i=0; i<ny; ++i) {
		for(j=0; j<nx; ++j) {
			f[i][j]=fun(x[j],y[i]);
			w[i][j]=1.e-5;
		}
	}

/*	Use both bspfit2 and bspfitw2 to calculate the fit to same data */

	i=bspfit2(nx,ny,x,y,&f[0][0],k,ax,ay,np,vx,vy,id,sx,sy,&a[0][0],xf,
		yf,mx,my,&fy[0][0],reps,rlm,ide,&chi);
	printf(" no. of points = %d %d, no. of knots = %d %d, order of B-spline = %d \n", nx,ny,mx,my,k);
	printf("  Regularisation parameter = %e    order of derivative = %d\n", rlm,ide);
	printf("bspfit2 : ier = %d   chisq = %e \n", i,chi);

/*	dimension ia should be > (mx+k-2)**2  */
	ia=200;
	i=bspfitw2(nx,ny,x,y,&f[0][0],&w[0][0],k,ax,ia,vx,ia,sx,&c[0][0],np,xf,
		yf,mx,my,&fy[0][0],reps,rlm,ide,&chi);
	printf("bspfitw2 : ier = %d   chisq = %e \n", i,chi);


	ide=2;
	for(i1=0; i1<99; ++i1) {
		printf("type xx,yy = point at which function is to be evaluated\n");
		printf("                            (quits when xx<-100)\n");
		scanf(" %le %le", &xx,&yy);
		if(xx<-100) return 0;

		fx=bspev2(mx,my,xf,yf,k,ide,&a[0][0],np,xx,yy,&dfx,&dfy,
				&dfxx,&dfxy,&dfyy,&ier);
		printf(" x = %e %e   f =  %e   f' = %e  %e \n",xx,yy,fx,dfx,dfy);
		printf(" f'' = %e %e %e \n",dfxx,dfxy,dfyy);

		fx=bspev2(mx,my,xf,yf,k,ide,&c[0][0],np,xx,yy,&dfx,&dfy,
				&dfxx,&dfxy,&dfyy,&ier);
		printf(" using coef from bspfitw2 :    f =  %e   f' = %e  %e \n",fx,dfx,dfy);
		printf(" f'' = %e %e %e \n",dfxx,dfxy,dfyy);

	}
	return;
}

double fun(double x, double y)

{
	double fx;
	fx=sin(x)*sin(y) +1.e-5*rangau(&seed); 
	return fx;

}



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


/*	To calculate function value using B-spline expansion

	N : (input) Number of knots to define B-splines
	X : (input) Array of length N containing the knots.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	NDERIV : (input) Number of derivatives required
		For NDERIV<=0 only function value is calculated
		For NDERIV=1 first derivative is also calculated
		For NDERIV>1 both first and second derivatives are calculated
	WT : (input) Coefficients of B-spline expansion
	X0 : (input) The point at which expansion has to be evaluated
	DF : (output) First derivative of function at X0
	DDF : (output) Second derivative of function at X0
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values of IER may be set by BSPLIN which is called

	BSPEVL = SUM_{i=1}^{N+K-2} WT(I) \phi_i(X0)
	where \phi_i(x) are B-spline basis functions on knots X

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier)

{
	int i,left;
	double f;
	double *wk, *wk1, *wk2;

	wk=(double *) calloc((size_t) (n+k), sizeof(double));
	wk1=(double *) calloc((size_t) (n+k), sizeof(double));
	wk2=(double *) calloc((size_t) (n+k), sizeof(double));
	*ier=bsplin(x,n,k,x0,nderiv,wk,wk1,wk2,&left);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	f=0.0; *df=0.0; *ddf=0.0;
	for(i=left; i<=left+k-1; ++i) {
		f = f + wt[i]*wk[i];
		*df = *df + wt[i]*wk1[i];
		*ddf = *ddf + wt[i]*wk2[i];
	}
	free(wk2); free(wk1); free(wk);
	return f;
}


/*	To calculate linear least squares fit to B-spline basis functions in 1 dimension

	N : (input) Number of data points to be fitted
	X : (input) Array of length N containing the coordinates
		of points at which function values is available
	F : (input) Array of length N containing the function values
		F[I] should be function value at X[I]
	EF : (input) Array of length N containing the estimated error in F[I]. 
	K : (input) Order of B-splines required, K=4 gives cubic B-splines
	A : (output) Array of length LA*2N containing the matrix
		U of SVD of the design matrix
	LA : (input) Second dimension of A in the calling function (LA>=NO+K-2)
	V : (output) Array of length IV*(NO+K-2) containing the matrix
		V of SVD of the design matrix
	IV : (input) Second dimension of V, COV in the calling function (IV>=NO+K-2)
	SIGMA : (output) Array of length NO+K-2 containing the singular
		values of the design matrix
	C : (output) Array of length 2N containing the fitted coefficients
		Note that although the number of coefficients is NO+K-2, the
		rest of array is used as scratch space
	XF : (input) Array of size NO, containing
		the knots used for defining B-spline basis functions.
		The knots must be distinct and in ascending order.
	NO : (input) Number of knots for B-splines, the number of basis
		functions would be NO+K-2
	Y : (output) Array of length N containing the values of fitted
		function at each of the tabular points
	IFLG : (input/output) Integer specifying the type of calculation required
		IFLG=0 The matrix will be calculated and solved for coefficients
			the fitted values Y and CHISQ are also calculated
		IFLG=1 The matrix will be calculated and SVD
			is obtained, but coefficients are not calculated
		IFLG=2 The SVD of matrix is assumed to be available in
			arrays A, V, SIGMA and coefficients C are calculated
		IFLG=3 The SVD of matrix is assumed to be available in arrays
			A, V, SIGMA and coefficients C are calculated and in
			addition fitted values Y and CHISQ are also calculated
	REPS : (input) Required accuracy for solution of equations using SVD
		singular values less than REPS times maximum will be set to zero
	RLM : (input) Parameter lambda for smoothing. If RLM<=0 no smoothing
		is applied
	IDE : (input) Order of derivative to be used for smoothing
		This is used only when RLM>0. IDE=1 for first derivative
		and IDE=2 for second derivative smoothing
	CHISQ : (output) The value of Chi square at minimum
	COV : (output) Array of length IV*M containing the covariance
                matrix of the fitted parameters. COV(I,I) will be the
                variance in C(I).

	Error status is returned by the value of the function BSPFIT.
		0 value implies successful execution
		608 implies that NO+K-2>N or K<2
		609 implies that RLM>0 and IDE is not acceptable
		610 implies that EF[I]<=0 for some I
		No calculations are done in all these cases
		Other values may be set by SVD or BSPLIN

	Required functions : BSPLIN, BSPEVL, SVD, SVDEVL

	THE ARGUMENTS OF THIS FUNCTION HAVE CHANGED FROM THE EARLIER VERSION.
	NOW THERE IS AN ADDITIONAL ARGUMENT COV TO CALCULATE THE COVARIANCE
	MATRIX.

*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier);

int bspfit(int n, double x[], double f[], double ef[], int k, double *a,
	int la, double *v, int iv, double sigma[], double c[], double xf[],
	int no, double y[], int *iflg, double reps, double rlm, int ide,
	double *chisq, double *cov)

{
	int i,j,ik,n1,m,ndb,nb,nderiv,left,ier;
	double xb,r1,df,ddf,sigmax;
	double *wk;

	if(n<no+k-2 || k<2) return 608;

	if(rlm>0.0 && (ide<1 || ide>2)) return 609;

/*	N1 is the number of equations to be solved */
	n1=n;
	if(rlm>0.0) n1=2*n;

	if(*iflg<2) {
/*	Set up the matrix equation and obtain SVD
	M is the number of coefficients to be determined */
		m=no+k-2;
		ndb=m+1;
		nderiv=0;
		if(rlm>0.0) {
			nderiv=ide;
			nb=ndb; if(ide == 2) nb=2*ndb;
		}

		wk=(double *) calloc((size_t) (3*ndb), sizeof(double));
/*	Set up the matrix for equations */
		for(i=0; i<n; ++i) {
			xb=x[i];
			if(ef[i]<= 0.0) {free(wk); return 610;}
			ier=bsplin(xf,no,k,xb,nderiv,wk,&wk[ndb],&wk[ndb*2],&left);
			if(ier>100) {free(wk); return ier;}
			for(j=0; j<m; ++j) {
				a[j+i*la]=wk[j]/ef[i];
				if(rlm>0) a[j+(i+n)*la]=rlm*wk[nb+j];
			}
		}
		free(wk);
		ier=svd(m,n1,a,v,sigma,la,iv);
		if(ier>100) return ier;

		if(*iflg==1) {*iflg=2; return 0;}

	}

/*	Setup the RHS and solve the equations */
	for(i=0; i<n; ++i) {
		c[i]=f[i]/ef[i];
		if(rlm>0.0) c[i+n]=0.0;
	}

	ier=svdevl(m,n1,a,v,sigma,la,iv,c,reps);
	if(*iflg==2) return ier;
	*iflg=2;

/*	Calculate the \chi^2 */
	*chisq=0.0;
	nderiv=0;
	for(i=0; i<n; ++i) {
		y[i]=bspevl(no,xf,k,nderiv,c,x[i],&df,&ddf,&ier);
		r1=(f[i]-y[i])/ef[i];
		*chisq=(*chisq)+r1*r1;
	}
	 
/*	Computing the covariance matrix for fitted coefficients */
	sigmax=0.0;
	for(i=0; i<m; ++i)
		if(sigma[i]>sigmax) sigmax=sigma[i];
	for(i=0; i<m; ++i) {
		for(j=0; j<=i; ++j) {
			cov[j+i*iv]=0.0;
			for(ik=0; ik<m; ++ik) if(sigma[ik]>reps*sigmax)
			cov[j+i*iv]=cov[j+i*iv]+v[ik+j*iv]*v[ik+i*iv]/(sigma[ik]*sigma[ik]);
			cov[i+j*iv]=cov[j+i*iv];
		}
	}


	return 0;
}


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


/*	To calculate the B-spline basis functions at a specified point

	X : (input) Array of length NX containing the knots.
		The knots must be distinct and in ascending order.
	NX : (input) Number of knots
	K : (input) Order of B-spline, 0< K, K=4 gives cubic B-splines
	XB : (input) The point at which B-spline basis functions are to be evaluated
	NDERIV : (input) Number of derivatives required
		NDERIV<=0 only B-splines are calculated
		NDERIV=1 first derivative is also calculated
		NDERIV>1 first and second derivatives are also calculated
	B : (output) Array of length NX+K-2 containing the value of
		B-spline basis functions
	DB : (output) Array of length NX+K-2 containing the value of
		the first derivative of B-spline basis functions (if NDERIV>0)
	DDB : (output) Array of length NX+K-2 containing the value of
		the second derivative of B-spline basis functions (if NDERIV>1)
	LEFT : (output) XB is located between X[LEFT] and X[LEFT+1]
		
	Error status is returned by the value of the function BSPLIN.
		0 value implies successful execution
		26 implies XB > X[NX-1]
		27 implies XB < X[0]
		203 implies NX<2, K<1

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left)

{
	int i,j,igh, mid,nigh,lx,ier;
	double t1,t2,t3,p1,p2;
	double *wk, *dr, *dl;
	static int low = -1;

	if(nx <= 1 || k<1 ) return 203;
	ier=0;

/*	If the previous value of LOW is inadmissible, set the range to (0,N-1) */
	if(low<0 || low>=nx-1) {low=0; igh=nx-1;}
	else igh=low+1;

	while((xb<x[low] && xb<x[igh]) || (xb>x[low] && xb>x[igh])) {
/*	Extend the range */
		if( xb>x[low] ) {
/*	Extend the range on higher side */
			if(igh >= nx-1) {ier=26; low=nx-2; break;}
			else {
				nigh=igh+2*(igh-low); if(nx-1 < nigh) nigh=nx-1;
				low=igh; igh=nigh;
			}
		}

		else {
/*	Extend the range on lower side */
			if(low <= 0) {ier=27; igh=low+1; break;}
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

 
/*	Evaluate the B-spline basis functions

	Define the extra knots on either side of table
	Note that the function assumes knots from -K+2 to NX+K-1
	and the B-splines B_{i,k}, i ranges from 0 to NX+K-3 
	The knots are stored in scratch array wk. */

	wk=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));
	dr=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));
	dl=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));

	for(i=0; i<nx; ++i) wk[i+k]=x[i];
	for(i=1; i<=k; ++i) {
		wk[k-i]=x[0];
		wk[nx+i-1+k]=x[nx-1];
	}

	for(i=0; i<nx+k-2; ++i) {b[i]=0.0; db[i]=0.0; ddb[i]=0.0;}
	*left=low;
	lx=low-1;
	b[lx+1]=1;

/*	The recurrence relation for B-splines */
	for(j=1; j<=k-1; ++j) {
		dr[j] = wk[low+j+k] - xb;
		dl[j] = xb - wk[low+1-j+k];
		t1=0.0;
		for(i=1; i<=j; ++i) {
			t2=b[lx+i]/(dr[i]+dl[j+1-i]);
			b[lx+i]=t1+t2*dr[i];
			t1=t2*dl[j+1-i];
		}
		b[lx+j+1]=t1;
			
/*	Calculate the first derivative using recurrence relations */
		if(j == k-2 && nderiv > 0) {
			t1=0.0;
			for(i=1; i<=j+1; ++i) {
				t2=b[lx+i]/(wk[low+i+k]-wk[low+i+1]);
				db[lx+i]=(k-1)*(t1-t2);
				t1=t2;
			}
			db[lx+j+2]=(k-1)*t1;
		}
 
/*	Calculate the second derivative using recurrence relations */
		if(j == k-3 && nderiv>1) {
			t2=0.0; p1=0.0;
			for(i=1; i<=j+1; ++i) {
				t3=b[lx+i]/(wk[low+i+k]-wk[low+i+2]);
				p2=(t2-t3)/(wk[low+i+k]-wk[low+i+1]);
				ddb[lx+i]=(k-2)*(k-1)*(p1-p2);
				t2=t3; p1=p2;
			}
			p2=t2/(wk[low+j+2+k]-wk[low+j+3]);
			ddb[lx+j+2]=(k-2)*(k-1)*(p1-p2);
			ddb[lx+j+3]=(k-2)*(k-1)*p2;
		}
	}

/*	For K=2 the first derivative has to be calculated outside the loop */
	if(k == 2 && nderiv > 0) {
		t2=1./(wk[low+1+k]-wk[low+2]);
		db[lx+1]=-t2;
		db[lx+2]=t2;
	}

/*	For K=3 the second derivative has to be calculated outside the loop */
	if(k == 3 && nderiv > 1) {
		t3=1./(wk[low+1+k]-wk[low+3]);
		p2=-t3/(wk[low+1+k]-wk[low+2]);
		ddb[lx+1]=-2.*p2;
		p1=p2;
		p2=t3/(wk[low+2+k]-wk[low+3]);
		ddb[lx+2]=2.*(p1-p2);
		ddb[lx+3]=2.*p2;
	}
	free(dl); free(dr); free(wk);
	return ier;
}



/*     Solution of a system of linear equations using Gaussian elimination
     	for a band matrix

	N : (input) Number of equations to be solved
	KB : (input) Bandwidth of matrix A[i][j]=0 if abs(I-J)>KB
	NUM : (input) Number of different sets (each with N equations) of
		equations to be solved
	A : (input/output) The matrix of coefficients of size LJ*(3*KB+1)
		A[J-I+KB][I] is the coefficient of x_J in Ith equation
		at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
		X[j][i] is the ith element of jth right hand side
		at output it will contain the solutions
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer variable used as a flag to specify the type
		of computation required
     	If IFLG=-1, both elimination and solution are calculated
     		without pivoting and IFLG is set to 2
	If IFLG=0, both elimination and solution are computed
     		with partial pivoting and IFLG is set to 2
     	If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
     	If IFLG>=2 only solution is calculated, the triangular
     		decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUBND.
		0 value implies successful execution
		104 implies N<=0 or N>LJ or KB>N
		124 implies some pivot turned out to be zero and hence
	     		matrix must be nearly singular

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg)

{
	int i,j,k,km,l,kb1,m1,m2;
	double r1, t1;
	double *wk;

	if(n<=0 || n>lj || kb>n) return 104;

	kb1=kb+1;
	if(*iflg < 2) {
/*     Perform elimination */
		for(i=0; i<n ; ++i) {
			for(j=2*kb+1; j<=3*kb; ++j) a[j*lj+i] = 0.0;
		}

		*det=1.0; *idet=0;
		wk=(double *) calloc((size_t) (3*kb+1), sizeof(double));
		for(k=0; k<n-1; ++k) {
/*     Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			if(*iflg >= 0) {
				m1=k+kb; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) {
					if(fabs(a[l+(k-l+kb)*lj]) > r1) {
						r1=fabs(a[l+(k-l+kb)*lj]);
						km=l;
					}
				}
			}

			inc[k]=km;
			if(km != k) {
/*     Interchange the rows if needed  */
				m1=2*kb+k; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) wk[l-k]=a[k+(l-k+kb)*lj];
				for(l=k; l<=m1; ++l) {
					a[k+(l-k+kb)*lj]=a[km+(l-km+kb)*lj];
					a[km+(l-km+kb)*lj]=wk[l-k];
				}
				*det= -(*det);
			}

			*det = (*det)*a[kb*lj+k];
			if( a[kb*lj+k] == 0.0) {free(wk); return 124;}
			if(*det != 0.0) {
/*     Scale the value of the determinant   */
				while(fabs(*det) > 32.0) {
					*det = (*det)*0.03125e0; *idet = *idet + 5;
				}

				while(fabs(*det) < 0.03125e0) {
					*det = (*det)*32.0; *idet = *idet - 5;
				}
			}

			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<= m1; ++l) {
				a[l+(k-l+kb)*lj] = a[l+(k-l+kb)*lj]/a[kb*lj+k];
				m2=k+2*kb; if(n-1 < m2) m2=n-1;
				for(i=k+1; i<=m2; ++i)
					a[l+(i-l+kb)*lj]=a[l+(i-l+kb)*lj]-a[l+(k-l+kb)*lj]*a[k+(i-k+kb)*lj];
			}
		}

		free(wk);
		*det = (*det)*a[(n-1)+kb*lj];
		inc[n-1]=n-1;
		if(a[(n-1)+kb*lj] == 0.0) return 124;

		if(*iflg==1) {*iflg=2; return 0;}
		*iflg=2;
	}
		
/*     Solution for the NUM different right-hand sides */
	for(j=0; j<num; ++j) {
/*     Forward substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k]) {
				t1= x[j*lj+k];
				x[j*lj+k] = x[j*lj+inc[k]];
				x[j*lj+inc[k]] = t1;
			}
			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<=m1; ++l)
				x[j*lj+l]=x[j*lj+l]-a[l+(k-l+kb)*lj]*x[j*lj+k];
		}

/*     back-substitution  */
		x[j*lj+n-1] = x[j*lj+n-1]/a[(n-1)+kb*lj];
		for(k=n-2; k>=0; --k) {
			m1=k+2*kb; if(n-1 < m1) m1=n-1;
			for(l=m1; l>=k+1; --l)
				x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k+(l-k+kb)*lj];
			x[j*lj+k] = x[j*lj+k]/a[kb*lj+k];
		}
	}
	return 0;
}



/*	To calculate the Singular Value Decomposition of a matrix A=U D Vtranspose

	N : (input) Number of variables
	M : (input) Number of equations
	A : (input/output) Matrix of coefficients of size LA*M
		After execution it will contain the matrix U
	V : (output) The matrix V of size LV*N
	SIGMA : (output) Array of length N, containing the singular values
	LA : (input) Actual value of second dimension of A in the calling function
	LV : (input) Actual value of second dimension of V in the calling function
		
	Error status is returned by the value of the function SVD.
		0 value implies successful execution
		12 QR iteration failed to converge to required accuracy
		105 implies N<=0, N>LV, M<=0, N>LA, N>M

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv)

{
	int i,j,k,l,itr,ier, itmax=30;
	double f, g, h, rmax, s, r1, r2, c, x, y, z, aeps, eps=1.e-16;
	double *e;

	if(n>m || n<=0 || m<=0 || n>la || n>lv) return 105;
	ier=0;

/*	Reduction to Bidiagonal form using Householder transformations */
	g=0.0; rmax=0.0;
	e=(double *) calloc((size_t) n, sizeof(double));

	for(i=0; i<n; ++i) {
/*	Off-diagonal elements of bidiagonal form  */
		e[i]=g;
		s=0.0;
		for(j=i; j<m; ++j) s=s+a[i+j*la]*a[i+j*la];
		if(s <= 0.0) {
/*	transformation not required */
			g=0.0;
		}
		else {
			f= a[i+i*la];
			g=sqrt(s);
			if(f>=0.0) g=-g;
			h=f*g-s;
			a[i+i*la] = f-g;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]= a[j+k*la]+f*a[i+k*la];
			}
		}

/*	Diagonal elements of bidiagonal form  */
		sigma[i]=g;
		s=0.0;
		for(j=i+1; j<n; ++j) s=s+a[j+i*la]*a[j+i*la];

		if(s<= 0.0) g=0.0;
		else {
			f= a[i*la+(i+1)];
			g=sqrt(s);
			if(f>= 0.0) g=-g;
			h=f*g-s;
			a[i*la+(i+1)]=f-g;
			for(j=i+1; j<n; ++j) e[j]=a[j+i*la]/h;

			for(j=i+1; j<m; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+j*la]*a[k+i*la];
				for(k=i+1; k<n; ++k) a[k+j*la] = a[k+j*la]+s*e[k];
			}
		}
		r1=fabs(sigma[i])+fabs(e[i]);
		if(r1 > rmax) rmax=r1;
	}

/*	Accumulation of right hand transformation in array V */
	for(i=n-1; i>=0; --i) {
		if(g != 0.0) {
			h=g*a[i*la+(i+1)];
			for(j=i+1; j<n; ++j) v[i+j*lv]=a[j+i*la]/h;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+i*la]*v[j+k*lv];
				for(k=i+1; k<n; ++k) v[j+k*lv]=v[j+k*lv]+s*v[i+k*lv];
			}
		}

		for(j=i+1; j<n; ++j) {
			v[j+i*lv]=0.0; v[i+j*lv]=0.0;
		}
		v[i+i*lv]=1;
		g= e[i];
	}

/*	Accumulation of left hand transformation overwritten on matrix A */
	for(i=n-1; i>=0; --i) {
		g=sigma[i];
		for(j=i+1; j<n; ++j) a[j+i*la]=0.0;
		if(g != 0.0) {
			h=g*a[i+i*la];

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]=a[j+k*la]+f*a[i+k*la];
			}

			for(j=i; j<m; ++j) a[i+j*la]=a[i+j*la]/g;
		}
		else {
			for(j=i; j<m; ++j) a[i+j*la]=0.0;
		}
		a[i+i*la] = a[i+i*la]+1;
	}

/*	Diagonalisation of the bidiagonal form */
	aeps=eps*rmax;
/*	Loop over the singular values */
	for(k=n-1; k>=0; --k) {
/*	The QR transformation */
		for(itr=1; itr<=itmax; ++itr) {

/*	Test for splitting */
			for(l=k; l>=0; --l) {
				if(fabs(e[l]) < aeps) goto split;
				if(fabs(sigma[l-1]) < aeps) break;
			}

/*	cancellation of E[L] if L>1  */
			c=0.0; s=1.0;
			for(i=l; i<=k; ++i) {
				f=s*e[i];
				e[i] = c*e[i];
				if(fabs(f) < aeps) goto split;
				g=sigma[i];
				sigma[i]=sqrt(f*f+g*g);
				c=g/sigma[i];
				s=-f/sigma[i];

				for(j=0; j<m; ++j) {
					r1= a[j*la+(l-1)];
					r2= a[i+j*la];
					a[j*la+(l-1)]=r1*c+r2*s;
					a[i+j*la]=c*r2-s*r1;
				}
			}

split:			z=sigma[k];
			if(l == k) {
/*	QR iteration has converged */
				if(z < 0.0) {
					sigma[k] = -z;
					for(j=0; j<n; ++j) v[k+j*lv]=-v[k+j*lv];
				}
				break;
			}

			if(itr==itmax) {ier=12; break;}
					
/*	calculating shift from bottom 2x2 minor */
			x=sigma[l];
			y=sigma[k-1];
			g=e[k-1];
			h=e[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g=sqrt(1.+f*f);
			if(f < 0.0) g=-g;
			f=((x-z)*(x+z)+h*(y/(f+g)-h))/x;

/*	next QR transformation */
			c=1.0; s=1.0;
/*	Given's rotation  */
			for(i=l+1; i<=k; ++i) {
				g=e[i];
				y=sigma[i];
				h=s*g;
				g=c*g;
				e[i-1]=sqrt(f*f+h*h);
				c=f/e[i-1];
				s=h/e[i-1];
				f=c*x+s*g;
				g=c*g-s*x;
				h=s*y;
				y=c*y;

				for(j=0; j<n; ++j) {
					x=v[j*lv+(i-1)];
					z=v[i+j*lv];
					v[j*lv+(i-1)]=c*x+s*z;
					v[i+j*lv]=c*z-s*x;
				}

				sigma[i-1]=sqrt(f*f+h*h);
				if(sigma[i-1] != 0.0) {
					c=f/sigma[i-1];
					s=h/sigma[i-1];
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for(j=0; j<m; ++j) {
					y= a[j*la+(i-1)];
					z= a[i+j*la];
					a[j*la+(i-1)] = c*y+s*z;
					a[i+j*la] = c*z-s*y;
				}
			}

			e[l]=0.0;
			e[k]=f;
			sigma[k]=x;
		}
	}
	free(e);
	return ier;
}



/*	To evaluate the solution of a system of linear equations using SVD

	N : (input) Number of variables
	M : (input) Number of equations
	U : (input) Array of size LU*M containing the left-hand transformation
	V : (input) Array of size LV*N containing the right-hand transformation
	SIGMA : (input) Array of size N containing the singular values
	LU : (input) Second dimension of array U in the calling function
	LV : (input) Second dimension of array V in the calling function
	B : (input/output) Array of length M containing the RHS
		after execution it will contain the solution
	REPS : (input) Relative accuracy. All singular values < REPS*(Max of singular values)
		will be reduced to zero

	The returned value is always zero

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps)

{
	int i,j;
	double smax, aeps, s;
	double *wk;

/*	Finding the largest singular value */
	smax=0.0;
	for(i=0; i<n; ++i)
		if(sigma[i] > smax) smax=sigma[i];

	aeps=smax*reps;
	wk=(double *)calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i) {
		s=0.0;
/*	Only SIGMA[I] > AEPS contribute to the solution */
		if(sigma[i] > aeps) {
			for(j=0; j<m; ++j) s=s+b[j]*u[i+j*lu];
			s=s/sigma[i];
		}
		wk[i]=s;
	}

	for(i=0; i<n; ++i) {
		s=0.0;
		for(j=0; j<n; ++j) s=s+v[j+i*lv]*wk[j];
		b[i]=s;
	}
	free(wk);
	return 0;
}


/*	To generate random numbers with Gaussian probability distribution
	It generates random numbers with zero mean and variance of 1.
	
	SEED : (input/output) real seed, it should be positive and
		less than AN. It is updated by the function and should
		not be modified between two calls, unless a fresh
		sequence is required

	Required functions : None
	
	THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
	VERSION AS THE SEED IS NOW DOUBLE INSTEAD OF INT.

*/

#include <math.h>


double rangau(double *seed)

{
	int n2;
	double am=2147483648.0, a=45875.0, ac=453816693.0, an=2147483647.0, r1, rn;

	rn=a*(*seed)+ac; n2=rn/am; r1=rn-am*n2;
	if(*seed==0.0) *seed=0.1;
	rn=sqrt(2.0*log(an/(*seed)))*cos(2.0*M_PI*r1/an);
	*seed=r1;
	return rn;
}

