/*	Least squares fit using B-spline basis functions in n dimensions */

#include <stdio.h>
#include <math.h>

double fun(double x, double y, double z);
double rangau(double *iseed);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int bspfit(int n, double x[], double f[], double ef[], int k, double *a,
	int la, double *v, int iv, double sigma[], double c[], double xf[],
	int no, double y[], int *iflg, double reps, double rlm, int ide,
	double *chisq, double *cov);
double bspevn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double x0[], int *ier);

int bspfitn(int n, int nk[], double *x, double *f, int k, double *a, int la,
	double *v, int iv, double *sigma, double *c, double *xf, int mk[],
	double *fy, double reps, double rlm, double ide, double *chisq);
int bspfitwn(int n, int nk[], double *x, double *f, double *ef, int k,
	double *a, int la, double *v, int iv, double sigma[], double *c,
	double *xf, int mk[], double *fy, double reps, double rlm, int ide, double *chisq);

double seed;

main()
{
	int i,i1,j,nx,ny,mx,my,ide,id,k, iflg, ier,np,nd, nk[20], mk[20];
	double hh, x[100], f[20][10][10], w[20][10][10], ax[950000],
	dfxx,dfyy,dfxy,dfx,dfy, fx, xx,yy, chi, a[900],y[200],fy[20][10][10];
	double xy[5][100],x0[9],df[9],d2f[2][2],axy[900],c[10000];
	double reps,rlm,vx[290000],sx[590],xyf[5][50];

/*	EXERCISE  10.70  */

	seed=2; np=100; id=50;
/*	No. of points is fixed by dimensions 
	No. of knots should be such that dimensions of ax and vx are
	sufficient for bspfitwn */
	nx=10;
	nd=3;
	printf("type k=order of B-spline,  m=no. of knots,\n");
	printf("    ide=order of derivative for regularisation,  rlm=regularisation parameter\n");
	scanf(" %d %d %d %le", &k,&mx,&ide,&rlm);

/*	Set up the table of values
	For simplicity assume the no. of knots to be the same along each axis */
	hh=1.0/(nx-1.0);  reps=1.e-7;
	for(i=0; i<nd; ++i) {
		nk[i]=nx;
		mk[i]=mx;
	}
	for(i=0; i<nx; ++i) {
		for(j=0; j<nd; ++j) xy[j][i]=i*hh;
	}
	hh=1./(mx-1.0);
	for(i=0; i<mx; ++i) {
		for(j=0; j<nd; ++j) xyf[j][i]=i*hh;
	}
	for(i=0; i<nx; ++i) {
		for(j=0; j<nx; ++j) {
			for(i1=0; i1<nx; ++i1) {
				f[i1][i][j]=fun(xy[0][j],xy[1][i],xy[2][i1]);
				w[i1][i][j]=1.e-3;
			}
		}
	}


/*	Use both bspfitn and bspfitwn to calculate the fit to same data */

	i=bspfitn(nd,nk,&xy[0][0],&f[0][0][0],k,ax,np,vx,id,sx,a,&xyf[0][0],
		mk,&fy[0][0][0],reps,rlm,ide,&chi);
	printf("   no. of pts = %d    no. of knots = %d   order of B-spline = %d \n", nx,mx,k);
	printf("bspfitn:      ier = %d       Chi square = %e \n",i,chi);

	i=bspfitwn(nd,nk,&xy[0][0],&f[0][0][0],&w[0][0][0],k,ax,np,vx,id,sx,c,&xyf[0][0],
		mk,&fy[0][0][0],reps,rlm,ide,&chi);
	printf("bspfitwn:     ier = %d           Chi square = %e \n", i,chi);

	for(i1=0; i1<99; ++i1) {
		printf("type x value at which function is required \n");
		printf("                       (quits when x0[0]<-100)\n");
		scanf(" %le %le %le", &x0[0],&x0[1],&x0[2]);
		if(x0[0]<-100) return 0;

		fx=bspevn(nd,mk,&xyf[0][0],id,k,a,x0,&ier);
		printf(" ier = %d  x= %e %e %e    f =  %e  \n",ier,x0[0],x0[1],x0[2],fx);

/*	Use the coefficients computed by bspfitwn */
		fx=bspevn(nd,mk,&xyf[0][0],id,k,c,x0,&ier);
		printf("using coefs from bspfitwn :  ier = %d    f =  %e  \n",ier,fx);

	}
	return;
}

double fun(double x, double y, double z)

{
	double fx;
	fx=sin(x)*sin(y)*sin(z) + 1.e-3*rangau(&seed);
	return fx;

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




/*	To calculate function value using B-spline expansion in n-dimensions

	N : (input) Number of dimensions
	NK : (input) Array of length N containing the number of knots
		along each dimension
	X : (input) Array of length NXD*N containing the knots
		X[j][i] is the ith knot along jth dimension
	NXD : (input) Second dimension of array X in calling function.
		NXD must be greater than or equal to maximum of NK[I].
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Coefficients of B-spline expansion, the coefficients
		are assumed to be stored in natural FORTRAN order with
		no gaps in data. The calling function should have
	 	dimension  WT[NK[N-1]+K-2]...[NK[1]+K-2][NK[0]+K-2]
	X0 : (input) Array of length N, containing the coordinates
		of the point at which function needs to be evaluated
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values of IER may be set by BSPLIN which is called

	BSPEVN will give the value of the B-spline expansion at X0.
	This function does not calculate the derivatives, for that
	use BSPEVN1 or BSPEVN2.

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

double bspevn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double x0[], int *ier)

{
	int i,j,j1,n1,n2,ndp,index,nderiv;
	double f,term;
	int *iwk;
	double *wk, *wk1;
 
	nderiv=0;
	iwk=(int *) calloc((size_t) (2*n),sizeof(int));
	wk=(double *) calloc((size_t) (n*nxd), sizeof(double));
	wk1=(double *) calloc((size_t) (n*nxd), sizeof(double));

	for(i=0; i<n; ++i) {
		n1=nxd*i;
		*ier=bsplin((x+i*nxd),nk[i],k,x0[i],nderiv,&wk[n1],&wk1[n1],
			&wk1[n1],&iwk[i]);
		if(*ier>100) {free(wk1); free(wk); free(iwk); return 0.0;}
		iwk[i+n]=0;
	}

/*	Calculate the sum over all points */
	f=0.0;

	do {
		index=iwk[0]+iwk[n];
		term=wk[index];
		ndp=nk[0]+k-2;

		for(i=1; i<n; ++i) {
			n1=i*nxd;
			n2=iwk[i]+iwk[i+n];
			term=term*wk[n1+n2];
			index=index+n2*ndp;
			ndp=ndp*(nk[i]+k-2);
		}
		f=f+term*wt[index];

		j=0;
		while(j <= n-1) {
			if(iwk[j+n] >= k-1) {
				iwk[j+n]=0;
				j=j+1;
			}
			else {
				iwk[j+n]=iwk[j+n]+1;
				break;
			}
		}
	} while(j<n);

	free(wk1); free(wk); free(iwk);
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


/*	To calculate linear least squares fit to B-spline basis functions in N dimension
	Weights are assumed to be unity

	N : (input) Number of dimensions
	NK : (input) Integer array of length N containing the number of
		data points along each dimension
	X : (input) Array of length LA*N containing the coordinates
		of points at which function values is available
		X[j][i] contains the ith point along jth axis
	F : (input) Array of length NK[0]NK[1]...NK[N-1] containing the
		function values. The dimension of F in the calling function
		must match the size along each dimension, F[NK[N-1]]...[NK[1]]
	K : (input) Order of B-splines required, K=4 gives cubic B-splines
	A : (output) Array of length LA*IV*N containing the matrix
		U of SVD of the design matrix for fit along each axis
	LA : (input) Second dimension of array X as declared
		in the calling function (LA >= max(NK[I]))
	V : (output) Array of length IV*IV*N containing the matrix
		V of SVD of the design matrix for fit along each axis
	IV : (input) Second dimension of A, V, XF, SIGMA in the calling function
		IV >= MAX(MK[I])+K-2
	SIGMA : (output) Array of length IV*N containing the singular
		values of the design matrix for fit along each axis
	C : (output) Array of length 2NK[0]NK[1]...NK[N-1] containing
		the fitted coefficients. Note that although the number of
		coefficients is (MK[0]+K-2)(MK[1]+K-2)...(MK[N-1]+K-2)
		the rest of array is used as scratch space
		Dimensions of array could be declared as
		C[NX][MK[N-2]+K-2]...[MK[1]+K-2][MK[0]+K-2]
		where the first dimension is increased suitably to accommodate
		the scratch space required. The first dimension NX should
		be chosen such that the total size of array exceeds the
		required value.
	XF : (input) Array of size IV*N, containing the knots
		along each axis used for defining B-spline basis functions.
		XF[j][i] should contain the ith knot along jth dimension
	MK : (input) Integer array of length N containing the number of
		knots for B-splines along each axis.
	FY : (output) Array of length NK[0]NK[1]...NK[N-1] containing
		the values of fitted function at each of the tabular points
		Dimensions of this array must match those of F.
	REPS : (input) Required accuracy for solution of equations using SVD
		singular values less than REPS times maximum will be set to zero
	RLM : (input) Parameter lambda for smoothing. If RLM<=0 no smoothing
		is applied
	IDE : (input) Order of derivative to be used for smoothing
		This is used only when RLM>0. IDE=1 for first derivative
		and IDE=2 for second derivative smoothing
	CHISQ : (output) The value of Chi square at minimum

	Error status is returned by the value of the function BSPFITN.
		0 value implies successful execution
		609 implies that RLM>0 and IDE is not acceptable
		No calculations are done in this case
		Other values of may be set by SVD or BSPFIT

	Required functions : BSPFIT, BSPLIN, BSPEVL, BSPEVN, SVD, SVDEVL
*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int bspfit(int n, double x[], double f[], double ef[], int k, double *a,
	int la, double *v, int iv, double sigma[], double c[], double xf[],
	int no, double y[], int *iflg, double reps, double rlm, int ide,
	double *chisq, double *cov);
double bspevn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double x0[], int *ier);

int bspfitn(int n, int nk[], double *x, double *f, int k, double *a, int la,
	double *v, int iv, double *sigma, double *c, double *xf, int mk[],
	double *fy, double reps, double rlm, double ide, double *chisq)

{
	int i,j,i1,i2,j1,iflg,ier,nx,ny,ju,n0,n1,m,lj,nx1,lj1,ny1,num,ni,ij,id;
	double r1;
	int *iwk;
	double *wk,*cov;

	if(rlm>0.0 && (ide<1 || ide>2)) return 609;

	nx=1;
	for(i=0; i<n-1; ++i) nx=nx*nk[i];
	n0=nx*nk[n-1];
	if(rlm>0.0) n0=2*n0;

	wk=(double *) calloc((size_t) (n0), sizeof(double));
	cov=(double *) calloc((size_t) (la*la), sizeof(double));
/*	Set the weights to one for fits along each dimension */
	for(i=0; i<la; ++i) wk[i]=1.0;
 
/*	Calculate the SVD of matrices for fit along each dimension */
	for(i=0; i<n; ++i) {
		iflg=1;
		ier=bspfit(nk[i],&x[la*i],f,wk,k,&a[i*la*iv],la,&v[i*iv*iv],iv,
			&sigma[i*iv],c,&xf[i*iv],mk[i],fy,&iflg,reps,rlm,ide,chisq,cov);
		if(ier>100) {free(wk); return ier;}
	}

	ny=1;
	ju=n-1;

	if(n-2*(n/2) == 1) {
/*	If N is odd then fit along the last dimension outside the loop */
		n1=nk[n-1]; if(rlm>0.0) n1=2*nk[n-1];
/*	Set up the RHS for fit */
		lj=n1;
		for(i=0; i<nx; ++i) {
			for(j=0; j<nk[n-1]; ++j) {
				wk[j+i*lj]=f[i+j*nx];
				if(rlm>0.0) wk[j+nk[n-1]+i*lj]=0.0;
			}
		}
		m=mk[n-1]+k-2;
		for(i=0; i<nx; ++i) {
			ier=svdevl(m,n1,&a[(n-1)*la*iv],&v[(n-1)*iv*iv],&sigma[(n-1)*iv],
					la,iv,&wk[i*lj],reps);
		}
 
/*	setup the RHS for fit along the N-2 th dimension */
		nx1=nx/nk[n-2];
		ny=m;
		lj1=nk[n-2];
		if(rlm>0.0) lj1=2*lj1;
		for(i1=0; i1<nk[n-2]; ++i1) {
			for(i=0; i<nx1; ++i) {
				for(j=0; j<ny; ++j) {
					c[i1+i*lj1+j*nx1*lj1]=wk[j+i*lj+i1*nx1*lj];
					if(rlm>0.0) c[i1+nk[n-2]+i*lj1+j*nx1*lj1]=0.0;
				}
			}
		}
		nx=nx1; ju=n-2;
	}

	else {
/*	setup the RHS for fit along the N-1 th dimension */
		lj=nk[n-1];
		if(rlm>0.0) lj=2*nk[n-1];
		for(i=0; i<nx; ++i) {
			for(j=0; j<nk[n-1]; ++j) {
				c[j+i*lj]=f[i+j*nx];
				if(rlm>0.0) c[j+nk[n-1]+i*lj]=0.0;
			}
		}
	}
 
/*	Loop for fit along each dimension, each pass fits along 2 dimensions */
	for(j1=ju; j1>=0; j1 -= 2) {
		lj=nk[j1]; if(rlm>0.0) lj=2*nk[j1];
		n1=lj;
		m=mk[j1]+k-2;
		for(i=0; i<nx*ny; ++i) {
			ier=svdevl(m,n1,&a[j1*la*iv],&v[j1*iv*iv],&sigma[j1*iv],la,iv,&c[i*lj],reps);
		}

/*	Set up the RHS for fit along next dimension */
		nx1=nx/nk[j1-1];
		ny1=ny*m;
		lj1=nk[j1-1]; if(rlm>0.0) lj1=2*nk[j1-1];
		for(i1=0; i1<ny; ++i1) {
			for(i2=0; i2<m; ++i2) {
				for(i=0; i<nk[j1-1]; ++i) {
					for(j=0; j<nx1; ++j) {
						wk[i+j*lj1+i2*lj1*nx1+i1*nx1*lj1*m]=c[i2+j*lj+i*nx1*lj+i1*nx*lj];
						if(rlm>0.0) wk[i+nk[j1-1]+j*lj1+i2*lj1*nx1+i1*nx1*lj1*m]=0.0;
					}
				}
			}
		}
		nx=nx1; ny=ny1;
		lj=nk[j1-1];
		if(rlm>0.0) lj=2*nk[j1-1];
		m=mk[j1-1]+k-2;
		n1=lj;
		for(i=0; i<nx*ny; ++i) {
			ier=svdevl(m,n1,&a[(j1-1)*la*iv],&v[(j1-1)*iv*iv],&sigma[(j1-1)*iv],
					la,iv,&wk[i*lj],reps);
		}

		if(j1==1) {
/*	Store the fitted coefficients in array C */
			n1=nk[0]; if(rlm>0.0) n1=2*nk[0];
			for(i=0; i<m; ++i) {
				for(j=0; j<ny; ++j) c[i+j*m]=wk[i+j*lj];
			}
		}

		else {
/*	Set up the RHS for fit along next dimension */
			lj1=nk[j1-2]; if(rlm>0.0) lj1=2*nk[j1-2];
			m=mk[j1-1]+k-2;
			nx1=nx/nk[j1-2];
			ny1=ny*m;
			for(i1=0; i1<m; ++i1) {
				for(i2=0; i2<nk[j1-2]; ++i2) {
					for(i=0; i<nx1; ++i) {
						for(j=0; j<ny; ++j) {
							c[i2+i*lj1+i1*lj1*nx1+j*lj1*nx1*m]=wk[i1+i*lj+i2*lj*nx1+j*lj*nx];
							if(rlm>0.0) c[i2+nk[j1-2]+i*lj1+i1*lj1*nx1+j*lj1*nx1*m]=0.0;
						}
					}
				}
			}
			nx=nx1; ny=ny1;
		}
	}
 
/*	Calculate the Chi Square */
	*chisq=0.0;
	iwk=(int *) calloc((size_t) n,sizeof(int));
	for(i=0; i<n; ++i) iwk[i]=0;
 
/*	Loop over all points */
	do {
		ij=iwk[0];
		wk[0]=x[ij];
		id=nk[0];
		for(i=1; i<n; ++i) {
			ij=ij+id*iwk[i];
			id=id*nk[i];
			wk[i]=x[iwk[i]+i*la];
		}

		fy[ij]=bspevn(n,mk,xf,iv,k,c,wk,&ier);
		r1=f[ij]-fy[ij];
		*chisq=(*chisq)+r1*r1;

/*	Choose the next point */
		j=0;
		while(j <= n-1) {
			if(iwk[j]>=nk[j]-1) {
				iwk[j]=0;
				j=j+1;
			}
			else {
				iwk[j]=iwk[j]+1;
				break;
			}
		}
	} while(j<n);

	free(iwk); free(wk);
	return 0;
}



/*	To calculate linear least squares fit to B-spline basis functions in N dimension
	version for general weights, would require long time to calculate

	N : (input) Number of dimensions
	NK : (input) Integer array of length N containing the number of
		data points along each dimension
	X : (input) Array of length LA*N containing the coordinates
		of points at which function values is available
		X[j][i] contains the ith point along jth axis
	F : (input) Array of length NK[0]NK[1]...NK[N-1] containing the
		function values. The dimension of F in the calling function
		must match the size along each dimension, F[NK[N-1]]...[NK[0]]
	EF : (input) Array of length NK[0]NK[1]...NK[N-1] containing the
		estimated errors in F
	K : (input) Order of B-splines required, K=4 gives cubic B-splines
	A : (output) Array of length 
		NK[0]NK[1]...NK[N-1] * (MK[0]+K-2)(MK[1]+K-2)...(MK[N-1]+K-2)
		containing the matrix U of SVD of the design matrix
		If RLM>0 the size should be (N+1) times this value
	LA : (input) Second dimension of arrays X as declared
		in the calling function LA >= max(NK[I])
	V : (output) Array of length 
		[(MK[0]+K-2)(MK[1]+K-2)...(MK[N-1]+K-2)]**2
		containing the matrix V of SVD of the design matrix
	IV : (input) Second dimension of XF in the calling function
		IV >= max(MK[I])
	SIGMA : (output) Array of length 
		(MK[0]+K-2)(MK[1]+K-2)...(MK[N-1]+K-2)
		containing the singular values of the design matrix
	C : (output) Array of length NK[0]NK[1]...NK[N-1] containing
		the fitted coefficients. If RLM>0 then the required size
		will be N+1 times this value. Note that although the number of
		coefficients is (MK[0]+K-2)(MK[1]+K-2)...(MK[N-1]+K-2)
		the rest of array is used as scratch space
		Dimensions of array could be declared as
		C[NX][MK[N-2]+K-2]...[MK[1]+K-2][MK[0]+K-2]
		where the first dimension is increased suitably to accommodate
		the scratch space required. The first dimension NX should
		be chosen such that the total size of array exceeds the
		required value.
	XF : (input) Array of size IV*N, containing the knots
		along each axis used for defining B-spline basis functions.
		XF[j][i] should contain the ith knot along jth dimension
	MK : (input) Integer array of length N containing the number of
		knots for B-splines along each axis.
	FY : (output) Array of length NK[0]NK[1]...NK[N-1] containing
		the values of fitted function at each of the tabular points
		Dimensions of this array must match those of F.
	REPS : (input) Required accuracy for solution of equations using SVD
		singular values less than REPS times maximum will be set to zero
	RLM : (input) Parameter lambda for smoothing. If RLM<=0 no smoothing
		is applied
	IDE : (input) Order of derivative to be used for smoothing
		This is used only when RLM>0. IDE=1 for first derivative
		and IDE=2 for second derivative smoothing
	CHISQ : (output) The value of Chi square at minimum

	Error status is returned by the value of the function BSPFITWN.
		0 value implies successful execution
		608 implies that LA or IV are not large enough
		609 implies that RLM>0 and IDE is not acceptable
		No calculations are done in these cases
		Other values may be set by SVD or BSPLIN

	Required functions : BSPLIN, BSPEVN, SVD, SVDEVL
*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspevn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double x0[], int *ier);

int bspfitwn(int n, int nk[], double *x, double *f, double *ef, int k,
	double *a, int la, double *v, int iv, double sigma[], double *c,
	double *xf, int mk[], double *fy, double reps, double rlm, int ide, double *chisq)

{
	int i,j,m,n0,n1,n2,ndb,nderiv,ij,id,ji,j1,ier,left;
	double r1,xb,wk1;
	int *iwk;
	double *wk, *wk2;

	if(rlm>0.0 && (ide<1 || ide>2)) return 609;

/*     N1 is the number of equations to be solved
     M is the number of coefficients to be determined */
	 n1=1; m=1;
	 for(i=0; i<n; ++i) {
		 if(nk[i]>la || mk[i]+2*k+2 > iv) return 608;
		 n1=n1*nk[i];
		 m=m*(mk[i]+k-2);
	 }

	 n2=n1;
	 if(rlm>0.0) n1=(n+1)*n1;
	 iwk=(int *) calloc((size_t) (3*n), sizeof(int));

	 ndb=la+4; nderiv=0;
	 if(rlm>0.0) {
		 nderiv=ide;
		 for(i=0; i<n; ++i) iwk[i]=ndb*3*i+ide*ndb;
	 }

	wk=(double *) calloc((size_t) (3*ndb*n), sizeof(double));
/*     Set up the matrix for equations */
	for(i=0; i<n; ++i) {
		xb=x[i*la];
		ier=bsplin(&xf[i*iv],mk[i],k,xb,nderiv,&wk[3*ndb*i],&wk[3*ndb*i+ndb],
				&wk[3*ndb*i+2*ndb],&left);
		if(ier>100) {free(wk); free(iwk); return ier;}
		iwk[i+n]=0;
	}

	wk2=(double *) calloc((size_t) (n+2), sizeof(double));
/*	Loop over each tabular point */
	do {
		ij=iwk[n];
		id=nk[0];
		for (i=1; i<n; ++i) {
			ij=ij+iwk[i+n]*id;
			id=id*nk[i];
		}
		c[ij]=f[ij]/ef[ij];

/*     set each coefficient in the equation */
		for(i=0; i<n; ++i) iwk[i+2*n]=0;
 
/*	Loop over the coefficients */
		do {
			ji=iwk[2*n];
			for(i=0; i<=n; ++i) wk2[i]=wk[iwk[2*n]];
			if(rlm>0.0) wk2[1]=wk[iwk[0]+iwk[2*n]];

			id=mk[0]+k-2;
			for(i=1; i<n; ++i) {
				ji=ji+id*iwk[i+2*n];
				id=id*(mk[i]+k-2);
				wk1=wk[3*i*ndb+iwk[i+2*n]];
				wk2[0]=wk2[0]*wk1;
				if(rlm>0.0) {
					for(j=0; j<n; ++j) {
						if(i != j) wk2[j+1]=wk2[j+1]*wk1;
						else wk2[j+1]=wk2[j+1]*wk[iwk[i]+iwk[i+2*n]];
					}
				}
			}

/*	setup the coefficient */
			a[ij*m+ji]=wk2[0]/ef[ij];
			if(rlm>0.0) {
				for(i=1; i<=n; ++i) a[ji+m*(ij+i*n2)]=rlm*wk2[i];
			}

/*	Go to the next coefficient */
			j1=0;
			while(j1 <= n-1) {
				if(iwk[2*n+j1]>= (mk[j1]+k-3)) {
					iwk[2*n+j1]=0;
					j1=j1+1;
				}
				else {
					iwk[2*n+j1]=iwk[2*n+j1]+1;
					break;
				}
			}
		} while(j1<n);

/*	If all coefficients are done go to the next tabular point */
		j=0;
		while(j <= n-1) {
			if(iwk[j+n]>=nk[j]-1) {
/*	If Jth dimension is exhausted go to the next */
				iwk[j+n]=0;
				xb=x[iwk[j+n]+j*la];
				ier=bsplin(&xf[j*iv],mk[j],k,xb,nderiv,&wk[j*3*ndb],
						&wk[j*3*ndb+ndb],&wk[j*3*ndb+2*ndb],&left);
				if(ier>100) {free(wk2); free(wk); free(iwk); return ier;}
				j=j+1;
			}
			else {
				iwk[j+n]=iwk[j+n]+1;
				xb=x[iwk[j+n]+j*la];
				ier=bsplin(&xf[j*iv],mk[j],k,xb,nderiv,&wk[j*3*ndb],
						&wk[j*3*ndb+ndb],&wk[j*3*ndb+2*ndb],&left);
				if(ier>100) {free(wk2); free(wk); free(iwk); return ier;}
				break;
			}
		}
	} while(j<n);

 
/*	If all points are done find the SVD of the matrix */
	ier=svd(m,n1,a,v,sigma,m,m);
	if(ier>100) {free(wk2); free(wk); free(iwk); return ier;}

/*     Setup the remaining part of RHS and solve the equations */
	for(i=n2; i<n1; ++i) c[i]=0.0;
	ier=svdevl(m,n1,a,v,sigma,m,m,c,reps);
 
/*     Calculate the \chi^2 */
	*chisq=0.0;
	for(i=0; i<n; ++i) iwk[i]=0;
 
/*	loop over all points */
	do {
		ij=iwk[0];
		id=nk[0];
		wk[0]=x[ij];
		for(i=1; i<n; ++i) {
			ij=ij+id*iwk[i];
			id=id*nk[i];
			wk[i]=x[iwk[i]+i*la];
		}
		fy[ij]=bspevn(n,mk,xf,iv,k,c,wk,&ier);
		r1=(f[ij]-fy[ij])/ef[ij];
		*chisq=(*chisq)+r1*r1;
 
/*	Go to the next point */
		j=0;
		while(j<=n-1) {
			if(iwk[j] >= nk[j]-1) {
				iwk[j]=0;
				j=j+1;
			}
			else {
				iwk[j]=iwk[j]+1;
				break;
			}
		}
	} while(j<n);

	free(wk2); free(wk); free(iwk);
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

