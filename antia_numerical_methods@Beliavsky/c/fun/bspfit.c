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
