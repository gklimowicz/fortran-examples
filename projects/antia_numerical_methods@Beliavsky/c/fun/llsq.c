/*	Linear least squares fit in K dimensions

	N : (input) Number of data points to be fitted
	M : (input) Number of basis functions to be used
	K : (input) Number of dimensions
	X : (input) Array of length IX*N containing the coordinates
		of points at which function values is available
		X[j][i] is the ith coordinate of jth point
		The points may have arbitrary distribution in K-space
	IX : (input) Second dimension of X in the calling function, IX >= K 
	F : (input) Array of length N containing the function values
		F[I] should be function value at X[I][0],...,X[I][K-1]
	EF : (input) Array of length N containing the estimated error in F[I]. 
	A : (output) Array of length N containing the fitted coefficients
		Note that although the number of coefficients is M, the
		rest of array is used as scratch space
	U : (output) Array of length IU*N containing the matrix U of SVD
		of the design matrix
	V : (output) Array of length IV*M containing the matrix V of SVD
		of the design matrix
	IU : (input) Second dimension of U in the calling function (IU>=M)
	IV : (input) Second dimension of V in the calling function (IV>=M)
	SIGMA : (output) Array of length M containing the singular values
		of the design matrix
	Y : (output) Array of length N containing the values of fitted
		function at each of the tabular points
	PHI : (input) Name of function to calculate the basis functions
		at any given point
	REPS : (input) Required accuracy for solution of equations using SVD
		singular values less than REPS times maximum will be set to zero
	CHISQ : (output) The value of Chi square at minimum
        COV : (output) Array of length IV*M containing the covariance
                matrix of the fitted parameters. COV(I,I) will be the
                variance in A(I).
		
	Error status is returned by the value of the function LLSQ.
		0 value implies successful execution
		606 implies that M>N, M<=0, N<=0, or K>IX
		607 implies that EF[I]<=0 for some I
		No calculations are done in both these cases
		Other values may be set by SVD

	Function PHI(M,X,Y) must be supplied by the user to calculate
	the required basis functions. M is the number of basis functions,
	X is an array of length K containing the coordinates of point
	at which the basis function needs to be calculated. Y is an
	array of length M containing the computed basis functions at X

	Required functions : SVD, SVDEVL, PHI

	THE ARGUMENTS OF THIS FUNCTION HAVE CHANGED FROM THE EARLIER VERSION.
	NOW THERE IS AN ADDITIONAL ARGUMENT COV TO CALCULATE THE COVARIANCE
	MATRIX.

*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);

int llsq(int n, int m, int k, double *x, int ix, double f[], double ef[],
	double a[], double *u, double *v, int iu, int iv, double sigma[],
	double y[], void (* phi) (int , double * , double * ), double reps,
	double *chisq, double *cov)

{
	int i,j,ik,ier;
	double s1,sigmax;
	double *wk;

	if(m>n || m<=0 || n<=0 || k>ix) return 606;
	
	wk=(double *) calloc((size_t) m, sizeof(double));
/*	Setting up the design matrix and the RHS */
	for(i=0; i<n; ++i) {
		if(ef[i]<=0.0) {free(wk); return 607;}
		a[i]=f[i]/ef[i];
		phi(m,&x[i*ix],wk);
		for(j=0; j<m; ++j) u[j+i*iu]=wk[j]/ef[i];
	}

	ier=svd(m,n,u,v,sigma,iu,iv);
	if(ier>100) {free(wk); return ier;}

/*	Calculate the least squares solution */
	ier=svdevl(m,n,u,v,sigma,iu,iv,a,reps);
 
/*	Computing the \chi^2 from fitted coefficients */
	*chisq=0.0;
	for(i=0; i<n; ++i) {
		phi(m,&x[i*ix],wk);
		s1=0.0;
		for(j=0; j<m; ++j) s1=s1+a[j]*wk[j];
		y[i]=s1;
		s1=(f[i]-y[i])/ef[i];
		*chisq=(*chisq)+s1*s1;
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


	free(wk);
	return 0;
}
