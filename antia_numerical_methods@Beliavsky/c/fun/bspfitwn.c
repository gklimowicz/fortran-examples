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
