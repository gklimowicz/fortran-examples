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
