/*	Least squares polynomial fit using orthogonal polynomials in n dimensions
	Weights are assumed to be equal for all points and points are
	assumed to be on a hyper-rectangular mesh

	N : (input) Number of dimensions
	NK : (input) Integer array of length N containing the number of
		data points along each direction
	X : (input) Array of length LA*N containing the coordinates
		of tabular points, X[j][i] contains the ith point along
		jth dimension
	F : (input) Array of length NK[0]*NK[1]*...*NK[N-1] containing
		the function values. The dimension of F in the calling function
		must match the size along each dimension, F[NK[N-1]]...[NK[0]]
	AX : (output) Array of length LA*(3*N+3) containing information about
		fit along each direction. AX[3*J][I], AX[3*J+1][I], AX[3*J+2][I] will
		respectively contain the coefficients, alpha, beta, gamma
		for fit along Jth dimension.
		The rest of the array is used as scratch space
	LA : (input) Second dimension of arrays X and AX as declared
		in the calling function. LA >= MAX(NK[I])
	C : (output) Array of length (MK[0]+1)(MK[1]+1)...(MK[N-1]+1) containing
		the fitted coefficients of product of orthogonal polynomials 
		The dimension of F in the calling function must match the size
		along each dimension, C[MK[N-1]+1]...[MK[1]+1][MK[0]+1]
	MK : (input) Integer array of length N containing the required
		degree of polynomial in each dimension
	FY : (output) Array of same size and shape as F containing the
		fitted function values at each tabular point	
	CHISQ : (output) the Chi square value for the fit
		
	Error status is returned by the value of the function POLFITN.
		0 value implies successful execution
		605 implies that LA < MAX(NK[I]) 
		In  this case no calculations are done
		other values may be set by POLFIT1
	
	The fitted polynomial can be calculated at any value of x using POLEVN

	Required functions : POLFIT1, POLEVN, POLORT
*/

#include <math.h>
#include <stdlib.h>

int polfit1(int n, int m, int num, double x[], double *f, double w[],
	double *a, double alp[], double beta[], double gam[]);
int polevn(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f);

int polfitn(int n, int nk[], double *x, double *f, double *ax, int la,
	double *c, int mk[], double *fy, double *chisq)

{
	int i,j,j1,i1,i2,lj,nmax,nx,ny,num,m,n1,ni,nx1,ny1,lj1,ij,id,ier;
	double r1;
	int *iwk;
	double *wk,*wk1;

	n1=3*n*la;
/*     Set the weights to one for fits along each dimension */
	for(i=0; i<la; ++i) ax[i+n1]=1.0;

	nx=1; nmax=nk[n-1];
	for(i=0; i<n-1; ++i) {
		nx=nx*nk[i];
		if(nk[i]>nmax) nmax=nk[i];
	}
	if(nmax>la) return 605;

	ny=1;
	wk=(double *) calloc((size_t) (nx*nk[n-1]),sizeof(double));
	wk1=(double *) calloc((size_t) (nx*nk[n-1]),sizeof(double));
 
/*     Set up the RHS for fit along Nth dimension */
	lj=nk[n-1];
	for(i=0; i<nx; ++i) {
		for(j=0; j<nk[n-1]; ++j) wk[j+i*lj]=f[i+j*nx];
	}

/*     Loop for fit along each dimension */
	for(j1=n-1; j1>=0; --j1) {
		num=nx*ny;
		lj=nk[j1];
		m=mk[j1]+1;
		ni=j1*3*la;
		ier=polfit1(lj,mk[j1],num,&x[j1*la],wk,&ax[n1],wk1,&ax[ni],&ax[ni+la],&ax[ni+2*la]);
		if(ier>0) {free(wk1); free(wk); return ier;}

		if(j1>0) {
/*     Set up the RHS for fit along next dimension */
			nx1=nx/nk[j1-1];
			ny1=ny*m;
			lj1=nk[j1-1];
			for(i1=0; i1<ny; ++i1) {
				for(i2=0; i2<m; ++i2) {
					for(i=0; i<nk[j1-1]; ++i) {
						for(j=0; j<nx1; ++j) wk[i+j*lj1+i2*lj1*nx1+i1*nx1*lj1*m] =
							wk1[i2+j*m+i*nx1*m+i1*nx*m];
					}
				}
			}

			nx=nx1; ny=ny1;
		}

		else {
/*     Store the fitted coefficients in array C */
			for(i=0; i<m; ++i) {
				for(j=0; j<ny; ++j) c[i+j*m]=wk1[i+j*m];
			}
		}
	}

/*     Calculate the Chi Square */
	*chisq=0.0;
	free(wk1);
	iwk=(int *) calloc((size_t) n, sizeof(int));

	for(i=0; i<n; ++i) iwk[i]=0;
 
/*     Loop over all points */
	do {
		ij=iwk[0];
		wk[0]=x[ij];
		id=nk[0];
		for(i=1; i<n; ++i) {
			ij=ij+id*iwk[i];
			id=id*nk[i];
			wk[i]=x[iwk[i]+i*la];
		}

		ier=polevn(n,mk,ax,la,c,wk,&fy[ij]);
		r1=f[ij]-fy[ij];
		*chisq=(*chisq)+r1*r1;

/*     Choose the next point */
		j=0;
		while(j<=n-1) {
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
