/*	To find eigenvalues and eigenvectors of a complex Hermitian matrix

	AR : (input) Array of length IA*N containing the real part of the matrix
	AI : (input) Array of length IA*N containing the imaginary part of the matrix
	N : (input) Order of the matrix
	IA : (input) Second dimension of arrays AR, AI as declared in the calling function
	EI : (output) Array of length N containing the eigenvalues which
		should be real.
	VR : (output) Array of length IV*N containing the real part of eigenvectors
		VR[i][j] is the ith component of jth eigenvector
	VI : (output) Array of length IV*N containing the imaginary part of
		eigenvectors VI[i][j] is the ith component of jth eigenvector
	IV : (input) Second dimension of array VR, VI as declared in the calling function
	REPS : (input) Required tolerance, should be of the order of machine
		accuracy
		
	Error status is returned by the value of the function HEREVP.
		0 value implies successful execution
		111 implies that N<=1 or N>IA or N>IZ
			in which case no calculations are done
		Other values may be set by TRED2 or TQL2

	Required functions : TRED2, TQL2
*/

#include <math.h>
#include <stdlib.h>

int tred2(double *a, int n, int ia, double d[], double e[]);
int tql2(double *z, int n, int iz, double d[], double e[], double reps);

int herevp(double *ar, double *ai, int n, int ia, double ei[], double *vr,
	double *vi, int iv, double reps)

{
	int i,j,i1,n2,ier;
	double *wk,*d,*e;

	if(n<=1 || n>ia || n>iv) return 111;

	wk=(double *) calloc((size_t) (4*n*n), sizeof(double));
	d=(double *) calloc((size_t) (2*n), sizeof(double));
	e=(double *) calloc((size_t) (2*n), sizeof(double));
	n2=2*n;

/*	Setup the 2N*2N real matrix */
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) {
			wk[j+i*n2]=ar[j+i*ia];
			wk[j+n+(i+n)*n2]=ar[j+i*ia];
			wk[j+(i+n)*n2]=ai[j+i*ia];
			wk[j+n+i*n2]=-ai[j+i*ia];
		}
	}

/*	To reduce the 2N*2N matrix to tridiagonal form */
	ier=tred2(wk,n2,n2,d,e);
	if(ier>100) {free(e); free(d); free(wk); return ier;}

/*	Find eigenvalues and eigenvectors of tridiagonal matrix */
	ier=tql2(wk,n2,n2,d,e,reps);
	if(ier>100) {free(e); free(d); free(wk); return ier;}
 
/*	Since all eigenvalues are repeated and sorted in ascending order
	pick alternate eigenvalues and eigenvectors */
	for(i=0; i<n; ++i) {
		i1=2*i;
		ei[i]=d[i1];
		for(j=0; j<n; ++j) {
			vr[i+j*iv]=wk[i1+j*n2];
			vi[i+j*iv]=wk[i1+(j+n)*n2];
		}
	}
	free(e); free(d); free(wk);
	return 0;
}
