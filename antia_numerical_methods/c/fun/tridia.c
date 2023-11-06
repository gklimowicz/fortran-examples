/*	To find specified eigenvalues and eigenvectors of a real symmetric
	tridiagonal matrix using Sturm sequence coupled with inverse iteration

	E : (input) Array of length N containing the off-diagonal elements
		of the tridiagonal matrix, E[i+1]=A(i,i+1)=A(i+1,i)
	D : (input) Array of length N containing the diagonal elements
		of the tridiagonal matrix, D[i]=A(i,i)
	N : (input) Order of the matrix
	M1 : (input) Serial number of lowest eigenvalue to be determined.
		The eigenvalues are sorted in increasing order
	M2 : (input) Serial number of highest eigenvalue to be determined.
		All eigenvalues from M1 to M2 are determined
	EI : (output) Array of length M2-M1+1 containing the calculated
		eigenvalues
	EPS1 : (input) Relative accuracy to which eigenvalues are located by bisection
		before using inverse iteration. If the inverse iteration
		does not converge to nearest eigenvalues, EPS1 can be reduced.
		A value of 0.1-0.01 times typical eigenvalue is generally sufficient.
	REPS : (input) Desired relative accuracy in eigenvalues and eigenvectors
	EV : (output) Array of length IV*N containing the
		eigenvectors. EV[i][j] should contain the ith component of
		the jth eigenvector
	IV : (input) The second dimension of array EV as declared in the
		calling function (IV>=M2-M1+1)
		
	Error status is returned by the value of the function TRIDIA.
		0 value implies successful execution
		109 implies that N<=1 or M2-M1+1>IV or M1<0 or M2>=N
		Other values may be set by TINVIT, only the last value is
		returned.

	Required functions : STURM, TINVIT, RAN1
*/

#include <math.h>
#include <stdlib.h>

int sturm(double e[], double d[], int n, int m1, int m2, double el[],
	double eu[], int *num, double reps);
int tinvit(double e[], double d[], int n, double el, double eu, double *ei,
	double ev[], double reps, int iflg, int *num);

int tridia(double e[], double d[], int n, int m1, int m2, double ei[],
	double eps1, double reps, double *ev, int iv)

{
	int i,j,i1,iflg,ier,num;
	double *wl,*wu,*v;

	if(n<=1 || m2-m1+1>iv || m1<0 || m2>=n) return 109;

	if(m1>m2) return 0;

	wl=(double *) calloc((size_t) n,sizeof(double));
	wu=(double *) calloc((size_t) n,sizeof(double));
	v=(double *) calloc((size_t) n,sizeof(double));
/*	Locate the eigenvalues */
	ier=sturm(e,d,n,m1,m2,wl,wu,&num,eps1);
	ier=0;

/*	Loop for finding individual eigenvalues and eigenvectors */
	for(i=m1; i<=m2; ++i) {
		iflg=0;
		i1=i-m1;

		if(i>m1) {
			if(fabs(wl[i]-wl[i-1])<3.*reps*fabs(wl[i])) {
/*	Set the flag for close eigenvalues */
				iflg=1;
				for(j=0; j<n; ++j) v[j]=ev[i1-1+j*iv];
			}
		}
		j=tinvit(e,d,n,wl[i],wu[i],&ei[i1],v,reps,iflg,&num);
		if(j>0) ier=j;
		for(j=0; j<n; ++j) ev[i1+j*iv]=v[j];
	}
	free(v); free(wu); free(wl);
	return ier;
}
