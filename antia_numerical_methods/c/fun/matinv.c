/*	To calculate inverse of a square matrix

	N : (input) Order of matrix
	IA : (input) The second dimension of arrays A and AI as specified
		in the calling function
	A : (input) Array of length IA*N containing the matrix
	AI : (output) Array of length IA*N which will contain the
		calculated inverse of A
		
	Error status is returned by the value of the function MATINV.
		0 value implies successful execution
		Nonzero values may be set by function GAUELM

	It is possible to use CROUT instead of GAUELM for calculating
	the triangular decomposition.

	Required functions : GAUELM (or CROUT)
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double *a, double *x, double *det,
       int inc[], int lj, int * iflg);
int crout(int n, int num, double *a, double *x, double *det,
       double *idet, int inc[], int lj, int * iflg);

int matinv(int n, int ia, double *a, double *ai)

{
	int i,j,iflg,ier,idet;
	double det,r;
	int *iwk;

	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) ai[j+i*ia]=0.0;
		ai[i+i*ia]=1.0;
	}

	iwk=(int *) calloc((size_t) n, sizeof(int));
	iflg=0;
	ier=gauelm(n,n,a,ai,&det,iwk,ia,&iflg);
/*	ier=crout(n,n,a,ai,&det,&idet,iwk,ia,&iflg);  */
	free(iwk);
/*	transpose the matrix to get the correct inverse in normal form */
	for(i=0; i<n; ++i) {
		for(j=0; j<i; ++j) {
			r=ai[j+i*ia];
			ai[j+i*ia]=ai[i+j*ia];
			ai[i+j*ia]=r;
		}
	}
	return ier;
}
