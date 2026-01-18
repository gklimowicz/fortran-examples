/*	To reduce a general real matrix to Hessenberg form using stabilised
	elementary transformations
	It is advisable to balance the matrix before applying this transformations

	A : (input/output) Array of length IA*N containing the matrix
		After execution the reduced matrix will be overwritten on
		the same array
	N : (input) Order of matrix
	IA : (input) Second dimension of array A as declared in the calling
		function
	LOW, IGH : (input) After balancing only rows and columns from
		LOW to IGH need to be considered, since other rows or
		columns contain isolated eigenvalues. If the matrix is
		not balanced use LOW=0, IGH=N-1
	INC : (output) Integer array of length N, which will contain information
		about row and column interchanges during reduction
		
	Error status is returned by the value of the function ELMHES.
		0 value implies successful execution
		113 implies that N<=1 or N>IA

	Required functions : None
*/

#include <math.h>

int elmhes(double *a, int n, int ia, int low, int igh, int inc[])

{
	int i,j,m;
	double amax,t;

	if(n<=1 || n>ia) return 113;

	if(low>igh-2) return 0;

	for(i=0; i<n; ++i) inc[i]=i;
	for(m=low+1; m<=igh-1; ++m) {
		i=m;
/*	Find the pivot */
		amax=0.0;
		for(j=m; j<=igh; ++j) {
			if(fabs(a[m-1+j*ia]) > fabs(amax)) {amax=a[m-1+j*ia]; i=j;}
		}
		inc[m]=i;

		if(i != m) {
/*	Interchange the corresponding rows and columns */
			for(j=m-1; j<n; ++j) {
				t=a[j+i*ia]; a[j+i*ia]=a[j+m*ia]; a[j+m*ia]=t;
			}
			for(j=0; j<=igh; ++j) {
				t=a[i+j*ia]; a[i+j*ia]=a[m+j*ia]; a[m+j*ia]=t;
			}
		}

		if(amax != 0.0) {
/*	Perform Gaussian elimination */
			for(i=m+1; i<=igh; ++i) {
				t=a[m-1+i*ia];
				if(t != 0.0) {
					t=t/amax;
					a[m-1+i*ia]=t;
					for(j=m; j<n; ++j) a[j+i*ia]=a[j+i*ia]-t*a[j+m*ia];
					for(j=0; j<=igh; ++j) a[m+j*ia]=a[m+j*ia]+t*a[i+j*ia];
				}
			}
		}
	}
	return 0;
}
