/*	To perform back-transformation on eigenvectors of reduced tridiagonal matrix
	to obtain eigenvectors of original real symmetric matrix reduced by TRED2
	This function is not required if eigenvectors are computed using TQL2

	A : (input/output) Array of length IA*N containing the transformation
		matrix. The last column of this array is used as scratch
		space and will be destroyed during execution
	IA : (input) The second dimension of array A as declared in the
		calling function (IA>=N)
	N : (input) Order of the matrix
	Z : (input/output) Array of length IZ*N containing the
		eigenvectors of tridiagonal matrix. After execution it will
		be overwritten by eigenvectors of original matrix
	IZ : (input) Second dimension of array Z as declared in the calling
		function (IZ>=NZ)
	NZ : (input) Number of eigenvectors to be transformed

	The returned value is always zero.

	Required functions : None
*/

#include <math.h>

int trbak(double *a, int ia, int n, double *z, int iz, int nz)

{
	int i,j,k;
	double s;

/*	Loop on eigenvectors */
	for(i=0; i<nz; ++i) {
/*	The matrix multiplication */
		for(j=0; j<n; ++j) {
			s=0.0;
			for(k=0; k<n-1; ++k) s=s+a[k+j*ia]*z[i+k*iz];
/*	To take care of the last column of A, which is overwritten */
			if(j==n-1) s=s+z[i+(n-1)*iz];
			a[n-1+j*ia]=s;
		}

		for(j=0; j<n; ++j) z[i+j*iz]=a[n-1+j*ia];
	}
	return 0;
}
