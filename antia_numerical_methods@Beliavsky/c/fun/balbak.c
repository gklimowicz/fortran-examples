/*	Perform back-transformation of a set of right eigenvectors from
	those of balanced matrix to that for original matrix
	Only real eigenvectors are handled by this function
	Since transformation is real and linear, the real and imaginary
	parts of complex eigenvectors can be transformed by two separate calls
	to this function.

	N : (input) Order of matrix
	LOW, IGH : (input) After balancing only rows and columns from
		LOW to IGH need to be considered, since other rows or
		columns contain isolated eigenvalues
	CZ : (input/output) Array of length IZ*M containing the
		eigenvectors of balanced matrix. After execution it will
		be overwritten by eigenvectors of the original matrix.
	M : (input) Number of eigenvectors to be balanced
	IZ : (input) Second dimension of array CZ as declared in the calling
		function
	D : (input) Array of length N, containing information about
		transformation used for balancing.

	Returned value is always 0

	Required functions : None
*/

#include <math.h>

int balbak(int n, int low, int igh, double *cz, int m, int iz, double d[])

{
	int i,j,k;
	double s,cs;

	for(i=low; i<=igh; ++i) {
		s=d[i];
/*	For back-transforming left eigenvectors use the following
	statement instead of the preceding one
       s=1.0/d[i];  */

	   for(j=0; j<m; ++j) cz[j+i*iz]=cz[j+i*iz]*s;
	}

	for(i=low-1; i>=0; --i) {
		k=d[i];
		if(k != i) {
/*	Exchange the corresponding rows */
			for(j=0; j<m; ++j) {
				cs=cz[j+i*iz]; cz[j+i*iz]=cz[j+k*iz]; cz[j+k*iz]=cs;
			}
		}
	}

	for(i=igh+1; i<n; ++i) {
		k=d[i];
		if(k != i) {
/*	Exchange the corresponding rows */
			for(j=0; j<m; ++j) {
				cs=cz[j+i*iz]; cz[j+i*iz]=cz[j+k*iz]; cz[j+k*iz]=cs;
			}
		}
	}
	return 0;
}
