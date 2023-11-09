/*	To solve a linear programming problem in the standard form arising
		 in L_1 minimisation problems using the simplex method

	A : (input/output) Array of length IA*(N-M+1) containing
		the tableau of simplex algorithm
		A[I+1][0]=c_i, the cost coefficients
		Columns 1 to M contain constraints with A[0][j]=b_j and A[i+1][j]=a_i
	IA : (input) Second dimension of array A as declared in the calling
		function. IA >= M+2
	N : (input) Number of variables, each is constraint to be >=0
	M : (input) Number of constraints of form a^T X = b_i >= 0
	ID : (input/output) Integer array of length M+1 which contains
		information about interchange of variables on LHS
	IV : (input/output) Integer array of length N-M+1 which contains
		information about interchange of variables on RHS
	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
		assumed to be zero
		
	Error status is returned by the value of the function SIMPL1.
		0 value implies successful execution
		63 implies that the objective function is unbounded from below
		635 implies that the simplex algorithm failed to find
			the optimal feasible vector

	Required functions : None
*/

#include <math.h>

int simpl1(double *a, int ia, int n, int m, int id[], int iv[], double aeps)

{
	int i,j,k,k1,l,l1,jf,m1,n1,it, nit=20;
	double rmin,r1,aj;

	jf=0;
	m1=m+1; n1=n-m+1;

/*	The simplex iteration */
	for(it=1; it<=nit*(n+m); ++it) {
/*	Finding the minimum reduced cost coefficient */
		rmin=0.0; k=0;
		for(j=1; j<n1; ++j) {
			if(a[jf+j*ia]<rmin) {rmin=a[jf+j*ia]; k=j;}
			else if(iv[j]<=m-1 && 2-a[jf+j*ia]< rmin) {rmin=2-a[jf+j*ia]; k=-j;}
		}
		if(rmin>= -aeps) return 0;

/*	Finding the pivot element */
		k1=abs(k);
		rmin=0.0;
		l=0;
		for(j=1; j<=m; ++j) {
			aj=a[j+k1*ia];
			if(k<0) aj=-aj;
			if(aj > aeps) {
				r1=a[j]/aj;
				if(r1<rmin || l==0) {rmin=r1; l=j;}
			}
		}
		if(l==0) return 63;  /*	The objective function is unbounded from below */

		if(k<0) {
			a[jf+k1*ia]=-2+a[jf+k1*ia];
			iv[k1]=-iv[k1];
		}

/*	Exchange the variables */
		l1=id[l]; id[l]=iv[k1]; iv[k1]=l1;
		for(j=0; j<n1; ++j) {
			if(j != k1) {
				r1=a[l+j*ia]/a[l+k1*ia];
				for(i=0; i<m1; ++i) {
					if(i != l) a[i+j*ia]=a[i+j*ia]-a[i+k1*ia]*r1;
				}
			}
		}

		r1=fabs(a[l+k1*ia]);
		for(j=0; j<n1; ++j) {
			if(j != k1) a[l+j*ia]=a[l+j*ia]/r1;
		}
		for(i=0; i<m1; ++i) {
			if(i != l) a[i+k1*ia]= -a[i+k1*ia]/a[l+k1*ia];
		}
		a[l+k1*ia]=1./r1;
	}

	return 635;
}
