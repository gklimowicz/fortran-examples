/*	To solve a linear programming problem in the standard form
		 using the simplex method

	A : (input/output) Array of length IA*(N-M+1) containing
		the tableau of simplex algorithm
		A[I+1][0]=c_i, the cost coefficients
		Columns 1 to M contain constraints with A[0][j+1]=b_j and A[i+1][j+1]=a_i
		Column M+1 contains the cost coefficients for auxiliary problem
		when QF!=1.
	IA : (input) Second dimension of array A as declared in the calling
		function. IA >= M+2
	N : (input) Number of variables, each is constrained to be >=0
	M : (input) Number of constraints of form a^T X = b_i >= 0
	NV : (input) Number of variables, excluding the artificial variables
	QF : (input) Parameter to decide which objective function
		is to be minimised.
		If QF=1  the main objective function specified by
			the first column of A is to be minimised.
		Otherwise the auxiliary objective function specified by
			the last column of A is to be minimised.
	ID : (input/output) Integer array of length M+1 which contains
		information about interchange of variables on LHS
	IV : (input/output) Integer array of length N-M+1 which contains
		information about interchange of variables on RHS
	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
		assumed to be zero
		
	Error status is returned by the value of the function SIMPX.
		0 value implies successful execution
		57 implies that the objective function is unbounded from below
		531 implies that the simplex algorithm failed to find
			the optimal feasible vector

	Required functions : None
*/

#include <math.h>

int simpx(double *a, int ia, int n, int m, int nv, int qf, int id[],
	int iv[], double aeps)

{
	int i,j,k,jf,m1,n1,l,l1,it,nit=20;
	double rmin,rmax,r1;

/*	Minimise the objective function in the first column */
	if(qf==1) {jf=0; m1=m+1;}
/*	Minimise the objective function in the last column */
	else {jf=m+1; m1=jf+1;}
	n1=n-m+1;

	for(it=1; it<= nit*(n+m); ++it) {
/*	Find the minimum of the reduced cost coefficients */
		rmin=a[jf+ia]; k=1;
		for(j=2; j<n1; ++j) {
			if(a[jf+j*ia] < rmin) {rmin=a[jf+j*ia]; k=j;}
		}

		if(rmin>=0.0) {
/*	The objective function cannot be reduced further */
			if(qf==1 || a[jf]<-aeps) return 0;
/*	Check if any artificial variable is on the LHS */
			for(i=1; i<=m; ++i) {
				if(id[i]>=nv && fabs(a[i])<=aeps) {
					a[i]=0.0; l=i; k=0;
					rmax=0.0;
					for(j=1; j<n1; ++j) {
						if(fabs(a[i+j*ia])>rmax && iv[j]<=nv) {
							rmax=fabs(a[i+j*ia]);
							k=j;
						}
					}
/*	To exchange the artificial variable */
					if(k>0) goto exch;
				}
			}
			return 0;
		}

/*	Finding allowed change in the Kth variable */
		rmin=0.0; l=0;
		for(j=1; j<=m; ++j) {
			if(a[j+k*ia]>aeps) {
				r1=a[j]/a[j+k*ia];
				if(r1<rmin || l==0) {rmin=r1; l=j;}
			}
		}
		if(l==0) return 57;

/*	exchange the variables */
exch:	l1=id[l]; id[l]=iv[k]; iv[k]=l1;
		for(j=0; j<n1; ++j) {
			if(j!=k) {
				for(i=0; i<m1; ++i) {
					if(i != l) a[i+j*ia]=a[i+j*ia]-a[i+k*ia]*a[l+j*ia]/a[l+k*ia];
				}
			}
		}

		for(j=0; j<n1; ++j) {
			if(j != k) a[l+j*ia]=a[l+j*ia]/a[l+k*ia];
		}
		for(i=0; i<m1; ++i) {
			if(i != l) a[i+k*ia]= -a[i+k*ia]/a[l+k*ia];
		}
		a[l+k*ia]=1.0/a[l+k*ia];
	}

	return 531;
}
