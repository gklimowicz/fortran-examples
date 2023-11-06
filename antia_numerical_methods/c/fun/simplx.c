/*	To solve a linear programming problem using the simplex method

	A : (input/output) Array of length IA*(N+M2+1) containing
		the tableau of simplex algorithm
		A[I+1][0]=c_i, the cost coefficients
		columns 1 to M1 contain constraint of first (see below) type
		columns M1+1 to M1+M2 contain constraint of second type
		columns M1+M2+1 to M1+M2+M3 contain constraint of third type
		For all constraints A[0][j+1]=b_j and A[i+1][j+1]=a_i
	IA : (input) Second dimension of array A as declared in the calling
		function. IA >= M1+M2+M3+2
	N : (input) Number of variables, each is constrained to be >=0
	M1 : (input) Number of constraints of form a^T X <= b_i >= 0
	M2 : (input) Number of constraints of form a^T X >= b_i >= 0
	M3 : (input) Number of constraints of form a^T X = b_i >= 0
	X : (output) Array of length N, containing the optimal feasible
		vector
	F : (output) Optimum value of the objective function
	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
		assumed to be zero
		
	Error status is returned by the value of the function SIMPLX.
		0 value implies successful execution
		57 implies that the objective function is unbounded from below
		58 implies that there is no basic feasible vector
		505 implies that N,M1,M2 or M3 is negative or 
			IA < M1+M2+M3+2, no calculations are done
		506 implies that some coefficient b_i for constraints
			is negative, no calculations are done
		531 implies that the simplex algorithm failed to find
			the optimal feasible vector
		532 implies that a basic feasible vector to start the
			calculations could not be found

	Required functions : SIMPX
*/

#include <math.h>
#include <stdlib.h>

int simpx(double *a, int ia, int n, int m, int nv, int qf, int id[],
	int iv[], double aeps);

int simplx(double *a, int ia, int n, int m1, int m2, int m3, double x[],
	double *f, double aeps)

{
	int i,j,k,m,n1,nm,nv,is,qf,ier;
	double s;
	int *iwk;

	if(m1<0 || m2<0 || m3<0 || n<0 || ia<(m1+m2+m3+2)) return 505;

/*	Testing for negative b_i */
	m=m1+m2+m3;
	for(i=1; i<=m; ++i) {
		if(a[i]<0.0) return 506;
	}

	a[0]=0.0;
	if(m2>0) {
/*	Add extra artificial variables */
		for(j=n+1; j<=n+m2; ++j) {
			for(i=0; i<=m; ++i) a[i+j*ia]=0.0;
			a[m1+j-n+j*ia]=-1;
		}
	}
	n1=n+m+m2;
	nm=n+m2+1;
	nv=n+m1+m2;
	iwk=(int *) calloc((size_t) (n1+1),sizeof(int));
	for(i=0; i<nm-1; ++i) iwk[i+m+1]=i;
	for(i=0; i<m; ++i) iwk[i+1]=nm-1+i;

	if(m2+m3>0) {
/*	Cost coefficients for the auxiliary objective function */
		for(i=0; i<nm; ++i) {
			s=0.0;
			for(j=m1+1; j<=m; ++j) s=s+a[j+i*ia];
			a[m+1+i*ia]=-s;
		}
		qf=0;

/*	solve the auxiliary problem to get a basic feasible vector */
		ier=simpx(a,ia,n1,m,nv,qf,iwk,&iwk[m],aeps);
		if(ier>0) {free(iwk); return ier;}
/*	If the minimum is positive, then no basic feasible vector exists */
		if(a[m+1]<-aeps) {free(iwk); return 58;}

/*	Remove the artificial variables from the tableau */
		for(i=1; i<nm; ++i) {
			if(iwk[i+m]>=nv && i<=n1-m) {
				is=1;
				for(j=i+m+1; j<n1; ++j) {
					if(iwk[j]<nv) break;
					is=is+1;
				}
				n1=n1-is;
				if(n1-m>=i) {
					for(j=i; j<=n1-m; ++j) {
						iwk[j+m]=iwk[j+m+is];
						for(k=0; k<=m; ++k) a[k+j*ia]=a[k+j*ia+is*ia];
					}
				}
			}
		}

/*	If artificial variables are not removed, then quit */
		if(n1 != nv) {free(iwk); return 532;}

	}

	qf=1;
	ier=simpx(a,ia,n1,m,nv,qf,iwk,&iwk[m],aeps);
/*	The minimum value of the objective function */
	*f=-a[0];
	for(i=1; i<=m; ++i) {
		if(iwk[i]<=n-1) x[iwk[i]]=a[i];
	}
	for(i=m+1; i<n1; ++i) {
		if(iwk[i]<=n-1) x[iwk[i]]=0.0;
	}
	free(iwk);
	return 0;
}
