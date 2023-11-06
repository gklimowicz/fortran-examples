/*	To calculate coefficients of Pade approximation using the known
	coefficients of Maclaurin series

	M : (input) Degree of polynomial in the numerator
	K : (input) Degree of polynomial in the denominator
	A : (output) Array of length M+K+1 containing the coefficients
		of Pade approximation. A[I-1] is the coefficient of x**I in
		the denominator, the constant term being 1. A[K+J] is
		the coefficient of x**J in the numerator
	C : (input) Array of length M+K+1 containing the coefficients
		of Maclaurin series. C[I] is the coefficient of x**I
		
	Error status is returned by the value of the function PADE.
		0 value implies successful execution
		612 implies that M<0 or K<0, no calculations are done
		Other values may be set by GAUELM

	Required functions : GAUELM
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int pade(int m, int k, double a[], double c[])

{
	int i,j,ier,iflg,k1;
	double ai,det;
	int *iwk;
	double *wk;

	if(m<0 || k<0) return 612;

	if(k>0) {
		wk=(double *) calloc((size_t) (k*k), sizeof(double));
		iwk=(int *) calloc((size_t) k, sizeof(double));
/*	Setting up the coefficients matrix of linear equations to
	calculate the coefficients in the denominator */
		for(i=0; i<k; ++i) {
			for(j=0; j<k; ++j) {
				wk[i+j*k]=0.0;
				if(m+j-i >=0) wk[i+j*k]=c[m+j-i];
			}
/*	The right-hand side vector */
			a[i]=-c[m+i+1];
		}

		iflg=0;
/*	Solve the system of linear equations */
		ier=gauelm(k,1,wk,a,&det,iwk,k,&iflg);
		free(iwk); free(wk);
		if(ier>0) return ier;
	}

/*	Calculate the coefficients in the numerator */
	a[k]=c[0];
	if(m>0) {
		for(i=1; i<=m; ++i) {
			ai=c[i];
			if(k>0) {
				k1=k-1; if(i<k) k1=i-1;
				for(j=0; j<=k1; ++j) ai=ai+c[i-j-1]*a[j];
			}
			a[k+i]=ai;
		}
	}
	return 0;
}
