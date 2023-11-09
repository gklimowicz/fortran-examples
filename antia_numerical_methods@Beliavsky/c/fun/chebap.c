/*	To calculate coefficients of rational function Chebyshev approximation
	using the known coefficients of Chebyshev expansion

	M : (input) Degree of polynomial in the numerator
	K : (input) Degree of polynomial in the denominator
	A : (output) Array of length M+K+1 containing the coefficients
		of rational function approximation. A[I-1] is the coefficient
		of T_I(x) in the denominator, the constant term being 1.
		A[K+J] is the coefficient of T_J(x) in the numerator
	C : (input) Array of length M+2K+1 containing the coefficients
		of Chebyshev expansion. C[I] is the coefficient of T_I(x).
		Coefficient of T_0(x) is C[0]/2
		
	Error status is returned by the value of the function CHEBAP.
		0 value implies successful execution
		612 implies that M<0 or K<0, no calculations are done
		Other values may be set by GAUELM

	Required functions : GAUELM
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int chebap(int m, int k, double a[], double c[])

{
	int i,j,iflg,ier;
	double ai,det;
	int *iwk;
	double *wk;

	if(m<0 || k<0) return 612;

	if(k>0) {
/*	Setting up the coefficients of matrix of linear equations
	to calculate the coefficients in the denominator */
	wk=(double *) calloc((size_t) (k*k), sizeof(double));
	iwk=(int *) calloc((size_t) k, sizeof(int));
	for(i=0; i<k; ++i) {
		for(j=0; j<k; ++j) wk[i+j*k]=0.5*(c[m+2+i+j]+c[abs(m+j-i)]);

/*	The right-hand side vector for linear equations */
		a[i]=-c[m+i+1];
	}

	iflg=0;
/*	Solve the system of linear equations */
	ier=gauelm(k,1,wk,a,&det,iwk,k,&iflg);
	free(iwk); free(wk);
	if(ier>0) return ier;
	}

/*	Calculating the coefficients in the numerator */
	for(i=0; i<=m; ++i) {
		ai=c[i];
		for(j=0; j<k; ++j) ai=ai+0.5*a[j]*(c[i+j+1]+c[abs(i-1-j)]);
		a[k+i]=ai;
	}
	a[k]=0.5*a[k];
	return 0;
}
