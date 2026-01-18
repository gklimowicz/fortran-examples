/*     Solution of a system of linear equations with real symmetric
            positive definite matrix using Cholesky decomposition

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
        	equations to be solved
	A : (input/output) The matrix of coefficient of size ND*N
        	A[i][j] is the coefficient of x_j in ith equation
		at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size ND*NUM)
        	X[j][i] is the ith element of jth right hand side
		at output it will contain the solutions
	DET : (output) The determinant of the matrix
	ND : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter which determines the type of computation
		required.
        	If IFLG<=0, both elimination and solution are calculated
	    		and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
			decomposition should have been calculated earlier
		
	Error status is returned by the value of the function CHOLSK.
		0 value implies successful execution
		103 implies (N<=0 or N>ND)
		123 implies some pivot turned out to be zero

	Required functions : None
*/
 
#include <math.h>

int cholsk(int n, int num, double * a, double *x, double *det, int nd,
	int *iflg)

{
	int i,j,k,l;
	double sum;

	if(n<=0 || n>nd) return 103;

	if(*iflg < 2) {
/*     Perform elimination  */
		*det = 1.0;
		for(k=0; k<n; ++k) {
			for(i=0; i<=k-1; ++i) {
				if( a[i*nd+i] == 0.0) return 123;
				sum= a[k*nd+i];
				for(j=0; j<=i-1; ++j) sum=sum-a[i*nd+j]*a[k*nd+j];
				a[k*nd+i]=sum/a[i*nd+i];
			}

			sum = a[k*nd+k];
			for(j=0; j<=k-1; ++j) sum=sum-a[k*nd+j]*a[k*nd+j];
			if(sum<= 0.0) return 123;
			a[k*nd+k]=sqrt(sum);
			*det = (*det)*sum;
		}

		if(*iflg == 1) {*iflg = 2; return 0;}
		*iflg = 2;
	}
/*     Solution for the NUM different right-hand sides */
	for(j=0; j<num; ++j) {
/*     Forward substitution  */
		x[j*nd]= x[j*nd]/a[0];
		for(k=1; k<n; ++k) {
			for(l=0; l<=k-1; ++l) x[j*nd+k]=x[j*nd+k]-a[k*nd+l]*x[j*nd+l];
			x[j*nd+k]=x[j*nd+k]/a[k*nd+k];
		}
/*     back-substitution */
		x[j*nd+n-1]=x[j*nd+n-1]/a[(n-1)*nd+n-1];
		for(k=n-2; k>=0; --k) {
			for(l=n-1; l>=k+1; --l) x[j*nd+k]=x[j*nd+k]-x[j*nd+l]*a[l*nd+k];
			x[j*nd+k] = x[j*nd+k]/a[k*nd+k];
		}
	}
	return 0;
}

