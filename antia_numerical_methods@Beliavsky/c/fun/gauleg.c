/*     To calculate weights and abscissas of a Gauss-Legendre quadrature formula

	N : (input) Number of points in the required quadrature formula
	W : (output) Array of length N, which will contain the weights
	A : (output) Array of length N containing the abscissas
		
	Error status is returned by the value of the function GAULEG.
		0 value implies successful execution
		302 implies N<=0
		321 implies that some coefficient becomes imaginary
			during calculations.
		In both these cases calculations are abandoned.
		Other values may be set by TQL2

	Required functions : GAUSRC, TQL2
*/

#include <math.h>
#include <stdlib.h>

int gausrc(int n, double w[], double a[], double *cof, double ri0);

int gauleg(int n, double w[], double a[])

{
	int i,ier;
	double ri0;
	double *wk;

	ri0=2.0;
	wk=(double *) calloc((size_t) (3*(n+1)), sizeof(double));

/*     Coefficients of recurrence relation */
	for(i=0; i<=n; ++i) {
		wk[i*3]=(2*i+1.0)/(i+1.0);
		wk[1+i*3]=0.0;
		wk[2+i*3]=i/(i+1.0);
	}
	
	ier=gausrc(n,w,a,wk,ri0);
	free(wk);
	return ier;
}
