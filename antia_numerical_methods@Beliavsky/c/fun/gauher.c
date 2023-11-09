/*     To calculate weights and abscissas of a Gauss-Hermite quadrature formula

	N : (input) Number of points in the required quadrature formula
	W : (output) Array of length N, which will contain the weights
	A : (output) Array of length N containing the abscissas
		
	Error status is returned by the value of the function GAUHER.
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

int gauher(int n, double w[], double a[])

{
	int i,ier;
	double ri0, pis=1.772453850905516;
	double *wk;

	ri0=pis;
	wk=(double *) calloc((size_t) (3*(n+1)), sizeof(double));

/*     Coefficients of recurrence relation */
	for(i=0; i<=n; ++i) {
		wk[i*3]=2.0;
		wk[1+i*3]=0.0;
		wk[2+i*3]=i*2.0;
	}
	
	ier=gausrc(n,w,a,wk,ri0);
	free(wk);
	return ier;
}
 
