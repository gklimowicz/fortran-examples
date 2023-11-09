/*     To calculate weights and abscissas of a Gauss-Laguerre quadrature formula

	N : (input) Number of points in the required quadrature formula
	ALP : (input) Exponent of x in weight function, w(x)=x**ALP*EXP(-x)
	W : (output) Array of length N, which will contain the weights
	A : (output) Array of length N containing the abscissas
		
	Error status is returned by the value of the function LAGURW.
		0 value implies successful execution
		302 implies N<= 0
		313 implies ALP<= -1
		321 implies that some coefficient becomes imaginary
			during calculations.
		In both these cases calculations are abandoned.
		Other values may be set by TQL2

	Required functions : GAUSRC, TQL2, GAMMA
*/

#include <math.h>
#include <stdlib.h>

int gausrc(int n, double w[], double a[], double *cof, double ri0);

int lagurw(int n, double alp, double w[], double a[])

{
	int i,ier;
	double ri0;
	double *wk;

	if(alp<= -1) return 313;

	ri0=gamma(alp+1.0);
	wk=(double *) calloc((size_t) (3*(n+1)), sizeof(double));

/*     Coefficients of recurrence relation */
	for(i=0; i<=n; ++i) {
		wk[i*3]=-(alp+i+1.0)/(i+1.0);
		wk[1+i*3]=(2*i+1.0+alp)*(alp+i+1.0)/(i+1.0);
		wk[2+i*3]=(i+alp)*(i+alp)*(i+1+alp)/(i+1.0);
	}
	
	ier=gausrc(n,w,a,wk,ri0);
	free(wk);
	return ier;
}
