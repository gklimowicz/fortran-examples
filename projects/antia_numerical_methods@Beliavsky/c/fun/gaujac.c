/*     To calculate weights and abscissas of a Gauss-Jacobi quadrature formula

	N : (input) Number of points in the required quadrature formula
	ALP : (input) Exponent of 1-x in weight function
	BETA : (input) Exponent of 1+x in weight function
		w(x)=(1-x)**ALP*(1+x)**BETA
	W : (output) Array of length N, which will contain the weights
	A : (output) Array of length N containing the abscissas
		
	Error status is returned by the value of the function GAUJAC.
		0 value implies successful execution
		302 implies N<=0
		313 implies ALP<= -1 or BETA<= -1
		321 implies that some coefficient becomes imaginary
			during calculations.
		In both these cases calculations are abandoned.
		Other values may be set by TQL2

	Required functions : GAUSRC, TQL2, GAMMA
*/

#include <math.h>
#include <stdlib.h>

int gausrc(int n, double w[], double a[], double *cof, double ri0);

int gaujac(int n, double alp, double beta, double w[], double a[])

{
	int i,ier;
	double ri0,a1;
	double *wk;

	if(alp<= -1 || beta <= -1) return 313;

	ri0=pow(2.0,alp+beta+1.0)*gamma(alp+1.0)*gamma(beta+1.0)/gamma(alp+beta+2.0);
	wk=(double *) calloc((size_t) (3*(n+1)), sizeof(double));

/*     Coefficients of recurrence relation */
	wk[0]=1.0+(alp+beta)/2.0;
	wk[1]=(alp-beta)/2.0;
/*	wk(2) is not required and would give 0/0 form in some cases
	Hence do not calculate it. */

	for(i=1; i<=n; ++i) {
		a1=2.0*(i+1)*(i+1+alp+beta)*(2*i+alp+beta);
		wk[i*3]=(2*i+alp+beta)*(alp+beta+2*i+1)*(alp+beta+2*i+2)/a1;
		wk[1+i*3]=(2*i+1.0+alp+beta)*(alp*alp-beta*beta)/a1;
		wk[2+i*3]=2.0*(i+alp)*(i+beta)*(alp+beta+2*i+2)/a1;
	}
	
	ier=gausrc(n,w,a,wk,ri0);
	free(wk);
	return ier;
}
