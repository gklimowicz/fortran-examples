/*	To calculate the Chi square function for a nonlinear
	least squares fit, for use with function NMINF.
	Version of NLLSQ for use with NMINF, when derivatives are not available

	N : (input) Number of parameters to be fitted
	A : (input) Array of length N containing the parameter values
		at which the function needs to be evaluated
	F : (output) The function value at specified parameter values

	The data points are passed through global variables NNLSQ, FXLSQ,
	XLSQ, EFLSQ, FX1LSQ, which must be initialised in the calling
	function before calling NMINF

	NNLSQ : (input) The number of data points in the table to be fitted
	FXLSQ : (input) Array of length NPL containing the function values
	XLSQ : (input) Array of length NPL containing the abscissas at which
		function values are available
	EFLSQ : (input) Array of length NPL containing the estimated errors in 
		function values, for use in deciding the weights for each point
	FX1LSQ : (output) Array of length NPL containing the fitted value
		of the function at each tabular point

	The parameter NPL must be equal to the dimension of arrays as
	declared in the calling function.

	This function requires function FCN to calculate the required
	function which has to be fitted. There is no provision to pass on
	the name of this function and hence it must be changed explicitly
	to the required name. Function FCN(N,A,X,F) must be supplied
	by the user. Here N is the number of parameters, A is an array
	of length N containing the values of parameters and X is the value
	of independent variable where the fitting function needs to be evaluated.
	F is the calculated function value.

	Required functions : FCN
*/

#include <math.h>


/* To pass on these parameters to nllsq_f, this must be included before
   the calling NMINF */

void fcn(int n, double a[], double x, double *fa);

#define NPL 100
int NNLSQ;
double FXLSQ[NPL], XLSQ[NPL], EFLSQ[NPL], FX1LSQ[NPL];

void nllsq_f(int n, double a[], double *f)

{
	int i,j;
	double fa,r1;

/*     GENERATING THE FUNCTION FOR MINIMISATION */

	*f=0.0;

	for(i=0; i<NNLSQ; ++i) {
		fcn(n,a,XLSQ[i],&fa);
		FX1LSQ[i]=fa;

/*     SUM OF SQUARED DIFFERENCE */

		r1=(fa-FXLSQ[i])/EFLSQ[i];
		*f=(*f)+r1*r1;
	}
	return;
}
