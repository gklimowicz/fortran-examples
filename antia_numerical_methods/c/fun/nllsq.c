/*	Function to calculate the Chi square function for a nonlinear
	least squares fit, for use with function BFGS.

	N : (input) Number of parameters to be fitted
	A : (input) Array of length N containing the parameter values
		at which the function needs to be evaluated
	F : (output) The function value at specified parameter values
	G : (output) Array of length N containing the gradient vector
		G[I] will contain dF/dA[I]

	The data points are passed through global variables NNLSQ, FXLSQ,
	XLSQ, EFLSQ, FX1LSQ, which must be initialised in the calling
	function before calling BFGS

	NNLSQ : (input) The number of data points in the table to be fitted
	FXLSQ : (input) Array of length NPL containing the function values
	XLSQ : (input) Array of length NPL containing the abscissas at which
		function values are available
	EFLSQ : (input) Array of length NPL containing the estimated errors in 
		function values, for use in deciding the weights for each point
		Although EFLSQ should contain the error estimate, but in many
		cases it is found that multiplying all values by a suitable
		constant can improve the convergence of BFGS dramatically
	FX1LSQ : (output) Array of length NPL containing the fitted value
		of the function at each tabular point

	The parameter NPL must be equal to the dimension of arrays as
	declared in the calling function.

	This function requires function FCN to calculate the required
	function which has to be fitted. There is no provision to pass on
	the name of this function and hence it must be changed explicitly
	to the required name. Function FCN(N,A,X,F,DF) must be supplied
	by the user. Here N is the number of parameters, A is an array
	of length N containing the values of parameters and X is the value
	of independent variable where the fitting function needs to be evaluated.
	F is the calculated function value and DF is an array of length N
	containing the calculated derivatives.

	Required functions : FCN
*/

#include <math.h>
#include <stdlib.h>


/* To pass on these parameters to nllsq, this must be included before
   the calling BFGS */

void fcn(int n, double a[], double x, double *fa, double df[]);

#define NPL 100
int NNLSQ;
double FXLSQ[NPL], XLSQ[NPL], EFLSQ[NPL], FX1LSQ[NPL];

void nllsq(int n, double a[], double *f, double g[])

{
	int i,j;
	double fa,r1;
	double *df;

/*     GENERATING THE FUNCTION FOR MINIMISATION */

	*f=0.0;
	for(i=0; i<n; ++i) g[i]=0.0;

	df=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<NNLSQ; ++i) {
		fcn(n,a,XLSQ[i],&fa,df);
		FX1LSQ[i]=fa;

/*     SUM OF SQUARED DIFFERENCE AND ITS GRADIENT */

		r1=(fa-FXLSQ[i])/EFLSQ[i];
		*f=(*f)+r1*r1;
		for(j=0; j<n; ++j) g[j]=g[j]+2.*r1*df[j]/EFLSQ[i];
	}
	free(df);
	return;
}
