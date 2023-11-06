/*	To calculate the coefficients of Chebyshev expansion of a given
	function, which can be evaluated at any required point
	This function uses orthogonality relation over discrete points
	to compute the coefficients.

	N : (input) Number of coefficients required, this number should
		be much larger than the actual number required. The
		accuracy of computed coefficients increases with N.
	C : (output) Array of length N containing the required
		coefficients. The constant term is C[0]/2, while C[I]
		is the coefficient of T_I(x).
	FUN : (input) Name of function to compute the given function
		
	Error status is returned by the value of the function CHEBEX.
		0 value implies successful execution
		613 implies N<10, in which case no calculations are done
	
	There is no test for accuracy and this has to be ascertained by
	doing another calculation with larger N (say 2N).

	Function FUN(X) must be supplied by the user

	Required functions : FUN
*/

#include <math.h>

#define PI 3.14159265358979324

int chebex(int n, double c[], double fun(double ))

{
	int i,j;
	double x,t0,t1,t2,fx;

	if(n<10) return 613;

	for(i=0; i<n; ++i) c[i]=0.0;
 
/*     evaluate the sums over n points */
	for(i=0; i<n; ++i) {
		x=cos((2*i+1)*PI/(2*n));
		fx=fun(x);
		c[0]=c[0]+fx; t0=1.0;
		c[1]=c[1]+fx*x; t1=x;
		for(j=2; j<n; ++j) {
			t2=2*x*t1-t0;
			c[j]=c[j]+fx*t2;
			t0=t1; t1=t2;
		}
	}

	for(i=0; i<n; ++i) c[i]=2.0*c[i]/n;
	return 0;
}
