/*	To find coefficients of Chebyshev expansion */

#include <stdio.h>
#include <math.h>

int chebex(int n, double c[], double fun(double ));
double fun(double x);

main()
{
	int i,i1,j,n;
	double hh, rm, c[2000],a[40];

/*	Chebyshev coefficients for Arc Tan(x) */

/*	The number must be somewhat larger than the required number */
	for(i1=0; i1<99; ++i1) {
		printf("type n = No. of coefficients to be calculated \n");
		printf("                       (quits when n<=0)\n");
		scanf(" %d", &n);
		if(n<=0) return 0;

/*	Calculate the exact coefficients of expansion */
		rm=sqrt(2.0)-1.0;
		hh=1.0;
		for(i=0; i<=20; i += 2) {
			a[i]=0.0;
			a[i+1]=2.*hh*pow(rm, (double) (i+1))/(i+1.0);
			hh=-hh;
		}

		i=chebex(n,c,fun);
/*	Only the first 20 coefficients are printed out */
		printf(" ier= %d    n= %d \n    exact coef.     calculated coef.\n",i,n);
		for(i=0; i<20; ++i) printf(" %d  %e    %e \n",i,a[i],c[i]);
	}
	return;
}

double fun(double x)

{
	return atan(x);
}



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
