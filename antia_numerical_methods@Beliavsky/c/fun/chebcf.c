/*	To calculate the coefficients of Chebyshev expansion from those of
	power series or vice versa.

	N : (input) The degree of polynomial
	C : (input/output) Array of length N+2 containing the coefficients
		of Chebyshev series, C[I] is the coefficient of T_I(x) in
		expansion. These coefficients would be calculated if IFLG=0,
		otherwise they must be supplied
		In the latter case the contents of C are destroyed
	P : (input/output) Array of length N+2 containing the coefficients
		of power series. P[I] is the coefficient of x**I.
		These coefficients must be supplied if IFLG=0,
		otherwise they would be calculated.
	IFLG : (input) Flag to decide which coefficients need to be calculated
		If IFLG=0 then coefficients in Chebyshev expansion are
		calculated. In this case P[I] must be supplied.
		otherwise coefficients in power series expansion are
		calculated. In this case C[I] must be supplied and contents
		of array C are destroyed during calculations

	The returned value is always zero.

	Required functions : None
*/

#include <math.h>

int chebcf(int n, double c[], double p[], int iflg)

{
	int i,j;
	double cj,cj1;

	if(iflg==0) {
/*	Calculate coefficients of Chebyshev expansion */
		for(i=1; i<n+2; ++i) c[i]=0.0;
		c[0]=p[n];
/*	Nested multiplication */
		for(i=0; i<n; ++i) {
			cj=c[1];
			c[1]=c[0]+0.5*c[2];
			c[0]=p[n-1-i]+0.5*cj;
			for(j=2; j<=i+1; ++j) {
				cj1=0.5*(cj+c[j+1]);
				cj=c[j];
				c[j]=cj1;
			}
		}
	}

	else {
/*	Calculate coefficients in power series */
		if(n<=1) {
			p[0]=c[0];
			if(n==1) p[1]=c[1];
			return 0;
		}
		for(i=0; i<n-2; ++i) {
			cj=c[n-2-i];
			c[n-2-i]=2.*c[n-1-i];
			c[n-1-i]=2.*c[n-i];
			for(j=n-3-i; j>=1; --j) {
				cj1=2.*cj-c[j+2];
				cj=c[j];
				c[j]=cj1;
			}
			cj1=cj-0.5*c[2];
			cj=c[0];
			c[0]=cj1;
			p[i]=cj-0.5*c[1];
		}

		cj=c[0];
		c[0]=c[1];
		c[1]=2.*c[2];
		p[n-2]=cj-0.5*c[1];
		p[n-1]=c[0];
		p[n]=c[1];
	}
	return 0;
}
