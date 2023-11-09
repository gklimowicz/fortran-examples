/*	To compute integral of a tabulated function using cubic splines

	A : (input) Lower limit of integration
	B : (input) Upper limit of integration (B > A)
	SINT : (output) Computed value of integral using cubic splines
	TINT : (output) Computed value of integral using trapezoidal rule
	N : (input) Number of tabular points
	X : (input) Array of length N containing abscissas (in ascending order)
	F : (input) Array of length N containing the function values
	C : (input) Array of length 3*N containing the coefficients of splines
		
	Error status is returned by the value of the function SPLINT.
		0 value implies successful execution
		31 implies lower limit A is outside the table
		32 implies upper limit B is outside the table
		301 implies A>B or X[0]>X[N-1] and no calculations are done
		Other values may be set by SPLEVL

	Required functions : SPLEVL
*/

#include <math.h>


double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier);

int splint(double a, double b, double *sint, double *tint, int n, double x[],
	double f[], double c[][3])

{
	int i,ier,ier1;
	double a1,f1,f2,dfb,ddfb,b1,d1,d2;

	*sint=0.0; *tint=0.0;
	if(a == b) return 0;
	if(a>b || x[0]>=x[n-1]) return 301;

	a1=a;
/*	Evaluate the function at x=A for trapezoidal rule */
	f1=splevl(a1,n,x,f,c,&dfb,&ddfb,&ier);
	if(ier>100) return ier;
	if(a<x[0] || a>x[n-1]) ier=31;
	if(b<x[0] || b>x[n-1]) ier=32;

/*	Integrating over the n-1 subintervals */
	for(i=0; i<n-1; ++i) {
		if(a1 < x[i+1] || i == n-2) {
			b1=x[i+1]; if(b<b1) b1=b;
			if(i == n-2) b1=b;
			d1=a1-x[i];
			d2=b1-x[i];

/*	Add integral of cubic spline over [A1,B1] */
			*sint = *sint+f[i]*(d2-d1)+c[i][0]*(d2*d2-d1*d1)/2.+
				c[i][1]*(d2*d2*d2-d1*d1*d1)/3.+
				c[i][2]*(d2*d2*d2*d2-d1*d1*d1*d1)/4.;
			f2=f[i+1];
			if(b1 != x[i+1]) f2=splevl(b1,n,x,f,c,&dfb,&ddfb,&ier1);

/*	Trapezoidal rule approximation to integral */
			*tint = *tint+0.5*(b1-a1)*(f1+f2);
			if(b<=x[i+1]) return ier;
			a1=b1;
			f1=f2;
		}
	}
	return ier;
}
