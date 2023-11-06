/* Utility functions for complex arithmetic, complex numbers are
	represented as double array of length 2, with first element as
	real part and second element as imaginary part.

*/

/*	To find square root of a complex number

	CX : (input) Array of length 2 containing the complex number
		whose square root is required
	CR : (output) Array of length 2 containing the complex square root
	CR can be same as CX.

	Required functions : None

*/

#include <math.h>

void csqrt(double cx[2], double cr[2])

{
	double r,s;

	if(cx[0]==0.0 && cx[1]==0.0) {cr[0]=0.0; cr[1]=0.0; return;}
	else if(fabs(cx[0]) >= fabs(cx[1])) {
		r=cx[1]/cx[0];
		s=sqrt(fabs(cx[0]))*sqrt((1.0+sqrt(1.0+r*r))/2.0);
	}
	else {
		r=cx[0]/cx[1];
		s=sqrt(fabs(cx[1]))*sqrt((fabs(r)+sqrt(1.0+r*r))/2.0);
	}

	if(cx[0]>= 0.0) {cr[0]=s; cr[1]=cx[1]/(2.0*s);}
	else if(cx[1] >= 0.0) {cr[0]=fabs(cx[1])/(2.0*s); cr[1]=s;}
	else {cr[0]=fabs(cx[1])/(2.0*s); cr[1]=-s;}
	return;
}

