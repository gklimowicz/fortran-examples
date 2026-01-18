/*	To integrate a function over a finite interval using Gauss-Kronrod formula
	For use with ADPINT

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	DIF : (output) estimated (absolute) error achieved by the function
	N : (output) Number of function evaluations used
	F : (input) Name of the function to calculate the integrand

	Function F(X) must be supplied by the user

	The returned value is always zero.

	Required functions : F
*/

#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ))

{
	int k;
	double at,bt,fbt,r1,f1,f2;

/*	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
	WK7 are the weights for these points in Kronrod formula
	WK15 and AK15 are the weights and abscissas for the remaining points
	in Kronrod formula.
	Because of symmetry only half the points are given.
*/

	double w7[4] = {0.12948496616886969327e0, 0.27970539148927666790e0,
                        0.38183005050511894495e0, 0.41795918367346938775e0};
	double a7[4] = {0.94910791234275852452e0, 0.74153118559939443986e0,
                        0.40584515137739716690e0, 0.0};
	double wk7[4] = {0.06309209262997855329e0, 0.14065325971552591874e0,
                         0.19035057806478540991e0, 0.20948214108472782801e0};
	double wk15[4] = {0.02293532201052922496e0, 0.10479001032225018383e0,
                          0.16900472663926790282e0, 0.20443294007529889241e0};
	double ak15[4] = {0.99145537112081263920e0, 0.86486442335976907278e0,
                          0.58608723546769113029e0, 0.20778495500789846760e0};

	at=(b-a)/2.;
	bt=(b+a)/2.;
	fbt=f(bt);
	r1=w7[3]*fbt;
	*ri=wk7[3]*fbt;
	for(k=0; k<3; ++k) {
		f1=f(bt+at*a7[k]);
		f2=f(bt-at*a7[k]);
/*	7-point Gauss-Legendre formula */
		r1=r1+w7[k]*(f1+f2);
/*	15-point Kronrod formula */
		*ri=(*ri)+wk7[k]*(f1+f2);
	}

	for(k=0; k<4; ++k) *ri=(*ri)+wk15[k]*(f(bt+at*ak15[k]) + f(bt-at*ak15[k]));

	*ri=(*ri)*at;
	r1=r1*at;
	*dif=fabs(*ri-r1);
	*n=15;
	return 0;
}
