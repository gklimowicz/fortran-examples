/*	To integrate a function over a finite interval using 16 point
	Gauss-Legendre formula, for use with ADPINT

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

int gaus16(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ))

{
	int k;
	double at,bt,fbt,r1,f1,f2;

/*	W8 and A8 are the weights and abscissas for the 8-point Gauss formula
	W16 and A16 are the weights and abscissas for the 16-point Gauss formula
	Because of symmetry only half the points are given.
*/

	double w8[4] = {0.10122853629037625915e0, 0.22238103445337447054e0,
                        0.31370664587788728734e0, 0.36268378337836198297e0};
	double a8[4] = {0.96028985649753623168e0, 0.79666647741362673959e0,
                        0.52553240991632898582e0, 0.18343464249564980494e0};
	double w16[8] ={0.02715245941175409485e0, 0.06225352393864789286e0,
                        0.09515851168249278481e0, 0.12462897125553387205e0,
                        0.14959598881657673208e0, 0.16915651939500253819e0,
                        0.18260341504492358887e0, 0.18945061045506849629e0};
	double a16[8] ={0.98940093499164993260e0, 0.94457502307323257608e0,
                        0.86563120238783174388e0, 0.75540440835500303390e0,
                        0.61787624440264374845e0, 0.45801677765722738634e0,
                        0.28160355077925891323e0, 0.09501250983763744019e0};

	at=(b-a)/2.;
	bt=(b+a)/2.;
	r1=0.0;
/*	8-point Gauss-Legendre formula */
	for(k=0; k<4; ++k) r1=r1+w8[k]*(f(bt+at*a8[k])+f(bt-at*a8[k]));

	*ri=0.0;
/*	16-point Gauss-Legendre formula */
	for(k=0; k<8; ++k) *ri=(*ri)+w16[k]*(f(bt+at*a16[k])+f(bt-at*a16[k]));

	*ri=(*ri)*at;
	r1=r1*at;
	*dif=fabs(*ri-r1);
	*n=24;
	return 0;
}
