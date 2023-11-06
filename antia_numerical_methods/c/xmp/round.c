/*	To round a number to specified number of digits using the given base */

#include <stdio.h>
#include <math.h>

double round(double x, int n, double b);

main()
{
	int i,j,n;
	double b,x,xr;

	for(i=0; i<99; ++i) {
		printf("type x = no. to be rounded, n = no. of digits, b = base\n");
		printf("                        (quits when n<=0)\n");
		scanf("%lf %d %lf", &x,&n,&b);

		if(n<=0) return 0;
		xr=round(x,n,b);
		printf(" %f rounded to %f \n", x,xr);
	}
	return;
}

 

/*	To round a number to N digits using base B.

	X : (input) The number to be rounded
	N : (input) The number of digits required in rounded number
	B : (input) Base of number system to be used for rounding
	The returned value of the function should be the rounded value of X

	Required routines : None
*/
#include <math.h>

double round(double x, int n, double b)

{
	int lgx, nx;
	double rx, xa, fx;

	if(x==0.0) return 0.0;

	xa=fabs(x);
	lgx=log(xa)/log(b);
	if(xa<1.0) lgx=lgx-1;
	fx=pow(b,(double) (n-lgx-1));
	nx=xa*fx+0.5;
	rx=nx/fx;
	if(x<0) rx=-rx;
	return rx;
}

