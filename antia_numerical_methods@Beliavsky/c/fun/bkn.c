/*	To calculate the modified Bessel function of second kind of
		integral order for real argument
		The function value is not defined for x<=0 and no calculations
		are done in that case but no warning is issued

	N : (input) Order of Bessel function required, N must be positive
	X : (input) Argument at which the value is required
	BK : (output) Array of length at least abs(N)+1
		which will contain the value of Bessel function of order
		0,1,...,N. BK[i] will contain Bessel function of order i

	Required functions : BK0, BK1
*/

#include <math.h>

double bk0(double x);
double bk1(double x);

void bkn(int n, double x, double bk[])

{
	int i;
 
	bk[0]=bk0(x);
	bk[1]=bk1(x);
	if(x<= 0.0) return;
/*	Use the recurrence relation in Forward direction */
	for(i=2; i<=abs(n); ++i) bk[i]=2*(i-1)*bk[i-1]/x+bk[i-2];
	return;
}
 
