/*	To calculate the Bessel function of second kind of integral order
		for real argument

	N : (input) Order of Bessel function required, N may be negative
		or positive. For N=0,1 use BY0 and BY1 respectively.
	X : (input) Argument at which the value is required
	BY : (output) Array of length at least abs(N)+1
		which will contain the value of Bessel function of order
		0,1,...,ABS(N). BY[I] will contain Bessel function of order I
		or -I (if N<0).

	Required functions : BY0, BY1
*/

#include <math.h>

double by0(double x);
double by1(double x);

void byn(int n, double x, double by[])

{
	int i;
 
	by[0]=by0(x);
	by[1]=by1(x);
	if(x<= 0.0) return;
/*	Use the recurrence relation in Forward direction */
	for(i=2; i<=abs(n); ++i) by[i]=2*(i-1)*by[i-1]/x-by[i-2];

	if(n<0) {
		for(i=1; i<=abs(n); i +=2) by[i]= -by[i];
	}
	return;
}
