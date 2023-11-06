/*	To calculate the associated Legendre functions P_lm(X)

	L,M : (input) Order of function (L>=0), abs(M)<=L
	X : (input) Argument at which the value of polynomial is required
	P : (output) Array of length L+1, which will contain the
		calculated values of polynomials. P[j] will contain
		P_jM(X) for j>=M

	Required functions : None
*/

#include <math.h>

void plm(int l, int m, double x, double p[])

{
	int i,n,mm;
	double pm,rm;

	if(l<0 || abs(m)>l ) return;
 
/*	First compute P_MM(x) */
	mm=abs(m);
	pm=1.0;
	for(i=1; i<2*mm; i += 2) pm=pm*i;

	if(m<0) {
/*	Modify the normalisation factor */
		for(i=1; i<=2*mm; ++i) pm=pm/i;
		if(mm - 2*(mm/2) == 1) pm=-pm;
	}

	rm=mm/2.0;
	if(mm == 0) p[mm]=1.0;
	else if(fabs(x)<1.0) p[mm]=pm*pow(1.0-x*x, rm);
	else p[mm]=0.0;

/*	Use the recurrence relation to compute P_nM(x) */
	p[mm+1]=x*p[mm]*(2*mm+1)/(mm+1-m);
	for(n=mm+2; n<=l; ++n) p[n]=((2*n-1)*x*p[n-1]-(n-1+m)*p[n-2])/(n-m);

	return;
}
