/*	To calculate the probability that two uncorrelated sequences
	of length (n+2) will give a correlation coefficient exceeding |x|
	For large even values of n the series can give large roundoff
	in which case the distribution is approximated by Normal distribution

	N  : (input) Length of data sets should be N+2
        XX : (input) The function will calculate the probability for
                correlation exceeding |x| for uncorrelated sequences

	Required routines : GAMMALN, ERF
*/

#include <math.h>

double gammaln(double xg);
double erf(double x);

#	define PI 3.14159265358979323846

double pcor(int n, double xx)

{
	double f,p,ss,pl,pmax,x;
/*	pis is sqrt(pi) */
	double pis=1.7724538509055144;
	int l,k;

	x=fabs(xx);
	if(x>1.0) return 0.0;
	if(2*(n/2)==n)
	{
/*	if n is even */
		l=(n-2)/2; ss=x; p=x; 
		pmax=p;
		for(k=1; k<=l; k++)
		{
			p=-(l-k+1)*p*x*x/k; ss=ss+p/(2*k+1.0);
			if(fabs(p)>pmax) pmax=fabs(p);
		}
		pl=gammaln((n+1.0)/2)-gammaln(n/2.0);
		f= 2.*exp(pl)*ss/pis;
		if(pmax>1.e5||f>1.0) f=erf(x*sqrt(n/2.0));
	}
	else
	{
/*	if n is odd */
		l=(n-3)/2; ss=sqrt(1-x*x); p=ss; 
		if(n==1) ss=0.0;
		for(k=1; k<=l; k++)
		{
			p=p*2*k*(1-x*x)/(2*k+1.0); ss=ss+p;
		}
		f=(asin(x)+x*ss)*2/PI;
	}
	return 1-f;
}
