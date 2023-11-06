/*	To calculate the incomplete Beta function I_x(a,b) using
	the infinite series, called by BETAP

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln
*/


#include <math.h>
#include <stdio.h>

double gammaln(double x);

double betser(double a, double b, double x)
{
	double t,t1,tmax,b1,s3;
	int i;

	s3=1; t=a; tmax=a; i=-1;
	while(i<500 && fabs(t/(a+i+1))>1.e-15)
	{
		i=i+2; t1=t*x*(i-b)/i; t=t1*x*(i+1-b)/(i+1);
		s3=s3+t1/(a+i)+t/(a+i+1);
		if(fabs(t/(a+i+1))>tmax) tmax=fabs(t/(a+i+1));
	}
	if(s3<1.e-16*tmax) printf("  Roundoff error while evaluating the infinite series for incomplete Beta function at\n a= %e  b= %e  x= %e sum= %e  Max term= %e\n",a,b,x,s3,tmax);
	b1=a*log(x)+log(s3)-log(a);
	b1=b1+gammaln(a+b)-gammaln(a)-gammaln(b);
	return exp(b1);
}

