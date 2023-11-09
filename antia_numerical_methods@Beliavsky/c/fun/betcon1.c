/*	To calculate the incomplete Beta function I_x(a,b) using
	the continued fraction, called by betap

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln
*/


#include <math.h>
#include <stdio.h>

double gammaln(double x);

double betcon1(double a, double b, double x)
{
	double c[500],d[500],a1,b1,c1,c2;
	int i,m;

	c[0]=1; d[0]=1; c[1]=1; d[1]=1-x*(a+b)/(a+1);
	c2=c[1]/d[1]; c1=0.0; i=0;
	while(i<500 && fabs(c1-c2)>1.e-12)
	{
		i=i+2; m=i/2; c1=c2;
		c[i]=c[i-1]+m*(b-m)*x*c[i-2]/((a+i-1)*(a+i));
		d[i]=d[i-1]+m*(b-m)*x*d[i-2]/((a+i-1)*(a+i));
		c[i+1]=c[i]-(a+m)*(a+b+m)*x*c[i-1]/((a+i+1)*(a+i));
		d[i+1]=d[i]-(a+m)*(a+b+m)*x*d[i-1]/((a+i+1)*(a+i));
		if(fabs(c[i])>1.e200)
/*	Scale the numerator and denominator to prevent underflow/overflow */
		{
			c[i]=c[i]/1e200; c[i+1]=c[i+1]/1e200;
			d[i]=d[i]/1e200; d[i+1]=d[i+1]/1e200;
		}
		if(fabs(c[i])<1.e-200)
		{
			c[i]=c[i]*1e200; c[i+1]=c[i+1]*1e200;
			d[i]=d[i]*1e200; d[i+1]=d[i+1]*1e200;
		}
		c2=c[i+1]/d[i+1];
	}
	if(c2<0.0) printf(" ** Roundoff error while evaluating the continued fraction for incomplete Beta function at\n a= %e  b=%e  x= %e cont. frac.= %e\n",a,b,x,c2);
	b1=a*log(x)+b*log(1-x)+log(c2)-log(a);
	b1=b1+gammaln(a+b)-gammaln(a)-gammaln(b);
	return exp(b1);
}

