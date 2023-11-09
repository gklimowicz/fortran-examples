/*	To calculate the incomplete Beta function I_x(a,b) using
	the continued fraction (modified form), called by betap

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln
*/


#include <math.h>
#include <stdio.h>

double gammaln(double x);

double betcon(double a, double b, double x)
{
	double c[500],d[500],a1,b1,c1,c2,d1,d2,d3;
	int i;

	c[0]=a; d[0]=a*(1-x*(a+b)/(a+1));
	b1=1+x*(b-1)/(a+1)-x*(a+1)*(a+b+1)/(a+3);
	a1=a*x*x*(b-1)*(a+b)/pow(a+1,2.0);
	c[1]=a*b1; d[1]=b1*d[0]+a1; c2=c[1]/d[1]; c1=0.0; i=1;

	while(i<500 && fabs(c1-c2)>1.e-12)
	{
		i=i+1; c1=c2;
		d1=-x*(a+i-1)*(a+b+i-1)/((a+2*i-2)*(a+2*i-1));
		d3=-x*(a+i)*(a+b+i)/((a+2*i)*(a+2*i+1));
		d2=x*i*(b-i)/((a+2*i-1)*(a+2*i));
		c[i]=c[i-1]*(a+2*i)*(1+d2+d3)-c[i-2]*(a+2*i)*(a+2*i-2)*d2*d1;
		d[i]=d[i-1]*(a+2*i)*(1+d2+d3)-d[i-2]*(a+2*i)*(a+2*i-2)*d2*d1;
		if(fabs(c[i])>1.e200)
/*	Scale the numerator and denominator to prevent underflow/overflow */
		{
			c[i]=c[i]/1e200; c[i-1]=c[i-1]/1e200;
			d[i]=d[i]/1e200; d[i-1]=d[i-1]/1e200;
		}
		if(fabs(c[i])<1.e-200)
		{
			c[i]=c[i]*1e200; c[i-1]=c[i-1]*1e200;
			d[i]=d[i]*1e200; d[i-1]=d[i-1]*1e200;
		}
		c2=c[i]/d[i];
	}
	if(c2<0.0) printf(" ** Roundoff error while evaluating the modified continued fraction for incomplete Beta function at\n a= %e  b=%e  x= %e cont. frac.= %e\n",a,b,x,c2);
	b1=a*log(x)+b*log(1-x)+log(c2)-log(a);
	b1=b1+gammaln(a+b)-gammaln(a)-gammaln(b);
	return exp(b1);
}

