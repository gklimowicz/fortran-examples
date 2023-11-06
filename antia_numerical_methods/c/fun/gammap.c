/*	To calculate the incomplete Gamma function P(a,x) for
	positive real arguments
	It returns a value of -1 for a.le.0 or x<0

	a : (input) Argument for the complete Gamma function
	x : (input) Upper limit of integration defining the incomplete
			Gamma function

	Required routines : gamma, gammaln
*/
 

#include <math.h>

double gamma(double xg);
double gammaln(double xg);

double gammap(double a, double x)
{
	double c[500],b[500],f,s1,s2,g1,t,c1,c2;
	int i;

	if(a<=0 || x<0) return -1.0; 

	if(x<3)
	{
/*	Use power series */
        	s1=1; t=a; i=0;
		while(i<500 && fabs(t/(a+i))>1.e-14)
		{
			i=i+1; t=-t*x/i;
			s1=s1+t/(a+i);
		}
		if(a<140) f=pow(x,a)*s1/gamma(a+1);
		else f=0.0;
	}
	else if(a<1.2*x)
	{
/*	Use continued fraction */
        	c[0]=1; b[0]=x; c[1]=1; b[1]=x+1-a;
		c2=c[1]/b[1]; c1=0.0; i=0;
		while(i<500 && fabs(c1-c2)>1.e-12)
		{
			i=i+2; c1=c2;
			c[i]=x*c[i-1]+(i/2)*c[i-2];
			b[i]=x*b[i-1]+(i/2)*b[i-2];
			c[i+1]=c[i]+(i/2+1-a)*c[i-1];
			b[i+1]=b[i]+(i/2+1-a)*b[i-1];
			if(fabs(b[i+1])>1.e200)
			{
				c[i]=c[i]/1.e200; c[i+1]=c[i+1]/1.e200;
				b[i]=b[i]/1.e200; b[i+1]=b[i+1]/1.e200;
			}
			if(fabs(b[i+1])<1.e-200)
			{
				c[i]=c[i]*1.e200; c[i+1]=c[i+1]*1.e200;
				b[i]=b[i]*1.e200; b[i+1]=b[i+1]*1.e200;
			}
			c2=c[i+1]/b[i+1];
		}
		g1=-x+a*log(x)+log(c2)-gammaln(a); f=1-exp(g1);
	}
	else
	{
/*	Use the power series for a>x */
		s2=1; t=1; i=0;
		while(i<500 && fabs(t)>1.e-14)
		{
			i=i+1; t=t*x/(a+i); s2=s2+t;
		}
		g1=-x+a*log(x)+log(s2)-gammaln(a+1); f=exp(g1);
	}
	return f;
}

