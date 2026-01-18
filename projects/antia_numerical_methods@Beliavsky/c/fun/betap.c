/*	To calculate the incomplete Beta function I_x(a,b) for
	positive real arguments
	It returns a value of -1 for a<=0 or b<=0 or x<0 or x>1

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln, betai, betser, betcon, betcon1,
			    adpint, kronrd, fbeta
*/

#include <math.h>
#include <stdio.h>

	double betai(double a, double b, double x);
	double betser(double a, double b, double x);
	double betcon(double a, double b, double x);
	double betcon1(double a, double b, double x);
	double fbeta(double x);
	int adpint(double *rint, double xl, double xu, double reps, double aeps,
		double *dif, double f1(double x), int *n, int *nmax);
	int kronrd(double *ri, double a, double b, double *dif, int *n,
		double f1(double x));


double betap(double a, double b, double x)
{
	double amax,amin,betal,f;
	
	if(a<=0 || b<=0 || x<0 || x>1) return -1.0;
	amax=a; if(b>a) amax=b;
	amin=a; if(b<a) amin=b;
	betal=(a+b)*log(a+b)-a*log(a)-b*log(b);
	if(amax<=30)
	{
		if(x<=0.5) f=betser(a,b,x);
		else f=1-betser(b,a,1-x);
	}
	else if(b<=20 && x<=0.71) f=betser(a,b,x);
	else if(a<=20 && x>=0.3) f=1-betser(b,a,1-x);
	else if(b<=50 && x<=0.35) f=betser(a,b,x);
	else if(a<=50 && x>=0.65) f=1-betser(b,a,1-x);
	else if(b<=100 && x<=0.18) f=betser(a,b,x);
	else if(a<=100 && x>=0.82) f=1-betser(b,a,1-x);
	else if(b<=180 && x<=0.1) f=betser(a,b,x);
	else if(a<=180 && x>=0.9) f=1-betser(a,b,x);
	else if(x<0.5)
	{
		if(a<2) f=betcon(a,b,x);
		else if(betal>700) f=betcon1(a,b,x);
		else f=betai(a,b,x);
	}
	else
	{
		if(b<2) f=1-betcon(b,a,1-x);
		else if(betal>700) f=1-betcon1(b,a,1-x);
		else f=1-betai(b,a,1-x);
	}
	if(f<0 || f>1) printf("error in evaluating incomplete beta function at a = %e  b= %e  x= %e  betap= %e\n",a,b,x,f);
	return f;
}

