/*	To calculate the incomplete Beta function I_x(a,b) using
	the integral, called by betap

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln, adpint, kronrd, fbeta
*/

#include <math.h>

	double fbeta(double x);
	int adpint(double *rint, double xl, double xu, double reps, double aeps,
		double *dif, double f1(double x), int *n, int *nmax);
	int kronrd(double *ri, double a, double b, double *dif, int *n,
		double f1(double x));
	double AA, BB;

double betai(double a, double b, double x)
{
	double xl,xu,reps,aeps,b1,rint,dif,f;
	int i,npt,nmax;

	AA=a; BB=b; xl=0; xu=x; reps=1e-14; aeps=1e-280;
	i=adpint(&rint,xl,xu,reps,aeps,&dif,fbeta,&npt,&nmax);
	b1=log(rint)+gammaln(a+b)-gammaln(a)-gammaln(b);
	f=exp(b1);
	if(f>1) f=1;
	return f;
}
	
