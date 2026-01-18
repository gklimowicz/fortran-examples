/*	To calculate the modified Bessel function of first kind of
		integral order for real argument

	N : (input) Order of Bessel function required, N must be positive
	XB : (input) Argument at which the value is required
	BI : (output) Array of length at least N+16+MAX(25, 5*sqrt(N))
		which will contain the value of Bessel function of order
		0,1,...,n. BI[i] will contain Bessel function of order i
		Remaining elements of array are used as scratch space

	Required functions : BI0, BI1
*/
#include <stdio.h>
#include <math.h>

double bi0(double x);
double bi1(double x);

void bin(int n, double xb, double bi[])

{
	int i,j,n1,na;
/*	REPS should be less than the machine accuracy */
	double x,s,t,t0,xa, reps=1.e-17;

	x=fabs(xb);
	na=abs(n);
	if(xb==0.0) {
		bi[0]=1.0;
		for(i=1; i<=na; ++i) bi[i]=0.0;
		return;
	}

	else if(na<x-10 || na<2) {
/*	Use the recurrence relation in the forward direction  */
		bi[0]=bi0(x);
		bi[1]=bi1(x);
		for(i=2; i<=na; ++i) bi[i]=-2*(i-1)*bi[i-1]/x+bi[i-2];
	}

	else if(x<=4.0) {
/*	Use series expansion to calculate  bi[na], bi[na-1] */
		xa=x*x/4.0;
		t0=x/2.0;
		for(i=2; i<=na; ++i) {
			t0=t0*x/(2.0*i);
			if(i>=na-1) {
				t=t0; s=t0;
				for(j=1; j<51; ++j) {
					t=t*xa/(j*(j+i));
					s=s+t;
					if(fabs(t)<fabs(s*reps)) break;
				}
				bi[i]=s;
			}
		}
		for(i=na-1; i>=1; --i) bi[i-1]=2.*i*bi[i]/x+bi[i+1];
	}


	else {
/*	Use the recurrence relation in the backward direction  */
		n1=5.*sqrt((double) na); if(n1<25) n1=25;
		if(x<na/6.0) n1=n1/log(na*0.5/x);
		n1=n1+na;
		if(na<x+15.0) n1=n1+15;
		if(n1-na<15) n1=na+15;
		bi[n1]=0.0; bi[n1-1]=1.0;
		for(i=n1-1; i>=1; --i) bi[i-1]=2.*i*bi[i]/x+bi[i+1];

		s=bi[0]/bi0(x);
		for(i=0; i<=na; ++i) bi[i]=bi[i]/s;
/*	If fabs(bi[na+1])<1/reps, then the required accuracy may not be achieved
	hence printout an error message */
		if(fabs(bi[na+1]*reps)< 1.0) printf(" bin failed at n = %d  x = %e  s = %e \n", n,x,bi[na+1]);
	}

	if(xb<0.0) {
		for(i=1; i<=na; i += 2) bi[i]= -bi[i];
	}
	return;
}

 
