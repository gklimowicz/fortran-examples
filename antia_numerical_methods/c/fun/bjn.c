/*	To calculate the Bessel function of integral order for real argument

	N : (input) Order of Bessel function required, N may be negative
		or positive. For N=0,1 use BJ0 and BJ1 respectively.
	XB : (input) Argument at which the value is required
	BJ : (output) Array of length at least 
		abs(N)+16+MAX(25,5*sqrt(N))
		which will contain the value of Bessel function of order
		0,1,...,abs(N). BJ[i] will contain Bessel function of order i
		or -i (if N<0)
		Remaining elements of array are used as scratch space

	Required functions : BJ0, BJ1
*/

#include <math.h>

#include <stdio.h>

double bj0(double x);
double bj1(double x);

void bjn(int n, double xb, double bj[])

{
	int i,j,n1,na;
/*	REPS should be less than the machine accuracy */
	double x,s,xa,t0,t, reps=1.e-17;

	x=fabs(xb);
	na=abs(n);
	if(xb==0.0) {
		bj[0]=1.0;
		for(i=1; i<=na; ++i) bj[i]=0.0;
		return;
	}

	else if(na<x || na<2) {
/*	Use the recurrence relation in the forward direction  */
		bj[0]=bj0(x);
		bj[1]=bj1(x);
		for(i=2; i<=na; ++i) bj[i]=2*(i-1)*bj[i-1]/x-bj[i-2];
	}

	else if(x<=4.0) {
/*	Use series expansion to calculate  bj[na], bj[na-1] */
		xa=x*x/4.0;
		t0=x/2.0;
		for(i=2; i<=na; ++i) {
			t0=t0*x/(2.0*i);
			if(i>=na-1) {
				t=t0; s=t0;
				for(j=1; j<51; ++j) {
					t=-t*xa/(j*(j+i));
					s=s+t;
					if(fabs(t)<fabs(s*reps)) break;
				}
				bj[i]=s;
			}
		}
		for(i=na-1; i>=1; --i) bj[i-1]=2.*i*bj[i]/x-bj[i+1];
	}

	else {
/*	Use the recurrence relation in the backward direction  */
		n1=5.*sqrt((double) na); if(n1<25) n1=25;
		if(x<na/6.0) n1=n1/log(na*0.5/x);
		n1=n1+na;
		if(na<x+15.0) n1=n1+15;
		if(n1-na<15) n1=na+15;
		bj[n1]=0.0; bj[n1-1]=1.0;
		for(i=n1-1; i>=1; --i) bj[i-1]=2.*i*bj[i]/x-bj[i+1];

		s=bj[0];
		for(i=2; i<n1; i += 2) s=s+2.*bj[i];

		for(i=0; i<=na; ++i) bj[i]=bj[i]/s;
/*	If fabs(bj(na+1))<1/reps, then the required accuracy may not be achieved
	hence printout an error message */
		if(fabs(bj[na+1]*reps)< 1.0) printf(" bjn failed at n = %d  x = %e s = %e \n", n,x,bj[na+1]);
	}

	if((n<0) != (xb<0.0)) {
		for(i=1; i<=na; i += 2) bj[i]= -bj[i];
	}
	return;
}
