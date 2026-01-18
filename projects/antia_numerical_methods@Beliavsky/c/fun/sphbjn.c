/*	To calculate the spherical Bessel function of integral order
	(j_n(x)=Sqrt(PI/(2x))*J_{n+1/2}(x)) for a real argument

	N : (input) Order of Bessel function required, N may be negative
		or positive. 
	XB : (input) Argument at which the value is required
	BJ : (output) Array of length at least 
		ABS(N)+16+MAX(25,5*SQRT(N))
		which will contain the value of Bessel function of order
		0,1,...,ABS(N). BJ[I] will contain Bessel function of order I
		or -I (if N<0)
		Remaining elements of array are used as scratch space

	Required functions : None
*/

#include <math.h>
#include <stdio.h>

void sphbjn(int n, double xb, double bj[])

{
	int i,j,n1;
/*	REPS should be less than the machine accuracy */
	double s,x,xa,t,t0, reps=1.e-17; 

	x=fabs(xb);
	if(xb==0.0) {
		bj[0]=1.0;
		for(i=1; i<=abs(n); ++i) bj[i]=0.0;
		return;
	}
	if(n>0) {
		if(n<x || n<2) {
/*	Use the recurrence relation in the forward direction  */
 
			bj[0]=sin(x)/x;
			bj[1]=sin(x)/(x*x)-cos(x)/x;
			for(i=2; i<=n; ++i) bj[i]=(2*i-1)*bj[i-1]/x-bj[i-2];
		}

		else if(x<=4.0) {
/*	Use series expansion to calculate  bj[n], bj[n-1] */
			xa=x*x/4.0;
			t0=1.0;
			for(i=1; i<=n; ++i) {
				t0=t0*x/(2.0*i+1);
				if(i>=n-1) {
					t=t0; s=t0;
					for(j=1; j<51; ++j) {
						t=-t*xa/(j*(j+i+0.5));
						s=s+t;
						if(fabs(t)<fabs(s*reps)) break;
					}
					bj[i]=s;
				}
			}
			for(i=n-1; i>=1; --i) bj[i-1]=(2.*i+1)*bj[i]/x-bj[i+1];
		}

		else {
/*	Use the recurrence relation in the backward direction */
			n1=5.*sqrt((double) n); if(n1<25) n1=25;
			if(x<n/6.0) n1=n1/log(n*0.5/x);
			n1=n1+n;
			if(n<x+15.0) n1=n1+15;
			if(n1-n<15) n1=n+15;
			bj[n1]=0.0;
			bj[n1-1]=1.0;
			for(i=n1-1; i>=1; --i) bj[i-1]=(2.0*i+1)*bj[i]/x-bj[i+1];

			s=bj[0]*x/sin(x);
			for(i=0; i<=n; ++i) bj[i]=bj[i]/s;

/*	If ABS(bj[n+2])<1./reps, then the required accuracy may not be achieved */
			if(fabs(bj[n+2]*reps)<1.0) printf(" sphbjn failed at n = %d  x = %e s = %e \n",n,x,bj[n+2]);
		}
	}

	else {
/*	For negative N use the recurrence relation in the forward direction  */
		bj[0]=sin(x)/x;
		bj[1]=cos(x)/x;
		for(i=2; i<=abs(n); ++i) bj[i]=(3-2*i)*bj[i-1]/x-bj[i-2];
	}

	if(xb<0.0) {
		for(i=1; i<=abs(n); i +=2) bj[i]=-bj[i];
	}
	return;
}
