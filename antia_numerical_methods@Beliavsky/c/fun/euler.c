/*	Summation of alternating series using Euler transform

	N : (input) Number of terms to be summed.
		If N<1 it is assumed to be infinite
	M1 : (input/output) Number of terms in the beginning to be summed
		separately. If M1+M2+2*NMAX >= N then M1 is set to N and
		the sum is evaluated by direct summation
	M2 : (input) Number of terms at the end to be summed separately.
		The Euler transform is applied to terms from M1 to N-M2-1
	A0 : (input) The sign of first term. The sum assuming the first
		term to be positive is multiplied by A0
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(SUM))
	DIF : (output) Estimated (absolute) error achieved by the function
	N1 : (output) The number of terms actually used at the beginning
	N2 : (output) The number of terms actually used at the end
		N2 is relevant only for finite series
	SUM : (output) The calculated value of the sum
	TERM : (input) Name of the function to calculate the terms of
		series, function TERM(I) should calculate the Ith term
		without sign. The Ith term of the series will be
		A0*TERM(I)*(-1)**I. Here I is expected to range from 0 to N-1
		
	Error status is returned by the value of the function EULER.
		0 value implies successful execution
		31 implies M1+M2+2*NMAX >= N, in which case direct
			sum is calculated
		32 implies that the transformed series did not converge to
			specified accuracy at lower end
		34 implies that the transformed series did not converge to
			specified accuracy at upper end
		36 implies that the transformed series did not converge to
			specified accuracy at both ends

	FUNCTION TERM(I) must be supplied by the user

	Required functions : TERM
*/

#include <math.h>


int euler(int n, int *m1, int *m2, double a0, double reps, double aeps,
	double *dif, int *n1, int *n2, double *sum, double (*term) (int i))

{
	int i,j,k,ier, nmax=20;
	double a,t1,ts,r1,s1,s2,d[20][20];

	ier=0; *dif=0.0; *n2=0;
	if(((*m1)+(*m2)+2*nmax >= n) && (n>0)) {*m1=n; *n1=n; ier=31;}

/*	Sum the first M1 terms separately */
	s1=0.0; a=a0;
	for(i=0; i<(*m1); ++i) {
		s1=s1+a*term(i);
		a=-a;
	}
	*sum=s1;
	if((*m1>=n) && (n>0)) return ier;
 
/*	Sum the remaining terms using Euler transform */
	s2=0.0; t1=-1;
	for(i=(*m1); i<(*m1)+nmax; ++i) {
		d[i-(*m1)][0]=term(i);
/*	The differences */
		for(j=1; j<=i-(*m1); ++j) d[i-(*m1)][j]=d[i-(*m1)][j-1]-d[i-(*m1)-1][j-1];
		t1=-t1/2.;
		ts=t1*d[i-(*m1)][i-(*m1)];
		s2=s2+ts;
		r1=reps*fabs(s1+s2); if(aeps>r1) r1=aeps;
		if(fabs(ts) < r1) break;
	}
	if(i>=(*m1)+nmax) {ier=32; i=(*m1)+nmax-1;}

	*n1=i+1;

/*	Sum of the infinite series */
	*sum=s1+a*s2;
	*dif=fabs(ts);
	if(n<=0) return ier;
 
/*	For finite series calculate the contribution from the other end also */
	s1=0.0;
	a=a0; if(n-2*(n/2) == 0) a=-a0;
/*	Sum the last M2 terms separately */
	for(i=0; i<(*m2); ++i) {s1=s1+a*term(n-i-1); a= -a;}
 
/*	Sum the remaining terms using Euler transform */
	s2=0.0;
	t1=-1;
	for(i=(*m2); i<(*m2)+nmax; ++i) {
		d[i-(*m2)][0]=term(n-i-1);
		for(j=1; j<=i-(*m2); ++j) d[i-(*m2)][j]=d[i-(*m2)][j-1]-d[i-(*m2)-1][j-1];
		t1=-t1/2.;
		ts=t1*d[i-(*m2)][i-(*m2)];
		s2=s2+ts;
		r1=reps*fabs(*sum+s1+s2); if(aeps>r1) r1=aeps;
		if(fabs(ts) < r1) break;
	}

	if(i>=(*m2)+nmax) {ier=ier+4; i=(*m2)+nmax-1;}

	*n2=i+1;
/*	Add the two parts to get the sum of finite series */
	*dif=(*dif)+fabs(ts);
	*sum=(*sum)+s1+a*s2;
	if(ier>0 && ier<30) ier=ier+30;
	return ier;
}
