/*	Complex roots of a nonlinear equation using Muller's method with deflation */

#include <stdio.h>
#include <math.h>

void fun(double x[2], double f[2]);
int muller(double cx1[], double cx2[], double cx3[], double reps, double aeps,
	void cf(double * , double * ), int nz, double czero[][2], double rmax);
int zroot(int n, double cx[][2], double czero[][2], int *nz, double reps,
	double aeps, double rmax, void cf(double * , double * ));

main()
{
	int i,i1,j, id, iflg, ier,n,nz;
	double xl, xu, yl, rmx, reps, aeps, cx[50][2], czero[50][2];

/*	Example 7.6 */

	aeps=1.e-7; reps=1.e-7; rmx=100.0; nz=0;

/*	The zeros will keep accumulating in array czero */
	for(i1=0; i1<99; ++i1) {
		printf("type n = no. of zeros to be tried     (quits when n<=0)\n");
		scanf(" %d", &n);
		if(n<=0) return 0;

		printf("type n  complex starting values \n");
		for(i=0; i<n; ++i) {
			scanf(" %le %le", &cx[i][0],&cx[i][1]);
			printf(" %e %e \n", cx[i][0],cx[i][1]);
		}

		i=zroot(n,cx,czero,&nz,reps,aeps,rmx,fun);
		printf(" ier = %d  no. of zeros found = %d \n", i,nz);
		for(i=0; i<nz; ++i) printf(" %d  %e  %e \n",i,czero[i][0],czero[i][1]);

	}
	return;
}

void fun(double x[2], double f[2])

{
	f[0]=x[0]+sin(x[0])*cosh(x[1]);
	f[1]=x[1]+cos(x[0])*sinh(x[1]);
	return;
}




/*	Complex zero of a given function using Muller iteration with deflation

	CX1,CX2,CX3 : (input/output) Complex starting values for iteration
		These will be updated during execution and CX3 should be
		the best estimate for the zero. These are treated as arrays of
		length 2 containing the real and imaginary part.
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(CX3))
	CF : (input) Name of the function routine to calculate the function value
	NZ : (input) Number of zeros already known (for deflation)
	CZERO : (input) Array of length NZ*2 containing the known
		zeros for deflation. It should contain the real and imaginary
		parts of known zeros.
	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
		when ABS(CX) > RMAX
		
	Error status is returned by the value of the function MULLER.
		0 value implies successful execution
		42 implies that roundoff errors appear to be dominating
			and calculations are terminated.
		43 implies that iteration failed to converge to specified accuracy
			but a weaker convergence criterion is satisfied
		404 implies 2 of the 3 starting values are identical and
			no calculations are done
		431 implies that iteration goes outside the specified limits
		432 implies that iteration failed to converge to specified accuracy
		433 implies that denominator for calculating the iteration
			vanishes and it is not possible to continue further
	
	Function CF(Z,F) must be supplied by the user. Here Z and F are
		both arrays of length 2 containing the complex values of the
		z and f(z) respectively.

	Required functions : CF, CABS, CDIV, CSQRT
*/

#include <math.h>

void csqrt(double cx[2], double cr[2]);
double cabs(double cx[2]);
void cdiv(double c1[2], double c2[2], double cr[2]);

int muller(double cx1[], double cx2[], double cx3[], double reps, double aeps,
	void cf(double * , double * ), int nz, double czero[][2], double rmax)

{
	int i,j, nit=50;
	double cf1[2],cf2[2],cf3[2],ch1[2],ch2[2],c1[2],c2[2],c3[2],cg[2],ci[2];
	double cd[2],cd1[2],cdx[2],dx,dx1,r, reps0=1.0e-4;

	if((cx1[0]==cx2[0]&&cx1[1]==cx2[1]) || (cx3[0]==cx2[0]&&cx3[1]==cx2[1]) ||
		(cx1[0]==cx3[0]&&cx1[1]==cx3[1]) ) return 404;

	cf(cx1,cf1); cf(cx2,cf2);

/*	Perform deflation */
	for(j=0; j<nz; ++j) {
		c1[0]=cx1[0]-czero[j][0]; c1[1]=cx1[1]-czero[j][1];
		cdiv(cf1,c1,cf1);
		c1[0]=cx2[0]-czero[j][0]; c1[1]=cx2[1]-czero[j][1];
		cdiv(cf2,c1,cf2);
	}

/*	The divided difference f[CX1,CX2] */
	c1[0]=cf2[0]-cf1[0]; c1[1]=cf2[1]-cf1[1];
	c2[0]=cx2[0]-cx1[0]; c2[1]=cx2[1]-cx1[1];
	cdiv(c1,c2,ch1);
	c1[0]=cx3[0]-cx2[0]; c2[1]=cx3[1]-cx2[1];
	dx1=cabs(c1);
		
/*	Loop for Muller iteration */
	for(i=1; i<=nit; ++i) {
		cf(cx3,cf3);
/*	Perform deflation */
		for(j=0; j<nz; ++j) {
			c1[0]=cx3[0]-czero[j][0]; c1[1]=cx3[1]-czero[j][1];
			cdiv(cf3,c1,cf3);
		}

		if((cx3[0]==cx1[0]&&cx3[1]==cx1[1]) || (cx3[0]==cx2[0]&&cx3[1]==cx2[1]) ) return 0;
		if(cf3[0]==0.0 && cf3[1]==0.0) return 0;
		c1[0]=cf3[0]-cf2[0]; c1[1]=cf3[1]-cf2[1];
		c3[0]=cx3[0]-cx2[0]; c3[1]=cx3[1]-cx2[1];
		cdiv(c1,c3,ch2);
		c1[0]=ch2[0]-ch1[0]; c1[1]=ch2[1]-ch1[1];
		c2[0]=cx3[0]-cx1[0]; c2[1]=cx3[1]-cx1[1];
		cdiv(c1,c2,c2);
		cg[0]=ch2[0]+c3[0]*c2[0]-c3[1]*c2[1];
		cg[1]=ch2[1]+c3[0]*c2[1]+c3[1]*c2[0];
		c1[0]=cg[0]*cg[0]-cg[1]*cg[1]-4.*cf3[0]*c2[0]+4.*cf3[1]*c2[1];
		c1[1]=2.*cg[0]*cg[1]-4.*cf3[0]*c2[1]-4.*cf3[1]*c2[0];
		csqrt(c1,ci);
		cd[0]=cg[0]+ci[0]; cd[1]=cg[1]+ci[1];
		cd1[0]=cg[0]-ci[0]; cd1[1]=cg[1]-ci[1];
		if(cabs(cd1) > cabs(cd)) {cd[0]=cd1[0]; cd[1]=cd1[1];}

		if(cd[0]==0.0 && cd[1]==0.0) return 433;

		cdiv(cf3,cd,cdx); cdx[0]=-2.*cdx[0]; cdx[1]=-2.*cdx[1];
		cx1[0]=cx2[0]; cx1[1]=cx2[1];
		cx2[0]=cx3[0]; cx2[1]=cx3[1];
/*	The new approximation to zero */
		cx3[0]=cx3[0]+cdx[0]; cx3[1]=cx3[1]+cdx[1];
		cf1[0]=cf2[0]; cf1[1]=cf2[1];
		cf2[0]=cf3[0]; cf2[1]=cf3[1];
		ch1[0]=ch2[0]; ch1[1]=ch2[1];

		dx=cabs(cdx);
		r=reps*cabs(cx3); if(aeps>r) r=aeps;
		if(i>2 && dx<r) return 0;
		r=reps0*cabs(cx3); if(aeps>r) r=aeps;
		if(i>5 && dx<r && dx>dx1) return 42;

		dx1=dx;
		if(cabs(cx3)>rmax) return 431;
	}

	r=reps0*cabs(cx3); if(aeps>r) r=aeps;
	if(dx<r) return 43;
	else return 432;
}




/*	Complex zeros of a given function using Muller iteration with deflation

	N : (input) Number of zeros to be determined
	CX : (input) Array of length N*2 containing the initial guesses
		for the complex zeros
	CZERO : (input/output) Array of length (N+NZ)*2 containing the
		computed values of the complex zeros. The first NZ zeros which are
		already known must be supplied as input while other zeros will be added
	NZ : (input/output) Number of known zeros. At input it should
		contain the number already known which are included in
		array CZERO. This number will be incremented as more zeros
		are determined successfully.
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))
	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
		when ABS(CZ) > RMAX
	CF : (input) Name of the function routine to calculate the function
		Function CF(Z,F) must be supplied by the user, where F and Z
		are arrays of length 2 containing the complex values. F should
		be the calculated value of the complex function at Z.
		For use with MULER2 the function is calculated
		in the form F*2**JF and Function CF(Z,F,JF) should be provided.
		
	Error status is returned by the value of the function ZROOT.
		0 value implies successful execution
		41 implies that iteration did not converge to satisfactory accuracy for
			at least one zero, but value is acceptable at lower accuracy.
		429 implies that iteration failed to find at least one zero.
			The number of zeros found can be checked from the value of NZ.

	Required functions : MULLER (or MULER2), CF, CABS, CDIV, CSQRT
*/

#include <math.h>

int muller(double cx1[], double cx2[], double cx3[], double reps, double aeps,
	void cf(double * , double * ), int nz, double czero[][2], double rmax);

int zroot(int n, double cx[][2], double czero[][2], int *nz, double reps,
	double aeps, double rmax, void cf(double * , double * ))

{
	int i,ier,ier1;
	double cdx[2],cx1[2],cx2[2],cx3[2];

	if(n<=0) return 0;
	if(*nz<0) *nz=0;
	ier=0;

	for(i=0; i<n; ++i) {
/*	The starting values for Muller's iteration */
		cx3[0]=cx[i][0]; cx3[1]=cx[i][1];
		cdx[0]=0.01*cx3[0]; cdx[1]=0.01*cx3[1];
		if(i<n-1 && cdx[0]==0.0 && cdx[1]==0.0) {
			cdx[0]=0.1*(cx[i+1][0]-cx[i][0]);
			cdx[1]=0.1*(cx[i+1][1]-cx[i][1]);
		}
		if(i>0 && cdx[0]==0.0 && cdx[1]==0.0) {
			cdx[0]=0.1*(cx[i-1][0]-cx[i][0]);
			cdx[1]=0.1*(cx[i-1][1]-cx[i][1]);
		}
		if(cdx[0]==0.0 && cdx[1]==0.0) {cdx[0]=1.0e-3; cdx[1]=0.0;}
/*	These values may need to be changed in some cases */
		cx2[0]=cx3[0]+cdx[0]; cx2[1]=cx3[1]+cdx[1];
		cx1[0]=cx3[0]-cdx[0]; cx1[1]=cx3[1]-cdx[1];

/*	Find the next zero */
		ier1=muller(cx1,cx2,cx3,reps,aeps,cf,*nz,czero,rmax);
		if(ier1<100) {
/*	The zero is accepted */
			czero[*nz][0]=cx3[0]; czero[*nz][1]=cx3[1];
			*nz=(*nz)+1;
			if(ier==0 && ier1 != 0) ier=41;
		}
		else ier=429;
	}
	return ier;
}




/* Utility functions for complex arithmetic, complex numbers are
	represented as double array of length 2, with first element as
	real part and second element as imaginary part.

*/

/*	For finding absolute value of a complex number

	CX : (input) Array of length 2 containing the complex number
	CABS will be the absolute value of CX

	Required functions : None

*/

#include <math.h>

double cabs(double cx[2])

{
	double r;

	if(fabs(cx[0])>fabs(cx[1])) {
		r=cx[1]/cx[0];
		return fabs(cx[0])*sqrt(1+r*r);
	}
	else if(cx[1] != 0.0) {
		r=cx[0]/cx[1];
		return fabs(cx[1])*sqrt(1+r*r);
	}
	else return 0.0;

}




/* Utility functions for complex arithmetic, complex numbers are
	represented as double array of length 2, with first element as
	real part and second element as imaginary part.

*/

/*	Complex division

	C1 : (input) Array of length 2 containing the numerator 
	C2 : (input) Array of length 2 containing the denominator
	CR : (output) Array of length 2 containing the result (C1/C2)
		CR can be same as  C1 or C2

	Required functions : None

*/

#include <math.h>

void cdiv(double c1[2], double c2[2], double cr[2])

{
	double r,r1,den;

	if(fabs(c2[0])>fabs(c2[1])) {
		r=c2[1]/c2[0];
		den=c2[0]+c2[1]*r;
/*	To avoid overwriting on c1 if c1 and cr are same */
		r1=(c1[0]+c1[1]*r)/den;
		cr[1]=(c1[1]-c1[0]*r)/den;
		cr[0]=r1;
	}
	else {
		r=c2[0]/c2[1];
		den=c2[0]*r+c2[1];
		r1=(c1[0]*r+c1[1])/den;
		cr[1]=(c1[1]*r-c1[0])/den;
		cr[0]=r1;
	}
	return;
}






/* Utility functions for complex arithmetic, complex numbers are
	represented as double array of length 2, with first element as
	real part and second element as imaginary part.

*/

/*	To find square root of a complex number

	CX : (input) Array of length 2 containing the complex number
		whose square root is required
	CR : (output) Array of length 2 containing the complex square root
	CR can be same as CX.

	Required functions : None

*/

#include <math.h>

void csqrt(double cx[2], double cr[2])

{
	double r,s;

	if(cx[0]==0.0 && cx[1]==0.0) {cr[0]=0.0; cr[1]=0.0; return;}
	else if(fabs(cx[0]) >= fabs(cx[1])) {
		r=cx[1]/cx[0];
		s=sqrt(fabs(cx[0]))*sqrt((1.0+sqrt(1.0+r*r))/2.0);
	}
	else {
		r=cx[0]/cx[1];
		s=sqrt(fabs(cx[1]))*sqrt((fabs(r)+sqrt(1.0+r*r))/2.0);
	}

	if(cx[0]>= 0.0) {cr[0]=s; cr[1]=cx[1]/(2.0*s);}
	else if(cx[1] >= 0.0) {cr[0]=fabs(cx[1])/(2.0*s); cr[1]=s;}
	else {cr[0]=fabs(cx[1])/(2.0*s); cr[1]=-s;}
	return;
}

