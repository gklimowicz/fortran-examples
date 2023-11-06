/*	Complex roots of a nonlinear equation using Muller's method with deflation */

#include <stdio.h>
#include <math.h>

void fun(double x[2], double f[2], int *ix);
void muler2(double cx1[], double cx2[], double cx3[], double reps, double aeps,
	int *ier, double cf[], double cx[], int ix, int nz, double czero[][2], double rmax);
int zroot2(int n, double cx[][2], double czero[][2], int *nz, double reps,
	double aeps, double rmax, void cf(double * , double * , int * ));

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

		i=zroot2(n,cx,czero,&nz,reps,aeps,rmx,fun);
		printf(" ier = %d  no. of zeros found = %d \n", i,nz);
		for(i=0; i<nz; ++i) printf(" %d  %e  %e \n",i,czero[i][0],czero[i][1]);

	}
	return;
}

void fun(double x[2], double f[2], int *ix)

{
	f[0]=x[0]+sin(x[0])*cosh(x[1]);
	f[1]=x[1]+cos(x[0])*sinh(x[1]);
	ix=0;
	return;
}




/*	Complex zero of a given function using Muller iteration with deflation
	Function is assumed to be calculated as CF*2**IX
	This function uses reverse communication to calculate function
	values. If IER<0 the function should be evaluated and MULER2
	should be called back with new function value. Calculation of
	function value should not change any other variables in the
	call statement. First call should be with IER=0

	CX1,CX2,CX3 : (input/output) Complex starting values for iteration
		These will be updated during execution and CX3 should be
		the best estimate for the zero. These are treated as arrays of
		length 2 containing the real and imaginary part.
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(CX3))
	IER : (input/output) Error parameter, IER=0 implies successful execution
		For the first call IER should be set to 0.
		IER<0 implies that execution is not over and the function
			needs a new function evaluation at z=CX. After calculating
			the function value MULER2 should be called back.
		IER=42 implies that roundoff errors appear to be dominating
			and calculations are terminated.
		IER=43 implies that iteration failed to converge to specified accuracy
			but a weaker convergence criterion is satisfied
		IER=404 implies 2 of the 3 starting values are identical and
			no calculations are done
		IER=431 implies that iteration goes outside the specified limits
		IER=432 implies that iteration failed to converge to specified accuracy
		IER=433 implies that denominator for calculating the iteration
			vanishes and it is not possible to continue further
	CF : (input) Array of length 2 used to transmit the calculated
		value of the complex function at z=CX
		If the function exits with IER<0, then the calling function
		should calculate the function value at CX and call MULER2
		with this value stored in CF and IX. Other variables should
		not be changed.
	CX : (output) Array of length 2 specifying the value of z at
		which the function evaluation is needed (if IER <0)
	IX : (input) The exponent of function value, The function value 
		should be CF*2**IX
	NZ : (input) Number of zeros already known (for deflation)
	CZERO : (input) Array of length NZ*2 containing the known
		complex zeros for deflation. It should contain the real and
		imaginary parts of the zeros.
	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
		when ABS(CX) > RMAX
	
	Required functions : CABS, CSQRT, CDIV
*/

#include <math.h>

void csqrt(double cx[2], double cr[2]);
double cabs(double cx[2]);
void cdiv(double c1[2], double c2[2], double cr[2]);

void muler2(double cx1[], double cx2[], double cx3[], double reps, double aeps,
	int *ier, double cf[], double cx[], int ix, int nz, double czero[][2], double rmax)

{
	static int i,j,ifl,if1,if2,if3, nit=50;
	static double cf1[2],cf2[2],cf3[2],ch1[2],ch2[2],c1[2],c2[2],c3[2],cg[2],ci[2];
	static double cd[2],cd1[2],cdx[2],cf1i[2],cf2a[2],ch1a[2],dx,dx1,r, reps0=1.0e-4;

	ifl=-(*ier);
	if(ifl<=0 || ifl>3) {
/*	Initial call to the  function, start from beginning */

		if((cx1[0]==cx2[0]&&cx1[1]==cx2[1]) || (cx3[0]==cx2[0]&&cx3[1]==cx2[1]) ||
			(cx1[0]==cx3[0]&&cx1[1]==cx3[1]) ) {*ier=404; return;}

		*ier=-1;
		cx[0]=cx1[0]; cx[1]=cx1[1];
/*	To evaluate the function value at CX1 */
		return;
	}

	else if(ifl==1) {
		cf1[0]=cf[0]; cf1[1]=cf[1]; if1=ix;
/*	Perform deflation */
		for(j=0; j<nz; ++j) {
			c1[0]=cx1[0]-czero[j][0]; c1[1]=cx1[1]-czero[j][1];
			cdiv(cf1,c1,cf1);
		}

		*ier=-2;
		cx[0]=cx2[0]; cx[1]=cx2[1];
/*	To evaluate the function value at CX2 */
		return;
	}

	else if(ifl==2) {
		cf2[0]=cf[0]; cf2[1]=cf[1]; if2=ix;
/*	Perform deflation */
		for(j=0; j<nz; ++j) {
			c1[0]=cx2[0]-czero[j][0]; c1[1]=cx2[1]-czero[j][1];
			cdiv(cf2,c1,cf2);
		}

		r=pow(2.0, (double ) if1-if2);
		cf1i[0]=cf1[0]*r; cf1i[1]=cf1[1]*r;
		c1[0]=cf2[0]-cf1i[0]; c1[1]=cf2[1]-cf1i[1];
		c2[0]=cx2[0]-cx1[0]; c2[1]=cx2[1]-cx1[1];
/*	The divided difference f[CX1,CX2] */
		cdiv(c1,c2,ch1);
		c1[0]=cx3[0]-cx2[0]; c2[1]=cx3[1]-cx2[1];
		dx1=cabs(c1);

		i=0;
		*ier=-3;
		cx[0]=cx3[0]; cx[1]=cx3[1];
/*	To evaluate the function value at CX3 */
		return;
	}

	else if(ifl==3) {
/*	Loop for Muller iteration */
		i=i+1;
		cf3[0]=cf[0]; cf3[1]=cf[1]; if3=ix;
/*	Perform deflation */
		for(j=0; j<nz; ++j) {
			c1[0]=cx3[0]-czero[j][0]; c1[1]=cx3[1]-czero[j][1];
			cdiv(cf3,c1,cf3);
		}

		*ier=0;
		if((cx3[0]==cx1[0]&&cx3[1]==cx1[1]) || (cx3[0]==cx2[0]&&cx3[1]==cx2[1]) ) return;
		if(cf3[0]==0.0 && cf3[1]==0.0) return;
		r=pow(2.0, (double ) if2-if3);
		cf2a[0]=cf2[0]*r; cf2a[1]=cf2[1]*r;

		c1[0]=cf3[0]-cf2a[0]; c1[1]=cf3[1]-cf2a[1];
		c3[0]=cx3[0]-cx2[0]; c3[1]=cx3[1]-cx2[1];
		cdiv(c1,c3,ch2);
		ch1a[0]=ch1[0]*r; ch1a[1]=ch1[1]*r;
		c1[0]=ch2[0]-ch1a[0]; c1[1]=ch2[1]-ch1a[1];
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

		if(cd[0]==0.0 && cd[1]==0.0) {*ier=433; return;}

		cdiv(cf3,cd,cdx); cdx[0]=-2.*cdx[0]; cdx[1]=-2.*cdx[1];
		cx1[0]=cx2[0]; cx1[1]=cx2[1];
		cx2[0]=cx3[0]; cx2[1]=cx3[1];
/*	The new approximation to zero */
		cx3[0]=cx3[0]+cdx[0]; cx3[1]=cx3[1]+cdx[1];
		cf1[0]=cf2[0]; cf1[1]=cf2[1]; if1=if2;
		cf2[0]=cf3[0]; cf2[1]=cf3[1]; if2=if3;
		ch1[0]=ch2[0]; ch1[1]=ch2[1];

		dx=cabs(cdx);
		r=reps*cabs(cx3); if(aeps>r) r=aeps;
		if(i>2 && dx<r) return;
		r=reps0*cabs(cx3); if(aeps>r) r=aeps;
		if(i>5 && dx<r && dx>dx1) {*ier=42; return;}

		dx1=dx;
		if(cabs(cx3)>rmax) {*ier=431; return;}

		if(i<nit) {
			*ier=-3;
			cx[0]=cx3[0]; cx[1]=cx3[1];
/*	To evaluate the function value at CX3 and continue the loop */
			return;
		}

		r=reps0*cabs(cx3); if(aeps>r) r=aeps;
		if(dx<r) {*ier=43; return;}
		else {*ier=432; return;}

	}
}




/*	Complex zeros of a given function using Muller iteration with deflation
	For use with MULER2 which requires function in form CF*2**IDET

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
		Function CF(Z,F,JF) must be supplied by the user, where F and Z
		are arrays of length 2 containing the complex values. F*2**JF should
		be the calculated value of the complex function at Z.
		
	Error status is returned by the value of the function ZROOT2.
		0 value implies successful execution
		41 implies that iteration did not converge to satisfactory accuracy for
			at least one zero, but value is acceptable at lower accuracy.
		429 implies that iteration failed to find at least one zero.
			The number of zeros found can be checked from the value of NZ.

	Required functions : MULER2, CF, CABS, CDIV, CSQRT
*/

#include <math.h>

void muler2(double cx1[], double cx2[], double cx3[], double reps, double aeps,
	int *ier, double cf[], double cx[], int ix, int nz, double czero[][2], double rmax);

int zroot2(int n, double cx[][2], double czero[][2], int *nz, double reps,
	double aeps, double rmax, void cf(double * , double * , int * ))

{
	int i,ier,ier1,ix;
	double cdx[2],cx1[2],cx2[2],cx3[2],cf0[2],cx0[2];

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
		ier1=0;
		do {
			muler2(cx1,cx2,cx3,reps,aeps,&ier1,cf0,cx0,ix,*nz,czero,rmax);
			if(ier1>=0) break;
			cf(cx0,cf0,&ix);
		} while(ier1<0);

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

