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
