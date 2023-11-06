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
