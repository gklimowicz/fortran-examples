/*	Root of a polynomial with real coefficients using Laguerre iteration

	N : (input) The degree of the polynomial
	A : (input) Array of length N+1 containing the coefficients of
		the polynomial. A[0] is the constant term and A[N] is the
		coefficient of X**N
	CXI : (input/output) Array of length 2 containing the real and
		imaginary parts of the initial guess,
		after execution it will contain the computed root
		
	Error status is returned by the value of the function LAGITR.
		0 value implies successful execution
		438 implies that denominator is zero and iteration cannot
			be continued further
		439 implies that iteration has failed to converge

	Required functions : CABS, CSQRT, CDIV
*/

#include <math.h>

double cabs(double cx[2]);
void csqrt(double cx[2], double cs[2]);
void cdiv(double c1[2], double c2[2], double c[2]);

int lagitr(int n, double a[], double cxi[2])

{
	int ic,i,j,qc, itmax=50;
	double r,cdx1,cf[2],cfp[2],cfpp[2],ch[2],cden[2],cdx[2],cr[2];
	double reps=1.0e-4, aeps=1.0e-6;

	cdx1=cabs(cxi)+1.0;
	qc=0;
	ic=itmax;

/*	The Laguerre iteration */
	for(i=1; i<=itmax; ++i) {
/*	Calculate the polynomial and its derivatives */
		cf[0]=a[n]; cf[1]=0.0;
		cfp[0]=0.0; cfp[1]=0.0;
		cfpp[0]=0.0; cfpp[1]=0.0;
		for(j=n-1; j>=0; --j) {
			r=2.*cfp[0]+cxi[0]*cfpp[0]-cxi[1]*cfpp[1];
			cfpp[1]=2.*cfp[1]+cxi[0]*cfpp[1]+cxi[1]*cfpp[0];
			cfpp[0]=r;
			r=cf[0]+cxi[0]*cfp[0]-cxi[1]*cfp[1];
			cfp[1]=cf[1]+cxi[0]*cfp[1]+cxi[1]*cfp[0];
			cfp[0]=r;
			r=a[j]+cxi[0]*cf[0]-cxi[1]*cf[1];
			cf[1]=cxi[0]*cf[1]+cxi[1]*cf[0];
			cf[0]=r;
		}

		cr[0]=(n-1)*((n-1)*(cfp[0]*cfp[0]-cfp[1]*cfp[1])-
				n*(cf[0]*cfpp[0]-cf[1]*cfpp[1]));
		cr[1]=(n-1)*((n-1)*2.0*cfp[0]*cfp[1]-n*(cf[0]*cfpp[1]+cf[1]*cfpp[0]));
		csqrt(cr,ch);
		cden[0]=cfp[0]+ch[0]; cden[1]=cfp[1]+ch[1];
		cr[0]=cfp[0]-ch[0]; cr[1]=cfp[1]-ch[1];
		if(cabs(cr)>cabs(cden)) {cden[0]=cr[0]; cden[1]=cr[1];}

		if(cden[0] != 0 || cden[1] != 0) {
/*	Laguerre's iteration */
			cdiv(cf,cden,cr); cdx[0]=-n*cr[0]; cdx[1]=-n*cr[1];
			r=reps*cabs(cxi); if(aeps>r) r=aeps;
			if(cabs(cdx)<r && i>1 && qc==0) {qc=1; ic=i;}
			if(qc==1 && cabs(cdx)/cdx1>0.99) return 0;
			cdx1=cabs(cdx);
			if(cdx1==0.0) return 0;
			cxi[0]=cxi[0]+cdx[0]; cxi[1]=cxi[1]+cdx[1];
		}
		else {
			if(cf[0]==0.0 && cf[1]==0.0) return 0;
			else return 438;
		}
	}

	if(qc==1) return 0;
	else return 439;
}
