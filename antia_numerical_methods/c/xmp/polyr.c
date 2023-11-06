/*	Roots of a polynomial with real coefficients using Laguerre's method */

#include <stdio.h>
#include <math.h>

int lagitr(int n, double a[], double cxi[2]);
int polyr(int n, double a[], double rr[], double ri[], int qrefin);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double xl, xu, rr[20], ri[20], x0, a[20];

	iflg=0;
	for(i1=0; i1<99; ++i1) {
		printf("type np = degree of polynomial     (quits when np<=0)\n");
		scanf(" %d", &np);
		if(np<=0) return 0;

		printf("type the coefficients starting with highest degree term\n");
		for(i=np; i>=0; --i) scanf(" %le",&a[i]);

		i=polyr(np,a,rr,ri,iflg);
		printf(" ier = %d  np = %d    roots: \n", i,np);
		for(i=0; i<np; ++i) printf(" %e %e \n",rr[i],ri[i]);
		printf("\n");

	}
	return;
}




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



/*	Roots of a polynomial with real coefficients using Laguerre iteration

	N : (input) The degree of the polynomial
	A : (input) Array of length N+1 containing the coefficients of
		the polynomial. A[0] is the constant term and A[N] is the
		coefficient of X**N
	RR : (output) Array of length N, containing the real part of computed roots
	RI : (output) Array of length N, containing the imaginary part of computed roots
	QREFIN : (input) Parameter to decide if roots need to be refined
		If QREFIN=1 the roots are refined using original polynomial
		otherwise no refinement is done.
		
	Error status is returned by the value of the function POLYR.
		0 value implies successful execution
		k*11 implies that iteration for refining the roots failed
			to converge for k of the roots
		406 implies that N<=0 and no calculations are done
		408 implies that A[N]=0 and no calculations are done
		430 implies that iteration failed to converge for some root

	Required functions : LAGITR, CABS, CSQRT, CDIV
*/

#include <math.h>
#include <stdlib.h>

int lagitr(int n, double a[], double cxi[2]);

int polyr(int n, double a[], double rr[], double ri[], int qrefin)

{
	int i,ier,i1,j,np;
	double cxr[2],xrt,bn,an,xr,xi,p,s,bn1, hcross=1.e-16;
	double *wk;

	if(n<=0) return 406;
	if(a[n]==0.0) return 408;

	wk=(double *) calloc((size_t) (n+1), sizeof(double));
/*	Make a copy of coefficients to be used for deflation */
	for(i=0; i<=n; ++i) wk[i]=a[i];
	np=n;
/*	Starting value for iteration to allow the roots to be found
	in ascending order of magnitude */
	cxr[0]=0.0; cxr[1]=0.0;

	do {
		ier=lagitr(np,wk,cxr);
		if(ier != 0) {
/*	If the iteration fails to converge, then try once more */
			cxr[0]=1.123456; cxr[1]=0.0;
			ier=lagitr(np,wk,cxr);
			if(ier != 0) {free(wk); return 430;}
		}

		if(fabs(cxr[1]) <=10.0*hcross*fabs(cxr[0])) {
/*	Perform deflation for a real root */
			xrt=cxr[0];
			if(np>1) {
				bn=wk[np-1];
				wk[np-1]=wk[np];
				for(i=np-2; i>=0; --i) {
					an=xrt*wk[i+1]+bn;
					bn=wk[i];
					wk[i]=an;
				}
			}
			rr[n-np]=cxr[0]; ri[n-np]=cxr[1];
			np=np-1;
		}

		else {
/*	Perform deflation for a pair of complex conjugate roots */
			xr=cxr[0]; xi=fabs(cxr[1]);
			if(np>2) {
				p=2.*xr;
				s=xr*xr+xi*xi;
				bn=wk[np-2];
				bn1=wk[np-3];
				wk[np-2]=wk[np];
				wk[np-3]=wk[np-1]+wk[np-2]*p;
				for(i=np-4; i>=0; --i) {
					an=bn+p*wk[i+1]-wk[i+2]*s;
					bn=bn1;
					bn1=wk[i];
					wk[i]=an;
				}
			}
			np=np-2;
			rr[n-np-2]=xr; ri[n-np-2]=xi;
			rr[n-np-1]=xr; ri[n-np-1]=-xi;
		}
	} while(np>0);
	free(wk);

	ier=0;
	if(qrefin==1) {
/*	Refine the roots by iterating on original polynomial */
		for(i=1; i<n; ++i) {
			cxr[0]=rr[i]; cxr[1]=ri[i];
			ier=lagitr(n,a,cxr);
			if(ier==0) {rr[i]=cxr[0]; ri[i]=cxr[1];}
			else ier=ier+11;
		}
	}

/*	Sort the roots in ascending order of real part */
	for(i=1; i<n; ++i) {
		cxr[0]=rr[i]; cxr[1]=ri[i];
		for(j=i-1; j>=0; --j) {
			if(rr[j]<=cxr[0]) break;
			rr[j+1]=rr[j]; ri[j+1]=ri[j];
		}
		if(j<-1) j=-1;
		rr[j+1]=cxr[0]; ri[j+1]=cxr[1];
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

