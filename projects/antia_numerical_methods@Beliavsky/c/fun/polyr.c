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
