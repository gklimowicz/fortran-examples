/*	Abscissas and weights of Gaussian quadrature formulas with a
	specified weight function */

#include <stdio.h>
#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int polyr(int n, double a[], double cr[], double ci[], int qrefin);

int gauswt(int n, double w[], double ab[], double (*fmom) (int ), int qgaus);

double fmom(int i);
double f2(double x);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np;
	double xl, xu, rint, reps, aeps, df, a[50], w[50];

/*	Exercise 6.17  */

	aeps=1.e-18; reps=1.e-13;
	for(i1=0; i1<99; ++i1) {
		printf("type  nuse = No. of points in quadrature formula,  iflg=0/1\n");
		printf("   iflg =0/1  for Newton-Cotes / Gaussian formula \n");
		printf("                       (quits when  nuse <= 0)\n");
		scanf(" %d %d", &nuse, &iflg);
		if(nuse<=0) return 0;

		for(i=0; i<nuse; ++i) a[i]=i/(nuse-1.0);
		i=gauswt(nuse,w,a,fmom,iflg);
		printf(" ier = %d  n = %d  iflg = %d  \n", i,nuse,iflg);
		printf("   Abscissas             Weights  \n");
		for(i=0; i<nuse; ++i) printf(" %20.13e %20.13e \n",a[i],w[i]);
	}
	return;
}

/*	to calculate integrals of  w(x)*x**N over the required interval */
double fmom(int i)

{
	return 1.0/((1.0+i)*(1.0+i));
}



/*	Solution of a system of linear equations using Gaussian elimination
	with partial pivoting

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input/output) The matrix of coefficient of size LJ*N
	        A[i][j] is the coefficient of x_j in ith equation
	     	at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
	     	at output it will contain the solutions
	DET : (output) The determinant of the matrix
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter to specify the type of computation required
		If IFLG<=0, both elimination and solution are
			done and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUELM.
		0 value implies successful execution
		101 implies (N<=0 or N>LJ) 
		121 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg)

{
int i,j,k,km,l;
double r1,t1;

	if(n<=0 || n>lj) return 101;
 
	if((*iflg)<2) {
/*	Perform elimination  */
		*det=1.0;
		for(k=0; k<n-1; ++k) {
/*	Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			for(l=k; l<n; ++l)
				if(fabs(a[l*lj+k])>r1) {r1=fabs(a[l*lj+k]); km=l;}

			inc[k]=km;
			if(km != k) {
/*	Interchange the rows if needed  */
				for(l=k; l<n; ++l) 
				{t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;}
				*det=-(*det);
			}

			*det=(*det)*a[k*lj+k];
			if(a[k*lj+k]==0) return 121;

			for(l=k+1; l<n; ++l) {
				a[l*lj+k]=a[l*lj+k]/a[k*lj+k];
				for(i=k+1; i<n; ++i) a[l*lj+i]=a[l*lj+i]-a[l*lj+k]*a[k*lj+i];
			}
		}
		*det=(*det)*a[(n-1)*lj+n-1];
		inc[n-1]=n-1;
		if(a[(n-1)*lj+n-1]==0) return 121;

		if((*iflg)==1) {*iflg=2; return 0;}
		*iflg=2;
	}

/*	Solution for the num different right-hand sides  */
	for(j=0; j<num; ++j) {
/*	forward-substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k])
			{t1=x[j*lj+k]; x[j*lj+k]=x[j*lj+inc[k]]; x[j*lj+inc[k]]=t1;}
			for(l=k+1; l<n; ++l) x[j*lj+l]=x[j*lj+l]-a[l*lj+k]*x[j*lj+k];
		}

/*	back-substitution  */

		x[j*lj+n-1]=x[j*lj+n-1]/a[(n-1)*lj+n-1];
		for(k=n-2; k>=0; --k) {
			for(l=n-1; l>=k+1; --l) x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=x[j*lj+k]/a[k*lj+k];
		}
	}
	return 0;
}



/*	To calculate weights and abscissas of a quadrature formula with
	specified weight function.

	N : (input) Number of points in the required quadrature formula
	W : (output) Array of length N, which will contain the weights
	AB : (input/output) Array of length N containing the abscissas
		For Gaussian formulas (QGAUS=1) AB will be calculated, 
		while for other values of QGAUS abscissas must be supplied
	FMOM : (input) Name of the function to calculate the moments
		Function FMOM(I) should calculate integral of w(x)x**I
	QGAUS : (input) Parameter to decide type of formula to be obtained
		If QGAUS=1 a Gaussian formula is calculated. In this
		case both abscissas and weights are calculated.
		Otherwise an interpolatory formula is calculated.
		In this case only weights are calculated, while abscissas
		must be supplied.
		
	Error status is returned by the value of the function GAUSWT.
		0 value implies successful execution
		303 implies N<=0 or N>=NPMAX
		322 implies GAUELM failed to find coefficients of polynomial
		323 implies POLYR failed to find roots of polynomial
		324 implies GAUELM failed to find weights 

	Required functions : GAUELM, POLYR, LAGITR, FMOM, CABS, CSQRT, CDIV
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int polyr(int n, double a[], double cr[], double ci[], int qrefin);

int gauswt(int n, double w[], double ab[], double (*fmom) (int ), int qgaus)

{
	int i,j,k,lj,iflg,ier,qrefin;
	double det;
	double *a, *cof, *zero;
	int *inc;

	if(n <= 0) return 303;
	lj=n;
	a=(double *)calloc((size_t) (n*n), sizeof(double));
	inc=(int *) calloc((size_t) n, sizeof(int));

	if(qgaus == 1) {
		cof=(double *)calloc((size_t) (n+2), sizeof(double));
		zero=(double *)calloc((size_t) (2*n), sizeof(double));

/*	Calculating the coefficients of the orthogonal polynomial */
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) a[i*n+j]=fmom(i+j);
			cof[i]=-fmom(n+i);
		}
		iflg=0;
		ier=gauelm(n,1,a,cof,&det,inc,lj,&iflg);
		if(ier >100) {free(zero); free(cof); free(inc); free(a); return 322;}

/*	Find the roots of polynomial, which will be the abscissas */
		cof[n]=1.0;
		qrefin=1;
		polyr(n,cof,zero,&zero[n],qrefin);
		if(ier >100) {free(zero); free(cof); free(inc); free(a); return 323;}
		for(i=0; i<n; ++i) ab[i]=zero[i];
		free(zero); free(cof);
	}


/*	Calculate the weights */
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) {
			if(i==0) a[i*n+j]=1.0;
			else a[i*n+j]=ab[j]*a[i*n+j-n];
		}
		w[i]=fmom(i);
	}
	iflg=0;
	ier=gauelm(n,1,a,w,&det,inc,lj,&iflg);
	free(inc); free(a);
	if(ier >100) return 324;
	return 0;
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

