/*	Integration of tabulated function in 2 dimensions using B-splines */

#include <stdio.h>
#include <math.h>
int bspint2(int nx, int ny, double x[], double y[], double *f, int k, double *ax,
	double *ay, int la, double *c, double *xf, double *yf, int *mx, int *my,
	int *iflg, int intx[], int inty[]);
int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);
double bspqd2(int nx, int ny, double *x, double *y, int k, double *wt,
	int iw, double xl, double xu, double yl, double yu, int *ier);

main()
{
	int i,i1,j,k,ndim, np1,np2, n1,n2, iflg, if1, intx[20], inty[20];
	int m1,m2,nderiv,ier;
	double xl,h,xu, yl,yu, aeps, fb, ax[20][20],ay[20][20],cof[20][20];
	double xf[50],yf[50], dfxx, dfyy, dfxy, x1[20], x2[20],f[20][20];

	aeps=1.e-6;
	n1=n2=5;
	ndim=20;
	nderiv=0;
	iflg=0;
	printf("type  n1 = No. of points in each direction,  k = order of B-spline \n");
	i=scanf(" %d %d",&n1,&k);
	n2=n1;

/*	Setting up the abscissas */
	h=1.0/(n1-1.0);
	for(i=0; i<n1; ++i) {
		x1[i]=i*h;
		x2[i]=i*h;
	}

/*	Setting up the table of values using a known function */
	for(i=0; i<n1; ++i) {
		for(j=0; j<n2; ++j) f[j][i]=sin(x1[i])*sin(x2[j]);
	}
	i=bspint2(n1,n2,x1,x2,&f[0][0],k,&ax[0][0],&ay[0][0],ndim,&cof[0][0],xf,yf,
		   &m1,&m2,&iflg,intx,inty);
	printf(" ier = %d  k = %d   n1 = %d   n2 = %d   m1 = %d  m2 = %d\n", i,k,n1,n2,m1,m2);

	for(i1=0; i1<99; ++i1) {
		printf("type xl, xu, yl, yu  = limits for integration \n");
		printf("                       (quits when xl=xu)\n");
		scanf("%le %le %le %le", &xl, &xu, &yl, &yu);
		if(xl==xu) return 0;

		fb=bspqd2(m1,m2,xf,yf,k,&cof[0][0],ndim,xl,xu,yl,yu,&ier);
		printf(" ier = %d    x-limits = %e,%e   y-limits = %e %e\n", ier,xl,xu,yl,yu);
		printf(" integral = %e  \n", fb);

	}
	return;
}



/*	To calculate coefficients for B-spline interpolation

	N : (input) Number of entries in the table
	X : (input) Array of length N containing the abscissas
	F : (input) Array of length N containing the function values
		F[I] is the tabulated function value at X[I].
	K : (input) Order of B-spline required. K=4 gives cubic B-splines
	A : (input/output) Array of length LA*3K containing the
		triangular decomposition of equation matrix in band form
		For IFLG=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	LA : (input) The second dimension of A as specified in calling function
		LA>=N
	C : (output) Coefficients of expansion, which will be calculated
		provided IFLG!=1
	XF : (input/output) Array of size NO, containing
		the knots used for B-spline calculations.
		The knots must be distinct and in ascending order.
		For IFLG=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	NO : (input/output) Number of knots for B-splines
		For IFLG=2, this number must be supplied, for other values
		of IFLG it is calculated by the function
	IFLG : (input/output) Integer specifying the type of calculation required
		IFLG=0 The matrix will be calculated and solved for coefficients
		IFLG=1 The matrix will be calculated and triangular decomposition
			is obtained, but coefficients are not calculated
		IFLG=2 The triangular decomposition of matrix is assumed
			to be available in A and coefficients C are calculated
		IFLG=-1 same as 0, except that no pivoting will be used
	INC : (input/output) Integer array containing information about
		pivoting during solution of system of linear equations
		For IFLG=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
		
	Error status is returned by the value of the function BSPINT.
		0 value implies successful execution
		204 implies N<K or K<2
		other values may be set by BSPLIN or GAUBND

	Required functions : BSPLIN, GAUBND
*/

#include <math.h>
#include <stdlib.h>

int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[])

{
	int i,j,i1,i2,kl,ku,ndb,nderiv,kb,left,num,idet;
	double xb,det;
	double *wk;

	if(n<k || k<2) return 204;

	if(*iflg < 2) {
/*	set up the knots for B-splines by dropping points near the ends */
		xf[0] = x[0];
		kl=(k-2)/2;
		ku=(k-1)/2;
		for(i=1+kl; i<=n-2-ku; ++i) xf[i-kl] = x[i];
		xf[n-kl-ku-1] = x[n-1];
		*no=n-kl-ku;
		ndb=n+1;
		nderiv=0;
		wk=(double *) calloc((size_t) (2*n+5), sizeof(double));

/*	Set up the equation matrix for calculating coefficients of expansion
	The matrix is in band form A_{i,j} is stored in A[J-I+K-1][I] */
		for(i=0; i<n; ++i) {
			xb=x[i];
			j=bsplin(xf,*no,k,xb,nderiv,wk,(wk+ndb),(wk+ndb+2),&left);
			if(j>100) {free(wk); return j;}
			i1=i-k+1; if(i1<0) i1=0;
			i2=i+k-1; if(i2>n-1) i2=n-1;
			for(j=i1; j<=i2; ++j) a[la*(j-i+k-1)+i] = wk[j];
		}
		free(wk);
	}
	
/*	Solve the system of equations for a band matrix */
	num=1;
	kb=k-1;
	for(i=0; i<n; ++i) c[i]=f[i];
	i=gaubnd(n,kb,num,a,c,&det,&idet,inc,la,iflg);
	return i;
}



/*	To calculate coefficients for B-spline interpolation in 2 dimensions

	NX, NY : (input) Number of entries in the table along X, Y directions
	X, Y : (input) Array of length NX, NY containing the abscissas
	F : (input) Array of length LA*NY containing the function values
	K : (input) Order of B-spline required. K=4 gives cubic B-splines
	AX : (input/output) Array of length LA*3K containing the
		triangular decomposition of equation matrix in band form
		for interpolation along X.
		For IFLG>=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	AY : (input/output) Array of length LA*3K containing the
		triangular decomposition of equation matrix in band form
		for interpolation along Y.
		For IFLG>=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	LA : (input) the second dimension of arrays AX, AY, F, C as specified
			 in calling function, LA >= max(NX,NY)
	C : (output) Array of length LA*NY containing coefficients
		of expansion, which will be calculated by the function
	XF : (input/output) Array of size MX, containing
		the knots used for B-spline calculations along x.
		The knots must be distinct and in ascending order.
		For IFLG>1, this array must be supplied, for other values
		of IFLG it is calculated by the function
	YF : (input/output) Array of size MY, containing
		the knots used for B-spline calculations along y.
		The knots must be distinct and in ascending order.
		For IFLG>1, this array must be supplied, for other values
		of IFLG it is calculated by the function
	MX : (input/output) Number of knots for B-splines along X
		For IFLG>1, this value must be supplied, for other values
		of IFLG it is calculated by the function
	MY : (input/output) Number of knots for B-splines along Y
		For IFLG>1, this value must be supplied, for other values
		of IFLG it is calculated by the function
	IFLG : (input/output) Integer specifying the type of calculation required
		IFLG<2 The matrix will be calculated and solved for coefficients
		IFLG>1 The triangular decomposition of matrix is assumed
			to be available in AX, AY and coefficients C are calculated
	INTX, INTY : (input/output) Integer arrays of length NX, NY
		containing information about
		pivoting during solution of system of linear equations
		For IFLG>1, these arrays must be supplied, for other values
		of IFLG they are calculated by the function
		
	Error status is returned by the value of the function BSPINT2.
		0 value implies successful execution
		nonzero values may be set by BSPLIN, BSPINT or GAUBND

	Required functions : BSPINT, BSPLIN, GAUBND
*/

#include <math.h>
#include <stdlib.h>

int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);

int bspint2(int nx, int ny, double x[], double y[], double *f, int k, double *ax,
	double *ay, int la, double *c, double *xf, double *yf, int *mx, int *my,
	int *iflg, int intx[], int inty[])

{
	int i,j,kb,iflg1,idet;
	double det;
	double *wk;

	if(*iflg < 2) {
		iflg1=1;
		i=bspint(nx,x,f,k,ax,la,c,xf,mx,&iflg1,intx);
		if(i>100) return i;
		iflg1=1;
		i=bspint(ny,y,f,k,ay,la,c,yf,my,&iflg1,inty);
		if(i>100) return i;
		*iflg=2;
	}

	wk=(double *) calloc((size_t) (la*nx), sizeof(double));
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) wk[j+i*la]=f[i+j*la];
	}
	iflg1=2;
	kb=k-1;
	i=gaubnd(ny,kb,nx,ay,wk,&det,&idet,inty,la,&iflg1);
	if(i>100) {free(wk); return i;}
 
	for(j=0; j<ny; ++j) {
		for(i=0; i<nx; ++i) c[i+j*la]=wk[j+i*la];
	}
	iflg1=2;
	i=gaubnd(nx,kb,ny,ax,c,&det,&idet,intx,la,&iflg1);
	free(wk);
	return i;
}



/*	To calculate the B-spline basis functions at a specified point

	X : (input) Array of length NX containing the knots.
		The knots must be distinct and in ascending order.
	NX : (input) Number of knots
	K : (input) Order of B-spline, 0< K, K=4 gives cubic B-splines
	XB : (input) The point at which B-spline basis functions are to be evaluated
	NDERIV : (input) Number of derivatives required
		NDERIV<=0 only B-splines are calculated
		NDERIV=1 first derivative is also calculated
		NDERIV>1 first and second derivatives are also calculated
	B : (output) Array of length NX+K-2 containing the value of
		B-spline basis functions
	DB : (output) Array of length NX+K-2 containing the value of
		the first derivative of B-spline basis functions (if NDERIV>0)
	DDB : (output) Array of length NX+K-2 containing the value of
		the second derivative of B-spline basis functions (if NDERIV>1)
	LEFT : (output) XB is located between X[LEFT] and X[LEFT+1]
		
	Error status is returned by the value of the function BSPLIN.
		0 value implies successful execution
		26 implies XB > X[NX-1]
		27 implies XB < X[0]
		203 implies NX<2, K<1

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left)

{
	int i,j,igh, mid,nigh,lx,ier;
	double t1,t2,t3,p1,p2;
	double *wk, *dr, *dl;
	static int low = -1;

	if(nx <= 1 || k<1 ) return 203;
	ier=0;

/*	If the previous value of LOW is inadmissible, set the range to (0,N-1) */
	if(low<0 || low>=nx-1) {low=0; igh=nx-1;}
	else igh=low+1;

	while((xb<x[low] && xb<x[igh]) || (xb>x[low] && xb>x[igh])) {
/*	Extend the range */
		if( xb>x[low] ) {
/*	Extend the range on higher side */
			if(igh >= nx-1) {ier=26; low=nx-2; break;}
			else {
				nigh=igh+2*(igh-low); if(nx-1 < nigh) nigh=nx-1;
				low=igh; igh=nigh;
			}
		}

		else {
/*	Extend the range on lower side */
			if(low <= 0) {ier=27; igh=low+1; break;}
			else {
				nigh=low;
				low=low-2*(igh-low); if(low<0) low=0;
				igh=nigh;
			}
		}
	}


/*	Once the point is bracketed between two tabular points locate it by bisection */
	while((igh-low > 1) && (xb != x[low])) {
		mid=(low+igh)/2;
		if((xb<= x[mid]) == (xb<= x[low])) low=mid;
		else igh=mid;
	}

 
/*	Evaluate the B-spline basis functions

	Define the extra knots on either side of table
	Note that the function assumes knots from -K+2 to NX+K-1
	and the B-splines B_{i,k}, i ranges from 0 to NX+K-3 
	The knots are stored in scratch array wk. */

	wk=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));
	dr=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));
	dl=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));

	for(i=0; i<nx; ++i) wk[i+k]=x[i];
	for(i=1; i<=k; ++i) {
		wk[k-i]=x[0];
		wk[nx+i-1+k]=x[nx-1];
	}

	for(i=0; i<nx+k-2; ++i) {b[i]=0.0; db[i]=0.0; ddb[i]=0.0;}
	*left=low;
	lx=low-1;
	b[lx+1]=1;

/*	The recurrence relation for B-splines */
	for(j=1; j<=k-1; ++j) {
		dr[j] = wk[low+j+k] - xb;
		dl[j] = xb - wk[low+1-j+k];
		t1=0.0;
		for(i=1; i<=j; ++i) {
			t2=b[lx+i]/(dr[i]+dl[j+1-i]);
			b[lx+i]=t1+t2*dr[i];
			t1=t2*dl[j+1-i];
		}
		b[lx+j+1]=t1;
			
/*	Calculate the first derivative using recurrence relations */
		if(j == k-2 && nderiv > 0) {
			t1=0.0;
			for(i=1; i<=j+1; ++i) {
				t2=b[lx+i]/(wk[low+i+k]-wk[low+i+1]);
				db[lx+i]=(k-1)*(t1-t2);
				t1=t2;
			}
			db[lx+j+2]=(k-1)*t1;
		}
 
/*	Calculate the second derivative using recurrence relations */
		if(j == k-3 && nderiv>1) {
			t2=0.0; p1=0.0;
			for(i=1; i<=j+1; ++i) {
				t3=b[lx+i]/(wk[low+i+k]-wk[low+i+2]);
				p2=(t2-t3)/(wk[low+i+k]-wk[low+i+1]);
				ddb[lx+i]=(k-2)*(k-1)*(p1-p2);
				t2=t3; p1=p2;
			}
			p2=t2/(wk[low+j+2+k]-wk[low+j+3]);
			ddb[lx+j+2]=(k-2)*(k-1)*(p1-p2);
			ddb[lx+j+3]=(k-2)*(k-1)*p2;
		}
	}

/*	For K=2 the first derivative has to be calculated outside the loop */
	if(k == 2 && nderiv > 0) {
		t2=1./(wk[low+1+k]-wk[low+2]);
		db[lx+1]=-t2;
		db[lx+2]=t2;
	}

/*	For K=3 the second derivative has to be calculated outside the loop */
	if(k == 3 && nderiv > 1) {
		t3=1./(wk[low+1+k]-wk[low+3]);
		p2=-t3/(wk[low+1+k]-wk[low+2]);
		ddb[lx+1]=-2.*p2;
		p1=p2;
		p2=t3/(wk[low+2+k]-wk[low+3]);
		ddb[lx+2]=2.*(p1-p2);
		ddb[lx+3]=2.*p2;
	}
	free(dl); free(dr); free(wk);
	return ier;
}



/*	To calculate the integral of a function defined by B-spline expansion

	N : (input) Number of knots for B-spline expansion
	X : (input) Array of length N containing the knots.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Array of length N+K-2 containing the coefficients
		of B-spline expansion
	XL : (input) Lower limit of integration
	XU : (input) Upper limit of integration
	IER : (output) Error parameter, IER=0 implies successful execution
		IER=31 implies that XL is outside the range of table
		IER=32 implies that XU is outside the range of table

	BSPQD = Integral of \sum WT[I]\phi_I(x) over [XL,XU]

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
                double db[], double ddb[], int *left);

double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier)

{
	int i,j,k1,jk,lf1,lf2,nderiv;
	double f,f1,f2;
	double *wk1, *wk2, *wk;

	*ier=0;
	if(xl == xu) return 0.0;
 
/*	Calculate the B-spline basis of order K+1 at XL and XU */
	k1=k+1;
	nderiv=0;
	wk=(double *)calloc((size_t) (n+3*k), sizeof(double));
	wk1=(double *)calloc((size_t) (n+k), sizeof(double));
	wk2=(double *)calloc((size_t) (n+k), sizeof(double));

	*ier=bsplin(x,n,k1,xl,nderiv,wk1,wk,wk,&lf1);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	*ier=bsplin(x,n,k1,xu,nderiv,wk2,wk,wk,&lf2);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	if(xl<x[0] || xl>x[n-1]) *ier=31;
	if(xu<x[0] || xu>x[n-1]) *ier=32;

	for(i=0; i<n; ++i) wk[k+i]=x[i];
	for(i=1; i<=k; ++i) {
		wk[k-i]=x[0];
		wk[n+i-1+k]=x[n-1];
	}
 
/*	The sum for x=XL */
	f1=0.0;
	for(i=lf1; i<=lf1+k; ++i) {
		f=0.0;
		for(j=1; j<=i; ++j) f=f+wt[j-1]*(wk[k+j]-wk[j]);
		f1=f1+f*wk1[i]/k;
	}

/*	The sum for x=XU */
	f2=0.0;
	for(i=lf2; i<=lf2+k; ++i) {
		f=0.0;
		for(j=1; j<=i; ++j) f=f+wt[j-1]*(wk[k+j]-wk[j]);
		f2=f2+f*wk2[i]/k;
	}

	free(wk2); free(wk1); free(wk);
	return f2-f1;
}



/*	To calculate the integral of a function defined by B-spline expansion
	over 2 dimensions

	NX : (input) Number of knots for B-spline expansion along x direction
	NY : (input) Number of knots for B-spline expansion along y direction
	X : (input) Array of length NX containing the knots along x direction
		The knots must be distinct and in ascending order.
	Y : (input) Array of length NX containing the knots along y direction
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Array of length IW*(NY+K-2) containing coefficients
		of B-spline expansion. 
	IW : (input) The second dimension of WT as declared in calling function
	XL : (input) Lower limit of integration along x
	XU : (input) Upper limit of integration along x
	YL : (input) Lower limit of integration along y
	YU : (input) Upper limit of integration along y
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values may be set by BSPQD which is called

	BSPQD2 = Integral of \sum WT[J][I]\phi_I(x)\psi_J(y) over [XL,XU] x [YL,YU]

	Required functions : BSPLIN, BSPQD
*/

#include <math.h>
#include <stdlib.h>

double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);

double bspqd2(int nx, int ny, double *x, double *y, int k, double *wt,
	int iw, double xl, double xu, double yl, double yu, int *ier)

{
	int i;
	double ri;
	double *wk;

	wk=(double *) calloc((size_t) (ny+k), sizeof(double));

	for(i=0; i<ny+k-2; ++i) {
/*	Calculate integral along x */
		wk[i]=bspqd(nx,x,k,&wt[i*iw],xl,xu,ier);
		if(*ier>100) {free(wk); return 0.0;}
	}
		
/*	Calculate integral along y */
	ri=bspqd(ny,y,k,wk,yl,yu,ier);
	free(wk);
	return ri;
}



/*     Solution of a system of linear equations using Gaussian elimination
     	for a band matrix

	N : (input) Number of equations to be solved
	KB : (input) Bandwidth of matrix A[i][j]=0 if abs(I-J)>KB
	NUM : (input) Number of different sets (each with N equations) of
		equations to be solved
	A : (input/output) The matrix of coefficients of size LJ*(3*KB+1)
		A[J-I+KB][I] is the coefficient of x_J in Ith equation
		at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
		X[j][i] is the ith element of jth right hand side
		at output it will contain the solutions
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer variable used as a flag to specify the type
		of computation required
     	If IFLG=-1, both elimination and solution are calculated
     		without pivoting and IFLG is set to 2
	If IFLG=0, both elimination and solution are computed
     		with partial pivoting and IFLG is set to 2
     	If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
     	If IFLG>=2 only solution is calculated, the triangular
     		decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUBND.
		0 value implies successful execution
		104 implies N<=0 or N>LJ or KB>N
		124 implies some pivot turned out to be zero and hence
	     		matrix must be nearly singular

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg)

{
	int i,j,k,km,l,kb1,m1,m2;
	double r1, t1;
	double *wk;

	if(n<=0 || n>lj || kb>n) return 104;

	kb1=kb+1;
	if(*iflg < 2) {
/*     Perform elimination */
		for(i=0; i<n ; ++i) {
			for(j=2*kb+1; j<=3*kb; ++j) a[j*lj+i] = 0.0;
		}

		*det=1.0; *idet=0;
		wk=(double *) calloc((size_t) (3*kb+1), sizeof(double));
		for(k=0; k<n-1; ++k) {
/*     Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			if(*iflg >= 0) {
				m1=k+kb; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) {
					if(fabs(a[l+(k-l+kb)*lj]) > r1) {
						r1=fabs(a[l+(k-l+kb)*lj]);
						km=l;
					}
				}
			}

			inc[k]=km;
			if(km != k) {
/*     Interchange the rows if needed  */
				m1=2*kb+k; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) wk[l-k]=a[k+(l-k+kb)*lj];
				for(l=k; l<=m1; ++l) {
					a[k+(l-k+kb)*lj]=a[km+(l-km+kb)*lj];
					a[km+(l-km+kb)*lj]=wk[l-k];
				}
				*det= -(*det);
			}

			*det = (*det)*a[kb*lj+k];
			if( a[kb*lj+k] == 0.0) {free(wk); return 124;}
			if(*det != 0.0) {
/*     Scale the value of the determinant   */
				while(fabs(*det) > 32.0) {
					*det = (*det)*0.03125e0; *idet = *idet + 5;
				}

				while(fabs(*det) < 0.03125e0) {
					*det = (*det)*32.0; *idet = *idet - 5;
				}
			}

			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<= m1; ++l) {
				a[l+(k-l+kb)*lj] = a[l+(k-l+kb)*lj]/a[kb*lj+k];
				m2=k+2*kb; if(n-1 < m2) m2=n-1;
				for(i=k+1; i<=m2; ++i)
					a[l+(i-l+kb)*lj]=a[l+(i-l+kb)*lj]-a[l+(k-l+kb)*lj]*a[k+(i-k+kb)*lj];
			}
		}

		free(wk);
		*det = (*det)*a[(n-1)+kb*lj];
		inc[n-1]=n-1;
		if(a[(n-1)+kb*lj] == 0.0) return 124;

		if(*iflg==1) {*iflg=2; return 0;}
		*iflg=2;
	}
		
/*     Solution for the NUM different right-hand sides */
	for(j=0; j<num; ++j) {
/*     Forward substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k]) {
				t1= x[j*lj+k];
				x[j*lj+k] = x[j*lj+inc[k]];
				x[j*lj+inc[k]] = t1;
			}
			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<=m1; ++l)
				x[j*lj+l]=x[j*lj+l]-a[l+(k-l+kb)*lj]*x[j*lj+k];
		}

/*     back-substitution  */
		x[j*lj+n-1] = x[j*lj+n-1]/a[(n-1)+kb*lj];
		for(k=n-2; k>=0; --k) {
			m1=k+2*kb; if(n-1 < m1) m1=n-1;
			for(l=m1; l>=k+1; --l)
				x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k+(l-k+kb)*lj];
			x[j*lj+k] = x[j*lj+k]/a[kb*lj+k];
		}
	}
	return 0;
}

