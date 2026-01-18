/*	To integrate a tabulated function using cubic spline / B-spline */

#include <stdio.h>
#include <math.h>

int spline(double x[], double f[], int n, double c[][3]);
double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier);
int splint(double a, double b, double *sint, double *tint, int n, double x[],
	double f[], double c[][3]);
double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);

main()
{
	int i,i1,j,k,lj,nuse, ntab, iflg, ier,no, inc[20];
	double xl, xu, dfb, sint, tint, fi, ax[1000], c[20][3], x1[100],cf[30];

/*	Exercise 6.4 :  specified table of values */
	double x[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	double f[11] = {0.0, 0.3162, 0.4472, 0.5477, 0.6325, 0.7071,
		0.7746, 0.8367, 0.8944, 0.9487, 1.0};

	ntab=11;
	lj=20;
	i=spline(x, f,  ntab, c);
	printf(" ier = %d  ntab = %d   \n", i,ntab);

	printf("type k = order of B-spline\n");
	scanf(" %d",&k);
	iflg=0;
	i=bspint(ntab,x,f,k,ax,lj,cf,x1,&no,&iflg,inc);

	printf(" ier = %d  no. of knots = %d  \n", i,no);
	for(i1=0; i1<99; ++i1) {
		printf("type xl = lower limit,  xu = upper limit \n");
		printf("                    (quits when xl = xu)\n");
		scanf(" %le %le",  &xl, &xu);
		if(xl==xu) return 0;
/* The exact value of integral */
		fi=(pow(xu,1.5)-pow(xl,1.5))/1.5;

		i=splint(xl,xu,&sint,&tint, ntab, x, f, c);
		printf(" ier = %d   lower limit = %e   upper limit = %e \n", i,xl,xu);
		printf("  integral using  cubic spline = %e,   trapezoidal rule = %e\n",sint,tint);
		printf("  exact value = %e \n",fi);

		sint=bspqd(no,x1,k,cf,xl,xu,&ier);
		printf(" ier = %d   integral using B-spline = %e \n", ier,sint);

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




/*	To evaluate the cubic spline interpolant at a specified point

	XB : (input) point at which interpolation is required
	N : (input) Number of points in the table
	X : (input) Array of length N, containing the abscissas
	F : (input) Array of length N, containing the function values at X[I]
	C : (input) Array of length 3*N containing the spline coefficients
		which should have been calculated using SPLINE
	DFB : (output) First derivative of spline at x=XB
	DDFB : (output) Second derivative of spline at x=XB
	IER : (output) error parameter, IER=0 if execution is successful
		IER=24 implies XB is outside the range of table on higher side
		IER=25 implies XB is outside the range of table on lower side
		IER=201 implies N<2
	SPLEVL will be the interpolated value at x=XB

	Required functions : None
*/

#include <math.h>

double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier)

{
	int i,j,igh,nigh, mid;
	double dx, r1;
	static int low=-1;
	
	if(n<2) {*ier=201; return 0.0;}

	*ier=0;
/*	If the previous value of LOW is inadmissible, set the range to (0,N-1) */
	if(low<0 || low>=n-1) {low=0; igh=n-1;}
	else igh=low+1;

	while((xb<x[low] && xb<x[igh]) || (xb>x[low] && xb>x[igh])) {
/*	Extend the range */
		if((xb>x[low]) == (x[n-1]>x[0])) {
/*	Extend the range on higher side */
			if(igh >= n-1) {*ier=24; low=n-2; break;}
			else {
				nigh=igh+2*(igh-low); if(n-1 < nigh) nigh=n-1;
				low=igh; igh=nigh;
			}
		}

		else {
/*	Extend the range on lower side */
			if(low <= 0) {*ier=25; igh=low+1; break;}
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

	dx=xb-x[low];
	r1=((c[low][2]*dx+c[low][1])*dx+c[low][0])*dx+f[low];
	*dfb=(3.0*c[low][2]*dx+2.*c[low][1])*dx+c[low][0];
	*ddfb=6.*c[low][2]*dx+2.*c[low][1];
	return r1;
}



/*	To calculate coefficients of cubic spline interpolation with
		not-a-knot boundary conditions

	X : (input) Array of length N containing x values
	F : (input) Array of length N containing values of function at X[I]
		F[I] is the tabulated function value at X[I].
	N : (input) Length of table X, F
	C : (output) Array of length 3*N containing the spline coefficients
		
	Error status is returned by the value of the function SPLINE.
		0 value implies successful execution
		201 implies that N<2

	Required functions : None
*/

#include <math.h>

double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier);

int spline(double x[], double f[], int n, double c[][3])

{
	int i,j;
	double g, c1, cn, div12, div01;

	if(n<2) return 201;
	else if(n == 2) {
/*	Use linear interpolation */
		c[0][0]=(f[1]-f[0])/(x[1]-x[0]);
		c[0][1]=0.0;
		c[0][2]=0.0;
		return 0;
	}
	else if(n == 3) {
/*	Use quadratic interpolation */
		div01=(f[1]-f[0])/(x[1]-x[0]);
		div12=(f[2]-f[1])/(x[2]-x[1]);
		c[0][2]=0.0;
		c[1][2]=0.0;
		c[0][1]=(div12-div01)/(x[2]-x[0]);
		c[1][1]=c[0][1];
		c[0][0]=div01+c[0][1]*(x[0]-x[1]);
		c[1][0]=div12+c[0][1]*(x[1]-x[2]);
	        return 0;
	}
	else {
/*	Use cubic splines 

	Setting up the coefficients of tridiagonal matrix */
		c[n-1][2]=(f[n-1]-f[n-2])/(x[n-1]-x[n-2]);
		for(i=n-2; i>=1; --i) {
			c[i][2]=(f[i]-f[i-1])/(x[i]-x[i-1]);
			c[i][1]=2.*(x[i+1]-x[i-1]);
/*	The right hand sides */
			c[i][0]=3.*(c[i][2]*(x[i+1]-x[i])+c[i+1][2]*(x[i]-x[i-1]));
		}

/*	The not-a-knot boundary conditions */
		c1=x[2]-x[0];
		c[0][1]=x[2]-x[1];
		c[0][0]=c[1][2]*c[0][1]*(2.*c1+x[1]-x[0])+c[2][2]*(x[1]-x[0])*(x[1]-x[0]);
		c[0][0]=c[0][0]/c1;
		cn=x[n-1]-x[n-3];
		c[n-1][1]=x[n-2]-x[n-3];
		c[n-1][0]=c[n-1][2]*c[n-1][1]*(2.*cn+x[n-1]-x[n-2]);
		c[n-1][0]=(c[n-1][0]+c[n-2][2]*(x[n-1]-x[n-2])*(x[n-1]-x[n-2]))/cn;
/*	Solving the equation by Gaussian elimination */
		g=(x[2]-x[1])/c[0][1];
		c[1][1]=c[1][1]-g*c1;
		c[1][0]=c[1][0]-g*c[0][0];
		for(j=1; j<n-2; ++j) {
			g=(x[j+2]-x[j+1])/c[j][1];
			c[j+1][1]=c[j+1][1]-g*(x[j]-x[j-1]);
			c[j+1][0]=c[j+1][0]-g*c[j][0];
		}
		g=cn/c[n-2][1];
		c[n-1][1]=c[n-1][1]-g*(x[n-2]-x[n-3]);
		c[n-1][0]=c[n-1][0]-g*c[n-2][0];


/*	The back-substitution */
		c[n-1][0]=c[n-1][0]/c[n-1][1];
		for(i=n-2; i>=1; --i) c[i][0]=(c[i][0]-c[i+1][0]*(x[i]-x[i-1]))/c[i][1];
		c[0][0]=(c[0][0]-c[1][0]*c1)/c[0][1];

/*	Calculating the coefficients of cubic spline */
		for(i=0; i<n-1; ++i) {
			c[i][1]=(3.*c[i+1][2]-2.*c[i][0]-c[i+1][0])/(x[i+1]-x[i]);
			c[i][2]=(c[i][0]+c[i+1][0]-2.*c[i+1][2])/((x[i+1]-x[i])*(x[i+1]-x[i]));
		}
/*	Set the coefficients for interval beyond X(N) using continuity
	of second derivative, although they may not be used. */
		c[n-1][1]=c[n-1][1]+3*(x[n-1]-x[n-2])*c[n-2][2];
		c[n-1][2]=0.0;
		return 0;
	}
}



/*	To compute integral of a tabulated function using cubic splines

	A : (input) Lower limit of integration
	B : (input) Upper limit of integration (B > A)
	SINT : (output) Computed value of integral using cubic splines
	TINT : (output) Computed value of integral using trapezoidal rule
	N : (input) Number of tabular points
	X : (input) Array of length N containing abscissas (in ascending order)
	F : (input) Array of length N containing the function values
	C : (input) Array of length 3*N containing the coefficients of splines
		
	Error status is returned by the value of the function SPLINT.
		0 value implies successful execution
		31 implies lower limit A is outside the table
		32 implies upper limit B is outside the table
		301 implies A>B or X[0]>X[N-1] and no calculations are done
		Other values may be set by SPLEVL

	Required functions : SPLEVL
*/

#include <math.h>


double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier);

int splint(double a, double b, double *sint, double *tint, int n, double x[],
	double f[], double c[][3])

{
	int i,ier,ier1;
	double a1,f1,f2,dfb,ddfb,b1,d1,d2;

	*sint=0.0; *tint=0.0;
	if(a == b) return 0;
	if(a>b || x[0]>=x[n-1]) return 301;

	a1=a;
/*	Evaluate the function at x=A for trapezoidal rule */
	f1=splevl(a1,n,x,f,c,&dfb,&ddfb,&ier);
	if(ier>100) return ier;
	if(a<x[0] || a>x[n-1]) ier=31;
	if(b<x[0] || b>x[n-1]) ier=32;

/*	Integrating over the n-1 subintervals */
	for(i=0; i<n-1; ++i) {
		if(a1 < x[i+1] || i == n-2) {
			b1=x[i+1]; if(b<b1) b1=b;
			if(i == n-2) b1=b;
			d1=a1-x[i];
			d2=b1-x[i];

/*	Add integral of cubic spline over [A1,B1] */
			*sint = *sint+f[i]*(d2-d1)+c[i][0]*(d2*d2-d1*d1)/2.+
				c[i][1]*(d2*d2*d2-d1*d1*d1)/3.+
				c[i][2]*(d2*d2*d2*d2-d1*d1*d1*d1)/4.;
			f2=f[i+1];
			if(b1 != x[i+1]) f2=splevl(b1,n,x,f,c,&dfb,&ddfb,&ier1);

/*	Trapezoidal rule approximation to integral */
			*tint = *tint+0.5*(b1-a1)*(f1+f2);
			if(b<=x[i+1]) return ier;
			a1=b1;
			f1=f2;
		}
	}
	return ier;
}
