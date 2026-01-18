/*	Polynomial interpolation in two dimensions */

#include <stdio.h>
#include <math.h>

int nearst(double xb, double x[], int ntab);
int divdif0(double xb, double x[], double f[], int *nuse, int ntab,
		double fb[], double aeps, int iflg, int *if1);
int poly2(double xb1, double xb2, double x1[], double x2[], double *f, int ndim,
	int n1, int n2, int np1, int np2, double *fb);
int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);

int bspint2(int nx, int ny, double x[], double y[], double *f, int k, double *ax,
	double *ay, int la, double *c, double *xf, double *yf, int *mx, int *my,
	int *iflg, int intx[], int inty[]);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspev2(int nx, int ny, double *x, double *y, int k, int nderiv,
	double *wt, int iw, double x0, double y0, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy, int *ier);

main()
{
	int i,i1,j,ndim,k,nx,ny, m1,m2, n1,n2, iflg, if1, la,intx[20],inty[20];
	double xb1, xb2, dfx, dfy, aeps, fb, ax[20][5], ay[20][5], b[20][5];
	double xf[20],yf[20], dfxx, dfxy, dfyy;

/*	Exercise 4.31 :  specified table of values */
	double x1[5] = {50.0, 60.0, 70.0, 80.0, 90.0};
	double x2[5] = {50.0, 52.0, 54.0, 56.0, 58.0};
	double f[5][5] = {0.9401,0.9647,0.9876,1.0044,1.0107,
		0.9835,1.0118,1.0387,1.0587,1.0662,
		1.0277,1.0602,1.0915,1.1152,1.1242,
		1.0725,1.1097,1.1462,1.1743,1.1851,
		1.1180,1.1605,1.2030,1.2362,1.2492};

	aeps=1.e-6;
	n1=n2=5;
	ndim=5;
	printf("type k = order of B-spline \n");
	scanf("%d", &k);
	iflg=0;
	la=5;

/*	Calculate the coefficients of B-spline expansion */
	i=bspint2(n1,n2,x1,x2,&f[0][0],k,&ax[0][0],&ay[0][0],la,&b[0][0],xf,yf,&nx,&ny,&iflg,intx,inty);
	printf("ier = %d   k = %d   No. of knots = %d %d \n",i,k,nx,ny);
	printf(" Coefficients : \n");
	for(i=0; i<n2; ++i){
		for(j=0; j<n1; ++j) printf(" %e ",b[i][j]);
		printf("\n");
	}

	for(i1=0; i1<99; ++i1) {
		printf("type m1, m2 = No. of points in x,y\n");
		printf("            (quits when m1<=0  or m2<=0)\n");
		scanf("%d %d", &m1, &m2);
		if(m1<=0 || m2<=0) return 0;
		printf("type xb1, xb2 = coordinates of point where interpolation is required \n");
		scanf("%le %le", &xb1, &xb2);

/*	polynomial interpolation in 2 dimensions */
		i=poly2(xb1, xb2, x1, x2, &f[0][0], ndim, n1, n2, m1, m2, &fb);
		printf(" ier = %d  m1 = %d, m2 = %d,  x = %e,%e  f = %e\n", i,m1,m2,xb1,xb2,fb);

/* Evaluate the B-spline interpolation in 2 dimensions */
		if1=2;
		fb=bspev2(nx,ny,xf,yf,k,if1,&b[0][0],la,xb1,xb2,&dfx,&dfy,&dfxx,&dfxy,&dfyy,&i);
		printf(" B-spline interpolation :  f = %e   f' = %e  %e \n   f'' = %e  %e  %e \n",
				fb,dfx,dfy,dfxx,dfxy,dfyy);

	}
	return;
}

 

/*	To evaluate a B-spline expansion in 2 dimensions

	NX : (input) Number of knots to define B-splines along 1st dimension
	NY : (input) Number of knots to define B-splines along 2nd dimension
	X : (input) Array of length NX containing the knots.
		The knots must be distinct and in ascending order.
	Y : (input) Array of length NY containing the knots.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	NDERIV : (input) Number of derivatives required
		For NDERIV<=0 only function value is calculated
		For NDERIV=1 first derivative is also calculated
		For NDERIV>1 both first and second derivatives are calculated
	WT : (input) Array of length IW*(NY+K-2) containing coefficients
		of B-spline expansion,
	IW : (input) Second dimension of array WT as defined in the calling function
	X0,Y0 : (input) The point at which expansion has to be evaluated
	DFX : (output) First derivative of function w.r.t. X at X0, Y0
	DFY : (output) First derivative of function w.r.t. Y at X0, Y0
	DFXX : (output) Second derivative of function w.r.t X,X at X0, Y0
	DFXY : (output) Second derivative of function w.r.t X,Y at X0, Y0
	DFYY : (output) Second derivative of function w.r.t Y,Y at X0, Y0
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values of IER may be set by BSPLIN which is called

	BSPEV2 =SUM_{i=1}^{NX+K-2} SUM_{j=1}^{NY+K-2} WT[j][i]\phi_i(X0)\psi_j(Y0)
	where \phi_i(x) are B-spline basis functions on knots X
	and \psi_j(y) are B-spline basis functions on knots Y

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

double bspev2(int nx, int ny, double *x, double *y, int k, int nderiv,
	double *wt, int iw, double x0, double y0, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy, int *ier)

{
	int i,j,lx,ly;
	double f;
	double *wx, *wx1, *wx2, *wy, *wy1, *wy2;

	wx=(double *) calloc((size_t) (nx+k),sizeof(double));
	wy=(double *) calloc((size_t) (ny+k),sizeof(double));
	wx1=(double *) calloc((size_t) (nx+k),sizeof(double));
	wy1=(double *) calloc((size_t) (ny+k),sizeof(double));
	wx2=(double *) calloc((size_t) (nx+k),sizeof(double));
	wy2=(double *) calloc((size_t) (ny+k),sizeof(double));
	*ier=bsplin(x,nx,k,x0,nderiv,wx,wx1,wx2,&lx);
	if(*ier > 100) {
		free(wy2); free(wx2); free(wy1); free(wx1); free(wy); free(wx);
		return 0.0;
	}
	*ier=bsplin(y,ny,k,y0,nderiv,wy,wy1,wy2,&ly);
	if(*ier > 100) {
		free(wy2); free(wx2); free(wy1); free(wx1); free(wy); free(wx);
		return 0.0;
	}

	f=0.0; *dfx=0.0; *dfy=0.0;
	*dfxx=0.0; *dfyy=0.0; *dfxy=0.0;

	for(i=lx; i<=lx+k-1; ++i) {
		for(j=ly; j<=ly+k-1; ++j) {
			f=f+wt[i+j*iw]*wx[i]*wy[j];
			*dfx = *dfx + wt[i+j*iw]*wx1[i]*wy[j];
			*dfy = *dfy + wt[i+j*iw]*wx[i]*wy1[j];
			*dfxx = *dfxx + wt[i+j*iw]*wx2[i]*wy[j];
			*dfyy = *dfyy + wt[i+j*iw]*wx[i]*wy2[j];
			*dfxy = *dfxy + wt[i+j*iw]*wx1[i]*wy1[j];
		}
	}
	free(wy2); free(wx2); free(wy1); free(wx1); free(wy); free(wx);
	return f;
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


/*	Interpolation using Newton's divided difference formula
	simplified version of DIVDIF without derivative calculation

	XB : (input) Value of x at which interpolation is required
	X : (input) Array of length NTAB containing x values
	F : (input) Array of length NTAB containing function values
		F[I] is the tabulated function value at X[I].
	NUSE : (input/output) Number of points to be used for interpolation
		After execution it will contain the number actually used
	NTAB : (input) Number of points in the table
	FB : (output) Array containing interpolated values
	    FB[I] should contain interpolation using I points
	    FB[NUSE] should be the final value
	AEPS : (input) Required accuracy
	IFLG : (input) Flag to decide whether nearest point has to be found
		If IFLG=0 find the nearest point to XB to start interpolation
		otherwise if IF1 is admissible use X[IF1] as the first point
	IF1 : (input/output) The first point to be used for interpolation
		when IFLG!=0
		If IFLG=0 then IF1 is set to the index of nearest point in X
		
	Error status is returned by the value of the function DIVDIF0.
		0 value implies successful execution
		21 implies NUSE<1, in which case it is set to MIN(6,NTAB)
		22 implies NUSE>NTAB or NMAX, in which case it is reduced
		23 implies interpolation has not converged to specified accuracy

	Required functions : NEARST
*/

#include <math.h>

int nearst(double xb, double x[], int ntab);

int divdif0(double xb, double x[], double f[], int *nuse, int ntab,
		double fb[], double aeps, int iflg, int *if1)


{
	int i,j,k,next,in,ip,nit,ier, nmax=10;
	double err,px,xn[11],xd[11];

/*	Find the nearest point */

	if(iflg == 0 || *if1 < 0 || *if1 >= ntab) {
		next=nearst(xb,x,ntab);
		*if1=next;
	}
	else next=*if1;

	fb[1]=f[next];
	xd[1]=f[next];
	xn[1]=x[next];
	ier=0;
	px=1.0;

/*	Points between IN and IP are used for interpolation */

	ip=next; in=next;

/*	Maximum number of points to be used for interpolation */
	nit=*nuse; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
	if(*nuse>nmax || *nuse>ntab) ier=22;
	if(*nuse<1) {
		ier=21;
		nit=6; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
	}
	*nuse=1;
		
/*	Calculate successive interpolation polynomial */
	for(j=2; j<=nit; ++j) {

/*	Choose the next nearest point to XB */
		if(in<=0 ) {
			ip=ip+1; next=ip;
		}
		else if(ip >= ntab-1) {
			in=in-1; next=in;
		}
		else if(fabs(xb-x[ip+1]) < fabs(xb-x[in-1]) ) {
			ip=ip+1; next=ip;
		}
		else {
			in=in-1; next=in;
		}

/*	Calculating the divided differences */
		xd[j]=f[next];
		xn[j]=x[next];
		for(k=j-1; k>=1; --k) xd[k]=(xd[k+1]-xd[k])/(xn[j]-xn[k]);

		px=px*(xb-xn[j-1]);
		err=xd[1]*px;
		fb[j]=fb[j-1]+err;
		*nuse=j;

		if(fabs(err) < aeps) return ier;
	}
	return 23;
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





/*	To locate the nearest point in an ordered table using bisection

	XB : (input) given value of x for which nearest point is needed
	X : (input) array of length NTAB containing table of values
	NTAB : (input) length of table
	After execution X[NEARST] is the tabular point closest to XB 

	Required functions : None
*/

#include <math.h>

int nearst(double xb, double x[], int ntab)

{
	int low, igh, mid;

	low=0; igh=ntab-1;
	if((xb < x[low]) != (xb < x[igh]) ) {

/*	If the point is within the range of table, then locate it by bisection */

		while(igh-low > 1) {
			mid=(low+igh)/2;
			if((xb < x[mid]) == (xb < x[low])) low=mid;
			else igh=mid;
		}
	}

	if(fabs(xb-x[low]) < fabs(xb-x[igh])) return low;
	else return igh;
}
				
				


/*	To calculate polynomial interpolation in 2 dimensions on a
	rectangular mesh of points

	(XB1,XB2) : (input) is the point at which interpolation is required
	X1 : (input) Array of length N1 containing the abscissas
	X2 : (input) Array of length N2 containing the abscissas
	F : (input) Array of length NDIM*N2 containing the function values
		F[J][I]=f(X1(I),X2(J))
	NDIM : (input) Second dimension of array F as specified in calling function
	N1 : (input) Length of array X1, i.e Number of points along first dimension
	N2 : (input) Length of array X2, i.e Number of points along second dimension
	NP1 : (input) Number of points to be used for interpolation along X1
	NP2 : (input) Number of points to be used for interpolation along X2
	FB : (output) Interpolated value of function at (XB1,XB2)
		
	Error status is returned by the value of the function POLY2.
		0 value implies successful execution
		206 implies N1 > NDIM, in which case no calculations
			are done

	Required functions : DIVDIF0, NEARST
*/

#include <math.h>

int nearst(double xb, double x[], int ntab);
int divdif0(double xb, double x[], double f[], int *nuse, int ntab,
		double fb[], double aeps, int iflg, int *if1);

int poly2(double xb1, double xb2, double x1[], double x2[], double *f, int ndim,
	int n1, int n2, int np1, int np2, double *fb)

{
	int i,j,k,next,in,ip,nuse,iflg,if1,nit,ier, nmax=10;
	double err,px,reps, xn[11], xd[11], fb1[11];
	
	if(n1>ndim) return 206;
	ier=0;
	reps=0.0;
/*	Find the point nearest to XB1 in X1 */
	next=nearst(xb2,x2,n2);

	nuse=np1; iflg=0;
	i=divdif0(xb1,x1,(f+next*ndim),&nuse,n1,fb1,reps,iflg,&if1);

/*	Set IFLG=1 so that next time DIVDIF0 does not try to locate
	the point again in the table. */
	iflg=1;
	*fb=fb1[nuse];
	xd[1]=*fb;
	xn[1]=x2[next];
	px=1.0;

	ip=next; in=next;

/*	The number of points to be used along X1 */
	nit=np2; if(nmax<nit) nit=nmax; if(n2<nit) nit=n2;
	if(np2<1) {
		nit=4; if(nmax<nit) nit=nmax; if(n2<nit) nit=n2;
	}

/*	Calculate the successive interpolation polynomials */
	for(j=2; j<=nit; ++j) {
	
/*	Find the next nearest point in X1 */
		if(in<=0 ) {
			ip=ip+1; next=ip;
		}
		else if(ip >= n2-1) {
			in=in-1; next=in;
		}
		else if(fabs(xb2-x2[ip+1]) < fabs(xb2-x2[in-1]) ) {
			ip=ip+1; next=ip;
		}
		else {
			in=in-1; next=in;
		}

/*	Interpolate along X2 to calculate function value at (x1[next], xb2) */
		nuse=np1;
		i=divdif0(xb1,x1,(f+next*ndim),&nuse,n1,fb1,reps,iflg,&if1);
		xd[j]=fb1[nuse];
		xn[j]=x2[next];

/*	Calculate the divided differences for interpolation in X1 */
		for(k=j-1; k>=1; --k) xd[k]=(xd[k+1]-xd[k])/(xn[j]-xn[k]);

		px=px*(xb2-xn[j-1]);
		err=xd[1]*px;
		*fb=*fb+err;

	}
	return 0;
}
