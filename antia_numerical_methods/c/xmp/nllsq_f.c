/*	Nonlinear least squares fit using direction set method */

#include <stdio.h>
#include <math.h>

void fcn(int n, double a[], double x, double *f);
void nllsq_f(int n, double a[], double *f);
double fln(void fcn(int , double * , double * ), double x, double v[],
	double x0[], int n, int *num);
double rangau(double *seed);

int linmnf(double *x0, double *x1, double *f0, double reps, double aeps,
	void (*f) (int , double * , double * ), double v[], double xi[],
	int n, int *num);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int nminf(int n, double x[], double *f, int *num, double reps, double aeps,
	void (*fcn) (int , double * , double * ));

#define NPL 100
int NNLSQ;
double FXLSQ[NPL], XLSQ[NPL], EFLSQ[NPL], FX1LSQ[NPL];

/*	Example 10.4 */

/*	With most starting values the iteration does not converge to
	satisfactory accuracy, but the final result may still be acceptable,
	The sequence of iterates also depend on roundoff errors and it is
	unlikely that the results given in output file will be reproduced */

main()
{
	int i,i1,j,nuse, iflg, ier,np,nmax;
	double x0[20], seed, reps, aeps, hh, fx;

	aeps=1.e-8; reps=1.e-7;
	np=6; seed=2;
/*	The number of parameters np should be even */

	for(i1=0; i1<99; ++i1) {
		printf("type np=No. of parameters,    NNLSQ=No. of data pts\n");
		printf("                   (quits when np<=1)\n");
		scanf(" %d %d",&np,&NNLSQ);
		if(np<=1) return 0;

/*	Generating input data with random error added */
		hh=1.0/(NNLSQ-1.);
		for(i=0; i<NNLSQ; ++i) {
			XLSQ[i]=i*hh;
			EFLSQ[i]=1.e-2;
			FXLSQ[i]=exp(-XLSQ[i])+2.*exp(-2.*XLSQ[i])+3.0*exp(-3.*XLSQ[i])+1.e-6*rangau(&seed);
		}
	
		printf("type initial guess for parameters \n");
		for(i=0; i<np; ++i) scanf(" %le", &x0[i]);
		printf(" Initial guess : ");
		for(i=0; i<np; ++i) printf(" %e", x0[i]);
		printf("\n");

		i=nminf(np,x0,&fx,&nuse,reps,aeps,nllsq_f);
		printf(" ier = %d   No. of data pts = %d    No. of parameters %d  no. of function evaluations = %d  \n x = ", i,NNLSQ,np,nuse);
		printf("  Squared difference = %e   coefficients : \n",fx);
		for(i=0; i<np; ++i) printf(" %e", x0[i]);
		printf(" \n  ");
	}
	return;
}



void fcn(int n, double a[], double x, double *f)

{
	int i;

	*f=0.0;
	for(i=0; i<n; i += 2) *f=(*f)+a[i]*exp(-x*a[i+1]);

	return;
}




/*	Function to calculate the function value as required for
	line search without derivatives

	FCN : (input) Name of function to calculate the required function
	X : (input) Parameter along the line to specifying the point where
		function evaluation is required
	V : (input/output) Array of length 2N, first N elements specify the
		direction of line search. After execution next N elements will contain
		the coordinates of the point at which function is evaluated.
	X0 : (input) Array of length N, containing the coordinates
		of the starting point for line search
	N : (input) Number of variables in the function to be minimised
	NUM : (input/output) Integer variable to keep count of function evaluations

	Function FCN(N,X,FX) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, FX is the
		function value at X. X is an array of length N.

	Required functions : FCN
*/

double fln(void fcn(int , double * , double * ), double x, double v[],
	double x0[], int n, int *num)

{
	int i;
	double fl;

	*num=(*num)+1;
/*	coordinates of the required points */
	for(i=0; i<n; ++i) x0[n+i]=x0[i]+v[i]*x;
	fcn(n,&x0[n],&fl);
	return fl;
}




/*	To perform a line search for minimum of several variables as required
	by direction set method
	This function should not be used for any other purpose

	X0 : (input/output) Starting value for the line search. After execution
		it should contain the distance to minimiser along the line
	X1 : (input/output) Initial estimate for the minimum along the line.
		 This value will be modified by the function.
	F0 : (input/output) The function value at X1, this value must be supplied
	REPS : (input) Required relative accuracy 
	AEPS : (input) Required absolute accuracy,
		This criterion is only used to terminate line search under
		certain conditions and not generally applicable to the line
		search. These values should be same as what is used in function NMINF
	F : (input) Name of the function routine to calculate the function value
	V : (input/output) Array of length 2N. First N element specify
		the direction in which minimisation is required. The next
		N elements will contain the coordinates of the minimiser
		found by LINMNF.
	XI : (input) Array of length N containing the coordinates of
		starting point for line search
	N : (input) Number of variables in the function to be minimised
	NUM : (output) Integer variable to keep count of the number of
		function evaluations used so far
		
	Error status is returned by the value of the function LINMNF.
		0 value implies successful execution
		56 implies that the function failed to find acceptable point
			in this case the starting point itself is accepted

	Function F(N,X,FX) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, FX is the
		function value at X. X is an array of length N.

	Required functions :  FLN, F
*/

#include <math.h>

double fln(void fcn(int , double * , double * ), double x, double v[],
	double x0[], int n, int *num);

int linmnf(double *x0, double *x1, double *f0, double reps, double aeps,
	void (*f) (int , double * , double * ), double v[], double xi[],
	int n, int *num)

{
	int i,nit=25;
	double f1,x2,f2,r,r1,r2,t,p,xp,fp;
	double t1=9.0, gr=1.618034, gc=0.381966;

	f1=fln(f,*x1,v,xi,n,num);
	if(f1<(*f0)) x2=(*x0)+gr*(*x1-(*x0));
	else x2=(*x0)-gr*(*x1-(*x0));
	f2=fln(f,x2,v,xi,n,num);

	for(i=1; i<=nit; ++i) {
/*	Parabolic interpolation */
		r=(x2-(*x1))*(f2-(*f0));
		t=(x2-(*x0))*(f2-f1);
		p=(x2-(*x0))*t-(x2-(*x1))*r;
		t=2.*(t-r);
		if(t>0.0) p=-p;
		t=fabs(t);
		xp=(*x0);
		if(fabs(p)<fabs(t1*t*(x2-(*x1)))) {
/*	Try the interpolated value */
			xp=x2+p/t;
			fp=fln(f,xp,v,xi,n,num);
			if(fp<(*f0)) {
/*	accept the point */
				*x0=xp; *f0=fp;
				return 0;
			}
		}

		if(f1<(*f0)) {*x0=(*x1); *f0=f1; return 0;}
		else if(f2<(*f0)) {*x0=x2; *f0=f2; return 0;}
		else if(xp != (*x0)) {
/*	subdivide the interval */
			if( (xp-(*x0) > 0.0) == (*x1-(*x0) > 0.0) ) {
				if(fabs(xp-(*x0))<0.5*fabs(*x1-(*x0))) {*x1=xp; f1=fp;}
				else {
/*	use golden section */
					*x1=(*x0)+gc*(*x1-(*x0));
					f1=fln(f,*x1,v,xi,n,num);
				}
			} else {
				if(fabs(xp-(*x0))<0.5*fabs(x2-(*x0))) {x2=xp; f2=fp;}
				else {
/*	use golden section */
					x2=(*x0)+gc*(x2-(*x0));
					f2=fln(f,x2,v,xi,n,num);
				}
			}

/*	If interpolated point is not acceptable use golden section */
		} else if(fabs(x2-(*x0))>fabs(*x1-(*x0))) {
			x2=(*x0)+gc*(x2-(*x0));
			f2=fln(f,x2,v,xi,n,num);
		} else {
			*x1=(*x0)+gc*(*x1-(*x0));
			f1=fln(f,*x1,v,xi,n,num);
		}

/*	If the change in function value is too small, then quit */
		r1=reps*fabs(*f0); if(aeps>r1) r1=aeps;
		r2=f1-(*f0); if(f2-(*f0) < r2) r2=f2-(*f0);
		if(r2<r1) return 0;
	}

	return 56;
}




/*	To minimise a function of several variables using direction set method

	N : (input) Number of variables
	X : (input/output) Array of length N containing the initial
		guess for the minimum.
		After execution it should contain the coordinates of minimiser
	F : (output) The function value at X
	NUM : (output) Number of function evaluations used by NMINF
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy, iteration will stop when
		change in function value is less than max(AEPS, REPS*fabs(F))
	FCN : (input) Name of the function to calculate the function value
		
	Error status is returned by the value of the function NMINF.
		0 value implies successful execution
		504 implies that N < 2, in which case no calculations are done
		529 implies that iteration failed to converge to specified accuracy
		Other values may be set by LINMNF

	Function FCN(N,X,F) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, F is the
		function value at X. X is an array of length N.

	Required functions : LINMNF, FLN, SVD, FCN
*/

#include <math.h>
#include <stdlib.h>

int linmnf(double *x0, double *x1, double *f0, double reps, double aeps,
	void (*f) (int , double * , double * ), double v[], double xi[],
	int n, int *num);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);

int nminf(int n, double x[], double *f, int *num, double reps, double aeps,
	void (*fcn) (int , double * , double * ))

{
	int i,j,it,k,kmax,ier,ier1, nit=200;
	double f0,fi,dfmax,x1,x2,fun,r1;
	double *wk, *wk1, *u;

	if(n<2) return 504;
	wk=(double *) calloc((size_t) (n*n),sizeof(double));
	wk1=(double *) calloc((size_t) (4*n),sizeof(double));

/*	Initialise the direction set matrix to identity matrix */
	for(i=0; i<n; ++i) {
		wk1[i+2*n]=0.0;
		wk1[i+3*n]=0.0;
		for(j=0; j<n; ++j) wk[j+i*n]=0.0;
		wk[i+i*n]=1.0;
	}

	fcn(n,x,f);
	*num=1;
	ier=0;

/*	The main iteration loop */
	for(it=1; it<=nit; ++it) {
		fi=(*f);
		for(k=0; k<n; ++k) {

/*	The starting point for line search */
			for(i=0; i<n; ++i) wk1[i]=x[i];
			f0=(*f);
			dfmax=0.0; kmax=0;

			for(i=0; i<n; ++i) {
				x1=0.0;
/*	Use previous value as initial approximation to minimum */
				x2=wk1[i+2*n];
				fun=f0;
				if(x2==0.0) x2=1.0;
				ier1=linmnf(&x1,&x2,&f0,reps,aeps,fcn,&wk[i*n],wk1,n,num);
				if(ier1>0) ier=ier1;

				wk1[i+2*n]=x1;
/*	Estimate of second derivative along the line */
				if(x1 != 0.0) wk1[i+3*n]=fabs((fun-f0)/(x1*x1));
/*	The new starting point */
				for(j=0; j<n; ++j) wk1[j]=wk1[j]+x1*wk[j+i*n];
				if(fun-f0>=dfmax && i<=n-k-1) {dfmax=fun-f0; kmax=i;}
			}

/*	Remove the KMAX th direction */
			for(i=kmax; i<n-1; ++i) {
				wk1[i+2*n]=wk1[i+1+2*n];
				wk1[i+3*n]=wk1[i+1+3*n];
				for(j=0; j<n; ++j) wk[j+i*n]=wk[j+i*n+n];
			}
/*	Add a new direction */
			for(i=0; i<n; ++i) wk[i+n*n-n]=wk1[i]-x[i];

			x1=0.0; x2=1.0;
			fun=(*f);
			wk1[4*n-1]=0.0;
/*	Starting point for the final line search in the loop */
			for(i=0; i<n; ++i) wk1[i]=x[i];
			ier1=linmnf(&x1,&x2,f,reps,aeps,fcn,&wk[n*n-n],wk1,n,num);
			if(ier1>0) ier=ier1;

			wk1[3*n-1]=x1;
			if(x1!=0.0) wk1[4*n-1]=fabs((fun-(*f))/(x1*x1));
			for(j=0; j<n; ++j) x[j]=x[j]+x1*wk[j+n*n-n];
		}

		r1=reps*fabs(*f); if(aeps>r1) r1=aeps;
		if(fabs(*f-fi)<r1) {free(wk1); free(wk); return ier;}

/*	The matrix V for SVD */
		for(j=0; j<n; ++j) {
			if(wk1[j+3*n]>0.0) {
				for(i=0; i<n; ++i) wk[i+j*n]=wk[i+j*n]/sqrt(wk1[j+3*n]);
			}
		}
/*	Take the transpose of matrix as the array is stored in different order */
		for(i=0; i<n; ++i){
			for(j=i+1; j<n; ++j)
				{r1=wk[i*n+j]; wk[i*n+j]=wk[i+j*n]; wk[i+j*n]=r1;}
		}

		u=(double *)calloc((size_t) (n*n),sizeof(double));
		ier=svd(n,n,wk,u,wk1,n,n);
		free(u);
		if(ier>0) {free(wk1); free(wk); return ier;}
/*	Transpose the matrix back to earlier order */
		for(i=0; i<n; ++i){
			for(j=i+1; j<n; ++j) {r1=wk[i*n+j]; wk[i*n+j]=wk[i+j*n]; wk[i+j*n]=r1;}
		}

		for(i=0; i<n; ++i) {
			wk1[i+2*n]=sqrt(fun-(*f))*wk1[i];
			wk1[i+3*n]=0.0;
			if(wk1[i]!=0.0) wk1[i+3*n]=1./(wk1[i]*wk1[i]);
		}
	}

	free(wk1); free(wk);
	return 529;
}




/*	To calculate the Chi square function for a nonlinear
	least squares fit, for use with function NMINF.
	Version of NLLSQ for use with NMINF, when derivatives are not available

	N : (input) Number of parameters to be fitted
	A : (input) Array of length N containing the parameter values
		at which the function needs to be evaluated
	F : (output) The function value at specified parameter values

	The data points are passed through global variables NNLSQ, FXLSQ,
	XLSQ, EFLSQ, FX1LSQ, which must be initialised in the calling
	function before calling NMINF

	NNLSQ : (input) The number of data points in the table to be fitted
	FXLSQ : (input) Array of length NPL containing the function values
	XLSQ : (input) Array of length NPL containing the abscissas at which
		function values are available
	EFLSQ : (input) Array of length NPL containing the estimated errors in 
		function values, for use in deciding the weights for each point
	FX1LSQ : (output) Array of length NPL containing the fitted value
		of the function at each tabular point

	The parameter NPL must be equal to the dimension of arrays as
	declared in the calling function.

	This function requires function FCN to calculate the required
	function which has to be fitted. There is no provision to pass on
	the name of this function and hence it must be changed explicitly
	to the required name. Function FCN(N,A,X,F) must be supplied
	by the user. Here N is the number of parameters, A is an array
	of length N containing the values of parameters and X is the value
	of independent variable where the fitting function needs to be evaluated.
	F is the calculated function value.

	Required functions : FCN
*/

#include <math.h>


/* To pass on these parameters to nllsq_f, this must be included before
   the calling NMINF */

void fcn(int n, double a[], double x, double *fa);

#define NPL 100
int NNLSQ;
double FXLSQ[NPL], XLSQ[NPL], EFLSQ[NPL], FX1LSQ[NPL];

void nllsq_f(int n, double a[], double *f)

{
	int i,j;
	double fa,r1;

/*     GENERATING THE FUNCTION FOR MINIMISATION */

	*f=0.0;

	for(i=0; i<NNLSQ; ++i) {
		fcn(n,a,XLSQ[i],&fa);
		FX1LSQ[i]=fa;

/*     SUM OF SQUARED DIFFERENCE */

		r1=(fa-FXLSQ[i])/EFLSQ[i];
		*f=(*f)+r1*r1;
	}
	return;
}




/*	To calculate the Singular Value Decomposition of a matrix A=U D Vtranspose

	N : (input) Number of variables
	M : (input) Number of equations
	A : (input/output) Matrix of coefficients of size LA*M
		After execution it will contain the matrix U
	V : (output) The matrix V of size LV*N
	SIGMA : (output) Array of length N, containing the singular values
	LA : (input) Actual value of second dimension of A in the calling function
	LV : (input) Actual value of second dimension of V in the calling function
		
	Error status is returned by the value of the function SVD.
		0 value implies successful execution
		12 QR iteration failed to converge to required accuracy
		105 implies N<=0, N>LV, M<=0, N>LA, N>M

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv)

{
	int i,j,k,l,itr,ier, itmax=30;
	double f, g, h, rmax, s, r1, r2, c, x, y, z, aeps, eps=1.e-16;
	double *e;

	if(n>m || n<=0 || m<=0 || n>la || n>lv) return 105;
	ier=0;

/*	Reduction to Bidiagonal form using Householder transformations */
	g=0.0; rmax=0.0;
	e=(double *) calloc((size_t) n, sizeof(double));

	for(i=0; i<n; ++i) {
/*	Off-diagonal elements of bidiagonal form  */
		e[i]=g;
		s=0.0;
		for(j=i; j<m; ++j) s=s+a[i+j*la]*a[i+j*la];
		if(s <= 0.0) {
/*	transformation not required */
			g=0.0;
		}
		else {
			f= a[i+i*la];
			g=sqrt(s);
			if(f>=0.0) g=-g;
			h=f*g-s;
			a[i+i*la] = f-g;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]= a[j+k*la]+f*a[i+k*la];
			}
		}

/*	Diagonal elements of bidiagonal form  */
		sigma[i]=g;
		s=0.0;
		for(j=i+1; j<n; ++j) s=s+a[j+i*la]*a[j+i*la];

		if(s<= 0.0) g=0.0;
		else {
			f= a[i*la+(i+1)];
			g=sqrt(s);
			if(f>= 0.0) g=-g;
			h=f*g-s;
			a[i*la+(i+1)]=f-g;
			for(j=i+1; j<n; ++j) e[j]=a[j+i*la]/h;

			for(j=i+1; j<m; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+j*la]*a[k+i*la];
				for(k=i+1; k<n; ++k) a[k+j*la] = a[k+j*la]+s*e[k];
			}
		}
		r1=fabs(sigma[i])+fabs(e[i]);
		if(r1 > rmax) rmax=r1;
	}

/*	Accumulation of right hand transformation in array V */
	for(i=n-1; i>=0; --i) {
		if(g != 0.0) {
			h=g*a[i*la+(i+1)];
			for(j=i+1; j<n; ++j) v[i+j*lv]=a[j+i*la]/h;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+i*la]*v[j+k*lv];
				for(k=i+1; k<n; ++k) v[j+k*lv]=v[j+k*lv]+s*v[i+k*lv];
			}
		}

		for(j=i+1; j<n; ++j) {
			v[j+i*lv]=0.0; v[i+j*lv]=0.0;
		}
		v[i+i*lv]=1;
		g= e[i];
	}

/*	Accumulation of left hand transformation overwritten on matrix A */
	for(i=n-1; i>=0; --i) {
		g=sigma[i];
		for(j=i+1; j<n; ++j) a[j+i*la]=0.0;
		if(g != 0.0) {
			h=g*a[i+i*la];

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]=a[j+k*la]+f*a[i+k*la];
			}

			for(j=i; j<m; ++j) a[i+j*la]=a[i+j*la]/g;
		}
		else {
			for(j=i; j<m; ++j) a[i+j*la]=0.0;
		}
		a[i+i*la] = a[i+i*la]+1;
	}

/*	Diagonalisation of the bidiagonal form */
	aeps=eps*rmax;
/*	Loop over the singular values */
	for(k=n-1; k>=0; --k) {
/*	The QR transformation */
		for(itr=1; itr<=itmax; ++itr) {

/*	Test for splitting */
			for(l=k; l>=0; --l) {
				if(fabs(e[l]) < aeps) goto split;
				if(fabs(sigma[l-1]) < aeps) break;
			}

/*	cancellation of E[L] if L>1  */
			c=0.0; s=1.0;
			for(i=l; i<=k; ++i) {
				f=s*e[i];
				e[i] = c*e[i];
				if(fabs(f) < aeps) goto split;
				g=sigma[i];
				sigma[i]=sqrt(f*f+g*g);
				c=g/sigma[i];
				s=-f/sigma[i];

				for(j=0; j<m; ++j) {
					r1= a[j*la+(l-1)];
					r2= a[i+j*la];
					a[j*la+(l-1)]=r1*c+r2*s;
					a[i+j*la]=c*r2-s*r1;
				}
			}

split:			z=sigma[k];
			if(l == k) {
/*	QR iteration has converged */
				if(z < 0.0) {
					sigma[k] = -z;
					for(j=0; j<n; ++j) v[k+j*lv]=-v[k+j*lv];
				}
				break;
			}

			if(itr==itmax) {ier=12; break;}
					
/*	calculating shift from bottom 2x2 minor */
			x=sigma[l];
			y=sigma[k-1];
			g=e[k-1];
			h=e[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g=sqrt(1.+f*f);
			if(f < 0.0) g=-g;
			f=((x-z)*(x+z)+h*(y/(f+g)-h))/x;

/*	next QR transformation */
			c=1.0; s=1.0;
/*	Given's rotation  */
			for(i=l+1; i<=k; ++i) {
				g=e[i];
				y=sigma[i];
				h=s*g;
				g=c*g;
				e[i-1]=sqrt(f*f+h*h);
				c=f/e[i-1];
				s=h/e[i-1];
				f=c*x+s*g;
				g=c*g-s*x;
				h=s*y;
				y=c*y;

				for(j=0; j<n; ++j) {
					x=v[j*lv+(i-1)];
					z=v[i+j*lv];
					v[j*lv+(i-1)]=c*x+s*z;
					v[i+j*lv]=c*z-s*x;
				}

				sigma[i-1]=sqrt(f*f+h*h);
				if(sigma[i-1] != 0.0) {
					c=f/sigma[i-1];
					s=h/sigma[i-1];
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for(j=0; j<m; ++j) {
					y= a[j*la+(i-1)];
					z= a[i+j*la];
					a[j*la+(i-1)] = c*y+s*z;
					a[i+j*la] = c*z-s*y;
				}
			}

			e[l]=0.0;
			e[k]=f;
			sigma[k]=x;
		}
	}
	free(e);
	return ier;
}




/*	To generate random numbers with Gaussian probability distribution
	It generates random numbers with zero mean and variance of 1.
	
	SEED : (input/output) real seed, it should be positive and
		less than AN. It is updated by the function and should
		not be modified between two calls, unless a fresh
		sequence is required

	Required functions : None
	
	THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
	VERSION AS THE SEED IS NOW DOUBLE INSTEAD OF INT.

*/

#include <math.h>


double rangau(double *seed)

{
	int n2;
	double am=2147483648.0, a=45875.0, ac=453816693.0, an=2147483647.0, r1, rn;

	rn=a*(*seed)+ac; n2=rn/am; r1=rn-am*n2;
	if(*seed==0.0) *seed=0.1;
	rn=sqrt(2.0*log(an/(*seed)))*cos(2.0*M_PI*r1/an);
	*seed=r1;
	return rn;
}

