/*	To calculate rational function minimax approximation using Remes algorithm */

#include <stdio.h>
#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int brentm(double *a, double *b, double *x, double *fx, double reps,
	double aeps, double (*f) (double ));
double fm(double x);
double fun(double x);
double fund(double x);
int remes(int m, int k, int *n, double xl, double xu, double a[], double x[],
	double f[], double ex[], int *ie, double *emax, double eps, double epsm,
	int iflg);

main()
{
	int i,i1,j,k,m,n, iflg, ier,np,nmax;
	double hh, rm,xl,xu,eps,epsm,emax,a[20],e[20],x[500],f[500];

/*	Example 10.10 : Relative minimax approximation to Arc Tan(x) */

	n=100; eps=1.e-8; epsm=1.e-5;
	for(i1=0; i1<99; ++i1) {
		printf("type m=degree of numerator,  k=degree of denominator\n");
		printf("                          (quits when m+k<=0)\n");
		scanf(" %d %d", &m,&k);
		if(m+k<=0) return 0;

		printf("type iflg=0/1/2,   xl=lower limit,    xu=upper limit\n");
		scanf(" %d %le %le", &iflg,&xl,&xu);

		np=m+k+2;
		printf(" m= %d    k = %d    iflg = %d \n",m,k,iflg);
		printf("  lower limit = %e    upper limit = %e \n",xl,xu);

		if(iflg==0) {
			printf("type initial guess for coefficients \n");
			for(i=0; i<=m+k; ++i) scanf(" %le",&a[i]);
			printf("initial guess for coefficients :\n");
			for(i=0; i<=m+k; ++i) printf(" %e",a[i]);
			printf("\n");
		}
		else if(iflg==2) {
			printf("type initial guess for extrema \n");
			for(i=0; i<=m+k+1; ++i) scanf(" %le",&e[i]);
			printf("initial guess for extrema :\n");
			for(i=0; i<=m+k+1; ++i) printf(" %e",e[i]);
			printf("\n");
		}
			
		i=remes(m,k,&n,xl,xu,a,x,f,e,&np,&emax,eps,epsm,iflg);
		printf(" ier= %d    no. of extrema = %d    maximum error =%e \n  ",i,np,emax);
		printf(" Coefficients in numerator : ");
		for(i=k; i<=m+k; ++i) printf(" %e",a[i]);
		printf(" Coefficients in denominator : ");
		for(i=0; i<k; ++i) printf(" %e",a[i]);
		printf("\n  Extrema: ");
		for(i=0; i<np; ++i) printf(" %e",e[i]);
		printf("\n");
	}
	return;
}

double fun(double x)

{
	return 1.0;
}

double fund(double x)

{
	double x1;

	if(x==0.0) return 1.0;
	else {
		x1=sqrt(x);
		return x1/atan(x1);
	}
}



/*	To minimise a function in one dimension using Brent's method

	A,B,X : (input/output) Triplet which brackets the minimum.
		After execution these will contain the final triplet
		with X giving the best estimate for minimiser.
		X must be between A and B; and F(X)<min(F(A),F(B))
	FX : (output) The function value at X.
	reps : (input) Required relative accuracy
	aeps : (input) Required absolute accuracy
		The bracketing interval is subdivided until
		fabs(B-A) < MAX(AEPS, REPS*fabs(X))
	F : (input) Name of the function to calculate the function
		which is to be minimised

	Error status is returned by the value of the function BRENTM.
		0 value implies successful execution
		51 implies that function failed to reduce the
			bracketing interval to required level
		523 implies that initial values of A,X,B do not bracket
			the minimum

	Function F(X) to calculate the required function, must be supplied
		by the user.

	Required functions : F
*/

#include <math.h>

int brentm(double *a, double *b, double *x, double *fx, double reps,
	double aeps, double (*f) (double ))

{
	int i, nit=75;
	double fa,fb,t,e,v,w,fv,fw,xm,eps,eps2,p,r,r1,d,u,fu;
	double gr=0.381966e0;

	if(*a>(*b)) {t=(*a); *a=(*b); *b=t;}
	fa=f(*a); fb=f(*b); *fx=f(*x);
	if(fa<(*fx) || fb<(*fx) || (*x)<=(*a) || (*x)>=(*b) ) return 523;

	e=0.0; v=(*x); w=(*x);
	fv=(*fx); fw=(*fx);

/*	Loop for iteration */
	for(i=1; i<=nit; ++i) {
		xm=0.5*(*a+(*b));
		eps2=reps*fabs(*x); if(aeps>eps2) eps2=aeps;
		eps=0.5*eps2;

/*	The convergence test */
		if(fabs(*x-xm)<eps2-0.5*(*b-(*a))) return 0;

		p=0.0; t=0.0; r=0.0;
		if(fabs(e)>eps) {
/*	Parabolic interpolation */
			r=(*x-w)*(*fx-fv);
			t=(*x-v)*(*fx-fw);
			p=(*x-v)*t-(*x-w)*r;
			t=2.*(t-r);
			if(t>0) p=-p;
			t=fabs(t);
			r=e;
			e=d;
		}

		if(fabs(p)<fabs(0.5*t*r) && p>t*(*a-(*x)) && p<t*(*b-(*x))) {
/*	accept the interpolated point */
			d=p/t; u=(*x)+d;
			if(u-(*a)<eps2 || (*b)-u<eps2) {
/*	If it is too close to end points shift it by eps at least */
				d=eps;
				if((*x)>=xm) d=-eps;
			}
		}
		else {
/*	Perform golden section search */
			e=(*b)-(*x);
			if((*x)>=xm) e=(*a)-(*x);
			d=gr*e;
		}
		if(fabs(d)>=eps) u=(*x)+d;
		else {
/*	Shift the point by at least eps */
			if(d>=0) u=(*x)+eps;
			else u=(*x)-eps;
		}
		fu=f(u);

/*	Updating the bracketing triplet */
		if(fu<=(*fx)) {
/*	(a, u, x) is the triplet */
			if(u<(*x)) *b=(*x);
/*	(x, u, b) is the triplet */
			else (*a)=(*x);
			v=w; fv=fw;
			w=(*x); fw=(*fx);
			*x=u; *fx=fu;
		}
		else {
/*	(u, x, b) is the triplet */
			if(u<(*x)) *a=u;
/*	(a, x, u) is the triplet */
			else *b=u;
			if(fu<=fw || w==(*x)) {
				v=w; fv=fw;
				w=u; fw=fu;
			}
			else if(fu<=fv || v==(*x) || v==w) {v=u; fv=fw;}
		}
	}

	return 51;
}



/*	To calculate coefficients of Rational function minimax approximation
	for a given function over a finite interval using Remes algorithm

	M : (input) Required degree of polynomial in the numerator
	K : (input) Required degree of polynomial in the denominator
	N : (input/output) Number of points to be used for scanning the
		extrema of error curves. If N<3(M+K+2) it will be set to
		this value.
	XL : (input) Lower limit of the interval over which approximation is required
	XU : (input) Upper limit of the interval over which approximation is required
	A : (input/output) Array of length M+K+1 containing the coefficients
		of approximation. A[I-1] is the coefficient of x**I in
		the denominator, the constant term being 1. A[K+J] is
		the coefficient of x**J in the numerator. If IFLG=0 the
		initial guess to coefficients must be supplied.
	X : (output) Array of length N containing the points at which
		function value is calculated for scanning the extrema in
		error curve
	F : (output) Array of length N containing the function values at X[I]
	EX : (input/output) Array of length M+K+5 containing the
		extrema of error curve. If IFLG=2, then the initial guess
		for these extrema must be supplied
	IE : (input/output) Number of extrema in error curve. This value
		must be supplied if IFLG=2
	EMAX : (output) Maximum error in the calculated approximation
	EPS : (input) Required accuracy, the iteration for calculating the
		coefficients of approximations is continued until
		difference in maximum error is less than EPS
	EPSM : (input) Required accuracy to which extrema in error curve
		are determined.
	IFLG : (input) Integer parameter specifying the nature of initial
		approximation that is supplied. 
		If IFLG=0 then iteration is started from initial guess for
			coefficients supplied in array A
		If IFLG=1 then no initial guess is required. Iteration is
			started by assuming that the extrema of error curve
			coincides with those of T_{M+K+1}(x).
		If IFLG=2 then iteration is started from initial guess for
			position of extrema in error curve. These must be
			supplied in array EX and IE should be set to the number
			of extrema
		
	Error status is returned by the value of the function REMES.
		0 value implies successful execution
		614 implies that M<0, K<0, XU<=XL or M+K+2>NMAXR (=50)
			in this case no calculations are done
		632 implies that the Remes iteration did not converge
			to the specified accuracy
		633 implies that at some stage the error curve does not
			have the required number of extrema.
		Other values may be set by GAUELM or BRENTM

	Functions FUN(X) and FUND(X) must be supplied by the user.
	Here X, FUN and FUND are real variables and we seek approximation
	of the form FUN(x) ~ FUND(X)*R_mk(x)
	To obtain minimax approximation with respect to absolute error
	set FUND(X)=1 and FUN(X) to required function
	To obtain minimax approximation with respect to relative error
	set FUN(X)=1 and FUND(X) to reciprocal of required function.

	Required functions : GAUELM, BRENTM, FM, FUN, FUND
*/

#include <math.h>
#include <stdlib.h>

/* To pass on these parameters to FM put these statements before that */
#define NMAXR 50
#define PI 3.14159265358979324
int MM, KK;
double AA[NMAXR],SI;

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int brentm(double *a, double *b, double *x, double *fx, double reps,
	double aeps, double (*f) (double ));
double fm(double x);

int remes(int m, int k, int *n, double xl, double xu, double a[], double x[],
	double f[], double ex[], int *ie, double *emax, double eps, double epsm,
	int iflg)

{
	int i,j,ier,nk,it,jt,ifl, nit=30, njt=10;
	double h,reps,e0,e1,e2,emin,ei,ai,b,xi,fx,r1,fi,fd,dif,det;
	int *iwk;
	double *wk;

	nk=m+k+2;
	if(m<0 || k<0 || nk>NMAXR || xu <= xl) return 614;
	
/*	Copy the value in the global variables for use by function FM */
	MM=m; KK=k;
	if(*n < 3*nk) *n=3*nk;

/*	Generating the mesh points for crude scan of extrema */
	h=(xu-xl)/(*n-1.0);
	for(i=0; i<(*n); ++i) {
		x[i]=xl+i*h;
		f[i]=fun(x[i]);
	}

	reps=0.0;
	wk=(double *) calloc((size_t) (nk*nk), sizeof(double));
	iwk=(int *) calloc((size_t) nk, sizeof(int));

/*	Loop for Remes iteration */
	for(it=0; it<nit; ++it) {
		e1=0.0;
/*	Copy the coefficients to global variables for function FM */
		for(i=0; i<=m+k; ++i) AA[i]=a[i];

		if(it==0 && iflg >= 1) {
			if(iflg==1) {
/*	Use the extrema of Chebyshev polynomial T_{M+K+1} as initial guess */
				ex[0]=x[0];
				*ie=1;
				for(i=1; i<nk-1; ++i) ex[i]=0.5*(xl+xu)+0.5*(xu-xl)*cos((nk-i-1)*PI/(nk-1));
				ex[nk-1]=xu;
			}
			ei=0.0;
			*emax=0.0;
			emin=0.0;
		}
		else {

/*	Locate the extrema of the error curve */
			for(i=0; i<(*n); ++i) {
/*	Flag to avoid calculating the function repeatedly */
				SI=22.0;
				ei=f[i]+fm(x[i]);
				if(i>1) {
					if((e1 >= e2) == (e1 >= ei)) {
/*	An extrema is bracketed */
						SI=1.0;
/*	To convert maxima to minima */
						if(e1>e2 || e1>ei) SI=-1.0;
						ai=x[i-2]; b=x[i]; xi=x[i-1];
						ier=brentm(&ai,&b,&xi,&fx,reps,epsm,fm);
						if(ier > 0) {free(iwk); free(iwk); return ier;}

						ex[*ie]=xi;
						wk[*ie]=fx*SI;
						*ie=(*ie)+1;
						if(fabs(fx)>(*emax)) (*emax)=fabs(fx);
						if(fabs(fx)<emin) emin=fabs(fx);
/*	To ensure that dimensions of EX are not exceeded */
						if(*ie> m+k+4) goto ext;
					}
				}
				else if(i==0) {
/*	The end point is always included in the list of extrema */
					*ie=1;
					ex[0]=x[0];
					*emax=fabs(ei);
					emin=(*emax);
					wk[i]=ei;
					e0=ei;
				}
				e2=e1; e1=ei;
			}

/*	The end point is always included in the list of extrema */
			ex[*ie]=x[*n-1];
			wk[*ie]=ei; *ie =(*ie)+1;
			if(fabs(ei)>(*emax)) (*emax)=fabs(ei);
			if(fabs(ei)<emin) emin=fabs(ei);

/*	If the number of extrema is less than required then quit */
			if(*ie < nk) {free(iwk); free(wk); return 633;}

ext:		while(*ie > nk) {
/*	Remove one extrema from the list */
				*ie=(*ie)-1;
				if(fabs(wk[*ie]) > fabs(wk[0])) {
/*	remove the first extrema */
					emin=fabs(wk[1]);
					e0=wk[1];
					for(i=0; i<(*ie); ++i) {
						ex[i]=ex[i+1];
						if(fabs(wk[i+1])<emin) emin=fabs(wk[i+1]);
						wk[i]=wk[i+1];
					}
				}

				else {
/*	remove the last extrema */
					emin=fabs(wk[0]);
					for(i=1; i<(*ie); ++i) {
						if(fabs(wk[i])<emin) emin=fabs(wk[i]);
					}
				}
			}

/*	Convergence check, quit if difference between various extrema is less than 1% */
			if(*emax-emin < 0.01*(*emax)) {free(iwk); free(wk); return 0;}
			ei=0.5*(*emax+emin); if(ei>2.*emin) ei=2.*emin;
			if(e0<0.0) ei=-ei;
		}

/*	Iteration to determine the coefficients */
		for(jt=0; jt<njt; ++jt) {
			ai=1.0;
/*	Setting up the equation matrix */
			for(i=0; i<nk; ++i) {
				fi=fun(ex[i]);
				fd=fund(ex[i]);
				r1=1.0;
				for(j=0; j<k; ++j) {
					r1=r1*ex[i];
					wk[j+i*nk]=-(fi-ai*ei)*r1;
				}
				wk[k+i*nk]=fd;
				r1=1.0;
				for(j=0; j<m; ++j) {
					r1=r1*ex[i];
					wk[j+k+1+i*nk]=fd*r1;
				}
				wk[nk-1+i*nk]=ai;
				a[i]=fi;
				ai=-ai;
			}

			ifl=0;
			ier=gauelm(nk,1,wk,a,&det,iwk,nk,&ifl);
			if(ier>0) {free(iwk); free(wk); return ier;}
			dif=fabs(a[nk-1]-ei);
			ei=a[nk-1];
/*	convergence check */
			if(dif<eps || dif<0.3*(*emax-emin)) break;
			if(it==0 && iflg>=1) *emax=0.3*fabs(a[nk-1]);
		}

/*	Even if iteration for calculating coefficients fails continue further */

	}

	free(iwk); free(wk);
	return 632;
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



/*	To calculate the error in rational function approximation
	For use with function REMES. It is called by BRENTM to find
	extrema of error curve. It may not be used for any other purpose.

	X : (input) the value at which error is to be calculated
	FM = (FUN(X) - FUND(X)*R_mk(X))*SI  is the calculated error

	The parameters for rational function are passed through global variables
	MM,KK are degree of numerator and denominator
	AA is an array containing the coefficient of rational function approximation
	SI is -1 for maximum and +1 for minimum. The error is multiplied
		by SI to make it a minimum.
		For initial scan SI>10 and function is not evaluated.

	Functions FUN(x) and FUND(x) must be provided by the user to
		seek approximation of form FUN(X) = FUND(X)*RMK(X)

	Required functions : FUN, FUND
*/

/* To pass on these parameters to FM put these statements before that */
#define NMAXR 50
int MM, KK;
double AA[NMAXR],SI;

double fm(double x)

{
	int j,nk;
	double fd,fn,f;

	nk=MM+KK;
	f=0.0;
	if(SI<10.0) f=fun(x);

/*	Calculate the numerator using nested multiplication */
	fn=AA[KK+MM];
	for(j=1; j<=MM; ++j) fn=fn*x+AA[nk-j];

	if(KK>0) {
/*	Calculate the denominator using nested multiplication */
		fd=AA[KK-1];
		for(j=2; j<=KK; ++j) fd=fd*x+AA[KK-j];
		fd=fd*x+1.0;
	}
	else fd=1.0;

	f=f-fund(x)*fn/fd;
	if(SI<0.0) f=-f;
	return f;
}
