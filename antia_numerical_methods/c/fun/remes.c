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
