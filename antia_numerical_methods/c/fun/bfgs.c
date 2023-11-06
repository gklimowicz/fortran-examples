/*	To minimise a function of several variables using quasi-Newton method

	N : (input) Number of variables
	X : (input/output) Array of length N containing the initial
		guess for the minimum.
		After execution it should contain the coordinates of minimiser
	F : (output) The function value at X
	G : (output) Array of length N containing the gradient vector at X
	H : (output) Array of length N*N containing the estimated
		Hessian matrix at X. The second dimension of H should be N.
	NUM : (output) Number of function evaluations used by BFGS
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy, all components of the
		Minimiser will be calculated with accuracy MAX(AEPS, REPS*ABS(X[I]))
	FCN : (input) Name of the function routine to calculate the function value
		and its derivatives

	Error status is returned by the value of the function BFGS.
		0 value implies successful execution
		53 implies that Hessian matrix at final point is probably singular
			The iteration may have converged to a saddle point
		503 implies that N < 1, in which case no calculations are done
		526 implies that iteration failed to converge to specified accuracy
		Other values may be set by LINMIN

	Function FCN(N,X,F,G) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, F is the
		function value at X and G is the gradient vector. X and G
		are arrays of length N. F and G must be calculated by FCN.

	Required functions : LINMIN, FLNM, FCN
*/

#include <math.h>
#include <stdlib.h>

int linmin(double *x1, double *x2, double *f1, double *df1, double reps,
	double aeps, void (*f) (int , double * , double * , double * ),
	double v[], double xi[], int n, int *num);

int bfgs(int n, double x[], double *f, double g[], double *h, int *num,
	double reps, double aeps, void fcn(int , double * , double * , double * ))

{
	int i,j,it,n2,qc,ier,ier1, nit=200;
	double df,df1,h1,h2,x1,x2,f1,r1,gi,dg,ghg;
	double *wk;

	if(n<1) return 503;

	for(i=0; i<n; ++i) {
/*	Initialise the Hessian matrix to unit matrix */
		for(j=0; j<n; ++j) h[j*n+i]=0.0;
		h[i*n+i]=1.0;
	}
/*	If some variable needs to be kept fixed the corresponding
	diagonal elements, h[i][i] should be set to zero. */

	fcn(n,x,f,g);
	df=fabs(*f);
	if(df==0.0) df=1.0;
	n2=2*n;
	*num=1;
	ier=0;
	h2=1.0;
	wk=(double *) calloc((size_t) (3*n), sizeof(double));

/*	The iteration loop */
	for(it=1; it<=nit; ++it) {
		df1=0.0; h1=h2; h2=0.0;
/*	Calculating the search direction wk =S^(k) */
		for(i=0; i<n; ++i) {
			wk[i]=0.0;
			for(j=0; j<n; ++j) {
				h2=h2+fabs(h[i*n+j]);
				wk[i]=wk[i]-h[i*n+j]*g[j];
			}
			df1=df1+wk[i]*g[i];
		}

		if(df1==0.0) {
/*	If gradient vanishes, then quit
	If Hessian matrix appears to be singular, set the error flag */
			if(fabs(h2/h1)>1.3) {free(wk); return 53;}
			else {free(wk); return 0;}
		}

/*	Initial guess for line search */
		x1=0.0;
		x2=-2.*df/df1; if(x2>1.0) x2=1.0;
		f1=(*f);
		if(x2<=0.0) x2=1.0;
		ier1=linmin(&x1,&x2,&f1,&df1,reps,aeps,fcn,wk,x,n,num);
		if(ier1>0) ier=ier1;
/*	If line search fails, then quit */
		if(ier>100) {free(wk); return ier;}

/*	The convergence test */
		qc=1;
		for(i=0; i<n; ++i) {
			x[i]=wk[n+i];
			wk[n+i]=x1*wk[i];
			r1=reps*fabs(x[i]); if(aeps>r1) r1=aeps;
			if(fabs(wk[n+i])>r1) qc=0;
			gi=wk[n2+i];
			wk[n2+i]=wk[n2+i]-g[i];
			g[i]=gi;
		}
/*	It is possible to apply convergence check on Function value using df
	instead of x[i] */
		df=(*f)-f1;
		*f=f1;
		if(qc==1) {
			if(fabs(h2/h1)>1.3) {free(wk); return 53;}
			else {free(wk); return ier;}
		}

/*	Update the matrix using BFGS formula */
		for(i=0; i<n; ++i) {
			wk[i]=0.0;
			for(j=0; j<n; ++j) wk[i]=wk[i]+h[i*n+j]*wk[n2+j];
		}
		ghg=0.0; dg=0.0;
		for(i=0; i<n; ++i) {
			dg=dg+wk[n2+i]*wk[n+i];
			ghg=ghg+wk[i]*wk[n2+i];
		}
		r1=(1.+ghg/dg)/dg;
		for(j=0; j<n; ++j) {
			for(i=0; i<n; ++i) h[i*n+j]=h[i*n+j]+r1*wk[n+i]*wk[n+j]-
				(wk[n+i]*wk[j]+wk[i]*wk[n+j])/dg;
		}
	}

	free(wk);
	return 526;
}
