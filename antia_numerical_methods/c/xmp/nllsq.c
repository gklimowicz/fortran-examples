/*	Nonlinear least squares fit using quasi newton method */

#include <stdio.h>
#include <math.h>

double rangau(double *seed);
void fcn(int n, double a[], double x, double *f, double g[]);
double flnm(void (*fcn) (int , double * , double * , double * ), double x,
	double *df, double v[], double x0[], int n, int *num);

int linmin(double *x1, double *x2, double *f1, double *df1, double reps,
	double aeps, void (*f) (int , double * , double * , double * ),
	double v[], double xi[], int n, int *num);
int bfgs(int n, double x[], double *f, double g[], double *h, int *num,
	double reps, double aeps, void fcn(int , double * , double * , double * ));
void nllsq(int n, double a[], double *f, double g[]);

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
	double x0[20],g[20], seed, reps, aeps, hh,h[6][6], fx;

	aeps=1.e-6; reps=1.e-5;
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

		i=bfgs(np,x0,&fx,g,&h[0][0],&nuse,reps,aeps,nllsq);
		printf(" ier = %d   No. of data pts = %d    No. of parameters %d  no. of function evaluations = %d  \n x = ", i,NNLSQ,np,nuse);
		printf("  Squared difference = %e   coefficients : \n",fx);
		for(i=0; i<np; ++i) printf(" %e", x0[i]);
		printf(" \n  gradient = ");
		for(i=0; i<np; ++i) printf(" %e", g[i]);
		printf(" \n ");
	}
	return;
}

/*	The fitting function with gradient */

void fcn(int n, double a[], double x, double *f, double g[])

{
	int i;

	*f=0.0;
	for(i=0; i<n; i += 2) *f=(*f)+a[i]*exp(-x*a[i+1]);

	for(i=0; i<n; i += 2) {
		g[i]=exp(-x*a[i+1]);
		g[i+1]=-x*a[i]*exp(-x*a[i+1]);
	}
	return;
}




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



/*	To perform a line search for minimum of several variables as required
	by Quasi-Newton methods
	This function should not be used for any other purpose

	X1 : (input/output) Starting value for the line search. After execution
		it should contain the distance to minimiser along the line
	X2 : (input/output) Initial estimate for the minimum along the line.
		 This value will be modified by the function.
	F1 : (input/output) The function value at X1, this value must be supplied
	DF1 : (output) The first derivative along the search direction at X1
	REPS : (input) Required relative accuracy 
	AEPS : (input) Required absolute accuracy,
		These criterion is only used to terminate line search under
		certain conditions and not generally applicable to the line
		search. These values should be same as what is used in function BFGS
	F : (input) Name of the function routine to calculate the function value
		and its derivatives
	V : (input/output) Array of length 3N. First N element specify
		the direction in which minimisation is required. The next
		N elements will contain the coordinates of the minimiser
		found by LINMIN. The last 3N elements will contain the gradient
		vector at the minimiser.
	XI : (input) Array of length N containing the coordinates of
		starting point for line search
	N : (input) Number of variables in the function to minimised
	NUM : (output) Integer variable to keep count of the number of
		function evaluations used so far
		
	Error status is returned by the value of the function LINMIN.
		0 value implies successful execution
		54 implies that LINMIN failed to find acceptable point
			but the function value at last point is less than
			that at beginning.
		55 implies that LINMIN failed to find acceptable point
			even though the interval has been reduced to required accuracy
			but the function value at last point is less than
			that at beginning.
		527 implies that iteration failed to find any point where
			the function value is less than the starting value
		528 implies that iteration failed to find any point where
			the function value is less than the starting value
			even though the interval has been reduced to required accuracy

	Function F(N,X,FX,G) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, FX is the
		function value at X and G is the gradient vector. X and G
		are arrays of length N.

	Required functions :  FLNM, F
*/

#include <math.h>

double flnm(void (*fcn) (int , double * , double * , double * ), double x,
	double *df, double v[], double x0[], int n, int *num);

int linmin(double *x1, double *x2, double *f1, double *df1, double reps,
	double aeps, void (*f) (int , double * , double * , double * ),
	double v[], double xi[], int n, int *num)

{
	int i,qb, nit=15;
	double f2,df2,x0,f0,dx1,df0,fc,df12,dx,r,r1,r2,xa,fa,dfa;
	double rho=0.01, sigma=0.1, t1=9.0, t2=0.9, t3=0.5;

/*	Select the bracketing phase */
	qb=0;
	f2=flnm(f,*x2,&df2,v,xi,n,num);
	dx1=(*x2)-(*x1);
	f0=(*f1); df0=(*df1); x0=(*x1);

	for(i=1; i<=nit; ++i) {
		fc=f0+df0*rho*(*x2-x0);
		if(fabs(df2)<=(-sigma*df0) && f2<=fc) {
/*	Found an acceptable point */
			*x1=(*x2); *f1=f2; *df1=df2;
			return 0;
		}

/*	Test for bracketing */
		if(qb==0) {
			if(f2>fc || f2>(*f1) || df2>=0.0) qb=1;
		}

/*	Hermite cubic interpolation */
		df12=(f2-(*f1))/dx1;
		r=2.*df2+(*df1)-3.*df12;
		r1=3.*df12-df2-(*df1); r1=r1*r1-(*df1)*df2;
		dx=0.0;
		if(r1>0.0) {
			if(r>=0.0) r1=sqrt(r1);
			else r1=-sqrt(r1);
			r=r+r1;
			dx=-df2*dx1/r;
		} else {
/*	try parabolic interpolation */
			r2=2.*(df12-df2);
			if(r2==0.0) dx=dx1*df2/r2;
		}

		if(qb==1) {
/*	Minimum is bracketed and hence improve on the bracket */
			if(dx<(-t2*dx1)) dx=-t2*dx1;
			if(dx>(-t3*dx1)) dx=-t3*dx1;
			xa=(*x2)+dx;
			fa=flnm(f,xa,&dfa,v,xi,n,num);
			fc=f0+df0*rho*(xa-x0);

			if(fabs(dfa)<=(-sigma*df0) && fa<=fc) {
/*	The new point is acceptable */
				*x1=xa; *f1=fa; *df1=dfa;
				return 0;
			} else if(fa>fc || fa>=(*f1) || dfa>0.0) {
				*x2=xa; f2=fa; df2=dfa;
			} else {
				*x1=xa; *f1=fa; *df1=dfa;
			}

			dx1=(*x2)-(*x1);
			r1=reps*fabs(*x2); if(aeps>r1) r1=aeps;
			if(fabs(dx1)<r1) {
/*	If the interval is too small, then quit */
				if(f2<=f0) {
/*	Accept the last point in any case */
					*x1=(*x2); *f1=f2; *df1=df2;
					return 55;
				} else return 528;
			}
		} else {
/*	Minimum hasn't been bracketed, choose the point further down. */
			if(dx< (*x2)-(*x1) && dx>0.0) dx=(*x2)-(*x1);
			if(dx> t1*(*x2-(*x1)) || dx<= 0.0) dx=t1*(*x2-(*x1));
			*x1=(*x2);
			*x2=(*x2)+dx;
			dx1=dx;
			*f1=f2;
			*df1=df2;
			f2=flnm(f,*x2,&df2,v,xi,n,num);
		}
	}

/*	Iteration has failed to find an acceptable point */
	if(f2<=f0) {
/*	accept this point if function value is smaller */
		*f1=f2; *x1=(*x2); *df1=df2;
		return 54;
	}
	else if((*f1)<=f0 && (*x1) != x0) return 54;
	else return 527;
}




/*	Function to calculate the function value and derivative
	as required for line search

	FCN : (input) Name of function to calculate the required function
	X : (input) Parameter along the line to specify the point where
		function evaluation is required
	DF : (output) First derivative of function along the line at X
	V : (input/output) Array of length 3N, first N elements specify the
		direction of line search. Next N elements will contain the
		coordinates of the point at which function is evaluated,
		while the last N elements contain the gradient vector at the point
	X0 : (input) Array of length N, containing the coordinates
		of starting point for line search
	N : (input) Number of variables in the function to be minimised
	NUM : (input/output) Integer variable to keep count of function evaluations

	Function FCN(N,X,FX,G) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, FX is the
		function value at X and G is the gradient vector. X and G
		are arrays of length N.

	Required functions : FCN
*/

double flnm(void (*fcn) (int , double * , double * , double * ), double x,
	double *df, double v[], double x0[], int n, int *num)

{
	int i,n2;
	double fl;

	*num=(*num)+1;
	n2=2*n;

/*	The coordinates of the required point */
	for(i=0; i<n; ++i) v[n+i]=x0[i]+v[i]*x;

	fcn(n,&v[n],&fl,&v[n2]);

/*	The first derivative along the search direction */
	*df=0.0;
	for(i=0; i<n; ++i) *df=(*df)+v[i]*v[n2+i];
	return fl;
}




/*	Function to calculate the Chi square function for a nonlinear
	least squares fit, for use with function BFGS.

	N : (input) Number of parameters to be fitted
	A : (input) Array of length N containing the parameter values
		at which the function needs to be evaluated
	F : (output) The function value at specified parameter values
	G : (output) Array of length N containing the gradient vector
		G[I] will contain dF/dA[I]

	The data points are passed through global variables NNLSQ, FXLSQ,
	XLSQ, EFLSQ, FX1LSQ, which must be initialised in the calling
	function before calling BFGS

	NNLSQ : (input) The number of data points in the table to be fitted
	FXLSQ : (input) Array of length NPL containing the function values
	XLSQ : (input) Array of length NPL containing the abscissas at which
		function values are available
	EFLSQ : (input) Array of length NPL containing the estimated errors in 
		function values, for use in deciding the weights for each point
		Although EFLSQ should contain the error estimate, but in many
		cases it is found that multiplying all values by a suitable
		constant can improve the convergence of BFGS dramatically
	FX1LSQ : (output) Array of length NPL containing the fitted value
		of the function at each tabular point

	The parameter NPL must be equal to the dimension of arrays as
	declared in the calling function.

	This function requires function FCN to calculate the required
	function which has to be fitted. There is no provision to pass on
	the name of this function and hence it must be changed explicitly
	to the required name. Function FCN(N,A,X,F,DF) must be supplied
	by the user. Here N is the number of parameters, A is an array
	of length N containing the values of parameters and X is the value
	of independent variable where the fitting function needs to be evaluated.
	F is the calculated function value and DF is an array of length N
	containing the calculated derivatives.

	Required functions : FCN
*/

#include <math.h>
#include <stdlib.h>


/* To pass on these parameters to nllsq, this must be included before
   the calling BFGS */

void fcn(int n, double a[], double x, double *fa, double df[]);

#define NPL 100
int NNLSQ;
double FXLSQ[NPL], XLSQ[NPL], EFLSQ[NPL], FX1LSQ[NPL];

void nllsq(int n, double a[], double *f, double g[])

{
	int i,j;
	double fa,r1;
	double *df;

/*     GENERATING THE FUNCTION FOR MINIMISATION */

	*f=0.0;
	for(i=0; i<n; ++i) g[i]=0.0;

	df=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<NNLSQ; ++i) {
		fcn(n,a,XLSQ[i],&fa,df);
		FX1LSQ[i]=fa;

/*     SUM OF SQUARED DIFFERENCE AND ITS GRADIENT */

		r1=(fa-FXLSQ[i])/EFLSQ[i];
		*f=(*f)+r1*r1;
		for(j=0; j<n; ++j) g[j]=g[j]+2.*r1*df[j]/EFLSQ[i];
	}
	free(df);
	return;
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

