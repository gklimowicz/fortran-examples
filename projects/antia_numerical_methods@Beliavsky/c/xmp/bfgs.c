/*	To minimise a function of several variables using quasi-Newton method
	with BFGS formula */

#include <stdio.h>
#include <math.h>

void fun(int n, double x[], double *f, double df[]);
double flnm(void (*fcn) (int , double * , double * , double * ), double x,
	double *df, double v[], double x0[], int n, int *num);

int linmin(double *x1, double *x2, double *f1, double *df1, double reps,
	double aeps, void (*f) (int , double * , double * , double * ),
	double v[], double xi[], int n, int *num);
int bfgs(int n, double x[], double *f, double g[], double *h, int *num,
	double reps, double aeps, void fcn(int , double * , double * , double * ));


main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double x0[20],g[20], reps, aeps, h[2][2], fx;

/*	Example 8.4 : Rosenbrock's function */

	aeps=1.e-13; reps=1.e-11;
	np=2;
	for(i1=0; i1<99; ++i1) {
		printf("type initial guesses for x0[0] to x0[np-1] \n");
		printf("                   (quits when x0[0]<-100)\n");
		for(i=0; i<np; ++i) scanf(" %le", &x0[i]);
		if(x0[0]<-100) return 0;
		printf(" Initial guess to minimiser : \n");
		for(i=0; i<np; ++i) printf(" %e", x0[i]);
		printf("\n");

		i=bfgs(np,x0,&fx,g,&h[0][0],&nuse,reps,aeps,fun);
		printf(" ier = %d  no. of function evaluations = %d , minimum = %e \n minimiser = ", i,nuse,fx);
		for(i=0; i<np; ++i) printf(" %e", x0[i]);
		printf(" \n  gradient = ");
		for(i=0; i<np; ++i) printf(" %e", g[i]);
		printf(" \n  Inverse of Hessian matrix = ");
		for(i=0; i<np; ++i) {
			for(j=0; j<np; ++j) printf(" %e", h[i][j]);
			printf("\n");
		}
	}
	return;
}

void fun(int n, double x[], double *f, double g[])

{
	*f=(x[1]-x[0]*x[0]);
	g[0]=-400.*(*f)*x[0]-2.*(1.-x[0]);
	g[1]=200.0*(*f);
	*f=100.*(*f)*(*f)+(1.0-x[0])*(1.0-x[0])+1.0;
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
