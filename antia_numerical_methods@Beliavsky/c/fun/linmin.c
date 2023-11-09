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
