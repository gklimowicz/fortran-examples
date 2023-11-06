/*	To generate the starting values for multistep methods in solution
	of initial value problem in ordinary differential equations.
	It is used by MSTEP.

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length 4N, the first N elements
		should contain the initial values of variables. During
		execution, the solution at the next three points will
		be stored in this array. 
	DY : (output) Array of length 4N, containing the derivatives
		of Y the points stored in Y.
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input/output) Initial guess for the step size. After execution
		it will contain the step size used by the function
	T : (input/output) Initial value of independent variable t, at
		which initial values are specified. After execution it will
		be set to T+3H if execution is successful.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls local truncation error and hence
		actual error could be larger
	IFLG : (input) Integer variable used as flag to decide the
		type of computation.
		If IFLG=0 the step size is adjusted to meet the specified accuracy.
		If IFLG=1 the step size is kept fixed.
	TSTEP : (input) The size of interval over which integration is requested
			It is used only for convergence check.
	NSTP : (input/output) Number of calls to DIF required since the starting
		of integration. The count is updated by the function.
	WK : Scratch array of length 2N used to pass on the value for modifier
		
	Error status is returned by the value of the function STRT4.
		0 value implies successful execution
		724 implies that step size becomes too small
		725 implies that iteration to generate starting values
			failed to converge

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : RK4, DIF
*/

#include <math.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ));

int strt4(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double *h, double *t, double reps, int iflg, double tstep, int *nstp,
	double *wk)

{
	int i,it, nit=10;
	double t0,h2,err,r1,r2, sfac=0.9, eps=1.0e-30;

	t0=(*t);
	for(it=1; it<=nit; ++it) {
		dif(t0,n,y,dy);
/*	Generate y_1 */
		rk4(n,t0,y,dy,*h,&y[n],dif);
		*t=t0+(*h);
		dif(*t,n,&y[n],&dy[n]);
/*	Generate y_2 */
		rk4(n,*t,&y[n],&dy[n],*h,&y[2*n],dif);

		h2=2.0*(*h);
/*	Calculate y_2 using double step */
		rk4(n,t0,y,dy,h2,&y[3*n],dif);
		*nstp=*nstp+11;

/*	Estimate the truncation error in y_2 */
		err=0.0;
		for(i=0; i<n; ++i) {
			r2=fabs(y[i])+fabs(y[i+n])+fabs(y[i+2*n])+eps;
			r1=fabs(y[i+2*n]-y[i+3*n])/r2;
			if(r1>err) err=r1;
		}
		err=err*tstep/(*h);

		if(err<=reps || iflg>0) {
/*	Accept the computed values */
			*t=(*t)+(*h);
			dif(*t,n,&y[2*n],&dy[2*n]);
/*		Generate y_3 */
			rk4(n,*t,&y[2*n],&dy[2*n],*h,&y[3*n],dif);
			*t=(*t)+(*h);
			dif(*t,n,&y[3*n],&dy[3*n]);
			*nstp=(*nstp)+5;
/*	Store Y_3 for use by modifier */
			for(i=0; i<n; ++i) wk[i+n]=y[i+3*n];
			return 0;
		}

		else {
/*	Reduce the step size */
			*h=sfac*(*h)*pow(reps/err,0.25);
			if(fabs(*h/tstep) < reps || (*t)==t0) return 724;
		}
	}

	return 725;
}
