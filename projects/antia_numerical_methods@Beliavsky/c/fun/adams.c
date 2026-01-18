/*	To perform one step of solution of initial value problems in ordinary
	differential equations using fourth order Adams-Bashforth-Moulton
	predictor corrector method. It is called by function MSTEP.

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length 7N containing the solution
		at last four points. After execution new point will be added.
	DY : (input/output) Array of length 7N containing the derivatives
		of Y at the points in Y
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input) Step length to be used.
	T : (input) Value of independent variable t, at which the solution
		is required.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls error in iteration on corrector.
	NSTP : (input/output) Number of calls to DIF required so far.
		The count is updated by the function.
	IJ : (input) The index j+1 in corrector formula
	IJM1 : (input) The index j in corrector formula
	IJM2 : (input) The index j-1 in corrector formula
	IJM3 : (input) The index j-2 in corrector formula
	IJM4 : (input) The index j-3 in corrector formula
	WK : Array of length 2N used to pass the predicted values

	Error status is returned by the value of the function ADAMS.
		0 value implies successful execution
		729 implies that iteration on corrector failed to converge

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : DIF
*/
	
#include <math.h>
int adams(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
	int ijm3, int ijm4, double *wk)

{
	int j,k, nit=10;
	double err,t1,t2,r1, cfac=0.2, eps=1.e-30;

	for(j=0; j<n; ++j) {
/*	The predictor */
		t1=y[j+ijm1*n]+h*(55*dy[j+ijm1*n]-59*dy[j+ijm2*n]+37*dy[j+ijm3*n]-9*dy[j+ijm4*n])/24.0;
/*	Modifier */
		y[j+ij*n]=t1+251.0*(y[j+ijm1*n]-wk[j+n])/270.0;
		wk[j]=t1;
	}

/*	Iteration on corrector */
	for(j=1; j<=nit; ++j) {
		*nstp=(*nstp)+1;
		dif(t,n,&y[ij*n],&dy[ij*n]);
		err=0.0;
		for(k=0; k<n; ++k) {
/*	The corrector */
			t2=y[k+ijm1*n]+h*(9*dy[k+ij*n]+19*dy[k+ijm1*n]-5*dy[k+ijm2*n]+dy[k+ijm3*n])/24.0;
			r1=fabs(t2)+fabs(t2-y[k+ijm1*n])+eps;
			t1=fabs((y[k+ij*n]-t2)/r1);
			if(t1>err) err=t1;
			y[k+ij*n]=t2;
		}

/*	The convergence test */
		if(err<cfac*reps) return 0;
	}

	return 729;
}
