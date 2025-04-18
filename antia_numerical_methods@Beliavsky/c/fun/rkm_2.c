/*	To solve initial value problems in ordinary differential equations
	using a second-order Runge Kutta method with adaptive
	step size control

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length N containing the initial
		values of variables. After execution it will contain the
		values at the last point where the integration
		has been successful.
	DY : (output) Array of length N containing the derivatives
		of Y at the last point
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input/output) Initial guess for the step size. After execution
		it will contain the step size used by the function
	T0 : (input/output) Initial value of independent variable t, at
		which initial values are specified. After execution it will
		be set to the point up to which integration has been successful
	TN : (input) The final value of t at which the solution is required.
		If integration is successful T0 will be set equal to TN.
		Intermediate values will not be preserved so if solution
		is required at intermediate points, TN must be set to first
		such value and multiple calls will be needed to calculate
		all required values. For each subsequent call only TN needs
		to be updated.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls local truncation error and hence
		actual error could be larger
	NSTEP : (output) Number of steps required to complete the integration
		Each step requires 4 or 5 calls to DIF (using RK2)
	NMAX : (input/output) Maximum number of steps to be used. If NMAX<=0
		it will be set to a default value of NMX=10000.
		
	Error status is returned by the value of the function RKM_2.
		0 value implies successful execution
		701 implies N<=0, no calculations are done
		721 implies that step-size has become smaller than REPS*|TN-T0|
		722 implies that step size is too small for arithmetic used
		723 implies that integration could not be completed in
			the specified number of steps.

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : RK2 (or RK4), DIF
*/	

#include <math.h>
#include <stdlib.h>

void rk2(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ));

int rkm_2(int n, double y[], double dy[], void dif(double , int , double * , double *),
	double *h, double *t0, double tn, double reps, int *nstep, int *nmax)

{
	int i, nmx=10000;
	double tstep,h2,t1,err,r1,r2;
/*	double e1=0.20, e2=0.25, ce=15.0, sfac=0.9, eps=1.e-30;
 For use with rk4 use th preceding statement instead of following one */
	double e1=0.33, e2=0.5, ce=3.0, sfac=0.9, eps=1.e-30; 
	double *wk3,*wk4,*wk5;

	if(n<=0) return 701;

	if(*nmax <=0) *nmax=nmx;
	*nstep=0;
	tstep=tn-(*t0);
	if(tstep==0.0) return 0;
/*	Adjust the initial value of H if needed */
	if(*h==0.0 || fabs(*h)>fabs(tstep)) *h=tstep;
	if((*h<0.0) == (tn>(*t0))) *h=-(*h);
/*	Calculate the derivatives at the initial point */
	dif(*t0,n,y,dy);
	wk3=(double *) calloc((size_t) n, sizeof(double));
	wk4=(double *) calloc((size_t) n, sizeof(double));
	wk5=(double *) calloc((size_t) n, sizeof(double));

/*	Loop for integration */
	while(*nstep<(*nmax)) {
/*	Use two steps of h/2 */
		++*nstep;
		h2=(*h)/2.0;
		rk2(n,*t0,y,dy,h2,wk4,dif);
		t1=(*t0)+h2;
		dif(t1,n,wk4,wk5);
		rk2(n,t1,wk4,wk5,h2,wk3,dif);

/*	Use single step of h */
		rk2(n,*t0,y,dy,*h,wk4,dif);

/*	Estimate the truncation error */
		err=0.0;
		for(i=0; i<n; ++i) {
			r2=fabs(y[i])+fabs(wk3[i]-y[i])+eps;
			r1=fabs(wk4[i]-wk3[i])/r2;
			if(r1>err) err=r1;
		}
		err=fabs(err*tstep/(ce*(*h)));

		if(t1==(*t0)) {free(wk5); free(wk4); free(wk3); return 722;}

		if(err<reps) {
/*	Integration at this step is successful, update T0 and Y */
			*t0=(*t0)+(*h);
			for(i=0; i<n; ++i) y[i]=wk3[i]-(wk4[i]-wk3[i])/ce;
			dif(*t0,n,y,dy);
			if(fabs((tn-(*t0))/tstep)<reps) {free(wk5); free(wk4); free(wk3); return 0;}

/*	Adjust the step size */
			if(err==0.0) *h=2.0*(*h);
			else *h=sfac*(*h)*pow(reps/err,e1);
			if((*t0+(*h)>tn) == ((*h)>0.0)) (*h)=tn-(*t0);
		}
		
/*	If the integration is not successful, try again with smaller step size */
		else  *h=sfac*(*h)*pow(reps/err,e2);
		
/*	Step size is too small then quit */
		if(fabs(*h/tstep)<reps) {free(wk5); free(wk4); free(wk3); return 721;}

	}

	free(wk5); free(wk4); free(wk3);
	return 723;
}
