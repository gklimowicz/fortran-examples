/*	To solve two point boundary value problem in ordinary differential equations
	using shooting technique */

#include <stdio.h>
#include <math.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ));
int rkm(int n, double y[], double dy[], void dif(double , int , double * , double *),
	double *h, double *t0, double tn, double reps, int *nstep, int *nmax);
void dif(double t, int n, double y[], double dy[]);
int brent(double *a, double *b, double *x, double reps, double aeps,
	double (*f) (double ));
double fs(double xs);


main()
{
	int i,i1,j,n,m, id, iflg, ier,np,nmax;
	double hh, x[100], dx[100],tn,x1,xl,xu,reps,aeps,t0;

/*	Example 12.11 */

	id=4; n=2; nmax=100000; reps=1.e-7, aeps=1.e-8;

/*	the initial value x for f'(0) is determined by solving a nonlinear
	equation defined by the second boundary condition
	Supply the expected limits on this value */

	for(i1=0; i1<99; ++i1) {
		printf("type xl=lower limit,   xu=upper limit on x\n");
		printf("                      (quits when xl=xu)\n");
		scanf(" %le %le", &xl, &xu);
		if(xl==xu) return 0;
		printf(" lower limit on x = %e    upper limit on x = %e\n",xl,xu);

		i=brent(&xl,&xu,&x1,reps,aeps,fs);
		printf(" ier = %d    initial slope =  %e \n",i,x1);

/*	Once the initial values are known generate the solution at the required points */

		t0=0; hh=0.01;
		x[0]=0; x[1]=x1;
		for(i=1; i<=5; ++i) {
			tn=i*0.2;
			j=rkm(n,x,dx,dif,&hh,&t0,tn,reps,&np,&nmax);
			printf("  ier = %d    t = %e   solution = ",j,tn);
			for(j=0; j<n; ++j) printf(" %e ",x[j]);
			printf(" \n");
		}

	}
	return;
}

void dif(double t, int n, double y[], double dy[])

{
	double al=5.0;

	dy[0]=y[1];
	dy[1]=al*sinh(al*y[0]);
	return;
}



/*	Real zero of a given function using Brent's method

	A : (input/output) Lower limit of interval where zero is located
	B : (input/output) Upper limit of interval where zero is located
		The limits A, B are updated by the function
	X : (output) Computed value of the zero
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))
	F : (input) Name of the function routine to calculate the function

	Error status is returned by the value of the function BRENT.
		0 value implies successful execution
		427 implies that function has same sign at x=a and b
			in which case no calculations are done
		428 implies that iteration failed to converge to specified
			accuracy

		Function F(X) must be supplied by the user.

	Required functions : F
*/

#include <math.h>

int brent(double *a, double *b, double *x, double reps, double aeps,
	double (*f) (double ))

{

	int i,nit=75;
	double fa,fb,c,fc,d,e,dx,r1,r2,r3,p,p1,eps;

	fa=f(*a); fb=f(*b);
/*	If the function has the same sign at both the end points, then quit */
	if( (fa>0.0) == (fb>0.0)) return 427;

	c=(*a); fc=fa;
	d=(*b)-(*a); e=d;

/*	The iteration loop */
	for(i=1; i<=nit; ++i) {
		if(fabs(fc) < fabs(fb)) {
/*	If f(c) < f(b), then interchange b and c */
			*a=(*b); *b=c; c=(*a);
			fa=fb; fb=fc; fc=fa;
		}

		eps=reps*fabs(*b); if(aeps>eps) eps=aeps;
		dx=0.5*(c-(*b));
		if(fabs(dx)<eps || fb==0.0) {*x=(*b); return 0;}

/*	If the previous step is too small or if f(a)<=f(b), perform bisection */
		if(fabs(e)<eps || fabs(fa)<=fabs(fb)) {d=dx; e=dx;}
		else {
			r3=fb/fa;
			if((*a)==c) {
/*	Try linear interpolation */
				p=2.*dx*r3;
				p1=1.-r3;
			}
			else {
/*	Try inverse quadratic interpolation */
				r1=fa/fc; r2=fb/fc;
				p=r3*(2.*dx*r1*(r1-r2)-(*b-(*a))*(r2-1.0));
				p1=(r1-1.0)*(r2-1.0)*(r3-1.0);
			}

			if(p>0.0) p1=-p1;
			else p=-p;

			e=d;
			if(2.0*p < (3.*dx*p1-fabs(eps*p1)) && p<fabs(0.5*e*p1)) d=p/p1;
			else {d=dx; e=dx;}
		}

		*a=(*b); fa=fb;
		if(fabs(d)>eps) *b=(*b)+d;
		else {
			if(dx>=0.0) *b=(*b)+eps;
			else *b=(*b)-eps;
		}
		fb=f(*b);
		if((fb>0.0) == (fc>0.0)) {
/*	If f(b) and f(c) have the same sign, then replace c by a */
			c=(*a); fc=fa;
			d=(*b)-(*a); e=d;
		}
	}

	*x=(*b);
	return 428;
}



/*	To perform one step of integration of ordinary differential equations
	using a fourth-order Runge Kutta method 

	N : (input) Number of first order differential equations to be solved
	T : (input) Initial value of independent variable t, at
		which initial values are specified. This value is not updated.
	Y0 : (input) Array of length N containing the initial
		values of variables.
	DY0 : (input) Array of length N containing the derivatives
		of Y at the initial point Y0
	H : (input) The step size to be used for integration
	Y1 : (output) Array of length N containing the solution at t=T+H
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : DIF
*/	

#include <math.h>
#include <stdlib.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ))

{
	int i;
	double h2,t1;
	double *wk1, *wk2;

	h2=h/2.;
	wk1=(double *) calloc((size_t) n, sizeof(double));
	wk2=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i) {
		y1[i]=h2*dy0[i];
/*	y_0+0.5k_1 */
		wk1[i]=y0[i]+y1[i];
	}

	t1=t+h2;
	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i) {
		y1[i]=y1[i]+h*wk2[i];
/*	y_0+0.5k_2 */
		wk1[i]=y0[i]+h2*wk2[i];
	}

	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i) {
		y1[i]=y1[i]+h*wk2[i];
/*	y_0+k_3 */
		wk1[i]=y0[i]+h*wk2[i];
	}

	t1=t+h;
	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i)  y1[i]=y0[i]+(y1[i]+h2*wk2[i])/3.0;
	free(wk2); free(wk1);
	return;
}
	
	


/*	To solve initial value problems in ordinary differential equations
	using a second or fourth-order Runge Kutta method with adaptive
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
		Each step requires 10 or 11 calls to DIF (using RK4)
	NMAX : (input/output) Maximum number of steps to be used. If NMAX<=0
		it will be set to a default value of NMX=10000.
		
	Error status is returned by the value of the function RKM.
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

	Required functions : RK4 (or RK2), DIF
*/	

#include <math.h>
#include <stdlib.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ));

int rkm(int n, double y[], double dy[], void dif(double , int , double * , double *),
	double *h, double *t0, double tn, double reps, int *nstep, int *nmax)

{
	int i, nmx=10000;
	double tstep,h2,t1,err,r1,r2;
	double e1=0.20, e2=0.25, ce=15.0, sfac=0.9, eps=1.e-30;
/* For use with rk2 use th following statement instead of preceding one 
	double e1=0.33, e2=0.5, ce=3.0, sfac=0.9, eps=1.e-30; */
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
		rk4(n,*t0,y,dy,h2,wk4,dif);
		t1=(*t0)+h2;
		dif(t1,n,wk4,wk5);
		rk4(n,t1,wk4,wk5,h2,wk3,dif);

/*	Use single step of h */
		rk4(n,*t0,y,dy,*h,wk4,dif);

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



/*	To specify the second boundary condition whose zero will give
	the required initial conditions */

double fs(double xs)

{
	int n,np,nmax,i;
	double t0,tn,x[100],dx[100],h,reps;

	t0=0; tn=1.0; h=0.01;
	x[0]=0; x[1]=xs;
	n=2; nmax=10000;
	reps=1.e-7;

	i=rkm(n,x,dx,dif,&h,&t0,tn,reps,&np,&nmax);

	if(i>0) return 1000.0;
	else return x[0]-1.0;
}

