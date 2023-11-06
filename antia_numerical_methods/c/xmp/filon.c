/*	Integrate an oscillatory functions using Filon's formula */

#include <stdio.h>
#include <math.h>

double fun(double x);
int filon(double *ri, double xl, double xu, double rk, int qsin, double reps,
	double aeps, double *dif, int *n, double (*fun) (double ));

main()
{
	int i,i1,j,nuse, id, iflg, ier,np;
	double xl, xu, rint, reps, aeps, df, b, rk;

/*	Exercise 6.21 : Integral I1 */

	aeps=1.e-18; reps=1.e-13;
	iflg=1;
	b=4.*acos(0.0);

	for(i1=0; i1<99; ++i1) {
		printf("type xl=lower limit,   rk=coef. in oscillatory factor\n");
		printf("                        (quits when rk=0)\n");
		scanf(" %le %le", &xl,  &rk);
		if(rk==0) return 0;

/*	Set the upper limit to 2*Pi */
		xu=b;
		i=filon(&rint,xl,xu,rk,iflg,reps,aeps,&df,&nuse,fun);
		printf(" ier = %d  no. of function evaluations = %d  interval = %e , %e \n", i,nuse,xl,xu);
		printf(" rk =  %e    integral = %e   error = %e  \n", rk,rint, df);
	}
	return;
}

/*  The integrand without the oscillatory factor cos(k*x) or sin(k*x) */
double fun(double x)

{
	return x*cos(x); 
}

 
/*	To calculate integrals with oscillatory integrand of the form
	FUN(x)*SIN(RK*x)  or  FUN(x)*COS(RK*x)

	RI : (output) Computed value of the integral
	XL : (input) The lower limit
	XU : (input) The upper limit
	RK : (input) Coefficient of x in the oscillatory term in the integrand
	QSIN : (input) Integer variable to specify the form of integrand
		If QSIN=1 the oscillatory factor is sin(RK*x)
		otherwise the oscillatory factor is cos(RK*x)
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) estimated (absolute) error achieved by the function
	N : (output) Number of function evaluations used
	FUN : (input) Name of the function to calculate the integrand
		excluding the oscillatory factor
		
	Error status is returned by the value of the function FILON.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy

		Function FUN(X) must be supplied by the user.
	
	Required functions : FUN
*/

#include <math.h>

int filon(double *ri, double xl, double xu, double rk, int qsin, double reps,
	double aeps, double *dif, int *n, double (*fun) (double ))
	
{
	int i,j, nmin=5, nmax=15;
	double r1,alpha,beta,gamma,t,h,fend,odd,even,x,x1,h2;
	double thc=0.005;


	if(qsin == 1) {
/*	For Sin(kx) */
		fend=fun(xl)*cos(rk*xl)-fun(xu)*cos(rk*xu);
		even=0.5*(fun(xl)*sin(rk*xl)+fun(xu)*sin(rk*xu));
	}
	else {
/*	For Cos(kx) */
		fend=fun(xu)*sin(rk*xu)-fun(xl)*sin(rk*xl);
		even=0.5*(fun(xl)*cos(rk*xl)+fun(xu)*cos(rk*xu));
	}

	odd=0.0; *ri=0.0; *dif=0.0;
	*n=2;
	h=(xu-xl);
	if(h == 0.0) return 0;

	for(i=0; i<nmax; ++i) {
		h=h/2.;
		even=even+odd;
		odd=0.0;
		x1=xl+h;
		h2=2.*h;

/*	Starting with 3 points subdivide the intervals into 2 until convergence */
		for(j=0; j<(*n)/2; ++j) {
			x=x1+h2*j;
			if(qsin==1) odd=odd+fun(x)*sin(rk*x);
			else odd=odd+fun(x)*cos(rk*x);
		}

		t=rk*h;
		if(fabs(t) > thc) {
/*	Use normal functions */
			alpha=(t*t+sin(t)*(t*cos(t)-2.*sin(t)))/(t*t*t);
			beta=2.*(t+cos(t)*(t*cos(t)-2.*sin(t)))/(t*t*t);
			gamma=4.*(sin(t)-t*cos(t))/(t*t*t);
		}
		else {
/*	Use Taylor series expansion */
			alpha=2.*t*t*t*(1.+t*t*(-1.+t*t/15.)/7.)/45.;
			beta=2.0/3.+2.*t*t*(1.0/15.+t*t*(-2.0/105.+t*t/567.));
			gamma=4.0/3.+t*t*(-2.0/15.+t*t*(1.0/210.-t*t/11340.));
		}
		r1=h*(alpha*fend+beta*even+gamma*odd);

		*dif=r1-(*ri);
		*ri=r1;
		if(i >= nmin) {
			t=reps*fabs(r1); if(aeps>t) t=aeps;
			if(fabs(*dif) < t) return 0;
		}
		*n=(*n)*2;
	}

	*n=(*n)/2;
	return 30;
}
