/*	Differentiation using h ---> 0  extrapolation */

#include <stdio.h>
#include <math.h>
double drvt(double a, int id, double hh0, double reps, double aeps,
	double (*f)(double ), int *ier);

double fun(double x);

main()
{
	int i,i1,j, id,  ier;
	double xb, h, reps, aeps, df, fb[20];

/*	Example 5.2 : Derivatives of Exp(x) */

	aeps=1.e-6; reps=1.e-5;
	for(i1=0; i1<99; ++i1) {
		printf("type  id = order of derivative, x, h = initial spacing\n");
		printf("                              (quits when id<=0)\n");
		scanf("%d %le %le", &id, &xb, &h);
		if(id<=0) return 0;

		df=drvt(xb, id, h, reps, aeps, fun, &ier);
		printf(" ier = %d   id = %d   x = %e  df = %e \n", ier,id,xb,df);
		printf("  initial spacing = %e    exact derivative = %e\n",h,exp(xb));
	}
	return;
}

/*  The given function  */

double fun(double x)

{
	return exp(x);
}

 


/*	To calculate derivative of a function

	A : (input) Value of x at which the derivative needs to be calculated
	ID : (input) Order of derivative required, ID may be 1,2,3 or 4
	HH0 : (input) Initial value of step length to be used.
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy.
		Calculations will continue till estimated error is less
		than max(AEPS, REPS*fabs(DRVT))
	F : (input) The name of function routine to calculate the given function
		Function F(X) must be supplied by the user.
	IER : (output) Error parameter, IER=0 for successful execution
		IER=28 implies required accuracy is not achieved
		IER=29 implies roundoff error appears to be dominating
		IER=208 implies ID<1 or ID>4, no calculations are done.
		
	DRVT is the calculated derivative at x=A

	Function F(X) must be supplied by the user

	Required functions : F
*/

#include <math.h>

double drvt(double a, int id, double hh0, double reps, double aeps,
	double (*f)(double ), int *ier)

{
	int i,j,ju, nmax=15, ncmax=6;
	double df, e1, err, err1, hh, h[15], t[15][6], hfac=1.5;

/*	Set spacing to the specified initial value */

	hh=hh0; err=0.0; *ier=0;

	for(i=0; i<nmax; ++i) {
		h[i]=hh;

		if(id == 1) t[i][0] = (f(a+hh)-f(a-hh))/(2.*hh);
		else if(id == 2) t[i][0] = (f(a-hh)-2.*f(a)+f(a+hh))/(hh*hh);
		else if(id == 3) t[i][0] = (f(a+2*hh)-f(a-2*hh)-2.*(f(a+hh)-f(a-hh)))/(2.*hh*hh*hh);
		else if(id == 4)
			t[i][0] = (f(a+2*hh)+f(a-2*hh)-4.*(f(a+hh)+f(a-hh))+6*f(a))/(hh*hh*hh*hh);
		else {*ier=208; return 0.0;}

/*	Reduce spacing by a factor of hfac */

		hh=hh/hfac;
		if(i>0) {
/*	Calculate at most ncmax columns of the T-table */
			ju=i; if(ju>ncmax-1) ju=ncmax-1;
			for(j=1; j<=ju; ++j) t[i][j]=t[i][j-1]+(t[i][j-1]-t[i-1][j-1])/
				(h[i-j]*h[i-j]/(h[i]*h[i])-1.);

			df=t[i][ju];
			err1=err;
			j=i-1; if(j>ncmax-1) j=ncmax-1;
			e1=fabs(df-t[i-1][j]);
			err=fabs(df-t[i][ju-1]); if(e1<err) err=e1;
			e1=fabs(df)*reps; if(e1<aeps) e1=aeps;
			if(err<e1) return df;

			if(i>4 && err>2.*err1) {
/*	Roundoff error dominates */
				*ier=29;
				return t[i-1][j];
			}
		}
	}
	*ier=28;
	return df;
}

