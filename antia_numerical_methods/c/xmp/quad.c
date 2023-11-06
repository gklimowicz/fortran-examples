/*	Integration over a finite interval in one dimension, using Simpson's rule,
	Romberg integration, epsilon-algorithm, Gauss-Legendre formula or
	Adaptive integration */

#include <stdio.h>
#include <math.h>

double fun(double x);
int simson(double *ri, double xl, double xu, double reps, double aeps,
	double *dif, int *n, double (*fun)(double ));
int rombrg(double *ri, double a, double b, double gi[], double reps,
	double aeps, double *dif, int *n, double (*fun)(double ));
int epsiln(double *ri, double a, double b, double reps, double aeps,
	double *dif, int *n, double (*fun)(double ));
int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ));
int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));
int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double xl, xu, rint, reps, aeps, dif, fb[20];

/*	Example 6.1 */

	aeps=1.e-18; reps=1.e-13;
	nmax=0;
	for(i=0; i<20; ++i) fb[i]=0.0;
	for(i1=0; i1<99; ++i1) {
/*	np is the no. of abscissas in Gaussian rule to be used */
		printf("type np = 2/4/8/16/32, xl = lower limit, xu = upper limit\n");
		printf("             quits when xl=xu)\n");
		scanf(" %d %le %le", &np, &xl, &xu);
		if(xl==xu) return 0;

		i=simson(&rint,xl,xu,reps,aeps,&dif,&nuse,fun);
		printf(" ier = %d  n = %d  interval = %e , %e \n", i,nuse,xl,xu);
		printf("simson:   integral = %e   error = %e  \n", rint, dif);

/*	rombrg is used with default error exponents, which are suitable
	for smooth integrands */
		nuse=0;
		i=rombrg(&rint,xl,xu,fb,reps,aeps,&dif,&nuse,fun);
		printf("rombrg:   ier = %d  n = %d  integral = %e  error = %e \n", i,nuse,rint,dif);
		nuse=0;
		i=epsiln(&rint,xl,xu,reps,aeps,&dif,&nuse,fun);
		printf("epsiln:    ier = %d  n = %d  integral = %e  error = %e \n", i,nuse,rint,dif);
		i=gauss(&rint,xl,xu,&np,reps,aeps,&dif,&nuse,fun);
		printf(" %d pt Gauss formula:   ier = %d  n = %d  integral = %e  error = %e \n", np,i,nuse,rint,dif);
		i=adpint(&rint,xl,xu,reps,aeps,&dif,fun,&nuse,&nmax);
		printf("adpint:   ier = %d  n = %d  integral = %e  error = %e \n", i,nuse,rint,dif);
	}
	return;
}

/*  Specify the integrand */
double fun(double x)

{
	return sqrt(x); 
}

 

/*	To integrate a function over finite interval using adaptive control
	of step size

	RINT : (output) Calculated value of the integral
	XL : (input) The lower limit
	XU : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	NPT : (output) Number of function evaluations
	NMAX : (input/output) Maximum number of function evaluations to be tried
		If NMAX<=0 it is set to MAXPT (=100000)

	Error status is returned by the value of the function ADPINT.
		0 value implies successful execution
		31 implies specified accuracy was not achieved on
			at least one subinterval
		32 implies that this failure occurred more than IFMAX (=5) times
		325 implies that function failed to attain required
			accuracy using NMAX function evaluations
		In all cases DIF will contain the estimated accuracy

		Function F(X) must be supplied by the user.

	Required functions : KRONRD (or GAUS16), F
*/

#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));
int gaus16(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));


int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax)

{
	int i,j,k,ier,ifail,iu,np,q, ipmax=100, ifmax=5, maxpt=100000;
	double aepsl,rl,ru,rm,fint,dif0,r1,xu1[100];
 
	*rint=0.0; *dif=0.0;
	if(xl == xu) return 0;
	ifail=0; ier=0;
	if(*nmax <= 0) *nmax = maxpt;
	aepsl= aeps;
	*npt=0;
	rl=xl; ru=xu; iu=0;

	do {
/*	To evaluate the integral over [rl,ru] */
		i=kronrd(&fint,rl,ru,&dif0,&np,f);
/*		i=gaus16(&fint,rl,ru,&dif0,&np,f);  */
		*npt = (*npt)+np;
		rm=0.5*(rl+ru);
/*	q=.TRUE. if the interval cannot be divided further */
		q = ((iu>=ipmax) || (rm == rl) || (rm == ru));
		r1=fabs(fint)*reps; if(aepsl>r1) r1=aepsl;

		if( (dif0 < r1) || q) {
/*	Accept the value of fint if adequate convergence or if the interval
	cannot be subdivided further */
			*rint=(*rint)+fint;
			*dif=(*dif)+dif0;
			r1=fabs(*rint)*reps; if(aepsl>r1) r1=aepsl;

			if(q && (dif0>r1)) {
/*	Integration fails to converge on this subinterval. Go to the next subinterval */
				ier=31; ifail=ifail+1;
/*	If failure is frequent then adjust the convergence criterion. */
				if(ifail > ifmax) {ier=32; aepsl=(*dif)*0.5;}
			}

/*	If all subintervals are exhausted then return */
			if(iu <= 0) return ier;

/*	otherwise try next subinterval */
			rl=ru; ru=xu1[iu]; iu=iu-1;
		}
		else {

/*	Subdivide the current interval and try again */
			iu=iu+1;
			xu1[iu]=ru;
			ru=rm;
		}
	} while(*npt < (*nmax));


/*	If the number of function evaluations has exceeded the limit then
	try a last call to estimate the integral over the remaining interval */
	ru=xu;
	i=kronrd(&fint,rl,ru,&dif0,&np,f);
/*	i=gaus16(&fint,rl,ru,&dif0,&np,f); */
	*npt=(*npt)+np;
	*rint=(*rint)+fint;
	*dif=(*dif)+dif0;
	return 325;
}




/*	To integrate a function over finite interval using Epsilon algorithm

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) estimated (absolute) error achieved by the function
	N : (input/output) On input it should contain the number of function
		evaluations to be used for first estimate. If N<2 or N>NPT it
		is set to 2. After execution it will contain the number of
		function evaluations actually used by the function
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function EPSILN.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy
		33 implies that N>NPT (=100) in which case it is set to 2
		34 implies that at some stage denominator vanished while
			calculating epsilon table. This value is ignored.
		35 implies that roundoff error appears to be dominating
		
	Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>

int epsiln(double *ri, double a, double b, double reps, double aeps,
	double *dif, int *n, double (*fun)(double ))

{
	int i,j,nd,ier, nmin=3, nmax=13, npt=100;
	double h,ri1,dif0,dif1,s1,den,r1,t[13][13];

	ier=0;
	if((*n)<=1) *n=2;
	if((*n)>npt) {ier=33; *n=2;}

/*	Sum of the end points for trapezoidal rule */
	s1=0.5*(fun(a)+fun(b));
	nd=1; *ri=0.0; t[0][0]=0.0;

	for(i=0; i<nmax-1; ++i) {
		h=(b-a)/(*n-1);

		for(j=1; j<=(*n)-2; j=j+nd) s1=s1+fun(a+j*h);
/*	The trapezoidal rule approximation */
		t[i][1]=s1*h;
		t[i+1][0]=0.0;
		ri1=*ri;
		if(i>0) {*dif=fabs(t[i][1]-t[i-1][1]); *ri=t[i][1];}

/*	Construct the Epsilon table */
		for(j=2; j<=i+1; ++j) {
			den=t[i-j+2][j-1]-t[i-j+1][j-1];

/*	If denominator is zero set the error flag */
			if(den != 0.0) t[i-j+1][j]=t[i-j+2][j-2]+1./den;
			else {ier=34; t[i-j+1][j]=t[i-j+2][j-2];}
		}

		if(i>3) {

/*	DIF is the minimum difference between two rows of epsilon table */
			for(j=3; j<=i-1; j=j+2) {
				dif1=fabs(t[i-j+1][j]-t[i-j][j]);
				if(dif1 < (*dif)) {*dif=dif1; *ri=t[i-j+1][j];}
			}
		}

		nd=2;
		if(i>5 && *dif>dif0) {
/*	Roundoff error appears to be dominating, retain the previous value of RI */
			*ri=ri1;
			return 35;
		}
		if(i>nmin) {
			dif0=*dif;
			r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
			if(*dif < r1) return ier;
		}
		*n=2*(*n)-1;
	}

	*n=(*n+1)/2;
	return 30;
}



/*	To integrate a function over finite interval using Gauss-Legendre formulas

	RINT : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	NP : (input/output) The Gauss-Legendre formula to be used, (NP=2,4,8,16,32)
		For other values of NP it will be set to 8.
		It will use composite NP-point formula
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	NPT : (output) Number of function evaluations used by the function
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function GAUSS.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy
		36 implies that NP was not 2,4,8,16 or 32. In which case
			it is set to 8.

		FUNCTION FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>

int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ))

{
	int i,j,k,n,no,np2,ier;
	double r1,dx,a1,s1,at,bt, nmax=9;

/*	Weights and abscissas for Gauss-Legendre quadrature.
	For N-point formula W[K]=W[N-K-1] and X[K]=-X[N-K-1]
		For K=0,1,...,N/2-1. Hence only half points are tabulated.
	For 2-point W[0]; 4-point W[1], W[2]; 8-point W[3],...,W[6];
	16-point W[7],...,W[14]; 32-point W[15],...,W[30] are the
	weights corresponding to abscissas X[I].
*/

	double w[31]={1.0e0,
      0.34785484513745385737e0, 0.65214515486254614263e0,
      0.10122853629037625915e0, 0.22238103445337447054e0,
      0.31370664587788728734e0, 0.36268378337836198297e0,
      0.02715245941175409485e0, 0.06225352393864789286e0,
      0.09515851168249278481e0, 0.12462897125553387205e0,
      0.14959598881657673208e0, 0.16915651939500253819e0,
      0.18260341504492358887e0, 0.18945061045506849629e0,
      0.00701861000947009660e0, 0.01627439473090567061e0,
      0.02539206530926205945e0, 0.03427386291302143310e0,
      0.04283589802222668066e0, 0.05099805926237617620e0,
      0.05868409347853554714e0, 0.06582222277636184684e0,
      0.07234579410884850623e0, 0.07819389578707030647e0,
      0.08331192422694675522e0, 0.08765209300440381114e0,
      0.09117387869576388471e0, 0.09384439908080456564e0,
      0.09563872007927485942e0, 0.09654008851472780057e0};

	double x[31] = {0.57735026918962576451e0,
      0.86113631159405257522e0, 0.33998104358485626480e0,
      0.96028985649753623168e0, 0.79666647741362673959e0,
      0.52553240991632898582e0, 0.18343464249564980494e0,
      0.98940093499164993260e0, 0.94457502307323257608e0,
      0.86563120238783174388e0, 0.75540440835500303390e0,
      0.61787624440264374845e0, 0.45801677765722738634e0,
      0.28160355077925891323e0, 0.09501250983763744019e0,
      0.99726386184948156354e0, 0.98561151154526833540e0,
      0.96476225558750643077e0, 0.93490607593773968917e0,
      0.89632115576605212397e0, 0.84936761373256997013e0,
      0.79448379596794240696e0, 0.73218211874028968039e0,
      0.66304426693021520098e0, 0.58771575724076232904e0,
      0.50689990893222939002e0, 0.42135127613063534536e0,
      0.33186860228212764978e0, 0.23928736225213707454e0,
      0.14447196158279649349e0, 0.04830766568773831623e0};

	n=1;
	dx=b-a;
	ier=0;
	*rint=0.0;
	*npt=0;

	if(*np==2) no=0;
	else if(*np==4) no=1;
	else if(*np==8) no=3;
	else if(*np==16) no=7;
	else if(*np==32) no=15;
	else {
/*	If NP-point formula is not available use NP=8 */
		*np=8; no=3; ier=36;
	}
/*	X[NO],...,X[NP2] are the abscissas for the formula */
	np2=no+(*np)/2-1;

/*	Subdivide the interval until convergence */
	for(i=0; i<nmax; ++i) {
		r1=0.0;
		for(j=0; j<n; ++j) {
			a1=a+j*dx;
			at=dx/2.;
			bt=a1+dx/2.;

/*	To reduce roundoff errors sum over each subinterval is evaluated separately */
			s1=0.0;
			for(k=no;  k<=np2; ++k) s1=s1+w[k]*(fun(at*x[k]+bt)+fun(bt-at*x[k]));
			r1=r1+s1;
		}
		r1=r1*dx/2.;

/*	convergence check */
		*dif=r1-(*rint);
		*rint=r1;
		*npt=(*npt)+n*(*np);
		s1=reps*fabs(*rint); if(aeps>s1) s1=aeps;
		if(i>0 && fabs(*dif) < s1) return ier;
		dx=dx/2.;
		n=n*2;
	}
	return 30;
}



/*	To integrate a function over a finite interval using Gauss-Kronrod formula
	For use with ADPINT

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	DIF : (output) estimated (absolute) error achieved by the function
	N : (output) Number of function evaluations used
	F : (input) Name of the function to calculate the integrand

	Function F(X) must be supplied by the user

	The returned value is always zero.

	Required functions : F
*/

#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ))

{
	int k;
	double at,bt,fbt,r1,f1,f2;

/*	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
	WK7 are the weights for these points in Kronrod formula
	WK15 and AK15 are the weights and abscissas for the remaining points
	in Kronrod formula.
	Because of symmetry only half the points are given.
*/

	double w7[4] = {0.12948496616886969327e0, 0.27970539148927666790e0,
                        0.38183005050511894495e0, 0.41795918367346938775e0};
	double a7[4] = {0.94910791234275852452e0, 0.74153118559939443986e0,
                        0.40584515137739716690e0, 0.0};
	double wk7[4] = {0.06309209262997855329e0, 0.14065325971552591874e0,
                         0.19035057806478540991e0, 0.20948214108472782801e0};
	double wk15[4] = {0.02293532201052922496e0, 0.10479001032225018383e0,
                          0.16900472663926790282e0, 0.20443294007529889241e0};
	double ak15[4] = {0.99145537112081263920e0, 0.86486442335976907278e0,
                          0.58608723546769113029e0, 0.20778495500789846760e0};

	at=(b-a)/2.;
	bt=(b+a)/2.;
	fbt=f(bt);
	r1=w7[3]*fbt;
	*ri=wk7[3]*fbt;
	for(k=0; k<3; ++k) {
		f1=f(bt+at*a7[k]);
		f2=f(bt-at*a7[k]);
/*	7-point Gauss-Legendre formula */
		r1=r1+w7[k]*(f1+f2);
/*	15-point Kronrod formula */
		*ri=(*ri)+wk7[k]*(f1+f2);
	}

	for(k=0; k<4; ++k) *ri=(*ri)+wk15[k]*(f(bt+at*ak15[k]) + f(bt-at*ak15[k]));

	*ri=(*ri)*at;
	r1=r1*at;
	*dif=fabs(*ri-r1);
	*n=15;
	return 0;
}


/*	To integrate a function over finite interval using Romberg integration

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	GI : (input/output) Array of length NMAX (=13), containing
		the expected values of exponents \gamma_i in error expansion
		If GI[I]<=0 it will be set to 2I+2, the correct value for
		a smooth function
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) Estimated (absolute) error achieved by the function
	N : (input/output) On input it should contain the number of function evaluations
		to be used for first estimate. If N<2 or N>NPT it is set to 2.
		After execution it will contain the number of function
		evaluations actually used by ROMBRG.
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function ROMBRG.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy
		33 implies that N>NPT (=100) in which case it is set to 2

	Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/		

#include <math.h>

int rombrg(double *ri, double a, double b, double gi[], double reps,
	double aeps, double *dif, int *n, double (*fun)(double ))

{
	int i,j,nd,ier, nmin=3, nmax=13, npt=100;
	double dif1,fj,r1,s1,h, t[13][13];

	for(i=0; i<nmax; ++i) {
		if(gi[i] <= 0.0) gi[i]=2*i+2;
	}

	ier=0;
	if((*n)<=1) *n=2;
	if((*n)>npt) {ier=33; *n=2;}

/*	Contribution from the end points */
	s1=0.5*(fun(a)+fun(b));
/*	First time use all points */
	nd=1; *dif=0.0;

	for(i=0; i<nmax; ++i) {
		h=(b-a)/(*n-1);

/*	Add new points to the sum */
		for(j=1; j<=(*n)-2; j=j+nd) s1=s1+fun(a+j*h);
/*	The trapezoidal rule approximation */
		t[i][0]=s1*h;

/*	The Richardson's extrapolation */
		for(j=0; j<i-1; ++j) {
			fj=pow(2.,gi[j]);
			t[i][j+1]=t[i][j]+(t[i][j]-t[i-1][j])/(fj-1.0);
			dif1=fabs(t[i][j]-t[i-1][j]);
/*	Find the minimum difference between the last two rows of T-table */
			if(dif1 < (*dif) || j == 0) {*dif=dif1; *ri=t[i][j];}
		}

/*	On second and subsequent pass add only new points to the sum */
		nd=2;
		if(i>=nmin) {
			r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
			if(*dif < r1) return ier;
		}
		*n=2*(*n)-1;
	}

	*n=((*n)+1)/2;
	return 30;
}




/*	To integrate a function over finite interval using Simpson's rule

	RI : (output) Calculated value of the integral
	XL : (input) The lower limit
	XU : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	DIF : (output) estimated (absolute) error achieved by the function
	N : (output) Number of function evaluations used by SIMSON
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function SIMSON.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy

	Function FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>


int simson(double *ri, double xl, double xu, double reps, double aeps,
	double *dif, int *n, double (*fun)(double ))

{
	int i,j,n2, nmin=3, nmax=13;
	double r1,r2,fend,odd,even,x,x1,h2,h;

	fend=fun(xl)+fun(xu);
	even=0.0; odd=0.0;
	*ri=0.0; *dif=0.0;

	*n=2;
/*	starting with 2+1 points, subdivide the intervals into 2 until convergence */
	h=(xu-xl);
	if(h==0) return 0;

	for(i=1; i<= nmax; ++i) {
		h=h/2.;
		even=even+odd;
		odd=0.0;
		x1=xl+h;
		n2=(*n)/2;
		h2=2.*h;

		for(j=0; j<n2; ++j) {
			x=x1+h2*j;
			odd=odd+fun(x);
		}

/*	Estimate for the integral */
		r1=(fend+4.*odd+2.*even)*h/3.;
		*dif=r1-*ri;
		*ri=r1;
/*	To avoid spurious convergence in first few trials skip the convergence test */
		if(i>nmin) {
			r2=reps*fabs(r1); if(aeps>r2) r2=aeps;
			if(fabs(*dif) < r2) return 0;
		}
		*n = (*n)*2;
	}

	*n = (*n)/2;
	return 30;
}
