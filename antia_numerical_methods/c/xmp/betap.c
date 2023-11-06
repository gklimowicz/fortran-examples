/*	Calculate the incomplete beta function and incomplete gamma function */


#include <stdio.h>
#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));
int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax);
double betai(double a, double b, double x);
double betser(double a, double b, double x);
double betcon(double a, double b, double x);
double betcon1(double a, double b, double x);
double fbeta(double x);
double betap(double a, double b, double x);
double gamma(double xg);
double gammaln(double xg);
double gammap(double a,double x);

	double AA, BB;

main()
{
	int i;
	double a,b,x,f1,f2;

	for(i=0; i<99; ++i){
		printf("Type x, a, b  (quits when x<0)\n");
		scanf("%le %le %le",&x,&a,&b);
		if(x<0) return;
		f1=gammap(a,x); f2=betap(a,b,x);
		printf("x = %le  a = %le  b = %le \n",x,a,b);
		printf("Pgamma(a,x) = %le    I_x(a,b) = %le\n",f1,f2);
	}
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


/*	To calculate the incomplete Beta function I_x(a,b) using
	the integral, called by betap

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln, adpint, kronrd, fbeta
*/

#include <math.h>

	double fbeta(double x);
	int adpint(double *rint, double xl, double xu, double reps, double aeps,
		double *dif, double f1(double x), int *n, int *nmax);
	int kronrd(double *ri, double a, double b, double *dif, int *n,
		double f1(double x));

double betai(double a, double b, double x)
{
	double xl,xu,reps,aeps,b1,rint,dif,f;
	int i,npt,nmax;

	AA=a; BB=b; xl=0; xu=x; reps=1e-14; aeps=1e-280;
	i=adpint(&rint,xl,xu,reps,aeps,&dif,fbeta,&npt,&nmax);
	b1=log(rint)+gammaln(a+b)-gammaln(a)-gammaln(b);
	f=exp(b1);
	if(f>1) f=1;
	return f;
}
	
/*	To calculate the incomplete Beta function I_x(a,b) for
	positive real arguments
	It returns a value of -1 for a<=0 or b<=0 or x<0 or x>1

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln, betai, betser, betcon, betcon1,
			    adpint, kronrd, fbeta
*/

#include <math.h>
#include <stdio.h>

	double betai(double a, double b, double x);
	double betser(double a, double b, double x);
	double betcon(double a, double b, double x);
	double betcon1(double a, double b, double x);
	double fbeta(double x);
	int adpint(double *rint, double xl, double xu, double reps, double aeps,
		double *dif, double f1(double x), int *n, int *nmax);
	int kronrd(double *ri, double a, double b, double *dif, int *n,
		double f1(double x));


double betap(double a, double b, double x)
{
	double amax,amin,betal,f;
	
	if(a<=0 || b<=0 || x<0 || x>1) return -1.0;
	amax=a; if(b>a) amax=b;
	amin=a; if(b<a) amin=b;
	betal=(a+b)*log(a+b)-a*log(a)-b*log(b);
	if(amax<=30)
	{
		if(x<=0.5) f=betser(a,b,x);
		else f=1-betser(b,a,1-x);
	}
	else if(b<=20 && x<=0.71) f=betser(a,b,x);
	else if(a<=20 && x>=0.3) f=1-betser(b,a,1-x);
	else if(b<=50 && x<=0.35) f=betser(a,b,x);
	else if(a<=50 && x>=0.65) f=1-betser(b,a,1-x);
	else if(b<=100 && x<=0.18) f=betser(a,b,x);
	else if(a<=100 && x>=0.82) f=1-betser(b,a,1-x);
	else if(b<=180 && x<=0.1) f=betser(a,b,x);
	else if(a<=180 && x>=0.9) f=1-betser(a,b,x);
	else if(x<0.5)
	{
		if(a<2) f=betcon(a,b,x);
		else if(betal>700) f=betcon1(a,b,x);
		else f=betai(a,b,x);
	}
	else
	{
		if(b<2) f=1-betcon(b,a,1-x);
		else if(betal>700) f=1-betcon1(b,a,1-x);
		else f=1-betai(b,a,1-x);
	}
	if(f<0 || f>1) printf("error in evaluating incomplete beta function at a = %e  b= %e  x= %e  betap= %e\n",a,b,x,f);
	return f;
}


/*	To calculate the incomplete Beta function I_x(a,b) using
	the continued fraction (modified form), called by betap

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln
*/


#include <math.h>
#include <stdio.h>

double gammaln(double x);

double betcon(double a, double b, double x)
{
	double c[500],d[500],a1,b1,c1,c2,d1,d2,d3;
	int i;

	c[0]=a; d[0]=a*(1-x*(a+b)/(a+1));
	b1=1+x*(b-1)/(a+1)-x*(a+1)*(a+b+1)/(a+3);
	a1=a*x*x*(b-1)*(a+b)/pow(a+1,2.0);
	c[1]=a*b1; d[1]=b1*d[0]+a1; c2=c[1]/d[1]; c1=0.0; i=1;

	while(i<500 && fabs(c1-c2)>1.e-12)
	{
		i=i+1; c1=c2;
		d1=-x*(a+i-1)*(a+b+i-1)/((a+2*i-2)*(a+2*i-1));
		d3=-x*(a+i)*(a+b+i)/((a+2*i)*(a+2*i+1));
		d2=x*i*(b-i)/((a+2*i-1)*(a+2*i));
		c[i]=c[i-1]*(a+2*i)*(1+d2+d3)-c[i-2]*(a+2*i)*(a+2*i-2)*d2*d1;
		d[i]=d[i-1]*(a+2*i)*(1+d2+d3)-d[i-2]*(a+2*i)*(a+2*i-2)*d2*d1;
		if(fabs(c[i])>1.e200)
/*	Scale the numerator and denominator to prevent underflow/overflow */
		{
			c[i]=c[i]/1e200; c[i-1]=c[i-1]/1e200;
			d[i]=d[i]/1e200; d[i-1]=d[i-1]/1e200;
		}
		if(fabs(c[i])<1.e-200)
		{
			c[i]=c[i]*1e200; c[i-1]=c[i-1]*1e200;
			d[i]=d[i]*1e200; d[i-1]=d[i-1]*1e200;
		}
		c2=c[i]/d[i];
	}
	if(c2<0.0) printf(" ** Roundoff error while evaluating the modified continued fraction for incomplete Beta function at\n a= %e  b=%e  x= %e cont. frac.= %e\n",a,b,x,c2);
	b1=a*log(x)+b*log(1-x)+log(c2)-log(a);
	b1=b1+gammaln(a+b)-gammaln(a)-gammaln(b);
	return exp(b1);
}

/*	To calculate the incomplete Beta function I_x(a,b) using
	the continued fraction, called by betap

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln
*/


#include <math.h>
#include <stdio.h>

double gammaln(double x);

double betcon1(double a, double b, double x)
{
	double c[500],d[500],a1,b1,c1,c2;
	int i,m;

	c[0]=1; d[0]=1; c[1]=1; d[1]=1-x*(a+b)/(a+1);
	c2=c[1]/d[1]; c1=0.0; i=0;
	while(i<500 && fabs(c1-c2)>1.e-12)
	{
		i=i+2; m=i/2; c1=c2;
		c[i]=c[i-1]+m*(b-m)*x*c[i-2]/((a+i-1)*(a+i));
		d[i]=d[i-1]+m*(b-m)*x*d[i-2]/((a+i-1)*(a+i));
		c[i+1]=c[i]-(a+m)*(a+b+m)*x*c[i-1]/((a+i+1)*(a+i));
		d[i+1]=d[i]-(a+m)*(a+b+m)*x*d[i-1]/((a+i+1)*(a+i));
		if(fabs(c[i])>1.e200)
/*	Scale the numerator and denominator to prevent underflow/overflow */
		{
			c[i]=c[i]/1e200; c[i+1]=c[i+1]/1e200;
			d[i]=d[i]/1e200; d[i+1]=d[i+1]/1e200;
		}
		if(fabs(c[i])<1.e-200)
		{
			c[i]=c[i]*1e200; c[i+1]=c[i+1]*1e200;
			d[i]=d[i]*1e200; d[i+1]=d[i+1]*1e200;
		}
		c2=c[i+1]/d[i+1];
	}
	if(c2<0.0) printf(" ** Roundoff error while evaluating the continued fraction for incomplete Beta function at\n a= %e  b=%e  x= %e cont. frac.= %e\n",a,b,x,c2);
	b1=a*log(x)+b*log(1-x)+log(c2)-log(a);
	b1=b1+gammaln(a+b)-gammaln(a)-gammaln(b);
	return exp(b1);
}

/*	To calculate the incomplete Beta function I_x(a,b) using
	the infinite series, called by BETAP

	a,b : (input) Arguments for the complete Beta function
	x : (input) Upper limit of integration defining the incomplete
			Beta function

	Required routines : gammaln
*/


#include <math.h>
#include <stdio.h>

double gammaln(double x);

double betser(double a, double b, double x)
{
	double t,t1,tmax,b1,s3;
	int i;

	s3=1; t=a; tmax=a; i=-1;
	while(i<500 && fabs(t/(a+i+1))>1.e-15)
	{
		i=i+2; t1=t*x*(i-b)/i; t=t1*x*(i+1-b)/(i+1);
		s3=s3+t1/(a+i)+t/(a+i+1);
		if(fabs(t/(a+i+1))>tmax) tmax=fabs(t/(a+i+1));
	}
	if(s3<1.e-16*tmax) printf("  Roundoff error while evaluating the infinite series for incomplete Beta function at\n a= %e  b= %e  x= %e sum= %e  Max term= %e\n",a,b,x,s3,tmax);
	b1=a*log(x)+log(s3)-log(a);
	b1=b1+gammaln(a+b)-gammaln(a)-gammaln(b);
	return exp(b1);
}


/*	To calculate the integrand for incomplete beta function
	The variable AA, BB need to be defined outside

	Required routines : none
*/

double AA, BB;

#include <math.h>

	double fbeta(double x)
{
	return pow(x,AA-1)*pow(1-x,BB-1);
}



/*	To calculate Gamma function for any real value of XG
	Use GAMMALN for calculating the logarithm of Gamma function
	which may be useful for large arguments or when argument is
	close to a negative integer.

	Required functions : None
*/

#include <math.h>

/*	PIS is SQRT(2*PI) */
#define PI 3.14159265358979323846
#define PIS 2.5066282746310005024

double gamma(double xg)

{
	double x,gx,f,y,rmk,f1,x1,fn,fd;
	
/*	The coefficients for rational function approximations */
	double a[2] ={1.767971449569122937e+00,  2.909421117928672645e-01};
	double b[3] ={8.333333333333231537e-02,  1.445531763554246280e-01,
                      2.012779361583001035e-02};
	double a1[6] ={3.905731686764559737e+03,  2.204952264401381785e+03,
                      -1.932467485468849660e+03,  4.643360871045442213e+02,
                      -4.818088806916028754e+01,  1.896853765546068169e+00};
	double b1[7] ={3.918055655523400310e+03, -1.088116266563809683e+02,
                       8.203258626193993149e+02, -9.289402000761705906e+01,
                       6.521113026294866877e+01, -6.090618615608719044e+00,
                       1.475909104740280784e+00};
 
	x=fabs(xg);
	if(x>1000.0) {
/*	Use asymptotic approximation (Stirling formula) */
        gx=(1+1.0/(12*x)+1./(288*x*x)-139.0/(51840*x*x*x)-571.0/(2488320e0*x*x*x*x));
		f=pow(x,x-0.5)*exp(-x)*PIS*gx;
	}

	else if(x>8.0) {
/*	Use rational function approximation for Log(Gamma)  */
		y=1.0/(x*x);
        rmk=((b[2]*y+b[1])*y+b[0])/((a[1]*y+a[0])*y+1);
		f=pow(x,x-0.5)*exp(-x)*PIS*exp(rmk/x);
	}
	else if(x>=2.0) {
/*	Use rational function approximation for (Gamma) over [2,3]
	after translating the range if necessary */
		f1=1.0; x1=x;
		while(x1>3.0) {
			f1=f1*(x1-1.0);
			x1=x1-1.0;
		}

		if(x1==3) f=f1*2.0;
		else if(x1==2) f=f1;

		fn=(((((b1[6]*x1+b1[5])*x1+b1[4])*x1+b1[3])*x1+b1[2])*x1+b1[1])*x1+b1[0];
		fd=(((((a1[5]*x1+a1[4])*x1+a1[3])*x1+a1[2])*x1+a1[1])*x1+a1[0])*x1+1;
		f=f1*fn/fd;
	}
	else if(x>0.0) {
/*	Use rational function approximation for (Gamma) over [2,3]
	after translating the range if necessary */
		f1=1./x;
		x1=x+1.0;
		if(x<1) {f1=f1/x1; x1=x1+1.0;}
		if(x1==2) f=f1;

		fn=(((((b1[6]*x1+b1[5])*x1+b1[4])*x1+b1[3])*x1+b1[2])*x1+b1[1])*x1+b1[0];
		fd=(((((a1[5]*x1+a1[4])*x1+a1[3])*x1+a1[2])*x1+a1[1])*x1+a1[0])*x1+1;
		f=f1*fn/fd;
	}

	if(xg>0.0) return f;
	f=PI/(xg*sin(PI*x)*f);
	return f;
}
 
/*	To calculate the Logarithm of Gamma function for a real argument
 	For negative values it give ln(abs(Gamma(x)))
 
 	Required routines : None
*/

#include <math.h>

#define PI 3.14159265358979323846

double gammaln(double xg)

{
	double x,gx,f,y,rmk,f1,x1,fn,fd;
/*	PIL is LOG(SQRT(2*PI)) */
	double PI2L=0.918938533204672741;
	
/*	The coefficients for rational function approximations */
	double a[4] = {1.013142782275024216e-2,  7.645657825398191944e-1,
             3.381172379819461227e-4,  1.595363637547538209e-2};
	double b[4] ={0.0, 8.333333333338911768e-2, 8.442856404442060242e-4,
                 6.093603832366013704e-2};
	double a1[6] ={4.681163846241230144,  3.208225429683256526,
              5.145525793448859216e-1, 1.581117883959157936e-2,
             -6.398416804905407512e-5, 5.264566254181773919e-7};
	double b1[6] ={2.938038561191284576,  1.489364948862436743,
             -5.466291543917642961,  1.972497734170110410e-1,
              7.830146473241555157e-1, 5.756753067834747499e-2};
 
	x=fabs(xg);
	if(x>10.0) {
/*	Use asymptotic approximation (Stirling formula) */
		y=1.0/x;
		fn=((b[3]*y+b[2])*y+b[1])*y+b[0];
		fd=(((a[3]*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=(x-0.5)*log(x)-x+PI2L+fn/fd;
	}

	else {
/*	Use rational function approximation for (Gamma) over [4,5]
	after translating the range if necessary */
		f1=1.0; x1=x;
		while(x1>5.0) {
			f1=f1*(x1-1.0);
			x1=x1-1.0;
		}
		while(x1<4.0) {
			f1=f1/(x1);
			x1=x1+1.0;
		}


		fn=((((b1[5]*x1+b1[4])*x1+b1[3])*x1+b1[2])*x1+b1[1])*x1+b1[0];
		fd=(((((a1[5]*x1+a1[4])*x1+a1[3])*x1+a1[2])*x1+a1[1])*x1+a1[0])*x1+1;
		f=log(f1)+fn/fd;
	}

	if(xg>0.0) return f;
	f=log(PI/(x*fabs(sin(PI*x))))-f;
	return f;
}
 
/*	To calculate the incomplete Gamma function P(a,x) for
	positive real arguments
	It returns a value of -1 for a.le.0 or x<0

	a : (input) Argument for the complete Gamma function
	x : (input) Upper limit of integration defining the incomplete
			Gamma function

	Required routines : gamma, gammaln
*/
 

#include <math.h>

double gamma(double xg);
double gammaln(double xg);

double gammap(double a,double x)
{
	double c[500],b[500],f,s1,s2,g1,t,c1,c2;
	int i;

	if(a<=0 || x<0) return -1.0; 

	if(x<3)
	{
/*	Use power series */
        	s1=1; t=a; i=0;
		while(i<500 && fabs(t/(a+i))>1.e-14)
		{
			i=i+1; t=-t*x/i;
			s1=s1+t/(a+i);
		}
		if(a<140) f=pow(x,a)*s1/gamma(a+1);
		else f=0.0;
	}
	else if(a<1.2*x)
	{
/*	Use continued fraction */
        	c[0]=1; b[0]=x; c[1]=1; b[1]=x+1-a;
		c2=c[1]/b[1]; c1=0.0; i=0;
		while(i<500 && fabs(c1-c2)>1.e-12)
		{
			i=i+2; c1=c2;
			c[i]=x*c[i-1]+(i/2)*c[i-2];
			b[i]=x*b[i-1]+(i/2)*b[i-2];
			c[i+1]=c[i]+(i/2+1-a)*c[i-1];
			b[i+1]=b[i]+(i/2+1-a)*b[i-1];
			if(fabs(b[i+1])>1.e200)
			{
				c[i]=c[i]/1.e200; c[i+1]=c[i+1]/1.e200;
				b[i]=b[i]/1.e200; b[i+1]=b[i+1]/1.e200;
			}
			if(fabs(b[i+1])<1.e-200)
			{
				c[i]=c[i]*1.e200; c[i+1]=c[i+1]*1.e200;
				b[i]=b[i]*1.e200; b[i+1]=b[i+1]*1.e200;
			}
			c2=c[i+1]/b[i+1];
		}
		g1=-x+a*log(x)+log(c2)-gammaln(a); f=1-exp(g1);
	}
	else
	{
/*	Use the power series for a>x */
		s2=1; t=1; i=0;
		while(i<500 && fabs(t)>1.e-14)
		{
			i=i+1; t=t*x/(a+i); s2=s2+t;
		}
		g1=-x+a*log(x)+log(s2)-gammaln(a+1); f=exp(g1);
	}
	return f;
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
