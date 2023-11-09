/*	To calculate Gamma function and its natural logarithm */

#include <stdio.h>
#include <math.h>

double gammaln(double xg);
double gamma(double xg);

main()
{
	int i,i1,j,k,n;
	double x,f1,f2;

	for(i1=0; i1<99; ++i1) {
		printf("type x,        (quits when x=-9999)\n");
		scanf(" %le", &x);
		if(x==-9999) return 0;

		f1=gamma(x);
		f2=gammaln(x);
		printf(" x= %e    Gamma(x)= %e    log(|Gamma(x)|) =%e \n",x,f1,f2);
	}
	return;
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
	int ix;
	
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
	ix=x;
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
	if(x>ix) f=PI/(xg*sin(PI*x)*f);
	else f=pow(-1.0,ix)/0.0;
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
	int ix;
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
	ix=x;
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
	if(x>ix) f=log(PI/(x*fabs(sin(PI*x))))-f;
	else f=pow(-1.0,ix)/0.0;
	return f;
}
 
