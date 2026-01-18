/*	Calculate the correlation coefficient and the probability of
        finding a value as high as that for uncorrelated dats sets */


#include <stdio.h>
#include <math.h>

double gammaln(double xg);
double erf(double x);
double pcor(int n, double xx);
double ran1(double *seed);


main()
{
	int i,n;
	double seed,r,xmean,ymean,sxx,syy,sxy,rcor,prob,rx,ry;
	double x[100],y[100];

	seed=9; n=100;
/*	Generate the data */
	for(i=0; i<n; ++i) {
		r=ran1(&seed); x[i]=i*0.01+0.001*(r-0.5);
		r=ran1(&seed); y[i]=i*1e-5+0.01*(r-0.5);
	}

/*	Calculate the correlation coefficient between X and Y */
	xmean=0.0; ymean=0.0;
	for(i=0; i<n; ++i){
		xmean=xmean+x[i]; ymean=ymean+y[i];}
	xmean=xmean/n; ymean=ymean/n;
	sxx=0.0; syy=0.0; sxy=0.0;
	for(i=0; i<n; ++i){
		rx=x[i]-xmean; ry=y[i]-ymean;
		sxx=sxx+rx*rx; syy=syy+ry*ry; sxy=sxy+rx*ry;
	}
	rcor=sxy/sqrt(sxx*syy);
	prob=pcor(n-2,fabs(rcor));
	printf("Correlation Coefficient = %f\n",rcor);
	printf("Probability of finding this value or higher for uncorrelated sequence =,%le\n",prob);
}

/*	To calculate the Error function for a real argument 

	Required functions : None
*/

#include <math.h>

double erf(double x)

{
	double fn,fd,f,y;

/*	The coefficients of rational function approximations */
	double a[7] = {4.891837297874833514e-01,  1.110175596090841322e-01,
                       1.526977188817169289e-02,  1.388143322498740953e-03,
                       8.446154421878829637e-05,  3.239915935631695227e-06,
                       6.200069065781009292e-08};
	double b[8] = {1.128379167095512575e+00,  1.758583405424390318e-01,
                       5.411290829613803886e-02,  3.805760982281134630e-03,
                       4.314620532020313106e-04,  1.423201530735308681e-05,
                       6.188698151904667110e-07,  2.293915472550064153e-09};
	double a1[7] = {3.040844316651787557e+01,  3.358927547920980388e+02,
                        1.703170048554673596e+03,  4.133195606956137156e+03,
                        4.555611776312694034e+03,  1.935778559646575488e+03,
                        2.051089518116697193e+02};
	double b1[8] = {5.641895835477562749e-01,  1.687403209467950089e+01,
                        1.813522721872712655e+02,  8.779664433851135986e+02,
                        1.965115658619443782e+03,  1.865558781108286245e+03,
                        5.828865035607128761e+02,  2.558559157228883880e+01};
 
	if(fabs(x)<2.0) {
/*	Use approximation for Erf(x)/x */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=x*fn/fd;
	}

	else {
/*	Use approximation for x*Exp(x^2)*Erfc(x) */
		y=1.0/(x*x);
		fn=((((((b1[7]*y+b1[6])*y+b1[5])*y+b1[4])*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((((a1[6]*y+a1[5])*y+a1[4])*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		f=1.-exp(-x*x)*fn/(fd*fabs(x));
		if(x<0.0) f=-f;
	}
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
 
/*	To calculate the probability that two uncorrelated sequences
	of length (n+2) will give a correlation coefficient exceeding |x|
	For large even values of n the series can give large roundoff
	in which case the distribution is approximated by Normal distribution

	N  : (input) Length of data sets should be N+2
        XX : (input) The function will calculate the probability for
                correlation exceeding |x| for uncorrelated sequences

	Required routines : GAMMALN, ERF
*/

#include <math.h>

double gammaln(double xg);
double erf(double x);

#	define PI 3.14159265358979323846

double pcor(int n, double xx)

{
	double f,p,ss,pl,pmax,x;
/*	pis is sqrt(pi) */
	double pis=1.7724538509055144;
	int l,k;

	x=fabs(xx);
	if(x>1.0) return 0.0;
	if(2*(n/2)==n)
	{
/*	if n is even */
		l=(n-2)/2; ss=x; p=x; 
		pmax=p;
		for(k=1; k<=l; k++)
		{
			p=-(l-k+1)*p*x*x/k; ss=ss+p/(2*k+1.0);
			if(fabs(p)>pmax) pmax=fabs(p);
		}
		pl=gammaln((n+1.0)/2)-gammaln(n/2.0);
		f= 2.*exp(pl)*ss/pis;
		if(pmax>1.e5||f>1.0) f=erf(x*sqrt(n/2.0));
	}
	else
	{
/*	if n is odd */
		l=(n-3)/2; ss=sqrt(1-x*x); p=ss; 
		if(n==1) ss=0.0;
		for(k=1; k<=l; k++)
		{
			p=p*2*k*(1-x*x)/(2*k+1.0); ss=ss+p;
		}
		f=(asin(x)+x*ss)*2/PI;
	}
	return 1-f;
}


/*	To generate uniformly distributed random numbers in interval (0,1)

	SEED : (input/output) is a real value used as the seed
		It should be positive during initial call and
		should not be modified between different calls

	Required functions : None
*/

#include <math.h>

double ran1(double *seed)

{
	double am=2147483648e0, a=45875e0, ac=453816693e0, an=2147483647e0;

	*seed=fmod((*seed)*a+ac,am);
	return (*seed)/an;
}
	
