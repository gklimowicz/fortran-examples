/*	To simulate a binomial distribution */

#include <stdio.h>
#include <math.h>

double gamma(double xg);
double gammaln(double xg);
double gammap(double a,double x);
int iranbin(double *seed, int n, double p, double c[]);

main()
{
	int i,n,ii,nu,num,ip[200];
	double p,p0,r,rn,seed,chi,prob,c[200];

	rn=1e4; num=rn;
	printf(" Type n = No. of trials (<200)\n");
	printf(" and p = probability of event in binomial distribution\n");
	scanf("%d %le",&n,&p);
	if(n>=200) return;
	seed=-1;
	for(i=0; i<200; ++i) ip[i]=0.0;
/*	Generate num random numbers with binomial distribution and
        place them in the appropriate bin */
	for(i=0; i<num; ++i){
		ii=iranbin(&seed,n,p,c);
		ip[ii]=ip[ii]+1;
	}

/*	Apply the Chi-square test to check the distribution
        p0 is the expected probability
        This will not work if initial P0 gives underflow	*/
	p0=pow(1.0-p,n); chi=0.0;
	printf("       I        Prob(I)\n");
	for(i=0; i<=n; ++i) {
		printf(" %d %e %e \n",i,ip[i]/rn,p0);
		r=p0*rn-ip[i]; chi=chi+r*r/(p0*rn);
		p0=p0*(n-i)*p/((i+1)*(1.0-p));
	}
	nu=n;
	printf(" Chi-square = %e    Degree of freedom = %d \n",chi,nu);
/*	If NUM is increased and approximation to Gaussian distribution
        is used in IRANBIN the test may fail if P is close to 1
        as the Gaussian distribution may not be a good approximation in
        that case. RMAX in IRANBIN may be increased to avoid use of
        Gaussian approximation */
	prob=1-gammap(nu/2.0,chi/2);
	printf(" Probability of finding a value of CHI or larger = %e\n",prob);
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


/*	To generate random numbers with Binomial distribution 
        For large values of mean it approximates the distribution by a
        normal distribution with the same mean and standard deviation.
 
        SEED : (input/output) real seed, it should be negative during
               the first call and should not be modified between two
               calls with the same N, P. If a new value of N or P is used
               it should be reset to a negative value.
        N   :  (input) Number of trials in the binomial distribution
        P   :  (input) probability of the event
        C   :  (output) Real array of length N used to store the
               cumulative probability table for use in subsequent
               calculations. This array should not be modified by the
               user. This is used only if the mean N*P<RMAX and there
               is no underflow while calculating (1-P)**N
 
        Required routines : None

	If P is close to 1 the approximation by Gaussian may not be very
        good and in that case increasing RMAX may help although it will
        take longer time to generate the numbers.
*/


#include <math.h>

int iranbin(double *seed, int n, double p, double c[])
{
	int i;
	double am=2147483648.0,a=45875.0,ac=453816693.0,an=2147483647.0;
	double pi=3.14159265358979324,rmax=100;
	double rmu,sig,r1,rangau,p0,r;

	rmu=n*p; sig=sqrt(rmu*(1-p)); p0=pow(1-p,(double) n);

/*	For large mu or if P0 gives underflow use normal distribution */
	if(rmu > rmax || p0==0) {
/*	Initialise the seed during first call */
		if(*seed<0.0) *seed=fmod(fabs(*seed),am);
		r1=fmod(a*(*seed)+ac,am);
		if(*seed==0.0) *seed=0.1;
		rangau=sqrt(2.0*log(an/(*seed)))*cos(2.0*pi*r1/an);
		*seed=r1;
		i=rmu+sig*rangau+0.5;
		if(i<0) i=0; if(i>n) i=n;
		return i;
	}
	else {
/*       Initialise the array during the first call */
		if(*seed<0.0) {
			c[0]=p0;
			for(i=1; i<=n; ++i) {
				p0=p0*p*(n-i+1)/((1-p)*i);
				c[i]=c[i-1]+p0;
			}
			*seed=fabs(*seed);
		}

/*	generate the random number */
		*seed=fmod((*seed)*a+ac,am); r=(*seed)/an;
		for(i=0; i<=n; ++i) {
			if(r<c[i]) return i;
		}
		return n;
	}
}
	





	

	

