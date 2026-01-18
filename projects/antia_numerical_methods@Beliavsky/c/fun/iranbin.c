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
	





	

	

