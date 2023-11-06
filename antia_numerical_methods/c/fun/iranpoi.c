/*	To generate random numbers with Binomial distribution 
        For large values of mean it approximates the distribution by a
        normal distribution with the same mean and standard deviation.
 
        SEED : (input/output) real seed, it should be negative during
               the first call and should not be modified between two
               calls with the same N, P. If a new value of N or P is used
               it should be reset to a negative value.
        RMU :  (input) Mean value of the Poisson distribution
        C   :  (output) Real array of length NMAX used to store the
               cumulative probability table for use in subsequent
               calculations. This array should not be modified by the
               user. This is used only if RMU is smaller than some
               critical value.
 
        Required routines : None

	Increasing NMAX may give underflow while calculating exp(-rmu)

*/


#include <math.h>

int iranpoi(double *seed, double rmu, double p[])
{
	int i,n;
	double am=2147483648.0,a=45875.0,ac=453816693.0,an=2147483647.0;
	double pi=3.14159265358979324,nmax=800;
	double sig,r1,rangau,p0,r;

	sig=sqrt(rmu); n=rmu+7*sig; if(n<20) n=20;

/*	For large rmu use normal distribution */
	if(n > nmax) {
/*	Initialise the seed during first call */
		if(*seed<0.0) *seed=fmod(fabs(*seed),am);
		r1=fmod(a*(*seed)+ac,am);
		if(*seed==0.0) *seed=0.1;
		rangau=sqrt(2.0*log(an/(*seed)))*cos(2.0*pi*r1/an);
		*seed=r1;
		i=rmu+sig*rangau+0.5;
		if(i<0) i=0;
		return i;
	}
	else {
/*       Initialise the array during the first call */
		if(*seed<0.0) {
			p0=exp(-rmu);
			p[0]=p0;
			for(i=1; i<=n; ++i) {
				p0=p0*rmu/i;
				p[i]=p[i-1]+p0;
			}
			*seed=fabs(*seed);
		}

/*	generate the random number */
		*seed=fmod((*seed)*a+ac,am); r=(*seed)/an;
		for(i=0; i<n; ++i) {
			if(r<p[i]) return i;
		}
		return n-1;
	}
}
