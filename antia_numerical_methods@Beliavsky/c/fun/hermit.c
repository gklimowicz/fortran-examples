/*     To integrate a function over infinite interval using Gauss-Hermite formulas

	RINT : (output) Calculated value of the integral
	AEPS : (input) The required absolute accuracy
	REPS : (input) The required relative accuracy
     		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the
     		integrand (multiplied by EXP(X**2))
	NPT : (output) Number of function evaluations used
		
	Error status is returned by the value of the function HERMIT.
		0 value implies successful execution
	     	30 implies specified accuracy was not achieved
     		DIF will contain the estimated accuracy

	Function F(X) must be supplied by the user

	Required functions : F
*/

#include <math.h>

int hermit(double *rint, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt)

{
	int i,j,n;
	double r1,r2,s1;
 
/*     Weights and abscissas for Gauss-Hermite quadrature
     Because of symmetry only positive abscissas are listed
     W[N/2-1],...,W[N-2], are the weights for N-point rule and
     X[N/2-1],...,X[N-2], the corresponding abscissas
     Weights and abscissas are available for N=2,4,8,16,32
*/
 
	double x[31] = {    0.707106781186548e0,  0.524647623275290e0,
     1.650680123885785e0,  0.3811869902073221e0, 1.1571937124467802e0,
     1.98165675669584293e0,2.93063742025724402e0,0.27348104613815245e0,
     0.82295144914465589e0,1.3802585391988808e0, 1.9517879909162540e0,
     2.54620215784748136e0,3.17699916197995603e0,3.8694479048601227e0,
     4.68873893930581836e0,0.19484074156939933e0,0.58497876543593245e0,
     0.97650046358968284e0,1.37037641095287184e0,1.76765410946320160e0,
     2.16949918360611217e0,2.57724953773231745e0,2.99249082500237421e0,
     3.41716749281857074e0,3.85375548547144464e0,4.30554795335119845e0,
     4.77716450350259639e0,5.27555098651588013e0,5.81222594951591383e0,
     6.409498149269660412e0,7.1258139098307275728e0};

	double w[31] = {    0.886226925452758e0,  0.804914090005513e0,
     8.13128354472452e-2,  0.6611470125582413e0, 0.207802325814892e0,
     1.707798300741347e-2, 1.9960407221136762e-4,0.50792947901661374e0,
     0.2806474585285342e0, 8.3810041398985349e-2,1.2880311535509972e-2,
     9.3228400862418052e-4,2.7118600925378815e-5,2.3209808448652106e-7,
     2.654807474011182e-10,0.37523835259280239e0,0.27745814230252993e0,
     0.15126973407664245e0,6.0458130955912614e-2,1.7553428831573430e-2,
     3.6548903266544280e-3,5.3626836552797204e-4,5.4165840618199826e-5,
     3.6505851295623761e-6,1.5741677925455940e-7,4.0988321647708966e-9,
     5.933291463396639e-11,4.215010211326448e-13,1.197344017092849e-15,
     9.2317365365182922e-19,7.3106764273841624e-23};
 
 
/*     The 2-point formula */
	r1=w[0]*(f(x[0])+f(-x[0]));
	*npt=2;
	n=2;
 
/*     Use higher order formula until convergence */
	for(j=2; j<=5; ++j) {
		n=n*2;
		r2=0.0;
		for(i=n/2-1; i<n-1; ++i) r2=r2+(f(x[i])+f(-x[i]))*w[i];

		*npt=(*npt)+n;
		*dif=r2-r1;
		*rint=r2;
		s1=reps*fabs(*rint); if(aeps>s1) s1=aeps;
		if(fabs(*dif)<s1) return 0;
		r1=r2;
	}

	return 30;
}
