/*	Multiple integration over a hyper-rectangle in n-dimensions
	using equidistributed sequences

	A : (input) Array of length N containing the lower limit
		along each dimension
	B : (input) Array of length N containing the upper limit
		along each dimension
	N : (input) The number of dimensions
	NPT : (input) Maximum number of function evaluations to be used
	F : (input) Name of the function to calculate the integrand
		Function F(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	S1 : (output) The calculated value of the integral
	S2 : (output) Another approximation to the value of the integral
		For smooth functions S2 is expected to be better approximation
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(S2))
	DIF : (output) estimated (absolute) error achieved by the function
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function EQUIDS.
		0 value implies successful execution
		39 implies specified accuracy was not achieved in
			which case DIF will contain the estimated accuracy
		312 implies N<1 or N>21 and no calculations are done

	Function F(N,X) must be supplied by the user
	
	Required functions :  F
*/

#include <math.h>


int equids(double a[], double b[], int n, int npt, double (*f) (int , double * ),
	double *s1, double *s2, double reps, double aeps, double *dif, int *np)

{
	int i,j,npt1,ia, nchk=100, nmax=21;
	double ri,hh,ri1,dif1,diff,ss1,a1,r1,h[21],xa[21],wt[21];
/*	The first 21 prime numbers */
	double at[21] ={2.,3.,5.,7.,11.,13.,17.,19.,23.,29.,31.,37.,41.,
             43.,47.,53.,59.,61.,67.,71.,73.};

	*s1=f(n,a);
	*s2=(*s1);
	*np=1;
	if(n>nmax || n<1) return 312;

	hh=1.0;
	for(i=0; i<n; ++i) {
		h[i]=b[i]-a[i];
		hh=hh*h[i];
/*	The irrational numbers for generating equidistributed sequences */
		wt[i]=sqrt(at[i]);
	}

	ri=0.0; ri1=0.0; *dif=0.0;
	npt1=nchk;
	for(i=1; i<=npt; ++i) {
/*	Generate the abscissas using equidistributed sequences */
		for(j=0; j<n; ++j) {
			a1=i*wt[j];
			ia=a1+0.5; a1=2.*fabs(a1-ia)*h[j];
			xa[j]=a[j]+a1;
		}
/*	Accumulate the sum */
		*s1=(*s1)+2.*(f(n,xa)-ri);
		*s2=(*s2)+(*s1);

		if(i-nchk*(i/nchk) == 0) {
/*	To control the roundoff error form partial sums */
			ss1=ri+(*s1)/(2*i+1);
			diff=(*s2)/((i+1.)*(i+1.));
			*s2=0.0;
			*s1=(*s1)-(2*i+1)*diff;
/*	The new approximation to the average value of function */
			ri=ri+diff;

			if(i==npt1) {
/*	Check for convergence */
				dif1=*dif;
				*dif=fabs(ri-ri1);
				ri1=ri;
				r1=reps*fabs(ri); if(aeps/hh > r1) r1=aeps/hh;
				if((*dif+dif1 < r1) && i > 5*nchk) {
					*s1=ss1*hh;
					*s2=ri*hh;
					*dif=(*dif+dif1)*hh;
					*np=i+1;
					return 0;
				}

				npt1=npt1*2;
			}
		}
	}

/*	Integral fails to converge */
	*s1=(ri+(*s1)/(2*npt+1))*hh;
	*s2=(ri+(*s2)/((npt+1.0)*(npt+1.0)))*hh;
	*dif=(*dif+dif1)*hh;
	*np=npt+1;
	return 39;
}
