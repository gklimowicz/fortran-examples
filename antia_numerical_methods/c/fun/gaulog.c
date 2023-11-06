/*     To integrate a function with logarithmic singularity Gaussian formulas

	RINT : (output) Calculated value of the integral
	A : (input) The upper limit
	AEPS : (input) The required absolute accuracy
	REPS : (input) The required relative accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the
		integrand (divided by LOG(A/X))
	NPT : (output) Number of function evaluations used
		
	Error status is returned by the value of the function GAULOG.
		0 value implies successful execution
     		30 implies specified accuracy was not achieved
     		DIF will contain the estimated accuracy

	Function F(X) must be supplied by the user
	Note that function calculates integral of F(X)*LOG(A/X)

	Required functions : F
*/

#include <math.h>

int gaulog(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt)

{
	int i,j,n;
	double r1,s1,r2;

/*     Weights and abscissas for Gaussian formula for logarithmic singularity
     W[N-2],...,W[2N-3], are the weights for N-point rule and
     X[N-2],...,X[2N-3], the corresponding abscissas
     Weights and abscissas are available for N=2,4,8,16
*/
 
	double x[30] = {1.120088061669761830e-1, 6.022769081187381028e-1,
			4.144848019938322080e-2, 2.452749143206022519e-1,
                        5.561654535602758372e-1, 8.489823945329851746e-1,
                        1.332024416089246501e-2, 7.975042901389493841e-2,
                        1.978710293261880538e-1, 3.541539943519094197e-1,
                        5.294585752349172777e-1, 7.018145299390999638e-1,
                        8.493793204411066760e-1, 9.533264500563597888e-1,
                        3.897834487115909095e-3, 2.302894561687320045e-2,
                        5.828039830624031972e-2, 1.086783650910538817e-1,
                        1.726094549098437244e-1, 2.479370544705782363e-1,
                        3.320945491299168705e-1, 4.221839105819483085e-1,
                        5.150824733814623250e-1, 6.075561204477284747e-1,
                        6.963756532282138523e-1, 7.784325658732652431e-1,
                        8.508502697153909688e-1, 9.110868572222718348e-1,
                        9.570255717035421226e-1, 9.870478002479844660e-1};

	double w[30] = {7.185393190303844407e-1,2.814606809696155593e-1,
                        3.834640681451351249e-1, 3.868753177747626273e-1,
                        1.904351269501424154e-1, 3.922548712995983245e-2,
                        1.644166047280028868e-1, 2.375256100233060205e-1,
                        2.268419844319191264e-1, 1.757540790060702450e-1,
                        1.129240302467590519e-1, 5.787221071778207240e-2,
                        2.097907374213297804e-2, 3.686407104027619013e-3,
                        6.079171004359114509e-2, 1.029156775175820228e-1,
                        1.223556620460090919e-1, 1.275692469370159323e-1,
                        1.230135746000709083e-1, 1.118472448554857552e-1,
                        9.659638515212439849e-2, 7.935666435147320573e-2,
                        6.185049458196527197e-2, 4.543524650772672381e-2,
                        3.109897475158184829e-2, 1.945976592736087029e-2,
                        1.077625496320554213e-2, 4.972542890087649610e-3,
                        1.678201110051197249e-3, 2.823537646684367889e-4};
 
/*     The 2-point formula */
	r1=(f(a*x[0])*w[0]+f(a*x[1])*w[1])*a;
	*npt=2;
	n=2;
 
/*     Use higher order formula until convergence */
	for(j=2; j<=4; ++j) {
		n=n*2;
		r2=0.0;
		for(i=n-2; i<2*n-2; ++i) r2=r2+f(x[i]*a)*w[i];
		r2=r2*a;
 
		*npt=(*npt)+n;
		*dif=r2-r1;
		*rint=r2;
		s1=reps*fabs(*rint); if(aeps>s1) s1=aeps;
		if(fabs(*dif) < s1) return 0;
		r1=r2;
	}
 
	return 30;
}
