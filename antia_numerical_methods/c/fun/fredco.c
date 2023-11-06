/*	To solve linear Fredholm equation of first or second kind using
	collocation method
	It can be used to solve Volterra equations by defining the kernel
	to be zero for t>x

	N : (input) Number of points to be used in the collocation method
	A : (input) Lower limit of the integral
	B : (input) Upper limit of the integral
	F : (output) Array of length N containing the calculated coefficients
		of expansion
		solution = SUM_I F[I]*PHI(I,X)
	X : (input) Array of length N containing the points to
		be used for collocation.
	REPS : (input) Required relative accuracy to which integrals are to
		be calculated
	AEPS : (input) Required absolute accuracy to which integrals are to
		be calculated
		REPS and AEPS are passed on to function ADPINT for calculating
		the integrals when IQ=0. Otherwise these variables are not used.
	IQ : (input) Integer variable to specify the treatment for integrals
		PSI(I,X)=Integral[FKER(X,T)*PHI(I,T) dT] over [A,B]
		If IQ=0, then the integrals are evaluated using function
			ADPINT, using function FUNK to calculate the
			integrand, which in turn requires, PHI(I,T) and FKER(X,T). 
		Otherwise the integrals are calculated using a user supplied
			function PSI(I,X).
	IT : (input) Integer variable to specify the type of integral equation
		If IT=1 Fredholm equation of the first kind is solved
		If IT=2 Fredholm equation of the second kind is solved
		
	Error status is returned by the value of the function FREDCO.
		0 value implies successful execution
		708 implies N<1, IT>2 or IT<=0, No calculations are done.
		Other values may be set by GAUELM and ADPINT

	Functions FG(X), PHI(I,X) and FKER(X,T) (for IQ=0)
	or function PSI(I,X) (for IQ!=0) must be supplied by the user.
	Names of these functions are fixed. FG(X) is the right
	hand side function g(x), FKER(X,T) is the kernel, PHI(I,X) calculates
	the basis functions phi_i(x), while PSI(I,X) calculates the integrals
	defined above. Global variables XFRED and IFRED are used to pass
	on the values of X and I to FUNK.

	Required functions : GAUELM, ADPINT, KRONRD, FUNK, FG, FKER, PHI, PSI
*/

#include <math.h>
#include <stdlib.h>

double funk(double x);
double fg(double x);
double fker(double x, double t);
double phi(int i, double x);
double psi(int i, double x);
int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax);

/* To pass arguments to FUNK */
double XFRED;
int IFRED;


int fredco(int n, double a, double b, double f[], double x[], double reps,
	double aeps, int iq, int it)

{
	int i,j,nmax,npt,npt1,ier,iflg;
	double det,xi,ri,ri1,dif,dif1;
	int *iwk;
	double *wk;

	if(n<1 || it>2 || it<=0) return 708;
	nmax=10000;
	iwk=(int *) calloc((size_t) n, sizeof(int));
	wk=(double *) calloc((size_t) (n*n), sizeof(double));

/*	Setting up the system of linear equations */
	for(j=0; j<n; ++j) {
		f[j]=fg(x[j]);
		for(i=0; i<n; ++i) {
			if(iq==0) {
/*	Evaluate the integrals numerically, split the range into two 
	to tackle possible discontinuity at t=X */
				XFRED=x[i];  IFRED=j;
				ier=adpint(&ri,a,x[i],reps,aeps,&dif,funk,&npt,&nmax);
				if(ier>100) {free(wk); free(iwk); return ier;}
				ier=adpint(&ri1,x[i],b,reps,aeps,&dif1,funk,&npt1,&nmax);
				if(ier>100) {free(wk); free(iwk); return ier;}
				wk[j+i*n]=ri+ri1;
			}

			else {
/*	Calculate the integrals PSI(I,X) using user supplied function */
				wk[j+i*n]=psi(j,x[i]);
			}
			if(it==2) wk[j+i*n]=wk[j+i*n]-phi(j,x[i]);
		}
	}

/*	Solve the resulting system of linear equations */
	iflg=0;
	ier=gauelm(n,1,wk,f,&det,iwk,n,&iflg);
	free(wk); free(iwk);
	return 0;
}
