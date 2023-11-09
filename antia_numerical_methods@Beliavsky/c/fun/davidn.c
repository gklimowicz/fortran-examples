/*	Real root of a system of nonlinear equations using Davidenko's method

	FCN : (input) Name of the function to calculate the vector function
		A sample function can be found at the end of this file.
	NP : (input) Number of variables and equations to be solved
	X : (input/output) Array of length NP containing the starting
		values. After execution it should contain the computed roots
	F : (output) Array of length NP containing the function values
		at the point specified by array X
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X(I)))
	THETA : which is read by the function is the parameter introduced
		into the equations. Normally THETA=1 to start with and
		initial guess should satisfy the equations with this parameter
		THETA=0 corresponds to the equations that need to be solved.
		The value of THETA should be gradually reduced from 1 to 0
		
	Error status is returned by the value of the function DAVIDN.
		0 value implies successful execution
		407 implies NP > NMAXD or NP < 1, in which case no
			calculations are done
		440 implies that THETA is nonzero even after 100 steps
		Other values may be set by functions NEWTON or BROYDN

	Function FCN(NP,X,F,DF) or FCN(NP,X,F)
	must be supplied by the user. Here X and F are arrays of length
	NP with F containing the value of vector function at point X. DF is
	an array of length NP*NP containing the Jacobian which is required
	by function NEWTON. If function BROYDN is used then DF is not
	required. The function FCN should parameterise the equations in
	terms of THETA as explained in section 7.16. A sample function is
	appended at the end of this file.
	DAVIDN_B is the version to be used with BROYDN

	Required functions : NEWTON (or BROYDN), GAUELM, FCN
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NMAXD 200

/*	To pass on parameters to FCN this must be included before fcn */
double THETA, X0[NMAXD], F0[NMAXD];
int newton(void fcn(int , double * ,double *, double * ), int np,
	double x[], double f[], double reps, double aeps);

int davidn(void (*fcn) (int ,double * , double * ,double *), int np, double x[],
	double f[], double reps, double aeps)

{
	int i,j,ier;
	double reps0,aeps0;
	double *wk;


	if(np>NMAXD || np<1) return 407;

	printf(" Type in the initial guess \n");
	for(i=0; i<np; ++i) scanf(" %le",&x[i]); 

/*	To initialise the parameterised function */
      THETA=0.0;
	  wk=(double *) calloc((size_t) (np*np),sizeof(double));
	  fcn(np,x,F0,wk);
	  free(wk);

	  for(i=0; i<np; ++i) X0[i]=x[i];
	  reps0=reps; if(0.001>reps0) reps0=0.001;
	  aeps0=aeps; if(1.0e-4>aeps0) aeps0=1.0e-4;

	  for(i=1; i<=100; ++i) {
		  printf(" Type THETA \n");
		  scanf(" %le",&THETA);
		  printf("  THETA = %e \n",THETA);
		  if(THETA == 0.0) {reps0=reps; aeps0=aeps;}

		  ier=newton(fcn,np,x,f,reps0,aeps0);
		  if(ier>0)  return ier;
		  printf(" x =");
		  for(j=0; j<np; ++j) printf(" %e",x[j]);
		  printf("\n");
		  if(THETA == 0) return 0;
	  }

	  return 440;
}
/*
	double THETA, X0[NMAXD], F0[NMAXD];
	void fcn(int np, double x[], double f[], double *df)
	{

	fun(np,x,f,df);

	for(i=0; i<np; ++i) f[i]=f[i]-THETA*F0[i];
	}

*/

