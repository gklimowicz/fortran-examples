/*	Real root of a system of nonlinear equations using Newton's method

	FCN : (input) Name of the function to calculate the vector function
	NP : (input) Number of variables and equations to be solved
	X : (input/output) Array of length NP containing the starting
		values. After execution it should contain the computed roots
	F : (output) Array of length NP containing the function values
		at the point specified by array X
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X[I]))
		
	Error status is returned by the value of the function NEWTON.
		0 value implies successful execution
		441 implies that GAUELM failed to solve the required
			system of linear equations
		442 implies that iteration did not converge to required
			accuracy

	Function FCN(NP,X,F,DF) must be supplied by the user.
	Here X and F are arrays of length NP containing the value of
	vector function at point X. DF is an array of length NP*NP 
	containing the Jacobian. The second dimension of DF in FCN should be NP.
	The system of equations to be solved are F[I]=0

	Required functions : GAUELM, FCN
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int newton(void fcn(int , double * ,double *, double * ), int np,
	double x[], double f[], double reps, double aeps)

{
	int i,j,lj,ipas,num,iflg,ier,qchk;
	double det,r1;
	double *wk, *wk1;
	int *intc;

	lj=np; num=1;
	wk=(double *) calloc((size_t) (np*np), sizeof(double));
	wk1=(double *) calloc((size_t) np, sizeof(double));
	intc=(int *) calloc((size_t) np, sizeof(int));

/*	Loop for Newton's iteration */
	for(ipas=1; ipas<=200; ++ipas) {
		fcn(np,x,f,wk);
/*	The right hand side of system of equations for Newton's iteration */
		for(i=0; i<np; ++i) wk1[i]=f[i];

		iflg=0;
		ier=gauelm(np,num,wk,wk1,&det,intc,lj,&iflg);
		if(ier>0) {free(intc); free(wk1); free(wk); return 441;}

/*	Convergence check */
		qchk=1;
		for(j=0; j<np; ++j) {
			r1=reps*fabs(x[j]); if(aeps>r1) r1=aeps;
			if(fabs(wk1[j]) > r1) qchk=0;
			x[j]=x[j]-wk1[j];
		}
		if(qchk==1) {free(intc); free(wk1); free(wk); return 0;}
	}

	free(intc); free(wk1); free(wk);
	return 442;
}
