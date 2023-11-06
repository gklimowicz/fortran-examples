/*	Real root of a system of nonlinear equations using Broydon's method

	FCN : (input) Name of the function to calculate the vector function
	NP : (input) Number of variables and equations to be solved
	X : (input/output) Array of length np containing the starting
		values. After execution it should contain the computed roots
	F : (output) Array of length np containing the function values
		at the point specified by array X
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))

	Error status is returned by the value of the function BROYDN.
		0 value implies successful execution
		441 implies that GAUELM failed to solve the required
			system of linear equations
		442 implies that iteration did not converge to required
			accuracy

	Function FCN(NP,X,F) must be supplied by the user.
	Here X and F are arrays of length NP containing the value of
	vector function at point X.
	The system of equations to be solved are F[I]=0

	Required functions : GAUELM, FCN
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int broydn(void fcn(int , double * , double * ), int np, double x[],
	double f[], double reps, double aeps)

{
	int i,j,k,lj,ip,ipas,num,iflg,ier,qchk;
	double det,r1,x1,dk,s1,s2,ss;
	double *wk, *wk1;
	int *intc;

	lj=np;
	fcn(np,x,f);
	wk=(double *) calloc((size_t) (np*np),sizeof(double));
	ip=np*np; if(4>np) ip=4*np;
	wk1=(double *) calloc((size_t) (ip),sizeof(double));

/*	Generating initial approximation to the Jacobian */
	for(ip=0; ip<np; ++ip) {
		x1=x[ip];
		dk=0.01*fabs(x1);
		if(dk==0.0) {dk=aeps*100.0; if(1.0e-6 > dk) dk=1.0e-6;}
		x[ip]=x1+dk;
		fcn(np,x,wk);
		for(j=0; j<np; ++j) wk1[j*np+ip]=(wk[j]-f[j])/dk;
		x[ip]=x1;
	}

	num=np; iflg=0;
	intc=(int *) calloc((size_t) (np),sizeof(int));
	for(i=0; i<np; ++i) {
		for(j=0; j<np; ++j) wk[j*np+i]=0.0;
		wk[i*np+i]=1.0;
	}

/*	To calculate H = B**(-1) */
	ier=gauelm(np,num,wk1,wk,&det,intc,lj,&iflg);
	if(ier>0) {free(intc); free(wk1); free(wk); return 441;}

	for(i=0; i<np; ++i) {
		wk1[i]=f[i];
		wk1[i+np]=x[i];
	}

/*	Broyden's iteration */
	for(ipas=1; ipas<=200; ++ipas) {

/*	Convergence check */
		qchk=1;
		for(j=0; j<np; ++j) {
			s1=0.0;
			for(k=0; k<np; ++k) s1=s1+wk[k*np+j]*wk1[k];
			r1=reps*fabs(wk1[j+np]); if(aeps>r1) r1=aeps;
			if(fabs(s1) > r1) qchk=0;
			x[j]=wk1[j+np]-s1;
			wk1[j+2*np]=-s1;
		}
		if(qchk==1) {free(intc); free(wk1); free(wk); return 0;}

		fcn(np,x,f);

/*	Updating the inverse matrix */
		for(j=0; j<np; ++j) {
			wk1[j+3*np]=f[j]-wk1[j];
			wk1[j]=f[j];
			wk1[j+np]=x[j];
		}

		ss=0.0;
		for(j=0; j<np; ++j) {
			s1=0.0; s2=0.0;
			for(k=0; k<np; ++k) {
				s1=s1+wk1[k+2*np]*wk[j*np+k];
				s2=s2+wk[k*np+j]*wk1[k+3*np];
			}
			f[j]=s1;
			x[j]=s2-wk1[j+2*np];
			ss=ss+s1*wk1[j+3*np];
		}

		for(j=0; j<np; ++j) {
			for(k=0; k<np; ++k) wk[j*np+k]=wk[j*np+k]-x[k]*f[j]/ss;
		}
	}
	free(intc); free(wk1); free(wk);
	return 442;
}
