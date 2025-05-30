/*	Solution of a system of nonlinear equations using Broyden's method */

#include <stdio.h>
#include <math.h>

void fun(int n, double x[], double f[]);
int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int davidn_b(void (*fcn) (int ,double * , double * ), int np, double x[],
	double f[], double reps, double aeps);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double f[20],reps, aeps, x0;
	double x[8] = {-0.25,0.1,-0.25,0.35,-0.25,0.6,-0.25,0.85};

/*	Example 7.15 */

	aeps=1.e-9; reps=1.e-7;
	np=8;
	i=davidn_b(fun,np,x,f,reps,aeps);
	printf(" ier = %d  f=", i);
	for(i=0; i<np; ++i) printf(" %e",f[i]);
	printf("\n");

	return;
}



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



/*	Real root of a system of nonlinear equations using Davidenko's method
	For use with function BROYDN

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
		
	Error status is returned by the value of the function DAVIDN_B.
		0 value implies successful execution
		407 implies NP > NMAXD or NP < 1, in which case no
			calculations are done
		440 implies that THETA is nonzero even after 100 steps
		Other values may be set by functions NEWTON or BROYDN

	Function FCN(NP,X,F) or FCN(NP,X,F,DF)
	must be supplied by the user. Here X and F are arrays of length
	NP with F containing the value of vector function at point X. DF is
	an array of length NP*NP containing the Jacobian which is required
	by function NEWTON. If function BROYDN is used then DF is not
	required. The function FCN should parameterise the equations in
	terms of THETA as explained in section 7.16. A sample function is
	appended at the end of this file.
	DAVIDN is the version to be used with NEWTON.

	Required functions : BROYDN (or NEWTON), GAUELM, FCN
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NMAXD 200

/*	To pass on parameters to FCN this must be included before fcn */
double THETA, X0[NMAXD], F0[NMAXD];
int broydn(void fcn(int , double * ,double * ), int np,
	double x[], double f[], double reps, double aeps);

int davidn_b(void (*fcn) (int ,double * , double * ), int np, double x[],
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
/*	  wk=(double *) calloc((size_t) (np*np),sizeof(double)); 
	  fcn(np,x,F0,wk);
	  free(wk);
*/
	  fcn(np,x,F0);

	  for(i=0; i<np; ++i) X0[i]=x[i];
	  reps0=reps; if(0.001>reps0) reps0=0.001;
	  aeps0=aeps; if(1.0e-4>aeps0) aeps0=1.0e-4;

	  for(i=1; i<=100; ++i) {
		  printf(" Type THETA \n");
		  scanf(" %le",&THETA);
		  printf("  THETA = %e \n",THETA);
		  if(THETA == 0.0) {reps0=reps; aeps0=aeps;}

/*		  ier=newton(fcn,np,x,f,reps0,aeps0); */
		  ier=broydn(fcn,np,x,f,reps0,aeps0);
		  if(ier>0) return ier;
		  printf(" x =");
		  for(j=0; j<np; ++j) printf(" %e",x[j]);
		  printf("\n");
		  if(THETA == 0) return 0;
	  }

	  return 440;
}
/*
	double THETA, X0[NMAXD], F0[NMAXD];
	void fcn(int np, double x[], double f[])
	{

	fun(np,x,f);

	for(i=0; i<np; ++i) f[i]=f[i]-THETA*F0[i];
	}

*/




/*	Solution of a system of linear equations using Gaussian elimination
	with partial pivoting

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input/output) The matrix of coefficient of size LJ*N
	        A[i][j] is the coefficient of x_j in ith equation
	     	at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
	     	at output it will contain the solutions
	DET : (output) The determinant of the matrix
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter to specify the type of computation required
		If IFLG<=0, both elimination and solution are
			done and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUELM.
		0 value implies successful execution
		101 implies (N<=0 or N>LJ) 
		121 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg)

{
int i,j,k,km,l;
double r1,t1;

	if(n<=0 || n>lj) return 101;
 
	if((*iflg)<2) {
/*	Perform elimination  */
		*det=1.0;
		for(k=0; k<n-1; ++k) {
/*	Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			for(l=k; l<n; ++l)
				if(fabs(a[l*lj+k])>r1) {r1=fabs(a[l*lj+k]); km=l;}

			inc[k]=km;
			if(km != k) {
/*	Interchange the rows if needed  */
				for(l=k; l<n; ++l) 
				{t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;}
				*det=-(*det);
			}

			*det=(*det)*a[k*lj+k];
			if(a[k*lj+k]==0) return 121;

			for(l=k+1; l<n; ++l) {
				a[l*lj+k]=a[l*lj+k]/a[k*lj+k];
				for(i=k+1; i<n; ++i) a[l*lj+i]=a[l*lj+i]-a[l*lj+k]*a[k*lj+i];
			}
		}
		*det=(*det)*a[(n-1)*lj+n-1];
		inc[n-1]=n-1;
		if(a[(n-1)*lj+n-1]==0) return 121;

		if((*iflg)==1) {*iflg=2; return 0;}
		*iflg=2;
	}

/*	Solution for the num different right-hand sides  */
	for(j=0; j<num; ++j) {
/*	forward-substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k])
			{t1=x[j*lj+k]; x[j*lj+k]=x[j*lj+inc[k]]; x[j*lj+inc[k]]=t1;}
			for(l=k+1; l<n; ++l) x[j*lj+l]=x[j*lj+l]-a[l*lj+k]*x[j*lj+k];
		}

/*	back-substitution  */

		x[j*lj+n-1]=x[j*lj+n-1]/a[(n-1)*lj+n-1];
		for(k=n-2; k>=0; --k) {
			for(l=n-1; l>=k+1; --l) x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=x[j*lj+k]/a[k*lj+k];
		}
	}
	return 0;
}


void fun(int n, double x[], double f[])

{
	int i,j;
	double fi;

	for(i=0; i<n; ++i) {
		fi=1.0/((i+1)*(i+1));
		for(j=0; j<n; j=j+2) {
			if(i == 0) {
				fi=fi+x[j];
			}
			else {
				fi=fi+x[j]*pow(x[j+1], (double) i);
			}
		}

/*	Parameterisation for Davidenko's method */
		f[i]=fi-THETA*F0[i];
	}
	return;
}

