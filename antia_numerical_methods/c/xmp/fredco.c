/*	Solution of Fredholm equation using Collocation method */

#include <stdio.h>
#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax);
int fredco(int n, double a, double b, double f[], double x[], double reps,
	double aeps, int iq, int it);

double funk(double x);
double fg(double x);
double fker(double x, double t);
double phi(int i, double x);
double psi(int i, double x);

double al;

main()
{
	int i,i1,j,n,m, id, iflg, ier,np,nmax;
	double hh, x[100], wt[100],f[100],fc[100],tn,erc,reps,t0,fi,aeps;

/*	Example 13.1 */

	reps=1.e-7; aeps=1.e-8;
	t0=0.0; tn=1.0; iflg=2; id=0;

/*	id =0 for evaluating the integrals using adpint
	id=1  for evaluating the integrals using user supplied function psi */

	for(i1=0; i1<99; ++i1) {
		printf("type n=no. of points    id = 0/1  (quits when n<=0)\n");
		scanf(" %d %d",&n, &id);
		if(n<=0) return 0;

/*	Set up the collocation points */
		for(i=0; i<=n; ++i) x[i]=(i+1.0)/n;

		i=fredco(n,t0,tn,f,x,reps,aeps,id,iflg);
		printf(" ier = %d     no. of points =  %d     id = %d\n",i,n,id);
		printf("    coefficients :");
		for(i=0; i<n; ++i) printf(" %e ",f[i]);
		printf("\n");

/*	Calculate the solution at required points using the coefficients of expansion */

		printf("    x          solution \n");
		for(i=0; i<6; ++i) {
			hh=i*0.2; fi=0.0;
			for(j=0; j<n; ++j) fi=fi+phi(j,hh)*f[j];
			printf(" %e  %e \n",hh,fi);
		}

	}
	return;
}

/*	The right hand side function */

double fg(double x)

{
	return -x*x*x; 
}

/*	The kernel */

double fker(double x, double t)

{
	if(t>x) return -x*(1.0-t);
	else return -t*(1.0-x); 
}

/*	The integrals involving basis functions */
double psi(int i, double x)

{
		return -x*(1.0-pow(x, (double ) 2*i+2))/((2*i+2)*(2*i+3.0));
}

/*	The basis functions,  x**(2i+1), only odd powers are used */

double phi(int i, double x)

{
	if(x==0.0) return 0.0;
	else return pow(x, (double) 2*i+1);
}



/*	To integrate a function over finite interval using adaptive control
	of step size

	RINT : (output) Calculated value of the integral
	XL : (input) The lower limit
	XU : (input) The upper limit
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	NPT : (output) Number of function evaluations
	NMAX : (input/output) Maximum number of function evaluations to be tried
		If NMAX<=0 it is set to MAXPT (=100000)

	Error status is returned by the value of the function ADPINT.
		0 value implies successful execution
		31 implies specified accuracy was not achieved on
			at least one subinterval
		32 implies that this failure occurred more than IFMAX (=5) times
		325 implies that function failed to attain required
			accuracy using NMAX function evaluations
		In all cases DIF will contain the estimated accuracy

		Function F(X) must be supplied by the user.

	Required functions : KRONRD (or GAUS16), F
*/

#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));
int gaus16(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ));


int adpint(double *rint, double xl, double xu, double reps, double aeps,
	double *dif, double (*f) (double ), int *npt, int *nmax)

{
	int i,j,k,ier,ifail,iu,np,q, ipmax=100, ifmax=5, maxpt=100000;
	double aepsl,rl,ru,rm,fint,dif0,r1,xu1[100];
 
	*rint=0.0; *dif=0.0;
	if(xl == xu) return 0;
	ifail=0; ier=0;
	if(*nmax <= 0) *nmax = maxpt;
	aepsl= aeps;
	*npt=0;
	rl=xl; ru=xu; iu=0;

	do {
/*	To evaluate the integral over [rl,ru] */
		i=kronrd(&fint,rl,ru,&dif0,&np,f);
/*		i=gaus16(&fint,rl,ru,&dif0,&np,f);  */
		*npt = (*npt)+np;
		rm=0.5*(rl+ru);
/*	q=.TRUE. if the interval cannot be divided further */
		q = ((iu>=ipmax) || (rm == rl) || (rm == ru));
		r1=fabs(fint)*reps; if(aepsl>r1) r1=aepsl;

		if( (dif0 < r1) || q) {
/*	Accept the value of fint if adequate convergence or if the interval
	cannot be subdivided further */
			*rint=(*rint)+fint;
			*dif=(*dif)+dif0;
			r1=fabs(*rint)*reps; if(aepsl>r1) r1=aepsl;

			if(q && (dif0>r1)) {
/*	Integration fails to converge on this subinterval. Go to the next subinterval */
				ier=31; ifail=ifail+1;
/*	If failure is frequent then adjust the convergence criterion. */
				if(ifail > ifmax) {ier=32; aepsl=(*dif)*0.5;}
			}

/*	If all subintervals are exhausted then return */
			if(iu <= 0) return ier;

/*	otherwise try next subinterval */
			rl=ru; ru=xu1[iu]; iu=iu-1;
		}
		else {

/*	Subdivide the current interval and try again */
			iu=iu+1;
			xu1[iu]=ru;
			ru=rm;
		}
	} while(*npt < (*nmax));


/*	If the number of function evaluations has exceeded the limit then
	try a last call to estimate the integral over the remaining interval */
	ru=xu;
	i=kronrd(&fint,rl,ru,&dif0,&np,f);
/*	i=gaus16(&fint,rl,ru,&dif0,&np,f); */
	*npt=(*npt)+np;
	*rint=(*rint)+fint;
	*dif=(*dif)+dif0;
	return 325;
}



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



/*	To integrate a function over a finite interval using Gauss-Kronrod formula
	For use with ADPINT

	RI : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	DIF : (output) estimated (absolute) error achieved by the function
	N : (output) Number of function evaluations used
	F : (input) Name of the function to calculate the integrand

	Function F(X) must be supplied by the user

	The returned value is always zero.

	Required functions : F
*/

#include <math.h>

int kronrd(double *ri, double a, double b, double *dif, int *n,
	double (*f)(double ))

{
	int k;
	double at,bt,fbt,r1,f1,f2;

/*	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
	WK7 are the weights for these points in Kronrod formula
	WK15 and AK15 are the weights and abscissas for the remaining points
	in Kronrod formula.
	Because of symmetry only half the points are given.
*/

	double w7[4] = {0.12948496616886969327e0, 0.27970539148927666790e0,
                        0.38183005050511894495e0, 0.41795918367346938775e0};
	double a7[4] = {0.94910791234275852452e0, 0.74153118559939443986e0,
                        0.40584515137739716690e0, 0.0};
	double wk7[4] = {0.06309209262997855329e0, 0.14065325971552591874e0,
                         0.19035057806478540991e0, 0.20948214108472782801e0};
	double wk15[4] = {0.02293532201052922496e0, 0.10479001032225018383e0,
                          0.16900472663926790282e0, 0.20443294007529889241e0};
	double ak15[4] = {0.99145537112081263920e0, 0.86486442335976907278e0,
                          0.58608723546769113029e0, 0.20778495500789846760e0};

	at=(b-a)/2.;
	bt=(b+a)/2.;
	fbt=f(bt);
	r1=w7[3]*fbt;
	*ri=wk7[3]*fbt;
	for(k=0; k<3; ++k) {
		f1=f(bt+at*a7[k]);
		f2=f(bt-at*a7[k]);
/*	7-point Gauss-Legendre formula */
		r1=r1+w7[k]*(f1+f2);
/*	15-point Kronrod formula */
		*ri=(*ri)+wk7[k]*(f1+f2);
	}

	for(k=0; k<4; ++k) *ri=(*ri)+wk15[k]*(f(bt+at*ak15[k]) + f(bt-at*ak15[k]));

	*ri=(*ri)*at;
	r1=r1*at;
	*dif=fabs(*ri-r1);
	*n=15;
	return 0;
}


/*	Function to calculate the integrand for calculating
	PSI(I,X) as required by function FREDCO.

	Function FKER(X,T) is the kernel K(x,t) and PHI(I,T) is the Ith
	basis function, phi_i(t). The argument X and I are passed through
	global variables XFRED and IFRED respectively.

	Functions FKER(X,T) and PHI(I,T) must be supplied by the user

	Required functions : FKER, PHI
*/

/*	To pass parameters from function FREDCO put these statements before fredco */
double XFRED;
int IFRED;

double funk(double t)

{
	return fker(XFRED,t)*phi(IFRED,t);
}
