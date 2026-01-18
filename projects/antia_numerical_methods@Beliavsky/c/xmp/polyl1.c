/*	to calculate polynomial L1-approximation of a tabulated function */

#include <stdio.h>
#include <math.h>

double fun(double x);
double rangau(double *seed);
int simpl1(double *a, int ia, int n, int m, int id[], int iv[], double aeps);
int polyl1(int m, int n, double a[], double x[], double f[], double eps,
	double *esum);
int linl1(int m, int n, double a[], double f[], double *g, int ig,
	double eps, double *esum);

double seed;

main()
{
	int i,i1,j,n,m, id, iflg, ier,np,m1;
	double hh, x[100], f[100], w[100][100], reps, fx, xx, a[20],y[200],esum;

/*	Example 10.13 */

	seed=2; id=100; reps=1.e-9;
	for(i1=0; i1<99; ++i1) {
		printf("type  m = degree,    n = no. of data points \n");
		printf("                       (quits when n<=0)\n");
		scanf(" %d %d", &m,&n);
		if(n<=0) return 0;

/*	Generating the input Data set */
		hh=1.0/(n-1.0);
		for(i=0; i<n; ++i) {
			x[i]=i*hh;
			f[i]=fun(x[i]);

/*	Basis functions for use with linl1 */

			w[i][0]=1.0;
			for(j=1; j<=m; ++j) w[i][j]=pow(x[i],(double) j);
		}

	i=polyl1(m,n,a,x,f,reps,&esum);
	printf(" ier = %d    no. of pts = %d    degree = %d   Error sum = %e\n", i,n,m,esum);
	printf(" coefficients : ");
	for(i=0; i<=m; ++i) printf(" %e",a[i]);
	printf("\n");

	m1=m+1;
	i=linl1(m1,n,a,f,&w[0][0],id,reps,&esum);
	printf(" ier = %d  no. of pts = %d   no. of basis functions = %d    Error sum = %e\n", i,n,m1,esum);
	printf(" coefficients : ");
	for(i=0; i<=m; ++i) printf(" %e",a[i]);
	printf("\n");
	}
	return;
}

double fun(double x)

{
	double fx;
	fx=((231*x*x-315)*x*x+105)*x*x-5+rangau(&seed)*1.e-3;
	return fx;

}



/*	To calculate coefficients of linear L_1 approximation in terms of
	specified basis functions for a tabulated function

	M : (input) Number of basis functions in the required approximation
	N : (input) Number of points in the table of values. N> M
	A : (output) Array of length M+1 containing the coefficients
		of approximation. A[I] is the coefficient of Phi_I(x) in the
		approximation
	F : (input) Array of length N containing the function values
	G : (input) Array of length IG*N, containing the values of
		basis functions at each point in table. G[j][i] is the
		value of Phi_i(x) at jth tabular point.
	IG : (input) Second dimension of array G as declared in the calling
		function, IG >= M
	EPS : (input) Estimate of roundoff error, passed on to SIMPL1
	ESUM : (output) L_1 norm of the residual at the calculated approximation
		
	Error status is returned by the value of the function LINL1.
		0 value implies successful execution
		616 implies that M<1  or N<M+1
			in this case no calculations are done
		Other values may be set by SIMPL1

	Required functions : SIMPL1
*/
	
#include <math.h>
#include <stdlib.h>

int simpl1(double *a, int ia, int n, int m, int id[], int iv[], double aeps);

int linl1(int m, int n, double a[], double f[], double *g, int ig,
	double eps, double *esum)

{
	int i,j,lj,nv,m3,nc,ier;
	double fi,si,s,ti;
	int *iwk;
	double *wk;

	if(m<=0 || n<m+1) return 616;

	lj=n+2; nv=m+1;
	m3=n;		/*	The number of constraints */

	iwk=(int *) calloc((size_t) (n+m+3), sizeof(int));
	wk=(double *) calloc((size_t) (lj*(nv+1)), sizeof(double));

/*	Setting up the tableau for simplex algorithm */
	for(j=0; j<m+2; ++j) wk[j*lj]=0.0;
	for(i=0; i<n; ++i) {
		fi=f[i];
		si=1.0; if(fi<0.0) si=-1.0;
		iwk[i+1]=i*si;
		wk[i+1]=fi*si;
		wk[0]=wk[0]-fi*si;
		s=0.0;

		for(j=0; j<m; ++j) {
			ti=g[j+i*ig]*si;
			s=s+ti;
			wk[i+1+(j+1)*lj]=ti;
			wk[(j+1)*lj]=wk[(j+1)*lj]-ti;
		}
		wk[i+1+(m+1)*lj]=-s;
		wk[(m+1)*lj]=wk[(m+1)*lj]+s;
	}

	for(j=0; j<m+1; ++j) iwk[n+1+j]=n+j;

	nc=nv+m3;
	ier=simpl1(wk,lj,nc,m3,iwk,&iwk[m3],eps);
	
/*	L_1 norm of the residual */
	*esum=-wk[0];
/*	Finding the coefficients from the tableau */
	for(i=m3+1; i<=nc; ++i) {
		if(iwk[i]>=n) a[iwk[i]-n]=0.0;
	}
	for(i=1; i<=m3; ++i) {
		if(iwk[i]>=n) a[iwk[i]-n]=wk[i];
	}
	for(i=0; i<m; ++i) a[i]=a[i]-a[m];
	
	free(wk); free(iwk);
	return ier;
}




/*	To calculate coefficients of polynomial L_1 approximation
	for a tabulated function

	M : (input) Required degree of polynomial
	N : (input) Number of points in the table of values. N> M+1
	A : (output) Array of length M+2 containing the coefficients
		of approximation. A[I] is the coefficient of x**I
	X : (input) Array of length N containing the points at which
		function value is available
	F : (input) Array of length N containing the function values at X[I]
	EPS : (input) Estimate of roundoff error, passed on to SIMPL1
	ESUM : (output) L_1 norm of the residual at the calculated approximation
		
	Error status is returned by the value of the function POLYL1.
		0 value implies successful execution
		616 implies that M<0 or N<M+2
			in this case no calculations are done
		Other values may be set by SIMPL1

	Required functions : SIMPL1
*/

#include <math.h>
#include <stdlib.h>

int simpl1(double *a, int ia, int n, int m, int id[], int iv[], double aeps);

int polyl1(int m, int n, double a[], double x[], double f[], double eps,
	double *esum)

{
	int i,j,lj,nv,m3,nc,ier;
	double fi,si,s,ti;
	int *iwk;
	double *wk;

	if(m<0 || n<m+2) return 616;

	lj=n+2; nv=m+2;
	m3=n;		/*	The number of constraints */

	iwk=(int *) calloc((size_t) (n+m+3), sizeof(int));
	wk=(double *) calloc((size_t) (lj*(nv+1)), sizeof(double));

/*	Setting up the tableau for simplex algorithm */
	for(j=0; j<m+3; ++j) wk[j*lj]=0.0;
	for(i=0; i<n; ++i) {
		fi=f[i];
		si=1.0; if(fi<0.0) si=-1.0;
		iwk[i+1]=i*si;
		wk[i+1]=fi*si;
		wk[0]=wk[0]-fi*si;
		wk[i+1+lj]=si;
		wk[lj]=wk[lj]-si;
		s=si; ti=si;

		for(j=0; j<m; ++j) {
			ti=ti*x[i];
			s=s+ti;
			wk[i+1+(j+2)*lj]=ti;
			wk[(j+2)*lj]=wk[(j+2)*lj]-ti;
		}
		wk[i+1+(m+2)*lj]=-s;
		wk[(m+2)*lj]=wk[(m+2)*lj]+s;
	}

	for(j=0; j<m+2; ++j) iwk[n+1+j]=n+j;

	nc=nv+m3;
	ier=simpl1(wk,lj,nc,m3,iwk,&iwk[m3],eps);
	
/*	L_1 norm of the residual */
	*esum=-wk[0];
/*	Finding the coefficients from the tableau */
	for(i=m3+1; i<=nc; ++i) {
		if(iwk[i]>=n) a[iwk[i]-n]=0.0;
	}
	for(i=1; i<=m3; ++i) {
		if(iwk[i]>=n) a[iwk[i]-n]=wk[i];
	}
	for(i=0; i<=m; ++i) a[i]=a[i]-a[m+1];
	
	free(wk); free(iwk);
	return ier;
}



/*	To solve a linear programming problem in the standard form arising
		 in L_1 minimisation problems using the simplex method

	A : (input/output) Array of length IA*(N-M+1) containing
		the tableau of simplex algorithm
		A[I+1][0]=c_i, the cost coefficients
		Columns 1 to M contain constraints with A[0][j]=b_j and A[i+1][j]=a_i
	IA : (input) Second dimension of array A as declared in the calling
		function. IA >= M+2
	N : (input) Number of variables, each is constraint to be >=0
	M : (input) Number of constraints of form a^T X = b_i >= 0
	ID : (input/output) Integer array of length M+1 which contains
		information about interchange of variables on LHS
	IV : (input/output) Integer array of length N-M+1 which contains
		information about interchange of variables on RHS
	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
		assumed to be zero
		
	Error status is returned by the value of the function SIMPL1.
		0 value implies successful execution
		63 implies that the objective function is unbounded from below
		635 implies that the simplex algorithm failed to find
			the optimal feasible vector

	Required functions : None
*/

#include <math.h>

int simpl1(double *a, int ia, int n, int m, int id[], int iv[], double aeps)

{
	int i,j,k,k1,l,l1,jf,m1,n1,it, nit=20;
	double rmin,r1,aj;

	jf=0;
	m1=m+1; n1=n-m+1;

/*	The simplex iteration */
	for(it=1; it<=nit*(n+m); ++it) {
/*	Finding the minimum reduced cost coefficient */
		rmin=0.0; k=0;
		for(j=1; j<n1; ++j) {
			if(a[jf+j*ia]<rmin) {rmin=a[jf+j*ia]; k=j;}
			else if(iv[j]<=m-1 && 2-a[jf+j*ia]< rmin) {rmin=2-a[jf+j*ia]; k=-j;}
		}
		if(rmin>= -aeps) return 0;

/*	Finding the pivot element */
		k1=abs(k);
		rmin=0.0;
		l=0;
		for(j=1; j<=m; ++j) {
			aj=a[j+k1*ia];
			if(k<0) aj=-aj;
			if(aj > aeps) {
				r1=a[j]/aj;
				if(r1<rmin || l==0) {rmin=r1; l=j;}
			}
		}
		if(l==0) return 63;  /*	The objective function is unbounded from below */

		if(k<0) {
			a[jf+k1*ia]=-2+a[jf+k1*ia];
			iv[k1]=-iv[k1];
		}

/*	Exchange the variables */
		l1=id[l]; id[l]=iv[k1]; iv[k1]=l1;
		for(j=0; j<n1; ++j) {
			if(j != k1) {
				r1=a[l+j*ia]/a[l+k1*ia];
				for(i=0; i<m1; ++i) {
					if(i != l) a[i+j*ia]=a[i+j*ia]-a[i+k1*ia]*r1;
				}
			}
		}

		r1=fabs(a[l+k1*ia]);
		for(j=0; j<n1; ++j) {
			if(j != k1) a[l+j*ia]=a[l+j*ia]/r1;
		}
		for(i=0; i<m1; ++i) {
			if(i != l) a[i+k1*ia]= -a[i+k1*ia]/a[l+k1*ia];
		}
		a[l+k1*ia]=1./r1;
	}

	return 635;
}



/*	To generate random numbers with Gaussian probability distribution
	It generates random numbers with zero mean and variance of 1.
	
	SEED : (input/output) real seed, it should be positive and
		less than AN. It is updated by the function and should
		not be modified between two calls, unless a fresh
		sequence is required

	Required functions : None
	
	THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
	VERSION AS THE SEED IS NOW DOUBLE INSTEAD OF INT.

*/

#include <math.h>


double rangau(double *seed)

{
	int n2;
	double am=2147483648.0, a=45875.0, ac=453816693.0, an=2147483647.0, r1, rn;

	rn=a*(*seed)+ac; n2=rn/am; r1=rn-am*n2;
	if(*seed==0.0) *seed=0.1;
	rn=sqrt(2.0*log(an/(*seed)))*cos(2.0*M_PI*r1/an);
	*seed=r1;
	return rn;
}

