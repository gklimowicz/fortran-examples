/*	To calculate rational function minimax approximation of tabulated functions */

#include <stdio.h>
#include <math.h>

int simpx(double *a, int ia, int n, int m, int nv, int qf, int id[],
	int iv[], double aeps);
int minmax(int m, int k, int n, double a[], double x[], double f[],
	double eps, double *emax);

main()
{
	int i,i1,j,k,m,n, iflg, ier,np,nmax;
	double hh, rm,xl,xu,eps,epsm,emax,a[20],e[20],x[500],f[500];

/*	Exercise 10.57 */

	n=100; eps=1.e-6; epsm=1.e-5;
	for(i1=0; i1<99; ++i1) {
		printf("type m=degree of numerator,  k=degree of denominator,  n=no. of data points\n");
		printf("                               (quits when n<=0)\n");
		scanf(" %d %d %d", &m,&k,&n);
		if(n<=0) return 0;

		np=m+k+1;
		printf(" m= %d     k = %d     n = %d \n",m,k,n);
		printf("type initial guess for coefficients \n");
		for(i=0; i<=m+k; ++i) scanf(" %le",&a[i]);
		printf(" Initial guess for coefficients :\n");
		for(i=0; i<=m+k; ++i) printf(" %e",a[i]);
		printf("\n");

/*	Generating the input data set */
		hh=2./(n-1.0);
		for(i=0; i<n; ++i) {
			x[i]=-1+i*hh;
			f[i]=atan(x[i]);
		}
			
		i=minmax(m,k,n,a,x,f,eps,&emax);
		printf(" ier = %d   Maximum error = %e \n coefficients in numerator : \n ",i,emax);
		for(i=k; i<=m+k; ++i) printf(" %e",a[i]);
		printf(" \n   coefficients in denominator :\n");
		for(i=0; i<k; ++i) printf(" %e",a[i]);
		printf(" \n");
	}
	return;
}




/*	To calculate coefficients of Rational function minimax approximation
	for a tabulated function

	M : (input) Required degree of polynomial in the numerator
	K : (input) Required degree of polynomial in the denominator
	N : (input) Number of points in the table of values, N> M+K+1
	A : (input/output) Array of length M+K+2 containing the coefficients
		of approximation. A[I-1] is the coefficient of x**I in
		the denominator, the constant term being 1. A[K+J] is
		the coefficient of x**J in the numerator. The 
		initial guess for coefficients must be supplied.
	X : (input) Array of length N containing the points at which
		function value is available
	F : (input) Array of length N containing the function values at X[I]
	EPS : (input) Required accuracy, the iteration is continued until
		change in EMAX is less than EPS
	EMAX : (output) Maximum error in the calculated approximation
		
	Error status is returned by the value of the function MINMAX.
		0 value implies successful execution
		615 implies that M<0, K<0 or N<M+K+2
			in this case no calculations are done
		634 implies that the iteration did not converge
			to the specified accuracy
		Other values may be set by SIMPX

	Required functions : SIMPX
*/

#include <math.h>
#include <stdlib.h>

int simpx(double *a, int ia, int n, int m, int nv, int qf, int id[],
	int iv[], double aeps);

int minmax(int m, int k, int n, double a[], double x[], double f[],
	double eps, double *emax)

{
	int i,j,nk,lj,qf,nc,n2,it,m3,nv,ier, nit=30;
	double xi,fd,fn,ei,fi,ti,s1,s2,s3,dif;
	int *iwk;
	double *wk;

	nk=m+k+2;
	if(m<0 || k<0 || n<nk) return 615;

	lj=nk+3; nv=3*n;
	m3=nk+1;    /*	The number of constraints */
	qf=1;
	nc=nv+m3; n2=2*n;

	iwk=(int *) calloc((size_t) nc, sizeof(int));
	wk=(double *) calloc((size_t) (lj*(nv+1)), sizeof(double));

/*	Loop for iteration */
	for(it=1; it<=nit; ++it) {

/*	Finding the maximum error in approximation */
		*emax=0.0;
		for(i=0; i<n; ++i) {
			xi=x[i];
			fd=0.0;
			for(j=k-1; j>=0; --j) fd=fd*xi+a[j];
			fd=fd*xi+1.0;
			wk[lj-1+i*lj]=fd;

			fn=a[m+k];
			for(j=m+k-1; j>=k; --j) fn=fn*xi+a[j];
			ei=f[i]-fn/fd;
			if(fabs(ei)>(*emax)) *emax=fabs(ei);
		}

/*	Setting up the tableau for Simplex algorithm */
		for(i=0; i<n; ++i) {
			xi=x[i];
			fd=wk[lj-1+i*lj];
			fi=f[i];

			iwk[i+m3+1]=i;
			iwk[i+n+m3+1]=i+n;
			iwk[i+n2+m3+1]=i+n2;
			wk[(i+1)*lj]=-fi+(*emax)-(*emax)*fd;
			wk[(i+n+1)*lj]=fi+(*emax)-(*emax)*fd;
			wk[(i+n2+1)*lj]=1.0;
			wk[1+(i+1)*lj]=fd;
			wk[1+(i+n+1)*lj]=fd;
			wk[1+(i+n2+1)*lj]=0.0;
			s3=0.0; ti=1.0;
			for(j=0; j<k; ++j) {
				ti=ti*xi;
				s3=s3+ti;
				wk[2+j+(i+1)*lj]=ti*((*emax)-fi);
				wk[2+j+(i+n+1)*lj]=ti*((*emax)+fi);
				wk[2+j+(i+n2+1)*lj]=ti;
			}
			s1=s3*((*emax)-fi)+1.0;
			s2=s3*((*emax)+fi)-1.0;

			wk[2+k+(i+1)*lj]=1.0;
			wk[2+k+(i+n+1)*lj]=-1.0;
			wk[2+k+(i+n2+1)*lj]=0.0;
			ti=1.0;
			for(j=0; j<m; ++j) {
				ti=ti*xi;
				s1=s1+ti;
				s2=s2-ti;
				wk[3+k+j+(i+1)*lj]=ti;
				wk[3+k+j+(i+n+1)*lj]=-ti;
				wk[3+k+j+(i+n2+1)*lj]=0.0;
			}
			wk[3+k+m+(i+1)*lj]=-s1;
			wk[3+k+m+(i+n+1)*lj]=-s2;
			wk[3+k+m+(i+n2+1)*lj]=-s3;
		}

		wk[0]=0.0; wk[1]=1.0;
		iwk[1]=3*n;
		for(j=2; j<m+k+4; ++j) {
			iwk[j]=3*n+j-1;
			wk[j]=0.0;
		}

		ier=simpx(wk,lj,nc,m3,nv,qf,iwk,&iwk[m3],eps);
		if(ier >0) {free(wk); free(iwk); return ier;}

/*	The maximum error */
		ei=wk[0];
/*	Obtaining the coefficients from the tableau */
		for(i=m3+1; i<=nc; ++i) {
			if(iwk[i]>nv) a[iwk[i]-nv-1]=wk[(i-m3)*lj];
		}
		for(i=1; i<=m3; ++i) {
			if(iwk[i]>nv) a[iwk[i]-nv-1]=0.0;
		}
		for(i=0; i<=m+k; ++i) a[i]=a[i]-a[m+k+1];
		dif=(*emax)-ei;
		*emax=ei;

/*	The convergence test */
		if(dif<eps || k==0) {free(wk); free(iwk); return 0;}
	}

	free(wk); free(iwk);
	return 634;
}



/*	To solve a linear programming problem in the standard form
		 using the simplex method

	A : (input/output) Array of length IA*(N-M+1) containing
		the tableau of simplex algorithm
		A[I+1][0]=c_i, the cost coefficients
		Columns 1 to M contain constraints with A[0][j+1]=b_j and A[i+1][j+1]=a_i
		Column M+1 contains the cost coefficients for auxiliary problem
		when QF!=1.
	IA : (input) Second dimension of array A as declared in the calling
		function. IA >= M+2
	N : (input) Number of variables, each is constrained to be >=0
	M : (input) Number of constraints of form a^T X = b_i >= 0
	NV : (input) Number of variables, excluding the artificial variables
	QF : (input) Parameter to decide which objective function
		is to be minimised.
		If QF=1  the main objective function specified by
			the first column of A is to be minimised.
		Otherwise the auxiliary objective function specified by
			the last column of A is to be minimised.
	ID : (input/output) Integer array of length M+1 which contains
		information about interchange of variables on LHS
	IV : (input/output) Integer array of length N-M+1 which contains
		information about interchange of variables on RHS
	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
		assumed to be zero
		
	Error status is returned by the value of the function SIMPX.
		0 value implies successful execution
		57 implies that the objective function is unbounded from below
		531 implies that the simplex algorithm failed to find
			the optimal feasible vector

	Required functions : None
*/

#include <math.h>

int simpx(double *a, int ia, int n, int m, int nv, int qf, int id[],
	int iv[], double aeps)

{
	int i,j,k,jf,m1,n1,l,l1,it,nit=20;
	double rmin,rmax,r1;

/*	Minimise the objective function in the first column */
	if(qf==1) {jf=0; m1=m+1;}
/*	Minimise the objective function in the last column */
	else {jf=m+1; m1=jf+1;}
	n1=n-m+1;

	for(it=1; it<= nit*(n+m); ++it) {
/*	Find the minimum of the reduced cost coefficients */
		rmin=a[jf+ia]; k=1;
		for(j=2; j<n1; ++j) {
			if(a[jf+j*ia] < rmin) {rmin=a[jf+j*ia]; k=j;}
		}

		if(rmin>=0.0) {
/*	The objective function cannot be reduced further */
			if(qf==1 || a[jf]<-aeps) return 0;
/*	Check if any artificial variable is on the LHS */
			for(i=1; i<=m; ++i) {
				if(id[i]>=nv && fabs(a[i])<=aeps) {
					a[i]=0.0; l=i; k=0;
					rmax=0.0;
					for(j=1; j<n1; ++j) {
						if(fabs(a[i+j*ia])>rmax && iv[j]<=nv) {
							rmax=fabs(a[i+j*ia]);
							k=j;
						}
					}
/*	To exchange the artificial variable */
					if(k>0) goto exch;
				}
			}
			return 0;
		}

/*	Finding allowed change in the Kth variable */
		rmin=0.0; l=0;
		for(j=1; j<=m; ++j) {
			if(a[j+k*ia]>aeps) {
				r1=a[j]/a[j+k*ia];
				if(r1<rmin || l==0) {rmin=r1; l=j;}
			}
		}
		if(l==0) return 57;

/*	exchange the variables */
exch:	l1=id[l]; id[l]=iv[k]; iv[k]=l1;
		for(j=0; j<n1; ++j) {
			if(j!=k) {
				for(i=0; i<m1; ++i) {
					if(i != l) a[i+j*ia]=a[i+j*ia]-a[i+k*ia]*a[l+j*ia]/a[l+k*ia];
				}
			}
		}

		for(j=0; j<n1; ++j) {
			if(j != k) a[l+j*ia]=a[l+j*ia]/a[l+k*ia];
		}
		for(i=0; i<m1; ++i) {
			if(i != l) a[i+k*ia]= -a[i+k*ia]/a[l+k*ia];
		}
		a[l+k*ia]=1.0/a[l+k*ia];
	}

	return 531;
}
