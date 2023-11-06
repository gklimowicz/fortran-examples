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
