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
