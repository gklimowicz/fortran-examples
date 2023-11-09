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
