/*	To calculate weights and abscissas of a quadrature formula with
	specified weight function when recurrence relation for orthogonal
	polynomial is known.

	N : (input) Number of points in the required quadrature formula
	W : (output) Array of length N, which will contain the weights
	AB : (output) Array of length N containing the abscissas
	COF : (input) Array of length 3*N containing the coefficients of
		the recurrence relation for orthogonal polynomials
		P_i(x)=(COF[3*i]*x+COF[3*i+1])P_{i-1}(x) - COF[3*i+2]*P_{i-2}(x)
	RI0 : (input) The integral of weight function over the required interval.
		
	Error status is returned by the value of the function GAUSRC.
		0 value implies successful execution
		302 implies N<=0
		321 implies that some coefficient becomes imaginary
			during calculations.
		In both these cases calculations are abandoned.
		Other values may be set by TQL2

     Required functions : TQL2
*/
 
#include <math.h>
#include <stdlib.h>

int tql2(double *z, int n, int iz, double d[], double e[], double reps);

int gausrc(int n, double w[], double ab[], double *cof, double ri0)

{
	int i,j,ier;
	double r1,reps=1.0e-15;
	double *wk,*d,*e;

	if(n<=0) return 302;
	wk=(double *) calloc((size_t) (n*n), sizeof(double));
	d=(double *) calloc((size_t) n, sizeof(double));
	e=(double *) calloc((size_t) n, sizeof(double));

/*     Calculate the coefficients of symmetric tridiagonal matrix */
	for(i=0; i<n; ++i) {
		d[i]=-cof[3*i+1]/cof[3*i];
		if(i<n-1) {
			r1=cof[3*i+5]/(cof[3*i]*cof[3*i+3]);
			if(r1>=0.0) e[i+1]=sqrt(r1);
			else {free(e); free(d); free(wk); return 321;}
		}
		for(j=0; j<n; ++j) wk[j+i*n]=0.0;
		wk[i+i*n]=1.0;
	}

/*     Find eigenvalues and eigenvectors of the tridiagonal matrix */
	ier=tql2(wk,n,n,d,e,reps);
	if(ier>0) {free(e); free(d); free(wk); return ier;}

/*     Calculate the abscissas and weights */
	for(i=0; i<n; ++i) {
		ab[i]=d[i];
		w[i]=wk[i]*wk[i]*ri0;
	}
	free(e); free(d); free(wk);
	return 0;
}
