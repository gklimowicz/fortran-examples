/*	To calculate weights and abscissas of a quadrature formula with
	specified weight function.

	N : (input) Number of points in the required quadrature formula
	W : (output) Array of length N, which will contain the weights
	AB : (input/output) Array of length N containing the abscissas
		For Gaussian formulas (QGAUS=1) AB will be calculated, 
		while for other values of QGAUS abscissas must be supplied
	FMOM : (input) Name of the function to calculate the moments
		Function FMOM(I) should calculate integral of w(x)x**I
	QGAUS : (input) Parameter to decide type of formula to be obtained
		If QGAUS=1 a Gaussian formula is calculated. In this
		case both abscissas and weights are calculated.
		Otherwise an interpolatory formula is calculated.
		In this case only weights are calculated, while abscissas
		must be supplied.
		
	Error status is returned by the value of the function GAUSWT.
		0 value implies successful execution
		303 implies N<=0 or N>=NPMAX
		322 implies GAUELM failed to find coefficients of polynomial
		323 implies POLYR failed to find roots of polynomial
		324 implies GAUELM failed to find weights 

	Required functions : GAUELM, POLYR, LAGITR, FMOM, CABS, CSQRT, CDIV
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int polyr(int n, double a[], double cr[], double ci[], int qrefin);

int gauswt(int n, double w[], double ab[], double (*fmom) (int ), int qgaus)

{
	int i,j,k,lj,iflg,ier,qrefin;
	double det;
	double *a, *cof, *zero;
	int *inc;

	if(n <= 0) return 303;
	lj=n;
	a=(double *)calloc((size_t) (n*n), sizeof(double));
	inc=(int *) calloc((size_t) n, sizeof(int));

	if(qgaus == 1) {
		cof=(double *)calloc((size_t) (n+2), sizeof(double));
		zero=(double *)calloc((size_t) (2*n), sizeof(double));

/*	Calculating the coefficients of the orthogonal polynomial */
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) a[i*n+j]=fmom(i+j);
			cof[i]=-fmom(n+i);
		}
		iflg=0;
		ier=gauelm(n,1,a,cof,&det,inc,lj,&iflg);
		if(ier >100) {free(zero); free(cof); free(inc); free(a); return 322;}

/*	Find the roots of polynomial, which will be the abscissas */
		cof[n]=1.0;
		qrefin=1;
		polyr(n,cof,zero,&zero[n],qrefin);
		if(ier >100) {free(zero); free(cof); free(inc); free(a); return 323;}
		for(i=0; i<n; ++i) ab[i]=zero[i];
		free(zero); free(cof);
	}


/*	Calculate the weights */
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) {
			if(i==0) a[i*n+j]=1.0;
			else a[i*n+j]=ab[j]*a[i*n+j-n];
		}
		w[i]=fmom(i);
	}
	iflg=0;
	ier=gauelm(n,1,a,w,&det,inc,lj,&iflg);
	free(inc); free(a);
	if(ier >100) return 324;
	return 0;
}
