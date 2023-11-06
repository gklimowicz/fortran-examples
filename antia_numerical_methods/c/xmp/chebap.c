/*	To find coefficients of rational function approximation in terms
	of Chebyshev polynomials using the known Chebyshev expansion */

#include <stdio.h>
#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int chebap(int m, int k, double a[], double c[]);

main()
{
	int i,i1,j,k,m, iflg, ier,np,nmax;
	double hh, rm, a[20],c[20];

/*	Example 10.9 */

	for(i1=0; i1<99; ++i1) {
		printf("type m=degree of numerator,  k=degree of denominator\n");
		printf("                        (quits when m+k<=0)\n");
		scanf(" %d %d", &m,&k);
		if(m+k<=0) return 0;

/*	Coefficients of Chebyshev expansion of Arc Tan(x) */
		rm=sqrt(2.0)-1.0;
		hh=1.0;
		for(i=0; i<=m+2*k; i += 2) {
			c[i]=0.0;
			c[i+1]=2.*hh*pow(rm, (double) (i+1))/(i+1.0);
			hh=-hh;
		}

		i=chebap(m,k,a,c);
		printf(" ier= %d    m= %d   k =%d\n  Coefficients in numerator :\n ",i,m,k);
		for(i=k; i<=m+k; ++i) printf(" %e",a[i]);
		printf("\n  Coefficients in denominator : \n");
		for(i=0; i<k; ++i) printf(" %e",a[i]);
		printf(" \n");
	}
	return;
}



/*	To calculate coefficients of rational function Chebyshev approximation
	using the known coefficients of Chebyshev expansion

	M : (input) Degree of polynomial in the numerator
	K : (input) Degree of polynomial in the denominator
	A : (output) Array of length M+K+1 containing the coefficients
		of rational function approximation. A[I-1] is the coefficient
		of T_I(x) in the denominator, the constant term being 1.
		A[K+J] is the coefficient of T_J(x) in the numerator
	C : (input) Array of length M+2K+1 containing the coefficients
		of Chebyshev expansion. C[I] is the coefficient of T_I(x).
		Coefficient of T_0(x) is C[0]/2
		
	Error status is returned by the value of the function CHEBAP.
		0 value implies successful execution
		612 implies that M<0 or K<0, no calculations are done
		Other values may be set by GAUELM

	Required functions : GAUELM
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int chebap(int m, int k, double a[], double c[])

{
	int i,j,iflg,ier;
	double ai,det;
	int *iwk;
	double *wk;

	if(m<0 || k<0) return 612;

	if(k>0) {
/*	Setting up the coefficients of matrix of linear equations
	to calculate the coefficients in the denominator */
	wk=(double *) calloc((size_t) (k*k), sizeof(double));
	iwk=(int *) calloc((size_t) k, sizeof(int));
	for(i=0; i<k; ++i) {
		for(j=0; j<k; ++j) wk[i+j*k]=0.5*(c[m+2+i+j]+c[abs(m+j-i)]);

/*	The right-hand side vector for linear equations */
		a[i]=-c[m+i+1];
	}

	iflg=0;
/*	Solve the system of linear equations */
	ier=gauelm(k,1,wk,a,&det,iwk,k,&iflg);
	free(iwk); free(wk);
	if(ier>0) return ier;
	}

/*	Calculating the coefficients in the numerator */
	for(i=0; i<=m; ++i) {
		ai=c[i];
		for(j=0; j<k; ++j) ai=ai+0.5*a[j]*(c[i+j+1]+c[abs(i-1-j)]);
		a[k+i]=ai;
	}
	a[k]=0.5*a[k];
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
