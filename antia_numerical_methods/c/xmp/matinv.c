/*	To calculate inverse of a square matrix using Gaussian elimination */

#include <stdio.h>
#include <math.h>

int gauelm(int n, int num, double *a, double *x, double * det,
       int inc[], int lj, int * iflg);
int matinv(int n, int ia, double *a, double *ai);

main()
{
	int i,i1,j,n, num, iflg, lj, inc[20];
	double r1, det, x[2][20], wk[60], am[20][20], bm[20][20];

	lj=20;
	for(i1=0; i1<99; ++i1) {
		printf("type n = order of matrix    (quits when n<=0)\n");
		scanf("%d", &n);
		if(n<=0) return 0;
		printf("type the matrix row-wise \n");

/*	Read the matrix */
		for(i=0; i<n; ++i) {
			printf(" %d th row \n",i);
			for(j=0; j<n; ++j) scanf(" %le",&am[i][j]);
		}
		i=matinv(n,lj,&am[0][0],&bm[0][0]);
		printf(" ier = %d   n = %d \n The inverse matrix is \n", i,n);
		for(j=0; j<n; j++) {
			for(i=0; i<n; ++i) printf(" %e ",bm[j][i]);
			printf("\n");
		}
	}
	return;
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



/*	To calculate inverse of a square matrix

	N : (input) Order of matrix
	IA : (input) The second dimension of arrays A and AI as specified
		in the calling function
	A : (input) Array of length IA*N containing the matrix
	AI : (output) Array of length IA*N which will contain the
		calculated inverse of A
		
	Error status is returned by the value of the function MATINV.
		0 value implies successful execution
		Nonzero values may be set by function GAUELM

	It is possible to use CROUT instead of GAUELM for calculating
	the triangular decomposition.

	Required functions : GAUELM (or CROUT)
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double *a, double *x, double *det,
       int inc[], int lj, int * iflg);
int crout(int n, int num, double *a, double *x, double *det,
       double *idet, int inc[], int lj, int * iflg);

int matinv(int n, int ia, double *a, double *ai)

{
	int i,j,iflg,ier,idet;
	double det,r;
	int *iwk;

	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) ai[j+i*ia]=0.0;
		ai[i+i*ia]=1.0;
	}

	iwk=(int *) calloc((size_t) n, sizeof(int));
	iflg=0;
	ier=gauelm(n,n,a,ai,&det,iwk,ia,&iflg);
/*	ier=crout(n,n,a,ai,&det,&idet,iwk,ia,&iflg);  */
	free(iwk);
/*	transpose the matrix to get the correct inverse in normal form */
	for(i=0; i<n; ++i) {
		for(j=0; j<i; ++j) {
			r=ai[j+i*ia];
			ai[j+i*ia]=ai[i+j*ia];
			ai[i+j*ia]=r;
		}
	}
	return ier;
}
