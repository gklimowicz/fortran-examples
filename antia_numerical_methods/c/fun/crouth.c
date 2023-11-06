/*	Solution of a system of linear equations using Crout's algorithm 
            with iterative refinement

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input) The matrix of coefficient of size LJ*N
	         A[i][j] is the coefficient of x_j in ith equation
	B : (output) Array of size LJ*N containing triangular decomposition of matrix A
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
			at output it will contain the solutions
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
			interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	REPS : (input) Required relative precision in solution
	IFLG : (input) Integer parameter which determines the type of computation
		required.
		If IFLG<=0, both elimination and solution are calculated
			and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
	WK : (output) Array of length NUM containing the estimated errors.
		WK[I] will contain estimated error in the Ith solution.
		
	Error status is returned by the value of the function CROUTH.
		0 value implies successful execution
		11 implies that iterative refinement did not converge
		102 implies (N<=0 or N>LJ) 
		122 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : CROUT
*/

#include <math.h>
#include <stdlib.h>

int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg);

int crouth(int n, int num, double *a, double *b, double *x, double *det,
	int *idet, int inc[], int lj, double reps, int *iflg, double wk[])

{
	int i,j,k,it,iflg1,ier,num1, nitr=10;
	double r1,r2,rp1,rp;
/* If long double is supported, use that for residuals */
	long double d1,d2;
	double *wk1, *wk2;

	if(n<=0 || n>lj) return 102;

	if(*iflg<2) {
/*	Preserving the matrix for calculating the residuals */
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) b[j*lj+i]= a[j*lj+i];
		}
		iflg1=1;
/*	Perform LU decomposition using crout */
		ier=crout(n,num,b,x,det,idet,inc,lj,&iflg1);
		if(ier>100) return ier;
	}

	if(*iflg==1) {*iflg=2; return 0;}
	*iflg=2;
	num1=1;
	wk1=(double *) calloc((size_t) n, sizeof(double));
	wk2=(double *) calloc((size_t) n, sizeof(double));

/*	Solving the systems with NUM different right-hand sides */
	for(j=0; j<num; ++j) {
		for(i=0; i<n; ++i) {
			wk1[i]= x[j*lj+i];
/*	Preserving the RHS for calculating residuals */
			wk2[i]= x[j*lj+i];
			x[j*lj+i]=0.0;
		}

		rp1=0.0;
/*	The iterative refinement */
		for(it=1; it<= nitr ; ++it) {
			ier=crout(n,num1,b,wk1,det,idet,inc,lj,&iflg1);
			if(ier>100) {free(wk2); free(wk1); return ier;}
			r1=0.0; r2=0.0;
			for(i=0; i<n; ++i) {
				if(fabs(wk1[i])>r1) r1=fabs(wk1[i]);
				x[j*lj+i]=x[j*lj+i]+wk1[i];
				if(fabs(x[j*lj+i])>r2) r2=fabs(x[j*lj+i]);
			}

/*	The error estimate  */
			if(it==2) wk[j]=r1;
			if(r2==0.0) goto nextit;
			rp=r1/r2;
			if(rp<reps) goto nextit;
			if(rp>rp1 && it>1) {ier=11; goto nextit;}
			rp1=rp;

/*	Calculating the residue  */
			for(i=0; i<n; ++i) {
				d1= wk2[i];
				for(k=0; k<n; ++k) {
					d2=a[i*lj+k];
					d1=d1-d2*x[j*lj+k];
				}
				wk1[i]=d1;
			}
		}
		ier=11;
nextit: ;
	}
	free(wk2); free(wk1);
	return ier;
}

