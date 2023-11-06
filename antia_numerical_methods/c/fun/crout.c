/*	Solution of a system of linear equations using Crout's algorithm 
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
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
			interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter which determines the type of computation
		required.
		If IFLG<=0, both elimination and solution are calculated
			and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function CROUT
		0 value implies successful execution
		102 implies (N<0 or N>LJ) 
		122 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg)

{
	int i,j,k,km,l;
	double r1,t1,r2,d1;
	double *wk;

	if(n<=0 || n>lj) return 102;

	if((*iflg)<2)
/*	Perform LU decomposition  */
	{wk=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i){
		r1=0.0;
		for(j=0; j<n; ++j){
			if(fabs(a[i*lj+j]) > r1) r1=fabs(a[i*lj+j]);}
		wk[i]=r1;
		if(r1==0.0) {free(wk); return 122;}
	}

	*det=1.0; *idet=0;
	for(k=0; k<n; ++k){
		r1=0.0; km=k;
/*	Generate the Kth column of L  */
		for (l=k; l<n; ++l){
			d1=a[l*lj+k];
			for(i=0; i<=k-1; ++i) d1=d1-a[l*lj+i]*a[i*lj+k];
			a[l*lj+k]=d1;
/*	Finding the pivot  */
			r2=fabs(d1/wk[l]);
			if(r2>r1) {r1=r2; km=l;}
		}

		inc[k]=km;
/*	Interchange the rows if needed  */
		if(km != k) {
			*det=-(*det);
			for(l=0; l<n; ++l) {
				t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;
			}
			t1=wk[k]; wk[k]=wk[km]; wk[km]=t1;
		}

		*det=*det*a[k*lj+k];
		if(a[k*lj+k]==0.0) {free(wk); return 122;}

		if(*det != 0.0) {
/*	Scale the value of the determinant DET */
			while(fabs(*det)>32.0) {
				*det=(*det)*0.03125; *idet=*idet+5;}

			while(fabs(*det)<0.03125) {
				*det=(*det)*32.0; *idet=*idet-5;}
		}

/*	Generate the Kth row of U */
		for(l=k+1; l<n; ++l) {
			d1=a[k*lj+l];
			for(i=0; i<=k-1; ++i) d1=d1-a[k*lj+i]*a[i*lj+l];
			a[k*lj+l]=d1/a[k*lj+k];
		}
	}
	free(wk);

	if(*iflg==1) {*iflg=2; return 0;}
	*iflg=2;
	}

/*	Solution for NUM different right-hand sides  */
	for(j=0; j<num; ++j){
/*	Forward substitution  */
		for(k=0; k<n; ++k) {
			if(k != inc[k]) {
				t1=x[j*lj+k];
				x[j*lj+k]=x[j*lj+inc[k]];
				x[j*lj+inc[k]]=t1;
			}

			d1=x[j*lj+k];
			for(l=0; l<=k-1; ++l) d1=d1-a[k*lj+l]*x[j*lj+l];
			x[j*lj+k]=d1/a[k*lj+k];
		}

/*	Back-substitution  */
		for(k=n-2; k>=0; --k){
			d1=x[j*lj+k];
			for(l=n-1; l>=k+1; --l) d1=d1-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=d1;
		}
	}
	return 0;
}
