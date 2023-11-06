/*	Evaluating the fitted polynomial at any point using known
	coefficients of orthogonal polynomials in N dimensions.
	Should be used to evaluate the polynomial using coefficients calculated
	by POLFITN. It does not calculate the derivatives, for derivatives
	use POLEVN1 or POLEVN2

	N : (input) Number of dimensions
	NK : (input) Integer array of length N containing the degree of
		polynomial in each direction
	AX : (input) Array of length LA*(3*N+3) containing the coefficients
		alpha and beta for orthogonal polynomials in each direction
		AX[3*J][I] contains alpha and AX[3*J+1][I] contains beta
		for polynomials in Jth dimension
	LA : (input) Second dimension of array AX in the calling function.
		It must be same as what was used in call to POLFITN while
		calculating the coefficients.
	WT : (input) Array of length (MK[0]+1)(MK[1]+1)...(MK[N-1]+1)
		containing the coefficients of the fit. The dimensions of WT in
		the calling function must match the size along each dimension,
		WT[MK[N-1]+1]...[MK[1]+1][MK[0]+1]
	X0 : (input) Array of length N containing the coordinates of
		the point at which polynomial needs to be evaluated
	F : (output) Calculated value of the fitted polynomial at X0
	
	The returned value is always zero.

	Required functions : POLORT
*/

#include <math.h>
#include <stdlib.h>

int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[]);

int polevn(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f)

{
	int i,j,n1,nj,ndp,index;
	double xb,term;
	int *iwk;
	double *wk, *wk1, *wk2;

	wk=(double *) calloc((size_t) (n*la), sizeof(double));
	wk1=(double *) calloc((size_t) (n*la), sizeof(double));
	wk2=(double *) calloc((size_t) (n*la), sizeof(double));
	iwk=(int *) calloc((size_t) n, sizeof(int));

/*     Calculate the orthogonal polynomials along each dimension */
	for(i=0; i<n; ++i) {
		n1=la*i; nj=3*i*la;
		xb=x0[i];
		j=polort(nk[i],&ax[nj],&ax[nj+la],xb,&wk[n1],&wk1[n1],&wk2[n1]);
		iwk[i]=0;
	}
 
/*     Calculate the summation over n dimensions */
	*f=0.0;
	do {
		index=iwk[0];
		term=wk[index];
		ndp=nk[0]+1;
		for(i=1; i<n; ++i) {
			n1=i*la;
			term=term*wk[n1+iwk[i]];
			index=index+iwk[i]*ndp;
			ndp=ndp*(nk[i]+1);
		}
		*f=(*f)+term*wt[index];

/*     Choose the next point */
		j=0;
		while(j<=n-1) {
			if(iwk[j]>=nk[j]) {
				iwk[j]=0;
				j=j+1;
			}
			else {
				iwk[j]=iwk[j]+1;
				break;
			}
		}
	} while(j<n);
	
	free(iwk); free(wk2); free(wk1); free(wk);
	return 0;
}
