/*	Evaluating the fitted polynomial and its derivatives at any point
	using known coefficients of orthogonal polynomials in N dimensions.
	Should be used to evaluate the polynomial using coefficients calculated
	by POLFITN. If second derivative is not required use POLEVN1, if
	no derivatives are required then use POLEVN

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
	DF : (output) Array of length N containing the first derivatives
		of F at X0
	DDF : (output) Array of length N*N containing the second derivatives
		of F at X0, DDF[I][J]=d^2F/(dX[I] dX[J])
	
	The returned value is always zero.

	Required functions : POLORT
*/

#include <math.h>
#include <stdlib.h>

int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[]);

int polevn2(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f, double df[], double *ddf)

{
	int i,j,j1,n1,nj,ndp,index;
	double xb,term;
	int *iwk;
	double *wk, *wk1, *wk2, *wd, *wd2;

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
	for(i=0; i<n; ++i)  df[i]=0.0;
	for(i=0; i<n*n; ++i) ddf[i]=0.0;
	wd=(double *) calloc((size_t) n, sizeof(double));
	wd2=(double *) calloc((size_t) (n*n), sizeof(double));

	do {
		index=iwk[0];
		term=wk[index];
		wd[0]=wk1[index];
		wd2[0]=wk2[index];
		for (i=1; i<n; ++i) {wd[i]=wk[index]; wd2[i]=wk1[index];}
		for (i=n; i<n*n; ++i) wd2[i]=wk[index];

		ndp=nk[0]+1;
		for(i=1; i<n; ++i) {
			n1=i*la;
			term=term*wk[n1+iwk[i]];
			for (j=0; j<n; ++j) {
				if(j==i) wd[j]=wd[j]*wk1[n1+iwk[i]];
				else wd[j]=wd[j]*wk[n1+iwk[i]];

				for (j1=0; j1<=j; ++j1) {
					if(j==i && j1==i) wd2[j+j1*n]=wd2[j+j1*n]*wk2[n1+iwk[i]];
					else if(j==i || j1==i) wd2[j+j1*n]=wd2[j+j1*n]*wk1[n1+iwk[i]];
					else wd2[j+j1*n]=wd2[j+j1*n]*wk[n1+iwk[i]];
				}
			}

			index=index+iwk[i]*ndp;
			ndp=ndp*(nk[i]+1);
		}

		*f=(*f)+term*wt[index];
		for(i=0; i<n; ++i) {
			df[i]=df[i]+wd[i]*wt[index];
			for(j=0; j<=i; ++j) ddf[i+j*n]=ddf[i+j*n]+wd2[i+j*n]*wt[index];
		}

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
	
	for(i=0; i<n; ++i) {
		for(j=i+1; j<n; ++j) ddf[i+j*n]=ddf[j+i*n];
	}

	free(wd2); free(wd); free(iwk); free(wk2); free(wk1); free(wk);
	return 0;
}
