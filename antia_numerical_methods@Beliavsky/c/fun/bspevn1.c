/*	To calculate function value and first derivative using B-spline
		expansion in n-dimensions

	N : (input) Number of dimensions
	NK : (input) Array of length N containing the number of knots
		along each dimension
	X : (input) Array of length NXD*N containing the knots
		X[j][i] is the ith knot along jth dimension
	NXD : (input) Second dimension of array X in calling function.
		NXD must be greater than or equal to maximum of NK[I].
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Coefficients of B-spline expansion, the coefficients
		are assumed to be stored in natural FORTRAN order with
		no gaps in data. The calling function should have
	 	dimension  WT[NK[N-1]+K-2]...[NK[1]+K-2][NK[0]+K-2]
	X0 : (input) Array of length N, containing the coordinates
		of the point at which function needs to be evaluated
	DF : (output) Array of length N which will contain the
		calculated values of first derivatives.
		DF[i] will be derivative w.r.t. X_i
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values of IER may be set by BSPLIN which is called

	BSPEVN1 will give the value of the B-spline expansion at X0.
	This function does not calculate the second derivative, for that
	use BSPEVN2. To improve efficiency if derivatives are not
	required then use BSPEVN.

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

double bspevn1(int n, int nk[], double *x, int nxd, int k, double *wt,
	double x0[], double df[], int *ier)

{
	int i,j,j1,n1,n2,ndp,index,nderiv;
	double f,term;
	int *iwk;
	double *wk, *wk1, *wk2, *wd1;
 
	nderiv=1;
	iwk=(int *) calloc((size_t) (2*n),sizeof(int));
	wk=(double *) calloc((size_t) (n*nxd), sizeof(double));
	wk1=(double *) calloc((size_t) (n*nxd), sizeof(double));
	wk2=(double *) calloc((size_t) (n*nxd), sizeof(double));

	for(i=0; i<n; ++i) {
		n1=nxd*i;
		*ier=bsplin((x+i*nxd),nk[i],k,x0[i],nderiv,&wk[n1],&wk1[n1],
			&wk2[n1],&iwk[i]);
		if(*ier>100) {free(wk2); free(wk1); free(wk); free(iwk); return 0.0;}
		iwk[i+n]=0;
	}

/*	Calculate the sum over all points */
	f=0.0;
	for(i=0; i<n; ++i) df[i]=0.0;
	wd1=(double *) calloc((size_t) n, sizeof(double));

	do {
		index=iwk[0]+iwk[n];
		term=wk[index];
		wd1[0]=wk1[index];
		for(i=1; i<n; ++i) wd1[i]=wk[index];
		ndp=nk[0]+k-2;

		for(i=1; i<n; ++i) {
			n1=i*nxd;
			n2=iwk[i]+iwk[i+n];
			term=term*wk[n1+n2];

			for(j=0; j<n; ++j) {
				if(j == i) wd1[j]=wd1[j]*wk1[n1+n2];
				else wd1[j]=wd1[j]*wk[n1+n2];
			}
			index=index+n2*ndp;
			ndp=ndp*(nk[i]+k-2);
		}
		f=f+term*wt[index];
		for(i=0; i<n; ++i) df[i]=df[i]+wd1[i]*wt[index];

		j=0;
		while(j <= n-1) {
			if(iwk[j+n] >= k-1) {
				iwk[j+n]=0;
				j=j+1;
			}
			else {
				iwk[j+n]=iwk[j+n]+1;
				break;
			}
		}
	} while(j<n);

	free(wd1); free(wk2); free(wk1); free(wk); free(iwk);
	return f;
}

