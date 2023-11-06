/*	Least squares polynomial fit in n dimensions using orthogonal polynomials */

#include <stdio.h>
#include <math.h>

double fun(double x, double y, double z);
double rangau(double *seed);
int polevl(int m, double a[], double alp[], double beta[], double x,
	double *f, double *df, double *ddf);
int polfit1(int n, int m, int num, double x[], double *f, double w[],
	double *a, double alp[], double beta[], double gam[]);
int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[]);
int polevn2(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f, double df[], double *ddf);
int polevn(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f);
int polevn1(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f, double df[]);
int polfitn(int n, int nk[], double *x, double *f, double *ax, int la,
	double *c, int mk[], double *fy, double *chisq);

double seed;

main()
{
	int i,i1,j,nx,ny,mx,my, id, iflg, ier,np,n, nk[20], mk[20];
	double hh, f[20][20][20], fx,  chi, a[20][20],y[200],fy[20][20][20];
	double xy[5][100],x0[9],df[9],d2f[3][3],axy[1900],c[20000];

	seed=2; np=100; id=20;
/*	No. of data points is fixed by dimension */
	n=3; nx=20;
	printf("type mk[0],mk[1],mk[2] = degree of polynomial in 3 dimension \n");
	scanf(" %d %d %d", &mk[0],&mk[1],&mk[2]);

/*	Generating the data with random error */
	hh=1.0/(nx-1.0); 
	nk[0]=nx; nk[1]=nx; nk[2]=nx;
	for(i=0; i<nx; ++i) {
		xy[0][i]=i*hh;
		xy[1][i]=i*hh;
		xy[2][i]=i*hh;
	}
	for(i=0; i<nx; ++i) {
		for(j=0; j<nx; ++j) {
			for(i1=0; i1<nx; ++i1) f[i1][i][j]=fun(xy[0][j],xy[1][i],xy[2][i1]);
		}
	}
	i=polfitn(n,nk,&xy[0][0],&f[0][0][0],axy,np,&c[0],mk,&fy[0][0][0],&chi);
	printf(" ier = %d, no. of pts = %d %d %d , degree = %d %d %d  , chisq = %e \n", i,nk[0],nk[1],nk[2],mk[0],mk[1],mk[2],chi);

	for(i1=0; i1<99; ++i1) {
		printf("type x0[0], x0[1], x0[2] = coordinates of point where fitted value is required \n");
		printf("                             (quits when x0[0]<-100)\n");
		scanf(" %le %le %le", &x0[0],&x0[1],&x0[2]);
		if(x0[0]<-100) return 0;

/*	Evaluate the fitted polynomial and derivatives using all three versions of polevn */
		i=polevn2(n,mk,axy,np,c,x0,&fx,df,&d2f[0][0]);
		printf(" x= %e %e %e   f =  %e    f' = %e  %e %e \n  f'' : ",x0[0],x0[1],x0[2],fx,df[0],df[1],df[2]);
		for(i=0; i<n; ++i) printf(" %e %e %e \n",d2f[i][0],d2f[i][1],d2f[i][2]);

		i=polevn(n,mk,axy,np,c,x0,&fx);
		printf("  f =  %e  \n",fx);

		i=polevn1(n,mk,axy,np,c,x0,&fx,df);
		printf(" f =  %e    f' = %e  %e %e\n",fx,df[0],df[1],df[2]);
	}
	return;
}



double fun(double x, double y, double z)

{
	double fx;
	fx=((231*x*y-315)*z*x+105)*z*y-5+rangau(&seed)*1.e-3;
	return fx;

}




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




/*	Evaluating the fitted polynomial and its first derivative at any point
	using known coefficients of orthogonal polynomials in N dimensions.
	Should be used to evaluate the polynomial using coefficients calculated
	by POLFITN. It does not calculate the second derivatives, for that
	use POLEVN2, for no derivatives use POLEVN

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
	
	The returned value is always zero.

	Required functions : POLORT
*/

#include <math.h>
#include <stdlib.h>

int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[]);

int polevn1(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f, double df[])

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
	wd=(double *) calloc((size_t) n, sizeof(double));

	do {
		index=iwk[0];
		term=wk[index];
		wd[0]=wk1[index];
		for (i=1; i<n; ++i) wd[i]=wk[index];

		ndp=nk[0]+1;
		for(i=1; i<n; ++i) {
			n1=i*la;
			term=term*wk[n1+iwk[i]];

			for (j=0; j<n; ++j) {
				if(j==i) wd[j]=wd[j]*wk1[n1+iwk[i]];
				else wd[j]=wd[j]*wk[n1+iwk[i]];
			}

			index=index+iwk[i]*ndp;
			ndp=ndp*(nk[i]+1);
		}

		*f=(*f)+term*wt[index];
		for(i=0; i<n; ++i)  df[i]=df[i]+wd[i]*wt[index];

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
	
	free(wd); free(iwk); free(wk2); free(wk1); free(wk);
	return 0;
}




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




/*	Least squares polynomial fit using orthogonal polynomials in 1 dimension
	Modified version of POLFIT to fit multiple sets of function values
	This function is meant to be used for fit in multiple dimensions

	N : (input) Number of data points
	M : (input) Required degree of polynomial
	NUM : (input) Number of different RHS (function values) to be fitted
		Each set must be defined over the same abscissas and
		with same weights.
	X : (input) Array of length N containing the abscissas
	F : (input) Array of length N*NUM containing the function values
		F[J][I] is the value at X[I] in Jth data set
		The second dimension of F is assumed to be exactly equal
		to N to minimise storage requirement.
	SIG : (input) Array of length N containing the errors associated
		with each point. Errors are same for all data sets.
		SIG[I] is the error for F[J][I].
	A : (output) Array of length (M+1)*NUM containing the coefficients
		for the fit for each RHS. The second dimension of A
		is assumed to be M+1 to minimise storage requirements
		A[j][i] is the ith coefficient for jth RHS.
	ALP, BETA : (output) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
	GAM : (output) Array of length M+1, containing the quantities
		\gamma_i for the orthogonal polynomials
		
	Error status is returned by the value of the function POLFIT1.
		0 value implies successful execution
		601 implies that N<M+1 or M<0 or N < 1
		621 implies that GAM[I] vanishes at some I and calculations
			are abandoned
	
	The fitted polynomial can be calculated at any value of x using POLEVL

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int polfit1(int n, int m, int num, double x[], double *f, double sig[],
	double *a, double alp[], double beta[], double gam[])

{
	int i,j,i1,i2,it,j1,m1;
	double gam0,s;
	double *wk;

	if(m>=n || m<0 || n<1) return 601;

/*	Initialisation */
	i1=0; i2=1; m1=m+1;
	gam0=0.0;
	wk=(double *) calloc((size_t) (2*n), sizeof(double));
	for(i=0; i<n; ++i) {
		wk[i]=0.0; wk[i+n]=1.0;
		gam0=gam0+1/(sig[i]*sig[i]);
	}
	gam[0]=gam0;
	beta[0]=0.0;

/*	Loop over the degree of polynomial */
	for(j=0; j<=m; ++j) {
		if(gam[j] <=0.0) {free(wk); return 621;}

		for(j1=0; j1<num; ++j1) {
			s=0.0;
			for(i=0; i<n; ++i) s=s+f[i+j1*n]*wk[i+i2*n]/(sig[i]*sig[i]);

/*	The coefficient a_j */
			a[j+j1*m1]=s/gam[j];
		}
		if(j == m) {free(wk); return 0;}

		s=0.0;
		for(i=0; i<n; ++i) s=s+x[i]*wk[i+i2*n]*wk[i+i2*n]/(sig[i]*sig[i]);
/*	The coefficient \alpha_{j+1} */
		alp[j]=s/gam[j];
		gam0=0.0;
		for(i=0; i<n; ++i) {
			wk[i+i1*n]=(x[i]-alp[j])*wk[i+i2*n]-beta[j]*wk[i+i1*n];
			gam0=gam0+wk[i+i1*n]*wk[i+i1*n]/(sig[i]*sig[i]);
		}
/*	The coefficient \beta_{j+1} */
		beta[j+1]=gam0/gam[j];
/*	The coefficient \gamma_{j+1} */
		gam[j+1]=gam0;
/*	Interchange indices I1, I2 so that only last two columns of WK are stored */
		it=i1; i1=i2; i2=it;
	}
	free(wk);
	return 0;
}




/*	Least squares polynomial fit using orthogonal polynomials in n dimensions
	Weights are assumed to be equal for all points and points are
	assumed to be on a hyper-rectangular mesh

	N : (input) Number of dimensions
	NK : (input) Integer array of length N containing the number of
		data points along each direction
	X : (input) Array of length LA*N containing the coordinates
		of tabular points, X[j][i] contains the ith point along
		jth dimension
	F : (input) Array of length NK[0]*NK[1]*...*NK[N-1] containing
		the function values. The dimension of F in the calling function
		must match the size along each dimension, F[NK[N-1]]...[NK[0]]
	AX : (output) Array of length LA*(3*N+3) containing information about
		fit along each direction. AX[3*J][I], AX[3*J+1][I], AX[3*J+2][I] will
		respectively contain the coefficients, alpha, beta, gamma
		for fit along Jth dimension.
		The rest of the array is used as scratch space
	LA : (input) Second dimension of arrays X and AX as declared
		in the calling function. LA >= MAX(NK[I])
	C : (output) Array of length (MK[0]+1)(MK[1]+1)...(MK[N-1]+1) containing
		the fitted coefficients of product of orthogonal polynomials 
		The dimension of F in the calling function must match the size
		along each dimension, C[MK[N-1]+1]...[MK[1]+1][MK[0]+1]
	MK : (input) Integer array of length N containing the required
		degree of polynomial in each dimension
	FY : (output) Array of same size and shape as F containing the
		fitted function values at each tabular point	
	CHISQ : (output) the Chi square value for the fit
		
	Error status is returned by the value of the function POLFITN.
		0 value implies successful execution
		605 implies that LA < MAX(NK[I]) 
		In  this case no calculations are done
		other values may be set by POLFIT1
	
	The fitted polynomial can be calculated at any value of x using POLEVN

	Required functions : POLFIT1, POLEVN, POLORT
*/

#include <math.h>
#include <stdlib.h>

int polfit1(int n, int m, int num, double x[], double *f, double w[],
	double *a, double alp[], double beta[], double gam[]);
int polevn(int n, int nk[], double *ax, int la, double *wt, double x0[],
	double *f);

int polfitn(int n, int nk[], double *x, double *f, double *ax, int la,
	double *c, int mk[], double *fy, double *chisq)

{
	int i,j,j1,i1,i2,lj,nmax,nx,ny,num,m,n1,ni,nx1,ny1,lj1,ij,id,ier;
	double r1;
	int *iwk;
	double *wk,*wk1;

	n1=3*n*la;
/*     Set the weights to one for fits along each dimension */
	for(i=0; i<la; ++i) ax[i+n1]=1.0;

	nx=1; nmax=nk[n-1];
	for(i=0; i<n-1; ++i) {
		nx=nx*nk[i];
		if(nk[i]>nmax) nmax=nk[i];
	}
	if(nmax>la) return 605;

	ny=1;
	wk=(double *) calloc((size_t) (nx*nk[n-1]),sizeof(double));
	wk1=(double *) calloc((size_t) (nx*nk[n-1]),sizeof(double));
 
/*     Set up the RHS for fit along Nth dimension */
	lj=nk[n-1];
	for(i=0; i<nx; ++i) {
		for(j=0; j<nk[n-1]; ++j) wk[j+i*lj]=f[i+j*nx];
	}

/*     Loop for fit along each dimension */
	for(j1=n-1; j1>=0; --j1) {
		num=nx*ny;
		lj=nk[j1];
		m=mk[j1]+1;
		ni=j1*3*la;
		ier=polfit1(lj,mk[j1],num,&x[j1*la],wk,&ax[n1],wk1,&ax[ni],&ax[ni+la],&ax[ni+2*la]);
		if(ier>0) {free(wk1); free(wk); return ier;}

		if(j1>0) {
/*     Set up the RHS for fit along next dimension */
			nx1=nx/nk[j1-1];
			ny1=ny*m;
			lj1=nk[j1-1];
			for(i1=0; i1<ny; ++i1) {
				for(i2=0; i2<m; ++i2) {
					for(i=0; i<nk[j1-1]; ++i) {
						for(j=0; j<nx1; ++j) wk[i+j*lj1+i2*lj1*nx1+i1*nx1*lj1*m] =
							wk1[i2+j*m+i*nx1*m+i1*nx*m];
					}
				}
			}

			nx=nx1; ny=ny1;
		}

		else {
/*     Store the fitted coefficients in array C */
			for(i=0; i<m; ++i) {
				for(j=0; j<ny; ++j) c[i+j*m]=wk1[i+j*m];
			}
		}
	}

/*     Calculate the Chi Square */
	*chisq=0.0;
	free(wk1);
	iwk=(int *) calloc((size_t) n, sizeof(int));

	for(i=0; i<n; ++i) iwk[i]=0;
 
/*     Loop over all points */
	do {
		ij=iwk[0];
		wk[0]=x[ij];
		id=nk[0];
		for(i=1; i<n; ++i) {
			ij=ij+id*iwk[i];
			id=id*nk[i];
			wk[i]=x[iwk[i]+i*la];
		}

		ier=polevn(n,mk,ax,la,c,wk,&fy[ij]);
		r1=f[ij]-fy[ij];
		*chisq=(*chisq)+r1*r1;

/*     Choose the next point */
		j=0;
		while(j<=n-1) {
			if(iwk[j]>=nk[j]-1) {
				iwk[j]=0;
				j=j+1;
			}
			else {
				iwk[j]=iwk[j]+1;
				break;
			}
		}
	} while(j<n);

	free(iwk); free(wk);
	return 0;
}




/*	Evaluating the orthogonal polynomial basis functions at any value
	of x using known coefficients
	Should be used to evaluate the basis using coefficients calculated
	by POLFIT or POLFIT1.

	M : (input) Degree of polynomial
	ALP, BETA : (input) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
		ALP, BETA could be calculated using POLFIT
	X : (input) Value of x at which polynomials needs to be evaluated
	F : (output) Array of length M+1 containing the value of
		orthogonal polynomials at X
	DF : (output) Array of length M+1 containing first derivative of F at X
	DDF : (output) Array of length M+1 containing second derivative of F at X

	Returned value POLORT is always zero.
	
	Required functions : None
*/	

#include <math.h>

int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[])

{
	int j;
	
	f[0]=1.0; f[1]=x-alp[0];
	df[0]=0.0; df[1]=1.0;
	ddf[0]=0.0; ddf[1]=0.0;

 
/*	The recurrence relations */
	for(j=2; j<=m; ++j) {
		ddf[j]=2.*df[j-1]+(x-alp[j-1])*ddf[j-1]-beta[j-1]*ddf[j-2];
		df[j]=f[j-1]+(x-alp[j-1])*df[j-1]-beta[j-1]*df[j-2];
		f[j]=(x-alp[j-1])*f[j-1]-beta[j-1]*f[j-2];
	}
	return 0;
}



/*	To generate random numbers with Gaussian probability distribution
	It generates random numbers with zero mean and variance of 1.
	
	SEED : (input/output) real seed, it should be positive and
		less than AN. It is updated by the function and should
		not be modified between two calls, unless a fresh
		sequence is required

	Required functions : None
	
	THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
	VERSION AS THE SEED IS NOW DOUBLE INSTEAD OF INT.

*/

#include <math.h>


double rangau(double *seed)

{
	int n2;
	double am=2147483648.0, a=45875.0, ac=453816693.0, an=2147483647.0, r1, rn;

	rn=a*(*seed)+ac; n2=rn/am; r1=rn-am*n2;
	if(*seed==0.0) *seed=0.1;
	rn=sqrt(2.0*log(an/(*seed)))*cos(2.0*M_PI*r1/an);
	*seed=r1;
	return rn;
}

