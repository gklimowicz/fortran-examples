/*	Least squares polynomial fit using orthogonal polynomials */

#include <stdio.h>
#include <math.h>

double fun(double x);
double rangau(double *seed);
int polevl(int m, double a[], double alp[], double beta[], double x,
	double *f, double *df, double *ddf);
int polfit(int n, int m, double x[], double f[], double w[], double a[],
	double alp[], double beta[], double y[], double h[], double gam[]);

double seed;

main()
{
	int i,i1,j,n,m, id, iflg, ier,np,nmax;
	double hh, x[100], f[100], w[100], alp[20],beta[20],gam[20],
	d2f,df, fx, xx, a[20],y[200],h[20];

/*	Example 10.3 */

	seed=2;
	printf("type  n = no. of data points,   m = degree of polynomial to be fitted\n");
	scanf(" %d %d", &n,&m);

/*	Generate artificial data using a known function with random errors added */
	hh=1.0/(n-1.0);
	for(i=0; i<n; ++i) {
		x[i]=i*hh;
		f[i]=fun(x[i]);
		w[i]=5.0e-3;
	}

	i=polfit(n,m,x,f,w,a,alp,beta,y,h,gam);
	printf(" ier = %d  no. of pts = %d , degree of polynomial = %d \n", i,n,m);
	printf("Degree  Coefficient    Stnd. dev     Alpha         Beta        Chi square \n");
	for(i=0; i<=m; ++i) printf(" %d    %e   %e  %e  %e  %e \n",i,a[i],1.0/sqrt(gam[i]),alp[i],beta[i],h[i]);

	for(i1=0; i1<99; ++i1) {
		printf("type  xx = point at which polynomial is to evaluated    (quits when xx<-100)\n");
		scanf(" %le", &xx);
		if(xx<-100) return 0;

		i=polevl(m,a,alp,beta,xx,&fx,&df,&d2f);
		printf(" x = %e    f(x) =  %e    f'(x) = %e    f''(x) = %e \n",xx,fx,df,d2f);
	}
	return;
}

double fun(double x)

{
	double fx;
	fx=((231*x*x-315)*x*x+105)*x*x-5+rangau(&seed)*5.e-3;
	return fx;

}




/*	Evaluating the fitted polynomial and its derivatives at any value
	of x using known coefficients of orthogonal polynomials
	Should be used to evaluate the polynomial using coefficients calculated
	by POLFIT.

	M : (input) Degree of polynomial
	A : (input) Array of length M+1 containing the coefficients
		of the fit
	ALP, BETA : (input) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
		A, ALP, BETA could be calculated using POLFIT
	X : (input) Value of x at which polynomial needs to be evaluated
	F : (output) Calculated value of the polynomial at X
	DF : (output) First derivative of F(X) at X
	DDF : (output) Second derivative of F(X) at X

	The returned value of POLEVL is always zero
	
	Required functions : None
*/

#include <math.h>

int polevl(int m, double a[], double alp[], double beta[], double x,
	double *f, double *df, double *ddf)

{
	int j;
	double f1,df1,ddf1,ff,d,dd;

	*f=a[m-1]+(x-alp[m-1])*a[m];
	f1=a[m];
	*df=a[m]; df1=0.0;
	*ddf=0.0; ddf1=0.0;
	if(m<=1) return 0;

/*	Clenshaw's recurrence for F, DF and DDF */
	for(j=m-2; j>=0; --j) {
		dd=2.*(*df)+(x-alp[j])*(*ddf)-beta[j+1]*ddf1;
		d=(*f)+(x-alp[j])*(*df)-beta[j+1]*df1;
		ff=a[j]+(x-alp[j])*(*f)-beta[j+1]*f1;
		f1=(*f); *f=ff;
		df1=(*df); (*df)=d;
		ddf1=(*ddf); *ddf=dd;
	}
	return 0;
}



/*	Least squares polynomial fit using orthogonal polynomials in 1 dimension

	N : (input) Number of data points
	M : (input) Required degree of polynomial
	X : (input) Array of length N containing the abscissas
	F : (input) Array of length N containing the function values
		F[I] is the value at X[I]
        SIG : (input) Real array of length N containing the estimated
	        error in the function values, SIG(I) is the error in F(I)
	A : (output) Array of length M+1 containing the coefficients for the fit
	ALP, BETA : (output) Arrays of length M+1, containing the coefficients
		required for defining the orthogonal polynomials
	Y : (output) Array of length N containing the fitted values at X[I]
	H : (output) Array of length M+1 containing the \chi^2 values
		for residuals using polynomial fit of degrees 0,...,M
	GAM : (output) Array of length M+1, containing the quantities
		\gamma_i for the orthogonal polynomials
		
	Error status is returned by the value of the function POLFIT.
		0 value implies successful execution
		601 implies that N<M+1 or M<0 or N < 1
		621 implies that GAM[I] vanishes at some I and calculations
			are abandoned
	
	The fitted polynomial can be calculated at any value of x using POLEVL

	Required functions : None

	THE ARGUMENTS OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
	VERSION AS THE 5TH ARGUMENT IS NOW ERROR INSTEAD OF WEIGHT

*/

#include <math.h>
#include <stdlib.h>

int polfit(int n, int m, double x[], double f[], double sig[], double a[],
	double alp[], double beta[], double y[], double h[], double gam[])

{
	int i,j,i1,i2,it;
	double gam0,h0,s;
	double *wk;

	if(m>=n || m<0 || n<1) return 601;

/*	Initialisation */
	i1=0; i2=1;
	gam0=0.0; h0=0.0;
	wk=(double *) calloc((size_t) (2*n), sizeof(double));
	for(i=0; i<n; ++i) {
		wk[i]=0.0; wk[i+n]=1.0;
		y[i]=0.0;
		gam0=gam0+1/(sig[i]*sig[i]);
		h0=h0+f[i]*f[i]/(sig[i]*sig[i]);
	}
	gam[0]=gam0;
	beta[0]=0.0;

/*	Loop over the degree of polynomial */
	for(j=0; j<=m; ++j) {
		s=0.0;
		for(i=0; i<n; ++i) s=s+f[i]*wk[i+i2*n]/(sig[i]*sig[i]);
		if(gam[j] <=0.0) {free(wk); return 621;}

/*	The coefficient a_j */
		a[j]=s/gam[j];
		h0=h0-a[j]*a[j]*gam[j];
		h[j]=h0;
		for(i=0; i<n; ++i) y[i]=y[i]+a[j]*wk[i+i2*n];
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

