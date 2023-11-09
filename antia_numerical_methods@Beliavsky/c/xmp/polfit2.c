/*	Least squares polynomial fit in 2 dimensions using orthogonal polynomials */

#include <stdio.h>
#include <math.h>

double fun(double x, double y);
double rangau(double *seed);
int polfit1(int n, int m, int num, double x[], double *f, double w[],
	double *a, double alp[], double beta[], double gam[]);
int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[]);
int polev2(int nx, int ny, double *ax, double *ay, int la, double *wt,
	double x0, double y0, double *f, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy);
int polfit2(int nx, int ny, double x[], double y[], double *f, double *ax,
	double *ay, int la, double *c, int ic, int mx, int my, double *fy,
	double *chisq);

double seed;

main()
{
	int i,i1,j,nx,ny,mx,my, id, iflg, ier,np,nmax;
	double hh, x[100], f[100][100], w[100], ax[900],ay[200],gam[20],
	dfxx,dfyy,dfxy,dfx,dfy, fx, xx,yy, chi, a[20][20],y[200],fy[100][100];

	seed=2; np=100; id=20;
	printf("type nx=no. of data pts along each axis,   mx,my=degree of polynomial in x,y\n");
	scanf(" %d %d %d", &nx,&mx,&my);

/*	Generating the data with a known function plus random error */

	hh=1.0/(nx-1.0); ny=nx;
	for(i=0; i<nx; ++i) {
		x[i]=i*hh;
		y[i]=i*hh;
	}

	for(i=0; i<ny; ++i) {
		for(j=0; j<nx; ++j) f[i][j]=fun(x[j],y[i]);
	}

	i=polfit2(nx,ny,x,y,&f[0][0],ax,ay,np,&a[0][0],id,mx,my,&fy[0][0],&chi);
	printf(" ier = %d   no. of pts = %d %d   degree = %d %d   chisq = %e \n", i,nx,ny,mx,my,chi);

	for(i1=0; i1<99; ++i1) {
		printf("type x,y values where fitted value is to be evaluated \n");
		printf("                          (quits when xx<-100)\n");
		scanf(" %le %le", &xx,&yy);
		if(xx<-100) return 0;

		i=polev2(mx,my,ax,ay,id,&a[0][0],xx,yy,&fx,&dfx,&dfy,
				&dfxx,&dfxy,&dfyy);
		printf(" x= %e, %e    f =  %e    f' = %e  %e \n",xx,yy,fx,dfx,dfy);
		printf(" f'' = %e %e %e \n",dfxx,dfxy,dfyy);

	}
	return;
}

double fun(double x, double y)

{
	double fx;
	fx=((231*x*y-315)*x*x+105)*x*y-5+rangau(&seed)*5.e-3;
	return fx;

}




/*	Evaluating the fitted polynomial and its derivatives at any value
	of x using known coefficients of orthogonal polynomials in 2 dimensions
	Should be used to evaluate the polynomial using coefficients calculated
	by POLFIT2.

	NX : (input) Degree of polynomial in X
	NY : (input) Degree of polynomial in Y
	AX : (input) Array of length LA*2 containing the coefficients
		alpha and beta for orthogonal polynomials in X
		AX[0][I] contains alpha and AX[1][I] contains beta
	AY : (input) Array of length LA*2 containing the coefficients
		alpha and beta for orthogonal polynomials in Y
		AY[0][I] contains alpha and AY[1][I] contains beta
		The arrays AX and AY can be calculated using POLFIT2
	LA : (input) Second dimension of arrays AX, AY and WT in the calling
		function. LA > MAX(NX,NY)
	WT : (input) Array of length LA*(MY+1) containing the coefficients
		of the fit. WT[J][I] is the coefficient of
		PHI_I(X)PSI_J(Y), where PHI_I and PSI_J are orthogonal
		polynomials in X and Y
	X0,Y0 : (input) Coordinates of the point at which polynomial
		needs to be evaluated
	F : (output) Calculated value of the fitted polynomial at (X0,Y0)
	DFX : (output) First derivative  dF/dX at X0,Y0
	DFY : (output) First derivative dF/dY at X0,Y0
	DFXX : (output) Second derivative d^2F/dXdX at X0,Y0
	DFXY : (output) Second derivative d^2F/dXdY at X0,Y0
	DFYY : (output) Second derivative d^2F/dYdY at X0,Y0
		
	Error status is returned by the value of the function POLEV2.
		0 value implies successful execution
		604 implies that LA<=MAX(NX,NY), in which case no
			calculations are done.
	
	Required functions : POLORT
*/

#include <math.h>
#include <stdlib.h>

int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[]);

int polev2(int nx, int ny, double *ax, double *ay, int la, double *wt,
	double x0, double y0, double *f, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy)

{
	int i,j,nk;
	double *wkx,*wky;

	nk=nx+1; if(ny>nx) nk=ny+1;
	if(nk-1>la) return 604;

	wkx=(double *) calloc((size_t) (3*nk),sizeof(double));
	wky=(double *) calloc((size_t) (3*nk),sizeof(double));

/*	Calculate the orthogonal polynomials along each dimension */
	i=polort(nx,ax,&ax[la],x0,wkx,&wkx[nk],&wkx[2*nk]);
	i=polort(ny,ay,&ay[la],y0,wky,&wky[nk],&wky[2*nk]);
 
/*	Calculate the fitted polynomial and its derivatives */
	*f=0.0; *dfx=0.0; *dfy=0.0;
	*dfxx=0.0; *dfxy=0.0; *dfyy=0.0;
	for(i=0; i<=nx; ++i) {
		for(j=0; j<=ny; ++j) {
			*f=(*f)+wt[i+j*la]*wkx[i]*wky[j];
			*dfx=(*dfx)+wt[i+j*la]*wkx[nk+i]*wky[j];
			*dfy=(*dfy)+wt[i+j*la]*wkx[i]*wky[nk+j];
			*dfxx=(*dfxx)+wt[i+j*la]*wkx[2*nk+i]*wky[j];
			*dfyy=(*dfyy)+wt[i+j*la]*wkx[i]*wky[2*nk+j];
			*dfxy=(*dfxy)+wt[i+j*la]*wkx[nk+i]*wky[nk+j];
		}
	}
	free(wky); free(wkx);
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



/*	Least squares polynomial fit using orthogonal polynomials in 2 dimensions
	Weights are assumed to be equal for all points

	NX : (input) Number of data points along X direction
	NY : (input) Number of data points along Y direction
	X,Y : (input) Arrays of length NX,NY containing the coordinates
		of tabular points
	F : (input) Array of length LA*NY containing the function values
		F[j][i] is the value at X[I],Y[J]
	AX : (output) Array of length IC*3 containing information about
		fit along X direction. AX[0][I], AX[1][I], AX[2][I] will
		respectively contain the coefficients, alpha, beta, gamma
	AY : (output) Array of length IC*3 containing information about
		fit along Y direction. AY[0][I], AY[1][I], AY[2][I] will
		respectively contain the coefficients, alpha, beta, gamma
	LA : (input) Second dimension of arrays  F, FY as declared
		in the calling function. LA >= MAX(NX,NY)
	C : (output) Array of length IC*(MY+1) containing the fitted
		coefficients of product of orthogonal polynomials in X & Y
	IC : (input) Second dimension of arrays C, AX, AY as declared in the calling
		function. IC > MAX(MX,MY)+1
	MX : (input) Required degree of polynomial in X
	MY : (input) Required degree of polynomial in Y
	FY : (output) Array of length LA*NY containing the values of
		fitted function at each tabular point	
	CHISQ : (output) the Chi square value for the fit
		
	Error status is returned by the value of the function POLFIT2.
		0 value implies successful execution
		602 implies that IC<MX+1 or IC<MY+1 or LA<NX or LA<NY
		603 implies that NX<MX+1 or MX<0 or NY < MY+1 or MY<0
		In both these cases no calculations are done
		Other values may be set by POLFIT1
	
	The fitted polynomial can be calculated at any value of x using POLEV2

	Required functions : POLFIT1, POLEV2, POLORT
*/

#include <math.h>
#include <stdlib.h>

int polfit1(int n, int m, int num, double x[], double *f, double w[],
	double *a, double alp[], double beta[], double gam[]);
int polev2(int nx, int ny, double *ax, double *ay, int la, double *wt,
	double x0, double y0, double *f, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy);

int polfit2(int nx, int ny, double x[], double y[], double *f, double *ax,
	double *ay, int la, double *c, int ic, int mx, int my, double *fy,
	double *chisq)

{
	int i,j,lj,ln1,m,num,ier;
	double dfx,dfy,dfxx,dfxy,dfyy,r1;
	double *wk,*wk1, *aw;

	if(nx>la || ny>la || ic<mx+1 || ic<my+1) return 602;

	if(mx<0 || my<0 || nx<=mx || ny<=my) return 603;

/*     Set the weights to 1 */
	lj=nx; if(lj<ny) lj=ny;
	aw=(double *) calloc((size_t) lj, sizeof(double));
	for(i=0; i<lj; ++i) aw[i]=1.0;

/*     Set up the RHS for calculating the fits along y-axis */
	wk=(double *) calloc((size_t) (nx*ny), sizeof(double));
	wk1=(double *) calloc((size_t) (nx*(my+1)), sizeof(double));
	lj=ny;
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) wk[j+i*lj]=f[i+j*la];
	}

	m=my+1;
	num=nx;
	ier=polfit1(ny,my,num,y,wk,aw,wk1,ay,&ay[ic],&ay[2*ic]);
	if(ier>100) { free(wk1); free(wk); free(aw); return ier;}

/*     Set up the RHS for calculating the fits along x-axis */
	for(j=0; j<m; ++j) {
		for(i=0; i<nx; ++i) wk[i+j*nx]=wk1[j+i*m];
	}
	num=m;
	ier=polfit1(nx,mx,num,x,wk,aw,wk1,ax,&ax[ic],&ax[2*ic]);
	if(ier>100) {free(wk1); free(wk); free(aw); return ier;}

/*	Store the calculated coefficients in array C */
	m=mx+1;
	for(i=0; i<=my; ++i) {
		for(j=0; j<m; ++j) c[j+ic*i]=wk1[j+i*m];
	}
 
/*     Calculate the CHI square */
	*chisq=0.0;
	for(i=0; i<ny; ++i) {
		for(j=0; j<nx; ++j) {
			ier=polev2(mx,my,ax,ay,ic,c,x[j],y[i],&fy[j+i*la],
				&dfx,&dfy,&dfxx,&dfxy,&dfyy);
			r1=f[j+i*la]-fy[j+i*la];
			*chisq=(*chisq)+r1*r1;
		}
	}
	free(wk1); free(wk); free(aw);
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

