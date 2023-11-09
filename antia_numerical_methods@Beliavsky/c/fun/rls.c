/*	To solve a linear inverse problem in one dimension using RLS
	technique with B-spline basis functions

	NK : (input) Number of knots for defining B-splines, the number
		of basis functions would be NK+K-2
	XO : (input) Array of length NK containing the knots
		used for defining B-spline basis functions.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines required, K=4 gives cubic B-splines
	NR : (input) Number of points used in defining the kernels
	R : (input) Array of length NR containing the coordinates
		of points at which kernels are available.
	RKER : (input) Array of length IK*NR containing the kernels
		for the inverse problem. RKER[J][I] should contain the
		value at R[J] for the Ith kernel. This array
		must be supplied if IFLG<2, otherwise it is not required
	IK : (input) Second dimension of arrays RKER, AC and A, as specified
		in the calling function. IK>=NM+NS
	AC : (input/output) Array of length IK*(NK+K-2) containing
		the coefficients of matrix defining the inverse problem
		If IFLG<2, these coefficients are calculating by integrating
		the kernels with appropriate weights. For IFLG=2,3 these
		coefficients must be supplied.
	NM : (input) Number of data points in the inverse problem
	NS : (input) Number of points to be used for applying regularisation
		The function chooses a uniform mesh covering the full interval
		for applying smoothing.
	ALP : (input) Regularisation parameter, ALP>0.
	IDE : (input) Order of derivative to be used for regularisation,
		IDE should be 1 or 2 for first or second derivative smoothing
	DI : (input) Array of length NM, containing the data points for inversion
	DE : (input) Array of length NM, containing the estimated error in DI. 
	DF : (output) Array of length NM, containing the normalised
		residuals (DI-DI(fit))/DE for each data point.
	F : (output) Array of length NR which will contain the
		calculated solution at each point in array R.
	B : (output) Array of length NM+NS containing the coefficients
		of basis functions in fitted solution. Although the
		the number of coefficients is only NK+K-2, the rest
		of array is used as scratch space
	IFLG : (input/output) Integer specifying the type of calculation required.
		IFLG=0 : The matrix coefficients will be calculated using
			the kernels and then the equations are solved to
			find the solution for given data points DI.
			IFLG is set to 4 after calculations.
		IFLG=1 : The matrix coefficients will be calculated using
			the kernels and then the SVD of the full matrix is
			computed, but the solution is not computed.
			IFLG is set to 4 after calculations.
		IFLG=2 : The matrix coefficients are assumed to be available
			in array AC and the matrix is setup and solved
			to find the solution for given data points DI.
			IFLG is set to 4 after calculations.
		IFLG=3 : The matrix coefficients are assumed to be available
			in array AC and the matrix is setup and the SVD
			is computed, but the solution is not computed.
			IFLG is set to 4 after calculations.
		IFLG=4 : The SVD of matrix is assumed to be available
			from previous calculations and only the solution
			for given DI is computed.
		Since IFLG is set to 4 every-time, it should be reset to
		0 or 2 before next call when the data or error estimates
		or smoothing are changed.
	REPS : (input) Required accuracy for solution of equations using
		SVD. Singular values less than REPS times maximum will be
		set to zero.
	CHISQ : (output) The computed value of Chi square for the solution
	SUMD : (output) The computed value of the smoothing term
	A : (input/output) Array of length IV*(NM+NS) containing
		the SVD of the matrix of equations. If IFLG<4 this matrix
		will be calculated, otherwise it must be supplied.
	AV : (input/output) Array of length IV*(NK+K-2) containing
		the matrix V or SVD of the matrix of equations. if IFLG<4
		this matrix will be calculated, otherwise it must be supplied.
	IV : (input) The second dimension of A, AV as declared in the calling
		function. IV>=NK+K-2
	SIGMA : (input/output) Array of length NK+K-2 containing the
		singular values of the matrix A. If IFLG<4 this array will
		be calculated, otherwise it must be supplied.
	NSIM : (input) Number of sets to be tried for simulations to
		calculate the error estimates. If NSIM<=1 error estimates
		are not calculated.
	FE : (output) Array of length NR containing the estimated
		error in F[I]. This is calculated only if NSIM>1.
		
	Error status is returned by the value of the function RLS.
		0 value implies successful execution
		709 implies that NM<=NK+K-2 or IK<NM+NS or IV<NK+K-2
		710 implies that ALP<0 or IDE<1 or IDE>2
		other values may be set by BSPLIN or SVD or BSPEVL
	

	Required functions : BSPLIN, BSPEVL, SVD, SVDEVL, RANGAU
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
double rangau(double *seed);


int rls(int nk, double xo[], int k, int nr, double r[], double *rker,
	int ik, double *ac, int nm, int ns, double alp, int ide, double di[],
	double de[], double df[], double f[], double b[], int *iflg,
	double reps, double *chisq, double *sumd, double *a, double *av,
	int iv, double sigma[], int nsim, double fe[])

{
	int i,j,ir,nv,ne,ier,nderiv,is,left;
	double a1,h,fa,d0[2],xi,fi,s1,ssd;
	double *wk,*wb;

	nv=nk+k-2; ne=nm+ns;
	if(nm<=nv || ik<ne || iv<nv) return 709;
	if(alp<0 || ide<1 || ide>2) return 710;

	if(*iflg<2) {
/*     Setting up the system of linear equations */
		nderiv=0;
		for(i=0; i<nm; ++i) {

			for(j=0; j<nv; ++j) sigma[j]=0.0;
			h=(r[1]-r[0])/2.0;
			for(ir=0; ir<nr; ++ir) {
				ier=bsplin(xo,nk,k,r[ir],nderiv,df,av,&av[nv],&left);
				if(ier>100) return ier;

				for(j=0; j<nv; ++j) sigma[j]=sigma[j]+h*df[j]*rker[i+ir*ik];
				h=(r[ir+2]-r[ir])/2.0;
				if(ir==nr-2) h=(r[ir+1]-r[ir])/2.0;
			}

			for(j=0; j<nv; ++j) {
				a[j+i*iv]=sigma[j]/de[i];
				ac[i+j*ik]=sigma[j];
			}
		}
	}
	else if(*iflg<4) {
/*	The coefficients of matrix are available */
		for(i=0; i<nm; ++i) {
			for(j=0; j<nv; ++j) a[j+i*iv]=ac[i+j*ik]/de[i];
		}
	}

	if(*iflg<4) {
/*     The equations arising from regularisation term */
		h=(r[nr-1]-r[0])/(ns-1);
		nderiv=ide;
		fa=alp*sqrt(h);

		for(i=0; i<ns; ++i) {
			xi=r[0]+h*i;
			ier=bsplin(xo,nk,k,xi,nderiv,df,av,&av[nv],&left);
			if(ier>100) return ier;

			for(j=0; j<nv; ++j) a[j+(i+nm)*iv]=av[j+(ide-1)*nv]*fa;
		}

		ier=svd(nv,ne,a,av,sigma,iv,iv);
		if(ier>0) return ier;
	}

	if(*iflg==1 || *iflg==3) {*iflg=4; return 0;}
	*iflg=4;

/*	Set up the RHS of equations */
	for(i=0; i<nm; ++i) b[i]=di[i]/de[i];
	for(i=nm; i<nm+ns; ++i) b[i]=0.0;

/*	Solve the system of equations using SVD */
	ier=svdevl(nv,ne,a,av,sigma,iv,iv,b,reps);
	nderiv=0;
	for(i=0; i<nr; ++i) f[i]=bspevl(nk,xo,k,nderiv,b,r[i],&d0[0],&d0[1],&ier);

/*	Calculate the smoothing term */
	*sumd=0.0;
	nderiv=ide;
	h=(r[nr-1]-r[0])/(ns-1);
	for(i=0; i<ns; ++i) {
		xi=r[0]+i*h;
		fi=bspevl(nk,xo,k,nderiv,b,xi,&d0[0],&d0[1],&ier);
		*sumd=(*sumd)+d0[ide-1]*d0[ide-1];
	}
	*sumd=(*sumd)*h;

/*	Calculate the chi square */
	*chisq=0.0;
	for(i=0; i<nm; ++i) {
		s1=0.0;
		for(j=0; j<nv; ++j) s1=s1+ac[i+j*ik]*b[j];
		df[i]=(di[i]-s1)/de[i];
		*chisq=(*chisq)+df[i]*df[i];
	}
	if(nsim<2) return ier;
/*	Calculate the error estimate in the solution */
	ssd=123;	/* the seed for random number generator */
	nderiv=0;
	for(i=0; i<nr; ++i) fe[i]=0.0;
	wk=(double *) calloc((size_t) (nr*nsim), sizeof(double));
	wb=(double *) calloc((size_t) ne, sizeof(double));

	for(is=0; is<nsim; ++is) {
		for(i=0; i<nm; ++i) wb[i]=di[i]/de[i]+rangau(&ssd);
		for(i=nm; i<nm+ns; ++i) wb[i]=0.0;

		ier=svdevl(nv,ne,a,av,sigma,iv,iv,wb,reps);

		for(i=0; i<nr; ++i) {
			wk[i+is*nr]=bspevl(nk,xo,k,nderiv,wb,r[i],&d0[0],&d0[1],&ier);
			fe[i]=fe[i]+wk[i+is*nr];
		}
	}

	for(i=0; i<nr; ++i) {
		a1=fe[i]/nsim;
		s1=0.0;
		for(is=0; is<nsim; ++is) s1=s1+(wk[i+is*nr]-a1)*(wk[i+is*nr]-a1);
		fe[i]=sqrt(s1/nsim);
	}
	free(wb); free(wk);
	return 0;
}

