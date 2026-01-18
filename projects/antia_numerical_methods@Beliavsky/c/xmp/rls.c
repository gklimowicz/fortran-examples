/*	To solve an inverse problem with specified kernels.
	It first generates artificial data using a prescribed solution using
	forw and then these data are used for inversion */

#include <stdio.h>
#include <math.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);


int forw(int np, int nm, double r[], double *rker, int ik, double di[],
	double f[], double fun(double ), int *iflg);
int rls(int nk, double xo[], int k, int nr, double r[], double *rker,
	int ik, double *ac, int nm, int ns, double alp, int ide, double di[],
	double de[], double df[], double f[], double b[], int *iflg,
	double reps, double *chisq, double *sumd, double *a, double *av,
	int iv, double sigma[], int nsim, double fe[]);

double fun(double x);
double rangau(double *seed);

double aker[700][2000];

main()
{
	int i,i1,j,n,k,nk,nm,nsim,id,iflg, ier,np,ik,iv;
	double hh,seed, xo[100], a[2000][100],av[100][100],f[1000],b[2000],erc,reps,aeps;
	double sigma[100],r[1000],ac[100][2000],di[2000],fe[1000],
			alp,chi,sum,df[2000],de[2000],f0[1000];
	FILE *fp,*fp1;

/*	Exercise 13.23 */

	reps=1.e-6; aeps=1.e-7; erc=0.2;
/*	Seed for random numbers */
	seed=1;
	iflg=0;

/*	Read in the kernels */
	fp=fopen("kernel","r");
	fscanf(fp," %d",&np);
	for(i=0; i<np; ++i) fscanf(fp," %le",&r[i]);
	nm=1387;
	for(i=0; i<nm; ++i) {
		fscanf(fp," %d",&j);
		for(j=0; j<np; ++j) fscanf(fp," %le",&aker[j][i]);
		de[i]=erc;
	}
	fclose(fp);
	printf(" no. of points in r = %d     no. of modes = %d\n",np,nm);

	iflg=0;
	ik=2000; iv=100;
/*	Calculate the data points using forw */
	ier=forw(np,nm,r,&aker[0][0],ik,di,f0,fun,&iflg);
	printf(" ier = %d\n",ier);

/*	Add random errors to calculated data points */
	for(i=0; i<nm; ++i) di[i]=di[i]+de[i]*rangau(&seed);

/*	No. of points for smoothing */
	n=400;
/*	No. of sets for simulation */
	nsim=20;

	for(i1=0; i1<99; ++i1) {
		printf("type k=order of B-splines,   nk=no. of knots, \n");
		printf("     id=order of derivative for smoothing,    iflg,\n");
		printf("     alp=regularisation parameter \n");
		printf("                       (quits when k<=1)\n");
		scanf(" %d %d %d %d %le",&k,&nk,&id,&iflg,&alp);
		if(k<=1) return 0;

/*	Set up the knots with uniform spacing */
		hh=(r[np-1]-r[0])/(nk-1);
		for(i=0; i<=nk; ++i) xo[i]=r[0]+hh*i;

		i=rls(nk,xo,k,np,r,&aker[0][0],ik,&ac[0][0],nm,n,alp,id,di,de,df,f,b,
			&iflg,reps,&chi,&sum,&a[0][0],&av[0][0],iv,sigma,nsim,fe);
		printf(" ier = %d    order of B-spline =  %d    no. of knots = %d\n",i,k,nk);
		printf(" order of derivative for smoothing = %d     regularisation parameter = %e \n",id,alp); 
		printf("  chi squre = %e    regularisation term = %e \n",chi,sum); 
		printf(" coefficients : ");
		for(i=0; i<nk+k-2; ++i) {
			printf(" %e ", b[i]);
			if(i - 5*(i/5) ==4) printf("\n");
		}
		printf("\n");

/*	Write the results in invers.out,  The first nm points give the
	residuals in fitted data points, subsequent points give the inverted
	function compared with exact value */

		fp1=fopen("invers.out","w");
		fprintf(fp1,"#S. No.  Input data     error         residuals \n");
		for(i=0; i<nm; ++i) fprintf(fp1," %d    %e    %e    %e \n",i,di[i],de[i],df[i]);

		fprintf(fp1,"#    R       Inverted soln    error      exact value \n");
		for(i=0; i<np; ++i) fprintf(fp1," %e %e %e %e\n",r[i],f[i],fe[i],f0[i]);
		fclose(fp1);


	}
	return;
}


/*	The function to generate artificial data */

double fun(double x)

{
	return 430+20*sin(30*x);
}



/*	To calculate function value using B-spline expansion

	N : (input) Number of knots to define B-splines
	X : (input) Array of length N containing the knots.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	NDERIV : (input) Number of derivatives required
		For NDERIV<=0 only function value is calculated
		For NDERIV=1 first derivative is also calculated
		For NDERIV>1 both first and second derivatives are calculated
	WT : (input) Coefficients of B-spline expansion
	X0 : (input) The point at which expansion has to be evaluated
	DF : (output) First derivative of function at X0
	DDF : (output) Second derivative of function at X0
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values of IER may be set by BSPLIN which is called

	BSPEVL = SUM_{i=1}^{N+K-2} WT(I) \phi_i(X0)
	where \phi_i(x) are B-spline basis functions on knots X

	Required functions : BSPLIN
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier)

{
	int i,left;
	double f;
	double *wk, *wk1, *wk2;

	wk=(double *) calloc((size_t) (n+k), sizeof(double));
	wk1=(double *) calloc((size_t) (n+k), sizeof(double));
	wk2=(double *) calloc((size_t) (n+k), sizeof(double));
	*ier=bsplin(x,n,k,x0,nderiv,wk,wk1,wk2,&left);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	f=0.0; *df=0.0; *ddf=0.0;
	for(i=left; i<=left+k-1; ++i) {
		f = f + wt[i]*wk[i];
		*df = *df + wt[i]*wk1[i];
		*ddf = *ddf + wt[i]*wk2[i];
	}
	free(wk2); free(wk1); free(wk);
	return f;
}



/*	To calculate the B-spline basis functions at a specified point

	X : (input) Array of length NX containing the knots.
		The knots must be distinct and in ascending order.
	NX : (input) Number of knots
	K : (input) Order of B-spline, 0< K, K=4 gives cubic B-splines
	XB : (input) The point at which B-spline basis functions are to be evaluated
	NDERIV : (input) Number of derivatives required
		NDERIV<=0 only B-splines are calculated
		NDERIV=1 first derivative is also calculated
		NDERIV>1 first and second derivatives are also calculated
	B : (output) Array of length NX+K-2 containing the value of
		B-spline basis functions
	DB : (output) Array of length NX+K-2 containing the value of
		the first derivative of B-spline basis functions (if NDERIV>0)
	DDB : (output) Array of length NX+K-2 containing the value of
		the second derivative of B-spline basis functions (if NDERIV>1)
	LEFT : (output) XB is located between X[LEFT] and X[LEFT+1]
		
	Error status is returned by the value of the function BSPLIN.
		0 value implies successful execution
		26 implies XB > X[NX-1]
		27 implies XB < X[0]
		203 implies NX<2, K<1

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left)

{
	int i,j,igh, mid,nigh,lx,ier;
	double t1,t2,t3,p1,p2;
	double *wk, *dr, *dl;
	static int low = -1;

	if(nx <= 1 || k<1 ) return 203;
	ier=0;

/*	If the previous value of LOW is inadmissible, set the range to (0,N-1) */
	if(low<0 || low>=nx-1) {low=0; igh=nx-1;}
	else igh=low+1;

	while((xb<x[low] && xb<x[igh]) || (xb>x[low] && xb>x[igh])) {
/*	Extend the range */
		if( xb>x[low] ) {
/*	Extend the range on higher side */
			if(igh >= nx-1) {ier=26; low=nx-2; break;}
			else {
				nigh=igh+2*(igh-low); if(nx-1 < nigh) nigh=nx-1;
				low=igh; igh=nigh;
			}
		}

		else {
/*	Extend the range on lower side */
			if(low <= 0) {ier=27; igh=low+1; break;}
			else {
				nigh=low;
				low=low-2*(igh-low); if(low<0) low=0;
				igh=nigh;
			}
		}
	}


/*	Once the point is bracketed between two tabular points locate it by bisection */
	while((igh-low > 1) && (xb != x[low])) {
		mid=(low+igh)/2;
		if((xb<= x[mid]) == (xb<= x[low])) low=mid;
		else igh=mid;
	}

 
/*	Evaluate the B-spline basis functions

	Define the extra knots on either side of table
	Note that the function assumes knots from -K+2 to NX+K-1
	and the B-splines B_{i,k}, i ranges from 0 to NX+K-3 
	The knots are stored in scratch array wk. */

	wk=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));
	dr=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));
	dl=(double *) calloc((size_t) (nx+2*k+2), sizeof(double));

	for(i=0; i<nx; ++i) wk[i+k]=x[i];
	for(i=1; i<=k; ++i) {
		wk[k-i]=x[0];
		wk[nx+i-1+k]=x[nx-1];
	}

	for(i=0; i<nx+k-2; ++i) {b[i]=0.0; db[i]=0.0; ddb[i]=0.0;}
	*left=low;
	lx=low-1;
	b[lx+1]=1;

/*	The recurrence relation for B-splines */
	for(j=1; j<=k-1; ++j) {
		dr[j] = wk[low+j+k] - xb;
		dl[j] = xb - wk[low+1-j+k];
		t1=0.0;
		for(i=1; i<=j; ++i) {
			t2=b[lx+i]/(dr[i]+dl[j+1-i]);
			b[lx+i]=t1+t2*dr[i];
			t1=t2*dl[j+1-i];
		}
		b[lx+j+1]=t1;
			
/*	Calculate the first derivative using recurrence relations */
		if(j == k-2 && nderiv > 0) {
			t1=0.0;
			for(i=1; i<=j+1; ++i) {
				t2=b[lx+i]/(wk[low+i+k]-wk[low+i+1]);
				db[lx+i]=(k-1)*(t1-t2);
				t1=t2;
			}
			db[lx+j+2]=(k-1)*t1;
		}
 
/*	Calculate the second derivative using recurrence relations */
		if(j == k-3 && nderiv>1) {
			t2=0.0; p1=0.0;
			for(i=1; i<=j+1; ++i) {
				t3=b[lx+i]/(wk[low+i+k]-wk[low+i+2]);
				p2=(t2-t3)/(wk[low+i+k]-wk[low+i+1]);
				ddb[lx+i]=(k-2)*(k-1)*(p1-p2);
				t2=t3; p1=p2;
			}
			p2=t2/(wk[low+j+2+k]-wk[low+j+3]);
			ddb[lx+j+2]=(k-2)*(k-1)*(p1-p2);
			ddb[lx+j+3]=(k-2)*(k-1)*p2;
		}
	}

/*	For K=2 the first derivative has to be calculated outside the loop */
	if(k == 2 && nderiv > 0) {
		t2=1./(wk[low+1+k]-wk[low+2]);
		db[lx+1]=-t2;
		db[lx+2]=t2;
	}

/*	For K=3 the second derivative has to be calculated outside the loop */
	if(k == 3 && nderiv > 1) {
		t3=1./(wk[low+1+k]-wk[low+3]);
		p2=-t3/(wk[low+1+k]-wk[low+2]);
		ddb[lx+1]=-2.*p2;
		p1=p2;
		p2=t3/(wk[low+2+k]-wk[low+3]);
		ddb[lx+2]=2.*(p1-p2);
		ddb[lx+3]=2.*p2;
	}
	free(dl); free(dr); free(wk);
	return ier;
}



/*	To solve the forward problem corresponding to a linear inverse
	problem. This function may be used to generate artificial data
	for testing inversion techniques.

	NP : (input) Number of points used in defining the kernels.
	NM : (input) Number of data points in the inverse problem
		which should be same as the number of kernels that are
		supplied in array RKER.
	R : (input) Array of length NP containing the coordinates
		of points at which kernels are available.
	RKER : (input) Array of length IK*NP containing the kernels
		for the inverse problem. RKER[j][i] should contain the
		value at R[j] for the ith kernel.
	IK : (input) The second dimension of RKER as declared in the calling
		function
	DI : (output) Array of length NM containing the calculated
		data points using the kernel.
	F : (input/output) Array of length NP containing the function
		value at points in R. F[I] should contain the function
		value at R[I]. If IFLG=0, the function values are
		calculated using user supplied function FUN,
		otherwise, these values must be supplied while calling
		the function.
	FUN : (input) Name of function routine to calculate the given
		function. This is used only if IFLG=0, otherwise the
		function values are to be supplied in array F.
	IFLG : (input/output) Integer parameter used as a flag to decide
		the type of computation required.
		If IFLG=0, then the function values are calculated using
			a user supplied function FUN. These values are stored
			in array F and IFLG is set to 1 so that next time
			the values need not be calculated.
		For other values of IFLG the function values must be
			supplied in array F.
		
	Error status is returned by the value of the function FORW.
		0 value implies successful execution
		711 implies that IK<NM and no calculations are done

	Function FUN(X) must be supplied by the user

	Required functions : FUN
*/

#include <math.h>

int forw(int np, int nm, double r[], double *rker, int ik, double di[],
	double f[], double fun(double ), int *iflg)

{
	int i,ir;
	double h,s1;

	if(ik<nm) return 711;

	if(*iflg==0) {
/*     Calculate the function value using supplied routine */
		for(i=0; i<np; ++i) f[i]=fun(r[i]);
		*iflg=1;
	}
 
/*     Calculate the integrals */
	for(i=0; i<nm; ++i) {
		s1=0.0;
		h=(r[1]-r[0])/2.0;
		for(ir=0; ir<np; ++ir) {
			s1=s1+h*f[ir]*rker[i+ir*ik];
			h=(r[ir+2]-r[ir])/2.0;
			if(ir==np-2) h=(r[ir+1]-r[ir])/2.0;
		}
		di[i]=s1;
	}
	return 0;
}




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




/*	To calculate the Singular Value Decomposition of a matrix A=U D Vtranspose

	N : (input) Number of variables
	M : (input) Number of equations
	A : (input/output) Matrix of coefficients of size LA*M
		After execution it will contain the matrix U
	V : (output) The matrix V of size LV*N
	SIGMA : (output) Array of length N, containing the singular values
	LA : (input) Actual value of second dimension of A in the calling function
	LV : (input) Actual value of second dimension of V in the calling function
		
	Error status is returned by the value of the function SVD.
		0 value implies successful execution
		12 QR iteration failed to converge to required accuracy
		105 implies N<=0, N>LV, M<=0, N>LA, N>M

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv)

{
	int i,j,k,l,itr,ier, itmax=30;
	double f, g, h, rmax, s, r1, r2, c, x, y, z, aeps, eps=1.e-16;
	double *e;

	if(n>m || n<=0 || m<=0 || n>la || n>lv) return 105;
	ier=0;

/*	Reduction to Bidiagonal form using Householder transformations */
	g=0.0; rmax=0.0;
	e=(double *) calloc((size_t) n, sizeof(double));

	for(i=0; i<n; ++i) {
/*	Off-diagonal elements of bidiagonal form  */
		e[i]=g;
		s=0.0;
		for(j=i; j<m; ++j) s=s+a[i+j*la]*a[i+j*la];
		if(s <= 0.0) {
/*	transformation not required */
			g=0.0;
		}
		else {
			f= a[i+i*la];
			g=sqrt(s);
			if(f>=0.0) g=-g;
			h=f*g-s;
			a[i+i*la] = f-g;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]= a[j+k*la]+f*a[i+k*la];
			}
		}

/*	Diagonal elements of bidiagonal form  */
		sigma[i]=g;
		s=0.0;
		for(j=i+1; j<n; ++j) s=s+a[j+i*la]*a[j+i*la];

		if(s<= 0.0) g=0.0;
		else {
			f= a[i*la+(i+1)];
			g=sqrt(s);
			if(f>= 0.0) g=-g;
			h=f*g-s;
			a[i*la+(i+1)]=f-g;
			for(j=i+1; j<n; ++j) e[j]=a[j+i*la]/h;

			for(j=i+1; j<m; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+j*la]*a[k+i*la];
				for(k=i+1; k<n; ++k) a[k+j*la] = a[k+j*la]+s*e[k];
			}
		}
		r1=fabs(sigma[i])+fabs(e[i]);
		if(r1 > rmax) rmax=r1;
	}

/*	Accumulation of right hand transformation in array V */
	for(i=n-1; i>=0; --i) {
		if(g != 0.0) {
			h=g*a[i*la+(i+1)];
			for(j=i+1; j<n; ++j) v[i+j*lv]=a[j+i*la]/h;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+i*la]*v[j+k*lv];
				for(k=i+1; k<n; ++k) v[j+k*lv]=v[j+k*lv]+s*v[i+k*lv];
			}
		}

		for(j=i+1; j<n; ++j) {
			v[j+i*lv]=0.0; v[i+j*lv]=0.0;
		}
		v[i+i*lv]=1;
		g= e[i];
	}

/*	Accumulation of left hand transformation overwritten on matrix A */
	for(i=n-1; i>=0; --i) {
		g=sigma[i];
		for(j=i+1; j<n; ++j) a[j+i*la]=0.0;
		if(g != 0.0) {
			h=g*a[i+i*la];

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]=a[j+k*la]+f*a[i+k*la];
			}

			for(j=i; j<m; ++j) a[i+j*la]=a[i+j*la]/g;
		}
		else {
			for(j=i; j<m; ++j) a[i+j*la]=0.0;
		}
		a[i+i*la] = a[i+i*la]+1;
	}

/*	Diagonalisation of the bidiagonal form */
	aeps=eps*rmax;
/*	Loop over the singular values */
	for(k=n-1; k>=0; --k) {
/*	The QR transformation */
		for(itr=1; itr<=itmax; ++itr) {

/*	Test for splitting */
			for(l=k; l>=0; --l) {
				if(fabs(e[l]) < aeps) goto split;
				if(fabs(sigma[l-1]) < aeps) break;
			}

/*	cancellation of E[L] if L>1  */
			c=0.0; s=1.0;
			for(i=l; i<=k; ++i) {
				f=s*e[i];
				e[i] = c*e[i];
				if(fabs(f) < aeps) goto split;
				g=sigma[i];
				sigma[i]=sqrt(f*f+g*g);
				c=g/sigma[i];
				s=-f/sigma[i];

				for(j=0; j<m; ++j) {
					r1= a[j*la+(l-1)];
					r2= a[i+j*la];
					a[j*la+(l-1)]=r1*c+r2*s;
					a[i+j*la]=c*r2-s*r1;
				}
			}

split:			z=sigma[k];
			if(l == k) {
/*	QR iteration has converged */
				if(z < 0.0) {
					sigma[k] = -z;
					for(j=0; j<n; ++j) v[k+j*lv]=-v[k+j*lv];
				}
				break;
			}

			if(itr==itmax) {ier=12; break;}
					
/*	calculating shift from bottom 2x2 minor */
			x=sigma[l];
			y=sigma[k-1];
			g=e[k-1];
			h=e[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g=sqrt(1.+f*f);
			if(f < 0.0) g=-g;
			f=((x-z)*(x+z)+h*(y/(f+g)-h))/x;

/*	next QR transformation */
			c=1.0; s=1.0;
/*	Given's rotation  */
			for(i=l+1; i<=k; ++i) {
				g=e[i];
				y=sigma[i];
				h=s*g;
				g=c*g;
				e[i-1]=sqrt(f*f+h*h);
				c=f/e[i-1];
				s=h/e[i-1];
				f=c*x+s*g;
				g=c*g-s*x;
				h=s*y;
				y=c*y;

				for(j=0; j<n; ++j) {
					x=v[j*lv+(i-1)];
					z=v[i+j*lv];
					v[j*lv+(i-1)]=c*x+s*z;
					v[i+j*lv]=c*z-s*x;
				}

				sigma[i-1]=sqrt(f*f+h*h);
				if(sigma[i-1] != 0.0) {
					c=f/sigma[i-1];
					s=h/sigma[i-1];
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for(j=0; j<m; ++j) {
					y= a[j*la+(i-1)];
					z= a[i+j*la];
					a[j*la+(i-1)] = c*y+s*z;
					a[i+j*la] = c*z-s*y;
				}
			}

			e[l]=0.0;
			e[k]=f;
			sigma[k]=x;
		}
	}
	free(e);
	return ier;
}



/*	To evaluate the solution of a system of linear equations using SVD

	N : (input) Number of variables
	M : (input) Number of equations
	U : (input) Array of size LU*M containing the left-hand transformation
	V : (input) Array of size LV*N containing the right-hand transformation
	SIGMA : (input) Array of size N containing the singular values
	LU : (input) Second dimension of array U in the calling function
	LV : (input) Second dimension of array V in the calling function
	B : (input/output) Array of length M containing the RHS
		after execution it will contain the solution
	REPS : (input) Relative accuracy. All singular values < REPS*(Max of singular values)
		will be reduced to zero

	The returned value is always zero

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps)

{
	int i,j;
	double smax, aeps, s;
	double *wk;

/*	Finding the largest singular value */
	smax=0.0;
	for(i=0; i<n; ++i)
		if(sigma[i] > smax) smax=sigma[i];

	aeps=smax*reps;
	wk=(double *)calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i) {
		s=0.0;
/*	Only SIGMA[I] > AEPS contribute to the solution */
		if(sigma[i] > aeps) {
			for(j=0; j<m; ++j) s=s+b[j]*u[i+j*lu];
			s=s/sigma[i];
		}
		wk[i]=s;
	}

	for(i=0; i<n; ++i) {
		s=0.0;
		for(j=0; j<n; ++j) s=s+v[j+i*lv]*wk[j];
		b[i]=s;
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

