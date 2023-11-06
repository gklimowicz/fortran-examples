/*	Two point boundary value problem using expansion methods with
	B-spline basis functions */

#include <stdio.h>
#include <math.h>

void eqn(int j, int m, int ml, double par[], double *a, double *b,
	double x[], double f[], double t);
void bcs(int m, int ml, double par[], double *bc, double g[], double t1,
	double tn, double x1[], double xn[]);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);

int bspode(int nk, int k, int m, int ml, double par[], double *x,
	double *a, double t[], int n, double tx[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
	int iflag, double reps);


main()
{
	int i,i1,j,k,n,m,ml,id, iflg, ier,np,nmax;
	double hh,t[500],x[500][2],tx[500],tn,par[20],a[500],reps,t0,f1,f2;
	double df,ddf;

/*	Example 12.13 */

	m=2; ml=1; reps=1.e-8;

/*	Pass the parameter lambda in equation via par[0] */

	for(i1=0; i1<99; ++i1) {
		printf("type n=no. of knots,   k=order of B-spline,   np=no. of points,   lambda\n");
		printf("                      (quits when n<=0)\n");
		scanf(" %d %d %d %le",&n, &k, &np, &par[0]);
		if(n<=0) return 0;

		iflg=0;
/*	Set up the knots and initial values of coefficients */
		hh=1.0/(n-1);
		for(i=0; i<n; ++i) t[i]=i*hh;
		for(i=0; i<n+k-2; ++i) {
			a[i]=1.0;
			a[i+n+k-2]=1.0;
		}

		i=bspode(n,k,m,ml,par,&x[0][0],a,t,np,tx,eqn,bcs,iflg,reps);
		printf(" ier = %d    no. of knots = %d    order of B-splines = %d\n",i,n,k);
		printf("  no. of points = %d    lambda = %e\n",np,par[0]);

/*	Calculate the solution at required points using bspevl
	Also calculate the derivative of x[0] and compare it with x[1] */

		printf("    t          solution                  f'(t) \n");
		nmax=1.0;
		for(i=0; i<6; ++i) {
			t0=0.2*i;
			f1=bspevl(n,t,k,nmax,a,t0,&df,&ddf,&ier);
			f2=bspevl(n,t,k,nmax,&a[n+k-2],t0,&ddf,&ddf,&ier);

			printf(" %e %e %e %e\n",t0,f1,f2,df);
		}

/*
		printf("\n  Coefficients : ");
		for(i=0; i<m*(n+k-2); ++i) printf(" %e",a[i]);
*/

	}
	return;
}

/*	The differential equation */

void eqn(int j, int m, int ml, double par[], double *a, double *b,
	double y[], double f[], double t)

{
	int i,k;

	for(i=0; i<m; ++i) {
		f[i]=0.0;
		for(k=0; k<m; ++k) {
			a[k+i*m]=0.0;
			b[k+i*m]=0.0;
		}
	}
	b[0]=1.0;
	a[m]=1.0;
	b[1+m]=1.0;
	a[1]=par[0]*par[0]*cosh(par[0]*y[0]);

	f[0]=y[1];
	f[1]=par[0]*sinh(par[0]*y[0]);
	return;
}

/*	The boundary conditions */

void bcs(int m, int ml, double par[], double *bc, double g[],
	double t1, double tn, double y1[], double yn[])

{
	int i,k;

	bc[0]=1.0;
	bc[m]=0.0;
	bc[1]=1.0;
	bc[1+m]=0.0;

	g[0]=y1[0];
	g[1]=yn[0]-1.0;
	return;
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


/*	To solve two-point boundary value problem in ordinary differential
	equations using expansion method with B-spline basis functions

	NK : (input) Number of knots to be used for calculating B-splines.
	K : (input) Order of B-splines to be used, K=4 for cubic B-splines
	M : (input) Number of first order differential equations in the system
	ML : (input) Number of boundary conditions at the first boundary t=T[0]
	PAR : (input) Array to be passed on to EQN and BCS for calculating
		the equations and boundary conditions. This array is not used
		by BSPODE, but is merely used to pass on any extra parameters
		that may be required to specify the equations
	X : (output) Array of length M*N containing the solution.
		X[j][i] is the ith component of solution at jth mesh point
		TX[j]. Second dimension of X in calling function must be M
	A : (input/output) Array of length (NK+K-2)*M containing the
		coefficients of expansion. At the time of calling it should
		contain the initial guess. After execution it will contain
		the final values. The second dimension of A must be NK+K-2.
	T : (input) Array of length NK containing the knots.
		These points must be in ascending order with T[0] and T[NK-1]
		as the two boundaries. 
	N : (input) The number of mesh points to be used for calculating
		the coefficients. The solution will be calculated at
		all these points. N>=NK+K-2
	TX : (input/output) Array of length N containing the mesh points
		For IFLAG>1 the mesh points must be supplied, while for
		other values of IFLAG the function calculates the mesh point
		assuming uniform spacing. 
	EQN : (input) Name of function to specify the differential equation
		to be solved.
	BCS : (input) Name of function to calculate the boundary conditions
		at t=T[0] and T[NK-1]
	IFLAG : (input) Integer variable used as a flag to decide the type of
		computation required.
		IFLAG=0 implies that equations are nonlinear and mesh points
			are not supplied in array TX. These would be calculated
			assuming uniform spacing.
		IFLAG=1 implies that equations are linear and mesh points
			are not supplied in array TX. These would be calculated
			assuming uniform spacing.
		IFLAG=2 implies that equations are nonlinear and mesh points
			are supplied in array TX.
		IFLAG=3 implies that equations are linear and mesh points
			are supplied in array TX.
	REPS : (input) Required accuracy. This is only used to check convergence
		of Newton's method for solving the resulting system of equations.
		The truncation error depends on NK and K and is not
		controlled in this function.
		
	Error status is returned by the value of the function BSPODE.
		0 value implies successful execution
		705 implies that NK<3, M<=ML, ML<=0, or N<NK+K-2
			in which case no calculations are done
		741 implies that Newton's iteration for solving the
			system of equations failed to converge.
		Other values may be set by BSPLIN or SVD

	Functions EQN and BCS must be supplied by the user.
	Function EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
		sides for differential equations y'_i=f_i(t,y,par)
		J is the serial number of mesh point at which calculation
		is required. M is the number of first order differential
		equations in the system, ML the number of boundary conditions
		at the first point, PAR is an array which can be used
		to pass on any required parameters to define the equations.
		A and B are arrays of length M*M defining the
		differential equation B Y'=f(T,PAR,Y) and
		A[K][I]=dF_I/dY(K) is the Jacobian matrix
		Y is an array of length M specifying the solution at t=T
		F is an array of length M containing the right hand
		sides of differential equations f_I as defined above.
		F, A and B must be calculated by the function.

	Function BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) calculates the boundary
		conditions at both boundaries. M is the number of differential
		equations, ML is the number of boundary condition at the
		first mesh point t=T1. PAR is an array which can be used
		to pass any required parameters. BC is an array of length
		M*M which should contain the coefficients of boundary
		conditions. First ML rows will specify the boundary conditions
		at t=T1, while remaining rows will specify those at t=TN.
		G is an array of length M specifying the boundary conditions
		G[I]=g_i(T1,PAR,Y1) (I<ML) are the boundary conditions
		at T1 (g_i=0), while G[I]=g_i(TN,PAR,YN) (I>=ML) are the
		boundary conditions at TN. BC is the Jacobian matrix dg_i/dY(K)
		Y1 and YN are arrays of length M specifying the
		solution at t=T1 and TN.

	Required functions : BSPLIN, BSPEVL, SVD, SVDEVL, EQN, BCS
*/

#include <math.h>
#include <stdlib.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
		double *df, double *ddf, int *ier);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);

int bspode(int nk, int k, int m, int ml, double par[], double *x,
	double *a, double t[], int n, double tx[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
	int iflag, double reps)

{
	int i,j,nb,ne,me,nderiv,it,j1,ji,jk,ier,left, nit=50;
	double h,xb,df,ddf,s,t1,t2,xj,rerr,re,r2;
	double eps=1.0e-30;
	double *wk, *v;

	nb=nk+k-2;
	if(nk<3 || m<=ml || ml<=0 || n<nb) return 705;

	if(iflag<=1) {
/*	Setup the mesh assuming uniform spacing */
		h=(t[nk-1]-t[0])/(n-1);
		for(i=0; i<n; ++i) tx[i]=t[0]+i*h;
	}
	ne=n*m+m; me=nb*m;
	wk=(double *) calloc((size_t) (ne*(me+7)), sizeof(double));

/*	The iteration loop */
	nderiv=1;
	for(it=1; it<=nit; ++it) {

/*	Setup the equation matrix */
		for(i=0; i<n; ++i) {
			xb=tx[i];
			ier=bsplin(t,nk,k,xb,nderiv,&wk[me*ne],&wk[me*ne+ne],
					&wk[(me+2)*ne],&left);
			if(ier>100) {free(wk); return ier;}
			for(j=0; j<m; ++j) {
				x[j+i*m]=bspevl(nk,t,k,nderiv,&a[nb*j],xb,&df,&ddf,&ier);
				if(ier>100) {free(wk); return ier;}
				x[j+(i+1)*m]=df;
			}
			eqn(i,m,ml,par,&wk[(me+2)*ne],&wk[(me+3)*ne],&x[i*m],
					&wk[(me+4)*ne],xb);

			for(j=0; j<m; ++j) {
				s=0.0;
				for(j1=0; j1<m; ++j1) s=s+wk[j+j1*m+(me+3)*ne]*x[j1+(i+1)*m];
				wk[j+i*m+(me+5)*ne]=-s+wk[j+(me+4)*ne];
				for(j1=0; j1<me; ++j1) {
					jk=j1/nb; ji=j1-jk*nb;
					wk[(j+i*m)*me+j1]=wk[j+jk*m+(me+3)*ne]*wk[ji+(me+1)*ne]-
						wk[j+jk*m+(me+2)*ne]*wk[ji+me*ne];
				}
			}
		}

		t1=t[0];
		ier=bsplin(t,nk,k,t1,nderiv,&wk[me*ne],&wk[(me+1)*ne],
				&wk[(me+2)*ne],&left);
		if(ier>100) {free(wk); return ier;}
		t2=t[nk-1];
		ier=bsplin(t,nk,k,t2,nderiv,&wk[me*ne+ne],&wk[(me+2)*ne],
				&wk[(me+3)*ne],&left);
		if(ier>100) {free(wk); return ier;}

		bcs(m,ml,par,&wk[(me+2)*ne],&wk[(me+3)*ne],t1,t2,x,&x[(n-1)*m]);

		for(i=0; i<m; ++i) {
			wk[i+n*m+(me+5)*ne]=-wk[i+(me+3)*ne];
			for(j1=0; j1<me; ++j1) {
				jk=j1/nb; ji=j1-jk*nb;
				if(i<ml) wk[me*(i+n*m)+j1]=wk[i+jk*m+(me+2)*ne]*wk[ji+me*ne];
				else wk[me*(i+n*m)+j1]=wk[i+jk*m+(me+2)*ne]*wk[ji+me*ne+ne];
			}
		}

/*	Solve the system of equations using SVD */
		v=(double *) calloc((size_t) (me*me), sizeof(double));
		ier=svd(me,ne,wk,v,&wk[me*ne],me,me);
		if(ier>100) {free(v); free(wk); return ier;}

		ier=svdevl(me,ne,wk,v,&wk[me*ne],me,me,&wk[(me+5)*ne],reps);
		free(v);

/*		The convergence check */
		rerr=0.0;
		for(j=0; j<nb; ++j) {
			j1=j+1; if(j==nb-1) j1=j-1;
			for(i=0; i<m; ++i) {
				jk=j+i*nb;
				xj=a[jk]+wk[jk+(me+5)*ne];
				r2=fabs(xj)+fabs(x[j+i*m]-x[j1+i*m])+eps;
				re=fabs(wk[jk+(me+5)*ne]/r2);
				if(re>rerr) rerr=re;
				a[jk]=xj;
			}
		}
		if(rerr<reps || iflag==1 || iflag==3) {free(wk); return 0;}
	}

	free(wk);
	return 741;
}

/*	--------------------------------

void eqn(int j, int m, int ml, double par[], double *a, double *b,
	double y[], double f[], double t)

{
	int i,k;

	for(i=0; i<m; ++i) {
		f[i]=f_i(t,par,y);
		for(k=0; k<m; ++k) {
			a[k+i*m]=df_k/dy[i];
			b[k+i*m]=b_{ki}(t,par);
		}
	}
	return;
}

void bcs(int m, int ml, double par[], double *bc, double g[],
	double t1, double tn, double y1[], double yn[])

{
	int i,k;

	for(i=0; i<m; ++i) {
		if(i<ml) g[i]=g_i(t1,par,y1);
		else g[i]=g_i(tn,par,yn);

		for(k=0; k<m; ++k) bc[i+k*m]=dg_i/dy[k];
	}
	return;
}

*/




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
