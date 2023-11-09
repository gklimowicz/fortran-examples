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

