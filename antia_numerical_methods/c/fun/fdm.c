/*	To solve two-point boundary value problem in ordinary differential
	equations using finite difference method

	N : (input) Number of mesh points to be used.
	M : (input) Number of first order differential equations in the system
	ML : (input) Number of boundary conditions at the first boundary t=T[0]
	PAR : (input) Array to be passed on to EQN and BCS for calculating
		the equations and boundary conditions. This array is not used
		by FDM, but is merely used to pass on any extra parameters
		that may be required to specify the equations
	X : (input/output) Array of length M*N containing the solution.
		It should contain the initial guess to solution at the time
		of calling. After execution it will contain the calculated
		solution. X[j][i] is the ith component of solution at jth mesh point.
		Second dimension of X in calling function must be M
	XC : (output) Array of length M*N containing the solution after
		applying the deferred correction. It is stored in the same
		format as X.
	T : (input) Array of length N containing the mesh points.
		These points must be in ascending or descending order.
		For calculating the deferred correction the mesh spacing
		must be uniform.
	EQN : (input) Name of function to specify the differential equation
		to be solved.
	BCS : (input) Name of function to calculate the boundary conditions
		at t=T[0] and T[N-1]
		Outline of a sample for EQN and BCS can be found
		at the end of this file
	IFLAG : (input) Integer variable used as a flag to decide the type of
		computation required.
		IFLAG=0 implies that equations are nonlinear and deferred
			correction is to be calculated. In this case X must
			contain the initial guess and mesh spacing must be uniform
		IFLAG=1 implies that equations are nonlinear and deferred
			correction is not required. In this case X must contain
			the initial guess and mesh spacing could be arbitrary
		IFLAG=2 implies that equations are linear and deferred correction
			is required. In this case initial guess is not required,
			but mesh spacing must be uniform
		IFLAG=3 implies that equations are linear and deferred correction
			is not required. In this case initial guess is not
			required and mesh spacing can be arbitrary
	REPS : (input) Required accuracy. This is only used to check convergence
		of Newton's method for solving finite difference equations.
		The truncation error depends on mesh spacing and is not
		controlled in this function.
		
	Error status is returned by the value of the function FDM.
		0 value implies successful execution
		704 implies that N<3, M<=ML, or ML<=0, in which case
			no calculations are done
		734 implies that N<5 and deferred correction is requested.
			In this case deferred correction is not calculated
		735 implies that the finite difference matrix is singular
		736 implies that mesh spacing is not uniform and
			deferred correction is not calculated
		737 implies that Newton's iteration for solving the
			finite difference equations failed to converge.

	Functions EQN and BCS must be supplied by the user.
	Function EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
		sides for differential equations By'_i=f_i(t,y,par)
		J is the serial number of mesh point at which calculation
		is required. M is the number of first order differential
		equations in the system, ML the number of boundary conditions
		at the first point, PAR is an array which can be used
		to pass on any required parameters to define the equations.
		A and B are arrays of length (M+ML)*M defining the
		differential equation B Y'=f(T,PAR,Y) and
		A[K][I]=dF_I/dY[K] is the Jacobian matrix
		Y is an array of length M specifying the solution at t=T
		F is an array of length M containing the right hand
		sides of differential equations f_I as defined above.
		F, A and B must be calculated by the function.

	Function BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) calculates the boundary
		conditions at both boundaries. M is the number of differential
		equations, ML is the number of boundary condition at the
		first mesh point t=T1. PAR is an array which can be used
		to pass any required parameters. BC is an array of length
		(M+ML)*M which should contain the coefficients of boundary
		conditions. First ML rows will specify the boundary conditions
		at t=T1, while remaining rows will specify those at t=TN.
		G is an array of length M specifying the boundary conditions
		G[I]=g_i(T1,PAR,Y1) (I<ML) are the boundary conditions
		at T1 (g_i=0), while G[I]=g_i(TN,PAR,YN) (I>=ML) are the
		boundary conditions at TN. BC is the Jacobian matrix dg_i/dY[K]
		Y1 and YN are arrays of length M specifying the
		solution at t=T1 and TN.

	Required functions : SETMAT, GAUBLK, EQN, BCS
*/

#include <math.h>
#include <stdlib.h>

int setmat(int n, int m, int ml, double *a, double *bc, double *x,
	double *xc, double t[], double par[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ));
int gaublk(int n, int m, int ml, double *a, int *iflg, double *det,
	int *idet, int *inc, double *x);

int fdm(int n, int m, int ml, double par[], double *x, double *xc, double t[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
	int iflag, double reps)

{
	int i,j,j1,k,m1,idet,ier,iflg, nit=20;
	double det,rerr,xj,r1,r2,re,h,hh,tj,xd, eps=1.0e-30;
	int *iwk;
	double *wk,*bc;

	if(n<3 || m<=ml || ml<=0) return 704;

	wk=(double *) calloc((size_t) ((m+ml)*2*m*n), sizeof(double));
	bc=(double *) calloc((size_t) ((m+ml)*2*m), sizeof(double));
	iwk=(int *) calloc((size_t) (m*n), sizeof(int));

	for(i=1; i<=nit; ++i) {
/*	Set up the matrix of finite difference equations */
		ier=setmat(n,m,ml,wk,bc,x,xc,t,par,eqn,bcs);
		iflg=0;
/*	Solve the system of linear equations */
		ier=gaublk(n,m,ml,wk,&iflg,&det,&idet,iwk,xc);
		if(ier>0) {free(iwk); free(bc); free(wk); return ier;}

/*	Checking for convergence */
		rerr=0.0;
		for(j=0; j<n; ++j) {
			j1=j+1; if(j==n-1) j1=j-1;
			for(k=0; k<m; ++k) {
				xj=x[k+j*m]+xc[k+j*m];
				r2=fabs(xj)+fabs(x[k+j*m]-x[k+j1*m])+eps;
				re=fabs(xc[k+j*m]/r2);
				if(re>rerr) rerr=re;
				x[k+j*m]=xj;
			}
		}
		if(rerr<reps || iflag>1) break;
	}

	if(i>nit) {free(iwk); free(bc); free(wk); return 737;}

	if(iflag==1 || iflag>2) {free(iwk); free(bc); free(wk); return 0;}

	if(n<=4) {free(iwk); free(bc); free(wk); return 734;}

/*	Calculate the deffered correction */
	hh=t[1]-t[0];
	m1=m+ml;
	for(j=0; j<n-1; ++j) {
		tj=0.5*(t[j]+t[j+1]);
		h=t[j+1]-t[j];
		if(fabs(h-hh)>1.e-4*fabs(hh)) {free(iwk); free(bc); free(wk); return 736;}

		eqn(j,m,ml,par,bc,&bc[m1*m],&x[j*m],&xc[(j+1)*m],tj);
		for(i=0; i<m; ++i) {
			if(j==0) {
				bc[m+i*m1]=(-2*x[i]+7*x[i+m]-9*x[i+2*m]+5*x[i+3*m]-x[i+4*m])/24.0;
				bc[m+(i+m)*m1]=(43*x[i]-112*x[i+m]+102*x[i+2*m]-40*x[i+3*m]+
								7*x[i+4*m])*h/192.0;
			}
			else if(j==n-2) {
				bc[m+i*m1]=(2*x[i+m*(n-1)]-7*x[i+m*(n-2)]+9*x[i+m*(n-3)]-
								5*x[i+m*(n-4)]+x[i+m*(n-5)])/24.0;
				bc[m+(i+m)*m1]=(43*x[i+(n-1)*m]-112*x[i+(n-2)*m]+
					102*x[i+(n-3)*m]-40*x[i+(n-4)*m]+7*x[i+(n-5)*m])*h/192.0;
			}
			else {
				bc[m+i*m1]=(-x[i+(j-1)*m]+3.*x[i+j*m]-3.*x[i+(j+1)*m]
							+x[i+(j+2)*m])/24.0;
				bc[m+(i+m)*m1]=(x[i+(j-1)*m]-x[i+j*m]-x[i+(j+1)*m]+
							x[i+(j+2)*m])*h/16.0;
			}
		}

/*	Set up the RHS for deferred correction */
		for(i=0; i<m; ++i) {
			xd=0.0;
			for(k=0; k<m; ++k) xd=xd+bc[i+(k+m)*m1]*bc[m+k*m1]-
				bc[i+k*m1]*bc[m+(k+m)*m1];
			xc[ml+i+j*m]=xd;
		}
	}

	for(i=0; i<ml; ++i) xc[i]=0.0;
	for(i=ml; i<m; ++i) xc[i+(n-1)*m]=0.0;

/*	Calculate deferred correction */
	ier=gaublk(n,m,ml,wk,&iflg,&det,&idet,iwk,xc);
	for(j=0; j<n; ++j) {
		for(i=0; i<m; ++i) xc[i+j*m]=xc[i+j*m]+x[i+j*m];
	}
	free(iwk); free(bc); free(wk);
	return 0;
}


/*	--------------------------------

void eqn(int j, int m, int ml, double par[], double *a, double *b,
	double y[], double f[], double t)

{
	int i,k,m1;

	m1=m+ml;
	for(i=0; i<m; ++i) {
		f[i]=f_i(t,par,y);
		for(k=0; k<m; ++k) {
			a[k+i*m1]=df_k/dy[i];
			b[k+i*m1]=b_{ki}(t,par);
		}
	}
	return;
}

void bcs(int m, int ml, double par[], double *bc, double g[],
	double t1, double tn, double y1[], double yn[])

{
	int i,k,m1;

	m1=m+ml;
	for(i=0; i<m; ++i) {
		if(i<ml) g[i]=g_i(t1,par,y1);
		else g[i]=g_i(tn,par,yn);

		for(k=0; k<m; ++k) bc[i+k*m1]=dg_i/dy[k];
	}
	return;
}

*/
