/*	To solve two point boundary value problem using finite difference method */

#include <stdio.h>
#include <math.h>

int setmat(int n, int m, int ml, double *a, double *bc, double *x,
	double *xc, double t[], double par[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ));
int gaublk(int n, int m, int ml, double *a, int *iflg, double *det,
	int *idet, int *inc, double *x);
int fdm(int n, int m, int ml, double par[], double *x, double *xc, double t[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
	int iflag, double reps);
void eqn(int j, int m, int ml, double par[], double *a, double *b,
	double x[], double f[], double t);
void bcs(int m, int ml, double par[], double *bc, double g[], double t1,
	double tn, double x1[], double xn[]);

double al;

main()
{
	int i,i1,j,n,m,ml,id, iflg, ier,np,nmax;
	double hh,t[500],x[500][2],xc[500][2],tn,par[20],reps,t0;

/*	Example 12.13 */

	m=2; ml=1; reps=1.e-8;

/*	pass the parameter lambda in the equation via array par */

	for(i1=0; i1<99; ++i1) {
		printf("type n = no. of points,  lambda     (quits when n<=0)\n");
		scanf(" %d %le",&n, &par[0]);
		if(n<=0) return 0;

/*	Set up the mesh and initial values */
		iflg=0;
		hh=1.0/(n-1);
		for(i=0; i<n; ++i) {
			t[i]=i*hh;
			x[i][0]=1.0;
			x[i][1]=1.0;
		}
		i=fdm(n,m,ml,par,&x[0][0],&xc[0][0],t,eqn,bcs,iflg,reps);
		printf(" ier = %d    no. of points = %d    lamdba = %e\n",i,n,par[0]);
		printf("       t           solution                    corrected solution \n");
		for(i=0; i<n; i+=10) printf(" %d %e   %e %e   %e %e \n",i,t[i],x[i][0],x[i][1],xc[i][0],xc[i][1]);

	}
	return;
}


/*	The differential equation */

void eqn(int j, int m, int ml, double par[], double *a, double *b,
	double y[], double f[], double t)

{
	int i,k,m1;

	m1=m+ml;
	for(i=0; i<m; ++i) {
		f[i]=0.0;
		for(k=0; k<m; ++k) {
			a[k+i*m1]=0.0;
			b[k+i*m1]=0.0;
		}
	}
	b[0]=1.0;
	a[m1]=1.0;
	b[1+m1]=1.0;
	a[1]=par[0]*par[0]*cosh(par[0]*y[0]);

	f[0]=y[1];
	f[1]=par[0]*sinh(par[0]*y[0]);
	return;
}

/*	The boundary conditions */

void bcs(int m, int ml, double par[], double *bc, double g[],
	double t1, double tn, double y1[], double yn[])

{
	int i,k,m1;

	m1=m+ml;
	bc[0]=1.0;
	bc[m1]=0.0;
	bc[1]=1.0;
	bc[1+m1]=0.0;

	g[0]=y1[0];
	g[1]=yn[0]-1.0;
	return;
}



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



/*	To solve a system of linear equations arising from finite difference
	approximation of ordinary differential equations

	N : (input) Number of mesh points used.
	M : (input) Number of first order differential equations in the system
	ML : (input) Number of boundary conditions at the first boundary t=T[0]
	A : (input/output) Array of length (M+ML)*2M*N containing
		the matrix of equations. After execution it will contain
		the triangular decomposition
	IFLG : (input/output) Integer variable used as a flag to decide
		the nature of computation.
		If IFLG=0 the triangular decomposition, determinant and
			solution of equations is computed. IFLG is set to 2
		If IFLG=1 only the triangular decomposition and determinant
			are calculated and IFLG is set to 2
		If IFLG=2 then it is assumed that triangular decomposition
			is already done and available in A and only the solution
			of system of equations is solved
	DET : (output) Scaled value of determinant of the matrix
	IDET : (output) exponent of determinant, the value of determinant
		is DET*2**IDET
	INC : (input/output) Integer array of length M*N containing the
		information about interchanges used by Gaussian elimination.
		It is calculated if IFLG=0 or 1, while for IFLG=2 it must
		be supplied from previous calculation.
	X : (input/output) Array of length M*N containing the right hand
		side of equation at input. After execution it will be overwritten
		by the solution if IFLG=0 or 2. For IFLG=1 it is not used.
		
	Error status is returned by the value of the function GAUBLK.
		0 value implies successful execution
		735 implies that the matrix is singular

	Required functions : None
*/

#include <math.h>

int gaublk(int n, int m, int ml, double *a, int *iflg, double *det,
	int *idet, int *inc, double *x)

{
	int i,j,k,mr,mc,kmax,ki,kj,kk,m1,m2;
	double rmax,r1,at,xt,d1;

	m1=m+ml;
	m2=m1*2*m;

	if(*iflg<=1) {
		*idet=0; *det=1.0;
		mr=m1;		/*	The number of rows in each block of the matrix */
		mc=2*m;		/*	The number of columns in each block of the matrix*/

		for(j=0; j<n; ++j) {
			if(j==n-1) {mr=m; mc=m;}

			kk=mr-1; if(m<kk) kk=m;
			for(k=0; k<kk; ++k) {
				rmax=fabs(a[k+k*m1+j*m2]);
				kmax=k;
/*		Find the pivot */
				for(ki=k+1; ki<mr; ++ki) {
					r1=fabs(a[ki+k*m1+j*m2]);
					if(r1>rmax) {rmax=r1; kmax=ki;}
				}
				inc[k+j*m]=kmax;

				if(kmax != k) {
/*	exchange rows K and KMAX */
					*det= -(*det);
					for(ki=k; ki<mc; ++ki) {
						at=a[k+ki*m1+j*m2];
						a[k+ki*m1+j*m2]=a[kmax+ki*m1+j*m2];
						a[kmax+ki*m1+j*m2]=at;
					}
				}

				*det=(*det)*a[k+k*m1+j*m2];
/*		If the pivot is zero, then quit */
				if(a[k+k*m1+j*m2]==0.0) return 735;

/*	Gaussian elimination */
				for(ki=k+1; ki<mr; ++ki) {
					a[ki+k*m1+j*m2]=a[ki+k*m1+j*m2]/a[k+k*m1+j*m2];
					for(kj=k+1; kj<mc; ++kj) a[ki+kj*m1+j*m2] =
						a[ki+kj*m1+j*m2]-a[ki+k*m1+j*m2]*a[k+kj*m1+j*m2];
				}
			}

			if(*det != 0.0) {
/*	Scale the determinant if necessary */
				while(fabs(*det)>32.0) {
					*det=(*det)*0.03125;
					*idet=(*idet)+5;
				}

				while(fabs(*det)<0.03125) {
					*det=(*det)*32.0;
					*idet=(*idet)-5;
				}
			}

/*	Copy the overlapping elements into the next block */
			if(j<n-1) {
				for(k=0; k<ml; ++k) {
					for(ki=0; ki<m; ++ki) a[k+ki*m1+(j+1)*m2]=a[k+m+(ki+m)*m1+j*m2];
				}
			}
		}

		inc[n*m-1]=m-1;
		*det=(*det)*a[m-1+(m-1)*m1+(n-1)*m2];
		if(a[m-1+(m-1)*m1+(n-1)*m2]==0.0) return 735;

		if(*iflg==1) {*iflg=2; return 0;}
		*iflg=2;
	}

/*	Solve the system of linear equations */
	mr=m1;
	for(j=0; j<n; ++j) {
		if(j==n-1) mr=m;
		kj=mr-1; if(m<kj) kj=m;
		for(k=0; k<kj; ++k) {
			kk=inc[k+j*m];
			if(k != kk) {
/*	exchange the corresponding elements of RHS */
				xt=x[k+j*m];
				x[k+j*m]=x[kk+j*m];
				x[kk+j*m]=xt;
			}
			
/*	Gaussian elimination */
			for(ki=k+1; ki<mr; ++ki) x[ki+j*m]=x[ki+j*m]-a[ki+k*m1+j*m2]*x[k+j*m];
		}
	}

/*	back-substitution */
	mc=m;
	for(j=n-1; j>=0; --j) {
		for(k=m-1; k>=0; --k) {
			d1=x[k+j*m];
			for(ki=mc-1; ki>=k+1; --ki) d1=d1-x[ki+j*m]*a[k+ki*m1+j*m2];
			x[k+j*m]=d1/a[k+k*m1+j*m2];
		}
		mc=2*m;
	}
	return 0;
}



/*	To setup the matrix of system of linear equations arising from
	finite difference approximation of ordinary differential equations
	This function is called by FDM or GEVP

	N : (input) Number of mesh points used.
	M : (input) Number of first order differential equations in the system
	ML : (input) Number of boundary conditions at the first boundary t=T[0]
	A : (output) Array of length (M+ML)*2M*N containing the matrix of equations.
	BC : (output) Array of length (M+ML)*(M+1) containing the
		coefficients of boundary conditions
	X : (input) Array of length M*N containing the current approximation
		to the solution.
	XC : (output) Array of length M*N containing the right hand
		side of finite difference equations calculated by the function
	T : (input) Array of length N containing the mesh points.
	PAR : (input) Array containing the parameters to be passed
		on to functions EQN and BCS. This array is not used by
		the function.
	EQN : (input) Name of function to calculate the equation matrix
	BCS : (input) Name of function to calculate the boundary conditions

	Functions EQN(J,M,ML,PAR,A,B,Y,F,T) and
	BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) must be supplied by the user.
	The form of these functions is described in documentation for FDM or GEVP

	The returned value is always zero.

	Required functions : EQN, BCS
*/

#include <math.h>

int setmat(int n, int m, int ml, double *a, double *bc, double *x,
	double *xc, double t[], double par[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ))

{
	int i,j,k,m1,m2;
	double tk,h,xk;

	m1=m+ml;
	m2=m1*2*m;
/*	Loop over the mesh points */
	for(k=0; k<n-1; ++k) {
		tk=0.5*(t[k]+t[k+1]);	/*	t_{k+1/2} */
		h=t[k+1]-t[k];
		for(i=0; i<m; ++i) bc[i]=0.5*(x[i+k*m]+x[i+(k+1)*m]);	/* y_{k+1/2} */

/*	Calculate the equation matrix at t_{k+1/2} */
		eqn(k,m,ml,par,&a[k*m2],&a[m*m1+k*m2],bc,&xc[ml+k*m],tk);

/*	Setup the RHS of finite difference equations */
		for(j=0; j<m; ++j) {
			xk=xc[ml+j+k*m]*h;
			for(i=0; i<m; ++i) xk=xk-a[j+(m+i)*m1+k*m2]*(x[i+(k+1)*m]-x[i+k*m]);
			xc[ml+j+k*m]=xk;
		}

/*	Setup the finite difference matrix */
		for(j=0; j<m; ++j) {
			for(i=m-1; i>=0; --i) {
				a[ml+i+j*m1+k*m2]= -a[i+(j+m)*m1+k*m2]-0.5*h*a[i+j*m1+k*m2];
				a[ml+i+(j+m)*m1+k*m2]= a[i+(j+m)*m1+k*m2]-0.5*h*a[i+j*m1+k*m2];
			}
			for(i=0; i<ml; ++i) a[i+(j+m)*m1+k*m2]=0.0;
		}
	}

/*	The boundary conditions */
	bcs(m,ml,par,bc,&bc[m*m1],t[0],t[n-1],x,&x[(n-1)*m]);

/*	Boundary conditions at the first boundary */
	for(i=0; i<ml; ++i) {
		xc[i]= -bc[i+m*m1];
		for(j=0; j<m; ++j) a[i+j*m1]=bc[i+j*m1];
	}

/*	Boundary conditions at the second boundary */
	for(i=ml; i<m; ++i) {
		xc[i+(n-1)*m]= -bc[i+m*m1];
		for(j=0; j<m; ++j) a[i+j*m1+(n-1)*m2]=bc[i+j*m1];
	}
	return 0;
}
