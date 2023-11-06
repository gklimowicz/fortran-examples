/*	To solve eigenvalue problem in ordinary differential equations */

#include <stdio.h>
#include <math.h>

int setmat(int n, int m, int ml, double *a, double *bc, double *x,
	double *xc, double t[], double par[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ));
int gaublk(int n, int m, int ml, double *a, int *iflg, double *det,
	int *idet, int *inc, double *x);
void secani(double x0, double xl, double xu, double *x, double f, int jf,
	double reps,double aeps, int *ier);
void eqn(int j, int m, int ml, double par[], double *a, double *b,
	double x[], double f[], double t);
void bcs(int m, int ml, double par[], double *bc, double g[], double t1,
	double tn, double x1[], double xn[]);
void eqnd(int j, int m, int ml, double par[], double *a, double *b, double t);
void bcsd(int m, int ml, double par[], double *bc, double t1, double tn);

int gevp(int n, int m, int ml, double par[], double *x, double *xc,
	double t[], double *e0,
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
	void eqnd(int , int , int , double * , double * , double * , double ),
	void bcsd(int , int , double * , double * ,  double , double ),
	int iflag, double reps, double el, double eu);


main()
{
	int i,i1,j,n,m,ml,id, iflg, ier,np,nmax;
	double hh,t[500],x[500][2],xc[500][2],tn,par[20],reps,e0,el,eu;

/*	Example 12.14 : Spheroidal harmonics */

	m=2; ml=1; reps=1.e-8;
/*	Parameters C**2 and M */
	par[1]=1.0; par[2]=2.0;

	for(i1=0; i1<99; ++i1) {
		printf("type n = no. of points,   e0 = initial guess for eigenvalue\n");
		printf("                    (quits when n<=0)\n");
		scanf(" %d %le",&n, &e0);
		if(n<=0) return 0;
		el=e0-2; eu=e0+2;

		iflg=2;
		hh=1.0/(n-1);
		for(i=0; i<n; ++i) {
			t[i]=i*hh;
			x[i][0]=1.0;
			x[i][1]=1.0;
		}
		i=gevp(n,m,ml,par,&x[0][0],&xc[0][0],t,&e0,eqn,bcs,eqnd,
			bcsd,iflg,reps,el,eu);
		printf(" ier = %d    no. of points = %d\n  eigenvalue = %e    corrected eigenvalue = %e\n",i,n,par[0],e0);
		printf("     T             Eigenfunction  \n");
		for(i=0; i<n; i+=10) printf(" %d %e   %e %e \n",i,t[i],x[i][0],x[i][1]);

	}
	return;
}

/*	The differential equations */

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
	b[1+m1]=1.0-t*t;
	a[1]=-(par[0]-par[1]*t*t);
	a[1+m1]=2.*(par[2]+1.0)*t;

	return;
}

/*	The boundary conditions */

void bcs(int m, int ml, double par[], double *bc, double g[],
	double t1, double tn, double y1[], double yn[])

{
	int i,k,m1;

	m1=m+ml;
	bc[0]=0.0;
	bc[m1]=1.0;
	bc[1]=-0.5*(par[0]-par[1])/(par[2]+1.0);
	bc[1+m1]=1.0;

	g[0]=0.0;
	g[1]=0.0;
	return;
}

/*	Derivative of equation matrix w.r.t. eigenvalue */

void eqnd(int j, int m, int ml, double par[], double *a, double *b, double t)

{
	int i,k,m1;

	m1=m+ml;
	for(i=0; i<m; ++i) {
		for(k=0; k<m; ++k) {
			a[k+i*m1]=0.0;
			b[k+i*m1]=0.0;
		}
	}
	a[1]=-1.0;

	return;
}

/*	Derivative of boundary condition matrix w.r.t. eigenvalue */

void bcsd(int m, int ml, double par[], double *bc, double t1, double tn)

{
	int i,k,m1;

	m1=m+ml;
	bc[0]=0.0;
	bc[m1]=0.0;
	bc[1]=-0.5/(par[2]+1.0);
	bc[1+m1]=0.0;

	return;
}



/*	To solve a generalised eigenvalue problem in ordinary differential
	equations using finite difference method

	N : (input) Number of mesh points to be used.
	M : (input) Number of first order differential equations in the system
	ML : (input) Number of boundary conditions at the first boundary t=T[0]
	PAR : (input/output) Array to be passed on to EQN, EQND, BCS and BCSD
		for specifying the equations and boundary conditions. PAR[0]
		is used for passing the eigenvalue and hence should not be
		used for any other variable. After execution the eigenvalue
		will be available in PAR[0]. Other elements of this array are
		not used by GEVP, but are merely used to pass on any extra parameters
		that may be required to specify the equations
	X : (output) Array of length M*N containing the eigenvector.
		X[j][i] is the ith component of solution at jth mesh point.
		Second dimension of X in calling function must be M
	XC : (output) Array of length M*N containing the left eigenvector.
		It is stored in the same format as X.
	T : (input) Array of length N containing the mesh points.
		These points must be in ascending or descending order.
		For calculating the deferred correction the mesh spacing
		must be uniform.
	E0 : (output) The calculated eigenvalue including deferred correction
		PAR[0] contains the uncorrected eigenvalue. If deferred
		correction is not applied E0 is not calculated.
	EQN : (input) Name of function to specify the differential equation
		to be solved.
	BCS : (input) Name of function to calculate the boundary conditions
		at t=T[0] and T[N-1]
	EQND : (input) Name of function to calculate derivative of equation
		matrix with respect to the eigenvalue.
	BCSD : (input) Name of function to calculate derivative of the boundary
		conditions with respect to the eigenvalue
		A sample outline for functions EQN, BCS, EQND, BCSD is
		included at the end of this file.
	IFLAG : (input) Integer variable used as a flag to decide the type of
		computation required.
		IFLAG=0 implies that only the eigenvalue is calculated.
		IFLAG=1 implies that both eigenvalue and eigenvector are
			calculated.
		IFLAG=2 implies that first order correction to eigenvalue
			is also calculated in addition to eigenvector
			In this case the mesh spacing must be uniform
	REPS : (input) Required accuracy. This is only used to check convergence
		of Secant method for finding zeros of the determinant
		The truncation error depends on mesh spacing and is not
		controlled in this function.
	EL : (input) Lower limit on the eigenvalue.
	EU : (input) Upper limit on the eigenvalue.
		These parameters are passed on to SECANI and iteration is
		terminated if it goes outside these limits.
		
	Error status is returned by the value of the function GEVP.
		0 value implies successful execution
		704 implies that N<3, M<=ML, or ML<=0, in which case
			no calculations are done
		734 implies that N<5 and deferred correction is requested.
			In this case deferred correction is not calculated
		735 implies that the finite difference matrix is singular
			while calculating the eigenvector. In this case
			eigenvector is not calculated, but eigenvalue is
			available in PAR[0].
		736 implies that mesh spacing is not uniform and
			deferred correction is not calculated
		738 implies that eigenvector vanishes.
		739 implies that inverse iteration for calculating the
			eigenvector failed to converge.
		740 implies that inverse iteration for calculating the
			left eigenvector failed to converge.

	Functions EQN, EQND, BCS and BCSD must be supplied by the user.
	Function EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
		sides for differential equations B y'=A y
		J is the serial number of mesh point at which calculation
		is required. M is the number of first order differential
		equations in the system, ML the number of boundary conditions
		at the first point, PAR is an array which can be used
		to pass on any required parameters to define the equations.
		PAR[0] is the eigenvalue.
		A and B are arrays of length (M+ML)*M defining the
		differential equation B Y'= A Y.
		Y is an array of length M specifying the solution at t=T
		F is an array of length M which should be set to zero.
		F, A and B must be calculated by the function. T is the
		value of t at which solution Y is specified.

	Function BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) calculates the boundary
		conditions at both boundaries. M is the number of differential
		equations, ML is the number of boundary condition at the
		first mesh point t=T1. PAR is an array which can be used
		to pass any required parameters. PAR[0] is the eigenvalue.
		BC is an array of length (M+ML)*M which should contain
		the coefficients of boundary conditions, BC Y =0.
		First ML rows will specify the boundary conditions
		at t=T1, while remaining rows will specify those at t=TN.
		G is an array of length M which should be set to zero.
		Y1 and YN are arrays of length M specifying the
		solution at t=T1 and TN.

	Function EQND(J,M,ML,PAR,A,B,T) calculates the derivatives of
		matrices A and B (calculated by EQN) with respect to the
		eigenvalue (PAR[0]). 
		J is the serial number of mesh point at which calculation
		is required. M is the number of first order differential
		equations in the system, ML the number of boundary conditions
		at the first point, PAR is an array which can be used
		to pass on any required parameters to define the equations.
		PAR[0] is the eigenvalue.
		A and B are arrays of length (M+ML)*M which should
		give the derivatives of arrays A and B calculated by EQN.
		T is the value of t at which the derivatives are required.

	Function BCSD(M,ML,PAR,BC,T1,TN) calculates the derivative of
		matrix BC (calculated by BCS) with respect to the eigenvalue PAR[0].
		M is the number of differential equations,
		ML is the number of boundary condition at the
		first mesh point t=T1. PAR is an array which can be used
		to pass any required parameters. PAR[0] is the eigenvalue.
		BC is an array of length (M+ML)*M which should contain
		the derivative of the matrix BC calculated by BCS.
		First ML rows will specify the boundary conditions
		at t=T1, while remaining rows will specify those at t=TN.

	Required functions : SETMAT, GAUBLK, SECANI (or MULER2), EQN, EQND, BCS, BCSD
*/

#include <math.h>
#include <stdlib.h>

int setmat(int n, int m, int ml, double *a, double *bc, double *x,
	double *xc, double t[], double par[],
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ));
int gaublk(int n, int m, int ml, double *a, int *iflg, double *det,
	int *idet, int *inc, double *x);
void secani(double x0, double xl, double xu, double *x, double f, int jf,
	double reps,double aeps, int *ier);


int gevp(int n, int m, int ml, double par[], double *x, double *xc,
	double t[], double *e0,
	void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
	void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
	void eqnd(int , int , int , double * , double * , double * , double ),
	void bcsd(int , int , double * , double * ,  double , double ),
	int iflag, double reps, double el, double eu)

{
	int i,j,k,ier,idet,iflg,it,i1,k0,k1,k2,jj,m0,m1,m2,mr,kk, nit=20;
	double raps,cx,det,rmax,xmax,xt,hh,h,enu,eden,tj,xn,xd;
	double eps=1.0e-30;
	int *iwk;
	double *wk,*bc;

	if(n<3 || m<=ml || ml<=0) return 704;

	m1=m+ml; m2=m1*2*m;
	wk=(double *) calloc((size_t) (m2*n), sizeof(double));
	bc=(double *) calloc((size_t) m2, sizeof(double));
	iwk=(int *) calloc((size_t) (m*n), sizeof(int));
	raps=0.1*reps;
	ier=0;
	do {
		secani(*e0,el,eu,&cx,det,idet,reps,raps,&ier);
		if(ier>=0) break;
/*	calculate the determinant */
		par[0]=cx;
		i=setmat(n,m,ml,wk,bc,x,xc,t,par,eqn,bcs);
		iflg=1;
		i=gaublk(n,m,ml,wk,&iflg,&det,&idet,iwk,xc);
	} while(ier<0);

/*	the eigenvalue */
	par[0]=cx;
	if(ier>100 || iflag==0) {free(iwk); free(bc); free(wk); return ier;}
	if(i>0) {free(iwk); free(bc); free(wk); return 735;}

/*	Initialise the vector for inverse iteration */
	for(i=0; i<n; ++i) {
		for(j=0; j<m; ++j) {x[j+i*m]=1.0; xc[j+i*m]=1.0;}
	}
	iflg=2;

/*	Loop for inverse iteration to calculate the eigenvector */
	for(i=1; i<=nit; ++i) {
		ier=gaublk(n,m,ml,wk,&iflg,&det,&idet,iwk,x);

/*	Normalising the eigenvector */
		rmax=0.0;
		for(j=0; j<n; ++j) {
			for(k=0; k<m; ++k) {
				if(fabs(x[k+j*m])>rmax) {
					rmax=fabs(x[k+j*m]);
					xmax=x[k+j*m];
				}
			}
		}
		if(rmax==0.0) {free(iwk); free(bc); free(wk); return 738;}

		for(j=0; j<n; ++j) {
			for(k=0; k<m; ++k) x[k+j*m]=x[k+j*m]/xmax;
		}

/*	convergence check */
		if(rmax*reps>1.0) break;
	}
	if(i>nit) {free(iwk); free(bc); free(wk); return 739;}

	if(iflag<=1) {free(iwk); free(bc); free(wk); return 0;}
	if(n<=4) {free(iwk); free(bc); free(wk); return 734;}

/*	Loop for inverse iteration for the left eigenvector */
	for(it=1; it<=nit; ++it) {

/*	Forward substitution for U^T */
		for(i=0; i<n; ++i) {
			for(j=0; j<m; ++j) {
				if(i>0) {
					for(k=0; k<m; ++k)
						xc[j+i*m]=xc[j+i*m]-xc[k+(i-1)*m]*wk[k+(j+m)*m1+(i-1)*m2];
				}
				for(k=0; k<=j-1; ++k) xc[j+i*m]=xc[j+i*m]-xc[k+i*m]*wk[k+j*m1+i*m2];
				xc[j+i*m]=xc[j+i*m]/wk[j+j*m1+i*m2];
			}
		}

/*	Back substitution for L^T */
		mr=m;
		for(i=n-1; i>=0; --i) {
			if(i==n-2) mr=m1;
			m0=mr-1; if(m<m0) m0=m;
			for(j=m0-1; j>=0; --j) {
				for(k=j+1; k<mr; ++k) {
					k0=k; i1=i;
/*	Checking for backward interchange */
					for(k1=j+1; k1<=k; ++k1) {
						if(k1>m-1) k0=k-m;
						if(iwk[k1+i*m] == k0) goto endk;
					}

					k1=k;
					while(i1<n) {
						if(k1>=m) {k1=k1-m; i1=i1+1;}
						if(iwk[k1+i1*m]>k1) {
							k0=iwk[k1+i1*m];
							k2=k0;
							for(kk=k1+1; kk<=k2; ++kk) {
								if(kk>=m) k0=k2-m;
								if(iwk[kk+i1*m]==k0) {
									k1=kk;
									goto endk;
								}
							}
							k1=k2;
						}
					}
/*	Hopefully, this statement will not be executed */
					exit(1);

endk:				xc[j+i*m]=xc[j+i*m]-xc[k1+i1*m]*wk[k+j*m1+i*m2];
				}
			}
		}
					        
/*	Exchanges due to interchange matrix N */
		for(i=n-1; i>=0; --i) {
			for(j=m-1; j>=0; --j) {
				jj=iwk[j+i*m];
				if(jj != j) {
					xt=xc[j+i*m];
					xc[j+i*m]=xc[jj+i*m];
					xc[jj+i*m]=xt;
				}
			}
		}

		rmax=0.0;
		for(i=0; i<n; ++i) {
			for(j=0; j<m; ++j) {
				if(rmax<fabs(xc[j+i*m])) rmax=fabs(xc[j+i*m]);
			}
		}
		for(i=0; i<n; ++i) {
			for(j=0; j<m; ++j) xc[j+i*m]=xc[j+i*m]/rmax;
		}

		if(rmax*reps>1.0) break;
	}
	if(it>nit) {free(iwk); free(bc); free(wk); return 740;}

/*	calculating the deferred correction */
	hh=t[1]-t[0];
	enu=0.0; eden=0.0;
	for(j=0; j<n-1; ++j) {
		tj=0.5*(t[j]+t[j+1]);
		h=t[j+1]-t[j];
		if(fabs(h-hh)>1.e-4*fabs(hh)) {free(iwk); free(bc); free(wk); return 736;}

		eqn(j,m,ml,par,&wk[m2],&wk[m2+m*m1],&x[j*m],&wk[2*m2],tj);
		eqnd(j,m,ml,par,&wk[2*m2],&wk[2*m2+m*m1],tj);

		for(i=0; i<m; ++i) {
			if(j==0) {
				wk[i]=(-2*x[i]+7*x[i+m]-9*x[i+2*m]+5*x[i+3*m]-x[i+4*m])/24.0;
				wk[i+m1]=(43*x[i]-112*x[i+m]+102*x[i+2*m]-40*x[i+3*m]+
						7*x[i+4*m])*h/192.0;
			}
			else if(j==n-2) {
				wk[i]=(2*x[i+(n-1)*m]-7*x[i+(n-2)*m]+9*x[i+(n-3)*m]-
						5*x[i+(n-4)*m]+x[i+(n-5)*m])/24.0;
				wk[i+m1]=(43*x[i+(n-1)*m]-112*x[i+(n-2)*m]+102*x[i+(n-3)*m]
						-40*x[i+(n-4)*m]+7*x[i+(n-5)*m])*h/192.0;
			}
			else {
				wk[i]=(-x[i+(j-1)*m]+3*x[i+j*m]-3*x[i+(j+1)*m]+x[i+(j+2)*m])/24.0;
				wk[i+m1]=(x[i+(j-1)*m]-x[i+j*m]-x[i+(j+1)*m]+x[i+(j+2)*m])*h/16.0;
			}
		}

		for(i=0; i<m; ++i) {
			xn=0.0; xd=0.0;
			for(k=0; k<m; ++k) {
				xd=xd+wk[i+(k+m)*m1+2*m2]*(x[k+(j+1)*m]-x[k+j*m])-
					0.5*h*wk[i+k*m1+2*m2]*(x[k+j*m]+x[k+(j+1)*m]);
				xn=xn+wk[i+(k+m)*m1+m2]*wk[k]-wk[i+k*m1+m2]*wk[k+m1];
			}
			enu=enu+xc[ml+i+j*m]*xn;
			eden=eden+xc[ml+i+j*m]*xd;
		}
	}

	bcsd(m,ml,par,wk,t[0],t[n-1]);
	for(i=0; i<ml; ++i) {
		xd=0.0;
		for(k=0; k<m; ++k) xd=xd+wk[i+k*m1]*x[k];
		eden=eden+xd*xc[i];
	}
	for(i=ml; i<m; ++i) {
		xd=0.0;
		for(k=0; k<m; ++k) xd=xd+wk[i+k*m1]*x[k+(n-1)*m];
		eden=eden+xd*xc[i+(n-1)*m];
	}

	*e0=par[0]+enu/eden;
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
		f[i]=0.0;
		for(k=0; k<m; ++k) {
			a[k+i*m1]=a_{ki}(t,par);
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
		g[i]=0.0;

		for(k=0; k<m; ++k) {
			if(i<ml) bc[i+k*m1]=bc_{ik}(t1,par);
			else bc[i+k*m1]=bc_{ik}(tn,par);
		}
	}
	return;
}

void eqnd(int j, int m, int ml, double par[], double *a, double *b, double t)

{
	int i,k,m1;

	m1=m+ml;
	for(i=0; i<m; ++i) {
		for(k=0; k<m; ++k) {
			a[k+i*m1]=da_{ki}/dpar[0] (t,par);
			b[k+i*m1]=db_{ki}/dpar[0] (t,par);
		}
	}
	return;
}

void bcsd(int m, int ml, double par[], double *bc, double t1, double tn)

{
	int i,k,m1;

	m1=m+ml;
	for(i=0; i<m; ++i) {

		for(k=0; k<m; ++k) {
			if(i<ml) bc[i+k*m1]=d bc_{ik}/dpar[0] (t1,par);
			else bc[i+k*m1]=d bc_{ik}/dpar[0] (tn,par);
		}
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



/*      Real zero of a given function using secant iteration
	Function is calculated as FX*2**JF
	This function uses reverse communication to calculate function
	values. If IER<0 the function should be evaluated and SECANI should
	be called back with new function value. Calculation of function
	value should not change any other variables in the call statement.

	X0 : (input) Initial guess for the zero
	XL : (input) Lower limit of interval where zero is expected
	XU : (input) Upper limit of interval where zero is expected
	X : (output) Value of x at which the function evaluation is required.
		If IER=0 then it will contain the final value of zero computed
		by the function.
	F : (input) Calculated value of the function at X.
		If function exits with IER<0, then the calling function should
		calculate the function value at X and call SECANI with this value
		stored in F and JF. Other variables should not be changed.
	JF : (input) The exponent of function value, the function value
		should be F*2**JF
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
     		The estimated error should be less than max(AEPS, REPS*fabs(X))
	IER : (input/output) Error parameter, IER=0 implies successful execution
		Before the first call IER should be set to zero
		IER<0 implies that execution is not over and the function needs
			a new function evaluation at X. After calculating the
			function value SECANI should be called back.
     		IER=40 implies that function value is equal at two points
     			and it is not possible to continue the iteration
     		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
     		IER=422 implies that iteration goes outside the specified limits
     		IER=423 implies that iteration failed to converge to specified accuracy


	Required functions : None (Function value is calculated by the calling program)
*/

#include <math.h>


void secani(double x0, double xl, double xu, double *x, double f, int jf,
	double reps,double aeps, int *ier)

{
	int nis=75;
	static int l,jf1;
	static double f1,r1,dx,dx1;

	if(*ier==-1) goto func;

	if(xl>x0 || xu<x0) {*ier=402; return;}

	*x=x0;
/*	Select the increment for the next point X+DX */
	dx=(xu-x0)/100.0;
	r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
	if(fabs(dx)<5.0*r1) dx=(xl-x0)/5.;
	r1=fabs(*x); if(r1<100.0*aeps) r1=100.0*aeps;
	if(fabs(dx)>0.1*r1) {
		if(dx>=0.0) dx=r1;
		else dx=-r1;
	}
	f1=0.0; jf1=0; l=0;

loop:	l=l+1; *ier=-1;
/*	To evaluate the function at x */
		return;

func:	dx1=dx;
		f1=f1*pow(2.0, (double) (jf1-jf));

		if(f1-f == 0.0) {
			if(f == 0.0) {*ier = 0; return ;}
/*	If F1=F and F!=0, then quit */
			else {*ier=40; return;}
		}

/*	The secant iteration */
		if(l>1) dx=dx1*f/(f1-f);
		*x=(*x)+dx;
		f1=f; jf1=jf;

		r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
		if(fabs(dx)<r1 && l>2) {*ier=0; return;}
		if(*x<xl || (*x)>xu) {*ier=422; return;}
	if(l<nis) goto loop;

	*ier=423; return;
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
