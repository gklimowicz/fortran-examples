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
