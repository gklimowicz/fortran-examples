/*	To solve a system of hyperbolic differential equation using the
	Lax-Wendroff difference scheme
	The differential equations are assumed to be of the form

	du_i/dt + df_i/dx =0,   where f_i=f_i(x,t,u)

	while boundary conditions are either Dirichlet form or

	A_i1(t) f_i+A_i2(t) df_i/dx = A_i3(t)

	N : (input) Number of equations in the system
	T : (input/output) Initial value of "time" where the initial conditions
		are specified. After execution, it will be replaced by the
		value of T at the last point where execution is successful
	DT : (input) The time step to be used for computations. It is kept fixed.
	X0 : (input) Lower limit on X where the solution is required
	XN : (input) Upper limit on X where the solution is required.
		Solution is computed in the interval (X0,XN)
	NT : (input) Number of time steps each of length DT  to be executed.
	NX : (input) Number of mesh points in the X direction.
	X : (output) Array of length NX containing the mesh points used
		in X direction. These are calculated by assuming uniform spacing.
	U : (input/output) Array of length IU*NX containing the solution
		at the current time step. It should contain the initial
		values at the time of calling, if IFLG!=0. Otherwise it
		is computed using function FIC. After execution it will contain
		the computed solution at t=T. U[j][i] is the ith component
		at jth mesh point in X.
	IU : (input) The second dimension of U as declared in the calling
		function, IU>=N
	FLUX : (input) Name of the function to calculate the "fluxes"
		f_i(x,t,u)
	BC : (input) Name of the function to calculate the coefficients
		in the boundary conditions
	FIC : (input) Name of the function to calculate the initial
		values when IFLG=0. For other values of IFLG this function
		is not used, but a dummy function may be required by the compiler.
	IFLG : (input/output) Integer variable used as a flag to denote how the
			initial values are calculated.
		If IFLG=0 then the initial values are calculated using the
			function FIC to be supplied by the user. IFLG is set to 1
		Otherwise the initial values must be supplied in array U.
		
	Error status is returned by the value of the function LAX.
		0 value implies successful execution
		715 implies that DT=0, XN=X0, NX<3, or IU<N
			in which case no calculations are done
		763 implies that the difference equations are singular
			and solution cannot be continued further.

	Functions FLUX(N,X,T,U,F), BC(IB,N,X,T,IW,A)
	and FIC(N,X,T,U) must be supplied by the user 
	Function FLUX should calculate f_i(x,t,u) as
	defined above for given values of X,T,U. Here X,T are the values
	of x and t where the calculations are required, N is the number of
	equations, while U and F are arrays 
	of length N containing the solution U and corresponding fluxes F.
	The function should calculate F[I] using the values of X,T,U.

	Function BC should calculate the coefficients for the boundary
	conditions for given values of X, T. Here IB is an integer variable
	used as a flag to denote the boundary at which the boundary conditions
	are required. IB=1 for boundary at x=X0 and IB=2 for x=XN.
	N is the number of equations in the system, X and T
	are the values of x and t where the boundary conditions are required.
	IW is an integer array which is used as a flag to transmit the
	type of boundary conditions. This flag also must be set by function
	BC. IW[I] should be set to zero if the boundary conditions are of
	Dirichlet form, any other value will imply boundary condition on
	flux as described above. A is an array of length N*3 which
	contains the coefficients required to define the boundary conditions.
	This array must be dimensioned A[3][N] in the function. For
	Dirichlet conditions the boundary values are given by
	A[0][I]U[I]=A[2][I], otherwise the boundary condition is assumed as
	A[0][I]F[I]+A[1][I]dF[I]/dx=A[2][I]
	The coefficient A[J][I] and IW[I] for each I, must be set by function BC.

	Function FIC is required only if IFLG=0, otherwise a dummy function
	with this name will suffice. If IFLG=0, function FIC must calculate
	the initial values at required X,T. Here N is the number of equations
	and U is an array of length N, which should contain the calculated initial
	values of all components.

	Required functions : FLUX, BC, FIC
*/

#include <math.h>
#include <stdlib.h>

int lax(int n, double *t, double dt, double x0, double xn, int nt, int nx,
	double x[], double *u, int iu,
	void flux(int , double , double , double * , double * ),
	void bc(int , int , double , double , int *, double * ),
	void fic(int , double , double , double * ), int *iflg)

{
	int i,j,k,j1,j2,jt;
	double dx,s,t1,tt,d1,f1;
	int *iwk;
	double *wk,*wb;

	if(dt==0.0 || xn==x0 || nx<=2 || iu<n) return 715;

	dx=(xn-x0)/(nx-1);
	s=dt/dx;
/*	Setting the initial values */
	for(j=0; j<nx; ++j) {
		x[j]=x0+j*dx;
		if(*iflg==0) fic(n,x[j],*t,&u[j*iu]);
	}
	*iflg=1;
	wk=(double *) calloc((size_t) (n*(nx+2)), sizeof(double));
	wb=(double *) calloc((size_t) (n*3), sizeof(double));
	iwk=(int *) calloc((size_t) n, sizeof(int));

/*	Loop over time steps */
	for(i=1; i<=nt; ++i) {
		t1=(*t)+dt/2.0; tt=(*t)+dt;

/*	First half-step of Lax-Wendroff scheme */
		j1=nx; j2=nx+1;
		flux(n,x[0],*t,u,&wk[j1*n]);
		for(j=1; j<nx; ++j) {
			flux(n,x[j],*t,&u[j*iu],&wk[j2*n]);
			for(k=0; k<n; ++k) wk[k+j*n]=0.5*(u[k+j*iu]+u[k+(j-1)*iu])
								-0.5*s*(wk[k+j2*n]-wk[k+j1*n]);
			jt=j1; j1=j2; j2=jt;
		}

/*	The boundary conditions at X0 */
		j1=nx; j2=nx+1;
		flux(n,x[1],t1,&wk[n],&wk[n*j1]);
		bc(1,n,x0,tt,iwk,wb);
		for(k=0; k<n; ++k) {
			if(iwk[k]==0) {
				if(wb[k]==0.0) {free(iwk); free(wb); free(wk); return 763;}
				u[k]=wb[k+2*n]/wb[k];
			}
			else {
				d1=0.5*wb[k]-wb[k+n]/dx;
				if(d1==0.0) {free(iwk); free(wb); free(wk); return 763;}
				f1=(wb[k+2*n]-wk[k+j1*n]*(0.5*wb[k]+wb[k+n]/dx))/d1;
				u[k]=u[k]-s*(wk[k+j1*n]-f1);
			}
		}

/*	Second half-step of Lax-Wendroff scheme */
		for(j=1; j<nx-1; ++j) {
			flux(n,x[j+1],t1,&wk[(j+1)*n],&wk[j2*n]);
			for(k=0; k<n; ++k) u[k+j*iu]=u[k+j*iu]-s*(wk[k+j2*n]-wk[k+j1*n]);
			jt=j1; j1=j2; j2=jt;
		}

/*	Boundary conditions at XN */
		bc(2,n,xn,tt,iwk,wb);
		for(k=0; k<n; ++k) {
			if(iwk[k]==0) {
				if(wb[k]==0.0) {free(iwk); free(wb); free(wk); return 763;}
				u[k+(nx-1)*iu]=wb[k+2*n]/wb[k];
			}
			else {
				d1=0.5*wb[k]+wb[k+n]/dx;
				if(d1==0.0) {free(iwk); free(wb); free(wk); return 763;}
				f1=(wb[k+2*n]-wk[k+j1*n]*(0.5*wb[k]-wb[k+n]/dx))/d1;
				u[k+(nx-1)*iu]=u[k+(nx-1)*iu]-s*(f1-wk[k+j1*n]);
			}
		}

		*t=tt;
	}
	free(iwk); free(wb); free(wk);
	return 0;
}
