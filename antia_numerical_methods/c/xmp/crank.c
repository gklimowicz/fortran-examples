/*	Solution of a linear second order parabolic equation using
	Crank-Nicolson method */

#include <stdio.h>
#include <math.h>

void cof(double x, double t, double *a, double *b, double *c, double *d);
void bc(double t, double x0, double xn, double *a0, double *b0,
		double *f0, double *an, double *bn, double *fn);
double fic(double x, double t);

int crank(double *t, double dt, double x0, double xn, int nt, int nx,
	double x[], double u[],
	void cof(double , double , double * , double * , double * , double * ),
	void bc(double , double , double , double * , double * , double * ,
		double * , double * , double * ),
	double fic(double , double ), int *iflg);

main()
{
	int i,i1,j,nx,nt,m, id, iflg, ier,np,nmax;
	double hh, x[100], wt[100],f[100],fc[100],tn,erc,reps,t0,fi;
	double x0,xn,dt;

/*	Example 14.1 */

	x0=1.0; xn=0.5; iflg=0;
	printf("type  nx=no. of points in x,   t0=initial time \n");
	scanf(" %d %le",&nx,&t0);

	for(i1=0; i1<99; ++i1) {
		printf("type nt=no. of time steps,    dt=step size\n");
		printf("                       (quits when nt<=0)\n");
		scanf(" %d %le",&nt,&dt);
		if(nt<=0) return 0;

		i=crank(&t0,dt,x0,xn,nt,nx,x,f,cof,bc,fic,&iflg);
		printf(" ier = %d    nx =  %d    nt = %d    t = %e    dt = %e \n",i,nx,nt,t0,dt);
		printf("  x           f(x,t)        exact value\n");
		for(i=0; i<nx; ++i) printf(" %e %e %e \n",x[i],f[i],fic(x[i],t0));

	}
	return;
}

/*	Coefficient for diffusion equation */

void cof(double x, double t, double *a, double *b, double *c,
		double *d)

{
	*a=1; *b=0.0; *c=0.0; *d=0.0;
	return;
}

/*	The boundary conditions */

void bc(double t, double x0, double xn, double *a0, double *b0,
		double *f0, double *an, double *bn, double *fn)

{
	*a0=1.0; *b0=0.0; *f0=0.0;
	*an=0.0; *bn=1.0; *fn=0.0;
	return;
}

/*	The exact solution used to generate boundary values */

double fic(double x, double t)

{
	double pi=3.1415926535l;

	return sin(pi*x)*exp(-pi*pi*t);
}



/*	To solve linear parabolic differential equation using the
	Crank-Nicolson difference scheme
	The differential equation is assumed to be of the form

	du/dt = A(x,t)d^2u/dx^2 + B(x,t)du/dx + C(x,t)u + D(x,t)

	with boundary conditions

	A0(t) u(x0,t)+B0(t) du(x0,t)/dx = F0(t)
	AN(t) u(xN,t)+BN(t) du(xN,t)/dx = FN(t)

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
	U : (input/output) Array of length NX containing the solution
		at the current time step. It should contain the initial
		values at the time of calling, if IFLG!=0. Otherwise it
		is computed using Function FIC. After execution it will contain
		the computed solution at t=T.
		It may be noted that solution at intermediate time steps
		is not preserved. Hence, if solution is required at number
		of time steps, multiple calls to CRANK will be needed.
	COF : (input) Name of the function to calculate the coefficients
		in the equation
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
		
	Error status is returned by the value of the function CRANK.
		0 value implies successful execution
		713 implies that DT=0, XN=X0 or NX<3
			in which case no calculations are done
		761 implies that the difference equations are singular
			and solution cannot be continued further.

	Functions COF(X,T,A,B,C,D), BC(T,X0,XN,A0,B0,F0,AN,BN,FN)
	and FIC(X,T) must be supplied by the user 
	Function COF should calculate the coefficients A, B, C, D as
	defined above for given values of X,T.
	Function BC should calculate the coefficients A0, B0, F0, AN, BN, FN
	for given values of T, X0 and XN.
	Function FIC is required only if IFLG=0, otherwise a dummy function
	with this name will suffice. If IFLG=0, function FIC must calculate
	the initial values at required X,T.

	Required functions : COF, BC, FIC
*/

#include <math.h>
#include <stdlib.h>

int crank(double *t, double dt, double x0, double xn, int nt, int nx,
	double x[], double u[],
	void cof(double , double , double * , double * , double * , double * ),
	void bc(double , double , double , double * , double * , double * ,
		double * , double * , double * ),
	double fic(double , double ), int *iflg)

{
	int i,j,jl,ju,it;
	double dx,r,r1,a0,b0,f0,an,bn,fn,u0,un,t1,um,up,a,b,c,d,rp;
	double *wk;

	if(dt==0.0 || xn==x0 || nx<=2) return 713;
	wk=(double *) calloc((size_t) (8*nx), sizeof(double));

/*	Set up the step size in X as well as the initial values */
	dx=(xn-x0)/(nx-1);
	for(i=0; i<nx; ++i) {
		x[i]=x0+i*dx;
		cof(x[i],*t,&wk[i+4*nx],&wk[i+5*nx],&wk[i+6*nx],&wk[i+7*nx]);
		if(*iflg==0) u[i]=fic(x[i],*t);
	}
	*iflg=1;
	r=0.5*dt/(dx*dx);
	r1=0.25*dt/dx;

/*	Boundary condition at X0 */
	bc(*t,x0,xn,&a0,&b0,&f0,&an,&bn,&fn);
	jl=1;
	if(b0 != 0.0) {jl=0; u0=u[1]+2.*dx*(a0*u[0]-f0)/b0;}
	
/*	boundary condition at XN */
	ju=nx-2;
	if(bn != 0.0) {ju=nx-1; un=u[nx-2]-2.*dx*(an*u[nx-1]-fn)/bn;}

/*	Loop for time steps */
	for(it=1; it<=nt; ++it) {
		t1=(*t)+dt;
		bc(t1,x0,xn,&a0,&b0,&f0,&an,&bn,&fn);

		for(j=jl; j<=ju; ++j) {
			if(j==0) um=u0;
			else um=u[j-1];

			if(j==nx-1) up=un;
			else up=u[j+1];

/*	The Crank-Nicolson difference scheme */
			cof(x[j],t1,&a,&b,&c,&d);
			wk[j]=-r*a+r1*b;
			wk[j+nx]=1+2*r*a-0.5*dt*c;
			wk[j+2*nx]=-r*a-r1*b;
			wk[j+3*nx]=u[j]+r*wk[j+4*nx]*(up-2*u[j]+um)+r1*wk[j+5*nx]*(up-um)
						+0.5*dt*(wk[j+6*nx]*u[j]+wk[j+7*nx]+d);
			wk[j+4*nx]=a;
			wk[j+5*nx]=b;
			wk[j+6*nx]=c;
			wk[j+7*nx]=d;
		}

/*	Boundary condition at X0 */
		if(jl==1) {
			if(a0==0.0) {free(wk); return 761;}
			u[0]=f0/a0;
			wk[1+3*nx]=wk[1+3*nx]-u[0]*wk[1];
		}
		else {
			wk[nx]=wk[nx]+2.0*wk[0]*a0*dx/b0;
			wk[2*nx]=wk[2*nx]+wk[0];
			wk[3*nx]=wk[3*nx]+2.*wk[0]*f0*dx/b0;
		}

/*	Boundary condition at XN */
		if(ju==nx-1) {
			wk[ju]=wk[ju]+wk[ju+2*nx];
			wk[ju+nx]=wk[ju+nx]-2.*wk[ju+2*nx]*an*dx/bn;
			wk[ju+3*nx]=wk[ju+3*nx]-2.*wk[ju+2*nx]*fn*dx/bn;
		}
		else {
			if(an==0.0) {free(wk); return 761;}
			u[nx-1]=fn/an;
			wk[nx-2+3*nx]=wk[4*nx-2]-u[nx-1]*wk[3*nx-2];
		}

/*	Gaussian elimination for tridiagonal matrix */
		for(j=jl+1; j<=ju; ++j) {
			if(wk[j-1+nx]==0.0) {free(wk); return 761;}
			rp=-wk[j]/wk[j-1+nx];
			wk[j+nx]=wk[j+nx]+rp*wk[j-1+2*nx];
			wk[j+3*nx]=wk[j+3*nx]+rp*wk[j-1+3*nx];
		}
		if(wk[ju+nx]==0.0) {free(wk); return 761;}

/*	Back-substitution for tridiagonal system */
		u[ju]=wk[ju+3*nx]/wk[ju+nx];
		for(j=ju-1; j>=jl; --j) u[j]=(wk[j+3*nx]-wk[j+2*nx]*u[j+1])/wk[j+nx];

		if(jl==0) u0=u[1]+2.*dx*(a0*u[0]-f0)/b0;
		if(ju==nx-1) un=u[nx-2]-2.*dx*(an*u[nx-1]-fn)/bn;
		*t=t1;
	}
	free(wk);
	return 0;
}
