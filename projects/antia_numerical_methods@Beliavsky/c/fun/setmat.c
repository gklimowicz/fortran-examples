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
