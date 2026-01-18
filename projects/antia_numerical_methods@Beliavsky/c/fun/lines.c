/*	To solve a system of nonlinear parabolic differential equation using
	the method of lines
	The differential equation is assumed to be of the form

	du/dt = f_i( d^2u/dx^2, du/dx, u, x, t),    i=0,1,...,NELIN-1

	with Dirichlet boundary conditions

	u(x0,t) = g1(t);  	u(xN,t) = g2(t)

	This function is used through function MSTEP or RKM for solving
	a system of ordinary differential equations. This function specifies
	the equivalent system of ordinary differential equations.
	The parameters for the partial differential equations are passed
	through global variables as defined below. These variables must
	be initialised before calling MSTEP or RKM.

	T : (input) Value of "time" where the derivatives need to be calculated.
	N : (input) Number of ordinary differential equations to be solved
		after applying the method of lines, N=NELIN*(NXLIN-2)
	U : (input) Array of length N containing the estimated
		solution at t.
	DU : (output) Array of length N which will contain the derivatives
		of U at T.

	The global variable used for passing information about equations are:

	DXLIN : (input) The step length in X to be used for computations.
		It is assumed to be uniform.
	X0LIN : (input) Lower limit on X where the solution is required
	XNLIN : (input) Upper limit on X where the solution is required.
		Solution is computed in the interval (X0LIN,XNLIN)
	XLIN : (input) Array of length NXLIN-2 containing the mesh points used
		in X direction. The end points X0LIN and XNLIN are not included.
		This array must be set before calling MSTEP or RKM.
	U0LIN : (output) Array of length NELIN containing the solution
		at the first boundary X0LIN. It is calculated using the
		boundary conditions calculated by function BC.
	UNLIN : (output) Array of length NELIN containing the solution
		at the second boundary XNLIN. It is calculated using the
		boundary conditions calculated by function BC.
	NELIN : (input) Number of differential equations in the parabolic system
	NXLIN : (input) Number of mesh points in the X direction.

	This function requires function FCN to calculate the derivatives
	defining the differential equations and function BC to calculate
	the boundary conditions. The names of these functions are fixed and
	cannot be passed on as arguments.

	Functions FCN(NELIN,X,T,U,UX,UXX,DU) and
	BC(NELIN,T,X0LIN,XNLIN,U0LIN,UNLIN) must be supplied by the user 
	Function FCN should calculate the derivatives DU = dU/dt at
	specified T, X, U, UX (=dU/dX), UXX (=d^2U/dX^2)
	Here NELIN is the number of parabolic equations in the system, T and X
	are the values of t and x, U, UX, UXX and DU are arrays of
	length NELIN containing the solution and its derivative as described
	above.

	Function BC should calculate the solution at boundaries at x=X0LIN
	and XNLIN. Here NELIN is the number of parabolic equations in the system,
	T is the value of t at which boundary conditions are required.
	X0LIN and XNLIN are the two boundaries in X where boundary conditions are
	applied. U0LIN and UNLIN are arrays of length NELIN containing the
	solution at X0LIN and XNLIN.
	U0LIN and UNLIN must be calculated by function BC.

	Required functions : FCN, BC
*/

#include <math.h>
#include <stdlib.h>

/* Include the following before the calling function to pass on
the parameters to lines */

#define NPL 101
#define NQL 5

double DXLIN, X0LIN, XNLIN, XLIN[NPL], U0LIN[NQL], UNLIN[NQL];
int NXLIN, NELIN;

void fcn(int ne, double x, double t, double u[], double ux[], double uxx[],
	double du[]);
void bc(int ne, double t, double x0, double xn, double u0[], double un[]);

void lines(double t, int n, double u[], double du[])

{
	int i,j,k;
	double um,up;
	double *ux, *uxx;

/*	Calculate the boundary values */
	bc(NELIN, t, X0LIN, XNLIN, U0LIN, UNLIN);
	ux=(double *) calloc((size_t) (NELIN+1), sizeof(double));
	uxx=(double *) calloc((size_t) (NELIN+1), sizeof(double)); 
	for(i=0; i<n; i += NELIN) {
		j=i/NELIN;

		for(k=0; k<NELIN; ++k) {
			if(i==0) um=U0LIN[k];
			else um=u[(j-1)*NELIN+k];

			if(i+NELIN >=n) up=UNLIN[k];
			else up=u[(j+1)*NELIN+k];

			ux[k]=(up-um)/(2.*DXLIN);
			uxx[k]=(up-2.0*u[i+k]+um)/(DXLIN*DXLIN);
		}
		
		fcn(NELIN, XLIN[j],t,&u[i],ux,uxx, &du[i]);
	}
	free(uxx); free(ux); 
	return;
}

