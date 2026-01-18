/*	To solve linear second order elliptic differential equation using the
	successive over-relaxation (SOR) method
	The differential equation is assumed to be of the form

	Axx(x,y)d^2u/dx^2 + Axy(x,y)d^2u/dydx + Ayy(x,y)d^2u/dy^2 +
		Ax(x,y)du/dx + Ay(x,y)du/dy + A0(x,y)u + F(x,y)=0

	with Dirichlet boundary conditions on a rectangular region

	X0 : (input) Lower limit on X where the solution is required
	XN : (input) Upper limit on X where the solution is required.
		Solution is computed in the interval (X0,XN)
	Y0 : (input) Lower limit on Y where the solution is required
	YN : (input) Upper limit on Y where the solution is required.
		Solution is computed in the interval (Y0,YN)
	NX : (input) Number of mesh points in the X direction.
	NY : (input) Number of mesh points in the Y direction.
	X : (output) Array of length NX containing the mesh points used
		in X direction. These are calculated by assuming uniform spacing.
	Y : (output) Array of length NY containing the mesh points used
		in Y direction. These are calculated by assuming uniform spacing.
	U : (input/output) Array of length IU*NY containing the solution
		It should contain the initial values at the time of calling.
		After execution it will contain	the computed solution.
		U[j][i] is the solution at (x_i,y_j)
	IU : (input) The second dimension of U as declared in the calling
		function, IU>=NX
	COF : (input) Name of the function to calculate the coefficients
		in the equation
	BC : (input) Name of the function to calculate the boundary conditions
	OMEGA : (input/output) Value of the relaxation parameter, 1<OMEGA<2
		If OMEGA <= 0 then the function sets it to the
		optimal value for Poisson's equation.
	AEPS : (input) Required absolute accuracy. The SOR iteration is
			continued until the change in all elements is less than AEPS
	NIT : (output) Number of SOR iterations required by the function.
		
	Error status is returned by the value of the function SOR.
		0 value implies successful execution
		716 implies that YN=Y0, XN=X0, NX<3, NY<3 or IU<NX
			in which case no calculations are done
		764 implies that the diagonal term in the difference
			equation vanishes and calculations have to be abandoned.
		765 implies that SOR iteration failed to converge to
			specified accuracy

	Functions COF(X,Y,AXX,AXY,AYY,AX,AY,A0,F) and BC(IB,X,Y)
	must be supplied by the user 
	Function COF should calculate the coefficients AXX, AXY, AYY, 
	AX, AY, A0, F as defined above for given values of X,Y.
	Function BC should calculate the Boundary values at each boundary.
	Here IB is an integer denoting which boundary is being considered.
	The boundary conditions are assumed to be
	u(X0,Y)=BC(1,X0,Y);	u(XN,Y)=BC(2,XN,Y);
	u(x,Y0)=BC(3,x,Y0);	u(x,YN)=BC(4,x,YN)

	Required functions : COF, BC
*/

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979324

int sor(double x0, double xn, double y0, double yn, int nx, int ny,
	double x[], double y[], double *u, int iu,
	void cof(double , double , double * , double * , double * ,
		double * , double * , double * , double *),
	double bc(int , double , double ), double *omega, double aeps, int *nit)

{
	int i,j,k,it,n1, maxit=1000;
	double dx,dy,dx2,dy2,dxy,rj,ad,axx,axy,ayy,ax,ay,a0,f,err,res;
	double *wk;

	if(xn==x0 || yn==y0 || nx<=2 || ny<=2 || iu<nx) return 716;

/*	Setting up the mesh points */
	dx=(xn-x0)/(nx-1);
	for(i=0; i<nx; ++i) x[i]=x0+i*dx;
	dy=(yn-y0)/(ny-1);
	for(i=0; i<ny; ++i) y[i]=y0+i*dy;

	dx2=dx*dx;
	dy2=dy*dy;
	dxy=4.*dx*dy;
	if(*omega<=0.0) {
/*	Estimate the optimal value of OMEGA */
		rj=(dy2*cos(PI/(nx-1.0))+dx2*cos(PI/(ny-1.0)))/(dx2+dy2);
		*omega=2.0/(1.0+sqrt(1.-rj*rj));
	}
	wk=(double *) calloc((size_t) (9*nx*ny), sizeof(double));
	n1=9*nx;

/*	Calculate the coefficients of the difference equations */
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) {
			cof(x[i],y[j],&axx,&axy,&ayy,&ax,&ay,&a0,&f);
			ad=2.0*axx/dx2+2.*ayy/dy2-a0;
			if(ad==0) {free(wk); return 764;}
			wk[i*9+j*n1]=-axy/(dxy*ad);
			wk[1+i*9+j*n1]=(-ayy/dy2+0.5*ay/dy)/ad;
			wk[2+i*9+j*n1]=-wk[i*9+j*n1];
			wk[3+i*9+j*n1]=(-axx/dx2+0.5*ax/dx)/ad;
			wk[4+i*9+j*n1]=(-axx/dx2-0.5*ax/dx)/ad;
			wk[5+i*9+j*n1]=wk[2+i*9+j*n1];
			wk[6+i*9+j*n1]=(-ayy/dy2-0.5*ay/dy)/ad;
			wk[7+i*9+j*n1]=wk[i*9+j*n1];
			wk[8+i*9+j*n1]=f/ad;
		}
	}

/*	Calculate the boundary values */
	for(k=0; k<ny; ++k) {
		u[k*iu]=bc(1,x0,y[k]);
		u[nx-1+k*iu]=bc(2,xn,y[k]);
	}
	for(j=0; j<nx; ++j) {
		u[j]=bc(3,x[j],y0);
		u[j+(ny-1)*iu]=bc(4,x[j],yn);
	}

/*	Loop for the SOR iteration */
	for(it=1; it<=maxit; ++it) {
		err=0.0;
		for(k=1; k<ny-1; ++k) {
			for(j=1; j<nx-1; ++j) {
				res=wk[8+j*9+k*n1]-u[j+k*iu]-wk[j*9+k*n1]*u[j-1+(k-1)*iu]
					-wk[1+j*9+k*n1]*u[j+(k-1)*iu]-wk[2+j*9+k*n1]*u[j+1+(k-1)*iu]
					-wk[3+j*9+k*n1]*u[j-1+k*iu]-wk[4+j*9+k*n1]*u[j+1+k*iu]
					-wk[5+j*9+k*n1]*u[j-1+(k+1)*iu]-wk[6+j*9+k*n1]*u[j+(k+1)*iu]
					-wk[7+j*9+k*n1]*u[j+1+(k+1)*iu];

				if(fabs(res)>err) err=fabs(res);
				u[j+k*iu]=u[j+k*iu]+(*omega)*res;
			}
		}

		if(err<aeps) {*nit=it; free(wk); return 0;}

	}

	*nit=maxit;
	free(wk);
	return 765;
}
