/*	To solve linear parabolic differential equation in 2 dimensions
	using the alternating direction method
	The differential equation is assumed to be of the form

	du/dt = Axx(x,y,t)d^2u/dx^2 + Ayy(x,y,t)d^2u/dy^2 + Ax(x,y,t)du/dx
		+Ay(x,y,t)du/dy + Au(x,y,t)u +A0(x,y,t)

	with Dirichlet boundary conditions

	u(X0,y,t)=BC(1,X0,y,t);		u(XN,y,t)=BC(2,XN,y,t)
	u(x,Y0,t)=BC(3,x,Y0,t);		u(x,YN,t)=BC(4,x,YN,t)

	T : (input/output) Initial value of "time" where the initial conditions
		are specified. After execution, it will be replaced by the
		value of T at the last point where execution is successful
	DT : (input) The time step to be used for computations. It is kept fixed.
	X0 : (input) Lower limit on X where the solution is required
	XN : (input) Upper limit on X where the solution is required.
		Solution is computed in the interval (X0,XN)
	Y0 : (input) Lower limit on Y where the solution is required
	YN : (input) Upper limit on Y where the solution is required.
		Solution is computed in the interval (Y0,YN)
	NT : (input) Number of time steps each of length DT  to be executed.
	NX : (input) Number of mesh points in the X direction.
	NY : (input) Number of mesh points in the Y direction.
	X : (output) Array of length NX containing the mesh points used
		in X direction. These are calculated by assuming uniform spacing.
	Y : (output) Array of length NY containing the mesh points used
		in Y direction. These are calculated by assuming uniform spacing.
	U : (input/output) Array of length IU*NY containing the solution
		at the current time step. It should contain the initial
		values at the time of calling, if IFLG!=0. Otherwise it
		is computed using Function FIC. After execution it will contain
		the computed solution at t=T. U[j][i] is the solution at (x_i,y_j)
	IU : (input) Second dimension of array U as declared in the calling function
	COF : (input) Name of the function to calculate the coefficients
		in the equation
	BC : (input) Name of the function to calculate the boundary values
	FIC : (input) Name of the function to calculate the initial
		values when IFLG=0. For other values of IFLG this function
		is not used, but a dummy function may be required by the compiler.
	IFLG : (input/output) Integer variable used as a flag to denote how the
			initial values are calculated.
		If IFLG=0 then the initial values are calculated using the
			function FIC to be supplied by the user. IFLG is set to 1
		Otherwise the initial values must be supplied in array U.

	Error status is returned by the value of the function ADM.
		0 value implies successful execution
		714 implies that DT=0, XN=X0, YN=Y0, NX<3, NY<3 or IU<NX
			in which case no calculations are done
		762 implies that the difference equations are singular
			and solution cannot be continued further.

	Functions COF(X,Y,T,AXX,AYY,AX,AY,AU,A0), BC(IB,X,Y,T)
	and FIC(X,Y,T) must be supplied by the user 
	Function COF should calculate the coefficients AXX, AYY, AX, AY
	AU, A0 defined above for given values of X, Y, T.
	Function BC should calculate the value of solution at the boundaries
	as described earlier. Integer IB specifies which boundary is required.
	Function FIC is required only if IFLG=0, otherwise a dummy function
	with this name will suffice. If IFLG=0, function FIC must calculate
	the initial values at required X,Y,T.

	Required functions : COF, BC, FIC
*/

#include <math.h>
#include <stdlib.h>

int adm(double *t, double dt, double x0, double xn, double y0, double yn,
	int nt, int nx, int ny, double x[], double y[], double *u, int iu,
	void cof(double , double , double , double * , double * , double * ,
		double * , double * , double * ),
	double  bc(int , double , double , double ),
	double fic(double , double , double ), int *iflg)

{
	int i,j,k,it,nn;
	double dx,dy,r,r1,s,s1,t1,t2,u1,rp,axx,ayy,ax,ay,au,a0,bxx,byy,bx,by,bu,b0,wk4;
	double *wk;

	if(dt==0.0 || xn==x0 || yn==y0 || nx<=2 || ny<=2 || iu<nx) return 714;

/*	Setting up the grid points along X and Y */
	dx=(xn-x0)/(nx-1);
	for(i=0; i<nx; ++i) x[i]=x0+i*dx;
	dy=(yn-y0)/(ny-1);
	for(i=0; i<ny; ++i) y[i]=y0+i*dy;

	if(*iflg==0) {
/*	Calculate the initial values */
		for(j=0; j<ny; ++j) {
			for(i=0; i<nx; ++i) u[i+j*iu]=fic(x[i],y[j],*t);
		}
	}
	*iflg=1;
	r=0.5*dt/(dx*dx);
	r1=0.25*dt/dx;
	s=0.5*dt/(dy*dy);
	s1=0.25*dt/dy;
	nn=nx; if(ny>nx) nn=ny;
	wk=(double *) calloc((size_t) (nn*(4+nn)), sizeof(double));

/*	Loop over time steps */
	for(it=1; it<=nt; ++it) {
		t1=(*t)+dt/2.; t2=(*t)+dt;

/*	Setting up the equations for the first half-step */
		for(k=1; k<ny-1; ++k) {
			for(j=1; j<nx-1; ++j) {
				cof(x[j],y[k],*t,&axx,&ayy,&ax,&ay,&au,&a0);
				cof(x[j],y[k],t1,&bxx,&byy,&bx,&by,&bu,&b0);
				wk[j]=-r*bxx+r1*bx;
				wk[j+nn]=1.0+2.0*r*bxx-0.25*dt*bu;
				wk[j+2*nn]=-r*bxx-r1*bx;
				wk[j+(4+k)*nn]=u[j+k*iu]+s*ayy*(u[j+(k+1)*iu]-2*u[j+k*iu]+u[j+(k-1)*iu])
							+s1*ay*(u[j+(k+1)*iu]-u[j+(k-1)*iu])
							+0.25*dt*(au*u[j+k*iu]+a0+b0);
			}

/*	The boundary conditions */
			u1=bc(1,x0,y[k],t1);
			wk[(4+k)*nn]=u1;
			wk[1+(4+k)*nn]=wk[1+(4+k)*nn]-u1*wk[1];
			u1=bc(2,xn,y[k],t1);
			wk[nx-1+(4+k)*nn]=u1;
			wk[nx-2+(4+k)*nn]=wk[nx-2+(4+k)*nn]-u1*wk[nx-2+2*nn];

/*	Gaussian elimination for the tridiagonal system */
			for(j=2; j<nx-1; ++j) {
				if(wk[j-1+nn]==0.0) {free(wk); return 762;}
				rp=-wk[j]/wk[j-1+nn];
				wk[j+nn]=wk[j+nn]+rp*wk[j-1+2*nn];
				wk[j+(4+k)*nn]=wk[j+(4+k)*nn]+rp*wk[j-1+(4+k)*nn];
			}
			if(wk[nx-2+nn]==0.0) {free(wk); return 762;}

/*	Back-substitution */
			wk[nx-2+(4+k)*nn]=wk[nx-2+(4+k)*nn]/wk[nx-2+nn];
			for(j=nx-3; j>=1; --j) wk[j+(4+k)*nn]=(wk[j+(4+k)*nn]-
								wk[j+2*nn]*wk[j+1+(4+k)*nn])/wk[j+nn];
		}

/*	Setting up the equations for the second half step */
		for(k=1; k<nx-1; ++k) {
			for(j=1; j<ny-1; ++j) {
				cof(x[k],y[j],t1,&axx,&ayy,&ax,&ay,&au,&a0);
				cof(x[k],y[j],t2,&bxx,&byy,&bx,&by,&bu,&b0);
				wk[j]=-s*byy+s1*by;
				wk[j+nn]=1.0+2.0*s*byy-0.25*dt*bu;
				wk[j+2*nn]=-s*byy-s1*by;
				wk4=wk[k+(j+4)*nn]+r*axx*(wk[k+1+(j+4)*nn]-2*wk[k+(j+4)*nn]+wk[k-1+(j+4)*nn])
						+r1*ax*(wk[k+1+(j+4)*nn]-wk[k-1+(j+4)*nn])
						+0.25*dt*(au*wk[k+(j+4)*nn]+a0+b0);
				if(k>1) wk[k-1+(j+4)*nn]=wk[j+3*nn];
				wk[j+3*nn]=wk4;
			}

/*	The boundary conditions */
			wk[k+4*nn]=bc(3,x[k],y0,t2);
			wk[1+3*nn]=wk[1+3*nn]-wk[k+4*nn]*wk[1];
			wk[k+(ny+3)*nn]=bc(4,x[k],yn,t2);
			wk[ny-2+3*nn]=wk[ny-2+3*nn]-wk[k+(ny+3)*nn]*wk[ny-2+2*nn];

/*	Gaussian elimination for the tridiagonal system */
			for(j=2; j<ny-1; ++j) {
				if(wk[j-1+nn]==0.0) {free(wk); return 762;}
				rp=-wk[j]/wk[j-1+nn];
				wk[j+nn]=wk[j+nn]+rp*wk[j-1+2*nn];
				wk[j+3*nn]=wk[j+3*nn]+rp*wk[j-1+3*nn];
			}
			if(wk[ny-2+nn]==0.0) {free(wk); return 762;}

/*	Back-substitution */
			wk[ny-2+3*nn]=wk[ny-2+3*nn]/wk[ny-2+nn];
			for(j=ny-3; j>=1; --j) wk[j+3*nn]=(wk[j+3*nn]-wk[j+2*nn]*wk[j+1+3*nn])/wk[j+nn];
		}

		for(j=1; j<ny-1; ++j) wk[nx-2+(j+4)*nn]=wk[j+3*nn];
		for(j=0; j<ny; ++j) {
			u[j*iu]=bc(1,x0,y[j],t2);
			u[nx-1+j*iu]=bc(2,xn,y[j],t2);
			for(k=1; k<nx-1; ++k) u[k+j*iu]=wk[k+(j+4)*nn];
		}

		*t=t2;
	}
	free(wk);
	return 0;
}

