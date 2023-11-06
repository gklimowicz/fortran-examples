/*	To solve linear second order elliptic equations using ADI method */

#include <stdio.h>
#include <math.h>

void cof(double x, double y, double *axx, double *ayy,
		double *ax, double *ay, double *a0, double *f);
void  bc(int ib, double x, double y, double *a0, double *an, double *f);

int adi(double x0, double xn, double y0, double yn, int *kn, int nx, int ny,
	double x[], double y[], double *u, int iu,
	void cof(double , double , double * , double * , double * , double * ,
		double * , double * ),
	void bc(int , double , double , double * , double * , double * ),
	double *el, double *eu, double aeps, int *nit);

main()
{
	int i,i1,j,nx,k,m,ny,id, iflg, ier,np,nmax;
	double hh, x[100], f[100][51],reps,el,eu,fi,x0,xn,y0,yn,y[100];

/*	Example 14.5 : Laplace's equation */

	id=51; reps=1.e-6;
	x0=0.0; xn=1.0; 
	y0=0.0; yn=1.0;

/*	k is the parameter in ADI method, 2**k should be close to the number
	of iterations required. */

	for(i1=0; i1<99; ++i1) {
		printf("type  nx,ny = no. of points along x,y;  k\n");
		printf("                      (quits when nx<=0)\n");
		scanf(" %d %d %d",&nx,&ny,&k);
		if(nx<=0) return 0;

/*	Set the initial values to zero */
		for(i=0; i<ny; ++i) {
			for(j=0; j<nx; ++j) f[i][j]=0.0;
		}
		el=0.0; eu=0.0;

		i=adi(x0,xn,y0,yn,&k,nx,ny,x,y,&f[0][0],id,cof,bc,&el,&eu,reps,&np);
		printf(" ier = %d   nx =  %d   ny = %d    k = %d    No. of iterations = %d \n",i,nx,ny,k,np);
		printf(" The solution :\n");
		for(i=0; i<ny; i += 1) {
			printf(" %e",y[i]);
			for(j=0; j<nx; j += 1) printf(" %e ",f[i][j]);
			printf(" \n");
		}

	}
	return;
}

/*	The coefficients of elliptic equations */

void cof(double x, double y, double *axx,  double *ayy,
		double *ax, double *ay, double *a0, double *f)

{
	*axx=1.0; *ayy=1.0;  *ax=0.0; *ay=0.0;
	*a0=0.0; *f=0.0;
	return;
}


/*	The boundary conditions */

#define PI 3.14159265358979324

void bc(int ib, double x, double y, double *a0, double *ax, double *f)

{
	*a0=1.0; *ax=0.0;
	if(ib<=3) *f=0.0;
	else *f=sin(PI*x);
	return;
}



/*	To solve linear second order elliptic differential equation using the
	Alternating direction implicit iterative (ADI) method
	The differential equation is assumed to be of the form

	Axx(x,y)d^2u/dx^2 + Ayy(x,y)d^2u/dy^2 + Ax(x,y)du/dx
		+ Ay(x,y)du/dy + A0(x,y)u + F(x,y)=0

	with following boundary conditions on a rectangular region

		A0*u+An*dun=F; where dun is the normal derivative of u

	X0 : (input) Lower limit on X where the solution is required
	XN : (input) Upper limit on X where the solution is required.
		Solution is computed in the interval (X0,XN)
	Y0 : (input) Lower limit on Y where the solution is required
	YN : (input) Upper limit on Y where the solution is required.
		Solution is computed in the interval (Y0,YN)
	KN : (input/output) Parameter k in ADI iteration. The function repeats
		a cycle of 2^k iteration, but convergence is checked after
		each iteration. This parameter may be adjusted if it is
		outside acceptable limits.
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
	EL : (input/output) Lower limit on the eigenvalues of the partitions
		S and Y of the finite difference matrix. If EL<=0, then
		the value for Poisson's equation will be used.
	EU : (input/output) Upper limit on the eigenvalues of the partitions
		S and Y of the finite difference matrix. If EU<=EL, then
		the value for Poisson's equation will be used.
	AEPS : (input) Required absolute accuracy. The ADI iteration is
			continued until the change in all elements is less
			than AEPS
	NIT : (output) Number of ADI iterations required by the function.

	Error status is returned by the value of the function ADI.
		0 value implies successful execution
		717 implies that YN=Y0, XN=X0, NX<3, NY<3 or IU<NX
			in which case no calculations are done
		766 implies that the matrix for ADI iteration is singular
			and calculations have to be abandoned.
		767 implies that ADI iteration failed to converge to
			specified accuracy

	Functions COF(X,Y,AXX,AYY,AX,AY,A0,F) and BC(IB,X,Y,A0,AN,F)
	must be supplied by the user 
	Function COF should calculate the coefficients AXX, AYY, 
	AX, AY, A0, F as defined above for given values of X,Y.
	Function BC should calculate the Boundary values at each boundary.
	Here IB is an integer denoting which boundary is being considered.
	IB=1 implies x=X0, IB=2 implies x=XN, IB=3 implies y=Y0, IB=4 implies y=YN

	Required functions : COF, BC
*/

#include <math.h>
#include <stdlib.h>

#define PI2 1.570796326794896619

int adi(double x0, double xn, double y0, double yn, int *kn, int nx, int ny,
	double x[], double y[], double *u, int iu,
	void cof(double , double , double * , double * , double * , double * ,
		double * , double * ),
	void bc(int , double , double , double * , double * , double * ),
	double *el, double *eu, double aeps, int *nit)

{
	int i,j,k,i1,j1,l,jl,ju,kl,ku,n0,n1,n2,it, km=6, nk=32, maxit=1000;
	double dx,dy,dx2,dy2,el1,el2,eu1,eu2,pr,disc,axx,ayy,ax,ay,a0,f,r,rp;
	double err,u1,alp[6],bet[6],sa[6][32];
	double *wk,*wka;

	if(xn==x0 || yn==y0 || nx<=2 || ny<=2 || iu<nx) return 717;

/*	setting up the mesh points */
	dx=(xn-x0)/(nx-1);
	for(i=0; i<nx; ++i) x[i]=x0+i*dx;
	dy=(yn-y0)/(ny-1);
	for(i=0; i<ny; ++i) y[i]=y0+i*dy;

	dx2=dx*dx; dy2=dy*dy;
	if(*el<=0.0) {
/*	Use the value for Poisson's equation */
		r=sin(PI2/(nx-1.0));
		el1=4.*r*r/dx2;
		r=sin(PI2/(ny-1.0));
		el2=4.*r*r/dy2;
		*el=el1; if(el2<el1) *el=el2;
	}

	if(*eu<=(*el)) {
/*	Use the value for Poisson's equation */
		r=cos(PI2/(nx-1.0));
		eu1=4.*r*r/dx2;
		r=cos(PI2/(ny-1.0));
		eu2=4.*r*r/dy2;
		*eu=eu1; if(eu2>eu1) *eu=eu2;
	}
	if(*kn>km-1) *kn=km-1;
	if(*kn<2) *kn=2;

/*	Calculating the parameters r_i for ADI iteration */
	alp[0]=(*el); bet[0]=(*eu);
	for(j=0; j<(*kn); ++j) {
		alp[j+1]=sqrt(alp[j]*bet[j]);
		bet[j+1]=0.5*(alp[j]+bet[j]);
	}
	sa[0][0]=sqrt(alp[*kn]*bet[*kn]);
	j1=1;
	for(j=0; j<(*kn); ++j) {
		pr=alp[*kn-j-1]*bet[*kn-j-1];
		for(l=0; l<j1; ++l) {
			disc=sa[j][l]*sa[j][l]-pr;
			if(disc<0.0) disc=0.0;
			sa[j+1][2*l]=sa[j][l]+sqrt(disc);
			sa[j+1][2*l+1]=pr/sa[j+1][2*l];
		}
		j1=j1*2;
	}

	wka=(double *) calloc((size_t) (7*nx*ny), sizeof(double));
	n1=7*nx;
/*	Setting up the finite difference equations */
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) {
			cof(x[i],y[j],&axx,&ayy,&ax,&ay,&a0,&f);
			wka[i*7+j*n1]=-axx/dx2+0.5*ax/dx;
			wka[1+i*7+j*n1]=2.0*axx/dx2-0.5*a0;
			wka[2+i*7+j*n1]=-axx/dx2-0.5*ax/dx;
			wka[3+i*7+j*n1]=-ayy/dy2+0.5*ay/dy;
			wka[4+i*7+j*n1]=2.0*ayy/dy2-0.5*a0;
			wka[5+i*7+j*n1]=-ayy/dy2-0.5*ay/dy;
			wka[6+i*7+j*n1]=f;
		}
	}

/*	The boundary conditions at x=x0 */
	jl=1; bc(1,x0,y0,&a0,&ax,&f);
	if(ax != 0.0) jl=0;
	for(k=0; k<ny; ++k) {
		bc(1,x0,y[k],&a0,&ax,&f);
		if(jl==0) {
			if(ax==0.0) {free(wka); return 766;}
			wka[1+k*n1]=wka[1+k*n1]+wka[k*n1]*2.*dx*a0/ax;
			wka[2+k*n1]=wka[2+k*n1]+wka[k*n1];
			wka[6+k*n1]=wka[6+k*n1]+wka[k*n1]*2.*dx*f/ax;
		}
		else {
			if(a0==0.0) {free(wka); return 766;}
			u[k*iu]=f/a0;
			wka[6+7+k*n1]=wka[13+k*n1]-wka[7+k*n1]*u[k*iu];
		}
	}

/*	The boundary conditions at x=xn */
	ju=nx-2; bc(2,xn,y0,&a0,&ax,&f);
	if(ax != 0.0) ju=nx-1;
	for(k=0; k<ny; ++k) {
		bc(2,xn,y[k],&a0,&ax,&f);
		if(ju==nx-1) {
			if(ax==0.0) {free(wka); return 766;}
			wka[1+7*(nx-1)+k*n1]=wka[7*nx-6+k*n1]-wka[2+(nx-1)*7+k*n1]*2.0*dx*a0/ax;
			wka[7*(nx-1)+k*n1]=wka[7*nx-7+k*n1]+wka[2+(nx-1)*7+k*n1];
			wka[6+7*(nx-1)+k*n1]=wka[7*nx-1+k*n1]-wka[2+(nx-1)*7+k*n1]*2.0*dx*f/ax;
		}
		else {
			if(a0==0.0) {free(wka); return 766;}
			u[nx-1+k*iu]=f/a0;
			wka[6+(nx-2)*7+k*n1]=wka[6+(nx-2)*7+k*n1]-wka[2+(nx-2)*7+k*n1]*u[nx-1+k*iu];
		}
	}

/*	The boundary conditions at y=y0 */
	kl=1; bc(3,x0,y0,&a0,&ax,&f);
	if(ax != 0.0) kl=0;
	for(j=jl; j<=ju; ++j) {
		bc(3,x[j],y0,&a0,&ax,&f);
		if(kl==0) {
			if(ax==0.0) {free(wka); return 766;}
			wka[4+j*7]=wka[4+j*7]+wka[3+j*7]*2.0*dy*a0/ax;
			wka[5+j*7]=wka[5+j*7]+wka[3+j*7];
			wka[6+j*7]=wka[6+j*7]+wka[3+j*7]*2.0*dy*f/ax;
		}
		else {
			if(a0==0.0) {free(wka); return 766;}
			u[j]=f/a0;
			wka[6+j*7+n1]=wka[6+j*7+n1]-wka[3+j*7+n1]*u[j];
		}
	}

/*	The boundary conditions at y=yn */
	ku=ny-2; bc(4,x0,yn,&a0,&ax,&f);
	if(ax != 0.0) ku=ny-1;
	for(j=jl; j<=ju; ++j) {
		bc(4,x[j],yn,&a0,&ax,&f);
		if(ku==ny-1) {
			if(ax==0.0) {free(wka); return 766;}
			wka[4+j*7+(ny-1)*n1]=wka[4+j*7+(ny-1)*n1]-wka[5+j*7+(ny-1)*n1]*2.*dy*a0/ax;
			wka[3+j*7+(ny-1)*n1]=wka[3+j*7+(ny-1)*n1]+wka[5+j*7+(ny-1)*n1];
			wka[6+j*7+(ny-1)*n1]=wka[6+j*7+(ny-1)*n1]-wka[5+j*7+(ny-1)*n1]*2.*dy*f/ax;
		}
		else {
			if(a0==0.0) {free(wka); return 766;}
			u[j+(ny-1)*iu]=f/a0;
			wka[6+j*7+(ny-2)*n1]=wka[6+j*7+(ny-2)*n1]-wka[5+j*7+(ny-2)*n1]*u[j+(ny-1)*iu];
		}
	}

/*	Loop for ADI iteration */
	n0=pow(2.0, (double ) *kn);
	n2=nx; if(ny>nx) n2=ny;
	wk=(double *) calloc((size_t) (n2*(2+n2)), sizeof(double));
	for(it=0; it<maxit; ++it) {
		i1=it-n0*(it/n0);
		r=sa[*kn][i1];

/*	Gaussian elimination for the first half step */
		for(k=kl; k<=ku; ++k) {
			wk[jl]=wka[1+jl*7+k*n1]+r;
			wk[jl+(2+k)*n2]=wka[6+jl*7+k*n1]-(wka[4+jl*7+k*n1]-r)*u[jl+k*iu];
			if(k>kl) wk[jl+(2+k)*n2]=wk[jl+(2+k)*n2]-wka[3+jl*7+k*n1]*u[jl+(k-1)*iu];
			if(k<ku) wk[jl+(2+k)*n2]=wk[jl+(2+k)*n2]-wka[5+jl*7+k*n1]*u[jl+(k+1)*iu];
			for(j=jl+1; j<=ju; ++j) {
				if(wk[j-1]==0.0) {free(wk); free(wka); return 766;}
				rp=-wka[j*7+k*n1]/wk[j-1];
				wk[j]=wka[1+j*7+k*n1]+r+rp*wka[2+(j-1)*7+k*n1];
				wk[j+(2+k)*n2]=wka[6+j*7+k*n1]-(wka[4+j*7+k*n1]-r)*u[j+k*iu]
								+rp*wk[j-1+(2+k)*n2];
				if(k>kl) wk[j+(2+k)*n2]=wk[j+(2+k)*n2]-wka[3+j*7+k*n1]*u[j+(k-1)*iu];
				if(k<ku) wk[j+(2+k)*n2]=wk[j+(2+k)*n2]-wka[5+j*7+k*n1]*u[j+(k+1)*iu];
			}
			if(wk[ju]==0.0) {free(wk); free(wka); return 766;}

/*	Back-substitution */
			wk[ju+(2+k)*n2]=wk[ju+(2+k)*n2]/wk[ju];
			for(j=ju-1; j>=jl; --j) wk[j+(2+k)*n2]=(wk[j+(2+k)*n2]-
							wka[2+j*7+k*n1]*wk[j+1+(2+k)*n2])/wk[j];
		}
		
/*	Gaussian elimination for the second half-step */
		err=0.0;
		for(k=jl; k<=ju; ++k) {
			wk[kl]=wka[4+k*7+kl*n1]+r;
			wk[kl+n2]=wka[6+k*7+kl*n1]-(wka[1+k*7+kl*n1]-r)*wk[k+(kl+2)*n2];
			if(k>jl) wk[kl+n2]=wk[kl+n2]-wka[k*7+kl*n1]*wk[k-1+(kl+2)*n2];
			if(k<ju) wk[kl+n2]=wk[kl+n2]-wka[2+k*7+kl*n1]*wk[k+1+(kl+2)*n2];
			for(j=kl+1; j<=ku; ++j) {
				if(wk[j-1]==0.0) {free(wk); free(wka); return 766;}
				rp=-wka[3+k*7+j*n1]/wk[j-1];
				wk[j]=wka[4+k*7+j*n1]+r+rp*wka[5+k*7+(j-1)*n1];
				wk[j+n2]=wka[6+k*7+j*n1]-(wka[1+k*7+j*n1]-r)*wk[k+(j+2)*n2]
							+rp*wk[j-1+n2];
				if(k>jl) wk[j+n2]=wk[j+n2]-wka[k*7+j*n1]*wk[k-1+(j+2)*n2];
				if(k<ju) wk[j+n2]=wk[j+n2]-wka[2+k*7+j*n1]*wk[k+1+(j+2)*n2];
			}
			if(wk[ku]==0.0) {free(wk); free(wka); return 766;}

/*	Back-substitution */
			u1=wk[ku+n2]/wk[ku];
			if(fabs(u1-u[k+ku*iu])>err) err=fabs(u1-u[k+ku*iu]);
			u[k+ku*iu]=u1;
			for(j=ku-1; j>=kl; --j) {
				u1=(wk[j+n2]-wka[5+k*7+j*n1]*u[k+(j+1)*iu])/wk[j];
				if(fabs(u1-u[k+j*iu])>err) err=fabs(u1-u[k+j*iu]);
				u[k+j*iu]=u1;
			}

		}
		if(err<aeps) {*nit=it+1; free(wk); free(wka); return 0;}
	}

	free(wk); free(wka);
	return 767;
}
