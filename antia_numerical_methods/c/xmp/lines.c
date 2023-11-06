/*	To solve nonlinear parabolic equations using method of lines */

#include <stdio.h>
#include <math.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ));
void lines(double t, int n, double y[], double dy[]);
int strt4(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double *h, double *t, double reps, int iflg, double tstep, int *nstp,
	double *wk);
int gear(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
	int ijm3, int ijm4, int *iflag, double *wk1, double *wk);
int adams(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
	int ijm3, int ijm4, double *wk);
int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int mstep(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double *h, double *t0, double tn, double yf[], double reps, int *nstp,
	int *nmax, int *iflg, int ist, double *wk);
void fcn(int ne, double x, double t, double u[], double ux[], double uxx[],
	double du[]);
void bc(int ne, double t, double x0, double xn, double u0[], double un[]);

double al;
/* Include the following before the calling program to pass on
the parameters to lines */

#define NPL 101
#define NQL 5

double DXLIN, X0LIN, XNLIN, XLIN[NPL], U0LIN[NQL], UNLIN[NQL];
int NXLIN, NELIN;

/*	Exercise 14.12 : Burger's equation */

main()
{
	int i,i1,j,n,m, id, iflg, ier,np,nmax;
	double hh, x[100], dx[100],y[20],tn,wk[500],reps,t0,t1,xt;

	NELIN=1; nmax=100000; reps=1.e-6; 
	X0LIN=0.0; XNLIN=1.0; hh=0.01; iflg=0;
	t0=0.0;
	printf("type NXLIN=no. of points in x,  id=0/1  for Adams/Gear,  hh=initial step size\n");
	scanf(" %d %d %le",&NXLIN,&id,&hh);
	n=NXLIN-2; DXLIN=(XNLIN-X0LIN)/(NXLIN-1);
	printf("  No. of points in x = %d    id = %d    initial step size = %e\n",NXLIN,id,hh);

/*	Set up the initial values */
	for(i=0; i<n; ++i) {
		XLIN[i]=X0LIN+(i+1)*DXLIN;
		x[i]=2./(1.0+exp(XLIN[i]));
	}
	for(i1=0; i1<99; ++i1) {
		printf("type tn= required time      (quits when tn<-20)\n");
		scanf(" %le", &tn);
		if(tn<-20) return 0;

		i=mstep(n,x,dx,lines,&hh,&t0,tn,y,reps,&np,&nmax,&iflg,id,&wk[0]);
		printf(" ier = %d    no. of function evaluations =  %d    t = %e    step size = %e \n",i,np,tn,hh);
		printf("   x          solution     exact value \n");
		for(i=0; i<n; ++i) printf(" %e %e %e \n",XLIN[i],y[i],2./(1.0+exp(XLIN[i]-tn)));
		printf(" \n");

	}
	return;
}

/*	The differential equation */

void fcn(int ne, double x, double t, double u[], double ux[], double uxx[],
	double du[])

{
	du[0]=uxx[0]-u[0]*ux[0];
	return;
}

/*	The boundary conditions */

void bc(int ne, double t, double x0, double xn, double u0[], double un[])

{
	u0[0]=2.0/(1.+exp(x0-t));
	un[0]=2.0/(1.+exp(xn-t));
	return;
}




/*	To perform one step of solution of initial value problems in ordinary
	differential equations using fourth order Adams-Bashforth-Moulton
	predictor corrector method. It is called by function MSTEP.

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length 7N containing the solution
		at last four points. After execution new point will be added.
	DY : (input/output) Array of length 7N containing the derivatives
		of Y at the points in Y
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input) Step length to be used.
	T : (input) Value of independent variable t, at which the solution
		is required.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls error in iteration on corrector.
	NSTP : (input/output) Number of calls to DIF required so far.
		The count is updated by the function.
	IJ : (input) The index j+1 in corrector formula
	IJM1 : (input) The index j in corrector formula
	IJM2 : (input) The index j-1 in corrector formula
	IJM3 : (input) The index j-2 in corrector formula
	IJM4 : (input) The index j-3 in corrector formula
	WK : Array of length 2N used to pass the predicted values

	Error status is returned by the value of the function ADAMS.
		0 value implies successful execution
		729 implies that iteration on corrector failed to converge

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : DIF
*/
	
#include <math.h>
int adams(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
	int ijm3, int ijm4, double *wk)

{
	int j,k, nit=10;
	double err,t1,t2,r1, cfac=0.2, eps=1.e-30;

	for(j=0; j<n; ++j) {
/*	The predictor */
		t1=y[j+ijm1*n]+h*(55*dy[j+ijm1*n]-59*dy[j+ijm2*n]+37*dy[j+ijm3*n]-9*dy[j+ijm4*n])/24.0;
/*	Modifier */
		y[j+ij*n]=t1+251.0*(y[j+ijm1*n]-wk[j+n])/270.0;
		wk[j]=t1;
	}

/*	Iteration on corrector */
	for(j=1; j<=nit; ++j) {
		*nstp=(*nstp)+1;
		dif(t,n,&y[ij*n],&dy[ij*n]);
		err=0.0;
		for(k=0; k<n; ++k) {
/*	The corrector */
			t2=y[k+ijm1*n]+h*(9*dy[k+ij*n]+19*dy[k+ijm1*n]-5*dy[k+ijm2*n]+dy[k+ijm3*n])/24.0;
			r1=fabs(t2)+fabs(t2-y[k+ijm1*n])+eps;
			t1=fabs((y[k+ij*n]-t2)/r1);
			if(t1>err) err=t1;
			y[k+ij*n]=t2;
		}

/*	The convergence test */
		if(err<cfac*reps) return 0;
	}

	return 729;
}




/*	Solution of a system of linear equations using Gaussian elimination
	with partial pivoting

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input/output) The matrix of coefficient of size LJ*N
	        A[i][j] is the coefficient of x_j in ith equation
	     	at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
	     	at output it will contain the solutions
	DET : (output) The determinant of the matrix
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter to specify the type of computation required
		If IFLG<=0, both elimination and solution are
			done and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUELM.
		0 value implies successful execution
		101 implies (N<=0 or N>LJ) 
		121 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg)

{
int i,j,k,km,l;
double r1,t1;

	if(n<=0 || n>lj) return 101;
 
	if((*iflg)<2) {
/*	Perform elimination  */
		*det=1.0;
		for(k=0; k<n-1; ++k) {
/*	Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			for(l=k; l<n; ++l)
				if(fabs(a[l*lj+k])>r1) {r1=fabs(a[l*lj+k]); km=l;}

			inc[k]=km;
			if(km != k) {
/*	Interchange the rows if needed  */
				for(l=k; l<n; ++l) 
				{t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;}
				*det=-(*det);
			}

			*det=(*det)*a[k*lj+k];
			if(a[k*lj+k]==0) return 121;

			for(l=k+1; l<n; ++l) {
				a[l*lj+k]=a[l*lj+k]/a[k*lj+k];
				for(i=k+1; i<n; ++i) a[l*lj+i]=a[l*lj+i]-a[l*lj+k]*a[k*lj+i];
			}
		}
		*det=(*det)*a[(n-1)*lj+n-1];
		inc[n-1]=n-1;
		if(a[(n-1)*lj+n-1]==0) return 121;

		if((*iflg)==1) {*iflg=2; return 0;}
		*iflg=2;
	}

/*	Solution for the num different right-hand sides  */
	for(j=0; j<num; ++j) {
/*	forward-substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k])
			{t1=x[j*lj+k]; x[j*lj+k]=x[j*lj+inc[k]]; x[j*lj+inc[k]]=t1;}
			for(l=k+1; l<n; ++l) x[j*lj+l]=x[j*lj+l]-a[l*lj+k]*x[j*lj+k];
		}

/*	back-substitution  */

		x[j*lj+n-1]=x[j*lj+n-1]/a[(n-1)*lj+n-1];
		for(k=n-2; k>=0; --k) {
			for(l=n-1; l>=k+1; --l) x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=x[j*lj+k]/a[k*lj+k];
		}
	}
	return 0;
}



/*	To perform one step of solution of initial value problems in ordinary
	differential equations using fourth order stiffly stable method.
	It is called by function MSTEP.

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length 7N containing the solution
		at last four points. After execution new point will be added.
	DY : (input/output) Array of length 7N containing the derivatives
		of Y at the points in Y
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input) Step length to be used.
	T : (input) Value of independent variable t, at which the solution
		is required.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls error in iteration on corrector.
	NSTP : (input/output) Number of calls to DIF required so far.
		The count is updated by the function.
	IJ : (input) The index j+1 in corrector formula
	IJM1 : (input) The index j in corrector formula
	IJM2 : (input) The index j-1 in corrector formula
	IJM3 : (input) The index j-2 in corrector formula
	IJM4 : (input) The index j-3 in corrector formula
	IFLAG : (input/output) Integer variable used as a flag
		If IFLAG=0 initial approximation to Jacobian is generated
			and IFLAG is set to 1
		Otherwise old approximation to J^{-1} is used.
	WK1 : (input/output) Array of length 3N used to pass the predicted values
	WK : (input/output) Array of length N*N used to store J^{-1},
			the inverse of Jacobian
		
	Error status is returned by the value of the function GEAR.
		0 value implies successful execution
		729 implies that iteration on corrector failed to converge
			or cannot be continued because estimated Jacobian is singular

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : GAUELM, DIF
*/
	
#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int gear(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
	int ijm3, int ijm4, int *iflag, double *wk1, double *wk)

{
	int i,j,k,ipas,ier,ip,iflg, nit=20;
	double det,t1,x1,dk,err,s1,s2,r2,ss, cfac=0.01, eps=1.e-30;
	int *iwk;
	double  *wm;

	for(j=0; j<n; ++j) {
/*	predictor */
		t1=y[j+ijm1*n]+h*(55*dy[j+ijm1*n]-59*dy[j+ijm2*n]+37*dy[j+ijm3*n]-9*dy[j+ijm4*n])/24.0;
/*	Modifier */
		y[j+ij*n]=t1+6275.0*(y[j+ijm1*n]-wk1[j+n])/8003.0;
		wk1[j]=t1;
	}

	dif(t,n,&y[ij*n],&dy[ij*n]); *nstp=(*nstp)+1;

/*	Residual in the corrector */
	for(j=0; j<n; ++j) wk1[j+2*n]=y[j+ij*n]-(48*y[j+ijm1*n]-36*y[j+ijm2*n]+
					16*y[j+ijm3*n]-3*y[j+ijm4*n]+12.*h*dy[j+ij*n])/25.0;

	ip=n*n; if(ip<4*n) ip=4*n;
	wm=(double *) calloc((size_t) ip, sizeof(double)); 

	if(*iflag==0) {
/*	Generate initial approximation to Jacobian J  */
		for(ip=0; ip<n; ++ip) {
			x1=y[ip+ij*n];
			dk=100.0*reps*fabs(x1);
			if(dk==0.0) dk=100.0*reps;
			y[ip+ij*n]=x1+dk;
			dif(t,n,&y[ij*n],&dy[ij*n]); *nstp=(*nstp)+1;
			for(j=0; j<n; ++j) {
				wk[j]=y[j+ij*n]-(48*y[j+ijm1*n]-36*y[j+ijm2*n]+
					16*y[j+ijm3*n]-3*y[j+ijm4*n]+12*h*dy[j+ij*n])/25.0;
				wm[n*j+ip]=(wk[j]-wk1[j+2*n])/dk;
			}
			y[ip+ij*n]=x1;
		}

		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) wk[j+i*n]=0.0;
			wk[i+i*n]=1.0;
		}
		iflg=0;
		iwk=(int *) calloc((size_t) n,sizeof(int));

/*	Calculate inverse of Jacobian */
		ier=gauelm(n,n,wm,wk,&det,iwk,n,&iflg);
		free(iwk);
		if(ier>0) {free(wm); return 729;}

/*	Reset the flag so that Jacobian is not calculated again */
		*iflag=1;
	}

	for(i=0; i<n; ++i) {
		wm[i]=wk1[i+2*n];
		wm[i+n]=y[i+ij*n];
	}

/*	Iteration on parameters using Broyden's method */
	for(ipas=1; ipas<=nit; ++ipas) {

/*	The convergence check */
		err=0.0;
		for(j=0; j<n; ++j) {
			s1=0.0;
/*	Correction to the Jth component */
			for(k=0; k<n; ++k) s1=s1+wk[j+k*n]*wm[k];
			r2=fabs(y[j+ij*n])+fabs(y[j+ij*n]-y[j+ijm1*n])+eps;
			if(fabs(s1/r2)>err) err=fabs(s1/r2);
			y[j+ij*n]=wm[j+n]-s1;
			wm[j+2*n]=-s1;
		}

		dif(t,n,&y[ij*n],&dy[ij*n]); *nstp=(*nstp)+1;
		if(err < cfac*reps) {free(wm); return 0;}

/*	calculate the residuals */
		for(j=0; j<n; ++j) wk1[j+2*n]=y[j+ij*n]-(48*y[j+ijm1*n]-36*y[j+ijm2*n]+
				16*y[j+ijm3*n]-3*y[j+ijm4*n]+12*h*dy[j+ij*n])/25.0;

/*	Update J inverse using Broyden's formula */
		for(j=0; j<n; ++j) {
			wm[j+3*n]=wk1[j+2*n]-wm[j];
			wm[j]=wk1[j+2*n];
			wm[j+n]=y[j+ij*n];
		}
		ss=0.0;
		for(j=0; j<n; ++j) {
			s1=0.0; s2=0.0;
			for(k=0; k<n; ++k) {
				s1=s1+wm[k+2*n]*wk[k+j*n];
				s2=s2+wk[j+k*n]*wm[k+3*n];
			}
			wk1[j+2*n]=s1;
			y[j+ij*n]=s2-wm[j+2*n];
			ss=ss+s1*wm[j+3*n];
		}

		if(ss==0.0) {free(wm); return 729;}

		for(j=0; j<n; ++j) {
			for(k=0; k<n; ++k) wk[k+j*n]-y[k+ij*n]*wk1[j+2*n]/ss;
		}
	}

	free(wm);
	return 729;
}



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





/*	To solve initial value problems in ordinary differential equations
	using a fourth-order multistep method with adaptive step size control
	It can be used with ADAMS (Adams-Bashforth-Moulton predictor corrector method)
	or with GEAR (Gear's stiffly stable method)

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length 7N, the first N elements
		should contain the initial values of variables. During
		execution, the solution at various intermediate points will
		be stored in this array. The contents of this array must
		be preserved between two calls to MSTEP if the integration
		needs to be continued further.
	DY : (output) Array of length 7N, containing the derivatives
		of Y the points stored in Y. This array should also be
		preserved between two calls to MSTEP.
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input/output) Initial guess for the step size. After execution
		it will contain the step size used by the function
	T0 : (input/output) Initial value of independent variable t, at
		which initial values are specified. After execution it will
		be set to the point up to which integration has been successful
	TN : (input) The final value of t at which the solution is required.
		Intermediate values will not be preserved so if solution
		is required at intermediate points, TN must be set to first
		such value and multiple calls will be needed to calculate
		all required values. For each subsequent call only TN needs
		to be updated. Other variables including the scratch arrays
		in the call statement must be preserved to their old values.
	YF : (output) Array of length N, containing the solution at
		the last point up to which integration is successful.
		If execution is successful it will contain the required
		solution at TN.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls local truncation error and hence
		actual error could be larger
	NSTP : (input/output) Number of calls to DIF required since the starting
		of integration. The count is updated by the function.
	NMAX : (input/output) Maximum number of calls to DIF to be tried.
		If NMAX<=0 it will be set to a default value of NMX=100000.
	IFLG : (input/output) Integer variable used as flag to decide the
		type of computation.
		If IFLG=0 the computation is started afresh with starting
			values being generated and the step size is adjusted
			to meet the specified accuracy. NSTP is also initialised
			to zero in the beginning. After execution IFLG is set to 2
		If IFLG=1 the computation is started afresh with starting
			values being generated but the step size is kept fixed.
			NSTP is also initialised to zero in the beginning.
			After execution IFLG is set to 3
		If IFLG=2 the integration is continued using already available
			solution at previous steps. The step size is adjusted
			to meet the specified accuracy.
		If IFLG=3 the integration is continued using already available
			solution at previous steps. The step size is kept fixed.
		IFLG should not be changed after the first call, unless
			fresh starting values need to be generated.
	IST : (input/output) Integer variable used as flag to decide the
		multistep method to be used:
		If IST=0 Adams method is used, which may be good for
			general equations
		If IST=1 Gear's method is used, which may be used for
			stiff equations
	WK : Array used as scratch space. For use with ADAMS the size
		of WK should be 2N, while for use with GEAR it should be
		3N+N*N. This scratch space must be preserved
		between two calls to MSTEP.
		
	Error status is returned by the value of the function MSTEP.
		0 value implies successful execution
		702 implies N<=0, no calculations are done
		724 or 725 implies STRT4 failed to generate starting values
		726 implies that step-size has become smaller than
			REPS*|TN-T0|
		727 implies that step size is too small for arithmetic used
		728 implies that integration could not be completed in
			the specified number of steps.

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : ADAMS, GEAR, STRT4, RK4, GAUELM, DIF
*/

#include <math.h>

int strt4(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double *h, double *t, double reps, int iflg, double tstep, int *nstp,
	double *wk);
int gear(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
	int ijm3, int ijm4, int *iflag, double *wk1, double *wk);
int adams(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
	int ijm3, int ijm4, double *wk);

int mstep(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double *h, double *t0, double tn, double yf[], double reps, int *nstp,
	int *nmax, int *iflg, int ist, double *wk)

{
	static int i,j,i2,j1,j2,ij,ijm1,ijm2,ijm3,ijm4,iflag,ns,ier;
	int nmx=100000;
	double t,err,reps1,r1,r2,y1,y2,tx,del[4],tstep, eps=1.e-30;

	if(n<=0) return 702;
	if(*nmax<=0) *nmax=nmx;
	tstep=tn-(*t0);

	if(*iflg<=1) {
/*	Generate the starting values */
		if(tstep==0.0) return 0;
		*nstp=0;
		if(*h==0.0) *h=tstep/8.0;
/*	Change the sign of H if needed */
		if( (*h<0.0) == (tn>(*t0))) *h=-(*h);

		t=(*t0);
/*	Generate the starting values */
		ier=strt4(n,y,dy,dif,h,&t,reps,*iflg,tstep,nstp,wk);
		if(ier>0) return ier;
	
/*	Initialising the array indices */
		ij=3; ijm1=2; ijm2=1; ijm3=0;
		*t0=t;
		iflag=0;		 /*	Set the flag for function GEAR */
		*iflg=(*iflg)+2;     /*	Update IFLG */
		ns=4;
	}

	if((tn<=(*t0)) == (*h>0.0) ) goto interp;
	if(tstep==0.0) goto interp;

	while(*nstp<(*nmax)) {
/*	Updating the array indices */
		ijm4=ijm3; ijm3=ijm2; ijm2=ijm1; ijm1=ij;
		++ij; if(ij>6) ij=0;
		t=(*t0)+(*h);
		reps1=fabs(reps*(*h)/tstep);

/*	To perform one step of integration using Adams method */
		if(ist==0) ier=adams(n,y,dy,dif,*h,t,reps1,nstp,ij,ijm1,ijm2,ijm3,ijm4,wk);

/*	To perform one step of integration using stiffly stable method */
		else ier=gear(n,y,dy,dif,*h,t,reps1,nstp,ij,ijm1,ijm2,ijm3,ijm4,&iflag,wk,&wk[3*n]);

/*	Estimating the truncation error */
		err=0.0;
		for(i=0; i<n; ++i) {
			r2=fabs(y[i+ij*n])+fabs(y[i+ij*n]-y[i+ijm1*n])+eps;
			r1=fabs(wk[i]-y[i+ij*n])/r2;
			if(r1>err) err=r1;
		}
		err=fabs(0.25*err*tstep/(*h));

		if(t==(*t0)) {
/*	Step size is too small for the arithmetic used */
			for(i=0; i<n; ++i) yf[i]=y[i+ij*n];
			return 727;
		}

		if(*iflg<=2) {
			if(err<reps && ier==0) {
/*	Integration is successful */
				*t0=(*t0)+(*h);
				for(i=0; i<n; ++i) wk[i+n]=wk[i];
				ns=ns+1;
				if((tn<=(*t0)) == (*h>0.0)) goto interp;
				if(fabs((tn-(*t0))/tstep) <= reps) goto interp;

				if(err<reps/32.0 && ns>7) {
/*	Double the step size */
					*h=2.*(*h);
					i2=ijm4-2; if(i2<0) i2=i2+7;
					for(i=0; i<n; ++i) {
						y[i+ijm1*n]=y[i+ijm2*n];
						y[i+ijm2*n]=y[i+ijm4*n];
						y[i+ijm3*n]=y[i+i2*n];
						dy[i+ijm1*n]=dy[i+ijm2*n];
						dy[i+ijm2*n]=dy[i+ijm4*n];
						dy[i+ijm3*n]=dy[i+i2*n];
					}
					ns=4;
				}
			}

			else {
/*	If integration has failed, then reduce the step size */
				if(ns>-4 && err<16*reps) {
/*	H is halved */
					for(i=0; i<n; ++i) {
						y1=(45.*y[i+ijm1*n]+72.*y[i+ijm2*n]+11.*y[i+ijm3*n])/128.+
							(*h)*(-9.*dy[i+ijm1*n]+36*dy[i+ijm2*n]+3*dy[i+ijm3*n])/128.;
						y2=(11.*y[i+ijm1*n]+72.*y[i+ijm2*n]+45.*y[i+ijm3*n])/128.+
							(*h)*(-3.*dy[i+ijm1*n]-36*dy[i+ijm2*n]+9*dy[i+ijm3*n])/128.;
						y[i+ij*n]=y[i+ijm1*n];
						y[i+ijm1*n]=y1;
						y[i+ijm3*n]=y2;
						dy[i+ij*n]=dy[i+ijm1*n];
					}
					t=(*t0)-(*h)/2;
					dif(t,n, &y[ijm1*n], &dy[ijm1*n]);
					t=t-(*h);
					dif(t,n, &y[ijm3*n], &dy[ijm3*n]);
					*nstp=(*nstp)+2;
					*h=(*h)/2.0;
				}

				else {
/*	If error is too large or the halving has failed once, then
	generate fresh starting values with smaller h */
					*iflg=(*iflg)-2;
					t=(*t0); *h=(*h)/8.;
					for(i=0; i<n; ++i) y[i]=y[i+ijm1*n];
					ier=strt4(n,y,dy,dif,h,&t,reps,*iflg,tstep,nstp,wk);
					if(ier>0) {
/*			If STRT4 fails, then quit */
						for(i=0; i<n; ++i) yf[i]=y[i+ijm1*n];
						return ier;
					}
					ij=3; ijm1=2; ijm2=1; ijm3=0;
					*t0=t; *iflg=(*iflg)+2;
					iflag=0;
				}

				ns=-4;
				if(fabs((*h)/tstep)<reps) {
/*	The step size is too small then quit */
					for(i=0; i<n; ++i) yf[i]=y[i+ijm1*n];
					return 726;
				}
			}
		}

		else {
			if(ier>0) {
/*	For integration with fixed H, quit if corrector fails */
				for(i=0; i<n; ++i) yf[i]=y[i+ijm1*n];
				return ier;
			}
			*t0=(*t0)+(*h);
			for(i=0; i<n; ++i) wk[i+n]=wk[i];
			if((tn<=(*t0)) == (*h>0.0)) goto interp;
			if(fabs((tn-(*t0))/tstep)<=reps) goto interp;
		}

	}

/*	Quit if the specified number of function evaluations are exceeded */
	for(i=0; i<n; ++i) yf[i]=y[i+ij*n];
	return 728;

interp:	tx=(tn-(*t0))/(*h);
/*	Use interpolation to calculate solution at TN */
		for(i=0; i<n; ++i) {
			del[0]=dy[i+ij*n]*(*h);
			del[1]=y[i+ij*n]-y[i+ijm1*n];
			del[2]=dy[i+ijm1*n]*(*h);
			del[3]=y[i+ijm1*n]-y[i+ijm2*n];
			for(j1=1; j1<3; ++j1) {
				for(j2=3; j2>=j1; --j2) del[j2]=del[j2-1]-del[j2];
			}
			del[3]=del[2]-0.5*del[3];
			yf[i]=y[i+ij*n]+tx*(del[0]+tx*(del[1]+(tx+1.0)*(del[2]+(tx+1.0)*del[3]/2.)));
		}
		return 0;
}


/*	To perform one step of integration of ordinary differential equations
	using a fourth-order Runge Kutta method 

	N : (input) Number of first order differential equations to be solved
	T : (input) Initial value of independent variable t, at
		which initial values are specified. This value is not updated.
	Y0 : (input) Array of length N containing the initial
		values of variables.
	DY0 : (input) Array of length N containing the derivatives
		of Y at the initial point Y0
	H : (input) The step size to be used for integration
	Y1 : (output) Array of length N containing the solution at t=T+H
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : DIF
*/	

#include <math.h>
#include <stdlib.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ))

{
	int i;
	double h2,t1;
	double *wk1, *wk2;

	h2=h/2.;
	wk1=(double *) calloc((size_t) n, sizeof(double));
	wk2=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i) {
		y1[i]=h2*dy0[i];
/*	y_0+0.5k_1 */
		wk1[i]=y0[i]+y1[i];
	}

	t1=t+h2;
	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i) {
		y1[i]=y1[i]+h*wk2[i];
/*	y_0+0.5k_2 */
		wk1[i]=y0[i]+h2*wk2[i];
	}

	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i) {
		y1[i]=y1[i]+h*wk2[i];
/*	y_0+k_3 */
		wk1[i]=y0[i]+h*wk2[i];
	}

	t1=t+h;
	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i)  y1[i]=y0[i]+(y1[i]+h2*wk2[i])/3.0;
	free(wk2); free(wk1);
	return;
}
	
	
	


/*	To generate the starting values for multistep methods in solution
	of initial value problem in ordinary differential equations.
	It is used by MSTEP.

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length 4N, the first N elements
		should contain the initial values of variables. During
		execution, the solution at the next three points will
		be stored in this array. 
	DY : (output) Array of length 4N, containing the derivatives
		of Y the points stored in Y.
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input/output) Initial guess for the step size. After execution
		it will contain the step size used by the function
	T : (input/output) Initial value of independent variable t, at
		which initial values are specified. After execution it will
		be set to T+3H if execution is successful.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls local truncation error and hence
		actual error could be larger
	IFLG : (input) Integer variable used as flag to decide the
		type of computation.
		If IFLG=0 the step size is adjusted to meet the specified accuracy.
		If IFLG=1 the step size is kept fixed.
	TSTEP : (input) The size of interval over which integration is requested
			It is used only for convergence check.
	NSTP : (input/output) Number of calls to DIF required since the starting
		of integration. The count is updated by the function.
	WK : Scratch array of length 2N used to pass on the value for modifier
		
	Error status is returned by the value of the function STRT4.
		0 value implies successful execution
		724 implies that step size becomes too small
		725 implies that iteration to generate starting values
			failed to converge

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : RK4, DIF
*/

#include <math.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ));

int strt4(int n, double *y, double *dy, void dif(double , int , double * , double * ),
	double *h, double *t, double reps, int iflg, double tstep, int *nstp,
	double *wk)

{
	int i,it, nit=10;
	double t0,h2,err,r1,r2, sfac=0.9, eps=1.0e-30;

	t0=(*t);
	for(it=1; it<=nit; ++it) {
		dif(t0,n,y,dy);
/*	Generate y_1 */
		rk4(n,t0,y,dy,*h,&y[n],dif);
		*t=t0+(*h);
		dif(*t,n,&y[n],&dy[n]);
/*	Generate y_2 */
		rk4(n,*t,&y[n],&dy[n],*h,&y[2*n],dif);

		h2=2.0*(*h);
/*	Calculate y_2 using double step */
		rk4(n,t0,y,dy,h2,&y[3*n],dif);
		*nstp=*nstp+11;

/*	Estimate the truncation error in y_2 */
		err=0.0;
		for(i=0; i<n; ++i) {
			r2=fabs(y[i])+fabs(y[i+n])+fabs(y[i+2*n])+eps;
			r1=fabs(y[i+2*n]-y[i+3*n])/r2;
			if(r1>err) err=r1;
		}
		err=err*tstep/(*h);

		if(err<=reps || iflg>0) {
/*	Accept the computed values */
			*t=(*t)+(*h);
			dif(*t,n,&y[2*n],&dy[2*n]);
/*		Generate y_3 */
			rk4(n,*t,&y[2*n],&dy[2*n],*h,&y[3*n],dif);
			*t=(*t)+(*h);
			dif(*t,n,&y[3*n],&dy[3*n]);
			*nstp=(*nstp)+5;
/*	Store Y_3 for use by modifier */
			for(i=0; i<n; ++i) wk[i+n]=y[i+3*n];
			return 0;
		}

		else {
/*	Reduce the step size */
			*h=sfac*(*h)*pow(reps/err,0.25);
			if(fabs(*h/tstep) < reps || (*t)==t0) return 724;
		}
	}

	return 725;
}
