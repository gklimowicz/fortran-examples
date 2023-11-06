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
