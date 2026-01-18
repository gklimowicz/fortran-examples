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
