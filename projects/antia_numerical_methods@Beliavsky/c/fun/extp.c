/*	To solve initial value problems in ordinary differential equations
	using extrapolation method

	N : (input) Number of first order differential equations to be solved
	Y : (input/output) Array of length N containing the initial
		values of variables. After execution it will contain the
		values of variable at the last point where the integration
		has been successful.
	DY : (output) Array of length N containing the derivatives
		of Y at the last point
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)
	H : (input/output) Initial guess for the step size. After execution
		it will contain the step size used by the function
	T0 : (input/output) Initial value of independent variable t, at
		which initial values are specified. After execution it will
		be set to the point up to which integration has been successful
	TN : (input) The final value of t at which the solution is required.
		If integration is successful T0 will be set equal to TN.
		Intermediate values will not be preserved so if solution
		is required at intermediate points, TN must be set to first
		such value and multiple calls will be needed to calculate
		all required values. For each subsequent call only TN needs
		to be changed.
	REPS : (input) Required accuracy in each component of the solution.
		The function only controls local truncation error and hence
		actual error could be larger
	NSTEP : (output) Number of calls to DIF made by the function to
		complete the integration during each call to EXTP.
	NMAX : (input/output) Maximum number of function evaluations to be used.
		If NMAX<=0 it will be set to a default value of NMX=100000.
	IFLG : (input) Integer variable used as a flag to decide the type
		of extrapolation to be used.
		If IFLG=0 polynomial extrapolation is used
		otherwise rational function extrapolation is used
		
	Error status is returned by the value of the function EXTP.
		0 value implies successful execution
		703 implies N<=0, no calculations are done
		730 implies that step-size has become smaller than REPS*|TN-T0|
		731 implies that step size is too small for arithmetic used
		732 implies that integration could not be completed in
			the specified number of steps.
		733 implies that denominator for evaluating rational
			function extrapolation vanished

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : DIF
*/	

#include <math.h>
#include <stdlib.h>

int extp(int n, double y[], double dy[], void dif(double , int , double * , double * ),
	double *h, double *t0, double tn, double reps, int *nstep, int *nmax,
	int iflg)

{
	int i,j,ie,ier, nmx=100000;
	int nseq[12] = {2,4,6,8,12,16,24,32,48,64,96,128};
	double tstep,h2,t1,c1,cn,den,err,r1,r2,hi[12];
	double sfac=0.80, eps=1.0e-30;
	double *wk;

	if(n<=0) return 703;
	if(*nmax<=0) *nmax=nmx;

	*nstep=0;
	tstep=tn-(*t0);
	if(tstep==0.0) return 0;
	if(*h==0.0 || fabs(*h)>fabs(tstep)) *h=tstep;
	if((*h<0.0) == (tn>(*t0))) *h=-(*h);
	wk=(double *) calloc((size_t) (39*n), sizeof(double));
	ier=0;

/*	Loop for integration */
	while((*nstep)<(*nmax)) {
		dif(*t0,n,y,dy); *nstep=(*nstep)+1;

/*	Loop for extrapolation */
		for(ie=0; ie<12; ++ie) {
/*	Step size for the midpoint method */
			h2=(*h)/nseq[ie];
			hi[ie]=h2*h2;
			for(i=0; i<n; ++i) {
				wk[i]=y[i];
				wk[i+n]=y[i]+h2*dy[i];
			}

			t1=(*t0);
			for(j=0; j<nseq[ie]-1; ++j) {
				t1=t1+h2;
				dif(t1,n,&wk[n],&wk[2*n]);
				for(i=0; i<n; ++i) {
/*		The midpoint rule */
					c1=wk[i]+2.0*h2*wk[i+2*n];
					wk[i]=wk[i+n];
					wk[i+n]=c1;
				}
			}
			t1=t1+h2;
			dif(t1,n,&wk[n],&wk[2*n]);
			*nstep=(*nstep)+nseq[ie];
			for(i=0; i<n; ++i) {
/*		Modified midpoint method */
				wk[i+(3+ie)*n]=0.5*(wk[i]+wk[i+n]+h2*wk[i+2*n]);
				wk[i+(15+ie)*n]=wk[i+(3+ie)*n];
				wk[i+(27+ie)*n]=wk[i+(3+ie)*n];
			}

			if(ie>0) {
/*	Perform extrapolation */
				for(j=ie-1; j>=0; --j) {
					for(i=0; i<n; ++i) {
						if(iflg==0) {
/*		Use polynomial extrapolation */
							cn=(wk[i+(4+j)*n]-wk[i+(15+j)*n])/(hi[j]-hi[ie]);
							wk[i+(15+j)*n]=hi[ie]*cn;	/* D[j,ie-j] */
							wk[i+(3+j)*n]=hi[j]*cn;		/* C[j,ie-j] */
						}

						else{
/*	Use rational function extrapolation */
							den=hi[j]*wk[i+(15+j)*n]/hi[ie]-wk[i+(4+j)*n];
							if(den==0.0) {ier=733; goto next;}
							cn=(wk[i+(4+j)*n]-wk[i+(15+j)*n])/den;
/*								C(J,IE-J) */
							wk[i+(3+j)*n]=hi[j]*wk[i+(15+j)*n]*cn/hi[ie];
/*								D(J,IE-J) */
							wk[i+(15+j)*n]=wk[i+(4+j)*n]*cn;
						}
						wk[i+(27+j)*n]=wk[i+(3+j)*n]+wk[i+(27+j)*n];
					}
				}

/*	Estimating the truncation error */
				err=0.0;
				for(i=0; i<n; ++i) {
					r2=fabs(y[i])+fabs(wk[i+27*n]-y[i])+eps;
					r1=fabs(wk[i+27*n]-wk[i+28*n])/r2;
					if(r1>err) err=r1;
				}
				err=fabs(err*tstep/(*h));
				if(err<reps) break;
			}
		}

next:	if(t1==(*t0)) {
/*	The step size is too small for the arithmetic used */
			if(ier==0) ier=731;
			free(wk); return ier;
		}

		if(err<reps) {
/*	Integration is successful */
			*t0=(*t0)+(*h);
			for(i=0; i<n; ++i) y[i]=wk[i+27*n];
			if(fabs((tn-(*t0))/tstep)<reps) {free(wk); return ier;}

			if(ie<6) *h=sfac*nseq[6]*(*h)/nseq[ie]; /*	Increase the step size */
			else if(ie > 6) *h=sfac*(*h); /*	Decrease the step size */
			if((*t0+(*h)>tn) == (*h>0.0)) *h=tn-(*t0);
		}

/*	If the integration has failed, then decrease the step size */
		else *h=(*h)/32.0;

		if(fabs((*h)/tstep)<reps) {
/*	The step size is too small then quit */
			if(ier==0) ier=730;
			free(wk); return ier;
		}
	}

	free(wk);
	return 732;
}
