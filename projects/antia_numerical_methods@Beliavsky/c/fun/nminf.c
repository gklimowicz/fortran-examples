/*	To minimise a function of several variables using direction set method

	N : (input) Number of variables
	X : (input/output) Array of length N containing the initial
		guess for the minimum.
		After execution it should contain the coordinates of minimiser
	F : (output) The function value at X
	NUM : (output) Number of function evaluations used by NMINF
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy, iteration will stop when
		change in function value is less than max(AEPS, REPS*fabs(F))
	FCN : (input) Name of the function to calculate the function value
		
	Error status is returned by the value of the function NMINF.
		0 value implies successful execution
		504 implies that N < 2, in which case no calculations are done
		529 implies that iteration failed to converge to specified accuracy
		Other values may be set by LINMNF

	Function FCN(N,X,F) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, F is the
		function value at X. X is an array of length N.

	Required functions : LINMNF, FLN, SVD, FCN
*/

#include <math.h>
#include <stdlib.h>

int linmnf(double *x0, double *x1, double *f0, double reps, double aeps,
	void (*f) (int , double * , double * ), double v[], double xi[],
	int n, int *num);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);

int nminf(int n, double x[], double *f, int *num, double reps, double aeps,
	void (*fcn) (int , double * , double * ))

{
	int i,j,it,k,kmax,ier,ier1, nit=200;
	double f0,fi,dfmax,x1,x2,fun,r1;
	double *wk, *wk1, *u;

	if(n<2) return 504;
	wk=(double *) calloc((size_t) (n*n),sizeof(double));
	wk1=(double *) calloc((size_t) (4*n),sizeof(double));

/*	Initialise the direction set matrix to identity matrix */
	for(i=0; i<n; ++i) {
		wk1[i+2*n]=0.0;
		wk1[i+3*n]=0.0;
		for(j=0; j<n; ++j) wk[j+i*n]=0.0;
		wk[i+i*n]=1.0;
	}

	fcn(n,x,f);
	*num=1;
	ier=0;

/*	The main iteration loop */
	for(it=1; it<=nit; ++it) {
		fi=(*f);
		for(k=0; k<n; ++k) {

/*	The starting point for line search */
			for(i=0; i<n; ++i) wk1[i]=x[i];
			f0=(*f);
			dfmax=0.0; kmax=0;

			for(i=0; i<n; ++i) {
				x1=0.0;
/*	Use previous value as initial approximation to minimum */
				x2=wk1[i+2*n];
				fun=f0;
				if(x2==0.0) x2=1.0;
				ier1=linmnf(&x1,&x2,&f0,reps,aeps,fcn,&wk[i*n],wk1,n,num);
				if(ier1>0) ier=ier1;

				wk1[i+2*n]=x1;
/*	Estimate of second derivative along the line */
				if(x1 != 0.0) wk1[i+3*n]=fabs((fun-f0)/(x1*x1));
/*	The new starting point */
				for(j=0; j<n; ++j) wk1[j]=wk1[j]+x1*wk[j+i*n];
				if(fun-f0>=dfmax && i<=n-k-1) {dfmax=fun-f0; kmax=i;}
			}

/*	Remove the KMAX th direction */
			for(i=kmax; i<n-1; ++i) {
				wk1[i+2*n]=wk1[i+1+2*n];
				wk1[i+3*n]=wk1[i+1+3*n];
				for(j=0; j<n; ++j) wk[j+i*n]=wk[j+i*n+n];
			}
/*	Add a new direction */
			for(i=0; i<n; ++i) wk[i+n*n-n]=wk1[i]-x[i];

			x1=0.0; x2=1.0;
			fun=(*f);
			wk1[4*n-1]=0.0;
/*	Starting point for the final line search in the loop */
			for(i=0; i<n; ++i) wk1[i]=x[i];
			ier1=linmnf(&x1,&x2,f,reps,aeps,fcn,&wk[n*n-n],wk1,n,num);
			if(ier1>0) ier=ier1;

			wk1[3*n-1]=x1;
			if(x1!=0.0) wk1[4*n-1]=fabs((fun-(*f))/(x1*x1));
			for(j=0; j<n; ++j) x[j]=x[j]+x1*wk[j+n*n-n];
		}

		r1=reps*fabs(*f); if(aeps>r1) r1=aeps;
		if(fabs(*f-fi)<r1) {free(wk1); free(wk); return ier;}

/*	The matrix V for SVD */
		for(j=0; j<n; ++j) {
			if(wk1[j+3*n]>0.0) {
				for(i=0; i<n; ++i) wk[i+j*n]=wk[i+j*n]/sqrt(wk1[j+3*n]);
			}
		}
/*	Take the transpose of matrix as the array is stored in different order */
		for(i=0; i<n; ++i){
			for(j=i+1; j<n; ++j)
				{r1=wk[i*n+j]; wk[i*n+j]=wk[i+j*n]; wk[i+j*n]=r1;}
		}

		u=(double *)calloc((size_t) (n*n),sizeof(double));
		ier=svd(n,n,wk,u,wk1,n,n);
		free(u);
		if(ier>0) {free(wk1); free(wk); return ier;}
/*	Transpose the matrix back to earlier order */
		for(i=0; i<n; ++i){
			for(j=i+1; j<n; ++j) {r1=wk[i*n+j]; wk[i*n+j]=wk[i+j*n]; wk[i+j*n]=r1;}
		}

		for(i=0; i<n; ++i) {
			wk1[i+2*n]=sqrt(fun-(*f))*wk1[i];
			wk1[i+3*n]=0.0;
			if(wk1[i]!=0.0) wk1[i+3*n]=1./(wk1[i]*wk1[i]);
		}
	}

	free(wk1); free(wk);
	return 529;
}
