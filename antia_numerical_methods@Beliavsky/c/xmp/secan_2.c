/*	Real roots of nonlinear equation,  f(x) = f * 2**jf =0 
	using secant iteration */

#include <stdio.h>
#include <math.h>

double det(double x, int * j);
int secan2(double x0, double xl, double xu, double *x, double reps,
	double aeps, double (*fun) (double , int *));
int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg);
void secani(double x0, double xl, double xu, double *x, double f, int jf,
	double reps,double aeps, int *ier);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np,nmax;
	double xl, xu, root, reps, aeps, x0, f0, fb[20];

/*	Exercise 7.44 : Eigenvalues of Hilbert matrix */

	aeps=1.e-27; reps=1.e-9;
	for(i1=0; i1<99; ++i1) {
		printf("type xl=lower limit,  xu=upper limit,   x0=initial guess \n");
		printf("                           (quits  xl=xu)\n");
		scanf(" %le %le %le", &xl, &xu, &x0);
		if(xl==xu) return 0;

		i=secan_2(x0,xl,xu,&root,reps,aeps,det);
		printf(" ier = %d   initial guess = %e   root = %e   \n", i,x0,root);

/*	Using secani with reverse communication technique */
		ier=0;
		do {
			secani(x0,xl,xu,&root,f0,id,reps,aeps,&ier);
			if(ier>=0) break;
			f0=det(root,&id);
		} while(ier<0);
		printf(" ier = %d   lower limit = %e  upper limit = %e  root = %e   \n", ier,xl,xu,root);

	}
	return;
}


/*	Function to calculate the determinant of Hilbert matrix */

double det(double x, int *jf)

{
	int i,j,n, num,lj,iflg,intc[20];
	double cd, a[20][20], xr[20];

	n=10; num=1; lj=20; iflg=1;

/*	set up the Hilbert matrix  H - x*I  */
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) a[i][j]=1.0/(i+j+1.0);
		a[i][i]=a[i][i]-x;
	}

	i=crout(n,num,&a[0][0],&xr[0],&cd,jf,intc,lj,&iflg);

	return cd;
}



/*	Solution of a system of linear equations using Crout's algorithm 
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
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
			interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter which determines the type of computation
		required.
		If IFLG<=0, both elimination and solution are calculated
			and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function CROUT
		0 value implies successful execution
		102 implies (N<0 or N>LJ) 
		122 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg)

{
	int i,j,k,km,l;
	double r1,t1,r2,d1;
	double *wk;

	if(n<=0 || n>lj) return 102;

	if((*iflg)<2)
/*	Perform LU decomposition  */
	{wk=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i){
		r1=0.0;
		for(j=0; j<n; ++j){
			if(fabs(a[i*lj+j]) > r1) r1=fabs(a[i*lj+j]);}
		wk[i]=r1;
		if(r1==0.0) {free(wk); return 122;}
	}

	*det=1.0; *idet=0;
	for(k=0; k<n; ++k){
		r1=0.0; km=k;
/*	Generate the Kth column of L  */
		for (l=k; l<n; ++l){
			d1=a[l*lj+k];
			for(i=0; i<=k-1; ++i) d1=d1-a[l*lj+i]*a[i*lj+k];
			a[l*lj+k]=d1;
/*	Finding the pivot  */
			r2=fabs(d1/wk[l]);
			if(r2>r1) {r1=r2; km=l;}
		}

		inc[k]=km;
/*	Interchange the rows if needed  */
		if(km != k) {
			*det=-(*det);
			for(l=0; l<n; ++l) {
				t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;
			}
			t1=wk[k]; wk[k]=wk[km]; wk[km]=t1;
		}

		*det=*det*a[k*lj+k];
		if(a[k*lj+k]==0.0) {free(wk); return 122;}

		if(*det != 0.0) {
/*	Scale the value of the determinant DET */
			while(fabs(*det)>32.0) {
				*det=(*det)*0.03125; *idet=*idet+5;}

			while(fabs(*det)<0.03125) {
				*det=(*det)*32.0; *idet=*idet-5;}
		}

/*	Generate the Kth row of U */
		for(l=k+1; l<n; ++l) {
			d1=a[k*lj+l];
			for(i=0; i<=k-1; ++i) d1=d1-a[k*lj+i]*a[i*lj+l];
			a[k*lj+l]=d1/a[k*lj+k];
		}
	}
	free(wk);

	if(*iflg==1) {*iflg=2; return 0;}
	*iflg=2;
	}

/*	Solution for NUM different right-hand sides  */
	for(j=0; j<num; ++j){
/*	Forward substitution  */
		for(k=0; k<n; ++k) {
			if(k != inc[k]) {
				t1=x[j*lj+k];
				x[j*lj+k]=x[j*lj+inc[k]];
				x[j*lj+inc[k]]=t1;
			}

			d1=x[j*lj+k];
			for(l=0; l<=k-1; ++l) d1=d1-a[k*lj+l]*x[j*lj+l];
			x[j*lj+k]=d1/a[k*lj+k];
		}

/*	Back-substitution  */
		for(k=n-2; k>=0; --k){
			d1=x[j*lj+k];
			for(l=n-1; l>=k+1; --l) d1=d1-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=d1;
		}
	}
	return 0;
}




/*	Real zero of a given function using secant iteration
	Function is calculated as FX*2**JF

	X0 : (input) Initial guess for the zero
	XL : (input) Lower limit of interval where zero is expected
	XU : (input) Upper limit of interval where zero is expected
	X : (output) Computed value of the zero
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))
	FUN : (input) Name of the function routine to calculate the function
		
	Error status is returned by the value of the function SECAN_2.
		0 value implies successful execution
		40 implies that function value is equal at two points
			and it is not possible to continue the iteration
		402 implies XL>X0 or XU<X0, in which case no calculations are done
		422 implies that iteration goes outside the specified limits
		423 implies that iteration failed to converge to specified accuracy

	Function FUN(X,JF) must be supplied by the user.
		The function value should be FUN*2**JF

	Required functions : FUN
*/

#include <math.h>

int secan_2(double x0, double xl, double xu, double *x, double reps,
	double aeps, double (*fun) (double , int * ))

{
	int l,jf,jf1, nis=75;
	double f,f1,r1,dx,dx1;

	if(xl>x0 || xu<x0) return 402;

	*x=x0;
/*	Select the increment for the next point X+DX */
	dx=(xu-x0)/100.0;
	r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
	if(fabs(dx)<5.0*r1) dx=(xl-x0)/5.;
	r1=fabs(*x); if(r1<100.0*aeps) r1=100.0*aeps;
	if(fabs(dx)>0.1*r1) {
		if(dx>=0.0) dx=r1;
		else dx=-r1;
	}
	f1=0.0; jf1=0;

	for(l=1; l<=nis; ++l) {
		f=fun(*x,&jf);
		dx1=dx;
		f1=f1*pow(2.0, (double) (jf1-jf));

		if(f1-f == 0.0) {
			if(f == 0.0) return 0;
/*	If F1=F and F!=0, then quit */
			else return 40;
		}

/*	The secant iteration */
		if(l>1) dx=dx1*f/(f1-f);
		*x=(*x)+dx;
		f1=f; jf1=jf;

		r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
		if(fabs(dx)<r1 && l>2) return 0;
		if(*x<xl || (*x)>xu) return 422;
	}

	return 423;
}



/*      Real zero of a given function using secant iteration
	Function is calculated as FX*2**JF
	This function uses reverse communication to calculate function
	values. If IER<0 the function should be evaluated and SECANI should
	be called back with new function value. Calculation of function
	value should not change any other variables in the call statement.

	X0 : (input) Initial guess for the zero
	XL : (input) Lower limit of interval where zero is expected
	XU : (input) Upper limit of interval where zero is expected
	X : (output) Value of x at which the function evaluation is required.
		If IER=0 then it will contain the final value of zero computed
		by the function.
	F : (input) Calculated value of the function at X.
		If function exits with IER<0, then the calling function should
		calculate the function value at X and call SECANI with this value
		stored in F and JF. Other variables should not be changed.
	JF : (input) The exponent of function value, the function value
		should be F*2**JF
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
     		The estimated error should be less than max(AEPS, REPS*fabs(X))
	IER : (input/output) Error parameter, IER=0 implies successful execution
		Before the first call IER should be set to zero
		IER<0 implies that execution is not over and the function needs
			a new function evaluation at X. After calculating the
			function value SECANI should be called back.
     		IER=40 implies that function value is equal at two points
     			and it is not possible to continue the iteration
     		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
     		IER=422 implies that iteration goes outside the specified limits
     		IER=423 implies that iteration failed to converge to specified accuracy


	Required functions : None (Function value is calculated by the calling program)
*/

#include <math.h>


void secani(double x0, double xl, double xu, double *x, double f, int jf,
	double reps,double aeps, int *ier)

{
	int nis=75;
	static int l,jf1;
	static double f1,r1,dx,dx1;

	if(*ier==-1) goto func;

	if(xl>x0 || xu<x0) {*ier=402; return;}

	*x=x0;
/*	Select the increment for the next point X+DX */
	dx=(xu-x0)/100.0;
	r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
	if(fabs(dx)<5.0*r1) dx=(xl-x0)/5.;
	r1=fabs(*x); if(r1<100.0*aeps) r1=100.0*aeps;
	if(fabs(dx)>0.1*r1) {
		if(dx>=0.0) dx=r1;
		else dx=-r1;
	}
	f1=0.0; jf1=0; l=0;

loop:	l=l+1; *ier=-1;
/*	To evaluate the function at x */
		return;

func:	dx1=dx;
		f1=f1*pow(2.0, (double) (jf1-jf));

		if(f1-f == 0.0) {
			if(f == 0.0) {*ier = 0; return ;}
/*	If F1=F and F!=0, then quit */
			else {*ier=40; return;}
		}

/*	The secant iteration */
		if(l>1) dx=dx1*f/(f1-f);
		*x=(*x)+dx;
		f1=f; jf1=jf;

		r1=reps*fabs(*x); if(aeps>r1) r1=aeps;
		if(fabs(dx)<r1 && l>2) {*ier=0; return;}
		if(*x<xl || (*x)>xu) {*ier=422; return;}
	if(l<nis) goto loop;

	*ier=423; return;
}
