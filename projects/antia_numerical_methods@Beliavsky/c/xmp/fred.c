/*	To solve linear Fredholm equation using quadrature method */

#include <stdio.h>
#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int invit(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit);
int fred(int m, double a, double b, double wt[], double x[], double f[],
	double fc[], double fg(double ), double fker(double , double ),
	double *ei, int iq, int it, double reps);

double fg(double x);
double fker(double x, double t);

double al;

main()
{
	int i,i1,j,n,m, id, iflg, ier,np,nmax;
	double hh, x[100], wt[100],f[100],fc[100],tn,erc,reps,t0,fi;

/*	Example 13.1 */

	reps=1.e-8;
	t0=0.0; tn=1.0; iflg=2;
	for(i1=0; i1<99; ++i1) {
		printf("type n = no. of points,  id=1/2/4/8/16/32   quadrature formula\n");
		printf("                           (quits when n<=0)\n");
		scanf(" %d %d",&n,&id);
		if(n<=0) return 0;

		i=fred(n,t0,tn,wt,x,f,fc,fg,fker,&erc,id,iflg,reps);
		printf(" ier = %d    no. of points =  %d   id = %d \n",i,n,id);
		printf("    x         solution    corrected solution\n");
		for(i=0; i<n; ++i) printf(" %e %e %e \n",x[i],f[i],fc[i]);

/*	For Fredholm equation of second kind use quadrature formula to
	compute the solution at required points */

		printf("\n   x          solution \n");
		for(i=0; i<6; ++i) {
			hh=i*0.2; fi=-fg(hh);
			for(j=0; j<n; ++j) fi=fi+wt[j]*fker(hh,x[j])*f[j];
			printf(" %e  %e \n",hh,fi);
		}

	}
	return;
}

double fg(double x)

{
	return -x*x*x;
}

double fker(double x, double t)

{
	if(t>x) return -x*(1.0-t);
	else return -t*(1.0-x); 
}



/*	To solve linear Fredholm equation using quadrature method

	M : (input) Number of abscissas to be used in quadrature formula
	A : (input) Lower limit of the integral
	B : (input) Upper limit of the integral
	WT : (input/output) Array of length M containing the
		weights used in quadrature formula.
		If IQ is negative the weights must be supplied by the user
		otherwise they are calculated by the function.
	X : (input/output) Array of length M containing the
		abscissas used in quadrature formula.
		If IQ is negative the abscissas must be supplied by the user
		otherwise they are calculated by the function.
	F : (output) Array of length M containing the calculated solution
		F[I] is the value of solution at X[I].
	FC : (output) Array of length M containing the calculated solution
		after applying the deferred correction. This will be
		relevant only if IQ=1 and IT=1,2. In other cases a dummy
		array of any length may be supplied.
	FG : (input) Name of the function used to calculate the right
		hand side g(x). For IT=3 g(x) is not required, but in that
		case this function is used to calculate an initial guess
		for the eigenfunctions. In most cases the inverse iteration
		converges from essentially arbitrary initial guess and it
		is enough to set FG(X) to any nonzero value.
	FKER : (input) Name of the function used to calculate the kernel K(x,t)
	EI : (input/output) Initial guess for the eigenvalue. After execution
		it will contain the computed eigenvalue. Used only for IT=3
	IQ : (input) Integer variable to specify the quadrature formula to be used.
		If IQ=1, then trapezoidal rule is used and deferred correction
			is calculated using Gregory's formula
		If IQ=2, the Simpson's 1/3 rule is used
		If IQ=4,8,16,32 a composite rule using IQ point
			Gauss-Legendre formula is used.
		In all these cases the weights and abscissas are calculated
		If IQ is negative then it is assumed that weights and abscissas
		are supplied in arrays WT and X.
		Other values of IQ will cause an error return.
	IT : (input) Integer variable to specify the type of integral equation
		If IT=1 Fredholm equation of the first kind is solved
		If IT=2 Fredholm equation of the second kind is solved
		If IT=3 Fredholm equation of the third kind (eigenvalue
			problem) is solved
	REPS : (input) Required relative accuracy in calculating eigenvalue
		and eigenvectors. It is not used for IT=1,2. It only specifies
		the accuracy to which inverse iteration converges. It does not
		control the truncation error.
		
	Error status is returned by the value of the function FRED.
		0 value implies successful execution
		-11 implies that the calculations are actually performed
			using a smaller number of abscissas than M
		706 implies M<3, IT>3, IT<=0 or M is not sufficient
			to apply the required quadrature formula. No calculations
			are done.
		707 implies that IQ is not acceptable and no calculations
			are done
		Other values of may be set by GAUELM and INVIT

	Functions FG(X) and FKER(X,T) must be supplied by the user.
		FG is the right hand side function g(x) and FKER(X,T) is 
		the kernel. The integral equation is specified in the form
		given by Eq.(12.1)

	Required functions : GAUELM, INVIT, FG, FKER
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int invit(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit);

int fred(int m, double a, double b, double wt[], double x[], double f[],
	double fc[], double fg(double ), double fker(double , double ),
	double *ei, int iq, int it, double reps)

{
	int i,j,n,no,i1,ier,iflg,nit;
	double h,a1,det,p,rc,u2,u2m,t1,t2;
	int *iwk;
	double *wk;

/*	Weights and abscissas for Gauss-Legendre quadrature.
	For N-point formula W[K]=W[N-K-1] and X[K]=-X[N-K-1]
		For K=0,1,...,N/2-1. Hence only half points are tabulated.
	For 2-point W[0]; 4-point W[1], W[2]; 8-point W[3],...,W[6];
	16-point W[7],...,W[14]; 32-point W[15],...,W[30] are the
	weights corresponding to abscissas X[I].
*/

	double w[31]={1.0e0,
      0.34785484513745385737e0, 0.65214515486254614263e0,
      0.10122853629037625915e0, 0.22238103445337447054e0,
      0.31370664587788728734e0, 0.36268378337836198297e0,
      0.02715245941175409485e0, 0.06225352393864789286e0,
      0.09515851168249278481e0, 0.12462897125553387205e0,
      0.14959598881657673208e0, 0.16915651939500253819e0,
      0.18260341504492358887e0, 0.18945061045506849629e0,
      0.00701861000947009660e0, 0.01627439473090567061e0,
      0.02539206530926205945e0, 0.03427386291302143310e0,
      0.04283589802222668066e0, 0.05099805926237617620e0,
      0.05868409347853554714e0, 0.06582222277636184684e0,
      0.07234579410884850623e0, 0.07819389578707030647e0,
      0.08331192422694675522e0, 0.08765209300440381114e0,
      0.09117387869576388471e0, 0.09384439908080456564e0,
      0.09563872007927485942e0, 0.09654008851472780057e0};

	double xa[31] = {0.57735026918962576451e0,
      0.86113631159405257522e0, 0.33998104358485626480e0,
      0.96028985649753623168e0, 0.79666647741362673959e0,
      0.52553240991632898582e0, 0.18343464249564980494e0,
      0.98940093499164993260e0, 0.94457502307323257608e0,
      0.86563120238783174388e0, 0.75540440835500303390e0,
      0.61787624440264374845e0, 0.45801677765722738634e0,
      0.28160355077925891323e0, 0.09501250983763744019e0,
      0.99726386184948156354e0, 0.98561151154526833540e0,
      0.96476225558750643077e0, 0.93490607593773968917e0,
      0.89632115576605212397e0, 0.84936761373256997013e0,
      0.79448379596794240696e0, 0.73218211874028968039e0,
      0.66304426693021520098e0, 0.58771575724076232904e0,
      0.50689990893222939002e0, 0.42135127613063534536e0,
      0.33186860228212764978e0, 0.23928736225213707454e0,
      0.14447196158279649349e0, 0.04830766568773831623e0};

	if(m<3 || it>3 || it<=0) return 706;
	n=m;
	ier=0;
	if(iq==1) {
/*	Use the trapezoidal rule */
		h=(b-a)/(n-1);
		wt[0]=0.5*h; wt[n-1]=wt[0];
		x[0]=a; x[n-1]=b;
		for(i=1; i<n-1; ++i) {x[i]=a+i*h; wt[i]=h;}
	}

	else if(iq==2) {
/*	Use the Simpson's 1/3 rule, if N is even, then reduce it by 1 */
		n=2*((n-1)/2)+1;
		h=(b-a)/(n-1);
		wt[0]=h/3.0; x[0]=a;
		for(i=1; i<n-1; i += 2) {
			x[i]=a+i*h; x[i+1]=a+(i+1)*h;
			wt[i]=4.0*h/3.0; wt[i+1]=2.0*h/3.0;
		}
		wt[n-1]=wt[0]; x[n-1]=b;
	}

	else if(iq>=0) {
/*	Try Gauss-Legendre formulas */
		if(iq==4) no=1;
		else if(iq==8) no=3;
		else if(iq==16) no=7;
		else if(iq==32) no=15;
		else return 707;	 /*	If IQ is not acceptable then quit */
		n=(n/iq)*iq;
		if(n<iq) return 707;

/*	Setup the weights and abscissas for Gauss-Legendre formula */
		h=iq*(b-a)/n;
		for(i=0; i<n; i += iq) {
			a1=a+i*h/iq+h/2.0;
			for(i1=0; i1<iq/2; ++i1) {
				wt[i+i1]=w[no+i1]*h/2.0;
				wt[i+iq-i1-1]=wt[i+i1];
				x[i+i1]=a1-xa[no+i1]*h/2.0;
				x[i+iq-i1-1]=a1+xa[no+i1]*h/2.0;
			}
		}
		if(m != n) ier=-11;
	}

/*	Setting up the equation matrix and the right hand side */
	wk=(double *) calloc((size_t) (n*n), sizeof(double));
	iwk=(int *) calloc((size_t) n, sizeof(int));
	for(i=0; i<n; ++i) {
		f[i]=fg(x[i]);
		for(j=0; j<n; ++j) wk[i+j*n]=wt[i]*fker(x[j],x[i]);
		if(it==2) wk[i+i*n]=wk[i+i*n]-1.0;
	}

	if(it<=2) {
/*	Solve the system of linear equations */
		iflg=0;
		ier=gauelm(n,1,wk,f,&det,iwk,n,&iflg);
		if(ier>0) {free(iwk); free(wk); return ier;}
	}
	else {
/*	Solve the eigenvalue problem */
		iflg=0;
		p=*ei;
		nit=0;
		ier=invit(wk,n,n,&p,f,iflg,ei,&rc,reps,&nit);
	}

	if(iq==1 && it!=3) {
/*	Apply the deferred correction */
		for(i=0; i<n; ++i) {
			u2=fker(x[i],x[1])*f[1];
			u2m=fker(x[i],x[n-2])*f[n-2];
			t1=fker(x[i],x[0])*f[0]-u2+fker(x[i],x[n-1])*f[n-1]-u2m;
			t2=t1+fker(x[i],x[2])*f[2]-u2+fker(x[i],x[n-3])*f[n-3]-u2m;
			fc[i]=h*t1/12.0+h*t2/24.0;
		}
		ier=gauelm(n,1,wk,fc,&det,iwk,n,&iflg);

		for(i=0; i<n; ++i) fc[i]=f[i]+fc[i];
	}
	free(iwk); free(wk);
	return ier;
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


/*	Real eigenvalues and eigenvector of a real matrix using inverse iteration

	A : (input) Array of length IA*M containing the matrix elements
	M : (input) Order of the matrix
	IA : (input) The second dimension of A as declared in the calling function, IA>=M
	P : (input/output) Initial value of the shift. This will be modified
		by the function if IFLG>0
	U : (input/output) Array of length M, which should specify the
		initial approximation to eigenvector. After execution it
		will contain the calculated eigenvector.
	IFLG : (input) Integer variable to specify the type of iteration required
		If IFLG=0 the shift P is kept fixed
		If IFLG=1 the shift P is varied using Rayleigh quotient
		If IFLG=2 the shift P is varied using max(V_s+1)
	EI : (output) Estimated eigenvalue using simple inverse iteration
	ERC : (output) Estimated eigenvalue using Rayleigh quotient
	REPS : (input) Required absolute accuracy. Iteration is terminated
		when all components of eigenvector and the eigenvalue have
		converged to REPS.
	NIT : (input/output) Number of iterations required. If it is
		zero or negative NIT is set to NIT0=100
		
	Error status is returned by the value of the function INVIT.
		0 value implies successful execution
		106 implies that M<=1 or M>IA, in which case no
			calculations are done
		141 implies that vector is zero at some stage and
			calculations are aborted
		142 implies that inverse iteration has failed to converge

	Required functions : GAUELM
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int invit(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit)

{
	int i,j,k,ifl,km,ier, nit0=100;
	double epi,r1,ukm,s1,s2,det;
	int *iwk;
	double *wk,*v;

	if(m<=1 || m> ia) return 106;

	iwk=(int *) calloc((size_t) m, sizeof(int));
	v=(double *) calloc((size_t) m, sizeof(double));
	wk=(double *) calloc((size_t) (m*m), sizeof(double));

/*	Copy the matrix to WK and apply the shift */
	for(i=0; i<m; ++i) {
		v[i]=u[i];
		for(j=0; j<m; ++j) wk[i+j*m]=a[i+j*m];
		wk[i+i*m]=a[i+i*m]-(*p);
	}

	ifl=1;
	if(*nit<=0) *nit=nit0;
/*	Perform Gaussian elimination on A-pI */
	ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
	if(ier>0) {free(wk); free(v); free(iwk); return ier;}

	epi=0.0;
/*	Loop for inverse iteration */
	for(j=1; j<=(*nit); ++j) {
		ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
		if(ier>0) {free(wk); free(v); free(iwk); return ier;}

/*	Normalising the vector U */
		r1=0.0; km=0;
		for(k=0; k<m; ++k) {
			if(r1<fabs(u[k])) {r1=fabs(u[k]); km=k;}
		}
		ukm=u[km];
/*	If the vector is zero, then quit */
		if(ukm==0.0) {free(wk); free(v); free(iwk); return 141;}

/*	The eigenvalue */
		*ei=v[km]/ukm+(*p);
		s1=0.0; s2=0.0;
/*	Calculating the Rayleigh quotient */
		for(k=0; k<m; ++k) {
			s1=s1+u[k]*v[k];
			s2=s2+fabs(u[k])*fabs(u[k]);
			u[k]=u[k]/ukm;
		}
		*erc=(*p)+s1/s2;

/*	Convergence check */
		r1=fabs(*ei-epi);
		for(i=0; i<m; ++i) {
			if(r1<fabs(v[i]-u[i])) r1=fabs(v[i]-u[i]);
			v[i]=u[i];
		}
		if(fabs(r1)<reps) {free(wk); free(v); free(iwk); return 0;}
		epi=(*ei);

		if(iflg>=1) {
/*	Update the shift */
			*p=(*erc);
			if(iflg==2) *p=(*ei);
/*	Setting up the new matrix A-pI */
			for(i=0; i<m; ++i) {
				for(k=0; k<m; ++k) wk[i+k*m]=a[i+k*m];
				wk[i+i*m]=a[i+i*m]-(*p);
			}
			ifl=1;
			ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
			if(ier>0) {free(wk); free(v); free(iwk); return ier;}
		}
	}

	free(wk); free(v); free(iwk);
	return 142;
}
