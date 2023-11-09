/*	Real eigenvalue and eigenvector of a general real matrix using inverse iteration */

#include <stdio.h>
#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int invit(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit);
int invit_l(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit);

main()
{
	int i,i1,j,n,m, id, iflg, ier,np,nmax;
	double hh, x[100], p,ei,erc,reps,p0;
/*	double a[4][4]={4,-5,0,3, 0,4,-3,-5, 5,-3,4,0, 3,0,5,4}; */
	double a[4][4]={-2,2,2,2, -3,3,2,2, -2,0,4,2, -1,0,0,5};

	id=4; m=n=4; np=0; reps=1.e-14;
	printf(" The matrix is :\n");
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) printf(" %e ",a[i][j]);
		printf("\n");
	}

	for(i1=0; i1<99; ++i1) {
		printf("type iflg=0/1/2,  p=initial shift    (quits when iflg<0)\n");
		scanf(" %d %le",  &iflg, &p);
		if(iflg<0) return 0;

		p0=p;
/*	Initial guess for eigenvector */
		for(i=0; i<m; ++i) x[i]=1.0;
		i=invit(&a[0][0],m,id,&p,x,iflg,&ei,&erc,reps,&np);
		printf(" ier = %d    m =  %d    initial shift = %e \n",i,m,p0);
		printf(" final shift = %e    eigenvalue = %e    Rayleigh Quotient = %e\n",p,ei,erc);
		printf(" eigenvector : ");
		for(i=0; i<m; ++i) printf(" %e ",x[i]);
		printf(" \n");

/*	Calculate the left eigenvector for the same matrix */
		p=p0;
		for(i=0; i<m; ++i) x[i]=1.0;
		i=invit_l(&a[0][0],m,id,&p,x,iflg,&ei,&erc,reps,&np);
		printf(" ier = %d   final shift = %e    eigenvalue = %e    Rayleigh Quotient = %e\n",i,p,ei,erc);
		printf(" left eigenvector : ");
		for(i=0; i<m; ++i) printf(" %e ",x[i]);
		printf(" \n");
	}
	return;
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



/*	Real eigenvalues and left eigenvector of a real matrix using inverse iteration

	A : (input) Array of length IA*M containing the matrix elements
	M : (input) Order of the matrix
	IA : (input) The second dimension of A as declared in the calling function, IA>=M
	P : (input/output) Initial value of the shift. This will be modified
		by the function if IFLG>0
	U : (input/output) Array of length M, which should specify the
		initial approximation to eigenvector. After execution it
		will contain the calculated left eigenvector.
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
		
	Error status is returned by the value of the function INVIT_L.
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

int invit_l(double *a, int m, int ia, double *p, double u[], int iflg,
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
		for(j=0; j<m; ++j) wk[j+i*m]=a[i+j*m];
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
				for(k=0; k<m; ++k) wk[k+i*m]=a[i+k*m];
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
